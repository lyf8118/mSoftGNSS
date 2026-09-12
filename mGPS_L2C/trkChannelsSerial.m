function [trackResults, channel]= trkChannelsSerial(fid, channel, settings)
% Performs GPS L2C code and carrier tracking in channel-serial mode.
%
%[trackResults, channel] = trkChannelsSerial(fid, channel, settings)
%
% This implementation processes one channel at a time. For every active
% PRN it reads one integration interval from the IF file, calls the selected
% correlator backend, and updates the carrier PLL and code DLL before
% moving to the next interval. The serial form is useful as the reference
% implementation for the SIMD/GPU accelerated correlators.
%
%   Inputs:
%       fid             - File identifier of the IF signal record.
%       channel         - Channel state prepared from acquisition results.
%                       Each active element contains PRN, acquired carrier
%                       frequency, CM code phase, CL code phase and status.
%       settings        - Receiver settings, including sampling frequency,
%                       integration time, loop bandwidths, correlator type,
%                       file format and L2C code definitions.
%   Outputs:
%       trackResults    - Tracking results for each channel. The structure
%                       stores absolute sample positions, carrier/code
%                       frequencies, E/P/L correlator outputs, residual
%                       phases, DLL/PLL discriminator values, C/No and PLL
%                       lock detector values once per integration interval.
%       channel         - Input channel structure with the final tracking
%                       status copied into the corresponding result.

%--------------------------------------------------------------------------
%                         CU Multi-GNSS SDR
% (C) Developed for GPS L2C SDR by Yafeng Li, Nagaraj C. Shivaramaiah
% and Dennis M. Akos.
% Based on the original framework for GPS C/A SDR by Darius Plausinaitis,
% Peter Rinder, Nicolaj Bertelsen and Dennis M. Akos
%--------------------------------------------------------------------------

%This program is free software; you can redistribute it and/or
%modify it under the terms of the GNU General Public License
%as published by the Free Software Foundation; either version 2
%of the License, or (at your option) any later version.
%
%This program is distributed in the hope that it will be useful,
%but WITHOUT ANY WARRANTY; without even the implied warranty of
%MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%GNU General Public License for more details.
%
%You should have received a copy of the GNU General Public License
%along with this program; if not, write to the Free Software
%Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
%USA.
%--------------------------------------------------------------------------

%CVS record:
%$Id: tracking.m,v 1.14.2.31 2006/08/14 11:38:22 dpl Exp $

%% Initialize result structure ============================================
% Default channel status. It is replaced after successful tracking.
trackResults.status         = '-';      % No tracked signal, or lost lock
% Number of coherent integrations to process.
NumToProcess =  round(settings.msToProcess/1000/settings.intTime);
% Absolute sample index of each tracked L2C code epoch in the IF record.
trackResults.absoluteSample = zeros(1, NumToProcess);
% Code NCO frequency used for each integration interval.
trackResults.codeFreq       = inf(1, NumToProcess);
% Carrier NCO frequency used for each integration interval.
trackResults.carrFreq       = inf(1, NumToProcess);
% In-phase Early/Prompt/Late correlator outputs from the data branch.
trackResults.I_P            = zeros(1, NumToProcess);
trackResults.I_E            = zeros(1, NumToProcess);
trackResults.I_L            = zeros(1, NumToProcess);
% Quadrature Early/Prompt/Late correlator outputs from the data branch.
trackResults.Q_E            = zeros(1, NumToProcess);
trackResults.Q_P            = zeros(1, NumToProcess);
trackResults.Q_L            = zeros(1, NumToProcess);
if settings.pilotTRKflag
    % Prompt correlator outputs from the pilot branch.
    trackResults.Pilot_I_P            = zeros(1, NumToProcess);
    trackResults.Pilot_Q_P            = zeros(1, NumToProcess);
end

% Raw and filtered discriminator outputs from DLL and PLL.
trackResults.dllDiscr       = inf(1, NumToProcess);
trackResults.dllDiscrFilt   = inf(1, NumToProcess);
trackResults.pllDiscr       = inf(1, NumToProcess);
trackResults.pllDiscrFilt   = inf(1, NumToProcess);
% Residual code/carrier phases saved for later navigation processing.
trackResults.remCodePhase   = inf(1, NumToProcess);
trackResults.remCarrPhase   = inf(1, NumToProcess);

% C/No and PLL lock detector histories for the data branch.
trackResults.DataCNo  = zeros(1,floor(NumToProcess/settings.CNoInterval));
trackResults.DataPLD  = zeros(1,floor(NumToProcess/settings.CNoInterval));

% C/No and PLL lock detector histories for the pilot and combined L2C
% branches. These are updated once per CNoInterval integrations.
if (settings.pilotTRKflag == 1)
    trackResults.PilotCNo = zeros(1,floor(NumToProcess/settings.CNoInterval));
    trackResults.PilotPLD  = zeros(1,floor(NumToProcess/settings.CNoInterval));
    trackResults.L2C_CNo  = zeros(1,floor(NumToProcess/settings.CNoInterval));
end

%--- Allocate one result structure per configured channel ------------------
trackResults = repmat(trackResults, 1, settings.numberOfChannels);

%% Initialize tracking variables ==========================================
%--- DLL variables --------------------------------------------------------
% Code length of one coherent integration in equivalent L2C chips.
codeLength = settings.codeLength;
% Coherent integration time used by the code loop filter.
PDIcode = settings.intTime;
% Calculate second-order DLL filter coefficients.
[tau1code, tau2code] = calcLoopCoef(settings.dllNoiseBandwidth, ...
    settings.dllDampingRatio, 1.0);

%--- PLL variables --------------------------------------------------------
% Coherent integration time used by the carrier loop filter.
PDIcarr = settings.intTime;
% Calculate second-order PLL filter coefficients. For L2C the loop is
% updated every 20 ms, so use unit loop gain; the 1 ms QPSK receivers'
% k = 0.25 setting makes the discrete carrier loop unstable here.
[tau1carr, tau2carr] = calcLoopCoef(settings.pllNoiseBandwidth, ...
    settings.pllDampingRatio, 1.0);
% -------- Number of acquired signals --------------------------------------
% Count active channels for GUI progress reporting.
TrackedNr = nnz([channel.status]== 'T');

% Start waitbar.
hwb = waitbar(0,'Tracking...');

% Adjust the waitbar height to insert C/No and PRN status text.
CNoPos = get(hwb,'Position');
set(hwb,'Position',[CNoPos(1),CNoPos(2),CNoPos(3),90],'Visible','on');

% Real sample files consume one value per sample; interleaved I/Q files
% consume two values per complex sample.
if (settings.fileType==1)
    dataAdaptCoeff=1;
else
    dataAdaptCoeff=2;
end

if settings.correlatorType == 0
    % MATLAB correlator processes double-precision samples.
    dataConverter = strcat(settings.dataType,'=>double');
elseif (settings.correlatorType == 1 || settings.correlatorType == 2)
    % SIMD and GPU correlator processes int16 samples.
    dataConverter = strcat(settings.dataType,'=>int16');
end
% The GUI bar update period depends on the correlator backend.
if (settings.correlatorType == 0)
    barUprate = 3;
else
    barUprate = 100;
end
%% Start processing channels ==============================================
for channelNr = 1:settings.numberOfChannels

    % Only process a channel when acquisition assigned a PRN.
    if (channel(channelNr).PRN ~= 0)
        % Save the PRN being tracked by this channel.
        trackResults(channelNr).PRN     = channel(channelNr).PRN;

        % Seek to the acquisition code phase for this PRN, after applying
        % any user configured file offset. int16 samples use two bytes per
        % real value, and interleaved I/Q files are scaled by dataAdaptCoeff.
        if strcmp(settings.dataType,'int16')
            fseek(fid, dataAdaptCoeff*2*(settings.skipNumberOfSamples + ...
                channel(channelNr).codePhase-1), 'bof');
        else
            fseek(fid, dataAdaptCoeff*(settings.skipNumberOfSamples + ...
                channel(channelNr).codePhase-1), 'bof');
        end

        % Generate the CM data code sampled at one value per equivalent chip.
        cmCode = generateCMcode(channel(channelNr).PRN,settings);
        % Repeat CM to the same length as CL. This keeps the L2C data and
        % pilot code tables compatible with the common QPSK MEX interface,
        % so the same correlator code can process both branches.
        cmCode = repmat(cmCode,1,75);
        % Add one chip of wraparound on both sides for early/late indexing.
        cmCode = [cmCode(settings.CLCodeLength) cmCode cmCode(1)];  %#ok<AGROW>

        % CL code phase corresponding to the acquired CM code phase.
        CLCodePhase = channel(channelNr).CLCodePhase;
        % Generate the CL pilot code sampled at one value per equivalent chip.
        CLCode = generateCLcode(channel(channelNr).PRN,settings);
        % Add one chip of wraparound on both sides for early/late indexing.
        CLCode = [CLCode(settings.CLCodeLength) CLCode CLCode(1)]; %#ok<AGROW>

        % Place equal-length CM and CL tables contiguously. The correlator
        % selects the data or pilot half using the L2C code phase passed
        % below, matching the shared QPSK MEX input layout.
        L2CCodeTable = [cmCode CLCode];
        if settings.correlatorType == 1
            % SIMD serial correlator expects an int32 local code table.
            L2CCodeTable = int32(L2CCodeTable);
        elseif settings.correlatorType == 2
            % GPU serial correlator expects an int8 local code table.
            L2CCodeTable = int8(L2CCodeTable);
        end

        %--- Perform various initializations ------------------------------
        % Initial code NCO frequency.
        codeFreq      = channel(channelNr).codeFreq;
        codeFreqBasis = channel(channelNr).codeFreq;
        % Residual code phase carried from one integration to the next.
        remCodePhase  = 0.0;
        % Carrier NCO starts from the acquisition Doppler estimate.
        carrFreq      = channel(channelNr).acquiredFreq;
        % Keep the acquisition estimate as the fixed carrier frequency basis.
        carrFreqBasis = channel(channelNr).acquiredFreq;
        % Residual carrier phase carried from one integration to the next.
        remCarrPhase  = 0.0;
        % Previous DLL state for the recursive loop filter.
        oldCodeNco   = 0.0;
        oldCodeError = 0.0;
        % Previous PLL state for the second-order carrier loop filter.
        oldCarrNco   = 0.0;
        oldCarrError = 0.0;
        % C/No values: [data, pilot, combined L2C].
        CNoValue = zeros(1,3);
        tempCNoValue = zeros(1,3);

        %=== Process the configured number of coherent integrations ========
        for loopCnt =  1:NumToProcess

            %% GUI update -------------------------------------------------------------
            % Update the GUI periodically so MATLAB remains responsive
            % without spending too much time repainting the waitbar.
            if (rem(loopCnt, barUprate) == 0)

                Ln=newline;
                trackingStatus=['Tracking: Ch ', int2str(channelNr), ...
                    ' of ', int2str(TrackedNr),Ln ...
                    'PRN: ', int2str(channel(channelNr).PRN),Ln ...
                    'Completed ',int2str(loopCnt*20), ...
                    ' of ', int2str(NumToProcess*20), ' msec',Ln...
                    'Data C/No: ',int2str(CNoValue(1)),' (dB-Hz);',...
                    '   Pilot C/No: ',int2str(CNoValue(2)),' (dB-Hz)'];
                try
                    waitbar(loopCnt/NumToProcess, hwb, trackingStatus);
                catch %#ok<CTCH>
                    % The progress bar was closed. It is used as a signal
                    % to stop, "cancel" processing. Exit.
                    disp('Progress bar closed, exiting...');
                    return
                end
            end

            %% Read next block of data ------------------------------------------------
            % Record the absolute sample index of the current code epoch.
            if strcmp(settings.dataType,'int16')
                trackResults(channelNr).absoluteSample(loopCnt) = ...
                    (ftell(fid))/dataAdaptCoeff/2;
            else
                trackResults(channelNr).absoluteSample(loopCnt) = ...
                    (ftell(fid))/dataAdaptCoeff;
            end
            % Code phase step in chips/sample from the current code NCO.
            codePhaseStep = codeFreq / settings.samplingFreq;
            % Number of whole IF samples needed to complete this integration.
            blksize = ceil((codeLength-remCodePhase) / codePhaseStep);

            % Read the exact sample block needed for this integration.
            [rawSignal, samplesRead] = fread(fid, dataAdaptCoeff * blksize, dataConverter);

            % Stop tracking cleanly if the file does not contain enough data.
            if (samplesRead ~= dataAdaptCoeff*blksize)
                disp('Not able to read the specified number of samples  for tracking, exiting!')
                delete(hwb);
                return
            end

            %% Correlator implementation  --------------------------------------
            % Carrier phase step in rad/sample based on the current
            % carrier NCO and the sampling frequency.
            carrPhaseStep = carrFreq * 2.0 * pi / settings.samplingFreq;
            % Convert the residual phase within the current CM interval into
            % an absolute phase within the 75-period CL sequence.
            remCodePhaseL2C = remCodePhase + codeLength *(CLCodePhase-1);
            if settings.correlatorType == 0 % Using MATLAB correlator
                correValues = corrMatlabSerialTMBPSK(settings, rawSignal, ...
                    L2CCodeTable, remCarrPhase, carrPhaseStep, remCodePhaseL2C,...
                    codePhaseStep);
            elseif settings.correlatorType == 1 % Using SIMD of CPU
                % Compile the MEX file with: mex corrSIMDSerialQPSK.cpp
                correValues = corrSIMDSerialQPSK(settings,rawSignal, ...
                    L2CCodeTable, remCarrPhase,carrPhaseStep,remCodePhaseL2C, ...
                    codePhaseStep);
            elseif settings.correlatorType == 2 % Using GPU
                % Compile the MEX file with: mexcuda corrGPUSerialQPSK.cu
                correValues = corrGPUSerialQPSK(settings,rawSignal, ...
                    L2CCodeTable, remCarrPhase,carrPhaseStep,remCodePhaseL2C, ...
                    codePhaseStep, channel(channelNr).PRN);
            end

            % Unpack data-branch Early/Prompt/Late correlator outputs.
            I_E = correValues(1); Q_E = correValues(2);
            I_P = correValues(3); Q_P = correValues(4);
            I_L = correValues(5); Q_L = correValues(6);

            if (settings.pilotTRKflag == 1)
                % Unpack pilot-branch Early/Prompt/Late correlator outputs.
                I_ECL = correValues(7); Q_ECL = correValues(8);
                I_PCL = correValues(9); Q_PCL = correValues(10);
                I_LCL = correValues(11); Q_LCL = correValues(12);
            end

            % Save residual phases used by the current correlation.
            trackResults(channelNr).remCodePhase(loopCnt) = remCodePhase;
            trackResults(channelNr).remCarrPhase(loopCnt) = remCarrPhase;
            % Advance residual code/carrier phase for the next integration.
            remCodePhase = rem(blksize*codePhaseStep + remCodePhaseL2C,codeLength);
            remCarrPhase = rem(carrPhaseStep * blksize + remCarrPhase, 2 * pi);

            %% Find PLL error and update carrier NCO ----------------------

            % Carrier phase discriminator from the data prompt correlator.
            carrError = atan(Q_P / I_P) / (2.0 * pi);

            % Combine carrier phase error estimates from data and pilot
            % branches when pilot tracking is enabled.
            if (settings.pilotTRKflag == 1)
                % atan is not affected by the L2C pilot overlay modulation.
                carrErrorCL = atan(Q_PCL/I_PCL)/ (2.0 * pi);

                % Data and pilot branches are treated with equal weight.
                carrError = (carrError + carrErrorCL)/2;
            end

            % Second-order carrier loop filter and NCO command.
            carrNco = oldCarrNco + (tau2carr/tau1carr) * ...
                (carrError - oldCarrError) + carrError * (PDIcarr/tau1carr);
            oldCarrNco   = carrNco;
            oldCarrError = carrError;

            % Save carrier frequency for current correlation
            trackResults(channelNr).carrFreq(loopCnt) = carrFreq;

            % Apply PLL correction around the acquisition carrier frequency.
            carrFreq = carrFreqBasis + carrNco;

            %% Find DLL error and update code NCO -------------------------------------
            % Non-coherent early-minus-late envelope discriminator for DLL.
            codeError = (sqrt(I_E^2 + Q_E^2) - sqrt(I_L^2 + Q_L^2)) / ...
                (sqrt(I_E^2 + Q_E^2) + sqrt(I_L^2 + Q_L^2));

            % Combine code tracking errors from data and pilot branches.
            if (settings.pilotTRKflag == 1)
                codeErrorCL = (sqrt(I_ECL^2 + Q_ECL^2) - sqrt(I_LCL^2 + Q_LCL^2)) / ...
                    (sqrt(I_ECL^2 + Q_ECL^2) + sqrt(I_LCL^2 + Q_LCL^2));
                codeError = (codeError + codeErrorCL)/2;
            end

            % Advance the CL period index and wrap after the 75th CM period.
            CLCodePhase = CLCodePhase + 1;
            if (CLCodePhase >= 76)
                CLCodePhase = 1;
            end

            % Code loop filter and NCO command.
            codeNco = oldCodeNco + (tau2code/tau1code) * ...
                (codeError - oldCodeError) + codeError * (PDIcode/tau1code);
            oldCodeNco   = codeNco;
            oldCodeError = codeError;

            % Save code frequency for current correlation
            trackResults(channelNr).codeFreq(loopCnt) = codeFreq;

            % Apply DLL correction around the acquired L2C code frequency.
            codeFreq = codeFreqBasis - codeNco;

            %% Record values for postprocessing and diagnostics -----------------------
            trackResults(channelNr).dllDiscr(loopCnt)       = codeError;
            trackResults(channelNr).dllDiscrFilt(loopCnt)   = codeNco;
            trackResults(channelNr).pllDiscr(loopCnt)       = carrError;
            trackResults(channelNr).pllDiscrFilt(loopCnt)   = carrNco;

            % Data-branch E/P/L correlator outputs.
            trackResults(channelNr).I_E(loopCnt) = I_E;
            trackResults(channelNr).I_P(loopCnt) = I_P;
            trackResults(channelNr).I_L(loopCnt) = I_L;
            trackResults(channelNr).Q_E(loopCnt) = Q_E;
            trackResults(channelNr).Q_P(loopCnt) = Q_P;
            trackResults(channelNr).Q_L(loopCnt) = Q_L;

            % Pilot-branch prompt correlator outputs.
            if settings.pilotTRKflag
                trackResults(channelNr).Pilot_I_P(loopCnt) = I_PCL;
                trackResults(channelNr).Pilot_Q_P(loopCnt) = Q_PCL;
            end
            %% CNo calculation --------------------------------------------------------

            if (rem(loopCnt,settings.CNoInterval)==0)
                % Compute C/No and PLL lock detector outputs over the latest
                % CNoInterval window.
                [CNoValue, PllDetector]= ...
                    Calc_CNo_PLD(trackResults(channelNr),settings,loopCnt);

                CNoCnt = loopCnt/settings.CNoInterval;

                % Save C/No for data branch. A 0.5/0.5 smoother is used to
                % smooth the results
                trackResults(channelNr).DataCNo(CNoCnt) = ...
                    CNoValue(1) * 0.5 + tempCNoValue(1) * 0.5;
                % Save PLL lock detector output for data channel
                trackResults(channelNr).DataPLD(CNoCnt) = PllDetector(1);

                % Save C/No and PLL lock detector output for the pilot and
                % combined L2C branches.
                if (settings.pilotTRKflag == 1)
                    trackResults(channelNr).PilotCNo(CNoCnt) = ...
                        CNoValue(2) * 0.5 + tempCNoValue(2) * 0.5;
                    trackResults(channelNr).L2C_CNo(CNoCnt) = ...
                        CNoValue(3) * 0.5 + tempCNoValue(3) * 0.5;
                    trackResults(channelNr).PilotPLD(CNoCnt) = PllDetector(2);
                end
            end
            tempCNoValue = CNoValue;
        end % for loopCnt

        % If tracking reached the end of the requested interval, copy the
        % channel status. A future lock detector can replace this decision.
        trackResults(channelNr).status  = channel(channelNr).status;

    end % if a PRN is assigned
end % for channelNr

% Close the waitbar
close(hwb)
clear corrSIMDSerialQPSK corrGPUSerialQPSK;
