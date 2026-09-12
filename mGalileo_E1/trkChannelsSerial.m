function [trackResults, channel]= trkChannelsSerial(fid, channel, settings)
% Performs serial Galileo E1 code and carrier tracking for all enabled
% channels. Each channel is processed independently with the selected
% MATLAB, SIMD, or GPU correlator backend.
%
% [trackResults, channel] = trkChannelsSerial(fid, channel, settings)
%
%   Inputs:
%       fid             - File identifier of the signal record.
%       channel         - PRN, carrier frequencies and code phases of all
%                         satellites to be tracked, normally prepared from
%                         acquisition results.
%       settings        - Receiver settings.
%   Outputs:
%       trackResults    - Tracking results, including correlator outputs,
%                         loop discriminator values, code/carrier NCO
%                         states, residual phases, C/No estimates, and
%                         absolute sample positions for each integration.

%--------------------------------------------------------------------------
%                         CU Multi-GNSS SDR
% (C) Developed for Galileo E1 SDR by Yafeng Li, Nagaraj C. Shivaramaiah
% and Dennis M. Akos.
% Based on the original framework for GPS C/A by Darius Plausinaitis,
% Peter Rinder, Nicolaj Bertelsen and Dennis M. Akos
%
% Reference: Adapted within the CU Multi-GNSS SDR receiver framework for Galileo E1.
% implementation of an open-source Galileo E1 SDR receiver.
%
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

%% Initialize result structure ============================================

% Channel status: '-' means no tracked signal or lost lock.
trackResults.status         = '-';      % No tracked signal, or lost lock
% Number of coherent integration intervals to process.
NumToProcess =  round(settings.msToProcess/1000/settings.intTime);
% Absolute sample index at the start of each tracked E1 code epoch.
trackResults.absoluteSample = zeros(1, NumToProcess);
% Code NCO frequency for each integration interval.
trackResults.codeFreq       = inf(1, NumToProcess);
% Carrier NCO frequency for each integration interval.
trackResults.carrFreq       = inf(1, NumToProcess);
% Data-branch early, prompt, and late in-phase correlator outputs.
trackResults.I_P            = zeros(1, NumToProcess);
trackResults.I_E            = zeros(1, NumToProcess);
trackResults.I_L            = zeros(1, NumToProcess);
% Data-branch early, prompt, and late quadrature correlator outputs.
trackResults.Q_E            = zeros(1, NumToProcess);
trackResults.Q_P            = zeros(1, NumToProcess);
trackResults.Q_L            = zeros(1, NumToProcess);
% Pilot-branch prompt correlator outputs, saved when pilot tracking is on.
if (settings.pilotTRKflag == 1)
    trackResults.Pilot_I_P  = zeros(1, NumToProcess);
    trackResults.Pilot_Q_P  = zeros(1, NumToProcess);
end
% Loop discriminators
trackResults.dllDiscr       = inf(1, NumToProcess);
trackResults.dllDiscrFilt   = inf(1, NumToProcess);
trackResults.pllDiscr       = inf(1, NumToProcess);
trackResults.pllDiscrFilt   = inf(1, NumToProcess);
% Residual code and carrier phase carried into each tracking update.
trackResults.remCodePhase   = inf(1, NumToProcess);
trackResults.remCarrPhase   = inf(1, NumToProcess);
% Data-branch C/No and PLL lock detector values.
trackResults.DataCNo  = zeros(1,floor(NumToProcess/settings.CNoInterval));
trackResults.DataPLD  = zeros(1,floor(NumToProcess/settings.CNoInterval));
% Pilot and combined L1C C/No values, plus pilot PLL lock detector values.
if (settings.pilotTRKflag == 1)
    trackResults.PilotCNo = zeros(1,floor(NumToProcess/settings.CNoInterval));
    trackResults.PilotPLD  = zeros(1,floor(NumToProcess/settings.CNoInterval));
    trackResults.E1_CNo  = zeros(1,floor(NumToProcess/settings.CNoInterval));
end
%--- Copy initial settings for all channels -------------------------------
trackResults = repmat(trackResults, 1, settings.numberOfChannels);

%% Initialize tracking variables ==========================================
% Primary E1 code length in chips. The local code table below is expanded
% to the 12-sample-per-chip grid used by the TMBOC/QPSK correlator.
codeLength = settings.codeLength;
% Early-late spacing in primary-code chips, used by the DLL discriminator.
earlyLateSpc = settings.dllCorrelatorSpacing;
% Convert the early-late spacing to the 12x local-code grid expected by the
% shared QPSK correlator interface.
settings.dllCorrelatorSpacing = settings.dllCorrelatorSpacing*12;
% Coherent integration time for the DLL.
PDIcode = settings.intTime;

% DLL loop filter coefficients.
[tau1code, tau2code] = calcLoopCoef(settings.dllNoiseBandwidth, ...
    settings.dllDampingRatio, 1.0);

% Coherent integration time for the carrier loop.
PDIcarr = settings.intTime;
% Second-order PLL loop filter coefficients.
[tau1carr, tau2carr] = calcLoopCoef(settings.pllNoiseBandwidth, ...
    settings.pllDampingRatio, 1.0);

% Number of acquired signals marked for tracking.
TrackedNr = nnz([channel.status]== 'T');

% Start waitbar
hwb = waitbar(0,'Tracking...');

% Adjust the waitbar size so the status text can show PRN and C/No values.
CNoPos=get(hwb,'Position');
set(hwb,'Position',[CNoPos(1),CNoPos(2),CNoPos(3),90],'Visible','on');

if (settings.fileType==1)
    dataAdaptCoeff=1;
else
    dataAdaptCoeff=2;
end

if settings.correlatorType == 0
    % MATLAB correlator processes double-precision samples.
    dataConverter = strcat(settings.dataType,'=>double');
elseif (settings.correlatorType == 1 || settings.correlatorType == 2)
    % SIMD and GPU correlators process int16 samples.
    dataConverter = strcat(settings.dataType,'=>int16');
end
% The GUI bar update period depends on the correlator backend.
if (settings.correlatorType == 0)
    barUprate = 10;
else
    barUprate = 200;
end
%% Start processing channels ==============================================
for channelNr = 1:settings.numberOfChannels

    % Only process channels with a valid PRN from acquisition.
    if (channel(channelNr).PRN ~= 0)
        % Save the PRN assigned to this tracking channel.
        trackResults(channelNr).PRN     = channel(channelNr).PRN;

        % Seek to the acquisition code phase for this PRN. The coefficient
        % accounts for real/complex input, and int16 samples need two bytes
        % per I or Q value.

        if strcmp(settings.dataType,'int16')
            fseek(fid, dataAdaptCoeff*2*(settings.skipNumberOfSamples + ...
                channel(channelNr).codePhase-1), 'bof');
        else
            fseek(fid, dataAdaptCoeff*(settings.skipNumberOfSamples + ...
                channel(channelNr).codePhase-1), 'bof');
        end
        % Generate the L1C data waveform. generateDataBOC11 returns the
        % BOC(1,1) waveform on a 2-sample-per-chip grid, so repeat by 6 to
        % align it with the 12-sample-per-chip TMBOC/QPSK correlator grid.
        E1BDataCode = generateE1Bcode(channel(channelNr).PRN);
        % Add guard samples at both ends for early/late indexing around the
        % code boundary.
        E1BDataCode = [E1BDataCode(codeLength*12) E1BDataCode E1BDataCode(1)]; %#ok<AGROW>

        if (settings.pilotTRKflag == 1)
            % Generate the L1C pilot TMBOC(6,1,4/33) waveform. It is already
            % produced on the 12-sample-per-chip local-code grid.
            E1CPilotCode = generateE1Ccode(channel(channelNr).PRN);
            % Add guard samples for early/late indexing.
            E1CPilotCode = [E1CPilotCode(codeLength*12) E1CPilotCode E1CPilotCode(1)];  %#ok<AGROW>
        end

        % Store data and pilot branches in the same table layout used by
        % the shared QPSK serial correlators.
        E1CodeTable = [E1BDataCode E1CPilotCode];
        if settings.correlatorType == 1
            % SIMD serial correlator expects an int32 local code table.
            E1CodeTable = int32(E1CodeTable);
        elseif settings.correlatorType == 2
            % GPU serial correlator expects an int8 local code table.
            E1CodeTable = int8(E1CodeTable);
        end

        %--- Perform various initializations ------------------------------
        % Initial code NCO frequency from acquisition.
        codeFreq      = channel(channelNr).codeFreq;
        % Residual code phase in primary-code chips.
        remCodePhase  = 0.0;
        % Carrier NCO starts from the acquired carrier frequency. The basis
        % remains fixed while the loop filter output provides corrections.
        carrFreq      = channel(channelNr).acquiredFreq;
        carrFreqBasis = channel(channelNr).acquiredFreq;
        % Residual carrier phase in radians.
        remCarrPhase  = 0.0;
        % DLL loop memory.
        oldCodeNco   = 0.0;
        oldCodeError = 0.0;
        % PLL loop memory.
        oldCarrNco   = 0.0;
        oldCarrError = 0.0;
        % Current and previous C/No estimates for display smoothing.
        CNoValue = zeros(1,3);
        tempCNoValue = zeros(1,3);

        for loopCnt =  1:NumToProcess

            %% GUI update -------------------------------------------------------------
            % Update the GUI periodically so MATLAB remains responsive
            % without spending too much time repainting the waitbar.
            if (rem(loopCnt, barUprate) == 0)

                Ln = newline;
                trackingStatus=['Tracking: Ch ', int2str(channelNr), ...
                    ' of ', int2str(TrackedNr),Ln ...
                    'PRN: ', int2str(channel(channelNr).PRN),Ln ...
                    'Completed ',int2str(loopCnt*10), ...
                    ' of ', int2str(NumToProcess*10), ' msec',Ln...
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

            %% Read next block of data ------------------------------------
            % Record the absolute sample index of the current code epoch.
            if strcmp(settings.dataType,'int16')
                trackResults(channelNr).absoluteSample(loopCnt) = ...
                    (ftell(fid))/dataAdaptCoeff/2;
            else
                trackResults(channelNr).absoluteSample(loopCnt) = ...
                    (ftell(fid))/dataAdaptCoeff;
            end
            % Code phase step in primary-code chips per input sample.
            codePhaseStep = codeFreq / settings.samplingFreq;

            % Number of input samples needed to complete the current L1C
            % primary-code epoch.
            blksize = ceil((codeLength-remCodePhase) / codePhaseStep);

            % Read the samples needed for this coherent integration.
            [rawSignal, samplesRead] = fread(fid, dataAdaptCoeff * blksize, dataConverter);

            % Stop if the input file does not contain enough samples.
            if (samplesRead ~= dataAdaptCoeff*blksize)
                disp('Not able to read the specified number of samples  for tracking, exiting!')
                delete(hwb);
                return
            end

            %% Correlator implementation  ---------------------------------
            % Carrier phase step in rad/sample based on the current
            % carrier NCO and the sampling frequency.
            carrPhaseStep = carrFreq * 2.0 * pi / settings.samplingFreq;

            if settings.correlatorType == 0 % Using MATLAB correlator
                % remCodePhase and codePhaseStep are converted from primary
                % chips to the 12x local-code grid before correlation.
                correValues = corrMatlabSerialE1(settings, rawSignal, ...
                    E1CodeTable, remCarrPhase, carrPhaseStep, remCodePhase*12,...
                    codePhaseStep*12);
            elseif settings.correlatorType == 1 % Using SIMD of CPU
                % Compile with: mex corrSIMDSerialQPSK.cpp
                % Uses the same QPSK interface and 12x code phase units.
                correValues = corrSIMDSerialQPSK(settings,rawSignal, ...
                    E1CodeTable, remCarrPhase,carrPhaseStep,remCodePhase*12, ...
                    codePhaseStep*12);
            elseif settings.correlatorType == 2 % Using GPU
                % Compile with: mexcuda corrGPUSerialQPSK.cu
                % Uses the same QPSK interface and 12x code phase units.
                correValues = corrGPUSerialQPSK(settings,rawSignal, ...
                    E1CodeTable, remCarrPhase,carrPhaseStep,remCodePhase*12, ...
                    codePhaseStep*12, channel(channelNr).PRN);
            end

            % Data branch output order: early, prompt, late; each with I/Q.
            I_E = correValues(1); Q_E = correValues(2);
            I_P = correValues(3); Q_P = correValues(4);
            I_L = correValues(5); Q_L = correValues(6);
            % Pilot branch uses the same early, prompt, late I/Q order.
            if (settings.pilotTRKflag == 1)
                p_I_E = correValues(7); p_Q_E = correValues(8);
                p_I_P = correValues(9); p_Q_P = correValues(10);
                p_I_L = correValues(11); p_Q_L = correValues(12);
            end

            % Save residual phases used by this correlation interval.
            trackResults(channelNr).remCodePhase(loopCnt) = remCodePhase;
            trackResults(channelNr).remCarrPhase(loopCnt) = remCarrPhase;
            % Carry residual phases into the next integration interval.
            remCodePhase = blksize * codePhaseStep + remCodePhase - codeLength;
            remCarrPhase = rem(carrPhaseStep * blksize + remCarrPhase, 2 * pi);

            %% Find PLL error and update carrier NCO ----------------------
            % Data prompt phase discriminator, normalized to cycles.
            carrError = atan(Q_P /I_P) / (2.0 * pi);
            if (settings.pilotTRKflag == 1)
                % Pilot prompt phase discriminator. The composite PLL error
                % gives the pilot branch three times the data-branch weight.
                pCarrError = atan(p_Q_P/p_I_P)/ (2.0 * pi);
                carrError = ( carrError + pCarrError)/2;
            end

            % Second-order carrier loop filter and NCO correction.
            carrNco = oldCarrNco + (tau2carr/tau1carr) * ...
                (carrError - oldCarrError) + carrError * (PDIcarr/tau1carr);
            oldCarrNco   = carrNco;
            oldCarrError = carrError;

            % Save the carrier frequency used for this interval, then apply
            % the latest NCO correction for the next interval.
            trackResults(channelNr).carrFreq(loopCnt) = carrFreq;
            carrFreq = carrFreqBasis + carrNco;

            %% Find DLL error and update code NCO -------------------------------------
            % Noncoherent early-minus-late envelope discriminator for the
            % data branch.
            codeError = (sqrt(I_E ^2 + Q_E ^2) - sqrt(I_L ^2 + Q_L ^2)) / ...
                (sqrt(I_E ^2 + Q_E ^2) + sqrt(I_L ^2 + Q_L ^2))* (1-earlyLateSpc);
            if (settings.pilotTRKflag == 1)
                % Pilot discriminator and weighted composite DLL error.
                pCodeError = (sqrt(p_I_E ^2 + p_Q_E ^2) - ...
                    sqrt(p_I_L ^2 + p_Q_L ^2)) / ...
                    (sqrt(p_I_E ^2 + p_Q_E ^2) + ...
                    sqrt(p_I_L ^2 + p_Q_L ^2))* (1-earlyLateSpc);
                codeError = (codeError + pCodeError)/2;
            end

            % DLL loop filter and code NCO correction.
            codeNco = oldCodeNco + (tau2code/tau1code) * ...
                (codeError - oldCodeError) + codeError * (PDIcode/tau1code);
            oldCodeNco   = codeNco;
            oldCodeError = codeError;

            % Save the code frequency used for this interval, then apply the
            % latest NCO correction for the next interval.
            trackResults(channelNr).codeFreq(loopCnt) = codeFreq;
            codeFreq = channel(channelNr).codeFreq - codeNco;

            %% Record various measures to show in postprocessing ----------------------
            trackResults(channelNr).dllDiscr(loopCnt)       = codeError;
            trackResults(channelNr).dllDiscrFilt(loopCnt)   = codeNco;
            trackResults(channelNr).pllDiscr(loopCnt)       = carrError;
            trackResults(channelNr).pllDiscrFilt(loopCnt)   = carrNco;

            trackResults(channelNr).I_E(loopCnt) = I_E;
            trackResults(channelNr).I_P(loopCnt) = I_P;
            trackResults(channelNr).I_L(loopCnt) = I_L;
            trackResults(channelNr).Q_E(loopCnt) = Q_E;
            trackResults(channelNr).Q_P(loopCnt) = Q_P;
            trackResults(channelNr).Q_L(loopCnt) = Q_L;

            if (settings.pilotTRKflag == 1)
                trackResults(channelNr).Pilot_I_P(loopCnt) = p_I_P ;
                trackResults(channelNr).Pilot_Q_P(loopCnt) = p_Q_P;
            end

            %% CNo calculation --------------------------------------------------------

            if (rem(loopCnt,settings.CNoInterval)==0)
                % Compute C/No and PLL lock detector over the configured
                % accumulation interval.
                [CNoValue, PllDetector]= ...
                    Calc_CNo_PLD(trackResults(channelNr),settings,loopCnt);
                CNoCnt = loopCnt/settings.CNoInterval;

                % Smooth displayed C/No with a 0.5/0.5 one-step average.
                trackResults(channelNr).DataCNo(CNoCnt) = ...
                    CNoValue(1) * 0.5 + tempCNoValue(1) * 0.5;
                % Save PLL lock detector output for the data branch.
                trackResults(channelNr).DataPLD(CNoCnt) = PllDetector(1);

                % Save pilot-branch and combined L1C C/No estimates.
                if (settings.pilotTRKflag == 1)
                    trackResults(channelNr).PilotCNo(CNoCnt) = ...
                        CNoValue(2) * 0.5 + tempCNoValue(2) * 0.5;
                    trackResults(channelNr).E1_CNo(CNoCnt) = ...
                        CNoValue(3) * 0.5 + tempCNoValue(3) * 0.5;
                    trackResults(channelNr).PilotPLD(CNoCnt) = PllDetector(2);
                end
            end
            tempCNoValue = CNoValue;

        end % for loopCnt

        % Tracking completed for this channel. The current implementation
        % copies the acquisition status; a lock detector can update it later.
        trackResults(channelNr).status  = channel(channelNr).status;

    end % if a PRN is assigned
end % for channelNr

% Close the waitbar
close(hwb)
clear corrSIMDSerialQPSK corrGPUSerialQPSK;
