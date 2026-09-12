function [trackResults, channel]= trkChannelsParallel(fid, channel, settings)
% Performs parallel GPS L1C code and carrier tracking for all enabled
% channels. The IF data block is read once, and the selected MATLAB, SIMD,
% or GPU correlator backend processes all active channels together.
%
% [trackResults, channel] = trkChannelsParallel(fid, channel, settings)
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
% (C) Developed for GPS L1C SDR by Yafeng Li, Nagaraj C. Shivaramaiah
% and Dennis M. Akos.
% Based on the original SoftGNSS SDR framework by Darius Plausinaitis,
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

%% Initialize result structure ====================================

% Channel status: '-' means no tracked signal or lost lock.
trackResults.status         = '-';      % No tracked signal, or lost lock
% Number of coherent integration intervals to process.
NumToProcess =  round(settings.msToProcess/1000/settings.intTime);
% Absolute sample index at the start of each tracked L1C code epoch.
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
    trackResults.L1C_CNo  = zeros(1,floor(NumToProcess/settings.CNoInterval));
end

%--- Copy initial settings for all channels -------------------------------
trackResults = repmat(trackResults, 1, settings.numberOfChannels);

%% Construct variables for tracking loops =========================
% Number of acquired channels marked for parallel tracking.
channelCnt = nnz([channel.status]== 'T');
if channelCnt == 0
    return
end
% DLL loop memory for each active channel.
oldCodeNco   = zeros(1,channelCnt);
oldCodeError = zeros(1,channelCnt);
% PLL loop memory for each active channel.
oldCarrNco   = zeros(1,channelCnt);
oldCarrError = zeros(1,channelCnt);
% Initial code NCO frequencies from acquisition.
codeFreq      = [channel(1:channelCnt).codeFreq];
% Residual code phase in primary-code chips.
remCodePhase  = zeros(1,channelCnt);
% Residual carrier phase in radians.
remCarrPhase  = zeros(1,channelCnt);
% Carrier NCO starts from the acquired carrier frequency. The basis remains
% fixed while each channel's loop filter output provides corrections.
carrFreq      = [channel(1:channelCnt).acquiredFreq];
carrFreqBasis = [channel(1:channelCnt).acquiredFreq];
% Absolute sample index of each channel's current L1C code epoch.
absoluteSample = settings.skipNumberOfSamples + ...
    [channel(1:channelCnt).codePhase] - 1;

%% Initialize tracking variables ==================================
%--- DLL variables --------------------------------------------------------
% Coherent integration time for the DLL.
PDIcode = settings.intTime;
% DLL loop filter coefficients.
[tau1code, tau2code] = calcLoopCoef(settings.dllNoiseBandwidth, ...
    settings.dllDampingRatio,1.0);
%--- PLL variables --------------------------------------------------------
% Coherent integration time for the carrier loop.
PDIcarr = settings.intTime;
% Second-order PLL loop filter coefficients.
[tau1carr, tau2carr] = calcLoopCoef(settings.pllNoiseBandwidth, ...
    settings.pllDampingRatio, 1.0);

% Weighting factor for the composite data/pilot DLL discriminator.
factor = CalcWeighingFactor(settings);

% Primary L1C code length in chips. The local code table below is expanded
% to the 12-sample-per-chip grid used by the TMBOC/QPSK correlator.
codeLength = settings.codeLength;
% Early-late spacing in primary-code chips, used by the DLL discriminator.
earlyLateSpc = settings.dllCorrelatorSpacing;
% Convert the early-late spacing to the 12x local-code grid expected by the
% shared QPSK correlator interface.
settings.dllCorrelatorSpacing = settings.dllCorrelatorSpacing*12;

%% Generate the L1C code table for acquired signals =======================
for channelNr = 1:channelCnt
    % Generate the L1C data waveform. generateDataBOC11 returns the
    % BOC(1,1) waveform on a 2-sample-per-chip grid, so repeat by 6 to
    % align it with the 12-sample-per-chip TMBOC/QPSK correlator grid.
    L1CCodeD = generateDataBOC11(settings,channel(channelNr).PRN);
    L1CCodeD = repelem(L1CCodeD, 6);
    % Add guard samples at both ends for early/late indexing around the
    % code boundary.
    L1CCodeDTable(channelNr,:)  = [L1CCodeD(codeLength*12) L1CCodeD L1CCodeD(1)]; %#ok<AGROW>

    % Generate the L1C pilot TMBOC(6,1,4/33) waveform. It is already
    % produced on the 12-sample-per-chip local-code grid.
    L1CCodeP = generatePilotTMBOC61(settings,channel(channelNr).PRN);
    % Add guard samples for early/late indexing.
    L1CCodePTable(channelNr,:)  = [L1CCodeP(codeLength*12) L1CCodeP L1CCodeP(1)]; %#ok<AGROW>
end
% Data and pilot branches are stored in the shared QPSK table layout for
% all active channels.
L1CCodeTable = [L1CCodeDTable L1CCodePTable];
if settings.correlatorType == 1
    % SIMD backend expects one int32 code table shared by all active channels.
    L1CCodeTable = int32(L1CCodeTable');
elseif settings.correlatorType == 2
    % GPU backend expects one int8 code table shared by all active channels.
    L1CCodeTable = int8(L1CCodeTable');
end

%% Initialize waitbar and variables for data reading ==============
% -------- Initialize waitbar ---------------------------------------------
hwb = waitbar(0,'Tracking...');
% Adjust the waitbar size so the status text can show PRN and C/No values.
CNoPos = get(hwb,'Position');
set(hwb,'Position',[CNoPos(1),CNoPos(2),CNoPos(3),90],'Visible','on');

% -------- Variables for IF signal reading --------------------------------
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

% Last absolute sample index covered by the current IF data buffer.
blkEndIdx = 0;
% Nominal number of input samples per L1C primary-code period.
samplesPerCode = round(settings.samplingFreq / ...
    (settings.codeFreqBasis / settings.codeLength));
% Number of L1C primary-code periods in one second.
codePeriodsPerSec = round(settings.codeFreqBasis / settings.codeLength);
% Read one second of data at a time for the parallel correlator.
samplesPerSec = samplesPerCode * codePeriodsPerSec;

%% ============= Tracking processing for all channels =============
%  1)  GUI update
%  2)  Read next block of data
%  3)  SIMD/GPU Correlator implementation
%  4)  Update tracking-loop variables
%  5)  CNo calculation
% =========================================================================
% The GUI bar update period depends on the correlator backend.
if (settings.correlatorType == 0)
    barUprate = 5;
else
    barUprate = 100;
end
% For C/No display in the GUI. The waitbar only reports the first tracked
% channel, while cnoSmoothState keeps each channel's smoothing history.
CNoValue = zeros(1,3);
cnoSmoothState = zeros(channelCnt,3);
% Process the requested number of coherent integration intervals.
for loopCnt =  1:NumToProcess
    %% GUI update -------------------------------------------------
    % Update the GUI periodically so MATLAB remains responsive without
    % spending too much time repainting the waitbar.
    if (rem(loopCnt, barUprate) == 0)
        Ln = newline;
        trackingStatus = ['Tracking ', int2str(channelCnt),' PRNs:',Ln ...
            '1st PRN: ', int2str(channel(1).PRN),Ln ...
            'Completed ',int2str(loopCnt*10), ...
            ' of ', int2str(NumToProcess*10), ' msec',Ln...
            '1st Data C/No: ',int2str(CNoValue(1)),' (dB-Hz);',...
            '   1st Pilot C/No: ',int2str(CNoValue(2)),' (dB-Hz)'];
        try
            waitbar(loopCnt/NumToProcess,hwb,trackingStatus);
        catch
            % The progress bar was closed. It is used as a signal
            % to stop, "cancel" processing. Exit.
            disp('Progress bar closed, exiting...');
            return
        end
    end

    %% Read next block of data ------------------------------------
    % Code phase step in primary-code chips per input sample for each
    % active channel.
    codePhaseStep = codeFreq / settings.samplingFreq;
    % Number of input samples needed to complete the current L1C
    % primary-code epoch for each active channel.
    chSampSize = ceil((settings.codeLength - remCodePhase) ./ codePhaseStep);

    % Marks whether a new IF data buffer is loaded during this iteration.
    isDataRead = 0;
    % Load a new shared buffer when any channel's next integration interval
    % would exceed the current buffer boundary.
    if blkEndIdx <= max(absoluteSample + chSampSize)
        % Inform the MEX correlator that the input buffer has changed.
        isDataRead = 1;
        % Read from the earliest active channel so one buffer can serve all
        % channel start positions.
        blkStartIdx = min(absoluteSample);
        % Seek to the start of the next shared IF data block.
        if strcmp(settings.dataType,'int16')
            fseek(fid, dataAdaptCoeff * blkStartIdx * 2,'bof');
        else
            fseek(fid, dataAdaptCoeff * blkStartIdx,'bof');
        end

        % Read one second of IF samples for the next shared data block.
        [rawSignal, samplesRead] = fread(fid,...
            dataAdaptCoeff * samplesPerSec, dataConverter);

        % Stop if the input file does not contain enough samples.
        if (samplesRead ~= dataAdaptCoeff * samplesPerSec)
            disp('Not able to read the specified number of samples for tracking, exiting!')
            fclose(fid);
            return
        end
        % Update the absolute sample range covered by the current buffer.
        blkEndIdx = blkStartIdx + samplesPerSec - 1;
    end

    %% SIMD/GPU Correlator implementation -------------------------
    % Carrier phase step in rad/sample based on the current carrier NCO and
    % the sampling frequency.
    carrPhaseStep = carrFreq * 2.0 * pi / settings.samplingFreq;
    % Start index of each channel block within the current rawSignal buffer.
    startIdx = absoluteSample - blkStartIdx + 1;

    if settings.correlatorType == 0 % Using MATLAB correlator
        % remCodePhase and codePhaseStep are converted from primary chips to
        % the 12x local-code grid before correlation.
        correValues = corrMatlabParallelL1C(settings,rawSignal,L1CCodeTable, ...
            remCarrPhase,carrPhaseStep,remCodePhase*12,codePhaseStep*12,...
            startIdx,chSampSize,isDataRead);
    elseif settings.correlatorType == 1 % Using SIMD of CPU
        % Compile with: mex corrSIMDParallelQPSK.cpp
        % Uses the same QPSK interface and 12x code phase units.
        correValues = corrSIMDParallelQPSK(settings,rawSignal,L1CCodeTable, ...
            remCarrPhase,carrPhaseStep,remCodePhase*12,codePhaseStep*12,...
            int32(startIdx),int32(chSampSize),isDataRead);
    elseif settings.correlatorType == 2 % Using GPU
        % Compile with: mexcuda corrGPUParallelQPSK.cu
        % Uses the same QPSK interface and 12x code phase units.
        correValues = corrGPUParallelQPSK(settings,rawSignal,L1CCodeTable, ...
            remCarrPhase,carrPhaseStep,remCodePhase*12,codePhaseStep*12,...
            int32(startIdx),int32(chSampSize),isDataRead);
    end

    %% Update tracking-loop variables for all channels ------------
    for channelNr = 1:channelCnt
        % Extract Early / Prompt / Late correlator outputs for the current
        % channel.
        I_E = correValues(1,channelNr); Q_E = correValues(2,channelNr);
        I_P = correValues(3,channelNr); Q_P = correValues(4,channelNr);
        I_L = correValues(5,channelNr); Q_L = correValues(6,channelNr);
        % Pilot branch uses the same early, prompt, late I/Q order.
        if (settings.pilotTRKflag == 1)
            p_I_E = correValues(7,channelNr); p_Q_E = correValues(8,channelNr);
            p_I_P = correValues(9,channelNr); p_Q_P = correValues(10,channelNr);
            p_I_L = correValues(11,channelNr); p_Q_L = correValues(12,channelNr);
        end

        % Save residual phases used by this correlation interval.
        trackResults(channelNr).remCodePhase(loopCnt) = remCodePhase(channelNr);
        trackResults(channelNr).remCarrPhase(loopCnt) = remCarrPhase(channelNr);

        % Carry residual phases into this channel's next integration
        % interval.
        remCodePhase(channelNr) = chSampSize(channelNr) * codePhaseStep(channelNr)...
            + remCodePhase(channelNr) - codeLength;
        remCarrPhase(channelNr) = rem(carrPhaseStep(channelNr) * ...
            chSampSize(channelNr) + remCarrPhase(channelNr), 2 * pi);

        % Find PLL error and update carrier NCO ---------------------------
        % Data prompt phase discriminator, normalized to cycles.
        carrError = atan(Q_P / I_P) / (2.0 * pi);

        if (settings.pilotTRKflag == 1)
            % Pilot prompt phase discriminator. The composite PLL error
            % gives the pilot branch three times the data-branch weight.
            pCarrError = atan(p_Q_P/p_I_P)/ (2.0 * pi);
            carrError = ( carrError + pCarrError*3)/4;
        end

        % Second-order carrier loop filter and NCO correction.
        carrNco = oldCarrNco(channelNr) + (tau2carr/tau1carr) * ...
            (carrError - oldCarrError(channelNr)) + carrError * (PDIcarr/tau1carr);
        oldCarrNco(channelNr)   = carrNco;
        oldCarrError(channelNr) = carrError;
        % Save the carrier frequency used for this interval, then apply the
        % latest NCO correction for the next interval.
        trackResults(channelNr).carrFreq(loopCnt) = carrFreq(channelNr);
        carrFreq(channelNr) = carrFreqBasis(channelNr) + carrNco;

        %% Find DLL error and update code NCO -------------------------------------
        % Noncoherent early-minus-late envelope discriminator for the data
        % branch.
        codeError = (sqrt(I_E ^2 + Q_E ^2) - sqrt(I_L ^2 + Q_L ^2)) / ...
            (sqrt(I_E ^2 + Q_E ^2) + sqrt(I_L ^2 + Q_L ^2))* (1-earlyLateSpc);
        if (settings.pilotTRKflag == 1)
            % Pilot discriminator and weighted composite DLL error.
            pCodeError = (sqrt(p_I_E ^2 + p_Q_E ^2) - ...
                sqrt(p_I_L ^2 + p_Q_L ^2)) / ...
                (sqrt(p_I_E ^2 + p_Q_E ^2) + ...
                sqrt(p_I_L ^2 + p_Q_L ^2))* (1-earlyLateSpc);
            codeError = codeError * factor + pCodeError * (1 - factor);
        end
        % DLL loop filter and code NCO correction.
        codeNco = oldCodeNco(channelNr) + (tau2code/tau1code) * ...
            (codeError - oldCodeError(channelNr)) + codeError * (PDIcode/tau1code);
        oldCodeNco(channelNr)   = codeNco;
        oldCodeError(channelNr) = codeError;
        % Save the code frequency used for this interval, then apply the
        % latest NCO correction for the next interval.
        trackResults(channelNr).codeFreq(loopCnt) = codeFreq(channelNr);
        codeFreq(channelNr) = channel(channelNr).codeFreq - codeNco;

        % Record the absolute sample index of this channel's current code
        % epoch.
        trackResults(channelNr).absoluteSample(loopCnt) = absoluteSample(channelNr);
        % Advance this channel's absolute sample index for the next epoch.
        absoluteSample(channelNr) = absoluteSample(channelNr) + chSampSize(channelNr);

        % Record various measures to show in postprocessing ---------------
        trackResults(channelNr).dllDiscr(loopCnt)     = codeError;
        trackResults(channelNr).dllDiscrFilt(loopCnt) = codeNco;
        trackResults(channelNr).pllDiscr(loopCnt)     = carrError;
        trackResults(channelNr).pllDiscrFilt(loopCnt) = carrNco;
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

        % Save additional information before tracking is over
        if loopCnt == NumToProcess
            % Save this channel's tracked PRN.
            trackResults(channelNr).PRN     = channel(channelNr).PRN;
            % Tracking completed for this channel. The current implementation
            % copies the acquisition status; a lock detector can update it later.
            trackResults(channelNr).status  = channel(channelNr).status;
        end

        %% CNo calculation --------------------------------------------------------

        if (rem(loopCnt,settings.CNoInterval)==0)
            % Compute C/No and PLL lock detector output for this channel.
            [currentCNoValue, PllDetector]= ...
                Calc_CNo_PLD(trackResults(channelNr),settings,loopCnt);

            CNoCnt = loopCnt/settings.CNoInterval;

            % Save C/No for data channel. Use this channel's own previous
            % value so parallel channels do not contaminate one another.
            trackResults(channelNr).DataCNo(CNoCnt) = ...
                currentCNoValue(1) * 0.5 + cnoSmoothState(channelNr,1) * 0.5;
            % Save PLL lock detector output for the data branch.
            trackResults(channelNr).DataPLD(CNoCnt) = PllDetector(1);

            % Save pilot-branch and combined L1C C/No estimates.
            if (settings.pilotTRKflag == 1)
                trackResults(channelNr).PilotCNo(CNoCnt) = ...
                    currentCNoValue(2) * 0.5 + cnoSmoothState(channelNr,2) * 0.5;
                trackResults(channelNr).L1C_CNo(CNoCnt) = ...
                    currentCNoValue(3) * 0.5 + cnoSmoothState(channelNr,3) * 0.5;
                trackResults(channelNr).PilotPLD(CNoCnt) = PllDetector(2);
            end
            if channelNr == 1
                % Display the first channel's current C/No values in the GUI.
                CNoValue = currentCNoValue;
            end
            % Store the unsmoothed values for this channel's next update.
            cnoSmoothState(channelNr,:) = currentCNoValue;
        end

    end % for channelNr
end % for loopCnt
% Close the waitbar
close(hwb)
clear corrSIMDParallelQPSK corrGPUParallelQPSK;
