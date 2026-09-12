function [trackResults, channel]= trkChannelsParallel(fid, channel, settings)
% Performs GLONASS L2OC code and carrier tracking in channel-parallel mode.
%
%   [trackResults, channel] = trkChannelsParallel(fid, channel, settings)
%
% This implementation tracks all active PRNs together. Each loop reads a
% shared IF data buffer large enough for the active channel set, calls the
% selected parallel correlator backend, and then updates the per-channel
% DLL/PLL states. It is intended to match the tracking behavior of
% trkChannelsSerial while reducing repeated file I/O and correlator calls.
%
%   Inputs:
%       fid             - File identifier of the IF signal record.
%       channel         - Channel state prepared from acquisition results.
%                       Each active element contains PRN, acquired carrier
%                       frequency, code phase and status.
%       settings        - Receiver settings, including sampling frequency,
%                       integration time, loop bandwidths, correlator type,
%                       file format and L2OC code definitions.
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
% (C) Developed for GLONASS L2OC SDR by Yafeng Li, Nagaraj C. Shivaramaiah
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

%% Initialize result structure ====================================
% Default channel status. It is replaced after successful tracking.
trackResults.status         = '-';      % No tracked signal, or lost lock
% Number of coherent integrations to process.
NumToProcess =  round(settings.msToProcess/1000/settings.intTime);
% Absolute sample index of each tracked L2OC code epoch in the IF record.
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
% Prompt correlator outputs from the pilot branch.
% Raw and filtered discriminator outputs from DLL and PLL.
trackResults.dllDiscr       = inf(1, NumToProcess);
trackResults.dllDiscrFilt   = inf(1, NumToProcess);
trackResults.pllDiscr       = inf(1, NumToProcess);
trackResults.pllDiscrFilt   = inf(1, NumToProcess);
% Residual code/carrier phases saved for later navigation processing.
trackResults.remCodePhase   = inf(1, NumToProcess);
trackResults.remCarrPhase   = inf(1, NumToProcess);
% C/No and PLL lock detector histories for the data branch.
tempCnt = floor(NumToProcess/settings.CNoInterval);
trackResults.DataCNo  = zeros(1,tempCnt);
trackResults.DataPLD  = zeros(1,tempCnt);
% C/No and PLL lock detector histories for the pilot and combined L2OC
% branches. These are updated once per CNoInterval integrations.

%--- Allocate one result structure per configured channel ------------------
trackResults = repmat(trackResults, 1, settings.numberOfChannels);

%% Construct variables for tracking loops =========================
% Count active channels that acquisition marked for tracking.
channelCnt = nnz([channel.status]== 'T');
% Per-channel DLL state for the recursive loop filter.
oldCodeNco   = zeros(1,channelCnt);
oldCodeError = zeros(1,channelCnt);
% Per-channel PLL state for the second-order carrier loop filter.
oldCarrNco   = zeros(1,channelCnt);
oldCarrError = zeros(1,channelCnt);
% Initial code NCO frequency for each active channel.
codeFreq     = [channel(1:channelCnt).codeFreq];
codeFreqBasis = codeFreq;
% Residual code phase carried between integrations for each channel.
remCodePhase  = zeros(1,channelCnt);
% Residual carrier phase carried between integrations for each channel.
remCarrPhase  = zeros(1,channelCnt);
% Carrier NCO starts from each channel's acquisition Doppler estimate.
carrFreq      = [channel(1:channelCnt).acquiredFreq];
% Keep acquisition carrier estimates as fixed frequency bases.
carrFreqBasis = [channel(1:channelCnt).acquiredFreq];
% Absolute IF sample index of each channel's next L2OC code epoch.
absoluteSample = settings.skipNumberOfSamples + ...
    [channel(1:channelCnt).codePhase] - 1;
fineFactor = settings.L2OCFineFactor;
% Generate the local code table for all acquired channels.
for channelNr = 1:channelCnt
    % Get a vector with the B2a data code sampled 1x/chip
        codeLength = settings.codeLength * fineFactor;
        L2OCpCode = generateL2OcpBOCCode(channel(channelNr).PRN);
        % Then make it possible to do early and late versions
        L2OCpCodeTable(channelNr,:) = [L2OCpCode(end) L2OCpCode L2OCpCode(1)];
		L2CSICode = zeros(1, codeLength);
		L2CSICodeTable(channelNr,:) = [L2CSICode(end) L2CSICode L2CSICode(1)];
end
		
		
		L2OCCodeTable = [L2OCpCodeTable L2CSICodeTable];

        if settings.correlatorType == 1
            % SIMD serial correlator expects an int32 local code table.
            L2OCCodeTable = int32(L2OCCodeTable');
        elseif settings.correlatorType == 2
            % GPU serial correlator expects an int8 local code table.
            L2OCCodeTable = int8(L2OCCodeTable');
        end

%% Initialize tracking variables ==================================
%--- DLL variables --------------------------------------------------------
% Calculate second-order DLL filter coefficients.
[tau1code, tau2code] = calcLoopCoef(settings.dllNoiseBandwidth, ...
    settings.dllDampingRatio,1.0);
%--- PLL variables --------------------------------------------------------
% Coherent integration time used by the carrier loop filter.
PDIcarr = settings.intTime;
% Calculate second-order PLL filter coefficients. For L2OC the loop is
% updated every 20 ms, so use unit loop gain; the 1 ms QPSK receivers'
% k = 0.25 setting makes the discrete carrier loop unstable here.
[tau1carr, tau2carr] = calcLoopCoef(settings.pllNoiseBandwidth, ...
    settings.pllDampingRatio, 1.0);

%% Initialize waitbar and variables for data reading ==============
% -------- Initialize waitbar ---------------------------------------------
hwb = waitbar(0,'Tracking...');
% Adjust the waitbar height to insert C/No and PRN status text.
CNoPos = get(hwb,'Position');
set(hwb,'Position',[CNoPos(1),CNoPos(2),CNoPos(3),90],'Visible','on');

% -------- Variables for IF signal reading --------------------------------
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

% Absolute sample index of the end of the current shared IF buffer.
blkEndIdx = 0;
% Nominal samples per coherent integration at the basis code frequency.
samplesPerCode = round(settings.samplingFreq / ...
    (settings.codeFreqBasis / settings.codeLength));
% Shared buffer length. One second keeps enough look-ahead for all channels
% while avoiding one fread per channel per integration.
samplesPerSec = samplesPerCode * (1/settings.intTime);

%% ============= Tracking processing for all channels =============
%  1)  GUI update
%  2)  Read next block of data
%  3)  MATLAB/SIMD/GPU correlator implementation
%  4)  Update tracking-loop variables
%  5)  CNo calculation
% =========================================================================
% The GUI bar update period depends on the correlator backend.
if (settings.correlatorType == 0)
    barUprate = 3;
else
    barUprate = 100;
end
% For C/No display in the GUI. The waitbar only reports the first tracked
% channel, while cnoSmoothState keeps each channel's smoothing history.
CNoValue = zeros(1,3);
cnoSmoothState = zeros(channelCnt,3);
% Process the configured number of coherent integrations.
for loopCnt =  1:NumToProcess
    %% GUI update -------------------------------------------------
    % Update the GUI periodically so MATLAB remains responsive without
    % spending too much time repainting the waitbar.
    if (rem(loopCnt, barUprate) == 0)
        Ln = newline;
        trackingStatus = ['Tracking ', int2str(channelCnt),' PRNs:',Ln ...
            '1st PRN: ', int2str(channel(1).PRN),Ln ...
            'Completed ',int2str(loopCnt*20), ...
            ' of ', int2str(NumToProcess*20), ' msec',Ln...
            '1st Data C/No: ',int2str(CNoValue(1)),' (dB-Hz);'];
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
    % Code phase step in chips/sample from each channel's current code NCO.
    codePhaseStep = codeFreq / settings.samplingFreq;
    % Find the size of the next code period in whole samples for all
    % active channels.
    chSampSize = ceil((settings.codeLength - remCodePhase) ./ codePhaseStep);

    % Flag passed to the parallel MEX backend so it can refresh any
    % persistent input-buffer state only when a new IF block is read.
    isDataRead = 0;
    % If the maximum index of the starting point for the next block exceeds
    % the current high boundary of the data segment, then a new block will
    % be read.
    if blkEndIdx <= max(absoluteSample + chSampSize)
        % Update data reading flag
        isDataRead = 1;
        % The shared buffer starts at the earliest channel code epoch.
        blkStartIdx = min(absoluteSample);
        % Seek the starting point of next data block to be processed.
        if strcmp(settings.dataType,'int16')
            fseek(fid, dataAdaptCoeff * blkStartIdx * 2,'bof');
        else
            fseek(fid, dataAdaptCoeff * blkStartIdx,'bof');
        end

        % Read one second of samples for reuse across channel updates.
        [rawSignal, samplesRead] = fread(fid,...
            dataAdaptCoeff * samplesPerSec, dataConverter);

        % Stop tracking cleanly if the file does not contain enough data.
        if (samplesRead ~= dataAdaptCoeff * samplesPerSec)
            disp('Not able to read the specified number of samples for tracking, exiting!')
            fclose(fid);
            return
        end
        % Update the last absolute sample covered by the shared buffer.
        blkEndIdx = blkStartIdx + samplesPerSec - 1;
    end

    %% SIMD/GPU Correlator implementation -------------------------
    % Carrier phase step in rad/sample based on the current carrier NCO and
    % the sampling frequency.
    carrPhaseStep = carrFreq * 2.0 * pi / settings.samplingFreq;
    remCodePhaseL2OC = remCodePhase;
    % Start index of each channel block within the current rawSignal buffer.
    startIdx = absoluteSample - blkStartIdx + 1;
    settingsCorr = settings;
    settingsCorr.codeLength           = settings.codeLength * fineFactor;
    settingsCorr.codeFreqBasis        = settings.codeFreqBasis * fineFactor;
    settingsCorr.dllCorrelatorSpacing = settings.dllCorrelatorSpacing * fineFactor;

    if settings.correlatorType == 0 % Using MATLAB correlator
        correValues = corrMatlabParallelTMBPSK(settingsCorr,rawSignal,L2OCCodeTable, ...
            remCarrPhase,carrPhaseStep,remCodePhaseL2OC * fineFactor,codePhaseStep * fineFactor,...
            startIdx,chSampSize,isDataRead);
    elseif settings.correlatorType == 1 % Using SIMD of CPU
        % Compile the MEX file with: mex corrSIMDParallelQPSK.cpp
        correValues = corrSIMDParallelQPSK(settingsCorr,rawSignal,L2OCCodeTable, ...
            remCarrPhase,carrPhaseStep,remCodePhaseL2OC * fineFactor,codePhaseStep * fineFactor,...
            int32(startIdx),int32(chSampSize),isDataRead);
    elseif settings.correlatorType == 2 % Using GPU
        % Compile the MEX file with: mexcuda corrGPUParallelQPSK.cu
        correValues = corrGPUParallelQPSK(settingsCorr,rawSignal,L2OCCodeTable, ...
            remCarrPhase,carrPhaseStep,remCodePhaseL2OC * fineFactor,codePhaseStep * fineFactor,...
            int32(startIdx),int32(chSampSize),isDataRead);
    end
    
    %% Update tracking-loop variables for all channels ------------
    for channelNr = 1:channelCnt
        % Extract Early / Prompt / Late correlator outputs for the current
        % channel.
        I_E = correValues(1,channelNr); Q_E = correValues(2,channelNr);
        I_P = correValues(3,channelNr); Q_P = correValues(4,channelNr);
        I_L = correValues(5,channelNr); Q_L = correValues(6,channelNr);
        % Save residual phases used by the current correlation.
        trackResults(channelNr).remCodePhase(loopCnt) = remCodePhase(channelNr);
        % Save remCarrPhase for current correlation
        trackResults(channelNr).remCarrPhase(loopCnt) = remCarrPhase(channelNr);
        % Advance residual code/carrier phase for this channel.
        remCodePhase(channelNr) = rem(chSampSize(channelNr) * ...
            codePhaseStep(channelNr) + remCodePhase(channelNr), settings.codeLength);
        remCarrPhase(channelNr) = rem(carrPhaseStep(channelNr) * ...
            chSampSize(channelNr) + remCarrPhase(channelNr), 2 * pi);
        
        % Find PLL error and update carrier NCO ---------------------------
        % Carrier phase discriminator from the data prompt correlator.
        carrError = atan(Q_P / I_P) / (2.0 * pi);
        % Implement carrier loop filter and generate NCO command
        carrNco = oldCarrNco(channelNr) + (tau2carr/tau1carr) * ...
            (carrError - oldCarrError(channelNr)) + carrError * (PDIcarr/tau1carr);
        oldCarrNco(channelNr)   = carrNco;
        oldCarrError(channelNr) = carrError;
        % Save carrier frequency for current correlation
        trackResults(channelNr).carrFreq(loopCnt) = carrFreq(channelNr);
        % Apply PLL correction around the acquisition carrier frequency.
        carrFreq(channelNr) = carrFreqBasis(channelNr) + carrNco;
        
        % Find DLL error and update code NCO ------------------------------
        % Non-coherent early-minus-late envelope discriminator for DLL.
        codeError = (sqrt(I_E^2 + Q_E^2) - sqrt(I_L^2 + Q_L^2)) / ...
                (sqrt(I_E^2 + Q_E^2) + sqrt(I_L^2 + Q_L^2));
        % Implement code loop filter and generate NCO command
        codeNco = oldCodeNco(channelNr) + (tau2code/tau1code) * ...
            (codeError - oldCodeError(channelNr)) + codeError * (settings.intTime/tau1code);
        oldCodeNco(channelNr)   = codeNco;
        oldCodeError(channelNr) = codeError;
        % Save code frequency for current correlation
        trackResults(channelNr).codeFreq(loopCnt) = codeFreq(channelNr);
        % Apply DLL correction around the acquired L2OC code frequency.
        codeFreq(channelNr) = codeFreqBasis(channelNr) - codeNco;
        
        % Record this channel's absolute sample index for the current epoch.
        trackResults(channelNr).absoluteSample(loopCnt) = absoluteSample(channelNr);
        % Advance the next start sample by the channel-specific block size.
        absoluteSample(channelNr) = absoluteSample(channelNr) + chSampSize(channelNr);

        % Record values for postprocessing and diagnostics.
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
        
        % Save additional information when the requested tracking interval
        % has completed.
        if loopCnt == NumToProcess
            % Each channel's tracked PRN
            trackResults(channelNr).PRN     = channel(channelNr).PRN;
            % If tracking reached the end of the requested interval, copy
            % the channel status. A future lock detector can replace this.
            trackResults(channelNr).status  = channel(channelNr).status;
        end

        %% CNo calculation --------------------------------------------------------

        if (rem(loopCnt,settings.CNoInterval)==0)
            % Compute C/No and PLL lock detector outputs over the latest
            % CNoInterval window.
            [currentCNoValue, PllDetector]= ...
                Calc_CNo_PLD(trackResults(channelNr),settings,loopCnt);

            CNoCnt = loopCnt/settings.CNoInterval;

            % Save C/No for data branch. A 0.5/0.5 smoother is used to
            % smooth the results
            trackResults(channelNr).DataCNo(CNoCnt) = ...
                currentCNoValue(1) * 0.5 + cnoSmoothState(channelNr,1) * 0.5;
            % Save PLL lock detector output for data channel
            trackResults(channelNr).DataPLD(CNoCnt) = PllDetector(1);
            if channelNr == 1
                % Display the first channel's current C/No values in the GUI.
                CNoValue = currentCNoValue;
            end

            % Store the unsmoothed values for the next smoothing update.
            cnoSmoothState(channelNr,:) = currentCNoValue;
        end
    end % for channelNr
end % for loopCnt
% Close the waitbar
close(hwb)
clear corrSIMDParallelQPSK corrGPUParallelQPSK;
