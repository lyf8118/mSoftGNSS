function [trackResults, channel]= trkChannelsParallel(fid, channel, settings)
% Performs code and carrier tracking for all channels.
%
%[trackResults, channel] = trkChannelsParallel(fid, channel, settings)
%
%   Inputs:
%       fid             - file identifier of the signal record.
%       channel         - PRN, carrier frequencies and code phases of all
%                       satellites to be tracked (prepared by preRun.m from
%                       acquisition results).
%       settings        - receiver settings.
%   Outputs:
%       trackResults    - tracking results (structure array). Contains
%                       in-phase prompt outputs and absolute spreading
%                       code's starting positions, together with other
%                       observation data from the tracking loops. All are
%                       saved every millisecond.

%--------------------------------------------------------------------------
%                         CU Multi-GNSS SDR
% (C) Developed for BDS-3 B2a SDR by Yafeng Li, Nagaraj C. Shivaramaiah
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

%CVS record:
%$Id: tracking.m,v 1.14.2.31 2006/08/14 11:38:22 dpl Exp $

%% Initialize result structure ====================================
% Channel status
trackResults.status         = '-';      % No tracked signal, or lost lock
% The absolute sample in the record of the B2a code start:
trackResults.absoluteSample = zeros(1, settings.msToProcess);
% Freq of the PRN code:
trackResults.codeFreq       = inf(1, settings.msToProcess);
% Frequency of the tracked carrier wave:
trackResults.carrFreq       = inf(1, settings.msToProcess);
% Outputs from the correlators (In-phase):
trackResults.I_P            = zeros(1, settings.msToProcess);
trackResults.I_E            = zeros(1, settings.msToProcess);
trackResults.I_L            = zeros(1, settings.msToProcess);
% Outputs from the correlators (Quadrature-phase):
trackResults.Q_E            = zeros(1, settings.msToProcess);
trackResults.Q_P            = zeros(1, settings.msToProcess);
trackResults.Q_L            = zeros(1, settings.msToProcess);
% for pilot signal
if (settings.pilotTRKflag == 1)
    trackResults.Pilot_I_P      = zeros(1, settings.msToProcess);
    trackResults.Pilot_Q_P      = zeros(1, settings.msToProcess);
end
% Loop discriminators
trackResults.dllDiscr       = inf(1, settings.msToProcess);
trackResults.dllDiscrFilt   = inf(1, settings.msToProcess);
trackResults.pllDiscr       = inf(1, settings.msToProcess);
trackResults.pllDiscrFilt   = inf(1, settings.msToProcess);
% Remain code and carrier phase
trackResults.remCodePhase   = inf(1, settings.msToProcess);
trackResults.remCarrPhase   = inf(1, settings.msToProcess);
% C/No and PLL lock detector of data channel
trackResults.DataCNo  = zeros(1,floor(settings.msToProcess/settings.CNoInterval));
trackResults.DataPLD  = zeros(1,floor(settings.msToProcess/settings.CNoInterval));

% C/No and PLL lock detector of pilot channel
if (settings.pilotTRKflag == 1)
    trackResults.PilotCNo = zeros(1,floor(settings.msToProcess/settings.CNoInterval));
    trackResults.PilotPLD  = zeros(1,floor(settings.msToProcess/settings.CNoInterval));
    trackResults.B2a_CNo  = zeros(1,floor(settings.msToProcess/settings.CNoInterval));
end

%--- Copy result variables for all channels -------------------------------
% Construct result structure for all channels to be tracked
trackResults = repmat(trackResults, 1, settings.numberOfChannels);

%% Construct variables for tracking loops =========================
% Count of channels to be tracked
channelCnt = nnz([channel.status]== 'T');
% Code tracking loop parameters
oldCodeNco   = zeros(1,channelCnt);
oldCodeError = zeros(1,channelCnt);
% Carrier/Costas loop parameters
oldCarrNco   = zeros(1,channelCnt);
oldCarrError = zeros(1,channelCnt);
% Define initial code frequency basis of NCO
codeFreq      = [channel(1:channelCnt).codeFreq];
% Define residual code phase (in chips)
remCodePhase  = zeros(1,channelCnt);
% Define residual carrier phase
remCarrPhase  = zeros(1,channelCnt);
% Define carrier frequency which is used over whole tracking period
carrFreq      = [channel(1:channelCnt).acquiredFreq];
carrFreqBasis = [channel(1:channelCnt).acquiredFreq];
% Generate the B2a code table for acquired signals
for channelNr = 1:channelCnt
    % Get a vector with the B2a data code sampled 1x/chip
    B2aCodeD = generateB2aDataCode(channel(channelNr).PRN,settings);
    % Then make it possible to do early and late versions
    B2aCodeDTable(channelNr,:)  = [B2aCodeD(end) B2aCodeD B2aCodeD(1)]; %#ok<AGROW>
    % Get a vector with the B2a pilot code sampled 1x/chip
    B2aCodeP = generateB2aPilotCode(channel(channelNr).PRN,settings);
    % Then make it possible to do early and late versions
    B2aCodePTable(channelNr,:)  = [B2aCodeP(end) B2aCodeP B2aCodeP(1)]; %#ok<AGROW>
end
B2aCodeTable = [B2aCodeDTable B2aCodePTable];
if settings.correlatorType == 1
    % SIMD backend expects one int32 code table shared by all active channels.
    B2aCodeTable = int32(B2aCodeTable');
elseif settings.correlatorType == 2
    % GPU backend expects one int8 code table shared by all active channels.
    B2aCodeTable = int8(B2aCodeTable');
end

%% Initialize tracking variables ==================================
% Code periods to be processed
codePeriods = settings.msToProcess;     % For BDS one B2a code is one ms
%--- DLL variables --------------------------------------------------------
% Summation interval
PDIcode = settings.intTime;
% Calculate filter coefficient values
[tau1code, tau2code] = calcLoopCoef(settings.dllNoiseBandwidth, ...
    settings.dllDampingRatio,1.0);
%--- PLL variables --------------------------------------------------------
% Summation interval
PDIcarr = settings.intTime;
% Calculate filter coefficient values
[tau1carr, tau2carr] = calcLoopCoef(settings.pllNoiseBandwidth, ...
    settings.pllDampingRatio, 1.0);

%% Initialize waitbar and variables for data reading ==============
% -------- Initialize waitbar ---------------------------------------------
hwb = waitbar(0,'Tracking...');
%Adjust the size of the waitbar to insert text
CNoPos = get(hwb,'Position');
set(hwb,'Position',[CNoPos(1),CNoPos(2),CNoPos(3),90],'Visible','on');

% -------- Variables for IF signal reading --------------------------------
if (settings.fileType==1)
    dataAdaptCoeff=1;
else
    dataAdaptCoeff=2;
end
if strcmp(settings.dataType,'int16')
    bytesPerScalar = 2;
else
    bytesPerScalar = 1;
end
% The absolute sample in the record of the B2a code start for all channels.
% In addition skip through that data file to start at the appropriate
% sampling points. The configured skip and code phase are both expressed 
% in logical samples.
absoluteSample = settings.skipNumberOfSamples + ...
    [channel(1:channelCnt).codePhase] - 1;

if settings.correlatorType == 0
    % Matlab correlator processes double-precision samples.
    dataConverter = strcat(settings.dataType,'=>double');
elseif (settings.correlatorType == 1 || settings.correlatorType == 2)
    % SIMD and GPU correlator processes int16 samples.
    dataConverter = strcat(settings.dataType,'=>int16');
end

% Initialize the absolute sample indexe for the end of the current signal
% block read from the IF file:
blkEndIdx = 0;
% Find number of samples per spreading code
samplesPerCode = round(settings.samplingFreq / ...
    (settings.codeFreqBasis / settings.codeLength));
% Number of samples per second
samplesPerSec = samplesPerCode * 1000;

%% ============= Tracking processing for all channels =============
%  1)  GUI update
%  2)  Read next block of data
%  3)  SIMD/GPU Correlator implementation
%  4)  Update tracking-loop variables
%  5)  CNo calculation
% =========================================================================
% The GUI bar update period depends on the correlator backend.
if (settings.correlatorType == 0)
    barUprate = 50;
else
    barUprate = 2000;
end
% For C/No display in the GUI. The waitbar only reports the first tracked
% channel, while cnoSmoothState keeps each channel's smoothing history.
CNoValue = zeros(1,3);
cnoSmoothState = zeros(channelCnt,3);
% Process the number of specified code periods
for loopCnt =  1:codePeriods
    %% GUI update -------------------------------------------------
    % Update the GUI periodically so Matlab remains responsive without
    % spending too much time repainting the waitbar.
    if (rem(loopCnt, barUprate) == 0)
        Ln = newline;
        trackingStatus = ['Tracking ', int2str(channelCnt),' PRNs:',Ln ...
            '1st PRN: ', int2str(channel(1).PRN),Ln ...
            'Completed ',int2str(loopCnt), ...
            ' of ', int2str(codePeriods), ' msec',Ln...
            '1st Data C/No: ',int2str(CNoValue(1)),' (dB-Hz);',...
            '   1st Pilot C/No: ',int2str(CNoValue(2)),' (dB-Hz)'];
        try
            waitbar(loopCnt/codePeriods,hwb,trackingStatus);
        catch
            % The progress bar was closed. It is used as a signal
            % to stop, "cancel" processing. Exit.
            disp('Progress bar closed, exiting...');
            return
        end
    end

    %% Read next block of data ------------------------------------
    % Code phase step based on code freq. and sampling frequency (fixed)
    codePhaseStep = codeFreq / settings.samplingFreq;
    % Find the size of the next code period in whole samples for all
    % active channels.
    chSampSize = ceil((settings.codeLength - remCodePhase) ./ codePhaseStep);

    % Data reading flag
    isDataRead = 0;
    % If the maximum index of the starting point for the next block exceeds
    % the current high boundary of the data segment, then a new block will
    % be read.
    if blkEndIdx <= max(absoluteSample + chSampSize)
        % Update data reading flag
        isDataRead = 1;
        % The starting point is the minimum index for all channels
        blkStartIdx = min(absoluteSample);
        % Seek the starting point of next data block to be processed.
        fseek(fid, dataAdaptCoeff * bytesPerScalar * blkStartIdx, 'bof');

        % Read 1 second number of samples for the next data block
        [rawSignal, samplesRead] = fread(fid,...
            dataAdaptCoeff * samplesPerSec, dataConverter);

        % If did not read in enough samples, then could be out of
        % data - better exit
        if (samplesRead ~= dataAdaptCoeff * samplesPerSec)
            disp('Not able to read the specified number of samples for tracking, exiting!')
            fclose(fid);
            return
        end
        % Update end sample indexes of current block
        blkEndIdx = blkStartIdx + samplesPerSec - 1;
    end

    %% SIMD/GPU Correlator implementation -------------------------
    % Carrier phase step in rad/sample based on the current carrier NCO and
    % the sampling frequency.
    carrPhaseStep = carrFreq * 2.0 * pi / settings.samplingFreq;
    % Start index of each channel block within the current rawSignal buffer.
    startIdx = absoluteSample - blkStartIdx + 1;

    if settings.correlatorType == 0 % Using Matlab correlator
        correValues = corrMatlabParallelQPSK(settings,rawSignal,B2aCodeTable, ...
            remCarrPhase,carrPhaseStep,remCodePhase,codePhaseStep,...
            startIdx,chSampSize,isDataRead);
    elseif settings.correlatorType == 1 % Using SIMD of CPU
        % Compile the MEX file with: mex corrSIMDParallelQPSK.cpp
        correValues = corrSIMDParallelQPSK(settings,rawSignal,B2aCodeTable, ...
            remCarrPhase,carrPhaseStep,remCodePhase,codePhaseStep,...
            int32(startIdx),int32(chSampSize),isDataRead);
    elseif settings.correlatorType == 2 % Using GPU
        % Compile the MEX file with: mexcuda corrGPUParallelQPSK.cu
        correValues = corrGPUParallelQPSK(settings,rawSignal,B2aCodeTable, ...
            remCarrPhase,carrPhaseStep,remCodePhase,codePhaseStep,...
            int32(startIdx),int32(chSampSize),isDataRead);
    end

    %% Update tracking-loop variables for all channels ------------
    for channelNr = 1:channelCnt
        % Extract Early / Prompt / Late correlator outputs for the current
        % channel.
        I_E = correValues(1,channelNr); Q_E = correValues(2,channelNr);
        I_P = correValues(3,channelNr); Q_P = correValues(4,channelNr);
        I_L = correValues(5,channelNr); Q_L = correValues(6,channelNr);
        % For pilot channel signal tracking
        if (settings.pilotTRKflag == 1)
            % Now get early, late, and prompt values for pilot branch
            pilot_I_E = correValues(7,channelNr); pilot_Q_E = correValues(8,channelNr);
            pilot_I_P = correValues(9,channelNr); pilot_Q_P = correValues(10,channelNr);
            pilot_I_L = correValues(11,channelNr); pilot_Q_L = correValues(12,channelNr);
        end

        % Save remCodePhase for current correlation
        trackResults(channelNr).remCodePhase(loopCnt) = remCodePhase(channelNr);
        % Save remCarrPhase for current correlation
        trackResults(channelNr).remCarrPhase(loopCnt) = remCarrPhase(channelNr);
        % Update the remCodePhase and remCarrPhase
        remCodePhase(channelNr) = chSampSize(channelNr) * codePhaseStep(channelNr)...
            + remCodePhase(channelNr) - settings.codeLength;
        remCarrPhase(channelNr) = rem(carrPhaseStep(channelNr) * ...
            chSampSize(channelNr) + remCarrPhase(channelNr), 2 * pi);

        % Find PLL error and update carrier NCO ---------------------------
        % Implement carrier loop discriminator of data channel
        carrError = atan(Q_P / I_P) / (2.0 * pi);

        % Combined code tracking error estimation using data and pilot
        % chaannel signals
        if (settings.pilotTRKflag == 1)
            % B2a pilot channel carrier phase is pi/2 rad ahead of the
            % data channel carrier phase. Here we rotate the pilot
            % channel phase pi/2 back to the data channel phase.
            QI = (pilot_I_P + 1i * pilot_Q_P) * exp(-1i * pi/2);
            % atan is not affectede by the NH code modulation
            carrErrorQ = atan(imag(QI)/real(QI)) / (2.0 * pi);
            % As the data and pilot power is the same, so a simple
            % avergae is used as the carrier phase error estimate
            carrError = (carrError + carrErrorQ)/2;
        end

        % Implement carrier loop filter and generate NCO command
        carrNco = oldCarrNco(channelNr) + (tau2carr/tau1carr) * ...
            (carrError - oldCarrError(channelNr)) + carrError * (PDIcarr/tau1carr);
        oldCarrNco(channelNr)   = carrNco;
        oldCarrError(channelNr) = carrError;
        % Save carrier frequency for current correlation
        trackResults(channelNr).carrFreq(loopCnt) = carrFreq(channelNr);
        % Modify carrier freq based on NCO command
        carrFreq(channelNr) = carrFreqBasis(channelNr) + carrNco;

        %% Find DLL error and update code NCO -------------------------------------
        codeError = (sqrt(I_E^2 + Q_E^2) - sqrt(I_L^2 + Q_L^2)) / ...
            (sqrt(I_E^2 + Q_E^2) + sqrt(I_L^2 + Q_L^2));
        % Combined code tracking error estimation using data and pilot
        % chaannel signals
        if (settings.pilotTRKflag == 1)
            codeErrorQ = (sqrt(pilot_I_E^2 + pilot_Q_E^2) - ...
                sqrt(pilot_I_L^2 + pilot_Q_L^2)) / ...
                (sqrt(pilot_I_E^2 + pilot_Q_E^2) + ...
                sqrt(pilot_I_L^2 + pilot_Q_L^2));
            % Combined code tracking error
            codeError = (codeError + codeErrorQ)/2;
        end
        % Implement code loop filter and generate NCO command
        codeNco = oldCodeNco(channelNr) + (tau2code/tau1code) * ...
            (codeError - oldCodeError(channelNr)) + codeError * (PDIcode/tau1code);
        oldCodeNco(channelNr)   = codeNco;
        oldCodeError(channelNr) = codeError;
        % Save code frequency for current correlation
        trackResults(channelNr).codeFreq(loopCnt) = codeFreq(channelNr);
        % Modify code freq based on NCO command
        codeFreq(channelNr) = channel(channelNr).codeFreq - codeNco;

        % Record sample number based on the settings.dataType
        trackResults(channelNr).absoluteSample(loopCnt) = absoluteSample(channelNr);
        % Update sample index for next "loopCnt"
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
            trackResults(channelNr).Pilot_I_P(loopCnt) = pilot_I_P ;
            trackResults(channelNr).Pilot_Q_P(loopCnt) = pilot_Q_P;
        end

        % Save additional information before tracking is over
        if loopCnt == codePeriods
            % Each channel's tracked PRN
            trackResults(channelNr).PRN     = channel(channelNr).PRN;
            % If we got so far, this means that the tracking was successful
            % Now we only copy status, but it can be update by a lock detector
            % if implemented
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
            % Save PLL lock detector output for data channel
            trackResults(channelNr).DataPLD(CNoCnt) = PllDetector(1);
            
            % Save C/No and PLL lock detector output for pilot channel
            if (settings.pilotTRKflag == 1)
                trackResults(channelNr).PilotCNo(CNoCnt) = ... 
                    currentCNoValue(2) * 0.5 + cnoSmoothState(channelNr,2) * 0.5;
                trackResults(channelNr).B2a_CNo(CNoCnt) = ... 
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
