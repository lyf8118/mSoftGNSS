function [trackResults, channel]= trkChannelsParallel(fid, channel, settings)
%TRKCHANNELSPARALLEL Performs code and carrier tracking for all active
%channels using channel-parallel tracking mode.
%
%   [trackResults, channel] = trkChannelsParallel(fid, channel, settings)
%
%   Inputs:                                                        
%       fid             - file identifier of the signal record.
%       channel         - PRN, carrier frequencies and code phases of all
%                       satellites to be tracked (prepared by preRun.m from
%                       acquisition results).
%       settings        - receiver settings.
%   Outputs:
%       trackResults    - tracking results (structure array). Contains
%                       in-phase prompt outputs, absolute spreading-code
%                       starting positions, and other observation data
%                       from the tracking loops. All are saved every
%                       millisecond.

%--------------------------------------------------------------------------
%                         CU Multi-GNSS SDR
% (C) Developed by Yafeng Li, Nagaraj C. Shivaramaiah and Dennis M. Akos
% Based on the original framework by Darius Plausinaitis,Peter Rinder,
% Nicolaj Bertelsen and Dennis M. Akos
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
% The absolute sample in the record of the C/A code start:
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
% Loop discriminators
trackResults.dllDiscr       = inf(1, settings.msToProcess);
trackResults.dllDiscrFilt   = inf(1, settings.msToProcess);
trackResults.pllDiscr       = inf(1, settings.msToProcess);
trackResults.pllDiscrFilt   = inf(1, settings.msToProcess);
% Remain code and carrier phase
trackResults.remCodePhase   = inf(1, settings.msToProcess);
trackResults.remCarrPhase   = inf(1, settings.msToProcess);
%C/No
trackResults.CNo.VSMValue = ...
    zeros(1,floor(settings.msToProcess/settings.CNo.VSMinterval));
trackResults.CNo.VSMIndex = ...
    zeros(1,floor(settings.msToProcess/settings.CNo.VSMinterval));

%--- Copy result variables for all channels -------------------------------
% Construct result structure for all channels to be tracked
trackResults = repmat(trackResults, 1, settings.numberOfChannels);

%% Construct variables for tracking loops =========================
% Number of active channels
channelCnt = nnz([channel.status]== 'T');
% Code tracking loop parameters
oldCodeNco   = zeros(1,channelCnt);
oldCodeError = zeros(1,channelCnt);
% Carrier/Costas loop parameters
oldCarrNco   = zeros(1,channelCnt);
oldCarrError = zeros(1,channelCnt);
% Define initial code frequency basis of NCO
codeFreq      = [channel(1:channelCnt).codeFreq];
codeFreqBasis = codeFreq;
% C/No computation point count
vsmCnt = zeros(1,channelCnt);
% C/No values for all channles
CNo = zeros(1,channelCnt);
% Define residual code phase (in chips)
remCodePhase  = zeros(1,channelCnt);
% Define residual carrier phase
remCarrPhase  = zeros(1,channelCnt);
% Define carrier frequency which is used over whole tracking period
carrFreq      = [channel(1:channelCnt).acquiredFreq];
carrFreqBasis = [channel(1:channelCnt).acquiredFreq];
% The absolute sample in the record of the C/A code start for all channels.
% In addition skip through that data file to start at the appropriate
% sampling points.
absoluteSample = settings.skipNumberOfSamples + ...
    [channel(1:channelCnt).codePhase] - 1;

% Generate the local code table for all acquired channels.
for channelNr = 1:channelCnt
    % Get a vector with the C/A code sampled 1x/chip
    caCode = generateCAcode(channel(channelNr).PRN);
    % Then make it possible to do early and late versions
    caCodeTable(channelNr,:) = [caCode(end) caCode caCode(1)];  %#ok<AGROW>
end

if settings.correlatorType == 1
    % SIMD backend expects one int32 code table shared by all active channels.
    caCodeTable = int32(caCodeTable');
elseif settings.correlatorType == 2
    % GPU backend expects one int8 code table shared by all active channels.
    caCodeTable = int8(caCodeTable');
end

%% Initialize tracking variables ==================================
% Code periods to be processed
codePeriods = settings.msToProcess;     % For GPS one C/A code is one ms
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
    settings.pllDampingRatio, 0.25);

%% Initialize waitbar and variables for data reading ==============
% -------- Initialize waitbar ---------------------------------------------
hwb = waitbar(0,'Tracking...');
%Adjust the size of the waitbar to insert text
CNoPos = get(hwb,'Position');
set(hwb,'Position',[CNoPos(1),CNoPos(2),CNoPos(3),90],'Visible','on');

% -------- Variables for IF signal reading --------------------------------
if (settings.fileType==1)
    dataAdaptCoeff = 1;
else
    dataAdaptCoeff = 2;
end

if settings.correlatorType == 0
    % Matlab correlator processes double-precision samples.
    dataConverter = strcat(settings.dataType,'=>double');
elseif (settings.correlatorType == 1 || settings.correlatorType == 2)
    % SIMD and GPU correlator processes int16 samples.
    dataConverter = strcat(settings.dataType,'=>int16');
end

% Initialize the absolute sample index for the end of the current signal
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
%  3)  Correlator implementation for all active channels
%  4)  Update tracking-loop variables
%  5)  CNo calculation
% =========================================================================
% The GUI bar update period depends on the correlator backend.
if (settings.correlatorType == 0)
    barUprate = 50;
else
    barUprate = 2000;
end

% Process the number of specified code periods
for loopCnt =  1:codePeriods
    %% GUI update -------------------------------------------------
    % Update the GUI periodically so Matlab remains responsive without
    % spending too much time repainting the waitbar.
    if (rem(loopCnt, barUprate) == 0)
        Ln = newline;
        trackingStatus=['Tracking ', int2str(channelCnt),' PRNs:',Ln ...
            '1st PRN: ', int2str(channel(1).PRN),Ln ...
            'Completed ',int2str(loopCnt), ...
            ' of ', int2str(codePeriods), ' msec',Ln...
            '1st C/No: ',int2str(CNo(1)),' (dB-Hz)'];
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
    % Code phase step based on the current code NCO and the fixed
    % sampling frequency.
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
        if strcmp(settings.dataType,'int16')
            fseek(fid, dataAdaptCoeff * blkStartIdx * 2,'bof');
        else
            fseek(fid, dataAdaptCoeff * blkStartIdx,'bof');
        end
        
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
    
    %% Correlator implementation ------------------------------------------
    % Carrier phase step in rad/sample based on the current carrier NCO and
    % the sampling frequency.
    carrPhaseStep = carrFreq * 2.0 * pi / settings.samplingFreq;
    % Start index of each channel block within the current rawSignal buffer.
    startIdx = absoluteSample - blkStartIdx + 1;
    
    if settings.correlatorType == 0 % Using Matlab correlator
        correValues = corrMatlabParallelBPSK(settings,rawSignal,caCodeTable, ...
            remCarrPhase,carrPhaseStep,remCodePhase,codePhaseStep,...
            startIdx,chSampSize,isDataRead);
    elseif settings.correlatorType == 1 % Using SIMD of CPU
        % Compile the MEX file with: mex corrSIMDParallelBPSK.cpp
        correValues = corrSIMDParallelBPSK(settings,rawSignal,caCodeTable, ...
            remCarrPhase,carrPhaseStep,remCodePhase,codePhaseStep,...
            int32(startIdx),int32(chSampSize),isDataRead);
    elseif settings.correlatorType == 2 % Using GPU
        % Compile the MEX file with: mexcuda corrGPUParallelBPSK.cu
        correValues = corrGPUParallelFusedBPSK(settings,rawSignal,caCodeTable, ...
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
        % Implement carrier loop discriminator (phase detector)
        carrError = atan(Q_P / I_P) / (2.0 * pi);
        % Implement carrier loop filter and generate NCO command
        carrNco = oldCarrNco(channelNr) + (tau2carr/tau1carr) * ...
            (carrError - oldCarrError(channelNr)) + carrError * (PDIcarr/tau1carr);
        oldCarrNco(channelNr)   = carrNco;
        oldCarrError(channelNr) = carrError;
        % Save carrier frequency for current correlation
        trackResults(channelNr).carrFreq(loopCnt) = carrFreq(channelNr);
        % Modify carrier freq based on NCO command
        carrFreq(channelNr) = carrFreqBasis(channelNr) + carrNco;
        
        % Find DLL error and update code NCO ------------------------------
        codeError = (sqrt(I_E^2 + Q_E^2) - sqrt(I_L^2 + Q_L^2)) / ...
                (sqrt(I_E^2 + Q_E^2) + sqrt(I_L^2 + Q_L^2));
        % Implement code loop filter and generate NCO command
        codeNco = oldCodeNco(channelNr) + (tau2code/tau1code) * ...
            (codeError - oldCodeError(channelNr)) + codeError * (PDIcode/tau1code);
        oldCodeNco(channelNr)   = codeNco;
        oldCodeError(channelNr) = codeError;
        % Save code frequency for current correlation
        trackResults(channelNr).codeFreq(loopCnt) = codeFreq(channelNr);
        % Modify code freq based on NCO command
        codeFreq(channelNr) = codeFreqBasis(channelNr) - codeNco;
        
        % Record the absolute sample index of the current code epoch.
        trackResults(channelNr).absoluteSample(loopCnt) = absoluteSample(channelNr);
        % Update the absolute sample index for the next millisecond.
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
        
        % Save additional information before tracking is over
        if loopCnt == codePeriods
            % Each channel's tracked PRN
            trackResults(channelNr).PRN     = channel(channelNr).PRN;
            % If we got so far, this means that the tracking was successful
            % Now we only copy status, but it can be update by a lock detector
            % if implemented
            trackResults(channelNr).status  = channel(channelNr).status;
        end
        
        %% CNo calculation ----------------------------------------
        if (rem(loopCnt,settings.CNo.VSMinterval) == 0)
            vsmCnt(channelNr) = vsmCnt(channelNr) + 1;
            CNoValue = CNoVSM(trackResults(channelNr).I_P(loopCnt - ...
                settings.CNo.VSMinterval+1:loopCnt),...
                trackResults(channelNr).Q_P(loopCnt-settings.CNo.VSMinterval + ...
                1:loopCnt),settings.CNo.accTime);
            trackResults(channelNr).CNo.VSMValue(vsmCnt(channelNr)) = CNoValue;
            trackResults(channelNr).CNo.VSMIndex(vsmCnt(channelNr)) = loopCnt;
            CNo(channelNr) = CNoValue;
        end
    end % for channelNr
end % for loopCnt
% Close the waitbar
close(hwb)
clear corrSIMDParallelBPSK corrGPUParallelBPSK corrGPUParallelFusedBPSK;
