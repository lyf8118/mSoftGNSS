function [trackResults, channel]= trkChannelsSerial(fid, channel, settings)
%TRKCHANNELSSERIAL Performs code and carrier tracking for all acquired
%channels using channel-serial tracking mode.
%
%   [trackResults, channel] = trkChannelsSerial(fid, channel, settings)
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

%% Initialize result structure ============================================

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
trackResults.remCodePhase       = inf(1, settings.msToProcess);
trackResults.remCarrPhase       = inf(1, settings.msToProcess);

%C/No
trackResults.CNo.VSMValue = ...
    zeros(1,floor(settings.msToProcess/settings.CNo.VSMinterval));
trackResults.CNo.VSMIndex = ...
    zeros(1,floor(settings.msToProcess/settings.CNo.VSMinterval));

%--- Copy initial settings for all channels -------------------------------
trackResults = repmat(trackResults, 1, settings.numberOfChannels);

%% Initialize tracking variables ==========================================
% Signal period to be processed
codePeriods = settings.msToProcess;     % For GPS one C/A code is one ms

%--- DLL variables --------------------------------------------------------
% Summation interval
PDIcode = settings.intTime;

% Calculate filter coefficient values
[tau1code, tau2code] = calcLoopCoef(settings.dllNoiseBandwidth, ...
    settings.dllDampingRatio, 1.0);

%--- PLL variables --------------------------------------------------------
% Summation interval
PDIcarr = settings.intTime;

% Calculate filter coefficient values
[tau1carr, tau2carr] = calcLoopCoef(settings.pllNoiseBandwidth, ...
                                    settings.pllDampingRatio, 0.25);
% -------- Number of acquired signals ------------------------------------
% Number of acquired signals
TrackedNr = nnz([channel.status]== 'T');

% Start waitbar
hwb = waitbar(0,'Tracking...');

%Adjust the size of the waitbar to insert text
CNoPos=get(hwb,'Position');
set(hwb,'Position',[CNoPos(1),CNoPos(2),CNoPos(3),90],'Visible','on');

if (settings.fileType==1)
    dataAdaptCoeff=1;
else
    dataAdaptCoeff=2;
end

if settings.correlatorType == 0
    % Matlab correlator processes double-precision samples.
    dataConverter = strcat(settings.dataType,'=>double');
elseif (settings.correlatorType == 1 || settings.correlatorType == 2)
    % SIMD and GPU correlator processes int16 samples.
    dataConverter = strcat(settings.dataType,'=>int16');
end

% The GUI bar update period depends on the correlator backend.
if (settings.correlatorType == 0)
    barUprate = 50;
else
    barUprate = 2000;
end

%% Start processing channels ==============================================
for channelNr = 1:settings.numberOfChannels
    % Only process if PRN is non zero (acquisition was successful)
    if (channel(channelNr).PRN ~= 0)
        % Save additional information - each channel's tracked PRN
        trackResults(channelNr).PRN     = channel(channelNr).PRN;

        % Move the file pointer to the start of the current channel. This
        % allows tracking to begin at the acquired code phase even when the
        % data file stores either real or interleaved complex samples.
        if strcmp(settings.dataType,'int16')
            fseek(fid, dataAdaptCoeff*2*(settings.skipNumberOfSamples + ...
                channel(channelNr).codePhase-1), 'bof');
        else
            fseek(fid, dataAdaptCoeff*(settings.skipNumberOfSamples + ...
                channel(channelNr).codePhase-1), 'bof');
        end

        % Get a vector with the C/A code sampled 1x/chip
        caCode = generateCAcode(channel(channelNr).PRN);
        % Then make it possible to do early and late versions
        caCode = [caCode(1023) caCode caCode(1)];

        if settings.correlatorType == 1
            % SIMD serial correlator expects an int32 local code table.
            caCode = int32(caCode);
        elseif settings.correlatorType == 2
            % GPU serial correlator expects an int8 local code table.
            caCode = int8(caCode);
        end

        %--- Perform various initializations ------------------------------
        % define initial code frequency basis of NCO
        codeFreq      = channel(channelNr).codeFreq;
        codeFreqBasis = channel(channelNr).codeFreq;
        % define residual code phase (in chips)
        remCodePhase  = 0.0;
        % define carrier frequency which is used over whole tracking period
        carrFreq      = channel(channelNr).acquiredFreq;
        carrFreqBasis = channel(channelNr).acquiredFreq;
        % define residual carrier phase
        remCarrPhase  = 0.0;
        %code tracking loop parameters
        oldCodeNco   = 0.0;
        oldCodeError = 0.0;
        %carrier/Costas loop parameters
        oldCarrNco   = 0.0;
        oldCarrError = 0.0;
        %C/No computation
        vsmCnt  = 0;CNo = 0;

        %=== Process the number of specified code periods =================
        for loopCnt =  1:codePeriods
            %% GUI update -------------------------------------------------------------
            % Update the GUI periodically so Matlab remains responsive
            % without spending too much time repainting the waitbar.
            if (rem(loopCnt, barUprate) == 0)

                Ln = newline;
                trackingStatus=['Tracking: Ch ', int2str(channelNr), ...
                    ' of ', int2str(TrackedNr),Ln ...
                    'PRN: ', int2str(channel(channelNr).PRN),Ln ...
                    'Completed ',int2str(loopCnt), ...
                    ' of ', int2str(codePeriods), ' msec',Ln...
                    'C/No: ',CNo,' (dB-Hz)'];

                try
                    waitbar(loopCnt/codePeriods,hwb,trackingStatus);
                catch
                    % The progress bar was closed. It is used as a signal
                    % to stop, "cancel" processing. Exit.
                    disp('Progress bar closed, exiting...');
                    return
                end
            end

            %% Read next block of data ------------------------------------------------
            % Record the absolute sample index of the current code epoch.
            if strcmp(settings.dataType,'int16')
                trackResults(channelNr).absoluteSample(loopCnt) =(ftell(fid))/dataAdaptCoeff/2;
            else
                trackResults(channelNr).absoluteSample(loopCnt) =(ftell(fid))/dataAdaptCoeff;
            end
            % Update the code phase step based on the current code NCO and
            % the fixed sampling frequency.
            codePhaseStep = codeFreq / settings.samplingFreq;
            
            % Find the size of the next code period in whole samples.
            blksize = ceil((settings.codeLength-remCodePhase) / codePhaseStep);

            % Read the samples needed for the current millisecond.
            [rawSignal, samplesRead] = fread(fid, dataAdaptCoeff * blksize, dataConverter);

            % If did not read in enough samples, then could be out of
            % data - better exit
            if (samplesRead ~= dataAdaptCoeff*blksize)
                disp('Not able to read the specified number of samples  for tracking, exiting!')
				delete(hwb);
                return
            end

            %% Correlator implementation  --------------------------------------
            % Carrier phase step in rad/sample based on the current
            % carrier NCO and the sampling frequency.
            carrPhaseStep = carrFreq * 2.0 * pi / settings.samplingFreq;

            if settings.correlatorType == 0 % Using Matlab correlator
                correValues = corrMatlabSerialBPSK(settings, rawSignal, caCode, ...
                    remCarrPhase, carrPhaseStep, remCodePhase, codePhaseStep);
            elseif settings.correlatorType == 1 % Using SIMD of CPU
                % Compile the MEX file with: mex corrSIMDSerialBPSK.cpp
                correValues = corrSIMDSerialBPSK(settings,rawSignal,caCode, ...
                    remCarrPhase,carrPhaseStep,remCodePhase,codePhaseStep);
            elseif settings.correlatorType == 2 % Using GPU
                % Compile the MEX file with: mexcuda corrGPUSerialBPSK.cu
                correValues = corrGPUSerialBPSK(settings,rawSignal,caCode, ...
                    remCarrPhase,carrPhaseStep,remCodePhase,codePhaseStep,...
                    channel(channelNr).PRN);
            end

            % Extract Early / Prompt / Late correlator outputs.
            I_E = correValues(1); Q_E = correValues(2);
            I_P = correValues(3); Q_P = correValues(4);
            I_L = correValues(5); Q_L = correValues(6);
           
            % Save remCodePhase for current correlation
            trackResults(channelNr).remCodePhase(loopCnt) = remCodePhase;
            % Save remCarrPhase for current correlation
            trackResults(channelNr).remCarrPhase(loopCnt) = remCarrPhase;
            % Update the remCodePhase and remCarrPhase
            remCodePhase = blksize * codePhaseStep + remCodePhase - settings.codeLength;
            remCarrPhase = rem(carrPhaseStep * blksize + remCarrPhase, 2 * pi);
            %% Find PLL error and update carrier NCO ----------------------------------

            % Implement carrier loop discriminator (phase detector)
            carrError = atan(Q_P / I_P) / (2.0 * pi);

            % Implement carrier loop filter and generate NCO command
            carrNco = oldCarrNco + (tau2carr/tau1carr) * ...
                (carrError - oldCarrError) + carrError * (PDIcarr/tau1carr);
            oldCarrNco   = carrNco;
            oldCarrError = carrError;

            % Save carrier frequency for current correlation
            trackResults(channelNr).carrFreq(loopCnt) = carrFreq;

            % Modify carrier freq based on NCO command
            carrFreq = carrFreqBasis + carrNco;
       
            %% Find DLL error and update code NCO -------------------------------------
            codeError = (sqrt(I_E * I_E + Q_E * Q_E) - sqrt(I_L * I_L + Q_L * Q_L)) / ...
                (sqrt(I_E * I_E + Q_E * Q_E) + sqrt(I_L * I_L + Q_L * Q_L));

            % Implement code loop filter and generate NCO command
            codeNco = oldCodeNco + (tau2code/tau1code) * ...
                (codeError - oldCodeError) + codeError * (PDIcode/tau1code);
            oldCodeNco   = codeNco;
            oldCodeError = codeError;
            
            % Save code frequency for current correlation
            trackResults(channelNr).codeFreq(loopCnt) = codeFreq;

            % Modify code freq based on NCO command
            codeFreq = codeFreqBasis - codeNco;

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

            %% CNo calculation --------------------------------------
            if (rem(loopCnt,settings.CNo.VSMinterval)==0)
                vsmCnt=vsmCnt+1;
                CNoValue=CNoVSM(trackResults(channelNr).I_P(loopCnt-settings.CNo.VSMinterval+1:loopCnt),...
                    trackResults(channelNr).Q_P(loopCnt-settings.CNo.VSMinterval+1:loopCnt),settings.CNo.accTime);
                trackResults(channelNr).CNo.VSMValue(vsmCnt)=CNoValue;
                trackResults(channelNr).CNo.VSMIndex(vsmCnt)=loopCnt;
                CNo=int2str(CNoValue);
            end

        end % for loopCnt

        % If we got so far, this means that the tracking was successful
        % Now we only copy status, but it can be update by a lock detector
        % if implemented
        trackResults(channelNr).status  = channel(channelNr).status;

    end % if a PRN is assigned
end % for channelNr

% Close the waitbar
close(hwb);
clear corrSIMDSerialBPSK corrGPUSerialBPSK;
