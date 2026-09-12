function [trackResults, channel]= trkChannelsSerial(fid, channel, settings)
% Performs code and carrier tracking for E5b signals of all channels.
%
%[trackResults, channel] = tracking(fid, channel, settings)
%
%   Inputs:
%       fid             - file identifier of the signal record.
%       channel         - PRN, carrier frequencies and code phases of all
%                       satellites to be tracked (prepared by preRum.m from
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
% (C) Developed for Galileo E5b SDR by Yafeng Li, Nagaraj C. Shivaramaiah 
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

% Channel status
trackResults.status         = '-';      % No tracked signal, or lost lock

% The absolute sample in the record of the E1B code start:
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
    trackResults.Pilot_I_P  = zeros(1, settings.msToProcess);
    trackResults.Pilot_Q_P  = zeros(1, settings.msToProcess);
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
    trackResults.E5b_CNo  = zeros(1,floor(settings.msToProcess/settings.CNoInterval));
end

%--- Copy initial settings for all channels -------------------------------
trackResults = repmat(trackResults, 1, settings.numberOfChannels);

%% Initialize tracking variables ==========================================

codePeriods = settings.msToProcess;     % For GAl one E5b code is one ms

%--- DLL variables --------------------------------------------------------
% Summation interval
PDIcode = settings.intTime;

% Calculate filter coefficient values
[tau1code, tau2code] = calcLoopCoef(settings.dllNoiseBandwidth, ...
    settings.dllDampingRatio, 1.0);

%--- PLL variables --------------------------------------------------------
% Calculate filter coefficient values
[pf3,pf2,pf1] = calcLoopCoefCarr(settings);

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
        
        % Move the starting point of processing. Can be used to start the
        % signal processing at any point in the data record (e.g. for long
        % records). In addition skip through that data file to start at the
        % appropriate sample (corresponding to code phase). Assumes sample
        % type is schar (or 1 byte per sample)
        if strcmp(settings.dataType,'int16')
            fseek(fid, dataAdaptCoeff*2*(settings.skipNumberOfSamples + ...
                channel(channelNr).codePhase-1), 'bof');
        else
            fseek(fid, dataAdaptCoeff*(settings.skipNumberOfSamples + ...
                channel(channelNr).codePhase-1), 'bof');
        end
        
        % Get a vector with the E5bI tiered code sampled 1x/chip
        E5bICode = generateE5bIcode(channel(channelNr).PRN,1);
        % Then make it possible to do early and late versions
        E5bICode = [E5bICode(settings.codeLength) E5bICode E5bICode(1)];  %#ok<AGROW>
        
            % Get a vector with the E5bQ primary code sampled 1x/chip
            E5bQcode = generateE5bQcode(channel(channelNr).PRN,1);
            E5bQcode = [E5bQcode(end) E5bQcode E5bQcode(1)]; %#ok<AGROW>
        E5bCodeTable = [E5bICode E5bQcode];
        if settings.correlatorType == 1
            % SIMD serial correlator expects an int32 local code table.
            E5bCodeTable = int32(E5bCodeTable);
        elseif settings.correlatorType == 2
            % GPU serial correlator expects an int8 local code table.
            E5bCodeTable = int8(E5bCodeTable);
        end
        
        %--- Perform various initializations ------------------------------
        % define initial code frequency basis of NCO
        codeFreq      = channel(channelNr).codeFreq;
        % Define residual code phase (in chips)
        remCodePhase  = 0.0;
        % Define carrier frequency which is used over whole tracking period
        carrFreq      = channel(channelNr).acquiredFreq;
        carrFreqBasis = channel(channelNr).acquiredFreq;
        % Define residual carrier phase
        remCarrPhase  = 0.0;
        
        %code tracking loop parameters
        oldCodeNco   = 0.0;
        oldCodeError = 0.0;
        
        % Carrier/Costas loop parameters
        d2CarrError  = 0.0;
        dCarrError   = 0.0;
        
        % For C/No computation
        CNoValue = zeros(1,3);
        tempCNoValue = zeros(1,3);
        
        
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
                    'Data C/No: ',int2str(CNoValue(1)),' (dB-Hz);',...
                    '   Pilot C/No: ',int2str(CNoValue(2)),' (dB-Hz)'];
                try
                    waitbar(loopCnt/codePeriods, hwb, trackingStatus);
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
                correValues = corrMatlabSerialQPSK(settings, rawSignal, ...
                    E5bCodeTable, remCarrPhase, carrPhaseStep, remCodePhase,...
                    codePhaseStep);
            elseif settings.correlatorType == 1 % Using SIMD of CPU
                % Compile the MEX file with: mex corrSIMDSerialQPSK.cpp
                correValues = corrSIMDSerialQPSK(settings,rawSignal, ...
                    E5bCodeTable, remCarrPhase,carrPhaseStep,remCodePhase, ...
                    codePhaseStep);
            elseif settings.correlatorType == 2 % Using GPU
                % Compile the MEX file with: mexcuda corrGPUSerialQPSK.cu
                correValues = corrGPUSerialQPSK(settings,rawSignal, ...
                    E5bCodeTable, remCarrPhase,carrPhaseStep,remCodePhase, ...
                    codePhaseStep, channel(channelNr).PRN);
            end

            % Now get early, late, and prompt values for data branch
            I_E = correValues(1); Q_E = correValues(2);
            I_P = correValues(3); Q_P = correValues(4);
            I_L = correValues(5); Q_L = correValues(6);
            % For pilot channel signal tracking
            if (settings.pilotTRKflag == 1)
                % Now get early, late, and prompt values for pilot branch
                pilot_I_E = correValues(7); pilot_Q_E = correValues(8);
                pilot_I_P = correValues(9); pilot_Q_P = correValues(10);
                pilot_I_L = correValues(11); pilot_Q_L = correValues(12);
            end

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
            % Combined code tracking error estimation using data and pilot
            % chaannel signals
            if (settings.pilotTRKflag == 1)                
                % Pilot channel(E5bQ) carrier phase is pi/2 rad ahead of the
                % data channel (E5bI) carrier phase. Here we rotate the E5bQ
                % phase pi/2 back to the E5bI phase.
                QI = (pilot_I_P + 1i * pilot_Q_P) * exp(-1i * pi/2);
                
                % atan is not affectede by the NH code modulation
                carrErrorQ = atan(imag(QI)/real(QI)) / (2.0 * pi);
                
                % As the data and pilot power is the same, so a simple 
                % avergae is used as the carrier phase error estimate
                carrError = (carrError + carrErrorQ)/2;
            end
            
            % Implement carrier loop filter and generate NCO command
            d2CarrError = d2CarrError + carrError * pf3;
            dCarrError  = d2CarrError + carrError * pf2 + dCarrError;
            carrNco     = dCarrError + carrError * pf1;
            
            % Save carrier frequency for current correlation
            trackResults(channelNr).carrFreq(loopCnt) = carrFreq;
            % Modify carrier freq based on NCO command
            carrFreq = carrFreqBasis + carrNco;
            
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
            codeNco = oldCodeNco + (tau2code/tau1code) * ...
                (codeError - oldCodeError) + codeError * (PDIcode/tau1code);
            oldCodeNco   = codeNco;
            oldCodeError = codeError;
            
            % Save code frequency for current correlation
            trackResults(channelNr).codeFreq(loopCnt) = codeFreq;
            % Modify code freq based on NCO command
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
                trackResults(channelNr).Pilot_I_P(loopCnt) = pilot_I_P ;
                trackResults(channelNr).Pilot_Q_P(loopCnt) = pilot_Q_P;
            end
            %% CNo calculation --------------------------------------------------------
            if (rem(loopCnt,settings.CNoInterval)==0) 
                % Computation of CNo and PLL detector output
                [CNoValue, PllDetector]= ...
                    Calc_CNo_PLD(trackResults(channelNr),settings,loopCnt);
                
                CNoCnt = loopCnt/settings.CNoInterval;
                
                % Save C/No for data channel: a o.5-0.5 filter is used to
                % smooth the results
                trackResults(channelNr).DataCNo(CNoCnt) = ...
                    CNoValue(1) * 0.5 + tempCNoValue(1) * 0.5;
                % Save PLL lock detector output for data channel
                trackResults(channelNr).DataPLD(CNoCnt) = PllDetector(1);
                
                % Save C/No and PLL lock detector output for pilot channel
                if (settings.pilotTRKflag == 1)
                trackResults(channelNr).PilotCNo(CNoCnt) = ...
                    CNoValue(2) * 0.5 + tempCNoValue(2) * 0.5;
                trackResults(channelNr).E5b_CNo(CNoCnt) = ...
                    CNoValue(3) * 0.5 + tempCNoValue(3) * 0.5;
                trackResults(channelNr).PilotPLD(CNoCnt) = PllDetector(2);
                end
            end
            tempCNoValue = CNoValue;
            
        end % for loopCnt
        
        % If we got so far, this means that the tracking was successful
        % Now we only copy status, but it can be update by a lock detector
        % if implemented
        trackResults(channelNr).status  = channel(channelNr).status;
        
    end % if a PRN is assigned
end % for channelNr

% Close the waitbar
close(hwb)
clear corrSIMDSerialQPSK corrGPUSerialQPSK;
