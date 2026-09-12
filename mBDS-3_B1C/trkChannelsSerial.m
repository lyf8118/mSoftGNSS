function [trackResults, channel]= trkChannelsSerial(fid, channel, settings)
% Performs BDS-3 B1C code and carrier tracking for all channels using
% wideband correlating approach.
%
%[trackResults, channel] = trkChannelsSerial(fid, channel, settings)
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
% (C) Developed for BDS-3 B1C SDR by Yafeng Li, Nagaraj C. Shivaramaiah
% and Dennis M. Akos.
% Based on the original SoftGNSS SDR framework by Darius Plausinaitis,
% Peter Rinder, Nicolaj Bertelsen and Dennis M. Akos
%
% Reference: Li, Y., Shivaramaiah, N.C. & Akos, D.M. Design and
% implementation of an open-source BDS-3 B1C SDR receiver.
% GPS Solut (2019) 23: 60. https://doi.org/10.1007/s10291-019-0853-z
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
%$Id: fullbandTracking.m,v 1.14.2.31 2006/08/14 11:38:22 dpl Exp $

%% Initialize result structure ============================================

% Channel status
trackResults.status         = '-';      % No tracked signal, or lost lock

% Number of tracking loop updating
NumToProcess =  round(settings.msToProcess/1000/settings.intTime);

% The absolute sample in the record of the B1C code start:
trackResults.absoluteSample = zeros(1, NumToProcess);
% Freq of the PRN code:
trackResults.codeFreq       = inf(1, NumToProcess);
% Frequency of the tracked carrier wave:
trackResults.carrFreq       = inf(1, NumToProcess);

% Outputs from the correlators (In-phase):
trackResults.I_P            = zeros(1, NumToProcess);
trackResults.I_E            = zeros(1, NumToProcess);
trackResults.I_L            = zeros(1, NumToProcess);

% Outputs from the correlators (Quadrature-phase):
trackResults.Q_E            = zeros(1, NumToProcess);
trackResults.Q_P            = zeros(1, NumToProcess);
trackResults.Q_L            = zeros(1, NumToProcess);

% for pilot signal Outputs from the correlators (In-phase):
trackResults.Pilot_I_P  = zeros(1, NumToProcess);
trackResults.Pilot_Q_P  = zeros(1, NumToProcess);

% Loop discriminators
trackResults.dllDiscr       = inf(1, NumToProcess);
trackResults.dllDiscrFilt   = inf(1, NumToProcess);
trackResults.pllDiscr       = inf(1, NumToProcess);
trackResults.pllDiscrFilt   = inf(1, NumToProcess);

% Remaining code and carrier phase for each tracking update
trackResults.remCodePhase   = inf(1, NumToProcess);
trackResults.remCarrPhase   = inf(1, NumToProcess);

% C/No and PLL lock detector of data channel
trackResults.DataCNo  = zeros(1,floor(NumToProcess/settings.CNoInterval));
trackResults.DataPLD  = zeros(1,floor(NumToProcess/settings.CNoInterval));

% C/No and PLL lock detector of pilot channel
trackResults.PilotCNo = zeros(1,floor(NumToProcess/settings.CNoInterval));
trackResults.PilotPLD = zeros(1,floor(NumToProcess/settings.CNoInterval));
trackResults.B1C_CNo  = zeros(1,floor(NumToProcess/settings.CNoInterval));

%--- Copy initial settings for all channels -------------------------------
trackResults = repmat(trackResults, 1, settings.numberOfChannels);

%% Initialize tracking variables ==========================================
%--- DLL variables
% Consider the subcarrier of the BOC(1,1) modulation
codeLength = settings.codeLength;
% Early-late spacing in primary-code chips, used by the DLL discriminator.
earlyLateSpc = settings.dllCorrelatorSpacing;
% Convert the early-late spacing to the 12x local-code grid expected by the
% shared QMBOC correlator interface.
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

% Weighting factor
% For constructing of DLL composite code tracking error
if (settings.fullBandEn == 1)
    codeFactor = CalcWeighingFactor(settings);
    carrFactor = 1/4;
else
    codeFactor = 11/40;
    carrFactor = 11/40;
end

% -------- Number of acquired signals
TrackedNr = nnz([channel.status]== 'T');

% Start waitbar
hwb = waitbar(0,'Tracking...');

% Adjust the waitbar size so the status text can show PRN and C/No values.
CNoPos = get(hwb,'Position');
set(hwb,'Position',[CNoPos(1),CNoPos(2),CNoPos(3),90],'Visible','on');

if (settings.fileType == 1)
    dataAdaptCoeff = 1;
else
    dataAdaptCoeff = 2;
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
    barUprate = 5;
else
    barUprate = 200;
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

        % Get a vector of the B1C data-channel code with BOC(1,1)
        % modulation at rate of codeLength*2
        B1CData = generateDataBOC11(settings,channel(channelNr).PRN);
        % Then make it possible to do early and late versions
        B1CData = [B1CData(codeLength*12) B1CData B1CData(1)]; %#ok<AGROW>

        % Get a vector with the pilot BOC(1,1) spreading waveform
        pilotBOC11 = generatePilotBOC11(settings,channel(channelNr).PRN);
        % Then make it possible to do early and late versions
        pilotBOC11 = [pilotBOC11(codeLength*12) pilotBOC11 pilotBOC11(1)];

        % Get a vector with the pilot BOC(6,1) spreading waveform
        pilotBOC61 = generatePilotBOC61(settings,channel(channelNr).PRN);
        % Then make it possible to do early and late versions
        pilotBOC61 = [pilotBOC61(codeLength*12) pilotBOC61 pilotBOC61(1)];

        % Store data and pilot branches in the same table layout used by
        % the shared QMBOC serial correlators.
        B1CCodeTable = [B1CData pilotBOC11 pilotBOC61];
        if settings.correlatorType == 1
            % SIMD serial correlator expects an int32 local code table.
            B1CCodeTable = int32(B1CCodeTable);
        elseif settings.correlatorType == 2
            % GPU serial correlator expects an int8 local code table.
            B1CCodeTable = int8(B1CCodeTable);
        end

        %--- Perform various initializations ------------------------------

        % define initial code frequency basis of NCO,
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
        oldCarrNco   = 0.0;
        oldCarrError = 0.0;

        % For C/No computation
        CNoValue = zeros(1,3);
        tempCNoValue = zeros(1,3);

        %=== Process the number of specified code periods =================
        for loopCnt =  1:NumToProcess
            %% GUI update -------------------------------------------------
            % The GUI is updated every 200 ms. This way Matlab GUI is still
            % responsive enough. At the same time Matlab is not occupied
            % all the time with GUI task.
            if (rem(loopCnt, barUprate) == 0)
                Ln = newline;
                trackingStatus = ['Tracking: Ch ', int2str(channelNr), ...
                    ' of ', int2str(TrackedNr),Ln ...
                    'PRN: ', int2str(channel(channelNr).PRN),Ln ...
                    'Completed ',int2str(loopCnt*10), ...
                    ' of ', int2str(NumToProcess*10), ' msec',Ln...
                    'Data C/No: ',int2str(CNoValue(1)),' (dB-Hz);',...
                    '   Pilot C/No: ',int2str(CNoValue(2)),' (dB-Hz)'];

                try
                    waitbar(loopCnt/NumToProcess, ...
                        hwb, ...
                        trackingStatus);
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

            % Number of input samples needed to complete the current B1C
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

            %% Correlator implementation ==================================
            % Carrier phase step in rad/sample based on the current
            % carrier NCO and the sampling frequency.
            carrPhaseStep = carrFreq * 2.0 * pi / settings.samplingFreq;

            if settings.correlatorType == 0 % Using MATLAB correlator
                % remCodePhase and codePhaseStep are converted from primary
                % chips to the 12x local-code grid before correlation.
                correValues = corrMatlabSerialB1C(settings, rawSignal, ...
                    B1CCodeTable, remCarrPhase, carrPhaseStep, remCodePhase*12,...
                    codePhaseStep*12);
            elseif settings.correlatorType == 1 % Using SIMD of CPU
                % Compile with: mex corrSIMDSerialQMBOC.cpp
                % Uses the same QMBOC interface and 12x code phase units.
                correValues = corrSIMDSerialQMBOC(settings,rawSignal, ...
                    B1CCodeTable, remCarrPhase,carrPhaseStep,remCodePhase*12, ...
                    codePhaseStep*12);
            elseif settings.correlatorType == 2 % Using GPU
                % Compile with: mexcuda corrGPUSerialQMBOC.cu
                % Uses the same QMBOC interface and 12x code phase units.
                correValues = corrGPUSerialQMBOC(settings,rawSignal, ...
                    B1CCodeTable, remCarrPhase,carrPhaseStep,remCodePhase*12, ...
                    codePhaseStep*12, channel(channelNr).PRN);
            end

            % Data branch output order: early, prompt, late; each with I/Q.
            I_E = correValues(1); Q_E = correValues(2);
            I_P = correValues(3); Q_P = correValues(4);
            I_L = correValues(5); Q_L = correValues(6);

            % Correlation values for pilot BOC(1,1) spreading waveform
            p11_I_E = correValues(7); p11_Q_E = correValues(8);
            p11_I_P = correValues(9); p11_Q_P = correValues(10);
            p11_I_L = correValues(11); p11_Q_L = correValues(12);

            if (settings.fullBandEn == 1)
                % Correlation values for pilot BOC(6,1) spreading waveform
                p61_I_E = correValues(13); p61_Q_E = correValues(14);
                p61_I_P = correValues(15); p61_Q_P = correValues(16);
                p61_I_L = correValues(17); p61_Q_L = correValues(18);

                % Composite correlation values for the whole pilot channel
                p_I_E = -sqrt(4/33) * p61_I_E + sqrt(29/33)* p11_Q_E;
                p_Q_E = -sqrt(4/33) * p61_Q_E - sqrt(29/33)* p11_I_E;
                p_I_P = -sqrt(4/33) * p61_I_P + sqrt(29/33)* p11_Q_P;
                p_Q_P = -sqrt(4/33) * p61_Q_P - sqrt(29/33)* p11_I_P;
                p_I_L = -sqrt(4/33) * p61_I_L + sqrt(29/33)* p11_Q_L;
                p_Q_L = -sqrt(4/33) * p61_Q_L - sqrt(29/33)* p11_I_L;
            end

            % Save residual phases used by this correlation interval.
            trackResults(channelNr).remCodePhase(loopCnt) = remCodePhase;
            trackResults(channelNr).remCarrPhase(loopCnt) = remCarrPhase;
            % Carry residual phases into the next integration interval.
            remCodePhase = blksize * codePhaseStep + remCodePhase - codeLength;
            remCarrPhase = rem(carrPhaseStep * blksize + remCarrPhase, 2 * pi);

            %% Find PLL error and update carrier NCO ======================
            if (settings.fullBandEn == 0)
                % B1C pilot channel carrier phase is pi/2 rad ahead of the
                % data channel carrier phase. Here we rotate the pilot
                % channel phase pi/2 back to the data channel phase.
                p_QI = (p11_I_P + 1i * p11_Q_P) * exp(-1i * pi/2);
                % atan is not affectede by the NH code modulation
                p_carrError = atan(imag(p_QI)/real(p_QI)) / (2.0 * pi);
            elseif (settings.fullBandEn == 1)
                % Combined code tracking error estimation using data and pilot
                % chaannel signals
                % Pilot channel carrier tracking error
                p_carrError = atan(p_Q_P/p_I_P)/ (2.0 * pi);
            end

            % Implement carrier loop discriminator (phase detector)
            carrError = atan(Q_P / I_P) / (2.0 * pi);
            % Composite carrier tracking error
            carrError = carrError * carrFactor + p_carrError * (1 - carrFactor);

            % Second-order carrier loop filter and NCO correction.
            carrNco = oldCarrNco + (tau2carr/tau1carr) * ...
                (carrError - oldCarrError) + carrError * (PDIcarr/tau1carr);
            oldCarrNco   = carrNco;
            oldCarrError = carrError;

            % Save the carrier frequency used for this interval, then apply
            % the latest NCO correction for the next interval.
            trackResults(channelNr).carrFreq(loopCnt) = carrFreq;
            % Modify carrier freq based on NCO command
            carrFreq = carrFreqBasis + carrNco;

            %% Find DLL error and update code NCO =========================
            if (settings.fullBandEn == 0)
                % Combined code tracking error estimation using data and pilot
                % chaannel signals
                p_codeError = (sqrt(p11_I_E ^2 + p11_Q_E ^2) - ...
                    sqrt(p11_I_L ^2 + p11_Q_L ^2)) / ...
                    (sqrt(p11_I_E ^2 + p11_Q_E ^2) + ...
                    sqrt(p11_I_L ^2 + p11_Q_L ^2))* (1-earlyLateSpc);
            elseif (settings.fullBandEn == 1)
                % Combined code tracking error estimation using data and pilot
                % chaannel signals
                p_codeError = (sqrt(p_I_E ^2 + p_Q_E ^2) - ...
                    sqrt(p_I_L ^2 + p_Q_L ^2)) / ...
                    (sqrt(p_I_E ^2 + p_Q_E ^2) + ...
                    sqrt(p_I_L ^2 + p_Q_L ^2))* (1-earlyLateSpc);
            end

            codeError = (sqrt(I_E ^2 + Q_E ^2) - sqrt(I_L ^2 + Q_L ^2)) / ...
                (sqrt(I_E ^2 + Q_E ^2) + sqrt(I_L ^2 + Q_L ^2)) * (1-earlyLateSpc);

            % Composite code tracking error
            codeError = codeError * codeFactor + p_codeError * (1 - codeFactor);

            % Implement code loop filter and generate NCO command
            codeNco = oldCodeNco + (tau2code/tau1code) * ...
                (codeError - oldCodeError) + codeError * (PDIcode/tau1code);
            oldCodeNco   = codeNco;
            oldCodeError = codeError;

            % Save code frequency for current correlation
            trackResults(channelNr).codeFreq(loopCnt) = codeFreq;
            % Modify code freq based on NCO command
            codeFreq = channel(channelNr).codeFreq - codeNco;

            %% Record various measures to show in postprocessing ----------
            % Tracking parameters
            trackResults(channelNr).dllDiscr(loopCnt)       = codeError;
            trackResults(channelNr).dllDiscrFilt(loopCnt)   = codeNco;
            trackResults(channelNr).pllDiscr(loopCnt)       = carrError;
            trackResults(channelNr).pllDiscrFilt(loopCnt)   = carrNco;

            % Data channel correlation values
            trackResults(channelNr).I_E(loopCnt) = I_E;
            trackResults(channelNr).I_P(loopCnt) = I_P;
            trackResults(channelNr).I_L(loopCnt) = I_L;
            trackResults(channelNr).Q_E(loopCnt) = Q_E;
            trackResults(channelNr).Q_P(loopCnt) = Q_P;
            trackResults(channelNr).Q_L(loopCnt) = Q_L;

            % Pilot channel correlation values
            if (settings.fullBandEn == 0)
                trackResults(channelNr).Pilot_I_P(loopCnt) = p11_I_P ;
                trackResults(channelNr).Pilot_Q_P(loopCnt) = p11_Q_P;
            elseif (settings.fullBandEn == 1)
                trackResults(channelNr).Pilot_I_P(loopCnt) = p_I_P;
                trackResults(channelNr).Pilot_Q_P(loopCnt) = p_Q_P;
            end
            %% CNo calculation --------------------------------------------

            if (rem(loopCnt,settings.CNoInterval)==0)
                % Computation of CNo and PLL detector output
                [CNoValue, PllDetector]= ...
                    Calc_CNo_PLD_QMBOC(trackResults(channelNr),settings,loopCnt);

                CNoCnt = loopCnt/settings.CNoInterval;
                % Save C/No for data channel: a o.5-0.5 filter is used to
                % smooth the results
                trackResults(channelNr).DataCNo(CNoCnt) = ...
                    CNoValue(1) * 0.5 + tempCNoValue(1) * 0.5;
                % Save PLL lock detector output for data channel
                trackResults(channelNr).DataPLD(CNoCnt) = PllDetector(1);

                % Save C/No and PLL lock detector output for pilot channel
                trackResults(channelNr).PilotCNo(CNoCnt) = ...
                    CNoValue(2) * 0.5 + tempCNoValue(2) * 0.5;
                trackResults(channelNr).B1C_CNo(CNoCnt) = ...
                    CNoValue(3) * 0.5 + tempCNoValue(3) * 0.5;
                trackResults(channelNr).PilotPLD(CNoCnt) = PllDetector(2);
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
clear corrSIMDSerialQMBOC corrGPUSerialQMBOC;
