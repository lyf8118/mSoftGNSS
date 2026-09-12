  function acqResults = acquisitionCPU(fid,settings)
%Function performs cold start acquisition on the collected "data". It
%searches for E6 signals of all satellites, which are listed in field
%"acqSatelliteList" in the settings structure. Function saves code phase
%and frequency of the detected signals in the "acqResults" structure.
%
%acqResults = acquisition(fid,settings)
%
%   Inputs:
%       fid             - file identifier of the signal record.
%       settings      - Receiver settings. Provides information about
%                       sampling and intermediate frequencies and other
%                       parameters including the list of the satellites to
%                       be acquired.
%   Outputs:
%       acqResults    - Function saves code phases and frequencies of the
%                       detected signals in the "acqResults" structure. The
%                       field "carrFreq" is set to 0 if the signal is not
%                       detected for the given PRN number.

%--------------------------------------------------------------------------
%                         CU Multi-GNSS SDR
% (C) Developed for Galileo E6B/E6C SDR by Yafeng Li, Nagaraj C. Shivaramaiah
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

%% Read data for acquisition ======================================
%Initialize the multiplier to adjust for the data type
if (settings.fileType == 1)
    dataAdaptCoeff = 1;
else
    dataAdaptCoeff = 2;
end
% Move the starting point of processing. Can be used to start the
% signal processing at any point in the data record (e.g. good for long
% records or for signal processing in blocks).
if strcmp(settings.dataType,'int16')
    fseek(fid, dataAdaptCoeff*settings.skipNumberOfSamples * 2,'bof'); 
else
    fseek(fid, dataAdaptCoeff*settings.skipNumberOfSamples,'bof');
end
% Find number of samples per spreading code/chip
samplesPerCode = round(settings.samplingFreq / ...
    (settings.codeFreqBasis / settings.codeLength));
samplesPerChip   = round(settings.samplingFreq / settings.codeFreqBasis);

% At least 202ms of signal are needed for fine frequency estimation
codeLen = max(202,settings.acqNonCohTime+2); 
% Read data for acquisition.
longSignal  = fread(fid, dataAdaptCoeff*codeLen*samplesPerCode, settings.dataType)';
% Convert complex data into normal representation.
if (dataAdaptCoeff == 2) 
    longSignal = longSignal - mean(longSignal);
    longSignal = longSignal(1:2:end) + 1i .* longSignal(2:2:end);
end

%% Initialization =================================================
%--- Variables for coarse acquisition -------------------------------------
% Find sampling period
ts = 1 / settings.samplingFreq;
% Find phase points of 2ms local carrier wave (1ms for local duplicate,
% the other 1ms for zero padding)
phasePoints = (0 : (samplesPerCode * 2 -1)) * 2 * pi * ts;

% Number of the frequency bins for the specified search band
numberOfFreqBins = round(settings.acqSearchBand * 2 / settings.acqSearchStep) + 1;

%--- Initialize acqResults ------------------------------------------------
% Carrier frequencies of detected signals
acqResults.carrFreq     = zeros(1, 50);
% Code frequencies of detected signals
acqResults.codeFreq     = zeros(1, 50);
% E6b code phases of detected signals
acqResults.codePhase    = zeros(1, 50);
% Correlation peak ratios of the detected signals
acqResults.peakMetric   = zeros(1, 50);

%--- Variables for fine acquisition ---------------------------------------
% Phase points of the local carrier wave
finePhasePoints = (0 : (200*samplesPerCode-1)) * 2 * pi * ts;

% Perform search for all listed PRN numbers ...
fprintf('(');
for PRN = settings.acqSatelliteList
    
    %% Coarse acquisition ===========================================
    
    % Generate all E6B+E6C primary codes and sample them according to the
    % sampling freq.
    E6bCodesTable = makeE6BTable(PRN,settings); 
    E6cCodesTable = makeE6CTable(PRN,settings);
    % generate local code duplicate to do correlate
    localE6bCode = [E6bCodesTable, zeros(1,samplesPerCode)];
    localE6cCode = [E6cCodesTable, zeros(1,samplesPerCode)];
    % Search results of all frequency bins and code shifts (for one satellite)
    
    %--- Perform DFT of PRN code ------------------------------------------
    E6bCodeFreqDom = conj(fft(localE6bCode));
    E6cCodeFreqDom = conj(fft(localE6cCode));
    codePhaseMax = 0; freqMax = 0; peakMax = 0;
    %--- Make the correlation for all frequency bins
    for freqBinIndex = 1:numberOfFreqBins
        % Search results of one frequency bin and all code shifts
        results = zeros(1, samplesPerCode*2);
        % Generate carrier wave frequency grid
        coarseFreqBin = settings.IF - settings.acqSearchBand + ...
            settings.acqSearchStep * (freqBinIndex - 1);
        % Generate local sine and cosine
        sigCarr = exp(-1i * coarseFreqBin * phasePoints);
        %--- Do non-coherent integration ----------------------------------
        for nonCohIndex = 1: settings.acqNonCohTime
            % Take 2ms vectors of input data to do correlation
            signal = longSignal((nonCohIndex - 1) * samplesPerCode + ...
                1 : (nonCohIndex + 1) * samplesPerCode);
            % "Remove carrier" from the signal and convert the baseband 
            % signal to frequency domain
            IQfreqDom = fft(sigCarr .* signal);
            % Multiplication in the frequency domain (correlation in
            %domain)
            convE6b = IQfreqDom .* E6bCodeFreqDom;
            convE6c = IQfreqDom .* E6cCodeFreqDom;
            % Perform inverse DFT and non-coherent integration
            results = results + abs(ifft(convE6b)) + abs(ifft(convE6c));
        end
        %--- Look for correlation peaks -----------------------------------
        % Find the max correlation peak and corresponding code phase
        [maxPeakTemp,maxIndexTemp] = max(results);
        if maxPeakTemp > peakMax
            peakMax = maxPeakTemp;
            codePhaseMax = maxIndexTemp;
            freqMax = coarseFreqBin;
        end
    end % frqBinIndex = 1:numberOfFreqBins
    
    %--- Find 1 chip wide E6 code phase exclude range around the peak ----
    excludeIndex1 = codePhaseMax - samplesPerChip;
    excludeIndex2 = codePhaseMax + samplesPerChip;
    %--- Correct PRN code phase exclude range if the range includes array
    %boundaries
    if excludeIndex1 < 1
        codePhaseRange = excludeIndex2 : ...
                         (samplesPerCode + excludeIndex1);     
    elseif excludeIndex2 >= samplesPerCode
        codePhaseRange = (excludeIndex2 - samplesPerCode) : ...
                         excludeIndex1;
    else
        codePhaseRange = [1:excludeIndex1, ...
                          excludeIndex2 : samplesPerCode];
    end

    %--- Find the second highest correlation peak -------------------------
    secondPeak = max(results(codePhaseRange));
    % Save code phase acquisition result
    acqResults.codePhase(PRN) = codePhaseMax;
    % Store GLRT statistic
    acqResults.peakMetric(PRN) = peakMax/secondPeak;
    
    % If the result is above threshold, then there is a signal ...
    %% Fine carrier frequency search ==============================
    if acqResults.peakMetric(PRN) > settings.acqThreshold
        % Indicate PRN number of the detected signal
        fprintf('%02d ', PRN);
        
        %--- Prepare 200ms code, carrier and input signals ----------------
        % E6 data and pilot codes
        E6bCode = generateE6Bcode(PRN);
        E6cCode = generateE6Ccode(PRN,1);
        
        % E5a codes sample index
        codeValueIndex = floor( (0 : 200*samplesPerCode -1) * ts ...
            * settings.codeFreqBasis);
        
        % Sampled data and pilot codes
        E6bCode200ms = E6bCode((rem(codeValueIndex, settings.codeLength) + 1));
        E6cCode200ms = E6cCode((rem(codeValueIndex, settings.codeLength) + 1));
        
        % Take 200cm incoming signal for fine acquisition
         sig200ms = longSignal(codePhaseMax:codePhaseMax + 200*samplesPerCode -1);
        % Local carrier signal
        localCarr200cm = exp(-1i * freqMax * finePhasePoints);
        
        %--- Integration for each of the 200 codes ------------------------
        % Wipe off code and carrier from incoming signals
        basebandData = sig200ms .* E6bCode200ms .* localCarr200cm;
        basebandPilot = sig200ms .* E6cCode200ms .* localCarr200cm;
        % Integration for each code
        sumPerCodeData = sum(reshape(basebandData,[],200));
        sumPerCodePilot = sum(reshape(basebandPilot,[],200));
        
        %--- Find the fine carrier freq. ----------------------------------
        % Index of the max power
        [~,maxPowerIndex]  = max(abs(fft(sumPerCodeData.^2)) + ...
            abs(fft(sumPerCodePilot.^2)));  
        % FFT shift angle
        shiftAngle = angle(exp(-2*pi*1i*(maxPowerIndex-1)/200 ))/2;
        acqResults.carrFreq(PRN) = freqMax - shiftAngle/settings.intTime/2/pi;
        acqResults.codeFreq(PRN) = settings.codeFreqBasis + ...
            (acqResults.carrFreq(PRN) - settings.IF)/...
            settings.carrFreqBasis * settings.codeFreqBasis;
        % Signal found, if IF =0 just change to 1 Hz to allow processing
        if(acqResults.carrFreq(PRN) == 0)
            acqResults.carrFreq(PRN) = 1;
        end
    else
        %--- No signal with this PRN --------------------------------------
        fprintf('. ');
    end   % if (peakSize/secondPeakSize) > settings.acqThreshold
    
end    % for PRN = satelliteList

%=== Acquisition is over ==================================================
fprintf(')\n');
