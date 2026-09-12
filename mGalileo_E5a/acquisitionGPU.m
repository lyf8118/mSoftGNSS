  function acqResults = acquisitionGPU(fid,settings)
%Function performs cold start acquisition on the collected "data". It
%searches for E5a signals of all satellites, which are listed in field
%"acqSatelliteList" in the settings structure. Function saves code phase
%and frequency of the detected signals in the "acqResults" structure.
%
%acqResults = GPU_acquisition(fid,settings)
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
% (C) Developed for Galileo E5a SDR by Yafeng Li, Nagaraj C. Shivaramaiah
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
% Transfer the IF data into GPU memory
gpuSignal = gpuArray(single(longSignal));
% Find sampling period
ts = 1 / settings.samplingFreq;
% Find phase points of 2ms local carrier wave (1ms for local duplicate,
% the other 1ms for zero padding)
phasePoints = gpuArray(single((0 : (samplesPerCode * 2 -1)) * 2 * pi * ts));
% Number of the frequency bins for the specified search band
numberOfFreqBins = round(settings.acqSearchBand * 2 / settings.acqSearchStep) + 1;

%--- Initialize acqResults ------------------------------------------------
% Carrier frequencies of detected signals
acqResults.carrFreq     = zeros(1, 50);
% Code frequencies of detected signals
acqResults.codeFreq     = zeros(1, 50);
% E5a code phases of detected signals
acqResults.codePhase    = zeros(1, 50);
% Correlation peak ratios of the detected signals
acqResults.peakMetric   = zeros(1, 50);

%--- Variables for fine acquisition ---------------------------------------
% Phase points of the 1 ms local carrier wave
finePhasePoints = gpuArray(single((0 : samplesPerCode-1) * 2 * pi * ts));

% Perform search for all listed PRN numbers ...
fprintf('(');
for PRN = settings.acqSatelliteList
    %% Coarse acquisition ===========================================
    
    % Generate all E5aI+E5aQ primary codes and sample them
    E5aICodesTable = makeE5aITable(PRN,settings);
    E5aQCodesTable = makeE5aQTable(PRN,settings);
    % generate local code duplicate to do correlate
    E5aICode2ms = gpuArray(single([E5aICodesTable, zeros(1,samplesPerCode)]));
    E5aQCode2ms = gpuArray(single([E5aQCodesTable, zeros(1,samplesPerCode)]));
    
    %--- Perform DFT of PRN code ------------------------------------------
    E5aICodeFreqDom = conj(fft(E5aICode2ms));
    E5aQCodeFreqDom = conj(fft(E5aQCode2ms));
    codePhaseMax = 0; freqMax = 0; peakMax = 0;
    %--- Make the correlation for all frequency bins
    for freqBinIndex = 1:numberOfFreqBins
        % Search results of one frequency bin and all code shifts
        results = zeros(1, samplesPerCode*2,'single','gpuArray');
        % Generate carrier wave frequency grid
        coarseFreqBin = settings.IF - settings.acqSearchBand + ...
            settings.acqSearchStep * (freqBinIndex - 1);
        % Generate local sine and cosine
        sigCarr = exp(-1i * coarseFreqBin * phasePoints);
        
        %--- Do non-coherent integration ----------------------------------
        for nonCohIndex = 1: settings.acqNonCohTime
            % Take 2ms vectors of input data to do correlation
            signal = gpuSignal((nonCohIndex - 1) * samplesPerCode + ...
                1 : (nonCohIndex + 1) * samplesPerCode); 
            
            % "Remove carrier" from the signal and convert the baseband 
            % signal to frequency domain
            IQfreqDom = fft(sigCarr .* signal);
            
            % Multiplication in the frequency domain (correlation in
            % time domain)
            convE5aI = IQfreqDom .* E5aICodeFreqDom;
            convE5aQ = IQfreqDom .* E5aQCodeFreqDom;
            
            % Perform inverse DFT and non-coherent integration
            results = results + abs(ifft(convE5aI)) + abs(ifft(convE5aQ));
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
    
    %--- Find 1 chip wide E5a code phase exclude range around the peak ----
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
    acqResults.codePhase(PRN) = gather(codePhaseMax);
    % Store GLRT statistic
    acqResults.peakMetric(PRN) = gather(peakMax/secondPeak);
    
    % If the result is above threshold, then there is a signal ...
    %% Fine carrier frequency search ==============================
    if acqResults.peakMetric(PRN) > settings.acqThreshold
        
        % Indicate PRN number of the detected signal
        fprintf('%02d ', PRN);
        %--- Prepare 200ms code, carrier and input signals ----------------
        % E5a data and pilot codes
        E5aICode = gpuArray(single(generateE5aIcode(PRN,1)));
        E5aQCode = gpuArray(single(generateE5aQcode(PRN,1)));
        
        % E5a codes sample index
        codeValueIndex = gpuArray(single(floor((ts * ...
            (0 : 200*samplesPerCode -1)) / (1/settings.codeFreqBasis))));
        
        % Sampled data and pilot codes
        E5aICode200ms = E5aICode((rem(codeValueIndex, settings.codeLength) + 1));
        E5aQCode200ms = E5aQCode((rem(codeValueIndex, settings.codeLength) + 1));
        
        % Take 200cm incoming signal for fine acquisition
        sig200ms = gpuSignal(codePhaseMax:codePhaseMax + 200*samplesPerCode -1);
        % One-millisecond local carrier and the double-precision starting
        % phase of each millisecond segment.
        localCarr1ms = exp(-1i * freqMax * finePhasePoints).';
        initialPhase = rem(2*pi*freqMax*(0:199)*samplesPerCode*ts, 2*pi);
        initialCarr = gpuArray(single(exp(-1i * initialPhase)));
        
        %--- Integration for each of the 200 codes ------------------------
        % Wipe off code and carrier from incoming signals
        basebandData = reshape(sig200ms .* E5aICode200ms, ...
            samplesPerCode, 200) .* localCarr1ms;
        basebandPilot = reshape(sig200ms .* E5aQCode200ms, ...
            samplesPerCode, 200) .* localCarr1ms;
        % Integration for each code
        sumPerCodeData = sum(basebandData) .* initialCarr;
        sumPerCodePilot = sum(basebandPilot) .* initialCarr;
        
        %--- Find the fine carrier freq. ----------------------------------
        % Index of the max power
        [~,maxPowerIndex]  = max(abs(fft(sumPerCodeData.^2)) + ...
            abs(fft(sumPerCodePilot.^2)));  
        % FFT shift angle
        shiftAngle = angle(exp(-2*pi*1i*(gather(maxPowerIndex)-1)/200 ))/2;
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
