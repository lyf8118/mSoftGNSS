function acqResults = acquisitionGPU(fid,settings)
% Performs cold-start acquisition for GLONASS L2OC signals using MATLAB GPU
% arrays. The function searches all PRNs listed in settings.acqSatelliteList
% and stores the detected code phase,coarse/fine carrier
% frequency, and peak metric in the acqResults structure.
%
%acqResults = acquisitionGPU(fid,settings)
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

%% Read data for acquisition ======================================
% dataAdaptCoeff converts between real IF samples and complex interleaved
% samples. For complex data, two file samples represent one complex sample.
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

% Find number of samples per full L2OC code period and per chip. In this
% receiver, settings.codeLength corresponds to one 20 ms  code period.
samplesPerCode = round(settings.samplingFreq / ...
    (settings.codeFreqBasis / settings.codeLength));
samplesPerChip   = round(settings.samplingFreq / settings.codeFreqBasis);
% Read data for acquisition. At least 200ms of signal are needed for
% fine frequency estimation
longSignal  = fread(fid, dataAdaptCoeff*12*samplesPerCode, settings.dataType)';

% Convert complex interleaved data [I0,Q0,I1,Q1,...] into MATLAB complex
% samples. The DC component is removed before forming the complex sequence.
if (dataAdaptCoeff == 2)
    longSignal = longSignal - mean(longSignal);
    longSignal = longSignal(1:2:end) + 1i .* longSignal(2:2:end);
end

%% Initialization =================================================
%--- Variables for coarse acquisition -------------------------------------
% Transfer the IF data to GPU memory as single precision to reduce memory
% traffic and match the GPU FFT / vector operations used below.
gpuSignal = gpuArray(single(longSignal));
% Sampling period.
ts = 1 / settings.samplingFreq;
% Find phase points of the local carrier wave
phasePoints = gpuArray(single((0 : (samplesPerCode * 10-1)) * 2 * pi * ts));
% Phase points of the 1 ms local carrier used in fine acquisition.
samplesPer1ms = round(settings.samplingFreq / 1000);
finePhasePoints = gpuArray(single((0 : samplesPer1ms-1) * 2 * pi * ts));
% Number of the frequency bins for the specified search band
numberOfFreqBins = round(settings.acqSearchBand * 2 / settings.acqSearchStep) + 1;
% Make index array to read L2OCp code values
fineFactor = settings.L2OCFineFactor;
tc = 1/(settings.codeFreqBasis * fineFactor);
codeValueIndex = floor( gpuArray(single(ts * (0: 10*samplesPerCode-1)/tc)) );
codeValueIndex = rem(codeValueIndex, settings.codeLength * fineFactor) + 1;


%--- Initialize acqResults ------------------------------------------------
% Carrier frequencies of detected signals
acqResults.carrFreq     = zeros(1, 32);
% Initial code NCO frequencies of detected signals
acqResults.codeFreq     = zeros(1, 32);
% L2OCp code phases of detected signals
acqResults.codePhase    = zeros(1, 32);
% Correlation peak ratios of the detected signals
acqResults.peakMetric   = zeros(1, 32);


% Perform acquisition for all requested PRN numbers. A detected PRN is
% printed as its PRN number; an undetected PRN is printed as a dot.
fprintf('(');
for PRN = settings.acqSatelliteList
    
    %% Coarse acquisition =========================================
    % Generate all codes and sample them according to the sampling freq.
    cmCodesTable = makeL2OCpTable(settings,PRN);
    % Generate local code duplicate to do correlate
    localCmCode = gpuArray(single([cmCodesTable, zeros(1,samplesPerCode)]));
    
    %--- Perform DFT of L2OCp code ------------------------------------------
    cmCodeFreqDom = conj(fft(localCmCode));
    
    codePhaseMax = 0; freqMax = 0; peakMax = 0;
    %--- Make the correlation for all frequency bins ----------------------
    for freqBinIndex = 1:numberOfFreqBins
        % Current coarse carrier frequency hypothesis.
        coarseFreqBin = settings.IF - settings.acqSearchBand + ...
            settings.acqSearchStep * (freqBinIndex - 1);
        % Local carrier used to wipe off this frequency hypothesis.
        sigCarr = exp(-1i * coarseFreqBin * phasePoints(1:2*samplesPerCode));

        % "Remove carrier" from the signal and convert the baseband signal
        % to frequency domain
        IQfreqDom = fft(sigCarr .* gpuSignal(1:samplesPerCode*2));
        
        % Multiplication in the frequency domain is correlation in the time
        % domain. The largest IFFT magnitude gives the best code phase for
        % this frequency bin.
        results = abs(ifft(IQfreqDom .* cmCodeFreqDom));
        % Keep the global maximum across all coarse frequency bins.
        [maxPeakTemp,maxIndexTemp] = max(results);
        if maxPeakTemp > peakMax
            peakMax = maxPeakTemp;
            codePhaseMax = maxIndexTemp;
            freqMax = coarseFreqBin;
        end
    end
    % To prevent index from exceeding matrix dimensions in the fine
    % acquisition, move to previous code start position.
    if codePhaseMax >= samplesPerCode
        codePhaseMax = codePhaseMax - samplesPerCode;
    end

    %--- Find 1-chip-wide L2OC code phase exclude range around the peak ---
    % The second peak is searched outside this exclusion window so the peak
    % metric is not biased by the main peak's immediate sidelobes.
    excludeIndex1 = codePhaseMax - samplesPerChip;
    excludeIndex2 = codePhaseMax + samplesPerChip;

    % Correct the search range when the exclusion window crosses the circular
    % code-phase boundary.
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
    if acqResults.peakMetric(PRN) > settings.acqThreshold
        
        % Indicate PRN number of the detected signal
        fprintf('%02d ', PRN);
        %% Fine carrier frequency search ==========================
        %--- Prepare 200ms code, carrier and input signals ----------------
        % L2OCp codes
        L2CmCode = gpuArray(single(generateL2OcpBOCCode(PRN)));
        
        % Sampled data and pilot codes
        L2CmCode200ms = L2CmCode(codeValueIndex);
        
        % Take 200cm incoming signal for fine acquisition
        sig200ms = gpuSignal(codePhaseMax:codePhaseMax + 10*samplesPerCode -1);
        % One-millisecond local carrier and the double-precision starting
        % phase of each millisecond segment.
        localCarr1ms = exp(-1i * freqMax * finePhasePoints).';
        initialPhase = rem(2*pi*freqMax*(0:199)*samplesPer1ms*ts, 2*pi);
        initialCarr = gpuArray(single(exp(-1i * initialPhase)));
        
        %--- Integration for each of the 200 codes ------------------------
        % Wipe off code and carrier from incoming signals
        basebandSig = reshape(sig200ms(1:samplesPer1ms*200) .* ...
            L2CmCode200ms(1:samplesPer1ms*200), samplesPer1ms, 200) ...
            .* localCarr1ms;
        % Integration for each 1ms code
        sumPerCode = sum(basebandSig) .* initialCarr;
        
        %--- Find the fine carrier freq. ----------------------------------
        % Index of the strongest spectral component after data-bit removal.
        [~,maxPowerIndex]  = max(abs(fft(sumPerCode.^2)));
        % Convert FFT-bin phase to carrier frequency correction. The division
        % by two compensates for the squaring operation above.
        shiftAngle = angle(exp(-2*pi*1i*(gather(maxPowerIndex)-1)/200))/2;
        acqResults.carrFreq(PRN) = freqMax - shiftAngle/0.001/2/pi;
        % Initialize the code NCO using the carrier Doppler estimate.
        acqResults.codeFreq(PRN) = settings.codeFreqBasis + ...
            (acqResults.carrFreq(PRN) - settings.IF) / ...
            settings.carrFreqBasis * settings.codeFreqBasis;

        % Signal found. If IF is zero, use 1 Hz to keep downstream tracking
        % logic from treating zero frequency as "not acquired".
        if(acqResults.carrFreq(PRN) == 0)
            acqResults.carrFreq(PRN) = 1;
        end
        
    else
        %--- No signal with this PRN --------------------------------------
        fprintf('. ');
    end   % if acqResults.peakMetric(PRN) > settings.acqThreshold

end    % for PRN = satelliteList

%=== Acquisition is over ==================================================
fprintf(')\n');
