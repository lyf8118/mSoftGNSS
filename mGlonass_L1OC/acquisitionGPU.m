function acqResults = acquisitionGPU(fid,settings)
% Performs cold-start acquisition for GLONASS L1OC signals using CPU-based
% MATLAB arrays. The function searches all PRNs listed in
% settings.acqSatelliteList and stores the detected L1OCd code phase, L1OCp code
% phase, coarse/fine carrier frequency, and peak metric in the acqResults
% structure.
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
% (C) Developed for GLONASS L1OC SDR by Yafeng Li, Nagaraj C. Shivaramaiah
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
% Find number of samples per spreading code
samplesPerCode = round(settings.samplingFreq / ...
    (settings.codeFreqBasis / settings.codeLength));
samplesPerChip   = round(settings.samplingFreq / settings.codeFreqBasis);

% Read a signal block for acquisition. The coarse search uses the beginning
% of this block, while fine frequency estimation below uses 10 L1OCd periods
% (200 ms). A little extra data is read to protect later indexing.
longSignal  = fread(fid, dataAdaptCoeff*102*samplesPerCode, settings.dataType)';

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
phasePoints = gpuArray(single((0 : (samplesPerCode * 100-1)) * 2 * pi * ts));
% Phase points of the 1 ms local carrier used in fine acquisition.
samplesPer1ms = round(settings.samplingFreq / 1000);
finePhasePoints = gpuArray(single((0 : samplesPer1ms-1) * 2 * pi * ts));
% Number of frequency bins in the coarse Doppler search grid.
numberOfFreqBins = round(settings.acqSearchBand * 2 / settings.acqSearchStep) + 1;
% Make index array to read L1OCd code values
tc = 1/(settings.codeFreqBasis);
codeValueIndex = floor( gpuArray(single(ts * (0: 100*samplesPerCode-1)/tc)) );
codeValueIndex = rem(codeValueIndex, settings.codeLength) + 1;

% Scratch buffer for the 4 possible L1OCp code phase positions.
powerArray = zeros(1,4,'single','gpuArray');

%--- Initialize acqResults ------------------------------------------------
% Carrier frequencies of detected signals
acqResults.carrFreq     = zeros(1, 32);
% Initial code NCO frequencies of detected signals
acqResults.codeFreq     = zeros(1, 32);
% L1OCd code phases of detected signals
acqResults.codePhase    = zeros(1, 32);
% Correlation peak ratios of the detected signals
acqResults.peakMetric   = zeros(1, 32);
% L1OCp code phase
acqResults.L1OCpCodePhase =  zeros(1, 32);

% Perform acquisition for all requested PRN numbers. A detected PRN is
% printed as its PRN number; an undetected PRN is printed as a dot.
fprintf('(');
for PRN = settings.acqSatelliteList
    
    %% Coarse acquisition =========================================
    % Generate all L1OCd codes and sample them according to the sampling freq.
    L1OCdCodesTable = makeL1OCdTable(settings,PRN);
    % can cover the complete circular code-phase search range.
    localL1OCdCode = gpuArray(single([L1OCdCodesTable, zeros(1,samplesPerCode)]));
    
    %--- Perform DFT of L1OCd code ------------------------------------------
    % Conjugating the local code spectrum implements correlation rather than
    % convolution after multiplication with the incoming signal spectrum.
    cmCodeFreqDom = conj(fft(localL1OCdCode));
    
    codePhaseMax = 0; freqMax = 0; peakMax = 0;
    %--- Make the correlation for all frequency bins ----------------------
    for freqBinIndex = 1:numberOfFreqBins
        % Current coarse carrier frequency hypothesis.
        coarseFreqBin = settings.IF - settings.acqSearchBand + ...
            settings.acqSearchStep * (freqBinIndex - 1);
        % Local carrier used to wipe off this frequency hypothesis.
        sigCarr = exp(-1i * coarseFreqBin * phasePoints(1:2*samplesPerCode));

        % Wipe off the carrier hypothesis and convert the first two L1OCd code
        % periods of the input signal to the frequency domain.
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
    if codePhaseMax > samplesPerCode
        codePhaseMax = codePhaseMax - samplesPerCode;
    end

    %--- Find 1-chip-wide L1OCd code phase exclude range around the peak ---
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
    % Save L1OCd code phase acquisition result.
    acqResults.codePhase(PRN) = gather(codePhaseMax);
    % Store peak-ratio metric used for the acquisition decision.
    acqResults.peakMetric(PRN) = gather(peakMax/secondPeak);

    % If the peak metric exceeds the acquisition threshold, refine carrier
    % frequency and determine the L1OCp code phase.
    if acqResults.peakMetric(PRN) > settings.acqThreshold
        
        % Indicate PRN number of the detected signal
        fprintf('%02d ', PRN);
        %% Fine carrier frequency search ==========================
        %--- Prepare 200ms code, carrier and input signals ----------------
        % L1OCd codes
        L1OCdCode = gpuArray(single(generateL1OCdCode(PRN,settings)));
        
        % Sampled data and pilot codes
        L1OCdCode200ms = L1OCdCode(codeValueIndex);
        
        % Take 200cm incoming signal for fine acquisition
        sig200ms = gpuSignal(codePhaseMax:codePhaseMax + 100*samplesPerCode -1);
        % One-millisecond local carrier and the double-precision starting
        % phase of each millisecond segment.
        localCarr1ms = exp(-1i * freqMax * finePhasePoints).';
        initialPhase = rem(2*pi*freqMax*(0:199)*samplesPer1ms*ts, 2*pi);
        initialCarr = gpuArray(single(exp(-1i * initialPhase)));
        
        %--- Integration for each of the 200 codes ------------------------
        % Wipe off code and carrier from incoming signals
        basebandSig = reshape(sig200ms(1:samplesPer1ms*200) .* ...
            L1OCdCode200ms(1:samplesPer1ms*200), samplesPer1ms, 200) ...
            .* localCarr1ms;
        % Sum each 1 ms segment; the resulting 200-point sequence is used for
        % fine carrier estimation. Squaring removes the 180-degree data sign
        % transitions before the FFT-based frequency estimate.
        sumPerCode = sum(basebandSig) .* initialCarr;
        
        %--- Find the fine carrier freq. ----------------------------------
        % Index of the max power
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

        %% ========== Find the L1OCp code phase ====================
        % Wipe off the refined carrier over one L1OCd period and test all 4 L1OCp
        % phase hypotheses. The strongest coherent sum selects CLCodePhase.
        sigCarr = exp(-1i * acqResults.carrFreq(PRN) * phasePoints(1:samplesPerCode));

        % ---------- Fine-grid L1OCp BOC phase search ----------
        fineFactor = settings.L1OCFineFactor;

        tcFine = 1 / (settings.codeFreqBasis * fineFactor);

        codeValueIndexFine = floor(ts * (0:samplesPerCode-1) / tcFine);
        codeValueIndexFine = rem(codeValueIndexFine, settings.codeLength * fineFactor) + 1;

        L1OCpCodeFine = gpuArray(single(generateL1OCpBOCCode(PRN,settings)));

        % Search the 4 possible CL-code alignments relative to the L1OCd period.
        for ind = 1:4
            L1OCpCodeSample = L1OCpCodeFine(codeValueIndexFine + ...
                settings.codeLength * fineFactor * (ind - 1));
            powerArray(1,ind) = abs(sum(sig200ms(1:samplesPerCode) ...
                .* L1OCpCodeSample .*sigCarr));
        end

        % Store the selected L1OCp code phase index.
        [~,tempPhase]  = max(powerArray);
        acqResults.L1OCpCodePhase(PRN) = gather(tempPhase);
    else
        %--- No signal with this PRN --------------------------------------
        fprintf('. ');
    end   % if acqResults.peakMetric(PRN) > settings.acqThreshold

end    % for PRN = satelliteList

%=== Acquisition is over ==================================================
fprintf(')\n');
