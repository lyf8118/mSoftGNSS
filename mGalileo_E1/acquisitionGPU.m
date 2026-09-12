function acqResults = acquisitionGPU(fid, settings)
%Function performs cold start acquisition on the collected "data". It
%searches for Galileo E1 signals of all satellites, which are listed in field
%"acqSatelliteList" in the settings structure. Function saves code phase
%and frequency of the detected signals in the "acqResults" structure.
%
%acqResults = acquisition(longSignal, settings)
%
%   Inputs:
%       longSignal    - 20 ms of raw IF signal from the front-end.
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
% (C) Developed for Galileo E1 SDR by Yafeng Li, Nagaraj C. Shivaramaiah
% and Dennis M. Akos.
% Based on the original framework for GPS C/A SDR by Darius Plausinaitis,
% Peter Rinder, Nicolaj Bertelsen and Dennis M. Akos
%
% Reference: Adapted within the CU Multi-GNSS SDR receiver framework for Galileo E1.
% Signal-specific comments in this file refer to Galileo E1.
%
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
%$Id: acquisition.m,v 1.1.2.12 2006/08/14 12:08:03 dpl Exp $

%% Read data for acquisition ==============================================

% Find number of samples per spreading code
samplesPerCode = round(settings.samplingFreq / ...
    (settings.codeFreqBasis / settings.codeLength));
samplesPerChip   = round(settings.samplingFreq / settings.codeFreqBasis);
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
% Read data for acquisition: use 50 code length data to do fine
% acquisition; To ensure at least 50 complete E1 code periods are
% included, 51 code periods of data are needed
longSignal  = fread(fid, dataAdaptCoeff*51*samplesPerCode, settings.dataType)';
% Convert complex data into normal representation.
if (dataAdaptCoeff == 2)
    longSignal = longSignal(1:2:end) + 1i .* longSignal(2:2:end);
end
% Transfer the IF data into GPU memory
gpuSignal = gpuArray(single(longSignal));

%% Acquisition initialization =============================================
% Find sampling period
ts = 1 / settings.samplingFreq;
% Find phase points of the local carrier wave
phasePoints = gpuArray(single((0 : (samplesPerCode*2-1)) * 2 * pi * ts));

% Number of the frequency bins for the given acquisition band
numberOfFrqBins = round(settings.acqSearchBand * 2 / settings.acqStep) + 1;

%--- Initialize acqResults and related variables --------------------------
% Carrier frequencies of detected signals
acqResults.carrFreq     = zeros(1, max(settings.acqSatelliteList));
% Code frequencies of detected signals
acqResults.codeFreq     = zeros(1, max(settings.acqSatelliteList));
% PRN code phases of detected signals
acqResults.codePhase    = zeros(1, max(settings.acqSatelliteList));
% Correlation peak ratios of the detected signals
acqResults.peakMetric   = zeros(1, max(settings.acqSatelliteList));

%--- Variables for fine acquisition ---------------------------------------
% Code resampling index used by both E1C code fine acquisition and CL phase search.
% codeValueIndex maps each input sample to the corresponding local code chip.
tc = 1/(settings.codeFreqBasis*12);
codeValueIndex = floor( gpuArray(single(ts * (0: 50*samplesPerCode-1)/tc)) );
codeValueIndex = rem(codeValueIndex, settings.codeLength*12) + 1;

% Phase points of the 1 ms local carrier wave
samplesPer1ms = round(settings.samplingFreq / 1000);
finePhasePoints = gpuArray(single((0 : samplesPer1ms-1) * 2 * pi * ts));

% Perform search for all listed PRN numbers ...
fprintf('(');
for PRN = settings.acqSatelliteList
    %% Coarse acquisition =================================================

    % Generate E1B and E1C codes and sample them according to the sampling freq.
    E1bTable = makeE1BTable(settings,PRN);
    E1cTable = makeE1CTable(settings,PRN);
    % Add zero-padding samples
    E1bCodeSample = gpuArray(single([E1bTable zeros(1,samplesPerCode)]));
    E1cCodeSample = gpuArray(single([E1cTable zeros(1,samplesPerCode)]));

    %--- Perform DFT of E1 code ------------------------------------------
    E1bCodeFreqDom = conj(fft(E1bCodeSample));
    E1cCodeFreqDom = conj(fft(E1cCodeSample));

    codePhaseMax = 0; freqMax = 0; peakMax = 0;
    %--- Make the correlation for whole frequency band (for all freq. bins)
    for frqBinIndex = 1:numberOfFrqBins
        %--- Generate carrier wave frequency grid  -----------------------
        frqBins = settings.IF - settings.acqSearchBand + ...
            settings.acqStep * (frqBinIndex - 1);
        % Generate local sine and cosine
        sigCarr = exp(-1i * frqBins * phasePoints);
        % "Remove carrier" from the signal
        I1      = real(sigCarr .* gpuSignal(1:samplesPerCode*2));
        Q1      = imag(sigCarr .* gpuSignal(1:samplesPerCode*2));
        % Convert the baseband signal to frequency domain
        IQfreqDom = fft(I1 + 1i*Q1);
        % Multiplication in the frequency domain (correlation in time domain)
        convCodeIQ1 = IQfreqDom .* E1bCodeFreqDom;
        % Perform inverse DFT and store correlation results
        results = abs(ifft(convCodeIQ1));

        % Pilot signal components
        convCodeIQ1 = IQfreqDom .* E1cCodeFreqDom;
        % Non-coherent combining of data and pilot results
        results = results + abs(ifft(convCodeIQ1));

        %--- Look for correlation peaks -----------------------------------
        % Find the max correlation peak and corresponding code phase
        [maxPeakTemp,maxIndexTemp] = max(results);
        if maxPeakTemp > peakMax
            peakMax = maxPeakTemp;
            codePhaseMax = maxIndexTemp;
            freqMax = frqBins;
        end

    end % frqBinIndex = 1:numberOfFrqBins

    %--- Find 1 chip wide E1 code phase exclude range around the peak ----
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

    % To prevent index from exceeding matrix dimensions in the fine
    % acquisition, move to previous code start position.
    if (codePhaseMax + 50*samplesPerCode -1) > length(gpuSignal)
        codePhaseMax = codePhaseMax - samplesPerCode;
    end

    % If the result is above threshold, then there is a signal ...
    if acqResults.peakMetric(PRN) > settings.acqThreshold
        %% Fine resolution frequency search =========================
        % Indicate PRN number of the detected signal
        fprintf('%02d ', PRN);

        % Generate one unresampled E1C code period.
        E1CPilotCode = gpuArray(single(generateE1Ccode(PRN)));

        % Resample the local pilot code over 200 ms
        E1CPilotCode200ms = E1CPilotCode(codeValueIndex);

        % Extract 10 E1C code periods (200 ms) of incoming signal starting from the
        % detected E1C code phase.
        sig200ms = gpuSignal(codePhaseMax:codePhaseMax + 50*samplesPerCode -1);

        % One-millisecond local carrier and the double-precision starting
        % phase of each millisecond segment.
        localCarr1ms = exp(-1i * freqMax * finePhasePoints).';
        initialPhase = rem(2*pi*freqMax*(0:199)*samplesPer1ms*ts, 2*pi);
        initialCarr = gpuArray(single(exp(-1i * initialPhase)));

        %--- Integration for each 1 ms segment over 200 ms -----------------
        % Wipe off the E1C code and coarse carrier from the incoming signal.
        basebandSig = reshape(sig200ms .* E1CPilotCode200ms, ...
            samplesPer1ms, 200) .* localCarr1ms;
        % Sum each 1 ms segment; the resulting 200-point sequence is used for
        % fine carrier estimation. Squaring removes the 180-degree data sign
        % transitions before the FFT-based frequency estimate.
        sumPerCode = sum(basebandSig) .* initialCarr;

        %--- Find the fine carrier freq. ----------------------------------
        % Index of the strongest spectral component after data-bit removal.
        [~,maxPowerIndex]  = max(abs(fft(sumPerCode.^2)));
        % Convert FFT-bin phase to carrier frequency correction. The division
        % by two compensates for the squaring operation above.
        shiftAngle = angle(exp(-2*pi*1i*(gather(maxPowerIndex)-1)/200))/2;
        acqResults.carrFreq(PRN) = freqMax - shiftAngle/0.001/2/pi;
        acqResults.codeFreq(PRN) = settings.codeFreqBasis + ...
            (acqResults.carrFreq(PRN) - settings.IF)/...
            settings.carrFreqBasis * settings.codeFreqBasis;

        %signal found, if IF =0 just change to 1 Hz to allow processing
        if(acqResults.carrFreq(PRN) == 0)
            acqResults.carrFreq(PRN) = 1;
        end
        acqResults.codePhase(PRN) = gather(codePhaseMax);

    else
        %--- No signal with this PRN --------------------------------------
        fprintf('. ');
    end   % if (peakSize/secondPeakSize) > settings.acqThreshold

end    % for PRN = satelliteList

%=== Acquisition is over ==================================================
fprintf(')\n');
