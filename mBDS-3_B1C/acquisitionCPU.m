function acqResults = acquisitionCPU(fid, settings)
%Function performs cold start acquisition on the collected "data". It
%searches for BDS-3 B1C signals of all satellites, which are listed in field
%"acqSatelliteList" in the settings structure. Function saves code phase
%and frequency of the detected signals in the "acqResults" structure.
%
%acqResults = acquisitionCPU(fid, settings)
%
%   Inputs:
%       fid           - File identifier of the raw IF signal.
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
% Read data for acquisition: use 20 code length data to do fine
% acquisition; To ensure at least 20 complete B1C code periods are
% included, 21 code periods of data are needed
longSignal  = fread(fid, dataAdaptCoeff*21*samplesPerCode, settings.dataType)';
% Convert complex data into normal representation.
if (dataAdaptCoeff == 2)
    longSignal = longSignal(1:2:end) + 1i .* longSignal(2:2:end);
end

%% Acquisition initialization =============================================
% Find sampling period
ts = 1 / settings.samplingFreq;
% Find phase points of the local carrier wave
phasePoints = (0 : (samplesPerCode*2-1)) * 2 * pi * ts;

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
% Code resampling index used by B1C pilot-code fine acquisition.
% codeValueIndex maps each input sample to the corresponding local code chip.
tc = 1/(settings.codeFreqBasis*12);
codeValueIndex = floor(ts * (0: 20*samplesPerCode-1)/tc);
codeValueIndex = rem(codeValueIndex, settings.codeLength*12) + 1;

% Phase points of the local carrier wave
finePhasePoints = (0 : 20*samplesPerCode-1) * 2 * pi * ts;
% Perform search for all listed PRN numbers ...
fprintf('(');
for PRN = settings.acqSatelliteList
    %% Coarse acquisition ===========================================

    % Generate B1C data codes and sample them according to the sampling freq.
    DataPriTable = makeDataTable(settings,PRN);
    % Generate a zero-padded local code replica for correlation.
    localData = [DataPriTable(1:samplesPerCode), zeros(1,samplesPerCode)];

    % Perform DFT of B1C data code
    DataPriFreqDom = conj(fft(localData));

    % Use pilot signal power only when pilot-assisted acquisition is enabled.
    if settings.pilotACQflag
        PilotPriTable = makePilotTable(settings,PRN);
        localPilot = [PilotPriTable(1:samplesPerCode) ...
            zeros(1,samplesPerCode)];
        PilotPriFreqDom = conj(fft(localPilot));
    end

    codePhaseMax = 0; freqMax = 0; peakMax = 0;
    %--- Make the correlation for whole frequency band (for all freq. bins)
    for frqBinIndex = 1:numberOfFrqBins
        %--- Generate carrier wave frequency grid  -----------------------
        frqBins = settings.IF - settings.acqSearchBand + ...
            settings.acqStep * (frqBinIndex - 1);
        % Generate local sine and cosine
        sigCarr = exp(-1i * frqBins * phasePoints);
        % "Remove carrier" from the signal
        I1      = real(sigCarr .* longSignal(1:samplesPerCode*2));
        Q1      = imag(sigCarr .* longSignal(1:samplesPerCode*2));
        % Convert the baseband signal to frequency domain
        IQfreqDom = fft(I1 + 1i*Q1);
        % Multiplication in the frequency domain (correlation in time domain)
        convCodeIQ1 = IQfreqDom .* DataPriFreqDom;
        % Perform inverse DFT and store correlation results
        results = abs(ifft(convCodeIQ1));

        if settings.pilotACQflag
            % Pilot signal components; combine data and pilot power 1:3.
            convCodeIQ1 = IQfreqDom .* PilotPriFreqDom;
            results = (results + 3*abs(ifft(convCodeIQ1)))/4;
        end

        %--- Look for correlation peaks -----------------------------------
        % Find the max correlation peak and corresponding code phase
        [maxPeakTemp,maxIndexTemp] = max(results);
        if maxPeakTemp > peakMax
            peakMax = maxPeakTemp;
            codePhaseMax = maxIndexTemp;
            freqMax = frqBins;
        end

    end % frqBinIndex = 1:numberOfFrqBins

    %--- Find 1 chip wide B1C code phase exclude range around the peak ----
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

    % To prevent index from exceeding matrix dimensions in the fine
    % acquisition, move to previous code start position.
    if (codePhaseMax + 20*samplesPerCode -1) > length(longSignal)
        codePhaseMax = codePhaseMax - samplesPerCode;
    end

    % If the result is above threshold, then there is a signal ...
    if acqResults.peakMetric(PRN) > settings.acqThreshold
        %% Fine resolution frequency search =========================
        % Indicate PRN number of the detected signal
        fprintf('%02d ', PRN);

        % Use the same branch selected for coarse acquisition.
        if settings.pilotACQflag
            B1CCode = generatePilotBOC11(settings,PRN);
        else
            B1CCode = generateDataBOC11(settings,PRN);
        end
        B1CCode200ms = B1CCode(codeValueIndex);

        % Extract 10 B1C pilot code periods (200 ms) of incoming signal starting from the
        % detected B1C pilot code phase.
        sig200ms = longSignal(codePhaseMax:codePhaseMax + 20*samplesPerCode -1);

        % Coarse-frequency local carrier used before fine frequency estimate.
        localCarr200cm = exp(-1i * freqMax * finePhasePoints);

        %--- Integration for each 1 ms segment over 200 ms -----------------
        % Wipe off the selected B1C code and coarse carrier.
        basebandSig = sig200ms .* B1CCode200ms .* localCarr200cm;
        % Sum each 1 ms segment; the resulting 200-point sequence is used for
        % fine carrier estimation. Squaring removes the 180-degree data sign
        % transitions before the FFT-based frequency estimate.
        sumPerCode = sum(reshape(basebandSig,[],200));

        %--- Find the fine carrier freq. ----------------------------------
        % Index of the strongest spectral component after data-bit removal.
        [~,maxPowerIndex]  = max(abs(fft(sumPerCode.^2)));
        % Convert FFT-bin phase to carrier frequency correction. The division
        % by two compensates for the squaring operation above.
        shiftAngle = angle(exp(-2*pi*1i*(maxPowerIndex-1)/200))/2;
        acqResults.carrFreq(PRN) = freqMax - shiftAngle/0.001/2/pi;
        acqResults.codeFreq(PRN) = settings.codeFreqBasis + ...
            (acqResults.carrFreq(PRN) - settings.IF)/...
            settings.carrFreqBasis * settings.codeFreqBasis;

        %signal found, if IF =0 just change to 1 Hz to allow processing
        if(acqResults.carrFreq(PRN) == 0)
            acqResults.carrFreq(PRN) = 1;
        end
        acqResults.codePhase(PRN) = codePhaseMax;

    else
        %--- No signal with this PRN --------------------------------------
        fprintf('. ');
    end   % if (peakSize/secondPeakSize) > settings.acqThreshold

end    % for PRN = satelliteList

%=== Acquisition is over ==================================================
fprintf(')\n');
