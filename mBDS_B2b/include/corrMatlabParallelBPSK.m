function correValues = corrMatlabParallelBPSK(settings, rawSignal, caCodeTable, ...
    remCarrPhase, carrPhaseStep, remCodePhase, codePhaseStep, ...
    startIdx, chSampSize, isDataRead)
%CORRBPSKMATLABPARALLEL Matlab implementation of the tracking correlator.
%   This function keeps the same calling interface as the SIMD / GPU
%   correlators so it can be used directly in trkChannelsParallel.m when
%   settings.correlatorType == 0.
%
%   Input/output layout follows the MEX correlators:
%       correValues = [I_E; Q_E; I_P; Q_P; I_L; Q_L]
%   where each column corresponds to one tracking channel.

persistent rawSignalI rawSignalQ rawSignalT

% For complex data, split the interleaved buffer [I0 Q0 I1 Q1 ...] into I
% and Q row vectors once for each newly read data block. If the current
% call reuses the same rawSignal block, keep the cached I/Q buffers.
if (settings.fileType == 1) && (isDataRead || isempty(rawSignalT))
    rawSignalT = rawSignal';
elseif (settings.fileType == 2) && (isDataRead || isempty(rawSignalI) || isempty(rawSignalQ))
    rawSignalI = rawSignal(1:2:end)';
    rawSignalQ = rawSignal(2:2:end)';
end

% Number of active tracking channels in this millisecond.
channelCnt = numel(startIdx);

% Output order per channel: [IE QE IP QP IL QL]^T.
correValues = zeros(6, channelCnt);

% Define early-late offset (in chips).
earlyLateSpc = settings.dllCorrelatorSpacing;

for channelNr = 1:channelCnt
    % Find the size of the current code period in whole samples for the
    % current channel, and the start index of this channel block within the
    % 1-second rawSignal buffer.
    blksize = chSampSize(channelNr);
    codeStartIdx = startIdx(channelNr);

    % Get a vector with the local code for the current channel. The table
    % already contains the guard chips required by ceil(tcode) + 1.
    caCode = caCodeTable(channelNr, :);

    % Extract the current signal block to be processed by this channel.
    if settings.fileType == 1
        rawSignalBlock = rawSignalT(codeStartIdx : codeStartIdx + blksize - 1);
    else
        rawSignalBlockI = rawSignalI(codeStartIdx : codeStartIdx + blksize - 1);
        rawSignalBlockQ = rawSignalQ(codeStartIdx : codeStartIdx + blksize - 1);
    end

    %% Set up all the code phase tracking information ---------------------
    sampleIndex = (0:blksize-1);
    codePhase = remCodePhase(channelNr) + codePhaseStep(channelNr) .* sampleIndex;

    % Define index into early code vector
    tcode2 = ceil(codePhase - earlyLateSpc) + 1;
    earlyCode = caCode(tcode2);

    % Define index into late code vector
    tcode2 = ceil(codePhase + earlyLateSpc) + 1;
    lateCode = caCode(tcode2);

    % Define index into prompt code vector
    tcode2 = ceil(codePhase) + 1;
    promptCode = caCode(tcode2);

    %% Generate the carrier frequency to mix the signal to baseband -------
    % carrPhaseStep is already the carrier phase step in radians per
    % sample, so the local carrier phase is formed directly with sample
    % indices instead of a time vector in seconds.
    trigarg = carrPhaseStep(channelNr) * sampleIndex + remCarrPhase(channelNr);
    carrCos = cos(trigarg);
    carrSin = sin(trigarg);

    %% Do correlation to generate the six standard accumulated values -----
    % First mix to baseband
    if settings.fileType == 1
        iBasebandSignal = rawSignalBlock .* carrCos;
        qBasebandSignal = -rawSignalBlock .* carrSin;
    else
        iBasebandSignal = rawSignalBlockI .* carrCos + rawSignalBlockQ .* carrSin;
        qBasebandSignal = rawSignalBlockQ .* carrCos - rawSignalBlockI .* carrSin;
    end

    % Now get early, late, and prompt values for each
    I_E = sum(earlyCode  .* iBasebandSignal);
    Q_E = sum(earlyCode  .* qBasebandSignal);
    I_P = sum(promptCode .* iBasebandSignal);
    Q_P = sum(promptCode .* qBasebandSignal);
    I_L = sum(lateCode   .* iBasebandSignal);
    Q_L = sum(lateCode   .* qBasebandSignal);

    % Output order matches the SIMD / GPU correlator interface used by
    % trkChannelsParallel.m.
    correValues(:, channelNr) = [I_E; Q_E; I_P; Q_P; I_L; Q_L];
end
