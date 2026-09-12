function correValues = corrMatlabParallelTMBPSK(settings, rawSignal, L2CodeTable, ...
    remCarrPhase, carrPhaseStep, remCodePhase, codePhaseStep, ...
    startIdx, chSampSize, isDataRead)
%CORRMATLABPARALLELTMBPSK MATLAB channel-parallel GLONASS L1OC tracking correlator.
%   This function keeps the same calling interface as the SIMD / GPU
%   correlators so it can be used directly in trkChannelsParallel.m when
%   settings.correlatorType == 0.
%
%   Input code-table layout for MATLAB mode:
%       L1OCCodeTable(channel,:) = [L1OCd code with guard chips, ...
%                                 L1OCp code with guard chips]
%   The first half is the repeated L1OCd data branch, and the second half is
%   the L1OCp pilot branch. Each branch already includes one wraparound chip on
%   both sides, so ceil(codePhase) + 1 can be used safely for E/P/L indexing.
%
%   Input/output layout follows the SIMD / GPU correlators:
%       correValues = [I_E; Q_E; I_P; Q_P; I_L; Q_L;
%                      pilot_I_E; pilot_Q_E; pilot_I_P;
%                      pilot_Q_P; pilot_I_L; pilot_Q_L]
%   where each column corresponds to one active tracking channel.

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

% Number of active tracking channels in this coherent integration.
channelCnt = numel(startIdx);

% Output order per channel: [IE QE IP QP IL QL pilotIE pilotQE ...]^T.
correValues = zeros(12, channelCnt);

% Define early-late offset (in chips).
earlyLateSpc = settings.dllCorrelatorSpacing;

for channelNr = 1:channelCnt
    % Find the size of the current L1OC code period in whole samples for the
    % current channel, and the start index of this channel block within the
    % shared rawSignal buffer.
    blksize = chSampSize(channelNr);
    codeStartIdx = startIdx(channelNr);

    % Get the local L1OC code table for the current channel. The first half
    % is the L1OCd data branch; the second half is the L1OCp pilot branch.
    l2cCode = L2CodeTable(channelNr, :);
    codeLen = length(l2cCode) / 2;
    cmCode = l2cCode(1:codeLen);
    if settings.pilotTRKflag == 1
        clCode = l2cCode(1+codeLen:end);
    end

    % Extract the current signal block to be processed by this channel.
    if settings.fileType == 1
        rawSignalBlock = rawSignalT(codeStartIdx : codeStartIdx + blksize - 1);
    else
        rawSignalBlockI = rawSignalI(codeStartIdx : codeStartIdx + blksize - 1);
        rawSignalBlockQ = rawSignalQ(codeStartIdx : codeStartIdx + blksize - 1);
    end

    %% Set up all the code phase tracking information ---------------------
    sampleIndex = 0:blksize-1;
    codePhase = remCodePhase(channelNr) + codePhaseStep(channelNr) .* sampleIndex;

    % Define index into early CM / L1OCp code vectors.
    tcode = ceil(codePhase - earlyLateSpc) + 1;
    earlyCode = cmCode(tcode);
    if settings.pilotTRKflag == 1
        earlyCodeCL = clCode(tcode);
    end

    % Define index into late CM / L1OCp code vectors.
    tcode = ceil(codePhase + earlyLateSpc) + 1;
    lateCode = cmCode(tcode);
    if settings.pilotTRKflag == 1
        lateCodeCL = clCode(tcode);
    end

    % Define index into prompt CM / L1OCp code vectors.
    tcode = ceil(codePhase) + 1;
    promptCode = cmCode(tcode);
    if settings.pilotTRKflag == 1
        promptCodeCL = clCode(tcode);
    end

    %% Generate the carrier frequency to mix the signal to baseband -------
    % carrPhaseStep is already the carrier phase step in radians per sample,
    % so the local carrier phase is formed directly with sample indices.
    trigarg = carrPhaseStep(channelNr) * sampleIndex + remCarrPhase(channelNr);
    carrCos = cos(trigarg);
    carrSin = sin(trigarg);

    %% Do correlation to generate the twelve standard accumulated values ---
    % First mix to baseband. This is equivalent to multiplying the raw signal
    % by exp(-j*trigarg), but keeps the real and complex input paths explicit.
    if settings.fileType == 1
        iBasebandSignal = rawSignalBlock .* carrCos;
        qBasebandSignal = -rawSignalBlock .* carrSin;
    else
        iBasebandSignal = rawSignalBlockI .* carrCos + rawSignalBlockQ .* carrSin;
        qBasebandSignal = rawSignalBlockQ .* carrCos - rawSignalBlockI .* carrSin;
    end

    % Now get early, late, and prompt values for the L1OCd data branch.
    I_E = sum(earlyCode  .* iBasebandSignal);
    Q_E = sum(earlyCode  .* qBasebandSignal);
    I_P = sum(promptCode .* iBasebandSignal);
    Q_P = sum(promptCode .* qBasebandSignal);
    I_L = sum(lateCode   .* iBasebandSignal);
    Q_L = sum(lateCode   .* qBasebandSignal);

    % Output order matches the SIMD / GPU correlator interface used by
    % trkChannelsParallel.m.
    correValues(1:6, channelNr) = [I_E; Q_E; I_P; Q_P; I_L; Q_L];

    if settings.pilotTRKflag == 1
        % Now get early, late, and prompt values for the L1OCp pilot branch.
        I_ECL = sum(earlyCodeCL  .* iBasebandSignal);
        Q_ECL = sum(earlyCodeCL  .* qBasebandSignal);
        I_PCL = sum(promptCodeCL .* iBasebandSignal);
        Q_PCL = sum(promptCodeCL .* qBasebandSignal);
        I_LCL = sum(lateCodeCL   .* iBasebandSignal);
        Q_LCL = sum(lateCodeCL   .* qBasebandSignal);
        correValues(7:12, channelNr) = [I_ECL; Q_ECL; I_PCL; Q_PCL; I_LCL; Q_LCL];
    end
end
