function [normIp, normSigma2] = noiseEst(I_P, Q_P)
% =========================================================================
% Function Name: noiseEst
% Description:
%   Performs blind Signal-to-Noise Ratio (SNR) estimation using the
%   Moment Method (M2M4) based on the squared magnitude (Z = I^2 + Q^2),
%   and normalizes the prompt correlation values for soft-decision decoding.
% =========================================================================

%% BDS-3 B2b LDPC   
colH = 162;           
bitsPerSymbol = 6;

tiny = 1e-12;
I = I_P(:);
Q = Q_P(:);

Z     = I.^2 + Q.^2;
meanZ = mean(Z);
varZ  = var(Z, 1);          % population variance, matches moment definition
qVar  = var(Q, 1);

if qVar > tiny
    % Complex-noise / IQ-noise case:
    % meanZ = A^2 + 2*sigma^2,  varZ = 4*A^2*sigma^2 + 4*sigma^4
    Pav = sqrt(max(meanZ^2 - varZ, 0));          % Pav = A^2
    estNoiseVar  = 0.5 * max(meanZ - Pav, 0);    % sigma^2 (single-branch)
else
    % Q degenerates to 0 (real case):
    % meanZ = A^2 + sigma^2,  varZ = 4*A^2*sigma^2 + 2*sigma^4
    Pav = sqrt(max(meanZ^2 - varZ/2, 0));        % Pav = A^2
    estNoiseVar  = 1.0 * max(meanZ - Pav, 0);    % sigma^2 (single-branch)
end

estSignalAmp = sqrt(max(Pav, 0));                % A
estSignalAmp = max(estSignalAmp, tiny);          % avoid divide-by-zero

%% Normalize for decoder===============================================
normIp = I_P / estSignalAmp;
normSigma2  = estNoiseVar / (estSignalAmp^2);

%% Format Output for Decoder===========================================
% Reshape from 1x972 to 162x6 (Rows = symbols, Cols = bits)
normIp = reshape(normIp, bitsPerSymbol, colH).';
end
