function [normIp, normSigma2] = noiseEst(I_P, Q_P)
%NOISEEST Estimate normalized prompt-I samples and noise variance.
%   Uses the M2M4 moment estimator on I/Q prompt correlator outputs. The
%   returned normIp and normSigma2 are intended for soft-decision LLR
%   initialization.

tiny = 1e-12;

I = I_P(:);
Q = Q_P(:);

Z     = I.^2 + Q.^2;
meanZ = mean(Z);
varZ  = var(Z, 1);
qVar  = var(Q, 1);

if qVar > tiny
    Pav = sqrt(max(meanZ^2 - varZ, 0));
    estNoiseVar = 0.5 * max(meanZ - Pav, 0);
else
    Pav = sqrt(max(meanZ^2 - varZ / 2, 0));
    estNoiseVar = max(meanZ - Pav, 0);
end

estSignalAmp = sqrt(max(Pav, 0));
estSignalAmp = max(estSignalAmp, tiny);

normIp = I_P / estSignalAmp;
normSigma2 = estNoiseVar / (estSignalAmp^2);
normSigma2 = max(normSigma2, tiny);
end
