function correValues = corrMatlabSerialQPSK(settings, rawSignal, L5CCodeTable, ...
    remCarrPhase, carrPhaseStep, remCodePhase, codePhaseStep)

codeLen = length(L5CCodeTable)/2;
L5CCodeD = L5CCodeTable(1:codeLen);
L5CCodeP = L5CCodeTable(1+codeLen:end);
earlyLateSpc = settings.dllCorrelatorSpacing;

rawSignal = rawSignal';

% For complex data
if (settings.fileType == 2)
    rawSignal1 = rawSignal(1:2:end);
    rawSignal2 = rawSignal(2:2:end);
    rawSignal = rawSignal1 + 1i .* rawSignal2;  % transpose vector
end

blksize = length(rawSignal);

% Define index into early code vector
tcode       = (remCodePhase-earlyLateSpc) : ...
    codePhaseStep : ...
    ((blksize-1)*codePhaseStep+remCodePhase-earlyLateSpc);
tcode2      = ceil(tcode) + 1;
earlyCode   = L5CCodeD(tcode2);

% For pilot channel signal tracking
if (settings.pilotTRKflag == 1)
    earlyCodeQ   = L5CCodeP(tcode2);
end
% Define index into late code vector
tcode       = (remCodePhase+earlyLateSpc) : ...
    codePhaseStep : ...
    ((blksize-1)*codePhaseStep+remCodePhase+earlyLateSpc);
tcode2      = ceil(tcode) + 1;
lateCode    = L5CCodeD(tcode2);

% For pilot channel signal tracking
if (settings.pilotTRKflag == 1)
    lateCodeQ   = L5CCodeP(tcode2);
end
% Define index into prompt code vector
tcode       = remCodePhase : ...
    codePhaseStep : ...
    ((blksize-1)*codePhaseStep+remCodePhase);
tcode2      = ceil(tcode) + 1;
promptCode  = L5CCodeD(tcode2);

% For pilot channel signal tracking
if (settings.pilotTRKflag == 1)
    promptCodeQ   = L5CCodeP(tcode2);
end

%% Generate the carrier frequency to mix the signal to baseband -----------
% Get the argument to sin/cos functions
trigarg = carrPhaseStep .* (0:blksize) + remCarrPhase;

% Finally compute the signal to mix the collected data to
% bandband
carrsig = exp(-1i .* trigarg(1:blksize));

%% Do correlation to Generate the six standard accumulated values -----------
% First mix to baseband
iBasebandSignal = real(carrsig .* rawSignal);
qBasebandSignal = imag(carrsig .* rawSignal);

% Now get early, late, and prompt values for data
I_E = sum(earlyCode  .* iBasebandSignal);
Q_E = sum(earlyCode  .* qBasebandSignal);
I_P = sum(promptCode .* iBasebandSignal);
Q_P = sum(promptCode .* qBasebandSignal);
I_L = sum(lateCode   .* iBasebandSignal);
Q_L = sum(lateCode   .* qBasebandSignal);

correValues = [I_E,Q_E,I_P,Q_P,I_L,Q_L];

% For pilot channel signal tracking
if (settings.pilotTRKflag == 1)
    pilot_I_E = sum(earlyCodeQ  .* iBasebandSignal);
    pilot_Q_E = sum(earlyCodeQ  .* qBasebandSignal);
    pilot_I_P = sum(promptCodeQ .* iBasebandSignal);
    pilot_Q_P = sum(promptCodeQ .* qBasebandSignal);
    pilot_I_L = sum(lateCodeQ   .* iBasebandSignal);
    pilot_Q_L = sum(lateCodeQ   .* qBasebandSignal);
    correValuesP = [pilot_I_E,pilot_Q_E,pilot_I_P,pilot_Q_P,pilot_I_L,pilot_Q_L];
    correValues = [correValues, correValuesP];
end








