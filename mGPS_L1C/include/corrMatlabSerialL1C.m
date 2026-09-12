function correValues = corrMatlabSerialL1C(settings, rawSignal, L1CCodeTable, ...
    remCarrPhase, carrPhaseStep, remCodePhase, codePhaseStep)

codeLen = length(L1CCodeTable)/2;
L1CDataCode = L1CCodeTable(1:codeLen);
L1CPilotCode = L1CCodeTable(1+codeLen:end);
% Define early-late offset (in chips)
earlyLateSpc = settings.dllCorrelatorSpacing;

rawSignal = rawSignal';

correValues = zeros(1,12);

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
earlyCode   = L1CDataCode(tcode2);
% For pilot channel signal tracking
if (settings.pilotTRKflag == 1)
    pilotEarlyCode   = L1CPilotCode(tcode2);
end

% Define index into late code vector
tcode       = (remCodePhase+earlyLateSpc) : ...
    codePhaseStep : ...
    ((blksize-1)*codePhaseStep+remCodePhase+earlyLateSpc);
tcode2      = ceil(tcode) + 1;
lateCode    = L1CDataCode(tcode2);
% For pilot channel signal tracking
if (settings.pilotTRKflag == 1)
    pilotLateCode = L1CPilotCode(tcode2);
end

% Define index into prompt code vector
tcode       = remCodePhase : ...
    codePhaseStep : ...
    ((blksize-1)*codePhaseStep+remCodePhase);
tcode2      = ceil(tcode) + 1;
promptCode  = L1CDataCode(tcode2);
% For pilot channel signal tracking
if (settings.pilotTRKflag == 1)
    pilotPromptCode = L1CPilotCode(tcode2);
end
%% Generate the carrier frequency to mix the signal to baseband -----------
% Get the argument to sin/cos functions
trigarg = carrPhaseStep .* (0:blksize) + remCarrPhase;

% Finally compute the signal to mix the collected data to
% baseband
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

correValues(1:6) = [I_E,Q_E,I_P,Q_P,I_L,Q_L];

% For pilot channel signal tracking
if (settings.pilotTRKflag == 1)
    pilot_I_E = sum(pilotEarlyCode  .* iBasebandSignal);
    pilot_Q_E = sum(pilotEarlyCode  .* qBasebandSignal);
    pilot_I_P = sum(pilotPromptCode .* iBasebandSignal);
    pilot_Q_P = sum(pilotPromptCode .* qBasebandSignal);
    pilot_I_L = sum(pilotLateCode   .* iBasebandSignal);
    pilot_Q_L = sum(pilotLateCode   .* qBasebandSignal);
    correValues(7:12) = [pilot_I_E,pilot_Q_E,pilot_I_P,pilot_Q_P,pilot_I_L,pilot_Q_L];
end
