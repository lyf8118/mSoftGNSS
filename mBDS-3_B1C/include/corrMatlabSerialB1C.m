function correValues = corrMatlabSerialB1C(settings, rawSignal, B1CCodeTable, ...
    remCarrPhase, carrPhaseStep, remCodePhase, codePhaseStep)

codeLen       = length(B1CCodeTable)/3;
B1CDataCode   = B1CCodeTable(1:codeLen);
B1CPilotCode  = B1CCodeTable(1+codeLen:codeLen*2);
pilotBOC61    = B1CCodeTable(1+codeLen*2:end);
% Define early-late offset (in chips)
earlyLateSpc = settings.dllCorrelatorSpacing;

rawSignal = rawSignal';

correValues = zeros(1,18);

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
earlyCode   = B1CDataCode(tcode2);
% For pilot channel signal tracking
pilotEarlyCode   = B1CPilotCode(tcode2);
p61_earlyCode    = pilotBOC61(tcode2);

% Define index into late code vector
tcode       = (remCodePhase+earlyLateSpc) : ...
    codePhaseStep : ...
    ((blksize-1)*codePhaseStep+remCodePhase+earlyLateSpc);
tcode2      = ceil(tcode) + 1;
lateCode    = B1CDataCode(tcode2);
% For pilot channel signal tracking
pilotLateCode = B1CPilotCode(tcode2);
p61_lateCode   = pilotBOC61(tcode2);

% Define index into prompt code vector
tcode       = remCodePhase : ...
    codePhaseStep : ...
    ((blksize-1)*codePhaseStep+remCodePhase);
tcode2      = ceil(tcode) + 1;
promptCode  = B1CDataCode(tcode2);
% For pilot channel signal tracking
pilotPromptCode = B1CPilotCode(tcode2);
p61_promptCode  = pilotBOC61(tcode2);

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

pilot_I_E = sum(pilotEarlyCode  .* iBasebandSignal);
pilot_Q_E = sum(pilotEarlyCode  .* qBasebandSignal);
pilot_I_P = sum(pilotPromptCode .* iBasebandSignal);
pilot_Q_P = sum(pilotPromptCode .* qBasebandSignal);
pilot_I_L = sum(pilotLateCode   .* iBasebandSignal);
pilot_Q_L = sum(pilotLateCode   .* qBasebandSignal);
correValues(7:12) = [pilot_I_E,pilot_Q_E,pilot_I_P,pilot_Q_P,pilot_I_L,pilot_Q_L];
if (settings.fullBandEn == 1)
    % Correlation values for pilot BOC(6,1) spreading waveform
    p61_I_E = sum(p61_earlyCode  .* iBasebandSignal);
    p61_Q_E = sum(p61_earlyCode  .* qBasebandSignal);
    p61_I_P = sum(p61_promptCode .* iBasebandSignal);
    p61_Q_P = sum(p61_promptCode .* qBasebandSignal);
    p61_I_L = sum(p61_lateCode   .* iBasebandSignal);
    p61_Q_L = sum(p61_lateCode   .* qBasebandSignal);
    correValues(13:18) = [p61_I_E, p61_Q_E, p61_I_P, p61_Q_P, p61_I_L, p61_Q_L];
end
