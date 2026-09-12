function correValues = corrMatlabSerialBPSK(settings, rawSignal, caCode, ...
    remCarrPhase, carrPhaseStep, remCodePhase, codePhaseStep)

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
earlyCode   = caCode(tcode2);

% Define index into late code vector
tcode       = (remCodePhase+earlyLateSpc) : ...
    codePhaseStep : ...
    ((blksize-1)*codePhaseStep+remCodePhase+earlyLateSpc);
tcode2      = ceil(tcode) + 1;
lateCode    = caCode(tcode2);

% Define index into prompt code vector
tcode       = remCodePhase : ...
    codePhaseStep : ...
    ((blksize-1)*codePhaseStep+remCodePhase);
tcode2      = ceil(tcode) + 1;
promptCode  = caCode(tcode2);


%% Generate the carrier frequency to mix the signal to baseband -----------
% Get the argument to sin/cos functions
% time    = (0:blksize) ./ settings.samplingFreq;
% trigarg = ((carrFreq * 2.0 * pi) .* time) + remCarrPhase;
trigarg = carrPhaseStep .* (0:blksize) + remCarrPhase;
% Remaining carrier phase for each tracking update

% Finally compute the signal to mix the collected data to
% bandband
carrsig = exp(-1i .* trigarg(1:blksize));

%% Do correlation to Generate the six standard accumulated values -----------
% First mix to baseband
iBasebandSignal = real(carrsig .* rawSignal);
qBasebandSignal = imag(carrsig .* rawSignal);

% Now get early, late, and prompt values for each
I_E = sum(earlyCode  .* iBasebandSignal);
Q_E = sum(earlyCode  .* qBasebandSignal);
I_P = sum(promptCode .* iBasebandSignal);
Q_P = sum(promptCode .* qBasebandSignal);
I_L = sum(lateCode   .* iBasebandSignal);
Q_L = sum(lateCode   .* qBasebandSignal);

correValues = [I_E,Q_E,I_P,Q_P,I_L,Q_L];
