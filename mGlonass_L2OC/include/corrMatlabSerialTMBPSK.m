function correValues = corrMatlabSerialTMBPSK(settings, rawSignal, L2CCodeTable, ...
    remCarrPhase, carrPhaseStep, remCodePhaseL2C, codePhaseStep)

codeLen = length(L2CCodeTable)/2;
cmCode = L2CCodeTable(1:codeLen);
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
tcode       = (remCodePhaseL2C - earlyLateSpc) : ...
    codePhaseStep : ...
    ((blksize-1)*codePhaseStep+remCodePhaseL2C-earlyLateSpc);
tcode2      = ceil(tcode) + 1;
earlyCode   = cmCode(tcode2);



% Define index into late code vector
tcode       = (remCodePhaseL2C + earlyLateSpc) : ...
    codePhaseStep : ...
    ((blksize-1)*codePhaseStep + remCodePhaseL2C + earlyLateSpc);
tcode2      = ceil(tcode) + 1;
lateCode    = cmCode(tcode2);




% Define index into prompt code vector
tcode       = remCodePhaseL2C : ...
    codePhaseStep : ...
    ((blksize-1)*codePhaseStep + remCodePhaseL2C);
tcode2      = ceil(tcode) + 1;
promptCode  = cmCode(tcode2);

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

% Now get early, late, and prompt values for each
I_E = sum(earlyCode  .* iBasebandSignal);
Q_E = sum(earlyCode  .* qBasebandSignal);
I_P = sum(promptCode .* iBasebandSignal);
Q_P = sum(promptCode .* qBasebandSignal);
I_L = sum(lateCode   .* iBasebandSignal);
Q_L = sum(lateCode   .* qBasebandSignal);

correValues(1:6) = [I_E,Q_E,I_P,Q_P,I_L,Q_L];
