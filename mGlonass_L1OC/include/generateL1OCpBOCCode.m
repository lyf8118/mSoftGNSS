function L1OcpBOCCode = generateL1OCpBOCCode(PRN, settings)
% generateL1OCpBOCCode.m
% Generate GLONASS L1OCp PRN with MS(0101)/BOC(1,1) on fine TDM grid.
%
% Output fine-grid structure:
%   [0, 0, p, -p]
%
% Length:
%   raw L1OCp length = 4092
%   fine-grid length = 4092 * 4 = 16368

if nargin < 1
    error('generateL1OCpBOCCode:MissingInput', 'PRN/SV ID is required.');
end

if PRN ~= floor(PRN) || PRN < 0 || PRN > 63
    error('generateL1OCpBOCCode:InvalidPRN', ...
        'GLONASS L1OCp PRN/SV ID must be an integer in the range 0...63.');
end

%--- ICD parameters --------------------------------------------------------
CodeLength = 4092;

% DC1: 12 stages, feedback taps 6, 8, 11 and 12, IS1 = 000011000101.
reg1 = [0 0 0 0 1 1 0 0 0 1 0 1];

% DC2: 6 stages, feedback taps 1 and 6, IS2 = PRN/SV ID.
reg2 = dec2bin(PRN, 6) - '0';

%--- Generate raw binary PRN chips ----------------------------------------
rawCode = zeros(1, CodeLength);

for index = 1:CodeLength
    rawCode(index) = xor(reg1(12), reg2(6));

    fb1 = xor(xor(xor(reg1(6), reg1(8)), reg1(11)), reg1(12));
    fb2 = xor(reg2(1), reg2(6));

    reg1 = [fb1 reg1(1:11)];
    reg2 = [fb2 reg2(1:5)];
end

%--- Convert binary chips to bipolar chips: 0 -> +1, 1 -> -1 --------------
L1OcpCode = 1 - 2 * rawCode;

%--- L1OCp has MS(0101), forming BOC(1,1) -------------------------------
% Fine-grid structure:
%   [0, 0, p, -p]
tmp = zeros(1, 4 * length(L1OcpCode));

tmp(1:4:end) = 0;
tmp(2:4:end) = 0;
tmp(3:4:end) =  L1OcpCode;
tmp(4:4:end) = -L1OcpCode;

L1OcpBOCCode = tmp;

end