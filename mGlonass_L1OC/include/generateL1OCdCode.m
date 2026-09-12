function L1OcdCode = generateL1OCdCode(PRN, settings)
% generateL1OcdCode.m generates one GLONASS L1OCd PRN code.
%
% L1OcdCode = generateL1OcdCode(PRN, settings)
%
%   Inputs:
%       PRN       - GLONASS L1OC SV ID number, 0...63. PRN 0 is reserved
%                   in the ICD, but is kept here for test-vector checking.
%       settings  - receiver settings structure. It is kept for the same
%                   calling style as the L1OC code generator. This function
%                   does not require any field from settings.
%
%   Outputs:
%       L1OcdCode - return-zero/TDM-form code sequence. The non-zero chips
%                   are +1/-1 values generated from the L1OCd PRN; one zero
%                   is inserted after each PRN chip, consistent with the
%                   L1OC local-code table generation.
%
% Notes:
%   1) This function generates the L1OCd PRN replica only. The L1OCd data,
%      convolution encoder symbols and OC1 overlay are not applied here.
%   2) Raw L1OCd PRN: length 1023 chips, period 2 ms, chip rate 0.5115 MHz.
%--------------------------------------------------------------------------

%#ok<*INUSD>  % settings is intentionally unused, kept for interface style

%--- Check PRN/SV ID -------------------------------------------------------
if nargin < 1
    error('generateL1OcdCode:MissingInput', 'PRN/SV ID is required.');
end

if PRN ~= floor(PRN) || PRN < 0 || PRN > 63
    error('generateL1OcdCode:InvalidPRN', ...
        'GLONASS L1OCd PRN/SV ID must be an integer in the range 0...63.');
end

%--- ICD parameters --------------------------------------------------------
CodeLength = 1023;

% DC1: 10 stages, feedback taps 7 and 10, IS1 = 0011001000.
reg1 = [0 0 1 1 0 0 1 0 0 0];

% DC2: 10 stages, feedback taps 3, 7, 9 and 10, IS2 = PRN/SV ID.
% The least significant bit enters the last register stage, so dec2bin(PRN,10)
% can be loaded from left to right into stages 1...10.
reg2 = dec2bin(PRN, 10) - '0';

%--- Generate raw binary PRN chips ----------------------------------------
rawCode = zeros(1, CodeLength);
for index = 1:CodeLength
    % PRN output is modulo-2 sum of the last stages of DC1 and DC2.
    rawCode(index) = xor(reg1(10), reg2(10));

    fb1 = xor(reg1(7), reg1(10));
    fb2 = xor(xor(xor(reg2(3), reg2(7)), reg2(9)), reg2(10));

    % Shift direction is from lower trigger number to higher trigger number.
    reg1 = [fb1 reg1(1:9)];
    reg2 = [fb2 reg2(1:9)];
end

%--- Convert binary chips to bipolar chips: 0 -> +1, 1 -> -1 --------------
L1OcdCode = 1 - 2 * rawCode;
L1OcdCode = [L1OcdCode; zeros(1, length(L1OcdCode))];
L1OcdCode = reshape(L1OcdCode, 1, []);

end
