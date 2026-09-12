function Secondary = generate2ndCode(PRN)
% Generate the GPS L1C pilot secondary code in bipolar format (-1, +1).
% This implementation supports PRNs 1 to 63, which use the single-LFSR form.
%
% Secondary = generate2ndCode(PRN)
%
%   Inputs:
%       PRN        - Satellite PRN number.
%
%   Outputs:
%       Secondary  - 1800-chip L1C secondary-code sequence.

%--------------------------------------------------------------------------
%                         CU Multi-GNSS SDR
% Original authors: Yafeng Li, Nagaraj C. Shivaramaiah and Dennis M. Akos
% Reference: Adapted within the CU Multi-GNSS SDR receiver framework for GPS L1C.
%--------------------------------------------------------------------------

%% Check the supported PRN range.
if PRN < 1 || PRN > 63
    error('This function supports PRNs 1 to 63 only; PRN 64+ requires the dual-LFSR structure.');
end

%% L1C secondary-code parameter table.
% Each PRN maps to {S1 polynomial, initial state}, both stored as octal strings.
% Parameters are from IS-GPS-800 Table 3.2-3.
data_map = containers.Map('KeyType','int32','ValueType','any');

% PRNs 1 to 21.
data_map(1) = {'5111', '3266'};  data_map(2) = {'5421', '2040'};
data_map(3) = {'5501', '1527'};  data_map(4) = {'5403', '3307'};
data_map(5) = {'6417', '3756'};  data_map(6) = {'6141', '3026'};
data_map(7) = {'6351', '0562'};  data_map(8) = {'6501', '0420'};
data_map(9) = {'6205', '3415'};  data_map(10)= {'6235', '0337'};
data_map(11)= {'7751', '0265'};  data_map(12)= {'6623', '1230'};
data_map(13)= {'6733', '2204'};  data_map(14)= {'7627', '1440'};
data_map(15)= {'5667', '2412'};  data_map(16)= {'5051', '3516'};
data_map(17)= {'7665', '2761'};  data_map(18)= {'6325', '3750'};
data_map(19)= {'4365', '2701'};  data_map(20)= {'4745', '1206'};
data_map(21)= {'7633', '1544'};

% PRNs 22 to 42.
data_map(22)= {'6747', '1774'};  data_map(23)= {'4475', '0546'};
data_map(24)= {'4225', '2213'};  data_map(25)= {'7063', '3707'};
data_map(26)= {'4423', '2051'};  data_map(27)= {'6651', '3650'};
data_map(28)= {'4161', '1777'};  data_map(29)= {'7237', '3203'};
data_map(30)= {'4473', '1762'};  data_map(31)= {'5477', '2100'};
data_map(32)= {'6163', '0571'};  data_map(33)= {'7223', '3710'};
data_map(34)= {'6323', '3535'};  data_map(35)= {'7125', '3110'};
data_map(36)= {'7035', '1426'};  data_map(37)= {'4341', '0255'};
data_map(38)= {'4353', '0321'};  data_map(39)= {'4107', '3124'};
data_map(40)= {'5735', '0572'};  data_map(41)= {'6741', '1736'};
data_map(42)= {'7071', '3306'};

% PRNs 43 to 63.
data_map(43)= {'4563', '1307'};  data_map(44)= {'5755', '3763'};
data_map(45)= {'6127', '1604'};  data_map(46)= {'4671', '1021'};
data_map(47)= {'4511', '2624'};  data_map(48)= {'4533', '0406'};
data_map(49)= {'5357', '0114'};  data_map(50)= {'5607', '0077'};
data_map(51)= {'6673', '3477'};  data_map(52)= {'6153', '1000'};
data_map(53)= {'7565', '3460'};  data_map(54)= {'7107', '2607'};
data_map(55)= {'6211', '2057'};  data_map(56)= {'4321', '3467'};
data_map(57)= {'7201', '0706'};  data_map(58)= {'4451', '2032'};
data_map(59)= {'5411', '1464'};  data_map(60)= {'5141', '0520'};
data_map(61)= {'7041', '1766'};  data_map(62)= {'6637', '3270'};
data_map(63)= {'4577', '0341'};

%% Decode the LFSR parameters for the requested PRN.
params = data_map(PRN);
poly_oct = params{1};
init_oct = params{2};

% The polynomial is represented with 12 bits. Drop the two fixed end terms
% to obtain the 10 feedback taps aligned with state(2:11).
full_poly_bits = oct2bin_vec(poly_oct, 12);
taps = full_poly_bits(2:11);

% Per the ICD note, drop the leading zero bit from the 12-bit initial state.
full_init_bits = oct2bin_vec(init_oct, 12);
state = full_init_bits(2:12);

%% Generate the 1800-chip secondary sequence.
code_len = 1800;
binary_seq = zeros(1, code_len);

for k = 1:code_len
    output_bit = state(1);
    binary_seq(k) = output_bit;

    feedback_val = output_bit;
    intermediate_feedback = bitand(state(2:11), taps);
    feedback_val = xor(feedback_val, mod(sum(intermediate_feedback), 2));

    state(1:10) = state(2:11);
    state(11) = feedback_val;
end

% Convert to bipolar format: 0 -> +1, 1 -> -1.
Secondary = 1 - 2 * binary_seq;

end

function bin = oct2bin_vec(inputVal, numBits)
% Convert an octal string or number into a numeric binary vector.

if ischar(inputVal) || isstring(inputVal)
    octStr = char(inputVal);
else
    octStr = num2str(inputVal, '%04d');
end

val = base2dec(octStr, 8);
binStr = dec2bin(val, numBits);

bin = zeros(1, numBits);
for k = 1:numBits
    bin(k) = str2double(binStr(k));
end
end
