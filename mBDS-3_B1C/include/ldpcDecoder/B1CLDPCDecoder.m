function decodedNavBits = B1CLDPCDecoder(navBits_soft_I, navBits_soft_Q)
% B1CLDPCDecoder decodes one BeiDou B1C B-CNAV1 frame.
%
% Inputs:
%   navBits_soft_I : 1800-symbol B-CNAV1 soft sequence from data channel
%   navBits_soft_Q : 1800-symbol Q-branch soft/noise reference sequence
%
% Output:
%   decodedNavBits : 878-bit decoded B-CNAV1 navigation message
%
% B-CNAV1 frame structure:
%   SF1:
%       1  - 21   : BCH(21,6)
%       22 - 72   : BCH(51,8)
%
%   SF2/SF3:
%       73 - 1800 : interleaved LDPC-coded symbols
%
%   Decoded output:
%       decodedNavBits(1:6)     : SF1 part 1
%       decodedNavBits(7:14)    : SF1 part 2
%       decodedNavBits(15:614)  : SF2, 600 bits
%       decodedNavBits(615:878) : SF3, 264 bits

% -------------------------------------------------------------------------
decodedNavBits = [];

navBits_soft_I = navBits_soft_I(:).';
navBits_soft_Q = navBits_soft_Q(:).';

if length(navBits_soft_I) < 1800
    return;
end

if length(navBits_soft_Q) < 1800
    return;
end

navBits_soft_I = navBits_soft_I(1:1800);
navBits_soft_Q = navBits_soft_Q(1:1800);

% -------------------------------------------------------------------------
% Hard decision for BCH decoding and polarity determination
% -------------------------------------------------------------------------
bits_hard = zeros(size(navBits_soft_I));
bits_hard(navBits_soft_I > 0)  = 1;
bits_hard(navBits_soft_I <= 0) = 0;

polarity = 1;

% ==================== [SF1 BCH21_6 Decoding] =============================
checkBits = 1 - 2 * bits_hard(1:21);
[flag, decodedBits] = BCH21_6Decoding(checkBits);

% If BCH21_6 fails, try inverse polarity
if flag == 0
    bits_hard = 1 - bits_hard;
    polarity = -1;

    checkBits = 1 - 2 * bits_hard(1:21);
    [flag, decodedBits] = BCH21_6Decoding(checkBits);

    if flag == 0
        decodedNavBits = [];
        return;
    end
end

% Apply polarity correction to soft information
navBits_soft_I = navBits_soft_I * polarity;
navBits_soft_Q = navBits_soft_Q * polarity;

decodedNavBits = zeros(878, 1);
decodedNavBits(1:6) = decodedBits(:);

% ==================== [SF1 BCH51_8 Decoding] =============================
checkBits = 1 - 2 * bits_hard(22:72);
[flag, decodedBits] = BCH51_8Decoding(checkBits);

if flag == 0
    decodedNavBits = [];
    return;
end

decodedNavBits(7:14) = decodedBits(:);

% ==================== [Deinterleaving for SF2 and SF3] ===================
% Original B1C code:
%   temp_Bits_soft = reshape(bits_soft(73:end), [36, 48]);
%   Frame3Colum = 3:3:35;
%   Frame2Colum = setdiff(1:36, Frame3Colum);
%
% SF2 length: 1200 bits = 200 GF(64) symbols
% SF3 length: 528 bits  = 88 GF(64) symbols

tempBits_I = reshape(navBits_soft_I(73:end), [36, 48]);
tempBits_Q = reshape(navBits_soft_Q(73:end), [36, 48]);

Frame3Colum = 3:3:35;
Frame2Colum = setdiff(1:36, Frame3Colum);

Frame2_soft_I = reshape(tempBits_I(Frame2Colum, :)', 1, 1200);
Frame2_soft_Q = reshape(tempBits_Q(Frame2Colum, :)', 1, 1200);

Frame3_soft_I = reshape(tempBits_I(Frame3Colum, :)', 1, 528);
Frame3_soft_Q = reshape(tempBits_Q(Frame3Colum, :)', 1, 528);

% Keep the deinterleaved systematic hard decisions as a safe fallback. The
% LDPC code is systematic, so the first information bits must remain usable
% when EMS does not converge on a valid codeword.
Frame2_hard = Frame2_soft_I > 0;
Frame3_hard = Frame3_soft_I > 0;

% ==================== [SF2 LLR Soft Information Generation] ===============

[normIp_SF2, normSigma2_SF2] = noiseEstSF2(Frame2_soft_I, Frame2_soft_Q);
Symbol_LLR_SF2 = initLlrSF2(normIp_SF2, normSigma2_SF2);

% Optional dimension correction/check
if size(Symbol_LLR_SF2, 1) ~= 200 && size(Symbol_LLR_SF2, 2) == 200
    Symbol_LLR_SF2 = Symbol_LLR_SF2.';
end

if size(Symbol_LLR_SF2, 1) ~= 200
    error('B1CLDPCDecoder:LLRSizeError', ...
          'For B1C SF2, Symbol_LLR_SF2 should be a 200 x 64 matrix.');
end

% ==================== [SF2 EMS Decoding] =================================
decoded_syms_SF2 = B1C_EMS_decode_SF2(Symbol_LLR_SF2);
decoded_syms_SF2 = decoded_syms_SF2(:);
[~, ~, H_SF2] = getH_informationSF2();
sf2Syndrome = gf(decoded_syms_SF2.', 6) * gf(H_SF2.', 6);
sf2DecodeOK = all(sf2Syndrome == 0);

% First 100 GF(64) symbols are information symbols
if sf2DecodeOK
    decoded_syms_SF2_info = decoded_syms_SF2(1:100);
    bin_mat_SF2 = dec2bin(decoded_syms_SF2_info, 6) - '0';
    decodedNav_SF2 = reshape(bin_mat_SF2', 600, 1);
else
    decodedNav_SF2 = double(Frame2_hard(1:600).');
end

% ==================== [SF3 LLR Soft Information Generation] ===============
[normIp_SF3, normSigma2_SF3] = noiseEstSF3(Frame3_soft_I, Frame3_soft_Q);
Symbol_LLR_SF3 = initLlrSF3(normIp_SF3, normSigma2_SF3);

% Optional dimension correction/check
if size(Symbol_LLR_SF3, 1) ~= 88 && size(Symbol_LLR_SF3, 2) == 88
    Symbol_LLR_SF3 = Symbol_LLR_SF3.';
end

if size(Symbol_LLR_SF3, 1) ~= 88
    error('B1CLDPCDecoder:LLRSizeError', ...
          'For B1C SF3, Symbol_LLR_SF3 should be an 88 x 64 matrix.');
end

% ==================== [SF3 EMS Decoding] =================================
decoded_syms_SF3 = B1C_EMS_decode_SF3(Symbol_LLR_SF3);
decoded_syms_SF3 = decoded_syms_SF3(:);
[~, ~, H_SF3] = getH_informationSF3();
sf3Syndrome = gf(decoded_syms_SF3.', 6) * gf(H_SF3.', 6);
sf3DecodeOK = all(sf3Syndrome == 0);

% First 44 GF(64) symbols are information symbols
if sf3DecodeOK
    decoded_syms_SF3_info = decoded_syms_SF3(1:44);
    bin_mat_SF3 = dec2bin(decoded_syms_SF3_info, 6) - '0';
    decodedNav_SF3 = reshape(bin_mat_SF3', 264, 1);
else
    decodedNav_SF3 = double(Frame3_hard(1:264).');
end

% ==================== [Output Decoded B-CNAV1 Message] ===================
decodedNavBits(15:614)  = decodedNav_SF2;
decodedNavBits(615:878) = decodedNav_SF3;

end
