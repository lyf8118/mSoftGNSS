function decodedNavBits = B2aLDPCDecoder(navBits_soft_I, navBits_soft_Q)
% B2aLDPCDecoder decodes one B-CNAV2 frame of BeiDou B2a signal.
%
% Inputs:
%   navBits_soft_I : 600-bit B-CNAV2 soft sequence after secondary-code despreading
%   navBits_soft_Q : 600-bit Q-branch soft/noise reference sequence
%
% Output:
%   decodedNavBits : 288 decoded navigation bits for CRC and ephemeris decoding


% ==================== [LLR Soft Information Matrix Generation] ===========
% B-CNAV2 frame length before LDPC decoding: 600 bits
%   1  - 24   : preamble
%   25 - 600  : LDPC-coded payload, 576 bits = 96 GF(64) symbols
%
% Therefore, the LDPC payload section starts from the 25th bit.

payload_I = navBits_soft_I(25:end);
payload_Q = navBits_soft_Q(25:end);

% Make sure the payload length is exactly 576 bits
payload_I = payload_I(1:576);
payload_Q = payload_Q(1:576);

% Blind SNR estimation and soft information normalization
[normIp, normSigma2] = noiseEst(payload_I, payload_Q);

% Generate the LLR cost matrix required for non-binary LDPC decoding
% For B2a B-CNAV2:
%   576 coded bits = 96 GF(64) symbols
% Therefore, Symbol_LLR should be a 96 x 64 matrix.
Symbol_LLR = initLlr(normIp, normSigma2);


% ==================== [EMS Decoding] =====================================
decoded_symbols_96 = B2a_EMS_decode(Symbol_LLR);
decoded_symbols_96 = decoded_symbols_96(:);

% B-CNAV2 information part:
%   first 48 GF(64) symbols = 48 * 6 = 288 information bits
decoded_symbols_48 = decoded_symbols_96(1:48);

% ==================== [Symbol Inverse Mapping to Bitstream] ==============
bin_mat_48x6 = dec2bin(decoded_symbols_48, 6) - '0';

decodedNavBits = reshape(bin_mat_48x6', 288, 1);

end
