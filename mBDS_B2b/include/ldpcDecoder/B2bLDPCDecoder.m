function decodedNavBits = B2bLDPCDecoder(navBits_soft_I, navBits_soft_Q)
% ==================== [LLR Soft Information Matrix Generation] ===========
% 1. Extract soft information sequence of the LDPC payload section (from 
% the 29th bit to the end, totaling 972 floating-point numbers)
payload_I = navBits_soft_I(29:end);
payload_Q = navBits_soft_Q(29:end);

% 2. Blind SNR estimation and soft information normalization (M2M4 algorithm)
[normIp, normSigma2] = noiseEst(payload_I, payload_Q);

% 3. Generate the LLR cost matrix required for non-binary LDPC decoding [162 x 64]
Symbol_LLR = initLlr(normIp, normSigma2);

% ==================== [EMS Decoding] ====================
decoded_symbols_162 =  B2b_EMS_decode(Symbol_LLR);
decoded_symbols_162 = decoded_symbols_162(:);
decoded_symbols_81  = decoded_symbols_162(1:81);

% Symbol inverse mapping to bitstream
bin_mat_81x6 = dec2bin(decoded_symbols_81, 6) - '0';
decodedNavBits = reshape(bin_mat_81x6', 486, 1);