function [flag, decodedBits]= BCH51_8Decoding(checkBits)

% This function decodes the BCH(51,8) data symbols of the first subframe of
% in the CNAV-2 message. A total of 2^8-1 hypotheses are repeated to find
% the 8 original symbols.
%
% [flag,outputBits] = BCH51_8Decoding(bits)
%
%   Inputs:
%       bits         - Input row vector of message symbols to be decoded
%                      with 51 bits in bipolar (-1, +1)
%   Outputs:
%       flag         - Indicator for  decoding of the BCH(51,8) data 
%                      symbols, 1:success, 0: fail.
%       decodedBits  - Decoded message symbols for the BCH(51,8)

%--------------------------------------------------------------------------
%                           SoftGNSS v3.0
%
% Written by Yafeng Li
%--------------------------------------------------------------------------
%This program is free software; you can redistribute it and/or
%modify it under the terms of the GNU General Public License
%as published by the Free Software Foundation; either version 2
%of the License, or (at your option) any later version.
%
%This program is distributed in the hope that it will be useful,
%but WITHOUT ANY WARRANTY; without even the implied warranty of
%MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%GNU General Public License for more details.
%
%You should have received a copy of the GNU General Public License
%along with this program; if not, write to the Free Software
%Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
%USA.
%--------------------------------------------------------------------------

%CVS record:
%$Id: BCH51_8Decoding.m,v 1.1.2.5 2018/07/24 22:00:00 dpl Exp $

% --- Parameters ---
% Threshold can be adjusted based on noise level. 
% Max correlation is 52. 
% A threshold of 40 allows for approx 6 bit errors.
threshold = 40; 

% Ensure input is a row vector
checkBits = checkBits(:).';

% Check input length
if length(checkBits) ~= 52
    error('Input checkBits must be length 52 for L1C TOI decoding.');
end

% Total 2^8 = 256 hypotheses for the 8 LSBs
numHypotheses = 2^8;
hypoBits8 = zeros(numHypotheses, 8);
correValue = zeros(1, numHypotheses);

%% 1. Generate Hypotheses for 8 LSBs (Loop over all 256 possibilities)
% Reference: IS-GPS-800J Figure 3.2-4 / Text Section 3.2.3.2


for i = 0 : numHypotheses - 1
    idx = i + 1;
    
    % Get 8 bits for current hypothesis (b8...b1)
    % Using bitget to extract bits. b8 is MSB of this 8-bit chunk.
    % bitget(i, 8) is b8, bitget(i, 1) is b1.
    % The register is loaded with b8...b1.
    current8Bits = bitget(i, 8:-1:1); 
    hypoBits8(idx, :) = current8Bits;
    
    % --- Simulate LFSR (BCH 51,8 Generator) ---
    % Initial State: Loaded with the 8 data bits
    % Stage 1 (Left) ... Stage 8 (Right)
    % Diagram shows shift direction Left -> Right.
    % Register mapping: reg(1) corresponds to left-most, reg(8) to right-most.
    reg = fliplr(current8Bits); 
    
    generatedSeq = zeros(1, 51);
    
    for k = 1:51
        % Per image_8c9b55.jpg: "The last stage (8th) is selected as encoding output"
        outBit = reg(8);
        generatedSeq(k) = outBit;
        
        % Feedback Calculation (Polynomial 763 octal / 1+x+x^4+x^5+x^6+x^7+x^8)
        % Taps are at outputs of stages: 1, 4, 5, 6, 7, 8.
        % These are XORed (sum modulo 2) and fed back to Input (Stage 1).
        feedback = xor(reg(1), reg(4));
        feedback = xor(feedback, reg(5));
        feedback = xor(feedback, reg(6));
        feedback = xor(feedback, reg(7));
        feedback = xor(feedback, reg(8));
        
        % Shift Register: New bit enters at 1, others shift right
        reg = [feedback, reg(1:7)];
    end
    
    % --- Construct the Full 52-bit Hypothesis (Assuming b9 = 0) ---
    % Reference: image_8c9b55.jpg step (3) and (4)
    % The 9th bit (b9) is MSB of TOI.
    % Construction: [b9, (generatedSeq XOR b9)]
    % We assume b9 = 0 for the base hypothesis.
    % If b9 = 0: Sequence is [0, generatedSeq]
    
    hypoSequenceBin = [0, generatedSeq];
    
    % Convert to Bipolar for Correlation (+1 for 0, -1 for 1)
    hypoSyms = zeros(1, 52);
    hypoSyms(hypoSequenceBin == 0) = 1;
    hypoSyms(hypoSequenceBin == 1) = -1;
    
    % --- Correlation ---
    % Correlate received signal with the b9=0 hypothesis
    correValue(idx) = sum(hypoSyms .* checkBits);
end

%% 2. Decision Logic (Reference: image_8c9b39.jpg)
% "Find the maximum absolute value among the 256 correlation values."
[maxAbsVal, pos] = max(abs(correValue));

% Retrieve the raw correlation value (signed) to determine b9
rawCorr = correValue(pos);

% Check threshold
if maxAbsVal >= threshold
    flag = 1;
    
    % Extract the 8 LSBs corresponding to the best match
    best8Bits = hypoBits8(pos, :);
    
    % Determine the 9th bit (MSB / b9)
    % If max correlation is Positive: b9 was 0 (matches our hypothesis)
    % If max correlation is Negative: b9 was 1 (input was inverted)
    if rawCorr > 0
        b9 = 0;
    else
        b9 = 1;
    end
    
    % Combine to form full 9-bit TOI message
    decodedBits = [b9, best8Bits];
else
    flag = 0;
    decodedBits = [];
end

end
 