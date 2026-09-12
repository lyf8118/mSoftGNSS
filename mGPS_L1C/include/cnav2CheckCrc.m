function hasError = cnav2CheckCrc(bits)
%Checks the CRC-24Q of a decoded GPS L1C CNAV-2 subframe.
%
%The input contains the complete subframe, including its 24 transmitted CRC
%bits: 576 protected bits plus 24 CRC bits for Subframe 2, or 250 protected
%bits plus 24 CRC bits for Subframe 3. A zero remainder indicates that the
%CRC check has passed.
%
%hasError = cnav2CheckCrc(bits)
%
%   Inputs:
%       bits        - decoded CNAV-2 subframe bits represented by numeric
%                   values 0 or 1. Subframe 2 contains 600 bits and
%                   Subframe 3 contains 274 bits.
%
%   Outputs:
%       hasError    - CRC status flag. It is false if the CRC check passes
%                   and true if the CRC check fails.

%--- Initialize the CRC-24Q register and polynomial ----------------------
% The complete generator polynomial is 0x1864CFB. The leading x^24 term
% is implicit, so only its lower 24 bits are used in the CRC register.
crc = uint32(0);
poly = uint32(hex2dec('864CFB'));
mask = uint32(hex2dec('FFFFFF'));

%--- Compute the CRC remainder over the complete subframe ----------------
% Process the received bits in transmission order, most-significant bit
% first. The transmitted 24-bit CRC is included in this division.
for bitNr = 1:numel(bits)
    % Combine the input bit with the most-significant CRC-register bit.
    feedback = bitxor(bitget(crc, 24), uint32(bits(bitNr)));

    % Shift the CRC register and retain only its lower 24 bits.
    crc = bitand(bitshift(crc, 1), mask);

    % Apply the generator polynomial when the feedback bit is one.
    if feedback
        crc = bitxor(crc, poly);
    end
end

%--- A nonzero remainder indicates a CRC error ---------------------------
hasError = crc ~= 0;
end
