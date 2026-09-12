function [eph, firstSubFrame, TOD] = NAVdecoding(I_P_InputBits)
% findPreambles finds the first preamble occurrence in the bit stream of
% each channel. The preamble is verified by check of the spacing between
% preambles (6sec) and parity checking of the first two words in a
% subframe. At the same time function returns list of channels, that are in
% tracking state and with valid preambles in the nav data stream.
%
%[eph, firstSubFrame,TOD] = CNAVdecoding(I_P_InputBits)
%
%   Inputs:
%       I_P_InputBits   - output from the tracking function
%
%   Outputs:
%       firstSubframe   - Starting positions of the first message in the 
%                       input bit stream I_P_InputBits in each channel. 
%                       The position is CNAV bit(20ms before convolutional decoding) 
%                       count since start of tracking. Corresponding value will
%                       be set to inf if no valid preambles were detected in
%                       the channel.
%       TOD             - Time Of Week (TOW) of the first message(in seconds).
%                       Corresponding value will be set to inf if no valid preambles
%                       were detected in the channel.
%       eph             - SV ephemeris. 

%--------------------------------------------------------------------------
%                         CU Multi-GNSS SDR  
% (C) Written by Yafeng Li, Jakob Almqvist, Nagaraj C. Shivaramaiah and Dennis M. Akos
%--------------------------------------------------------------------------
%
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
%--- Initialize ephemeris structure  --------------------------------------
% This is in order to make sure variable 'eph' for each SV has a similar 
% structure when only one or even none of the three requisite messages
% is decoded for a given PRN.
    try
        eph = eph_structure_init(); 
    catch
        eph = struct(); 
    end
    firstSubFrame = inf;
    TOD = inf;
    got_types = [false, false, false]; 

    %% Secondary Code Stripping (De-Barker) ===================================
    % The L3OCd signal has a 5-chip Barker code secondary overlay. We must 
    % detect the correct phase and strip this code to recover the underlying 
    % symbol stream before Viterbi decoding.
    BarkerPattern = [-1 -1 -1  1 -1]; 
    BcLen = 5;
    maxEnergy = 0;
    bestPhase = 1;
    SymbolStream = [];
    %--- Search for best Barker phase -------------------------------------
    for phase = 1 : BcLen
        currentStream = I_P_InputBits(phase:end);
        nSym = floor(length(currentStream) / BcLen);
        if nSym < 600, continue; end
        reshaped = reshape(currentStream(1 : nSym*BcLen), BcLen, nSym);
        collapsedSymbols = BarkerPattern * reshaped;
        % Calculate energy to find the optimal phase
        currentEnergy = sum(abs(collapsedSymbols));
        if currentEnergy > maxEnergy
            maxEnergy = currentEnergy;
            bestPhase = phase;
            SymbolStream = collapsedSymbols; 
        end
    end
    
    if isempty(SymbolStream), return; end
    % Remove DC offset from the symbol stream
    SymbolStream = SymbolStream - mean(SymbolStream);

   %% Viterbi Decoding ====================================================
   %--- Take even number of input bits to do decoing ---------------------
    evenLen = length(SymbolStream) - rem(length(SymbolStream),2);
    encodedSymbols = SymbolStream(1:evenLen);
    dataBitsInput = (encodedSymbols < 0); 
    % Convert code polynomials to trellis description for rate 1/2 code
    trellis = poly2trellis(7, [133 171]);
    tblen = 105; 

    %--- Generate the preamble pattern ------------------------------------
    preamble_bin  = [0 0 0 0 0 1 0 0 1 0 0 1 0 1 0 0 1 1 1 0]; 
    preamble_corr = 1 - 2 * preamble_bin; 
    frameLenBits = 300; 
    %% NAV data decoding =====================================================
    % The first bit in dataBits may be G1 or G2 ouput, the bit
    % stream to be decoded must start at bit of G1 output. So we must
    % search the first two bits to find the right one corresponding to G1.
    for G1orG2 = 1:2
        decodedBits = vitdec(dataBitsInput(G1orG2:end-(G1orG2-1)), ...
                             trellis, tblen, 'trunc', 'hard');
        antipodalDecoded = 1 - 2 * decodedBits;
        
        tlmXcorrResult = xcorr(antipodalDecoded, preamble_corr);
        % Find all starting points of all preamble-like patterns
        xcorrCenter = (length(tlmXcorrResult) + 1) / 2;
        %--- Find at what index the preambles start -----------------------
        % Threshold is set to 18 (allowing up to 2 bit errors in preamble)
        index = find(abs(tlmXcorrResult(xcorrCenter : end)) >= 18)';
        
        for i = 1:length(index)
            startIdx = index(i);
            if (startIdx + frameLenBits - 1) > length(decodedBits), continue; end
            
            navBits = decodedBits(startIdx : startIdx + frameLenBits - 1);
            
            % Polarity Correction
            if sum(xor(navBits(1:20), preamble_bin)) > 10 
                 navBits = 1 - navBits;
            end
            if sum(xor(navBits(1:20), preamble_bin)) > 2, continue; end

           % The generation polynomial coefficients for CRC-24Q
             polyvec = [24 23 18 17 14 11 10 7 6 5 4 3 1 0];
             crcDet = comm.CRCDetector(polyvec);
             % Detect errors in input data using CRC
             [~, frmError] = step(crcDet, navBits(:));
            
             if (~frmError)
                typeStr = num2str(navBits(21:26));
                msgType = bin2dec(typeStr);
                % Filter for Ephemeris Types (10, 11, 12)
                if ~ismember(msgType, [10, 11, 12])
                     fprintf('  > Skipped Frame: Type %d (Not Ephemeris)\n', msgType);
                    continue; 
                end
                %---Save for first message --------------------------
                if isinf(firstSubFrame)
                    firstSubFrame = (startIdx - 1) * 2 * 5 + bestPhase + (G1orG2 - 1) * 5;
                    tsStr = num2str(navBits(27:41));
                    tsVal = bin2dec(tsStr);
                    TOD = tsVal * 3; 
                    eph.TOD = TOD;
                    eph.TOW = TOD; 
                end
                %--- Ephemeris decoding -----------------------------------
                try
                    [eph, ~] = ephemeris(navBits, eph);
                    % Track which message types we have successfully collected
                    if msgType == 10, got_types(1) = true; end
                    if msgType == 11, got_types(2) = true; end
                    if msgType == 12, got_types(3) = true; end
                    
                    fprintf('Decoded L3OCd Frame: Type %d, SV_ID %d\n', msgType, eph.SV_ID);
                    
                catch ME
                    warning(ME.identifier, 'Ephemeris parsing failed: %s', ME.message);
                end
                
                
                if all(got_types)
                    % disp('All 3 ephemeris strings (10, 11, 12) collected!');
                    return; 
                end
             end
        end
        % If all types collected, exit the G1/G2 loop as well
        if all(got_types), return; end
    end
end