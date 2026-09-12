function [eph, firstSubFrame,TOW] = NAVdecoding(I_P)

% findPreambles finds the first preamble occurrence in the bit stream of
% each channel. The preamble is verified by check of the spacing between
% preambles (6sec) and parity checking of the first two words in a
% subframe. At the same time function returns list of channels, that are in
% tracking state and with valid preambles in the nav data stream.
%
%[eph, subFrameStart,SOW] = CNAVdecoding(I_P_InputBits)
%
%   Inputs:
%       I_P_InputBits   - output from the tracking function
%
%   Outputs:
%       firstSubFrame   - Starting positions of the first message in the
%                       input bit stream I_P_InputBits in each channel.
%                       The position is CNAV bit(20ms before convolutional decoding)
%                       count since start of tracking. Corresponding value will
%                       be set to inf if no valid preambles were detected in
%                       the channel.
%       SOW             - Time Of Week (SOW) of the first message(in seconds).
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
eph = eph_structure_init();

% Preamble search can be delayed to a later point in the tracking results
% to avoid noise due to tracking loop transients
searchStartOffset = 0;

% Starting positions of the first message in the input bit stream
firstSubFrame = inf;

% TOW of the first message
TOW = inf;

%Creates a cyclic redundancy code (CRC) detector System object
crcDet = comm.CRCDetector([24 23 18 17 14 11 10 7 6 5 4 3 1 0]);

% Convert CNAV-producing convolutional code polynomials to trellis description
% Note that the difference from GPS is that the second branch G2 is
% inverted at the end (see ICD)
trellis = poly2trellis(7,[171 133]);

% Viterbi traceback depth for vitdec(function)
tblen = 35;

%--- Generate the sync pattern --------------------------------------------
% Secondary code is "842E9"
secondCode = [-1 1 1 1   1 -1 1 1   1 1 -1 1   -1 -1 -1 1   -1 1 1 -1];

% The preamble is [1 0 1 1 0 1 1 1 0 0 0 0], and
% the antipodal form of the preamble pattern is:
preamble_bits = [-1 1 -1 -1 1 -1 -1 -1 1 1 1 1];
sync_bits = preamble_bits <0;

% "Upsample" the preamble - make 25 values per one bit. The preamble must be
% found with precision of a sample.
preamble_ms = kron(preamble_bits, secondCode);

% Use the prompt correlator as the symbol stream, skip start of record if
%   set in initSettings to avoid tracking loop transients
bits = I_P(1 + searchStartOffset : end);

% Now threshold the output and convert it to -1 and +1
bits(bits > 0)  =  1;
bits(bits <= 0) = -1;

% Correlate tracking output with the preamble
tlmXcorrResult = xcorr(bits, preamble_ms);

% Find all starting points of all preamble-like patterns -----------------
clear index
xcorrLength = (length(tlmXcorrResult) +  1) /2;

% Find at what index/ms the preambles start
index = find(abs(tlmXcorrResult(xcorrLength : xcorrLength * 2 - 1)) > 239.99)';

% Analyze detected preamble-like patterns ================================
for ind = 1:numel(index) % For each occurrence
    
    % Check distance to all other possible sync patterns for current start
    index2 = index - index(ind);
    
    % Check spacing, need two 250 symbol pages (even and odd) and at least
    % one subframe (2500 nav bits)
    if (~isempty(find(index2 == 500*20,1))) && ...
            ( length(bits(index(ind):end))> 2500*20 )
        
        %=== Read bit values for CRC-24Q check and ephemeris decoding ======
        % Search every possible preamble pattern.
        temp_bits = bits(index(ind):index(ind)+ 2500*20 -1);
        
        % Group every 10 I_P to a row for corresponding bit
        I_P_group = reshape(temp_bits,20,[])';
        
        % Wipe off the 2nd code and form data bits
        navBits = sum(I_P_group .* repmat(secondCode,2500,1),2)';
        navBits = navBits <0;
        
        %--- Correct polarity of the all data bits according to preamble bits
        if(~isequal(navBits(1:12),sync_bits))
            navBits = not(navBits);
        end   
        
        % Decoding the even page part -------------------------------------
        % Pull out implied pages from the detected preamble
        pageSymInt = navBits(13:500);
        
        % De-interleave symbols
        symMat = reshape(pageSymInt,61,8)';
        pageSym = reshape(symMat,1,[])';
        
        % Restore the inverted G2 convolutional-code branch.
        pageSym(2:2:end) = ~pageSym(2:2:end);
        % Remove convolutional encoding from implied pages
        decBits = vitdec(pageSym,trellis,tblen,'trunc','hard');
        
        % Remove tail bits
        decBits = decBits(1:238);
        
        pageBits = decBits';
        
        % Detect errors in input data using CRC ---------------------------
        [~,frmError] = step(crcDet,pageBits');
        
        % CRC-24Q check was OK. Then to decode ephemeris message by message.
        if (~frmError)
            %--- Just save for first message ------------------------------
            % firstSubFrame is the starting positions of the first message
            % in the input bit stream dataBits
            firstSubFrame = index(ind) + searchStartOffset;
            
            %--- Ephemeris decoding ---------------------------------------
            eph = ephemeris(navBits',eph);
            TOW = eph.TOW;
            break
        end % if CRC is OK ...
        
    end % Check spacing, need two 250 symbol pages (even and odd)
end % for ind = 1:size(index)
