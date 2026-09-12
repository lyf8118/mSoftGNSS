function [eph, TOD] = ephemeris(navBits, eph)
%Function decodes ephemerides and TOD from the given bit stream. The stream
%(array) in the parameter BITS must contain 1275 bits. The first element in
%the array must be the first bit of a string. The string ID of the
%first string in the array is not important.
%
%
%[eph, TOD] = ephemeris(navBits, eph)
%
%   Inputs:
%       bits        - bits of the navigation messages (15 strings).
%                   Type is character array and it must contain only
%                   characters '0' or '1'.
%   Outputs:
%       TOD         - Time Of Day (TOD) of the first string in the bit
%                   stream (in seconds)
%       eph         - SV ephemeris
%
%--------------------------------------------------------------------------
%                           CU Multi-GNSS SDR  
%
% Copyright (C) Darius Plausinaitis and Kristin Larson
% Written by Darius Plausinaitis and Kristin Larson
%
% GLONASS modification by Jakob Almqvist
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
%$Id: ephemeris.m,v 1.1.2.7 2006/08/14 11:38:22 dpl Exp $

    %% Initialization============================================
    if nargin < 2 || isempty(eph)
        % Initialize empty structure if not provided
        eph = struct('SV_ID', [], 'TOD', [], 'Tb', [], ...
                     'N4', [], 'NT', [], ...
                     'Tau', [], 'Gamma', [], 'Tau_c', [], ...
                     'X', [], 'Y', [], 'Z', [], ...
                     'dX', [], 'dY', [], 'dZ', [], ...
                     'ddX', [], 'ddY', [], 'ddZ', []);
    end
    
    % Ensure input is a column vector
    navBits = navBits(:);
    
    % Check length 
    if length(navBits) ~= 300
        error('GLONASS L3OCd string must be exactly 300 bits long.');
    end
    
    %% Parse Common Service Fields ============================================
    % Type 
    typeVal = bin2dec_unsigned(navBits(21:26));
    % TS
    tsVal = bin2dec_unsigned(navBits(27:41));
    TOD   = tsVal * 3; 
    eph.TOD = TOD;
    %Satellite ID
    eph.SV_ID = bin2dec_unsigned(navBits(42:47));
    % L3Oc health
    eph.Health = navBits(48);
    %Data Validity Attribute
    eph.DataValid = navBits(49);
    %% 3. Parse Data Fields based on String Type ============================================
    switch typeVal
        case 10  %--- It is Message Type 10 ---------------------------------------------
             % It contains first part of ephemeris parameters
            % N4 
            eph.N4 = bin2dec_unsigned(navBits(58:62));
            % NT 
            eph.NT = bin2dec_unsigned(navBits(63:73));
            % tb,LSB: 90 s
            eph.Tb = bin2dec_unsigned(navBits(83:92)) * 90;
            % Tau,Scale: 2^-38 s
            eph.Tau = bin2dec_twos(navBits(123:154)) * (2^-38);  
            % Gamma,Scale: 2^-48 
            eph.Gamma = bin2dec_twos(navBits(155:173)) * (2^-48); 
            % Tau_c ,Scale: 2^-31 s 
            eph.Tau_c = bin2dec_twos(navBits(189:228)) * (2^-31); 
            
        case 11 %--- It is Message Type 11 ---------------------------------------------
           % It contains second part of ephemeris parameter
            % X coordinates 
            eph.X = bin2dec_twos(navBits(58:97)) * (2^-20);
            % Y coordinates 
            eph.Y = bin2dec_twos(navBits(98:137)) * (2^-20);
            % Z coordinates 
            eph.Z = bin2dec_twos(navBits(138:177)) * (2^-20);
            % dX Velocity (35 bits): 178 to 212
            eph.dX = bin2dec_twos(navBits(178:212)) * (2^-30);
            % dY Velocity (35 bits): 213 to 247
            eph.dY = bin2dec_twos(navBits(213:247)) * (2^-30);
            
        case 12 %--- It is Message Type 12 ---------------------------------------------
            % It contains third part of ephemeris parameter
            % dZ Velocity 
            eph.dZ = bin2dec_twos(navBits(58:92)) * (2^-30);
            % ddX Acceleration
            eph.ddX = bin2dec_twos(navBits(93:107)) * (2^-39);
            % ddY Acceleration 
            eph.ddY = bin2dec_twos(navBits(108:122)) * (2^-39);
            % ddZ Acceleration
            eph.ddZ = bin2dec_twos(navBits(123:137)) * (2^-39);
    end

end

%% Helper Functions for GLONASS Data Format ============================================
function val = bin2dec_unsigned(bits)
    % Convert logical binary array to unsigned decimal 
    val = 0;
    len = length(bits);
    for k = 1:len
        if bits(k) == 1
            val = val + 2^(len - k);
        end
    end
end
function val = bin2dec_twos(bits)
    % Convert logical binary array to signed decimal using SIGN-MAGNITUDE.
    signBit = bits(1);
    magnitudeBits = bits(2:end);
    mag = bin2dec_unsigned(magnitudeBits);
    
    if signBit == 0
        val = mag;  % Positive
    else
        val = -mag; % Negative
    end
end