function CAcode = generateB12Icode(PRN)
% generateB12Icode.m generates one of the 37 Beidou satellite B1I codes.
%
% B12Icode = generateB12IcodeBDS(PRN)
%
%   Inputs:
%       PRN         - PRN number of the sequence.
%
%   Outputs:
%       B12Icode      - a vector containing the desired B1I/B2I code sequence
%                   (chips).

%--------------------------------------------------------------------------
%                           CU Multi-GNSS SDR
% (C) Developed for BDS B1I/B2I SDR by Yafeng Li, Daehee Won,
% Nagaraj C. Shivaramaiah and Dennis M. Akos.
% Based on the original framework for GPS C/A SDR by Darius Plausinaitis,
% Peter Rinder, Nicolaj Bertelsen and Dennis M. Akos
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
%$Id: generateB12Icode.m,v 1.1.2.5 2006/08/14 11:38:22 dpl Exp $

%--------------------------------------------------------------------------
% Modified for Beidou by Daehee Won
% Final update: Sep. 26, 2013

%--- Generate G1 code
% Initialize G1 output
g1 = zeros(1, 2046);    % Chip length of Beidou is 2046
% Load shift register (11 bits)
reg = -1*[-1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1];

% Generate G1 signal chips based on the G1 generator polynomial
for i=1:2046
    g1(i)     = reg(11);
    saveBit   = reg(1)*reg(7)*reg(8)*reg(9)*reg(10)*reg(11);
    reg(2:11) = reg(1:10);
    reg(1)    = saveBit;
end

%--- Generate G2 code
% Initialize G2 output
g2 = zeros(1, 2046);    % Chip length of Beidou is 2046
% Load shift register (11 bits)
reg = -1*[-1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1];
% Phase assignment of G2 sequence for PRN 1-37
g2s1 = [1, 1, 1, 1, 1, 1, 1, 1, 2, 3, ...
    3, 3, 3, 3, 3, 3, 4, 4, 4, 4, ...
    4, 4, 5, 5, 5, 5, 5, 6, 6, 6, ...
    6, 8, 8, 8, 9, 9, 10];
g2s2 = [ 3, 4, 5, 6, 8, 9,10,11, 7, 4, ...
    5, 6, 8, 9,10,11, 5, 6, 8, 9, ...
    10,11, 6, 8, 9,10,11, 8, 9,10, ...
    11, 9,10,11,10,11,11];

% Phase assignment of G2 sequence for PRN 38-63 (BDS3 adds)
g2s1_3 = [1, 1, 1, 1,  1,  1, 1, 1, 1, 1,  1, 1,  1, 1, 1,  ...
    1, 2, 2, 2, 3, 3, 3, 3,  3,  3, 3];
g2s2_3 = [2, 3, 3, 3,  3,  3, 4, 4, 5, 5,  5, 5,  6, 8, 9,  ...
    9, 3, 5, 7, 4, 4, 5, 5,  5,  5, 6];
g2s3_3 = [7, 4, 6, 8, 10, 11, 5, 9, 6, 8, 10, 11, 9, 9, 10,  ...
    11, 7, 7, 9, 5, 9, 6, 8, 10, 11, 9];
% Generate G2 signal chips based on the G2 generator polynomial
if PRN <= 37
    for i=1:2046
        g2(i)       = reg(g2s1(PRN))*reg(g2s2(PRN));
        saveBit     = reg(1)*reg(2)*reg(3)*reg(4)*reg(5)*reg(8)*reg(9)*reg(11);
        reg(2:11)   = reg(1:10);
        reg(1)      = saveBit;
    end
elseif PRN >= 38 && PRN <= 63
    for i=1:2046
        g2(i)       = reg(g2s1_3(PRN-37))*reg(g2s2_3(PRN-37))* reg(g2s3_3(PRN-37));
        saveBit     = reg(1)*reg(2)*reg(3)*reg(4)*reg(5)*reg(8)*reg(9)*reg(11);
        reg(2:11)   = reg(1:10);
        reg(1)      = saveBit;
    end
end

%--- Form single sample B1I/B2I code by multiplying G1 and G2 -----------------
CAcode = -(g1 .* g2);
