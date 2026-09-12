function L3OcpCode = generateL3OcpCode(PRN)
% generateL3OcpCode.m generates one L3OCp code.
%
% L3OcpCode = generateL3OcpCode(PRN)
%
%   Inputs:
%       PRN         - PRN number of the sequence.
%
%   Outputs:
%       L3OcpCode      - a vector containing the desired L3OCp code sequence 
%                   (chips).  

%--------------------------------------------------------------------------
%                         CU Multi-GNSS SDR
% (C) Developed for GLONASS L3OC SDR by Yafeng Li, Nagaraj C. Shivaramaiah
% and Dennis M. Akos.
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
code_length = 10230;

persistent codes;
if isempty(codes)
    codes = containers.Map('KeyType', 'double', 'ValueType', 'any');
end

    function result = g2_shift(x)  
        result = [bitxor(bitxor(bitxor(x(14), x(13)), x(8)), x(4)) x(1:13)];
    end

    function result = g1_shift(x) 
        result = [bitxor(x(7), x(6)) x(1:6)];
    end

    function s = seq(n)  
        s = zeros(1, 7);
        for i = 1:7
            s(i) = bitand(bitshift(n, -(6 - (i-1))), 1);  
        end
    end

    function x = make_l3ocp(n)   
        g1 = seq(n+64);
        g2 = [0, 0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0];
        x = zeros(1, code_length);  
        
        for i = 1:code_length
            x(i) = bitxor(g1(7), g2(14));
            g1 = g1_shift(g1);
            g2 = g2_shift(g2);
        end
    end

    function c = l3ocp_code(n)
        if ~codes.isKey(n)
            codes(n) = make_l3ocp(n);
        end
        c = codes(n);
    end

    function result_array = L3Ocp(n) 
        c = l3ocp_code(n);
        result_array = zeros(1, code_length);
        for i = 1:code_length
            result_array(i) = c(i);
        end
    end

 L3OcpCode = L3Ocp(PRN);
 L3OcpCode = 1 - 2 * L3OcpCode; 
 end
