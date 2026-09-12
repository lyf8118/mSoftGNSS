function B1CPilot = generatePilotBOC61(settings,PRN)
% This function generates the BDS-3 B1C pilot BOC(6,1) code in bipolar
% format (-1, +1). PRNs 1 to 63 are supported by the parameter table below.
%
% B1CPilot = generatePilotBOC61(settings, PRN)
%
%   Inputs:
%       PRN         - PRN number of the sequence.
%
%   Outputs:
%       B1CPilot    - vector containing the B1C pilot BOC(6,1) waveform.

%--------------------------------------------------------------------------
%                         CU Multi-GNSS SDR  
% (C) Written by Yafeng Li, Nagaraj C. Shivaramaiah and Dennis M. Akos

% Reference: Li, Y., Shivaramaiah, N.C. & Akos, D.M. Design and 
% implementation of an open-source BDS-3 B1C SDR receiver. 
% GPS Solut (2019) 23: 60. https://doi.org/10.1007/s10291-019-0853-z
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
%$Id: generatePilotBOC61.m,v 1.1.2.5 2017/11/27 22:00:00 dpl Exp $

wp_pilot = [796    7575;   156    2369;   4198   5688;   3941   539;...
      1374   2270;   1338   7306;   1833   6457;   2521   6254;...
      3175   5644;   168    7119;   2715   1402;   4408   5557;...
      3160   5764;   2796   1073;   459    7001;   3594   5910;...
      4813   10060;  586    2710;   1428   1546;   2371   6887;...
      2285   1883;   3377   5613;   4965   5062;   3779   1038;...
      4547   10170;  1646   6484;   1430   1718;   607    2535;...
      2118   1158;   4709   526;    1149   7331;   3283   5844;...
      2473   6423;   1006   6968;   3670   1280;   1817   1838;...
      771    1989;   2173   6468;   740    2091;   1433   1581;...
      2458   1453;   3459   6252;   2155   7122;   1205   7711;...
      413    7216;   874    2113;   2463   1095;   1106   1628;...
      1590   1713;   3873   6102;   4026   6123;   4272   6070;...
      3556   1115;   128    8047;   1200   6795;   130    2575;...
      4494   53;     1871   1729;   3073   6388;   4386   682;...
      4098   5565;   1923   7160;   1176   2277];
        
        
N = 10243;

legendre = zeros(1,N);
for ind = 1:N-1
    legendre(ind+1) = JacobiSymbol(ind,N);
end

legendre(legendre==-1) = 0;

% The length of the primary code of B1C is 10230
Primary = zeros(1,settings.codeLength);

p = wp_pilot(PRN,2);
w = wp_pilot(PRN,1);

for ind=0:10229

    k = mod((ind+p-1), N);
    
    Primary(ind+1) = xor( legendre(k+1), legendre(mod(k+w,N)+1) );
end


% Change to bipolar format (-1, +1)
Primary = 1-2*Primary;

% %--- Add BOC(6,1) subcarrier  ---------------------------------------------
B1CPilot = zeros(1, settings.codeLength*12);

for jj = 1:settings.codeLength
    kk = (jj-1)*12;
    for ii = 1 : 12
        B1CPilot(kk+ii)   =  (-1)^ii * Primary(jj);
    end
end

%% Code-generation test placeholder





