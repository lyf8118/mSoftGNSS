function L2OcpCodesTable = makeL2OCpTable(settings,PRN)
%Function generates L2OCp primary codes for specified satellites based on the settings
%provided in the structure "settings". The codes are digitized at the
%sampling frequency specified in the settings structure.
%One row in the "L2OcpCodesTable" is one L2OCp primary code. The row number is the PRN
%number of the L2OCp code.
%
%L2OcpCodesTable = makeL2OCpTable(settings,PRN)
%
%   Inputs:
%       settings          - receiver settings
%       PRN               - PRN number of the sequence.
%   Outputs:
%       L2OcpCodesTable   - an array of arrays (matrix) containing L2OCp codes
%                       for the specified PRN

%--------------------------------------------------------------------------
%                         CU Multi-GNSS SDR  
% (C) Updated by Yafeng Li, Nagaraj C. Shivaramaiah and Dennis M. Akos
% Based on the original work by Darius Plausinaitis,Peter Rinder, 
% Nicolaj Bertelsen and Dennis M. Akos
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
%$Id: makeL2OCpTable.m,v 1.1.2.6 2006/08/14 11:38:22 dpl Exp $

%--- Find number of samples per spreading code ----------------------------
samplesPerCode = round(settings.samplingFreq / ...
    (settings.codeFreqBasis / settings.codeLength));

%--- Find time constants --------------------------------------------------
ts = 1/settings.samplingFreq;   % Sampling period in sec
fineFactor = settings.L2OCFineFactor;
tc = 1 / (settings.codeFreqBasis * fineFactor);   % PRN chip period in sec

%--- Generate L2OCp primary code for given PRN -----------------------------------
L2OcpCode = generateL2OcpBOCCode(PRN);  
%L2OcpCode = [L2OcpCode(end) L2OcpCode L2OcpCode(1)];

%=== Digitizing =======================================================

%--- Make index array to read L2OC code values -------------------------
% The length of the index array depends on the sampling frequency -
% number of samples per millisecond (because one primary code period is one
% millisecond).
codeValueIndex = ceil((ts * (0:samplesPerCode-1)) / tc);

%--- Correct the last index (due to number rounding issues) -----------
codeValueIndex(end) = settings.codeLength;
codeValueIndex(1) = 1;
codeValueIndex(end) = settings.codeLength * fineFactor;

%--- Make the digitized version of the L2OCp code -----------------------
% The "upsampled" code is made by selecting values from the L2OCp code
% chip array for the time instances of each sample.
L2OcpCodesTable = L2OcpCode(codeValueIndex);


