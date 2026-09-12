function L5QCodesTable = makeL5QTable(PRN,settings)
%Function generates L5Q codes for the specified PRN based on the settings
%provided in the structure "settings". The codes are digitized at the
%sampling frequency specified in the settings structure.
%One row in the "L5QCodesTable" is one L5Q code. The row number is the PRN
%number of the L5Q code.
%
%L5QCodesTable = makeL5QTable(settings, PRN)
%
%   Inputs:
%       settings        - receiver settings
%   Outputs:
%       L5QCodesTable   - an array of arrays (matrix) containing L5Q codes
%                       for all satellite PRNs

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
%$Id: makeL5QTable.m,v 1.1.2.6 2006/08/14 11:38:22 dpl Exp $

%--- Find number of samples per spreading code ----------------------------
samplesPerCode = round(settings.samplingFreq / ...
    (settings.codeFreqBasis / settings.codeLength));

%--- Find time constants --------------------------------------------------
ts = 1/settings.samplingFreq;   % Sampling period in sec
tc = 1/settings.codeFreqBasis;  % C/A chip period in sec


%--- Generate L5 code for given PRN -----------------------------------
L5QCode = generateL5Qcode(PRN,settings);

%=== Digitizing =======================================================

%--- Make index array to read L5Q code values -------------------------
% The length of the index array depends on the sampling frequency -
% number of samples per millisecond (because one L5Q code period is one
% millisecond).
codeValueIndex = ceil((ts * (1:samplesPerCode)) / tc);

%--- Correct the last index (due to number rounding issues) -----------
codeValueIndex(end) = settings.codeLength;

%--- Make the digitized version of the L5Q code -----------------------
% The "upsampled" code is made by selecting values from the L5 code
% chip array (L5Q code) for the time instances of each sample.
L5QCodesTable = L5QCode(codeValueIndex);

