function B2bITable = makeB2bITable(PRN,settings)
%Function generates B2b I codes for the specified PRN based on the settings
%provided in the structure "settings". The codes are digitized at the
%sampling frequency specified in the settings structure.
%One row in the "B2bICodesTable" is one B2b I code. The row number is the PRN
%number of the B2b I code.
%
%B2bICodesTable = makeB2bITable(settings, PRN)
%
%   Inputs:
%       settings        - receiver settings
%   Outputs:
%       B2bICodesTable    - an array of arrays (matrix) containing B2b I codes
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
%$Id: makeB2bITable.m,v 1.1.2.6 2006/08/14 11:38:22 dpl Exp $

%--- Find number of samples per spreading code ----------------------------
samplesPerCode = round(settings.samplingFreq / ...
    (settings.codeFreqBasis / settings.codeLength));

%--- Find time constants --------------------------------------------------
ts = 1/settings.samplingFreq;   % Sampling period in sec
tc = 1/settings.codeFreqBasis;  % B2a chip period in sec


%--- Generate B2bI code for given PRN -----------------------------------
B2bICode = generateB2bICode(PRN,settings);

%=== Digitizing =======================================================

%--- Make index array to read B2bI code values -------------------------
% The length of the index array depends on the sampling frequency -
% number of samples per millisecond (because one B2bI code period is one
% millisecond).
codeValueIndex = ceil((ts * (1:samplesPerCode)) / tc);

%--- Correct the last index (due to number rounding issues) -----------
codeValueIndex(end) = settings.codeLength;

%--- Make the digitized version of the B2bI code -----------------------
% The "upsampled" code is made by selecting values from the B2bI code
% chip array (B2b I code) for the time instances of each sample.
B2bITable = B2bICode(codeValueIndex);
