function E6CCodesTable = makeE6CTable(PRN,settings)
%Function generates E6C codes for the specified PRN based on the settings
%provided in the structure "settings". The codes are digitized at the
%sampling frequency specified in the settings structure.
%One row in the "E6CCodesTable" is one E6C code. The row number is the PRN
%number of the E6C code.
%
%E6CCodesTable = makeE6CTable(settings, PRN)
%
%   Inputs:
%       settings        - receiver settings
%   Outputs:
%       E6CCodesTable   - an array of arrays (matrix) containing E6C codes
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

%--- Find number of samples per spreading code ----------------------------
samplesPerCode = round(settings.samplingFreq / ...
    (settings.codeFreqBasis / settings.codeLength));

%--- Find time constants --------------------------------------------------
ts = 1/settings.samplingFreq;   % Sampling period in sec
tc = 1/settings.codeFreqBasis;  % E6C chip period in sec

%--- Generate E6C primary code for given PRN -----------------------------------
E6CCode = generateE6Ccode(PRN,2);

%=== Digitizing =======================================================

%--- Make index array to read E6C code values -------------------------
% The length of the index array depends on the sampling frequency -
% number of samples per millisecond (because one E6C code period is one
% millisecond).
codeValueIndex = ceil((ts * (1:samplesPerCode)) / tc);

%--- Correct the last index (due to number rounding issues) -----------
codeValueIndex(end) = 5115;

%--- Make the digitized version of the E6C code -----------------------
% The "upsampled" code is made by selecting values from the E6 code
% chip array (E6C code) for the time instances of each sample.
E6CCodesTable = E6CCode(codeValueIndex);