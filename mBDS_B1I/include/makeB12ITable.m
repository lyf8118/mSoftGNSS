function B12ICodesTable = makeB12ITable(PRN,settings)
%Function generates B1I/B2I codes for all specified PRN based on the settings
%provided in the structure "settings". The codes are digitized at the
%sampling frequency specified in the settings structure.
%
%B12ICodesTable = makeB12ITable(PRN,settings)
%
%   Inputs:
%       PRN             - specified PRN for B1I/B2I code
%       settings        - receiver settings
%   Outputs:
%       B12ICodesTable   - a vector containing the sampled B1I/B2I codes
%                       for specified PRN

%--------------------------------------------------------------------------
%                         CU Multi-GNSS SDR
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
%$Id: makeB12ITable.m,v 1.1.2.6 2020/01/16 11:38:22 dpl Exp $

%--------------------------------------------------------------------------
% Modified for Beidou by Yafeng Li
% Final update: 2020/01/20

%--- Find number of samples per spreading code ----------------------------
samplesPerCode = round(settings.samplingFreq / ...
    (settings.codeFreqBasis / settings.codeLength));

%--- Find time constants --------------------------------------------------
ts = 1/settings.samplingFreq;   % Sampling period in sec
tc = 1/settings.codeFreqBasis;  % B1I/B2I chip period in sec

%--- Generate B1I/B2I code for given PRN --------------------------------------
B12ICode = generateB12Icode(PRN);

%=== Digitizing ===========================================================
%--- Make index vector to read B1I/B2I code values ----------------------------
% The length of the index vector depends on the sampling frequency -
% number of samples per millisecond (because one B1I/B2I code period is one
% millisecond).
codeValueIndex = floor((ts * (0:samplesPerCode-1)) / tc)+1;

%--- Correct the last index (due to number rounding issues) ---------------
codeValueIndex(end) = 2046;  % B1I/B2I chip length is 2046

%--- Make the digitized version of the B1I/B2I code ---------------------------
% The "upsampled" code is made by selecting values from the B1I/B2I code
% chip vector (B12ICode) for the time instances of each sample.
B12ICodesTable = B12ICode(codeValueIndex);
