function DataCodesTable = makeDataTable(settings,PRN)
%Function generates GPS L1C data BOC(1,1) codes for specified PRN based
% on the settings provided in the structure "settings". The codes are 
%digitized at the sampling frequency specified in the settings structure.
%
%DataCodesTable = makeDataTable(settings,PRN)
%
%   Inputs:
%       settings        - receiver settings
%   Outputs:
%       DataCodesTable  - sampled GPS L1C data BOC(1,1) spreading waveform
%                         

%--------------------------------------------------------------------------
%                         CU Multi-GNSS SDR  
% (C) Developed for GPS L1C SDR by Yafeng Li, Nagaraj C. Shivaramaiah 
% and Dennis M. Akos. 
% Based on the original SoftGNSS SDR framework by Darius Plausinaitis,
% Peter Rinder, Nicolaj Bertelsen and Dennis M. Akos

% Reference: Adapted within the CU Multi-GNSS SDR receiver framework for GPS L1C.
% Signal-specific comments in this file refer to GPS L1C.
%
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
ts = 1/settings.samplingFreq;           % Sampling period in sec
tc = 1/settings.codeFreqBasis/2;        % L1C data BOC(1,1) chip period in sec

%--- Generate L1C data BOC(1,1) code for given PRN -----------------------------
DataCode= generateDataBOC11(settings,PRN);

%=== Digitizing ===========================================================

%--- Make index array to read L1C code values -----------------------------
% The length of the index array depends on the sampling frequency
codeValueIndex = ceil((ts * (1:samplesPerCode)) / tc);

%--- Correct the last index (due to number rounding issues) ---------------
codeValueIndex(end) = settings.codeLength*2;
codeValueIndex(1) = 1;

%--- Make the digitized version of the L1C code ---------------------------
% The "upsampled" code is made by selecting values from the L1C code
% chip array for the time instances of each sample.
DataCodesTable = DataCode(codeValueIndex);
