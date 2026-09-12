function settings = initSettings()
%Functions initializes and saves settings. Settings can be edited inside of
%the function or updated from the command line.
%Edit this file to configure receiver parameters.  
%
%All settings are described inside function code.
%
%settings = initSettings()
%
%   Inputs: none
%
%   Outputs:
%       settings     - Receiver settings (a structure). 
%
%--------------------------------------------------------------------------
%                         CU Multi-GNSS SDR
% (C) Developed for Galileo E6B/E6C SDR by Yafeng Li, Nagaraj C. Shivaramaiah 
% and Dennis M. Akos. 
% Based on the original framework for GPS C/A SDR by Darius Plausinaitis,
% Peter Rinder, Nicolaj Bertelsen and Dennis M. Akos
%---------------------------------------------------------------------------
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

% CVS record:
% $Id: initSettings.m,v 1.9.2.31 2006/08/18 11:41:57 dpl Exp $

%% Processing settings ============================================
% Number of milliseconds to be processed used 50000 + any transients (see
% below - in Nav parameters) to ensure nav subframes are provided
settings.msToProcess        = 50000;        %[ms]
% Number of channels to be used for signal processing
settings.numberOfChannels   = 12;
% Move the starting point of processing. Can be used to start signal
% processing at any point in a long data record. The offset is expressed
% in logical samples: one real value or one interleaved complex I/Q pair
% counts as one sample.
settings.skipNumberOfSamples     = 0;

%% Raw signal file name and other parameter ===============================
% This is a "default" name of the data file (signal record) to be used in
% the post-processing mode
settings.fileName           ="../IF_Data_Set/GAL_E6.bin";  
% Data type used to store one sample
settings.dataType           = 'int16';
% File Types
%1 - 8 bit real samples S0,S1,S2,...
%2 - 8 bit I/Q samples I0,Q0,I1,Q1,I2,Q2,...                      
settings.fileType           = 2;
% Intermediate, sampling and code frequencies
settings.IF                 = 0.0e6;   % [Hz] 
settings.samplingFreq       = 30e6;             % [Hz] 
settings.codeFreqBasis      = 5.115e6;              % [Hz]
% Define number of chips in a E6B tiered code period
settings.codeLength          = 5115;

%% Acquisition settings ===========================================
% Enable/disable use of GPU (MATLAB function) to accelerate acquisition
settings.gpuACQflag         = 1;             % 0 - Off; 1 - On
% Skips acquisition in the script postProcessing.m if set to 1
settings.skipAcquisition    = 0;
% List of satellites to look for. Some satellites can be excluded to speed
% up acquisition. 
settings.acqSatelliteList   = 1:50; 
% Band around IF to search for satellite signal. Depends on max Doppler.
% It is single sideband, so the whole search band is twice of it.
settings.acqSearchBand          = 5000;           % [Hz]
% Non-coherent integration times after 1ms coherent integration
settings.acqNonCohTime      = 10;             %[ms] 
% Threshold for the signal presence decision rule
settings.acqThreshold       = 3;
% Frequency search step for coarse acquisition
settings.acqSearchStep      = 500;                % [Hz]

%% Tracking loops settings ========================================
% Enable/disable use of SIMD/GPU MEX to accelerate tracking
settings.trkMode                 = 1;   % 0 - Channel-serial tracking; 
                                        % 1 - Channel-parallel tracking;
% Enable/disable use of SIMD/GPU MEX to accelerate tracking
settings.correlatorType          = 2;   % 0 - Matlab correlator; 
                                        % 1 - SIMD correlator;
                                        % 2 - GPU correlator;
% Number of right-shift bits for IF data to prevent SIMD correlator 
% overflow when ADC valid bits occupy the high bits
if settings.correlatorType == 1
    settings.rShiftBits              = 4;   % 0 - No shift;
end

% Code tracking loop parameters
settings.dllDampingRatio         = 0.7;
settings.dllNoiseBandwidth       = 2;          %[Hz]
settings.dllCorrelatorSpacing    = 0.5;        % [chips]
% Carrier tracking loop parameters
settings.pllDampingRatio         = 0.7;
settings.pllNoiseBandwidth       = 15;         %[Hz]
% Integration time for DLL and PLL
settings.intTime                 = 0.001;      %[s]
% Enable/disable use of pilot channel for tracking
settings.pilotTRKflag            = 1;        % 0 - Off; 1 - On

%% Navigation solution settings ===========================================
% Period for calculating pseudoranges and position
settings.navSolPeriod       = 500;          % [ms]
% Elevation mask to exclude signals from satellites at low elevation
settings.elevationMask      = 5;            % [degrees 0 - 90]
% Enable/disable use of tropospheric correction
settings.useTropCorr        = 1;            % 0 - Off; 1 - On
% True position of the antenna in UTM system (if known). Otherwise enter
% all NaN's and mean position will be used as a reference .
settings.truePosition.E     = nan;
settings.truePosition.N     = nan;
settings.truePosition.U     = nan;

%% Plot settings ==========================================================
% Enable/disable plotting of the tracking results for each channel
settings.plotTracking       = 1;            % 0 - Off; 1 - On

%% Constants ==============================================================
settings.c                  = 299792458;    % The speed of light, [m/s]
settings.startOffset        = 68.802;       %[ms] Initial sign. travel time

%% CNo Settings====================================================
% Number of correlation values used to compute each C/No point 
  

settings.CNoInterval = 200;

%% E6 carrier frequency ===================================================
settings.carrFreqBasis    = 1278.75e6;       %[Hz]
