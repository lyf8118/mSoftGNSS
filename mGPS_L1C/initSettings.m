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

%--------------------------------------------------------------------------
%                         CU Multi-GNSS SDR
% (C) Developed for GPS L1C SDR by Yafeng Li.
% Based on the original SoftGNSS SDR framework by Darius Plausinaitis,
% Peter Rinder, Nicolaj Bertelsen and Dennis M. Akos.
%
% Reference: Adapted within the CU Multi-GNSS SDR receiver framework for GPS L1C.
% Signal-specific settings in this file refer to GPS L1C.
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

% CVS record:
% $Id: initSettings.m,v 1.9.2.31 2006/08/18 11:41:57 dpl Exp $

%% Important operation modes configuration ==========================
% Number of milliseconds to be processed used 36000 + any transients (see
% below - in Nav parameters) to ensure nav subframes are provided
settings.msToProcess        = 40000;        %[ms]
% Number of channels to be used for signal processing
settings.numberOfChannels   = 15;
% Move the starting point of processing. Can be used to start signal
% processing at any point in a long data record. The offset is expressed
% in logical samples: one real value or one interleaved complex I/Q pair
% counts as one sample.
settings.skipNumberOfSamples  = 0e6;

%% Raw signal file name and related parameter =======================
settings.fileName           = "../IF_Data_Set/GPS L1C.iq";
% Data type used to store one sample

settings.dataType           = 'int8';
% File Types
%1 - 8 bit real samples S0,S1,S2,...
%2 - 8 bit I/Q samples I0,Q0,I1,Q1,I2,Q2,...
settings.fileType           = 2;
%Intermediate, sampling, code freqency and code length
settings.IF                 = 0e6;         %[Hz]
settings.samplingFreq       = 25e6;        %[Hz]
% Front-end bandwidth for calculation of weighting factor in wideband (WB)
% tracking mode (this is required only for WB tracking mode)
settings.FEBW   = 27e6;                    % [Hz]

%% L1C signal parameters ============================================
% Define number of chips in a L1C primary code period.
settings.codeLength    = 10230;        % [Chips]
% Nominal rate of L1C primary code.
settings.codeFreqBasis = 1.023e6;      % [Hz]
% Carrier center frequency
settings.carrFreqBasis = 1575.42e6;    % [Hz]

%% CNAV-2 LDPC decoder settings =====================================
% Min-Sum decoder options. Set ldpcMsOffset to 0 for standard Min-Sum.
settings.ldpcMaxIterations = 50;
settings.ldpcMsAlpha       = 1;
settings.ldpcMsOffset      = 0.2;

%% Acquisition settings =============================================
% Skips acquisition in the script postProcessing.m if set to 1
settings.skipAcquisition    = 0;
% Enable use of GPU acceleration for acquisition
settings.gpuACQflag         = 1;             % 0 - Off; 1 - On
% List of satellites to look for.
settings.acqSatelliteList   = 1:32;
% Enable use of pilot channel to do acquisition
settings.pilotACQflag       = 1;             % 0 - Off; 1 - On
% It is single sideband, so the whole search band is twice of it.
settings.acqSearchBand      = 5000;        % [Hz]
% Frequency search step
settings.acqStep            = 50;          % [Hz]
% Threshold for the signal presence decision rule
settings.acqThreshold       = 2;

%% Tracking loops settings ==========================================
% Enable/disable use of SIMD/GPU MEX to accelerate tracking
settings.trkMode                 = 1;    % 0 - Channel-serial tracking;
                                         % 1 - Channel-parallel tracking;
% Enable/disable use of SIMD/GPU MEX to accelerate tracking
settings.correlatorType          = 2;    % 0 - Matlab correlator;
                                         % 1 - SIMD correlator;
                                         % 2 - GPU correlator;
% Enable/disable use of pilot channel for tracking
settings.pilotTRKflag            = 1;          % 0 - Off; 1 - On

% Number of right-shift bits for IF data to prevent SIMD correlator
% overflow when ADC valid bits occupy the high bits
if settings.correlatorType == 1
    settings.rShiftBits          = 5;   % 0 - No shift;
end
% Code tracking loop parameters
settings.dllDampingRatio         = 0.7;
settings.dllNoiseBandwidth       = 1;       %[Hz]
% Half of the E-L spacing should be kept below approximately 0.0765 chip
% for full-band tracking mode (see Betz 2016)
settings.dllCorrelatorSpacing    =  0.07;     % [primary code chip]
% Carrier tracking loop parameters
settings.pllDampingRatio         = 0.7;
settings.pllNoiseBandwidth       = 10;        % [Hz] 
% Integration time for DLL and PLL
settings.intTime = 0.01;                      % [s]

%% Navigation solution settings =====================================
% Period for calculating pseudoranges and position
settings.navSolPeriod       = 200;          %[ms]  200
% Elevation mask to exclude signals from satellites at low elevation
settings.elevationMask      = 5;           %[degrees 0 - 90]
% Enable/disable use of tropospheric correction
settings.useTropCorr        = 1;            % 0 - Off; 1 - On
% True position of the antenna in UTM system (if known). Otherwise enter
% all NaN's and mean position will be used as a reference .
settings.truePosition.E     = nan;
settings.truePosition.N     = nan;
settings.truePosition.U     = nan;

%% Plot settings ==========================================================
% Enable/disable plotting of the tracking results for each channel
settings.plotTracking    = 1;            % 0 - Off; 1 - On
settings.plotNavigation  = 1;

%% Constants ========================================================
settings.c                  = 299792458;    % The speed of light, [m/s]
settings.startOffset        = 68.802;       % [ms] Initial sign. travel time

%% CNo Settings======================================================
% Number of correlation values used to compute each C/No point
settings.CNoInterval = 20;
