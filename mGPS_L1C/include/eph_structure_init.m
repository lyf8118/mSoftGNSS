function eph= eph_structure_init()
% This is in order to make sure variable 'eph' for each SV has a similar
% structure when only one or even none of the three requisite sub-frames
% is decoded for a given PRN.
%--------------------------------------------------------------------------
%                         CU Multi-GNSS SDR  
% (C) Written by Yafeng Li, Nagaraj C. Shivaramaiah and Dennis M. Akos

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

%CVS record:
%$Id: eph_structure_init.m,v 1.1.2.7 2018/03/06 11:38:22 dpl Exp $

%% ===== Initialization & Flags =====================================
    eph.flag = [];           % Flag: 1 if a full set of ephemeris is valid/decoded
    eph.TOW  = [];          % Computed Time of Week (seconds)

    %% ===== Subframe 1: Time of Interval (TOI) =========================
    % Bits 1-9 in Subframe 1 (implied 9 bits in your decoder)
    % Note: Your decoder treats TOI as a local variable but uses it for TOW calculation.
    % We don't strictly need to store it in 'eph' unless you want to persist it.
    
    %% ===== Subframe 2: Ephemeris & Clock (Non-Variable) ===============
    % --- General Data Fields ---
    eph.WN = [];            % Week Number (13 bits)
    eph.ITOW = [];          % Interval Time of Week (8 bits)
    eph.t_op = [];          % Data Sequence Propagation Time (11 bits)
    eph.L1CHealth = [];     % L1C Signal Health (1 bit)
    eph.URAEDIndex = [];    % Elevation Dependent URA Index (5 bits)
    eph.t_oe = [];          % Ephemeris Reference Time (11 bits)
    
    % --- Orbit Parameters (Table 3.5-1) --------------------------------
    eph.deltaA = [];        % Semi-major axis difference (26 bits)
    eph.ADot = [];          % Change rate in semi-major axis (25 bits)
    eph.delta_n_0 = [];     % Mean Motion difference (17 bits)
    eph.delta_n_0Dot = [];  % Rate of Mean Motion difference (23 bits)
    eph.M_0 = [];           % Mean anomaly at reference time (33 bits)
    eph.e = [];             % Eccentricity (33 bits)
    eph.omega = [];         % Argument of perigee (33 bits)
    eph.omega_0 = [];       % Longitude of Ascending Node (33 bits)
    eph.i_0 = [];           % Inclination angle at reference time (33 bits)
    eph.omegaDot = [];      % Rate of Right Ascension (17 bits)
    eph.IDOT = [];          % Rate of inclination angle (15 bits)
    
    % --- Harmonic Correction Terms -------------------------------------
    eph.C_is = [];          % Amplitude of sine harmonic correction to inclination (16 bits)
    eph.C_ic = [];          % Amplitude of cosine harmonic correction to inclination (16 bits)
    eph.C_rs = [];          % Amplitude of sine correction term to orbit radius (24 bits)
    eph.C_rc = [];          % Amplitude of cosine correction term to orbit radius (24 bits)
    eph.C_us = [];          % Amplitude of sine harmonic correction to arg of latitude (21 bits)
    eph.C_uc = [];          % Amplitude of cosine harmonic correction to arg of latitude (21 bits)
    
    % --- Accuracy & Clock Parameters -----------------------------------
    eph.URANED0Index = [];  % URA NED Index 0 (5 bits)
    eph.URANED1Index = [];  % URA NED Index 1 (3 bits)
    eph.URANED2Index = [];  % URA NED Index 2 (3 bits)
    
    eph.a_0 = [];           % SV Clock Bias Correction (af0) (26 bits)
    eph.a_1 = [];           % SV Clock Drift Correction (af1) (20 bits)
    eph.a_2 = [];           % SV Clock Drift Rate Correction (af2) (10 bits)
    eph.T_GD = [];          % Group Delay Differential (13 bits)
    
    % --- Inter-Signal Corrections (ISC) --------------------------------
    eph.ISC_L1Cp = [];      % ISC for L1C pilot (13 bits)
    eph.ISC_L1Cd = [];      % ISC for L1C data (13 bits)
    
    eph.ISF = [];           % Integrity Status Flag (1 bit)
    eph.WN_op = [];         % CEI Data Sequence Propagation Week Number (8 bits)

    %% ===== Subframe 3: Variable Data (Common Header) ==================
    % Bits 1-8
    eph.PRN = [];           % PRN Number (8 bits) - Critical for L1C identity
    % Bits 9-14
    eph.PageID1 = [];        % Page 1 Decoded Flag
    eph.PageID2 = [];        % Page 2 Decoded Flag

    %% ===== Subframe 3 - Page 1: UTC & Ionosphere ======================
    % --- UTC Parameters ---
    eph.A0_n = [];          % UTC Bias coefficient (16 bits)
    eph.A1_n = [];          % UTC Drift coefficient (13 bits)
    eph.A2_n = [];          % UTC Drift rate coefficient (7 bits)
    eph.delta_t_LS = [];    % Current Leap Seconds (8 bits)
    eph.t_ot = [];          % UTC Reference Time (16 bits)
    eph.WN_ot = [];         % UTC Reference Week (13 bits)
    eph.WN_LSF = [];        % Future Leap Second Week (13 bits)
    eph.DN = [];            % Day Number (4 bits)
    eph.delta_t_LSF = [];   % Future Leap Seconds (8 bits)
    
    % --- Ionospheric Parameters (Klobuchar) ---
    eph.alpha0 = [];        % Alpha 0 (8 bits)
    eph.alpha1 = [];        % Alpha 1 (8 bits)
    eph.alpha2 = [];        % Alpha 2 (8 bits)
    eph.alpha3 = [];        % Alpha 3 (8 bits)
    eph.beta0 = [];         % Beta 0 (8 bits)
    eph.beta1 = [];         % Beta 1 (8 bits)
    eph.beta2 = [];         % Beta 2 (8 bits)
    eph.beta3 = [];         % Beta 3 (8 bits)
    
    % --- ISC Parameters (Page 1) ---
    eph.ISC_L1CA = [];      % ISC L1C/A (13 bits)
    eph.ISC_L2C = [];       % ISC L2C (13 bits)
    eph.ISC_L5I5 = [];      % ISC L5I5 (13 bits)
    eph.ISC_L5Q5 = [];      % ISC L5Q5 (13 bits)

    %% ===== Subframe 3 - Page 2: GGTO & EOP ============================
    % --- GGTO Parameters ---
    eph.GGTO.GNSS_ID = [];  % GNSS ID (3 bits)
    eph.GGTO.t_GGTO = [];   % GGTO Ref Time (16 bits)
    eph.GGTO.WN_GGTO = [];  % GGTO Ref Week (13 bits)
    eph.GGTO.A0_GGTO = [];  % Time Bias (16 bits)
    eph.GGTO.A1_GGTO = [];  % Time Drift (13 bits)
    eph.GGTO.A2_GGTO = [];  % Time Drift Rate (7 bits)
    
    % --- EOP Parameters ---
    eph.EOP.t_EOP = [];     % EOP Ref Time (16 bits)
    eph.EOP.PM_X = [];      % Polar Motion X (21 bits)
    eph.EOP.PM_X_dot = [];  % Polar Motion X Rate (15 bits)
    eph.EOP.PM_Y = [];      % Polar Motion Y (21 bits)
    eph.EOP.PM_Y_dot = [];  % Polar Motion Y Rate (15 bits)
    eph.EOP.Delta_UTGPS = [];     % UT1-UTC Diff (31 bits)
    eph.EOP.Delta_UTGPS_dot = []; % UT1-UTC Diff Rate (19 bits)

end
