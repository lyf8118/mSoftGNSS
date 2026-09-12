function [eph] = ephemeris(navBitsBin,eph)
%Function decodes GPS L1C CNAV-2 ephemeris and TOW from the decoded frame
%bit stream. navBitsBin must contain one 883-bit sequence:
%TOI(9) + Subframe 2(600) + Subframe 3(274), represented as characters '0'
%and '1'.
%
%Function assumes BCH, LDPC, and CRC checks were completed before call.
%
%[eph] = ephemeris(navBitsBin,eph)
%
%   Inputs:
%       navBitsBin  - bits of the navigation messages.Type is character array
%                   and it must contain only characters '0' or '1'.
%       eph         - The ephemeris for each PRN is decoded message by message.
%                   To prevent loss of previous decoded messages, the eph structure
%                   must be passed onto this function.
%   Outputs:
%       TOW         - Time Of Week (TOW) of the first decoded CNAV-2 frame
%                   in the bit stream (in seconds)
%       eph         - SV ephemeris

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
%$Id: ephemeris.m,v 1.1.2.7 2006/08/14 11:38:22 dpl Exp $

% For more details on message contents please refer to IS-GPS-800.


%% Check if the parameters are strings ==============================
if ~ischar(navBitsBin)
    error('The parameter BITS must be a character array!');
end
% 'bits' should be row vector for 'bin2dec' function.
[a, b] = size(navBitsBin);
if a > b
    navBitsBin = navBitsBin';
end

% Check if we have enough bits for a full L1C frame 
if length(navBitsBin) < 883
    return;
end

% Pi used in the GPS coordinate system 
gpsPi = 3.1415926535898;

%% ===== Initialization =============================================
% Note: Unlike legacy GPS LNAV, the L1C ephemeris message (Subframe 2) 
% does not contain the PRN number. The PRN is typically known by the 
% receiver channel or decoded from the almanac in Subframe 3.
% This function assumes eph is initialized.

if isfield(eph, 'flag') && isempty(eph.flag)
    
    %% ===== Decode Subframe 1: Time of Interval (TOI) ==================
    % Bits 1-9 (9 bits)
    % Represents the number of 18-second epochs since the start of the 
    % current 2-hour ITOW period.
    TOI = bin2dec(navBitsBin(1:9));
    
    %% ===== Decode Subframe 2: Ephemeris & Clock (Non-Variable) ========
    % Bits 10-609 (600 bits)
    sf2 = navBitsBin(10:609);
    
    % --- General Data Fields ---
    % Week Number (WN): Bits 1-13
    eph.WN = bin2dec(sf2(1:13));
    
    % Interval Time of Week (ITOW): Bits 14-21
    % The number of 2-hour epochs in the current week.
    eph.ITOW = bin2dec(sf2(14:21));
    
    % Data Sequence Propagation Time of Week (top): Bits 22-32
    eph.t_op = bin2dec(sf2(22:32)) * 300;
    
    % L1C Signal Health: Bit 33
    eph.L1CHealth = bin2dec(sf2(33));
    
    % URA_ED Index (Elevation Dependent User Range Accuracy): Bits 34-38
    eph.URAEDIndex = twosComp2dec(sf2(34:38));
    
    % Ephemeris Reference Time (toe): Bits 39-49
    % Unit: 300 seconds
    eph.t_oe = bin2dec(sf2(39:49)) * 300;
    % For L1C, t_oc (clock reference) is typically same as t_oe
    % eph.t_oc = eph.t_oe; 
    
    % --- Orbit Parameters (Ref: IS-GPS-800J Table 3.5-1) ---------------
    % Semi-major axis difference (Delta A): Bits 50-75
    eph.deltaA = twosComp2dec(sf2(50:75)) * 2^(-9);
    
    % Change rate in semi-major axis (Adot): Bits 76-100
    eph.ADot = twosComp2dec(sf2(76:100)) * 2^(-21);
    
    % Mean Motion difference (Delta n0): Bits 101-117
    % Unit: semi-circles/sec -> multiply by gpsPi to get rad/sec
    eph.delta_n_0 = twosComp2dec(sf2(101:117)) * 2^(-44) * gpsPi;
    
    % Rate of Mean Motion difference (Delta n0 dot): Bits 118-140
    eph.delta_n_0Dot = twosComp2dec(sf2(118:140)) * 2^(-57) * gpsPi;
    
    % Mean anomaly at reference time (M0): Bits 141-173
    eph.M_0 = twosComp2dec(sf2(141:173)) * 2^(-32) * gpsPi;
    
    % Eccentricity (e): Bits 174-206 (Unsigned)
    eph.e = bin2dec(sf2(174:206)) * 2^(-34);
    
    % Argument of perigee (omega): Bits 207-239
    eph.omega = twosComp2dec(sf2(207:239)) * 2^(-32) * gpsPi;
    
    % Inclination angle at reference time: Bits 240-272
    eph.omega_0 = twosComp2dec(sf2(240:272)) * 2^(-32) * gpsPi;
    
    % Rate of Right Ascension (OmegaDot): Bits 273-305
    eph.i_0 = twosComp2dec(sf2(273:305)) * 2^(-32) * gpsPi;
    
    % Rate of right ascension difference: Bits 306-322
    eph.omegaDot = twosComp2dec(sf2(306:322)) * 2^(-44) * gpsPi;
    
    % Rate of inclination angle : Bits 323-337
    eph.IDOT = twosComp2dec(sf2(323:337)) * 2^(-44) * gpsPi;
    
    % --- Harmonic Correction Terms -------------------------------------
    % NOTE: For L1C, these are in meters or radians directly. DO NOT multiply by PI.
    
    % Amplitude of sine harmonic correction to angle of inclination
    eph.C_is = twosComp2dec(sf2(338:353)) * 2^(-30); % Radians
    % Amplitude of cosine harmonic correction to angle of inclination
    eph.C_ic = twosComp2dec(sf2(354:369)) * 2^(-30); % Radians
    % Amplitude of sine correction term to the orbit radius
    eph.C_rs = twosComp2dec(sf2(370:393)) * 2^(-8);  % Meters
    % Amplitude of cosine correction term to the orbit radius
    eph.C_rc = twosComp2dec(sf2(394:417)) * 2^(-8);  % Meters
    % Amplitude of sine harmonic correction to argument of latitude
    eph.C_us = twosComp2dec(sf2(418:438)) * 2^(-30); % Radians
    % Amplitude of cosine harmonic correction to argument of latitude
    eph.C_uc = twosComp2dec(sf2(439:459)) * 2^(-30); % Radians
    
    % --- Accuracy & Clock Parameters -----------------------------------
    % URA NED (Non-Elevation Dependent) Indices
    eph.URANED0Index = twosComp2dec(sf2(460:464));
    eph.URANED1Index = bin2dec(sf2(465:467));
    eph.URANED2Index = bin2dec(sf2(468:470));
    
    % SV Clock Bias Correction Coefficient (af0): Bits 471-496
    eph.a_0 = twosComp2dec(sf2(471:496)) * 2^(-35);
    
    % SV Clock Drift Correction Coefficient (af1): Bits 497-516
    eph.a_1 = twosComp2dec(sf2(497:516)) * 2^(-48);
    
    % SV Clock Bias Correction Coefficient (af2): Bits 517-526
    eph.a_2 = twosComp2dec(sf2(517:526)) * 2^(-60);
    
    % Group Delay Differential (TGD): Bits 527-539
    eph.T_GD = twosComp2dec(sf2(527:539)) * 2^(-35);
    
    % Inter-Signal Corrections (ISC)
    % ISC for L1C pilot component
    eph.ISC_L1Cp = twosComp2dec(sf2(540:552)) * 2^(-35);
    % ISC for L1C data component
    eph.ISC_L1Cd = twosComp2dec(sf2(553:565)) * 2^(-35);
    
    % Integrity Status Flag (ISF)
    eph.ISF = bin2dec(sf2(566));
    
    % CEI Data Sequence Propagation Week Number (WN_op)
    eph.WN_op = bin2dec(sf2(567:574));
        
    %% ===== Decode Subframe 3: Variable Data ===========================
    % Bits 610-883 (274 bits)
    sf3 = navBitsBin(610:883);
    
    % Bits 1-8 of SF3: PRN (Useful if almanac, otherwise may be reserved)
    eph.PRN = bin2dec(sf3(1:8)); 
    
    % Page ID: Bits 9-14
    PageID = bin2dec(sf3(9:14));
    
    % Check Page ID to decode specific messages
 if PageID == 1 
    %% === Page 1: UTC & Ionosphere ===
    
    eph.PageID1 = 1;
    
    % --- UTC Parameters  ---
    % Bits 15-30: A0-n (16 bits) - UTC Bias coefficient
    eph.A0_n = twosComp2dec(sf3(15:30)) * 2^(-35); 
    
    % Bits 31-43: A1-n (13 bits) - UTC Drift coefficient
    eph.A1_n = twosComp2dec(sf3(31:43)) * 2^(-51); 
    
    % Bits 44-50: A2-n (7 bits) - UTC Drift rate coefficient
    eph.A2_n = twosComp2dec(sf3(44:50)) * 2^(-68); 
    
    % Bits 51-58: Delta t_LS (8 bits) 
    eph.delta_t_LS = twosComp2dec(sf3(51:58)); 
    
    % Bits 59-74: t_ot (16 bits) - UTC 
    eph.t_ot = bin2dec(sf3(59:74)) * 16; 
    
    % Bits 75-87: WN_ot (13 bits) - UTC 
    eph.WN_ot = bin2dec(sf3(75:87)); 
    
    % Bits 88-100: WN_LSF (13 bits) - 
    eph.WN_LSF = bin2dec(sf3(88:100)); 
    
    % Bits 101-104: DN (4 bits) - 
    eph.DN = bin2dec(sf3(101:104)); 
    
    % Bits 105-112: Delta t_LSF (8 bits) - 
    eph.delta_t_LSF = twosComp2dec(sf3(105:112)); 
    
    % --- Ionospheric Parameters  ---
    % Bits 113-120: Alpha 0 (8 bits)
    eph.alpha0 = twosComp2dec(sf3(113:120)) * 2^(-30);
    
    % Bits 121-128: Alpha 1 (8 bits)
    eph.alpha1 = twosComp2dec(sf3(121:128)) * 2^(-27);
    
    % Bits 129-136: Alpha 2 (8 bits)
    eph.alpha2 = twosComp2dec(sf3(129:136)) * 2^(-24);
    
    % Bits 137-144: Alpha 3 (8 bits)
    eph.alpha3 = twosComp2dec(sf3(137:144)) * 2^(-24);
    
    % Bits 145-152: Beta 0 (8 bits)
    eph.beta0 = twosComp2dec(sf3(145:152)) * 2^(11);
    
    % Bits 153-160: Beta 1 (8 bits)
    eph.beta1 = twosComp2dec(sf3(153:160)) * 2^(14);
    
    % Bits 161-168: Beta 2 (8 bits)
    eph.beta2 = twosComp2dec(sf3(161:168)) * 2^(16);
    
    % Bits 169-176: Beta 3 (8 bits)
    eph.beta3 = twosComp2dec(sf3(169:176)) * 2^(16);
    
    % --- ISC Parameters (Inter-Signal Correction) ---
    % Bits 177-189: ISC L1C/A (13 bits)
    eph.ISC_L1CA = twosComp2dec(sf3(177:189)) * 2^(-35);
    
    % Bits 190-202: ISC L2C (13 bits)
    eph.ISC_L2C = twosComp2dec(sf3(190:202)) * 2^(-35);
    
    % Bits 203-215: ISC L5I5 (13 bits)
    eph.ISC_L5I5 = twosComp2dec(sf3(203:215)) * 2^(-35);
    
    % Bits 216-228: ISC L5Q5 (13 bits)
    eph.ISC_L5Q5 = twosComp2dec(sf3(216:228)) * 2^(-35);
    
    % Bits 229-250: Reserved (22 bits) 
    
 elseif PageID == 2 
    %% === Page 2: GGTO & EOP 
    eph.PageID2 = 2;
    
    % --- GGTO (GPS/GNSS Time Offset) ---
    eph.GGTO.GNSS_ID = bin2dec(sf3(15:17));
    eph.GGTO.t_GGTO  = bin2dec(sf3(18:33)) * 16;       
    eph.GGTO.WN_GGTO = bin2dec(sf3(34:46));
    eph.GGTO.A0_GGTO = twosComp2dec(sf3(47:62)) * 2^(-35); 
    eph.GGTO.A1_GGTO = twosComp2dec(sf3(63:75)) * 2^(-51); 
    eph.GGTO.A2_GGTO = twosComp2dec(sf3(76:82)) * 2^(-68); 
    
    % --- EOP (Earth Orientation Parameters) ---
    eph.EOP.t_EOP = bin2dec(sf3(83:98)) * 16;          
    
    %  Bits 99-119 (21 bits)
    eph.EOP.PM_X = twosComp2dec(sf3(99:119)) * 2^(-20); 
    
    % PM_X_dot Bits 120-134 (15 bits)
    eph.EOP.PM_X_dot = twosComp2dec(sf3(120:134)) * 2^(-21);
    
    % PM_Y Bits 135-155 (21 bits)
    eph.EOP.PM_Y = twosComp2dec(sf3(135:155)) * 2^(-20);
    
    % PM_Y_dot Bits 156-170 (15 bits)
    eph.EOP.PM_Y_dot = twosComp2dec(sf3(156:170)) * 2^(-21);
    
    eph.EOP.Delta_UTGPS     = twosComp2dec(sf3(171:201)) * 2^(-24); 
    
    eph.EOP.Delta_UTGPS_dot = twosComp2dec(sf3(202:220)) * 2^(-25);
    
    % Bits 220-250: Reserved
end
if isempty(eph.flag)
    eph.TOW = eph.ITOW * 7200 + TOI * 18;
end
    % Set flag indicating successful decode
eph.flag = 1;
end
end
