function [satPositions, satClkCorr] = satpos(transmitTime, prnList, eph)
%SATPOS Calculation of X,Y,Z satellites coordinates at TRANSMITTIME for
%given ephemeris EPH. Coordinates are calculated for each satellite in the
%list PRNLIST.
%[satPositions, satClkCorr] = satpos(transmitTime, prnList, eph, settings);
%
%   Inputs:
%       transmitTime  - transmission time
%       prnList       - list of PRNs to be processed
%       eph           - ephemerides of satellites
%       settings      - receiver settings
%
%   Outputs:
%       satPositions  - positions of satellites (in ECEF system [X; Y; Z;])
%       satClkCorr    - correction of satellites clocks

%--------------------------------------------------------------------------
%                           SoftGNSS v3.0
%--------------------------------------------------------------------------
%Based on Kai Borre 04-09-96
%Copyright (c) by Kai Borre
%Updated by Darius Plausinaitis, Peter Rinder and Nicolaj Bertelsen
% Modified for GPS L1C by Yafeng Li
%
% CVS record:
% $Id: satpos.m,v 1.1.2.15 2006/08/22 13:45:59 dpl Exp $

%% Initialize constants ===================================================
numOfSatellites = size(prnList, 2);

% GPS Constants (IS-GPS-800J Table 3.2-1)
gpsPi        = 3.1415926535898;        % Pi
OmegaE       = 7.2921151467e-5;        % Earth rotation rate, [rad/s]
mu           = 3.986005e14;            % Earth's gravitational constant [m^3/s^2]
F            = -4.442807633e-10;       % Relativistic correction constant [s/m^1/2]

% Semi-major axis reference for GPS L1C (IS-GPS-800J Table 3.5-1)
A_ref_L1C    = 26559710;               % [m]%%%%%%%%%%
Omega_dot_ref = -2.6e-9 * gpsPi;       % [rad/s]
%% Initialize results ===============================================
satClkCorr   = zeros(1, numOfSatellites);
satPositions = zeros(3, numOfSatellites);

%% Process each satellite =================================================
for satNr = 1 : numOfSatellites
    
    prn = prnList(satNr);
    
    %% Find initial satellite clock correction --------------------------------
    % Note: For L1C CNAV-2, toe is used as the reference time for clock (toc=toe)
    
    %--- Find time difference ---------------------------------------------
    % check_t handles week crossovers (adjusts result to [-302400, 302400])
    dt = check_t(transmitTime(satNr) - eph(prn).t_oe);
    
    %--- Calculate clock correction ---------------------------------------
    % Formula: dt_sv = af0 + af1*dt + af2*dt^2 - T_GD + ISC
    % Assuming tracking L1C Pilot component: use ISC_L1Cp
    % If tracking Data component, change to eph(prn).ISC_L1Cd
    
    satClkCorr(satNr) = (eph(prn).a_2 * dt + eph(prn).a_1) * dt + ...
        eph(prn).a_0 - ...
        eph(prn).T_GD + eph(prn).ISC_L1Cp;%%%%%%%%%
    
    time = transmitTime(satNr) - satClkCorr(satNr);
    
    %% Find satellite's position ----------------------------------------------
    % Time correction relative to Ephemeris Reference Time
    tk = check_t(time - eph(prn).t_oe);
    
    % --- Semi-major axis (A) ---
    % A0 = A_ref + deltaA
    A0   = A_ref_L1C + eph(prn).deltaA;%%%%%%%%%
    % A = A0 + Adot * tk
    A    = A0 + eph(prn).ADot * tk;
    
    % --- Mean Motion (n) ---
    n0  = sqrt(mu / (A0^3));
    
    % Corrected mean motion (including rate of change term for CNAV-2)
    delta_n = eph(prn).delta_n_0 + 0.5 * eph(prn).delta_n_0Dot * tk;
    n = n0 + delta_n;
    
    % --- Mean anomaly (M) ---
    M = eph(prn).M_0 + n * tk;
    
    % Reduce mean anomaly to between 0 and 2*pi
    M = rem(M + 2*gpsPi, 2*gpsPi);
    
    % --- Eccentric Anomaly (E) by Iteration ------------------------------
    E = M;
    for ii = 1:10
        E_old   = E;
        E       = M + eph(prn).e * sin(E);
        dE      = rem(E - E_old, 2*gpsPi);
        
        if abs(dE) < 1.e-12
            break;
        end
    end
    E = rem(E + 2*gpsPi, 2*gpsPi);
    
    % --- Relativistic correction term ------------------------------------
    % IS-GPS-800J uses sqrt(A) in the relativistic formula
    dtr = F * eph(prn).e * sqrt(A) * sin(E);%%%%%%%%%%%
    
    % --- True Anomaly (nu) -----------------------------------------------
    nu   = atan2(sqrt(1 - eph(prn).e^2) * sin(E), cos(E)-eph(prn).e);
    
    % --- Argument of Latitude (Phi) --------------------------------------
    Phi = nu + eph(prn).omega;
    Phi = rem(Phi, 2*gpsPi);
    
    % --- Second Harmonic Perturbations -----------------------------------
    sin2Phi = sin(2 * Phi);
    cos2Phi = cos(2 * Phi);
    
    % Argument of Latitude Correction
    u = Phi + ...
        eph(prn).C_uc * cos2Phi + ...
        eph(prn).C_us * sin2Phi;
    % Radius Correction
    r = A * (1 - eph(prn).e*cos(E)) + ...
        eph(prn).C_rc * cos2Phi + ...
        eph(prn).C_rs * sin2Phi;
    % Inclination Correction
    i = eph(prn).i_0 + eph(prn).IDOT * tk + ...
        eph(prn).C_ic * cos2Phi + ...
        eph(prn).C_is * sin2Phi;
    
    % --- Longitude of Ascending Node (Omega) -----------------------------

    OmegaDot_Corrected = Omega_dot_ref + eph(prn).omegaDot;
    
    Omega = eph(prn).omega_0 + (OmegaDot_Corrected - OmegaE)*tk - ...
        OmegaE * eph(prn).t_oe;
        
    Omega = rem(Omega + 2*gpsPi, 2*gpsPi);  %%%%%%%%%
    
    % --- Compute satellite coordinates (ECEF) ----------------------------
    xp = r * cos(u);
    yp = r * sin(u);
    
    satPositions(1, satNr) = xp * cos(Omega) - yp * cos(i) * sin(Omega);
    satPositions(2, satNr) = xp * sin(Omega) + yp * cos(i) * cos(Omega);
    satPositions(3, satNr) = yp * sin(i);
    
    %% Final Clock Correction ---------------------------------------------
    % Include relativistic correction
    satClkCorr(satNr) = (eph(prn).a_2 * dt + eph(prn).a_1) * dt + ...
                         eph(prn).a_0 - eph(prn).T_GD + eph(prn).ISC_L1Cp + dtr;
    
end % for satNr
end