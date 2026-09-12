function [satPositions, satClkCorr] = satpos(transmitTime, prnList, eph)
%SATPOS_L3OC Calculation of GLONASS L3OCd satellite coordinates.
%   Inputs:
%       transmitTime  - transmission time (GLONASS Time)
%       prnList       - list of PRNs
%       eph           - ephemeris structure
%   Outputs:
%       satPositions  - [3 x N] ECEF coordinates (X; Y; Z) in METERS
%       satClkCorr    - [1 x N] Clock corrections in SECONDS

%% Initialize constants ===================================================
    numOfSatellites = length(prnList);
    satPositions = zeros(3, numOfSatellites);
    satClkCorr   = zeros(1, numOfSatellites);
    
 %---GLONASS Constants-----------------------------------------------------
    ae = 6378136.0;         % Semi-major axis [m]
    mu = 398600.4418e9;     % Gravitational constant [m^3/s^2]
    J2 = 1082.6257e-6;      % Second zonal harmonic
    omega = 7.292115e-5;    % Earth rotation rate [rad/s]
    DaySec = 86400;

%% Process each satellite =================================================
    for i = 1 : numOfSatellites
        prn = prnList(i);
        
        if isempty(eph(prn).X) || isempty(eph(prn).Tb)
            continue; 
        end
      
        %Calculate Time Difference
        tb = eph(prn).Tb;     % In Moscow Time
        t  = transmitTime(i); % In GLONASS Time
        
       % Adjust for day rollover
        dt = t - tb;
        if dt > DaySec/2
            dt = dt - DaySec;
        elseif dt < -DaySec/2
            dt = dt + DaySec;
        end
        
        % Clock Correction
        satClkCorr(i) = -eph(prn).Tau + eph(prn).Gamma * dt;
        
        %State Vector Integration
        r0 = [eph(prn).X; eph(prn).Y; eph(prn).Z] * 1000;
        v0 = [eph(prn).dX; eph(prn).dY; eph(prn).dZ] * 1000;
        
        % Perturbing accelerations 
        ls_acc = [eph(prn).ddX; eph(prn).ddY; eph(prn).ddZ] * 1000;
        
        % Integration step size
        step = 30; 
        if abs(dt) < step, step = dt; end
        
        % Current state
        r = r0;
        v = v0;
        t_current = 0;
        
        % Numerical Integration Loop
        while abs(t_current) < abs(dt)
            if (abs(dt) - abs(t_current)) < abs(step)
                h = sign(dt) * (abs(dt) - abs(t_current));
            else
                h = sign(dt) * abs(step);
            end
            
            [k1_r, k1_v] = diff_eq(r, v, ls_acc, mu, ae, J2, omega);
            [k2_r, k2_v] = diff_eq(r + 0.5*h*k1_r, v + 0.5*h*k1_v, ls_acc, mu, ae, J2, omega);
            [k3_r, k3_v] = diff_eq(r + 0.5*h*k2_r, v + 0.5*h*k2_v, ls_acc, mu, ae, J2, omega);
            [k4_r, k4_v] = diff_eq(r + h*k3_r, v + h*k3_v, ls_acc, mu, ae, J2, omega);
            
            r = r + (h/6) * (k1_r + 2*k2_r + 2*k3_r + k4_r);
            v = v + (h/6) * (k1_v + 2*k2_v + 2*k3_v + k4_v);
            
            t_current = t_current + h;
        end
        
        % Satellite Coordinates at Signal Transmit Time
        satPositions(:, i) = r;
    end
end

%% GLONASS Differential Equations =================================================
function [dr, dv] = diff_eq(r, v, a_ls, mu, ae, J2, omega)
    x = r(1); y = r(2); z = r(3);
    vx = v(1); vy = v(2); vz = v(3);
    
    r_sq = x^2 + y^2 + z^2;
    r_mag = sqrt(r_sq);
    
    % Central Body + J2 Term factors
    c1 = -mu / r_sq / r_mag;
    c2 = 1.5 * J2 * mu * (ae^2) / (r_sq * r_mag * r_sq);
    z_rat = 5 * z^2 / r_sq;
    
    ax = c1*x - c2*x*(1 - z_rat) + omega^2*x + 2*omega*vy + a_ls(1);
    ay = c1*y - c2*y*(1 - z_rat) + omega^2*y - 2*omega*vx + a_ls(2);
    az = c1*z - c2*z*(3 - z_rat) + a_ls(3);
    
    dr = [vx; vy; vz];
    dv = [ax; ay; az];
end
