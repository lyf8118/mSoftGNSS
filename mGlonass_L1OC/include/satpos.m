function [satPositions, satClkCorr] = satpos(transmitTime, prnList, eph)
%SATPOS Calculate GLONASS L1OC satellite positions and clock corrections.
%
%[satPositions, satClkCorr] = satpos(transmitTime, prnList, eph)
%
%   Inputs:
%       transmitTime  - signal transmission time, in seconds
%       prnList       - list of PRNs to be processed
%       eph           - ephemerides of satellites
%
%   Outputs:
%       satPositions  - satellite ECEF positions [X; Y; Z], in meters
%       satClkCorr    - satellite clock corrections, in seconds

%% Initialize constants ===================================================
numOfSatellites = length(prnList);

%--- Constants for satellite position calculation -------------------------
omega = 7.2921151467e-5;  % Earth rotation rate, [rad/s]
my    = 3.986004418e14;   % Earth's gravitational parameter, [m^3/s^2]
a     = 6.378136e6;       % Semi-major axis of Earth, [m]
J02   = 1.0826257e-3;     % Second zonal harmonic of the geopotential
daySec = 86400;           % Number of seconds in one day

%% Initialize results =====================================================
satClkCorr   = zeros(1, numOfSatellites);
satPositions = zeros(3, numOfSatellites);

%% Process each satellite =================================================
for satNr = 1 : numOfSatellites
    PRN = prnList(satNr);

    %%% Find the integration time -----------------------------------------
    % Transform the received signal time to the interval from the broadcast
    % ephemeris epoch according to Appendix D of the GLONASS CDMA ICD.
    timeDifference = transmitTime(satNr) + eph(PRN).Tau + ...
        eph(PRN).Tau_c - eph(PRN).Tb;
    dayNumber = sign(timeDifference) * ...
        floor(abs(timeDifference / daySec) + 0.5);
    timeDifference = timeDifference - ...
        dayNumber * daySec;
    deltaTb = timeDifference / ...
        (1 + eph(PRN).Gamma - eph(PRN).d_Tau_c);

    %%% Calculate the satellite clock correction --------------------------
    % Tau_c and d_Tau_c convert GLONASS system time to MT and therefore do
    % not form part of the satellite clock correction returned to PVT.
    satClkCorr(satNr) = -eph(PRN).Tau + eph(PRN).Gamma * deltaTb + ...
        eph(PRN).Beta * deltaTb^2;

    %%% Prepare the broadcast state vector --------------------------------
    % Convert the broadcast position, velocity and acceleration to SI units.
    state = [eph(PRN).X;  eph(PRN).Y;  eph(PRN).Z; ...
             eph(PRN).dX; eph(PRN).dY; eph(PRN).dZ] * 1e3;
    acceleration = [eph(PRN).ddX; eph(PRN).ddY; eph(PRN).ddZ] * 1e3;

    %%% Propagate the state vector with fourth-order Runge-Kutta -----------
    integratedTime = 0;
    while abs(deltaTb - integratedTime) > 1e-12
        step = sign(deltaTb - integratedTime) * ...
            min(30, abs(deltaTb - integratedTime));

        D1 = stateDerivative(state, acceleration, my, a, J02, omega);
        D2 = stateDerivative(state + 0.5 * step * D1, acceleration, ...
            my, a, J02, omega);
        D3 = stateDerivative(state + 0.5 * step * D2, acceleration, ...
            my, a, J02, omega);
        D4 = stateDerivative(state + step * D3, acceleration, ...
            my, a, J02, omega);

        state = state + step / 6 * (D1 + 2 * D2 + 2 * D3 + D4);
        integratedTime = integratedTime + step;
    end

    satPositions(:, satNr) = state(1:3);
end
end

%% GLONASS differential equations ========================================
function derivative = stateDerivative(state, acceleration, my, a, J02, ...
                                      omega)
%STATEDERIVATIVE Return the state derivative in the rotating ECEF frame.

x  = state(1);
y  = state(2);
z  = state(3);
Vx = state(4);
Vy = state(5);
Vz = state(6);
Ax = acceleration(1);
Ay = acceleration(2);
Az = acceleration(3);

radius = sqrt(x^2 + y^2 + z^2);
common = 1.5 * J02 * my * a^2 / radius^5;
dVx = -my * x / radius^3 - common * x * ...
    (1 - 5 * z^2 / radius^2) + omega^2 * x + 2 * omega * Vy + Ax;
dVy = -my * y / radius^3 - common * y * ...
    (1 - 5 * z^2 / radius^2) + omega^2 * y - 2 * omega * Vx + Ay;
dVz = -my * z / radius^3 - common * z * ...
    (3 - 5 * z^2 / radius^2) + Az;

derivative = [Vx; Vy; Vz; dVx; dVy; dVz];
end
