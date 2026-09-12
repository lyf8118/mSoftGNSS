function [navSolutions, eph] = postNavigation(trackResults, settings)
%Function calculates navigation solutions for the receiver (pseudoranges,
%positions). At the end it converts coordinates from the WGS84 system to
%the UTM, geocentric or any additional coordinate system.
%
%[navSolutions, eph] = postNavigation(trackResults, settings)
%
%   Inputs:
%       trackResults    - results from the tracking function (structure
%                       array).
%       settings        - receiver settings.
%   Outputs:
%       navSolutions    - contains measured pseudoranges, receiver
%                       clock error, receiver coordinates in several
%                       coordinate systems (at least ECEF and UTM).
%       eph             - received ephemerides of all SV (structure array).

%--------------------------------------------------------------------------
%                         CU Multi-GNSS SDR
% (C) Developed for GPS L1C SDR by Yafeng Li, Nagaraj C. Shivaramaiah 
% and Dennis M. Akos. 
% Based on the original SoftGNSS SDR framework by Darius Plausinaitis,
% Peter Rinder, Nicolaj Bertelsen and Dennis M. Akos
%
% Reference: Adapted within the CU Multi-GNSS SDR receiver framework for GPS L1C.
% implementation of an open-source L1C SDR receiver. 
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
%$Id: postNavigation.m,v 1.1.2.22 2006/08/09 17:20:11 dpl Exp $

%% Check is there enough data to obtain any navigation solution ===========
% GPS L1C frame length is 18 seconds (18000 ms). 
% 36000 ms ensures at least 2 full frames to guarantee one complete ephemeris.
if (settings.msToProcess < 36000)
    % Show the error message and exit
    disp('Record is too short . Exiting!');
    navSolutions = [];
    eph          = [];
    return
end

%% Pre-allocate space =======================================================
% Starting positions of the first message in the input bit stream.
% For L1C, this points to the start of the 1800-symbol frame (Subframe 1).
subFrameStart  = inf(1, settings.numberOfChannels);

% Time Of Week (TOW) of the first message(in seconds). 
TOW  = inf(1, settings.numberOfChannels);

%--- Make a list of channels excluding not tracking channels ---------------
activeChnList = find([trackResults.status] ~= '-');

for channelNr = activeChnList
    
    % Get PRN of current channel
    PRN = trackResults(channelNr).PRN;
    if settings.plotNavigation
        fprintf('Decoding CNAV2 for PRN %02d -------------------- \n', PRN);
    end
    
    %=== Decode ephemerides and TOW of the first sub-frame ==================
    % NOTE: Ensure 'CNAV2decoding' has been replaced with the L1C version!
    [eph(PRN), subFrameStart(channelNr), TOW(channelNr)] = ...
                                       CNAV2decoding(trackResults, channelNr, settings); %#ok<AGROW>
    
    %--- Exclude satellite if it does not have the necessary nav data -----
    if (eph(PRN).flag == 1)
        fprintf('    The requisite messages for PRN %02d all decoded!\n', PRN);
    else
        %--- Exclude channel from the list (from further processing) ------
        activeChnList = setdiff(activeChnList, channelNr);
        fprintf('    Ephemeris decoding fails for PRN %02d!!!\n', PRN);
    end
end

%% Check if the number of satellites is still above 3 =====================
if (isempty(activeChnList) || (size(activeChnList, 2) < 4))
    % Show error message and exit
    disp('Too few satellites with ephemeris data for position calculations. Exiting!');
    navSolutions = [];
    if ~exist('eph','var')
        eph          = [];
    end
    return
end

%% Set measurement-time point and step  =====================================
% Find start and end of measurement point locations in IF signal stream
sampleStart = zeros(1, settings.numberOfChannels);
sampleEnd = inf(1, settings.numberOfChannels);

for channelNr = activeChnList
    % 'subFrameStart' is the index in the tracking results
    sampleStart(channelNr) = ...
        trackResults(channelNr).absoluteSample(subFrameStart(channelNr));
    
    sampleEnd(channelNr) = trackResults(channelNr).absoluteSample(end);
end

% Add margin to avoid index exceeding dimensions
sampleStart = max(sampleStart) + 1;
sampleEnd = min(sampleEnd) - 1;

%--- Measurement step in unit of IF samples -------------------------------
measSampleStep = fix(settings.samplingFreq * settings.navSolPeriod/1000);

%--- Number of measurement points -----------------------------------------
measNrSum = fix((sampleEnd - sampleStart) / measSampleStep);

%% Initialization =========================================================
% Set satellite elevations to INF for the first iteration
satElev  = inf(1, settings.numberOfChannels);

% Set local time to inf for first calculation.
localTime = inf;

%##########################################################################
%#       Do the satellite and receiver position calculations              #
%##########################################################################
fprintf('Positions are being computed. Please wait... \n');

for currMeasNr = 1:measNrSum
    
    fprintf('Fix: Processing %02d of %02d \n', currMeasNr, measNrSum);
    
    %% Initialization of current measurement ==============================
    % Exclude satellites below elevation mask
    activeChnList = intersect(find(satElev >= settings.elevationMask), ...
        activeChnList);
    
    % Save list of satellites used
    navSolutions.PRN(activeChnList, currMeasNr) = ...
        [trackResults(activeChnList).PRN];
    
    % Initialize Azimuth/Elevation arrays with NaN
    navSolutions.el(:, currMeasNr) = NaN(settings.numberOfChannels, 1);
    navSolutions.az(:, currMeasNr) = NaN(settings.numberOfChannels, 1);
    
    navSolutions.transmitTime(:, currMeasNr) = ...
        NaN(settings.numberOfChannels, 1);
    navSolutions.satClkCorr(:, currMeasNr) = ...
        NaN(settings.numberOfChannels, 1);
    
    % Current sample index
    currMeasSample = sampleStart + measSampleStep*(currMeasNr-1);
    
    %% Find pseudoranges ======================================================
    % Raw pseudorange = (localTime - transmitTime) * c
    [navSolutions.rawP(:, currMeasNr), transmitTime, localTime]=  ...
        calculatePseudoranges(trackResults, subFrameStart, TOW, ...
        currMeasSample, localTime, activeChnList, settings);
    
    % Save transmitTime
    navSolutions.transmitTime(activeChnList, currMeasNr) = ...
        transmitTime(activeChnList);
    
    %% Find satellites positions and clocks corrections =======================
    % Ensure your 'satpos' function supports L1C ephemeris structure!
    % Specifically, check handling of 'deltaA' and L1C reference A.
    [satPositions, satClkCorr] = satpos(transmitTime(activeChnList), ...
                                        [trackResults(activeChnList).PRN], eph); 
                                                                      
    % Save satClkCorr
    navSolutions.satClkCorr(activeChnList, currMeasNr) = satClkCorr;
    
    %% Find receiver position =================================================
    if size(activeChnList, 2) > 3
        
        %=== Calculate receiver position ==================================
        % Correct pseudorange for SV clock error
        clkCorrRawP = navSolutions.rawP(activeChnList, currMeasNr)' + ...
            satClkCorr * settings.c;
        
        % Calculate receiver position (Least Squares)
        [xyzdt, navSolutions.el(activeChnList, currMeasNr), ...
            navSolutions.az(activeChnList, currMeasNr), ...
            navSolutions.DOP(:, currMeasNr)] = ...
            leastSquarePos(satPositions, clkCorrRawP, settings);
        
        %=== Save results =================================================
        navSolutions.X(currMeasNr)  = xyzdt(1);
        navSolutions.Y(currMeasNr)  = xyzdt(2);
        navSolutions.Z(currMeasNr)  = xyzdt(3);
        
        % Clock bias update
        if (currMeasNr == 1)
            navSolutions.dt(currMeasNr) = 0; 
        else
            navSolutions.dt(currMeasNr) = xyzdt(4);
        end
        
        %=== Correct local time by clock error estimation =================
        localTime = localTime - xyzdt(4)/settings.c;
        navSolutions.localTime(currMeasNr) = localTime;
        
        navSolutions.currMeasSample(currMeasNr) = currMeasSample;
        
        % Update satellite elevations
        satElev = navSolutions.el(:, currMeasNr)';
        
        %=== Correct pseudorange measurements =============================
        navSolutions.correctedP(activeChnList, currMeasNr) = ...
            navSolutions.rawP(activeChnList, currMeasNr) + ...
            satClkCorr' * settings.c - xyzdt(4);
        
        %% Coordinate conversion ==========================================
        
        %=== Convert to geodetic coordinates (WGS84) ======================
        [navSolutions.latitude(currMeasNr), ...
            navSolutions.longitude(currMeasNr), ...
            navSolutions.height(currMeasNr)] = cart2geo(...
            navSolutions.X(currMeasNr), ...
            navSolutions.Y(currMeasNr), ...
            navSolutions.Z(currMeasNr), ...
            5); % 5 usually denotes WGS84 ellipsoid
        
        %=== Convert to UTM ===============================================
        navSolutions.utmZone = findUtmZone(navSolutions.latitude(currMeasNr), ...
            navSolutions.longitude(currMeasNr));
        
        [navSolutions.E(currMeasNr), ...
            navSolutions.N(currMeasNr), ...
            navSolutions.U(currMeasNr)] = cart2utm(xyzdt(1), xyzdt(2), ...
            xyzdt(3), ...
            navSolutions.utmZone);
        
    else
        %--- Not enough satellites ----------------------------------------
        disp(['   Measurement No. ', num2str(currMeasNr), ...
            ': Not enough information for position solution.']);
        
        navSolutions.X(currMeasNr)           = NaN;
        navSolutions.Y(currMeasNr)           = NaN;
        navSolutions.Z(currMeasNr)           = NaN;
        navSolutions.dt(currMeasNr)          = NaN;
        navSolutions.DOP(:, currMeasNr)      = zeros(5, 1);
        navSolutions.latitude(currMeasNr)    = NaN;
        navSolutions.longitude(currMeasNr)   = NaN;
        navSolutions.height(currMeasNr)      = NaN;
        navSolutions.E(currMeasNr)           = NaN;
        navSolutions.N(currMeasNr)           = NaN;
        navSolutions.U(currMeasNr)           = NaN;
        
        navSolutions.az(activeChnList, currMeasNr) = ...
            NaN(1, length(activeChnList));
        navSolutions.el(activeChnList, currMeasNr) = ...
            NaN(1, length(activeChnList));
        
        disp('   Exit Program')
        return
        
    end 
    
    %=== Update local time for next step ==================================
    localTime = localTime + measSampleStep/settings.samplingFreq ;
    
end 
end
