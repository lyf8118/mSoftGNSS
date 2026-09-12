function eph = eph_structure_init()
% This is in order to make sure variable 'eph' for each SV has a similar 
% structure when only one or even none of the three requisite messages
% is decoded for a given PRN.
%--------------------------------------------------------------------------
%                         CU Multi-GNSS SDR  
% (C) Written by Yafeng Li, Nagaraj C. Shivaramaiah and Dennis M. Akos
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
    %% Identification & Time  =============================================
    % Satellite ID 
    eph.SV_ID   = []; 
    % Time of Day
    eph.TOD     = inf; 
    % Time of Week 
    eph.TOW     = inf;      
    % String Type identifier
    eph.Type    = [];       
    
    %% Message type 10 ====================================================
    % 4-year interval number
    eph.N4      = [];      
    % Day number within 4-year interval
    eph.NT      = [];       
    % SV modification flag
    eph.M       = [];       
    % Index time [s]
    eph.Tb      = [];       
    % SV Clock Bias [s]
    eph.Tau     = [];       
    % Frequency Offset
    eph.Gamma   = [];       
    % Frequency Drift Rate  [s^-1]
    eph.Beta    = [];       
    % GLONASS to UTC(SU) time correction [s]
    eph.Tau_c   = [];       
    % Rate of change of Tau_c [dimensionless]
    eph.d_Tau_c = [];       
    % Age of Ephemeris data
    eph.E_E     = [];      
    % Age of Time/Clock data
    eph.E_T     = [];       
    % User Range Accuracy (URA) index for Ephemeris
    eph.F_E     = [];       
    % URA index for Time
    eph.F_T     = [];       
    % Health flag (H^j)
    eph.Health  = [];       
    % Data Validity flag (l^j)
    eph.DataValid = [];     
    % Flag P1 (Service field)
    eph.P1      = [];       
    % Flag P2 (Sun-pointing / Maneuver)
    eph.P2      = [];       

    %% Message type 11&12 =================================================
    % Position at instant Tb 
    % X coordinate [km]
    eph.X       = []; 
    % Y coordinate [km]
    eph.Y       = [];      
    % Z coordinate [km]
    eph.Z       = [];       
    % Velocity at instant Tb
    % Velocity X component [km/s]
    eph.dX      = [];       
    % Velocity Y component [km/s]
    eph.dY      = [];       
    % Velocity Z component [km/s]
    eph.dZ      = [];       
    % Luni-Solar Acceleration at instant Tb
    % Acceleration X component [km/s^2]
    eph.ddX     = [];       
    % Acceleration Y component [km/s^2]
    eph.ddY     = [];       
    % Acceleration Z component [km/s^2]
    eph.ddZ     = [];       
    
    %% Message type 12 ====================================================
    % Antenna Phase Center offsets
    eph.Delta_X_pc = [];   
    eph.Delta_Y_pc = [];
    eph.Delta_Z_pc = [];
    % L3OCp to L3OCd time offset
    eph.Delta_Tau_L3 = []; 
    % GPS to GLONASS time offset correction
    eph.Tau_GPS      = [];  
    %% L1OCd service fields ===============================================
    eph.KP      = [];
    eph.A       = [];

    %% L1OCd additional immediate-data fields =============================
    eph.PS      = [];
    eph.R_E     = [];
    eph.R_T     = [];
    eph.DeltaXpc = [];
    eph.DeltaYpc = [];
    eph.DeltaZpc = [];
    eph.DeltaTauL2 = [];
    eph.TauGPS     = [];

    %% Decoding Status ====================================================
    eph.flag    = 0;        % 1 if fully decoded
    
end