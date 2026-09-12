function eph = eph_structure_init()
% Initialize the Ephemeris data structure
%--------------------------------------------------------------------------
%                         CU Multi-GNSS SDR  
%--------------------------------------------------------------------------

% Initialize with NaN (Not a Number) or 0 to allow matrix operations
% and prevent "dimension mismatch" errors in satpos.m

%--- SVID , Clock correction, SISA, Ionospheric correction ----------------
eph.SVID      = NaN;      % 6 bits
eph.IODnav1   = NaN;      % 10 bits
eph.t_oc      = NaN;      % 14 bits in s
eph.a_f0      = 0;        % 31 bit in s (Initialize to 0 to prevent crash)
eph.a_f1      = 0;        % 21 bit in s/s
eph.a_f2      = 0;        % 6 bit in s/(s^2)
eph.a_i0      = 0;        % 11 bit
eph.a_i1      = 0;        % 11 bit
eph.a_i2      = 0;        % 14 bit
eph.iono_SF1  = 0;        % 1 bit
eph.iono_SF2  = 0;        % 1 bit
eph.iono_SF3  = 0;        % 1 bit
eph.iono_SF4  = 0;        % 1 bit
eph.iono_SF5  = 0;        % 1 bit
eph.BGD_E1E5a = 0;        % 10 bit
eph.E5a_HS    = NaN;      % 2 bit
eph.WN        = NaN;      % 12 bit
eph.TOW       = NaN;      % 20 bit
eph.E5a_DVS   = NaN;      % 1 bit

%--- Ephemeris (1/3) and GST -------------------------------------
eph.IODnav2   = NaN;      % 10 bit
eph.M_0       = 0;        % 32 bit in rad
eph.OmegaDot  = 0;        % 24 bit in rad
eph.e         = 0;        % 32 bit
eph.sqrtA     = 0;        % 32 bit in m^0.2
eph.Omega_0   = 0;        % 32 bit in rad
eph.iDot      = 0;        % 14 bit in rad

%--- Ephemeris (2/3) and GST -------------------------------
eph.IODnav3   = NaN;      % 10 bit
eph.i_0       = 0;        % 32 bit in rad
eph.omega     = 0;        % 32 bit in rad
eph.deltan    = 0;        % 16 bit in rad
eph.CUC       = 0;        % 16 bit in rad
eph.CUS       = 0;        % 16 bit in rad
eph.CRC       = 0;        % 16 bit in m
eph.CRS       = 0;        % 16 bit in m
eph.t_oe      = NaN;      % 14 bit in s

%--- Ephemeris (3/3) ----------------------------------------
eph.IODnav4   = NaN;      % 10 bit
eph.CIC       = 0;        % 16 bit in rad
eph.CIS       = 0;        % 16 bit in rad
eph.A0        = 0;        % 32 bit
eph.A1        = 0;        % 24 bit
eph.delt_LS   = 0;        % 8 bit
eph.t_ot      = 0;        % 8 bit
eph.WN_ot     = 0;        % 8 bit
eph.WN_LSF    = 0;        % 8 bit
eph.DN        = 0;        % 3 bit
eph.delt_LSF  = 0;        % 8 bit
eph.t_og      = 0;        % 8 bit
eph.A0_G      = 0;        % 16 bit
eph.A1_G      = 0;        % 12 bit
eph.WN_og     = 0;        % 6 bit

%--- HAS Data Buffer (E6-B Specific) -------------------------
eph.HAS = [];             % Storage for decoded HAS data
eph.HAS_Buffer = struct();

% ephemeris decoded flag --------------------------------------------------
eph.flag         = 0;