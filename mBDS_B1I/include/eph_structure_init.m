function eph = eph_structure_init ()
%--------------------------------------------------------------------------
%                         CU Multi-GNSS SDR
% (C) Developed for BDS B1I/B2I SDR by Yafeng Li, Daehee Won, 
% Nagaraj C. Shivaramaiah and Dennis M. Akos. 
% Based on the original framework for GPS C/A SDR by Darius Plausinaitis,
% Peter Rinder, Nicolaj Bertelsen and Dennis M. Akos
%--------------------------------------------------------------------------
% Initialize the Ephemeris data before decoding.

eph.SatH1  = [];
eph.IODC   = [];
eph.URAI   = [];
eph.WN     = [];

eph.t_oc   = [];
eph.T_GD_1 = [];
eph.T_GD_2 = [];

% Ionospheric Delay Model Parameters (alpha, beta)
eph.alpha0 = [];
eph.alpha1 = [];
eph.alpha2 = [];
eph.alpha3 = [];

eph.beta0  = [];
eph.beta1  = [];
eph.beta2  = [];
eph.beta3  = [];

% Clock Correction Parameters
eph.a2     = [];
eph.a0     = [];
eph.a1     = [];

% Issue of Data, Ephemeris (IODE)
eph.IODE   = [];

% Ephemeris Parameters
eph.deltan = [];
eph.C_uc   = [];
eph.M_0    = [];
eph.e      = [];

eph.C_us   = [];
eph.C_rc   = [];
eph.C_rs   = [];

eph.sqrtA  = [];
eph.i_0      = [];
eph.C_ic     = [];
eph.omegaDot = [];
eph.C_is     = [];
eph.iDot     = [];
eph.omega_0  = [];
eph.omega    = [];
eph.t_oe    = [];
% Tow of first decoded subframe
eph.SOW         = [];
% ephemeris decoded flag
eph.flag         = 0;
