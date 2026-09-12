function eph = eph_structure_init()
% Initialize GLONASS L2OC ephemeris/status structure.
%
% Note:
%   The uploaded ICD "GLONASS CDMA L2, Edition 1.0, 2016" defines
%   L2OCp signal structure, PRN generator, OC2 and L2CSI/L2OCp TDM.
%   It does not define the L2 CSI navigation message format or ephemeris
%   bit fields. Therefore this structure is only a placeholder for
%   tracking-only L2OCp processing.

%--------------------------------------------------------------------------
% GLONASS L2OC signal metadata
%--------------------------------------------------------------------------
eph.system          = 'GLONASS';
eph.signal          = 'L2OCp';
eph.component       = 'pilot';
eph.navComponent    = 'L2CSI';
eph.navDecoded      = 0;

% Nominal L2 CDMA carrier frequency from ICD: 1248.06 MHz
eph.carrierFreq     = 1248.06e6;

% L2OCp PRN parameters from ICD
eph.L2OCpCodeLength = 10230;      % PRN_L2OCp chips
eph.L2OCpCodeRate   = 0.5115e6;   % chips/s
eph.L2OCpPeriod     = 0.020;      % 20 ms

% Equivalent TDM parameters used in your receiver implementation
% because L2CSI slot is set to zero and L2OCp/L2CSI are represented
% as an equivalent 20460-chip sequence.
eph.equivCodeLength = 20460;
eph.equivCodeRate   = 1.023e6;
eph.equivCodePeriod = 0.020;

%--------------------------------------------------------------------------
% GLONASS navigation/ephemeris placeholders
%--------------------------------------------------------------------------
% These fields are intentionally left empty because the current L2 ICD does
% not define the navigation message layout for L2 CSI.
eph.SVID        = [];
eph.PRN         = [];
eph.slot        = [];
eph.freqChannel = [];

% Time-related placeholders
eph.tow         = [];
eph.tk          = [];
eph.tb          = [];
eph.N4          = [];
eph.NT          = [];

% Clock correction placeholders
eph.tau_n       = [];
eph.gamma_n     = [];
eph.tau_c       = [];
eph.tau_gps     = [];

% Orbit state placeholders, GLONASS-style state vector form
eph.x           = [];
eph.y           = [];
eph.z           = [];
eph.xDot        = [];
eph.yDot        = [];
eph.zDot        = [];
eph.xDotDot     = [];
eph.yDotDot     = [];
eph.zDotDot     = [];

% Health/status placeholders
eph.health      = [];
eph.P           = [];
eph.P1          = [];
eph.P2          = [];
eph.P3          = [];
eph.P4          = [];

% Raw data/debug fields
eph.rawNavBits      = [];
eph.rawNavBitsLength = 0;
eph.decodeStatus    = 'L2OCp pilot tracking only; L2 CSI navigation message is not decoded.';

% Decoded flag
% 0: ephemeris not available
% 1: ephemeris decoded
eph.flag        = 0;

end