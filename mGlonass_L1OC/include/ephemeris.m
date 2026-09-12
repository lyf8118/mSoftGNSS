function [eph, TOD] = ephemeris(navBits, eph)
%% Initialization ==========================================================
if nargin < 2 || isempty(eph)
    eph = eph_structure_init();
end

% Ensure input is a column vector.
navBits = navBits(:);

% Check length.
if length(navBits) ~= 250
    error('GLONASS L1OCd string must be exactly 250 bits long.');
end

%% Parse Common Service Fields ============================================
% Type.
typeVal = bin2dec_unsigned(navBits(13:18));

% TS. The L1OCd timestamp unit is 2 seconds.
tsVal = bin2dec_unsigned(navBits(35:50));
TOD = tsVal * 2;
eph.TOD = TOD;
eph.TOW = TOD;

% Satellite ID.
eph.SV_ID = bin2dec_unsigned(navBits(19:24));

% Health and data validity.
eph.Health = navBits(25);
eph.DataValid = navBits(26);

% Service fields.
eph.Type = typeVal;
eph.P1   = bin2dec_unsigned(navBits(27:30));
eph.P2   = navBits(31);
eph.KP   = bin2dec_unsigned(navBits(32:33));
eph.A    = navBits(34);

%% Parse Data Fields based on String Type =================================
switch typeVal

    case 10
        %--- It is String Type 10 -----------------------------------------
        eph.N4      = bin2dec_unsigned(navBits(51:55));
        eph.NT      = bin2dec_unsigned(navBits(56:66));
        eph.M       = bin2dec_unsigned(navBits(67:69));
        eph.PS      = bin2dec_unsigned(navBits(70:75));

        % tb, LSB: 90 s.
        eph.Tb      = bin2dec_unsigned(navBits(76:85)) * 90;

        eph.E_E     = bin2dec_unsigned(navBits(86:93));
        eph.E_T     = bin2dec_unsigned(navBits(94:101));
        eph.R_E     = bin2dec_unsigned(navBits(102:103));
        eph.R_T     = bin2dec_unsigned(navBits(104:105));
        eph.F_E     = bin2dec_signed(navBits(106:110));
        eph.F_T     = bin2dec_signed(navBits(111:115));

        % Clock and time correction parameters.
        eph.Tau     = bin2dec_signed(navBits(116:147)) * 2^(-38);
        eph.Gamma   = bin2dec_signed(navBits(148:166)) * 2^(-48);
        eph.Beta    = bin2dec_signed(navBits(167:181)) * 2^(-57);
        eph.Tau_c   = bin2dec_signed(navBits(182:221)) * 2^(-31);
        eph.d_Tau_c = bin2dec_signed(navBits(222:234)) * 2^(-49);

    case 11
        %--- It is String Type 11 -----------------------------------------
        % Position at instant Tb.
        eph.X       = bin2dec_signed(navBits(51:90))  * 2^(-20);
        eph.Y       = bin2dec_signed(navBits(91:130)) * 2^(-20);
        eph.Z       = bin2dec_signed(navBits(131:170)) * 2^(-20);

        % Velocity X component at instant Tb.
        eph.dX      = bin2dec_signed(navBits(171:205)) * 2^(-30);

        % Phase-center corrections.
        eph.DeltaXpc = bin2dec_signed(navBits(206:218)) * 2^(-10);
        eph.DeltaYpc = bin2dec_signed(navBits(219:231)) * 2^(-10);

    case 12
        %--- It is String Type 12 -----------------------------------------
        % Phase-center correction.
        eph.DeltaZpc = bin2dec_signed(navBits(51:63)) * 2^(-10);

        % Velocity Y and Z components at instant Tb.
        eph.dY      = bin2dec_signed(navBits(64:98))  * 2^(-30);
        eph.dZ      = bin2dec_signed(navBits(99:133)) * 2^(-30);

        % Luni-solar acceleration components at instant Tb.
        eph.ddX     = bin2dec_signed(navBits(134:148)) * 2^(-39);
        eph.ddY     = bin2dec_signed(navBits(149:163)) * 2^(-39);
        eph.ddZ     = bin2dec_signed(navBits(164:178)) * 2^(-39);

        % Inter-frequency and GPS time correction parameters.
        eph.DeltaTauL2 = bin2dec_signed(navBits(179:196)) * 2^(-38);
        eph.TauGPS     = bin2dec_signed(navBits(197:226)) * 2^(-38);

    otherwise
        % Strings other than Type 10, 11 and 12 are ignored here.
end

end

%=== Helper functions ======================================================
function val = bin2dec_unsigned(bits)
    bits = bits(:)';
    val = bin2dec(char(bits + '0'));
end

function val = bin2dec_signed(bits)
    bits = bits(:)';

    if bits(1) == 0
        val = bin2dec(char(bits + '0'));
    else
        val = -bin2dec(char(bits(2:end) + '0'));
    end
end