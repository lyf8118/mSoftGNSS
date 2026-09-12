function [eph] = ephemeris(navBits,eph)

%--------------------------------------------------------------------------
%                         CU Multi-GNSS SDR
% (C) Written by Yafeng Li, Nagaraj C. Shivaramaiah and Dennis M. Akos
%--------------------------------------------------------------------------

%% Prepare related variables ========================================

% This is 15 full pages (one subframe)
if length(navBits) < 2500
    disp('The ephemeris bitstream must contain at least 2500 symbols!');
    return;
end

% Check polarity of the data bits
sync_bits = [1 0 1 1 0 1 1 1 0 0 0 0];

% Convert CNAV-producing convolutional code polynomials to trellis description.
% Note that the difference from GPS is that the second branch G2 is
% inverted at the end (see ICD)
trellis = poly2trellis(7,[171 133]);

%Creates a cyclic redundancy code (CRC) detector System object
crcDet = comm.CRCDetector([24 23 18 17 14 11 10 7 6 5 4 3 1 0]);

% Viterbi traceback depth for vitdec(function)
tblen = 35;

%Indicator that the requisite messages are all decoded
word_valid = inf(1,4);
%% Define required constants ==============================================
% Pi used in the Galileo coordinate system (same as for GPS)
galPi = 3.1415926535898;
% Initialize the eph structure
eph = eph_structure_init();

%% Decode Messages ========================================================
for ii = 1:5
    
    %--- Pull out bits for both page parts --------------------------------
    pagePart = navBits(500*(ii-1)+1 : 500* ii   )';
    
    %--- Correct polarity of the all data bits according to preamble bits
    if(~isequal(pagePart(1:12),sync_bits))
        pageSymInt = not(pagePart(13:500));
    else
        pageSymInt = pagePart(13:500);
    end
    
    %--- De-interleave page parts -----------------------------------------
    symMat = reshape(pageSymInt,61,8)';
    pageSym = reshape(symMat,1,[])';
    
    %--- Restore the inverted G2 convolutional-code branch ----------------
    pageSym(2:2:end) = ~pageSym(2:2:end);
    %--- Remove convolutional encoding from page parts --------------------
    decBits = vitdec(pageSym,trellis,tblen,'trunc','hard');
    
    %--- Reconstruct full page --------------------------------------------
    page = decBits(1:238)';
    
    %--- Check the CRC ----------------------------------------------------
    [~,frmError] = step(crcDet,page');
    if (frmError)
        continue
    end
    
    %--- Pull out navigation data word ------------------------------------
    navWordDec = page;
    
    %--- Convert navigation data word to binary ---------------------------
    navWord = dec2bin(navWordDec)';
    
    %--- Decode the message type ------------------------------------------
    wordType = bin2dec(navWord(1:6));
    
    %--- Decode sub-frame based on the message type -----------------------
    switch wordType
        case 1   %--- SVID , Clock correction, SISA, Ionospheric correction,
            % BGD, GST, Signal health and Data validity status
            eph.SVID      = bin2dec(navWord(7:12));                     % 6 bits
            eph.IODnav1   = bin2dec(navWord(13:22));                    % 10 bits
            eph.t_oc      = bin2dec(navWord(23:36))        * 60;        % 14 bits in s
            eph.a_f0      = twosComp2dec(navWord(37:67))   * 2^(-34);   % 31 bit in s
            eph.a_f1      = twosComp2dec(navWord(68:88)) * 2^(-46);     % 21 bit in s/s
            eph.a_f2      = twosComp2dec(navWord(89:94)) * 2^(-59);     % 6 bit in s/(s^2)
            %             eph.SISA_E1_E5a = bin2dec(navWord(95:102));                 % 8 bit           eph.SISA      = navWord(121:128) % content not currently defined
            eph.a_i0      = bin2dec(navWord(103:113))         * 2^(-2);    % 11 bit
            eph.a_i1      = twosComp2dec(navWord(114:124))   * 2^(-8);    % 11 bit
            eph.a_i2      = twosComp2dec(navWord(125:138))   * 2^(-15);   % 14 bit
            eph.iono_SF1  = bin2dec(navWord(139));      % 1 bit
            eph.iono_SF2  = bin2dec(navWord(140));      % 1 bit
            eph.iono_SF3  = bin2dec(navWord(141));      % 1 bit
            eph.iono_SF4  = bin2dec(navWord(142));      % 1 bit
            eph.iono_SF5  = bin2dec(navWord(143));      % 1 bit
            eph.BGD_E1E5a = twosComp2dec(navWord(144:153))   * 2^(-32);  % 10 bit
            eph.E5a_HS    = bin2dec(navWord(154:155));   % 2 bit  E5a signal Health Status (0-OK)
            eph.WN        = bin2dec(navWord(156:167));   % 12 bit
            eph.TOW       = bin2dec(navWord(168:187));  % 20 bit
            eph.E5a_DVS   = bin2dec(navWord(188));      % 1 bit E5a Data validity status (0-valid)
            % Correct TOW to time for first page part
            eph.TOW = eph.TOW - (ii-1)*10;
            if ~isequal(eph.E5a_HS,1)   % && isequal(eph.E5a_DVS,0)
                word_valid(1) = 1;
            end
            
        case 2   %--- Ephemeris (1/3) and GST -------------------------------------
            eph.IODnav2   = bin2dec(navWord(7:16));                           % 10 bit
            eph.M_0       = twosComp2dec(navWord(17:48))   * 2^(-31) * galPi; % 32 bit in rad
            eph.OmegaDot  = twosComp2dec(navWord(49:72))   * 2^(-43) * galPi; % 24 bit in rad
            eph.e         = bin2dec(navWord(73:104))        * 2^(-33);         % 32 bit
            eph.sqrtA     = bin2dec(navWord(105:136))       * 2^(-19);         % 32 bit in m^0.2
            eph.Omega_0   = twosComp2dec(navWord(137:168))   * 2^(-31) * galPi; % 32 bit in rad
            eph.iDot      = twosComp2dec(navWord(169:182)) * 2^(-43) * galPi; % 14 bit in rad
            
            word_valid(2) = 1;
            
        case 3   %--- Ephemeris (2/3) and GST -------------------------------
            eph.IODnav3   = bin2dec(navWord(7:16));                           % 10 bit
            eph.i_0       = twosComp2dec(navWord(17:48))   * 2^(-31) * galPi; % 32 bit in rad
            eph.omega     = twosComp2dec(navWord(49:80))  * 2^(-31) * galPi; % 32 bit in rad
            eph.deltan    = twosComp2dec(navWord(81:96))   * 2^(-43) * galPi; % 16 bit in rad
            eph.CUC       = twosComp2dec(navWord(97:112))   * 2^(-29);         % 16 bit in rad
            eph.CUS       = twosComp2dec(navWord(113:128))   * 2^(-29);         % 16 bit in rad
            eph.CRC       = twosComp2dec(navWord(129:144))  * 2^(-5);          % 16 bit in m
            eph.CRS       = twosComp2dec(navWord(145:160)) * 2^(-5);          % 16 bit in m
            eph.t_oe      = bin2dec(navWord(161:174))        * 60;              % 14 bit in s
            
            word_valid(3) = 1;
            
        case 4   %--- Ephemeris (3/3), GST-UTC conversion, GST-GPS conversion and TO W
            eph.IODnav4   = bin2dec(navWord(7:16));                    % 10 bit
            eph.CIC       = twosComp2dec(navWord(17:32))   * 2^(-29);  % 16 bit in rad
            eph.CIS       = twosComp2dec(navWord(33:48))   * 2^(-29);  % 16 bit in rad
            eph.A0        = twosComp2dec(navWord(49:80))    * 2^(-30);  % 32 bit
            eph.A1        = twosComp2dec(navWord(81:104))   * 2^(-50);  % 24 bit
            eph.delt_LS   = twosComp2dec(navWord(105:112));              % 8 bit
            eph.t_ot      = bin2dec(navWord(113:120))        * 3600;     % 8 bit
            eph.WN_ot     = bin2dec(navWord(121:128));                   % 8 bit
            eph.WN_LSF    = bin2dec(navWord(129:136));                   % 8 bit
            eph.DN        = bin2dec(navWord(137:139));                   % 3 bit
            eph.delt_LSF  = twosComp2dec(navWord(140:147));             % 8 bit
            eph.t_og      = bin2dec(navWord(148:155))      * 3600;     % 8 bit
            eph.A0_G      = twosComp2dec(navWord(156:171))  * 2^(-35);  % 16 bit
            eph.A1_G      = twosComp2dec(navWord(172:183)) * 2^(-51);  % 12 bit
            eph.WN_og     = bin2dec(navWord(184:189));                 % 6 bit
            
            word_valid(4) = 1;
            
            %case 5   %--- Almanac -------------------------------------
            
            %case 6   %--- Almanac -------------------------------------
            
            
    end % switch word type
    if sum(word_valid) == 4
        % Check if the Issue of Data (IOD ) values are the same
        IODnav = [eph.IODnav1 eph.IODnav2 eph.IODnav3 eph.IODnav4];
        if (length(unique(IODnav)) == 1)
            eph.flag = 1;
        else
            disp('    The IOD values of each pages are different!');
            eph.flag = 0;
        end
        
        break
    end
    
end % all pages
