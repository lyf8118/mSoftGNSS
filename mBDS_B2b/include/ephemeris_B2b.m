function [eph] = ephemeris_B2b(navBitsBin,eph)
%Function decodes ephemerides and TOW from the given bit stream. The stream
%(array) in the parameter BITS must contain 486 bits. The first element in
%the array must be the first bit of a subframe.
%
%
%[eph] = ephemeris(navBitsBin,eph)
%
%   Inputs:
%       navBitsBin  - bits of the navigation messages.Type is character array
%                   and it must contain only characters '0' or '1'.
%       eph         - The ephemeris for each PRN is decoded message by message.
%                   To prevent loss of previous decoded messages, the eph structure
%                   must be passed onto this function.
%   Outputs:
%       eph         - SV ephemeris

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

%CVS record:
%$Id: ephemeris.m,v 1.1.2.7 2006/08/14 11:38:22 dpl Exp $

%% Preparation for data message decoding ============================
if length(navBitsBin) < 486
    error('The parameter BITS must contain 1000 bits!');
end

% Check if the parameters are strings
if ~ischar(navBitsBin)
    error('The parameter BITS must be a character array!');
end

% 'bits' should be row vector for 'bin2dec' function.
[a, b] = size(navBitsBin);
if a > b
    navBitsBin = navBitsBin';
end

% Pi used in the GPS coordinate system
gpsPi = 3.1415926535898;

% Decode the message id
MesType = bin2dec(navBitsBin(1:6));

%%  Decode messages based on the message id =========================
% The task is to select the necessary bits and convert them to decimal
% numbers. For more details on message contents please refer to BDS-3
% ICD (BDS-SIS-ICD-B2b-1.0).
switch MesType
    % Message type 10  provides users
    % the requisite data to calculate SV position.
    case 10  %--- It is Message Type 10 -----------------------------------
        % It contains first part of ephemeris parameters
        eph.idValid(1) = 10;
        % SOW
        if isempty(eph.SOW)
            eph.SOW  = bin2dec(navBitsBin(7:26));%%比例因子为1
        end
        
        %星历I
        % Ephemeris data reference time of week
        eph.t_oe        = bin2dec(navBitsBin(31:41)) * 300;
        % Satellite type
        SatType     = bin2dec(navBitsBin(42:43));
        if (SatType == 1)
            eph.SatType = 'GEO';
        elseif (SatType == 2)
            eph.SatType = 'IGSO';
        elseif (SatType == 3)
            eph.SatType = 'MEO';
        end
        % Semi-major axis difference at reference time
        eph.deltaA      = twosComp2dec(navBitsBin(44:69)) * 2^(-9) ;
        % Change rate in semi-major axis
        eph.ADot        = twosComp2dec(navBitsBin(70:94)) * 2^(-21);
        % Mean Motion difference from computed value at reference time
        eph.delta_n_0   = twosComp2dec(navBitsBin(95:111)) * 2^(-44)* gpsPi;
        % Rate of mean motion difference from computed value
        eph.delta_n_0Dot= twosComp2dec(navBitsBin(112:134)) * 2^(-57)* gpsPi;
        % Mean anomaly at reference time
        eph.M_0         = twosComp2dec(navBitsBin(135:167)) * 2^(-32) * gpsPi;
        % Eccentricity
        eph.e           = bin2dec(navBitsBin(168:200))* 2^(-34);
        % Argument of perigee
        eph.omega       = twosComp2dec(navBitsBin(201:233))* 2^(-32) * gpsPi;
        
        %星历II
        % Longitude of Ascending Node of Orbit Plane at Weekly Epoch
        eph.omega_0     = twosComp2dec(navBitsBin(234:266))* 2^(-32) * gpsPi;
        % Inclination angle at reference time
        eph.i_0         = twosComp2dec(navBitsBin(267:299))* 2^(-32) * gpsPi;
        % Rate of right ascension difference
        eph.omegaDot  = twosComp2dec(navBitsBin(300:318)) * 2^(-44) * gpsPi;
        % Rate of inclination angle
        eph.i_0Dot      = twosComp2dec(navBitsBin(319:333)) * 2^(-44) * gpsPi;
        % Amplitude of the sine harmonic correction term to the angle of inclination
        eph.C_is        = twosComp2dec(navBitsBin(334:349)) * 2^(-30);
        % Amplitude of the cosine harmonic correction term to the angle of inclination
        eph.C_ic        = twosComp2dec(navBitsBin(350:365)) * 2^(-30);
        % Amplitude of the sine correction term to the orbit radius
        eph.C_rs        = twosComp2dec(navBitsBin(366:389)) * 2^(-8);
        % Amplitude of the cosine correction term to the orbit radius
        eph.C_rc        = twosComp2dec(navBitsBin(390:413)) * 2^(-8);
        % Amplitude of the sine harmonic correction term to the argument of latitude
        eph.C_us        = twosComp2dec(navBitsBin(414:434)) * 2^(-30);
        % Amplitude of the cosine harmonic correction term to the argument of latitude
        eph.C_uc        = twosComp2dec(navBitsBin(435:455)) * 2^(-30);
        
        % DIF
        eph.DIFI  = bin2dec(navBitsBin(456));
        % SIF
        eph.SIFI  = bin2dec(navBitsBin(457));
        % AIF
        eph.AIFI  = bin2dec(navBitsBin(458));    
        
    case 30  %--- It is Message Type 30 -----------------------------------
        % It contains second part of ephemeris parameters
        eph.idValid(2)  = 30;
        % SOW
        if isempty(eph.SOW)
            eph.SOW  = bin2dec(navBitsBin(7:26));
        end
        
        % Week No.
        eph.WN  = bin2dec(navBitsBin(27:39));    
        
        % Clock data reference time of Week
        eph.t_oc        = bin2dec(navBitsBin(44:54)) * 300;
        % SV Clock Bias Correction Coefficient
        eph.a_0        = twosComp2dec(navBitsBin(55:79)) * 2^(-34);
        % SV Clock Drift Correction Coefficient
        eph.a_1        = twosComp2dec(navBitsBin(80:101)) * 2^(-50);
        % SV Clock Drift Rate Correction Coefficient
        eph.a_2        = twosComp2dec(navBitsBin(102:112)) * 2^(-66);
        
        % Group delay differential of the B2bI component
          eph.T_GDB2bI=twosComp2dec(navBitsBin(113:124)) * 2^(-34);
          
        % The ionospheric parameters
        eph.alpha1      = bin2dec(navBitsBin(125:134)) * 2^(-3);
        eph.alpha2      = twosComp2dec(navBitsBin(135:142)) * 2^(-3);
        eph.alpha3      = bin2dec(navBitsBin(143:150)) * 2^(-3);
        eph.alpha4      = bin2dec(navBitsBin(151:158)) * 2^(-3);
        eph.alpha5       = bin2dec(navBitsBin(159:166)) * 2^(-3);
        eph.alpha6       = twosComp2dec(navBitsBin(167:174)) * 2^(-3);
        eph.alpha7       = twosComp2dec(navBitsBin(175:182)) * 2^(-3);
        eph.alpha8       = twosComp2dec(navBitsBin(183:190)) * 2^(-3);
        eph.alpha9       = twosComp2dec(navBitsBin(191:198)) * 2^(-3);
        
        
        % BDT-UTC ------------------------------------------
        eph.A_0UTC        = twosComp2dec(navBitsBin(199:214)) * 2^(-35);
        eph.A_1UTC        = twosComp2dec(navBitsBin(215:227)) * 2^(-51);
        eph.A_2UTC        = twosComp2dec(navBitsBin(228:234)) * 2^(-68);
        eph.delta_t_LS    = twosComp2dec(navBitsBin(235:242));
        eph.t_ot          = bin2dec(navBitsBin(243:258)) * 2^(4);
        eph.WN_ot         = bin2dec(navBitsBin(259:271));
        eph.WN_LSF        = bin2dec(navBitsBin(272:284));
        eph.DN            = bin2dec(navBitsBin(285:287));
        eph.delta_t_LSF   = twosComp2dec(navBitsBin(288:295));
      
        % HS
        eph.HS  = bin2dec(navBitsBin(461:462));
         
    case 40 %--- It is Message Type 40 ------------------------------------
        eph.idValid(3) = 40;       
        % SOW        
        if isempty(eph.SOW)
            eph.SOW  = bin2dec(navBitsBin(7:26));
        end       
        
        % BGTO --------------------------------------
        % GNSS ID
        eph.GNSS_ID   = bin2dec(navBitsBin(27:29));
        eph.WN_0BGTO   = bin2dec(navBitsBin(30:42));
        eph.t_0BGTO   = bin2dec(navBitsBin(43:58))* 2^(4);
        eph.A_0BGTO   = twosComp2dec(navBitsBin(59:74))* 2^(-35);
        eph.A_1BGTO   = twosComp2dec(navBitsBin(75:87))* 2^(-51);
        eph.A_2BGTO   = twosComp2dec(navBitsBin(88:94))* 2^(-68);
        % Other terms not decoded at the moment...
    otherwise
        eph.idValid(4) = MesType;
        eph.SOW  = bin2dec(navBitsBin(7:26));%%比例因子为1
end
end % switch MesType ...
