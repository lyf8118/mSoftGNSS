function [eph, SOW] = ephemeris(bits, PRN)
%Function decodes ephemerides and SOW from the given bit stream. D1 decoding
%uses 1500-2100 bits (five to seven subframes), while D2 decoding uses 15000
%bits (fifty pages). The first element in the array must be the first bit of
%a subframe or page. The first supplied subframe/page ID is not important.
%
%Function does not check parity!
%
%[eph, SOW] = ephemeris(bits, PRN)
%
%   Inputs:
%       bits        - bits of the navigation messages.
%                   Type is character array and it must contain only
%                   characters '0' or '1'.
%       PRN         - PRN number to separate the decoding process.
%                   GEO & MEO/IGSO have different message structure.
%
%   Outputs:
%       SOW         - Second Of Week (SOW) of the first sub-frame in the bit
%                   stream (in seconds)
%       eph         - SV ephemeris
%--------------------------------------------------------------------------
%                         CU Multi-GNSS SDR  
% (C) Written by Daehee Won, Yafeng Li, Nagaraj C. Shivaramaiah and Dennis M. Akos
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

% Modified by Daehee Won.
% Final revision: Sep. 20, 2012.


%--- Decode Navigation mesaage ----------------------
% The task is to select the necessary bits and convert them to decimal
% numbers. For more details on sub-frame contents please refer to
% Beidou ICD (Version 1.0, December 2012).


eph = eph_structure_init();
SOW = inf;
%% Check if there is enough data ==========================================
if length(bits) < 1500
    error('The parameter BITS must contain at least 1500 bits!');
end

%% Check if the parameters are strings ====================================

% Pi used in the GPS coordinate system
BeidouPi = 3.1415926535898;
% Preamble for Beidou
preamble = [1 1 1 0 0 0 1 0 0 1 0];
Ipreamble = [0 0 0 1 1 1 0 1 1 0 1];

% 'bits' should be row vector for 'bin2dec' function.
[a, b] = size(bits);
if a > b
    bits = bits';
end

if isequal(bits(1:11),preamble) || isequal(bits(1:11),Ipreamble)
    %--- Correct polarity of the data bits in all 10 words ----------------
    if isequal(bits(1:11),Ipreamble)
        bits = bits*-1 +1;
    end
else
    disp(['Preamble does NOT match in for PRN ',num2str(PRN),'!']);
    return;
end

%% Ephemeris for GEO ======================================================
if ((1 <= PRN) && (PRN <=5 )) || ((59 <= PRN) && (PRN <= 63))
    
    % GEO has 10 words in subframe 1. Each word has 150 navigation bits
    % i.e., 300 raw bits (300 ms)
    
    % BCH decoding & interleaving should be apllied!!!
    
    % Keep every ten-page D2 cycle separate. The cycle index is based on
    % cycleKey = i-Pnum1, with an offset that makes it a MATLAB index.
    cycleEph   = cell(59,1);
    cycleParts = cell(59,1);
    cycleValid = false(59,10);
    
    for i = 1:50
        
        %--- "Cut" one sub-frame's 300 bits -------------------------------
        subframe = bits(300*(i-1)+1 : 300*i);
        
        % Deinterleaving---------------------------------------------------
        % The first word
        DeSubframe = subframe(1:30);
        % Other 4 words
        for k = 1:4 % 150 LSBs are reserved
            DeSubframe = [DeSubframe, ...
                subframe(k*30+1:2:k*30+22), ...
                subframe(k*30+2:2:k*30+22), ...
                subframe(k*30+23:2:k*30+30), ...
                subframe(k*30+24:2:k*30+30)]; %#ok<AGROW>
        end
        subframe = DeSubframe;
        
        %--- Decode the sub-frame id ------------------------------------------
        % Decode the received signal in code using an [15,11] BCH decoder
        [decoded,cnumerr] = bchdec( gf(subframe(16:30),1),15,11 ); 
        if cnumerr ~= -1
            subframe(16:26) = double(decoded.x);
            IDbin = num2str(subframe(16:18));
            subframeID = bin2dec(IDbin);
        else
            disp(['ID of subframe ',num2str(i),'for PRN# ',...
                num2str(PRN),'can not be decoded!']);
            continue
        end
        
        % For more details on sub-frame contents please refer to Beidou IS.
        % Only Subframe 1 includes the navigation message.
        if subframeID == 1
            % BCH decode
            bchValid = true;
            for ind = 1:4
                codeword = [subframe(30*ind+1:30*ind+11), subframe(30*ind+23:30*ind+26);...
                    subframe(30*ind+12:30*ind+22), subframe(30*ind+27:30*ind+30)];
                [decoded,cnumerr] = bchdec( gf(codeword,1),15,11 );
                decodeBCH = double(decoded.x);
                if all(cnumerr ~= -1)
                    subframe(30*ind+1:30*ind+11) = decodeBCH(1,1:11);
                    subframe(30*ind+12:30*ind+22) = decodeBCH(2,1:11);
                else
                    disp(['BCH decoding for PRN#',num2str(PRN),'fails!']);
                    bchValid = false;
                    break
                end
            end
            if ~bchValid
                continue
            end
            
            subframe = num2str(subframe')';
            Pnum1 = bin2dec(subframe(43:46));   % range: 1-10
            if Pnum1 < 1 || Pnum1 > 10
                continue
            end
            
            if SOW == inf
                SOW = bin2dec(subframe([19:26, 31:42])) - 0.6*(i-1);
                eph.SOW = SOW;
            end

            cycleIndex = i-Pnum1+10;
            if isempty(cycleEph{cycleIndex})
                candidate = eph_structure_init();
                parts = struct('a1_msb',[],'a1_lsb',[], ...
                    'C_uc_msb',[],'C_uc_lsb',[], ...
                    'e_msb',[],'e_lsb',[], ...
                    'C_ic_msb',[],'C_ic_lsb',[], ...
                    'i_0_msb',[],'i_0_lsb',[], ...
                    'omegaDot_msb',[],'omegaDot_lsb',[], ...
                    'omega_msb',[],'omega_lsb',[]);
            else
                candidate = cycleEph{cycleIndex};
                parts = cycleParts{cycleIndex};
            end
            
            switch Pnum1
                case 1
                    candidate.SatH1 = bin2dec( subframe(47) );  % Autonomous Satellite Health flag: 0, 1
                    candidate.IODC  = bin2dec( subframe(48:52) ); % Issue of Data, Clock
                    candidate.URAI  = bin2dec( subframe(61:64) ); % User Range Accuracy Index: 0~15
                    candidate.WN    = bin2dec( subframe(65:77) ); % Week Number: 0~8191
                     
                    candidate.t_oc = bin2dec( subframe([78:82, 91:102]) ) *2^(3); % Clock Correction Parameters
                    % B3I is the BDS clock-reference signal, so no
                    % inter-signal group-delay correction is applied.
                    candidate.T_GD_1 = 0;
                     
                case 2  % Ionospheric Delay Model Parameters (alpha, beta)
                    candidate.alpha0 = twosComp2dec( subframe([47:52, 61:62]) ) *2^(-30);  % [s]
                    candidate.alpha1 = twosComp2dec( subframe(63:70) ) *2^(-27);           % [s/pi]
                    candidate.alpha2 = twosComp2dec( subframe(71:78) ) *2^(-24);           % [s/pi^2]
                    candidate.alpha3 = twosComp2dec( subframe([79:82, 91:94]) ) *2^(-24);  % [s/pi^3]
                     
                    candidate.beta0 = twosComp2dec( subframe(95:102) ) *2^(11);           % [s]
                    candidate.beta1 = twosComp2dec( subframe(103:110) ) *2^(14);          % [s/pi]
                    candidate.beta2 = twosComp2dec( subframe([111:112, 121:126]) ) *2^(16); % [s/pi^2]
                    candidate.beta3 = twosComp2dec( subframe(127:134) ) *2^(16);          % [s/pi^3]
                     
                case 3  % Clock Correction Parameters
                    candidate.a0 = twosComp2dec( subframe([101:112, 121:132]) ) *2^(-33); % [s]
                    parts.a1_msb = subframe(133:136);
                     
                case 4
                    % Clock Correction Parameters (Cont.)
                    parts.a1_lsb = subframe([47:52, 61:72]);
                    candidate.a2 = twosComp2dec( subframe([73:82, 91]) ) *2^(-66); % [s/s^2]
                     
                    % Issue of Data, Ephemeris (IODE)
                    candidate.IODE = bin2dec( subframe(92:96) );
                     
                    % Ephemeris Parameters
                    candidate.deltan = twosComp2dec( subframe(97:112) ) * 2^(-43) * BeidouPi; % [pi/s]
                    parts.C_uc_msb = subframe(121:134);
                     
                case 5  % Ephemeris Parameters (Cont.)
                    parts.C_uc_lsb = subframe(47:50);
                    candidate.M_0 = twosComp2dec( subframe([51:52, 61:82, 91:98]) ) * 2^(-31) * BeidouPi; % [pi]
                    candidate.C_us = twosComp2dec( subframe([99:112, 121:124]) ) * 2^(-31);
                    parts.e_msb = subframe(125:134);
                     
                case 6  % Ephemeris Parameters (Cont.)
                    parts.e_lsb = subframe([47:52, 61:76]);
                    candidate.sqrtA = bin2dec( subframe([77:82, 91:112, 121:124]) ) * 2^(-19);
                    parts.C_ic_msb = subframe(125:134);
                     
                case 7  % Ephemeris Parameters (Cont.)
                    parts.C_ic_lsb = subframe([47:52, 61:62]);
                    candidate.C_is = twosComp2dec( subframe(63:80) ) * 2^(-31);
                    candidate.t_oe = bin2dec( subframe([81:82, 91:105]) ) * 2^3;
                    parts.i_0_msb = subframe([106:112, 121:134]);
                     
                case 8  % Ephemeris Parameters (Cont.)
                    parts.i_0_lsb = subframe([47:52, 61:65]);
                    candidate.C_rc = twosComp2dec( subframe([66:82, 91]) ) * 2^(-6);
                    candidate.C_rs = twosComp2dec( subframe(92:109) ) * 2^(-6);
                    parts.omegaDot_msb = subframe([110:112, 121:136]);
                     
                case 9  % Ephemeris Parameters (Cont.)
                    parts.omegaDot_lsb = subframe(47:51);
                    candidate.omega_0 = twosComp2dec( subframe([52, 61:82, 91:99]) ) * 2^(-31) * BeidouPi;
                    parts.omega_msb = subframe([100:112, 121:134]);
                case 10  % Ephemeris Parameters (Cont.)
                    parts.omega_lsb = subframe(47:51);
                    candidate.iDot = twosComp2dec( subframe([52, 61:73]) ) * 2^(-43) * BeidouPi;
            end
            cycleValid(cycleIndex,Pnum1) = true;
        
        %% MSB & LSB combination
        if length([parts.a1_msb parts.a1_lsb]) == 22
            candidate.a1 = twosComp2dec( [parts.a1_msb, parts.a1_lsb] ) *2^(-50);
        end
        
        if length([parts.C_uc_msb, parts.C_uc_lsb]) == 18
            candidate.C_uc = twosComp2dec( [parts.C_uc_msb, parts.C_uc_lsb] ) * 2^(-31);
        end
        
        if length([parts.e_msb, parts.e_lsb]) == 32
            candidate.e = bin2dec( [parts.e_msb, parts.e_lsb] ) * 2^(-33);
        end
        
        if length([parts.C_ic_msb, parts.C_ic_lsb]) == 18
            candidate.C_ic = twosComp2dec( [parts.C_ic_msb, parts.C_ic_lsb] ) * 2^(-31);
        end
        
        if length([parts.i_0_msb, parts.i_0_lsb]) == 32
            candidate.i_0 = twosComp2dec( [parts.i_0_msb, parts.i_0_lsb] ) * 2^(-31) * BeidouPi;
        end
        
        if length([parts.omegaDot_msb, parts.omegaDot_lsb]) == 24
            candidate.omegaDot = twosComp2dec( [parts.omegaDot_msb, parts.omegaDot_lsb] ) * 2^(-43) * BeidouPi;
        end
        
        if length([parts.omega_msb, parts.omega_lsb]) == 32
            candidate.omega = twosComp2dec( [parts.omega_msb, parts.omega_lsb] ) * 2^(-31) * BeidouPi;
        end
        cycleEph{cycleIndex} = candidate;
        cycleParts{cycleIndex} = parts;
        end % if subframeID == 1
    end % for i = 1:50
    
    completeCycle = find(all(cycleValid,2),1);
    if ~isempty(completeCycle)
        eph = cycleEph{completeCycle};
        eph.SOW = SOW;
        eph.flag = 1;
    end
    
    % Compute the second of week (SOW) of the first sub-frames in the array ====
    % Also correct the SOW. The transmitted SOW is actual SOW of the next
    % subframe and we need the SOW of the first subframe in this data block
    % (the variable subframe at this point contains bits of the last subframe).
    % D2 subframe is 3 seconds long.
    %     subframe = num2str(subframe')';
    %     SOW = bin2dec(subframe([19:26, 31:42]))-27;
    
    %% Ephemeris for MEO/IGSO  ================================================
elseif (6 <= PRN) && (PRN <= 58)
    % Keep subframes 1, 2 and 3 within one five-subframe cycle. The input
    % can contain five to seven subframes, depending on available data.
    cycleEph   = cell(9,1);
    cycleToe   = cell(9,2);
    cycleValid = false(9,3);
    frameCount = floor(length(bits)/300);
    
    % BCH decoding & interleaving should be applied!!!
    for i = 1:frameCount
        
        %--- "Cut" one sub-frame's bits ---------------------------------------
        subframe = bits(300*(i-1)+1 : 300*i);
        
        % Deinterleaving --------------------------------------------------
        DeSubframe = subframe(1:30);
        for k = 1:9
            DeSubframe = [DeSubframe, ...
                subframe(k*30+1:2:k*30+22), ...
                subframe(k*30+2:2:k*30+22), ...
                subframe(k*30+23:2:k*30+30), ...
                subframe(k*30+24:2:k*30+30)]; %#ok<AGROW>
        end
        subframe = DeSubframe;
        
        %--- Decode the sub-frame id ------------------------------------------
        % For more details on sub-frame contents please refer to GPS IS.
        %--- Decode the sub-frame id ------------------------------------------
        [decoded,cnumerr] = bchdec( gf(subframe(16:30),1),15,11 );
        if cnumerr ~= -1
            subframe(16:26) = double(decoded.x);
            IDbin = num2str(subframe(16:18));
            subframeID = bin2dec(IDbin);
        else
            disp(['ID of subframe ',num2str(i),'for PRN# ',...
                num2str(PRN),'can not be decoded!']);
            continue
        end
        
        % do BCH decoding
        if subframeID == 1 || subframeID == 2 || subframeID == 3
            bchValid = true;
            for ind = 1:4
                codeword = [subframe(30*ind+1:30*ind+11), subframe(30*ind+23:30*ind+26);...
                    subframe(30*ind+12:30*ind+22), subframe(30*ind+27:30*ind+30)];
                [decoded,cnumerr] = bchdec( gf(codeword,1),15,11 );
                decodeBCH = double(decoded.x);
                if all(cnumerr ~= -1)
                    subframe(30*ind+1:30*ind+11) = decodeBCH(1,1:11);
                    subframe(30*ind+12:30*ind+22) = decodeBCH(2,1:11);
                else
                    disp(['BCH decoding for PRN#',num2str(PRN),'fails!']);
                    bchValid = false;
                    break
                end
            end
            if ~bchValid
                continue
            end
            subframe = num2str(subframe')';
            
            if SOW==inf
                SOW = bin2dec(subframe([19:26, 31:42])) - (i-1)*6;
                eph.SOW = SOW;
            end
        else
            continue
        end

        % frameKey = (i-1)-(subframeID-1) = i-subframeID. The offset of
        % three maps every possible key for seven input subframes to 1:9.
        cycleIndex = i-subframeID+3;
        if isempty(cycleEph{cycleIndex})
            candidate = eph_structure_init();
        else
            candidate = cycleEph{cycleIndex};
        end
        
        switch subframeID
            case 1  %--- It is subframe 1 -------------------------------------
                
                candidate.SatH1 = bin2dec( subframe(43) );    % Autonomous Satellite Health flag: 0, 1
                candidate.IODC  = bin2dec( subframe(44:48) ); % Issue of Data, Clock
                candidate.URAI  = bin2dec( subframe(49:52) ); % User Range Accuracy Index: 0~15
                candidate.WN    = bin2dec( subframe(61:73) ); % Week Number: 0~8191
                
                candidate.t_oc = bin2dec( subframe([74:82, 91:98]) ) * 2^(3); % Clock Correction Parameters
                % B3I is the BDS clock-reference signal, so no inter-signal
                % group-delay correction is applied to B3I observations.
                candidate.T_GD_1 = 0;
                
                % Ionospheric Delay Model Parameters (alpha, beta)
                candidate.alpha0 = twosComp2dec( subframe(127:134) ) *2^(-30); % [s]
                candidate.alpha1 = twosComp2dec( subframe(135:142) ) *2^(-27); % [s/pi]
                candidate.alpha2 = twosComp2dec( subframe(151:158) ) *2^(-24); % [s/pi^2]
                candidate.alpha3 = twosComp2dec( subframe(159:166) ) *2^(-24); % [s/pi^3]
                
                candidate.beta0 = twosComp2dec( subframe([167:172,181:182]) ) *2^(11); % [s]
                candidate.beta1 = twosComp2dec( subframe(183:190) ) *2^(14);           % [s/pi]
                candidate.beta2 = twosComp2dec( subframe(191:198) ) *2^(16);           % [s/pi^2]
                candidate.beta3 = twosComp2dec( subframe([199:202, 211:214]) ) *2^(16);% [s/pi^3]
                
                % Clock Correction Parameters
                candidate.a2 = twosComp2dec( subframe(215:225) ) *2^(-66);            % [s/s^2]
                candidate.a0 = twosComp2dec( subframe([226:232, 241:257]) ) *2^(-33); % [s]
                candidate.a1 = twosComp2dec( subframe([258:262, 271:287]) ) *2^(-50); % [s/s]
                
                % Issue of Data, Ephemeris (IODE)
                candidate.IODE = bin2dec( subframe(288:292) );
                cycleValid(cycleIndex,1) = true;
                
            case 2  %--- It is subframe 2 -------------------------------------
                % Ephemeris Parameters
                candidate.deltan = twosComp2dec( subframe([43:52, 61:66])) * 2^(-43) * BeidouPi;
                candidate.C_uc = twosComp2dec( subframe([67:82, 91:92])) * 2^(-31);
                candidate.M_0 = twosComp2dec( subframe([93:112, 121:132])) * 2^(-31) * BeidouPi;
                candidate.e = bin2dec( subframe([133:142, 151:172])) * 2^(-33);
                
                candidate.C_us = twosComp2dec( subframe(181:198)) * 2^(-31);
                candidate.C_rc = twosComp2dec( subframe([199:202, 211:224])) * 2^(-6);
                candidate.C_rs = twosComp2dec( subframe([225:232, 241:250])) * 2^(-6);
                
                candidate.sqrtA = bin2dec( subframe([251:262, 271:290])) * 2^(-19);
                cycleToe{cycleIndex,1} = subframe(291:292);
                cycleValid(cycleIndex,2) = true;
                
            case 3  %--- It is subframe 3 -------------------------------------
                cycleToe{cycleIndex,2} = subframe([43:52, 61:65]);
                candidate.i_0 = twosComp2dec( subframe([66:82, 91:105])) * 2^(-31) * BeidouPi;
                candidate.C_ic = twosComp2dec( subframe([106:112, 121:131])) * 2^(-31);
                candidate.omegaDot = twosComp2dec( subframe([132:142, 151:163])) * 2^(-43) * BeidouPi;
                candidate.C_is = twosComp2dec( subframe([164:172, 181:189])) * 2^(-31);
                candidate.iDot = twosComp2dec( subframe([190:202, 211])) * 2^(-43) * BeidouPi;
                candidate.omega_0 = twosComp2dec( subframe([212:232, 241:251])) * 2^(-31) * BeidouPi;
                candidate.omega = twosComp2dec( subframe([252:262, 271:291])) * 2^(-31) * BeidouPi;
                cycleValid(cycleIndex,3) = true;
                
        end % switch subframeID ...
        cycleEph{cycleIndex} = candidate;
    end % for all available subframes ...
    
    completeCycle = find(all(cycleValid,2),1);
    if ~isempty(completeCycle)
        eph = cycleEph{completeCycle};
        toeBits = [cycleToe{completeCycle,1},cycleToe{completeCycle,2}];
        if length(toeBits) == 17
            eph.t_oe = bin2dec(toeBits) * 2^3;
        else
            eph.t_oe = nan;
        end
        eph.SOW = SOW;
        eph.flag = 1;
    end
    % Compute the second of week (SOW) of the first sub-frames in the array ====
    % Also correct the SOW. The transmitted SOW is actual SOW of the next
    % subframe and we need the SOW of the first subframe in this data block
    % (the variable subframe at this point contains bits of the last subframe).
    % D1 subframe is 6 seconds long.
else
    disp('PRN is NOT in range between 1-63 ');
    
end
