function [eph, firstSubFrame, TOW] = BCNAV2decoding(I_P_InputBits, Q_P_InputBits, PRNin)
% BCNAV2decoding decodes BeiDou B2a B-CNAV2 navigation message.
%
% This version follows the structure of BCNAV3decoding:
%   decodedNavBits = B2aLDPCDecoder(navBits_soft_I, navBits_soft_Q);
%
% Inputs:
%   I_P_InputBits   - Prompt I-branch tracking output, soft bits
%   Q_P_InputBits   - Prompt Q-branch tracking output, soft/noise reference
%   PRNin           - Current PRN number, optional
%
% Outputs:
%   eph             - SV ephemeris
%   firstSubFrame   - Starting position of the first valid B-CNAV2 frame
%   TOW             - Time of week / seconds of week

% -------------------------------------------------------------------------
% Initialize ephemeris structure
% -------------------------------------------------------------------------
eph = eph_structure_init();

firstSubFrame = inf;
TOW = inf;

%% Bit and frame synchronization ==========================================
searchStartOffset = 0;

% B2a secondary code, length = 5
secondCode = [1 1 1 -1 1];

% B-CNAV2 preamble:
% [1 1 1 0 0 0 1 0 0 1 0 0 1 1 0 1 1 1 1 0 1 0 0 0]
% Antipodal form:
preamble_bits = [-1 -1 -1 1 1 1 -1 1 1 -1 1 1 ...
                 -1 -1 1 -1 -1 -1 -1 1 -1 1 1 1];

% Upsampled preamble with secondary code
preamble_ms = kron(preamble_bits, secondCode);

% -------------------------------------------------------------------------
% Extract raw soft information
% -------------------------------------------------------------------------
raw_soft_bits_I = I_P_InputBits(1 + searchStartOffset : end);

if nargin >= 2 && ~isempty(Q_P_InputBits)
    raw_soft_bits_Q = Q_P_InputBits(1 + searchStartOffset : end);
else
    raw_soft_bits_Q = zeros(size(raw_soft_bits_I));
end

% Make sure I/Q branches have the same length
minLen = min(length(raw_soft_bits_I), length(raw_soft_bits_Q));
raw_soft_bits_I = raw_soft_bits_I(1:minLen);
raw_soft_bits_Q = raw_soft_bits_Q(1:minLen);

% Generate hard decision bits only for preamble synchronization
bits_hard = raw_soft_bits_I;
bits_hard(bits_hard > 0)  =  1;
bits_hard(bits_hard <= 0) = -1;

% -------------------------------------------------------------------------
% Correlate tracking output with the preamble
% -------------------------------------------------------------------------
tlmXcorrResult = xcorr(bits_hard, preamble_ms);

xcorrLength = (length(tlmXcorrResult) + 1) / 2;

% Ideal correlation peak is 24*5 = 120
index = find(abs(tlmXcorrResult(xcorrLength : xcorrLength * 2 - 1)) > 115)';

% CRC-24Q detector
crcDet = comm.CRCDetector([24 23 18 17 14 11 10 7 6 5 4 3 1 0]);
% LDPC/EMS decoder path
addpath ('include/ldpcDecoder/') 
%% B-CNAV2 decoding ========================================================
for i = 1:length(index)

    % One B-CNAV2 frame:
    % 600 navigation bits before LDPC decoding
    % Each bit is spread by 5-chip secondary code
    % Total length = 600*5 = 3000 samples
    if ((length(bits_hard) - index(i) + 1) >= 600 * 5)

        % -----------------------------------------------------------------
        % Extract one complete 3000-sample B-CNAV2 frame
        % -----------------------------------------------------------------
        frame_soft_I = raw_soft_bits_I(index(i) : index(i) + 600*5 - 1);
        frame_soft_Q = raw_soft_bits_Q(index(i) : index(i) + 600*5 - 1);
        frame_hard   = bits_hard(index(i)      : index(i) + 600*5 - 1);

        % -----------------------------------------------------------------
        % Despread with secondary code
        % -----------------------------------------------------------------
        I_group_soft = reshape(frame_soft_I, 5, [])';
        Q_group_soft = reshape(frame_soft_Q, 5, [])';
        hard_group   = reshape(frame_hard,   5, [])';

        navBits_soft_I = sum(I_group_soft .* repmat(secondCode, 600, 1), 2)';
        navBits_soft_Q = sum(Q_group_soft .* repmat(secondCode, 600, 1), 2)';

        navBits_hard = sign(sum(hard_group .* repmat(secondCode, 600, 1), 2))';
        navBits_hard(navBits_hard == 0) = -1;

        % -----------------------------------------------------------------
        % Polarity correction
        % -----------------------------------------------------------------
        if ~isequal(navBits_hard(1:24), preamble_bits)
            navBits_hard   = -navBits_hard;
            navBits_soft_I = -navBits_soft_I;
            navBits_soft_Q = -navBits_soft_Q;
        end

        % Check preamble after polarity correction
        if ~isequal(navBits_hard(1:24), preamble_bits)
            continue;
        end

        % -----------------------------------------------------------------
        % B2a B-CNAV2 LDPC / EMS decoder
        % -----------------------------------------------------------------
        % Decode the 576-bit LDPC payload into the 288-bit B-CNAV2 data block.
        decodedNavBits = B2aLDPCDecoder(navBits_soft_I, navBits_soft_Q);

        % -----------------------------------------------------------------
        % CRC-24Q check
        % -----------------------------------------------------------------
        NavDataLogic = (decodedNavBits == 1);

        [~, frmError] = step(crcDet, NavDataLogic);

        % -----------------------------------------------------------------
        % Decode ephemeris message by message
        % -----------------------------------------------------------------
        if ~frmError

            navBitsChar = char(NavDataLogic + '0');

            eph_tmp = ephemeris(navBitsChar', eph);

            % Optional PRN check
            if nargin >= 3 && ~isempty(PRNin)
                if isfield(eph_tmp, 'PRN')
                    if ~isempty(eph_tmp.PRN) && ~isnan(eph_tmp.PRN)
                        if eph_tmp.PRN ~= PRNin
                            continue;
                        end
                    end
                end
            end

            eph = eph_tmp;

            if isinf(firstSubFrame)
                firstSubFrame = index(i) + searchStartOffset;

                if isfield(eph, 'SOW')
                    TOW = eph.SOW;
                else
                    TOW = inf;
                end
            end
        end
    end
end

end
