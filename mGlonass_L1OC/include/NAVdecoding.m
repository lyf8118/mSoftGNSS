function [eph, firstSubFrame, TOD] = NAVdecoding(I_P)
%NAVdecoding Decode GLONASS L1OCd navigation data.
%
%   [eph, firstSubFrame, TOD] = NAVdecoding(I_P)
%
%   Input:
%       I_P           Prompt correlator output, one value per 2 ms.
%
%   Output:
%       eph           Ephemeris structure.
%       firstSubFrame Start index of the first valid L1OCd string in I_P.
%                     The unit is 2 ms prompt samples.
%       TOD           Time of day of the first valid string, in seconds.
%
%   ICD basis:
%       L1OCd data rate              : 125 bps
%       convolutional encoder         : rate 1/2, K = 7, (133,171)
%       CE symbol rate                : 250 sps
%       OC1                           : 01, 500 sps
%       normal L1OCd string length    : 250 bits, 2 s
%       preamble                      : 010111110001
%       CRC                           : CRC(250,234)

%--------------------------------------------------------------------------
% Modified for GLONASS L1OCd navigation message decoding.
%--------------------------------------------------------------------------

%% Initialization
try
    eph = eph_structure_init();
catch
    eph = struct();
end

I_P = I_P(:).';
I_P = I_P - mean(I_P);

firstSubFrame    = inf;
TOD              = inf;
searchStartOffset = 0;
gotTypes         = false(1, 3);     % Type 10, 11, 12

if numel(I_P) < 1200
    return;
end

%% ICD constants
PREAMBLE_BITS = logical([0 1 0 1 1 1 1 1 0 0 0 1]);
STR_LEN_BITS  = 250;

IDX.PREAMBLE = 1:12;
IDX.TYPE     = 13:18;
IDX.SVID     = 19:24;
IDX.TS       = 35:50;

% OC1 = 01. With antipodal mapping 0 -> +1, 1 -> -1.
OC1_PATTERN = [1 -1];
OC1_LEN     = 2;

% ICD convolutional encoder: K = 7, G = (133,171).
trellis = poly2trellis(7, [133 171]);
tblen   = 105;

%% 1. Remove OC1
% One I_P sample corresponds to 2 ms. Two adjacent samples correspond to
% one CE symbol after OC1 stripping.
[ceSymbols, bestOC1Phase] = removeOC1(I_P, OC1_PATTERN, OC1_LEN, STR_LEN_BITS);

if isempty(ceSymbols)
    return;
end

ceSymbols = ceSymbols - mean(ceSymbols);
ceSymbols = ceSymbols(1 : end - rem(numel(ceSymbols), 2));

%% 2. Viterbi decoding and string search
% Two ambiguities are kept:
%   1) 180-degree carrier phase ambiguity: CE symbols may be inverted.
%   2) CE branch alignment ambiguity: the first available CE symbol may be
%      the first or second encoder output of one data bit.
for polarity = [1 -1]

    hardCESymbols = logical((polarity * ceSymbols) < 0);

    for ceOffset = 0:1

        if ceOffset == 0
            vitInput = hardCESymbols;
            ipOffset = 0;
        else
            vitInput = hardCESymbols(2:end);
            ipOffset = 2;       % one CE symbol = two 2-ms I_P samples
        end

        vitInput = vitInput(1 : end - rem(numel(vitInput), 2));

        if numel(vitInput) < 2 * STR_LEN_BITS
            continue;
        end

        decodedBits = vitdec(vitInput, trellis, tblen, 'trunc', 'hard');
        decodedBits = logical(decodedBits(:).');

        candidateStarts = findPreambleCandidates(decodedBits, PREAMBLE_BITS, STR_LEN_BITS);

        for k = 1:numel(candidateStarts)

            startIdx = candidateStarts(k);
            navBits  = decodedBits(startIdx : startIdx + STR_LEN_BITS - 1);

            if ~checkL1OCdCRC250(navBits)
                continue;
            end

            msgType = bin2dec_unsigned(navBits(IDX.TYPE));
            svID    = bin2dec_unsigned(navBits(IDX.SVID));

            if ~ismember(msgType, [10 11 12])
                fprintf('  > Skipped L1OCd String: Type %d, SV_ID %d, CRC OK, not ephemeris.\n', ...
                        msgType, svID);
                continue;
            end

            if isinf(firstSubFrame)
                firstSubFrame = searchStartOffset + bestOC1Phase + ipOffset + ...
                                (startIdx - 1) * 4;

                TOD = bin2dec_unsigned(navBits(IDX.TS)) * 2;
                eph.TOD = TOD;
                eph.TOW = TOD;
            end

            try
                [eph, ~] = ephemeris(navBits, eph);

                switch msgType
                    case 10
                        gotTypes(1) = true;
                    case 11
                        gotTypes(2) = true;
                    case 12
                        gotTypes(3) = true;
                end

                if isfield(eph, 'SV_ID') && ~isempty(eph.SV_ID)
                    fprintf('Decoded L1OCd String: Type %d, SV_ID %d\n', ...
                            msgType, eph.SV_ID);
                else
                    fprintf('Decoded L1OCd String: Type %d, SV_ID %d\n', ...
                            msgType, svID);
                end

            catch ME
                warning('NAVdecoding:EphemerisParsingFailed', ...
                        'Ephemeris parsing failed for Type %d, SV_ID %d: %s', ...
                        msgType, svID, ME.message);
            end

            if all(gotTypes)
                return;
            end
        end
    end
end

end

%% ========================================================================
function [bestSymbols, bestPhase] = removeOC1(I_P, ocPattern, ocLen, strLenBits)
%REMOVEOC1 Search the best OC1 phase and collapse 2-ms samples to CE symbols.

bestMetric  = -inf;
bestPhase   = 1;
bestSymbols = [];

minCESymbols = 2 * strLenBits;

for phase = 1:ocLen

    x = I_P(phase:end);
    nPair = floor(numel(x) / ocLen);

    if nPair < minCESymbols
        continue;
    end

    x = reshape(x(1 : nPair * ocLen), ocLen, nPair);
    ce = ocPattern * x;

    metric = sum(abs(ce));

    if metric > bestMetric
        bestMetric  = metric;
        bestPhase   = phase;
        bestSymbols = ce;
    end
end

end

%% ========================================================================
function candidateStarts = findPreambleCandidates(bits, preambleBits, strLenBits)
%FINDPREAMBLECANDIDATES Find exact preamble matches with enough following bits.

nSearch = numel(bits) - strLenBits + 1;

if nSearch <= 0
    candidateStarts = [];
    return;
end

candidateStarts = [];

for idx = 1:nSearch
    if isequal(bits(idx : idx + numel(preambleBits) - 1), preambleBits)
        candidateStarts(end + 1) = idx; %#ok<AGROW>
    end
end

end

%% ========================================================================
function crcOK = checkL1OCdCRC250(navBits)
%CHECKL1OCDCRC250 Check normal 250-bit L1OCd string by CRC syndrome.

navBits = logical(navBits(:).');
crcOK = false;

if numel(navBits) ~= 250
    return;
end

% g(x) = x^16 + x^14 + x^13 + x^11 + x^10 + x^9 + x^8
%      + x^6 + x^5 + x + 1
poly = logical([1 0 1 1 0 1 1 1 1 0 1 1 0 0 0 1 1]);

remainder = crcRemainderMSB(navBits, poly);
crcOK = ~any(remainder);

end

%% ========================================================================
function remainder = crcRemainderMSB(codeBits, poly)
%CRCREMAINDERMSB MSB-first polynomial division over GF(2).

codeBits = logical(codeBits(:).');
poly     = logical(poly(:).');

crcLen = numel(poly) - 1;
work   = codeBits;

for k = 1:(numel(codeBits) - crcLen)
    if work(k)
        work(k : k + crcLen) = xor(work(k : k + crcLen), poly);
    end
end

remainder = work(end - crcLen + 1 : end);

end

%% ========================================================================
function val = bin2dec_unsigned(bits)
%BIN2DEC_UNSIGNED Convert an MSB-first binary field to unsigned decimal.

bits = logical(bits(:).');

if isempty(bits)
    val = 0;
else
    val = bin2dec(char(double(bits) + '0'));
end

end