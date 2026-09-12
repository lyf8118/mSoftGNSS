function [eph, message] = ephemeris_E6B(pageRecord, eph)
%Assemble receiver-wide Galileo E6B HAS messages from decoded C/NAV pages.
%
%[eph, message] = ephemeris_E6B(pageRecord, eph)
%
%   Inputs:
%       pageRecord - one CRC-valid page returned by NAVdecoding.
%       eph        - receiver-wide HAS decoder state shared by all channels.
%
%   Outputs:
%       eph        - updated HAS decoder state.
%       message    - one completed MT1 message, or [] when incomplete.

message = [];
if nargin < 2 || isempty(eph)
    eph = struct();
end
if ~isfield(eph, 'flag'), eph.flag = 0; end
if ~isfield(eph, 'TOH'), eph.TOH = inf; end
if ~isfield(eph, 'referenceGST'), eph.referenceGST = inf; end
if ~isfield(eph, 'HAS'), eph.HAS = []; end
if ~isfield(eph, 'HASS'), eph.HASS = []; end
if ~isfield(eph, 'HAS_Buffer'), eph.HAS_Buffer = empty_HAS_buffer(); end
if ~isfield(eph, 'completedMessages'), eph.completedMessages = struct([]); end

%% Parse the 24-bit HAS page header ========================================
pageBits = double(pageRecord.pageBits(:)');
hasPage = pageBits(15:462);
header = hasPage(1:24);
if bits2uint(header) == hex2dec('AF3BC3')
    return;
end

HASS = bits2uint(header(1:2));
% HASS=3 means that all previously received HAS messages shall be discarded.
if HASS == 3
    eph.HAS_Buffer = empty_HAS_buffer();
    eph.completedMessages = struct([]);
    eph.HAS = [];
    eph.flag = 0;
    eph.TOH = inf;
    eph.referenceGST = inf;
    eph.HASS = HASS;
    return;
elseif HASS == 2
    return;
end

MT = bits2uint(header(5:6));
MID = bits2uint(header(7:11));
MSraw = bits2uint(header(12:16));
PID = bits2uint(header(17:24));
if PID < 1 || PID > 255
    return;
end

% Do not combine test-mode and operational-mode encoded pages.
if ~isempty(eph.HASS) && any(eph.HASS == [0 1]) && eph.HASS ~= HASS
    eph.HAS_Buffer = empty_HAS_buffer();
    eph.HAS = [];
    eph.flag = 0;
    eph.TOH = inf;
    eph.referenceGST = inf;
end
eph.HASS = HASS;
if MT ~= 1
    return;
end

%% Expire and select the receiver-wide message buffer =====================
if ~isempty(eph.HAS_Buffer)
    firstRxTime = [eph.HAS_Buffer.firstRxTime];
    eph.HAS_Buffer(pageRecord.rxTime-firstRxTime >= 150.0) = [];
end
if ~isempty(eph.HAS_Buffer)
    conflict = [eph.HAS_Buffer.MT] == MT & ...
        [eph.HAS_Buffer.MID] == MID & ...
        [eph.HAS_Buffer.HASS] == HASS & ...
        [eph.HAS_Buffer.MSraw] ~= MSraw;
    eph.HAS_Buffer(conflict) = [];
end

bufferIndex = find([eph.HAS_Buffer.MT] == MT & ...
    [eph.HAS_Buffer.MID] == MID & ...
    [eph.HAS_Buffer.MSraw] == MSraw & ...
    [eph.HAS_Buffer.HASS] == HASS, 1);
if isempty(bufferIndex)
    bufferIndex = length(eph.HAS_Buffer)+1;
    eph.HAS_Buffer(bufferIndex).MT = MT;
    eph.HAS_Buffer(bufferIndex).MID = MID;
    eph.HAS_Buffer(bufferIndex).MSraw = MSraw;
    eph.HAS_Buffer(bufferIndex).HASS = HASS;
    eph.HAS_Buffer(bufferIndex).K = MSraw+1;
    eph.HAS_Buffer(bufferIndex).firstRxTime = pageRecord.rxTime;
    eph.HAS_Buffer(bufferIndex).lastRxTime = pageRecord.rxTime;
    eph.HAS_Buffer(bufferIndex).receivedMask = false(1, 255);
    eph.HAS_Buffer(bufferIndex).pages = zeros(255, 53);
    eph.HAS_Buffer(bufferIndex).sources = [];
end

encodedBits = hasPage(25:end);
encodedBytes = zeros(1, 53);
for byteNr = 1:53
    encodedBytes(byteNr) = bits2uint( ...
        encodedBits((byteNr-1)*8+1:byteNr*8));
end
if ~eph.HAS_Buffer(bufferIndex).receivedMask(PID)
    eph.HAS_Buffer(bufferIndex).pages(PID, :) = encodedBytes;
    eph.HAS_Buffer(bufferIndex).receivedMask(PID) = true;
end
eph.HAS_Buffer(bufferIndex).lastRxTime = pageRecord.rxTime;
eph.HAS_Buffer(bufferIndex).sources = unique( ...
    [eph.HAS_Buffer(bufferIndex).sources pageRecord.PRN]);

%% Recover the message after K distinct Page IDs have arrived =============
K = eph.HAS_Buffer(bufferIndex).K;
if nnz(eph.HAS_Buffer(bufferIndex).receivedMask) < K
    return;
end
validPID = find(eph.HAS_Buffer(bufferIndex).receivedMask);
usePIDs = validPID(1:K);
rxBlock = eph.HAS_Buffer(bufferIndex).pages(usePIDs, :);
msgBytes = rs_erasure_decode(rxBlock, usePIDs, K, ...
    load_RS_Generator_Matrix(), init_gf256());

% MATLAB column-major order requires this transpose to preserve page order.
flatBytes = msgBytes';
flatBytes = flatBytes(:);
bitsStr = dec2bin(flatBytes, 8)';
fullMsgBits = bitsStr(:)'-'0';

% Parse the fixed MT1 header. Dynamic correction blocks remain raw and are
% not applied until a matching Mask ID/IOD Set ID context is available.
message.TOH = bits2uint(fullMsgBits(1:12));
message.Flags = struct( ...
    'Mask', logical(fullMsgBits(13)), ...
    'Orbit', logical(fullMsgBits(14)), ...
    'ClockFull', logical(fullMsgBits(15)), ...
    'ClockSubset', logical(fullMsgBits(16)), ...
    'CodeBias', logical(fullMsgBits(17)), ...
    'PhaseBias', logical(fullMsgBits(18)));
message.Reserved = bits2uint(fullMsgBits(19:22));
message.MaskID = bits2uint(fullMsgBits(23:27));
message.IODSetID = bits2uint(fullMsgBits(28:32));
message.RawBody = fullMsgBits(33:end);
message.BlocksDecoded = false;
message.HASS = HASS;
message.Operational = HASS == 1;
message.MT = MT;
message.MID = MID;
message.MSraw = MSraw;
message.K = K;
message.PIDs = usePIDs;
message.Sources = eph.HAS_Buffer(bufferIndex).sources;
message.FirstRxTime = eph.HAS_Buffer(bufferIndex).firstRxTime;
message.LastRxTime = eph.HAS_Buffer(bufferIndex).lastRxTime;
message.Data = fullMsgBits;

eph.TOH = double(message.TOH);
eph.HAS = message;
eph.flag = 1;
if isempty(eph.completedMessages)
    eph.completedMessages = message;
else
    eph.completedMessages(end+1) = message;
end
eph.HAS_Buffer(bufferIndex) = [];
end

function buffer = empty_HAS_buffer()
buffer = struct('MT', {}, 'MID', {}, 'MSraw', {}, 'HASS', {}, 'K', {}, ...
    'firstRxTime', {}, 'lastRxTime', {}, 'receivedMask', {}, 'pages', {}, ...
    'sources', {});
end

function value = bits2uint(bits)
value = 0;
for bitNr = 1:length(bits)
    value = 2*value+double(bits(bitNr));
end
end

function T = init_gf256()
persistent tables
if isempty(tables)
    tables.exp = zeros(1, 512);
    tables.log = zeros(1, 256);
    value = 1;
    for exponent = 0:254
        tables.exp(exponent+1) = value;
        tables.log(value+1) = exponent;
        value = bitshift(value, 1);
        if bitand(value, 256)
            value = bitxor(value, 285);
        end
    end
    tables.exp(256:511) = tables.exp(1:256);
end
T = tables;
end

function G = load_RS_Generator_Matrix()
persistent generator
if isempty(generator)
    matrixName = ['Galileo-HAS-SIS-ICD_1.0_Annex_B_' ...
        'Reed_Solomon_Generator_Matrix.txt'];
    matrixPath = fullfile(fileparts(mfilename('fullpath')), matrixName);
    generator = readmatrix(matrixPath);
    if ~isequal(size(generator), [255 32]) || ...
            any(generator(:) < 0 | generator(:) > 255) || ...
            ~isequal(generator(1:32, :), eye(32))
        error('Invalid Galileo HAS RS generator matrix.');
    end
end
G = generator;
end

function value = gf_mul(a, b, T)
if a == 0 || b == 0
    value = 0;
else
    value = T.exp(T.log(a+1)+T.log(b+1)+1);
end
end

function value = gf_inv(a, T)
if a == 0
    error('GF(256) division by zero.');
end
value = T.exp(255-T.log(a+1)+1);
end

function decodedMsg = rs_erasure_decode(rxBlock, pids, K, G, T)
D = G(pids, 1:K);
D_inv = eye(K);
for row = 1:K
    if D(row, row) == 0
        swapRow = find(D(row+1:end, row) ~= 0, 1)+row;
        if isempty(swapRow)
            error('Singular matrix in RS decode.');
        end
        D([row swapRow], :) = D([swapRow row], :);
        D_inv([row swapRow], :) = D_inv([swapRow row], :);
    end
    invPivot = gf_inv(D(row, row), T);
    for column = 1:K
        D(row, column) = gf_mul(D(row, column), invPivot, T);
        D_inv(row, column) = gf_mul(D_inv(row, column), invPivot, T);
    end
    for otherRow = 1:K
        if otherRow == row
            continue;
        end
        factor = D(otherRow, row);
        if factor ~= 0
            for column = 1:K
                D(otherRow, column) = bitxor(D(otherRow, column), ...
                    gf_mul(factor, D(row, column), T));
                D_inv(otherRow, column) = bitxor(D_inv(otherRow, column), ...
                    gf_mul(factor, D_inv(row, column), T));
            end
        end
    end
end

decodedMsg = zeros(K, 53);
for row = 1:K
    for column = 1:53
        value = 0;
        for index = 1:K
            value = bitxor(value, gf_mul(D_inv(row, index), ...
                rxBlock(index, column), T));
        end
        decodedMsg(row, column) = value;
    end
end
end
