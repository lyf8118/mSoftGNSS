function [decodedCodeword, realIter] = msDecode(llr, hMatrix, settings)
%MSDECODE Min-Sum LDPC decoder for GPS L1C CNAV-2.
%   decodedCodeword = MSDECODE(llr, hMatrix, settings) decodes the input
%   LLR vector using the supplied parity-check matrix. This keeps the SDR
%   receiver on its existing H2/H3 matrices while using the Min-Sum message
%   update logic from the standalone GPS L1C decoder.

if nargin == 2
    settings = hMatrix;
    if ~isfield(settings, 'hMatrix')
        error('msDecode:MissingHMatrix', ...
            'Call msDecode with hMatrix or provide settings.hMatrix.');
    end
    hMatrix = settings.hMatrix;
end

llr = llr(:);
hMatrix = sparse(logical(hMatrix));
[numRows, numCols] = size(hMatrix);

if length(llr) ~= numCols
    error('msDecode:InputSize', ...
        'LLR length (%d) does not match LDPC matrix columns (%d).', ...
        length(llr), numCols);
end

maxIterations = getSetting(settings, 'maxIterations', 50);
alpha = getSetting(settings, 'alpha', 1);
offset = getSetting(settings, 'offset', 0);

useOffset = offset > 0;
useAlpha = alpha > 0 && alpha < 1;

persistent graphCache
if isempty(graphCache)
    graphCache = containers.Map();
end

messageType = getSetting(settings, 'messageType', 'CNAV2');
cacheKey = sprintf('%s_%dx%d_%d', messageType, numRows, numCols, nnz(hMatrix));

if isKey(graphCache, cacheKey)
    graph = graphCache(cacheKey);
else
    [rowIndices, colIndices] = find(hMatrix);
    numEdges = length(rowIndices);
    graph.hMatrix = hMatrix;
    graph.colIndices = colIndices;
    graph.cnConnections = accumarray(rowIndices, (1:numEdges).', ...
        [numRows 1], @(y){y});
    graph.numEdges = numEdges;
    graphCache(cacheKey) = graph;
end

colIndices = graph.colIndices;
cnConnections = graph.cnConnections;
numEdges = graph.numEdges;
hMatrix = graph.hMatrix;

if maxIterations == 0
    decodedCodeword = double(llr < 0);
    realIter = 0;
    return;
end

chkToVarMsg = zeros(numEdges, 1);
varToChkMsg = llr(colIndices);
decodedCodeword = double(llr < 0);
realIter = maxIterations;

for iter = 1:maxIterations
    for cn = 1:numRows
        edgeIndices = cnConnections{cn};
        degree = length(edgeIndices);

        if degree < 2
            chkToVarMsg(edgeIndices) = 0;
            continue;
        end

        qMessages = varToChkMsg(edgeIndices);

        signVec = sign(qMessages);
        signVec(signVec == 0) = 1;
        signProduct = prod(signVec);

        absQ = abs(qMessages);
        [min1, minIdx] = min(absQ);
        absQ(minIdx) = inf;
        min2 = min(absQ);

        if useOffset
            min1 = max(min1 - offset, 0);
            min2 = max(min2 - offset, 0);
        elseif useAlpha
            min1 = min1 * alpha;
            min2 = min2 * alpha;
        end

        magnitudes = min1 * ones(degree, 1);
        magnitudes(minIdx) = min2;
        chkToVarMsg(edgeIndices) = (signProduct .* signVec) .* magnitudes;
    end

    sumChkToVar = accumarray(colIndices, chkToVarMsg, [numCols 1], @sum, 0);
    posteriorLlr = llr + sumChkToVar;
    decodedCodeword = double(posteriorLlr < 0);

    if ~any(mod(hMatrix * decodedCodeword, 2))
        realIter = iter;
        return;
    end

    varToChkMsg = posteriorLlr(colIndices) - chkToVarMsg;
end
end

function value = getSetting(settings, fieldName, defaultValue)
if isstruct(settings) && isfield(settings, fieldName)
    value = settings.(fieldName);
else
    value = defaultValue;
end
end
