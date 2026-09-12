function [pageRecords, firstSubFrame] = NAVdecoding(I_P, absoluteSample, ...
    samplingFreq, PRN)
%Find all CRC-valid Galileo E6B C/NAV pages in one tracking channel.
%
%[pageRecords, firstSubFrame] = NAVdecoding(I_P, absoluteSample, ...
%    samplingFreq, PRN)
%
%   Inputs:
%       I_P             - prompt correlator output from tracking.
%       absoluteSample  - absolute IF-sample position of each tracking epoch.
%       samplingFreq    - IF sampling frequency in Hz.
%       PRN             - source Galileo satellite number.
%
%   Outputs:
%       pageRecords     - CRC-valid decoded C/NAV pages with comparable
%                       receiver time and source PRN.
%       firstSubFrame   - tracking-epoch index of the first valid page, or
%                       inf when no valid page is found.

% Preamble search can be delayed to avoid tracking-loop transients.
searchStartOffset = 0;
firstSubFrame = inf;
pageRecords = struct('rxTime', {}, 'absoluteSample', {}, 'PRN', {}, ...
    'pageBits', {});

% CRC-24Q detector and the non-inverted convolutional-code trellis.
crcDet = comm.CRCDetector([24 23 18 17 14 11 10 7 6 5 4 3 1 0]);
trellis = poly2trellis(7, [171 133]);
tblen = 35;

% Antipodal form of the 16-symbol E6B preamble.
preambleBits = [-1 1 -1 -1 1 -1 -1 -1 1 -1 -1 -1 1 1 1 1];
bits = I_P(1 + searchStartOffset:end);
bits(bits > 0) = 1;
bits(bits <= 0) = -1;

% Find preamble candidates and verify the expected one-second spacing.
tlmXcorrResult = xcorr(bits, preambleBits);
xcorrLength = (length(tlmXcorrResult) + 1)/2;
index = find(abs(tlmXcorrResult(xcorrLength:xcorrLength*2-1)) > 15)';

for ind = 1:length(index)
    if ~any(index-index(ind) == 1000)
        continue;
    end
    if index(ind)+999 > length(bits)
        continue;
    end

    navBits = bits(index(ind):index(ind)+999);
    if ~isequal(navBits(1:16), preambleBits)
        navBits = -navBits;
    end

    % Fill interleaver columns, read rows, then invert the G2 branch.
    pageSym = reshape(reshape(navBits(17:1000), 123, 8)', 1, [])';
    pageSym = (1-pageSym)/2;
    pageSym(2:2:end) = 1-pageSym(2:2:end);

    % Remove the rate-1/2 convolutional code and its six zero-tail bits.
    decBits = vitdec(pageSym, trellis, tblen, 'trunc', 'hard');
    decBits = decBits(1:486);
    [~, frmError] = step(crcDet, decBits);
    if frmError
        continue;
    end

    pageIndex = index(ind)+searchStartOffset;
    pageSample = double(absoluteSample(pageIndex));
    pageRecords(end+1) = struct( ...
        'rxTime', pageSample/double(samplingFreq), ...
        'absoluteSample', pageSample, 'PRN', PRN, ...
        'pageBits', decBits(:)'); %#ok<AGROW>
    if isinf(firstSubFrame)
        firstSubFrame = pageIndex;
    end
end
end
