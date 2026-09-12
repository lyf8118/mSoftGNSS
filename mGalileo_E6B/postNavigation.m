function [navSolutions, eph] = postNavigation(trackResults, settings)
%Decode Galileo E6B HAS pages without attempting an independent PVT fix.
%
%[navSolutions, eph] = postNavigation(trackResults, settings)
%
%   Inputs:
%       trackResults - results from the tracking function.
%       settings     - receiver settings.
%
%   Outputs:
%       navSolutions - receiver-wide decoded HAS messages.
%       eph          - receiver-wide HAS assembler state.
%
% Galileo E6B supplies HAS correction pages and Time Of Hour. It does not
% supply the independent broadcast ephemeris and GST synchronization needed
% to treat TOH as TOW or to form a standalone position solution.

%% Collect physical pages from every locked tracking channel ==============
pageRecords = struct('rxTime', {}, 'absoluteSample', {}, 'PRN', {}, ...
    'pageBits', {});
activeChnList = find([trackResults.status] ~= '-');
for channelNr = activeChnList
    PRN = trackResults(channelNr).PRN;
    fprintf('Decoding HAS pages for PRN %02d --------------------\n', PRN);
    [channelPages, ~] = NAVdecoding(trackResults(channelNr).I_P, ...
        trackResults(channelNr).absoluteSample, settings.samplingFreq, PRN);
    pageRecords = [pageRecords channelPages]; %#ok<AGROW>
    fprintf('    %d CRC-valid HAS pages found.\n', length(channelPages));
end

%% Assemble messages in true receiver-time order ==========================
navSolutions = struct('HAS', struct([]), 'pageCount', length(pageRecords), ...
    'positionAvailable', false);
eph = [];
if isempty(pageRecords)
    fprintf('No CRC-valid Galileo E6B HAS pages were found.\n');
    return;
end
[~, pageOrder] = sort([pageRecords.absoluteSample]);
pageRecords = pageRecords(pageOrder);

for pageNr = 1:length(pageRecords)
    [eph, message] = ephemeris_E6B(pageRecords(pageNr), eph);
    if isempty(message)
        continue;
    end
    fprintf(['    HAS MT%d MID %02d decoded from %d satellites; ' ...
        'TOH=%d s.\n'], message.MT, message.MID, ...
        length(message.Sources), message.TOH);
end

% HASS=3 clears completedMessages in the assembler. Build public output only
% after all time-ordered page records and status changes have been handled.
if isfield(eph, 'completedMessages')
    navSolutions.HAS = eph.completedMessages;
end

if isempty(navSolutions.HAS)
    fprintf('No complete Galileo HAS message was assembled.\n');
else
    fprintf('%d Galileo HAS messages were assembled.\n', ...
        length(navSolutions.HAS));
end
fprintf(['E6B HAS does not provide standalone broadcast ephemeris/GST; ' ...
    'position calculation was not performed.\n']);
end
