function [navSolutions, eph] = postNavigation(trackResults, settings)
%POSTNAVIGATION GLONASS L2OC tracking-only post processing.
%
% Current receiver tracks GLONASS L2OCp pilot component only.
% L2 CSI navigation message is not decoded, so pseudorange-based
% navigation solution is not available in this mode.

activeChnList = find([trackResults.status] ~= '-');

% Initialize ephemeris placeholder array
eph = repmat(eph_structure_init(), 1, settings.numberOfChannels);

for channelNr = activeChnList

    PRN = trackResults(channelNr).PRN;

    eph(PRN) = eph_structure_init();
    eph(PRN).PRN = PRN;
    eph(PRN).SVID = PRN;

    % L2OCp tracking-only status
    eph(PRN).navDecoded = 0;
    eph(PRN).flag = 0;
    eph(PRN).decodeStatus = ...
        'L2OCp pilot tracking only; L2 CSI navigation message is not decoded.';

end

navSolutions = [];
disp('GLONASS L2OCp tracking-only mode: no L2 CSI ephemeris decoded, navigation solution skipped.');

return