function L1CPilot = generatePilotTMBOC61(settings,PRN)
% Generate the GPS L1C pilot-channel TMBOC(6,1,4/33) code in bipolar
% format (-1, +1). PRNs 1 to 63 are supported.
%
% L1CPilot = generatePilotTMBOC61(settings,PRN)
%
%   Inputs:
%       PRN       - Satellite PRN number.
%
%   Outputs:
%       L1CPilot  - L1C pilot-channel code after TMBOC modulation.
%--------------------------------------------------------------------------
%                         CU Multi-GNSS SDR
% Original authors: Yafeng Li, Nagaraj C. Shivaramaiah and Dennis M. Akos
% Reference: Adapted within the CU Multi-GNSS SDR receiver framework for GPS L1C.
%--------------------------------------------------------------------------

% L1C pilot primary-code table. Column 1 is the Weil phase difference w;
% column 2 is the 7-chip extension insertion point p.
wp_pilot = [5111    412;   5109    161;   5108   1;   5106   303;...
      5103   207;   5101   4971;   5100   4496;   5098   5;...
      5095   4557;   5094    485;   5093   253;   5091   4676;...
      5090   1;   5081   66;   5080    4485;   5069   282;...
      5068   193;  5054    5211;   5044   729;   5027   4848;...
      5026   982;   5014   5955;   5004   9805;   4980   670;...
      4915   464;  4909   29;   4893   429;   4885    394;...
      4832   616;   4824   9457;    4591   4429;   3706   4771;...
      5092   365;   4986   9705;   4965   9489;   4920   4193;...
      4917    9947;   4858   824;   4847    864;   4790   347;...
      4770   677;   4318   6544;   4126   6312;   3961   9804;...
      3790    278;   4911    9461;   4881   444;   4827   4839;...
      4795   4144;   4789   9875;   4725   197;   4675   1156;...
      4539   4674;   4535    10035;   4458   4504;   4197    5;...
      4096   9937;     3484   430;   3481   5;   3393   355;...
      3175   909;   2360   1622;   1852   6284];

% Compute the Legendre/Jacobi sequence used as the base sequence for Weil codes.
N = 10223;
legendre = zeros(1,N);
for ind = 1:N-1
    legendre(ind+1) = JacobiSymbol(ind,N);
end

% Convert {-1,+1} to {0,1} so xor can be used for Weil-code generation.
legendre(legendre==-1) = 0;
Primary = zeros(1,10223);

% Read the PRN-specific phase difference and insertion point.
p = wp_pilot(PRN,2);
w = wp_pilot(PRN,1);

% Pilot-channel Weil code: W(n) = L(n) xor L(n+w).
for ind = 0:10222
    k = ind;
    Primary(ind+1) = xor(legendre(k+1), legendre(mod(k+w,N)+1));
end

% Convert to bipolar format and insert the 7-chip extension sequence.
Primary = 1 - 2*Primary;
extendedSequence = [0, 1, 1, 0, 1, 0, 0];
extendedSequence = 1 - 2*extendedSequence;
L1CWithExtension = [Primary(1:p-1) extendedSequence Primary(p:end)];

%% Apply TMBOC(6,1,4/33) modulation.
% Use 12 samples per chip so BOC(1,1) and BOC(6,1) can share one grid.
samplesPerChip = 12;
numChips = length(L1CWithExtension);
totalLen = numChips * samplesPerChip;
L1CPilot = zeros(1, totalLen);

% BOC(1,1): first half-chip is +1 and second half-chip is -1.
subcarrierBOC11 = [ones(1, 6), -ones(1, 6)];

% BOC(6,1): six alternating +1/-1 subcarrier cycles per chip.
subcarrierBOC61 = repmat([1, -1], 1, 6);

% ICD TMBOC pattern: 4 of every 33 chips use BOC(6,1); the rest use BOC(1,1).
% These MATLAB indices correspond to 0-based positions [0,4,6,29].
tm_indices = [1, 5, 7, 30];

for jj = 1:numChips
    posInCycle = mod(jj-1, 33) + 1;

    if ismember(posInCycle, tm_indices)
        currentSub = subcarrierBOC61;
    else
        currentSub = subcarrierBOC11;
    end

    idx_start = (jj-1) * samplesPerChip + 1;
    idx_end   = jj * samplesPerChip;
    L1CPilot(idx_start:idx_end) = L1CWithExtension(jj) * currentSub;
end

end
