function L1CData = generateDataBOC11(settings,PRN)
% Generate the GPS L1C data-channel BOC(1,1) primary code in bipolar
% format (-1, +1). PRNs 1 to 63 are supported.
%
% L1CData = generateDataBOC11(settings,PRN)
%
%   Inputs:
%       PRN       - Satellite PRN number.
%
%   Outputs:
%       L1CData   - L1C data-channel code after extension insertion and
%                   BOC(1,1) subcarrier modulation.

%--------------------------------------------------------------------------
%                         CU Multi-GNSS SDR
% Original authors: Yafeng Li, Nagaraj C. Shivaramaiah and Dennis M. Akos
% Reference: Adapted within the CU Multi-GNSS SDR receiver framework for GPS L1C.
%--------------------------------------------------------------------------

% L1C data primary-code table. Column 1 is the Weil phase difference w;
% column 2 is the 7-chip extension insertion point p.
wp_data = [5097   181;    5110   359;    5079    72;   4403    1110;...
      4121   1480;    5043   5034;   5042    4622;   5104   1;...
      4940   4547;   5035   826;   4372   6284;   5064   4195;...
      5084   368;   5048    1;   4950   4796;     5019   523;...
      5076   151;   3736   713;   4993    9850;   5060   5734;...
      5061   34;   5096   6142;   4983    190;   4783   644;...
      4991   467;   4815   5384;   4443   801;   4769   594;...
      4879   4450;    4894   9437;   4985   4307;   5056   5906;...
      4921   378;   5036   9448;   4812   9432;   4838    5849;...
      4855   5547;   4904    9546;   4753   9132;   4483    403;...
      4942   3766;   4813     3;   4957   684;    4618    9711;...
      4669   333;   4969   6124;   5031   10216;   5038   4251;...
      4740   9893;   4073    9884;   4843    4627;   4979    4449;...
      4867    9798;   4964   985;    5025   4272;   4579   126;...
      4390   10024;    4763   434;   4612   1029;   4784    561;...
      3716   289;   4703    638;   4851   4353];

% Compute the Legendre/Jacobi sequence used as the base sequence for Weil codes.
N = 10223;
legendre = zeros(1,N);
for ind = 1:N-1
    legendre(ind+1) = JacobiSymbol(ind,N);
end

% Convert {-1,+1} to {0,1} so xor can be used for Weil-code generation.
legendre(legendre==-1) = 0;
Primary  = zeros(1,10223);

% Read the PRN-specific phase difference and insertion point.
p = wp_data(PRN,2);
w = wp_data(PRN,1);

% Data-channel Weil code: W(n) = L(n) xor L(n+w).
for ind = 0:10222
    k = ind;
    Primary(ind+1) = xor(legendre(k+1), legendre(mod(k+w,N)+1));
end

% Convert to bipolar format: 0 -> +1, 1 -> -1.
Primary = 1 - 2*Primary;

%% Insert the 7-chip extension sequence and apply the BOC(1,1) subcarrier.
extendedSequence = [0, 1, 1, 0, 1, 0, 0];
extendedSequence = 1 - 2*extendedSequence;
L1CWithExtension = [Primary(1:p-1) extendedSequence Primary(p:end)];

L1CData = zeros(1, length(L1CWithExtension)*2);
jj = 1;
for ii = 1:2:length(L1CData)-1
    L1CData(ii)   =  L1CWithExtension(jj);
    L1CData(ii+1) = -L1CWithExtension(jj);
    jj = jj + 1;
end
