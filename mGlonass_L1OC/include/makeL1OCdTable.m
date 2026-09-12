function L1OcdCodesTable = makeL1OCdTable(settings, PRN)
%Function generates GLONASS L1OCd local codes for the specified satellite 
%based on the settings provided in the structure "settings". The codes are 
%digitized at the sampling frequency specified in the settings structure.
%One row or one vector in "L1OcdCodesTable" is one sampled L1OCd code.
%
%L1OcdCodesTable = makeL1OCdTable(settings, PRN)
%
%   Inputs:
%       settings        - receiver settings structure.
%                         The following fields are used:
%                         settings.samplingFreq   - sampling frequency [Hz]
%                         settings.codeFreqBasis  - code frequency [Hz]
%                         settings.codeLength     - code length after
%                                                   return-zero/TDM expansion
%       PRN             - GLONASS L1OC satellite ID number.
%                         According to ICD, the valid range is 0...63.
%                         PRN 0 is reserved, normally use 1...63.
%
%   Outputs:
%       L1OcdCodesTable - an array containing the sampled L1OCd local code
%                         for the specified satellite PRN.
%
%--------------------------------------------------------------------------
%                         CU Multi-GNSS SDR
% Modified for GLONASS L1OCd signal based on the L1OC local-code table structure.
%--------------------------------------------------------------------------
%
% Notes:
%   1) GLONASS L1OCd original PRN code length is 1023 chips.
%   2) Original L1OCd PRN chip rate is 0.5115 MHz and period is 2 ms.
%   3) In this receiver, the L1OCd local code is generated in return-zero /
%      TDM form, so one zero is inserted after each L1OCd PRN chip.
%   4) Therefore, the equivalent local code length becomes 1023*2 = 2046,
%      and the equivalent code frequency becomes 0.5115e6*2 = 1.023 MHz.
%   5) This function only samples the local PRN code. L1OCd navigation data,
%      convolution encoder symbols and OC1 overlay code are not generated here.
%--------------------------------------------------------------------------

%--- Find number of samples per spreading code ----------------------------
% The number of samples in one complete L1OCd code period is determined by
% the sampling frequency and the code period.
%
% For the current return-zero/TDM setting:
%   settings.codeLength    = 2046 chips
%   settings.codeFreqBasis = 1.023e6 Hz
%
% Therefore:
%   code period = settings.codeLength / settings.codeFreqBasis = 2 ms
samplesPerCode = round(settings.samplingFreq / ...
                      (settings.codeFreqBasis / settings.codeLength));

%--- Find time constants --------------------------------------------------
ts = 1 / settings.samplingFreq;      % Sampling period in seconds
tc = 1 / settings.codeFreqBasis;     % Code chip period in seconds

%--- Generate L1OCd code for given PRN ------------------------------------
% generateL1OcdCode returns the L1OCd local code in return-zero/TDM form.
% Its length should be equal to settings.codeLength, namely 2046.
L1OcdCode = generateL1OCdCode(PRN, settings);

%=== Digitizing ===========================================================

%--- Make index array to read L1OCd code values ---------------------------
% The length of the index array depends on the sampling frequency and the
% duration of one L1OCd code period.
%
% For each sampling instant, this array gives the corresponding L1OCd code
% chip index. The sampled local code is then obtained by selecting values
% from the original L1OCd code sequence.
codeValueIndex = ceil((ts * (0:samplesPerCode-1)) / tc);

%--- Correct the first and last indexes -----------------------------------
% Due to rounding, the first or last index may be slightly inaccurate.
% Force the first sample to use the first code chip, and the last sample to
% use the last code chip of one complete L1OCd period.
codeValueIndex(1)   = 1;
codeValueIndex(end) = settings.codeLength;

%--- Make the digitized version of the L1OCd code -------------------------
% The "upsampled" local code is made by selecting values from the L1OCd code
% chip array according to the sampling instants.
L1OcdCodesTable = L1OcdCode(codeValueIndex);

end