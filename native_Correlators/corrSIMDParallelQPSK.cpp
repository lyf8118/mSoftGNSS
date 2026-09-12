/*================================================================================
 * Filename: corrSIMDParallelQPSK.cpp
 * Description: SIMD-assisted channel-parallel QPSK tracking correlator for
 *              trkChannelsParallel.m.
 *
 * Authors: Yafeng Li (School of Automation, Beijing Information Science and Technology University)
 * Time: Aug, 1, 2026
 *
 * The core processing stages are intentionally organized as:
 *   1) mixCarrReal() / mixCarrComplex()
 *   2) localCodeGen()
 *   3) correlator()
 *
 * The MEX gateway framework in this file processes all tracking channels
 * for one code period per call, matching trkChannelsParallel.m.
 *==============================================================================*/
 
#include "mex.h"
#include <math.h>
#include <string.h>
/* ------- SIMD ------- */
#include <emmintrin.h>
#include <tmmintrin.h>
#include <immintrin.h>

/* Settings struct -------------------------------------------------------------
 * This structure stores the receiver settings required by the SIMD correlator.
 * trkMode        : 0 - channel-serial tracking; 1 - channel-parallel tracking
 * fileType       : 1 - real samples; 2 - complex samples
 * correlatorType : 0 - Matlab correlator; 1 - SIMD correlator; 2 - GPU correlator
 * rShiftBits     : number of right-shift bits for IF data to prevent SIMD
 *                  correlator overflow when ADC valid bits occupy the high bits
 * earlyLateSpc   : half of the early - late code correlation spacing(chips)
 *------------------------------------------------------------------------------*/
typedef struct
{
	int  trkMode;        /* 0 - Channel-serial tracking; 1 - Channel-parallel tracking */
	int  fileType;       /* 1 - real samples; 2 - complex samples */
	int  correlatorType;  /* 0 - Matlab correlator; 1 - SIMD correlator; 2 - GPU correlator */ 
	int  rShiftBits;     /* Right-shift bits for IF data before SIMD processing */
	double earlyLateSpc; /* half of the early - late code correlation spacing(chips) */
} Settings;

/* Function declaration ---------------------------------------------------------
 * mixCarrReal()        : mix local carrier with real IF data
 * mixCarrComplex()     : mix local carrier with complex IF data
 * localCodeGen()       : generate data/pilot early/prompt/late local code sequences
 * correlator()         : compute one branch of IE/QE/IP/QP/IL/QL correlation outputs
 * cleanup()            : Release all persistent host buffers when the MEX function is cleared
 * checkedMalloc()      : Allocate non-aligned host memory and stop execution if allocation fails
 * checkedAlignedMalloc(): Allocate 32-byte aligned host memory for SIMD buffers
 * gatherCode16()       : Gather 16 local-code samples and repack them to int16 order
 * parseSettings()      : Extract the required settings fields from the MATLAB settings structure
 * getStructIntField()  : Read a scalar numeric field from the MATLAB settings struct
 * getStructDoubleField():Read a scalar numeric field from the MATLAB settings struct
 * checkInputs()        : Validate the number, type and basic form of the MEX input arguments
 *------------------------------------------------------------------------------*/
void mixCarrReal(const short*, double, double, mwSize, short*, short*);
void mixCarrComplex(const short*, const short*, double, double, mwSize, short*, short*);
void localCodeGen(const int*, mwSize, double, double, double, mwSize, short**);
void correlator(short*, short*, short**, mwSize, double*);
void cleanup(void);
static void* checkedMalloc(size_t nbytes, const char* name);
static void* checkedAlignedMalloc(size_t nbytes, const char* name);
static __m256i gatherCode16(const int* caCode, __m256i Index_reg1, __m256i Index_reg2);
static Settings parseSettings(const mxArray* s);
static int getStructIntField(const mxArray* s, const char* name);
static double getStructDoubleField(const mxArray* s, const char* name);
static void checkInputs(int nlhs, int nrhs, const mxArray* prhs[]);

/* Constants -------------------------------------------------------------------- */
#define PI          3.1415926535897932  // pi 
#define DPI         (2.0*PI)            // 2*pi 
#define CSCALE      (1.0/16.0)          // carrier lookup table scale (LSB) 

/* ------------- To preserve variables between MEX function calls  -------------
 * initialized     : initialization flag for persistent buffers
 * iBasebandSignal : in-phase baseband samples after carrier wipeoff
 * qBasebandSignal : quadrature baseband samples after carrier wipeoff
 * rawSignalI/Q    : split I/Q branches for complex input data
 * localCode       : local code buffers for data/pilot Early, Prompt and Late replicas
 * settings        : cached receiver settings parsed from MATLAB struct
 *------------------------------------------------------------------------------*/ 
static int initialized = 0;
static short *iBasebandSignal = NULL, *qBasebandSignal = NULL, *rawSignalI = NULL, *rawSignalQ = NULL;
static short ** localCode = NULL;
static Settings settings;

 /* The gateway function --------------------------------------------------------
 * Args   : Settings settings         I   receiver settings used by the SIMD correlator
 *          short  *rawSignal         I   Input IF signal
 *          int    *caCodeTable       I   QPSK PRN code lookup table arranged as
 *                                     [dataCodeTable; pilotCodeTable]
 *          double remCarrPhase       I   initial phase (rad)
 *          double carrPhaseStep      I   carrier sampling interval (cycle: Hz*t)
 *          double remCodePhase       I   initial code phase(chip)
 *          double codePhaseStep      I   code sampling interval(chip)
 *          int    *startIndex        I   Sample start index of each channel within the 'rawSignal'
 *          int    *chSampleSize      I   Sample sizes of code periods for all channel
 *          int    isDataRead         I   data-read flag
 * Return : double *corrValues: outputs arranged as
 *                               [I_E Q_E I_P Q_P I_L Q_L
 *                                pilot_I_E pilot_Q_E pilot_I_P
 *                                pilot_Q_P pilot_I_L pilot_Q_L]
 *------------------------------------------------------------------------------*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	/* --------------- Declare all variables --------------- */
	short* rawSignal;
	int* caCodeTable;
	int *startIndex, * chSampleSize;
	double* remCarrPhase, * carrPhaseStep, * remCodePhase, * codePhaseStep;
	double *corrValues;
	static size_t channelCnt;
	static size_t codeLen;
	int isDataRead;
	static size_t sampleIQSize; /* Half of the total sample count if rawSignal is complex IQ data */

	if (!initialized)
	{
		/* Validate the MEX interface before dereferencing prhs[k]. Keep this
		   check in the first-time initialization path so later calls stay on the
		   existing fast path, but still guarantee the first call is checked before
		   any input pointers are unpacked. */
		checkInputs(nlhs, nrhs, prhs);

		// Basic setting 
		settings = parseSettings(prhs[0]);

		/* --------------- Find array dimensions --------------- */
		// Size of one branch code table ([lastChip code firstChip]).
		codeLen = mxGetM(prhs[2]) / 2;
		//Channel count
		channelCnt = mxGetN(prhs[2]);
		// Sample size for the 1st channel
		chSampleSize = (int*)mxGetData(prhs[8]);
		size_t blksize = chSampleSize[0];

		if (settings.fileType == 2) {
			sampleIQSize = mxGetNumberOfElements(prhs[1])/2;
			rawSignalI = (short*)checkedAlignedMalloc(sizeof(short) * sampleIQSize, "rawSignalI");
			rawSignalQ = (short*)checkedAlignedMalloc(sizeof(short) * sampleIQSize, "rawSignalQ");
		}

		/* --------------------- Allocate memory --------------------- */
		/* I/Q baseband signals with carrier wiped off */
		iBasebandSignal = (short*)checkedAlignedMalloc(sizeof(short) * (blksize + 100), "iBasebandSignal");
		qBasebandSignal = (short*)checkedAlignedMalloc(sizeof(short) * (blksize + 100), "qBasebandSignal");

		/* For local sampling code buffers:
		   data Early/Prompt/Late plus pilot Early/Prompt/Late. */
		localCode = (short**)checkedMalloc(sizeof(short*) * 6, "localCode");
		for (int i = 0; i < 6; i++) {
			localCode[i] = (short*)checkedAlignedMalloc(sizeof(short) * (blksize + 100), "localCode[i]");
		}

		mexAtExit(cleanup);
		initialized = 1;
		mexPrintf("   MEX-file initialized ...\n");
		mexPrintf("   Channel-parallel tracking with the SIMD-accelerated correlator developed by Yafeng Li.\n");
	}

    /* ----------------- Get input arrays -------------------
     * Input order:
     *   prhs[0]  - settings
     *   prhs[1]  - rawSignal
     *   prhs[2]  - caCodeTable
     *   prhs[3]  - remCarrPhase
     *   prhs[4]  - carrPhaseStep
     *   prhs[5]  - remCodePhase
     *   prhs[6]  - codePhaseStep
     *   prhs[7]  - startIndex
     *   prhs[8]  - chSampleSize
     *   prhs[9] - isDataRead
     *--------------------------------------------------------- */
	rawSignal = (short *)mxGetData(prhs[1]);
	caCodeTable = (int *)mxGetData(prhs[2]);
	remCarrPhase = mxGetPr(prhs[3]);
	carrPhaseStep = mxGetPr(prhs[4]);
	remCodePhase = mxGetPr(prhs[5]);
	codePhaseStep = mxGetPr(prhs[6]);
	startIndex = (int*)mxGetData(prhs[7]);
	chSampleSize = (int*)mxGetData(prhs[8]);
	isDataRead = (int)mxGetScalar(prhs[9]);

	/* --------------- Initialize output array ---------------
	   Rows 1..6  : data branch   [I_E Q_E I_P Q_P I_L Q_L]
	   Rows 7..12 : pilot branch  [I_E Q_E I_P Q_P I_L Q_L] */
	plhs[0] = mxCreateDoubleMatrix((mwSize)12, (mwSize)channelCnt, mxREAL);
	corrValues = (double*)mxGetData(plhs[0]);

	/* ----- Separate "rawSignal" into I and Q branches for complex signals ------ */
	if (isDataRead == 1 && settings.fileType == 2) 	{
		for (size_t k = 0; k < sampleIQSize; ++k) {
			rawSignalI[k] = rawSignal[2 * k];
			rawSignalQ[k] = rawSignal[2 * k + 1];
		}
	}

	mwIndex dims[2] = { 0, 0 };
	for (int chInd = 0; chInd < channelCnt; ++chInd)
	{
		dims[1] = chInd;
		size_t index = mxCalcSingleSubscript(prhs[2], (size_t)2, dims);
		int* caCode = &caCodeTable[index];

		/* -------- Carrier wiping off, code generation and correlating -------- */
		if (settings.fileType == 1)
		{ // for real data
			mixCarrReal(rawSignal + startIndex[chInd] - 1, carrPhaseStep[chInd], \
				remCarrPhase[chInd], chSampleSize[chInd], iBasebandSignal, qBasebandSignal);
		}
		else if (settings.fileType == 2)
		{ // for complex data
			mixCarrComplex(rawSignalI + startIndex[chInd] - 1, rawSignalQ + startIndex[chInd] - 1, carrPhaseStep[chInd], \
				remCarrPhase[chInd], chSampleSize[chInd], iBasebandSignal, qBasebandSignal);
		}
		
		// Generate data/pilot local code sequences
		localCodeGen(caCode, (mwSize)codeLen, remCodePhase[chInd],    \
			codePhaseStep[chInd], settings.earlyLateSpc, chSampleSize[chInd], localCode);
		// Data branch correlator
		correlator(iBasebandSignal, qBasebandSignal, localCode, chSampleSize[chInd], corrValues + 12 * chInd);
		// Pilot branch correlator
		correlator(iBasebandSignal, qBasebandSignal, localCode+3, chSampleSize[chInd], corrValues+ 12 * chInd +6);
	}
}

/* Mix local carrier for real data -----------------------------------------------
 * Mix local carrier to input signal
 * Args   : short   *rawSignal         I   Input IF signal
 *          double carrPhaseStep       I   carrier sampling interval (cycle: rad*t)
 *          double remCarrPhase        I   initial carrier phase (rad)
 *          mwSize blksize             I   number of samples to be generated
 *          short  *iBasebandSignal    O   I component of input signal with carrier wiped off 
 *          short  *qBasebandSignal    O   Q component of input signal with carrier wiped off
 * Return : None
 *------------------------------------------------------------------------------*/
void mixCarrReal(const short *rawSignal, double carrPhaseStep, double remCarrPhase, mwSize blksize,
	         short *iBasebandSignal, short *qBasebandSignal)
{
	
	/* ---------- Inphase and quadrature carrier lookup table --------- */
	static char cost[16] = { 0 }, sint[16] = { 0 };  
	//Carrier lookup table for exp(-j*carrPhase)
	if (!cost[0]) {
        for (int i = 0; i < 16; ++i) {
            cost[i] = (char)floor(cos(PI / 8.0 * i) / CSCALE + 0.5);
            sint[i] = (char)floor(-sin(PI / 8.0 * i) / CSCALE + 0.5);
		}
	}

	/* Convert radians to a 16-state phase accumulator used by the LUT.
	   One integer step corresponds to pi/8, so phase indexes are simply
	   obtained by rounding the scaled phase and masking with 0x0F. */
    remCarrPhase = remCarrPhase * 16.0 / DPI;
	/* Carrier phase step */
    const double phaseStep = carrPhaseStep * 16.0 / DPI;

	__m256i Index_mm256, Index_reg1, Index_reg2, cos_reg_mm256, sin_reg_mm256, IF_reg, iBaseband_reg, qBaseband_reg;
	__m256 carrPhase_reg1, carrPhase_reg2;
	__m128i Index_reg, cos_reg_mm128, sin_reg_mm128;
	/* The seemingly odd order compensates for the lane-local packing order of
	   _mm256_packs_epi32. After the later pack/shuffle/permute sequence, the
	   resulting 16 carrier samples still map to chronological sample order
	   0, 1, 2, ..., 15. */
	__m256 stepNumber_reg = _mm256_setr_ps(0.0f, 1.0f, 2.0f, 3.0f, 8.0f, 9.0f, 10.0f, 11.0f);
	const __m256 remCarrPhase_reg = _mm256_set1_ps((float)remCarrPhase);
	const __m256 carrPhaseStep_reg = _mm256_set1_ps((float)phaseStep);
	const __m256i mask4 = _mm256_set1_epi32(15);
	const __m256 sixteenStep = _mm256_set1_ps(16.0f);
	const __m256 fourPhaseStep = _mm256_set1_ps((float)(phaseStep * 4.0));
	
	const __m128i local_cos = _mm_loadu_si128((__m128i *)cost);   
	const __m128i local_sin = _mm_loadu_si128((__m128i *)sint);

    const short *IFData = rawSignal;
	for ( ; IFData <= rawSignal + blksize - 16; IFData+=16,iBasebandSignal+=16,qBasebandSignal+=16) {
		/* ------------------ Inphase and quadrature carriers generation ------------------ */
		/* Generate the carrier phase of 16 consecutive samples. Each AVX register
		   holds 8 phase values; together carrPhase_reg1/carrPhase_reg2 cover the
		   16-sample block that will later be packed back into time order. */
		carrPhase_reg1 = _mm256_fmadd_ps(carrPhaseStep_reg, stepNumber_reg, remCarrPhase_reg);
		carrPhase_reg2 = _mm256_add_ps(carrPhase_reg1, fourPhaseStep);
		
		/* Convert phase to LUT indexes. cvtps rounds to nearest integer, which
		   reduces carrier quantization error compared with truncation. */
		Index_reg1 = _mm256_cvtps_epi32(carrPhase_reg1);   // using _mm256_cvttps_epi32 will result in larger error    
		Index_reg1 = _mm256_and_si256(Index_reg1, mask4);   
		Index_reg2 = _mm256_cvtps_epi32(carrPhase_reg2); 
		Index_reg2 = _mm256_and_si256(Index_reg2, mask4);
		
		/* Pack 8+8 int32 indexes to 16 int16 indexes. The chosen step-number order
		   makes the low/high 128-bit halves already correspond to samples 0..7
		   and 8..15 respectively. */
		Index_mm256 = _mm256_packs_epi32(Index_reg1, Index_reg2);
		
		/* The LUT index range is already 0..15, so each packed int16 index has a
		   zero high byte. Extract the two 128-bit halves and pack them directly
		   to 16 contiguous bytes for the byte-wise LUT shuffle. */
		Index_reg = _mm_packus_epi16(_mm256_castsi256_si128(Index_mm256),
			_mm256_extracti128_si256(Index_mm256, 1));

		/* Look up 16 quantized carrier samples from the 16-entry cosine/sine LUTs.
		   Each byte in Index_reg selects one byte from local_cos/local_sin. */
		cos_reg_mm128 = _mm_shuffle_epi8(local_cos, Index_reg);
		sin_reg_mm128 = _mm_shuffle_epi8(local_sin, Index_reg);

		/* Expand carrier samples to int16 so they can be multiplied directly with
		   the int16 IF samples. Carrier amplitude is quantized by CSCALE = 1/16. */
		cos_reg_mm256 = _mm256_cvtepi8_epi16(cos_reg_mm128);
		sin_reg_mm256 = _mm256_cvtepi8_epi16(sin_reg_mm128);

		/* ----------------- Carrier wiping off from the input IF signals ----------------- */
		/* Load 16 IF samples. The persistent buffers were 32-byte aligned, so the
		   later stores can use aligned writes even though the input may be unaligned. */
		IF_reg = _mm256_loadu_si256((__m256i *)(IFData));
        /* Right-shift IF samples when ADC valid bits occupy the high bits. */
        if (settings.rShiftBits != 0) {
            IF_reg = _mm256_srai_epi16(IF_reg, settings.rShiftBits);
        }
		
		/* Real IF data multiplied by exp(-j*phase):
		     I = x * cos(phase)
		     Q = x * (-sin(phase))
		   The minus sign is already baked into sint[] in mixCarrReal(). */
		iBaseband_reg = _mm256_mullo_epi16(IF_reg, cos_reg_mm256);
		qBaseband_reg = _mm256_mullo_epi16(IF_reg, sin_reg_mm256);

		/* Store the de-rotated baseband block. localCodeGen() and correlator()
		   will reuse these buffers immediately after this stage. */
		_mm256_store_si256((__m256i*)iBasebandSignal, iBaseband_reg);
		_mm256_store_si256((__m256i*)qBasebandSignal, qBaseband_reg);

		/* Advance the synthetic sample index by 16 samples for the next SIMD block. */
		stepNumber_reg = _mm256_add_ps(stepNumber_reg, sixteenStep);
	}	

	/* ----------------- For sampling points outside the SIMD iteration  ----------------- */
	/* Convert the SIMD-loop progress back into the scalar phase accumulator so the
	   tail loop continues from exactly the same quantized carrier state. */
	remCarrPhase = (blksize / 16) * 16 * phaseStep + remCarrPhase;
	
    for (; IFData < rawSignal + blksize; ++IFData, ++iBasebandSignal, ++qBasebandSignal) {
        short sample = IFData[0];
        if (settings.rShiftBits != 0) {
            sample = (short)(sample >> settings.rShiftBits);
        }
		/* Scalar fallback uses the same LUT and scaling as the SIMD path, so the
		   last 0..15 samples remain bit-consistent with the vectorized section. */
		int n = ((int)remCarrPhase) % 16;
        iBasebandSignal[0] = (short)(cost[n] * sample);
        qBasebandSignal[0] = (short)(sint[n] * sample);
        remCarrPhase += phaseStep;
	}
}

/* Mix local carrier for complex data --------------------------------------------
 * Mix local carrier to input signal
 * Args   : short   *rawSignalI        I   Real part of the input complex IF signal
 *          short   *rawSignalQ        I   Imaginary part of the input complex IF signal
 *          double carrPhaseStep       I   carrier sampling interval (cycle: rad*t)
 *          double remCarrPhase        I   Initial carrier phase (rad)
 *          mwSize blksize             I   Number of samples to be generated
 *          short  *iBasebandSignal    O   I component of input signal with carrier wiped off
 *          short  *qBasebandSignal    O   Q component of input signal with carrier wiped off
 * Return : None
 *------------------------------------------------------------------------------*/
void mixCarrComplex(const short* rawSignalI, const short* rawSignalQ, double carrPhaseStep, double remCarrPhase, mwSize blksize,
	short* iBasebandSignal, short* qBasebandSignal)
{
	/* Inphase and quadrature carrier lookup table */
	static char cost[16] = { 0 }, sint[16] = { 0 };

	//Carrier lookup table for exp(-j*carrPhase)
	if (!cost[0]) 	{
		for (int i = 0; i < 16; ++i) {
			cost[i] = (char)floor(cos(PI / 8.0 * i) / CSCALE + 0.5);
			sint[i] = (char)floor(sin(PI / 8.0 * i) / CSCALE + 0.5);
		}
	}

	/* Complex carrier wipeoff uses the same 16-state LUT idea as the real-data
	   version. Here sint[] is positive, and the signs are introduced later by
	   the complex multiply formulas for exp(-j*phase). */
    remCarrPhase = remCarrPhase * 16.0 / DPI;
	/* Carrier phase step */
    const double phaseStep = carrPhaseStep * 16.0 / DPI;

	__m256i Index_mm256, Index_reg1, Index_reg2, cos_reg_mm256, sin_reg_mm256, IF_regI, IF_regQ, \
		sin_regI, cos_regI, sin_regQ, cos_regQ, iBaseband_reg, qBaseband_reg;
	__m256 carrPhase_reg1, carrPhase_reg2;
	__m128i Index_reg, cos_reg_mm128, sin_reg_mm128;
	
	/* Same index-order trick as mixCarrReal(): setr_ps makes the lane contents
	   explicit as [0 1 2 3 | 8 9 10 11], which still compacts to samples
	   0, 1, 2, ..., 15 after the later pack-and-compact path. */
	__m256 stepNumber_reg = _mm256_setr_ps(0.0f, 1.0f, 2.0f, 3.0f, 8.0f, 9.0f, 10.0f, 11.0f);
	const __m256 remCarrPhase_reg = _mm256_set1_ps((float)remCarrPhase);
	const __m256 carrPhaseStep_reg = _mm256_set1_ps((float)phaseStep);
	const __m256i mask4 = _mm256_set1_epi32(15);
	const __m256 sixteenStep = _mm256_set1_ps(16.0f);
	const __m256 fourPhaseStep = _mm256_set1_ps((float)(phaseStep * 4.0));
	const __m128i local_cos = _mm_loadu_si128((__m128i*)cost);
	const __m128i local_sin = _mm_loadu_si128((__m128i*)sint);

    const short *IFDataI = rawSignalI;
    const short *IFDataQ = rawSignalQ;
    for (; IFDataI <= rawSignalI + blksize - 16; IFDataI += 16, IFDataQ += 16, iBasebandSignal += 16, qBasebandSignal += 16) 
	{
		/* ------------------ Inphase and quadrature carriers generation ------------------ */
		/* Build the carrier phase for 16 adjacent complex samples. */
		carrPhase_reg1 = _mm256_fmadd_ps(carrPhaseStep_reg, stepNumber_reg, remCarrPhase_reg);
		carrPhase_reg2 = _mm256_add_ps(carrPhase_reg1, fourPhaseStep);

		/* Convert those phases into 0..15 LUT indexes. */
		Index_reg1 = _mm256_cvtps_epi32(carrPhase_reg1);   // using _mm256_cvttps_epi32 will result in larger error    
		Index_reg1 = _mm256_and_si256(Index_reg1, mask4);
		Index_reg2 = _mm256_cvtps_epi32(carrPhase_reg2);
		Index_reg2 = _mm256_and_si256(Index_reg2, mask4);

		/* Pack and compact the LUT indexes so they can feed a byte-wise pshufb lookup. */
		Index_mm256 = _mm256_packs_epi32(Index_reg1, Index_reg2);

		/* The LUT index range is already 0..15, so each packed int16 index has a
		   zero high byte. Extract the two 128-bit halves and pack them directly
		   to 16 contiguous bytes for the byte-wise LUT shuffle. */
		Index_reg = _mm_packus_epi16(_mm256_castsi256_si128(Index_mm256),
			_mm256_extracti128_si256(Index_mm256, 1));

		/* Fetch the 16 cosine and sine coefficients from the quantized LUTs. */
		cos_reg_mm128 = _mm_shuffle_epi8(local_cos, Index_reg);
		sin_reg_mm128 = _mm_shuffle_epi8(local_sin, Index_reg);

		/* Expand LUT output to int16 for the complex multiply below. */
		cos_reg_mm256 = _mm256_cvtepi8_epi16(cos_reg_mm128);
		sin_reg_mm256 = _mm256_cvtepi8_epi16(sin_reg_mm128);

		/* ----------------- Carrier wiping off from the input IF signals ----------------- */
		/* Load 16 complex samples split as I and Q arrays. */
		IF_regI = _mm256_loadu_si256((__m256i*)IFDataI);
		IF_regQ = _mm256_loadu_si256((__m256i*)IFDataQ);
		
        if (settings.rShiftBits != 0) {
            IF_regI = _mm256_srai_epi16(IF_regI, settings.rShiftBits);
            IF_regQ = _mm256_srai_epi16(IF_regQ, settings.rShiftBits);
		}

		/* Compute the four products needed for:
		     (I + jQ) * (cos - j sin)
		   and then combine them into the real/imag parts of the baseband signal. */
		sin_regI = _mm256_mullo_epi16(IF_regI, sin_reg_mm256);
		cos_regI = _mm256_mullo_epi16(IF_regI, cos_reg_mm256);
		sin_regQ = _mm256_mullo_epi16(IF_regQ, sin_reg_mm256);
		cos_regQ = _mm256_mullo_epi16(IF_regQ, cos_reg_mm256);

		/* Combine the four components to real and imag part of the complex baseband */
		/* iBasebandSignal = real(exp(-j*carrPhase).*rawSignal);
		   qBasebandSignal = imag(exp(-j*carrPhase).*rawSignal); */
		iBaseband_reg = _mm256_add_epi16(cos_regI, sin_regQ);
		qBaseband_reg = _mm256_sub_epi16(cos_regQ, sin_regI);
		

		/* Save the de-rotated complex baseband samples for the code correlator. */
		_mm256_store_si256((__m256i*)iBasebandSignal, iBaseband_reg);
		_mm256_store_si256((__m256i*)qBasebandSignal, qBaseband_reg);

		stepNumber_reg = _mm256_add_ps(stepNumber_reg, sixteenStep);
	}

	/* ----------------- For sampling points outside the SIMD iteration  ----------------- */
	/* Reconstruct the scalar phase state that follows the vectorized section. */
	remCarrPhase = (blksize / 16) * 16 * phaseStep + remCarrPhase;

	for (; IFDataI < rawSignalI + blksize; ++IFDataI, ++IFDataQ, ++iBasebandSignal, ++qBasebandSignal)
	{
        short sampleI = IFDataI[0];
        short sampleQ = IFDataQ[0];
        if (settings.rShiftBits != 0) {
            sampleI = (short)(sampleI >> settings.rShiftBits);
            sampleQ = (short)(sampleQ >> settings.rShiftBits);
        }
		/* Scalar fallback uses the same complex multiply as the SIMD path:
		     Ibb =  I*cos + Q*sin
		     Qbb =  Q*cos - I*sin */
		int n = ((int)remCarrPhase) % 16;
        iBasebandSignal[0] = (short)(cost[n] * sampleI + sint[n] * sampleQ);
        qBasebandSignal[0] = (short)(cost[n] * sampleQ - sint[n] * sampleI);
		remCarrPhase += phaseStep;
	}
}

 /* Generate local code ----------------------------------------------------------
 * Generate local code sequences for both data and pilot branches.
 * Args   : int     *caCode          I   concatenated QPSK code table
 *                                     [dataCode pilotCode]
 *          mwSize codeLen           I   chip number of one branch code table
 *          double remCodePhase      I   initial code phase (chip)
 *          double codePhaseStep     I   code sampling interval (chip)
 *          double earlyLateSpc      I   half of the early-late code correlation spacing (chips)
 *          int    blksize           I   number of samples to be generated 
 *          short  *localCode        O   local code replica outputs for
 *                                     data/pilot Early, Prompt and Late
 * return : None 
 *------------------------------------------------------------------------------*/
void localCodeGen(const int *caCode, mwSize codeLen, double remCodePhase, double codePhaseStep,
	              double earlyLateSpc, mwSize blksize, short **localCode)
{
	/* Initial code phase of early code: 1.0 is due to the first code
	chip added to caCode for early replica generation */
	remCodePhase = remCodePhase - earlyLateSpc + 1.0;
	if (remCodePhase >= codeLen)
		remCodePhase -= floor(remCodePhase / codeLen)*codeLen;
	/* Keep the large integer part of the code phase out of AVX float
	   arithmetic. Long QPSK codes such as GPS L2C can pass million-chip
	   absolute phases; converting that whole value to float loses fractional
	   chip precision. */
	const int codeIndexBase = (int)floor(remCodePhase);
	remCodePhase -= (double)codeIndexBase;
	
	/* Each branch of caCode is laid out as [lastChip, code(1..N), firstChip].
	   Shifting the starting phase by -earlyLateSpc + 1.0 means the early replica can safely
	   address one chip before the nominal prompt position, and the late replica
	   can safely walk one chip beyond the end without needing modulo per sample. */
	short *eCodeD = localCode[0], *pCodeD = localCode[1], * lCodeD = localCode[2];
	short* eCodeP = localCode[3], * pCodeP = localCode[4], * lCodeP = localCode[5];

	__m256i Index_reg1, Index_reg2, EPL_reg;
	__m256 codePhase_reg1, codePhase_reg2;
	/* Here the natural 0..7 order is fine because code indexes are gathered as
	   two 8-lane AVX2 vectors and packed back to 16 x int16 code chips. */
	__m256 stepNumber_reg = _mm256_setr_ps(0.0f, 1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f);
	const __m256 remCodePhase_reg = _mm256_set1_ps((float)remCodePhase);
	const __m256 codePhaseStep_reg = _mm256_set1_ps((float)codePhaseStep);
	const __m256 earlyLateSpc_reg = _mm256_set1_ps((float)earlyLateSpc);
	const __m256 sixteenStep = _mm256_set1_ps(16.0f);
	const __m256 eightCodePhase = _mm256_set1_ps((float)(codePhaseStep * 8.0));
	const __m256i codeIndexBase_reg = _mm256_set1_epi32(codeIndexBase);

	for ( ; eCodeD <= localCode[0] + blksize - 16;
        eCodeD += 16, pCodeD += 16, lCodeD += 16, eCodeP += 16, pCodeP += 16, lCodeP += 16)
	{
		/* ----------------------- Early code generation ----------------------- */
		/* Generate the code phase of 16 successive samples. Because the phase is
		   always non-negative in this routine, truncation toward zero is equivalent
		   to floor(), matching the MATLAB integer indexing behavior. */
		codePhase_reg1 = _mm256_fmadd_ps(codePhaseStep_reg, stepNumber_reg, remCodePhase_reg);
		codePhase_reg2 = _mm256_add_ps(codePhase_reg1, eightCodePhase);
		Index_reg1 = _mm256_cvttps_epi32(codePhase_reg1);
		Index_reg2 = _mm256_cvttps_epi32(codePhase_reg2);
		Index_reg1 = _mm256_add_epi32(Index_reg1, codeIndexBase_reg);
		Index_reg2 = _mm256_add_epi32(Index_reg2, codeIndexBase_reg);
		/* Gather 16 local-code chips directly from the int32 code table,
		   then pack them back to int16 so correlator() can keep using madd_epi16. */
		EPL_reg = gatherCode16(caCode, Index_reg1, Index_reg2);
		_mm256_store_si256((__m256i*)eCodeD, EPL_reg);

		/* Pilot early branch uses the second half of the concatenated code table. */
		EPL_reg = gatherCode16(caCode + codeLen, Index_reg1, Index_reg2);
		_mm256_store_si256((__m256i*)eCodeP, EPL_reg);

		/* ----------------------- Prompt code generation ----------------------- */
		/* Prompt is just Early shifted by +earlyLateSpc chips. */
		codePhase_reg1 = _mm256_add_ps(codePhase_reg1,earlyLateSpc_reg);
		codePhase_reg2 = _mm256_add_ps(codePhase_reg2, earlyLateSpc_reg);
		Index_reg1 = _mm256_cvttps_epi32(codePhase_reg1);
		Index_reg2 = _mm256_cvttps_epi32(codePhase_reg2);
		Index_reg1 = _mm256_add_epi32(Index_reg1, codeIndexBase_reg);
		Index_reg2 = _mm256_add_epi32(Index_reg2, codeIndexBase_reg);
		EPL_reg = gatherCode16(caCode, Index_reg1, Index_reg2);
		_mm256_store_si256((__m256i*)pCodeD, EPL_reg);

		EPL_reg = gatherCode16(caCode + codeLen, Index_reg1, Index_reg2);
		_mm256_store_si256((__m256i*)pCodeP, EPL_reg);
		/* ----------------------- Late code generation ----------------------- */
		/* Late is Prompt shifted once more by +earlyLateSpc, i.e. Early + 2*spacing. */
		codePhase_reg1 = _mm256_add_ps(codePhase_reg1, earlyLateSpc_reg);
		codePhase_reg2 = _mm256_add_ps(codePhase_reg2, earlyLateSpc_reg);
		Index_reg1 = _mm256_cvttps_epi32(codePhase_reg1);
		Index_reg2 = _mm256_cvttps_epi32(codePhase_reg2);
		Index_reg1 = _mm256_add_epi32(Index_reg1, codeIndexBase_reg);
		Index_reg2 = _mm256_add_epi32(Index_reg2, codeIndexBase_reg);
		EPL_reg = gatherCode16(caCode, Index_reg1, Index_reg2);
		_mm256_store_si256((__m256i*)lCodeD, EPL_reg);

		EPL_reg = gatherCode16(caCode + codeLen, Index_reg1, Index_reg2);
		_mm256_store_si256((__m256i*)lCodeP, EPL_reg);

		/* Move to the next batch of 16 code samples. */
		stepNumber_reg = _mm256_add_ps(stepNumber_reg, sixteenStep);
	}

	/* ------------ For sampling points outside the SIMD iteration  ------------ */
	/* Recreate the code phase at the first leftover sample after the SIMD loop. */
	double tempCodePhase = (blksize / 16) * 16 * codePhaseStep + remCodePhase;
	double twoEarlyLateSpc = earlyLateSpc * 2;
	for (; eCodeD < localCode[0] + blksize;
        ++eCodeD, ++pCodeD, ++lCodeD, ++eCodeP, ++pCodeP, ++lCodeP)
	{
		/* Scalar tail keeps the same data/pilot Early/Prompt/Late geometry as the SIMD path. */
		eCodeD[0] = caCode[codeIndexBase + (int)tempCodePhase];
		pCodeD[0] = caCode[codeIndexBase + (int)(tempCodePhase + earlyLateSpc)];
		lCodeD[0] = caCode[codeIndexBase + (int)(tempCodePhase + twoEarlyLateSpc)];
		eCodeP[0] = caCode[codeLen + codeIndexBase + (int)tempCodePhase];
		pCodeP[0] = caCode[codeLen + codeIndexBase + (int)(tempCodePhase + earlyLateSpc)];
		lCodeP[0] = caCode[codeLen + codeIndexBase + (int)(tempCodePhase + twoEarlyLateSpc)];
		tempCodePhase += codePhaseStep;
	}
}

/* Correlating function ------------------------------------------------------------
 * Correlate one branch of local code replicas with the carrier-wiped baseband signal.
 * Args :   short  *iBasebandSignal    I   I component of input signal with carrier wiped off
 *          short  *qBasebandSignal    I   Q component of input signal with carrier wiped off
 *          short  **localCode         I   one three-replica branch {E, P, L}
 *          mwSize blksize             I   number of samples to be generated
 *          double *corrValues         O   outputs of one branch arranged as
 *                                        [I_E Q_E I_P Q_P I_L Q_L]
 * Return : None
 *------------------------------------------------------------------------------*/
void correlator(short *iBasebandSignal, short *qBasebandSignal, short **localCode, mwSize blksize, double *corrValues)
{
	int I_E_sum[8], I_P_sum[8], I_L_sum[8], Q_E_sum[8], Q_P_sum[8], Q_L_sum[8];
	long long I_E = 0, I_P = 0, I_L = 0, Q_E = 0, Q_P = 0, Q_L = 0;
	
	const short* pI = iBasebandSignal;
	const short* pQ = qBasebandSignal;
	short* eCode = localCode[0];
	short* pCode = localCode[1];
	short* lCode = localCode[2];

	/* Each local code sample is stored as int16 +/-1. _mm256_madd_epi16 performs:
	     [a0*b0 + a1*b1, a2*b2 + a3*b3, ...]
	   so one AVX instruction produces eight int32 partial sums from 16 samples. */
	__m256i I_reg, Q_reg, I_E_reg, I_P_reg, I_L_reg, Q_E_reg, Q_P_reg, Q_L_reg, EPL_reg, mul_reg;
	I_E_reg = _mm256_setzero_si256();
	I_P_reg = _mm256_setzero_si256();
	I_L_reg = _mm256_setzero_si256();
	Q_E_reg = _mm256_setzero_si256();
	Q_P_reg = _mm256_setzero_si256();
	Q_L_reg = _mm256_setzero_si256();
	
	for (; eCode <= localCode[0] + blksize - 16; eCode += 16, pCode += 16, lCode += 16, pI += 16, pQ += 16) {
		
		/* --------------- Load IF signal with carrier wiped off --------------- */
		/* These loads are aligned because iBasebandSignal/qBasebandSignal/localCode
		   were allocated with 32-byte alignment and advanced in 16-sample steps. */
		I_reg = _mm256_load_si256((__m256i *)pI);
		Q_reg = _mm256_load_si256((__m256i *)pQ);

		/* --------------- Compute the early correlation value for data channel --------------- */
		EPL_reg = _mm256_load_si256((__m256i *)eCode);
		/* Multiply 16 samples by 16 code chips and horizontally add adjacent pairs,
		   leaving eight int32 partial sums to accumulate in the AVX register. */
		mul_reg = _mm256_madd_epi16(I_reg, EPL_reg);  //Multiplies signed packed 16-bit integer data elements of two vectors. 
		I_E_reg = _mm256_add_epi32(I_E_reg, mul_reg);
		mul_reg = _mm256_madd_epi16(Q_reg, EPL_reg);
		Q_E_reg = _mm256_add_epi32(Q_E_reg, mul_reg);

		/* --------------- Compute the prompt correlation value for data channel --------------- */
		EPL_reg = _mm256_load_si256((__m256i *)pCode);
		mul_reg = _mm256_madd_epi16(I_reg, EPL_reg);
		I_P_reg = _mm256_add_epi32(I_P_reg, mul_reg);
		mul_reg = _mm256_madd_epi16(Q_reg, EPL_reg);
		Q_P_reg = _mm256_add_epi32(Q_P_reg, mul_reg);

		/* --------------- Compute the late correlation value for data channel --------------- */
		EPL_reg = _mm256_load_si256((__m256i *)lCode);
		mul_reg = _mm256_madd_epi16(I_reg, EPL_reg);
		I_L_reg = _mm256_add_epi32(I_L_reg, mul_reg);
		mul_reg = _mm256_madd_epi16(Q_reg, EPL_reg);
		Q_L_reg = _mm256_add_epi32(Q_L_reg, mul_reg);
	}

	/* --------------- Integration for SIMD register elements ---------------
 * The SIMD partial sums are first stored to scalar arrays and then reduced
 * to 64-bit accumulators to reduce the risk of overflow in long integrations.
 *---------------------------------------------------------------------- */
	_mm256_storeu_si256((__m256i*)I_E_sum, I_E_reg);
	_mm256_storeu_si256((__m256i*)Q_E_sum, Q_E_reg);
	_mm256_storeu_si256((__m256i*)I_P_sum, I_P_reg);
	_mm256_storeu_si256((__m256i*)Q_P_sum, Q_P_reg);
	_mm256_storeu_si256((__m256i*)I_L_sum, I_L_reg);
	_mm256_storeu_si256((__m256i*)Q_L_sum, Q_L_reg);
    
	for (int i = 0; i < 8; ++i) {
		I_E += I_E_sum[i]; Q_E += Q_E_sum[i];
		I_P += I_P_sum[i]; Q_P += Q_P_sum[i];
		I_L += I_L_sum[i]; Q_L += Q_L_sum[i];
	}

	/* ----------------- For sampling points outside the SIMD iteration  ----------------- */
	/* Scalar cleanup directly accumulates the remaining 0..15 samples. */
	for (; eCode < localCode[0] + blksize; ++pI, ++pQ, ++eCode, ++pCode, ++lCode) {
		I_E += (long long)pI[0] * (long long)eCode[0];
		Q_E += (long long)pQ[0] * (long long)eCode[0];
		I_P += (long long)pI[0] * (long long)pCode[0];
		Q_P += (long long)pQ[0] * (long long)pCode[0];
		I_L += (long long)pI[0] * (long long)lCode[0];
		Q_L += (long long)pQ[0] * (long long)lCode[0];
	}
	/* Convert the six correlator outputs to double precision and undo the 16x
	   carrier LUT amplitude scaling introduced during carrier wipeoff. */
    corrValues[0] = (double)I_E * CSCALE;
    corrValues[1] = (double)Q_E * CSCALE;
    corrValues[2] = (double)I_P * CSCALE;
    corrValues[3] = (double)Q_P * CSCALE;
    corrValues[4] = (double)I_L * CSCALE;
    corrValues[5] = (double)Q_L * CSCALE;
}

/* checkedMalloc ---------------------------------------------------------------
 * Allocate non-aligned host memory and stop execution if allocation fails.
 * Args   : size_t nbytes             I   number of bytes to allocate
 *          const char* name          I   buffer name for error reporting
 * Return : void*                        allocated host pointer
 *------------------------------------------------------------------------------*/
static void* checkedMalloc(size_t nbytes, const char* name)
{
	void* p = malloc(nbytes);
	if (p == NULL) {
		mexErrMsgIdAndTxt("simdCorr:malloc", "Memory allocation failed for %s.", name);
	}
	return p;
}

/* checkedAlignedMalloc --------------------------------------------------------
 * Allocate 32-byte aligned host memory for SIMD buffers.
 * Args   : size_t nbytes             I   number of bytes to allocate
 *          const char* name          I   buffer name for error reporting
 * Return : void*                        allocated aligned pointer
 *------------------------------------------------------------------------------*/
static void* checkedAlignedMalloc(size_t nbytes, const char* name)
{
	void* p = _mm_malloc(nbytes, 32);
	if (p == NULL) {
		mexErrMsgIdAndTxt("simdCorr:malloc", "Aligned memory allocation failed for %s.", name);
	}
	return p;
}

/* gatherCode16 ----------------------------------------------------------------
 * Gather 16 local-code chips from the int32 code table and repack them
 * into one AVX register of 16 x int16 samples for the existing correlator path.
 * Args   : const int* caCode         I   int32 code table for one PRN
 *          __m256i Index_reg1        I   code indexes for samples 0..7
 *          __m256i Index_reg2        I   code indexes for samples 8..15
 * Return : __m256i                      gathered code chips in chronological order
 *------------------------------------------------------------------------------*/
static __m256i gatherCode16(const int* caCode, __m256i Index_reg1, __m256i Index_reg2)
{
	__m256i EPL_reg1 = _mm256_i32gather_epi32(caCode, Index_reg1, 4);
	__m256i EPL_reg2 = _mm256_i32gather_epi32(caCode, Index_reg2, 4);
	__m256i EPL_reg = _mm256_packs_epi32(EPL_reg1, EPL_reg2);
	return _mm256_permute4x64_epi64(EPL_reg, 0xD8);
}

/* getStructIntField -----------------------------------------------------------
 * Read a scalar numeric field from the MATLAB settings struct.
 * Args   : const mxArray* s          I   settings structure
 *          const char* name          I   field name
 * Return : int                           field value converted to int
 *------------------------------------------------------------------------------*/
static int getStructIntField(const mxArray* s, const char* name)
{
	const mxArray* f = mxGetField(s, 0, name);
	if (f == NULL) {
        mexErrMsgIdAndTxt("simdCorr:field", "Missing field: settings.%s", name);
	}
	
	if (!mxIsDouble(f) || mxIsComplex(f) || mxGetNumberOfElements(f) != 1) {
        mexErrMsgIdAndTxt("simdCorr:type", "settings.%s must be a real scalar.", name);
	}

	return (int)mxGetScalar(f);
}

/* getStructDoubleField --------------------------------------------------------
 * Read a scalar numeric field from the MATLAB settings struct.
 * Args   : const mxArray* s          I   settings structure
 *          const char* name          I   field name
 * Return : double                        field value converted to double
 *------------------------------------------------------------------------------*/
static double getStructDoubleField(const mxArray* s, const char* name)
{
	const mxArray* f = mxGetField(s, 0, name);
	if (f == NULL) {
        mexErrMsgIdAndTxt("simdCorr:field", "Missing field: settings.%s", name);
	}

	if (!mxIsDouble(f) || mxIsComplex(f) || mxGetNumberOfElements(f) != 1) {
        mexErrMsgIdAndTxt("simdCorr:type", "settings.%s must be a real scalar.", name);
	}

	return mxGetScalar(f);
}

/* parseSettings ---------------------------------------------------------------
 * Extract the required settings fields from the MATLAB settings structure.
 * Args   : const mxArray* s          I   settings structure
 * Return : Settings                      parsed settings structure
 *------------------------------------------------------------------------------*/
static Settings parseSettings(const mxArray* s)
{
	Settings cfg;
	cfg.rShiftBits = 0;

	if (!mxIsStruct(s)) {
		mexErrMsgIdAndTxt("simdCorr:type", "The input settings must be a struct.");
	}

	cfg.trkMode = getStructIntField(s, "trkMode");
	cfg.fileType = getStructIntField(s, "fileType");
	/* dllCorrelatorSpacing is fractional in this receiver (default 0.5 chip),
	   so reading it as int collapses Early/Prompt/Late into the same code. */
	cfg.earlyLateSpc = getStructDoubleField(s, "dllCorrelatorSpacing");
	cfg.correlatorType = getStructIntField(s, "correlatorType");
	if (cfg.correlatorType == 1) // SIMD correlator
	{
		cfg.rShiftBits = getStructIntField(s, "rShiftBits");
	}

	return cfg;
}

/* checkInputs ----------------------------------------------------------------
 * Validate the number, type and basic form of the MEX input arguments.
 * Args   : int nlhs                  I   number of output arguments
 *          int nrhs                  I   number of input arguments
 *          const mxArray* prhs[]     I   input argument list
 * Return : None
 *------------------------------------------------------------------------------*/
static void checkInputs(int nlhs, int nrhs, const mxArray* prhs[])
{
	/* ----------------- Get input arrays -------------------
	* Input order:
 	*   prhs[0]  - settings
 	*   prhs[1]  - rawSignal
 	*   prhs[2]  - caCodeTable
 	*   prhs[3]  - remCarrPhase
 	*   prhs[4]  - carrPhaseStep
 	*   prhs[5]  - remCodePhase
 	*   prhs[6]  - codePhaseStep
 	*   prhs[7]  - startIndex
 	*   prhs[8]  - chSampleSize
 	*   prhs[9] - isDataRead
    *--------------------------------------------------------- */

	if (nrhs != 10) {
		mexErrMsgIdAndTxt("simdCorr:nrhs", "10 input arguments are required.");
	}
	if (nlhs > 1) {
		mexErrMsgIdAndTxt("simdCorr:nlhs", "Too many output arguments.");
	}

	if (!mxIsStruct(prhs[0])) {
		mexErrMsgIdAndTxt("simdCorr:type", "The 1st input must be a settings struct.");
	}
	if (!mxIsInt16(prhs[1]) || mxIsComplex(prhs[1])) {
		mexErrMsgIdAndTxt("simdCorr:type", "rawSignal must be a real int16 array.");
	}
	if (!mxIsInt32(prhs[2]) || mxIsComplex(prhs[2])) {
		mexErrMsgIdAndTxt("simdCorr:type", "caCodeTable must be a real int32 array.");
	}
	if ((mxGetM(prhs[2]) & 1) != 0) {
		mexErrMsgIdAndTxt("simdCorr:type",
			"caCodeTable must be arranged as [dataCodeTable; pilotCodeTable], so its row count must be even.");
	}
	for (int k = 3; k <= 6; ++k) {
		if (!mxIsDouble(prhs[k]) || mxIsComplex(prhs[k])) {
			mexErrMsgIdAndTxt("simdCorr:type", "Phase-related inputs must be real double arrays.");
		}
	}
	if (!mxIsInt32(prhs[7]) || mxIsComplex(prhs[7])) {
		mexErrMsgIdAndTxt("simdCorr:type", "startIndex must be an int32 array.");
	}
	if (!mxIsInt32(prhs[8]) || mxIsComplex(prhs[8])) {
		mexErrMsgIdAndTxt("simdCorr:type", "chSampleSize must be an int32 array.");
	}
	if (!mxIsDouble(prhs[9]) || mxIsComplex(prhs[9]) || mxGetNumberOfElements(prhs[9]) != 1) {
		mexErrMsgIdAndTxt("simdCorr:type", "isDataRead must be a real scalar.");
	}
}

/* cleanup --------------------------------------------------------------------
 * Release all persistent host buffers when the MEX function is cleared.
 * Return : None
 *------------------------------------------------------------------------------*/
void cleanup(void)
{
	mexPrintf("   MEX-file is terminating, destroying allocated memory ...\n");
	/* --------------- Free memory --------------- */
	if (localCode != NULL) {
		for (int i = 0; i < 6; ++i) {
			if (localCode[i] != NULL) {
				_mm_free(localCode[i]);
				localCode[i] = NULL;
			}
		}
		free(localCode);
		localCode = NULL;
	}

	if (iBasebandSignal != NULL) { _mm_free(iBasebandSignal); iBasebandSignal = NULL; }
	if (qBasebandSignal != NULL) { _mm_free(qBasebandSignal); qBasebandSignal = NULL; }
	if (rawSignalI != NULL) { _mm_free(rawSignalI); rawSignalI = NULL; }
	if (rawSignalQ != NULL) { _mm_free(rawSignalQ); rawSignalQ = NULL; }
	initialized = 0;
}


