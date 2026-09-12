/*================================================================================
 * Filename: corrGpuParallelFusedBPSK.cu
 * Description: CUDA MEX implementation of GPU-assisted GNSS tracking correlator.
 *
 * Authors: Yafeng Li (School of Automation, Beijing Information Science and Technology University)
 * Time: Feb, 20, 2026
 *
 * Notes:
 *   - Input/Output interface is kept compatible with the original implementation.
 *   - fileType = 1: real IF samples [I0 I1 I2 ...]
 *   - fileType = 2: complex interleaved IF samples [I0 Q0 I1 Q1 ...]
 *   - Output layout is 6 x channelCnt:
 *         [IE; QE; IP; QP; IL; QL]
 *================================================================================*/

#include "mex.h"
#include <math.h>
#include <string>
#include <complex>
#include <cuda.h>
#include <cuda_runtime.h>

// ============================ Configuration ====================================
#define THREADS_PER_BLOCK  256
#define CORR_NUMBER 3    // Early / Prompt / Late    
#define ACCUM_N 256

/* Settings struct -------------------------------------------------------------
 * trkMode        : 0 - channel-serial tracking; 1 - channel-parallel tracking
 * fileType       : 1 - real samples; 2 - complex samples
 * earlyLateSpc   : half of the early - late code correlation spacing (chips)
 *------------------------------------------------------------------------------*/
typedef struct {
    int trkMode;
    int fileType;
    double earlyLateSpc;
} Settings;

// ============================ CUDA helpers =====================================
#define CUDA_CHECK(call)                                                                          \
    do {                                                                                          \
        cudaError_t err__ = (call);                                                               \
        if (err__ != cudaSuccess) {                                                               \
            mexErrMsgIdAndTxt("gpuCorr:CUDA", "%s failed at %s:%d: %s",                           \
                              #call, __FILE__, __LINE__, cudaGetErrorString(err__));              \
        }                                                                                         \
    } while (0)

#define CUDA_KERNEL_CHECK()                                                                       \
    do {                                                                                          \
        cudaError_t err__ = cudaPeekAtLastError();                                                \
        if (err__ != cudaSuccess) {                                                               \
            mexErrMsgIdAndTxt("gpuCorr:Kernel", "Kernel launch failed at %s:%d: %s",              \
                              __FILE__, __LINE__, cudaGetErrorString(err__));                     \
        }                                                                                         \
    } while (0)

/* ========== GPU internal data types for complex numbers ========== */
struct GPU_Complex
{
	float r;
	float i;
	__host__ __device__ GPU_Complex() : r(0.0f), i(0.0f) {};

	// Constructor for a complex number with real amd imag parts of float
    __host__ __device__ GPU_Complex(float realPart, float imagPart) : r(realPart), i(imagPart) {}

	// Magnitude of a complex number
	__host__ __device__ float magnitude2() { return r * r + i * i; }

	__device__ GPU_Complex operator*(const GPU_Complex& a) {
		return GPU_Complex(__fmul_rn(r, a.r) - __fmul_rn(i, a.i), __fmul_rn(i, a.r) + __fmul_rn(r, a.i));
	}

	__device__ GPU_Complex operator*(const float& a) {
		return GPU_Complex(__fmul_rn(r, a), __fmul_rn(i, a));
	}

	__host__ __device__ GPU_Complex operator+(const GPU_Complex& a) {
		return GPU_Complex(r + a.r, i + a.i);
	}

	__host__ __device__ void operator+=(const GPU_Complex& a) {
		r += a.r;
		i += a.i;
	}

    // Accumulate a * b into *this, where b is complex.
	__device__ void multiply_acc(const GPU_Complex& a, const GPU_Complex& b) {
		//real part
		r = __fmaf_rn(a.r, b.r, r);
		r = __fmaf_rn(-a.i, b.i, r);
		//imag part
		i = __fmaf_rn(a.i, b.r, i);
		i = __fmaf_rn(a.r, b.i, i);
	}

    // Accumulate a * b into *this, where b is real.
	__device__ void multiply_acc(const GPU_Complex& a, const float& b) {
		//real part
		r = __fmaf_rn(a.r, b, r);
		//imag part
		i = __fmaf_rn(a.i, b, i);
	}
};

/* Function declaration ---------------------------------------------------------
 * fusedCorrReal()    : fused carrier wipeoff + Early/Prompt/Late correlation for real IF data
 * fusedCorrComplex() : fused carrier wipeoff + Early/Prompt/Late correlation for complex IF data
 * cleanup()          : release persistent GPU / host memory when MEX is cleared
 * checkInputs()   : validate the MEX input argument types and counts
 * parseSettings()    : extract required settings from MATLAB settings struct
 * getStructIntField(): read a scalar integer-like field from settings struct
 * getStructDoubleField(): read a scalar numeric field from settings struct
 *------------------------------------------------------------------------------*/
__global__ void fusedCorrReal(GPU_Complex*, const short*, const signed char* __restrict__, float, int, float, float, int, float, float);
__global__ void fusedCorrComplex(GPU_Complex*, const short*, const signed char* __restrict__, float, int, float, float, int, float, float);
void cleanup(void);

/* checkInputs -------------------------------------------------------------
 * Validate the number, type and basic form of the MEX input arguments.
 * Args   : int nlhs                  I   number of output arguments
 *          int nrhs                  I   number of input arguments
 *          const mxArray* prhs[]     I   input argument list
 * Return : None
 *------------------------------------------------------------------------------*/
static void checkInputs(int nrhs, int nlhs, const mxArray* prhs[]);
static int getStructIntField(const mxArray* s, const char* name);
static double getStructDoubleField(const mxArray* s, const char* name);
static Settings parseSettings(const mxArray* s);

/* Persistent state -------------------------------------------------------------
 * d_pRawSignal        : device buffer for input IF signal block (int16 samples)
 * d_ppCaCode          : per-channel device pointers to local code tables
 * d_corrValues        : device pointer mapped to host correlation buffer
 * h_corrValues        : host mapped correlation buffer (Early/Prompt/Late)
 * streams             : one CUDA stream for each tracking channel
 * channelCnt          : number of tracking channels
 * settings            : cached receiver settings parsed from MATLAB struct
 *------------------------------------------------------------------------------*/
static short * d_pRawSignal;
static signed char ** d_ppCaCode;
static GPU_Complex * d_corrValues;
static std::complex<float>* h_corrValues;
static cudaStream_t* streams;
static size_t channelCnt;
static Settings settings;

/* The gateway function --------------------------------------------------------
 * Args   : prhs[0]   Settings  settings       I   Receiver settings structure
 *          prhs[1]   int16*    rawSignal      I   Input IF signal samples
 *          prhs[2]   int8*     caCodeTable    I   Local code table for all channels
 *          prhs[3]   double*   remCarrPhase   I   Residual carrier phase for each channel (rad)
 *          prhs[4]   double*   carrPhaseStep  I   Carrier phase step for each sample
 *          prhs[5]   double*   remCodePhase   I   Residual code phase for each channel (chips)
 *          prhs[6]   double*   codePhaseStep  I   Code phase step for each sample (chips/sample)
 *          prhs[7]   int*      startIndex     I   Start sample index for each channel
 *          prhs[8]   int*      chSampSize     I   Number of samples to be processed for each channel
 *          prhs[9]   int       isDataRead     I   Data-read flag (1: read new data, 0: reuse previous data)
 *          plhs[0]   double*   corrValues     O   Correlator outputs for all channels
 * Return : plhs[0]   corrValues
 *          Output layout is 6 x channelCnt: [IE; QE; IP; QP; IL; QL]
 *------------------------------------------------------------------------------*/
void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
	/* --------------- Declare all variables --------------- */
	short* rawSignal;
	signed char* caCodeTable;
	int* startIndex, * chSampSize;
	double* remCarrPhase, * carrPhaseStep, * remCodePhase, * codePhaseStep;
	double* corrValues;
	static size_t codeLen;
	static int isDataRead;

	static int rawSignalLen;
	static int initialized = 0;
	/* ================ Allocate host/device memory and initialize variables ================= */
	if (!initialized)
	{
		checkInputs(nrhs, nlhs, prhs);
		settings = parseSettings(prhs[0]);

        /* ------ Get array dimensions ------
         * codeLen    : number of rows in caCodeTable
         * channelCnt : number of tracking channels
         * fileType   : 1 - real IF; 2 - complex interleaved IF
         *------------------------------------------------ */
		//Size of the PRN code chips
		codeLen = mxGetM(prhs[2]);
		//Channel count
		channelCnt = mxGetN(prhs[2]);
		// file type: real or complex
		if (settings.fileType != 1 && settings.fileType != 2) {
			mexErrMsgIdAndTxt("gpuCorr:fileType", "fileType must be 1 (real) or 2 (complex interleaved).");
		}

        /* --------- Copy PRN codes to per-channel device buffers ---------
         * Each channel owns one device-side local code table.
         * The carrier wipeoff result is no longer staged in a separate
         * baseband buffer; it is consumed directly inside the fused kernel.
         *------------------------------------------------------------------------ */
		// PRN codes memory
		caCodeTable = (signed char*)mxGetData(prhs[2]);
		d_ppCaCode = (signed char**)malloc(sizeof(signed char*) * channelCnt);
		
		if (!d_ppCaCode) {
			mexErrMsgIdAndTxt("gpuCorr:HostAlloc", "Failed to allocate host-side bookkeeping memory.");
		}
		
		mwIndex dims[2] = { 0, 0 };
		for (int chInd = 0; chInd < channelCnt; chInd++)
		{
			dims[1] = chInd;
			size_t index = mxCalcSingleSubscript(prhs[2], (size_t)2, dims);
			signed char* caCode = &caCodeTable[index];
			//malloc device global memory for caCodeTable
			CUDA_CHECK(cudaMalloc(&d_ppCaCode[chInd], sizeof(signed char) * codeLen));
			CUDA_CHECK(cudaMemcpy(d_ppCaCode[chInd], caCode, codeLen * sizeof(signed char), cudaMemcpyHostToDevice));
		}

		/* --------- For input IF signal --------- */
		// malloc device global memory for RawSignal
		rawSignalLen = (int)mxGetNumberOfElements(prhs[1]);
		CUDA_CHECK(cudaMalloc(&d_pRawSignal, rawSignalLen * sizeof(short)));

        /* --------- Correlation values ---------
         * Host mapped pinned memory is used so the GPU can write correlation
         * results directly to a host-visible buffer without an extra copy.
         *---------------------------------------- */
		//malloc device host mapped page-locked memory for corrValues
		CUDA_CHECK(cudaHostAlloc(&h_corrValues, CORR_NUMBER*channelCnt*sizeof(std::complex<float>), cudaHostAllocMapped));
		// Get device pointer from host memory
		CUDA_CHECK(cudaHostGetDevicePointer(&d_corrValues, h_corrValues, 0));

		/* ------- Alloctate CUDA streams ------- */
		streams = (cudaStream_t*)malloc(channelCnt * sizeof(cudaStream_t));
		if (!streams) {
			mexErrMsgIdAndTxt("gpuCorr:", "Failed to allocate streams memory.");
		}
		for (int chInd = 0; chInd < channelCnt; chInd++){
			CUDA_CHECK(cudaStreamCreate(&streams[chInd]));
		}

		/* ------ Free host/device memory ------- */
		mexAtExit(cleanup);
		initialized = 1;
		mexPrintf("   CUDA initialized ...\n");
		mexPrintf("   Channel-parallel tracking with the GPU-accelerated correlator developed by Yafeng Li.\n");
	} // for initialization

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
     *   prhs[8]  - chSampSize
     *   prhs[9]  - isDataRead
     *--------------------------------------------------------- */
	remCarrPhase = mxGetPr(prhs[3]);
	carrPhaseStep = mxGetPr(prhs[4]);
	remCodePhase = mxGetPr(prhs[5]);
	codePhaseStep = mxGetPr(prhs[6]);
	startIndex = (int*)mxGetData(prhs[7]);
	chSampSize = (int*)mxGetData(prhs[8]);
	isDataRead = (int)mxGetScalar(prhs[9]);

	/* -------- Output allocation --------- */
	plhs[0] = mxCreateDoubleMatrix((mwSize)6, (mwSize)channelCnt, mxREAL);
	corrValues = (double*)mxGetData(plhs[0]);

    /* ========================= GPU correlator implementation =========================
     * 1) Upload a new IF block when requested
     * 2) For each channel:
     *      - run one fused kernel that generates the carrier replica,
     *        forms the local code for Early/Prompt/Late, and accumulates
     *        the three correlations inside one block
     * 3) Synchronize all streams
     * 4) Pack results back to MATLAB output matrix
     *=============================================================================== */
	/* Upload IF data when needed */
	if (isDataRead == 1) 	{
		rawSignal = (short*)mxGetData(prhs[1]);
		CUDA_CHECK(cudaMemcpy(d_pRawSignal, rawSignal, rawSignalLen * sizeof(short), cudaMemcpyHostToDevice));
	}

	/* --------- Per-channel GPU processing --------- */
	for (int chInd = 0; chInd < channelCnt; chInd++)
	{
		float initCodePhase = 1.0f - (float)settings.earlyLateSpc + (float)remCodePhase[chInd];
		// 1) Run one fused block per channel.
		if (settings.fileType == 1)
		{ // for real data
			fusedCorrReal<<<1, THREADS_PER_BLOCK, 0, streams[chInd]>>>(
				d_corrValues + CORR_NUMBER * chInd,
				d_pRawSignal + startIndex[chInd] - 1,
				d_ppCaCode[chInd],
				(float)settings.earlyLateSpc,
				(int)codeLen,
				(float)codePhaseStep[chInd],
				initCodePhase,
				chSampSize[chInd],
				(float)remCarrPhase[chInd],
				(float)carrPhaseStep[chInd]);
		}
		else if (settings.fileType == 2)
		{ // for complex data
			fusedCorrComplex<<<1, THREADS_PER_BLOCK, 0, streams[chInd]>>>(
				d_corrValues + CORR_NUMBER * chInd,
				d_pRawSignal + (startIndex[chInd] - 1) * 2,
				d_ppCaCode[chInd],
				(float)settings.earlyLateSpc,
				(int)codeLen,
				(float)codePhaseStep[chInd],
				initCodePhase,
				chSampSize[chInd],
				(float)remCarrPhase[chInd],
				(float)carrPhaseStep[chInd]);
		}
		CUDA_KERNEL_CHECK();
	}

	/* --------- Wait for all streams to finish before reading mapped host memory --------- */
	CUDA_CHECK(cudaDeviceSynchronize());
	// -------------------- Pack output back to MATLAB ---------------------------
	// Output order per channel: [IE, QE, IP, QP, IL, QL]^T.
	for (int chInd = 0; chInd < channelCnt; chInd++) {
		int correIndex = CORR_NUMBER * chInd;
		for (int ind = 0; ind < CORR_NUMBER; ind++) {
			corrValues[6 * chInd + ind * 2] = (double)h_corrValues[correIndex + ind].real();
			corrValues[6 * chInd + ind * 2 + 1] = (double)h_corrValues[correIndex + ind].imag();}
	}
	// Display the last error for debug
	//CUDA_CHECK(cudaGetLastError());
}

// ============================ CUDA kernels =====================================
/*
 * fusedCorrReal
 * -------------
 * Fused carrier wipeoff and Early / Prompt / Late correlation for real IF data.
 * One block handles one tracking channel. Each thread generates the local carrier
 * and the three local code branches for its strided samples, then a shared-memory
 * reduction combines the per-thread partial sums into [E, P, L].
 */
__global__ void fusedCorrReal(GPU_Complex* d_corrValues, const short* d_pRawSignal, const signed char* __restrict__ d_caCodeTable,
	float earlyLateSpc, int codeLen, float codePhaseStep, float initCodePhase, int blksize,
	float remCarrPhase, float carrPhaseStep)
{
	__shared__ GPU_Complex accumResult[CORR_NUMBER][ACCUM_N];
	float initPhase[CORR_NUMBER];
	#pragma unroll
	for (int vec = 0; vec < CORR_NUMBER; ++vec) {
		initPhase[vec] = initCodePhase + earlyLateSpc * vec;
	}

	for (int iAccum = threadIdx.x; iAccum < ACCUM_N; iAccum += blockDim.x)
	{
		GPU_Complex sumIQ[CORR_NUMBER] = { GPU_Complex(0.0f, 0.0f), GPU_Complex(0.0f, 0.0f), GPU_Complex(0.0f, 0.0f) };

		for (int pos = iAccum; pos < blksize; pos += ACCUM_N) {
			float sinVal, cosVal;
			__sincosf(remCarrPhase + pos * carrPhaseStep, &sinVal, &cosVal);
			GPU_Complex baseband = GPU_Complex(cosVal, -sinVal) * (float)d_pRawSignal[pos];

			#pragma unroll
			for (int vec = 0; vec < CORR_NUMBER; ++vec) {
				int chipIndex = __float2int_rd(fmodf(codePhaseStep * __int2float_rd(pos) + initPhase[vec], codeLen));
				sumIQ[vec].multiply_acc(baseband, (float)__ldg(&d_caCodeTable[chipIndex]));
			}
		}

		#pragma unroll
		for (int vec = 0; vec < CORR_NUMBER; ++vec) {
			accumResult[vec][iAccum] = sumIQ[vec];
		}
	}

	// Tree reduction in shared memory. ACCUM_N must be a power of two.
	for (int stride = ACCUM_N / 2; stride > 0; stride >>= 1) {
		__syncthreads();
		for (int iAccum = threadIdx.x; iAccum < stride; iAccum += blockDim.x) {
			#pragma unroll
			for (int vec = 0; vec < CORR_NUMBER; ++vec) {
				accumResult[vec][iAccum] += accumResult[vec][stride + iAccum];
			}
		}
	}

	if (threadIdx.x == 0) {
		#pragma unroll
		for (int vec = 0; vec < CORR_NUMBER; ++vec) {
			d_corrValues[vec] = accumResult[vec][0];
		}
	}
}

/*
 * fusedCorrComplex
 * ----------------
 * Fused carrier wipeoff and Early / Prompt / Late correlation for complex IF data.
 * The input is interleaved as [I0, Q0, I1, Q1, ...]. The carrier and three code
 * branches are generated inside the same block, then reduced to one [E, P, L]
 * output tuple for the channel. Complex raw samples are loaded as short2 so I/Q
 * are fetched together with one vector load.
 */
__global__ void fusedCorrComplex(GPU_Complex* d_corrValues, const short* d_pRawSignal, const signed char* __restrict__ d_caCodeTable,
	float earlyLateSpc, int codeLen, float codePhaseStep, float initCodePhase, int blksize,
	float remCarrPhase, float carrPhaseStep)
{
	__shared__ GPU_Complex accumResult[CORR_NUMBER][ACCUM_N];
	const short2* d_pRawSignalIQ = reinterpret_cast<const short2*>(d_pRawSignal);
	float initPhase[CORR_NUMBER];
	#pragma unroll
	for (int vec = 0; vec < CORR_NUMBER; ++vec) {
		initPhase[vec] = initCodePhase + earlyLateSpc * vec;
	}

	for (int iAccum = threadIdx.x; iAccum < ACCUM_N; iAccum += blockDim.x)
	{
		GPU_Complex sumIQ[CORR_NUMBER] = { GPU_Complex(0.0f, 0.0f), GPU_Complex(0.0f, 0.0f), GPU_Complex(0.0f, 0.0f) };

		for (int pos = iAccum; pos < blksize; pos += ACCUM_N) {
			float sinVal, cosVal;
			__sincosf(remCarrPhase + pos * carrPhaseStep, &sinVal, &cosVal);
			short2 rawIQ = d_pRawSignalIQ[pos];
			GPU_Complex baseband = GPU_Complex(cosVal, -sinVal) *
				GPU_Complex((float)rawIQ.x, (float)rawIQ.y);

			#pragma unroll
			for (int vec = 0; vec < CORR_NUMBER; ++vec) {
				int chipIndex = __float2int_rd(fmodf(codePhaseStep * __int2float_rd(pos) + initPhase[vec], codeLen));
				sumIQ[vec].multiply_acc(baseband, (float)__ldg(&d_caCodeTable[chipIndex]));
			}
		}

		#pragma unroll
		for (int vec = 0; vec < CORR_NUMBER; ++vec) {
			accumResult[vec][iAccum] = sumIQ[vec];
		}
	}

	// Tree reduction in shared memory. ACCUM_N must be a power of two.
	for (int stride = ACCUM_N / 2; stride > 0; stride >>= 1) {
		__syncthreads();
		for (int iAccum = threadIdx.x; iAccum < stride; iAccum += blockDim.x) {
			#pragma unroll
			for (int vec = 0; vec < CORR_NUMBER; ++vec) {
				accumResult[vec][iAccum] += accumResult[vec][stride + iAccum];
			}
		}
	}

	if (threadIdx.x == 0) {
		#pragma unroll
		for (int vec = 0; vec < CORR_NUMBER; ++vec) {
			d_corrValues[vec] = accumResult[vec][0];
		}
	}
}

// ============================ MATLAB validation ================================
static int getStructIntField(const mxArray* s, const char* name)
{
	const mxArray* f = mxGetField(s, 0, name);
	if (f == NULL) {
		mexErrMsgIdAndTxt("gpuCorr:field", "Missing field: settings.%s", name);
	}
	if (!mxIsDouble(f) || mxIsComplex(f) || mxGetNumberOfElements(f) != 1) {
		mexErrMsgIdAndTxt("gpuCorr:type", "settings.%s must be a real scalar.", name);
	}
	return (int)mxGetScalar(f);
}

static double getStructDoubleField(const mxArray* s, const char* name)
{
	const mxArray* f = mxGetField(s, 0, name);
	if (f == NULL) {
		mexErrMsgIdAndTxt("gpuCorr:field", "Missing field: settings.%s", name);
	}
	if (!mxIsDouble(f) || mxIsComplex(f) || mxGetNumberOfElements(f) != 1) {
		mexErrMsgIdAndTxt("gpuCorr:type", "settings.%s must be a real scalar.", name);
	}
	return mxGetScalar(f);
}

static Settings parseSettings(const mxArray* s)
{
	Settings cfg;

    if (!mxIsStruct(s)) {
        mexErrMsgIdAndTxt("gpuCorr:type", "The input settings must be a struct.");
    }

    cfg.trkMode  = getStructIntField(s, "trkMode");
    cfg.fileType = getStructIntField(s, "fileType");
    cfg.earlyLateSpc = getStructDoubleField(s, "dllCorrelatorSpacing");

    return cfg;
}
static void checkInputs(int nrhs, int nlhs, const mxArray* prhs[])
{
	if (nlhs > 1) {
		mexErrMsgIdAndTxt("gpuCorr:Output", "One output is expected.");
	}

	if (nrhs != 10) {
		mexErrMsgIdAndTxt("gpuCorr:InputCount",
			"Expected 10 inputs: settings, rawSignal, caCodeTable, remCarrPhase, carrPhaseStep, remCodePhase, codePhaseStep, startIndex, chSampSize, isDataRead.");
	}

	if (!mxIsStruct(prhs[0])) {
		mexErrMsgIdAndTxt("gpuCorr:settings", "settings must be a struct.");
	}
	if (!mxIsInt16(prhs[1]) || mxIsComplex(prhs[1])) {
		mexErrMsgIdAndTxt("gpuCorr:rawSignal", "rawSignal must be a real int16 array.");
	}
	if (!mxIsInt8(prhs[2]) || mxIsComplex(prhs[2])) {
		mexErrMsgIdAndTxt("gpuCorr:caCode", "caCodeTable must be a real int8 matrix.");
	}

	for (int k = 3; k <= 6; ++k) {
		if (!mxIsDouble(prhs[k]) || mxIsComplex(prhs[k])) {
			mexErrMsgIdAndTxt("gpuCorr:Type", "Inputs 4-7 must be real double arrays.");
		}
	}

	if (!mxIsInt32(prhs[7]) || mxIsComplex(prhs[7])) {
		mexErrMsgIdAndTxt("gpuCorr:startIndex", "startIndex must be an int32 array.");
	}
	if (!mxIsInt32(prhs[8]) || mxIsComplex(prhs[8])) {
		mexErrMsgIdAndTxt("gpuCorr:chSampSize", "chSampSize must be an int32 array.");
	}
	if (!mxIsDouble(prhs[9]) || mxIsComplex(prhs[9]) || mxGetNumberOfElements(prhs[9]) != 1) {
		mexErrMsgIdAndTxt("gpuCorr:isDataRead", "isDataRead must be a real double scalar.");
	}
}

/* cleanup --------------------------------------------------------------------
 * Release all persistent device and host-side resources when the MEX file is
 * cleared by MATLAB.
 * Return : None
 *------------------------------------------------------------------------------*/
void cleanup(void)
{
	mexPrintf("   CUDA is terminating, destroying allocated memory ...\n");

	/* --------- Free device-side memory --------- */
	for (size_t chInd = 0; chInd < channelCnt; ++chInd) {
		if (d_ppCaCode != NULL && d_ppCaCode[chInd] != NULL) {
			CUDA_CHECK(cudaFree(d_ppCaCode[chInd]));
		}
		if (streams != NULL) {
			CUDA_CHECK(cudaStreamDestroy(streams[chInd]));
		}
	}

	if (d_pRawSignal != NULL) {
		CUDA_CHECK(cudaFree(d_pRawSignal));
		d_pRawSignal = NULL;
	}

	/* --------- Free host-side memory --------- */
	if (h_corrValues != NULL) {
		CUDA_CHECK(cudaFreeHost(h_corrValues));
		h_corrValues = NULL;
		d_corrValues = NULL;
	}

	if (streams != NULL)  { free(streams); streams = NULL; }
	if (d_ppCaCode != NULL)  { free(d_ppCaCode); d_ppCaCode = NULL; }
}


