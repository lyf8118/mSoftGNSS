/*================================================================================
 * Filename: corrGPUSerialQMBOC.cu
 * Description: GPU-assisted single-channel QMBOC tracking correlator for
 *              trkChannelsSerial.m.
 *
 * Authors: Yafeng Li (School of Automation, Beijing Information Science and Technology University)
 * Time: Apr, 21, 2026
 *
 * The core processing stages are intentionally kept aligned with the
 * GPU channel-parallel correlator:
 *   1) mixCarrReal() / mixCarrComplex()
 *   2) correlator()
 *
 * The difference is only in the MEX gateway framework: this file processes
 * one channel / one code period per call, matching trkChannelsSerial.m.
 *===============================================================================*/

#include "mex.h"
#include <math.h>
#include <complex>
#include <cuda.h>
#include <cuda_runtime.h>

// ============================ Configuration ====================================
#define THREADS_PER_BLOCK  256
#define CORR_NUMBER 9    // Data + pilot BOC(1,1) + pilot BOC(6,1), each Early/Prompt/Late
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
    __host__ __device__ GPU_Complex() : r(0.0f), i(0.0f) {}

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
 * mixCarrReal()      : mix local carrier with real IF data on GPU
 * mixCarrComplex()   : mix local carrier with complex interleaved IF data on GPU
 * correlator()       : compute data/pilot Early / Prompt / Late correlations on GPU
 * cleanup()          : release persistent GPU / host memory when MEX is cleared
 * checkInputs()      : validate the MEX input argument types and counts
 * parseSettings()    : extract required settings from MATLAB settings struct
 * getStructIntField(): read a scalar integer-like field from settings struct
 * getStructDoubleField(): read a scalar numeric field from settings struct
 *------------------------------------------------------------------------------*/
__global__ void mixCarrReal(const short*, GPU_Complex*, int, float, float);
__global__ void mixCarrComplex(const short*, GPU_Complex*, int, float, float);
__global__ void correlator(GPU_Complex*, GPU_Complex*, const signed char* __restrict__, float, int, float, int, int);
void cleanup(void);

/* checkInputs -------------------------------------------------------------
 * Validate the number, type and basic form of the MEX input arguments.
 * Args   : int nrhs                  I   number of input arguments
 *          int nlhs                  I   number of output arguments
 *          const mxArray* prhs[]     I   input argument list
 * Return : None
 *------------------------------------------------------------------------------*/
static void checkInputs(int nrhs, int nlhs, const mxArray* prhs[]);
static int getStructIntField(const mxArray* s, const char* name);
static double getStructDoubleField(const mxArray* s, const char* name);
static Settings parseSettings(const mxArray* s);

/* Persistent state -------------------------------------------------------------
 * d_pRawSignal     : device buffer for one IF signal block (int16 samples)
 * d_caCode         : device buffer for one local code table
 * d_BasebandSignal : device buffer for carrier-wiped baseband samples
 * d_corrValues     : device pointer mapped to host correlation buffer
 * h_corrValues     : host mapped correlation buffer for nine QMBOC branches
 * settings         : cached receiver settings parsed from MATLAB struct
 * blockPerGrid     : grid size used by the carrier wipeoff kernels
 *------------------------------------------------------------------------------*/
static short* d_pRawSignal = NULL;
static signed char* d_caCode = NULL;
static GPU_Complex* d_BasebandSignal = NULL;
static std::complex<float>* h_corrValues = NULL;
static GPU_Complex* d_corrValues = NULL;
static Settings settings;
__constant__ float d_initCodePhase[CORR_NUMBER];


/* The gateway function --------------------------------------------------------
 * Args   : prhs[0]   Settings    settings       I   receiver settings structure
 *          prhs[1]   int16*      rawSignal      I   one code-period IF signal block
 *          prhs[2]   int8*       caCode         I   one PRN local code table arranged as
 *                                              [B1CData pilotBOC11 pilotBOC61]
 *          prhs[3]   double      remCarrPhase   I   residual carrier phase (rad)
 *          prhs[4]   double      carrPhaseStep  I   carrier phase step (rad/sample)
 *          prhs[5]   double      remCodePhase   I   residual code phase (chips)
 *          prhs[6]   double      codePhaseStep  I   code phase step (chips/sample)
 *          prhs[7]   double      inPRN          I   PRN number of the current channel
 *          plhs[0]   double*     corrValues     O   [data E/P/L I/Q,
 *                                                 pilotBOC11 E/P/L I/Q,
 *                                                 pilotBOC61 E/P/L I/Q]
 * Return : plhs[0]   corrValues
 *------------------------------------------------------------------------------*/
void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    /* --------------- Declare all variables --------------- */
    short* rawSignal;
    signed char* caCode;
    double remCarrPhase, carrPhaseStep, remCodePhase, codePhaseStep, inPRN;
    double* corrValues;
    size_t rawSignalLen;
    int blksize;
    static size_t codeLen = 0;
    static size_t branchCodeLen;
    static double PRN = 0.0;
    static int initialized = 0;
    static size_t blockPerGrid = 0;

    if (!initialized)
    {
        /* Validate the MEX interface before dereferencing prhs[k]. Keep this
           check in the first-time initialization path so later calls stay on
           the existing fast path. */
        checkInputs(nrhs, nlhs, prhs);

        // Basic setting
        settings = parseSettings(prhs[0]);

        /* --------------- Find array dimensions ---------------
         * codeLen      : total number of chips in the concatenated
         *                [B1CData pilotBOC11 pilotBOC61] local code table
         * rawSignalLen : number of int16 entries in the first signal block
         * blksize      : number of complex samples processed by the kernels
         *------------------------------------------------------ */
        codeLen = mxGetNumberOfElements(prhs[2]);
        branchCodeLen = codeLen / 3;
        rawSignalLen = mxGetNumberOfElements(prhs[1]);
        blksize = (settings.fileType == 1) ? (int)rawSignalLen : (int)(rawSignalLen / 2);

        /* --------------------- Allocate memory ---------------------
         * Add a small safety margin because blksize can vary slightly
         * across tracking loops as code/carrier NCOs change.
         *------------------------------------------------------------ */
        rawSignalLen = rawSignalLen + ((settings.fileType == 1) ? 100 : 200);

        CUDA_CHECK(cudaSetDeviceFlags(cudaDeviceMapHost));

        /* -------------------- Device-side memory -------------------- */
        CUDA_CHECK(cudaMalloc(&d_pRawSignal, rawSignalLen * sizeof(short)));
        CUDA_CHECK(cudaMalloc(&d_caCode, codeLen * sizeof(signed char)));
        CUDA_CHECK(cudaMalloc(&d_BasebandSignal, (blksize + 100) * sizeof(GPU_Complex)));

        // Grid size for carrier wipeoff kernels
        blockPerGrid = (blksize + 100 + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;


        /* The combined QMBOC code table is arranged as:
             [B1CData(last,1..N,first), pilotBOC11(last,1..N,first),
              pilotBOC61(last,1..N,first)].
           Only the small Early/Prompt/Late phase shift is kept here. The data/pilot
           branch offset is added later as an integer index in the CUDA kernel,
           avoiding loss of fractional code-phase precision in float arithmetic.
           The shifts are used with ceil-style indexing to match the MATLAB
           reference correlator: index = ceil(codePhase + shift). */
        float initCodePhase[CORR_NUMBER] = {
            (float)(-settings.earlyLateSpc),
            0.0f,
            (float)(settings.earlyLateSpc),
            (float)(-settings.earlyLateSpc),
            0.0f,
            (float)(settings.earlyLateSpc),
            (float)(-settings.earlyLateSpc),
            0.0f,
            (float)(settings.earlyLateSpc)
        };

        CUDA_CHECK(cudaMemcpyToSymbol(d_initCodePhase, initCodePhase, sizeof(float) * CORR_NUMBER));


        /* ----------------- Correlation values ----------------
         * Host mapped page-locked memory is used so the GPU can
         * write correlation results directly to a host-visible buffer.
         *------------------------------------------------------ */
        CUDA_CHECK(cudaHostAlloc(&h_corrValues, CORR_NUMBER * sizeof(std::complex<float>), cudaHostAllocMapped));
        CUDA_CHECK(cudaHostGetDevicePointer(&d_corrValues, h_corrValues, 0));

        mexAtExit(cleanup);
        initialized = 1;
        mexPrintf("   CUDA initialized ...\n");
        mexPrintf("   Channel-serial tracking with the GPU-accelerated correlator developed by Yafeng Li.\n");
    }

    /* ----------------- Get input arrays -------------------
     * Input order:
     *   prhs[0] - settings
     *   prhs[1] - rawSignal
     *   prhs[2] - caCode
     *   prhs[3] - remCarrPhase
     *   prhs[4] - carrPhaseStep
     *   prhs[5] - remCodePhase
     *   prhs[6] - codePhaseStep
     *   prhs[7] - PRN
     *--------------------------------------------------------- */
    rawSignal = (short*)mxGetData(prhs[1]);
    remCarrPhase = mxGetScalar(prhs[3]);
    carrPhaseStep = mxGetScalar(prhs[4]);
    remCodePhase = mxGetScalar(prhs[5]);
    codePhaseStep = mxGetScalar(prhs[6]);
    inPRN = mxGetScalar(prhs[7]);

    rawSignalLen = mxGetNumberOfElements(prhs[1]);
    blksize = (settings.fileType == 1) ? (int)rawSignalLen : (int)(rawSignalLen / 2);
    

    /* --------------- Initialize output array --------------- */
    plhs[0] = mxCreateDoubleMatrix((mwSize)1, (mwSize)(2 * CORR_NUMBER), mxREAL);
    corrValues = (double*)mxGetData(plhs[0]);

    /* -------------------- Copy PRN codes ------------------
     * In channel-serial tracking, only one PRN is processed per call.
     * Reuse the device-side code table until the tracked PRN changes.
     *-------------------------------------------------------- */
    if (fabs(PRN - inPRN) > 1e-6) {
        caCode = (signed char*)mxGetData(prhs[2]);
        CUDA_CHECK(cudaMemcpy(d_caCode, caCode, codeLen * sizeof(signed char), cudaMemcpyHostToDevice));
        PRN = inPRN;
    }

    /* ---------- Copy input IF signals to device --------------- */
    CUDA_CHECK(cudaMemcpy(d_pRawSignal, rawSignal, rawSignalLen * sizeof(short), cudaMemcpyHostToDevice));

    /* ------------ CUDA correlator implementation ------------
     * 1) carrier wipeoff to generate complex baseband samples
     * 2) Early / Prompt / Late correlation with the local code
     *--------------------------------------------------------- */
    double remCodePhaseWrapped = fmod(remCodePhase, (double)branchCodeLen);
    int codePhaseBase = (int)floor(remCodePhaseWrapped);
    float remCodePhaseFrac = (float)(remCodePhaseWrapped - (double)codePhaseBase);

    if (settings.fileType == 1) {
        mixCarrReal<<<blockPerGrid, THREADS_PER_BLOCK>>>(
            d_pRawSignal, d_BasebandSignal, blksize, (float)remCarrPhase, (float)carrPhaseStep);
    } else {
        mixCarrComplex<<<blockPerGrid, THREADS_PER_BLOCK>>>(
            d_pRawSignal, d_BasebandSignal, blksize, (float)remCarrPhase, (float)carrPhaseStep);
    }
    CUDA_KERNEL_CHECK();

    correlator<<<CORR_NUMBER, THREADS_PER_BLOCK>>>(
        d_corrValues, d_BasebandSignal, d_caCode,
        (float)codePhaseStep, codePhaseBase, remCodePhaseFrac, blksize,
        (int)branchCodeLen);

    CUDA_KERNEL_CHECK();

    /* ------------ Retrieve results from device ----------------- */
    CUDA_CHECK(cudaDeviceSynchronize());
    /* Output order:
         [data E, data P, data L, pilotBOC11 E, pilotBOC11 P, pilotBOC11 L,
          pilotBOC61 E, pilotBOC61 P, pilotBOC61 L]
       and each branch contributes [I, Q]. */
    for (int ind = 0; ind < CORR_NUMBER; ++ind) {
        corrValues[ind * 2] = (double)h_corrValues[ind].real();
        corrValues[ind * 2 + 1] = (double)h_corrValues[ind].imag();
    }
}

// ============================ CUDA kernels =====================================
/*
 * mixCarrReal
 * ----------
 * Mix a real-valued IF sequence with a complex carrier replica:
 *     bb[n] = raw[n] * exp(-j*(remCarrPhase + n*carrPhaseStep))
 */
__global__ void mixCarrReal(const short* d_pRawSignal, GPU_Complex* d_BasebandSignal, int blksize, float remCarrPhase, float carrPhaseStep)
{
    // CUDA version of floating point NCO and vector dot product integrated

    for (int index = blockIdx.x * blockDim.x + threadIdx.x; index < blksize; index += blockDim.x * gridDim.x) {
        float sin, cos;
        __sincosf(remCarrPhase + index * carrPhaseStep, &sin, &cos);
        d_BasebandSignal[index] = GPU_Complex(cos, -sin) * (float)d_pRawSignal[index];
    }
}

/*
 * mixCarrComplex
 * -------------
 * Mix a complex interleaved IF sequence [I,Q,I,Q,...] with the same local carrier:
 *     bb[n] = (I[n] + jQ[n]) * exp(-j*(remCarrPhase + n*carrPhaseStep))
 */
__global__ void mixCarrComplex(const short* d_pRawSignal, GPU_Complex* d_BasebandSignal,
    int blksize, float remCarrPhase, float carrPhaseStep)
{
    // CUDA version of floating point NCO and vector dot product integrated
    const short2* d_pRawSignalIQ = reinterpret_cast<const short2*>(d_pRawSignal);
    for (int index = blockIdx.x * blockDim.x + threadIdx.x; index < blksize; index += blockDim.x * gridDim.x) {
        float sin, cos;
        short2 rawIQ = d_pRawSignalIQ[index];
        __sincosf(remCarrPhase + index * carrPhaseStep, &sin, &cos);
        d_BasebandSignal[index] = GPU_Complex(cos, -sin) *
            GPU_Complex((float)rawIQ.x, (float)rawIQ.y);
    }
}

/*
 * correlator
 * ----------
 * Compute Early / Prompt / Late correlations.
 * blockIdx.x selects the correlator branch:
 *   vec = 0..2 -> B1CData Early / Prompt / Late
 *   vec = 3..5 -> pilotBOC11 Early / Prompt / Late
 *   vec = 6..8 -> pilotBOC61 Early / Prompt / Late
 *
 * Each thread accumulates a strided partial sum, then a shared-memory tree reduction
 * combines them to obtain one complex correlation value per branch.
 */
__global__ void correlator(GPU_Complex* d_corrValues, GPU_Complex* d_BasebandSignal, const signed char* __restrict__ d_caCodeTable,
    float codePhaseStep, int codePhaseBase, float remCodePhase, int blksize, int branchCodeLen)
{
    /* Each block computes one QMBOC correlator branch:
         block 0..2 -> data Early / Prompt / Late
         block 3..5 -> pilot BOC(1,1) Early / Prompt / Late
         block 6..8 -> pilot BOC(6,1) Early / Prompt / Late */
    __shared__ GPU_Complex accumResult[ACCUM_N];

    for (int vec = blockIdx.x; vec < CORR_NUMBER; vec += gridDim.x)
    {
        int branchOffset = (vec / 3) * branchCodeLen;
        float initPhase = d_initCodePhase[vec] + remCodePhase;

        for (int iAccum = threadIdx.x; iAccum < ACCUM_N; iAccum += blockDim.x)
        {
            GPU_Complex sumIQ(0.0f, 0.0f);

            for (int pos = iAccum; pos < blksize; pos += ACCUM_N) {
                // 1.resample local code for the current shift
                float localPhase = codePhaseStep * __int2float_rd(pos) + initPhase;
                int chipIndex = branchOffset + codePhaseBase + __float2int_ru(localPhase);

                // 2.correlate
                sumIQ.multiply_acc(d_BasebandSignal[pos], (float)__ldg(&d_caCodeTable[chipIndex]));
            }
            accumResult[iAccum] = sumIQ;
        }    //for (int iAccum = threadIdx.x; iAccum < ACCUM_N; iAccum += blockDim.x)

        /* Perform tree-like reduction of accumulators' results.
        ACCUM_N has to be power of two at this stage */
        for (int stride = ACCUM_N / 2; stride > 0; stride >>= 1) {
            __syncthreads();
            for (int iAccum = threadIdx.x; iAccum < stride; iAccum += blockDim.x) {
                accumResult[iAccum] += accumResult[stride + iAccum];
            }
        }

        if (threadIdx.x == 0) {
            d_corrValues[vec] = accumResult[0];
        }

    } // for (int vec = blockIdx.x; vec < CORR_NUMBER; vec += gridDim.x)
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

    cfg.trkMode = getStructIntField(s, "trkMode");
    cfg.fileType = getStructIntField(s, "fileType");
    cfg.earlyLateSpc = getStructDoubleField(s, "dllCorrelatorSpacing");

    if (cfg.fileType != 1 && cfg.fileType != 2) {
        mexErrMsgIdAndTxt("gpuCorr:fileType", "fileType must be 1 (real) or 2 (complex interleaved).");
    }

    return cfg;
}

static void checkInputs(int nrhs, int nlhs, const mxArray* prhs[])
{
    if (nlhs > 1) {
        mexErrMsgIdAndTxt("gpuCorr:Output", "One output is expected.");
    }

    if (nrhs != 8) {
        mexErrMsgIdAndTxt("gpuCorr:InputCount",
            "Expected 8 inputs: settings, rawSignal, caCode, remCarrPhase, carrPhaseStep, remCodePhase, codePhaseStep, PRN.");
    }

    if (!mxIsStruct(prhs[0])) {
        mexErrMsgIdAndTxt("gpuCorr:settings", "settings must be a struct.");
    }
    if (!mxIsInt16(prhs[1]) || mxIsComplex(prhs[1])) {
        mexErrMsgIdAndTxt("gpuCorr:rawSignal", "rawSignal must be a real int16 array.");
    }
    if (!mxIsInt8(prhs[2]) || mxIsComplex(prhs[2])) {
        mexErrMsgIdAndTxt("gpuCorr:caCode", "caCode must be a real int8 vector.");
    }
    if ((mxGetNumberOfElements(prhs[2]) % 3) != 0) {
        mexErrMsgIdAndTxt("gpuCorr:caCodeLength",
            "caCode must contain concatenated [B1CData pilotBOC11 pilotBOC61], so its length must be divisible by 3.");
    }

    for (int k = 3; k <= 7; ++k) {
        if (!mxIsDouble(prhs[k]) || mxIsComplex(prhs[k]) || mxGetNumberOfElements(prhs[k]) != 1) {
            mexErrMsgIdAndTxt("gpuCorr:Type", "Inputs 4-8 must be real double scalars.");
        }
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
    if (d_pRawSignal != NULL) {
        CUDA_CHECK(cudaFree(d_pRawSignal));
        d_pRawSignal = NULL;
    }

    if (d_caCode != NULL) {
        CUDA_CHECK(cudaFree(d_caCode));
        d_caCode = NULL;
    }

    if (d_BasebandSignal != NULL) {
        CUDA_CHECK(cudaFree(d_BasebandSignal));
        d_BasebandSignal = NULL;
    }

    /* --------- Free host-side memory --------- */
    if (h_corrValues != NULL) {
        CUDA_CHECK(cudaFreeHost(h_corrValues));
        h_corrValues = NULL;
        d_corrValues = NULL;
    }
}
