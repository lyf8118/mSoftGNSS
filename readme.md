An Open-Source MATLAB GNSS SDR Toolbox with SIMD/GPU Acceleration: mSoftGNSS 
===============================================================================



Overview
-------------------------------------------------------------------------------
mSoftGNSS is an open-source MATLAB-based toolbox for post-processing recorded 
GNSS intermediate-frequency (IF) signals. Building on the SoftGNSS receiver 
architecture, it supports GPS, Galileo, GLONASS, and BeiDou signals within a 
configurable processing framework. Signal-dependent capabilities include 
CPU/GPU acquisition, data/pilot tracking, navigation-message decoding, 
pseudorange generation, and positioning. The toolbox provides channel-serial 
and channel-parallel tracking modes, each supporting MATLAB-based, 
SIMD-accelerated, and GPU-accelerated correlators. Channel-parallel tracking 
shares IF data buffers across channels and batches correlation tasks to reduce 
repeated data access. Computationally intensive correlators are implemented in 
C++ and CUDA C++ and integrated through MEX interfaces, while receiver scheduling, 
tracking-loop control, navigation decoding, and positioning remain in MATLAB. 
This separation combines the accessibility of high-level algorithm development 
with efficient native computation. Together with modernized-signal support and 
LDPC decoding, the toolbox provides an extensible platform for GNSS algorithm 
research, receiver prototyping, and reproducible evaluation using recorded IF data.



Authors
-------------------------------------------------------------------------------
Yafeng Li
E-Mail: <lyf8118@126.com>
Wechat: lyf8118521
QQ Group for Technical Discussions on GNSS Software Receivers: 147304049


Dennis Akos  
E-Mail: <dma@colorado.edu>
HP: <http://www.colorado.edu/aerospace/dennis-akos>



Features
-------------------------------------------------------------------------------
* GNSS signal processing functions written in MATLAB
    * Local code generation
    * CPU and GPU acquisition
    * Channel-serial and channel-parallel tracking
    * MATLAB, SIMD MEX, and GPU/CUDA MEX correlators
    * Data/pilot tracking for modernized signals
    * Navigation-message decoding (including LDPC decoder)
    * Pseudorange generation
    * Position/clock-bias calculation and result plotting where implemented
* SIMD/GPU acceleration and correlator options
    * settings.gpuACQflag selects CPU or MATLAB GPU-array acquisition.
    * settings.trkMode selects channel-serial or channel-parallel tracking.
    * settings.correlatorType = 0 uses MATLAB reference correlators.
    * settings.correlatorType = 1 uses CPU SIMD MEX correlators.
    * settings.correlatorType = 2 uses CUDA GPU MEX correlators.
    * MATLAB correlator variants include BPSK, QPSK, TMBPSK, and QMBOC
      implementations inside each receiver's include folder.
    * Shared MEX correlator backends include BPSK, QPSK, and QMBOC
      SIMD and CUDA variants.
* Supported signals
    * GPS L1 C/A
    * GPS L1C
    * GPS L2C (data + pilot)
    * GPS L5 (data + pilot)
    * Galileo E1 (data + pilot)
    * Galileo E5a (data + pilot)
    * Galileo E5b (data + pilot)
    * Galileo E6B (HAS page/message assembly; no standalone positioning)
    * GLONASS L1OF
    * GLONASS L2OF
    * GLONASS L1OC
    * GLONASS L2OC (pilot tracking only)
    * GLONASS L3OC
    * BeiDou B1I/B2I
    * BeiDou B3I
    * BDS-3 B1C (data + pilot)
    * BDS-3 B2a (data + pilot)
    * BeiDou B2b (data)
* RF binary-file post processing
    * Supports real-sample and I/Q-sample IF files configured by
      settings.fileType.
    * Supports int8 and int16 input sample formats configured by
      settings.dataType.
    * LimeSDR complex IF and NUT4NT real-sample examples are documented below;
      match the settings to the recording metadata, not just the file name.
    * The current SDR set has been maintained and tested with recent MATLAB
      releases on Windows, including MATLAB R2025a/R2025b-era workflows.
	


Directory and Files
-------------------------------------------------------------------------------
Root folders
    ./Doc                     Documentation, ICD material, papers, and receiver
                              summary documents.
    ./IF_Data_Set             Optional location for IF data files and metadata.
    ./native_Correlators      Shared SIMD and CUDA MEX correlator source files
                              and compiled .mexw64 binaries.

Receiver folders
    Each receiver folder is self-contained except for the shared ./native_Correlators,
    ./IF_Data_Set, and ./Doc folders. The current receiver folders are:

    ./mGPS_L1CA               GPS L1 C/A SDR receiver
    ./mGPS_L1C                GPS L1C SDR receiver
    ./mGPS_L2C                GPS L2C SDR receiver
    ./mGPS_L5                 GPS L5 SDR receiver
    ./mGalileo_E1             Galileo E1 SDR receiver
    ./mGalileo_E5a            Galileo E5a SDR receiver
    ./mGalileo_E5b            Galileo E5b SDR receiver
    ./mGalileo_E6B            Galileo E6B SDR receiver
    ./mGlonass_L1_L2          GLONASS L1/L2 SDR receiver
    ./mGlonass_L1OC           GLONASS L1OC SDR receiver
    ./mGlonass_L2OC           GLONASS L2OC SDR receiver
    ./mGlonass_L3OC           GLONASS L3OC SDR receiver
    ./mBDS_B1I                BeiDou B1I/B2I SDR receiver
    ./mBDS_B3I                BeiDou B3I SDR receiver
    ./mBDS-3_B1C              BDS-3 B1C SDR receiver
    ./mBDS-3_B2a              BDS-3 B2a SDR receiver
    ./mBDS_B2b                BeiDou B2b SDR receiver

Common layout inside each receiver folder
    ./init.m                  Receiver startup script. It sets paths, loads
                              settings, probes raw IF data, and starts processing.
    ./initSettings.m          Receiver-specific parameter configuration.
    ./postProcessing.m        Top-level processing script for acquisition,
                              channel initialization, tracking, navigation,
                              result saving, and plotting.
    ./acquisitionCPU.m        MATLAB/CPU acquisition implementation.
    ./acquisitionGPU.m        MATLAB GPU-array acquisition implementation.
    ./trkChannelsSerial.m     Channel-serial tracking implementation.
    ./trkChannelsParallel.m   Channel-parallel tracking implementation.
    ./postNavigation.m        Navigation-message decoding and position computation.
                              For ./mGlonass_L3OC this file is located under
                              ./include/postNavigation.m.
    ./include                 Signal-specific helper functions, local-code
                              generators, MATLAB correlators, navigation-message
                              decoders, plotting helpers, and channel setup.
    ./include/ldpcDecoder     LDPC decoder support files, present only in
                              receivers that need LDPC-based navigation decoding
                              such as mBDS-3_B1C, mBDS-3_B2a, mBDS_B2b, and mGPS_L1C.
    ./Common                  Receiver-local common math and positioning utilities such as
                              pseudorange calculation, coordinate transforms,
                              least-squares positioning, and loop coefficients.

Shared correlator backends in ./native_Correlators
    corrSIMDSerial*.cpp       SIMD serial tracking correlator source files.
    corrSIMDParallel*.cpp     SIMD channel-parallel tracking correlator source files.
    corrGPUSerial*.cu         CUDA serial tracking correlator source files.
    corrGPUParallel*.cu       CUDA channel-parallel tracking correlator source files.
    corrGPUParallelFused*.cu  CUDA fused channel-parallel tracking correlator
                              source files where available.
    *.mexw64                  Compiled MEX binaries used by settings.correlatorType.

Correlator selection
    settings.correlatorType = 0
                              MATLAB reference correlators in the receiver's
                              ./include folder. These are easiest to inspect
                              and debug.
    settings.correlatorType = 1
                              CPU SIMD MEX correlators from ./native_Correlators.
                              These accelerate tracking without requiring CUDA.
    settings.correlatorType = 2
                              CUDA GPU MEX correlators from ./native_Correlators.
                              These are intended for high-throughput serial or
                              channel-parallel tracking.



Software Dependencies
-------------------------------------------------------------------------------
* MATLAB is required. The current code base is maintained with recent MATLAB
  releases on Windows; the bundled MEX binaries are Windows .mexw64 files.
* MATLAB Signal Processing Toolbox is used by IF spectrum plotting (pwelch).
* MATLAB Communications Toolbox is needed by signal-dependent navigation decoding:
  -- Create Galois field array: gf()
  -- BCH decoder: bchdec()
  -- Detect errors in input data using CRC: comm.CRCDetector()
  -- Convert convolutional code polynomials to trellis description: poly2trellis()
  -- Convolutionally decode binary data using Viterbi algorithm: vitdec()
  -- Detect errors in input data using CRC: step()
* MATLAB Parallel Computing Toolbox and a CUDA-capable GPU are needed for
  GPU-array acquisition. Direct CUDA MEX execution requires a compatible
  NVIDIA GPU, driver, and CUDA runtime.
* An AVX2-capable CPU is needed for SIMD MEX correlators; a supported C/C++
  compiler is needed to rebuild them with mex.
  A supported CUDA toolkit/compiler setup is needed to rebuild CUDA MEX files
  with mexcuda (Parallel Computing Toolbox).
* CPU-only operation is possible by selecting MATLAB or SIMD correlators and
  disabling GPU acquisition in each receiver's initSettings.m.
  
  
 
How to use
-------------------------------------------------------------------------------
* Step 1: Enter the receiver folder for the target signal, for example
          "mBDS-3_B1C" or "mGPS_L1CA".
* Step 2: Put the IF data file under "IF_Data_Set", or set an absolute file
          path in settings.fileName.
* Step 3: Configure processing, IF-data, acquisition, tracking, correlator, and
          navigation parameters in "initSettings.m".
* Step 4: Select CPU/GPU acquisition with settings.gpuACQflag, channel-serial or
          channel-parallel tracking with settings.trkMode, and MATLAB/SIMD/GPU
          correlation with settings.correlatorType.
* Step 5: Start processing by running "init.m" from the receiver folder.
* Step 6: Review acquisition, tracking, navigation, position, and plotting outputs
          generated by "postProcessing.m" and "postNavigation.m".



Implementation details
-------------------------------------------------------------------------------
The receiver set follows the original SoftGNSS-style MATLAB SDR organization
and extends it with modern GNSS signals, CPU/GPU acquisition, channel-serial and
channel-parallel tracking, and shared SIMD/CUDA MEX correlator backends.

See the reference for the BDS-3 B1C/B2a receiver design:
Li, Y., Shivaramaiah, N.C. & Akos, D.M. Design and implementation of an open-source 
BDS-3 B1C/B2a SDR receiver. GPS Solut 23, 60 (2019). 
https://doi.org/10.1007/s10291-019-0853-z


	   
Test signal (collected by LimeSDR) and parameter configurations
-------------------------------------------------------------------------------
    * L1CA_E1_B1C.bin
       -- for GPS L1C/A, Galielo E1 and BDS B1C 
       -- dataType: int16
       -- fileType: 16 bit complex samples
       -- IF: 0 Hz
       -- sampling Frequency: 30e6 Hz
       -- link: https://pan.baidu.com/s/1NdZeVrEdEjDnJiV91Cg0qg?pwd=bvh5 password: bvh5
    * GPS L1C.iq
       -- for GPS L1C 
       -- dataType: int8
       -- fileType: 8 bit complex samples
       -- IF: 0 Hz
       -- sampling Frequency: 25e6 Hz
       -- link: https://pan.baidu.com/s/1hZ4S8R4VWB5OtPr9VyPCYQ?pwd=tv8b password: tv8b
    * GPS_L2.bin
       -- for GPS L2C
       -- dataType: int16
       -- fileType: 16 bit complex samples
       -- IF: 0 Hz
       -- sampling Frequency: 30e6 Hz
       -- link: https://pan.baidu.com/s/1KqdyWyuVPAgx0EoUdoLwag?pwd=5tya password: 5tya
    * L5_E5a_B2a.bin
       -- for GPS L5, Galielo E5a and BDS B2a
       -- dataType: int16
       -- fileType: 16 bit complex samples
       -- IF: 0 Hz
       -- sampling Frequency: 30e6 Hz
       -- link: https://pan.baidu.com/s/1rX2lguJgMmgjAdZIJK1FNg?pwd=bnzm password: bnzm
    
    * E5b_B2b.bin   
       -- for Galielo E5b, BDS B2b
       -- dataType: int16
       -- fileType: 16 bit complex samples
       -- IF: 0 Hz
       -- sampling Frequency: 30e6 Hz
       -- link: https://pan.baidu.com/s/128dwHIKUG9iM0MigS_7vJw?pwd=x4uq password: x4uq
    * GAL_E6.bin  
       -- for Galielo E6
       -- dataType: int16
       -- fileType: 16 bit complex samples
       -- IF: 0 Hz
       -- sampling Frequency: 30e6 Hz
       -- link: https://pan.baidu.com/s/1OrHtpNG-dQxSmmVc9rX6tA?pwd=qfx8 password: qfx8
    * GLO_L1OC_L1.bin 
       -- for GLONASS L1OF  
       -- dataType: int16
       -- fileType: 16 bit complex samples
       -- IF: 1.005e6 Hz
       -- sampling Frequency: 30e6 Hz
       -- GLONASS L1OC
       -- dataType: int16
       -- fileType: 16 bit complex samples
       -- IF: 0 Hz
       -- sampling Frequency: 30e6 Hz
       -- link: https://pan.baidu.com/s/1w-fYtmM3oW_pxR4ZH_HGyw?pwd=js93 password: js93
    * GLO_L2OC_L2.bin  
       -- for GLONASS L2OF 
       -- dataType: int16
       -- fileType: 16 bit complex samples
       -- IF: -2.06e6 Hz
       -- sampling Frequency: 30e6 Hz
    
       -- for GLONASS L2OC 
       -- dataType: int16
       -- fileType: 16 bit complex samples
       -- IF: 0 Hz
       -- sampling Frequency: 30e6 Hz
       -- link: https://pan.baidu.com/s/1jefXiCkeI9qHSqwZ_RXr2g?pwd=33w6 password: 33w6
    * GLO_L3OC.bin  
       -- for GLONASS L3OC	
       -- dataType: int16
       -- fileType: 16 bit complex samples
       -- IF: 0 Hz
       -- sampling Frequency: 30e6 Hz
       -- link: https://pan.baidu.com/s/1ySNHn3W66Q6n9fxERR7nOQ?pwd=289m password: 289m
    * BDS_B1I.bin
       -- for BDS B1I 
       -- dataType: int16
       -- fileType: 16 bit complex samples
       -- IF: 0 Hz
       -- sampling Frequency: 30e6 Hz
       -- link: https://pan.baidu.com/s/1qV-TCtROGHL8iJZlqw5mZA?pwd=39vd password: 39vd
    * BDS_B3I.bin       
       -- for BDS B3I    
       -- dataType: int16
       -- fileType: 16 bit complex samples
       -- IF: 0 Hz
       -- sampling Frequency: 30e6 Hz
       -- link: https://pan.baidu.com/s/1T_g1to1e8OW5OFtHlQsZeg?pwd=6w1u password: 6w1u 



Test signal (collected by NUT4NT) and parameter configurations
-------------------------------------------------------------------------------
    * L1CA_E1_B1C_real.bin
       -- for GPS L1 C/A, Galileo E1 and BDS B1C
       -- dataType: int8
       -- fileType: 8 bit real samples
       -- IF: -14.58e6 Hz
       -- sampling Frequency: 53e6 Hz
       -- link: https://pan.baidu.com/s/1cFGCBwcWA0vHBeiL7WU_2Q?pwd=r8au password: r8au
    * B1I_real.bin
       -- for BDS B1I
       -- dataType: int8
       -- fileType: 8 bit real samples
       -- IF: -28.902e6 Hz
       -- sampling Frequency: 99.375e6 Hz
       -- link: https://pan.baidu.com/s/1y4LfFZSO-Qm3Tvshvg_Gzg?pwd=86ii password: 86ii
    * GLO_L1OF_L1OC_real.bin
       -- for GLONASS L1OF
       -- dataType: int8
       -- fileType: 8 bit real samples
       -- IF: +12e6 Hz
       -- sampling Frequency: 53e6 Hz
    
       -- for GLONASS L1OC
       -- dataType: int8
       -- fileType: 8 bit real samples
       -- IF: +10.995e6 Hz
       -- sampling Frequency: 53e6 Hz
       -- link: https://pan.baidu.com/s/1-fxXe6VW36izl2ve3dMQsw?pwd=j7p9 password: j7p9
    * GLO_L2OF_L2OC_real.bin
       -- for GLONASS L2OF
       -- dataType: int8
       -- fileType: 8 bit real samples
       -- IF: -14e6 Hz
       -- sampling Frequency: 70e6 Hz
    
       -- for GLONASS L2OC
       -- dataType: int8
       -- fileType: 8 bit real samples
       -- IF: -11.94e6 Hz
       -- sampling Frequency: 70e6 Hz
       -- link: https://pan.baidu.com/s/1P-xWe-aAUu3esymnbsPd-w?pwd=pnr8 password: pnr8
    * L2C_real.bin
       -- for GPS L2C (L2CM and L2CL)
       -- dataType: int8
       -- fileType: 8 bit real samples
       -- IF: -7.4e6 Hz
       -- sampling Frequency: 53e6 Hz
       -- link: https://pan.baidu.com/s/16OZPphh6SAObSKf1t2ykPA?pwd=ss7k password: ss7k
    * E5b_B2b_GLO_L3OC_real.bin
       -- for Galileo E5b, BDS B2b
       -- dataType: int8
       -- fileType: 8 bit real samples
       -- IF: +17.14e6 Hz
       -- sampling Frequency: 70e6 Hz
    
       -- for GLONASS L3OC
       -- dataType: int8
       -- fileType: 8 bit real samples
       -- IF: +12.025e6 Hz
       -- sampling Frequency: 70e6 Hz
       -- link: https://pan.baidu.com/s/1Z-LK5UZ_DQgQk654cY0LgA?pwd=puzw password: puzw
    * L5_E5a_B2a_real.bin
       -- for GPS L5, Galileo E5a and BDS B2a
       -- dataType: int8
       -- fileType: 8 bit real samples
       -- IF: -13.55e6 Hz
       -- sampling Frequency: 99.375e6 Hz
       -- link: https://pan.baidu.com/s/1xlftENWf76QwEVzzmxSUDw?pwd=m74w password: m74w
    * B3I_E6_real.bin
       -- for BDS B3I
       -- dataType: int8
       -- fileType: 8 bit real samples
       -- IF: +8.52e6 Hz
       -- sampling Frequency: 70e6 Hz
    
       -- for Galileo E6B
       -- dataType: int8
       -- fileType: 8 bit real samples
       -- IF: +18.75e6 Hz
       -- sampling Frequency: 70e6 Hz
       -- link: https://pan.baidu.com/s/14q3K1MKNnAtxz64mUu9DKw?pwd=unn8 password: unn8
