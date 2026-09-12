An Open-Source MATLAB GNSS SDR Toolbox with SIMD/GPU Acceleration: mSoftGNSS 
===============================================================================



Overview
-------------------------------------------------------------------------------
mSoftGNSS 是一个基于 MATLAB 的开源工具箱，用于对录制的 GNSS 中频（IF）信号进行后处理。
该工具箱基于 SoftGNSS 接收机架构，在可配置的处理框架下支持 GPS、Galileo、GLONASS 和北斗信号。
根据各信号已实现的处理能力，提供 CPU/GPU 捕获、数据/导频跟踪、导航电文解码、伪距生成和定位功能。
工具箱提供通道串行和通道并行两种跟踪模式，每种模式均支持 MATLAB 相关器、SIMD 加速相关器和
GPU 加速相关器。通道并行跟踪通过跨通道共享 IF 数据缓冲和批量执行相关任务，减少重复数据访问。
计算密集型相关器由 C++ 和 CUDA C++ 实现，并通过 MEX 接口集成；接收机调度、跟踪环路控制、
导航解码和定位仍在 MATLAB 中完成。这种功能划分将高层算法开发的便利性与原生代码的计算效率相结合。
结合现代化信号支持和 LDPC 解码功能，工具箱为 GNSS 算法研究、接收机原型开发以及基于录制 IF 数据的
可重复评估提供了可扩展的平台。



Authors
-------------------------------------------------------------------------------
* Yafeng Li
    * E-Mail: <lyf8118@126.com>
    * Wechat: lyf8118521
    * GNSS软件接收机技术讨论QQ群：147304049


* Dennis Akos  
    * E-Mail: <dma@colorado.edu>
    * HP: <http://www.colorado.edu/aerospace/dennis-akos>




Features
-------------------------------------------------------------------------------
* 基于 MATLAB 实现的 GNSS 信号处理功能
    * 本地码生成
    * CPU 和 GPU 捕获
    * 通道串行和通道并行跟踪
    * MATLAB、SIMD MEX 和 CUDA MEX 相关器
    * 现代化信号的数据/导频跟踪
    * 导航电文解码（包括 LDPC 解码器）
    * 伪距生成
    * 按已实现的处理范围进行位置/钟差解算和结果绘图
* SIMD/GPU 加速和相关器选项
    * settings.gpuACQflag 用于选择 CPU 捕获或 MATLAB GPU-array 捕获。
    * settings.trkMode 用于选择通道串行跟踪或通道并行跟踪。
    * settings.correlatorType = 0 表示使用 MATLAB 参考相关器。
    * settings.correlatorType = 1 表示使用 CPU SIMD MEX 相关器。
    * settings.correlatorType = 2 表示使用 CUDA GPU MEX 相关器。
    * MATLAB 相关器包含 BPSK、QPSK、TMBPSK、QMBOC等变体，
      位于各接收机的 include 目录中。
    * 共用 MEX 相关器后端包含 BPSK、QPSK 和 QMBOC 的
      SIMD 和 CUDA 实现。
* 支持的信号
    * GPS L1 C/A
    * GPS L1C
    * GPS L2C (data + pilot)
    * GPS L5 (data + pilot)
    * Galileo E1 (data + pilot)
    * Galileo E5a (data + pilot)
    * Galileo E5b (data + pilot)
    * Galileo E6B（HAS 页/电文组帧，不独立定位）
    * GLONASS L1OF
    * GLONASS L2OF
    * GLONASS L1OC
    * GLONASS L2OC（仅导频跟踪）
    * GLONASS L3OC
    * BeiDou B1I/B2I
    * BeiDou B3I
    * BDS-3 B1C (data + pilot)
    * BDS-3 B2a (data + pilot)
    * BeiDou B2b (data)
* 支持 RF 二进制文件后处理
    * 通过 settings.fileType 配置实采样或 I/Q 复采样 IF 数据。
    * 通过 settings.dataType 配置 int8 或 int16 输入采样格式。
    * README 后面保留 LimeSDR 复数 IF 和 NUT4NT 实采样数据示例；
      参数应与录制数据的元数据一致，不能只修改文件名。
    * 当前接收机在较新的 Windows MATLAB 环境中维护和测试，
      包括 MATLAB R2025a/R2025b。



Directory and Files
-------------------------------------------------------------------------------
根目录
    ./Doc                     文档、ICD、论文和各接收机说明材料。
    ./IF_Data_Set             可放置待处理 IF 数据及对应 metadata 文件。
    ./native_Correlators      所有接收机共用的 SIMD/CUDA MEX 相关器源码和
                              已编译的 .mexw64 文件。

接收机目录
    每个接收机目录都是相对独立的软件接收机；除 ./native_Correlators、
    ./IF_Data_Set 和 ./Doc 外，各接收机代码互不依赖。当前接收机目录如下：

    ./mGPS_L1CA               GPS L1 C/A 软件接收机
    ./mGPS_L1C                GPS L1C 软件接收机
    ./mGPS_L2C                GPS L2C 软件接收机
    ./mGPS_L5                 GPS L5 软件接收机
    ./mGalileo_E1             Galileo E1 软件接收机
    ./mGalileo_E5a            Galileo E5a 软件接收机
    ./mGalileo_E5b            Galileo E5b 软件接收机
    ./mGalileo_E6B            Galileo E6B 软件接收机
    ./mGlonass_L1_L2          GLONASS L1/L2 软件接收机
    ./mGlonass_L1OC           GLONASS L1OC 软件接收机
    ./mGlonass_L2OC           GLONASS L2OC 软件接收机
    ./mGlonass_L3OC           GLONASS L3OC 软件接收机
    ./mBDS_B1I                北斗 B1I/B2I 软件接收机
    ./mBDS_B3I                北斗 B3I 软件接收机
    ./mBDS-3_B1C              北斗三号 B1C 软件接收机
    ./mBDS-3_B2a              北斗三号 B2a 软件接收机
    ./mBDS_B2b                北斗 B2b 软件接收机

每个接收机目录的通用结构
    ./init.m                  接收机启动脚本：设置路径、读取配置、探测原始 IF
                              数据并启动处理。
    ./initSettings.m          当前信号接收机的参数配置文件。
    ./postProcessing.m        顶层处理脚本：调度捕获、通道初始化、跟踪、导航解算、
                              结果保存和绘图。
    ./acquisitionCPU.m        CPU/MATLAB 捕获实现。
    ./acquisitionGPU.m        MATLAB GPU-array 捕获实现。
    ./trkChannelsSerial.m     通道串行跟踪实现。
    ./trkChannelsParallel.m   通道并行跟踪实现。
    ./postNavigation.m        导航电文解码和定位解算。
                              ./mGlonass_L3OC 的该文件位于 ./include/postNavigation.m。
    ./include                 与当前信号相关的辅助函数、本地码生成、MATLAB 相关器、
                              导航电文解码、绘图和通道初始化函数。
    ./include/ldpcDecoder     LDPC 解码相关函数；仅存在于需要 LDPC 导航解码的接收机，
                              例如 mBDS-3_B1C、mBDS-3_B2a、mBDS_B2b 和 mGPS_L1C。
    ./Common                  当前接收机本地公共数学/定位工具函数，例如伪距计算、
                              坐标转换、最小二乘定位和环路参数计算。

./native_Correlators 中的共用相关器后端
    corrSIMDSerial*.cpp       SIMD 通道串行跟踪相关器源码。
    corrSIMDParallel*.cpp     SIMD 通道并行跟踪相关器源码。
    corrGPUSerial*.cu         CUDA 通道串行跟踪相关器源码。
    corrGPUParallel*.cu       CUDA 通道并行跟踪相关器源码。
    corrGPUParallelFused*.cu  CUDA fused 通道并行跟踪相关器源码。
    *.mexw64                  已编译 MEX 文件，由 settings.correlatorType 选择调用。

相关器选择
    settings.correlatorType = 0
                              使用各接收机 ./include 目录中的 MATLAB 参考相关器，
                              便于阅读、调试和对比。
    settings.correlatorType = 1
                              使用 ./native_Correlators 中的 CPU SIMD MEX 相关器，
                              不依赖 CUDA 即可加速跟踪。
    settings.correlatorType = 2
                              使用 ./native_Correlators 中的 CUDA GPU MEX 相关器，
                              用于高吞吐的通道串行或通道并行跟踪。

Software Dependencies
-------------------------------------------------------------------------------
* 需要 MATLAB。当前工程包含 Windows .mexw64 文件，主要面向 Windows 下的
  较新 MATLAB 环境。
* IF 频谱绘图使用 MATLAB Signal Processing Toolbox 中的 pwelch。
* 导航电文解码按信号需要使用 MATLAB Communications Toolbox 中的部分函数：
  -- Create Galois field array: gf()
  -- BCH decoder: bchdec()
  -- Detect errors in input data using CRC: comm.CRCDetector()
  -- Convert convolutional code polynomials to trellis description: poly2trellis()
  -- Convolutionally decode binary data using Viterbi algorithm: vitdec()
  -- Detect errors in input data using CRC: step()
* GPU-array 捕获需要 MATLAB Parallel Computing Toolbox 和支持 CUDA 的 GPU。
  直接运行 CUDA MEX 相关器需要兼容的 NVIDIA GPU、驱动和 CUDA 运行库。
* SIMD MEX 需要支持 AVX2 的 CPU；重新编译需要 MATLAB 支持的 C/C++ 编译器。
  重新编译 CUDA MEX 需要可用的 CUDA toolkit/compiler 和 mexcuda 环境
  （Parallel Computing Toolbox）。
* 如果只需要 CPU 运行，可在各接收机 initSettings.m 中关闭 GPU 捕获，并选择
  MATLAB 相关器或 SIMD 相关器。
  
  
 
How to use
-------------------------------------------------------------------------------
* Step 1: 进入目标信号对应的接收机目录，例如 "mBDS-3_B1C" 或 "mGPS_L1CA"。
* Step 2: 将 IF 数据放入 "IF_Data_Set"，或在 settings.fileName 中设置数据文件
          的绝对路径。
* Step 3: 在 "initSettings.m" 中配置处理时长、IF 数据、捕获、跟踪、相关器和
          导航解算参数。
* Step 4: 通过 settings.gpuACQflag 选择 CPU/GPU 捕获，通过 settings.trkMode
          选择通道串行/并行跟踪，通过 settings.correlatorType 选择 MATLAB、
          SIMD MEX 或 GPU MEX 相关器。
* Step 5: 在当前接收机目录下运行 "init.m" 启动处理。
* Step 6: 通过 "postProcessing.m" 和 "postNavigation.m" 生成并查看捕获、跟踪、
          导航、定位和绘图结果。



Implementation details
-------------------------------------------------------------------------------
本软件接收机套件沿用了 SoftGNSS 风格的 MATLAB SDR 组织方式，并扩展了多种
现代化 GNSS 信号、CPU/GPU 捕获、通道串行/并行跟踪，以及共用 SIMD/CUDA MEX
相关器后端。

BDS-3 B1C/B2a 接收机设计可参考：
Li, Y., Shivaramaiah, N.C. & Akos, D.M. Design and implementation of an open-source
BDS-3 B1C/B2a SDR receiver. GPS Solut (2019) 23: 60.
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
