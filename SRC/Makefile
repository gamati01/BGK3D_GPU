#----------------------------------------------
# makefile BGK3D_GPU
# Giorgio Amati, CASPUR 2000-2012
# Giorgio Amati, CINECA 2013-202?

# targets:
#
#		<none> default ---> GPU code with do concurrent (nvfortran)
#		offload ---> GPU code with openmp offload
#		openacc ---> GPU code with openACC
#		serial  ---> Serial CPU code
#		mulitcore -> Multicore CPU code
#
#----------------------------------------------
# BUILD EXAMPLES (copy/paste)
#----------------------------------------------
# Notes:
#  - Vendor: NVIDIA=1 or AMD=1 (GNU is the default if neither is set).
#  - Precision defaults to SINGLE; override with DOUBLE=1 / MIXED1=1 / MIXED2=1.
#  - A test case is required: LDC (cavity), POF, TGV, VKS, KVX, HIT.
#  - Native GPU targets (cuda/hip) require the fused collision: FUSED=1.
#  - GPU arch:  NVIDIA -> NV_ARCH=sm_80 ;  AMD -> GFX_ARCH=gfx942.
#  - ENERGY=1    : GPU power/energy monitoring (NVML on NVIDIA, AMD SMI on AMD).
#  - PROFILING=1 : per-phase timers (adds device barriers -> slower; profiling only).
#
# 1) OFFLOAD  (OpenMP target offload)
#    NVIDIA:
#      make offload NVIDIA=1 LDC=1
#      make offload NVIDIA=1 LDC=1 ENERGY=1
#      make offload NVIDIA=1 LDC=1 PROFILING=1
#      make offload NVIDIA=1 LDC=1 ENERGY=1 PROFILING=1
#    AMD:
#      make offload AMD=1 GFX_ARCH=gfx942 LDC=1
#      make offload AMD=1 GFX_ARCH=gfx942 LDC=1 ENERGY=1
#      make offload AMD=1 GFX_ARCH=gfx942 LDC=1 PROFILING=1
#      make offload AMD=1 GFX_ARCH=gfx942 LDC=1 ENERGY=1 PROFILING=1
#
# 1b) OFFLOAD KERNEL SYNTAX  (OpenMP Offload Kernel Language)
#    Adds KSYNTAX=1 to the offload target. This selects the
#    "target teams parallel" kernel-language loops, the omp_get_*_dim
#    runtime API and the num_teams/num_threads dims(N) modifier.
#    Requirements / notes:
#      - Needs OpenMP 6.1, so the build injects -fopenmp-version=61.
#      - Needs the upstream LLVM flang/clang (module: aomp_ksyntax); the
#        Makefile automatically swaps amdflang/amdclang for flang/clang.
#      - The executable name is tagged (../RUN/bgk3d.offload_ksyntax.x)
#        so it does not clobber the plain offload build.
#      - Add KS3D=1 for the 3D thread-block (8x8x4) launch variant, which
#        instead builds ../RUN/bgk3d.offload_ksyntax_3d.x.
#    AMD:
#      make offload AMD=1 GFX_ARCH=gfx942 LDC=1 KSYNTAX=1
#      make offload AMD=1 GFX_ARCH=gfx942 LDC=1 KSYNTAX=1 KS3D=1
#      make offload AMD=1 GFX_ARCH=gfx942 LDC=1 KSYNTAX=1 ENERGY=1
#      make offload AMD=1 GFX_ARCH=gfx942 LDC=1 KSYNTAX=1 PROFILING=1
#    NVIDIA:
#      make offload NVIDIA=1 LDC=1 KSYNTAX=1
#      make offload NVIDIA=1 LDC=1 KSYNTAX=1 KS3D=1
#
# 2) NATIVE  (CUDA on NVIDIA, HIP on AMD; needs FUSED=1)
#    NVIDIA (cuda):
#      make cuda NVIDIA=1 FUSED=1 LDC=1 NV_ARCH=sm_80
#      make cuda NVIDIA=1 FUSED=1 LDC=1 NV_ARCH=sm_80 ENERGY=1
#      make cuda NVIDIA=1 FUSED=1 LDC=1 NV_ARCH=sm_80 PROFILING=1
#      make cuda NVIDIA=1 FUSED=1 LDC=1 NV_ARCH=sm_80 ENERGY=1 PROFILING=1
#    AMD (hip):
#      make hip AMD=1 FUSED=1 LDC=1 GFX_ARCH=gfx942
#      make hip AMD=1 FUSED=1 LDC=1 GFX_ARCH=gfx942 ENERGY=1
#      make hip AMD=1 FUSED=1 LDC=1 GFX_ARCH=gfx942 PROFILING=1
#      make hip AMD=1 FUSED=1 LDC=1 GFX_ARCH=gfx942 ENERGY=1 PROFILING=1
#
# 3) DO CONCURRENT  (stdpar GPU offload)
#    NVIDIA:
#      make doconcurrent NVIDIA=1 LDC=1
#      make doconcurrent NVIDIA=1 LDC=1 ENERGY=1
#      make doconcurrent NVIDIA=1 LDC=1 PROFILING=1
#      make doconcurrent NVIDIA=1 LDC=1 ENERGY=1 PROFILING=1
#    AMD:
#      NOT SUPPORTED: the AMD vendor block sets DC = "not supported"
#      (no stdpar GPU offload via amdflang here). Use the native HIP
#      target (section 2, AMD) instead.
#----------------------------------------------
DATE = $(shell date +%d%m%y)
FINE = "******** End Compilation step  *******"
CHECK1 = NVIDIA
CHECK2 = SINGLE
#----------------------------------------------
# Compilers & Flags
#----------------------------------------------
# flags
# default (GNU)
MULTI =   -fopenmp -DMULTICORE -ftree-parallelize-loops=4
OPENACC = -fopenacc -DNOMANAGED
OFFLOAD = -fopenmp  -foffload=nvptx-none -DOFFLOAD
DC = not supported
SFC = gfortran
SCC = gcc
FOPT =  -Ofast
COPT =  -O3
OPT =   -Ofast
DEBUG = -Wuninitialized -DDEBUG_1 -Wunused
#
CHECK1=GNU
CHECK2=SINGLE
CHECK3=ORIGINAL
CHECK4=LDC
#
# GPU energy/power monitoring wrapper (NVML by default, AMD SMI under AMD)
GPU_ENERGY_OBJ = nvml_wrapper.o
# Matching Fortran interface module (NVML by default, AMD SMI under AMD)
GPU_IFACE_MOD = nvml_interface_module.o
#
# Native GPU port of the fused collision (HIP / CUDA back-ends)
HIPCC     ?= hipcc
NVCC      ?= nvcc
GFX_ARCH  ?= gfx90a
NV_ARCH   ?= sm_80
ROCM_PATH ?= /opt/rocm
CUDA_PATH ?= /usr/local/cuda
#
ifdef NVIDIA
 OFFLOAD = -mp=gpu -DOFFLOAD  -Minfo=accel
 # G100
 OPENACC = -Minfo=acc -O2 -acc -DOPENACC -ta=tesla:cc70,managed
 # A100
 # OPENACC = -Minfo=acc -O2 -acc -DOPENACC -ta=tesla:cc80,managed
 # H100
 # OPENACC = -Minfo=acc -O2 -acc -DOPENACC -acc -gpu=cc90,mem:managed
 MULTI = -stdpar=multicore -Minfo=all -DMULTICORE
 DC  = -stdpar=gpu -DDC
 SCC = nvcc
 SFC = scorep-nvfortran  -DPGI
 SFC = nvfortran -DPGI
 FOPT = -O2 -Mnodepchk -Mcontiguous
 COPT = -O2
 OPT  = -O2 -Mnodepchk -Mcontiguous
 CHECK1=NVIDIA
 # fix for ENERGY/G100
 INC1=/g100_work/PROJECTS/spack/preprod/01/install/0.18/linux-centos8-skylake_avx512/gcc-8.4.1/nvhpc-22.3-jgtetvzbyjbc2o6qpmlcjcjxkp5vi7r5/Linux_x86_64/22.3/cuda/11.6/targets/x86_64-linux/include
 LIB1=/g100_work/PROJECTS/spack/preprod/01/install/0.18/linux-centos8-skylake_avx512/gcc-8.4.1/nvhpc-22.3-jgtetvzbyjbc2o6qpmlcjcjxkp5vi7r5/Linux_x86_64/22.3/cuda/11.6/targets/x86_64-linux/lib/stubs -lnvidia-ml
 # fix for ENERGY/Leonardo
 #INC1=/leonardo/prod/spack/5.2/install/0.21/linux-rhel8-icelake/gcc-8.5.0/nvhpc-24.3-v63z4inohb4ywjeggzhlhiuvuoejr2le/Linux_x86_64/24.3/cuda/12.3/targets/x86_64-linux/include
 #LIB1=/leonardo/prod/spack/5.2/install/0.21/linux-rhel8-icelake/gcc-8.5.0/nvhpc-24.3-v63z4inohb4ywjeggzhlhiuvuoejr2le/Linux_x86_64/24.3/cuda/12.3/targets/x86_64-linux/lib/stubs -lnvidia-ml
 # fix for ENERGY/Pitagora
 #INC1=/pitagora/prod/spack/6.1/install/0.22/linux-rhel9-zen4/gcc-11.4.1/nvhpc-25.11-ntshdsgl52b6ckb6iu7xazd6uvqi3wqi/Linux_x86_64/25.11/cuda/13.0/targets/x86_64-linux/include
 #LIB1=/pitagora/prod/spack/6.1/install/0.22/linux-rhel9-zen4/gcc-11.4.1/nvhpc-25.11-ntshdsgl52b6ckb6iu7xazd6uvqi3wqi/Linux_x86_64/25.11/cuda/13.0/targets/x86_64-linux/lib/stubs -lnvidia-ml
endif
#
ifdef AMD
 # AMD GPU target architecture; override from the command line, e.g.
 #   make AMD=1 GFX_ARCH=gfx942 ...
 GFX_ARCH ?= gfx90a
 MULTI = -DMULTICORE
 OFFLOAD = -fopenmp --offload-arch=$(GFX_ARCH) -DOFFLOAD
 OPENACC = not supported
 DC = not supported
 SFC = amdflang -DPGI -DAMDSMI
 SCC = amdclang -DPGI
 FOPT =  -O3
 COPT =  -O3
 OPT =   -O3
 CHECK1=AMD
 # use the AMD SMI wrapper / interface module instead of the NVML ones
 GPU_ENERGY_OBJ = amdsmi_wrapper.o
 GPU_IFACE_MOD = amdsmi_interface_module.o
 # fix for ENERGY/AMD (adjust the ROCm path to the installed version)
 INC1=$(ROCM_PATH)/include
 LIB1=$(ROCM_PATH)/lib -lamd_smi
endif
#
ifdef CRAY
 MULTI =
 OFFLOAD = -DPGI -h omp -DOFFLOAD -DCRAY
 OPENACC = -DPGI -h acc -DNOMANAGED -DOPENACC
 DC = not supported
 SCC = cc
 SFC = ftn
 FOPT = -O3
 COPT =
 OPT = -O3
 CHECK1=CRAY
endif
#
ifdef GNU
# do nothing (default)
endif
#
ifdef INTEL
 MULTI =  -qopenmp -DMULTICORE
 OPENACC = not supported
 OFFLOAD =
 DC =  -DDC
 SFC = ifx
 SCC = icx
 FOPT = -O3 -xCORE-AVX512 -mtune=skylake-avx512 -assume contiguous_pointer
 COPT = -O3
 OPT =  -O3
 OPENMP = -qopenmp
 INFO = -warn all
 CHECK1=INTEL
endif
#
#
# to remove
#ifdef PWR
# SCC = xlc_r
# SFC = xlf_r
# OFFLOAD = -qsmp=omp -qoffload -DOFFLOAD
# FOPT = -O3 -qassert=contig
# COPT = -O3 -qassert=contig
# OPT = -O3 -qassert=contig
# CHECK1=PWR
#endif
#
#----------------------------------------------
# GPU energy/power monitoring (optional, opt-in)
#----------------------------------------------
# Enable with ENERGY=1 on the make command line, e.g.
#   make offload ENERGY=1 AMD=1 FUSED=1 GFX_ARCH=gfx942
# Placed after the vendor blocks so the +=/= appends to the selected vendor.
ENERGY_DEF =
ifdef ENERGY
 ENERGY_DEF = -DENERGY
 FOPT += -DENERGY
 OPT  += -DENERGY
 COPT += -DENERGY
endif
#
#----------------------------------------------
# profiling / timing instrumentation (optional, opt-in)
#----------------------------------------------
# Enable with PROFILING=1 on the make command line, e.g.
#   make cuda NVIDIA=1 FUSED=1 LDC=1 PROFILING=1 NV_ARCH=sm_80
# When OFF (default) all per-phase timers are compiled out (zero overhead)
# and neither bgk.time.log nor bgk.perf are written.
# When ON, the phase timers run and, on the native GPU build, a device
# barrier is issued before each timer stop so the CPU-side timers measure
# real GPU kernel time instead of async launch overhead (serializes
# CPU/GPU -> use for profiling only).
ifdef PROFILING
 FOPT += -DPROFILING
 OPT  += -DPROFILING
endif
#
#----------------------------------------------
# precision used
#----------------------------------------------
#
# default
ifdef SINGLE
 PREC =
 CHECK2=SINGLE
endif
#
ifdef MIXED1
 PREC =-DHALF_P -DMIXEDPRECISION
 CHECK2=MIXED1
endif
#
ifdef MIXED2
 PREC =-DMIXEDPRECISION
 CHECK2=MIXED2
endif
#
ifdef DOUBLE
 PREC =-DDOUBLE_P
 CHECK2=DOUBLE
endif
#
#----------------------------------------------
# implementations
#----------------------------------------------
#
# default, ORIGINAL
ifdef ORIG
 IMPL =
 CHECK3=ORIG
endif
#
ifdef FUSED
 IMPL = -DFUSED
 CHECK3=FUSED
endif
#
ifdef KSYNTAX
 # The OpenMP Offload Kernel Language ("target teams parallel", the
 # omp_get_*_dim runtime API and the num_teams/num_threads dims(N) modifier)
 # needs OpenMP 6.1, so bump the version.
 IMPL += -DOFFLOAD_KERNEL_SYNTAX -fopenmp-version=61
 CHECK3 := $(CHECK3)+KSYNTAX
 # Tag the executable name so the kernel-syntax build is distinguishable
 # (e.g. ../RUN/bgk3d.offload_ksyntax.x).
 KSYNTAX_SUFFIX := _ksyntax
 # Kernel Syntax needs the upstream LLVM flang/clang (module: aomp_ksyntax),
 # which provides the OpenMP Offload Kernel Language; swap the AMD drivers
 # (amdflang/amdclang) for flang/clang while keeping the other flags.
 SFC := $(subst amdflang,flang,$(SFC))
 SCC := $(subst amdclang,clang,$(SCC))
 # Optional 3D thread block (8x8x4) variant of the kernel-syntax loops, for
 # benchmarking against the default 1D-flat (256x1x1) launch: make ... KSYNTAX=1 KS3D=1
 ifdef KS3D
  IMPL += -DKS_BLOCK_3D
  CHECK3 := $(CHECK3)+3D
  KSYNTAX_SUFFIX := _ksyntax_3d
 endif
endif
#
# Build the plain OpenMP-offload kernels with the upstream LLVM flang/clang
# (module: aomp_ksyntax) instead of the AMD amdflang/amdclang drivers, e.g.
# to compare flang vs amdflang performance on the SAME offload code:
#   make offload AMD=1 FLANG=1 FUSED=1 GFX_ARCH=gfx942 LDC=1 ENERGY=1
# (KSYNTAX already performs this driver swap; FLANG lets you do it for the
#  standard "target teams distribute parallel do simd" kernels too.)
ifdef FLANG
 SFC := $(subst amdflang,flang,$(SFC))
 SCC := $(subst amdclang,clang,$(SCC))
 CHECK1 := $(CHECK1)+FLANG
 # Tag the executable so the flang build is distinguishable and does not
 # clobber the amdflang one (e.g. ../RUN/bgk3d.offload_flang.x).
 FLANG_SUFFIX := _flang
endif
#
#
#----------------------------------------------
# test case
#----------------------------------------------
#
# default, Lid Driven Cavity
ifdef LDC
 PREPROC = -DLDC
 CHECK4=LDC
endif
#
# Poiseulle Flow
ifdef POF
 PREPROC = -DCHANNEL -DPOF
 CHECK4=POF
endif
#
# Taylor Green Vortex
ifdef TGV
 PREPROC = -DPERIODIC -DTGV
 CHECK4=TGV
endif
#
# Von Karman Streets  (SPHERE)
ifdef VKS
 PREPROC =-DINFLOW -DOBSTACLE -DDRAG -DVKS
 CHECK4=VKS
endif
#
# Kida Vortex
ifdef KVX
 PREPROC =-DPERIODIC -DKVX
 CHECK4=KVX
endif
#
# Homogenous Isotropic Turbulence  (not forced)
ifdef HIT
 PREPROC =-DPERIODIC -DTGV3D -DHIT
 PREPROC =-DPERIODIC -DTGV3D -DVTK3D -DHIT
 CHECK4=HIT
endif
#
#----------------------------------------------
# module & subroutines
#----------------------------------------------
MOD0 = get_mem_mod.o \

MOD1 = precision.o \
       $(GPU_IFACE_MOD) \
       storage.o \
       bcond_gpu_mod.o \
       mod_tools.o \
       unixtime_module.o \

# 0th level subroutine
OBJ0 = bgk3d.o \

# 1st level subroutine
OBJ1 = setup.o \
       initialize.o \
       diagnostic.o \
       boundaries.o \
       propagation.o \
       save.o\
       finalize.o \

# 2nd level subroutine
OBJ2 = check_case.o \
       input.o \
       alloca.o \
       bcond.driven.o \
       bcond.periodic.o \
       bcond.channel.o \
       bcond.inflow.o \
       bcond.obs.o \
       build.bcond.o \
       build.obs.o \
       restore_raw.o \
       init.o\
       hencol.o \
       outdat.o \
       $(GPU_ENERGY_OBJ) \
       vtk_obstacle.o \
       vtk_3d.binary.o \
       vtk_xy.o \
       vtk_xy.binary.o \
       vtk_xz.binary.o \
       vtk_yz.binary.o \
       vtk_om.xy.binary.o \
       varm.o \
       prof_k.o \
       prof_j.o \
       prof_i.o \
       diagno.o \
       dissipation.o \
       dissipation2.o \
       probe.o \
       movef.o\
       save_raw.o\
       save_vel.o\
       profile.o\
       git_info.o \
       make_info.o \
       draglift.o \
       dealloca.o

# 3rd level subroutine
OBJ3 = get_mem.o \
       collision.o \
       col_MC.o \
       col.o

OBJ = $(OBJ1) $(OBJ0) $(OBJ2) $(OBJ3)

#----------------------------------------------
# targets
#----------------------------------------------
#
#
serial: all
#
#----------------------------------------------
all: VER := serial
all: FC := $(SFC)
all: CC := $(SCC)
all: $(MOD0) $(MOD1) $(OBJ) $(OBJ4)
	$(FC) $(OPT) $(MOD1) $(OBJ) $(OBJ4) $(IMPL) $(LIB) $(FIX) $(if $(LIB1),-L$(LIB1)) -o ../RUN/bgk3d.$(VER).x
	@echo "-------------------------------------------------"
	@echo "options: "$(CHECK1) $(CHECK2) $(CHECK3) $(CHECK4)
	@echo "further option: "$(FIX)
	@echo $(FINE)
	@echo "-------------------------------------------------"
	rm make_info.F90
	rm git_info.F90
#
#----------------------------------------------
doconcurrent: VER := doconcurrent
doconcurrent: FOPT += $(DC)
doconcurrent: OPT  += $(DC)
doconcurrent: FC := $(SFC)
doconcurrent: CC := $(SCC)
doconcurrent: $(MOD0) $(MOD1) $(OBJ) $(OBJ4)
	$(FC) $(OPT) $(MOD1) $(OBJ) $(OBJ4) $(IMPL) $(LIB) $(FIX) -L$(LIB1) -o ../RUN/bgk3d.$(VER).x
	@echo "-------------------------------------------------"
	@echo "options: "$(CHECK1) $(CHECK2) $(CHECK3) $(CHECK4)
	@echo "further option: "$(FIX)
	@echo $(FINE)
	@echo "-------------------------------------------------"
	rm make_info.F90
	rm git_info.F90
#
#----------------------------------------------
offload: FOPT += $(OFFLOAD)
offload: OPT  += $(OFFLOAD)
offload: VER := offload$(KSYNTAX_SUFFIX)$(FLANG_SUFFIX)
offload: FC := $(SFC)
offload: CC := $(SCC)
offload:$(MOD0) $(MOD1) $(OBJ) $(OBJ4)
	$(FC) $(OPT) $(MOD1) $(OBJ) $(OBJ4) $(IMPL) $(FIX) $(LIB) -L$(LIB1) -o ../RUN/bgk3d.$(VER).x
	@echo "-------------------------------------------------"
	@echo "options: "$(CHECK1) $(CHECK2) $(CHECK3) $(CHECK4)
	@echo "further option: "$(FIX)
	@echo $(FINE)
	@echo "-------------------------------------------------"
	rm make_info.F90
	rm git_info.F90
#
#----------------------------------------------
# Native GPU collision via HIP.
# Data residency + boundary kernels stay on OpenMP offload; only the
# fused collision (col_MC) is replaced by a HIP kernel.
# Typical use (AMD):
#   make hip AMD=1 FUSED=1 LDC=1 GFX_ARCH=gfx942
# GPU power/energy monitoring (AMD SMI) is backend-agnostic: just add
# ENERGY=1 exactly as for the offload/openacc targets, e.g.
#   make hip AMD=1 FUSED=1 LDC=1 ENERGY=1 GFX_ARCH=gfx942
# Results (energy in J, mean power in W) are written to RUN/bgk.perf.
hip: FOPT += $(OFFLOAD) -DGPU_NATIVE
hip: OPT  += $(OFFLOAD) -DGPU_NATIVE
hip: VER := hip
hip: FC := $(SFC)
hip: CC := $(SCC)
hip: $(MOD0) $(MOD1) $(OBJ) $(OBJ4) col_MC_hip.o bcond_hip.o
	$(FC) $(OPT) $(MOD1) $(OBJ) $(OBJ4) col_MC_hip.o bcond_hip.o $(IMPL) $(FIX) $(LIB) -L$(LIB1) -L$(ROCM_PATH)/lib -lamdhip64 -lstdc++ -o ../RUN/bgk3d.$(VER).x
	@echo "-------------------------------------------------"
	@echo "options: "$(CHECK1) $(CHECK2) $(CHECK3) $(CHECK4)" HIP/"$(GFX_ARCH)
	@echo "further option: "$(FIX)
	@echo $(FINE)
	@echo "-------------------------------------------------"
	rm make_info.F90
	rm git_info.F90
#
#----------------------------------------------
# Native GPU collision via CUDA (NVIDIA). Mirror of the HIP target.
# On NVIDIA, data residency + boundary kernels are handled by OpenACC
# (managed memory, like the `openacc` target); only the fused collision
# (col_MC) is replaced by a native CUDA kernel.  Device pointers are
# handed to the CUDA launcher via `!$acc host_data use_device` (see
# col_MC.F90), so OpenMP offload is NOT used on the NVIDIA path.
# Typical use (NVIDIA):
#   make cuda NVIDIA=1 FUSED=1 LDC=1 NV_ARCH=sm_80
# GPU power/energy monitoring (NVML) is backend-agnostic: just add
# ENERGY=1 exactly as for the offload/openacc targets, e.g.
#   make cuda NVIDIA=1 FUSED=1 LDC=1 ENERGY=1 NV_ARCH=sm_80
# Results (energy in J, mean power in W) are written to RUN/bgk.perf.
cuda: FOPT += $(OPENACC) -DGPU_NATIVE
cuda: OPT  += $(OPENACC) -DGPU_NATIVE
cuda: VER := cuda
cuda: FC := $(SFC)
cuda: CC := $(SCC)
cuda: $(MOD0) $(MOD1) $(OBJ) $(OBJ4) col_MC_cuda.o bcond_cuda.o
	$(FC) $(OPT) $(MOD1) $(OBJ) $(OBJ4) col_MC_cuda.o bcond_cuda.o $(IMPL) $(FIX) $(LIB) -L$(LIB1) -L$(CUDA_PATH)/lib64 -lcudart -lstdc++ -o ../RUN/bgk3d.$(VER).x
	@echo "-------------------------------------------------"
	@echo "options: "$(CHECK1) $(CHECK2) $(CHECK3) $(CHECK4)" CUDA/"$(NV_ARCH)
	@echo "further option: "$(FIX)
	@echo $(FINE)
	@echo "-------------------------------------------------"
	rm make_info.F90
	rm git_info.F90
#
#----------------------------------------------
openacc: FOPT += $(OPENACC)
openacc: OPT  += $(OPENACC)
openacc: VER := openacc
openacc: FC := $(SFC)
openacc: CC := $(SCC)
openacc:$(MOD0) $(MOD1) $(OBJ) $(OBJ4)
	$(FC) $(OPT) $(MOD1) $(OBJ) $(OBJ4) $(IMPL) $(FIX) $(LIB) -L$(LIB1) -o ../RUN/bgk3d.$(VER).x
	@echo "-------------------------------------------------"
	@echo "options: "$(CHECK1) $(CHECK2) $(CHECK3) $(CHECK4)
	@echo "further option: "$(FIX)
	@echo $(FINE)
	@echo "-------------------------------------------------"
	rm make_info.F90
	rm git_info.F90
#
#----------------------------------------------
multicore: VER:= multicore
multicore: FC := $(SFC) $(MULTI)
multicore: CC := $(SCC)
multicore: $(MOD0) $(MOD1) $(OBJ) $(OBJ4)
	$(FC) $(OPT) $(MOD1) $(OBJ) $(OBJ4) $(IMPL) $(FIX) $(LIB) -o ../RUN/bgk3d.$(VER).x
	@echo "-------------------------------------------------"
	@echo "options: "$(CHECK1) $(CHECK2) $(CHECK3) $(CHECK4)
	@echo "further option: "$(FIX)
	@echo $(FINE)
	@echo "-------------------------------------------------"
	rm make_info.F90
	rm git_info.F90
#
#----------------------------------------------
# rules
#----------------------------------------------

nvml_wrapper.o: nvml_wrapper.c
	$(SCC) -O3 -I$(INC1) $(ENERGY_DEF) $(FIX) -c nvml_wrapper.c

amdsmi_wrapper.o: amdsmi_wrapper.c
	$(SCC) -O3 -I$(INC1) $(ENERGY_DEF) $(FIX) -c amdsmi_wrapper.c

# Native GPU collision back-ends (kernel body shared via col_MC_gpu.hpp).
# $(PREC) carries -DDOUBLE_P etc. so real_t matches the Fortran storage kind.
col_MC_hip.o: col_MC_hip.cpp col_MC_gpu.hpp
	$(HIPCC) -O3 --offload-arch=$(GFX_ARCH) $(PREC) $(FIX) -c col_MC_hip.cpp -o col_MC_hip.o

col_MC_cuda.o: col_MC_cuda.cu col_MC_gpu.hpp
	$(NVCC) -O3 -arch=$(NV_ARCH) $(PREC) $(FIX) -c col_MC_cuda.cu -o col_MC_cuda.o

# Native GPU boundary-condition back-ends (kernels shared via bcond_gpu.hpp).
bcond_hip.o: bcond_hip.cpp bcond_gpu.hpp
	$(HIPCC) -O3 --offload-arch=$(GFX_ARCH) $(PREC) $(FIX) -c bcond_hip.cpp -o bcond_hip.o

bcond_cuda.o: bcond_cuda.cu bcond_gpu.hpp
	$(NVCC) -O3 -arch=$(NV_ARCH) $(PREC) $(FIX) -c bcond_cuda.cu -o bcond_cuda.o

%.o %.mod: $(INC) %.f90
	$(FC) $(FOPT) $(PREC) $(PREPROC) $(IMPL) $(FIX) -c $<

%.o: $(INC) %.F90
	$(FC) $(FOPT) $(PREC) $(PREPROC) $(IMPL) $(FIX) -c $<

%.o: $(INC) %.f
	$(FC) $(FOPT) $(PREC) $(PREPROC) $(IMPL) $(FIX) -c $<

%.o: $(INC) %.c
	$(CC) $(COPT)                          -c $<

make_info.F90:
	echo "    subroutine make_info"                             >  make_info.F90
	echo "    write( 6,*) \""INFO: `echo $(CHECK1) ` \"         >> make_info.F90
	echo "    write( 6,*) \""INFO: `echo $(CHECK2) ` \"         >> make_info.F90
	echo "    write( 6,*) \""INFO: `echo $(CHECK3) ` \"         >> make_info.F90
	echo "    write( 6,*) \""INFO: `echo $(CHECK4) ` \"         >> make_info.F90
	echo "    write( 6,*) \""INFO: `echo $(FIX)    ` \"         >> make_info.F90
	echo "    write( 6,*) \""INFO: `echo $(FC)     ` \"         >> make_info.F90
	echo "    write( 6,*) \""INFO: `echo $(FOPT)   ` \"         >> make_info.F90
	echo "    write( 6,*) \""INFO: `echo $(PREPROC)` \"         >> make_info.F90
	echo "    end subroutine make_info"                         >> make_info.F90

git_info.F90:
	echo "    subroutine git_info"                              >> git_info.F90
	echo "    write( 6,*) \""INFO: `git log  | head  -n 1` \"   >> git_info.F90
	echo "    write(16,*) \""INFO: `git log  | head  -n 1` \"   >> git_info.F90
	echo "    end subroutine git_info"                          >> git_info.F90
#
#----------------------------------------------
# utilities
#----------------------------------------------
clean:
	rm -f *.o;
	rm -f *.mod;
	rm -f core;
#----------------------------------------------
info:  $(INC)
	@echo "Modules          =  "$(MOD0) $(MOD1);
	@echo "Objects          =  "$(OBJ);
	@echo "Compiler         =  "$(SFC);
	@echo "Compiler  flags  =  "$(OPT);
	@echo "Preproc   flags  =  "$(PREPROC);
	@echo "Fix       flags  =  "$(FIX);
	@echo "Precision flags  =  "$(CHECK2);
	@echo "Implemet. flags  =  "$(CHECK3);
	@echo "Version          =  "$(CHECK4);

#----------------------------------------------
test:
	make clean;
	make serial;
	make clean;
	make multicore;
	make clean;
	make doconcurrent;
	make clean;
	make openacc;
	make clean;
	make offload;
#----------------------------------------------
tar:
	make clean
	rm -rf ../RUN/*.vtk
	rm -rf ../RUN/*.bin
	rm -rf ../RUN/*.h5
#	cd ../TEST/
#	./clean.x
#	cd -
	tar -cvf ../../archive.bgk3d.$(DATE).tar ../../BGK_3D/*
#----------------------------------------------
help:
	@echo "different target available"
	@echo "do concurrent (default, GPU) "
	@echo "serial (CPU)"
	@echo "multicore (CPU)"
	@echo "offload (GPU)"
	@echo "openacc (GPU)"
#----------------------------------------------

