// ====================================================================
//  bcond_hip.cpp
//
//  HIP back-end for the native GPU port of the D3Q19 boundary kernels.
//  Compiled with hipcc and -DUSE_HIP.  All kernels, launchers and the
//  extern "C" entry points live in the shared header bcond_gpu.hpp so
//  the CUDA back-end (bcond_cuda.cu) reuses them verbatim.
//
//  The Fortran caller (bcond.*.F90, GPU_NATIVE path) passes the host
//  array pa(1:19) of DEVICE pointers (obtained via c_loc under the
//  OpenMP `target data use_device_addr` residency model), plus the
//  scalar parameters by value.
// ====================================================================
#define USE_HIP
#include "bcond_gpu.hpp"
