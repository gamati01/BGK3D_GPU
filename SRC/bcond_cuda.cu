// ====================================================================
//  bcond_cuda.cu
//
//  CUDA back-end for the native GPU port of the D3Q19 boundary kernels.
//  Compiled with nvcc and -DUSE_CUDA.  All kernels, launchers and the
//  extern "C" entry points live in the shared header bcond_gpu.hpp
//  (identical to the HIP back-end).
//
//  The Fortran caller (bcond.*.F90, GPU_NATIVE path) passes the host
//  array pa(1:19) of DEVICE pointers (obtained via c_loc inside an
//  OpenACC `host_data use_device` region, or directly under managed
//  memory), plus the scalar parameters by value.
// ====================================================================
#define USE_CUDA
#include "bcond_gpu.hpp"
