// ====================================================================
//  col_MC_gpu.hpp
//
//  Native GPU (HIP / CUDA) port of the fused move+collide D3Q19 BGK
//  kernel originally implemented in SRC/col_MC.F90 (FUSED path).
//
//  This header is shared by both back-ends:
//    - SRC/col_MC_hip.cpp  (compiled with hipcc, -DUSE_HIP)
//    - SRC/col_MC_cuda.cu  (compiled with nvcc,   -DUSE_CUDA)
//
//  The kernel reproduces the math of the Fortran collapse(3) loop
//  exactly.  Memory is NOT owned here: the caller (col_MC.F90) hands us
//  device pointers for the storage arrays.  Those addresses come from
//  the GPU data environment that also drives the boundary kernels:
//    - CUDA / NVIDIA : OpenACC  (`!$acc host_data use_device`, managed mem)
//    - HIP  / AMD    : OpenMP   (`!$omp target data use_device_addr`)
//
//  Currently covered: default (LDC) config, optional uniform forcing
//  (forcex/forcey/forcez), NOSHIFT handled through cte0/cte1.
//  NOT yet covered: LES (per-cell omega) and TGV per-cell forcing.
// ====================================================================
#ifndef COL_MC_GPU_HPP
#define COL_MC_GPU_HPP

#if defined(USE_HIP)
#  include <hip/hip_runtime.h>
#  define GPU_CHECK(call)  do {                                            \
       hipError_t _e = (call);                                             \
       if (_e != hipSuccess) {                                             \
           fprintf(stderr, "HIP error %s at %s:%d\n",                      \
                   hipGetErrorString(_e), __FILE__, __LINE__);             \
       }                                                                   \
   } while (0)
#  define GPU_DEVICE_SYNC() hipDeviceSynchronize()
#elif defined(USE_CUDA)
#  include <cuda_runtime.h>
#  define GPU_CHECK(call)  do {                                            \
       cudaError_t _e = (call);                                            \
       if (_e != cudaSuccess) {                                            \
           fprintf(stderr, "CUDA error %s at %s:%d\n",                     \
                   cudaGetErrorString(_e), __FILE__, __LINE__);            \
       }                                                                   \
   } while (0)
#  define GPU_DEVICE_SYNC() cudaDeviceSynchronize()
#else
#  error "Define USE_HIP or USE_CUDA when compiling col_MC_gpu.hpp"
#endif

#include <cstdio>
#include <cstdint>

// --------------------------------------------------------------------
//  Storage precision: must match `mystorage` in SRC/precision.F90.
//  Default build is SINGLE -> float.  Pass -DGPU_DOUBLE for DOUBLE_P.
// --------------------------------------------------------------------
#if defined(DOUBLE_P) || defined(GPU_DOUBLE)
typedef double real_t;
#else
typedef float  real_t;
#endif

// Scalar parameters passed by value from the Fortran side.  Keeping
// them in one struct avoids a 30-argument launcher signature.
struct ColParams {
    int    l, m, n;        // interior extents (loops run 1..l, 1..m, 1..n)
    int    nx, ny;         // leading dimensions = size(a,1), size(a,2)
    real_t omega;
    real_t cte0, cte1;     // NOSHIFT control (cte1 = 0 if NOSHIFT else 1)
    real_t p0, p1, p2;     // equilibrium weights
    real_t rf, qf, tre;    // linear / quadratic / "3" factors
    real_t forcex, forcey, forcez;
};

// Bundles of the 19 input (a*) and 18+1 output (b*/a19) device pointers.
// Order matches direction index 01..19.
struct ColArrays {
    const real_t* a[19];   // a01..a19  (read)
    real_t*       b[19];   // b01..b18 written; b[18] aliases a19 (in-place)
};

// --------------------------------------------------------------------
//  Device kernel: one thread per interior cell.
// --------------------------------------------------------------------
__global__ void col_mc_kernel(ColArrays arr, ColParams p)
{
    const long tid = (long)blockIdx.x * blockDim.x + threadIdx.x;
    const long ncell = (long)p.l * p.m * p.n;
    if (tid >= ncell) return;

    // Map linear id -> (i,j,k) in 1..l, 1..m, 1..n
    const int i = (int)(tid % p.l) + 1;
    const int j = (int)((tid / p.l) % p.m) + 1;
    const int k = (int)(tid / ((long)p.l * p.m)) + 1;

    const long nx  = p.nx;
    const long nxy = (long)p.nx * p.ny;
    const long base = (long)i + nx * (long)j + nxy * (long)k;

    // gather the 19 streamed populations (matches Fortran offsets)
    const real_t x01 = arr.a[0][base - 1 + nx];        // (i-1,j+1,k  )
    const real_t x02 = arr.a[1][base - 1 + nxy];       // (i-1,j  ,k+1)
    const real_t x03 = arr.a[2][base - 1 - nx];        // (i-1,j-1,k  )
    const real_t x04 = arr.a[3][base - 1 - nxy];       // (i-1,j  ,k-1)
    const real_t x05 = arr.a[4][base - 1];             // (i-1,j  ,k  )
    const real_t x06 = arr.a[5][base - nxy];           // (i  ,j  ,k-1)
    const real_t x07 = arr.a[6][base - nx - nxy];      // (i  ,j-1,k-1)
    const real_t x08 = arr.a[7][base - nx];            // (i  ,j-1,k  )
    const real_t x09 = arr.a[8][base - nx + nxy];      // (i  ,j-1,k+1)
    const real_t x10 = arr.a[9][base + 1 + nx];        // (i+1,j+1,k  )
    const real_t x11 = arr.a[10][base + 1 + nxy];      // (i+1,j  ,k+1)
    const real_t x12 = arr.a[11][base + 1 - nx];       // (i+1,j-1,k  )
    const real_t x13 = arr.a[12][base + 1 - nxy];      // (i+1,j  ,k-1)
    const real_t x14 = arr.a[13][base + 1];            // (i+1,j  ,k  )
    const real_t x15 = arr.a[14][base + nxy];          // (i  ,j  ,k+1)
    const real_t x16 = arr.a[15][base + nx + nxy];     // (i  ,j+1,k+1)
    const real_t x17 = arr.a[16][base + nx];           // (i  ,j+1,k  )
    const real_t x18 = arr.a[17][base + nx - nxy];     // (i  ,j+1,k-1)
    const real_t x19 = arr.a[18][base];                // (i  ,j  ,k  )

    const real_t cte0 = p.cte0;
    const real_t cte1 = p.cte1;
    const real_t tre  = p.tre;
    const real_t qf   = p.qf;
    const real_t rf   = p.rf;

    real_t rho = x01 + x02 + x03 + x04 + x05 + x06 + x07 + x08
               + x09 + x10 + x11 + x12 + x13 + x14 + x15 + x16
               + x17 + x18 + x19 + cte1;

    const real_t rhoinv = (real_t)1 / rho;

    real_t vx = (x01 + x02 + x03 + x04 + x05 - x10 - x11 - x12 - x13 - x14) * rhoinv;
    real_t vy = (x03 + x07 + x08 + x09 + x12 - x01 - x10 - x16 - x17 - x18) * rhoinv;
    real_t vz = (x04 + x06 + x07 + x13 + x18 - x02 - x09 - x11 - x15 - x16) * rhoinv;

    const real_t vx2 = vx * vx;
    const real_t vy2 = vy * vy;
    const real_t vz2 = vz * vz;
    const real_t vsq = vx2 + vy2 + vz2;

    real_t vxpy = vx + vy;
    const real_t qxpy = cte0 + qf * (tre * vxpy * vxpy - vsq);
    real_t vxmy = vx - vy;
    const real_t qxmy = cte0 + qf * (tre * vxmy * vxmy - vsq);
    real_t vxpz = vx + vz;
    const real_t qxpz = cte0 + qf * (tre * vxpz * vxpz - vsq);
    real_t vxmz = vx - vz;
    const real_t qxmz = cte0 + qf * (tre * vxmz * vxmz - vsq);
    real_t vypz = vy + vz;
    const real_t qypz = cte0 + qf * (tre * vypz * vypz - vsq);
    real_t vymz = vy - vz;
    const real_t qymz = cte0 + qf * (tre * vymz * vymz - vsq);
    const real_t qx = cte0 + qf * (tre * vx2 - vsq);
    const real_t qy = cte0 + qf * (tre * vy2 - vsq);
    const real_t qz = cte0 + qf * (tre * vz2 - vsq);
    const real_t q0 = cte0 + qf * (-vsq);

    vx   = rf * vx;
    vy   = rf * vy;
    vz   = rf * vz;
    vxpy = rf * vxpy;
    vxmy = rf * vxmy;
    vxpz = rf * vxpz;
    vxmz = rf * vxmz;
    vypz = rf * vypz;
    vymz = rf * vymz;

    const real_t rp0 = rho * p.p0;
    const real_t rp1 = rho * p.p1;
    const real_t rp2 = rho * p.p2;
    const real_t P0  = p.p0;
    const real_t P1  = p.p1;
    const real_t P2  = p.p2;

    const real_t e01 = rp2 * (+vxmy + qxmy) + cte1 * (rp2 - P2);
    const real_t e02 = rp2 * (+vxmz + qxmz) + cte1 * (rp2 - P2);
    const real_t e03 = rp2 * (+vxpy + qxpy) + cte1 * (rp2 - P2);
    const real_t e04 = rp2 * (+vxpz + qxpz) + cte1 * (rp2 - P2);
    const real_t e05 = rp1 * (+vx   + qx  ) + cte1 * (rp1 - P1);
    const real_t e06 = rp1 * (+vz   + qz  ) + cte1 * (rp1 - P1);
    const real_t e07 = rp2 * (+vypz + qypz) + cte1 * (rp2 - P2);
    const real_t e08 = rp1 * (+vy   + qy  ) + cte1 * (rp1 - P1);
    const real_t e09 = rp2 * (+vymz + qymz) + cte1 * (rp2 - P2);
    const real_t e10 = rp2 * (-vxpy + qxpy) + cte1 * (rp2 - P2);
    const real_t e11 = rp2 * (-vxpz + qxpz) + cte1 * (rp2 - P2);
    const real_t e12 = rp2 * (-vxmy + qxmy) + cte1 * (rp2 - P2);
    const real_t e13 = rp2 * (-vxmz + qxmz) + cte1 * (rp2 - P2);
    const real_t e14 = rp1 * (-vx   + qx  ) + cte1 * (rp1 - P1);
    const real_t e15 = rp1 * (-vz   + qz  ) + cte1 * (rp1 - P1);
    const real_t e16 = rp2 * (-vypz + qypz) + cte1 * (rp2 - P2);
    const real_t e17 = rp1 * (-vy   + qy  ) + cte1 * (rp1 - P1);
    const real_t e18 = rp2 * (-vymz + qymz) + cte1 * (rp2 - P2);
    const real_t e19 = rp0 * (       + q0 ) + cte1 * (rp0 - P0);

    const real_t omega  = p.omega;
    const real_t forcex = p.forcex;
    const real_t forcey = p.forcey;
    const real_t forcez = p.forcez;

    arr.b[0][base]  = x01 - omega * (x01 - e01) + forcex - forcey;
    arr.b[1][base]  = x02 - omega * (x02 - e02) + forcex          - forcez;
    arr.b[2][base]  = x03 - omega * (x03 - e03) + forcex + forcey;
    arr.b[3][base]  = x04 - omega * (x04 - e04) + forcex          + forcez;
    arr.b[4][base]  = x05 - omega * (x05 - e05) + forcex;
    arr.b[5][base]  = x06 - omega * (x06 - e06)                   + forcez;
    arr.b[6][base]  = x07 - omega * (x07 - e07)          + forcey + forcez;
    arr.b[7][base]  = x08 - omega * (x08 - e08)          + forcey;
    arr.b[8][base]  = x09 - omega * (x09 - e09)          + forcey - forcez;
    arr.b[9][base]  = x10 - omega * (x10 - e10) - forcex - forcey;
    arr.b[10][base] = x11 - omega * (x11 - e11) - forcex          - forcez;
    arr.b[11][base] = x12 - omega * (x12 - e12) - forcex + forcey;
    arr.b[12][base] = x13 - omega * (x13 - e13) - forcex          + forcez;
    arr.b[13][base] = x14 - omega * (x14 - e14) - forcex;
    arr.b[14][base] = x15 - omega * (x15 - e15)                   - forcez;
    arr.b[15][base] = x16 - omega * (x16 - e16)          - forcey - forcez;
    arr.b[16][base] = x17 - omega * (x17 - e17)          - forcey;
    arr.b[17][base] = x18 - omega * (x18 - e18)          - forcey + forcez;
    arr.b[18][base] = x19 - omega * (x19 - e19);   // b[18] aliases a19 (in place)
}

// --------------------------------------------------------------------
//  Host launcher shared by both back-ends.
// --------------------------------------------------------------------
static inline void col_mc_gpu_launch(const ColArrays& arr, const ColParams& p)
{
    const long ncell = (long)p.l * p.m * p.n;
    const int  block = 256;
    const int  grid  = (int)((ncell + block - 1) / block);

#if defined(USE_HIP)
    hipLaunchKernelGGL(col_mc_kernel, dim3(grid), dim3(block), 0, 0, arr, p);
#else
    col_mc_kernel<<<grid, block>>>(arr, p);
#endif

    // The boundary kernels (next iteration) run on the OpenACC (NVIDIA) or
    // OpenMP-offload (AMD) queue, which is distinct from this default GPU
    // stream.  Synchronize so the collision results are visible before they
    // execute, and so the Fortran-side collision timer measures the real
    // kernel cost.
    GPU_CHECK(GPU_DEVICE_SYNC());
}

#endif // COL_MC_GPU_HPP
