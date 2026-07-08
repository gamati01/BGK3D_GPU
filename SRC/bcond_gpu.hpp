// ====================================================================
//  bcond_gpu.hpp
//
//  Native GPU (HIP / CUDA) port of the D3Q19 boundary-condition kernels
//  originally implemented in SRC/bcond.*.F90 (OFFLOAD / OPENACC paths).
//
//  Shared by both back-ends:
//    - SRC/bcond_hip.cpp  (compiled with hipcc, -DUSE_HIP)
//    - SRC/bcond_cuda.cu  (compiled with nvcc,  -DUSE_CUDA)
//
//  These kernels reproduce the math of the Fortran boundary loops
//  EXACTLY, following the OFFLOAD/OPENACC branch conventions (which are
//  the ones that run on the GPU).  Memory is NOT owned here: the caller
//  (bcond.*.F90, GPU_NATIVE path) hands us device pointers for the a*
//  storage arrays (and the obstacle mask for bcond_obs).
//
//  Storage layout (see SRC/alloca.F90):
//    a01..a19 allocated as (0:l+1+ipad, 0:m+1+jpad, 0:n+1+kpad)
//    -> lower bound 0, leading dims nx = size(a,1), ny = size(a,2).
//    Element (i,j,k) lives at linear offset  i + nx*j + nxy*k  from the
//    base pointer c_loc(a(0,0,0)) -- identical convention to col_MC.
//    The obstacle mask obs is (1:l,1:m,1:n) -> lower bound 1, leading
//    dims obs_nx = l, obs_ny = m; element (i,j,k) at
//    (i-1) + obs_nx*(j-1) + obs_nx*obs_ny*(k-1).
//
//  All launches use the default stream (stream 0), the SAME stream as
//  the native col_MC kernel, so per-step cross-stream synchronization is
//  unnecessary: kernels serialize in issue order.  A device sync is only
//  required where the HOST reads device data (diagno / dissipation /
//  draglift / save), handled by gpu_device_sync() on the Fortran side.
// ====================================================================
#ifndef BCOND_GPU_HPP
#define BCOND_GPU_HPP

#if defined(USE_HIP)
#  include <hip/hip_runtime.h>
#elif defined(USE_CUDA)
#  include <cuda_runtime.h>
#else
#  error "Define USE_HIP or USE_CUDA when compiling bcond_gpu.hpp"
#endif

#include <cstdio>
#include <cstdint>

// Storage precision: must match `mystorage` in SRC/precision.F90.
#if defined(DOUBLE_P) || defined(GPU_DOUBLE)
typedef double real_t;
#else
typedef float  real_t;
#endif

// Scalar parameters shared by the boundary kernels (superset; each
// kernel uses only the fields it needs).
struct BcParams {
    int  l, m, n;          // interior extents
    int  l1, m1, n1;       // l+1, m+1, n+1
    long nx, nxy;          // leading dim and nx*ny (element strides)
    // obstacle mask geometry / loop bounds (bcond_obs)
    int  imin, imax, jmin, jmax, kmin, kmax;
    long obs_nx, obs_nxy;
    // lid-driven forcing (bcond_driven)
    real_t force;
    // inflow equilibrium (bcond_inflow)
    real_t u_inflow;
    real_t rf, qf, tre;
    real_t p0, p1, p2;
    real_t cte1;
};

// 19 storage-array device pointers, index 0..18 == a01..a19.
struct BcArrays {
    real_t* a[19];
};

// --------------------------------------------------------------------
//  Linear index helper (matches the Fortran column-major layout).
// --------------------------------------------------------------------
#define AIDX(i,j,k)  ((long)(i) + P.nx*(long)(j) + P.nxy*(long)(k))

// convenience: A.a[label-1] is aNN
#define AA(label,i,j,k)  A.a[(label)-1][AIDX(i,j,k)]

// ====================================================================
//  Launch helpers
// ====================================================================
static inline void gpu_launch_1d(void (*)(BcArrays, BcParams)) {}  // (unused placeholder)

static inline int grid_for(long count, int block) {
    return (int)((count + block - 1) / block);
}

// Non-blocking check of the last kernel LAUNCH error.  hipGetLastError /
// cudaGetLastError do NOT synchronize the device, so this is free to call
// after every launch and does not perturb async timing.  Execution-time
// faults still surface later at gpu_device_sync().  Flag failures loudly so
// a silently-skipped kernel can never masquerade as a fast, valid run.
#if defined(USE_HIP)
#  define GPU_CHECK_LAUNCH(name)  do {                                        \
       hipError_t _le = hipGetLastError();                                    \
       if (_le != hipSuccess) {                                               \
           fprintf(stderr, "ERROR: HIP launch failed for %s: %s at %s:%d\n",  \
                   (name), hipGetErrorString(_le), __FILE__, __LINE__);       \
           fflush(stderr);                                                    \
       }                                                                      \
   } while (0)
#else
#  define GPU_CHECK_LAUNCH(name)  do {                                        \
       cudaError_t _le = cudaGetLastError();                                  \
       if (_le != cudaSuccess) {                                             \
           fprintf(stderr, "ERROR: CUDA launch failed for %s: %s at %s:%d\n", \
                   (name), cudaGetErrorString(_le), __FILE__, __LINE__);      \
           fflush(stderr);                                                    \
       }                                                                      \
   } while (0)
#endif

// Macro to launch a 1D kernel over `count` cells on the default stream.
// Each launch is followed by a non-blocking launch-error check.
#if defined(USE_HIP)
#  define LAUNCH1D(KERN, COUNT, ARR, PAR)  do {                            \
     hipLaunchKernelGGL(KERN, dim3(grid_for((COUNT), 256)), dim3(256), 0, 0, (ARR), (PAR)); \
     GPU_CHECK_LAUNCH(#KERN);                                              \
   } while (0)
#else
#  define LAUNCH1D(KERN, COUNT, ARR, PAR)  do {                            \
     KERN<<<grid_for((COUNT), 256), 256>>>((ARR), (PAR));                  \
     GPU_CHECK_LAUNCH(#KERN);                                              \
   } while (0)
#endif

// ====================================================================
//  PERIODIC boundary  (bcond.periodic.F90)
// ====================================================================
// loop 1: x-faces, over j=0..m+1, k=0..n+1
__global__ void bc_periodic_x(BcArrays A, BcParams P)
{
    const long naj = (long)P.m + 2;      // j in 0..m+1
    const long nak = (long)P.n + 2;      // k in 0..n+1
    const long t = (long)blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= naj * nak) return;
    const int j = (int)(t % naj);
    const int k = (int)(t / naj);
    // front (x = l)
    AA(1,0,j,k) = AA(1,P.l,j,k);
    AA(2,0,j,k) = AA(2,P.l,j,k);
    AA(3,0,j,k) = AA(3,P.l,j,k);
    AA(4,0,j,k) = AA(4,P.l,j,k);
    AA(5,0,j,k) = AA(5,P.l,j,k);
    // rear (x = 1)
    AA(10,P.l1,j,k) = AA(10,1,j,k);
    AA(11,P.l1,j,k) = AA(11,1,j,k);
    AA(12,P.l1,j,k) = AA(12,1,j,k);
    AA(13,P.l1,j,k) = AA(13,1,j,k);
    AA(14,P.l1,j,k) = AA(14,1,j,k);
}

// loop 2: y-faces, over i=0..l+1, k=0..n+1
__global__ void bc_periodic_y(BcArrays A, BcParams P)
{
    const long nai = (long)P.l + 2;
    const long nak = (long)P.n + 2;
    const long t = (long)blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= nai * nak) return;
    const int i = (int)(t % nai);
    const int k = (int)(t / nai);
    // left (y = 1)
    AA(3,i,0,k) = AA(3,i,P.m,k);
    AA(7,i,0,k) = AA(7,i,P.m,k);
    AA(8,i,0,k) = AA(8,i,P.m,k);
    AA(9,i,0,k) = AA(9,i,P.m,k);
    AA(12,i,0,k) = AA(12,i,P.m,k);
    // right (y = m)
    AA(1,i,P.m1,k) = AA(1,i,1,k);
    AA(10,i,P.m1,k) = AA(10,i,1,k);
    AA(16,i,P.m1,k) = AA(16,i,1,k);
    AA(17,i,P.m1,k) = AA(17,i,1,k);
    AA(18,i,P.m1,k) = AA(18,i,1,k);
}

// loop 3: z-faces, over i=0..l+1, j=0..m+1
__global__ void bc_periodic_z(BcArrays A, BcParams P)
{
    const long nai = (long)P.l + 2;
    const long naj = (long)P.m + 2;
    const long t = (long)blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= nai * naj) return;
    const int i = (int)(t % nai);
    const int j = (int)(t / nai);
    // bottom (z = 1)
    AA(4,i,j,0) = AA(4,i,j,P.n);
    AA(6,i,j,0) = AA(6,i,j,P.n);
    AA(7,i,j,0) = AA(7,i,j,P.n);
    AA(13,i,j,0) = AA(13,i,j,P.n);
    AA(18,i,j,0) = AA(18,i,j,P.n);
    // up (z = n)
    AA(2,i,j,P.n1) = AA(2,i,j,1);
    AA(9,i,j,P.n1) = AA(9,i,j,1);
    AA(11,i,j,P.n1) = AA(11,i,j,1);
    AA(15,i,j,P.n1) = AA(15,i,j,1);
    AA(16,i,j,P.n1) = AA(16,i,j,1);
}

static inline void bcond_periodic_gpu_launch(const BcArrays& arr, const BcParams& p)
{
    LAUNCH1D(bc_periodic_x, ((long)p.m + 2) * ((long)p.n + 2), arr, p);
    LAUNCH1D(bc_periodic_y, ((long)p.l + 2) * ((long)p.n + 2), arr, p);
    LAUNCH1D(bc_periodic_z, ((long)p.l + 2) * ((long)p.m + 2), arr, p);
}

// ====================================================================
//  CHANNEL boundary  (bcond.channel.F90)
//    x/y periodic (full faces), z no-slip bounce-back (interior 1..l,1..m)
// ====================================================================
// loop 1: x-faces, over j=0..m+1, k=0..n+1
__global__ void bc_channel_x(BcArrays A, BcParams P)
{
    const long naj = (long)P.m + 2;
    const long nak = (long)P.n + 2;
    const long t = (long)blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= naj * nak) return;
    const int j = (int)(t % naj);
    const int k = (int)(t / naj);
    // front, periodic (x = 1)
    AA(10,P.l1,j,k) = AA(10,1,j,k);
    AA(11,P.l1,j,k) = AA(11,1,j,k);
    AA(12,P.l1,j,k) = AA(12,1,j,k);
    AA(13,P.l1,j,k) = AA(13,1,j,k);
    AA(14,P.l1,j,k) = AA(14,1,j,k);
    // rear, periodic (x = l)
    AA(1,0,j,k) = AA(1,P.l,j,k);
    AA(2,0,j,k) = AA(2,P.l,j,k);
    AA(3,0,j,k) = AA(3,P.l,j,k);
    AA(4,0,j,k) = AA(4,P.l,j,k);
    AA(5,0,j,k) = AA(5,P.l,j,k);
}

// loop 2: y-faces periodic, identical to periodic loop 2
// (reuse bc_periodic_y at launch time)

// loop 3: z no-slip, over i=1..l, j=1..m
__global__ void bc_channel_z(BcArrays A, BcParams P)
{
    const long t = (long)blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= (long)P.l * P.m) return;
    const int i = (int)(t % P.l) + 1;
    const int j = (int)(t / P.l) + 1;
    // bottom, no slip (z = 1)
    AA(6,i,  j,  0)  = AA(15,i,j,1);
    AA(4,i-1,j,  0)  = AA(11,i,j,1);
    AA(7,i,  j-1,0)  = AA(16,i,j,1);
    AA(13,i+1,j, 0)  = AA(2,i,j,1);
    AA(18,i,  j+1,0) = AA(9,i,j,1);
    // top, no slip (z = n)
    AA(15,i,  j,  P.n1) = AA(6,i,j,P.n);
    AA(2,i-1,j,  P.n1)  = AA(13,i,j,P.n);
    AA(9,i,  j-1,P.n1)  = AA(18,i,j,P.n);
    AA(11,i+1,j, P.n1)  = AA(4,i,j,P.n);
    AA(16,i,  j+1,P.n1) = AA(7,i,j,P.n);
}

static inline void bcond_channel_gpu_launch(const BcArrays& arr, const BcParams& p)
{
    LAUNCH1D(bc_channel_x, ((long)p.m + 2) * ((long)p.n + 2), arr, p);
    LAUNCH1D(bc_periodic_y, ((long)p.l + 2) * ((long)p.n + 2), arr, p);
    LAUNCH1D(bc_channel_z, (long)p.l * p.m, arr, p);
}

// ====================================================================
//  DRIVEN (lid-driven cavity) boundary  (bcond.driven.F90)
//    interior loops 1..m/1..n etc. (OFFLOAD convention)
// ====================================================================
// loop 1: x-faces, over j=1..m, k=1..n
__global__ void bc_driven_x(BcArrays A, BcParams P)
{
    const long t = (long)blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= (long)P.m * P.n) return;
    const int j = (int)(t % P.m) + 1;
    const int k = (int)(t / P.m) + 1;
    // front (x = l)
    AA(12,P.l1,j-1,k  ) = AA(1,P.l,j,k);
    AA(13,P.l1,j,  k-1) = AA(2,P.l,j,k);
    AA(10,P.l1,j+1,k  ) = AA(3,P.l,j,k);
    AA(11,P.l1,j,  k+1) = AA(4,P.l,j,k);
    AA(14,P.l1,j,  k  ) = AA(5,P.l,j,k);
    // rear (x = 1)
    AA(3,0,j-1,k  ) = AA(10,1,j,k);
    AA(4,0,j,  k-1) = AA(11,1,j,k);
    AA(1,0,j+1,k  ) = AA(12,1,j,k);
    AA(2,0,j,  k+1) = AA(13,1,j,k);
    AA(5,0,j,  k  ) = AA(14,1,j,k);
}

// loop 2: y-faces, over i=1..l, k=1..n
__global__ void bc_driven_y(BcArrays A, BcParams P)
{
    const long t = (long)blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= (long)P.l * P.n) return;
    const int i = (int)(t % P.l) + 1;
    const int k = (int)(t / P.l) + 1;
    // left (y = 1)
    AA(8,i,  0,k  ) = AA(17,i,1,k);
    AA(12,i+1,0,k ) = AA(1,i,1,k);
    AA(3,i-1,0,k  ) = AA(10,i,1,k);
    AA(7,i,  0,k-1) = AA(16,i,1,k);
    AA(9,i,  0,k+1) = AA(18,i,1,k);
    // right (y = m)
    AA(10,i+1,P.m1,k  ) = AA(3,i,P.m,k);
    AA(16,i,  P.m1,k+1) = AA(7,i,P.m,k);
    AA(17,i,  P.m1,k  ) = AA(8,i,P.m,k);
    AA(18,i,  P.m1,k-1) = AA(9,i,P.m,k);
    AA(1,i-1, P.m1,k  ) = AA(12,i,P.m,k);
}

// loop 3: z-faces, over i=1..l, j=1..m  (lid forcing at top)
__global__ void bc_driven_z(BcArrays A, BcParams P)
{
    const long t = (long)blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= (long)P.l * P.m) return;
    const int i = (int)(t % P.l) + 1;
    const int j = (int)(t / P.l) + 1;
    const real_t force = P.force;
    // bottom (z = 1)
    AA(6,i,  j,  0) = AA(15,i,j,1);
    AA(4,i-1,j,  0) = AA(11,i,j,1);
    AA(7,i,  j-1,0) = AA(16,i,j,1);
    AA(13,i+1,j, 0) = AA(2,i,j,1);
    AA(18,i,  j+1,0) = AA(9,i,j,1);
    // top (z = n) lid-wall
    AA(15,i,  j,  P.n1) = AA(6,i,j,P.n);
    AA(2,i-1,j,   P.n1) = AA(13,i,j,P.n) + force;
    AA(9,i,  j-1,P.n1)  = AA(18,i,j,P.n);
    AA(11,i+1,j, P.n1)  = AA(4,i,j,P.n) - force;
    AA(16,i,  j+1,P.n1) = AA(7,i,j,P.n);
}

static inline void bcond_driven_gpu_launch(const BcArrays& arr, const BcParams& p)
{
    LAUNCH1D(bc_driven_x, (long)p.m * p.n, arr, p);
    LAUNCH1D(bc_driven_y, (long)p.l * p.n, arr, p);
    LAUNCH1D(bc_driven_z, (long)p.l * p.m, arr, p);
}

// ====================================================================
//  INFLOW boundary  (bcond.inflow.F90)
//    x: outflow (front) + prescribed inflow (rear); y/z periodic
// ====================================================================
// loop 1: x-faces, over j=0..m+1, k=0..n+1 (equilibrium reconstruction)
__global__ void bc_inflow_x(BcArrays A, BcParams P)
{
    const long naj = (long)P.m + 2;
    const long nak = (long)P.n + 2;
    const long t = (long)blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= naj * nak) return;
    const int j = (int)(t % naj);
    const int k = (int)(t / naj);

    const real_t rf   = P.rf;
    const real_t qf   = P.qf;
    const real_t cte1 = P.cte1;
    const real_t p1   = P.p1;
    const real_t p2   = P.p2;
    const real_t crho = (real_t)1;
    const real_t rhoinv = (real_t)1;
    const real_t three  = (real_t)3;

    // front, outflow (x = l): reconstruct from moments at (l,j,k)
    real_t xj = ( AA(1,P.l,j,k) + AA(2,P.l,j,k) + AA(3,P.l,j,k)
                + AA(4,P.l,j,k) + AA(5,P.l,j,k)
                - AA(10,P.l,j,k) - AA(11,P.l,j,k) - AA(12,P.l,j,k)
                - AA(13,P.l,j,k) - AA(14,P.l,j,k) ) * rhoinv;
    real_t yj = ( AA(3,P.l,j,k) + AA(7,P.l,j,k) + AA(8,P.l,j,k)
                + AA(9,P.l,j,k) + AA(12,P.l,j,k)
                - AA(1,P.l,j,k) - AA(10,P.l,j,k) - AA(16,P.l,j,k)
                - AA(17,P.l,j,k) - AA(18,P.l,j,k) ) * rhoinv;
    real_t zj = ( AA(4,P.l,j,k) + AA(6,P.l,j,k) + AA(7,P.l,j,k)
                + AA(13,P.l,j,k) + AA(18,P.l,j,k)
                - AA(2,P.l,j,k) - AA(9,P.l,j,k) - AA(11,P.l,j,k)
                - AA(15,P.l,j,k) - AA(16,P.l,j,k) ) * rhoinv;
    real_t cvsq = xj*xj + yj*yj + zj*zj;

    real_t cx10 = rf*(-xj-yj    ) + qf*(three*(xj+yj)*(xj+yj) - cvsq);
    real_t cx11 = rf*(-xj    -zj) + qf*(three*(xj+zj)*(xj+zj) - cvsq);
    real_t cx12 = rf*(-xj+yj    ) + qf*(three*(xj-yj)*(xj-yj) - cvsq);
    real_t cx13 = rf*(-xj    +zj) + qf*(three*(xj-zj)*(xj-zj) - cvsq);
    real_t cx14 = rf*(-xj       ) + qf*(three*(xj   )*(xj   ) - cvsq);

    AA(10,P.l1,j,k) = crho*p2*(cte1 + cx10);
    AA(11,P.l1,j,k) = crho*p2*(cte1 + cx11);
    AA(12,P.l1,j,k) = crho*p2*(cte1 + cx12);
    AA(13,P.l1,j,k) = crho*p2*(cte1 + cx13);
    AA(14,P.l1,j,k) = crho*p1*(cte1 + cx14);

    // rear, inflow (x = 0): prescribed velocity u_inflow along x
    xj = P.u_inflow;
    yj = (real_t)0;
    zj = (real_t)0;
    cvsq = xj*xj + yj*yj + zj*zj;

    real_t cx01 = rf*(xj-yj    ) + qf*(three*(xj-yj)*(xj-yj) - cvsq);
    real_t cx02 = rf*(xj    -zj) + qf*(three*(xj-zj)*(xj-zj) - cvsq);
    real_t cx03 = rf*(xj+yj    ) + qf*(three*(xj+yj)*(xj+yj) - cvsq);
    real_t cx04 = rf*(xj    +zj) + qf*(three*(xj+zj)*(xj+zj) - cvsq);
    real_t cx05 = rf*(xj       ) + qf*(three*(xj   )*(xj   ) - cvsq);

    AA(1,0,j,k) = crho*p2*(cte1 + cx01);
    AA(2,0,j,k) = crho*p2*(cte1 + cx02);
    AA(3,0,j,k) = crho*p2*(cte1 + cx03);
    AA(4,0,j,k) = crho*p2*(cte1 + cx04);
    AA(5,0,j,k) = crho*p1*(cte1 + cx05);
}

static inline void bcond_inflow_gpu_launch(const BcArrays& arr, const BcParams& p)
{
    LAUNCH1D(bc_inflow_x, ((long)p.m + 2) * ((long)p.n + 2), arr, p);
    LAUNCH1D(bc_periodic_y, ((long)p.l + 2) * ((long)p.n + 2), arr, p);
    LAUNCH1D(bc_periodic_z, ((long)p.l + 2) * ((long)p.m + 2), arr, p);
}

// ====================================================================
//  OBSTACLE bounce-back  (bcond.obs.F90)
//    over i=imin..imax, j=jmin..jmax, k=kmin..kmax where obs(i,j,k)==1
// ====================================================================
__global__ void bc_obs_kernel(BcArrays A, BcParams P, const int* obs)
{
    const long ni = (long)(P.imax - P.imin + 1);
    const long nj = (long)(P.jmax - P.jmin + 1);
    const long nk = (long)(P.kmax - P.kmin + 1);
    const long t = (long)blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= ni * nj * nk) return;
    const int i = P.imin + (int)(t % ni);
    const int j = P.jmin + (int)((t / ni) % nj);
    const int k = P.kmin + (int)(t / (ni * nj));

    const long oidx = (long)(i - 1) + P.obs_nx * (long)(j - 1) + P.obs_nxy * (long)(k - 1);
    if (obs[oidx] != 1) return;

    AA(1,i,j,k)  = AA(12,i+1,j-1,k  );
    AA(2,i,j,k)  = AA(13,i+1,j,  k-1);
    AA(3,i,j,k)  = AA(10,i+1,j+1,k  );
    AA(4,i,j,k)  = AA(11,i+1,j,  k+1);
    AA(5,i,j,k)  = AA(14,i+1,j,  k  );
    AA(6,i,j,k)  = AA(15,i,  j,  k+1);
    AA(7,i,j,k)  = AA(16,i,  j+1,k+1);
    AA(8,i,j,k)  = AA(17,i,  j+1,k  );
    AA(9,i,j,k)  = AA(18,i,  j+1,k-1);
    AA(10,i,j,k) = AA(3,i-1,j-1,k  );
    AA(11,i,j,k) = AA(4,i-1,j,  k-1);
    AA(12,i,j,k) = AA(1,i-1,j+1,k  );
    AA(13,i,j,k) = AA(2,i-1,j,  k+1);
    AA(14,i,j,k) = AA(5,i-1,j,  k  );
    AA(15,i,j,k) = AA(6,i,  j,  k-1);
    AA(16,i,j,k) = AA(7,i,  j-1,k-1);
    AA(17,i,j,k) = AA(8,i,  j-1,k  );
    AA(18,i,j,k) = AA(9,i,  j-1,k+1);
}

static inline void bcond_obs_gpu_launch(const BcArrays& arr, const BcParams& p, const int* obs)
{
    const long count = (long)(p.imax - p.imin + 1)
                     * (long)(p.jmax - p.jmin + 1)
                     * (long)(p.kmax - p.kmin + 1);
#if defined(USE_HIP)
    hipLaunchKernelGGL(bc_obs_kernel, dim3(grid_for(count, 256)), dim3(256), 0, 0, arr, p, obs);
#else
    bc_obs_kernel<<<grid_for(count, 256), 256>>>(arr, p, obs);
#endif
    GPU_CHECK_LAUNCH("bc_obs_kernel");
}

// ====================================================================
//  extern "C" entry points called from the Fortran GPU_NATIVE branches.
//  `a` is the host array pa(1:19) of DEVICE pointers to a01..a19 (each
//  resolved with c_loc + the residency model, exactly as col_MC does).
//  Only ONE back-end translation unit (hip XOR cuda) includes this
//  header per build, so these definitions never clash.
// ====================================================================
static inline void fill_arrays(BcArrays& arr, real_t** a)
{
    for (int d = 0; d < 19; ++d) arr.a[d] = a[d];
}

static inline BcParams make_params(int l, int m, int n, int nx, int ny)
{
    BcParams p;
    p.l = l; p.m = m; p.n = n;
    p.l1 = l + 1; p.m1 = m + 1; p.n1 = n + 1;
    p.nx = (long)nx; p.nxy = (long)nx * (long)ny;
    p.imin = p.imax = p.jmin = p.jmax = p.kmin = p.kmax = 0;
    p.obs_nx = p.obs_nxy = 0;
    p.force = (real_t)0;
    p.u_inflow = (real_t)0;
    p.rf = p.qf = p.tre = (real_t)0;
    p.p0 = p.p1 = p.p2 = (real_t)0;
    p.cte1 = (real_t)0;
    return p;
}

extern "C" void bcond_periodic_gpu(real_t** a, int l, int m, int n, int nx, int ny)
{
    BcArrays arr; fill_arrays(arr, a);
    BcParams p = make_params(l, m, n, nx, ny);
    bcond_periodic_gpu_launch(arr, p);
}

extern "C" void bcond_channel_gpu(real_t** a, int l, int m, int n, int nx, int ny)
{
    BcArrays arr; fill_arrays(arr, a);
    BcParams p = make_params(l, m, n, nx, ny);
    bcond_channel_gpu_launch(arr, p);
}

extern "C" void bcond_driven_gpu(real_t** a, int l, int m, int n, int nx, int ny, real_t force)
{
    BcArrays arr; fill_arrays(arr, a);
    BcParams p = make_params(l, m, n, nx, ny);
    p.force = force;
    bcond_driven_gpu_launch(arr, p);
}

extern "C" void bcond_inflow_gpu(real_t** a, int l, int m, int n, int nx, int ny,
                                 real_t u_inflow, real_t rf, real_t qf, real_t tre,
                                 real_t p1, real_t p2, real_t cte1)
{
    BcArrays arr; fill_arrays(arr, a);
    BcParams p = make_params(l, m, n, nx, ny);
    p.u_inflow = u_inflow;
    p.rf = rf; p.qf = qf; p.tre = tre;
    p.p1 = p1; p.p2 = p2; p.cte1 = cte1;
    bcond_inflow_gpu_launch(arr, p);
}

extern "C" void bcond_obs_gpu(real_t** a, const int* obs, int l, int m, int n, int nx, int ny,
                              int imin, int imax, int jmin, int jmax, int kmin, int kmax,
                              int obs_nx, int obs_ny)
{
    BcArrays arr; fill_arrays(arr, a);
    BcParams p = make_params(l, m, n, nx, ny);
    p.imin = imin; p.imax = imax;
    p.jmin = jmin; p.jmax = jmax;
    p.kmin = kmin; p.kmax = kmax;
    p.obs_nx = (long)obs_nx; p.obs_nxy = (long)obs_nx * (long)obs_ny;
    bcond_obs_gpu_launch(arr, p, obs);
}

// Full-device barrier used by the Fortran side at host-read points only
// (diagno / dissipation / draglift / save), replacing the per-step sync.
extern "C" void gpu_device_sync()
{
    // A device-sync error means a kernel faulted during execution: every
    // result after this point (and the reported Mlups) is garbage.  Fail
    // loudly instead of continuing silently, so a broken run can never
    // masquerade as a fast, valid one.
#if defined(USE_HIP)
    hipError_t e = hipDeviceSynchronize();
    if (e != hipSuccess) {
        fprintf(stderr, "ERROR: HIP device sync failed: %s at %s:%d\n",
                hipGetErrorString(e), __FILE__, __LINE__);
        fflush(stderr);
        abort();
    }
#else
    cudaError_t e = cudaDeviceSynchronize();
    if (e != cudaSuccess) {
        fprintf(stderr, "ERROR: CUDA device sync failed: %s at %s:%d\n",
                cudaGetErrorString(e), __FILE__, __LINE__);
        fflush(stderr);
        abort();
    }
#endif
}

#endif // BCOND_GPU_HPP
