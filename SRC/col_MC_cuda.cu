// ====================================================================
//  col_MC_cuda.cu
//
//  CUDA back-end for the native GPU port of the fused collision kernel.
//  Compiled with nvcc and -DUSE_CUDA.  The kernel + launch logic live in
//  the shared header col_MC_gpu.hpp (identical to the HIP back-end).
//
//  The Fortran caller (col_MC.F90, GPU_NATIVE path) passes DEVICE
//  pointers obtained inside an OpenACC `host_data use_device` region
//  (residency handled by OpenACC managed memory on NVIDIA), plus the
//  scalar parameters by value.
// ====================================================================
#define USE_CUDA
#include "col_MC_gpu.hpp"

extern "C" void col_mc_gpu(
    const real_t* a01, const real_t* a02, const real_t* a03, const real_t* a04,
    const real_t* a05, const real_t* a06, const real_t* a07, const real_t* a08,
    const real_t* a09, const real_t* a10, const real_t* a11, const real_t* a12,
    const real_t* a13, const real_t* a14, const real_t* a15, const real_t* a16,
    const real_t* a17, const real_t* a18, const real_t* a19,
    real_t* b01, real_t* b02, real_t* b03, real_t* b04, real_t* b05, real_t* b06,
    real_t* b07, real_t* b08, real_t* b09, real_t* b10, real_t* b11, real_t* b12,
    real_t* b13, real_t* b14, real_t* b15, real_t* b16, real_t* b17, real_t* b18,
    int l, int m, int n, int nx, int ny,
    real_t omega, real_t cte0, real_t cte1,
    real_t p0, real_t p1, real_t p2,
    real_t rf, real_t qf, real_t tre,
    real_t forcex, real_t forcey, real_t forcez)
{
    ColArrays arr;
    arr.a[0]  = a01; arr.a[1]  = a02; arr.a[2]  = a03; arr.a[3]  = a04;
    arr.a[4]  = a05; arr.a[5]  = a06; arr.a[6]  = a07; arr.a[7]  = a08;
    arr.a[8]  = a09; arr.a[9]  = a10; arr.a[10] = a11; arr.a[11] = a12;
    arr.a[12] = a13; arr.a[13] = a14; arr.a[14] = a15; arr.a[15] = a16;
    arr.a[16] = a17; arr.a[17] = a18; arr.a[18] = a19;

    arr.b[0]  = b01; arr.b[1]  = b02; arr.b[2]  = b03; arr.b[3]  = b04;
    arr.b[4]  = b05; arr.b[5]  = b06; arr.b[6]  = b07; arr.b[7]  = b08;
    arr.b[8]  = b09; arr.b[9]  = b10; arr.b[10] = b11; arr.b[11] = b12;
    arr.b[12] = b13; arr.b[13] = b14; arr.b[14] = b15; arr.b[15] = b16;
    arr.b[16] = b17; arr.b[17] = b18;
    arr.b[18] = const_cast<real_t*>(a19);   // direction 19 collides in place into a19

    ColParams p;
    p.l = l; p.m = m; p.n = n; p.nx = nx; p.ny = ny;
    p.omega = omega; p.cte0 = cte0; p.cte1 = cte1;
    p.p0 = p0; p.p1 = p1; p.p2 = p2;
    p.rf = rf; p.qf = qf; p.tre = tre;
    p.forcex = forcex; p.forcey = forcey; p.forcez = forcez;

    col_mc_gpu_launch(arr, p);
}
