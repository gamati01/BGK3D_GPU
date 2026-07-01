!=====================================================================
!     ****** LBE/bcond_gpu_mod
!
!     DESCRIPTION
!       Fortran glue for the native GPU (HIP/CUDA) boundary kernels.
!       Provides the ISO-C interfaces to the extern "C" entry points in
!       bcond_gpu.hpp (bcond_hip.cpp / bcond_cuda.cu) and a helper that
!       resolves the CURRENT device base address of a01..a19 (following
!       the a<->b pointer swap), mirroring the logic in col_MC.F90.
!
!       Active only for the native GPU build (-DGPU_NATIVE); otherwise
!       this compiles to an empty module.
!=====================================================================
        module bcond_gpu_mod
#ifdef GPU_NATIVE
        use iso_c_binding
        use storage
#if !defined(OPENACC)
        use omp_lib
#endif
        implicit none
!
! Interoperable storage kind (must match mystorage in precision.F90).
#ifdef DOUBLE_P
        integer, parameter :: c_real = c_double
#else
        integer, parameter :: c_real = c_float
#endif
!
        interface
          subroutine bcond_periodic_gpu(a,l,m,n,nx,ny)                  &
                     bind(C, name="bcond_periodic_gpu")
            import :: c_ptr, c_int
            type(c_ptr)           :: a(19)
            integer(c_int), value :: l,m,n,nx,ny
          end subroutine bcond_periodic_gpu
!
          subroutine bcond_channel_gpu(a,l,m,n,nx,ny)                   &
                     bind(C, name="bcond_channel_gpu")
            import :: c_ptr, c_int
            type(c_ptr)           :: a(19)
            integer(c_int), value :: l,m,n,nx,ny
          end subroutine bcond_channel_gpu
!
          subroutine bcond_driven_gpu(a,l,m,n,nx,ny,force)              &
                     bind(C, name="bcond_driven_gpu")
            import :: c_ptr, c_int, c_real
            type(c_ptr)           :: a(19)
            integer(c_int), value :: l,m,n,nx,ny
            real(c_real),   value :: force
          end subroutine bcond_driven_gpu
!
          subroutine bcond_inflow_gpu(a,l,m,n,nx,ny,                    &
                     u_inflow,rf,qf,tre,p1,p2,cte1)                     &
                     bind(C, name="bcond_inflow_gpu")
            import :: c_ptr, c_int, c_real
            type(c_ptr)           :: a(19)
            integer(c_int), value :: l,m,n,nx,ny
            real(c_real),   value :: u_inflow,rf,qf,tre,p1,p2,cte1
          end subroutine bcond_inflow_gpu
!
#ifdef OBSTACLE
          subroutine bcond_obs_gpu(a,obs,l,m,n,nx,ny,                   &
                     imin,imax,jmin,jmax,kmin,kmax,obs_nx,obs_ny)       &
                     bind(C, name="bcond_obs_gpu")
            import :: c_ptr, c_int
            type(c_ptr)           :: a(19)
            type(c_ptr),    value :: obs
            integer(c_int), value :: l,m,n,nx,ny
            integer(c_int), value :: imin,imax,jmin,jmax,kmin,kmax
            integer(c_int), value :: obs_nx,obs_ny
          end subroutine bcond_obs_gpu
#endif
!
          subroutine gpu_device_sync() bind(C, name="gpu_device_sync")
          end subroutine gpu_device_sync
        end interface
!
        contains
!
! Resolve the CURRENT device base address of a01..a19 into pa(1:19),
! honouring the residency model (managed / NOMANAGED / OpenMP), exactly
! like the col_MC.F90 GPU_NATIVE path so the a<->b swap is tracked.
        subroutine resolve_a_devptrs(pa)
          type(c_ptr), intent(out) :: pa(19)
#if !defined(OPENACC)
          integer :: idev, ii
#endif
          pa( 1)=c_loc(a01); pa( 2)=c_loc(a02); pa( 3)=c_loc(a03)
          pa( 4)=c_loc(a04); pa( 5)=c_loc(a05); pa( 6)=c_loc(a06)
          pa( 7)=c_loc(a07); pa( 8)=c_loc(a08); pa( 9)=c_loc(a09)
          pa(10)=c_loc(a10); pa(11)=c_loc(a11); pa(12)=c_loc(a12)
          pa(13)=c_loc(a13); pa(14)=c_loc(a14); pa(15)=c_loc(a15)
          pa(16)=c_loc(a16); pa(17)=c_loc(a17); pa(18)=c_loc(a18)
          pa(19)=c_loc(a19)
#ifdef OPENACC
#ifdef NOMANAGED
!$acc host_data use_device(                                            &
!$acc&   a01,a02,a03,a04,a05,a06,a07,a08,a09,a10,                       &
!$acc&   a11,a12,a13,a14,a15,a16,a17,a18,a19)
          pa( 1)=c_loc(a01); pa( 2)=c_loc(a02); pa( 3)=c_loc(a03)
          pa( 4)=c_loc(a04); pa( 5)=c_loc(a05); pa( 6)=c_loc(a06)
          pa( 7)=c_loc(a07); pa( 8)=c_loc(a08); pa( 9)=c_loc(a09)
          pa(10)=c_loc(a10); pa(11)=c_loc(a11); pa(12)=c_loc(a12)
          pa(13)=c_loc(a13); pa(14)=c_loc(a14); pa(15)=c_loc(a15)
          pa(16)=c_loc(a16); pa(17)=c_loc(a17); pa(18)=c_loc(a18)
          pa(19)=c_loc(a19)
!$acc end host_data
#endif
          ! managed memory: host address IS the device address (no-op)
#else
          idev = omp_get_default_device()
          do ii=1,19
             pa(ii) = omp_get_mapped_ptr(pa(ii), idev)
          enddo
#endif
        end subroutine resolve_a_devptrs
!
#ifdef OBSTACLE
! Resolve the device base address of the obstacle mask.
        subroutine resolve_obs_devptr(pobs)
          type(c_ptr), intent(out) :: pobs
#if !defined(OPENACC)
          integer :: idev
#endif
          pobs = c_loc(obs)
#ifdef OPENACC
#ifdef NOMANAGED
!$acc host_data use_device(obs)
          pobs = c_loc(obs)
!$acc end host_data
#endif
#else
          idev = omp_get_default_device()
          pobs = omp_get_mapped_ptr(pobs, idev)
#endif
        end subroutine resolve_obs_devptr
#endif
!
#endif
        end module bcond_gpu_mod
