!=======================================================================
!     ****** LBE/bcond_obs
!
!     COPYRIGHT
!       (c) 2000-2011 by CASPUR/G.Amati
!       (c) 2013-20?? by CINECA/G.Amati
!     NAME
!       bcond_obs
!     DESCRIPTION
!       obstcle
!     INPUTS
!       
!     OUTPUT
!       none
!     TODO
!
!     NOTES
!     *****
!=======================================================================
!
        subroutine bcond_obs
!
        use timing
        use storage
#if defined(GPU_NATIVE) && defined(OBSTACLE)
        use iso_c_binding
        use bcond_gpu_mod
#endif
!
        implicit none
!
        integer:: i,j,k
#if defined(GPU_NATIVE) && defined(OBSTACLE)
        type(c_ptr) :: pa(19), pobs
#endif
!
#ifdef PROFILING
! start timing
        call SYSTEM_CLOCK(countO0, count_rate, count_max)
        call time(tcountO0)
#endif
!
#if defined(GPU_NATIVE) && defined(OBSTACLE)
        call resolve_a_devptrs(pa)
        call resolve_obs_devptr(pobs)
        call bcond_obs_gpu(pa,pobs,                                     &
             int(l,c_int),int(m,c_int),int(n,c_int),                    &
             int(size(a01,1),c_int),int(size(a01,2),c_int),             &
             int(imin,c_int),int(imax,c_int),int(jmin,c_int),           &
             int(jmax,c_int),int(kmin,c_int),int(kmax,c_int),           &
             int(l,c_int),int(m,c_int))
#else
#ifdef OFFLOAD
!$OMP target teams distribute parallel do simd collapse(3)
        do k=kmin,kmax
        do j=jmin,jmax
        do i=imin,imax
#elif OPENACC
!$acc parallel
!$acc loop independent collapse(3)
        do k=kmin,kmax
        do j=jmin,jmax
        do i=imin,imax
#else
        do concurrent (k=kmin:kmax, j=jmin:jmax, i=imin:imax)
#endif
            if(obs(i,j,k)==1) then
                   a01(i,j,k) = a12(i+1,j-1,k  )
                   a02(i,j,k) = a13(i+1,j  ,k-1)
                   a03(i,j,k) = a10(i+1,j+1,k  )
                   a04(i,j,k) = a11(i+1,j  ,k+1)
                   a05(i,j,k) = a14(i+1,j  ,k  )
                   a06(i,j,k) = a15(i  ,j  ,k+1)
                   a07(i,j,k) = a16(i  ,j+1,k+1)
                   a08(i,j,k) = a17(i  ,j+1,k  )
                   a09(i,j,k) = a18(i  ,j+1,k-1)
                   a10(i,j,k) = a03(i-1,j-1,k  )
                   a11(i,j,k) = a04(i-1,j  ,k-1)
                   a12(i,j,k) = a01(i-1,j+1,k  )
                   a13(i,j,k) = a02(i-1,j  ,k+1)
                   a14(i,j,k) = a05(i-1,j  ,k  )
                   a15(i,j,k) = a06(i  ,j  ,k-1)
                   a16(i,j,k) = a07(i  ,j-1,k-1)
                   a17(i,j,k) = a08(i  ,j-1,k  )
                   a18(i,j,k) = a09(i  ,j-1,k+1)
            endif
        end do
#ifdef OFFLOAD
        end do
        end do
!$OMP end target teams distribute parallel do simd
#elif OPENACC
        end do
        end do
!$acc end parallel
#endif
#endif
!
#ifdef PROFILING
! stop timing
#if defined(GPU_NATIVE) && defined(OBSTACLE)
! flush the async native kernels so the CPU timer captures real GPU time
        call gpu_device_sync()
#endif
        call time(tcountO1)
        call SYSTEM_CLOCK(countO1, count_rate, count_max)
        time_obs = time_obs + real(countO1-countO0)/(count_rate)
        time_obs1 = time_obs1 + (tcountO1-tcountO0)
#endif
!
#ifdef DEBUG_2
        if(myrank == 0) then
           write(6,*) "DEBUG2: Exiting from sub. bcond_obs", &
                               imin,imax,jmin,imax
        endif
#endif
!
        end subroutine bcond_obs
