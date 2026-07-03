!=====================================================================
!     ****** LBE/bcond_channel
!
!     COPYRIGHT
!       (c) 2000-2011 by CASPUR/G.Amati
!       (c) 2013-20?? by CINECA/G.Amati
!     NAME
!       bcond_channel: simple (periodic) channel flow with
!                      * periodic (rear/front)
!                      * periodic (left/right)
!                      * no slip (top/bottom)
!     DESCRIPTION
!       3D Channel
!     INPUTS
!       none
!     OUTPUT
!       none
!     TODO
!
!     NOTES
!       Order of upgrading bc
!       1) front (x = l)
!       2) rear  (x = 0)
!       3) left  (y = 0)
!       4) right (y = m)
!       5) right (z = 0)
!       6) right (z = n)
!
!     *****
!=====================================================================
!
        subroutine bcond_channel
!
        use timing
        use storage
#ifdef GPU_NATIVE
        use iso_c_binding
        use bcond_gpu_mod
#elif OFFLOAD_KERNEL_SYNTAX
        use omp_lib
#endif
!
        implicit none
!
        integer      :: i,j,k
        real(mykind) :: cte1
#ifdef GPU_NATIVE
        type(c_ptr) :: pa(19)
#elif OFFLOAD_KERNEL_SYNTAX
        integer :: na, nb
#  ifdef KS_BLOCK_3D
        integer :: ntx, nty, nbx, nby
#  else
        integer :: tid, ncell, nthrd, nblck
#  endif
#endif
!
#ifdef CHANNEL
!
# ifdef NOSHIFT
       cte1 = uno
# else
       cte1 = zero
# endif

!
#ifdef PROFILING
! start timing...
        call SYSTEM_CLOCK(countA0, count_rate, count_max)
        call time(tcountA0)
#endif
!
#ifdef GPU_NATIVE
        call resolve_a_devptrs(pa)
        call bcond_channel_gpu(pa,                                      &
             int(l,c_int),int(m,c_int),int(n,c_int),                    &
             int(size(a01,1),c_int),int(size(a01,2),c_int))
#else
!
! ----------------------------------------------
! loop foused for performance reason (for GPU)
! -------------------------------------------------------------
!
! ----------------------------------------------
! rear (x = 1)  ! front (x = n)
! ----------------------------------------------
!
#ifdef OFFLOAD_KERNEL_SYNTAX
        na = m + 2
        nb = n + 2
#  ifdef KS_BLOCK_3D
        ntx = 16
        nty = 16
        nbx = (na + ntx - 1)/ntx
        nby = (nb + nty - 1)/nty
!$omp target teams parallel thread_limit(dims(3):ntx,nty,1) num_teams(dims(3):nbx,nby,1)
        j = omp_get_thread_num_dim(0) + omp_get_team_num_dim(0)*omp_get_num_threads_dim(0)
        k = omp_get_thread_num_dim(1) + omp_get_team_num_dim(1)*omp_get_num_threads_dim(1)
        if (j <= m+1 .and. k <= n+1) then
#  else
        nthrd = 256
        ncell = na*nb
        nblck = (ncell + nthrd - 1)/nthrd
!$omp target teams parallel num_threads(dims(3):nthrd) num_teams(dims(3):nblck) thread_limit(nthrd)
        tid = omp_get_thread_num_dim(0) + omp_get_team_num_dim(0)*omp_get_num_threads_dim(0)
        if (tid < ncell) then
          j = mod(tid, na)
          k = tid/na
#  endif
#elif defined(OFFLOAD)
!$OMP target teams distribute parallel do simd collapse(2)
        do k=0,n+1
        do j=0,m+1
#elif OPENACC
!$acc parallel
!$acc loop independent collapse(2)
        do k=0,n+1
        do j=0,m+1
#else
        do concurrent (k=0:n+1,j=0:m+1)
#endif
!
! front, periodic bc (x = 1)
           a10(l1,j,k) = a10(1,j,k)
           a11(l1,j,k) = a11(1,j,k)
           a12(l1,j,k) = a12(1,j,k)
           a13(l1,j,k) = a13(1,j,k)
           a14(l1,j,k) = a14(1,j,k)
!
! rear, periodic bc (x = l)
           a01(0,j,k)  = a01(l,j,k)
           a02(0,j,k)  = a02(l,j,k)
           a03(0,j,k)  = a03(l,j,k)
           a04(0,j,k)  = a04(l,j,k)
           a05(0,j,k)  = a05(l,j,k)
#ifdef OFFLOAD_KERNEL_SYNTAX
        end if
!$omp end target teams parallel
#elif defined(OFFLOAD)
        end do
        end do
!$OMP end target teams distribute parallel do simd
#elif OPENACC
        end do
        end do
!$acc end parallel
#else
        end do
#endif
!
! ----------------------------------------------
! left (y = 1)  ! right (y = m)
! ----------------------------------------------
!
#ifdef OFFLOAD_KERNEL_SYNTAX
        na = l + 2
        nb = n + 2
#  ifdef KS_BLOCK_3D
        ntx = 16
        nty = 16
        nbx = (na + ntx - 1)/ntx
        nby = (nb + nty - 1)/nty
!$omp target teams parallel thread_limit(dims(3):ntx,nty,1) num_teams(dims(3):nbx,nby,1)
        i = omp_get_thread_num_dim(0) + omp_get_team_num_dim(0)*omp_get_num_threads_dim(0)
        k = omp_get_thread_num_dim(1) + omp_get_team_num_dim(1)*omp_get_num_threads_dim(1)
        if (i <= l+1 .and. k <= n+1) then
#  else
        nthrd = 256
        ncell = na*nb
        nblck = (ncell + nthrd - 1)/nthrd
!$omp target teams parallel num_threads(dims(3):nthrd) num_teams(dims(3):nblck) thread_limit(nthrd)
        tid = omp_get_thread_num_dim(0) + omp_get_team_num_dim(0)*omp_get_num_threads_dim(0)
        if (tid < ncell) then
          i = mod(tid, na)
          k = tid/na
#  endif
#elif defined(OFFLOAD)
!$OMP target teams distribute parallel do simd collapse(2)
        do k=0,n+1
        do i=0,l+1
#elif OPENACC
!$acc parallel
!$acc loop independent collapse(2)
        do k=0,n+1
        do i=0,l+1
#else
        do concurrent (k=0:n+1,i=0:l+1)
#endif
! left, periodic bc (y = 0)
           a03(i,0,k) = a03(i,m,k)
           a07(i,0,k) = a07(i,m,k)
           a08(i,0,k) = a08(i,m,k)
           a09(i,0,k) = a09(i,m,k)
           a12(i,0,k) = a12(i,m,k)
!
! right, periodic bc (y = m)
           a01(i,m1,k) = a01(i,1,k)
           a10(i,m1,k) = a10(i,1,k)
           a16(i,m1,k) = a16(i,1,k)
           a17(i,m1,k) = a17(i,1,k)
           a18(i,m1,k) = a18(i,1,k)
#ifdef OFFLOAD_KERNEL_SYNTAX
        end if
!$omp end target teams parallel
#elif defined(OFFLOAD)
        enddo
        enddo
!$OMP end target teams distribute parallel do simd
#elif OPENACC
        enddo
        enddo
!$acc end parallel
#else
        enddo
#endif
!
! ----------------------------------------------
! down (z = 1)  ! top (z = n)
! ----------------------------------------------
!
#ifdef OFFLOAD_KERNEL_SYNTAX
        na = l
        nb = m
#  ifdef KS_BLOCK_3D
        ntx = 16
        nty = 16
        nbx = (na + ntx - 1)/ntx
        nby = (nb + nty - 1)/nty
!$omp target teams parallel thread_limit(dims(3):ntx,nty,1) num_teams(dims(3):nbx,nby,1)
        i = 1 + omp_get_thread_num_dim(0) + omp_get_team_num_dim(0)*omp_get_num_threads_dim(0)
        j = 1 + omp_get_thread_num_dim(1) + omp_get_team_num_dim(1)*omp_get_num_threads_dim(1)
        if (i <= l .and. j <= m) then
#  else
        nthrd = 256
        ncell = na*nb
        nblck = (ncell + nthrd - 1)/nthrd
!$omp target teams parallel num_threads(dims(3):nthrd) num_teams(dims(3):nblck) thread_limit(nthrd)
        tid = omp_get_thread_num_dim(0) + omp_get_team_num_dim(0)*omp_get_num_threads_dim(0)
        if (tid < ncell) then
          i = 1 + mod(tid, na)
          j = 1 + tid/na
#  endif
#elif defined(OFFLOAD)
!$OMP target teams distribute parallel do simd collapse(2)
        do j=1,m
        do i=1,l
#elif OPENACC
!$acc parallel
!$acc loop independent collapse(2)
        do j=1,m
        do i=1,l
#else
        do concurrent (j=0:m+1,i=0:l+1)
#endif
! bottom, no slip bc (z = 1)
           a06(i  ,j  ,0)  = a15(i,j,1)
           a04(i-1,j  ,0)  = a11(i,j,1)
           a07(i  ,j-1,0)  = a16(i,j,1)
           a13(i+1,j  ,0)  = a02(i,j,1)
           a18(i  ,j+1,0)  = a09(i,j,1)
!
! top, no slip bc (y = m)
           a15(i  ,j  ,n1) = a06(i,j,n)
           a02(i-1,j  ,n1) = a13(i,j,n)
           a09(i  ,j-1,n1) = a18(i,j,n)
           a11(i+1,j  ,n1) = a04(i,j,n)
           a16(i  ,j+1,n1) = a07(i,j,n)
#ifdef OFFLOAD_KERNEL_SYNTAX
        end if
!$omp end target teams parallel
#elif defined(OFFLOAD)
        enddo
        enddo
!$OMP end target teams distribute parallel do simd
#elif OPENACC
        enddo
        enddo
!$acc end parallel
#else
        enddo
#endif
#endif

!
! ----------------------------------------------
#ifdef PROFILING
! stop timing
#ifdef GPU_NATIVE
! flush the async native kernels so the CPU timer captures real GPU time
        call gpu_device_sync()
#endif
        call time(tcountA1)
        call SYSTEM_CLOCK(countA1, count_rate, count_max)
        time_bc = time_bc + real(countA1-countA0)/(count_rate)
        time_bc1 = time_bc1 + (tcountA1-tcountA0)
#endif
!
#endif
!
#ifdef DEBUG_2
        if(myrank == 0) then
           write(6,*) "DEBUG2: Exiting from sub. bcond_channel"
        endif
#endif
        end subroutine bcond_channel
