!=====================================================================
!     ****** LBE/bcond_driven
!
!     COPYRIGHT
!       (c) 2000-2011 by CASPUR/G.Amati
!       (c) 2013-20?? by CINECA/G.Amati
!     NAME
!       bcond
!     DESCRIPTION
!
!     INPUTS
!       none
!     OUTPUT
!       none
!     TODO
!
!     NOTES
!       Order of upgrading bc
!       1) front  (x = l)
!       2) rear   (x = 0)
!       3) left   (y = 0)
!       4) right  (y = m)
!       5) bottom (z = 0)
!       6) top    (z = n)
!
!     *****
!=====================================================================
!
        subroutine bcond_driven
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
#ifdef UNIFIED_MEMORY
!$omp requires unified_shared_memory
#endif
!
        integer      :: i,j,k
        real(mykind) :: force
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
#ifdef PROFILING
! start timing...
        call SYSTEM_CLOCK(countA0, count_rate, count_max)
        call time(tcountA0)
#endif
!
        force =  u00/(6.0)
!
#ifdef GPU_NATIVE
        call resolve_a_devptrs(pa)
        call bcond_driven_gpu(pa,                                       &
             int(l,c_int),int(m,c_int),int(n,c_int),                    &
             int(size(a01,1),c_int),int(size(a01,2),c_int),             &
             real(force,c_real))
#else
!
! bc. along X direction
!
#ifdef OFFLOAD_KERNEL_SYNTAX
        na = m
        nb = n
#  ifdef KS_BLOCK_3D
        ntx = 16
        nty = 16
        nbx = (na + ntx - 1)/ntx
        nby = (nb + nty - 1)/nty
!$omp target teams parallel thread_limit(dims(3):ntx,nty,1) num_teams(dims(3):nbx,nby,1)
        j = 1 + omp_get_thread_num_dim(0) + omp_get_team_num_dim(0)*omp_get_num_threads_dim(0)
        k = 1 + omp_get_thread_num_dim(1) + omp_get_team_num_dim(1)*omp_get_num_threads_dim(1)
        if (j <= m .and. k <= n) then
#  else
        nthrd = 256
        ncell = na*nb
        nblck = (ncell + nthrd - 1)/nthrd
!$omp target teams parallel num_threads(dims(3):nthrd) num_teams(dims(3):nblck) thread_limit(nthrd)
        tid = omp_get_thread_num_dim(0) + omp_get_team_num_dim(0)*omp_get_num_threads_dim(0)
        if (tid < ncell) then
          j = 1 + mod(tid, na)
          k = 1 + tid/na
#  endif
#elif defined(OFFLOAD)
!$OMP target teams distribute parallel do simd collapse(2)
        do k=1,n
        do j=1,m
#elif OPENACC
!$acc parallel
!$acc loop independent collapse(2)
        do k=1,n
        do j=1,m
#else
        do concurrent (j=1:m,k=1:n)
#endif
! front (x = l)
           a12(l1,j-1,k  ) = a01(l,j,k)
           a13(l1,j  ,k-1) = a02(l,j,k)
           a10(l1,j+1,k  ) = a03(l,j,k)
           a11(l1,j  ,k+1) = a04(l,j,k)
           a14(l1,j  ,k  ) = a05(l,j,k)
!
! rear (x = 1)
           a03(0,j-1,k  ) = a10(1,j,k)
           a04(0,j  ,k-1) = a11(1,j,k)
           a01(0,j+1,k  ) = a12(1,j,k)
           a02(0,j  ,k+1) = a13(1,j,k)
           a05(0,j  ,k  ) = a14(1,j,k)
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
! bc. along y direction
!
#ifdef OFFLOAD_KERNEL_SYNTAX
        na = l
        nb = n
#  ifdef KS_BLOCK_3D
        ntx = 16
        nty = 16
        nbx = (na + ntx - 1)/ntx
        nby = (nb + nty - 1)/nty
!$omp target teams parallel thread_limit(dims(3):ntx,nty,1) num_teams(dims(3):nbx,nby,1)
        i = 1 + omp_get_thread_num_dim(0) + omp_get_team_num_dim(0)*omp_get_num_threads_dim(0)
        k = 1 + omp_get_thread_num_dim(1) + omp_get_team_num_dim(1)*omp_get_num_threads_dim(1)
        if (i <= l .and. k <= n) then
#  else
        nthrd = 256
        ncell = na*nb
        nblck = (ncell + nthrd - 1)/nthrd
!$omp target teams parallel num_threads(dims(3):nthrd) num_teams(dims(3):nblck) thread_limit(nthrd)
        tid = omp_get_thread_num_dim(0) + omp_get_team_num_dim(0)*omp_get_num_threads_dim(0)
        if (tid < ncell) then
          i = 1 + mod(tid, na)
          k = 1 + tid/na
#  endif
#elif defined(OFFLOAD)
!$OMP target teams distribute parallel do simd collapse(2)
        do k=1,n
        do i=1,l
#elif OPENACC
!$acc parallel
!$acc loop independent collapse(2)
        do k=1,n
        do i=1,l
#else
        do concurrent (k=1:n,i=1:l)
#endif
! left (y = 1)
            a08(i  ,0,k  )  = a17(i,1,k)
            a12(i+1,0,k  )  = a01(i,1,k)
            a03(i-1,0,k  )  = a10(i,1,k)
            a07(i  ,0,k-1)  = a16(i,1,k)
            a09(i  ,0,k+1)  = a18(i,1,k)

! right (y = m)
            a10(i+1,m1,k  ) = a03(i,m,k)
            a16(i  ,m1,k+1) = a07(i,m,k)
            a17(i  ,m1,k  ) = a08(i,m,k)
            a18(i  ,m1,k-1) = a09(i,m,k)
            a01(i-1,m1,k  ) = a12(i,m,k)
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
! bc. along z direction
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
        do concurrent (j=1:m,i=1:l)
#endif
! bottom (z = 1)
            a06(i  ,j  ,0)  = a15(i,j,1)
            a04(i-1,j  ,0)  = a11(i,j,1)
            a07(i  ,j-1,0)  = a16(i,j,1)
            a13(i+1,j  ,0)  = a02(i,j,1)
            a18(i  ,j+1,0)  = a09(i,j,1)
!
! top (z = n) lid-wall
            a15(i  ,j  ,n1) = a06(i,j,n)
            a02(i-1,j  ,n1) = a13(i,j,n) + force
            a09(i  ,j-1,n1) = a18(i,j,n)
            a11(i+1,j  ,n1) = a04(i,j,n) - force
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
#ifdef DEBUG_2
        if(myrank == 0) then
           write(6,*) "DEBUG2: Exiting from sub. bcond_driven", force
        endif
#endif
!
        end subroutine bcond_driven
