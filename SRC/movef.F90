!====================================================
!     ****** LBE/movef
!
!     COPYRIGHT
!       (c) 2000-2011 by CASPUR/G.Amati
!       (c) 2013-20?? by CINECA/G.Amati
!     NAME
!       movef
!     DESCRIPTION
!       "In place" streaming of populations
!     INPUTS
!       none
!     OUTPUT
!       none
!     TODO
!
!     NOTES
!       integer variables used: i,j
!
!     *****
!====================================================
!
        subroutine movef
!
        use storage
#ifdef OFFLOAD_KERNEL_SYNTAX
        use omp_lib
#endif
!
        implicit none

        integer :: i,j,k
#ifdef OFFLOAD_KERNEL_SYNTAX
#  ifdef KS_BLOCK_3D
        integer :: ntx, nty, ntz, nbx, nby, nbz
#  else
        integer :: tid, ncell, nthrd, nblck
#  endif
#endif
!
!------------------------------------------------------
! Best (?) decomposition for BW
!------------------------------------------------------
#ifdef OFFLOAD_KERNEL_SYNTAX
#  ifdef KS_BLOCK_3D
! 3D thread block (8x8x4) + 3D teams grid: one thread per cell, mapped to
! (i,j,k) via dims 0/1/2. The block dimensions are carried by thread_limit
! (num_threads only honors dim 0 on this toolchain); num_teams carries the grid.
        ntx = 8
        nty = 8
        ntz = 4
        nbx = (l + ntx - 1)/ntx
        nby = (m + nty - 1)/nty
        nbz = (n + ntz - 1)/ntz
!$omp target teams parallel thread_limit(dims(3):ntx,nty,ntz) num_teams(dims(3):nbx,nby,nbz)
        i = 1 + omp_get_thread_num_dim(0) + omp_get_team_num_dim(0)*omp_get_num_threads_dim(0)
        j = 1 + omp_get_thread_num_dim(1) + omp_get_team_num_dim(1)*omp_get_num_threads_dim(1)
        k = 1 + omp_get_thread_num_dim(2) + omp_get_team_num_dim(2)*omp_get_num_threads_dim(2)
        if (i <= l .and. j <= m .and. k <= n) then
#  else
! 1D flat thread block (256x1x1): global id decomposed into (i,j,k)
        nthrd = 256
        ncell = l*m*n
        nblck = (ncell + nthrd - 1)/nthrd
!$omp target teams parallel num_threads(dims(3):nthrd) num_teams(dims(3):nblck) thread_limit(nthrd)
        tid = omp_get_thread_num_dim(0) + omp_get_team_num_dim(0)*omp_get_num_threads_dim(0)
        if (tid < ncell) then
          i = mod(tid, l) + 1
          j = mod(tid/l, m) + 1
          k = tid/(l*m) + 1
#  endif
#elif defined(OFFLOAD)
!$OMP target teams distribute parallel do simd collapse(3)
        do k=1,n
        do j=1,m
        do i=1,l
#elif OPENACC
!$acc parallel 
!$acc loop independent collapse(3)
        do k=1,n
        do j=1,m
        do i=1,l
#else
        do concurrent (k=1:n,j=1:m,i=1:l)
#endif
                  b01(i,j,k) = a01(i-1,j+1,k  )
                  b02(i,j,k) = a02(i-1,j  ,k+1)
                  b03(i,j,k) = a03(i-1,j-1,k  )
                  b04(i,j,k) = a04(i-1,j  ,k-1)
                  b05(i,j,k) = a05(i-1,j  ,k  )
                  b06(i,j,k) = a06(i  ,j  ,k-1)
                  b07(i,j,k) = a07(i  ,j-1,k-1)
                  b08(i,j,k) = a08(i  ,j-1,k  )
                  b09(i,j,k) = a09(i  ,j-1,k+1)
                  b10(i,j,k) = a10(i+1,j+1,k  )
                  b11(i,j,k) = a11(i+1,j  ,k+1)
                  b12(i,j,k) = a12(i+1,j-1,k  )
                  b13(i,j,k) = a13(i+1,j,  k-1)
                  b14(i,j,k) = a14(i+1,j  ,k  )
                  b15(i,j,k) = a15(i  ,j  ,k+1)
                  b16(i,j,k) = a16(i  ,j+1,k+1)
                  b17(i,j,k) = a17(i  ,j+1,k  )
                  b18(i,j,k) = a18(i  ,j+1,k-1)
#ifdef OFFLOAD_KERNEL_SYNTAX
        end if
!$omp end target teams parallel
#elif defined(OFFLOAD)
        enddo
        enddo
        enddo
!$OMP end target teams distribute parallel do simd
#elif OPENACC
        enddo
        enddo
        enddo
!$acc end parallel
#else
        enddo
#endif
!
#ifdef DEBUG_2
        if(myrank == 0) then
           write(6,*) "DEBUG2: Exiting from sub. move", & 
                   b01(l/2+1,m,k) , a01(l/2,m1,k)
        endif
#endif
!
        end subroutine movef
