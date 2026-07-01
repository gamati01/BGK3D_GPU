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
#ifdef OFFLOAD_KERNEL_SYNTAX
        use omp_lib
#endif
!
        implicit none
!
        integer:: i,j,k
#ifdef OFFLOAD_KERNEL_SYNTAX
        integer :: ni, nj, nk
#  ifdef KS_BLOCK_3D
        integer :: ntx, nty, ntz, nbx, nby, nbz
#  else
        integer :: tid, ncell, nthrd, nblck
#  endif
#endif
!
        call SYSTEM_CLOCK(countO0, count_rate, count_max)
        call time(tcountO0)
!
#ifdef OFFLOAD_KERNEL_SYNTAX
        ni = imax - imin + 1
        nj = jmax - jmin + 1
        nk = kmax - kmin + 1
#  ifdef KS_BLOCK_3D
        ntx = 8
        nty = 8
        ntz = 4
        nbx = (ni + ntx - 1)/ntx
        nby = (nj + nty - 1)/nty
        nbz = (nk + ntz - 1)/ntz
!$omp target teams parallel thread_limit(dims(3):ntx,nty,ntz) num_teams(dims(3):nbx,nby,nbz)
        i = imin + omp_get_thread_num_dim(0) + omp_get_team_num_dim(0)*omp_get_num_threads_dim(0)
        j = jmin + omp_get_thread_num_dim(1) + omp_get_team_num_dim(1)*omp_get_num_threads_dim(1)
        k = kmin + omp_get_thread_num_dim(2) + omp_get_team_num_dim(2)*omp_get_num_threads_dim(2)
        if (i <= imax .and. j <= jmax .and. k <= kmax) then
#  else
        nthrd = 256
        ncell = ni*nj*nk
        nblck = (ncell + nthrd - 1)/nthrd
!$omp target teams parallel num_threads(dims(3):nthrd) num_teams(dims(3):nblck) thread_limit(nthrd)
        tid = omp_get_thread_num_dim(0) + omp_get_team_num_dim(0)*omp_get_num_threads_dim(0)
        if (tid < ncell) then
          i = imin + mod(tid, ni)
          j = jmin + mod(tid/ni, nj)
          k = kmin + tid/(ni*nj)
#  endif
#elif defined(OFFLOAD)
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
#ifdef OFFLOAD_KERNEL_SYNTAX
        end if
!$omp end target teams parallel
#elif defined(OFFLOAD)
        end do
        end do
        end do
!$OMP end target teams distribute parallel do simd
#elif OPENACC
        end do
        end do
        end do
!$acc end parallel
#else
        end do
#endif
!
! stop timing
        call time(tcountO1)
        call SYSTEM_CLOCK(countO1, count_rate, count_max)
        time_obs = time_obs + real(countO1-countO0)/(count_rate)
        time_obs1 = time_obs1 + (tcountO1-tcountO0)
!
#ifdef DEBUG_2
        if(myrank == 0) then
           write(6,*) "DEBUG2: Exiting from sub. bcond_obs", &
                               imin,imax,jmin,imax
        endif
#endif
!
        end subroutine bcond_obs
