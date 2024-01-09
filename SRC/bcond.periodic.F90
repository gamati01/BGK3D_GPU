!=====================================================================
!     ****** LBE/bcond_periodic
!
!     COPYRIGHT
!       (c) 2000-2011 by CASPUR/G.Amati
!       (c) 2013-20?? by CINECA/G.Amati
!     NAME
!       bcond
!     DESCRIPTION
!       3D periodic bc
!     INPUTS
!       none
!     OUTPUT
!       none
!     TODO
!       
!     NOTES
!       Order of upgrading bc
!       1) front  (x = l)
!       2) rear   (x = 1)
!       3) left   (y = 1)
!       4) right  (y = m)
!       5) bottom (z = 1)
!       6) top    (z = n)
!
!     *****
!=====================================================================
!
        subroutine bcond_periodic
!
        use timing
        use storage
!
        implicit none
!
        integer      :: i,j,k
!
#ifdef PERIODIC
!
! start timing...
        call SYSTEM_CLOCK(countA0, count_rate, count_max)
        call time(tcountA0)
!
! ----------------------------------------------
! front (x = l) ! rear (x = 1)
! ----------------------------------------------
!
# ifdef OFFLOAD
!$OMP target teams distribute parallel do simd
        do k=0,n+1
        do j=0,m+1
# elif OPENACC
!$acc parallel
!$acc loop independent
        do k=0,n+1
        do j=0,m+1
# else
        do concurrent (k=0:n+1,j=0:m+1)
# endif
! front (x = l)
           a01(0,j,k)  = a01(l,j,k)
           a02(0,j,k)  = a02(l,j,k)
           a03(0,j,k)  = a03(l,j,k)
           a04(0,j,k)  = a04(l,j,k)
           a05(0,j,k)  = a05(l,j,k)

! rear (x = 1)
           a10(l1,j,k) = a10(1,j,k)
           a11(l1,j,k) = a11(1,j,k)
           a12(l1,j,k) = a12(1,j,k)
           a13(l1,j,k) = a13(1,j,k)
           a14(l1,j,k) = a14(1,j,k)
        enddo
# ifdef OFFLOAD
        enddo
!$OMP end target teams distribute parallel do simd
# elif OPENACC
        enddo
!$acc end parallel
# endif
!
! ----------------------------------------------
! left (y = 1)  ! right (y = m) 
! ----------------------------------------------
!
# ifdef OFFLOAD
!$OMP target teams distribute parallel do simd 
        do k=0,n+1
        do i=0,l+1
# elif OPENACC
!$acc parallel
!$acc loop independent
        do k=0,n+1
        do i=0,l+1
# else
        do concurrent (k=0:n+1,i=0:l+1)
# endif
! left (y = 1)  
           a03(i,0,k) = a03(i,m,k)
           a07(i,0,k) = a07(i,m,k)
           a08(i,0,k) = a08(i,m,k)
           a09(i,0,k) = a09(i,m,k)
           a12(i,0,k) = a12(i,m,k)
!
! right (y = m) 
           a01(i,m1,k) = a01(i,1,k)
           a10(i,m1,k) = a10(i,1,k)
           a16(i,m1,k) = a16(i,1,k)
           a17(i,m1,k) = a17(i,1,k)
           a18(i,m1,k) = a18(i,1,k)
        enddo
# ifdef OFFLOAD
        enddo
!$OMP end target teams distribute parallel do simd
# elif OPENACC
        enddo
!$acc end parallel
# endif
!
!
! ----------------------------------------------
! bottom (z = 1)  ! up (z = m) 
! ----------------------------------------------
!
# ifdef OFFLOAD
!$OMP target teams distribute parallel do simd 
        do j=0,m+1
        do i=0,l+1
# elif OPENACC
!$acc parallel
!$acc loop independent
        do j=0,m+1
        do i=0,l+1
# else
        do concurrent (j=0:m+1,i=0:l+1)
# endif
! bottom (z = 1)  
           a04(i,j,0) = a04(i,j,n)
           a06(i,j,0) = a06(i,j,n)
           a07(i,j,0) = a07(i,j,n)
           a13(i,j,0) = a13(i,j,n)
           a18(i,j,0) = a18(i,j,n)
!
! up (z = n) 
           a02(i,j,n1)  = a02(i,j,1)
           a09(i,j,n1)  = a09(i,j,1)
           a11(i,j,n1)  = a11(i,j,1)
           a15(i,j,n1)  = a15(i,j,1)
           a16(i,j,n1)  = a16(i,j,1)
        enddo
# ifdef OFFLOAD
        enddo
!$OMP end target teams distribute parallel do simd
# elif OPENACC
        enddo
!$acc end parallel
# endif
!
! ----------------------------------------------
! stop timing
        call time(tcountA1)
        call SYSTEM_CLOCK(countA1, count_rate, count_max)
        time_bc = time_bc + real(countA1-countA0)/(count_rate)
        time_bc1 = time_bc1 + (tcountA1-tcountA0)
!
#endif
!
#ifdef DEBUG_2
        if(myrank == 0) then
           write(6,*) "DEBUG2: Exiting from sub. bcond_periodic"
        endif
#endif
!       
        end subroutine bcond_periodic
