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
!
        implicit none
!
        integer      :: i,j,k
        real(mykind) :: cte1
!
#ifdef CHANNEL
!
# ifdef NOSHIFT
       cte1 = uno
# else
       cte1 = zero
# endif

!
! start timing...
        call SYSTEM_CLOCK(countA0, count_rate, count_max)
        call time(tcountA0)
!
! ----------------------------------------------
! loop foused for performance reason (for GPU)
! -------------------------------------------------------------
!
! ----------------------------------------------
! rear (x = 1)  ! front (x = n)
! ----------------------------------------------
!
#ifdef OFFLOAD
!$OMP target teams distribute parallel do simd
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
        end do
#ifdef OFFLOAD
        end do
!$OMP end target teams distribute parallel do simd
#elif OPENACC
        end do
!$acc end parallel
#endif
!
! ----------------------------------------------
! left (y = 1)  ! right (y = m) 
! ----------------------------------------------
!
#ifdef OFFLOAD
!$OMP target teams distribute parallel do simd 
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
        enddo
#ifdef OFFLOAD
        enddo
!$OMP end target teams distribute parallel do simd
#elif OPENACC
        enddo
!$acc end parallel
#endif
!
! ----------------------------------------------
! down (z = 1)  ! top (z = n) 
! ----------------------------------------------
!
#ifdef OFFLOAD
!$OMP target teams distribute parallel do simd
        do j=0,m+1
        do i=0,l+1
#elif OPENACC
!$acc parallel
!$acc loop independent collapse(2)
        do j=0,m+1
        do i=0,l+1
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
        enddo
#ifdef OFFLOAD
        enddo
!$OMP end target teams distribute parallel do simd
#elif OPENACC
        enddo
!$acc end parallel
#endif

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
           write(6,*) "DEBUG2: Exiting from sub. bcond_channel"
        endif
#endif
        end subroutine bcond_channel
