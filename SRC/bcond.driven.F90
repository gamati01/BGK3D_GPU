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
!
        implicit none
!
        integer      :: i,j,k
        real(mykind) :: force
!
! start timing...
        call SYSTEM_CLOCK(countA0, count_rate, count_max)
        call time(tcountA0)
!
        force =  u00/(6.0)
!
! bc. along X direction
!        
#ifdef OFFLOAD
!$OMP target teams distribute parallel do simd 
        do k=1,n
        do j=1,m
#elif OPENACC
!$acc parallel
!$acc loop independent 
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
        end do
#ifdef OFFLOAD
        end do
!$OMP end target teams distribute parallel do simd
#elif OPENACC
        end do
        !$acc end parallel
#endif
!        
! bc. along y direction
!        
#ifdef OFFLOAD
!$OMP target teams distribute parallel do simd 
        do k=1,n
        do i=1,l
#elif OPENACC
!$acc parallel
!$acc loop independent
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
        enddo
#ifdef OFFLOAD
        enddo
!$OMP end target teams distribute parallel do simd
#elif OPENACC
        enddo
        !$acc end parallel
#endif
!
! bc. along z direction
!        
#ifdef OFFLOAD
!$OMP target teams distribute parallel do simd 
        do j=1,m
        do i=1,l
#elif OPENACC
!$acc parallel
!$acc loop independent
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
        enddo
#ifdef OFFLOAD
        enddo
!$OMP end target teams distribute parallel do simd
#elif OPENACC
        enddo
        !$acc end parallel
#endif
!
! stop timing
        call time(tcountA1)
        call SYSTEM_CLOCK(countA1, count_rate, count_max)
        time_bc = time_bc + real(countA1-countA0)/(count_rate)
        time_bc1 = time_bc1 + (tcountA1-tcountA0)
!
#ifdef DEBUG_2
        if(myrank == 0) then
           write(6,*) "DEBUG2: Exiting from sub. bcond_driven", force
        endif
#endif
!        
        end subroutine bcond_driven
