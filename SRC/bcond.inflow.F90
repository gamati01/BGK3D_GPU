!=====================================================================
!     ****** LBE/bcond_inflow
!
!     COPYRIGHT
!       (c) 2000-2011 by CASPUR/G.Amati
!       (c) 2013-20?? by CINECA/G.Amati
!     NAME
!       bcond_inflow: simple flow with
!                      * inflow (rear)
!                      * outflow (front)
!                      * periodic b.c. (left/right)
!     DESCRIPTION
!       2D periodic bc
!     INPUTS
!       none
!     OUTPUT
!       none
!     TODO
!       remove hardwritten inflow velocity       
!       
!     NOTES
!       Order of upgrading bc
!       1) front  (x = l)
!       2) rear   (x = 1)
!       3) left   (y = 1)
!       4) right  (y = m)
!       4) bottom (y = 1)
!       4) up     (y = n)
!
!       inflow velocity  is hardwritten : 0.1       
!
!     *****
!=====================================================================
!
        subroutine bcond_inflow
!
        use timing
        use storage
!
        implicit none
!
        integer      :: i,j,k
        real(mykind) :: xj,yj,zj
        real(mykind) :: cx01,cx02,cx03,cx04,cx05
        real(mykind) :: cx10,cx11,cx12,cx13,cx14
        real(mykind) :: cvsq,crho,rhoinv
        real(mykind) :: cte1
!
#ifdef INFLOW
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
        u_inflow=0.1
!        
! ----------------------------------------------
! loop fused for perfomance reason (on GPU)        
! ----------------------------------------------
! front, outflow  (x = l)
# ifdef OFFLOAD
!$OMP target teams distribute parallel do simd collapse(2)
        do k=0,n+1
        do j=0,m+1
# elif OPENACC
!$acc parallel 
!$acc loop independent collapse(2)
        do k=0,n+1
        do j=0,m+1
# else
        do concurrent (k=0:n+1,j=0:m+1)
# endif
           crho  =uno
           rhoinv=uno
!           
! front, outflow  (x = l)
           xj = +( a01(l,j,k)+a02(l,j,k)+a03(l,j,k) &
                  +a04(l,j,k)+a05(l,j,k)            &
                  -a10(l,j,k)-a11(l,j,k)-a12(l,j,k) &
                  -a13(l,j,k)-a14(l,j,k) )*rhoinv
!
           yj = +( a03(l,j,k)+a07(l,j,k)+a08(l,j,k) &
                  +a09(l,j,k)+a12(l,j,k)            &
                  -a01(l,j,k)-a10(l,j,k)-a16(l,j,k) &
                  -a17(l,j,k)-a18(l,j,k) )*rhoinv
!
           zj = +( a04(l,j,k)+a06(l,j,k)+a07(l,j,k) &
                  +a13(l,j,k)+a18(l,j,k)            &
                  -a02(l,j,k)-a09(l,j,k)-a11(l,j,k) &
                  -a15(l,j,k)-a16(l,j,k) )*rhoinv
!
           cvsq=xj*xj+yj*yj+zj*zj
!
           cx10 = rf*(-xj-yj   )+qf*(3.0*(xj+yj)*(xj+yj)-cvsq)
           cx11 = rf*(-xj   -zj)+qf*(3.0*(xj+zj)*(xj+zj)-cvsq)
           cx12 = rf*(-xj+yj   )+qf*(3.0*(xj-yj)*(xj-yj)-cvsq)
           cx13 = rf*(-xj   +zj)+qf*(3.0*(xj-zj)*(xj-zj)-cvsq)
           cx14 = rf*(-xj      )+qf*(3.0*(xj   )*(xj   )-cvsq)
!
           a10(l1,j,k) = crho*p2*(cte1+cx10)
           a11(l1,j,k) = crho*p2*(cte1+cx11)
           a12(l1,j,k) = crho*p2*(cte1+cx12)
           a13(l1,j,k) = crho*p2*(cte1+cx13)
           a14(l1,j,k) = crho*p1*(cte1+cx14)
!
! rear, inflow (x = 0)
           xj = u_inflow
           yj = zero
           zj = zero
!           
           cvsq=xj*xj+yj*yj+zj*zj
!
           cx01 = rf*(xj-yj   )+qf*(3.0*(xj-yj)*(xj-yj)-cvsq)
           cx02 = rf*(xj   -zj)+qf*(3.0*(xj-zj)*(xj-zj)-cvsq)
           cx03 = rf*(xj+yj   )+qf*(3.0*(xj+yj)*(xj+yj)-cvsq)
           cx04 = rf*(xj   +zj)+qf*(3.0*(xj+zj)*(xj+zj)-cvsq)
           cx05 = rf*(xj      )+qf*(3.0*(xj   )*(xj   )-cvsq)
!
           a01(0,j,k) = crho*p2*(cte1+cx01)
           a02(0,j,k) = crho*p2*(cte1+cx02)
           a03(0,j,k) = crho*p2*(cte1+cx03)
           a04(0,j,k) = crho*p2*(cte1+cx04)
           a05(0,j,k) = crho*p1*(cte1+cx05)
        end do
# ifdef OFFLOAD
        end do
!$OMP end target teams distribute parallel do simd
# elif OPENACC
        end do
!$acc end parallel
# endif
!
! ----------------------------------------------
! left (y = 1)  ! right (y = m)
! ----------------------------------------------
!
# ifdef OFFLOAD
!$OMP target teams distribute parallel do simd collapse(2)
        do k=0,n+1
        do i=0,l+1
# elif OPENACC
!$acc parallel 
!$acc loop independent collapse(2)
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
! ----------------------------------------------
! bottom (z = 1)  ! up (z = m)
! ----------------------------------------------
!
# ifdef OFFLOAD
!$OMP target teams distribute parallel do simd collapse(2)
        do j=0,m+1
        do i=0,l+1
# elif OPENACC
!$acc parallel 
!$acc loop independent collapse(2)
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
           write(6,*) "DEBUG2: Exiting from sub. bcond_inflow", & 
           u_inflow
        endif
#endif
!        
        end subroutine bcond_inflow
