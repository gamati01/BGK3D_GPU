!=======================================================================
!     ****** LBE/varm
!
!     COPYRIGHT
!       (c) 2000-2011 by CASPUR/G.Amati
!       (c) 2013-20?? by CINECA/G.Amati
!     NAME
!       varm
!     DESCRIPTION
!       Diagnostic subroutine:
!       mean profile of macroscopic variables (in z direction)
!       write on unit 62 (u_med.dat)
!     INPUTS
!       time --> timestep
!     OUTPUT
!       none
!     TODO
!
!     NOTES
!       integer variables used: i,k, itime
!       real variables used: rho, rhoinv, rvol
!                            u(1:n), v(1:n), den(1:n)
!
!     *****
!=======================================================================
!
      subroutine varm(itime)
!
      use storage
      use timing
!
      implicit none
!
      integer :: i,j,k
      integer, INTENT(in) :: itime
!
      real(mykind) :: u(1:n)     ! mean velocity profiles
      real(mykind) :: w(1:n)     ! mean velocity profiles
      real(mykind) :: v(1:n)     ! mean velocity profiles
      real(mykind) :: den(1:n)   ! mean density profile
      real(mykind) :: rho, rhoinv, rvol
      real(mykind) :: cte1
!
#ifdef NOSHIFT
       cte1 = zero
#else
       cte1 = uno
#endif
!
      do j = 1, n
         u(j)  = zero
         w(j)  = zero
         den(j)= zero
      enddo
!      
      do k = 1,n
      do j = 1,m
         do i = 1, l
            rho = +a01(i,j,k)+a02(i,j,k)+a03(i,j,k) &
                  +a04(i,j,k)+a05(i,j,k)+a06(i,j,k) &
                  +a07(i,j,k)+a08(i,j,k)+a09(i,j,k) &
                  +a10(i,j,k)+a11(i,j,k)+a12(i,j,k) &
                  +a13(i,j,k)+a14(i,j,k)+a15(i,j,k) &
                  +a16(i,j,k)+a17(i,j,k)+a18(i,j,k) &
                  +a19(i,j,k)+cte1

!
            rhoinv = uno/rho
!            
            den(j) = den(j) + rho 
!
                u(k) = u(k) &
                       +( a01(i,j,k)+a02(i,j,k)+a03(i,j,k) &
                         +a04(i,j,k)+a05(i,j,k) &
                         -a10(i,j,k)-a11(i,j,k)-a12(i,j,k) &
                         -a13(i,j,k)-a14(i,j,k) ) * rhoinv
!
                w(k) = w(k) &
                       +( a03(i,j,k)+a07(i,j,k)+a08(i,j,k) &
                         +a09(i,j,k)+a12(i,j,k) &
                         -a01(i,j,k)-a10(i,j,k)-a16(i,j,k) &
                         -a17(i,j,k)-a18(i,j,k) )*rhoinv
!
                v(k) = v(k) &
                       +( a04(i,j,k)+a06(i,j,k)+a07(i,j,k) &
                         +a13(i,j,k)+a18(i,j,k) &
                         -a02(i,j,k)-a09(i,j,k)-a11(i,j,k) &
                         -a15(i,j,k)-a16(i,j,k) )*rhoinv
!
         end do
         end do
      end do
!
      rvol = 1.0/(float(l))
!
      write(62,1005) itime
      do j=1,m
         write(62,1003) j+offset(2), & 
              u(j)*rvol ,w(j)*rvol, den(j)*rvol
      end do
      write(62,'(a1)')
      write(62,'(a1)')
!
#ifdef DEBUG_1
      if(myrank == 0) then
         write(6,*) "DEBUG1: Exiting from sub. varm"
      endif
#endif
!
# ifdef MEM_CHECK
      if(myrank == 0) then
         mem_stop = get_mem();
         write(6,*) "MEM_CHECK: after sub. varm mem =", mem_stop
      endif
# endif
!
1003    format(i5,3(e14.6,1x))
1005    format("# t=",i7)
!
       end subroutine varm
