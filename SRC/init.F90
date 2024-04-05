!=======================================================================
!     ****** LBE/init
!
!     COPYRIGHT
!       (c) 2000-2011 by CASPUR/G.Amati
!       (c) 2013-20?? by CINECA/G.Amati
!     NAME
!       init
!     DESCRIPTION
!       initial condition for 3d bgk simulation
!       write on unit 76 (fort.76) x & z coordinates for check
!     INPUTS
!      	rest flow  (default)
!      	Kida vortices
!      	TG vortex (2D)
!      	TG vortex (3D)
!     OUTPUT
!       none
!     TODO
!
!     NOTES
!       integer variables used: i,j.k,opt
!       real variables used: x,z,xj,zj,vsq,rho,pi
!                            x02,x04,x05,x06,x11,x13,x14,x15,x19
!       u = velocity along x direction (i, streamwise)
!       v = velocity along z direction (k, normal-to-wall)
!
!     *****
!=======================================================================
!
        subroutine init(opt)
!
        use storage
        implicit none
!
        integer             :: i,j,k
        integer, INTENT(in) :: opt
!
#ifdef HALF_P
        real(sp) ::  x,y,z,xj,yj,zj
        real(sp) ::  cvsq,crho,pi
        real(sp) ::  cx01,cx02,cx03,cx04,cx05,cx06
        real(sp) ::  cx07,cx08,cx09,cx10,cx11,cx12
        real(sp) ::  cx13,cx14,cx15,cx16,cx17,cx18
        real(sp) ::  cx19
        real(sp) ::  zstart, ystart, xstart
        real(sp) ::  cte1
#else
!
        real(mykind) ::  x,y,z,xj,yj,zj
        real(mykind) ::  cvsq,crho,pi
        real(mykind) ::  cx01,cx02,cx03,cx04,cx05,cx06
        real(mykind) ::  cx07,cx08,cx09,cx10,cx11,cx12
        real(mykind) ::  cx13,cx14,cx15,cx16,cx17,cx18
        real(mykind) ::  cx19
        real(mykind) ::  zstart, ystart, xstart
        real(mykind) ::  cte1
#endif
!
        integer            :: ll, mm, nn
        integer, parameter :: kappa=4       ! # armoniche
!
        parameter(pi=3.141592653589793238462643383279)
!
! check parameter opt
        if((opt.lt.0).or.(opt.gt.5)) then
           write(6,*) "Initial condition out of range[0,5]",opt
           stop
        endif
!
#ifdef NOSHIFT
       cte1=uno
#else
       cte1=zero
#endif
!
        ll = l
        mm = m
        nn = n
!
! builds the populations
        crho = 1.0_mykind
!
        do k = 0, n1
           z = (real(k+zstart) - 0.5 - 0.5*real(nn))/(0.5*real(nn))
!           write(76,*) k, z
        enddo
!        write(76,*) " "
!        write(76,*) " "
        do j = 0, m1
           y = (real(j,mykind)-0.5d0)/real(mm,mykind)  ! 0<x<1 (taylor)
!           write(76,*) j, y
        enddo
!        write(76,*) " "
!        write(76,*) " "
        do i = 0, l1
           x = (real(i,mykind)-0.5d0)/real(ll,mykind)! 0<x<1 (taylor)
!           write(76,*) i, x
        enddo
!        write(76,*) " "
!        write(76,*) " "
!
! rest flow 
        xj = 0.0
        yj = 0.0
        zj = 0.0
!
#ifdef TGV3D
        write(6,*) "INF0: kappa =", kappa
#endif
!        
        do k = 0, n1
#ifdef PERIODIC
           z = (real(k,mykind)-0.5d0)/real(nn,mykind)  ! 0<z<1 
#endif
           do j = 0, m1
#ifdef PERIODIC
              y = (real(j,mykind)-0.5d0)/real(mm,mykind)  ! 0<y<1 
#endif
              do i = 0, l1
#ifdef PERIODIC
                 x = (real(i,mykind)-0.5d0)/real(ll,mykind)! 0<x<1 
#endif
!                 
#ifdef PERIODIC
! default TaylorGreen (2D)
                 xj = 0.1d0*sin(real(2,mykind)*pi*x)*cos(real(2,mykind)*pi*z)
                 yj = zero
                 zj =-0.1d0*cos(real(2,mykind)*pi*x)*sin(real(2,mykind)*pi*z)
# ifdef KVX
! if Kida overwrite....                 
                 xj = 0.1d0*sin(2*pi*x)*    & 
                      (cos(6*pi*y)*cos(2*pi*z)-cos(2*pi*y)*cos(6*pi*z))
                 yj = 0.1d0*sin(2*pi*y)*    & 
                      (cos(6*pi*z)*cos(2*pi*x)-cos(2*pi*z)*cos(6*pi*x))
                 zj =-0.1d0*sin(2*pi*z)*    & 
                      (cos(6*pi*x)*cos(2*pi*y)-cos(2*pi*x)*cos(6*pi*y))
              stop
# elif TGV3D
! Article by Emerson et al.              
!                 xj = 0.1d0*sin(2*pi*x)*cos(2*pi*y)*cos(2*pi*z)
!                 yj =-0.1d0*cos(2*pi*x)*sin(2*pi*y)*cos(2*pi*z)
!                 zj = zero

!                 crho = uno + ((0.1*0.1)/16.0)* & 
!                         (cos(4*pi*x)+cos(4*pi*y))*(cos(4*pi*z)+2)
!
!
! Kajzer et al. (2014)              
                 xj = 0.10d0*cos(2*kappa*pi*x)*sin(2*kappa*pi*y)*sin(2*kappa*pi*z)
                 yj =-0.05d0*sin(2*kappa*pi*x)*cos(2*kappa*pi*y)*sin(2*kappa*pi*z)
                 zj =-0.05d0*sin(2*kappa*pi*x)*sin(2*kappa*pi*y)*cos(2*kappa*pi*z)
# endif
#endif
!
                 cvsq=xj*xj+yj*yj+zj*zj
!
                 cx01 = rf*( xj-yj   )+qf*(3.0*(xj-yj)*(xj-yj)-cvsq)
                 cx02 = rf*( xj   -zj)+qf*(3.0*(xj-zj)*(xj-zj)-cvsq)
                 cx03 = rf*( xj+yj   )+qf*(3.0*(xj+yj)*(xj+yj)-cvsq)
                 cx04 = rf*( xj   +zj)+qf*(3.0*(xj+zj)*(xj+zj)-cvsq)
                 cx05 = rf*( xj      )+qf*(3.0*(xj   )*(xj   )-cvsq)
                 cx06 = rf*(       zj)+qf*(3.0*(zj   )*(zj   )-cvsq)
                 cx07 = rf*(    yj+zj)+qf*(3.0*(yj+zj)*(yj+zj)-cvsq)
                 cx08 = rf*(    yj   )+qf*(3.0*(yj   )*(yj   )-cvsq)
                 cx09 = rf*(    yj-zj)+qf*(3.0*(yj-zj)*(yj-zj)-cvsq)
                 cx10 = rf*(-xj-yj   )+qf*(3.0*(xj+yj)*(xj+yj)-cvsq)
                 cx11 = rf*(-xj   -zj)+qf*(3.0*(xj+zj)*(xj+zj)-cvsq)
                 cx12 = rf*(-xj+yj   )+qf*(3.0*(xj-yj)*(xj-yj)-cvsq)
                 cx13 = rf*(-xj   +zj)+qf*(3.0*(xj-zj)*(xj-zj)-cvsq)
                 cx14 = rf*(-xj      )+qf*(3.0*(xj   )*(xj   )-cvsq)
                 cx15 = rf*(      -zj)+qf*(3.0*(zj   )*(zj   )-cvsq)
                 cx16 = rf*(   -yj-zj)+qf*(3.0*(yj+zj)*(yj+zj)-cvsq)
                 cx17 = rf*(   -yj   )+qf*(3.0*(yj   )*(yj   )-cvsq)
                 cx18 = rf*(   -yj+zj)+qf*(3.0*(yj-zj)*(yj-zj)-cvsq)
                 cx19 = rf*(   0.0   )+qf*(3.0*( 0.0 )*( 0.0 )-cvsq)

                 a01(i,j,k) = crho*p2*(cte1+cx01)
                 a02(i,j,k) = crho*p2*(cte1+cx02)
                 a03(i,j,k) = crho*p2*(cte1+cx03)
                 a04(i,j,k) = crho*p2*(cte1+cx04)
                 a05(i,j,k) = crho*p1*(cte1+cx05)
                 a06(i,j,k) = crho*p1*(cte1+cx06)
                 a07(i,j,k) = crho*p2*(cte1+cx07)
                 a08(i,j,k) = crho*p1*(cte1+cx08)
                 a09(i,j,k) = crho*p2*(cte1+cx09)
                 a10(i,j,k) = crho*p2*(cte1+cx10)
                 a11(i,j,k) = crho*p2*(cte1+cx11)
                 a12(i,j,k) = crho*p2*(cte1+cx12)
                 a13(i,j,k) = crho*p2*(cte1+cx13)
                 a14(i,j,k) = crho*p1*(cte1+cx14)
                 a15(i,j,k) = crho*p1*(cte1+cx15)
                 a16(i,j,k) = crho*p2*(cte1+cx16)
                 a17(i,j,k) = crho*p1*(cte1+cx17)
                 a18(i,j,k) = crho*p2*(cte1+cx18)
                 a19(i,j,k) = crho*p0*(cte1+cx19)
!
              end do
           end do
        end do
!
#ifdef DEBUG_1
        if(myrank == 0) then
           write(6,*) "DEBUG1: Exiting from sub. init"
        endif
#endif
!
        end subroutine init
