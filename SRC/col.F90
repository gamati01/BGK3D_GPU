!====================================================
!     ****** LBE/coll
!
!     COPYRIGHT
!       (c) 2000-2011 by CASPUR/G.Amati
!       (c) 2013-20?? by CINECA/G.Amati
!     NAME
!       coll
!     DESCRIPTION
!       Collision according to bgk style (\omega(f_f^(eq)))
!       forcing is included and is proportional to u_0
!     INPUTS
!       itime --> timestep 
!       none
!     TODO
!
!     NOTES
!       integer variables used: i,k,itime
!       real variables used: 
!
!
!     *****
!====================================================
!
        subroutine col(itime)
!
        use storage
        use real_kinds
        use timing
!
        implicit none
!
        integer, intent(IN) :: itime
        integer             :: i,j,k
!
        real(mykind) :: x01,x02,x03,x04,x05,x06,x07
        real(mykind) :: x08,x09,x10,x11,x12,x13
        real(mykind) :: x14,x15,x16,x17,x18,x19
        real(mykind) :: e01,e02,e03,e04,e05,e06,e07
        real(mykind) :: e08,e09,e10,e11,e12,e13
        real(mykind) :: e14,e15,e16,e17,e18,e19
        real(mykind) :: n01,n02,n03,n04,n05,n06,n07
        real(mykind) :: n08,n09,n10,n11,n12,n13
        real(mykind) :: n14,n15,n16,n17,n18,n19
        real(mykind) :: rho,rhoinv,vx,vy,vz
        real(mykind) :: vx2,vy2,vz2,vsq
        real(mykind) :: rp1,rp2,rp0
        real(mykind) :: vxpz,vxmz
        real(mykind) :: vxpy,vxmy
        real(mykind) :: vypz,vymz
        real(mykind) :: qx,qy,qz,q0
        real(mykind) :: qxpz,qxmz,qxpy,qxmy,qypz,qymz
        real(mykind) :: forcex, forcey, forcez
        real(mykind) :: cte1,cte0
        real(mykind) :: Pxx,Pyy,Pzz,Pxy,Pyx,Pxz,Pzx,Pyz,Pzy
        real(mykind) :: Ptotal,Ts,cteS
!
#ifdef DEBUG_3
        real(mykind) :: cte
        character(len=17) :: file_nameD
        file_nameD = 'debug.xxx.xxx.log'
        write(file_nameD(7:9),3300) itime
        write(file_nameD(11:13),3300) myrank
        open(41,file=file_nameD, status='unknown')        ! debug file
!
!        call probe(itime,(3*l/4),(m/2))
#endif
!
#ifdef NOSHIFT
        cte1 = zero
#else
        cte1 = uno
#endif
        cte0 = uno - cte1
!
! initialize constant.....        
        forcex = zero
        forcey = zero
        forcez = zero
!
#ifdef OFFLOAD
!$OMP target teams distribute parallel do simd collapse(3)
        do k = 1,n
        do j = 1,m
        do i = 1,l
#elif OPENACC
!$acc parallel
!$acc loop independent collapse(3)
        do k = 1,n
        do j = 1,m
        do i = 1,l
#else
        do concurrent (k=1:n,j=1:m,i=1:l)
#endif
!
           x01 = b01(i,j,k)
           x02 = b02(i,j,k)
           x03 = b03(i,j,k)
           x04 = b04(i,j,k)
           x05 = b05(i,j,k)
           x06 = b06(i,j,k)
           x07 = b07(i,j,k)
           x08 = b08(i,j,k)
           x09 = b09(i,j,k)
           x10 = b10(i,j,k)
           x11 = b11(i,j,k)
           x12 = b12(i,j,k)
           x13 = b13(i,j,k)
           x14 = b14(i,j,k)
           x15 = b15(i,j,k)
           x16 = b16(i,j,k)
           x17 = b17(i,j,k)
           x18 = b18(i,j,k)
           x19 = a19(i,j,k)
!
           rho = x01 + x02 + x03 + x04 + x05 + x06 + x07 + x08 &
                +x09 + x10 + x11 + x12 + x13 + x14 + x15 + x16 &
                +x17 + x18 + x19 + cte1
!
           rhoinv = uno/rho
!
           vx = (x01+x02+x03+x04+x05-x10-x11-x12-x13-x14)*rhoinv
           vy = (x03+x07+x08+x09+x12-x01-x10-x16-x17-x18)*rhoinv
           vz = (x04+x06+x07+x13+x18-x02-x09-x11-x15-x16)*rhoinv

#ifdef DEBUG_3
           write(41,3131) i+offset(1),j+offset(2), obs(i,j),   & 
                        vx,vy,rho
!                       x01,x03,x05,x08,x10,x12,x14,x17,x19
#endif
!
! Quadratic terms
           vx2 = vx*vx
           vy2 = vy*vy
           vz2 = vz*vz

           vsq = vx2+vy2+vz2

           vxpy = vx+vy
           qxpy = cte0+qf*(tre*vxpy*vxpy-vsq)
           vxmy = vx-vy
           qxmy = cte0+qf*(tre*vxmy*vxmy-vsq)
           vxpz = vx+vz
           qxpz = cte0+qf*(tre*vxpz*vxpz-vsq)
           vxmz = vx-vz
           qxmz = cte0+qf*(tre*vxmz*vxmz-vsq)
           vypz = vy+vz
           qypz = cte0+qf*(tre*vypz*vypz-vsq)
           vymz = vy-vz
           qymz = cte0+qf*(tre*vymz*vymz-vsq)
           qx   = cte0+qf*(tre*vx2      -vsq)
           qy   = cte0+qf*(tre*vy2      -vsq)
           qz   = cte0+qf*(tre*vz2      -vsq)
           q0   = cte0+qf*(             -vsq)
!
! linear terms
           vx   = rf*vx
           vy   = rf*vy
           vz   = rf*vz
           vxpy = rf*vxpy
           vxmy = rf*vxmy
           vxpz = rf*vxpz
           vxmz = rf*vxmz
           vypz = rf*vypz
           vymz = rf*vymz
!
! constant terms
           rp0 = rho*p0
           rp1 = rho*p1
           rp2 = rho*p2
!
! equilibrium distribution
           e01 = rp2*(+vxmy+qxmy) + cte1*(rp2-p2)
           e02 = rp2*(+vxmz+qxmz) + cte1*(rp2-p2)
           e03 = rp2*(+vxpy+qxpy) + cte1*(rp2-p2)
           e04 = rp2*(+vxpz+qxpz) + cte1*(rp2-p2)
           e05 = rp1*(+vx  +qx  ) + cte1*(rp1-p1)
           e06 = rp1*(+vz  +qz  ) + cte1*(rp1-p1)
           e07 = rp2*(+vypz+qypz) + cte1*(rp2-p2)
           e08 = rp1*(+vy  +qy  ) + cte1*(rp1-p1)
           e09 = rp2*(+vymz+qymz) + cte1*(rp2-p2)
           e10 = rp2*(-vxpy+qxpy) + cte1*(rp2-p2)
           e11 = rp2*(-vxpz+qxpz) + cte1*(rp2-p2)
           e12 = rp2*(-vxmy+qxmy) + cte1*(rp2-p2)
           e13 = rp2*(-vxmz+qxmz) + cte1*(rp2-p2)
           e14 = rp1*(-vx  +qx  ) + cte1*(rp1-p1)
           e15 = rp1*(-vz  +qz  ) + cte1*(rp1-p1)
           e16 = rp2*(-vypz+qypz) + cte1*(rp2-p2)
           e17 = rp1*(-vy  +qy  ) + cte1*(rp1-p1)
           e18 = rp2*(-vymz+qymz) + cte1*(rp2-p2)
           e19 = rp0*(     +q0  ) + cte1*(rp0-p0)
!
#ifdef LES
! compute les
!
!non-equilibrium distribution
           n01 = x01-e01
           n02 = x02-e02
           n03 = x03-e03
           n04 = x04-e04
           n05 = x05-e05
           n06 = x06-e06
           n07 = x07-e07
           n08 = x08-e08
           n09 = x09-e09
           n10 = x10-e10
           n11 = x11-e11
           n12 = x12-e12
           n13 = x13-e03
           n14 = x14-e14
           n15 = x15-e15
           n16 = x16-e16
           n17 = x17-e17
           n18 = x18-e18
           n19 = x19-e19
!
! compute Pij (six terms)           
           Pxx = cx(01)*cx(01)*n01 + &
                 cx(02)*cx(02)*n02 + &
                 cx(03)*cx(03)*n03 + &
                 cx(04)*cx(04)*n04 + &
                 cx(05)*cx(05)*n05 + &
                 cx(06)*cx(06)*n06 + &
                 cx(07)*cx(07)*n07 + &
                 cx(08)*cx(08)*n08 + &
                 cx(09)*cx(09)*n09 + &
                 cx(10)*cx(10)*n10 + &
                 cx(11)*cx(11)*n11 + &
                 cx(12)*cx(12)*n12 + &
                 cx(13)*cx(13)*n13 + &
                 cx(14)*cx(14)*n14 + &
                 cx(15)*cx(15)*n15 + &
                 cx(16)*cx(16)*n16 + &
                 cx(17)*cx(17)*n17 + &
                 cx(18)*cx(18)*n18 + &
                 cx(19)*cx(19)*n19 
!
           Pyy = cy(01)*cy(01)*n01 + &
                 cy(02)*cy(02)*n02 + &
                 cy(03)*cy(03)*n03 + &
                 cy(04)*cy(04)*n04 + &
                 cy(05)*cy(05)*n05 + &
                 cy(06)*cy(06)*n06 + &
                 cy(07)*cy(07)*n07 + &
                 cy(08)*cy(08)*n08 + &
                 cy(09)*cy(09)*n09 + &
                 cy(10)*cy(10)*n10 + &
                 cy(11)*cy(11)*n11 + &
                 cy(12)*cy(12)*n12 + &
                 cy(13)*cy(13)*n13 + &
                 cy(14)*cy(14)*n14 + &
                 cy(15)*cy(15)*n15 + &
                 cy(16)*cy(16)*n16 + &
                 cy(17)*cy(17)*n17 + &
                 cy(18)*cy(18)*n18 + &
                 cy(19)*cy(19)*n19 
!
           Pzz = cz(01)*cz(01)*n01 + &
                 cz(02)*cz(02)*n02 + &
                 cz(03)*cz(03)*n03 + &
                 cz(04)*cz(04)*n04 + &
                 cz(05)*cz(05)*n05 + &
                 cz(06)*cz(06)*n06 + &
                 cz(07)*cz(07)*n07 + &
                 cz(08)*cz(08)*n08 + &
                 cz(09)*cz(09)*n09 + &
                 cz(10)*cz(10)*n10 + &
                 cz(11)*cz(11)*n11 + &
                 cz(12)*cz(12)*n12 + &
                 cz(13)*cz(13)*n13 + &
                 cz(14)*cz(14)*n14 + &
                 cz(15)*cz(15)*n15 + &
                 cz(16)*cz(16)*n16 + &
                 cz(17)*cz(17)*n17 + &
                 cz(18)*cz(18)*n18 + &
                 cz(19)*cz(19)*n19 
!
           Pxz = cx(01)*cz(01)*n01 + &
                 cx(02)*cz(02)*n02 + &
                 cx(03)*cz(03)*n03 + &
                 cx(04)*cz(04)*n04 + &
                 cx(05)*cz(05)*n05 + &
                 cx(06)*cz(06)*n06 + &
                 cx(07)*cz(07)*n07 + &
                 cx(08)*cz(08)*n08 + &
                 cx(09)*cz(09)*n09 + &
                 cx(10)*cz(10)*n10 + &
                 cx(11)*cz(11)*n11 + &
                 cx(12)*cz(12)*n12 + &
                 cx(13)*cz(13)*n13 + &
                 cx(14)*cz(14)*n14 + &
                 cx(15)*cz(15)*n15 + &
                 cx(16)*cz(16)*n16 + &
                 cx(17)*cz(17)*n17 + &
                 cx(18)*cz(18)*n18 + &
                 cx(19)*cz(19)*n19 
!
           Pxy = cx(01)*cy(01)*n01 + &
                 cx(02)*cy(02)*n02 + &
                 cx(03)*cy(03)*n03 + &
                 cx(04)*cy(04)*n04 + &
                 cx(05)*cy(05)*n05 + &
                 cx(06)*cy(06)*n06 + &
                 cx(07)*cy(07)*n07 + &
                 cx(08)*cy(08)*n08 + &
                 cx(09)*cy(09)*n09 + &
                 cx(10)*cy(10)*n10 + &
                 cx(11)*cy(11)*n11 + &
                 cx(12)*cy(12)*n12 + &
                 cx(13)*cy(13)*n13 + &
                 cx(14)*cy(14)*n14 + &
                 cx(15)*cy(15)*n15 + &
                 cx(16)*cy(16)*n16 + &
                 cx(17)*cy(17)*n17 + &
                 cx(18)*cy(18)*n18 + &
                 cx(19)*cy(19)*n19 
!
           Pyz = cy(01)*cz(01)*n01 + &
                 cy(02)*cz(02)*n02 + &
                 cy(03)*cz(03)*n03 + &
                 cy(04)*cz(04)*n04 + &
                 cy(05)*cz(05)*n05 + &
                 cy(06)*cz(06)*n06 + &
                 cy(07)*cz(07)*n07 + &
                 cy(08)*cz(08)*n08 + &
                 cy(09)*cz(09)*n09 + &
                 cy(10)*cz(10)*n10 + &
                 cy(11)*cz(11)*n11 + &
                 cy(12)*cz(12)*n12 + &
                 cy(13)*cz(13)*n13 + &
                 cy(14)*cz(14)*n14 + &
                 cx(15)*cz(15)*n15 + &
                 cy(16)*cz(16)*n16 + &
                 cy(17)*cz(17)*n17 + &
                 cy(18)*cz(18)*n18 + &
                 cy(19)*cz(19)*n19 
!
                 Pyx = Pxy
                 Pzx = Pxz
                 Pzy = Pyz
!           
! calculate Pi total
                 Ptotal =sqrt((Pxx)**2 + (Pyy)**2 + (Pyy)**2 + &
                              (2.0*Pxy*Pyx)                  + &
                              (2.0*Pxz*Pzx)                  + &
                              (2.0*Pyz*Pzy)) 
                 cteS = 0.1
!           
! adding turbulent viscosity
                 Ts = 1/(2*omega1) + sqrt(18*(cteS)**2 *Ptotal+(1/omega1)**2)/2
                 omega = 1/Ts
!
#endif 
!
! forcing term
!
#ifdef CHANNEL
# ifdef FORCING_Y
        forcex = zero
        forcey = fgrad*rho
        forcez = zero
# elif FORCING_Z
        forcex = zero
        forcex = zero
        forcez = fgrad*rho
# else
        forcex = fgrad*rho
        forcey = zero
        forcez = zero
# endif
#endif
!
! loop on populations
        a01(i,j,k) = x01 - omega*(x01-e01) + forcex - forcey         
        a02(i,j,k) = x02 - omega*(x02-e02) + forcex          - forcez
        a03(i,j,k) = x03 - omega*(x03-e03) + forcex + forcey         
        a04(i,j,k) = x04 - omega*(x04-e04) + forcex          + forcez
        a05(i,j,k) = x05 - omega*(x05-e05) + forcex                  
        a06(i,j,k) = x06 - omega*(x06-e06)                   + forcez
        a07(i,j,k) = x07 - omega*(x07-e07)          + forcey + forcez
        a08(i,j,k) = x08 - omega*(x08-e08)          + forcey         
        a09(i,j,k) = x09 - omega*(x09-e09)          + forcey - forcez
        a10(i,j,k) = x10 - omega*(x10-e10) - forcex - forcey         
        a11(i,j,k) = x11 - omega*(x11-e11) - forcex          - forcez
        a12(i,j,k) = x12 - omega*(x12-e12) - forcex + forcey         
        a13(i,j,k) = x13 - omega*(x13-e13) - forcex          + forcez
        a14(i,j,k) = x14 - omega*(x14-e14) - forcex                  
        a15(i,j,k) = x15 - omega*(x15-e15)                   - forcez
        a16(i,j,k) = x16 - omega*(x16-e16)          - forcey - forcez
        a17(i,j,k) = x17 - omega*(x17-e17)          - forcey         
        a18(i,j,k) = x18 - omega*(x18-e18)          - forcey + forcez
        a19(i,j,k) = x19 - omega*(x19-e19)                           

        end do
#ifdef OFFLOAD
        end do
        end do
!$OMP end target teams distribute parallel do simd
#elif OPENACC
        end do
        end do
        !$acc end parallel
#endif
!
#ifdef DEBUG_3
!
! format
3300  format(i3.3)
!3131  format(3(i,1x),6(e14.6,1x))
3131  format(3(i,1x),3(e14.6,1x))
        close(41)
#endif
!
#ifdef DEBUG_2
        if(myrank == 0) then
           write(6,*) "DEBUG2: Exiting from sub. coll", &
                       forcex, forcey, forcez
        endif
#endif
        end subroutine col
