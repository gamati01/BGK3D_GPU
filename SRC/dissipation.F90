!====================================================
!     ****** LBE/dissipation
!
!     COPYRIGHT
!       (c) 2000-2011 by CASPUR/G.Amati
!       (c) 2013-20?? by CINECA/G.Amati
!     NAME
!       dissipation
!     DESCRIPTION
!       diagnostic subroutine:
!       computing dissipation, turbulent kinetic energy and othe stuff
!       for HIT test case
!       using LBM stress tensor
!     INPUTS
!       itime --> timestep 
!     TODO
!
!     NOTES
!
!     *****
!====================================================
!
        subroutine dissipation(itime)
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
        real(mykind) :: x,y,z
        real(mykind) :: pi
        real(mykind) :: x01,x02,x03,x04,x05,x06,x07
        real(mykind) :: x08,x09,x10,x11,x12,x13
        real(mykind) :: x14,x15,x16,x17,x18,x19
        real(mykind) :: e01,e02,e03,e04,e05,e06,e07
        real(mykind) :: e08,e09,e10,e11,e12,e13
        real(mykind) :: e14,e15,e16,e17,e18,e19
        real(mykind) :: n01,n02,n03,n04,n05,n06,n07
        real(mykind) :: n08,n09,n10,n11,n12,n13
        real(mykind) :: n14,n15,n16,n17,n18
        real(mykind) :: rho,rhoinv,vx,vy,vz
        real(mykind) :: vx2,vy2,vz2,vsq
        real(mykind) :: rp1,rp2,rp0
        real(mykind) :: vxpz,vxmz
        real(mykind) :: vxpy,vxmy
        real(mykind) :: vypz,vymz
        real(mykind) :: qx,qy,qz,q0
        real(mykind) :: qxpz,qxmz,qxpy,qxmy,qypz,qymz
        real(mykind) :: cte1,cte0
        real(mykind) :: Pxx,Pyy,Pzz,Pxy,Pyx,Pxz,Pzx,Pyz,Pzy
        real(mykind) :: Ptotal
        real(mykind) :: diss, tke, kolmog
        real(mykind) :: dt, dx
!        
!
#ifdef NOSHIFT
        cte1 = zero
#else
        cte1 = uno
#endif
        cte0 = uno - cte1
!
        dx = (2.0*3.141592)/float(l)
        dt = dx*u00
        diss = 0.0
        tke = 0.0
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
           x01 = a01(i-1,j+1,k  )
           x02 = a02(i-1,j  ,k+1)
           x03 = a03(i-1,j-1,k  )
           x04 = a04(i-1,j  ,k-1)
           x05 = a05(i-1,j  ,k  )
           x06 = a06(i  ,j  ,k-1)
           x07 = a07(i  ,j-1,k-1)
           x08 = a08(i  ,j-1,k  )
           x09 = a09(i  ,j-1,k+1)
           x10 = a10(i+1,j+1,k  )
           x11 = a11(i+1,j  ,k+1)
           x12 = a12(i+1,j-1,k  )
           x13 = a13(i+1,j  ,k-1)
           x14 = a14(i+1,j  ,k  )
           x15 = a15(i  ,j  ,k+1)
           x16 = a16(i  ,j+1,k+1)
           x17 = a17(i  ,j+1,k  )
           x18 = a18(i  ,j+1,k-1)
           x19 = a19(i  ,j  ,k  )

           rho = x01 + x02 + x03 + x04 + x05 + x06 + x07 + x08 &
                +x09 + x10 + x11 + x12 + x13 + x14 + x15 + x16 &
                +x17 + x18 + x19 + cte1
!
           rhoinv = uno/rho
!
           vx = (x01+x02+x03+x04+x05-x10-x11-x12-x13-x14)*rhoinv
           vy = (x03+x07+x08+x09+x12-x01-x10-x16-x17-x18)*rhoinv
           vz = (x04+x06+x07+x13+x18-x02-x09-x11-x15-x16)*rhoinv
!           
! Quadratic terms
           vx2 = vx*vx
           vy2 = vy*vy
           vz2 = vz*vz
!
           vsq = vx2+vy2+vz2
!
           tke = tke + 0.5*vsq           
!
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
           n13 = x13-e13
           n14 = x14-e14
           n15 = x15-e15
           n16 = x16-e16
           n17 = x17-e17
           n18 = x18-e18
!
! compute Pij (six terms)
           Pxx = n01 + &
                 n02 + &
                 n03 + &
                 n04 + &
                 n05 + &
                 n10 + &
                 n11 + &
                 n12 + &
                 n13 + &
                 n14 
!
           Pyy = n01 + &
                 n03 + &
                 n07 + &
                 n08 + &
                 n09 + &
                 n10 + &
                 n12 + &
                 n16 + &
                 n17 + &
                 n18 
        !
           Pzz = n02 + &
                 n04 + &
                 n06 + &
                 n07 + &
                 n09 + &
                 n11 + &
                 n13 + &
                 n15 + &
                 n16 + &
                 n18 
        !
           Pxz = -n02 &
                 +n04 &
                 +n11 &
                 -n13 
        !
           Pxy = -n01 &
                 +n03 &
                 +n10 &
                 -n12 
        !
           Pyz = +n07 &
                 -n09 &
                 +n16 &
                 -n18 
        !
           Pyx = Pxy
           Pzx = Pxz
           Pzy = Pyz
        !           
        ! calculate Pi total
           Ptotal =(Pxx)**2 + (Pyy)**2 + (Pzz)**2 + &
                   (2.0*Pxy*Pyx)                  + &
                   (2.0*Pxz*Pzx)                  + &
                   (2.0*Pyz*Pzy)
        !           
           diss = diss + 2.0*svisc*Ptotal
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
! write results (dissipation.dat)
! # timestep, time, tke, dissipation, kolmogorov lenght        
        kolmog = ((svisc**3)/diss/float(l)/float(m)/float(n))**(0.25)
        write(77,1004) itime                          &
                     , float(itime)*dt                &
                     , tke/float(l)/float(m)/float(n)/(u00*u00)  &
                     , diss                           &
                     , kolmog                          
!                            
        flush(77)      
!
1004    format(i8,4(e14.6,1x))
!
#ifdef DEBUG_2
        if(myrank == 0) then
           write(6,*) "DEBUG2: Exiting from sub. dissipation",  &
                      dt, dx, u0, u00
        endif
#endif
        end subroutine dissipation
