!=====================================================================
!     ****** LBE/dissipation
!
!     COPYRIGHT
!       (c) 2000-2011 by CASPUR/G.Amati
!       (c) 2013-20?? by CINECA/G.Amati
!     NAME
!       diagno
!     DESCRIPTION
!       diagnostic subroutine:
!     INPUTS
!       itime --> timestep
!     OUTPUT
!       none
!     TODO
!
!     NOTES
!       integer variables used: i,j
!       real variables used: rtot,xtot,stot
!                            xj, rho, rhoinv
!
!     *****
!=====================================================================
!
       subroutine dissipation(itime)
!
       use storage
       implicit none
!
       integer, INTENT(in) :: itime
       integer             :: i,j,k
!
       real(mystorage), dimension(0:l+1,0:m+1,0:n+1) :: density
       real(mystorage), dimension(0:l+1,0:m+1,0:n+1) :: vel_u
       real(mystorage), dimension(0:l+1,0:m+1,0:n+1) :: vel_v
       real(mystorage), dimension(0:l+1,0:m+1,0:n+1) :: vel_w
!       
       real(mykind):: diss,tke
       real(mykind):: rho,rhoinv
       real(mykind):: x01,x02,x03,x04,x05,x06,x07
       real(mykind):: x08,x09,x10,x11,x12,x13
       real(mykind):: x14,x15,x16,x17,x18,x19
       real(mykind):: cte1
       real(mykind):: u_y, u_z, v_x, v_z, w_x, w_y
       real(mykind):: mean_u, mean_v, mean_w
       real(mykind):: dissip, turbke, kolmog
!
#ifdef NOSHIFT
       cte1 = zero
#else
       cte1 = uno
#endif
!
       diss = zero
       tke  = zero
       u_y  = zero
       u_z  = zero
       v_x  = zero
       v_z  = zero
       w_x  = zero
       w_y  = zero
       mean_u = zero
       mean_v = zero
       mean_w = zero
!
#ifdef NOMANAGED
!$acc update self(a01,a02,a03,a04,a05,a06,a07,a08,a09,a10, &
!$acc&            a11,a12,a13,a14,a15,a16,a17,a18,a19)
#endif
!
#ifdef OFFLOAD
!$omp target update from(a01,a02,a03,a04,a05,a06,a07,a08,a09,a10, &
!$omp&                   a11,a12,a13,a14,a15,a16,a17,a18,a19)
#endif
       do k=0,n+1
         do j=0,m+1
            do i=0,l+1
               x01 = a01(i,j,k)
               x02 = a02(i,j,k)
               x03 = a03(i,j,k)
               x04 = a04(i,j,k)
               x05 = a05(i,j,k)
               x06 = a06(i,j,k)
               x07 = a07(i,j,k)
               x08 = a08(i,j,k)
               x09 = a09(i,j,k)
               x10 = a10(i,j,k)
               x11 = a11(i,j,k)
               x12 = a12(i,j,k)
               x13 = a13(i,j,k)
               x14 = a14(i,j,k)
               x15 = a15(i,j,k)
               x16 = a16(i,j,k)
               x17 = a17(i,j,k)
               x18 = a18(i,j,k)
               x19 = a19(i,j,k)

               density(i,j,k)  = x01 + x02 + x03 + x04 + x05 + x06 + x07 + x08 &
                    +x09 + x10 + x11 + x12 + x13 + x14 + x15 + x16 &
                    +x17 + x18 + x19 + cte1
!
               rhoinv = uno/density(i,j,k)
!
               vel_u(i,j,k) = (x01+x02+x03+x04+x05-x10-x11-x12-x13-x14)*rhoinv
               vel_v(i,j,k) = (x03+x07+x08+x09+x12-x01-x10-x16-x17-x18)*rhoinv
               vel_w(i,j,k) = (x04+x06+x07+x13+x18-x02-x09-x11-x15-x16)*rhoinv
!
               mean_u = mean_u + vel_u(i,j,k)
               mean_v = mean_v + vel_v(i,j,k)
               mean_w = mean_w + vel_w(i,j,k)
            enddo
         enddo
       enddo
!
       mean_u = mean_u/float(l+2)/float(m+2)/float(n+2)       
       mean_v = mean_v/float(l+2)/float(m+2)/float(n+2)       
       mean_w = mean_w/float(l+2)/float(m+2)/float(n+2)       
!
! compute energy dissipation / turbulent kinetic energy
! maybe not the most efficent way....      
       do k=1,n
         do j=1,m
            do i=1,l

! compute u_y
                u_y = 0.5*(vel_u(i,j+1,k)-vel_u(i,j-1,k))
! compute u_z
                u_z = 0.5*(vel_u(i,j,k+1)-vel_u(i,j,k-1))
! compute v_x
                v_x = 0.5*(vel_v(i+1,j,k)-vel_v(i-1,j,k))
! compute v_z
                v_z = 0.5*(vel_v(i,j,k+1)-vel_v(i,j,k-1))
! compute w_x
                w_x = 0.5*(vel_w(i+1,j,k)-vel_w(i-1,j,k))
! compute w_y
                w_y = 0.5*(vel_w(i,j+1,k)-vel_w(i,j-1,k))
!
! CHECK SIGN....
                diss = diss + (v_x + u_y)**2 +  & 
                              (u_z + w_x)**2 +  &
                              (w_y + v_z)**2 
!                      
                tke = tke + 0.5*((vel_u(i,j,k)-mean_u)**2 +  & 
                                 (vel_v(i,j,k)-mean_v)**2 +  &
                                 (vel_w(i,j,k)-mean_w)**2)     
             enddo
          enddo
       enddo
!
!# time, dissipation, turbulent kinetic energy
       dissip = 0.5*svisc*diss/float(l)/float(m)/float(n)       ! dissipation
       turbke = tke/float(l)/float(m)/float(n)                  ! turbulent kinetic energy
       kolmog = ((svisc**3)/dissip)**(0.25)                     ! kolmogorov microscale
!
       write(77,1004) itime,                                     &
                      dissip,                                    &
                      turbke,                                    &
                      kolmog
       flush(77)
!
#ifdef DEBUG_1
       write(6,*) "DEBUG1: Exiting from sub. dissipation"
#endif
!
! formats...
!
1001   format(" Timestep ",i8)
1002   format("       mean rho ",1(e14.6,1x))
1003   format("       mean vel ",4(e14.6,1x))
1004   format(i8,3(e14.6,1x))
!
       end subroutine dissipation
