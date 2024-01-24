!=====================================================================
!     ****** LBE/diagno
!
!     COPYRIGHT
!       (c) 2000-2011 by CASPUR/G.Amati
!       (c) 2013-20?? by CINECA/G.Amati
!     NAME
!       diagno
!     DESCRIPTION
!       diagnostic subroutine:
!       check of conserved quantities and mean value
!       write on unit 63 (diagno.dat)
!       write on unit 16 (bgk.log)
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
       subroutine diagno(itime)
!
       use storage
       implicit none
!
       integer, INTENT(in) :: itime
       integer             :: i,j,k
!
       real(mykind):: rtot,xtot,ytot,ztot,stot
       real(mykind):: tke
       real(mykind):: rho,rhoinv
       real(mykind):: x01,x02,x03,x04,x05,x06,x07
       real(mykind):: x08,x09,x10,x11,x12,x13
       real(mykind):: x14,x15,x16,x17,x18,x19
       real(mykind):: cte1
       real(mykind):: xj,yj,zj
!
#ifdef NOSHIFT
       cte1 = zero
#else
       cte1 = uno
#endif
!
       tke  = zero
       rtot = zero
       xtot = zero
       ytot = zero
       ztot = zero
       stot = zero
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
       do k=1,n
         do j=1,m
            do i=1,l
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

               rho = x01 + x02 + x03 + x04 + x05 + x06 + x07 + x08 &
                    +x09 + x10 + x11 + x12 + x13 + x14 + x15 + x16 &
                    +x17 + x18 + x19 + cte1
!
               rhoinv = uno/rho
!
               xj = (x01+x02+x03+x04+x05-x10-x11-x12-x13-x14)*rhoinv
               yj = (x03+x07+x08+x09+x12-x01-x10-x16-x17-x18)*rhoinv
               zj = (x04+x06+x07+x13+x18-x02-x09-x11-x15-x16)*rhoinv

!
               rtot = rtot+rho
               xtot = xtot+xj
               ytot = ytot+yj
               ztot = ytot+yj
               stot = stot+(xj*xj+yj*yj)
!
            enddo
         enddo
       enddo
!
       rtot = (rtot/float(l))/float(m)/float(n)
       xtot = (xtot/float(l))/float(m)/float(n)
       ytot = (ytot/float(l))/float(m)/float(n)
       ztot = (ytot/float(l))/float(m)/float(n)
       stot = (stot/float(l))/float(m)/float(n)
!
       write(16,1001) itime
       write(16,1002) rtot
       write(16,1003) xtot,ytot,ztot,stot
       flush(16)
!       
       write(63,1004) itime, xtot, ytot, ztot, rtot, stot
       flush(63)
!
! conmpute turbulent kinetic energy
       do k=1,n
         do j=1,m
            do i=1,l
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

               rho = x01 + x02 + x03 + x04 + x05 + x06 + x07 + x08 &
                    +x09 + x10 + x11 + x12 + x13 + x14 + x15 + x16 &
                    +x17 + x18 + x19 + cte1
!
               rhoinv = uno/rho
!
               xj = (x01+x02+x03+x04+x05-x10-x11-x12-x13-x14)*rhoinv
               yj = (x03+x07+x08+x09+x12-x01-x10-x16-x17-x18)*rhoinv
               zj = (x04+x06+x07+x13+x18-x02-x09-x11-x15-x16)*rhoinv

               tke = tke + 0.5*((xtot - xj)*(xtot - xj) + & 
                                (ytot - yj)*(ytot - yj) + & 
                                (ztot - zj)*(ztot - zj))  


           enddo
         enddo
       enddo

       write(666,*) itime, tke/(float(l)*float(m)*float(n))
       flush(666)

#ifdef DEBUG_1
       if (myrank == 0) then
          write(6,*) "DEBUG1: Exiting from sub. diagno"
       endif
#endif
!
! formats...
!
1001   format(" Timestep ",i8)
1002   format("       mean rho ",1(e14.6,1x))
1003   format("       mean vel ",4(e14.6,1x))
1004   format(i8,5(e14.6,1x))
!
       end subroutine diagno
