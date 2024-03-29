!=====================================================================
!     ****** LBE/probe
!
!     COPYRIGHT
!       (c) 2000-2011 by CASPUR/G.Amati
!       (c) 2013-20?? by CINECA/G.Amati
!     NAME
!       probe
!     DESCRIPTION
!       Diagnostic subroutine:
!       probe popolations and/or velocities in one single point	
!       write on unit 68 (probe.dat)
!     INPUTS
!       itime --> timestep
!       i0     --> x coordinate
!       j0     --> y coordinate
!       k0     --> z coordinate
!     OUTPUT
!       none
!     TODO
!
!     NOTES
!       integer variables used: i,j,k,itime
!       ifdef
!            * NOSHIFT
!            * DEBUG_1      
!
!     *****
!=====================================================================
!
        subroutine probe(itime,i0,j0,k0)
!
        use storage
        implicit none
!
        integer, INTENT(in) :: i0,j0,k0,itime
!
        real(mykind) :: rho, xj, yj, zj
        real(mykind) :: x01,x02,x03,x04,x05,x06,x07,x08,x09,x10
        real(mykind) :: x11,x12,x13,x14,x15,x16,x17,x18,x19
        real(mykind) :: cte1
!
#ifdef NOSHIFT
       cte1 = zero
#else
       cte1 = uno
#endif
!
        x01 = a01(i0,j0,k0)
        x02 = a02(i0,j0,k0)
        x03 = a03(i0,j0,k0)
        x04 = a04(i0,j0,k0)
        x05 = a05(i0,j0,k0)
        x06 = a06(i0,j0,k0)
        x07 = a07(i0,j0,k0)
        x08 = a08(i0,j0,k0)
        x09 = a09(i0,j0,k0)
        x10 = a10(i0,j0,k0)
        x11 = a11(i0,j0,k0)
        x12 = a12(i0,j0,k0)
        x13 = a13(i0,j0,k0)
        x14 = a14(i0,j0,k0)
        x15 = a15(i0,j0,k0)
        x16 = a16(i0,j0,k0)
        x17 = a17(i0,j0,k0)
        x18 = a18(i0,j0,k0)
        x19 = a19(i0,j0,k0)
!
        rho = x01 + x02 + x03 + x04 + x05 + x06 + x07 + x08 &
             +x09 + x10 + x11 + x12 + x13 + x14 + x15 + x16 &
             +x17 + x18 + x19 + cte1
!
        xj = (x01+x02+x03+x04+x05-x10-x11-x12-x13-x14)
        yj = (x03+x07+x08+x09+x12-x01-x10-x16-x17-x18)
        zj = (x04+x06+x07+x13+x18-x02-x09-x11-x15-x16)
!
      write(68,1002) itime, xj/rho, yj/rho, zj/rho, rho
!
#ifdef DEBUG_1
      if(myrank == 0) then
         write(6,*) "DEBUG1: Exiting from sub. probe", i0,j0,k0
      endif
#endif
!
! format
1002    format(i8,4(e14.6,1x))
!
       end subroutine probe
