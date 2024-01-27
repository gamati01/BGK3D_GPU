!=====================================================================
!     ****** LBE/build_bcond
!
!     COPYRIGHT
!       (c) 2000-2011 by CASPUR/G.Amati
!       (c) 2013-20?? by CINECA/G.Amati
!     NAME
!       bcond
!     DESCRIPTION
!       to remove?
!     INPUTS
!       none
!     OUTPUT
!       none
!     TODO
!       
!     NOTES
!
!     *****
!=====================================================================
!
      subroutine build_bcond
!
      use timing
      use storage
!
      implicit none
!
#ifdef CHANNEL
      u00 = 0.0            ! boundary condition
      u0  = u0             ! volume force
#else
      u00 = u0             ! boundary condition
      u0  = 0.0            ! volume force
#endif
!
      write(16,*) "INFO: reference velocities --->", u0, u00
!
#ifdef DEBUG_1
      write(6,*) "DEBUG1: Exiting from sub. build_bcond"
#endif
!
      end subroutine build_bcond
