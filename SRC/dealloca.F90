! =====================================================================
!     ****** LBE/dealloca
!
!     COPYRIGHT
!       (c) 2000-2011 by CASPUR/G.Amati
!       (c) 2013-20?? by CINECA/G.Amati
!     NAME
!       dealloca
!     DESCRIPTION
!       main program for LBM 3D
!     INPUTS
!       none
!     OUTPUT
!       none
!     TODO
!
!     NOTES
! =====================================================================

subroutine dealloca()
      use storage; 
      implicit none
!
      deallocate(a01)
      deallocate(a02)
      deallocate(a03)
      deallocate(a04)
      deallocate(a05)
      deallocate(a06)
      deallocate(a07)
      deallocate(a08)
      deallocate(a09)
      deallocate(a10)
      deallocate(a11)
      deallocate(a12)
      deallocate(a13)
      deallocate(a14)
      deallocate(a15)
      deallocate(a16)
      deallocate(a17)
      deallocate(a18)
      deallocate(a19)
!
      deallocate(b01)
      deallocate(b02)
      deallocate(b03)
      deallocate(b04)
      deallocate(b05)
      deallocate(b06)
      deallocate(b07)
      deallocate(b08)
      deallocate(b09)
      deallocate(b10)
      deallocate(b11)
      deallocate(b12)
      deallocate(b13)
      deallocate(b14)
      deallocate(b15)
      deallocate(b16)
      deallocate(b17)
      deallocate(b18)
      deallocate(b19)
!
#ifdef DEBUG_1
        if(myrank == 0) then
           write(6,*) "DEBUG1: Exiting from sub. dealloca"
        endif
#endif
!
end subroutine dealloca

