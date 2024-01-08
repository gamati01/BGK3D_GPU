!=======================================================================
!     ****** LBE/prof_j
!
!     COPYRIGHT
!       (c) 2000-2011 by CASPUR/G.Amati
!       (c) 2013-20?? by CINECA/G.Amati
!     NAME
!       prof_j
!     DESCRIPTION
!       Diagnostic subroutine:
!       istantaneous profile along y-direction (i fixed)
!       write on unit 64 (prof_j.dat)
!     INPUTS
!       itime   -->  timestep
!       icoord  -->  x coordinate
!     OUTPUT
!       none
!     TODO
!
!     NOTES
!       integer variables used: 
!       real variables used: 
!                            
!
!     *****
!=======================================================================
!
        subroutine prof_j(itime,icoord,kcoord)
!
        use storage
        implicit none
!
        integer             :: j
        integer, INTENT(in) :: itime
        integer, INTENT(in) :: icoord, kcoord
!
        real(mykind) :: u(1:m),v(1:m),w(1:m)   ! istantaneous velocity fields
        real(mykind) :: den(1:m)               ! istantaneous density field
        real(mykind) :: cte1
!
#ifdef NOSHIFT
       cte1 = zero
#else
       cte1 = uno
#endif
!
! density
        do concurrent (j=1:m)
           den(j) = +a01(icoord,j,kcoord)+a02(icoord,j,kcoord)+a03(icoord,j,kcoord) &
                    +a04(icoord,j,kcoord)+a05(icoord,j,kcoord)+a06(icoord,j,kcoord) &
                    +a07(icoord,j,kcoord)+a08(icoord,j,kcoord)+a09(icoord,j,kcoord) &
                    +a10(icoord,j,kcoord)+a11(icoord,j,kcoord)+a12(icoord,j,kcoord) &
                    +a13(icoord,j,kcoord)+a14(icoord,j,kcoord)+a15(icoord,j,kcoord) &
                    +a16(icoord,j,kcoord)+a17(icoord,j,kcoord)+a18(icoord,j,kcoord) &
                    +a19(icoord,j,kcoord)+cte1
        enddo
!        
        do j = 1,m         ! streamwise velocity
           u(j) = +( a01(icoord,j,kcoord)+a02(icoord,j,kcoord)+a03(icoord,j,kcoord) &
                    +a04(icoord,j,kcoord)+a05(icoord,j,kcoord) &
                    -a10(icoord,j,kcoord)-a11(icoord,j,kcoord)-a12(icoord,j,kcoord) &
                    -a13(icoord,j,kcoord)-a14(icoord,j,kcoord) ) / den(j)
        end do
!        
        do j = 1,m         ! spanwise velocity
           w(j) = +( a03(icoord,j,kcoord)+a07(icoord,j,kcoord)+a08(icoord,j,kcoord) &
                    +a09(icoord,j,kcoord)+a12(icoord,j,kcoord) &
                    -a01(icoord,j,kcoord)-a10(icoord,j,kcoord)-a16(icoord,j,kcoord) &
                    -a17(icoord,j,kcoord)-a18(icoord,j,kcoord) ) / den(j)
        end do
!        
        do j = 1,m         ! normal_to_wall velocity
           v(j) = +( a04(icoord,j,kcoord)+a06(icoord,j,kcoord)+a07(icoord,j,kcoord) &
                    +a13(icoord,j,kcoord)+a18(icoord,j,kcoord) &
                    -a02(icoord,j,kcoord)-a09(icoord,j,kcoord)-a11(icoord,j,kcoord) &
                    -a15(icoord,j,kcoord)-a16(icoord,j,kcoord) ) / den(j)
        end do
! 
        write(64,1005) itime
!
        do j=1,m
           write(64,1002) (j-0.5), u(j),w(j),v(j),den(j)
        end do
        write(64,'(a1)') 
        write(64,'(a1)') 
!
#ifdef DEBUG_1
        if(myrank == 0) then
           write(6,*) "DEBUG1: Exiting from sub. prof_j", icoord, kcoord
        endif
#endif
!
!	format
1002    format(5(e14.6,1x))
1005    format("# t=",i7)
!       
        end subroutine prof_j
