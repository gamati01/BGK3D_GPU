!=======================================================================
!     ****** LBE/prof_k
!
!     COPYRIGHT
!       (c) 2000-2011 by CASPUR/G.Amati
!       (c) 2013-20?? by CINECA/G.Amati
!     NAME
!       prof_k
!     DESCRIPTION
!       Diagnostic subroutine:
!       istantaneous profile along y-direction (i fixed)
!       write on unit 65 (prof_k.dat)
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
        subroutine prof_k(itime,icoord,jcoord)
!
        use storage
        implicit none
!
        integer             :: k
        integer, INTENT(in) :: itime
        integer, INTENT(in) :: icoord, jcoord
!
        real(mykind) :: u(1:n),v(1:n),w(1:n)   ! istantaneous velocity fields
        real(mykind) :: den(1:n)               ! istantaneous density field
        real(mykind) :: cte1
!
#ifdef NOSHIFT
       cte1 = zero
#else
       cte1 = uno
#endif
!
        do k = 1,n         ! density  
           den(k) = +a01(icoord,jcoord,k)+a02(icoord,jcoord,k)+a03(icoord,jcoord,k) &
                    +a04(icoord,jcoord,k)+a05(icoord,jcoord,k)+a06(icoord,jcoord,k) &
                    +a07(icoord,jcoord,k)+a08(icoord,jcoord,k)+a09(icoord,jcoord,k) &
                    +a10(icoord,jcoord,k)+a11(icoord,jcoord,k)+a12(icoord,jcoord,k) &
                    +a13(icoord,jcoord,k)+a14(icoord,jcoord,k)+a15(icoord,jcoord,k) &
                    +a16(icoord,jcoord,k)+a17(icoord,jcoord,k)+a18(icoord,jcoord,k) &
                    +a19(icoord,jcoord,k)+cte1

        enddo
!        
        do k = 1,n         ! streamwise velocity
           u(k) = +( a01(icoord,jcoord,k)+a02(icoord,jcoord,k)+a03(icoord,jcoord,k) &
                    +a04(icoord,jcoord,k)+a05(icoord,jcoord,k) &
                    -a10(icoord,jcoord,k)-a11(icoord,jcoord,k)-a12(icoord,jcoord,k) &
                    -a13(icoord,jcoord,k)-a14(icoord,jcoord,k) ) / den(k)
        end do
!        
        do k = 1,n         ! spanwise velocity
           w(k) = +( a03(icoord,jcoord,k)+a07(icoord,jcoord,k)+a08(icoord,jcoord,k) &
                    +a09(icoord,jcoord,k)+a12(icoord,jcoord,k) &
                    -a01(icoord,jcoord,k)-a10(icoord,jcoord,k)-a16(icoord,jcoord,k) &
                    -a17(icoord,jcoord,k)-a18(icoord,jcoord,k) ) / den(k)

        end do
!        
        do k = 1,n         ! normal_to_wall velocity
           v(k) = +( a04(icoord,jcoord,k)+a06(icoord,jcoord,k)+a07(icoord,jcoord,k) &
                    +a13(icoord,jcoord,k)+a18(icoord,jcoord,k) &
                    -a02(icoord,jcoord,k)-a09(icoord,jcoord,k)-a11(icoord,jcoord,k) &
                    -a15(icoord,jcoord,k)-a16(icoord,jcoord,k) ) / den(k)

        end do
! 
        write(65,1005) itime
!
        do k=1,n
           write(65,1002) (k-0.5), u(k),w(k),v(k),den(k)
        end do
        write(65,'(a1)') 
        write(65,'(a1)') 
!
#ifdef DEBUG_1
        if(myrank == 0) then
           write(6,*) "DEBUG1: Exiting from sub. prof_k", icoord, jcoord
        endif
#endif
!
!	format
1002    format(5(e14.6,1x))
1005    format("# t=",i7)
!       
        end subroutine prof_k
