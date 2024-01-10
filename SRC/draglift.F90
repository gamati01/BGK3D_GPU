!=====================================================================
!     ****** LBE/draglift
!
!     COPYRIGHT
!       (c) 2000-2011 by CASPUR/G.Amati
!       (c) 2013-20?? by CINECA/G.Amati
!     NAME
!       draglift
!     DESCRIPTION
!       computing total drag/lift  in "easy way" around a box with the obstacle
!     INPUTS
!       itime
!     OUTPUT
!     TODO
!
!     NOTES
!       drag/lift on file drag.lift (unit=66)
!       verify that the obstacle is inside the box       
!
!     *****
!=====================================================================
!
        subroutine draglift(itime)
!
        use storage
        use timing
        implicit none
!
        integer             :: i,j,k
        integer, INTENT(in) :: itime
        integer             :: istart,istop
        integer             :: jstart,jstop
        integer             :: kstart,kstop
        integer             :: icoord, jcoord, kcoord
        integer             :: delta
!
        real(mykind) :: norm
        real(mykind) :: forceX, forceY, forceZ
        real(mykind) :: force01, force02, force03, force04
        real(mykind) :: force05, force06, force07, force08
        real(mykind) :: force09, force10, force11, force12
        real(mykind) :: force13, force14, force15, force16
        real(mykind) :: force17, force18
!
! sphere center
        icoord = 2*l/5
        jcoord = m/2
        kcoord = n/2
! 
! border 
        delta = 3
!
! bounding box
        istart = icoord - radius - delta
        istop  = icoord + radius + delta
        jstart = jcoord - radius - delta
        jstop  = jcoord + radius + delta
        kstart = kcoord - radius - delta
        kstop  = kcoord + radius + delta
!          
! computing drag (force along x) 
        forceX = zero
!
#ifdef OFFLOAD
!$omp target update from(a01,a03,a05,a08,a10,a12,a14,a17,a19)
#endif
!
        do k = kstart, kstop
           do j = jstart, jstop
              do i = istart, istop
              if(obs(i,j,k)==0) then
                 force01 = 2.0*cx(01)*a01(i,j,k)*obs(i+icx(01),j+icy(01),k+icz(01))
                 force02 = 2.0*cx(02)*a02(i,j,k)*obs(i+icx(02),j+icy(02),k+icz(02))
                 force03 = 2.0*cx(03)*a03(i,j,k)*obs(i+icx(03),j+icy(03),k+icz(03))
                 force04 = 2.0*cx(04)*a04(i,j,k)*obs(i+icx(04),j+icy(04),k+icz(04))
                 force05 = 2.0*cx(05)*a05(i,j,k)*obs(i+icx(05),j+icy(05),k+icz(05))
                 force06 = 2.0*cx(06)*a06(i,j,k)*obs(i+icx(06),j+icy(06),k+icz(06))
                 force07 = 2.0*cx(07)*a07(i,j,k)*obs(i+icx(07),j+icy(07),k+icz(07))
                 force08 = 2.0*cx(08)*a08(i,j,k)*obs(i+icx(08),j+icy(08),k+icz(08))
                 force09 = 2.0*cx(09)*a09(i,j,k)*obs(i+icx(09),j+icy(09),k+icz(09))
                 force10 = 2.0*cx(10)*a10(i,j,k)*obs(i+icx(10),j+icy(10),k+icz(10))
                 force11 = 2.0*cx(11)*a11(i,j,k)*obs(i+icx(11),j+icy(11),k+icz(11))
                 force12 = 2.0*cx(12)*a12(i,j,k)*obs(i+icx(12),j+icy(12),k+icz(12))
                 force13 = 2.0*cx(13)*a13(i,j,k)*obs(i+icx(13),j+icy(13),k+icz(13))
                 force14 = 2.0*cx(14)*a14(i,j,k)*obs(i+icx(14),j+icy(14),k+icz(14))
                 force15 = 2.0*cx(15)*a15(i,j,k)*obs(i+icx(15),j+icy(15),k+icz(15))
                 force16 = 2.0*cx(16)*a16(i,j,k)*obs(i+icx(16),j+icy(16),k+icz(16))
                 force17 = 2.0*cx(17)*a17(i,j,k)*obs(i+icx(17),j+icy(17),k+icz(17))
                 force18 = 2.0*cx(18)*a18(i,j,k)*obs(i+icx(18),j+icy(18),k+icz(18))
!                                
                 forceX = forceX+(force01+force02+force03+force04+force05+ & 
                                  force06+force07+force08+force09+force10+ & 
                                  force11+force12+force13+force14+force15+ & 
                                  force16+force17+force18)
              endif 
              enddo
           enddo
        enddo
!
! computing lift (force along y) 
        forceY = zero
        do k = kstart, kstop
           do j = jstart, jstop
              do i = istart, istop
              if(obs(i,j,k)==0) then
                 force01 = 2.0*cy(01)*a01(i,j,k)*obs(i+icx(01),j+icy(01),k+icz(01))
                 force02 = 2.0*cy(02)*a02(i,j,k)*obs(i+icx(02),j+icy(02),k+icz(02))
                 force03 = 2.0*cy(03)*a03(i,j,k)*obs(i+icx(03),j+icy(03),k+icz(03))
                 force04 = 2.0*cy(04)*a04(i,j,k)*obs(i+icx(04),j+icy(04),k+icz(04))
                 force05 = 2.0*cy(05)*a05(i,j,k)*obs(i+icx(05),j+icy(05),k+icz(05))
                 force06 = 2.0*cy(06)*a06(i,j,k)*obs(i+icx(06),j+icy(06),k+icz(06))
                 force07 = 2.0*cy(07)*a07(i,j,k)*obs(i+icx(07),j+icy(07),k+icz(07))
                 force08 = 2.0*cy(08)*a08(i,j,k)*obs(i+icx(08),j+icy(08),k+icz(08))
                 force09 = 2.0*cy(09)*a09(i,j,k)*obs(i+icx(09),j+icy(09),k+icz(09))
                 force10 = 2.0*cy(10)*a10(i,j,k)*obs(i+icx(10),j+icy(10),k+icz(10))
                 force11 = 2.0*cy(11)*a11(i,j,k)*obs(i+icx(11),j+icy(11),k+icz(11))
                 force12 = 2.0*cy(12)*a12(i,j,k)*obs(i+icx(12),j+icy(12),k+icz(12))
                 force13 = 2.0*cy(13)*a13(i,j,k)*obs(i+icx(13),j+icy(13),k+icz(13))
                 force14 = 2.0*cy(14)*a14(i,j,k)*obs(i+icx(14),j+icy(14),k+icz(14))
                 force15 = 2.0*cy(15)*a15(i,j,k)*obs(i+icx(15),j+icy(15),k+icz(15))
                 force16 = 2.0*cy(16)*a16(i,j,k)*obs(i+icx(16),j+icy(16),k+icz(16))
                 force17 = 2.0*cy(17)*a17(i,j,k)*obs(i+icx(17),j+icy(17),k+icz(17))
                 force18 = 2.0*cy(18)*a18(i,j,k)*obs(i+icx(18),j+icy(18),k+icz(18))
!
                 forceY = forceY+(force01+force02+force03+force04+force05+ &
                                  force06+force07+force08+force09+force10+ &
                                  force11+force12+force13+force14+force15+ &
                                  force16+force17+force18)
              endif
              enddo
           enddo
        enddo
!
!
! computing lift (force along y) 
        forceY = zero
        do k = kstart, kstop
           do j = jstart, jstop
              do i = istart, istop
              if(obs(i,j,k)==0) then
                 force01 = 2.0*cz(01)*a01(i,j,k)*obs(i+icx(01),j+icy(01),k+icz(01))
                 force02 = 2.0*cz(02)*a02(i,j,k)*obs(i+icx(02),j+icy(02),k+icz(02))
                 force03 = 2.0*cz(03)*a03(i,j,k)*obs(i+icx(03),j+icy(03),k+icz(03))
                 force04 = 2.0*cz(04)*a04(i,j,k)*obs(i+icx(04),j+icy(04),k+icz(04))
                 force05 = 2.0*cz(05)*a05(i,j,k)*obs(i+icx(05),j+icy(05),k+icz(05))
                 force06 = 2.0*cz(06)*a06(i,j,k)*obs(i+icx(06),j+icy(06),k+icz(06))
                 force07 = 2.0*cz(07)*a07(i,j,k)*obs(i+icx(07),j+icy(07),k+icz(07))
                 force08 = 2.0*cz(08)*a08(i,j,k)*obs(i+icx(08),j+icy(08),k+icz(08))
                 force09 = 2.0*cz(09)*a09(i,j,k)*obs(i+icx(09),j+icy(09),k+icz(09))
                 force10 = 2.0*cz(10)*a10(i,j,k)*obs(i+icx(10),j+icy(10),k+icz(10))
                 force11 = 2.0*cz(11)*a11(i,j,k)*obs(i+icx(11),j+icy(11),k+icz(11))
                 force12 = 2.0*cz(12)*a12(i,j,k)*obs(i+icx(12),j+icy(12),k+icz(12))
                 force13 = 2.0*cz(13)*a13(i,j,k)*obs(i+icx(13),j+icy(13),k+icz(13))
                 force14 = 2.0*cz(14)*a14(i,j,k)*obs(i+icx(14),j+icy(14),k+icz(14))
                 force15 = 2.0*cz(15)*a15(i,j,k)*obs(i+icx(15),j+icy(15),k+icz(15))
                 force16 = 2.0*cz(16)*a16(i,j,k)*obs(i+icx(16),j+icy(16),k+icz(16))
                 force17 = 2.0*cz(17)*a17(i,j,k)*obs(i+icx(17),j+icy(17),k+icz(17))
                 force18 = 2.0*cz(18)*a18(i,j,k)*obs(i+icx(18),j+icy(18),k+icz(18))
!
                 forceZ = forceZ+(force01+force02+force03+force04+force05+ &
                                  force06+force07+force08+force09+force10+ &
                                  force11+force12+force13+force14+force15+ &
                                  force16+force17+force18)
              endif
              enddo
           enddo
        enddo
!
        norm = uno/(u_inflow*u_inflow*radius)
        write(66,*) itime, forceX*norm, forceY*norm, forceZ*norm
!
        call flush(66)            ! flush for drag/lift
!
#ifdef DEBUG_2
        if(myrank == 0) then
              write(6,*) "DEBUG2: Exiting from sub. drag2", itime
              write(6,*) "DEBUG2", forceX, forceY, norm
        endif
#endif
!
        end subroutine draglift
