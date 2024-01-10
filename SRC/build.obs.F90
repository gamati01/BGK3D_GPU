!=====================================================================
!     ****** LBE/build_obs
!
!     COPYRIGHT
!       (c) 2000-2011 by CASPUR/G.Amati
!       (c) 2013-20?? by CINECA/G.Amati
!     NAME
!       bcond
!     DESCRIPTION
!       build obstacles: setting 1 in obs(i,j) field 
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
      subroutine build_obs
!
      use timing
      use storage
!
      implicit none
!
      integer:: i, j, k
      integer:: icoord, jcoord, kcoord
      integer:: itime
      real(mykind):: d2, R2a, R2b, R
!
      imin = l
      jmin = m
      kmin = n
      imax = 0
      jmax = 0
      kmax = 0
!      
      nobs = 0
!
      myrank = 0
      itime  = 0
!
! creating  obstacle 
!
! creating obstacles
      write( 6,*) "INFO: creating obstacle (sphere)" 
      write(16,*) "INFO: creating obstacle (sphere)" 
!
      icoord = 2*l/5
      jcoord = m/2
      kcoord = n/2
!      
      if (m.gt.256) then
          radius = 32
      else
          radius = ly/8  
      endif
!      
      write(16,*) "WARNING: radius hardwired ", radius
      write( 6,*) "WARNING: radius hardwired ", radius
      write( 6,*) "INFO: Cyl radius    -->", radius
      write( 6,*) "INFO: Cyl icoord    -->", icoord
      write( 6,*) "INFO: Cyl jcoord    -->", jcoord
!
      R = radius*uno
      R2a = (R-2)*(R-2)   ! lower radius
      R2b = (R+2)*(R+2)   ! upper radius
!
      do k = 1, n
         do j = 1, m
            do i = 1, l
!
               d2 = (icoord-i)*(icoord-i)  &
                   +(jcoord-j)*(jcoord-j)  & 
                   +(kcoord-k)*(kcoord-k)  
!        
               if((d2.gt.R2a).and.(d2.lt.R2b)) then
!
                   obs(i,j,k) = 1
                   nobs = nobs + 1
!
                   imin = min(imin,i)
                   jmin = min(jmin,j)
                   kmin = min(kmin,k)
!
                   imax = max(imax,i)
                   jmax = max(jmax,j)
                   kmax = max(kmax,k)
!
                endif
            end do
         end do
      end do
!
      write(6,*) "INFO: num. obs         -->",nobs
      write(6,*) "INFO: ratio obs/size   -->",nobs/float(l*m*n)
      write(6,*) "INFO: obs (x)             ",imax, imin
      write(6,*) "INFO: obs (y)             ",jmax, jmin
      write(6,*) "INFO: obs (z)             ",kmax, kmin
!
!#ifdef DEBUG_1
      if(myrank == 0) then
         write(6,*) "DEBUG1: Exiting from sub. build_obs"
      endif
!#endif
!
      end subroutine build_obs
