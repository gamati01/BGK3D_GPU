!=======================================================================
!     ****** LBE/save_vel
!
!     COPYRIGHT
!       (c) 2000-2011 by CASPUR/G.Amati
!       (c) 2013-20?? by CINECA/G.Amati
!     NAME
!       save_raw
!     DESCRIPTION
!       save all veocities on save_vel.XXXXXXXX.vtk, where XXXXXXX is 
!       the timestep
!     INPUTS
!       itime --> timestep
!     OUTPUT
!       none
!     TODO
!	
!     NOTES
!       character*22 file_name
!       max time allowed  99'999'999
!       ifdef:
!              * NOSHIFT 
!              * DEBUG_1
!
!     *****
!=======================================================================
!
      subroutine save_vel(itime)
!
      use storage
      implicit none
!
      character(len=22) :: file_name
!
      integer, INTENT(in) :: itime
      integer             :: i,j,k
!
      real(mykind):: x01,x02,x03,x04,x05,x06,x07
      real(mykind):: x08,x09,x10,x11,x12,x13
      real(mykind):: x14,x15,x16,x17,x18,x19
      real(mykind):: cte1
      real(mykind):: rho,rhoinv
      real(mykind):: xj,yj,zj
!
      file_name = 'save_vel.xxxxxxxx.bin'
!
      write(file_name(10:17),4000) itime
      open(21,file=file_name,form="unformatted",status='unknown')
!
#ifdef NOSHIFT
       cte1 = zero
#else
       cte1 = uno
#endif
!
      if(myrank == 0) then
         write(16,*) 'task ', myrank, 'saving t=', itime, file_name
         write(6,*)  'task ', myrank, 'saving t=', itime, file_name
      endif
!
      write(21) itime
!
      do k = 1,n
         do j = 1,m
            do i = 1,l
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
!
               rho = x01 + x02 + x03 + x04 + x05 + x06 + x07 + x08 &
                    +x09 + x10 + x11 + x12 + x13 + x14 + x15 + x16 &
                    +x17 + x18 + x19 + cte1
!
               rhoinv = uno/rho
!
               xj = (x01+x02+x03+x04+x05-x10-x11-x12-x13-x14)*rhoinv
               yj = (x03+x07+x08+x09+x12-x01-x10-x16-x17-x18)*rhoinv
               zj = (x04+x06+x07+x13+x18-x02-x09-x11-x15-x16)*rhoinv

               write(21) xj, yj, zj
!
! check               
               if((i==l/2).AND.(j==m/2).AND.(k==n/2)) then
                  write(6,*) i, j, k, xj, yj, zj
               endif 
            enddo
         enddo
      enddo
      close(21) 
!
#ifdef DEBUG_1
      if(myrank == 0) then
         write(6,*) "DEBUG1: Exiting from sub. save_vel"
      endif
#endif
!
4000  format(i8.8)
!
      end subroutine save_vel
