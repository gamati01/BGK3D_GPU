!=======================================================================
!     ****** LBE/save_raw
!
!     COPYRIGHT
!       (c) 2000-2011 by CASPUR/G.Amati
!       (c) 2013-20?? by CINECA/G.Amati
!     NAME
!       save_raw
!     DESCRIPTION
!       save all microscopic variables (all populations)
!       write on unit 20 (save.xxxxxxxx.bin, unformatted, 
!                          xxxxxxxx is timestep)
!     INPUTS
!       itime --> timestep
!     OUTPUT
!       none
!     TODO
!	
!     NOTES
!       character*18 file_name
!       integer itime,i,k
!       max time allowed  99'999'999
!
!     *****
!=======================================================================
!
      subroutine save_raw(itime)
!
      use storage
      implicit none
!
      character(len=18) :: file_name
!
      integer, INTENT(in) :: itime
      integer             :: i,j,k
!
      file_name = 'save.xxxxxxxx.bin'
!
      write(file_name(6:13),4000) itime
      open(20,file=file_name,form="unformatted",status='unknown')
!
      if(myrank == 0) then
         write(16,*) 'task ', myrank, 'saving t=', itime, file_name
         write(6,*)  'task ', myrank, 'saving t=', itime, file_name
      endif
!
      write(20) itime
!
      write(20) (((a01(i,j,k),i=0,l1),j=0,m1),k=0,n1)
      write(20) (((a02(i,j,k),i=0,l1),j=0,m1),k=0,n1)
      write(20) (((a03(i,j,k),i=0,l1),j=0,m1),k=0,n1)
      write(20) (((a04(i,j,k),i=0,l1),j=0,m1),k=0,n1)
      write(20) (((a05(i,j,k),i=0,l1),j=0,m1),k=0,n1)
      write(20) (((a06(i,j,k),i=0,l1),j=0,m1),k=0,n1)
      write(20) (((a07(i,j,k),i=0,l1),j=0,m1),k=0,n1)
      write(20) (((a08(i,j,k),i=0,l1),j=0,m1),k=0,n1)
      write(20) (((a09(i,j,k),i=0,l1),j=0,m1),k=0,n1)
      write(20) (((a10(i,j,k),i=0,l1),j=0,m1),k=0,n1)
      write(20) (((a11(i,j,k),i=0,l1),j=0,m1),k=0,n1)
      write(20) (((a12(i,j,k),i=0,l1),j=0,m1),k=0,n1)
      write(20) (((a13(i,j,k),i=0,l1),j=0,m1),k=0,n1)
      write(20) (((a14(i,j,k),i=0,l1),j=0,m1),k=0,n1)
      write(20) (((a15(i,j,k),i=0,l1),j=0,m1),k=0,n1)
      write(20) (((a16(i,j,k),i=0,l1),j=0,m1),k=0,n1)
      write(20) (((a17(i,j,k),i=0,l1),j=0,m1),k=0,n1)
      write(20) (((a18(i,j,k),i=0,l1),j=0,m1),k=0,n1)
      write(20) (((a19(i,j,k),i=0,l1),j=0,m1),k=0,n1)
      flush(20)
      close(20)
!
#ifdef DEBUG_1
      if(myrank == 0) then
         write(6,*) "DEBUG1: Exiting from sub. save_raw"
      endif
#endif
!
4000  format(i8.8)
!
      end subroutine save_raw
