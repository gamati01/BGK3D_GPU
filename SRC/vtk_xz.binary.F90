!=======================================================================
!     ****** LBE/vtk_xz_bin
!
!     COPYRIGHT
!       (c) 2009 by CASPUR/G.Amati
!       (c) 2000-2011 by CASPUR/G.Amati
!       (c) 2013-20?? by CINECA/G.Amati
!     NAME
!       vtk_xz_bin
!     DESCRIPTION
!       Graphic subroutine:
!       write 2D binary output for VTK with velocity + pressure field
!       write on unit 53 (vtk_xy.xxxxxxx.dat) where 
!                                  xxxxxxx is the timestep
!       the file is closed at the end of the subroutine
!     INPUTS
!       itime   ---> timestep
!       j0      ---> k coordinate of the 2D slice
!     OUTPUT
!       none
!     TODO
!       
!     NOTES
!       character used:  file_name (19)
!       integer variables used: itime,i,k, x_scale, z_scale
!       max time allowed  99'999'999
!       single precision only, to saving space..
!
!     *****
!=======================================================================
!
        subroutine vtk_xz_bin(itime,j0)
!
        use storage
        implicit none
!
        integer :: i,k,itime
        integer :: j0
!	
        character(len=24) :: file_name
!
        real(sp) :: u,w,v,den
        real(sp) :: cte1
        real(sp), parameter :: zerol=0.0
!
        file_name = 'tec_xz.xxxxxxxx.vtk'
!
        myrank = 0
!
#ifdef NOSHIFT
       cte1 = zero
#else
       cte1 = uno
#endif
!
        write(file_name(8:15),4000) itime
!
! first write legal header (ASCII)
        open(53,file=file_name,status='unknown')
!
        write(53,'(A26)')'# vtk DataFile Version 2.0'
        write(53,'(A5)') 'Campo'
        write(53,'(A6)') 'BINARY'
        write(53,'(A25)')'DATASET STRUCTURED_POINTS'
        write(53,'(A11,I10,A1,I10,A1,I10)') 'DIMENSIONS ',l,' ',1,' ',n
        write(53,'(A7,I10,A1,I10,A1,I10)')  'ORIGIN ',1,' ',j0,' ',1  
        write(53,'(A8,I10,A1,I10,A1,I10)') 'SPACING ',1,' ',1,' ',1     
        write(53,'(A10,I10)')'POINT_DATA ',l*n*1
!
! vector
        write(53,'(A23)')'VECTORS velocity float'
        close(53)
!
! then write output (binary)
        open(53,file=file_name,status='old', position='append', & 
                form='unformatted',access='STREAM',CONVERT="BIG_ENDIAN")
!
        do k = 1,n
           do i = 1,l
              u = a01(i,j0,k)+a02(i,j0,k)+a03(i,j0,k)+a04(i,j0,k)+a05(i,j0,k) &
                 -a10(i,j0,k)-a11(i,j0,k)-a12(i,j0,k)-a13(i,j0,k)-a14(i,j0,k)
!
              w = a03(i,j0,k)+a07(i,j0,k)+a08(i,j0,k)+a09(i,j0,k)+a12(i,j0,k) &
                 -a01(i,j0,k)-a10(i,j0,k)-a16(i,j0,k)-a17(i,j0,k)-a18(i,j0,k)
!
              v = a04(i,j0,k)+a06(i,j0,k)+a07(i,j0,k)+a13(i,j0,k)+a18(i,j0,k) &
                 -a02(i,j0,k)-a09(i,j0,k)-a11(i,j0,k)-a15(i,j0,k)-a16(i,j0,k)
!
              den = a01(i,j0,k)+a02(i,j0,k)+a03(i,j0,k)+a04(i,j0,k) &
                   +a05(i,j0,k)+a06(i,j0,k)+a07(i,j0,k)+a08(i,j0,k) &
                   +a09(i,j0,k)+a10(i,j0,k)+a11(i,j0,k)+a12(i,j0,k) &
                   +a13(i,j0,k)+a14(i,j0,k)+a15(i,j0,k)+a16(i,j0,k) &
                   +a17(i,j0,k)+a18(i,j0,k)+a19(i,j0,k)+cte1
!
              write(53) u/den, w/den, v/den
           end do
        end do
        close(53)
!
        write(16,*) "I/O: plane xz (vtk,binary) done", j0
        write(6,*)  "I/O: plane xz (vtk,binary) done", j0
!
#ifdef DEBUG_1
        if(myrank == 0) then
           write(6,*) "DEBUG1: Exiting from sub. vtk_xz_bin"
        endif
#endif
!
4000    format(i8.8)
!
       end subroutine vtk_xz_bin
