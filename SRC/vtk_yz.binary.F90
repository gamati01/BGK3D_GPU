!=======================================================================
!     ****** LBE/vtk_yz_bin
!
!     COPYRIGHT
!       (c) 2009 by CASPUR/G.Amati
!       (c) 2000-2011 by CASPUR/G.Amati
!       (c) 2013-20?? by CINECA/G.Amati
!     NAME
!       vtk_yz_bin
!     DESCRIPTION
!       Graphic subroutine:
!       write 2D binary output for VTK with velocity + pressure field
!       write on unit 54 (vtk_yz.xxxxxxx.dat) where 
!                                xxxxxxx is the timestep
!       the file is closed at the end of the subroutine
!     INPUTS
!       itime   ---> timestep
!       io      ---> i coordinate of the 2D slice
!     OUTPUT
!       none
!     TODO
!       
!     NOTES
!       character used:  file_name (19)
!       integer variables used: itime,j,k
!       max time allowed  99'999'999
!       single precision only, to saving space..
!
!     *****
!=======================================================================
!
        subroutine vtk_yz_bin(itime,i0)
!
        use storage
        implicit none
!
        integer :: j,k,itime
        integer :: i0
!	
        character(len=24) :: file_name
!
        real(sp) :: u,w,v,den
        real(sp) :: cte1
        real(sp), parameter :: zerol=0.0
!
        file_name = 'tec_yz.xxxxxxxx.vtk'
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
        open(54,file=file_name,status='unknown')
!
        write(54,'(A26)')'# vtk DataFile Version 2.0'
        write(54,'(A5)') 'Campo'
        write(54,'(A6)') 'BINARY'
        write(54,'(A25)')'DATASET STRUCTURED_POINTS'
        write(54,'(A11,I10,A1,I10,A1,I10)') 'DIMENSIONS ',1,' ',m,' ',n
        write(54,'(A7,I10,A1,I10,A1,I10)')  'ORIGIN ',i0,' ',1,' ',1  
        write(54,'(A8,I10,A1,I10,A1,I10)') 'SPACING ',1,' ',1,' ',1     
        write(54,'(A10,I10)')'POINT_DATA ',m*n*1
!
! vector
        write(54,'(A23)')'VECTORS velocity float'
        close(54)
!
! then write output (binary)
        open(54,file=file_name,status='old', position='append', & 
                form='unformatted',access='STREAM',CONVERT="BIG_ENDIAN")
!
        do k = 1,n
           do j = 1,m
              u = a01(i0,j,k)+a02(i0,j,k)+a03(i0,j,k)+a04(i0,j,k)+a05(i0,j,k) &
                 -a10(i0,j,k)-a11(i0,j,k)-a12(i0,j,k)-a13(i0,j,k)-a14(i0,j,k)
!
              w = a03(i0,j,k)+a07(i0,j,k)+a08(i0,j,k)+a09(i0,j,k)+a12(i0,j,k) &
                 -a01(i0,j,k)-a10(i0,j,k)-a16(i0,j,k)-a17(i0,j,k)-a18(i0,j,k)
!
              v = a04(i0,j,k)+a06(i0,j,k)+a07(i0,j,k)+a13(i0,j,k)+a18(i0,j,k) &
                 -a02(i0,j,k)-a09(i0,j,k)-a11(i0,j,k)-a15(i0,j,k)-a16(i0,j,k)
!
              den = a01(i0,j,k)+a02(i0,j,k)+a03(i0,j,k)+a04(i0,j,k) &
                   +a05(i0,j,k)+a06(i0,j,k)+a07(i0,j,k)+a08(i0,j,k) &
                   +a09(i0,j,k)+a10(i0,j,k)+a11(i0,j,k)+a12(i0,j,k) &
                   +a13(i0,j,k)+a14(i0,j,k)+a15(i0,j,k)+a16(i0,j,k) &
                   +a17(i0,j,k)+a18(i0,j,k)+a19(i0,j,k)+cte1
!
              write(54) u/den, w/den, v/den
           end do
        end do
        close(54)
!
        write(16,*) "I/O: plane yz (vtk,binary) done", i0
        write(6,*)  "I/O: plane yz (vtk,binary) done", i0
!
#ifdef DEBUG_1
        if(myrank == 0) then
           write(6,*) "DEBUG1: Exiting from sub. vtk_yz_bin"
        endif
#endif
!
4000    format(i8.8)
!
       end subroutine vtk_yz_bin
