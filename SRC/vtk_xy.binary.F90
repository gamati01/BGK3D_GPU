!=======================================================================
!     ****** LBE/vtk_xy_bin
!
!     COPYRIGHT
!       (c) 2009 by CASPUR/G.Amati
!       (c) 2000-2011 by CASPUR/G.Amati
!       (c) 2013-20?? by CINECA/G.Amati
!     NAME
!       vtk_xy_bin
!     DESCRIPTION
!       Graphic subroutine:
!       write 2D binary output for VTK with velocity + pressure field
!       write on unit 52 (vtk_xy.xxxxxxx.dat) where 
!                                  xxxxxxx is the timestep
!       the file is closed at the end of the subroutine
!     INPUTS
!       itime   ---> timestep
!       k0      ---> k coordinate of the 2D slice
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
        subroutine vtk_xy_bin(itime,k0)
!
        use storage
        implicit none
!
        integer :: i,j,itime
        integer :: k0
!	
        character(len=24) :: file_name
!
        real(sp) :: u,w,v,den
        real(sp) :: cte1
        real(sp), parameter :: zerol=0.0
!
        file_name = 'tec_xy.xxxxxxxx.vtk'
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
        open(52,file=file_name,status='unknown')
!
        write(52,'(A26)')'# vtk DataFile Version 2.0'
        write(52,'(A5)') 'Campo'
        write(52,'(A6)') 'BINARY'
        write(52,'(A25)')'DATASET STRUCTURED_POINTS'
        write(52,'(A11,I10,A1,I10,A1,I10)') 'DIMENSIONS ',l,' ',m,' ',1
        write(52,'(A7,I10,A1,I10,A1,I10)')  'ORIGIN ',1,' ',1,' ',k0        
        write(52,'(A8,I10,A1,I10,A1,I10)') 'SPACING ',1,' ',1,' ',1     
        write(52,'(A10,I10)')'POINT_DATA ',l*m*1
!
! vector
        write(52,'(A23)')'VECTORS velocity float'
        close(52)
!
! then write output (binary)
        open(52,file=file_name,status='old', position='append', & 
                form='unformatted',access='STREAM',CONVERT="BIG_ENDIAN")
!
        do j = 1,m
           do i = 1,l
              u = a01(i,j,k0)+a02(i,j,k0)+a03(i,j,k0)+a04(i,j,k0)+a05(i,j,k0) &
                 -a10(i,j,k0)-a11(i,j,k0)-a12(i,j,k0)-a13(i,j,k0)-a14(i,j,k0)
!
              w = a03(i,j,k0)+a07(i,j,k0)+a08(i,j,k0)+a09(i,j,k0)+a12(i,j,k0) &
                 -a01(i,j,k0)-a10(i,j,k0)-a16(i,j,k0)-a17(i,j,k0)-a18(i,j,k0)
!
              v = a04(i,j,k0)+a06(i,j,k0)+a07(i,j,k0)+a13(i,j,k0)+a18(i,j,k0) &
                 -a02(i,j,k0)-a09(i,j,k0)-a11(i,j,k0)-a15(i,j,k0)-a16(i,j,k0)
!
              den = a01(i,j,k0)+a02(i,j,k0)+a03(i,j,k0)+a04(i,j,k0) &
                   +a05(i,j,k0)+a06(i,j,k0)+a07(i,j,k0)+a08(i,j,k0) &
                   +a09(i,j,k0)+a10(i,j,k0)+a11(i,j,k0)+a12(i,j,k0) &
                   +a13(i,j,k0)+a14(i,j,k0)+a15(i,j,k0)+a16(i,j,k0) &
                   +a17(i,j,k0)+a18(i,j,k0)+a19(i,j,k0)+cte1
!
              write(52) u/den, w/den, v/den
           end do
        end do
        close(52)
!
#ifdef QQQQQQ
!  w (scalar)
        open(52,file=file_name,status='old', position='append')
        write(52,'(A21)')'SCALARS w float'
        write(52,'(A20)')'LOOKUP_TABLE default'
        close(52)
!        
        open(52,file=file_name,status='old', position='append', &
                form='unformatted',access='STREAM',CONVERT="BIG_ENDIAN")
        do j = 1,m
           do i = 1,l
!
              w = a03(i,j)+a08(i,j)+a12(i,j) &
                 -a01(i,j)-a10(i,j)-a17(i,j)
!
              den =( a01(i,j)+a03(i,j)+a05(i,j)+a08(i,j) &
                    +a10(i,j)+a12(i,j)+a14(i,j)+a17(i,j)+a19(i,j))+cte1
!
              write(52) w/den
           end do
        end do
        close(52)
!
!  u (scalar)
        open(52,file=file_name,status='old', position='append')
        write(52,'(A21)')'SCALARS u float'
        write(52,'(A20)')'LOOKUP_TABLE default'
        close(52)
!        
        open(52,file=file_name,status='old', position='append', & 
                form='unformatted',access='STREAM',CONVERT="BIG_ENDIAN")
        do j = 1,m
           do i = 1,l

              u = a01(i,j)+a03(i,j)+a05(i,j) &
                 -a10(i,j)-a12(i,j)-a14(i,j)
!
              den =( a01(i,j)+a03(i,j)+a05(i,j)+a08(i,j) &
                    +a10(i,j)+a12(i,j)+a14(i,j)+a17(i,j)+a19(i,j))+cte1
!
              write(52) u/den
           end do
        end do
        close(52)
!        
!  den (scalar)
        open(52,file=file_name,status='old', position='append')
        write(52,'(A21)')'SCALARS rho float'
        write(52,'(A20)')'LOOKUP_TABLE default'
        close(52)
!
        open(52,file=file_name,status='old', position='append', &
                form='unformatted',access='STREAM',CONVERT="BIG_ENDIAN")
        do j = 1,m
           do i = 1,l

              den =( a01(i,j)+a03(i,j)+a05(i,j)+a08(i,j) &
                    +a10(i,j)+a12(i,j)+a14(i,j)+a17(i,j)+a19(i,j))+cte1
!
              write(52) den
           end do
        end do
        close(52)
#endif
!
        write(16,*) "I/O: plane xy (vtk,binary) done", k0
        write(6,*)  "I/O: plane xy (vtk,binary) done", k0
!
#ifdef DEBUG_1
        if(myrank == 0) then
           write(6,*) "DEBUG1: Exiting from sub. vtk_xy_bin"
        endif
#endif
!
4000    format(i8.8)
!
       end subroutine vtk_xy_bin
