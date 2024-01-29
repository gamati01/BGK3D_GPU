!=======================================================================
!     ****** LBE/vtk_om_xy_bin
!
!     COPYRIGHT
!       (c) 2009 by CASPUR/G.Amati
!       (c) 2000-2011 by CASPUR/G.Amati
!       (c) 2013-20?? by CINECA/G.Amati
!     NAME
!       vtk_om_bin
!     DESCRIPTION
!       Graphic subroutine:
!       write binary output for VTK with vorticity (stream field to do)
!       write on unit 55 (tec_om.xxxxxxx.dat) where 
!                                     xxxxxxx is the timestep 
!       the file is closed at the end of the subroutine
!     INPUTS
!       itime   ---> timestep
!     OUTPUT
!       none
!     TODO
!       
!     NOTES
!       character used:  file_name (19)
!       integer variables used: itime,i,k
!       max time allowed  99'999'999
!       max task allowed  9999
!       single precision only, to saving space..
!
!     *****
!=======================================================================
!
        subroutine vtk_om_xy_bin(itime,k0)
!
        use storage
        implicit none
!
        integer :: i,j,k
        integer, intent(in):: k0,itime
!	
        character(len=23) :: file_name
!
        real(sp) :: x01,x02,x03,x04,x05,x06,x07
        real(sp) :: x08,x09,x10,x11,x12,x13
        real(sp) :: x14,x15,x16,x17,x18,x19
!        
        real(sp), dimension(0:l+1,0:m+1,-1:+1) :: uu
        real(sp), dimension(0:l+1,0:m+1,-1:+1) :: vv
        real(sp), dimension(0:l+1,0:m+1,-1:+1) :: ww

        real(sp) :: rho, rhoinv
        real(sp) :: dudy, dudz
        real(sp) :: dvdx, dvdz
        real(sp) :: dwdx, dwdy
        real(sp) :: cte1
        real(sp) :: vorticity_x, vorticity_y, vorticity_z
!
        file_name = 'tec_om.xy.xxxxxxxx.vtk'
!
        myrank = 0
!
#ifdef NOSHIFT
       cte1 = zero
#else
       cte1 = uno
#endif
!
! precompute velocity
        do k = -1,1
           do j = 1,m
           do i = 1,l
!
           x01 = a01(i,j,k0+k)
           x02 = a02(i,j,k0+k)
           x03 = a03(i,j,k0+k)
           x04 = a04(i,j,k0+k)
           x05 = a05(i,j,k0+k)
           x06 = a06(i,j,k0+k)
           x07 = a07(i,j,k0+k)
           x08 = a08(i,j,k0+k)
           x09 = a09(i,j,k0+k)
           x10 = a10(i,j,k0+k)
           x11 = a11(i,j,k0+k)
           x12 = a12(i,j,k0+k)
           x13 = a13(i,j,k0+k)
           x14 = a14(i,j,k0+k)
           x15 = a15(i,j,k0+k)
           x16 = a16(i,j,k0+k)
           x17 = a17(i,j,k0+k)
           x18 = a18(i,j,k0+k)
           x19 = a19(i,j,k0+k)
!
           rho = x01 + x02 + x03 + x04 + x05 + x06 + x07 + x08 &
                +x09 + x10 + x11 + x12 + x13 + x14 + x15 + x16 &
                +x17 + x18 + x19 + cte1
!
           rhoinv = uno/rho
!
           uu(i,j,k) = (x01+x02+x03+x04+x05-x10-x11-x12-x13-x14)*rhoinv
           vv(i,j,k) = (x03+x07+x08+x09+x12-x01-x10-x16-x17-x18)*rhoinv
           ww(i,j,k) = (x04+x06+x07+x13+x18-x02-x09-x11-x15-x16)*rhoinv
!                 
!           write(6,*) u(i,j,k),v(i,j,k),w(i,j,k)
           enddo
           enddo
        enddo
!                 
        write(file_name(11:18),4000) itime
        open(55,file=file_name,status='unknown')
!
! vtk header
        write(55,'(A26)')'# vtk DataFile Version 2.0'
        write(55,'(A5)')'Campo'
        write(55,'(A6)')'BINARY'
        write(55,'(A25)')'DATASET STRUCTURED_POINTS'
        write(55,'(A11,I10,A1,I10,A1,I10)') &
                              'DIMENSIONS ',l-2,' ',m-2,' ',1
        write(55,'(A7,I10,A1,I10,A1,I10)')  'ORIGIN ',2,' ' &
                                                     ,2,' ' &
                                                     ,k0
        write(55,'(A8,I10,A1,I10,A1,I10)') 'SPACING ',1,' ',1,' ',1
        write(55,'(A10,I10)')'POINT_DATA ',(l-2)*(m-2)
!
! vorticity (vector)
        write(55,'(A24)')'VECTORS vorticity float'
        close(55)
!
! then write output (binary)
        open(55,file=file_name,status='old', position='append', &
                form='unformatted',access='STREAM',CONVERT="BIG_ENDIAN")
!
        k = 0
        do j = 2,m-1
           do i = 2,l-1
!
! delta_x = delta_y = delta_z = 2.0
!         
              dudy = 0.5*(uu(i,j+1,k)-uu(i,j-1,k))
              dudz = 0.5*(uu(i,j,k+1)-uu(i,j,k-1))
!              
              dvdx = 0.5*(vv(i+1,j,k)-vv(i-1,j,k))
              dvdz = 0.5*(vv(i,j,k+1)-vv(i,j,k-1))
!              
              dwdx = 0.5*(ww(i+1,j,k)-ww(i-1,j,k))
              dwdy = 0.5*(ww(i,j+1,k)-ww(i,j-1,k))
!              
              vorticity_x = dvdz - dwdy
              vorticity_y = dwdx - dudz
              vorticity_z = dudy - dvdx
!              
              write(55) vorticity_x, vorticity_y, vorticity_z
!              write(6,*) vorticity_x, vorticity_y, vorticity_z
!
           end do
        end do
        close(55)
!
        write(6,*)  "I/O: vorticity xy (vtk,binary) done", k0
        write(16,*) "I/O: vorticity xy (vtk,binary) done", k0
!
#ifdef DEBUG_1
        if(myrank == 0) then
           write(6,*) "DEBUG1: Exiting from sub. vtk_om_xy_bin"
        endif
#endif
!
4000    format(i8.8)
!
       end subroutine vtk_om_xy_bin

