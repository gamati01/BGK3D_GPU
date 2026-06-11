! =====================================================================
!     ****** LBE/bgk3D
!
!     COPYRIGHT
!       (c) 2021 by CINECA/G.Amati
!     NAME
!       bgk2d
!     DESCRIPTION
!       main program for LBM 3D
!     INPUTS
!       none
!     OUTPUT
!       none
!     TODO
!
!     NOTES
!
!       integer variables used: itfin, itstart, ivtim
!                               itime, itsave, icheck, itrestart
!                               isignal
!       real variables used: tempo1, tempo2     
!       open the following unit: 16 (bgk.log)
!                                60 (prof_k.dat)
!                                61 (prof_i.dat)
!                                62 (u_med.dat)
!                                63 (diagno.dat)
!                                68 (probe.dat)
!                                69 (bgk.perf)
!       velocity directions:
!        direction  1    unit vector = ( 1,-1, 0)   tipo 2
!        direction  2    unit vector = ( 1, 0,-1)   tipo 2
!        direction  3    unit vector = ( 1, 1, 0)   tipo 2
!        direction  4    unit vector = ( 1, 0, 1)   tipo 2
!        direction  5    unit vector = ( 1, 0, 0)   tipo 1
!        direction  6    unit vector = ( 0, 0, 1)   tipo 1
!        direction  7    unit vector = ( 0, 1, 1)   tipo 2
!        direction  8    unit vector = ( 0, 1, 0)   tipo 1
!        direction  9    unit vector = ( 0, 1,-1)   tipo 2
!        direction 10    unit vector = (-1,-1, 0)   tipo 2
!        direction 11    unit vector = (-1, 0,-1)   tipo 2
!        direction 12    unit vector = (-1, 1, 0)   tipo 2
!        direction 13    unit vector = (-1, 0, 1)   tipo 2
!        direction 14    unit vector = (-1, 0, 0)   tipo 1
!        direction 15    unit vector = ( 0, 0,-1)   tipo 1
!        direction 16    unit vector = ( 0,-1,-1)   tipo 2
!        direction 17    unit vector = ( 0,-1, 0)   tipo 1
!        direction 18    unit vector = ( 0,-1, 1)   tipo 2
!        direction 19    unit vector = ( 0, 0, 0)   tipo 0
!                              
!     *****
! =====================================================================
!
      program bgk2d
!
      use storage
      use timing
      use real_kinds
      use unix_time_mod
!
      implicit none
!
      INTEGER:: itfin, itstart, ivtim
      INTEGER:: itime, itsave, icheck, itrestart, init_v
      INTEGER:: isignal
!
      call SYSTEM_CLOCK(countH0, count_rate, count_max)
      call time(tcountH0)
!      
! set up the simulation...
      call setup(itfin,ivtim,isignal,itsave,icheck,itrestart, & 
                 init_v)
!      
! initialize the flow...
      call initialize(itrestart,init_v,itfin,itstart,ivtim,isignal, & 
                      itsave,icheck)
!
!
! timing just before the time loop
      utime1 = unix_time()
!
! get GPU energy at start
#ifdef ENERGY
!
      mydev0 = 0
      mydev1 = 1
      mydev2 = 2
      mydev3 = 3
!
      mydev0_c=int(mydev0,kind=c_int)
      mydev1_c=int(mydev1,kind=c_int)
      mydev2_c=int(mydev2,kind=c_int)
      mydev3_c=int(mydev3,kind=c_int)
!
      ierrc=get_gpu_energy_mJ_u64(mydev0_c,energy0_1)
      if (ierrc /= 0_c_int) then
         write(*,*) "NVML error reading energy for device 0, err=", ierrc
      endif
!
      ierrc=get_gpu_energy_mJ_u64(mydev1_c,energy1_1)
      if (ierrc /= 0_c_int) then
         write(*,*) "NVML error reading energy for device 1, err=", ierrc
      endif
!
      ierrc=get_gpu_energy_mJ_u64(mydev2_c,energy2_1)
      if (ierrc /= 0_c_int) then
         write(*,*) "NVML error reading energy for device 2, err=", ierrc
      endif
!
      ierrc=get_gpu_energy_mJ_u64(mydev3_c,energy3_1)
      if (ierrc /= 0_c_int) then
        write(*,*) "NVML error reading energy for device 3, err=", ierrc
      endif
!
      write( 6,*) "INFO: Start Energy GPU0 =", utime1, energy0_1
      write( 6,*) "INFO: Start Energy GPU1 =", utime1, energy1_1
      write( 6,*) "INFO: Start Energy GPU2 =", utime1, energy2_1
      write( 6,*) "INFO: Start Energy GPU3 =", utime1, energy3_1
      write(16,*) "INFO: Start Energy GPU0 =", utime1, energy0_1
      write(16,*) "INFO: Start Energy GPU1 =", utime1, energy1_1
      write(16,*) "INFO: Start Energy GPU2 =", utime1, energy2_1
      write(16,*) "INFO: Start Energy GPU3 =", utime1, energy3_1
#endif
!
#ifdef NOMANAGED
!$acc data copyin(a01,a02,a03,a04,a05,a06,a07,a08,a09,a10,   &
!$acc&            a11,a12,a13,a14,a15,a16,a17,a18,a19,       &
!$acc&            b01,b02,b03,b04,b05,b06,b07,b08,b09,b10,   &
!$acc&            b11,b12,b13,b14,b15,b16,b17,b18,b19,obs)   
#endif
!
      call SYSTEM_CLOCK(countH1, count_rate, count_max)
      call time(tcountH1)
      time_init =  real(countH1-countH0)/(count_rate)
      time_init1 = (tcountH1-tcountH0)
!      
      call SYSTEM_CLOCK(countE0, count_rate, count_max)
      call time(tcountE0)
!
      call SYSTEM_CLOCK(countD0, count_rate, count_max)
      call time(tcountD0)
!
! main loop starts here.....
!
      if(myrank==0) then 
         call system("date       >> time.log")
      endif
!
#ifdef OFFLOAD
!$OMP target data map(tofrom:a01,a02,a03,a04,a05,a06,a07,a08,a09,a10, &
!$OMP&                       a11,a12,a13,a14,a15,a16,a17,a18,a19,     &
!$OMP&                       b01,b02,b03,b04,b05,b06,b07,b08,b09,b10, &
!$OMP&                       b11,b12,b13,b14,b15,b16,b17,b18,b19,obs)
#endif
!      
      do itime=itstart+1,itfin
!
#ifdef DEBUG_2
         if(myrank == 0) then
            write(6,*) "DEBUG2: starting time step =", itime
         endif
#endif
!
         call boundaries         ! boundary conditions
         call propagation        ! propagation step
         call collision(itime)   ! collision step
!
! get macroscopic values
         call diagnostic(itime,ivtim,icheck,itsave)
!
! get timing/profiling values
         if (mod(itime,isignal).eq.0) then
            if (myrank == 0 ) then 
               call profile(itime,itfin,isignal) 
            endif
         endif
      enddo
!
!
! timing just after the time loop
      utime2 = unix_time()
!      
!     some global timings
      call SYSTEM_CLOCK(countE1, count_rate, count_max)
      call time(tcountE1)
      time_loop = real(countE1-countE0)/(count_rate)
      time_loop1 = tcountE1-tcountE0
!
! get GPU energy at end
#ifdef ENERGY
      ierrc=get_gpu_energy_mJ_u64(mydev0_c,energy0_2)
      ierrc=get_gpu_energy_mJ_u64(mydev1_c,energy1_2)
      ierrc=get_gpu_energy_mJ_u64(mydev2_c,energy2_2)
      ierrc=get_gpu_energy_mJ_u64(mydev3_c,energy3_2)
!
      write(6,*)  "INFO: End Energy GPU0 =", utime2, energy0_2
      write(6,*)  "INFO: End Energy GPU1 =", utime2, energy1_2
      write(6,*)  "INFO: End Energy GPU2 =", utime2, energy2_2
      write(6,*)  "INFO: End Energy GPU3 =", utime2, energy3_2
      write(16,*) "INFO: End Energy GPU0 =", utime2, energy0_2
      write(16,*) "INFO: End Energy GPU1 =", utime2, energy1_2
      write(16,*) "INFO: End Energy GPU2 =", utime2, energy2_2
      write(16,*) "INFO: End Energy GPU3 =", utime2, energy3_2
#endif
!
!     final diagnostic (for check)
      call diagno(itime-1)
      call varm(itime-1)
      call prof_i(itime-1,m/2,n/2)
      call prof_j(itime-1,l/2,n/2)
!      call vtk_3d_bin(itime-1)
!
#ifdef OFFLOAD
!$omp end target data
#endif
!      
#ifdef NOMANAGED
!$acc end data
#endif
!
      call finalize(itstart,itfin)     ! finalize all
!
      if(myrank==0) then
         call system("date       >> time.log")
         write(6,*) "That's all folks!!!!"
      endif
!
      end program bgk2d
