!=====================================================================
!     ****** LBE/collision
!
!     COPYRIGHT
!       (c) 2000-2011 by CASPUR/G.Amati
!       (c) 2013-20?? by CINECA/G.Amati
!     NAME
!       collision
!     DESCRIPTION
!       Simple wrapper for movef routine...
!       - none (if -DFUSED) 
!       - movef else
!     INPUTS
!       itime
!     OUTPUT
!       none
!     TODO
!       
!     NOTES
!
!     *****
!=====================================================================
!
        subroutine collision(itime)
!
        use timing
        use storage
#if defined(PROFILING) && defined(GPU_NATIVE)
        use bcond_gpu_mod, only : gpu_device_sync
#endif
!
        implicit none
!
        integer, INTENT(in) :: itime
!
#ifdef PROFILING
! start timing
        call SYSTEM_CLOCK(countC0, count_rate, count_max)
        call time(tcountC0)
#endif
!
#ifdef FUSED
        call col_MC(itime)
#else
        call col(itime)
#endif
!
#ifdef PROFILING
! stop timing
#ifdef GPU_NATIVE
! flush the async native kernels so the CPU timer captures real GPU time
        call gpu_device_sync()
#endif
        call time(tcountC1)
        call SYSTEM_CLOCK(countC1, count_rate, count_max)
        time_coll = time_coll + real(countC1-countC0)/(count_rate)
        time_coll1 = time_coll1 + (tcountC1-tcountC0)
#endif
!
#ifdef DEBUG_2
        if(myrank == 0) then
           write(6,*) "DEBUG2: Exiting from sub. collision"
        endif
#endif
!
        end subroutine collision
