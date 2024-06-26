! =====================================================================
!     ****** LBE/storage.f90
!
!     COPYRIGHT
!       (c) 2000-2008 by CASPUR/G.Amati
!     NAME
!       storage
!     DESCRIPTION
!       module for storage
!     INPUTS
!       none
!     OUTPUT
!       none
!     TODO
!
!     NOTES
!
!       integer variables defined:  l,n,l1,n1
!       real variables defined: lx, ly, dt, p0, p1, p2, rf, qf
!                               svisc, u0, omega,fgrad
!
!     *****
! =====================================================================
!
        module storage
!
        use real_kinds
!
        implicit none
!        
        integer:: lx, ly, lz          ! global size        
!
        integer:: l, m, n                    ! local (task) size
        integer:: l1, m1, n1
!
        integer, parameter::  mpid=3      ! mpi dimension
        integer, parameter::  zeroi=0      ! hint to help compiler...
        integer, parameter::  unoi=1       ! hint to help compiler...
!
!
! PGI doesn't support quad precision
#if defined PGI || defined ARM
        real(dp), parameter::  zero_qp=0.d0     ! hint to help compiler...
        real(dp), parameter::  uno_qp=1.d0      ! hint to help compiler...
        real(dp), parameter::  due_qp=2.d0      ! hint to help compiler...
        real(dp), parameter::  tre_qp=3.d0      ! hint to help compiler...
        real(dp), parameter::  qua_qp=4.d0      ! hint to help compiler...
!
        real(dp), parameter :: rf_qp = 3.d0
        real(dp), parameter :: qf_qp = 1.5d0
!
! 1/3
        real(dp), parameter :: p0_qp = 1.d0/3.d0
! 1/18
        real(dp), parameter :: p1_qp = 1.d0/18.d0
! 1/36
        real(dp), parameter :: p2_qp = uno_qp/(rf_qp*rf_qp*(uno_qp+tre_qp)) 
#else
        real(qp), parameter::  zero_qp=0.d0     ! hint to help compiler...
        real(qp), parameter::  uno_qp=1.d0      ! hint to help compiler...
        real(qp), parameter::  due_qp=2.d0      ! hint to help compiler...
        real(qp), parameter::  tre_qp=3.d0      ! hint to help compiler...
        real(qp), parameter::  qua_qp=4.d0      ! hint to help compiler...
!
        real(qp), parameter :: rf_qp = 3.d0
        real(qp), parameter :: qf_qp = 1.5d0
!
! 4/9
        real(qp), parameter :: p0_qp = 1.d0/3.d0
! 1/9
        real(qp), parameter :: p1_qp = 1.d0/18.d0
! 1/36
        real(qp), parameter :: p2_qp = uno_qp/(rf_qp*rf_qp*(uno_qp+tre_qp))
#endif
!
        integer:: myrank
        integer:: imax, imin                    ! obstacle stuff
        integer:: jmax, jmin                    ! obstacle stuff
        integer:: kmax, kmin                    ! obstacle stuff
        integer:: nobs                          ! #of obstacles per task
        integer:: offset(2)
        integer:: ipad,jpad,kpad
        integer:: flag1,flag2,flag3
!
        real(mykind), dimension(1:19) :: cx,cy,cz
        integer, dimension(1:19) :: icx,icy,icz
!
#ifdef FUSED
        real(mystorage), dimension(:,:,:), contiguous, pointer :: a01,a02,a03,a04,a05
        real(mystorage), dimension(:,:,:), contiguous, pointer :: a06,a07,a08,a09,a10
        real(mystorage), dimension(:,:,:), contiguous, pointer :: a11,a12,a13,a14,a15
        real(mystorage), dimension(:,:,:), contiguous, pointer :: a16,a17,a18,a19
!
        real(mystorage), dimension(:,:,:), contiguous, pointer :: b01,b02,b03,b04,b05
        real(mystorage), dimension(:,:,:), contiguous, pointer :: b06,b07,b08,b09,b10
        real(mystorage), dimension(:,:,:), contiguous, pointer :: b11,b12,b13,b14,b15
        real(mystorage), dimension(:,:,:), contiguous, pointer :: b16,b17,b18,b19
!
        real(mystorage), dimension(:,:,:), contiguous, pointer :: c01,c02,c03,c04,c05
        real(mystorage), dimension(:,:,:), contiguous, pointer :: c06,c07,c08,c09,c10
        real(mystorage), dimension(:,:,:), contiguous, pointer :: c11,c12,c13,c14,c15
        real(mystorage), dimension(:,:,:), contiguous, pointer :: c16,c17,c18,c19
!
#else
        !default (ORIGINAL)
        real(mystorage), dimension(:,:,:), allocatable :: a01,a02,a03,a04,a05
        real(mystorage), dimension(:,:,:), allocatable :: a06,a07,a08,a09,a10
        real(mystorage), dimension(:,:,:), allocatable :: a11,a12,a13,a14,a15
        real(mystorage), dimension(:,:,:), allocatable :: a16,a17,a18,a19
!
        real(mystorage), dimension(:,:,:), allocatable :: b01,b02,b03,b04,b05
        real(mystorage), dimension(:,:,:), allocatable :: b06,b07,b08,b09,b10
        real(mystorage), dimension(:,:,:), allocatable :: b11,b12,b13,b14,b15
        real(mystorage), dimension(:,:,:), allocatable :: b16,b17,b18,b19
#endif
!
        integer, dimension(:,:,:), allocatable :: obs
        integer:: radius
!
        real(mykind):: svisc, u0, u00, fgrad
        real(mykind):: u0x, u0y, u0z
        real(mykind):: u_inflow
        real(mykind):: cteS                ! Smagorinski constant  
        integer:: mydev, ndev              ! openacc variables
!
! correct casting
        real(mykind) :: omega
        real(mykind) :: omega1
        real(mykind), parameter :: zero = zero_qp
        real(mykind), parameter :: uno  = uno_qp
        real(mykind), parameter :: due  = due_qp
        real(mykind), parameter :: tre  = tre_qp
        real(mykind), parameter :: qua  = qua_qp
! 
        real(mykind), parameter :: rf = rf_qp
        real(mykind), parameter :: qf = qf_qp
        real(mykind), parameter :: p0 = p0_qp
        real(mykind), parameter :: p1 = p1_qp
        real(mykind), parameter :: p2 = p2_qp

        end module  storage
