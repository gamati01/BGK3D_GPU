!====================================================
!     ****** LBE/col_MC
!
!     COPYRIGHT
!       (c) 2000-2011 by CASPUR/G.Amati
!       (c) 2013-20?? by CINECA/G.Amati
!     NAME
!       col_MC
!     DESCRIPTION
!       Collision according to bgk style (\omega(f_f^(eq)))
!       forcing is included and is proportional to u_0
!       MOVE + COLLIDE version
!     INPUTS
!       itime --> timestep 
!     TODO
!
!     NOTES
!       integer variables used: i,j,k,itime
!       ifdefs
!         * #ifdef FUSED
!         * #ifdef NOSHIFT
!         * #ifdef OFFLOAD
!         * #ifdef OPENACC
!         * #ifdef TGVFORCING
!         * #ifdef LES
!         * #ifdef CHANNEL
!         * #ifdef FORCING_Y
!         * #ifdef FORCING_Z
!         * #ifdef DEBUG_2
!         * #ifdef DEBUG_3
!
!     *****
!====================================================
!
        subroutine col_MC(itime)
!
        use storage
        use real_kinds
        use timing
#ifdef GPU_NATIVE
        use iso_c_binding
#if !defined(OPENACC)
        use omp_lib
#endif
#elif defined(OFFLOAD_KERNEL_SYNTAX)
        use omp_lib
#endif
!
        implicit none
!
        integer, intent(IN) :: itime
        integer             :: i,j,k
#ifdef OFFLOAD_KERNEL_SYNTAX
#  ifdef KS_BLOCK_3D
        integer :: ntx, nty, ntz, nbx, nby, nbz
#  else
        integer :: tid, ncell, nthrd, nblck
#  endif
#endif
!
        real(mykind) :: x,y,z
        real(mykind) :: pi
        real(mykind) :: x01,x02,x03,x04,x05,x06,x07
        real(mykind) :: x08,x09,x10,x11,x12,x13
        real(mykind) :: x14,x15,x16,x17,x18,x19
        real(mykind) :: e01,e02,e03,e04,e05,e06,e07
        real(mykind) :: e08,e09,e10,e11,e12,e13
        real(mykind) :: e14,e15,e16,e17,e18,e19
        real(mykind) :: n01,n02,n03,n04,n05,n06,n07
        real(mykind) :: n08,n09,n10,n11,n12,n13
        real(mykind) :: n14,n15,n16,n17,n18
        real(mykind) :: rho,rhoinv,vx,vy,vz
        real(mykind) :: vx2,vy2,vz2,vsq
        real(mykind) :: rp1,rp2,rp0
        real(mykind) :: vxpz,vxmz
        real(mykind) :: vxpy,vxmy
        real(mykind) :: vypz,vymz
        real(mykind) :: qx,qy,qz,q0
        real(mykind) :: qxpz,qxmz,qxpy,qxmy,qypz,qymz
        real(mykind) :: forcex, forcey, forcez
        real(mykind) :: cte1,cte0
        real(mykind) :: Pxx,Pyy,Pzz,Pxy,Pyx,Pxz,Pzx,Pyz,Pzy
        real(mykind) :: Ptotal,Ts
!        
                parameter(pi=3.141592653589793238462643383279)
!
#ifdef GPU_NATIVE
! native GPU (HIP/CUDA) port of the fused collision: interoperable kind
! for the storage arrays (must match mystorage in precision.F90).
#ifdef DOUBLE_P
        integer, parameter :: c_real = c_double
#else
        integer, parameter :: c_real = c_float
#endif
        integer :: nx_c, ny_c
! Device base addresses for the 19 input (a*) and 18 output (b*) arrays.
! They are rebuilt every step from the CURRENT host pointers so the
! a<->b pointer swap (see end of this routine) is always honoured.
        type(c_ptr) :: pa(19), pb(18)
        integer     :: idev, ii
        interface
          subroutine col_mc_gpu(                                          &
              a01,a02,a03,a04,a05,a06,a07,a08,a09,a10,                    &
              a11,a12,a13,a14,a15,a16,a17,a18,a19,                        &
              b01,b02,b03,b04,b05,b06,b07,b08,b09,b10,                    &
              b11,b12,b13,b14,b15,b16,b17,b18,                            &
              l,m,n,nx,ny,                                               &
              omega,cte0,cte1,p0,p1,p2,rf,qf,tre,                         &
              forcex,forcey,forcez) bind(C, name="col_mc_gpu")
            import :: c_real, c_int, c_ptr
! Device pointers passed BY VALUE (a C `real_t*` argument).  Using c_ptr
! instead of `dimension(*)` lets us hand over addresses we resolved
! ourselves (omp_get_mapped_ptr / acc_deviceptr / managed host==device),
! which is robust against the Fortran pointer swap.
            type(c_ptr), value :: a01,a02,a03,a04,a05,a06,a07,a08,a09,a10
            type(c_ptr), value :: a11,a12,a13,a14,a15,a16,a17,a18,a19
            type(c_ptr), value :: b01,b02,b03,b04,b05,b06,b07,b08,b09,b10
            type(c_ptr), value :: b11,b12,b13,b14,b15,b16,b17,b18
            integer(c_int), value :: l,m,n,nx,ny
            real(c_real), value :: omega,cte0,cte1,p0,p1,p2,rf,qf,tre
            real(c_real), value :: forcex,forcey,forcez
          end subroutine col_mc_gpu
        end interface
#endif

#ifdef FUSED
!
#ifdef DEBUG_3
        real(mykind) :: cte
        character(len=17) :: file_nameD
        file_nameD = 'debug.xxx.xxx.log'
        write(file_nameD(7:9),3300) itime
        write(file_nameD(11:13),3300) myrank
        open(41,file=file_nameD, status='unknown')        ! debug file
!
        call probe(itime,(3*l/4),(m/2))
#endif
!
#ifdef NOSHIFT
        cte1 = zero
#else
        cte1 = uno
#endif
        cte0 = uno - cte1
!
! initialize constant.....        
        forcex = zero
        forcey = zero
        forcez = zero
!
#ifdef GPU_NATIVE
! ---- native GPU (HIP/CUDA) collision -------------------------------
! Resolve the DEVICE base address of every storage array from its
! CURRENT host pointer and hand those device pointers to the HIP/CUDA
! launcher.  Doing the translation ourselves (instead of relying on the
! `use_device_addr` / `host_data` substitution) is mandatory here because
! the FUSED scheme swaps the a*<->b* Fortran pointers at the end of every
! step: a directive-based substitution does NOT track that swap and ends
! up handing the SAME device buffer for both a and b on alternate steps,
! aliasing the kernel input and output.  c_loc() always reflects the
! current association, so the resolved pointers are correct every step.
        nx_c = size(a01,1)
        ny_c = size(a01,2)
!
! current host base addresses (follow the a<->b swap automatically)
        pa( 1)=c_loc(a01); pa( 2)=c_loc(a02); pa( 3)=c_loc(a03)
        pa( 4)=c_loc(a04); pa( 5)=c_loc(a05); pa( 6)=c_loc(a06)
        pa( 7)=c_loc(a07); pa( 8)=c_loc(a08); pa( 9)=c_loc(a09)
        pa(10)=c_loc(a10); pa(11)=c_loc(a11); pa(12)=c_loc(a12)
        pa(13)=c_loc(a13); pa(14)=c_loc(a14); pa(15)=c_loc(a15)
        pa(16)=c_loc(a16); pa(17)=c_loc(a17); pa(18)=c_loc(a18)
        pa(19)=c_loc(a19)
        pb( 1)=c_loc(b01); pb( 2)=c_loc(b02); pb( 3)=c_loc(b03)
        pb( 4)=c_loc(b04); pb( 5)=c_loc(b05); pb( 6)=c_loc(b06)
        pb( 7)=c_loc(b07); pb( 8)=c_loc(b08); pb( 9)=c_loc(b09)
        pb(10)=c_loc(b10); pb(11)=c_loc(b11); pb(12)=c_loc(b12)
        pb(13)=c_loc(b13); pb(14)=c_loc(b14); pb(15)=c_loc(b15)
        pb(16)=c_loc(b16); pb(17)=c_loc(b17); pb(18)=c_loc(b18)
!
! translate host -> device addresses according to the residency model:
!   * CUDA / NVIDIA managed : host address IS the device address.
!   * CUDA / NVIDIA NOMANAGED: explicit `!$acc data` region -> device
!                             address via `host_data use_device` + c_loc.
!   * HIP  / AMD  (OpenMP)  : OpenMP target data -> omp_get_mapped_ptr,
!                             keyed by the CURRENT host address (swap-safe).
#ifdef OPENACC
#ifdef NOMANAGED
!$acc host_data use_device(                                             &
!$acc&   a01,a02,a03,a04,a05,a06,a07,a08,a09,a10,                        &
!$acc&   a11,a12,a13,a14,a15,a16,a17,a18,a19,                            &
!$acc&   b01,b02,b03,b04,b05,b06,b07,b08,b09,b10,                        &
!$acc&   b11,b12,b13,b14,b15,b16,b17,b18)
        pa( 1)=c_loc(a01); pa( 2)=c_loc(a02); pa( 3)=c_loc(a03)
        pa( 4)=c_loc(a04); pa( 5)=c_loc(a05); pa( 6)=c_loc(a06)
        pa( 7)=c_loc(a07); pa( 8)=c_loc(a08); pa( 9)=c_loc(a09)
        pa(10)=c_loc(a10); pa(11)=c_loc(a11); pa(12)=c_loc(a12)
        pa(13)=c_loc(a13); pa(14)=c_loc(a14); pa(15)=c_loc(a15)
        pa(16)=c_loc(a16); pa(17)=c_loc(a17); pa(18)=c_loc(a18)
        pa(19)=c_loc(a19)
        pb( 1)=c_loc(b01); pb( 2)=c_loc(b02); pb( 3)=c_loc(b03)
        pb( 4)=c_loc(b04); pb( 5)=c_loc(b05); pb( 6)=c_loc(b06)
        pb( 7)=c_loc(b07); pb( 8)=c_loc(b08); pb( 9)=c_loc(b09)
        pb(10)=c_loc(b10); pb(11)=c_loc(b11); pb(12)=c_loc(b12)
        pb(13)=c_loc(b13); pb(14)=c_loc(b14); pb(15)=c_loc(b15)
        pb(16)=c_loc(b16); pb(17)=c_loc(b17); pb(18)=c_loc(b18)
!$acc end host_data
#endif
#else
        idev = omp_get_default_device()
        do ii=1,19
           pa(ii) = omp_get_mapped_ptr(pa(ii), idev)
        enddo
        do ii=1,18
           pb(ii) = omp_get_mapped_ptr(pb(ii), idev)
        enddo
#endif
        call col_mc_gpu(                                                 &
             pa( 1),pa( 2),pa( 3),pa( 4),pa( 5),pa( 6),pa( 7),pa( 8),    &
             pa( 9),pa(10),pa(11),pa(12),pa(13),pa(14),pa(15),pa(16),    &
             pa(17),pa(18),pa(19),                                       &
             pb( 1),pb( 2),pb( 3),pb( 4),pb( 5),pb( 6),pb( 7),pb( 8),    &
             pb( 9),pb(10),pb(11),pb(12),pb(13),pb(14),pb(15),pb(16),    &
             pb(17),pb(18),                                              &
             int(l,c_int),int(m,c_int),int(n,c_int),                     &
             int(nx_c,c_int),int(ny_c,c_int),                            &
             real(omega,c_real),real(cte0,c_real),real(cte1,c_real),     &
             real(p0,c_real),real(p1,c_real),real(p2,c_real),            &
             real(rf,c_real),real(qf,c_real),real(tre,c_real),           &
             real(forcex,c_real),real(forcey,c_real),real(forcez,c_real))
#else
#ifdef OFFLOAD_KERNEL_SYNTAX
#  ifdef KS_BLOCK_3D
        ntx = 8
        nty = 8
        ntz = 4
        nbx = (l + ntx - 1)/ntx
        nby = (m + nty - 1)/nty
        nbz = (n + ntz - 1)/ntz
!$omp target teams parallel thread_limit(dims(3):ntx,nty,ntz) num_teams(dims(3):nbx,nby,nbz)
        i = 1 + omp_get_thread_num_dim(0) + omp_get_team_num_dim(0)*omp_get_num_threads_dim(0)
        j = 1 + omp_get_thread_num_dim(1) + omp_get_team_num_dim(1)*omp_get_num_threads_dim(1)
        k = 1 + omp_get_thread_num_dim(2) + omp_get_team_num_dim(2)*omp_get_num_threads_dim(2)
        if (i <= l .and. j <= m .and. k <= n) then
#  else
        nthrd = 256
        ncell = l*m*n
        nblck = (ncell + nthrd - 1)/nthrd
!$omp target teams parallel num_threads(dims(3):nthrd) num_teams(dims(3):nblck) thread_limit(nthrd)
        tid = omp_get_thread_num_dim(0) + omp_get_team_num_dim(0)*omp_get_num_threads_dim(0)
        if (tid < ncell) then
          i = mod(tid, l) + 1
          j = mod(tid/l, m) + 1
          k = tid/(l*m) + 1
#  endif
#elif defined(OFFLOAD)
!$OMP target teams distribute parallel do simd collapse(3)
        do k = 1,n
        do j = 1,m
        do i = 1,l
#elif OPENACC
!$acc parallel
!$acc loop independent collapse(3)
        do k = 1,n
        do j = 1,m
        do i = 1,l
#else
        do concurrent (k=1:n,j=1:m,i=1:l)
#endif
!
           x01 = a01(i-1,j+1,k  )
           x02 = a02(i-1,j  ,k+1)
           x03 = a03(i-1,j-1,k  )
           x04 = a04(i-1,j  ,k-1)
           x05 = a05(i-1,j  ,k  )
           x06 = a06(i  ,j  ,k-1)
           x07 = a07(i  ,j-1,k-1)
           x08 = a08(i  ,j-1,k  )
           x09 = a09(i  ,j-1,k+1)
           x10 = a10(i+1,j+1,k  )
           x11 = a11(i+1,j  ,k+1)
           x12 = a12(i+1,j-1,k  )
           x13 = a13(i+1,j  ,k-1)
           x14 = a14(i+1,j  ,k  )
           x15 = a15(i  ,j  ,k+1)
           x16 = a16(i  ,j+1,k+1)
           x17 = a17(i  ,j+1,k  )
           x18 = a18(i  ,j+1,k-1)
           x19 = a19(i  ,j  ,k  )

           rho = x01 + x02 + x03 + x04 + x05 + x06 + x07 + x08 &
                +x09 + x10 + x11 + x12 + x13 + x14 + x15 + x16 &
                +x17 + x18 + x19 + cte1
!
           rhoinv = uno/rho
!
           vx = (x01+x02+x03+x04+x05-x10-x11-x12-x13-x14)*rhoinv
           vy = (x03+x07+x08+x09+x12-x01-x10-x16-x17-x18)*rhoinv
           vz = (x04+x06+x07+x13+x18-x02-x09-x11-x15-x16)*rhoinv
!           
#ifdef TGVFORCING
!shift equilibrium velocity           
           z = (real(k,mykind)-0.5d0)/real(n,mykind)  ! 0<z<1
           y = (real(j,mykind)-0.5d0)/real(m,mykind)  ! 0<y<1
           x = (real(i,mykind)-0.5d0)/real(l,mykind)  ! 0<x<1
           vx = vx + 0.0001*(u00*sin(due*pi*x)*cos(due*pi*y)*cos(due*pi*z))
           vy = vy - 0.0001*(u00*cos(due*pi*x)*sin(due*pi*y)*cos(due*pi*z))
#endif
!
#ifdef DEBUG_3
           write(41,3131) i+offset(1),j+offset(2), obs(i,j),   & 
                        vx,vy,rho
!                       x01,x03,x05,x08,x10,x12,x14,x17,x19
#endif
!
! Quadratic terms
           vx2 = vx*vx
           vy2 = vy*vy
           vz2 = vz*vz

           vsq = vx2+vy2+vz2

           vxpy = vx+vy
           qxpy = cte0+qf*(tre*vxpy*vxpy-vsq)
           vxmy = vx-vy
           qxmy = cte0+qf*(tre*vxmy*vxmy-vsq)
           vxpz = vx+vz
           qxpz = cte0+qf*(tre*vxpz*vxpz-vsq)
           vxmz = vx-vz
           qxmz = cte0+qf*(tre*vxmz*vxmz-vsq)
           vypz = vy+vz
           qypz = cte0+qf*(tre*vypz*vypz-vsq)
           vymz = vy-vz
           qymz = cte0+qf*(tre*vymz*vymz-vsq)
           qx   = cte0+qf*(tre*vx2      -vsq)
           qy   = cte0+qf*(tre*vy2      -vsq)
           qz   = cte0+qf*(tre*vz2      -vsq)
           q0   = cte0+qf*(             -vsq)
!
! linear terms
           vx   = rf*vx
           vy   = rf*vy
           vz   = rf*vz
           vxpy = rf*vxpy
           vxmy = rf*vxmy
           vxpz = rf*vxpz
           vxmz = rf*vxmz
           vypz = rf*vypz
           vymz = rf*vymz
!
! constant terms
           rp0 = rho*p0
           rp1 = rho*p1
           rp2 = rho*p2
!
! equilibrium distribution
           e01 = rp2*(+vxmy+qxmy) + cte1*(rp2-p2)
           e02 = rp2*(+vxmz+qxmz) + cte1*(rp2-p2)
           e03 = rp2*(+vxpy+qxpy) + cte1*(rp2-p2)
           e04 = rp2*(+vxpz+qxpz) + cte1*(rp2-p2)
           e05 = rp1*(+vx  +qx  ) + cte1*(rp1-p1)
           e06 = rp1*(+vz  +qz  ) + cte1*(rp1-p1)
           e07 = rp2*(+vypz+qypz) + cte1*(rp2-p2)
           e08 = rp1*(+vy  +qy  ) + cte1*(rp1-p1)
           e09 = rp2*(+vymz+qymz) + cte1*(rp2-p2)
           e10 = rp2*(-vxpy+qxpy) + cte1*(rp2-p2)
           e11 = rp2*(-vxpz+qxpz) + cte1*(rp2-p2)
           e12 = rp2*(-vxmy+qxmy) + cte1*(rp2-p2)
           e13 = rp2*(-vxmz+qxmz) + cte1*(rp2-p2)
           e14 = rp1*(-vx  +qx  ) + cte1*(rp1-p1)
           e15 = rp1*(-vz  +qz  ) + cte1*(rp1-p1)
           e16 = rp2*(-vypz+qypz) + cte1*(rp2-p2)
           e17 = rp1*(-vy  +qy  ) + cte1*(rp1-p1)
           e18 = rp2*(-vymz+qymz) + cte1*(rp2-p2)
           e19 = rp0*(     +q0  ) + cte1*(rp0-p0)
!
#ifdef LES
! compute les
!
!non-equilibrium distribution
           n01 = x01-e01
           n02 = x02-e02
           n03 = x03-e03
           n04 = x04-e04
           n05 = x05-e05
           n06 = x06-e06
           n07 = x07-e07
           n08 = x08-e08
           n09 = x09-e09
           n10 = x10-e10
           n11 = x11-e11
           n12 = x12-e12
           n13 = x13-e13
           n14 = x14-e14
           n15 = x15-e15
           n16 = x16-e16
           n17 = x17-e17
           n18 = x18-e18
!
!
! compute Pij (six terms)
           Pxx = n01 + &
                 n02 + &
                 n03 + &
                 n04 + &
                 n05 + &
                 n10 + &
                 n11 + &
                 n12 + &
                 n13 + &
                 n14 
!
           Pyy = n01 + &
                 n03 + &
                 n07 + &
                 n08 + &
                 n09 + &
                 n10 + &
                 n12 + &
                 n16 + &
                 n17 + &
                 n18 
!
           Pzz = n02 + &
                 n04 + &
                 n06 + &
                 n07 + &
                 n09 + &
                 n11 + &
                 n13 + &
                 n15 + &
                 n16 + &
                 n18 
!
           Pxz = -n02 &
                 +n04 &
                 +n11 &
                 -n13 
!
           Pxy = -n01 &
                 +n03 &
                 +n10 &
                 -n12 
!
           Pyz = +n07 &
                 -n09 &
                 +n16 &
                 -n18 
!
           Pyx = Pxy
           Pzx = Pxz
           Pzy = Pyz
!           
! calculate Pi total
           Ptotal =sqrt((Pxx)**2 + (Pyy)**2 + (Pzz)**2 + &
                        (2.0*Pxy*Pyx)                  + &
                        (2.0*Pxz*Pzx)                  + &
                        (2.0*Pyz*Pzy))
!           
! adding turbulent viscosity
           Ts = 1/(2*omega1) + sqrt(18*(cteS)**2 *Ptotal+(1/omega1)**2)/2
           omega = 1/Ts
!
#endif                  
!
! forcing term
!
#ifdef CHANNEL
# ifdef FORCING_Y
        forcex = zero
        forcey = fgrad*rho
        forcez = zero
# elif FORCING_Z
        forcex = zero
        forcex = zero
        forcez = fgrad*rho
# else
        forcex = fgrad*rho
        forcey = zero
        forcez = zero
# endif
#else
! skip this subroutine (ORIGINAL version)
#endif
!
#ifdef TGVFORCING_TRY
           z = (real(k,mykind)-0.5d0)/real(n,mykind)  ! 0<z<1
           y = (real(j,mykind)-0.5d0)/real(m,mykind)  ! 0<y<1
           x = (real(i,mykind)-0.5d0)/real(l,mykind)  ! 0<x<1
           forcex = +0.00005*(2.0*u00*sin((2*pi*x)+(2*pi*y)+(2*pi*z)))
           forcey = -0.00005*(    u00*sin((2*pi*x)+(2*pi*y)+(2*pi*z)))
           forcez = -0.00005*(    u00*sin((2*pi*x)+(2*pi*y)+(2*pi*z)))
!           write(6,*) forcex, forcey, forcez
#endif
!
! loop on populations
        b01(i,j,k) = x01 - omega*(x01-e01) + forcex - forcey         
        b02(i,j,k) = x02 - omega*(x02-e02) + forcex          - forcez
        b03(i,j,k) = x03 - omega*(x03-e03) + forcex + forcey         
        b04(i,j,k) = x04 - omega*(x04-e04) + forcex          + forcez
        b05(i,j,k) = x05 - omega*(x05-e05) + forcex                  
        b06(i,j,k) = x06 - omega*(x06-e06)                   + forcez
        b07(i,j,k) = x07 - omega*(x07-e07)          + forcey + forcez
        b08(i,j,k) = x08 - omega*(x08-e08)          + forcey         
        b09(i,j,k) = x09 - omega*(x09-e09)          + forcey - forcez
        b10(i,j,k) = x10 - omega*(x10-e10) - forcex - forcey         
        b11(i,j,k) = x11 - omega*(x11-e11) - forcex          - forcez
        b12(i,j,k) = x12 - omega*(x12-e12) - forcex + forcey         
        b13(i,j,k) = x13 - omega*(x13-e13) - forcex          + forcez
        b14(i,j,k) = x14 - omega*(x14-e14) - forcex                  
        b15(i,j,k) = x15 - omega*(x15-e15)                   - forcez
        b16(i,j,k) = x16 - omega*(x16-e16)          - forcey - forcez
        b17(i,j,k) = x17 - omega*(x17-e17)          - forcey         
        b18(i,j,k) = x18 - omega*(x18-e18)          - forcey + forcez
        a19(i,j,k) = x19 - omega*(x19-e19)                           

#ifdef OFFLOAD_KERNEL_SYNTAX
        end if
!$omp end target teams parallel
#elif defined(OFFLOAD)
        end do
        end do
        end do
!$OMP end target teams distribute parallel do simd
#elif OPENACC
        end do
        end do
        end do
        !$acc end parallel
#else
        end do
#endif
#endif
!
! fix: swap populations (pointers)
        c01 => a01
        c02 => a02
        c03 => a03
        c04 => a04
        c05 => a05
        c06 => a06
        c07 => a07
        c08 => a08
        c09 => a09
        c10 => a10
        c11 => a11
        c12 => a12
        c13 => a13
        c14 => a14
        c15 => a15
        c16 => a16
        c17 => a17
        c18 => a18
!
! new ---> current
        a01 => b01
        a02 => b02
        a03 => b03
        a04 => b04
        a05 => b05
        a06 => b06
        a07 => b07
        a08 => b08
        a09 => b09
        a10 => b10
        a11 => b11
        a12 => b12
        a13 => b13
        a14 => b14
        a15 => b15
        a16 => b16
        a17 => b17
        a18 => b18
!
        b01 => c01
        b02 => c02
        b03 => c03
        b04 => c04
        b05 => c05
        b06 => c06
        b07 => c07
        b08 => c08
        b09 => c09
        b10 => c10
        b11 => c11
        b12 => c12
        b13 => c13
        b14 => c14
        b15 => c15
        b16 => c16
        b17 => c17
        b18 => c18
!
#ifdef DEBUG_3
!
! format
3300  format(i3.3)
!3131  format(3(i,1x),6(e14.6,1x))
3131  format(3(i,1x),3(e14.6,1x))
        close(41)
#endif
#endif
!
#ifdef DEBUG_2
        if(myrank == 0) then
           write(6,*) "DEBUG2: Exiting from sub. col_MC", &
                       forcex, forcey, forcez
        endif
#endif
        end subroutine col_MC
