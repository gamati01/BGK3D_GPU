!=====================================================================
!     ****** LBE/alloca
!
!     COPYRIGHT
!       (c) 2021 by CINECA/G.Amati
!     NAME
!       alloca
!     DESCRIPTION
!
!     INPUTS
!       none
!     OUTPUT
!       none
!     TODO
!       
!     NOTES
!
!     *****
!=====================================================================
!
      subroutine alloca()
!            
      use storage 
      use timing
! 
      implicit none
!
!     define unit vector  (real/integer)
      cx(  1    ) =  1
      cx(  2    ) =  1
      cx(  3    ) =  1
      cx(  4    ) =  1
      cx(  5    ) =  1
      cx(  6    ) =  0
      cx(  7    ) =  0
      cx(  8    ) =  0
      cx(  9    ) =  0
      cx( 10    ) = -1
      cx( 11    ) = -1
      cx( 12    ) = -1
      cx( 13    ) = -1
      cx( 14    ) = -1
      cx( 15    ) =  0
      cx( 16    ) =  0
      cx( 17    ) =  0
      cx( 18    ) =  0
      cx( 19    ) =  0

      icx(  1    ) =  1
      icx(  2    ) =  1
      icx(  3    ) =  1
      icx(  4    ) =  1
      icx(  5    ) =  1
      icx(  6    ) =  0
      icx(  7    ) =  0
      icx(  8    ) =  0
      icx(  9    ) =  0
      icx( 10    ) = -1
      icx( 11    ) = -1
      icx( 12    ) = -1
      icx( 13    ) = -1
      icx( 14    ) = -1
      icx( 15    ) =  0
      icx( 16    ) =  0
      icx( 17    ) =  0
      icx( 18    ) =  0
      icx( 19    ) =  0
!
      cy(  1    ) = -1
      cy(  2    ) =  0
      cy(  3    ) =  1
      cy(  4    ) =  0
      cy(  5    ) =  0
      cy(  6    ) =  0
      cy(  7    ) =  1
      cy(  8    ) =  1
      cy(  9    ) =  1
      cy( 10    ) = -1
      cy( 11    ) =  0
      cy( 12    ) =  1
      cy( 13    ) =  0
      cy( 14    ) =  0
      cy( 15    ) =  0
      cy( 16    ) = -1
      cy( 17    ) = -1
      cy( 18    ) = -1
      cy( 19    ) =  0
!      
      icy(  1    ) = -1
      icy(  2    ) =  0
      icy(  3    ) =  1
      icy(  4    ) =  0
      icy(  5    ) =  0
      icy(  6    ) =  0
      icy(  7    ) =  1
      icy(  8    ) =  1
      icy(  9    ) =  1
      icy( 10    ) = -1
      icy( 11    ) =  0
      icy( 12    ) =  1
      icy( 13    ) =  0
      icy( 14    ) =  0
      icy( 15    ) =  0
      icy( 16    ) = -1
      icy( 17    ) = -1
      icy( 18    ) = -1
      icy( 19    ) =  0
!
      !
      cz(  1    ) =  0
      cz(  2    ) = -1
      cz(  3    ) =  0
      cz(  4    ) =  1
      cz(  5    ) =  0
      cz(  6    ) =  1
      cz(  7    ) =  1
      cz(  8    ) =  0
      cz(  9    ) = -1
      cz( 10    ) =  0
      cz( 11    ) = -1
      cz( 12    ) =  0
      cz( 13    ) =  1
      cz( 14    ) =  0
      cz( 15    ) = -1
      cz( 16    ) = -1
      cz( 17    ) =  0
      cz( 18    ) =  1
      cz( 19    ) =  0
!
      icz(  1    ) =  0
      icz(  2    ) = -1
      icz(  3    ) =  0
      icz(  4    ) =  1
      icz(  5    ) =  0
      icz(  6    ) =  1
      icz(  7    ) =  1
      icz(  8    ) =  0
      icz(  9    ) = -1
      icz( 10    ) =  0
      icz( 11    ) = -1
      icz( 12    ) =  0
      icz( 13    ) =  1
      icz( 14    ) =  0
      icz( 15    ) = -1
      icz( 16    ) = -1
      icz( 17    ) =  0
      icz( 18    ) =  1
      icz( 19    ) =  0
!
      allocate(obs(1:l,1:m,1:n))
      obs = 0
!
      allocate(a01(0:l+1+ipad,0:m+1+jpad,0:n+1+kpad))
      allocate(a02(0:l+1+ipad,0:m+1+jpad,0:n+1+kpad))
      allocate(a03(0:l+1+ipad,0:m+1+jpad,0:n+1+kpad))
      allocate(a04(0:l+1+ipad,0:m+1+jpad,0:n+1+kpad))
      allocate(a05(0:l+1+ipad,0:m+1+jpad,0:n+1+kpad))
      allocate(a06(0:l+1+ipad,0:m+1+jpad,0:n+1+kpad))
      allocate(a07(0:l+1+ipad,0:m+1+jpad,0:n+1+kpad))
      allocate(a08(0:l+1+ipad,0:m+1+jpad,0:n+1+kpad))
      allocate(a09(0:l+1+ipad,0:m+1+jpad,0:n+1+kpad))
      allocate(a10(0:l+1+ipad,0:m+1+jpad,0:n+1+kpad))
      allocate(a11(0:l+1+ipad,0:m+1+jpad,0:n+1+kpad))
      allocate(a12(0:l+1+ipad,0:m+1+jpad,0:n+1+kpad))
      allocate(a13(0:l+1+ipad,0:m+1+jpad,0:n+1+kpad))
      allocate(a14(0:l+1+ipad,0:m+1+jpad,0:n+1+kpad))
      allocate(a15(0:l+1+ipad,0:m+1+jpad,0:n+1+kpad))
      allocate(a16(0:l+1+ipad,0:m+1+jpad,0:n+1+kpad))
      allocate(a17(0:l+1+ipad,0:m+1+jpad,0:n+1+kpad))
      allocate(a18(0:l+1+ipad,0:m+1+jpad,0:n+1+kpad))
      allocate(a19(0:l+1+ipad,0:m+1+jpad,0:n+1+kpad))
!
      allocate(b01(0:l+1+ipad,0:m+1+jpad,0:n+1+kpad))
      allocate(b02(0:l+1+ipad,0:m+1+jpad,0:n+1+kpad))
      allocate(b03(0:l+1+ipad,0:m+1+jpad,0:n+1+kpad))
      allocate(b04(0:l+1+ipad,0:m+1+jpad,0:n+1+kpad))
      allocate(b05(0:l+1+ipad,0:m+1+jpad,0:n+1+kpad))
      allocate(b06(0:l+1+ipad,0:m+1+jpad,0:n+1+kpad))
      allocate(b07(0:l+1+ipad,0:m+1+jpad,0:n+1+kpad))
      allocate(b08(0:l+1+ipad,0:m+1+jpad,0:n+1+kpad))
      allocate(b09(0:l+1+ipad,0:m+1+jpad,0:n+1+kpad))
      allocate(b10(0:l+1+ipad,0:m+1+jpad,0:n+1+kpad))
      allocate(b11(0:l+1+ipad,0:m+1+jpad,0:n+1+kpad))
      allocate(b12(0:l+1+ipad,0:m+1+jpad,0:n+1+kpad))
      allocate(b13(0:l+1+ipad,0:m+1+jpad,0:n+1+kpad))
      allocate(b14(0:l+1+ipad,0:m+1+jpad,0:n+1+kpad))
      allocate(b15(0:l+1+ipad,0:m+1+jpad,0:n+1+kpad))
      allocate(b16(0:l+1+ipad,0:m+1+jpad,0:n+1+kpad))
      allocate(b17(0:l+1+ipad,0:m+1+jpad,0:n+1+kpad))
      allocate(b18(0:l+1+ipad,0:m+1+jpad,0:n+1+kpad))
      allocate(b19(0:l+1+ipad,0:m+1+jpad,0:n+1+kpad))
!
      a01 = huge(mystorage)
      a02 = huge(mystorage)
      a03 = huge(mystorage)
      a04 = huge(mystorage)
      a05 = huge(mystorage)
      a06 = huge(mystorage)
      a07 = huge(mystorage)
      a08 = huge(mystorage)
      a09 = huge(mystorage)
      a10 = huge(mystorage)
      a11 = huge(mystorage)
      a12 = huge(mystorage)
      a13 = huge(mystorage)
      a14 = huge(mystorage)
      a15 = huge(mystorage)
      a16 = huge(mystorage)
      a17 = huge(mystorage)
      a18 = huge(mystorage)
      a19 = huge(mystorage)
!
      b01 = huge(mystorage)
      b02 = huge(mystorage)
      b03 = huge(mystorage)
      b04 = huge(mystorage)
      b05 = huge(mystorage)
      b06 = huge(mystorage)
      b07 = huge(mystorage)
      b08 = huge(mystorage)
      b09 = huge(mystorage)
      b10 = huge(mystorage)
      b11 = huge(mystorage)
      b12 = huge(mystorage)
      b13 = huge(mystorage)
      b14 = huge(mystorage)
      b15 = huge(mystorage)
      b16 = huge(mystorage)
      b17 = huge(mystorage)
      b18 = huge(mystorage)
      b19 = huge(mystorage)
!
#ifdef FUSED
      c01 => null()
      c02 => null()
      c03 => null()
      c04 => null()
      c05 => null()
      c06 => null()
      c07 => null()
      c08 => null()
      c09 => null()
      c10 => null()
      c11 => null()
      c12 => null()
      c13 => null()
      c14 => null()
      c15 => null()
      c16 => null()
      c17 => null()
      c18 => null()
      c19 => null()
#endif
!
#ifdef PGI
! do nothing
#else
      mem_start = get_mem()
#endif

#ifdef DEBUG_1
      if(myrank == 0) then
         write(6,*) "DEBUG1: Exiting from sub. alloca"
      endif
#endif
!
# ifdef MEM_CHECK
      if(myrank == 0) then
         mem_stop = get_mem();
         write(6,*) "MEM_CHECK: after sub. alloca mem =", mem_stop
      endif
# endif

      end subroutine alloca

