c svn $Id$
c
c mex routine implementing rect.f as taken from Ives, D.C. and
c R. M. Zacharias "Conformal mapping and Orthogonal Grid
c Generation", aiaa-87-2057.
c
c
c 
c The calling sequence from MATLAB should be
c
c   >> zn = mexrect ( z, n, n1, n2, n3, n4 );
c
c So 
c    1)  nlhs = 1
c    2)  plhs = 
c    3)  nrhs = 6
c    4)  prhs(1) = z
c        prhs(2) = n
c        prhs(3) = n1
c        prhs(4) = n2
c        prhs(5) = n3
c        prhs(6) = n4
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine mexFunction ( nlhs, plhs, nrhs, prhs )

#if defined (__alpha__) || defined (__x86_64__)
      integer*8 plhs(*), prhs(*)
#else 
      integer plhs(*), prhs(*)
#endif

#if defined (__alpha__) || defined (__x86_64__)
      integer*8 mxcreatedoublematrix, mxgetpr
#endif


      integer nlhs, nrhs

c
c     Size of input arguments
      integer size_m, size_n

c
c     Grid outline to be transformed.
      Complex*16 z(1000)

c
c     Pointers to real and imaginary parts of z input argument,
c     which is prhs(1)
#if defined (__alpha__) || defined (__x86_64__)
      integer*8 zpr, zpi
#else 
      integer zpr, zpi
#endif

c
c     Pointers to real parts of n, n[1234] input arguments, which
c     are prhs(2) thru prhs(6)
#if defined (__alpha__) || defined (__x86_64__)
      integer*8 np, n1p, n2p, n3p, n4p
#else 
      integer np, n1p, n2p, n3p, n4p
#endif

c
c     Size of input z in terms of matlab dimensions, i.e.,
c     #rows x #columns
      integer zsize

c 
c     These point to the real and imaginary parts of the 
c     matlab structure plhs(1)
#if defined (__alpha__) || defined (__x86_64__)
      integer*8 zlpr, zlpi
#else 
      integer zlpr, zlpi
#endif


c     We apparently have to copy the matlab arguments to
c     real arrays first, then cast them to integers.
      real*8 nr(1), n1r(1), n2r(1), n3r(1), n4r(1)
      integer n, n1, n2, n3, n4



c     Check for proper number of arguments.
      if ( nrhs .ne. 6 ) then
         call mexErrMsgTxt('mexrect requires 6 arguments')
      elseif ( nlhs .ne. 1 ) then
         call mexErrMsgTxt ( 'Only one output argument allowed' )
      endif


c     Check to see that the first input was complex.
c     This is the "Z" argument.
      if ( mxIsComplex(prhs(1)) .ne. 1 ) then
        call mexErrMsgTxt 
     +    ( 'First argument to mexrect must be complex' )
      endif


c     Check to see that the 2nd thru 6th arguments are numeric.
c     This is "N", "N1", "N2", "N3", and "N4".
      if ( mxIsNumeric(prhs(2)) .ne. 1 ) then
        call mexErrMsgTxt 
     +    ( '2nd argument to mexrect must be numeric' )
      elseif ( mxIsNumeric(prhs(3)) .ne. 1 ) then
        call mexErrMsgTxt 
     +    ( '3rd argument to mexrect must be numeric' )
      elseif ( mxIsNumeric(prhs(4)) .ne. 1 ) then
        call mexErrMsgTxt 
     +    ( '4th argument to mexrect must be numeric' )
      elseif ( mxIsNumeric(prhs(5)) .ne. 1 ) then
        call mexErrMsgTxt 
     +    ( '5th argument to mexrect must be numeric' )
      elseif ( mxIsNumeric(prhs(6)) .ne. 1 ) then
        call mexErrMsgTxt 
     +    ( '6th argument to mexrect must be numeric' )
      endif

c
c     Check to see that the 2nd thru 6th arguments are scalars.
      size_m = mxGetM ( prhs(2) )
      size_n = mxGetN ( prhs(2) )
      if ( size_n .ne. 1 .or. size_m .ne. 1 ) then
        call mexErrMsgTxt ('2nd argument to mexrect must be a scalar' )
      endif

      size_m = mxGetM ( prhs(3) )
      size_n = mxGetN ( prhs(3) )
      if ( size_n .ne. 1 .or. size_m .ne. 1 ) then
        call mexErrMsgTxt ('3rd argument to mexrect must be a scalar' )
      endif

      size_m = mxGetM ( prhs(4) )
      size_n = mxGetN ( prhs(4) )
      if ( size_n .ne. 1 .or. size_m .ne. 1 ) then
        call mexErrMsgTxt ('4th argument to mexrect must be a scalar' )
      endif

      size_m = mxGetM ( prhs(5) )
      size_n = mxGetN ( prhs(5) )
      if ( size_n .ne. 1 .or. size_m .ne. 1 ) then
        call mexErrMsgTxt ('5th argument to mexrect must be a scalar' )
      endif

      size_m = mxGetM ( prhs(6) )
      size_n = mxGetN ( prhs(6) )
      if ( size_n .ne. 1 .or. size_m .ne. 1 ) then
        call mexErrMsgTxt ('6th argument to mexrect must be a scalar' )
      endif



c
c     Create matrix for the return argument.
      size_m = mxGetM ( prhs(1) )
      size_n = mxGetN ( prhs(1) )
#ifdef VERBOSE
      write (*,*) ' size_m is ', size_m 
      write (*,*) ' size_n is ', size_n
#endif
      zsize = size_m*size_n
#if defined ( MATLAB_R11 ) || defined ( MATLAB_R12 )
      plhs(1) = mxCreateFull(size_m, size_n, 1)
#else 
      plhs(1) = mxCreateDoubleMatrix(size_m, size_n, 1)
#endif
      zlpr = mxgetpr(plhs(1))
      zlpi = mxGetPi(plhs(1))


c
c     Load the data into Fortran matrices.
      zpr = mxGetPr(prhs(1))
      zpi = mxGetPi(prhs(1))
      call mxCopyPtrToComplex16(zpr,zpi,z,zsize)

      np = mxGetPr(prhs(2))
      call mxCopyPtrToReal8(np,nr,1)
      n = nr(1)

      n1p = mxGetPr(prhs(3))
      call mxCopyPtrToReal8(n1p,n1r,1)
      n1 = n1r(1)

      n2p = mxGetPr(prhs(4))
      call mxCopyPtrToReal8(n2p,n2r,1)
      n2 = n2r(1)

      n3p = mxGetPr(prhs(5))
      call mxCopyPtrToReal8(n3p,n3r,1)
      n3 = n3r(1)

      n4p = mxGetPr(prhs(6))
      call mxCopyPtrToReal8(n4p,n4r,1)
      n4 = n4r(1)



c
c     Call the computational subroutine.
      call rect ( z, n, n1, n2, n3, n4 )


c
c     Load the output into a MATLAB array.
      call mxCopyComplex16ToPtr(z,zlpr,zlpi,zsize)
      return
      end





c
      Subroutine RECT(Z,N,N1,N2,N3,N4)
c
c     this subroutine is taken from ives,d.c. and
c     r.m.zacharias "conformal mapping and orthogonal grid generation"
c     aiaa-87-2057.
c
c     The only modification is the addition of the "tracking_error"
c     warning message.  This allows matlab to gracefully exit a call
c     to RECT that is going badly.
c
      Implicit Double Precision (A-H,O-Z)
      Complex*16 Z(1), Z0, ZD
      Double Precision R(1000), T(1000)
      PI = 4.D0 * DATAN(1.D0)
      Do 70 I = 1, N
        IM = N - MOD(N-I+1,N)
        IP = 1 + MOD(I,N)
        ZD = (Z(IM)-Z(I)) / (Z(IP)-Z(I))
        ALPHA = DATAN2(DIMAG(ZD),DREAL(ZD))
        If (ALPHA.LT.0) ALPHA = ALPHA + PI + PI
        PWR = PI / ALPHA
        If (I.EQ.N1.OR.I.EQ.N2.OR.I.EQ.N3.OR.I.EQ.N4) Then
          ZD = (Z(IM)-Z(I)) / CDABS(Z(IM)-Z(I))
          Do 10 J = 1, N
            Z(J) = DCMPLX(0.D0,1.D0) * Z(J) / ZD
   10     Continue
          PWR = PWR / 2.
        End If
        H1IN = 100.
        H1AX = -100.
        TP = 0.
        Do 40 J = 2, N
          ZD = Z(MOD(J+I-2,N)+1) - Z(I)
          R(J) = CDABS(ZD)
          T(J) = DATAN2(DIMAG(ZD),DREAL(ZD)) - PI - PI - PI - PI - PI -
     *        PI
          Do 20 K = 1, 7
            If (DABS(T(J)-TP).LE.PI) Go To 30
            T(J) = T(J) + PI + PI
   20     Continue
         call mexErrMsgTxt ( 'warning - tracking error ' )
           return
c        pause ' warning - tracking error '
c        Call CRASH(TTYOUT,4,0)
   30     TP = T(J)
          H1AX = DMAX1(H1AX,T(J)*PWR)
          H1IN = DMIN1(H1IN,T(J)*PWR)
   40   Continue
        PWR = DMIN1(PWR,1.98D0*PI*PWR/(H1AX-H1IN))
        Z(I) = DCMPLX(0.D0,0.D0)
        Do 50 J = 2, N
          Z(MOD(J+I-2,N)+1) = R(J) ** PWR * 
     *        CDEXP(DCMPLX(0.D0,T(J)*PWR))
   50   Continue
        ZD = 1. / (Z(N2)-Z(N1))
        Z0 = Z(N1)
        Do 60 J = 1, N
          Z(J) = (Z(J)-Z0) * ZD
   60   Continue
   70 Continue
      Return
      End
