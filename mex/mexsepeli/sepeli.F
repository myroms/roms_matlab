c svn $Id$

      subroutine sepeli_setup ( l2_in, m2_in, seta_in, sxi_in, 
     &                          x_in, y_in, 
     &                          sepx, sepy, ierr )


c     input ranges for independent variables.
c     0 < x < l2, 0 < y < y2
      integer l2, m2

      integer neta, nxi
      integer ierr

#include "mexsepeli.inc"

      real*8 sxi_in(0:NX2)
      real*8 seta_in(0:NY2)
      real*8 x_in(0:NX2*NY2) 
      real*8 y_in(0:NX2*NY2)
      real*8 sepx(0:NX2*NY2) 
      real*8 sepy(0:NX2*NY2)
      real*8 pertrb
c
c l2
      l2 = l2_in
c
c m2
      m2 = m2_in

      neta = m2+1
      nxi = l2+1
c

c
c     x
      k = 0
      do 10 j = 0, nxi - 1
         do 5 i = 0, neta - 1
            x(j,i) = x_in(k)
            k = k + 1 
    5    continue
   10 continue  
           

c     y
      k = 0
      do 20 j = 0, nxi - 1
         do 15 i = 0, neta - 1
            y(j,i) = y_in(k)
            k = k + 1 
   15    continue
   20 continue  
           

c
c     seta
      do 30 j = 0, neta - 1
         seta(j) = seta_in(j)
   30 continue  



c
c     sxi
      do 40 j = 0, nxi - 1
         sxi(j) = sxi_in(j)
   40 continue  

c
c     Construct right hand side
      do 50 j = 0, neta - 1
        do 45 i = 0, nxi - 1
            rhs(j,i) = 0.
   45   continue
   50 continue

c
c     Call the computational subroutine.
      ewrk(1) = nwrk
      call sepeli(0,2,0.d0,dble(l2),l2,1,wrk,wrk,wrk,wrk,0.d0,
     *    dble(m2),m2,
     *    1,wrk,wrk,wrk,wrk,rhs,x,nx2+1,ewrk,pertrb,ierr)
      if ( ierr .ne. 0 ) then
          return
      end if

      ewrk(1) = nwrk
      call sepeli(0,2,0.d0,dble(l2),l2,1,wrk,wrk,wrk,wrk,0.d0,
     *    dble(m2),m2,
     *    1,wrk,wrk,wrk,wrk,rhs,y,nx2+1,ewrk,pertrb,ierr)
      if ( ierr .eq. 1 ) then
          return
      end if
c

c
c Copy into output matrices.
      k = 0
      do 70 i = 0, nxi - 1
         do 65 j = 0, neta - 1
            sepx(k) = x(i,j)
            k = k + 1 
   65    continue
   70 continue  

      k = 0
      do 80 i = 0, nxi - 1
         do 75 j = 0, neta - 1
            sepy(k) = y(i,j)
            k = k + 1 
   75    continue
   80 continue  

      return
      end


c      Subroutine SEPELI(INTL,IORDER,A,B,M,MBDCND,BDA,ALPHA,BDB,BETA,C,D,
c     *    N,NBDCND,BDC,GAMA,BDD,XNU,COFX,COFY,GRHS,USOL,IDMN,W,PERTRB,
c     *    IERROR)
      SUBROUTINE sepeli(INTL,IORDER,A,B,M,MBDCND,BDA,ALPHA,BDB,BETA,C,D,
     *    N,NBDCND,BDC,GAMA,BDD,XNU,GRHS,USOL,IDMN,W,PERTRB,
     *    IERROR)
c
c dimension of           bda(n+1), bdb(n+1), bdc(m+1), bdd(m+1),
c arguments              usol(idmn,n+1),grhs(idmn,n+1),
c                        w (see argument list)
c
c latest revision        march 1985
c
c purpose                sepeli solves for either the second-order
c                        finite difference approximation or a
c                        fourth-order approximation to a separable
c                        elliptic equation
c
c                          2    2
c                          af(x)*d u/dx + bf(x)*du/dx  + cf(x)*u +
c                          2    2
c                          df(y)*d u/dy  + ef(y)*du/dy + ff(y)*u
c
c                          = g(x,y)
c
c                        on a rectangle (x greater than or equal to a
c                        and less than or equal to b; y greater than
c                        or equal to c and less than or equal to d).
c                        any combination of periodic or mixed boundary
c                        conditions is allowed.
c
c                        the possible boundary conditions are:
c                        in the x-direction:
c                        (0) periodic, u(x+b-a,y)=u(x,y) for all
c                            y,x (1) u(a,y), u(b,y) are specified for
c                            all y
c                        (2) u(a,y), du(b,y)/dx+beta*u(b,y) are
c                            specified for all y
c                        (3) du(a,y)/dx+alpha*u(a,y),du(b,y)/dx+
c                            beta*u(b,y) are specified for all y
c                        (4) du(a,y)/dx+alpha*u(a,y),u(b,y) are
c                            specified for all y
c
c                        in the y-direction:
c                        (0) periodic, u(x,y+d-c)=u(x,y) for all x,y
c                        (1) u(x,c),u(x,d) are specified for all x
c                        (2) u(x,c),du(x,d)/dy+xnu*u(x,d) are
c                            specified for all x
c                        (3) du(x,c)/dy+gama*u(x,c),du(x,d)/dy+
c                            xnu*u(x,d) are specified for all x
c                        (4) du(x,c)/dy+gama*u(x,c),u(x,d) are
c                            specified for all x
c
c usage                  call sepeli (intl,iorder,a,b,m,mbdcnd,bda,
c                                     alpha,bdb,beta,c,d,n,nbdcnd,bdc,
c                                     gama,bdd,xnu,cofx,cofy,grhs,usol,
c                                     idmn,w,pertrb,ierror)
c
c arguments
c on input               intl
c                          = 0 on initial entry to sepeli or if any
c                              of the arguments c,d, n, nbdcnd, cofy
c                              are changed from a previous call
c                          = 1 if c, d, n, nbdcnd, cofy are unchanged
c                              from the previous call.
c
c                        iorder
c                          = 2 if a second-order approximation
c                              is sought
c                          = 4 if a fourth-order approximation
c                              is sought
c
c                        a,b
c                          the range of the x-independent variable,
c                          i.e., x is greater than or equal to a
c                          and less than or equal to b.  a must be
c                          less than b.
c
c                        m
c                          the number of panels into which the
c                          interval [a,b] is subdivided. hence,
c                          there will be m+1 grid points in the x-
c                          direction given by xi=a+(i-1)*dlx
c                          for i=1,2,...,m+1 where dlx=(b-a)/m is
c                          the panel width.  m must be less than
c                          idmn and greater than 5.
c
c                        mbdcnd
c                          indicates the type of boundary condition
c                          at x=a and x=b
c
c                          = 0 if the solution is periodic in x, i.e.,
c                              u(x+b-a,y)=u(x,y)  for all y,x
c                          = 1 if the solution is specified at x=a
c                              and x=b, i.e., u(a,y) and u(b,y) are
c                              specified for all y
c                          = 2 if the solution is specified at x=a and
c                              the boundary condition is mixed at x=b,
c                              i.e., u(a,y) and du(b,y)/dx+beta*u(b,y)
c                              are specified for all y
c                          = 3 if the boundary conditions at x=a and
c                              x=b are mixed, i.e.,
c                              du(a,y)/dx+alpha*u(a,y) and
c                              du(b,y)/dx+beta*u(b,y) are specified
c                              for all y
c                          = 4 if the boundary condition at x=a is
c                              mixed and the solution is specified
c                              at x=b, i.e., du(a,y)/dx+alpha*u(a,y)
c                              and u(b,y) are specified for all y
c
c                        bda
c                          a one-dimensional array of length n+1
c                          that specifies the values of
c                          du(a,y)/dx+ alpha*u(a,y) at x=a, when
c                          mbdcnd=3 or 4.
c                          bda(j) = du(a,yj)/dx+alpha*u(a,yj),
c                          j=1,2,...,n+1. when mbdcnd has any other
c                          other value, bda is a dummy parameter.
c
c                        alpha
c                          the scalar multiplying the solution in
c                          case of a mixed boundary condition at x=a
c                          (see argument bda).  if mbdcnd is not
c                          equal to 3 or 4 then alpha is a dummy
c                          parameter.
c
c                        bdb
c                          a one-dimensional array of length n+1
c                          that specifies the values of
c                          du(b,y)/dx+ beta*u(b,y) at x=b.
c                          when mbdcnd=2 or 3
c                          bdb(j) = du(b,yj)/dx+beta*u(b,yj),
c                          j=1,2,...,n+1. when mbdcnd has any other
c                          other value, bdb is a dummy parameter.
c
c                        beta
c                          the scalar multiplying the solution in
c                          case of a mixed boundary condition at
c                          x=b (see argument bdb).  if mbdcnd is
c                          not equal to 2 or 3 then beta is a dummy
c                          parameter.
c
c                        c,d
c                          the range of the y-independent variable,
c                          i.e., y is greater than or equal to c
c                          and less than or equal to d.  c must be
c                          less than d.
c
c                        n
c                          the number of panels into which the
c                          interval [c,d] is subdivided.
c                          hence, there will be n+1 grid points
c                          in the y-direction given by
c                          yj=c+(j-1)*dly for j=1,2,...,n+1 where
c                          dly=(d-c)/n is the panel width.
c                          in addition, n must be greater than 4.
c
c                        nbdcnd
c                          indicates the types of boundary conditions
c                          at y=c and y=d
c
c                          = 0 if the solution is periodic in y,
c                              i.e., u(x,y+d-c)=u(x,y)  for all x,y
c                          = 1 if the solution is specified at y=c
c                              and y = d, i.e., u(x,c) and u(x,d)
c                              are specified for all x
c                          = 2 if the solution is specified at y=c
c                              and the boundary condition is mixed
c                              at y=d, i.e., u(x,c) and
c                              du(x,d)/dy+xnu*u(x,d) are specified
c                              for all x
c                          = 3 if the boundary conditions are mixed
c                              at y=c and y=d, i.e.,
c                              du(x,d)/dy+gama*u(x,c) and
c                              du(x,d)/dy+xnu*u(x,d) are specified
c                              for all x
c                          = 4 if the boundary condition is mixed
c                              at y=c and the solution is specified
c                              at y=d, i.e. du(x,c)/dy+gama*u(x,c)
c                              and u(x,d) are specified for all x
c
c                        bdc
c                          a one-dimensional array of length m+1
c                          that specifies the value of
c                          du(x,c)/dy+gama*u(x,c) at y=c.
c                          when nbdcnd=3 or 4 bdc(i) = du(xi,c)/dy +
c                          gama*u(xi,c), i=1,2,...,m+1.
c                          when nbdcnd has any other value, bdc
c                          is a dummy parameter.
c
c                        gama
c                          the scalar multiplying the solution in
c                          case of a mixed boundary condition at
c                          y=c (see argument bdc).  if nbdcnd is
c                          not equal to 3 or 4 then gama is a dummy
c                          parameter.
c
c                        bdd
c                          a one-dimensional array of length m+1
c                          that specifies the value of
c                          du(x,d)/dy + xnu*u(x,d) at y=c.
c                          when nbdcnd=2 or 3 bdd(i) = du(xi,d)/dy +
c                          xnu*u(xi,d), i=1,2,...,m+1.
c                          when nbdcnd has any other value, bdd
c                          is a dummy parameter.
c
c                        xnu
c                          the scalar multiplying the solution in
c                          case of a mixed boundary condition at
c                          y=d (see argument bdd).  if nbdcnd is
c                          not equal to 2 or 3 then xnu is a
c                          dummy parameter.
c
c                        cofx
c                          a user-supplied subprogram with
c                          parameters x, afun, bfun, cfun which
c                          returns the values of the x-dependent
c                          coefficients af(x), bf(x), cf(x) in the
c                          elliptic equation at x.
c
c                        cofy
c                          a user-supplied subprogram with parameters
c                          y, dfun, efun, ffun which returns the
c                          values of the y-dependent coefficients
c                          df(y), ef(y), ff(y) in the elliptic
c                          equation at y.
c
c                          note:  cofx and cofy must be declared
c                          external in the calling routine.
c                          the values returned in afun and dfun
c                          must satisfy afun*dfun greater than 0
c                          for a less than x less than b, c less
c                          than y less than d (see ierror=10).
c                          the coefficients provided may lead to a
c                          matrix equation which is not diagonally
c                          dominant in which case solution may fail
c                          (see ierror=4).
c
c                        grhs
c                          a two-dimensional array that specifies the
c                          values of the right-hand side of the
c                          elliptic equation, i.e.,
c                          grhs(i,j)=g(xi,yi), for i=2,...,m,
c                          j=2,...,n.  at the boundaries, grhs is
c                          defined by
c
c                          mbdcnd   grhs(1,j)   grhs(m+1,j)
c                          ------   ---------   -----------
c                            0      g(a,yj)     g(b,yj)
c                            1         *           *
c                            2         *        g(b,yj)  j=1,2,...,n+1
c                            3      g(a,yj)     g(b,yj)
c                            4      g(a,yj)        *
c
c                          nbdcnd   grhs(i,1)   grhs(i,n+1)
c                          ------   ---------   -----------
c                            0      g(xi,c)     g(xi,d)
c                            1         *           *
c                            2         *        g(xi,d)  i=1,2,...,m+1
c                            3      g(xi,c)     g(xi,d)
c                            4      g(xi,c)        *
c
c                          where * means these quantities are not used.
c                          grhs should be dimensioned idmn by at least
c                          n+1 in the calling routine.
c
c                        usol
c                          a two-dimensional array that specifies the
c                          values of the solution along the boundaries.
c                          at the boundaries, usol is defined by
c
c                          mbdcnd   usol(1,j)   usol(m+1,j)
c                          ------   ---------   -----------
c                            0         *           *
c                            1      u(a,yj)     u(b,yj)
c                            2      u(a,yj)        *     j=1,2,...,n+1
c                            3         *           *
c                            4         *        u(b,yj)
c
c                          nbdcnd   usol(i,1)   usol(i,n+1)
c                          ------   ---------   -----------
c                            0         *           *
c                            1      u(xi,c)     u(xi,d)
c                            2      u(xi,c)        *     i=1,2,...,m+1
c                            3         *           *
c                            4         *        u(xi,d)
c
c                          where * means the quantities are not used
c                          in the solution.
c
c                          if iorder=2, the user may equivalence grhs
c                          and usol to save space.  note that in this
c                          case the tables specifying the boundaries
c                          of the grhs and usol arrays determine the
c                          boundaries uniquely except at the corners.
c                          if the tables call for both g(x,y) and
c                          u(x,y) at a corner then the solution must
c                          be chosen.  for example, if mbdcnd=2 and
c                          nbdcnd=4, then u(a,c), u(a,d), u(b,d) must
c                          be chosen at the corners in addition
c                          to g(b,c).
c
c                          if iorder=4, then the two arrays, usol and
c                          grhs, must be distinct.
c
c                          usol should be dimensioned idmn by at least
c                          n+1 in the calling routine.
c
c                        idmn
c                          the row (or first) dimension of the arrays
c                          grhs and usol as it appears in the program
c                          calling sepeli.  this parameter is used
c                          to specify the variable dimension of grhs
c                          and usol.  idmn must be at least 7 and
c                          greater than or equal to m+1.
c
c                        w
c                          a one-dimensional array that must be
c                          provided by the user for work space.
c                          let k=int(log2(n+1))+1 and set l=2**(k+1).
c                          then (k-2)*l+k+10*n+12*m+27 will suffice
c                          as a length of w.  the actual length of w
c                          in the calling routine must be set in w(1)
c                          (see ierror=11).
c
c on output              usol
c                          contains the approximate solution to the
c                          elliptic equation.
c                          usol(i,j) is the approximation to u(xi,yj)
c                          for i=1,2...,m+1 and j=1,2,...,n+1.
c                          the approximation has error
c                          o(dlx**2+dly**2) if called with iorder=2
c                          and o(dlx**4+dly**4) if called with
c                          iorder=4.
c
c                        w
c                          contains intermediate values that must not
c                          be destroyed if sepeli is called again with
c                          intl=1.  in addition w(1) contains the
c                          exact minimal length (in floating point)
c                          required for the work space (see ierror=11).
c
c                        pertrb
c                          if a combination of periodic or derivative
c                          boundary conditions
c                          (i.e., alpha=beta=0 if mbdcnd=3;
c                          gama=xnu=0 if nbdcnd=3) is specified
c                          and if the coefficients of u(x,y) in the
c                          separable elliptic equation are zero
c                          (i.e., cf(x)=0 for x greater than or equal
c                          to a and less than or equal to b;
c                          ff(y)=0 for y greater than or equal to c
c                          and less than or equal to d) then a
c                          solution may not exist.  pertrb is a
c                          constant calculated and subtracted from
c                          the right-hand side of the matrix equations
c                          generated by sepeli which insures that a
c                          solution exists. sepeli then computes this
c                          solution which is a weighted minimal least
c                          squares solution to the original problem.
c
c                        ierror
c                          an error flag that indicates invalid input
c                          parameters or failure to find a solution
c                          = 0 no error
c                          = 1 if a greater than b or c greater than d
c                          = 2 if mbdcnd less than 0 or mbdcnd greater
c                              than 4
c                          = 3 if nbdcnd less than 0 or nbdcnd greater
c                              than 4
c                          = 4 if attempt to find a solution fails.
c                              (the linear system generated is not
c                              diagonally dominant.)
c                          = 5 if idmn is too small
c                              (see discussion of idmn)
c                          = 6 if m is too small or too large
c                              (see discussion of m)
c                          = 7 if n is too small (see discussion of n)
c                          = 8 if iorder is not 2 or 4
c                          = 9 if intl is not 0 or 1
c                          = 10 if afun*dfun less than or equal to 0
c                               for some interior mesh point (xi,yj)
c                          = 11 if the work space length input in w(1)
c                               is less than the exact minimal work
c                               space length required output in w(1).
c
c                          note (concerning ierror=4):  for the
c                          coefficients input through cofx, cofy,
c                          the discretization may lead to a block
c                          tridiagonal linear system which is not
c                          diagonally dominant (for example, this
c                          happens if cfun=0 and bfun/(2.*dlx) greater
c                          than afun/dlx**2).  in this case solution
c                          may fail.  this cannot happen in the limit
c                          as dlx, dly approach zero.  hence, the
c                          condition may be remedied by taking larger
c                          values for m or n.
c
c special conditions     see cofx, cofy argument descriptions above.
c
c i/o                    none
c
c precision              single
c
c required library       blktri, comf, q8qst4, and sepaux, which are
c files                  loaded by default on ncar's cray machines.
c
c language               fortran
c
c history                developed at ncar during 1975-76 by
c                        john c. adams of the scientific computing
c                        division.  released on ncar's public software
c                        libraries in january 1980.
c
c portability            fortran 66
c
c algorithm              sepeli automatically discretizes the
c                        separable elliptic equation which is then
c                        solved by a generalized cyclic reduction
c                        algorithm in the subroutine, blktri.  the
c                        fourth-order solution is obtained using
c                        'deferred corrections' which is described
c                        and referenced in sections, references and
c                        method.
c
c timing                 the operational count is proportional to
c                        m*n*log2(n).
c
c accuracy               the following accuracy results were obtained
c                        on a cdc 7600.  note that the fourth-order
c                        accuracy is not realized until the mesh is
c                        sufficiently refined.
c
c                                     second-order  fourth-order
c                            m    n     error         error
c
c                             6    6    6.8e-1        1.2e0
c                            14   14    1.4e-1        1.8e-1
c                            30   30    3.2e-2        9.7e-3
c                            62   62    7.5e-3        3.0e-4
c                           126  126    1.8e-3        3.5e-6
c
c portability            fortran 66
c
c references             keller, h.b., numerical methods for two-point
c                        boundary-value problems, blaisdel (1968),
c                        waltham, mass.
c
c                        swarztrauber, p., and r. sweet (1975):
c                        efficient fortran subprograms for the
c                        solution of elliptic partial differential
c                        equations.  ncar technical note
c                        ncar-tn/ia-109, pp. 135-137.
c***********************************************************************
      IMPLICIT NONE
      integer IDMN, INTL, IORDER, MBDCND, NBDCND, M, N, IERROR
      integer L, LOGB2N, LL, K, LENGTH, LINPUT
      integer LOUTPT
      integer I1, I2, I3, I4, I5, I6, I7, I8, I9, I10, I11, I12, I13
      Dimension GRHS(IDMN,1), USOL(IDMN,IDMN)
      Dimension BDA(1), BDB(1), BDC(1), BDD(1), W(1)
      real*8 BDA
      real*8 ALPHA
      real*8 BDB
      real*8 BETA
      real*8 BDC
      real*8 GAMA
      real*8 BDD
      real*8 XNU
      real*8 GRHS
      real*8 USOL
      real*8 W
      real*8 A, B, C, D, PERTRB
c      External COFX, COFY
      Logical Q8Q4
      Save Q8Q4
      Data Q8Q4 /.TRUE./
c
      If (Q8Q4) Then
c         call q8qst4('loclib','sepeli','sepeli','version 01')
        Q8Q4 = .FALSE.
      End If
c
c     check input parameters
c
c      Call CHKPRM(INTL,IORDER,A,B,M,MBDCND,C,D,N,NBDCND,COFX,COFY,IDMN,
c     *    IERROR)
      Call CHKPRM(INTL,IORDER,A,B,M,MBDCND,C,D,N,NBDCND,IDMN,
     *    IERROR)
      If (IERROR.NE.0) Return
c
c     compute minimum work space and check work space length input
c
      L = N + 1
      If (NBDCND.EQ.0) L = N
      LOGB2N = INT(DLOG(DBLE(L)+0.d5)/DLOG(2.d0)) + 1
      LL = 2 ** (LOGB2N+1)
      K = M + 1
      L = N + 1
      LENGTH = (LOGB2N-2) * LL + LOGB2N + MAX0(2*L,6*K) + 5
      If (NBDCND.EQ.0) LENGTH = LENGTH + 2 * L
      IERROR = 11
      LINPUT = INT(W(1)+0.5)
      LOUTPT = LENGTH + 6 * (K+L) + 1
      W(1) = DBLE(LOUTPT)
      If (LOUTPT.GT.LINPUT) Return
      IERROR = 0
c
c     set work space indices
c
      I1 = LENGTH + 2
      I2 = I1 + L
      I3 = I2 + L
      I4 = I3 + L
      I5 = I4 + L
      I6 = I5 + L
      I7 = I6 + L
      I8 = I7 + K
      I9 = I8 + K
      I10 = I9 + K
      I11 = I10 + K
      I12 = I11 + K
      I13 = 2
c      Call SPELIP(INTL,IORDER,A,B,M,MBDCND,BDA,ALPHA,BDB,BETA,C,D,N,
c     *    NBDCND,BDC,GAMA,BDD,XNU,COFX,COFY,W(I1),W(I2),W(I3),W(I4),
c     *    W(I5),W(I6),W(I7),W(I8),W(I9),W(I10),W(I11),W(I12),GRHS,USOL,
c     *    IDMN,W(I13),PERTRB,IERROR)
      Call SPELIP(INTL,IORDER,A,B,M,MBDCND,BDA,ALPHA,BDB,BETA,C,D,N,
     *    NBDCND,BDC,GAMA,BDD,XNU,W(I1),W(I2),W(I3),W(I4),
     *    W(I5),W(I6),W(I7),W(I8),W(I9),W(I10),W(I11),W(I12),GRHS,USOL,
     *    IDMN,W(I13),PERTRB,IERROR)
      Return
      End













c      Subroutine SPELIP(INTL,IORDER,A,B,M,MBDCND,BDA,ALPHA,BDB,BETA,C,D,
c     *    N,NBDCND,BDC,GAMA,BDD,XNU,COFX,COFY,AN,BN,CN,DN,UN,ZN,AM,BM,
c     *    CM,DM,UM,ZM,GRHS,USOL,IDMN,W,PERTRB,IERROR)
      Subroutine SPELIP(INTL,IORDER,A,B,M,MBDCND,BDA,ALPHA,BDB,BETA,C,D,
     *    N,NBDCND,BDC,GAMA,BDD,XNU,AN,BN,CN,DN,UN,ZN,AM,BM,
     *    CM,DM,UM,ZM,GRHS,USOL,IDMN,W,PERTRB,IERROR)
c
c     spelip sets up vectors and arrays for input to blktri
c     and computes a second order solution in usol.  a return jump to
c     sepeli occurrs if iorder=2.  if iorder=4 a fourth order
c     solution is generated in usol.
c
      IMPLICIT NONE
      integer IDMN
      Dimension BDA(1), BDB(1), BDC(1), BDD(1), W(1)
      Dimension GRHS(IDMN,1), USOL(IDMN,IDMN)
      Dimension AN(1), BN(1), CN(1), DN(1), UN(1), ZN(1)
      Dimension AM(1), BM(1), CM(1), DM(1), UM(1), ZM(1)

      integer KSWX, KSWY, K, MBDCND, NBDCND 

      integer AIT, BIT, CIT, DIT, IS, MIT, NIT, L
      real*8 DLX, DLY, TDLX3, TDLY3, DLX4, DLY4
      integer MS, JS, NS
      Common /SPLP/ KSWX, KSWY, K, L, 
     *    DLX, DLY, TDLX3, TDLY3, DLX4, DLY4,
     *    AIT, BIT, CIT, DIT, MIT, NIT, IS, MS, JS, NS

      Logical SINGLR

      real*8 BDA, BDB, BDC, BDD
      real*8 ALPHA, BETA, GAMA
      real*8 A, B, C, D, W
      real*8 XNU
      real*8 GRHS, USOL
      real*8 AN, BN, CN, DN, UN, ZN
      real*8 AM, BM, CM, DM, UM, ZM
      real*8 YJ, DJ, EJ, FJ
      real*8 XI, AI, BI, CI
      real*8 PRTRB, PERTRB

      integer M, N, I, J, I1, MP, NP, IORDER, IORD, INTL, IERROR
      real*8 AXI, BXI, CXI, DYJ, EYJ, FYJ, AX1, CXM, DY1, FYN

c      External COFX, COFY
c
c     set parameters internally
c
      KSWX = MBDCND + 1
      KSWY = NBDCND + 1
      K = M + 1
      L = N + 1
      AIT = A
      BIT = B
      CIT = C
      DIT = D
c
c     set right hand side values from grhs in usol on the interior
c     and non-specified boundaries.
c
      Do 20 I = 2, M
        Do 10 J = 2, N
          USOL(I,J) = GRHS(I,J)
   10   Continue
   20 Continue
      If (KSWX.NE.2.AND.KSWX.NE.3) Then
        Do 30 J = 2, N
          USOL(1,J) = GRHS(1,J)
   30   Continue
      End If
      If (KSWX.NE.2.AND.KSWX.NE.5) Then
        Do 40 J = 2, N
          USOL(K,J) = GRHS(K,J)
   40   Continue
      End If
      If (KSWY.NE.2.AND.KSWY.NE.3) Then
        Do 50 I = 2, M
          USOL(I,1) = GRHS(I,1)
   50   Continue
      End If
      If (KSWY.NE.2.AND.KSWY.NE.5) Then
        Do 60 I = 2, M
          USOL(I,L) = GRHS(I,L)
   60   Continue
      End If
      If (KSWX.NE.2.AND.KSWX.NE.3.AND.KSWY.NE.2.AND.KSWY.NE.3) 
     *    USOL(1,1) = GRHS(1,1)
      If (KSWX.NE.2.AND.KSWX.NE.5.AND.KSWY.NE.2.AND.KSWY.NE.3) 
     *    USOL(K,1) = GRHS(K,1)
      If (KSWX.NE.2.AND.KSWX.NE.3.AND.KSWY.NE.2.AND.KSWY.NE.5) 
     *    USOL(1,L) = GRHS(1,L)
      If (KSWX.NE.2.AND.KSWX.NE.5.AND.KSWY.NE.2.AND.KSWY.NE.5) 
     *    USOL(K,L) = GRHS(K,L)
      I1 = 1
c
c     set switches for periodic or non-periodic boundaries
c
      MP = 1
      NP = 1
      If (KSWX.EQ.1) MP = 0
      If (KSWY.EQ.1) NP = 0
c
c     set dlx,dly and size of block tri-diagonal system generated
c     in nint,mint
c
      DLX = (BIT-AIT) / DBLE(M)
      MIT = K - 1
      If (KSWX.EQ.2) MIT = K - 2
      If (KSWX.EQ.4) MIT = K
      DLY = (DIT-CIT) / DBLE(N)
      NIT = L - 1
      If (KSWY.EQ.2) NIT = L - 2
      If (KSWY.EQ.4) NIT = L
      TDLX3 = 2.0 * DLX ** 3
      DLX4 = DLX ** 4
      TDLY3 = 2.0 * DLY ** 3
      DLY4 = DLY ** 4
c
c     set subscript limits for portion of array to input to blktri
c
      IS = 1
      JS = 1
      If (KSWX.EQ.2.OR.KSWX.EQ.3) IS = 2
      If (KSWY.EQ.2.OR.KSWY.EQ.3) JS = 2
      NS = NIT + JS - 1
      MS = MIT + IS - 1
c
c     set x - direction
c
      Do 70 I = 1, MIT
        XI = AIT + DBLE(IS+I-2) * DLX
        Call COFX(XI,AI,BI,CI)
        AXI = (AI/DLX-0.5*BI) / DLX
        BXI = -2. * AI / DLX ** 2 + CI
        CXI = (AI/DLX+0.5*BI) / DLX
        AM(I) = AXI
        BM(I) = BXI
        CM(I) = CXI
   70 Continue
c
c     set y direction
c
      Do 80 J = 1, NIT
        YJ = CIT + DBLE(JS+J-2) * DLY
        Call COFY(YJ,DJ,EJ,FJ)
        DYJ = (DJ/DLY-0.5*EJ) / DLY
        EYJ = (-2.*DJ/DLY**2+FJ)
        FYJ = (DJ/DLY+0.5*EJ) / DLY
        AN(J) = DYJ
        BN(J) = EYJ
        CN(J) = FYJ
   80 Continue
c
c     adjust edges in x direction unless periodic
c
      AX1 = AM(1)
      CXM = CM(MIT)
      Go To (130,90,110,120,100), KSWX
c
c     dirichlet-dirichlet in x direction
c
   90 AM(1) = 0.0
      CM(MIT) = 0.0
      Go To 130
c
c     mixed-dirichlet in x direction
c
  100 AM(1) = 0.0
      BM(1) = BM(1) + 2. * ALPHA * DLX * AX1
      CM(1) = CM(1) + AX1
      CM(MIT) = 0.0
      Go To 130
c
c     dirichlet-mixed in x direction
c
  110 AM(1) = 0.0
      AM(MIT) = AM(MIT) + CXM
      BM(MIT) = BM(MIT) - 2. * BETA * DLX * CXM
      CM(MIT) = 0.0
      Go To 130
c
c     mixed - mixed in x direction
c
  120 Continue
      AM(1) = 0.0
      BM(1) = BM(1) + 2. * DLX * ALPHA * AX1
      CM(1) = CM(1) + AX1
      AM(MIT) = AM(MIT) + CXM
      BM(MIT) = BM(MIT) - 2. * DLX * BETA * CXM
      CM(MIT) = 0.0
  130 Continue
c
c     adjust in y direction unless periodic
c
      DY1 = AN(1)
      FYN = CN(NIT)
      Go To (180,140,160,170,150), KSWY
c
c     dirichlet-dirichlet in y direction
c
  140 Continue
      AN(1) = 0.0
      CN(NIT) = 0.0
      Go To 180
c
c     mixed-dirichlet in y direction
c
  150 Continue
      AN(1) = 0.0
      BN(1) = BN(1) + 2. * DLY * GAMA * DY1
      CN(1) = CN(1) + DY1
      CN(NIT) = 0.0
      Go To 180
c
c     dirichlet-mixed in y direction
c
  160 AN(1) = 0.0
      AN(NIT) = AN(NIT) + FYN
      BN(NIT) = BN(NIT) - 2. * DLY * XNU * FYN
      CN(NIT) = 0.0
      Go To 180
c
c     mixed - mixed direction in y direction
c
  170 Continue
      AN(1) = 0.0
      BN(1) = BN(1) + 2. * DLY * GAMA * DY1
      CN(1) = CN(1) + DY1
      AN(NIT) = AN(NIT) + FYN
      BN(NIT) = BN(NIT) - 2.0 * DLY * XNU * FYN
      CN(NIT) = 0.0
  180 If (KSWX.NE.1) Then
c
c     adjust usol along x edge
c
        Do 190 J = JS, NS
          If (KSWX.EQ.2.OR.KSWX.EQ.3) Then
            USOL(IS,J) = USOL(IS,J) - AX1 * USOL(1,J)
          Else
            USOL(IS,J) = USOL(IS,J) + 2.0 * DLX * AX1 * BDA(J)
          End If
          If (KSWX.EQ.2.OR.KSWX.EQ.5) Then
            USOL(MS,J) = USOL(MS,J) - CXM * USOL(K,J)
          Else
            USOL(MS,J) = USOL(MS,J) - 2.0 * DLX * CXM * BDB(J)
          End If
  190   Continue
      End If
      If (KSWY.NE.1) Then
c
c     adjust usol along y edge
c
        Do 200 I = IS, MS
          If (KSWY.EQ.2.OR.KSWY.EQ.3) Then
            USOL(I,JS) = USOL(I,JS) - DY1 * USOL(I,1)
          Else
            USOL(I,JS) = USOL(I,JS) + 2.0 * DLY * DY1 * BDC(I)
          End If
          If (KSWY.EQ.2.OR.KSWY.EQ.5) Then
            USOL(I,NS) = USOL(I,NS) - FYN * USOL(I,L)
          Else
            USOL(I,NS) = USOL(I,NS) - 2.0 * DLY * FYN * BDD(I)
          End If
  200   Continue
      End If
c
c     save adjusted edges in grhs if iorder=4
c
      If (IORDER.EQ.4) Then
        Do 210 J = JS, NS
          GRHS(IS,J) = USOL(IS,J)
          GRHS(MS,J) = USOL(MS,J)
  210   Continue
        Do 220 I = IS, MS
          GRHS(I,JS) = USOL(I,JS)
          GRHS(I,NS) = USOL(I,NS)
  220   Continue
      End If
      IORD = IORDER
      PERTRB = 0.0
c
c     check if operator is singular
c
c      Call CHKSNG(MBDCND,NBDCND,ALPHA,BETA,GAMA,XNU,COFX,COFY,SINGLR)
      Call CHKSNG(MBDCND,NBDCND,ALPHA,BETA,GAMA,XNU,SINGLR)
c
c     compute non-zero eigenvector in null space of transpose
c     if singular
c
      If (SINGLR) Call SEPTRI(MIT,AM,BM,CM,DM,UM,ZM)
      If (SINGLR) Call SEPTRI(NIT,AN,BN,CN,DN,UN,ZN)
c
c     make initialization call to blktri
c
      If (INTL.EQ.0) Call BLKTRI(INTL,NP,NIT,AN,BN,CN,MP,MIT,AM,BM,CM,
     *    IDMN,USOL(IS,JS),IERROR,W)
      If (IERROR.NE.0) Return
c
c     adjust right hand side if necessary
c
  230 Continue
      If (SINGLR) Call SEPORT(USOL,IDMN,ZN,ZM,PERTRB)
c
c     compute solution
c
      Call BLKTRI(I1,NP,NIT,AN,BN,CN,MP,MIT,AM,BM,CM,IDMN,USOL(IS,JS),
     *    IERROR,W)
      If (IERROR.NE.0) Return
c
c     set periodic boundaries if necessary
c
      If (KSWX.EQ.1) Then
        Do 240 J = 1, L
          USOL(K,J) = USOL(1,J)
  240   Continue
      End If
      If (KSWY.EQ.1) Then
        Do 250 I = 1, K
          USOL(I,L) = USOL(I,1)
  250   Continue
      End If
c
c     minimize solution with respect to weighted least squares
c     norm if operator is singular
c
      If (SINGLR) Call SEPMIN(USOL,IDMN,ZN,ZM,PRTRB)
c
c     return if deferred corrections and a fourth order solution are
c     not flagged
c
      If (IORD.EQ.2) Return
      IORD = 2
c
c     compute new right hand side for fourth order solution
c
c      Call DEFER(COFX,COFY,IDMN,USOL,GRHS)
      Call DEFER(IDMN,USOL,GRHS)
      Go To 230
      End



















c      Subroutine CHKPRM(INTL,IORDER,A,B,M,MBDCND,C,D,N,NBDCND,COFX,COFY,
c     *    IDMN,IERROR)
      Subroutine CHKPRM(INTL,IORDER,A,B,M,MBDCND,C,D,N,NBDCND,
     *    IDMN,IERROR)
c
c     this program checks the input parameters for errors
c
c      External COFX, COFY
c
c     check definition of solution region
c
      IMPLICIT NONE
      integer IERROR, MBDCND, NBDCND, IDMN, M, N, IORDER, INTL
      real*8 A, B, C, D

      real*8 YJ, DJ, EJ, FJ
      real*8 XI, AI, BI, CI, DLX, DLY
      integer I, J




      IERROR = 1
      If (A.GE.B.OR.C.GE.D) Return
c
c     check boundary switches
c
      IERROR = 2
      If (MBDCND.LT.0.OR.MBDCND.GT.4) Return
      IERROR = 3
      If (NBDCND.LT.0.OR.NBDCND.GT.4) Return
c
c     check first dimension in calling routine
c
      IERROR = 5
      If (IDMN.LT.7) Return
c
c     check m
c
      IERROR = 6
      If (M.GT.(IDMN-1).OR.M.LT.6) Return
c
c     check n
c
      IERROR = 7
      If (N.LT.5) Return
c
c     check iorder
c
      IERROR = 8
      If (IORDER.NE.2.AND.IORDER.NE.4) Return
c
c     check intl
c
      IERROR = 9
      If (INTL.NE.0.AND.INTL.NE.1) Return
c
c     check that equation is elliptic
c
      DLX = (B-A) / DBLE(M)
      DLY = (D-C) / DBLE(N)
      Do 20 I = 2, M
        XI = A + DBLE(I-1) * DLX
        Call COFX(XI,AI,BI,CI)
        Do 10 J = 2, N
          YJ = C + DBLE(J-1) * DLY
          Call COFY(YJ,DJ,EJ,FJ)
          If (AI*DJ.LE.0.0) Then
            IERROR = 10
            Return
          End If
   10   Continue
   20 Continue
c
c     no error found
c
      IERROR = 0
      Return
      End















c      Subroutine CHKSNG(MBDCND,NBDCND,ALPHA,BETA,GAMA,XNU,COFX,COFY,
c     *    SINGLR)
      Subroutine CHKSNG(MBDCND,NBDCND,ALPHA,BETA,GAMA,XNU,
     *    SINGLR)
c
c     this subroutine checks if the pde   sepeli
c     must solve is a singular operator
c
      IMPLICIT NONE

      integer MBDCND, NBDCND

      integer KSWX, KSWY, K, L
      integer AIT, BIT, CIT, DIT, IS, MIT, NIT
      real*8 DLX, DLY, TDLX3, TDLY3, DLX4, DLY4
      integer MS, JS, NS
      Common /SPLP/ KSWX, KSWY, K, L, 
     *    DLX, DLY, TDLX3, TDLY3, DLX4, DLY4,
     *    AIT, BIT, CIT, DIT, MIT, NIT, IS, MS, JS, NS


      integer I, J

      Logical SINGLR

      real*8 ALPHA
      real*8 BETA
      real*8 GAMA
      real*8 XNU



      real*8 YJ, DJ, EJ, FJ
      real*8 XI, AI, BI, CI

      SINGLR = .FALSE.
c
c     check if the boundary conditions are
c     entirely periodic and/or mixed
c
      If ((MBDCND.NE.0.AND.MBDCND.NE.3).OR.(NBDCND.NE.0.AND.NBDCND.NE.3)
     *    ) Return
c
c     check that mixed conditions are pure neuman
c
      If (MBDCND.EQ.3) Then
        If (ALPHA.NE.0.0.OR.BETA.NE.0.0) Return
      End If
      If (NBDCND.EQ.3) Then
        If (GAMA.NE.0.0.OR.XNU.NE.0.0) Return
      End If
c
c     check that non-derivative coefficient functions
c     are zero
c
      Do 10 I = IS, MS
        XI = AIT + DBLE(I-1) * DLX
        Call COFX(XI,AI,BI,CI)
        If (CI.NE.0.0) Return
   10 Continue
      Do 20 J = JS, NS
        YJ = CIT + DBLE(J-1) * DLY
        Call COFY(YJ,DJ,EJ,FJ)
        If (FJ.NE.0.0) Return
   20 Continue
c
c     the operator must be singular if this point is reached
c
      SINGLR = .TRUE.
      Return
      End


















c      Subroutine DEFER(COFX,COFY,IDMN,USOL,GRHS)
      Subroutine DEFER(IDMN,USOL,GRHS)
c
c     this subroutine first approximates the truncation error given by
c     trun1(x,y)=dlx**2*tx+dly**2*ty where
c     tx=afun(x)*uxxxx/12.0+bfun(x)*uxxx/6.0 on the interior and
c     at the boundaries if periodic(here uxxx,uxxxx are the third
c     and fourth partial derivatives of u with respect to x).
c     tx is of the form afun(x)/3.0*(uxxxx/4.0+uxxx/dlx)
c     at x=a or x=b if the boundary condition there is mixed.
c     tx=0.0 along specified boundaries.  ty has symmetric form
c     in y with x,afun(x),bfun(x) replaced by y,dfun(y),efun(y).
c     the second order solution in usol is used to approximate
c     (via second order finite differencing) the truncation error
c     and the result is added to the right hand side in grhs
c     and then transferred to usol to be used as a new right
c     hand side when calling blktri for a fourth order solution.
c
      IMPLICIT NONE
      integer KSWX, KSWY, K 
      integer AIT, BIT, CIT, DIT, IS, MIT, NIT

      integer L
      real*8 DLX, DLY, TDLX3, TDLY3, DLX4, DLY4
      integer MS, JS, NS
      Common /SPLP/ KSWX, KSWY, K, L, 
     *    DLX, DLY, TDLX3, TDLY3, DLX4, DLY4,
     *    AIT, BIT, CIT, DIT, MIT, NIT, IS, MS, JS, NS


      integer IDMN
      Dimension GRHS(IDMN,1), USOL(IDMN,1)
      real*8 GRHS, USOL, UXXX, UXXXX, UYYY, UYYYY
c      External COFX, COFY

      real*8 YJ, DJ, EJ, FJ
      real*8 XI, AI, BI, CI, TX, TY

      integer I, J
c
c     compute truncation error approximation over the entire mesh
c
      Do 20 J = JS, NS
        YJ = CIT + DBLE(J-1) * DLY
        Call COFY(YJ,DJ,EJ,FJ)
        Do 10 I = IS, MS
          XI = AIT + DBLE(I-1) * DLX
          Call COFX(XI,AI,BI,CI)
c
c     compute partial derivative approximations at (xi,yj)
c
          Call SEPDX(USOL,IDMN,I,J,UXXX,UXXXX)
          Call SEPDY(USOL,IDMN,I,J,UYYY,UYYYY)
          TX = AI * UXXXX / 12.0 + BI * UXXX / 6.0
          TY = DJ * UYYYY / 12.0 + EJ * UYYY / 6.0
c
c     reset form of truncation if at boundary which is non-periodic
c
          If (.NOT.(KSWX.EQ.1.OR.(I.GT.1.AND.I.LT.K))) Then
            TX = AI / 3.0 * (UXXXX/4.0+UXXX/DLX)
          End If
          If (.NOT.(KSWY.EQ.1.OR.(J.GT.1.AND.J.LT.L))) Then
            TY = DJ / 3.0 * (UYYYY/4.0+UYYY/DLY)
          End If
          GRHS(I,J) = GRHS(I,J) + DLX ** 2 * TX + DLY ** 2 * TY
   10   Continue
   20 Continue
c
c     reset the right hand side in usol
c
      Do 40 I = IS, MS
        Do 30 J = JS, NS
          USOL(I,J) = GRHS(I,J)
   30   Continue
   40 Continue
      Return
c
c revision history---
c
c december 1979    first added to nssl
c-----------------------------------------------------------------------
      End











      Subroutine BLKTRI(IFLG,NP,N,AN,BN,CN,MP,M,AM,BM,CM,IDIMY,Y,IERROR,
     *    W)
c
c dimension of           an(n),bn(n),cn(n),am(m),bm(m),cm(m),y(idimy,n),
c arguments              w(see argument list)
c
c latest revision        january 1985
c
c usage                  call blktri (iflg,np,n,an,bn,cn,mp,m,am,bm,
c                                     cm,idimy,y,ierror,w)
c
c purpose                blktri solves a system of linear equations
c                        of the form
c
c                        an(j)*x(i,j-1) + am(i)*x(i-1,j) +
c                        (bn(j)+bm(i))*x(i,j) + cn(j)*x(i,j+1) +
c                        cm(i)*x(i+1,j) = y(i,j)
c
c                        for i = 1,2,...,m  and  j = 1,2,...,n.
c
c                        i+1 and i-1 are evaluated modulo m and
c                        j+1 and j-1 modulo n, i.e.,
c
c                        x(i,0) = x(i,n),  x(i,n+1) = x(i,1),
c                        x(0,j) = x(m,j),  x(m+1,j) = x(1,j).
c
c                        these equations usually result from the
c                        discretization of separable elliptic
c                        equations.  boundary conditions may be
c                        dirichlet, neumann, or periodic.
c
c arguments
c
c on input               iflg
c
c                          = 0  initialization only.
c                               certain quantities that depend on np,
c                               n, an, bn, and cn are computed and
c                               stored in the work array w.
c
c                          = 1  the quantities that were computed
c                               in the initialization are used
c                               to obtain the solution x(i,j).
c
c                               note:
c                               a call with iflg=0 takes
c                               approximately one half the time
c                               as a call with iflg = 1.
c                               however, the initialization does
c                               not have to be repeated unless np,
c                               n, an, bn, or cn change.
c
c                        np
c                          = 0  if an(1) and cn(n) are not zero,
c                               which corresponds to periodic
c                               bounary conditions.
c
c                          = 1  if an(1) and cn(n) are zero.
c
c                        n
c                          the number of unknowns in the j-direction.
c                          n must be greater than 4.
c                          the operation count is proportional to
c                          mnlog2(n), hence n should be selected
c                          less than or equal to m.
c
c                        an,bn,cn
c                          one-dimensional arrays of length n
c                          that specify the coefficients in the
c                          linear equations given above.
c
c                        mp
c                          = 0  if am(1) and cm(m) are not zero,
c                               which corresponds to periodic
c                               boundary conditions.
c
c                          = 1  if am(1) = cm(m) = 0  .
c
c                        m
c                          the number of unknowns in the i-direction.
c                           m must be greater than 4.
c
c                        am,bm,cm
c                          one-dimensional arrays of length m that
c                          specify the coefficients in the linear
c                          equations given above.
c
c                        idimy
c                          the row (or first) dimension of the
c                          two-dimensional array y as it appears
c                          in the program calling blktri.
c                          this parameter is used to specify the
c                          variable dimension of y.
c                          idimy must be at least m.
c
c                        y
c                          a two-dimensional array that specifies
c                          the values of the right side of the linear
c                          system of equations given above.
c                          y must be dimensioned at least m*n.
c
c                        w
c                          a one-dimensional array that must be
c                          provided by the user for work space.
c                          if np=1 define k=int(log2(n))+1 and
c                          set l=2**(k+1) then w must have dimension
c                          (k-2)*l+k+5+max(2n,6m)
c
c                          if np=0 define k=int(log2(n-1))+1 and
c                          set l=2**(k+1) then w must have dimension
c                          (k-2)*l+k+5+2n+max(2n,6m)
c
c                          **important**
c                          for purposes of checking, the required
c                          dimension of w is computed by blktri and
c                          stored in w(1) in floating point format.
c
c arguments
c
c on output              y
c                          contains the solution x.
c
c                        ierror
c                          an error flag that indicates invalid
c                          input parameters.  except for number zer0,
c                          a solution is not attempted.
c
c                        = 0  no error.
c                        = 1  m is less than 5
c                        = 2  n is less than 5
c                        = 3  idimy is less than m.
c                        = 4  blktri failed while computing results
c                             that depend on the coefficient arrays
c                             an, bn, cn.  check these arrays.
c                        = 5  an(j)*cn(j-1) is less than 0 for some j.
c
c                             possible reasons for this condition are
c                             1. the arrays an and cn are not correct
c                             2. too large a grid spacing was used
c                                in the discretization of the elliptic
c                                equation.
c                             3. the linear equations resulted from a
c                                partial differential equation which
c                                was not elliptic.
c
c                        w
c                           contains intermediate values that must
c                           not be destroyed if blktri will be called
c                           again with iflg=1. w(1) contains the
c                           number of locations required by w in
c                           floating point format.
c
c
c special conditions     the algorithm may fail if abs(bm(i)+bn(j))
c                        is less than abs(am(i))+abs(an(j))+
c                        abs(cm(i))+abs(cn(j))
c                        for some i and j. the algorithm will also
c                        fail if an(j)*cn(j-1) is less than zero for
c                        some j.
c                        see the description of the output parameter
c                        ierror.
c
c i/o                    none
c
c precision              single
c
c required library       comf and q8qst4, which are automatically loaded
c files                  ncar's cray machines.
c
c language               fortran
c
c history                written by paul swarztrauber at ncar in the
c                        early 1970's.  rewritten and released in
c                        january, 1980.
c
c algorithm              generalized cyclic reduction
c
c portability            fortran 66.  approximate machine accuracy
c                        is computed in function epmach.
c
c references             swarztrauber,p. and r. sweet, 'efficient
c                        fortran subprograms for the solution of
c                        elliptic equations'
c                        ncar tn/ia-109, july, 1975, 138 pp.
c
c                        swarztrauber p. n.,a direct method for
c                        the discrete solution of separable
c                        elliptic equations, s.i.a.m.
c                        j. numer. anal.,11(1974) pp. 1136-1150.
c***********************************************************************

      IMPLICIT NONE

      integer IDIMY
      Dimension AN(1), BN(1), CN(1), AM(1), BM(1), CM(1), Y(IDIMY,1) 
      Dimension W(2)
      External PROD, PRODP, CPROD, CPRODP

      integer M, N, NP, MP
      integer K, NM, NCMPLX, IK, NL
      real*8 EPS, NPP, CNV
      Common /CBLKT/ K, NM, NCMPLX, IK, NPP, CNV, EPS

      real*8 AN, BN, CN, AM, BM, CM, Y, W



      integer IERROR, NH, IWAH, IW1, IWBH, IW2, IW3, IWD, IWW, IWU
      integer IFLG

c
c the following call is for monitoring library use at ncar
c
      Logical Q8Q4
      Save Q8Q4
      Data Q8Q4 /.TRUE./
      If (Q8Q4) Then
c         call q8qst4('loclib','blktri','blktri','version 01')
        Q8Q4 = .FALSE.
      End If
c
c test m and n for the proper form
c
      NM = N
      IERROR = 0
      If (M.LT.5) Then
        IERROR = 1
      Else
        If (NM.LT.3) Then
          IERROR = 2
        Else
          If (IDIMY.LT.M) Then
            IERROR = 3
          Else
            NH = N
            NPP = NP
            If (NPP.NE.0) Then
              NH = NH + 1
            End If
            IK = 2
            K = 1
   10       IK = IK + IK
            K = K + 1
            If (NH.GT.IK) Go To 10
            NL = IK
            IK = IK + IK
            NL = NL - 1
            IWAH = (K-2) * IK + K + 6
            If (NPP.NE.0) Then
c
c     divide w into working sub arrays
c
              IW1 = IWAH
              IWBH = IW1 + NM
              W(1) = DBLE(IW1-1+DMAX1(2.d0*NM,6.d0*M))
            Else
              IWBH = IWAH + NM + NM
              IW1 = IWBH
              W(1) = DBLE(IW1-1+DMAX1(2.d0*NM,6.d0*M))
              NM = NM - 1
            End If
c
c subroutine comp b computes the roots of the b polynomials
c
            If (IERROR.EQ.0) Then
              IW2 = IW1 + M
              IW3 = IW2 + M
              IWD = IW3 + M
              IWW = IWD + M
              IWU = IWW + M
              If (IFLG.EQ.0) Then
                Call COMPB(NL,IERROR,AN,BN,CN,W(2),W(IWAH),W(IWBH))
              Else
                If (MP.NE.0) Then
c
c subroutine blktr1 solves the linear system
c
                  Call BLKTR1(NL,AN,BN,CN,M,AM,BM,CM,IDIMY,Y,W(2),
     *                W(IW1),W(IW2),W(IW3),W(IWD),W(IWW),W(IWU),1)
c                  Call BLKTR1(NL,AN,BN,CN,M,AM,BM,CM,IDIMY,Y,W(2),
c     *                W(IW1),W(IW2),W(IW3),W(IWD),W(IWW),W(IWU),PROD,
c     *                CPROD)
                Else
c                  Call BLKTR1(NL,AN,BN,CN,M,AM,BM,CM,IDIMY,Y,W(2),
c     *                W(IW1),W(IW2),W(IW3),W(IWD),W(IWW),W(IWU),PRODP,
c     *                CPRODP)
                  Call BLKTR1(NL,AN,BN,CN,M,AM,BM,CM,IDIMY,Y,W(2),
     *                W(IW1),W(IW2),W(IW3),W(IWD),W(IWW),W(IWU),2)
                End If
              End If
            End If
          End If
        End If
      End If
      Return
      End









      Subroutine BLKTR1(N,AN,BN,CN,M,AM,BM,CM,IDIMY,Y,B,W1,W2,W3,WD,WW,
c         WU,PRDCT,CPRDCT)
     *    WU, PFLAG )
c
c blktr1 solves the linear system
c
c b  contains the roots of all the b polynomials
c w1,w2,w3,wd,ww,wu  are all working arrays
c prdct  is either prodp or prod depending on whether the boundary
c conditions in the m direction are periodic or not
c cprdct is either cprodp or cprod which are the complex versions
c of prodp and prod. these are called in the event that some
c of the roots of the b sub p polynomial are complex
c
c If PFLAG==1, then use PROD and CPROD as callable functions. 
c Otherwise use PRODP and CPRODP
c

      IMPLICIT NONE


      integer IDIMY
      Dimension AN(1), BN(1), CN(1), AM(1), BM(1), CM(1), B(1), W1(1), 
     *    W2(1), W3(1), WD(1), WW(1), WU(1), Y(IDIMY,1)
      real*8 AN, BN, CN, AM, BM, CM, B, W1
      real*8 W2, W3, WD, WW, WU, Y
      real*8 DUM
      integer M, N
      integer PFLAG


      integer K, NM, NCMPLX, IK
      real*8 EPS, NPP, CNV
      Common /CBLKT/ K, NM, NCMPLX, IK, NPP, CNV, EPS


      integer  I2,IR,IM2,NM2, I1, IRM1, IM3, NM3
      integer I3, IM1, NM1, KDO, L, I4, IF,I, IPI1, IPI2, IPI3
      integer NC, IDXA, IP2, NA, NP2, IP1, NP1, J, IP3, NP3
      integer IZ, NZ, IDXC, IZR, LL, IFD, IP, NP, IMI1, IMI2

c
c Added this to avoid compiler warning on g77
c John Evans
      Complex cbip, cw1, cw3, cww


c
c begin reduction phase
c
      KDO = K - 1
      Do 40 L = 1, KDO
        IR = L - 1
        I2 = 2 ** IR
        I1 = I2 / 2
        I3 = I2 + I1
        I4 = I2 + I2
        IRM1 = IR - 1
        Call INDXB(I2,IR,IM2,NM2)
        Call INDXB(I1,IRM1,IM3,NM3)
        Call INDXB(I3,IRM1,IM1,NM1)

        If (PFLAG.eq.1) Then
            Call PROD(NM2,B(IM2),NM3,B(IM3),NM1,B(IM1),0,DUM,Y(1,I2),
     *                 W3, M, AM,BM,CM,WD,WW,WU)
        Else
            Call PRODP(NM2,B(IM2),NM3,B(IM3),NM1,B(IM1),0,DUM,Y(1,I2),
     *                 W3,M,AM,BM,CM,WD,WW,WU)
        End If
c
c        Call PRDCT(NM2,B(IM2),NM3,B(IM3),NM1,B(IM1),0,DUM,Y(1,I2),W3,M,
c     *      AM,BM,CM,WD,WW,WU)
c

        IF = 2 ** K
        Do 30 I = I4, IF, I4
          If (I.LE.NM) Then
            IPI1 = I + I1
            IPI2 = I + I2
            IPI3 = I + I3
            Call INDXC(I,IR,IDXC,NC)
            If (I.LT.IF) Then
              Call INDXA(I,IR,IDXA,NA)
              Call INDXB(I-I1,IRM1,IM1,NM1)
              Call INDXB(IPI2,IR,IP2,NP2)
              Call INDXB(IPI1,IRM1,IP1,NP1)
              Call INDXB(IPI3,IRM1,IP3,NP3)

c              Call PRDCT(NM1,B(IM1),0,DUM,0,DUM,NA,AN(IDXA),W3,W1,M,AM,
c     *            BM,CM,WD,WW,WU)

              If (PFLAG .eq.1) Then
                Call PROD(NM1,B(IM1),0,DUM,0,DUM,NA,AN(IDXA),W3,W1,
     *                   M,AM,BM,CM,WD,WW,WU)
              Else
                Call PRODP(NM1,B(IM1),0,DUM,0,DUM,NA,AN(IDXA),W3,W1,
     *                    M,AM,BM,CM,WD,WW,WU)
              End If


              If (IPI2.GT.NM) Then
                Do 10 J = 1, M
                  W3(J) = 0.
                  W2(J) = 0.
   10           Continue
              Else

c                Call PRDCT(NP2,B(IP2),NP1,B(IP1),NP3,B(IP3),0,DUM,
c     *              Y(1,IPI2),W3,M,AM,BM,CM,WD,WW,WU)
c                Call PRDCT(NP1,B(IP1),0,DUM,0,DUM,NC,CN(IDXC),W3,W2,M,
c     *              AM,BM,CM,WD,WW,WU)

                If (PFLAG.eq.1) Then
                    Call PROD(NP2,B(IP2),NP1,B(IP1),NP3,B(IP3),0,DUM,
     *                  Y(1,IPI2),W3,M,AM,BM,CM,WD,WW,WU)
                    Call PROD(NP1,B(IP1),0,DUM,0,DUM,NC,CN(IDXC),W3,
     *                   W2,M,AM,BM,CM,WD,WW,WU)
                Else
                    Call PRODP(NP2,B(IP2),NP1,B(IP1),NP3,B(IP3),0,DUM,
     *                  Y(1,IPI2),W3,M,AM,BM,CM,WD,WW,WU)
                    Call PRODP(NP1,B(IP1),0,DUM,0,DUM,NC,CN(IDXC),W3,
     *                   W2,M,AM,BM,CM,WD,WW,WU)
                End If

              End If
              Do 20 J = 1, M
                Y(J,I) = W1(J) + W2(J) + Y(J,I)
   20         Continue
            End If
          End If
   30   Continue
   40 Continue
      If (NPP.EQ.0) Then
c
c     the periodic case is treated using the capacitance matrix method
c
        IF = 2 ** K
        I = IF / 2
        I1 = I / 2
        Call INDXB(I-I1,K-2,IM1,NM1)
        Call INDXB(I+I1,K-2,IP1,NP1)
        Call INDXB(I,K-1,IZ,NZ)



c        Call PRDCT(NZ,B(IZ),NM1,B(IM1),NP1,B(IP1),0,DUM,Y(1,I),W1,M,AM,
c     *      BM,CM,WD,WW,WU)
        If (PFLAG.eq.1) Then
            Call PROD(NZ,B(IZ),NM1,B(IM1),NP1,B(IP1),0,DUM,Y(1,I),W1,
     *                M,AM,BM,CM,WD,WW,WU)
        Else
            Call PRODP(NZ,B(IZ),NM1,B(IM1),NP1,B(IP1),0,DUM,Y(1,I),
     *                W1,M,AM,BM,CM,WD,WW,WU)
        End If


        IZR = I
        Do 50 J = 1, M
          W2(J) = W1(J)
   50   Continue
        Do 70 LL = 2, K
          L = K - LL + 1
          IR = L - 1
          I2 = 2 ** IR
          I1 = I2 / 2
          I = I2
          Call INDXC(I,IR,IDXC,NC)
          Call INDXB(I,IR,IZ,NZ)
          Call INDXB(I-I1,IR-1,IM1,NM1)
          Call INDXB(I+I1,IR-1,IP1,NP1)



c          Call PRDCT(NP1,B(IP1),0,DUM,0,DUM,NC,CN(IDXC),W1,W1,M,AM,BM,
c     *        CM,WD,WW,WU)
          If (PFLAG.eq.1) Then
              Call PROD(NP1,B(IP1),0,DUM,0,DUM,NC,CN(IDXC),W1,W1,M,
     *                  AM,BM,CM,WD,WW,WU)
          Else
              Call PRODP(NP1,B(IP1),0,DUM,0,DUM,NC,CN(IDXC),W1,W1,
     *                  M,AM,BM,CM,WD,WW,WU)
          End If



          Do 60 J = 1, M
            W1(J) = Y(J,I) + W1(J)
   60     Continue



c          Call PRDCT(NZ,B(IZ),NM1,B(IM1),NP1,B(IP1),0,DUM,W1,W1,M,AM,BM,
c     *        CM,WD,WW,WU)
          If (PFLAG.eq.1) Then
              Call PROD(NZ,B(IZ),NM1,B(IM1),NP1,B(IP1),0,DUM,W1,W1,M,
     *                  AM,BM,CM,WD,WW,WU)
          Else
              Call PRODP(NZ,B(IZ),NM1,B(IM1),NP1,B(IP1),0,DUM,W1,W1,
     *                  M,AM,BM,CM,WD,WW,WU)
          End If



   70   Continue
        Do 100 LL = 2, K
          L = K - LL + 1
          IR = L - 1
          I2 = 2 ** IR
          I1 = I2 / 2
          I4 = I2 + I2
          IFD = IF - I2
          Do 90 I = I2, IFD, I4
            If (I-I2-IZR.EQ.0) Then
              If (I.GT.NM) Go To 100
              Call INDXA(I,IR,IDXA,NA)
              Call INDXB(I,IR,IZ,NZ)
              Call INDXB(I-I1,IR-1,IM1,NM1)
              Call INDXB(I+I1,IR-1,IP1,NP1)



c              Call PRDCT(NM1,B(IM1),0,DUM,0,DUM,NA,AN(IDXA),W2,W2,M,AM,
c     *            BM,CM,WD,WW,WU)
              If (PFLAG.eq.1) Then
                  Call PROD(NM1,B(IM1),0,DUM,0,DUM,NA,AN(IDXA),W2,W2,
     *                      M,AM,BM,CM,WD,WW,WU)
              Else
                  Call PRODP(NM1,B(IM1),0,DUM,0,DUM,NA,AN(IDXA),W2,
     *                      W2,M,AM,BM,CM,WD,WW,WU)
              End If



              Do 80 J = 1, M
                W2(J) = Y(J,I) + W2(J)
   80         Continue




c              Call PRDCT(NZ,B(IZ),NM1,B(IM1),NP1,B(IP1),0,DUM,W2,W2,M,
c     *            AM,BM,CM,WD,WW,WU)
              If (PFLAG.eq.1) Then
                  Call PROD(NZ,B(IZ),NM1,B(IM1),NP1,B(IP1),0,DUM,W2,
     *                      W2,M,AM,BM,CM,WD,WW,WU)
              Else
                  Call PRODP(NZ,B(IZ),NM1,B(IM1),NP1,B(IP1),0,DUM,
     *                       W2,W2,M,AM,BM,CM,WD,WW,WU)
              End If




              IZR = I
              If (I.EQ.NM) Go To 110
            End If
   90     Continue
  100   Continue
  110   Do 120 J = 1, M
          Y(J,NM+1) = Y(J,NM+1) - CN(NM+1) * W1(J) - AN(NM+1) * W2(J)
  120   Continue
        Call INDXB(IF/2,K-1,IM1,NM1)
        Call INDXB(IF,K-1,IP,NP)



        If (NCMPLX.NE.0) Then

          cbip = CMPLX ( B(IP) )
          cw1 = CMPLX ( W1(1) )
          cw3 = CMPLX ( W3(1) )
          cww = CMPLX ( WW(1) )

c          Call CPRDCT(NM+1,B(IP),NM1,B(IM1),0,DUM,0,DUM,Y(1,NM+1),
c     *        Y(1,NM+1),M,AM,BM,CM,W1,W3,WW)
          If (PFLAG.eq.1) Then
              Call CPROD(NM+1, cbip, NM1,B(IM1),0,DUM,0,DUM,
     *                   Y(1,NM+1),Y(1,NM+1),M,AM,BM,CM, cw1, cw3, cww)
          Else
              Call CPRODP(NM+1,cbip,NM1,B(IM1),0,DUM,0,DUM,
     *                    Y(1,NM+1),Y(1,NM+1),M,AM,BM,CM, cw1, cw3, cww)
          End If


        Else

c          Call PRDCT(NM+1,B(IP),NM1,B(IM1),0,DUM,0,DUM,Y(1,NM+1),
c     *        Y(1,NM+1),M,AM,BM,CM,WD,WW,WU)
          If (PFLAG.eq.1) Then
              Call PROD(NM+1,B(IP),NM1,B(IM1),0,DUM,0,DUM,Y(1,NM+1),
     *            Y(1,NM+1),M,AM,BM,CM,WD,WW,WU)
          Else
              Call PRODP(NM+1,B(IP),NM1,B(IM1),0,DUM,0,DUM,Y(1,NM+1),
     *            Y(1,NM+1),M,AM,BM,CM,WD,WW,WU)
          End If


        End If




        Do 130 J = 1, M
          W1(J) = AN(1) * Y(J,NM+1)
          W2(J) = CN(NM) * Y(J,NM+1)
          Y(J,1) = Y(J,1) - W1(J)
          Y(J,NM) = Y(J,NM) - W2(J)
  130   Continue
        Do 150 L = 1, KDO
          IR = L - 1
          I2 = 2 ** IR
          I4 = I2 + I2
          I1 = I2 / 2
          I = I4
          Call INDXA(I,IR,IDXA,NA)
          Call INDXB(I-I2,IR,IM2,NM2)
          Call INDXB(I-I2-I1,IR-1,IM3,NM3)
          Call INDXB(I-I1,IR-1,IM1,NM1)


c         Call PRDCT(NM2,B(IM2),NM3,B(IM3),NM1,B(IM1),0,DUM,W1,W1,M,AM,
c    *        BM,CM,WD,WW,WU)
c         Call PRDCT(NM1,B(IM1),0,DUM,0,DUM,NA,AN(IDXA),W1,W1,M,AM,BM,
c    *        CM,WD,WW,WU)
          If (PFLAG.eq.1) Then
              Call PROD(NM2,B(IM2),NM3,B(IM3),NM1,B(IM1),0,DUM,W1,W1,
     *                  M,AM,BM,CM,WD,WW,WU)
              Call PROD(NM1,B(IM1),0,DUM,0,DUM,NA,AN(IDXA),W1,W1,M,
     *                  AM,BM,CM,WD,WW,WU)
          Else
              Call PRODP(NM2,B(IM2),NM3,B(IM3),NM1,B(IM1),0,DUM,W1,
     *                   W1,M,AM,BM,CM,WD,WW,WU)
              Call PRODP(NM1,B(IM1),0,DUM,0,DUM,NA,AN(IDXA),W1,W1,
     *                   M,AM,BM,CM,WD,WW,WU)
          End If




          Do 140 J = 1, M
            Y(J,I) = Y(J,I) - W1(J)
  140     Continue
  150   Continue
c
        IZR = NM
        Do 180 L = 1, KDO
          IR = L - 1
          I2 = 2 ** IR
          I1 = I2 / 2
          I3 = I2 + I1
          I4 = I2 + I2
          IRM1 = IR - 1
          Do 170 I = I4, IF, I4
            IPI1 = I + I1
            IPI2 = I + I2
            IPI3 = I + I3
            If (IPI2.NE.IZR) Then
              If (I-IZR) 170, 180, 170
            End If
            Call INDXC(I,IR,IDXC,NC)
            Call INDXB(IPI2,IR,IP2,NP2)
            Call INDXB(IPI1,IRM1,IP1,NP1)
            Call INDXB(IPI3,IRM1,IP3,NP3)




c            Call PRDCT(NP2,B(IP2),NP1,B(IP1),NP3,B(IP3),0,DUM,W2,W2,M,
c     *          AM,BM,CM,WD,WW,WU)
c            Call PRDCT(NP1,B(IP1),0,DUM,0,DUM,NC,CN(IDXC),W2,W2,M,AM,BM,
c     *          CM,WD,WW,WU)
            If (PFLAG.eq.1) Then
                Call PROD(NP2,B(IP2),NP1,B(IP1),NP3,B(IP3),0,DUM,W2,
     *                    W2,M,AM,BM,CM,WD,WW,WU)
                Call PROD(NP1,B(IP1),0,DUM,0,DUM,NC,CN(IDXC),W2,W2,M,
     *                    AM,BM,CM,WD,WW,WU)
            Else
                Call PRODP(NP2,B(IP2),NP1,B(IP1),NP3,B(IP3),0,DUM,W2,
     *                     W2,M,AM,BM,CM,WD,WW,WU)
                Call PRODP(NP1,B(IP1),0,DUM,0,DUM,NC,CN(IDXC),W2,W2,M,
     *                     AM,BM,CM,WD,WW,WU)
            End If



            Do 160 J = 1, M
              Y(J,I) = Y(J,I) - W2(J)
  160       Continue
            IZR = I
            Go To 180
  170     Continue
  180   Continue
      End If
c
c begin back substitution phase
c
      Do 230 LL = 1, K
        L = K - LL + 1
        IR = L - 1
        IRM1 = IR - 1
        I2 = 2 ** IR
        I1 = I2 / 2
        I4 = I2 + I2
        IFD = IF - I2
        Do 220 I = I2, IFD, I4
          If (I.LE.NM) Then
            IMI1 = I - I1
            IMI2 = I - I2
            IPI1 = I + I1
            IPI2 = I + I2
            Call INDXA(I,IR,IDXA,NA)
            Call INDXC(I,IR,IDXC,NC)
            Call INDXB(I,IR,IZ,NZ)
            Call INDXB(IMI1,IRM1,IM1,NM1)
            Call INDXB(IPI1,IRM1,IP1,NP1)



            If (I.LE.I2) Then
              Do 190 J = 1, M
                W1(J) = 0.
  190         Continue

            Else

c              Call PRDCT(NM1,B(IM1),0,DUM,0,DUM,NA,AN(IDXA),Y(1,IMI2),
c     *            W1,M,AM,BM,CM,WD,WW,WU)
              If (PFLAG.eq.1) Then
                  Call PROD(NM1,B(IM1),0,DUM,0,DUM,NA,AN(IDXA),
     *                      Y(1,IMI2),W1,M,AM,BM,CM,WD,WW,WU)
              Else
                  Call PRODP(NM1,B(IM1),0,DUM,0,DUM,NA,AN(IDXA),
     *                       Y(1,IMI2),W1,M,AM,BM,CM,WD,WW,WU)
              End If

            End If



            If (IPI2.GT.NM) Then

              Do 200 J = 1, M
                W2(J) = 0.
  200         Continue

            Else

c              Call PRDCT(NP1,B(IP1),0,DUM,0,DUM,NC,CN(IDXC),Y(1,IPI2),
c     *            W2,M,AM,BM,CM,WD,WW,WU)
              If (PFLAG.eq.1) Then
                  Call PROD(NP1,B(IP1),0,DUM,0,DUM,NC,CN(IDXC),
     *                      Y(1,IPI2),W2,M,AM,BM,CM,WD,WW,WU)
              Else
                  Call PRODP(NP1,B(IP1),0,DUM,0,DUM,NC,CN(IDXC),
     *                       Y(1,IPI2),W2,M,AM,BM,CM,WD,WW,WU)
              End If


            End If



            Do 210 J = 1, M
              W1(J) = Y(J,I) + W1(J) + W2(J)
  210       Continue


c            Call PRDCT(NZ,B(IZ),NM1,B(IM1),NP1,B(IP1),0,DUM,W1,Y(1,I),M,
c     *          AM,BM,CM,WD,WW,WU)
            If (PFLAG.eq.1) Then
                Call PROD(NZ,B(IZ),NM1,B(IM1),NP1,B(IP1),0,DUM,W1,
     *                    Y(1,I),M,AM,BM,CM,WD,WW,WU)
            Else
                Call PRODP(NZ,B(IZ),NM1,B(IM1),NP1,B(IP1),0,DUM,W1,
     *                     Y(1,I),M,AM,BM,CM,WD,WW,WU)
            End If


          End If


  220   Continue
  230 Continue
      Return
      End

















      Subroutine INDXB(I,IR,IDX,IDP)
c
c b(idx) is the location of the first root of the b(i,ir) polynomial
c
      IMPLICIT NONE

      integer IDP, IDX, I, IR
      integer IZH, ID, IPL

      integer K, NM, NCMPLX, IK
      real*8 EPS, NPP, CNV
      Common /CBLKT/ K, NM, NCMPLX, IK, NPP, CNV, EPS
      IDP = 0
      If (IR) 30, 10, 20
   10 If (I.GT.NM) Go To 30
      IDX = I
      IDP = 1
      Return
   20 IZH = 2 ** IR
      ID = I - IZH - IZH
      IDX = ID + ID + (IR-1) * IK + IR + (IK-I) / IZH + 4
      IPL = IZH - 1
      IDP = IZH + IZH - 1
      If (I-IPL-NM.GT.0) Then
        IDP = 0
        Return
      End If
      If (I+IPL-NM.GT.0) Then
        IDP = NM + IPL - I + 1
      End If
   30 Return
      End








      Subroutine INDXA(I,IR,IDXA,NA)
      IMPLICIT NONE

      integer K, NM, NCMPLX, IK
      real*8 EPS, NPP, CNV
      Common /CBLKT/ K, NM, NCMPLX, IK, NPP, CNV, EPS

      integer NA, IDXA, I, IR



      NA = 2 ** IR
      IDXA = I - NA + 1
      If (I.GT.NM) Then
        NA = 0
      End If
      Return
      End





      Subroutine INDXC(I,IR,IDXC,NC)

      IMPLICIT NONE
      integer K, NM, NCMPLX, IK, NC, IR, I, IDXC
      real*8 EPS, NPP, CNV
      Common /CBLKT/ K, NM, NCMPLX, IK, NPP, CNV, EPS

      NC = 2 ** IR
      IDXC = I
      If (IDXC+NC-1-NM.GT.0) Then
        NC = 0
      End If
      Return
      End





      Subroutine PROD(ND,BD,NM1,BM1,NM2,BM2,NA,AA,X,Y,M,A,B,C,D,W,U)
c
c prod applies a sequence of matrix operations to the vector x and
c stores the result in y
c bd,bm1,bm2 are arrays containing roots of certian b polynomials
c nd,nm1,nm2 are the lengths of the arrays bd,bm1,bm2 respectively
c aa   array containing scalar multipliers of the vector x
c na is the length of the array aa
c x,y  the matrix operations are applied to x and the result is y
c a,b,c  are arrays which contain the tridiagonal matrix
c m  is the order of the matrix
c d,w,u are working arrays
c is  determines whether or not a change in sign is made
c
      IMPLICIT NONE

      Dimension A(1), B(1), C(1), X(1), Y(1), D(2), W(2), BD(1), 
     *    BM1(1), BM2(1), AA(1), U(1)

      real*8 bd,bm1, bm2, aa, a, b, c, d, u, w, x, y, rt
      real*8 den
      integer nd, nm1, nm2, na, m, j, mm, id, ibr, m1, m2, ia
      integer k



      Do 10 J = 1, M
        W(J) = X(J)
        Y(J) = W(J)
   10 Continue
      MM = M - 1
      ID = ND
      IBR = 0
      M1 = NM1
      M2 = NM2
      IA = NA
   20 If (IA.GT.0) Then
        RT = AA(IA)
        If (ND.EQ.0) RT = -RT
        IA = IA - 1
c
c scalar multiplication
c
        Do 30 J = 1, M
          Y(J) = RT * W(J)
   30   Continue
      End If
      If (ID.GT.0) Then
        RT = BD(ID)
        ID = ID - 1
        If (ID.EQ.0) IBR = 1
c
c begin solution to system
c
        D(M) = A(M) / (B(M)-RT)
        W(M) = Y(M) / (B(M)-RT)
        Do 40 J = 2, MM
          K = M - J
          DEN = B(K+1) - RT - C(K+1) * D(K+2)
          D(K+1) = A(K+1) / DEN
          W(K+1) = (Y(K+1)-C(K+1)*W(K+2)) / DEN
   40   Continue
        DEN = B(1) - RT - C(1) * D(2)
        W(1) = 1.
        If (DEN.NE.0) Then
          W(1) = (Y(1)-C(1)*W(2)) / DEN
        End If
        Do 50 J = 2, M
          W(J) = W(J) - D(J) * W(J-1)
   50   Continue
        If (NA) 80, 80, 20
   60   Do 70 J = 1, M
          Y(J) = W(J)
   70   Continue
        IBR = 1
        Go To 20
   80   If (M1.LE.0) Then
          If (M2) 60, 60, 90
        End If
        If (M2.GT.0) Then
          If (ABS(BM1(M1)).LE.ABS(BM2(M2))) Go To 90
        End If
        If (IBR.LE.0) Then
          If (ABS(BM1(M1)-BD(ID)).LT.ABS(BM1(M1)-RT)) Go To 60
        End If
        RT = RT - BM1(M1)
        M1 = M1 - 1
        Go To 100
   90   If (IBR.LE.0) Then
          If (ABS(BM2(M2)-BD(ID)).LT.ABS(BM2(M2)-RT)) Go To 60
        End If
        RT = RT - BM2(M2)
        M2 = M2 - 1
  100   Do 110 J = 1, M
          Y(J) = Y(J) + RT * W(J)
  110   Continue
        Go To 20
      End If
      Return
      End





      Subroutine PRODP(ND,BD,NM1,BM1,NM2,BM2,NA,AA,X,Y,M,A,B,C,D,U,W)
c
c prodp applies a sequence of matrix operations to the vector x and
c stores the result in y        periodic boundary conditions
c
c bd,bm1,bm2 are arrays containing roots of certian b polynomials
c nd,nm1,nm2 are the lengths of the arrays bd,bm1,bm2 respectively
c aa   array containing scalar multipliers of the vector x
c na is the length of the array aa
c x,y  the matrix operations are applied to x and the result is y
c a,b,c  are arrays which contain the tridiagonal matrix
c m  is the order of the matrix
c d,u,w are working arrays
c is  determines whether or not a change in sign is made
c
      IMPLICIT NONE



      Dimension A(1), B(1), C(1), X(1), Y(1), D(1), U(1), BD(1), 
     *    BM1(1), BM2(1), AA(1), W(1)


      real*8 bd,bm1, bm2, aa, a, b, c, d, u, w, x, y, rt, bh, ym, v
      real*8 den, am
      integer nd, nm1, nm2, na, m, j, mm, mm2, id, ibr, m1, m2, ia
      integer k


      Do 10 J = 1, M
        Y(J) = X(J)
        W(J) = Y(J)
   10 Continue
      MM = M - 1
      MM2 = M - 2
      ID = ND
      IBR = 0
      M1 = NM1
      M2 = NM2
      IA = NA
   20 If (IA.GT.0) Then
        RT = AA(IA)
        If (ND.EQ.0) RT = -RT
        IA = IA - 1
        Do 30 J = 1, M
          Y(J) = RT * W(J)
   30   Continue
      End If
      If (ID.GT.0) Then
        RT = BD(ID)
        ID = ID - 1
        If (ID.EQ.0) IBR = 1
c
c begin solution to system
c
        BH = B(M) - RT
        YM = Y(M)
        DEN = B(1) - RT
        D(1) = C(1) / DEN
        U(1) = A(1) / DEN
        W(1) = Y(1) / DEN
        V = C(M)
        If (MM2.GE.2) Then
          Do 40 J = 2, MM2
            DEN = B(J) - RT - A(J) * D(J-1)
            D(J) = C(J) / DEN
            U(J) = -A(J) * U(J-1) / DEN
            W(J) = (Y(J)-A(J)*W(J-1)) / DEN
            BH = BH - V * U(J-1)
            YM = YM - V * W(J-1)
            V = -V * D(J-1)
   40     Continue
        End If
        DEN = B(M-1) - RT - A(M-1) * D(M-2)
        D(M-1) = (C(M-1)-A(M-1)*U(M-2)) / DEN
        W(M-1) = (Y(M-1)-A(M-1)*W(M-2)) / DEN
        AM = A(M) - V * D(M-2)
        BH = BH - V * U(M-2)
        YM = YM - V * W(M-2)
        DEN = BH - AM * D(M-1)
        If (DEN.NE.0) Then
          W(M) = (YM-AM*W(M-1)) / DEN
        Else
          W(M) = 1.
        End If
        W(M-1) = W(M-1) - D(M-1) * W(M)
        Do 50 J = 2, MM
          K = M - J
          W(K) = W(K) - D(K) * W(K+1) - U(K) * W(M)
   50   Continue
        If (NA) 80, 80, 20
   60   Do 70 J = 1, M
          Y(J) = W(J)
   70   Continue
        IBR = 1
        Go To 20
   80   If (M1.LE.0) Then
          If (M2) 60, 60, 90
        End If
        If (M2.GT.0) Then
          If (ABS(BM1(M1)).LE.ABS(BM2(M2))) Go To 90
        End If
        If (IBR.LE.0) Then
          If (ABS(BM1(M1)-BD(ID)).LT.ABS(BM1(M1)-RT)) Go To 60
        End If
        RT = RT - BM1(M1)
        M1 = M1 - 1
        Go To 100
   90   If (IBR.LE.0) Then
          If (ABS(BM2(M2)-BD(ID)).LT.ABS(BM2(M2)-RT)) Go To 60
        End If
        RT = RT - BM2(M2)
        M2 = M2 - 1
  100   Do 110 J = 1, M
          Y(J) = Y(J) + RT * W(J)
  110   Continue
        Go To 20
      End If
      Return
      End

















      Subroutine CPROD(ND,BD,NM1,BM1,NM2,BM2,NA,AA,X,YY,M,A,B,C,D,W,Y)
c
c prod applies a sequence of matrix operations to the vector x and
c stores the result in yy           (complex case)
c aa   array containing scalar multipliers of the vector x
c nd,nm1,nm2 are the lengths of the arrays bd,bm1,bm2 respectively
c bd,bm1,bm2 are arrays containing roots of certian b polynomials
c na is the length of the array aa
c x,yy the matrix operations are applied to x and the result is yy
c a,b,c  are arrays which contain the tridiagonal matrix
c m  is the order of the matrix
c d,w,y are working arrays
c isgn  determines whether or not a change in sign is made
c
      IMPLICIT NONE
      Complex Y, D, W, BD, CRT, DEN, Y1, Y2
      Dimension A(1), B(1), C(1), X(1), Y(2), D(2), W(2), BD(1), 
     *    BM1(1), BM2(1), AA(1), YY(1)



      integer J, M, MM, ND, NM1, NM2, ID, M1, M2, IA
      integer NA, IFLG, K
      real*8 BM1, BM2, AA, X, YY, A, B, C, RT


      Do 10 J = 1, M
        Y(J) = DCMPLX(X(J),0.d0)
   10 Continue
      MM = M - 1
      ID = ND
      M1 = NM1
      M2 = NM2
      IA = NA
   20 IFLG = 0
      If (ID.GT.0) Then
        CRT = BD(ID)
        ID = ID - 1
c
c begin solution to system
c
        D(M) = A(M) / (B(M)-CRT)
        W(M) = Y(M) / (B(M)-CRT)
        Do 30 J = 2, MM
          K = M - J
          DEN = B(K+1) - CRT - C(K+1) * D(K+2)
          D(K+1) = A(K+1) / DEN
          W(K+1) = (Y(K+1)-C(K+1)*W(K+2)) / DEN
   30   Continue
        DEN = B(1) - CRT - C(1) * D(2)
        If (CABS(DEN).NE.0) Then
          Y(1) = (Y(1)-C(1)*W(2)) / DEN
        Else
          Y(1) = (1.,0.)
        End If
        Do 40 J = 2, M
          Y(J) = W(J) - D(J) * Y(J-1)
   40   Continue
      End If
      If (M1.LE.0) Then
        If (M2.LE.0) Go To 60
        RT = BM2(M2)
        M2 = M2 - 1
      Else
        If (M2.LE.0) Then
          RT = BM1(M1)
          M1 = M1 - 1
        Else
          If (ABS(BM1(M1)).GT.ABS(BM2(M2))) Then
            RT = BM1(M1)
            M1 = M1 - 1
          Else
            RT = BM2(M2)
            M2 = M2 - 1
          End If
        End If
      End If
      Y1 = (B(1)-RT) * Y(1) + C(1) * Y(2)
      If (MM.GE.2) Then
c
c matrix multiplication
c
        Do 50 J = 2, MM
          Y2 = A(J) * Y(J-1) + (B(J)-RT) * Y(J) + C(J) * Y(J+1)
          Y(J-1) = Y1
          Y1 = Y2
   50   Continue
      End If
      Y(M) = A(M) * Y(M-1) + (B(M)-RT) * Y(M)
      Y(M-1) = Y1
      IFLG = 1
      Go To 20
   60 If (IA.GT.0) Then
        RT = AA(IA)
        IA = IA - 1
        IFLG = 1
c
c scalar multiplication
c
        Do 70 J = 1, M
          Y(J) = RT * Y(J)
   70   Continue
      End If
      If (IFLG.GT.0) Go To 20
      Do 80 J = 1, M
        YY(J) = REAL(Y(J))
   80 Continue
      Return
      End











      Subroutine CPRODP(ND,BD,NM1,BM1,NM2,BM2,NA,AA,X,YY,M,A,B,C,D,U,Y)
c
c prodp applies a sequence of matrix operations to the vector x and
c stores the result in yy       periodic boundary conditions
c and  complex  case
c
c bd,bm1,bm2 are arrays containing roots of certian b polynomials
c nd,nm1,nm2 are the lengths of the arrays bd,bm1,bm2 respectively
c aa   array containing scalar multipliers of the vector x
c na is the length of the array aa
c x,yy the matrix operations are applied to x and the result is yy
c a,b,c  are arrays which contain the tridiagonal matrix
c m  is the order of the matrix
c d,u,y are working arrays
c isgn  determines whether or not a change in sign is made
c
      IMPLICIT NONE

      Complex Y, D, U, V, DEN, BH, YM, AM, Y1, Y2, YH, BD, CRT
      Dimension A(1), B(1), C(1), X(1), Y(2), D(1), U(1), BD(1), 
     *    BM1(1), BM2(1), AA(1), YY(1)


      integer J, M, MM, MM2, ND, NM1, NM2, ID, M1, M2, IA
      integer NA, IFLG, K
      real*8 BM1, BM2, AA, X, YY, A, B, C, RT



      Do 10 J = 1, M
        Y(J) = DCMPLX(X(J),0.d0)
   10 Continue
      MM = M - 1
      MM2 = M - 2
      ID = ND
      M1 = NM1
      M2 = NM2
      IA = NA
   20 IFLG = 0
      If (ID.GT.0) Then
        CRT = BD(ID)
        ID = ID - 1
        IFLG = 1
c
c begin solution to system
c
        BH = B(M) - CRT
        YM = Y(M)
        DEN = B(1) - CRT
        D(1) = C(1) / DEN
        U(1) = A(1) / DEN
        Y(1) = Y(1) / DEN
        V = DCMPLX(C(M),0.d0)
        If (MM2.GE.2) Then
          Do 30 J = 2, MM2
            DEN = B(J) - CRT - A(J) * D(J-1)
            D(J) = C(J) / DEN
            U(J) = -A(J) * U(J-1) / DEN
            Y(J) = (Y(J)-A(J)*Y(J-1)) / DEN
            BH = BH - V * U(J-1)
            YM = YM - V * Y(J-1)
            V = -V * D(J-1)
   30     Continue
        End If
        DEN = B(M-1) - CRT - A(M-1) * D(M-2)
        D(M-1) = (C(M-1)-A(M-1)*U(M-2)) / DEN
        Y(M-1) = (Y(M-1)-A(M-1)*Y(M-2)) / DEN
        AM = A(M) - V * D(M-2)
        BH = BH - V * U(M-2)
        YM = YM - V * Y(M-2)
        DEN = BH - AM * D(M-1)
        If (CABS(DEN).NE.0) Then
          Y(M) = (YM-AM*Y(M-1)) / DEN
        Else
          Y(M) = (1.,0.)
        End If
        Y(M-1) = Y(M-1) - D(M-1) * Y(M)
        Do 40 J = 2, MM
          K = M - J
          Y(K) = Y(K) - D(K) * Y(K+1) - U(K) * Y(M)
   40   Continue
      End If
      If (M1.LE.0) Then
        If (M2.LE.0) Go To 60
        RT = BM2(M2)
        M2 = M2 - 1
      Else
        If (M2.LE.0) Then
          RT = BM1(M1)
          M1 = M1 - 1
        Else
          If (DABS(BM1(M1)).GT.DABS(BM2(M2))) Then
            RT = BM1(M1)
            M1 = M1 - 1
          Else
            RT = BM2(M2)
            M2 = M2 - 1
          End If
        End If
      End If
c
c matrix multiplication
c
      YH = Y(1)
      Y1 = (B(1)-RT) * Y(1) + C(1) * Y(2) + A(1) * Y(M)
      If (MM.GE.2) Then
        Do 50 J = 2, MM
          Y2 = A(J) * Y(J-1) + (B(J)-RT) * Y(J) + C(J) * Y(J+1)
          Y(J-1) = Y1
          Y1 = Y2
   50   Continue
      End If
      Y(M) = A(M) * Y(M-1) + (B(M)-RT) * Y(M) + C(M) * YH
      Y(M-1) = Y1
      IFLG = 1
      Go To 20
   60 If (IA.GT.0) Then
        RT = AA(IA)
        IA = IA - 1
        IFLG = 1
c
c scalar multiplication
c
        Do 70 J = 1, M
          Y(J) = RT * Y(J)
   70   Continue
      End If
      If (IFLG.GT.0) Go To 20
      Do 80 J = 1, M
        YY(J) = REAL(Y(J))
   80 Continue
      Return
      End





      Subroutine PPADD(N,IERROR,A,C,CBP,BP,BH)
c
c     ppadd computes the eigenvalues of the periodic tridiagonal matrix
c     with coefficients an,bn,cn
c
c n is the order of the bh and bp polynomials
c     on output bp contians the eigenvalues
c cbp is the same as bp except type complex
c bh is used to temporarily store the roots of the b hat polynomial
c which enters through bp
c
      IMPLICIT NONE
      DOUBLE Complex CX, FSG, HSG, DD
      DOUBLE Complex F, FP, FPP, CDIS, R1, R2, R3, CBP
      Dimension A(1), C(1), BP(1), BH(3), CBP(1)
      real*8 A, C, BP, BH

      real*8 XR, XM



      real*8 SGN

      integer K, NM, NCMPLX, IK
      real*8 EPS, NPP, XL, CNV
      Common /CBLKT/ K, NM, NCMPLX, IK, NPP, CNV, EPS
      real*8 PSGF, PPSPF, PPSGF

      External PSGF, PPSPF, PPSGF

      real*8 SCNV, DB, BSRH, PSG
      integer IZ, N, IZM, IZM2, J, NT, MODIZ, IS, IF, IG, IT, ICV, I2
      integer I3, NHALF, IERROR


      SCNV = DSQRT(CNV)
      IZ = N
      IZM = IZ - 1
      IZM2 = IZ - 2
      If (BP(N)-BP(1)) 10, 240, 30
   10 Do 20 J = 1, N
        NT = N - J
        BH(J) = BP(NT+1)
   20 Continue
      Go To 50
   30 Do 40 J = 1, N
        BH(J) = BP(J)
   40 Continue
   50 NCMPLX = 0
      MODIZ = MOD(IZ,2)
      IS = 1
      If (MODIZ.NE.0) Then
        If (A(1)) 80, 240, 60
      End If
   60 XL = BH(1)
      DB = BH(3) - BH(1)
   70 XL = XL - DB
      If (PSGF(XL,IZ,C,A,BH).LE.0) Go To 70
      SGN = -1.
      CBP(1) = DCMPLX(BSRH(XL,BH(1),IZ,C,A,BH,PSGF,SGN),0.d0)
      IS = 2
   80 IF = IZ - 1
      If (MODIZ.NE.0) Then
        If (A(1)) 90, 240, 110
      End If
   90 XR = BH(IZ)
      DB = BH(IZ) - BH(IZ-2)
  100 XR = XR + DB
      If (PSGF(XR,IZ,C,A,BH).LT.0) Go To 100
      SGN = 1.
      CBP(IZ) = DCMPLX(BSRH(BH(IZ),XR,IZ,C,A,BH,PSGF,SGN),0.d0)
      IF = IZ - 2
  110 Do 180 IG = IS, IF, 2
        XL = BH(IG)
        XR = BH(IG+1)
        SGN = -1.
        XM = BSRH(XL,XR,IZ,C,A,BH,PPSPF,SGN)
        PSG = PSGF(XM,IZ,C,A,BH)
        If (ABS(PSG).GT.EPS) Then
          If (PSG*PPSGF(XM,IZ,C,A,BH)) 120, 130, 140
c
c     case of a real zero
c
  120     SGN = 1.
          CBP(IG) = DCMPLX(BSRH(BH(IG),XM,IZ,C,A,BH,PSGF,SGN),0.d0)
          SGN = -1.
          CBP(IG+1) = DCMPLX(BSRH(XM,BH(IG+1),IZ,C,A,BH,PSGF,SGN),0.d0)
          Go To 180
        End If
c
c     case of a multiple zero
c
  130   CBP(IG) = DCMPLX(XM,0.d0)
        CBP(IG+1) = DCMPLX(XM,0.d0)
        Go To 180
c
c     case of a complex zero
c
  140   IT = 0
        ICV = 0
        CX = DCMPLX(XM,0.d0)
  150   FSG = (1.,0.)
        HSG = (1.,0.)
        FP = (0.,0.)
        FPP = (0.,0.)
        Do 160 J = 1, IZ
          DD = 1. / (CX-BH(J))
          FSG = FSG * A(J) * DD
          HSG = HSG * C(J) * DD
          FP = FP + DD
          FPP = FPP - DD * DD
  160   Continue
        If (MODIZ.EQ.0) Then
          F = (1.,0.) - FSG - HSG
        Else
          F = (1.,0.) + FSG + HSG
        End If
        I3 = 0
        If (CDABS(FP).GT.0) Then
          I3 = 1
          R3 = -F / FP
        End If
        I2 = 0
        If (CDABS(FPP).GT.0) Then
          I2 = 1
          CDIS = CDSQRT(FP**2-2.*F*FPP)
          R1 = CDIS - FP
          R2 = -FP - CDIS
          If (CDABS(R1).GT.CDABS(R2)) Then
            R1 = R1 / FPP
          Else
            R1 = R2 / FPP
          End If
          R2 = 2. * F / FPP / R1
          If (CDABS(R2).LT.CDABS(R1)) R1 = R2
          If (I3.LE.0) Go To 170
          If (CDABS(R3).LT.CDABS(R1)) R1 = R3
        Else
          R1 = R3
        End If
  170   CX = CX + R1
        IT = IT + 1
        If (IT.GT.50) Go To 240
        If (CDABS(R1).GT.SCNV) Go To 150
        If (ICV.LE.0) Then
          ICV = 1
          Go To 150
        End If
        CBP(IG) = CX
        CBP(IG+1) = CONJG(CX)
  180 Continue
      If (CDABS(CBP(N))-CDABS(CBP(1))) 190, 240, 210
  190 NHALF = N / 2
      Do 200 J = 1, NHALF
        NT = N - J
        CX = CBP(J)
        CBP(J) = CBP(NT+1)
        CBP(NT+1) = CX
  200 Continue
  210 NCMPLX = 1
      Do 220 J = 2, IZ
        If (DIMAG(CBP(J)).NE.0) Go To 250
  220 Continue
      NCMPLX = 0
      Do 230 J = 2, IZ
        BP(J) = DREAL(CBP(J))
  230 Continue
      Go To 250
  240 IERROR = 4
  250 Continue
      Return
      End





      Function PSGF(X,IZ,C,A,BH)
      IMPLICIT NONE
      DOUBLE PRECISION PSGF
      Dimension A(1), C(1), BH(1)
      real*8 A, X, C, BH
      integer IZ

      real*8 FSG, HSG, DD
      integer J

      FSG = 1.
      HSG = 1.
      Do 10 J = 1, IZ
        DD = 1. / (X-BH(J))
        FSG = FSG * A(J) * DD
        HSG = HSG * C(J) * DD
   10 Continue
      If (MOD(IZ,2).EQ.0) Then
        PSGF = 1. - FSG - HSG
        Return
      End If
      PSGF = 1. + FSG + HSG
      Return
      End










      Function BSRH(XLL,XRR,IZ,C,A,BH,F,SGN)
      IMPLICIT NONE
      DOUBLE PRECISION BSRH
      DOUBLE PRECISION F
      Dimension A(1), C(1), BH(1)
      real*8 XLL, XRR, A, C, BH, SGN
      integer IZ

      integer K, NM, NCMPLX, IK
      real*8 EPS, NPP, XL, CNV
      Common /CBLKT/ K, NM, NCMPLX, IK, NPP, CNV, EPS

      real*8 DX, X, XR

      XL = XLL
      XR = XRR
      DX = .5 * ABS(XR-XL)
   10 X = .5 * (XL+XR)
      If (SGN*F(X,IZ,C,A,BH)) 30, 50, 20
   20 XR = X
      Go To 40
   30 XL = X
   40 DX = .5 * DX
      If (DX.GT.CNV) Go To 10
   50 BSRH = .5 * (XL+XR)
      Return
      End








      Function PPSGF(X,IZ,C,A,BH)
      IMPLICIT NONE
      real*8 PPSGF
      Dimension A(1), C(1), BH(1)
      real*8 A, C, BH, X
      integer J, IZ

      real*8 SUM
      SUM = 0.
      Do 10 J = 1, IZ
        SUM = SUM - 1. / (X-BH(J)) ** 2
   10 Continue
      PPSGF = SUM
      Return
      End







      Function PPSPF(X,IZ,C,A,BH)
      IMPLICIT NONE
      real*8 PPSPF
      Dimension A(1), C(1), BH(1)
      real*8 SUM, A, C, BH, X
      integer J, IZ
      SUM = 0.
      Do 10 J = 1, IZ
        SUM = SUM + 1. / (X-BH(J))
   10 Continue
      PPSPF = SUM
      Return
      End








      Subroutine COMPB(N,IERROR,AN,BN,CN,B,AH,BH)
c
c     compb computes the roots of the b polynomials using subroutine
c     tevls which is a modification the eispack program tqlrat.
c     ierror is set to 4 if either tevls fails or if a(j+1)*c(j) is
c     less than zero for some j.  ah,bh are temporary work arrays.
c
      IMPLICIT NONE
      Dimension AN(1), BN(1), CN(1), B(1), AH(1), BH(1)
      real*8 AN, BN, CN, B, AH, BH
      integer N, IERROR
      real*8 BNORM, EPMACH


 

      integer J, K, NM, NCMPLX, IK, DUM, IF, KDO, L, IR, I2, I4, IPL 
      integer IFD, I, IB, NB, JS, JF, LS, LH, NMP, L1, L2, J2, J1, N2M2
      real*8 EPS, NPP, CNV, ARG, D1, D2, D3
      Common /CBLKT/ K, NM, NCMPLX, IK, NPP, CNV, EPS

      EPS = EPMACH(DUM)
      BNORM = DABS(BN(1))
      Do 10 J = 2, NM
        BNORM = DMAX1(BNORM,DABS(BN(J)))
        ARG = AN(J) * CN(J-1)
        If (ARG.LT.0) Go To 110
        B(J) = SIGN(SQRT(ARG),AN(J))
   10 Continue
      CNV = EPS * BNORM
      IF = 2 ** K
      KDO = K - 1
      Do 50 L = 1, KDO
        IR = L - 1
        I2 = 2 ** IR
        I4 = I2 + I2
        IPL = I4 - 1
        IFD = IF - I4
        Do 40 I = I4, IFD, I4
          Call INDXB(I,L,IB,NB)
          If (NB.LE.0) Go To 50
          JS = I - IPL
          JF = JS + NB - 1
          LS = 0
          Do 20 J = JS, JF
            LS = LS + 1
            BH(LS) = BN(J)
            AH(LS) = B(J)
   20     Continue
          Call TEVLS(NB,BH,AH,IERROR)
          If (IERROR.NE.0) Go To 100
          LH = IB - 1
          Do 30 J = 1, NB
            LH = LH + 1
            B(LH) = -BH(J)
   30     Continue
   40   Continue
   50 Continue
      Do 60 J = 1, NM
        B(J) = -BN(J)
   60 Continue
      If (NPP.EQ.0) Then
        NMP = NM + 1
        NB = NM + NMP
        Do 70 J = 1, NB
          L1 = MOD(J-1,NMP) + 1
          L2 = MOD(J+NM-1,NMP) + 1
          ARG = AN(L1) * CN(L2)
          If (ARG.LT.0) Go To 110
          BH(J) = SIGN(SQRT(ARG),-AN(L1))
          AH(J) = -BN(L1)
   70   Continue
        Call TEVLS(NB,AH,BH,IERROR)
        If (IERROR.NE.0) Go To 100
        Call INDXB(IF,K-1,J2,LH)
        Call INDXB(IF/2,K-1,J1,LH)
        J2 = J2 + 1
        LH = J2
        N2M2 = J2 + NM + NM - 2
   80   D1 = DABS(B(J1)-B(J2-1))
        D2 = DABS(B(J1)-B(J2))
        D3 = DABS(B(J1)-B(J2+1))
        If (.NOT.((D2.LT.D1).AND.(D2.LT.D3))) Then
          B(LH) = B(J2)
          J2 = J2 + 1
          LH = LH + 1
          If (J2-N2M2) 80, 80, 90
        End If
        J2 = J2 + 1
        J1 = J1 + 1
        If (J2.LE.N2M2) Go To 80
   90   B(LH) = B(N2M2+1)
        Call INDXB(IF,K-1,J1,J2)
        J2 = J1 + NMP + NMP
        Call PPADD(NM+1,IERROR,AN,CN,B(J1),B(J1),B(J2))
      End If
      Return
  100 IERROR = 4
      Return
  110 IERROR = 5
      Return
      End








      Subroutine TEVLS(N,D,E2,IERR)
c
      IMPLICIT NONE
      Integer I, J, L, M, N, II, L1, MML, IERR
      Real*8 D(N), E2(N)
      Real*8 B, C, F, G, H, P, R, S
c
c     real sqrt,abs,sign
c

      integer K, NM, NCMPLX, IK, NHALF, NTOP
      real*8 MACHEP, NPP, CNV, DHOLD
      Common /CBLKT/ K, NM, NCMPLX, IK, NPP, CNV, MACHEP

c
c     this subroutine is a modification of the eispack subroutine tqlrat
c     algorithm 464, comm. acm 16, 689(1973) by reinsch.
c
c     this subroutine finds the eigenvalues of a symmetric
c     tridiagonal matrix by the rational ql method.
c
c     on input-
c
c        n is the order of the matrix,
c
c        d contains the diagonal elements of the input matrix,
c
c        e2 contains the                subdiagonal elements of the
c          input matrix in its last n-1 positions.  e2(1) is arbitrary.
c
c      on output-
c
c        d contains the eigenvalues in ascending order.  if an
c          error exit is made, the eigenvalues are correct and
c          ordered for indices 1,2,...ierr-1, but may not be
c          the smallest eigenvalues,
c
c        e2 has been destroyed,
c
c        ierr is set to
c          zero       for normal return,
c          j          if the j-th eigenvalue has not been
c                     determined after 30 iterations.
c
c     questions and comments should be directed to b. s. garbow,
c     applied mathematics division, argonne national laboratory
c
c
c     ********** machep is a machine dependent parameter specifying
c                the relative precision of floating point arithmetic.
c
c                **********
c
      IERR = 0
      If (N.NE.1) Then
c
        Do 10 I = 2, N
          E2(I-1) = E2(I) * E2(I)
   10   Continue
c
        F = 0.0
        B = 0.0
        E2(N) = 0.0
c
        Do 90 L = 1, N
          J = 0
          H = MACHEP * (DABS(D(L))+DSQRT(E2(L)))
          If (B.LE.H) Then
            B = H
            C = B * B
          End If
c
c     ********** look for small squared sub-diagonal element **********
c
          Do 20 M = L, N
            If (E2(M).LE.C) Go To 30
c
c     ********** e2(n) is always zero, so there is no exit
c                through the bottom of the loop **********
c
   20     Continue
c
   30     If (M.NE.L) Then
   40       If (J.EQ.30) Go To 110
            J = J + 1
c
c     ********** form shift **********
c
            L1 = L + 1
            S = SQRT(E2(L))
            G = D(L)
            P = (D(L1)-G) / (2.0*S)
            R = SQRT(P*P+1.0)
            D(L) = S / (P+SIGN(R,P))
            H = G - D(L)
c
            Do 50 I = L1, N
              D(I) = D(I) - H
   50       Continue
c
            F = F + H
c
c     ********** rational ql transformation **********
c
            G = D(M)
            If (G.EQ.0.0) G = B
            H = G
            S = 0.0
            MML = M - L
c
c     ********** for i=m-1 step -1 until l do -- **********
c
            Do 60 II = 1, MML
              I = M - II
              P = G * H
              R = P + E2(I)
              E2(I+1) = S * R
              S = E2(I) / R
              D(I+1) = H + S * (H+D(I))
              G = D(I) - E2(I) / G
              If (G.EQ.0.0) G = B
              H = G * P / R
   60       Continue
c
            E2(L) = S * G
            D(L) = H
c
c     ********** guard against underflowed h **********
c
            If (H.NE.0.0) Then
              If (ABS(E2(L)).GT.ABS(C/H)) Then
                E2(L) = H * E2(L)
                If (E2(L).NE.0.0) Go To 40
              End If
            End If
          End If
          P = D(L) + F
c
c     ********** order eigenvalues **********
c
          If (L.NE.1) Then
c
c     ********** for i=l step -1 until 2 do -- **********
c
            Do 70 II = 2, L
              I = L + 2 - II
              If (P.GE.D(I-1)) Go To 80
              D(I) = D(I-1)
   70       Continue
          End If
c
          I = 1
   80     D(I) = P
   90   Continue
c
        If (ABS(D(N)).GE.ABS(D(1))) Go To 120
        NHALF = N / 2
        Do 100 I = 1, NHALF
          NTOP = N - I
          DHOLD = D(I)
          D(I) = D(NTOP+1)
          D(NTOP+1) = DHOLD
  100   Continue
        Go To 120
c
c     ********** set error -- no convergence to an
c                eigenvalue after 30 iterations **********
c
  110   IERR = L
      End If
  120 Return
c
c     ********** last card of tqlrat **********
c
c
c revision history---
c
c december 1979    first added to nssl
c-----------------------------------------------------------------------
      End
c package comf           the entries in this package are lowlevel
c                        entries, supporting ulib packages blktri
c                        and cblktri. that is, these routines are
c                        not called directly by users, but rather
c                        by entries within blktri and cblktri.
c                        description of entries epmach and pimach
c                        follow below.
c
c latest revision        january 1985
c
c special conditions     none
c
c i/o                    none
c
c precision              single
c
c required library       none
c files
c
c language               fortran
c ********************************************************************
c
c function epmach (dum)
c
c purpose                to compute an approximate machine accuracy
c                        epsilon according to the following definition:
c                        epsilon is the smallest number such that
c                        (1.+epsilon).gt.1.)
c
c usage                  eps = epmach (dum)
c
c arguments
c on input               dum
c                          dummy value
c
c arguments
c on output              none
c
c history                the original version, written when the
c                        blktri package was converted from the
c                        cdc 7600 to run on the cray-1, calculated
c                        machine accuracy by successive divisions
c                        by 10.  use of this constant caused blktri
c                        to compute solutions on the cray-1 with four
c                        fewer places of accuracy than the version
c                        on the 7600.  it was found that computing
c                        machine accuracy by successive divisions
c                        of 2 produced a machine accuracy 29% less
c                        than the value generated by successive
c                        divisions by 10, and that use of this
c                        machine constant in the blktri package
c                        recovered the accuracy that appeared to
c                        be lost on conversion.
c
c algorithm              computes machine accuracy by successive
c                        divisions of two.
c
c portability            this code will execute on machines other
c                        than the cray1, but the returned value may
c                        be unsatisfactory.  see history above.
c ********************************************************************
c
c function pimach (dum)
c
c purpose                to supply the value of the constant pi
c                        correct to machine precision where
c                        pi=3.141592653589793238462643383279502884197
c                             1693993751058209749446
c
c usage                  pi = pimach (dum)
c
c arguments
c on input               dum
c                          dummy value
c
c arguments
c on output              none
c
c algorithm              the value of pi is set in a constant.
c
c portability            this entry is portable, but users should
c                        check to see whether greater accuracy is
c                        required.
c
c***********************************************************************
      Function EPMACH(DUM)
      IMPLICIT NONE
      real*8 EPMACH, V, EPS
      integer DUM
      Common /VALUE/ V
      EPS = 1.
   10 EPS = EPS / 2.
      Call STORE(EPS+1.)
      If (V.GT.1.) Go To 10
      EPMACH = 100. * EPS
      Return
      End
      Subroutine STORE(X)
      IMPLICIT NONE
      real*8 X, V
      Common /VALUE/ V
      V = X
      Return
      End







C      Function PIMACH(DUM)
C      IMPLICIT NONE
C      real*8 PIMACH
c     pi=3.1415926535897932384626433832795028841971693993751058209749446
c
C      PIMACH = 3.14159265358979
C      Return
C      End
c package sepaux         contains no user entry points.
c
c latest revision        march 1985
c
c purpose                this package contains auxiliary routines for
c                        ncar public software packages such as sepeli
c                        and sepx4.
c
c usage                  since this package contains no user entries,
c                        no usage instructions or argument descriptions
c                        are given here.
c
c special conditions     none
c
c i/o                    none
c
c precision              single
c
c required library       none
c files
c
c language               fortran
c
c history                developed in the late 1970's by john c. adams
c                        of ncar's scienttific computing division.
c
c portability            fortran 66
c **********************************************************************
      Subroutine SEPORT(USOL,IDMN,ZN,ZM,PERTRB)
c
c     this subroutine orthoganalizes the array usol with respect to
c     the constant array in a weighted least squares norm
c
      IMPLICIT NONE
      integer KSWX, KSWY, K , L
      integer AIT, BIT, CIT, DIT, IS, MIT, NIT
      real*8 DLX, DLY, TDLX3, TDLY3, DLX4, DLY4
      integer MS, JS, NS
      Common /SPLP/ KSWX, KSWY, K, L, 
     *    DLX, DLY, TDLX3, TDLY3, DLX4, DLY4,
     *    AIT, BIT, CIT, DIT, MIT, NIT, IS, MS, JS, NS
      integer IDMN
      Dimension USOL(IDMN,1), ZN(1), ZM(1)
      real*8 USOL, ZN, ZM, UTE, ETE, PERTRB
      integer ISTR, JSTR, IFNL, JFNL, I, II, J, JJ


      ISTR = IS
      IFNL = MS
      JSTR = JS
      JFNL = NS
c
c     compute weighted inner products
c
      UTE = 0.0
      ETE = 0.0
      Do 20 I = IS, MS
        II = I - IS + 1
        Do 10 J = JS, NS
          JJ = J - JS + 1
          ETE = ETE + ZM(II) * ZN(JJ)
          UTE = UTE + USOL(I,J) * ZM(II) * ZN(JJ)
   10   Continue
   20 Continue
c
c     set perturbation parameter
c
      PERTRB = UTE / ETE
c
c     subtract off constant pertrb
c
      Do 40 I = ISTR, IFNL
        Do 30 J = JSTR, JFNL
          USOL(I,J) = USOL(I,J) - PERTRB
   30   Continue
   40 Continue
      Return
      End








      Subroutine SEPMIN(USOL,IDMN,ZN,ZM,PERTB)
c
c     this subroutine orhtogonalizes the array usol with respect to
c     the constant array in a weighted least squares norm
c
      IMPLICIT NONE
      integer KSWX, KSWY, K, IDMN, ISTR, JSTR,L, IFNL, JFNL
      integer AIT, BIT, CIT, DIT, IS, MIT, NIT
      real*8 UTE, ETE, PERTB, PERTRB
      integer I, II, J, JJ
      real*8 DLX, DLY, TDLX3, TDLY3, DLX4, DLY4
      integer MS, JS, NS
      Common /SPLP/ KSWX, KSWY, K, L, 
     *    DLX, DLY, TDLX3, TDLY3, DLX4, DLY4,
     *    AIT, BIT, CIT, DIT, MIT, NIT, IS, MS, JS, NS
      Dimension USOL(IDMN,1), ZN(1), ZM(1)
      real*8 USOL, ZN, ZM
c
c     entry at sepmin occurrs when the final solution is
c     to be minimized with respect to the weighted
c     least squares norm
c
      ISTR = 1
      IFNL = K
      JSTR = 1
      JFNL = L
c
c     compute weighted inner products
c
      UTE = 0.0
      ETE = 0.0
      Do 20 I = IS, MS
        II = I - IS + 1
        Do 10 J = JS, NS
          JJ = J - JS + 1
          ETE = ETE + ZM(II) * ZN(JJ)
          UTE = UTE + USOL(I,J) * ZM(II) * ZN(JJ)
   10   Continue
   20 Continue
c
c     set perturbation parameter
c
      PERTRB = UTE / ETE
c
c     subtract off constant pertrb
c
      Do 40 I = ISTR, IFNL
        Do 30 J = JSTR, JFNL
          USOL(I,J) = USOL(I,J) - PERTRB
   30   Continue
   40 Continue
      Return
      End









      Subroutine SEPTRI(N,A,B,C,D,U,Z)
c
c     this subroutine solves for a non-zero eigenvector corresponding
c     to the zero eigenvalue of the transpose of the rank
c     deficient one matrix with subdiagonal a, diagonal b, and
c     superdiagonal c , with a(1) in the (1,n) position, with
c     c(n) in the (n,1) position, and all other elements zero.
c
      IMPLICIT NONE
      integer N
      Dimension A(N), B(N), C(N), D(N), U(N), Z(N)
      real*8 A, B, C, D, U, Z

      real*8 BN, V, DEN, AN
      integer K, NM1, NM2, J
      BN = B(N)
      D(1) = A(2) / B(1)
      V = A(1)
      U(1) = C(N) / B(1)
      NM2 = N - 2
      Do 10 J = 2, NM2
        DEN = B(J) - C(J-1) * D(J-1)
        D(J) = A(J+1) / DEN
        U(J) = -C(J-1) * U(J-1) / DEN
        BN = BN - V * U(J-1)
        V = -V * D(J-1)
   10 Continue
      DEN = B(N-1) - C(N-2) * D(N-2)
      D(N-1) = (A(N)-C(N-2)*U(N-2)) / DEN
      AN = C(N-1) - V * D(N-2)
      BN = BN - V * U(N-2)
      DEN = BN - AN * D(N-1)
c
c     set last component equal to one
c
      Z(N) = 1.0
      Z(N-1) = -D(N-1)
      NM1 = N - 1
      Do 20 J = 2, NM1
        K = N - J
        Z(K) = -D(K) * Z(K+1) - U(K) * Z(N)
   20 Continue
      Return
      End

























      Subroutine SEPDX(U,IDMN,I,J,UXXX,UXXXX)
c
c     this program computes second order finite difference
c     approximations to the third and fourth x
c     partial derivatives of u at the (i,j) mesh point
c
      IMPLICIT NONE
      integer KSWX, KSWY, K, L, I, J
      integer AIT, BIT, CIT, DIT, IS, MIT, NIT
      real*8 DLX, DLY, TDLX3, TDLY3, DLX4, DLY4
      integer MS, JS, NS
      Common /SPLP/ KSWX, KSWY, K, L, 
     *    DLX, DLY, TDLX3, TDLY3, DLX4, DLY4,
     *    AIT, BIT, CIT, DIT, MIT, NIT, IS, MS, JS, NS
      integer IDMN
      Dimension U(IDMN,1)
      real*8 U, UXXX, UXXXX
      If (I.LE.2.OR.I.GE.(K-1)) Then
        If (I.NE.1) Then
          If (I.EQ.2) Go To 10
          If (I.EQ.K-1) Go To 20
          If (I.EQ.K) Go To 30
        End If
c
c     compute partial derivative approximations at x=a
c
        If (KSWX.NE.1) Then
          UXXX = (-5.0*U(1,J)+18.0*U(2,J)-24.0*U(3,J)+14.0*U(4,J)-3.0*
     *        U(5,J)) / (TDLX3)
          UXXXX = (3.0*U(1,J)-14.0*U(2,J)+26.0*U(3,J)-24.0*U(4,J)+11.0*
     *        U(5,J)-2.0*U(6,J)) / DLX4
          Return
        End If
c
c     periodic at x=a
c
        UXXX = (-U(K-2,J)+2.0*U(K-1,J)-2.0*U(2,J)+U(3,J)) / (TDLX3)
        UXXXX = (U(K-2,J)-4.0*U(K-1,J)+6.0*U(1,J)-4.0*U(2,J)+U(3,J)) / 
     *      DLX4
        Return
c
c     compute partial derivative approximations at x=a+dlx
c
   10   If (KSWX.NE.1) Then
          UXXX = (-3.0*U(1,J)+10.0*U(2,J)-12.0*U(3,J)+6.0*U(4,J)-
     *        U(5,J)) / TDLX3
          UXXXX = (2.0*U(1,J)-9.0*U(2,J)+16.0*U(3,J)-14.0*U(4,J)+6.0*
     *        U(5,J)-U(6,J)) / DLX4
          Return
        End If
c
c     periodic at x=a+dlx
c
        UXXX = (-U(K-1,J)+2.0*U(1,J)-2.0*U(3,J)+U(4,J)) / (TDLX3)
        UXXXX = (U(K-1,J)-4.0*U(1,J)+6.0*U(2,J)-4.0*U(3,J)+U(4,J)) / 
     *      DLX4
        Return
      End If
c
c     compute partial derivative approximations on the interior
c
      UXXX = (-U(I-2,J)+2.0*U(I-1,J)-2.0*U(I+1,J)+U(I+2,J)) / TDLX3
      UXXXX = (U(I-2,J)-4.0*U(I-1,J)+6.0*U(I,J)-4.0*U(I+1,J)+U(I+2,J)) /
     *    DLX4
      Return
c
c     compute partial derivative approximations at x=b-dlx
c
   20 If (KSWX.NE.1) Then
        UXXX = (U(K-4,J)-6.0*U(K-3,J)+12.0*U(K-2,J)-10.0*U(K-1,J)+3.0*
     *      U(K,J)) / TDLX3
        UXXXX = (-U(K-5,J)+6.0*U(K-4,J)-14.0*U(K-3,J)+16.0*U(K-2,J)-9.0*
     *      U(K-1,J)+2.0*U(K,J)) / DLX4
        Return
      End If
c
c     periodic at x=b-dlx
c
      UXXX = (-U(K-3,J)+2.0*U(K-2,J)-2.0*U(1,J)+U(2,J)) / TDLX3
      UXXXX = (U(K-3,J)-4.0*U(K-2,J)+6.0*U(K-1,J)-4.0*U(1,J)+U(2,J)) / 
     *    DLX4
      Return
c
c     compute partial derivative approximations at x=b
c
   30 UXXX = -(3.0*U(K-4,J)-14.0*U(K-3,J)+24.0*U(K-2,J)-18.0*U(K-1,J)+
     *    5.0*U(K,J)) / TDLX3
      UXXXX = (-2.0*U(K-5,J)+11.0*U(K-4,J)-24.0*U(K-3,J)+26.0*U(K-2,J)-
     *    14.0*U(K-1,J)+3.0*U(K,J)) / DLX4
      Return
      End






















      Subroutine SEPDY(U,IDMN,I,J,UYYY,UYYYY)
      IMPLICIT NONE
c
c     this program computes second order finite difference
c     approximations to the third and fourth y
c     partial derivatives of u at the (i,j) mesh point
c
      integer KSWX, KSWY, K 
      integer AIT, BIT, CIT, DIT, IS, MIT, NIT
      real*8 DLX, DLY, TDLX3, TDLY3, DLX4, DLY4
      integer MS, JS, NS
      Common /SPLP/ KSWX, KSWY, K, L, 
     *    DLX, DLY, TDLX3, TDLY3, DLX4, DLY4,
     *    AIT, BIT, CIT, DIT, MIT, NIT, IS, MS, JS, NS
      integer IDMN, I, J, L
      Dimension U(IDMN,IDMN)
      real*8 U, UYYY, UYYYY

      If (J.LE.2.OR.J.GE.(L-1)) Then
        If (J.NE.1) Then
          If (J.EQ.2) Go To 10
          If (J.EQ.L-1) Go To 20
          If (J.EQ.L) Go To 30
        End If
c
c     compute partial derivative approximations at y=c
c
        If (KSWY.NE.1) Then
          UYYY = (-5.0*U(I,1)+18.0*U(I,2)-24.0*U(I,3)+14.0*U(I,4)-3.0*
     *        U(I,5)) / TDLY3
          UYYYY = (3.0*U(I,1)-14.0*U(I,2)+26.0*U(I,3)-24.0*U(I,4)+11.0*
     *        U(I,5)-2.0*U(I,6)) / DLY4
          Return
        End If
c
c     periodic at x=a
c
        UYYY = (-U(I,L-2)+2.0*U(I,L-1)-2.0*U(I,2)+U(I,3)) / TDLY3
        UYYYY = (U(I,L-2)-4.0*U(I,L-1)+6.0*U(I,1)-4.0*U(I,2)+U(I,3)) / 
     *      DLY4
        Return
c
c     compute partial derivative approximations at y=c+dly
c
   10   If (KSWY.NE.1) Then
          UYYY = (-3.0*U(I,1)+10.0*U(I,2)-12.0*U(I,3)+6.0*U(I,4)-
     *        U(I,5)) / TDLY3
          UYYYY = (2.0*U(I,1)-9.0*U(I,2)+16.0*U(I,3)-14.0*U(I,4)+6.0*
     *        U(I,5)-U(I,6)) / DLY4
          Return
        End If
c
c     periodic at y=c+dly
c
        UYYY = (-U(I,L-1)+2.0*U(I,1)-2.0*U(I,3)+U(I,4)) / TDLY3
        UYYYY = (U(I,L-1)-4.0*U(I,1)+6.0*U(I,2)-4.0*U(I,3)+U(I,4)) / 
     *      DLY4
        Return
      End If
c
c     compute partial derivative approximations on the interior
c
      UYYY = (-U(I,J-2)+2.0*U(I,J-1)-2.0*U(I,J+1)+U(I,J+2)) / TDLY3
      UYYYY = (U(I,J-2)-4.0*U(I,J-1)+6.0*U(I,J)-4.0*U(I,J+1)+U(I,J+2)) /
     *    DLY4
      Return
c
c     compute partial derivative approximations at y=d-dly
c
   20 If (KSWY.NE.1) Then
        UYYY = (U(I,L-4)-6.0*U(I,L-3)+12.0*U(I,L-2)-10.0*U(I,L-1)+3.0*
     *      U(I,L)) / TDLY3
        UYYYY = (-U(I,L-5)+6.0*U(I,L-4)-14.0*U(I,L-3)+16.0*U(I,L-2)-9.0*
     *      U(I,L-1)+2.0*U(I,L)) / DLY4
        Return
      End If
c
c     periodic at y=d-dly
c
      UYYY = (-U(I,L-3)+2.0*U(I,L-2)-2.0*U(I,1)+U(I,2)) / TDLY3
      UYYYY = (U(I,L-3)-4.0*U(I,L-2)+6.0*U(I,L-1)-4.0*U(I,1)+U(I,2)) / 
     *    DLY4
      Return
c
c     compute partial derivative approximations at y=d
c
   30 UYYY = -(3.0*U(I,L-4)-14.0*U(I,L-3)+24.0*U(I,L-2)-18.0*U(I,L-1)+
     *    5.0*U(I,L)) / TDLY3
      UYYYY = (-2.0*U(I,L-5)+11.0*U(I,L-4)-24.0*U(I,L-3)+26.0*U(I,L-2)-
     *    14.0*U(I,L-1)+3.0*U(I,L)) / DLY4
      Return
c
c revision history---
c
c december 1979    first added to nssl
c-----------------------------------------------------------------------
      End






      Subroutine COFX(XX,AFUN,BFUN,CFUN)
c
c     subroutine to compute the coefficients of the elliptic equation
c     solved to fill in the grid.  it is passed (along with cofy) to
c     the elliptic solver sepeli.
#include "mexsepeli.inc"
      integer I
      real*8 DXDI, AFUN, BFUN, CFUN, XX
      I = DNINT(XX)
      DXDI = 0.5 * (SXI(I+1)-SXI(I-1))
      AFUN = 1. / DXDI ** 2
      BFUN = (-SXI(I+1)+2.*SXI(I)-SXI(I-1)) / DXDI ** 3
      CFUN = 0.
      Return
      End




      Subroutine COFY(YY,DFUN,EFUN,FFUN)
#include "mexsepeli.inc"
      integer J
      real*8 DEDJ, DFUN, EFUN, FFUN, YY
      J = DNINT(YY)
      DEDJ = 0.5 * (SETA(J+1)-SETA(J-1))
      DFUN = 1. / DEDJ ** 2
      EFUN = (-SETA(J+1)+2.*SETA(J)-SETA(J-1)) / DEDJ ** 3
      FFUN = 0.
      Return
      End




