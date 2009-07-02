#include <fintrf.h>
c
c purpose       mexsepeli solves for the 2nd order finite 
c               difference approximation to the separable
c               elliptic equation arising from conformal
c               mapping and grid generation.
c
c
c
c 
c The calling sequence from MATLAB should be
c
c   >> [x, y, perturb, ierror] = mexsepeli (   x, y, ...
c                                              l2, m2, ...
c                                              seta, sxi, ...
c                                              nwrk, nx2 );
c
c
c Output Parameters
c     x,y:
c          represent "x" and "y" grid coordinates
c     perturb:
c          Optional.  Don't think this is needed for our 
c          application.
c     ierror:
c          Optional.  An error flag that indicates invalid 
c          input parameters or failure to find a solution.
c                          = 0 no error
c                          = 4 if attempt to find a solution fails.
c                              (the linear system generated is not
c                              diagonally dominant.)
c
c Input Parameters
c     x,y:
c          represent "x" and "y" grid coordinates
c     l2:
c          x independent variable ranges from 0 to l2
c     m2:
c          y independent variable ranges from 0 to m2
c     seta;
c          digitized points on boundaries 1,3
c     sxi:
c          digitized points on boundaries 2,4
c          
c
c So 
c    1)  nlhs = 2
c    2)  plhs(1) = x
c        plhs(2) = y
c    3)  nrhs = 6 
c    4)  prhs(1) = x 
c        prhs(2) = y 
c        prhs(3) = l2
c        prhs(4) = m2
c        prhs(5) = seta
c        prhs(6) = sxi
c


c
c svn $Id$
c Revision 1.2  2005/03/04 22:56:12  johnevans007
c A ton of minor changes.  Moved some common block variables around to make
c them fit nice.  Don't understand how the preprocessor was compiling before
c I changed the include statements.
c
c Revision 1.1.1.1  2004/11/09 21:14:02  johnevans007
c Initial import.
c
c Revision 1.3  2004/09/30 20:57:43  jevans
c Several automatic casts (from real*8 to complex) were being done that g77
c was complaining about.  Some new complex variables were added that made
c the cast explicit.  g77 no longer spits out any warnings.
c
c Revision 1.2  2004/09/30 20:32:39  jevans
c Removed all declarations of variables that the solaris compiler claimed
c were never used.  Seems to run fine.
c
c Revision 1.1.1.1  2004/09/30 15:22:15  jevans
c initial import
c
c Revision 1.1.1.1  2004/03/02 16:53:17  jevans
c MEXSEPELI
c
c Revision 1.4  2003/02/09 17:51:14  jevans
c Converted some more CMPLX statements to DCMPLX.  Linux compilation didn't
c catch it, but Alpha did.
c
c Revision 1.2  2003/01/15 23:23:04  jevans
c Lots of matrices redefined as for their sizes.  Who knows if it works?
c How the hell did it ever work before?  God how I luv fortran.
c
c Revision 1.1.1.1  2003/01/15 22:31:12  jevans
c Imported sources
c
c
c


c 
c     $Id$
c     Currently locked by $Locker:  $ (not locked if blank)
c     (Time in GMT, EST=GMT-5:00
c     $Log: mexsepeli.F,v $
c     Revision 1.2  2005/03/04 22:56:12  johnevans007
c     A ton of minor changes.  Moved some common block variables around to make
c     them fit nice.  Don't understand how the preprocessor was compiling before
c     I changed the include statements.
c
c     Revision 1.1.1.1  2004/11/09 21:14:02  johnevans007
c     Initial import.
c
c     Revision 1.3  2004/09/30 20:57:43  jevans
c     Several automatic casts (from real*8 to complex) were being done that g77
c     was complaining about.  Some new complex variables were added that made
c     the cast explicit.  g77 no longer spits out any warnings.
c
c     Revision 1.2  2004/09/30 20:32:39  jevans
c     Removed all declarations of variables that the solaris compiler claimed
c     were never used.  Seems to run fine.
c
c     Revision 1.1.1.1  2004/09/30 15:22:15  jevans
c     initial import
c
c     Revision 1.1.1.1  2004/03/02 16:53:17  jevans
c     MEXSEPELI
c
c     Revision 1.4  2003/02/09 17:51:14  jevans
c     Converted some more CMPLX statements to DCMPLX.  Linux compilation didn't
c     catch it, but Alpha did.
c
c     Revision 1.2  2003/01/15 23:23:04  jevans
c     Lots of matrices redefined as for their sizes.  Who knows if it works?
c     How the hell did it ever work before?  God how I luv fortran.
c
c     Revision 1.1.1.1  2003/01/15 22:31:12  jevans
c     Imported sources
c
c Revision 1.2  1997/10/03  19:43:04  jevans
c Had to change compile time options to include "-align dcommons"
c and "-r8".
c
c the "mexCopyReal8ToPtr" and related functions suffer from the
c fact that they assume all arrays to be a single column.  Since
c my x and y are dimensioned to be larger than their logical
c size (i.e. they are 800x800, but might logically be only 100x80
c or so), those arrays had to be columnified before passing them
c back to MATLAB, and "uncolumnified" upon passing them into
c MATLAB.
c
c Revision 1.1  1997/10/03  17:50:43  jevans
c Initial revision
c
c 
c  




      subroutine mexFunction ( nlhs, plhs, nrhs, prhs )
#ifdef __alpha__
      POINTER plhs(*), prhs(*)
#else
      integer plhs(*), prhs(*)
#endif
      integer nlhs, nrhs


#include "mexsepeli.inc"


c
c     These are the number of rows and columns (m by n) of the
c     matrices we work with.
      integer size_m, size_n, size_mn
      integer i, j, k


c
c     Need some temporary space to store the double-precision
c     matlab arguments.  The correct number of rows and columns
c     is not known at first.
	  real*8 tempx(0:nx2,0:ny2), tempy(0:nx2,0:ny2)


c
c     input ranges for independent variables.
c     0 < x < l2, 0 < y < y2
      integer l2, m2


c
c     Pointers to input arguments
#ifdef __alpha__
      POINTER p 
#else
      integer p 
#endif

      real*8 pa(1)
      real*8 PERTRB

c
c     Check for proper number of arguments.
c     There must be at least 1 argument on the left hand size and 
c     no more than 3.
      if ( nlhs .lt. 1 ) then
         call mexErrMsgTxt ( 'Not enough input args for mexsepeli.' )
      elseif ( nlhs .gt. 3 ) then
         call mexErrMsgTxt ( 'Too many input args for mexsepeli.' )
      endif


      if ( nrhs .ne. 6 ) then
         call mexErrMsgTxt('mexsepeli requires 6 input arguments')
      endif


c
c     Check that all input arguments were numeric.
      if ( mxIsNumeric(prhs(1)) .ne. 1 ) then
         call mexErrMsgTxt 
     +        ( 'x input array should be numeric' )
      endif
      if ( mxIsNumeric(prhs(2)) .ne. 1 ) then
         call mexErrMsgTxt 
     +        ( 'y input array should be numeric' )
      endif
      if ( mxIsNumeric(prhs(3)) .ne. 1 ) then
         call mexErrMsgTxt 
     +        ( 'l2 input should be numeric' )
      endif
      if ( mxIsNumeric(prhs(4)) .ne. 1 ) then
         call mexErrMsgTxt 
     +        ( 'm2 input should be numeric' )
      endif
      if ( mxIsNumeric(prhs(5)) .ne. 1 ) then
         call mexErrMsgTxt 
     +        ( 'seta input array should be numeric' )
      endif
      if ( mxIsNumeric(prhs(6)) .ne. 1 ) then
         call mexErrMsgTxt 
     +        ( 'sxi input array should be numeric' )
      endif


c
c     Check to see that the l2 and m2 arguments were scalars
      size_m = mxGetM ( prhs(3) )
      size_n = mxGetN ( prhs(3) )
	  if ( size_n .ne. 1 .or. size_m .ne. 1 ) then
         call mexErrMsgTxt
     +        ('l2 argument should be a scalar' )
      endif

      size_m = mxGetM ( prhs(4) )
      size_n = mxGetN ( prhs(4) )
	  if ( size_n .ne. 1 .or. size_m .ne. 1 ) then
         call mexErrMsgTxt
     +        ('m2 argument should be a scalar' )
      endif




c
c     Load the data into Fortran matrices.

c
c     x
      p = mxGetPr(prhs(1))
      size_m = mxGetM ( prhs(1) )
      size_n = mxGetN ( prhs(1) )
      size_mn = size_m*size_n
      call mxCopyPtrToReal8(p,tempx,size_mn)

      k = 0
      do 10 j = 0, size_n - 1
         do 5 i = 0, size_m - 1
            x(i,j) = tempx(k,0)
            k = k + 1 
    5    continue
   10 continue  
           


c
c     y
      p = mxGetPr(prhs(2))
      size_m = mxGetM ( prhs(2) )
      size_n = mxGetN ( prhs(2) )
      size_mn = size_m*size_n
      call mxCopyPtrToReal8(p,tempy,size_mn)

      k = 0
      do 20 j = 0, size_n - 1
         do 15 i = 0, size_m - 1
            y(i,j) = tempy(k,0)
            k = k + 1 
   15    continue
   20 continue  
           

c
c     l2
      p = mxGetPr(prhs(3))
      call mxCopyPtrToReal8(p,pa,1)
      l2 = pa(1)
      
c
c     m2
      p = mxGetPr(prhs(4))
      call mxCopyPtrToReal8(p,pa,1)
      m2 = pa(1)
      


c
c     seta
      p = mxGetPr(prhs(5))
      size_m = mxGetM ( prhs(5) )
      size_n = mxGetN ( prhs(5) )
      neta = size_m*size_n
      call mxCopyPtrToReal8(p,seta,neta)



c
c     sxi
      p = mxGetPr(prhs(6))
      size_m = mxGetM ( prhs(6) )
      size_n = mxGetN ( prhs(6) )
      nxi = size_m*size_n
      call mxCopyPtrToReal8(p,sxi,nxi)




c
c     Construct right hand side
      do 30 j = 0, neta - 1
        do 25 i = 0, nxi - 1
            rhs(i,j) = 0.
  25    continue
  30  continue


c
c     Call the computational subroutine.
      ewrk(1) = nwrk
      call sepeli(0,2,0.d0,dble(l2),l2,1,wrk,wrk,wrk,wrk,0.d0,
     *    dble(m2),m2,
     *    1,wrk,wrk,wrk,wrk,rhs,x,nx2+1,ewrk,pertrb,ierr)
      if ( ierr .eq. 1 ) then
         call mexErrMsgTxt
     +        ( 'range of independent variables is out of whack' )
      else if ( ierr .eq. 2 ) then
         call mexErrMsgTxt
     +        ( 'boundary condition mbdcnd wrongly specified' )
      else if ( ierr .eq. 3 ) then
         call mexErrMsgTxt
     +        ( 'boundary condition nbdcnd wrongly specified' )
      else if ( ierr .eq. 4 ) then
         call mexErrMsgTxt
     +        ( 'linear system generated is not diagonally dominant' )
      else if ( ierr .eq. 5 ) then
         call mexErrMsgTxt
     +        ( 'idmn is too small' )
      else if ( ierr .eq. 6 ) then
         call mexErrMsgTxt
     +        ( 'm is too small or too large' )
      else if ( ierr .eq. 7 ) then
         call mexErrMsgTxt
     +        ( 'n is too small or too large' )
      else if ( ierr .eq. 8 ) then
         call mexErrMsgTxt
     +        ( 'iorder is not 2 or 4' )
      else if ( ierr .eq. 9 ) then
         call mexErrMsgTxt
     +        ( 'intl is not 0 or 1' )
      else if ( ierr .eq. 10 ) then
         call mexErrMsgTxt
     +        ( 'afun*dfun less than or equal to 0' )
      else if ( ierr .eq. 11 ) then
         call mexErrMsgTxt
     +        ( 'work space length input in w(1) is not right' )
      end if
          



      ewrk(1) = nwrk
      call sepeli(0,2,0.d0,dble(l2),l2,1,wrk,wrk,wrk,wrk,0.d0,
     *    dble(m2),m2,
     *    1,wrk,wrk,wrk,wrk,rhs,y,nx2+1,ewrk,pertrb,ierr)
      if ( ierr .eq. 1 ) then
         call mexErrMsgTxt
     +        ( 'range of independent variables is out of whack' )
      else if ( ierr .eq. 2 ) then
         call mexErrMsgTxt
     +        ( 'boundary condition mbdcnd wrongly specified' )
      else if ( ierr .eq. 3 ) then
         call mexErrMsgTxt
     +        ( 'boundary condition nbdcnd wrongly specified' )
      else if ( ierr .eq. 4 ) then
         call mexErrMsgTxt
     +        ( 'linear system generated is not diagonally dominant' )
      else if ( ierr .eq. 5 ) then
         call mexErrMsgTxt
     +        ( 'idmn is too small' )
      else if ( ierr .eq. 6 ) then
         call mexErrMsgTxt
     +        ( 'm is too small or too large' )
      else if ( ierr .eq. 7 ) then
         call mexErrMsgTxt
     +        ( 'n is too small or too large' )
      else if ( ierr .eq. 8 ) then
         call mexErrMsgTxt
     +        ( 'iorder is not 2 or 4' )
      else if ( ierr .eq. 9 ) then
         call mexErrMsgTxt
     +        ( 'intl is not 0 or 1' )
      else if ( ierr .eq. 10 ) then
         call mexErrMsgTxt
     +        ( 'afun*dfun less than or equal to 0' )
      else if ( ierr .eq. 11 ) then
         call mexErrMsgTxt
     +        ( 'work space length input in w(1) is not right' )
      end if




c
c     Create matrices for output arguments
c     Since the mxCopyReal8ToPtr routine just assumes one long
c     column matrix, we have to "columnify" both x and y.

c
c     x
      size_m = mxGetM ( prhs(1) )
      size_n = mxGetN ( prhs(1) )
      size_mn = size_m*size_n
      plhs(1) = mxCreateFull(size_m, size_n, 0)
      p = mxGetPr(plhs(1))

      k = 0
      do 40 j = 0, size_n - 1
         do 35 i = 0, size_m - 1
            tempx(k,0) = x(i,j)
            k = k + 1 
   35    continue
   40 continue  
           
      call mxCopyReal8ToPtr ( tempx, p, size_mn )


c
c     y
      size_m = mxGetM ( prhs(2) )
      size_n = mxGetN ( prhs(2) )
      size_mn = size_m*size_n
      plhs(2) = mxCreateFull(size_m, size_n, 0)
      p = mxGetPr(plhs(2))

      k = 0
      do 50 j = 0, size_n - 1
         do 45 i = 0, size_m - 1
            tempy(k,0) = y(i,j)
            k = k + 1 
   45    continue
   50 continue  
           
      call mxCopyReal8ToPtr ( tempy, p, size_mn )



      return
      end





