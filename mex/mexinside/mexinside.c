/*
* Mex routine mexinside.
* Given a point x,y and the series xb(k),yb(k) (k=1...nb) defining
* vertices of a closed polygon.  ind is set to 1 if the point is in
* the polygon and 0 if outside.  each time a new set of bound points
* is introduced ind should be set to 999 on input.
* it is best to do a series of y for a single fixed x.
* method ... a count is made of the no. of times the boundary cuts
* the meridian thru (x,y) south of (x,y).   an odd count indicates
* the point is inside , even indicates outside.
* see a long way from euclid by constance reid  p 174 .
* oceanography emr   oct/69
*
* The calling sequence from MATLAB should be
*
*   >> inside_inds = mexinside ( x, y, xb, yb );
*
*
* So 
*   1)  nlhs = 1
*   2)  plhs = 
*   3)  nrhs = 4
*   4)  prhs(0) = x, 'x' coordinates of candidate points
*       prhs(1) = y, 'y' coordinates of candidate points
*       prhs(2) = xb, 'x' coordinates of polygon
*       prhs(3) = yb, 'y' coordinates of polygon
*/

/*
 *  svn $Id$
 */

#include "mex.h"
#include "math.h"

void inside1 ( double *x, double *y, int x_m, int x_n,
	       double *xb, double *yb, int xb_m, int xb_n, 
	       double *inds );

/*
 * The gateway function
 */
void mexFunction ( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {
      

    int
	x_m, x_n,     /* size of x input matrix */
	y_m, y_n,     /* size of y input matrix */
	xb_m, xb_n,   /* size of xb boundary matrix */
	yb_m, yb_n;   /* size of yb boundary matrix */

    double
	*x, *y,       /* pointers to candidate points coordinates */
	*xb, *yb,     /* pointers to boundary points */
	*inds;        /* pointer to return list */


    if ( nrhs != 4 ) {
	mexErrMsgTxt ( "mexinside requires 4 input arguments" );
    }

    if ( nlhs != 1 ) {
	mexErrMsgTxt ( "mexinside requires 1 output argument" );
    }


    /* 
     * Check to make sure that the input arguments were numeric.  
     */
    if ( !mxIsNumeric ( prhs[0] ) ||
         !mxIsNumeric ( prhs[1] ) ||
         !mxIsNumeric ( prhs[2] ) ||
         !mxIsNumeric ( prhs[3] ) ) {

	 mexErrMsgTxt ( "All input arguments must be numeric." );
    }

    /*
     * Get the sizes of the inputs to make sure that they
     * are consistent.
     */
    x_m = mxGetM ( prhs[0] );
    x_n = mxGetN ( prhs[0] );
    y_m = mxGetM ( prhs[1] );
    y_n = mxGetN ( prhs[1] );
    if ( (x_m != y_m) || (x_n != y_n) ) {
	mexErrMsgTxt ( "First two input arguments must have the same dimensions." );
    }

    xb_m = mxGetM ( prhs[2] );
    xb_n = mxGetN ( prhs[2] );
    yb_m = mxGetM ( prhs[3] );
    yb_n = mxGetN ( prhs[3] );
    if ( (x_m != y_m) || (x_n != y_n) ) {
	mexErrMsgTxt ( "Last two input arguments must have the same dimensions." );
    }


    /*
     * Set the size of the output matrix.  Since the output
     * is a zero is a point is outside and a 1 if the point
     * is inside, this matrix is the same size as x and y.
     */
    plhs[0] = mxCreateDoubleMatrix ( x_m, x_n, mxREAL );


    /*
     * Set pointers to the data (finally!)
     */
    x = mxGetPr ( prhs[0] );
    y = mxGetPr ( prhs[1] );
    xb = mxGetPr ( prhs[2] );
    yb = mxGetPr ( prhs[3] );
    inds = mxGetPr ( plhs[0] );

    /*
     * Computational stuff here.
     */

     /*
      * inside1
      */
    inside1 ( x, y, x_m, x_n, xb, yb, xb_m, yb_n, inds );

    /*
     * Last thing to do is to load the result.
     */
    mxSetPr ( plhs[0], inds );

}



void inside1 ( double *x, double *y, int x_m, int x_n,
	       double *xb, double *yb, int xb_m, int xb_n, 
	       double *inds ) {

    int
	nxy,          /* number of candidate points */
	nb,           /* number of boundary points */

	k;            /* loop index for each point
		       * x[k], y[k], k = 0..nxy
		       */

    int 
	i, j, c;      /* from the wustle graphics algorithm */

    nb = xb_m * xb_n;
    nxy = x_m * x_n;

    /*
     * For each point in x[] and y[], ...
     */
    for ( k = 0; k < nxy; ++k ) {
	inds[k] = 0;
	for ( i = 0, j = nb-1; i < nb; j = i++ ) {
	    if ((((yb[i]<=y[k]) && (y[k]<yb[j])) ||
		 ((yb[j]<=y[k]) && (y[k]<yb[i]))) &&
		  (x[k] < (xb[j]-xb[i])*(y[k]-yb[i])/(yb[j]-yb[i])+xb[i]))

		  inds[k] = !inds[k];
	}
    }

}
