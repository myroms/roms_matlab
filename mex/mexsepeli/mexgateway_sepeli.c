/**********************************************************************
 *
 * mexgateway.c
 *
 * This file functions as the mex-file entry point.  The intended mexnc
 * operation is gleaned from the first argument, and then we transfer
 * control to the source file that handles either the NetCDF-2 or
 * NetCDF-3 API.
 *
 *********************************************************************/

/*
 * svn $Id$ 
 * mexgateway_sepeli.c 2117 2006-12-31 17:24:19Z johnevans007
 *
 */

# include <stdio.h>
# include <stdlib.h>

# include "mex.h"

void mexFunction ( int nlhs, mxArray *plhs[], int nrhs, const mxArray	*prhs[] ) {

    char error_message[1000];

    /*
     * Number of elements in a matlab matrix.
     * */
    int num_elements;


    /*
     * Pointer to the data part of a matlab matrix.
     * */
    double *pr;

    /*
     * Pointers to the input x and y grids.  These will be transposed
     * before passing on to the fortran routines.
     * */
    double *xptr, *yptr;

    /*
     * Output grids from the fortran sepeli routines.  These will be
     * transposed before passing back to matlab.
     * */
    double *xsep, *ysep;

    /*
     * parameters to fortran sepeli routine.
     * */
    double *x, *y;
    int l2, m2;
    double *seta, *sxi;
    double *x_out, *y_out;
    int ierr;


    /*
     * Loop indices.
     * */
    int i, j, k;


    /*
     * Size of the output arrays.  
     * */
    int mxsize[2];


    /*
     * There must be two output arguments.
     * */
    if ( nlhs != 2 ) {
        mxErrMsgTxt ( "MEXSEPELI ERROR:  There must be two output arguments for mexsepeli.\n" );
    }

    /*
     * There must be six output arguments.
     * */
    if ( nrhs != 6 ) {
        mxErrMsgTxt ( "MEXSEPELI ERROR:  There must be six input arguments for mexsepeli.\n" );
    }


    /*
     * Check that all the input arguments were numeric.
     * */
    if ( mxIsNumeric(prhs[0]) != 1 ) {
        mxErrMsgTxt ( "MEXSEPELI ERROR:  x input array must be numeric.\n" );
    }
    if ( mxIsNumeric(prhs[1]) != 1 ) {
        mxErrMsgTxt ( "MEXSEPELI ERROR:  y input array must be numeric.\n" );
    }
    if ( mxIsNumeric(prhs[2]) != 1 ) {
        mxErrMsgTxt ( "MEXSEPELI ERROR:  l2 input array must be numeric.\n" );
    }
    if ( mxIsNumeric(prhs[3]) != 1 ) {
        mxErrMsgTxt ( "MEXSEPELI ERROR:  m2 input array must be numeric.\n" );
    }
    if ( mxIsNumeric(prhs[4]) != 1 ) {
        mxErrMsgTxt ( "MEXSEPELI ERROR:  seta input array must be numeric.\n" );
    }
    if ( mxIsNumeric(prhs[5]) != 1 ) {
        mxErrMsgTxt ( "MEXSEPELI ERROR:  sxi input array must be numeric.\n" );
    }


    /*
     * Make sure that X and Y are the same size
     * */
    if ( mxGetM(prhs[0]) != mxGetM(prhs[1]) ) {
        mxErrMsgTxt ( "MEXSEPELI ERROR:  X and Y must be the same size.\n" );
    }
    if ( mxGetN(prhs[0]) != mxGetN(prhs[1]) ) {
        mxErrMsgTxt ( "MEXSEPELI ERROR:  X and Y must be the same size.\n" );
    }



    /*
     * Check that the l2 and m2 inputs were scalars.
     * */
    if ( mxGetNumberOfElements(prhs[2]) != 1 ) {
        mxErrMsgTxt ( "MEXSEPELI ERROR:  l2 argument must be scalar.\n" );
    }
    if ( mxGetNumberOfElements(prhs[3]) != 1 ) {
        mxErrMsgTxt ( "MEXSEPELI ERROR:  m2 argument must be scalar.\n" );
    }


    /*
     * Extract L2 and M2
     * */
    pr = (double *) mxGetPr ( prhs[2] );
    l2 = (int)(pr[0]);

    pr = (double *) mxGetPr ( prhs[3] );
    m2 = (int)(pr[0]);

    /*
     * Take the handoff from matlab.  Extract X and Y.
     * */
    num_elements = mxGetNumberOfElements ( prhs[0] );
    xptr = (double *) mxGetData ( prhs[0] );

    num_elements = mxGetNumberOfElements ( prhs[0] );
    yptr = (double *) mxGetData ( prhs[1] );


    /*
     * Allocate space for the transposed input arrays.
     * */
    x = mxCalloc ( num_elements, sizeof(double) );
    y = mxCalloc ( num_elements, sizeof(double) );
    xsep = mxCalloc ( num_elements, sizeof(double) );
    ysep = mxCalloc ( num_elements, sizeof(double) );

    /*
     * Transpose them.
     * */
    k = 0;
    for ( j = 0; j <= m2; ++j ) { 
        for ( i = 0; i <= l2; ++i ) {
            x[i*(m2+1)+j] = xptr[k];
            y[i*(m2+1)+j] = yptr[k];
            k = k + 1;
        }
    }

           

    /*
     * Extract seta and sxi
     * */
    seta = mxGetPr ( prhs[4] );
    sxi = mxGetPr ( prhs[5] );

    /*
     * Create the output arrays.
     * */
    mxsize[0] = l2+1;
    mxsize[1] = m2+1;
    plhs[0] = mxCreateNumericArray ( 2, mxsize, mxDOUBLE_CLASS, mxREAL );
    plhs[1] = mxCreateNumericArray ( 2, mxsize, mxDOUBLE_CLASS, mxREAL );

    x_out = (double *) mxGetData ( plhs[0] );
    y_out = (double *) mxGetData ( plhs[1] );

    sepeli_setup_ ( &l2, &m2, seta, sxi, x, y, xsep, ysep, &ierr );

    switch ( ierr ) {
        case 1:
            mexErrMsgTxt ( "MEXSEPELI ERROR:  range of independent variables is out of whack" );

        case 2:
            mexErrMsgTxt ( "MEXSEPELI ERROR:  boundary condition mbdcnd wrongly specified" );

        case 3:
            mexErrMsgTxt ( "MEXSEPELI ERROR:  boundary condition nbdcnd wrongly specified" );

        case 4:
            mexErrMsgTxt ( "MEXSEPELI ERROR:  linear system generated is not diagonally dominant" );

        case 5:
            mexErrMsgTxt ( "MEXSEPELI ERROR:  idmn is too small" );

        case 6:
            mexErrMsgTxt ( "MEXSEPELI ERROR:  m is too small or too large" );

        case 7:
            mexErrMsgTxt ( "MEXSEPELI ERROR:  n is too small or too large" );

        case 8:
            mexErrMsgTxt ( "MEXSEPELI ERROR:  iorder is not 2 or 4" );

        case 9:
            mexErrMsgTxt ( "MEXSEPELI ERROR:  intl is not 0 or 1" );

        case 10:
            mexErrMsgTxt ( "MEXSEPELI ERROR:  afun*dfun less than or equal to 0" );

        case 11:
            mexErrMsgTxt ( "MEXSEPELI ERROR:  work space length input in w(1) is not right" );

    }


    /*
     * Ok, so transpose the output back.
     * */
    k = 0;
    for ( j = 0; j <= m2; ++j ) { 
        for ( i = 0; i <= l2; ++i ) {
            x_out[k] = xsep[i*(m2+1)+j];
            y_out[k] = ysep[i*(m2+1)+j];
            k = k + 1;
        }
    }

           


    return;

}




