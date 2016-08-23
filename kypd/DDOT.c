#include <math.h>
#include "mex.h"

#define DDOT ddot_
#define DGEMM dgemm_

double DDOT( const size_t * n, const double * x, const size_t * incx, const double * y, const size_t * incy );

static double multiply( const double * A,
			const double * B,
			const size_t n )
{
  size_t st_one;
  st_one = 1;
  return DDOT( &n, A, &st_one, B, &st_one );
}


void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray*prhs[] )
     
{ 
  const mxArray * mx_A;
  const mxArray * mx_B;
  mxArray * mx_M;

  size_t n;
  
  mx_A = prhs[ 0 ];
  mx_B = prhs[ 1 ];

  if( mxGetN( mx_A ) != 1 )
    {
      mexErrMsgTxt( "Left matrix must have exactly 1 column." ); 
    }
  if( mxGetN( mx_B ) != 1 )
    {
      mexErrMsgTxt( "Right matrix must have exactly 1 column." ); 
    }
  n = mxGetM( mx_A );
  if( mxGetM( mx_B ) != n )
    {
      mexErrMsgTxt( "Matrices must same length." ); 
    }

  mx_M = mxCreateDoubleMatrix( 1, 1, mxREAL ); 
  plhs[ 0 ] = mx_M;
  
  *mxGetPr( mx_M ) = multiply( mxGetPr( mx_A ), mxGetPr( mx_B ), n );

  return;
}
