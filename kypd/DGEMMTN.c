#include <math.h>
#include "mex.h"

#define DGEMV dgemv_
#define DGEMM dgemm_
#define DDOT ddot_

void displayArray( const double * pr, size_t m, size_t n )
{
  size_t r;
  size_t c;
  for( r = 0; r < m; ++r )
    {
      for( c = 0; c < n; ++c )
	{
	  mexPrintf( "%14.5e", *( pr + ( r + c * m ) ) );
	}
      mexPrintf( "\n" );
    }
}

static void multiply( double * M, 
		      const double * A,
		      const double * B,
		      const size_t m,
		      const size_t n,
		      const size_t k )
{
  double done;
  double dzero;
  done = 1;
  dzero = 0;
  DGEMM( "T", "N", &m, &n, &k, &done, A, &k, B, &k, &dzero, M, &m );
}

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray*prhs[] )
     
{ 
  const mxArray * mx_A;
  const mxArray * mx_B;
  mxArray * mx_M;

  size_t n;
  size_t m;
  size_t k;
  
  mx_A = prhs[ 0 ];
  mx_B = prhs[ 1 ];

  m = mxGetN( mx_A );
  n = mxGetN( mx_B );
  k = mxGetM( mx_A );
  if( mxGetM( mx_B ) != k )
    {
      mexErrMsgTxt( "Matrix dimensions don't match." ); 
    }

  mx_M = mxCreateDoubleMatrix( m, n, mxREAL ); 
  plhs[ 0 ] = mx_M;
  
  multiply( mxGetPr( mx_M ), mxGetPr( mx_A ), mxGetPr( mx_B ), m, n, k );

  return;
}
