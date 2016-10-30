#include <math.h>
#include "mex.h"

#define DGEMV dgemv_
#define DGEMM dgemm_
#define DDOT ddot_

double DDOT( const size_t * n, const double * x, const size_t * incx, const double * y, const size_t * incy );

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

static void compute_vW_sparse( double * vW, const double * W, const int * Wir, const int * Wjc,
			       const double * V,
			       const size_t n, const size_t rtot )
{
  double tmp;
  const double * srcW;
  const int * srcWir;
  const int * srcWjc;
  const int * srcWjcend;
  const double * srcbv1;
  const double * srcbvend;
  const double * srcWend;
  double * dstvW;

  srcWjcend = Wjc + ( 1 + n );
  srcbvend = V + rtot * n;
  dstvW = vW;
  for( srcbv1 = V; srcbv1 != srcbvend; srcbv1 += n )
    {
      srcW = W;
      srcWir = Wir;
      for( srcWjc = Wjc + 1; srcWjc != srcWjcend; ++srcWjc, ++dstvW )
	{
	  tmp = 0;
	  srcWend = W + (*srcWjc);
	  for( ; srcW != srcWend; ++srcW, ++srcWir )
	    {
	      tmp += (*(srcbv1 + (*srcWir))) * (*srcW);
	    }
	  *dstvW = tmp;
	}
    }
}

static void compute_vW_full( double * vW, const double * W,
			     const double * V,
			     const size_t n, const size_t rtot )
{
  double d_one;
  double d_zero;
  d_one = 1;
  d_zero = 0;
  DGEMM( "T", "N", &n, &rtot, &n, &d_one, W, &n, V, &n, &d_zero, vW, &n );
}

static void lowrankschur( double * M_res, const double * vW,
			  const size_t * rcum, const double * alpha, const double * V,
			  const size_t n, const size_t m )
{
  size_t rtot;
  double tmp;
  double tmp2;
  double tmp3;
  const double * srcv2;
  const double * srcv2begin;
  const double * srcv2end;
  const size_t * srcrcum1;
  const size_t * srcrcum2;
  const double * srca1;
  const double * srca2;
  const double * a1begin;
  const double * a2begin;
  const double * a1end;
  const double * a2end;
  const double * srcbvW;
  const double * srcb2vW;
  const double * srcvW;
  size_t i;
  size_t j;
  size_t st_one;
  st_one = 1;

  a1begin = alpha;
  for( i = 0, srcrcum1 = rcum; i < m; ++i, ++srcrcum1, a1begin = a1end )
    {
      srcbvW = vW + (*srcrcum1) * n;
      a1end = alpha + (*(srcrcum1+1));

      a2begin = a1begin;
      for( j = i, srcrcum2 = srcrcum1; j < m; ++j, ++srcrcum2, a2begin = a2end )
	{
	  a2end = alpha + (*(srcrcum2+1));
	  srcv2begin = V + ( (*srcrcum2) * n );
	  tmp = 0;
	  for( srca1 = a1begin, srcb2vW = srcbvW; srca1 != a1end; ++srca1, srcb2vW += n ) /* k = rcum( i ) : ( rcum( i + 1 ) - 1 ) */
	    {
	      tmp2 = (*srca1);
	      srcv2 = srcv2begin;
	      for( srca2 = a2begin; srca2 != a2end; ++srca2 ) /* l = rcum( j ) : ( rcum( j + 1 ) - 1 ) */
		{
		  tmp3 = DDOT( &n, srcb2vW, &st_one, srcv2, &st_one );
		  srcv2 += n;
		  tmp += tmp2 * (*srca2) * tmp3 * tmp3;
		}
	    }
	  
	  M_res[ i * m + j ] = tmp;
	  M_res[ j * m + i ] = tmp;
	}
    }
}

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray*prhs[] )
     
{ 
  const mxArray * mx_W;
  const mxArray * mx_aux;
  const mxArray * mx_rcum;
  size_t * rcum;
  const mxArray * mx_alpha;
  const mxArray * mx_V;
  mxArray * mx_M;

  size_t n;
  size_t m;
  size_t rtot;
  size_t subBegin;
  size_t subEnd;
  size_t * dst;
  size_t * end;
  const double * src;

  double * vW;
  
  mx_W = prhs[ 0 ];
  mx_aux = prhs[ 1 ];

  mx_rcum = mxGetCell( mx_aux, 0 );
  mx_alpha = mxGetCell( mx_aux, 1 );
  mx_V = mxGetCell( mx_aux, 2 );

  m = mxGetM( mx_rcum ) - 1;
  rcum = mxMalloc( ( m + 1 ) * sizeof( size_t ) );
  end = rcum + m + 1;
  for( dst = rcum, src = mxGetPr( mx_rcum ); dst != end; ++dst, ++src )
    {
      *dst = (size_t)(*src) - 1;
    }


  n = mxGetM( mx_W );
  if( mxGetM( mx_V ) != n )
    {
      mexErrMsgTxt( "The number of rows in V doesn't match the size of W." ); 
    }
  /*
  if( mxGetN( mx_V ) != rcum[ m ] )
    {
      mexErrMsgTxt( "The number of columns in V doesn't match rcum." ); 
    }
  if( mxGetM( mx_alpha ) != rcum[ m ] )
    {
      mexErrMsgTxt( "The number of rows in alpha doesn't match rcum." ); 
    }
  */
  /*
  if( mxGetN( mx_alpha ) != 1 )
    {
      mexErrMsgTxt( "The number of columns in alpha isn't 1." ); 
    }
  */
  /* Create a matrix for the return argument */ 
  mx_M = mxCreateDoubleMatrix( m, m, mxREAL ); 
  plhs[ 0 ] = mx_M;
  
  rtot = rcum[ m ]; /* rcum is zero-based, while mx_rcum is one-based */
  vW = malloc( rtot * n * sizeof( double ) ); /* vW is stored as W.' * V */

  /* Compute vW */
  if( mxIsSparse( mx_W ) )
    {
      compute_vW_sparse( vW, mxGetPr( mx_W ), mxGetIr( mx_W ), mxGetJc( mx_W ),
			 mxGetPr( mx_V ),
			 n, rtot );
    }
  else
    {
      compute_vW_full( vW, mxGetPr( mx_W ),
		       mxGetPr( mx_V ),
		       n, rtot );
    }

  /* Compute M */
  lowrankschur( mxGetPr( mx_M ), vW,
		rcum, mxGetPr( mx_alpha ), mxGetPr( mx_V ),
		n, m );

  mxFree( rcum );
  free( vW );

  return;
}
