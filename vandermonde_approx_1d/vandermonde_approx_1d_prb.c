# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <string.h>

# include "vandermonde_approx_1d.h"
# include "test_interp.h"
# include "qr_solve.h"
# include "r8lib.h"

int main ( );
void test01 ( int prob, int m );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for VANDERMONDE_APPROX_1D_PRB.

  Discussion:

    VANDERMONDE_APPROX_1D_PRB tests the VANDERMONDE_APPROX_1D library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    10 October 2012

  Author:

    John Burkardt
*/
{
  int j;
  int m;
  int m_test[8] = { 0, 1, 2, 3, 4, 5, 9, 12 };
  int m_test_num = 8;
  int prob;
  int prob_num;

  timestamp ( );
  printf ( "\n" );
  printf ( "VANDERMONDE_APPROX_1D_PRB:\n" );
  printf ( "  C version\n" );
  printf ( "  Test the VANDERMONDE_APPROX_1D library.\n" );
  printf ( "  The R8LIB library is needed.\n" );
  printf ( "  The QR_SOLVE library is needed.\n" );
  printf ( "  The test needs the CONDITION library.\n" );
  printf ( "  The test needs the TEST_INTERP libary.\n" );

  prob_num = p00_prob_num ( );
  for ( prob = 1; prob <= prob_num; prob++ )
  {
    for ( j = 0; j < m_test_num; j++ )
    {
      m = m_test[j];
      test01 ( prob, m );
    }
  }
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "VANDERMONDE_APPROX_1D_PRB:\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void test01 ( int prob, int m )

/******************************************************************************/
/*
  Purpose:

    TEST01 tests VANDERMONDE_APPROX_1D_MATRIX.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    10 October 2012

  Author:

    John Burkardt

  Parameters:

    Input, int PROB, the problem number.

    Input, int M, the polynomial degree.
*/
{
  double *a;
  double app_error;
  double *c;
  int debug = 0;
  int i;
  int j;
  double ld;
  double li;
  int nd;
  int ni;
  double *xd;
  double *xi;
  double xmax;
  double xmin;
  double *xy;
  double *yd;
  double *yi;
  double ymax;
  double ymin;

  printf ( "\n" );
  printf ( "TEST01:\n" );
  printf ( "  Approximate data from TEST_INTERP problem #%d\n", prob );

  nd = p00_data_num ( prob );
  printf ( "  Number of data points = %d\n", nd );

  xy = p00_data ( prob, 2, nd );
  
  if ( debug )
  {
    r8mat_transpose_print ( 2, nd, xy, "  Data array:" );
  }

  xd = ( double * ) malloc ( nd * sizeof ( double ) );
  yd = ( double * ) malloc ( nd * sizeof ( double ) );
  for ( i = 0; i < nd; i++ )
  {
    xd[i] = xy[0+i*2];
    yd[i] = xy[1+i*2];
  }
/*
  Compute the Vandermonde matrix.
*/
  printf ( "  Using polynomial approximant of degree %d\n", m );

  a = vandermonde_approx_1d_matrix ( nd, m, xd );
/*
  Solve linear system.
*/
  c = qr_solve ( nd, m + 1, a, yd );
/*
  #1:  Does approximant match function at data points?
*/
  ni = nd;
  xi = r8vec_copy_new ( ni, xd );
  yi = r8poly_value ( m, c, ni, xi );

  app_error = r8vec_norm_affine ( ni, yi, yd ) / ( double ) ( ni );

  printf ( "\n" );
  printf ( "  L2 data approximation error = %g\n", app_error );

  free ( xi );
  free ( yi );
/*
  #2: Compare estimated curve length to piecewise linear (minimal) curve length.
  Assume data is sorted, and normalize X and Y dimensions by (XMAX-XMIN) and
  (YMAX-YMIN).
*/
  xmin = r8vec_min ( nd, xd );
  xmax = r8vec_max ( nd, xd );
  ymin = r8vec_min ( nd, yd );
  ymax = r8vec_max ( nd, yd );

  ni = 501;
  xi = r8vec_linspace_new ( ni, xmin, xmax );

  yi = r8poly_value ( m, c, ni, xi );

  ld = 0.0;
  for ( i = 0; i < nd - 1; i++ )
  {
    ld = ld + sqrt ( pow ( ( xd[i+1] - xd[i] ) / ( xmax - xmin ), 2 )
                   + pow ( ( yd[i+1] - yd[i] ) / ( ymax - ymin ), 2 ) ); 
  }

  li = 0.0;
  for ( i = 0; i < ni - 1; i++ )
  {
    li = li + sqrt ( pow ( ( xi[i+1] - xi[i] ) / ( xmax - xmin ), 2 )
                   + pow ( ( yi[i+1] - yi[i] ) / ( ymax - ymin ), 2 ) );
  }

  printf ( "\n" );
  printf ( "  Normalized length of piecewise linear interpolant = %g\n", ld );
  printf ( "  Normalized length of polynomial interpolant       = %g\n", li );

  free ( a );
  free ( c );
  free ( xd );
  free ( xi );
  free ( xy );
  free ( yd );
  free ( yi );

  return;
}
