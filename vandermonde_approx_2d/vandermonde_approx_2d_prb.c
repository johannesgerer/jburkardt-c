# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <string.h>

# include "vandermonde_approx_2d.h"
# include "test_interp_2d.h"
# include "qr_solve.h"
# include "r8lib.h"

int main ( );
void test01 ( int prob, int grid, int m );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    VANDERMONDE_APPROX_2D_TEST tests VANDERMONDE_APPROX_2D.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 October 2012

  Author:

    John Burkardt
*/
{
  int j;
  int m;
  int m_test[5] = { 0, 1, 2, 4, 8 };
  int m_test_num = 5;
  int grid;
  int prob;
  int prob_num;

  timestamp ( );
  printf ( "\n" );
  printf ( "VANDERMONDE_APPROX_2D_TEST:\n" );
  printf ( "  C version\n" );
  printf ( "  Test the VANDERMONDE_APPROX_2D library.\n" );
  printf ( "  This test also needs the TEST_INTERP_2D library.\n" );

  prob_num = f00_num ( );

  for ( prob = 1; prob <= prob_num; prob++ )
  {
    grid = 1;
    for ( j = 0; j < m_test_num; j++ )
    {
      m = m_test[j];
      test01 ( prob, grid, m );
    }
  }
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "VANDERMONDE_APPROX_2D_TEST:\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void test01 ( int prob, int grd, int m )

/******************************************************************************/
/*
  Purpose:

    VANDERMONDE_APPROX_2D_TEST01 tests VANDERMONDE_APPROX_2D_MATRIX.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 October 2012

  Author:

    John Burkardt

  Parameters:

    Input, int PROB, the problem number.

    Input, int GRD, the grid number.
    (Can't use GRID as the name because that's also a plotting function.)

    Input, int M, the total polynomial degree.
*/
{
  double *a;
  double app_error;
  double *c;
  int nd;
  int ni;
  int tm;
  double *xd;
  double *xi;
  double *yd;
  double *yi;
  double *zd;
  double *zi;

  printf ( "\n" );
  printf ( "TEST01:\n" );
  printf ( "  Approximate data from TEST_INTERP_2D problem #%d\n", prob );
  printf ( "  Use grid from TEST_INTERP_2D with index #%d\n", grd );
  printf ( "  Using polynomial approximant of total degree %d\n", m );

  nd = g00_size ( grd );
  printf ( "  Number of data points = %d\n", nd );

  xd = ( double * ) malloc ( nd * sizeof ( double ) );
  yd = ( double * ) malloc ( nd * sizeof ( double ) );
  g00_xy ( grd, nd, xd, yd );

  zd = ( double * ) malloc ( nd * sizeof ( double ) );
  f00_f0 ( prob, nd, xd, yd, zd );

  if ( nd < 10 )
  {
    r8vec3_print ( nd, xd, yd, zd, "  X, Y, Z data:" );
  }
/*
  Compute the Vandermonde matrix.
*/
  tm = triangle_num ( m + 1 );
  a = vandermonde_approx_2d_matrix ( nd, m, tm, xd, yd );
/*
  Solve linear system.
*/
  c = qr_solve ( nd, tm, a, zd );
/*
  #1:  Does approximant match function at data points?
*/
  ni = nd;
  xi = r8vec_copy_new ( ni, xd );
  yi = r8vec_copy_new ( ni, yd );
  zi = r8poly_value_2d ( m, c, ni, xi, yi );

  app_error = r8vec_norm_affine ( ni, zi, zd ) / ( double ) ( ni );

  printf ( "\n" );
  printf ( "  L2 data approximation error = %g\n", app_error );

  free ( a );
  free ( c );
  free ( xd );
  free ( xi );
  free ( yd );
  free ( yi );
  free ( zd );
  free ( zi );

  return;
}
