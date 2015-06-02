# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <string.h>

# include "shepard_interp_2d.h"
# include "test_interp_2d.h"
# include "r8lib.h"

int main ( );
void test01 ( int prob, int g, double p );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for SHEPARD_INTERP_2D_PRB.

  Discussion:

    SHEPARD_INTERP_2D_PRB tests the SHEPARD_INTERP_2D library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 October 2012

  Author:

    John Burkardt
*/
{
  int g;
  int j;
  double p;
  double p_test[4] = { 1.0, 2.0, 4.0, 8.0 };
  int p_test_num = 4;
  int prob;
  int prob_num;

  timestamp ( );
  printf ( "\n" );
  printf ( "SHEPARD_INTERP_2D_PRB:\n" );
  printf ( "  C version\n" );
  printf ( "  Test the SHEPARD_INTERP_2D library.\n" );
  printf ( "  The R8LIB library is needed.\n" );
  printf ( "  This test also needs the TEST_INTERP_2D library.\n" );

  prob_num = f00_num ( );
  g = 1;

  for ( prob = 1; prob <= prob_num; prob++ )
  {
    for ( j = 0; j < p_test_num; j++ )
    {
      p = p_test[j];
      test01 ( prob, g, p );
    }
  }
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "SHEPARD_INTERP_2D_PRB:\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void test01 ( int prob, int g, double p )

/******************************************************************************/
/*
  Purpose:

    TEST01 tests SHEPARD_INTERP_2D.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 October 2012

  Author:

    John Burkardt

  Parameters:

    Input, int PROB, the problem number.

    Input, int G, the grid number.

    Input, double P, the power used in the distance weighting.
*/
{
  int debug = 0;
  double int_error;
  int nd;
  int ni;
  double *xd;
  double *xi;
  double *yd;
  double *yi;
  double *zd;
  double *zi;

  printf ( "\n" );
  printf ( "TEST01:\n" );
  printf ( "  Interpolate data from TEST_INTERP_2D problem #%d\n", prob );
  printf ( "  using grid #%d\n", g );
  printf ( "  using Shepard interpolation with P = %g\n", p );

  nd = g00_size ( g );
  printf ( "  Number of data points = %d\n", nd );

  xd = ( double * ) malloc ( nd * sizeof ( double ) );
  yd = ( double * ) malloc ( nd * sizeof ( double ) );
  g00_xy ( g, nd, xd, yd );
  
  zd = ( double * ) malloc ( nd * sizeof ( double ) );
  f00_f0 ( prob, nd, xd, yd, zd );

  if ( debug )
  {
    r8vec3_print ( nd, xd, yd, zd, "  X, Y, Z data:" );
  }
/*
  #1:  Does interpolant match function at interpolation points?
*/
  ni = nd;
  xi = r8vec_copy_new ( ni, xd );
  yi = r8vec_copy_new ( ni, yd );

  zi = shepard_interp_2d ( nd, xd, yd, zd, p, ni, xi, yi );

  int_error = r8vec_norm_affine ( ni, zi, zd ) / ( double ) ( ni );

  printf ( "\n" );
  printf ( "  L2 interpolation error averaged per interpolant node = %g\n",
    int_error );

  free ( xd );
  free ( xi );
  free ( yd );
  free ( yi );
  free ( zd );
  free ( zi );

  return;
}
