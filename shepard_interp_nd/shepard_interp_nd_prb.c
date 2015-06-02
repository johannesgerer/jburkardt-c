# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <string.h>

# include "shepard_interp_nd.h"
# include "test_interp_nd.h"
# include "r8lib.h"

int main ( );
void test01 ( int prob, double p, int m, int nd );
void test02 ( int prob, double p, int m, int n1d );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for SHEPARD_INTERP_ND_PRB.

  Discussion:

    SHEPARD_INTERP_ND_PRB tests the SHEPARD_INTERP_ND library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 October 2012

  Author:

    John Burkardt
*/
{
  int j;
  int m;
  int n1d;
  int nd;
  double p;
  double p_test[4] = { 1.0, 2.0, 4.0, 8.0 };
  int p_test_num = 4;
  int prob;
  int prob_num;

  timestamp ( );
  printf ( "\n" );
  printf ( "SHEPARD_INTERP_ND_PRB:\n" );
  printf ( "  C version\n" );
  printf ( "  Test the SHEPARD_INTERP_ND library.\n" );
  printf ( "  The R8LIB library is needed.\n" );
  printf ( "  This test also needs the TEST_INTERP_ND library.\n" );
/*
  Look at Shepard interpolant on an irregular grid.
*/
  nd = 25;

  prob_num = p00_prob_num ( );

  for ( prob = 1; prob <= prob_num; prob++ )
  {
    for ( m = 2; m <= 5; m = m + 3 )
    {
      for ( j = 0; j < p_test_num; j++ )
      {
        p = p_test[j];
        test01 ( prob, p, m, nd );
      }

    }
  }
/*
  Look at Shepard interpolant on a regular N1D^M grid.
*/
  n1d = 5;

  for ( prob = 1; prob <= prob_num; prob++ )
  {
    for ( m = 2; m <= 5; m = m + 3 )
    {
      for ( j = 0; j < p_test_num; j++ )
      {
        p = p_test[j];
        test02 ( prob, p, m, n1d );
      }
    }
  }
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "SHEPARD_INTERP_ND_PRB:\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void test01 ( int prob, double p, int m, int nd )

/******************************************************************************/
/*
  Purpose:

    TEST01 tests SHEPARD_INTERP on an irregular grid.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 October 2012

  Author:

    John Burkardt

  Parameters:

    Input, int PROB, the problem number.

    Input, double P, the power used in the distance weighting.

    Input, int M, the spatial dimension.

    Input, int ND, the number of data points.
*/
{
  double app_error;
  double *c;
  int i;
  double int_error;
  int j;
  int ni;
  int seed;
  double *w;
  double *xd;
  double *xi;
  double *zd;
  double *ze;
  double *zi;

  printf ( "\n" );
  printf ( "TEST01:\n" );
  printf ( "  Interpolate data from TEST_INTERP_ND problem #%d\n", prob );
  printf ( "  using Shepard interpolation with P = %g\n", p );
  printf ( "  spatial dimension M = %d\n", m );
  printf ( "  and an irregular grid of ND = %d data points.\n", nd );
/*
  Set problem parameters:
*/
  seed = 123456789;
  c = r8vec_uniform_01_new ( m, &seed );
  w = r8vec_uniform_01_new ( m, &seed );

  xd = r8mat_uniform_01_new ( m, nd, &seed );

  zd = p00_f ( prob, m, c, w, nd, xd );
/*
  #1:  Does interpolant match function at interpolation points?
*/
  ni = nd;
  xi = r8mat_copy_new ( m, ni, xd );
  zi = shepard_interp_nd ( m, nd, xd, zd, p, ni, xi );

  int_error = r8vec_norm_affine ( ni, zi, zd ) / ( double ) ( ni );

  printf ( "\n" );
  printf ( "  L2 interpolation error averaged per interpolant node = %g\n", int_error );

  free ( xi );
  free ( zi );
/*
  #2: Approximation test.  Estimate the integral (f-interp(f))^2.
*/
  ni = 1000;
  ni = 50;
  xi = r8mat_uniform_01_new ( m, ni, &seed );
  zi = shepard_interp_nd ( m, nd, xd, zd, p, ni, xi );
  ze = p00_f ( prob, m, c, w, ni, xi );

  app_error = r8vec_norm_affine ( ni, zi, ze ) / ( double ) ( ni );

  printf ( "  L2 approximation error averaged per 1000 samples =     %g\n", app_error );

  free ( c );
  free ( w );
  free ( xd );
  free ( xi );
  free ( zd );
  free ( ze );
  free ( zi );

  return;
}
/******************************************************************************/

void test02 ( int prob, double p, int m, int n1d )

/******************************************************************************/
/*
  Purpose:

    TEST02 tests SHEPARD_INTERP_ND on a regular N1D^M grid.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 October 2012

  Author:

    John Burkardt

  Parameters:

    Input, int PROB, the problem number.

    Input, double P, the power used in the distance weighting.

    Input, int M, the spatial dimension.

    Input, int N1D, the number of points in 1D.
*/
{
  double a;
  double app_error;
  double b;
  double *c;
  int i;
  double int_error;
  int nd;
  int ni;
  int seed;
  double *w;
  double *x1d;
  double *xd;
  double *xi;
  double *zd;
  double *ze;
  double *zi;
/*
  Set problem parameters:
*/
  seed = 123456789;
  c = r8vec_uniform_01_new ( m, &seed );
  w = r8vec_uniform_01_new ( m, &seed );

  nd = i4_power ( n1d, m );

  printf ( "\n" );
  printf ( "TEST02:\n" );
  printf ( "  Interpolate data from TEST_INTERP_ND problem #%d\n", prob );
  printf ( "  using Shepard interpolation with P = %g\n", p );
  printf ( "  spatial dimension M = %d\n", m );
  printf ( "  and a regular grid of N1D^M = %d data points.\n", nd );

  a = 0.0;
  b = 1.0;

  x1d = r8vec_linspace_new ( n1d, a, b );

  xd = ( double * ) malloc ( m * nd * sizeof ( double ) );
  for ( i = 0; i < m; i++ )
  {
    r8vec_direct_product ( i, n1d, x1d, m, nd, xd );
  }

  zd = p00_f ( prob, m, c, w, nd, xd );
/*
  #1:  Does interpolant match function at interpolation points?
*/
  ni = nd;
  xi = r8mat_copy_new ( m, nd, xd );
  zi = shepard_interp_nd ( m, nd, xd, zd, p, ni, xi );

  int_error = r8vec_norm_affine ( ni, zi, zd ) / ( double ) ( ni );

  printf ( "\n" );
  printf ( "  L2 interpolation error averaged per interpolant node = %g\n", int_error );

  free ( xi );
  free ( zi );
/*
  #2: Approximation test.  Estimate the integral (f-interp(f))^2.
*/
  ni = 1000;
  xi = r8mat_uniform_01_new ( m, ni, &seed );

  zi = shepard_interp_nd ( m, nd, xd, zd, p, ni, xi );

  ze = p00_f ( prob, m, c, w, ni, xi );

  app_error = r8vec_norm_affine ( ni, zi, ze ) / ( double ) ( ni );

  printf ( "  L2 approximation error averaged per 1000 samples =     %g\n", app_error );

  free ( c );
  free ( w );
  free ( xd );
  free ( xi );
  free ( zd );
  free ( ze );
  free ( zi );

  return;
}
