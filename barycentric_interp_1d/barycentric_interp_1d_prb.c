# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <string.h>

# include "barycentric_interp_1d.h"
# include "test_interp_1d.h"
# include "r8lib.h"

int main ( );
void test01 ( int prob, int nd );
void test02 ( int prob, int nd );
void test03 ( int prob, int nd );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for BARYCENTRIC_INTERP_1D_PRB.

  Discussion:

    BARYCENTRIC_INTERP_1D_PRB tests the BARYCENTRIC_INTERP_1D library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 October 2012

  Author:

    John Burkardt
*/
{
  int i;
  int nd;
  int nd_test[6] = { 4, 8, 16, 32, 64, 1000 };
  int nd_test_num = 6;
  int prob;
  int prob_num;

  timestamp ( );
  printf ( "\n" );
  printf ( "BARYCENTRIC_INTERP_1D_PRB:\n" );
  printf ( "  C version\n" );
  printf ( "  Test the BARYCENTRIC_INTERP_1D library.\n" );
  printf ( "  The R8LIB library is needed.\n" );
  printf ( "  These tests need the TEST_INTERP_1D library.\n" );

  prob_num = p00_prob_num ( );

  for ( prob = 1; prob <= prob_num; prob++ )
  {
    for ( i = 0; i < nd_test_num; i++ )
    {
      nd = nd_test[i];
      test01 ( prob, nd );
    }
  }

  for ( prob = 1; prob <= prob_num; prob++ )
  {
    for ( i = 0; i < nd_test_num; i++ )
    {
      nd = nd_test[i];
      test02 ( prob, nd );
    }
  }

  for ( prob = 1; prob <= prob_num; prob++ )
  {
    for ( i = 0; i < nd_test_num; i++ )
    {
      nd = nd_test[i];
      test03 ( prob, nd );
    }
  }
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "BARYCENTRIC_INTERP_1D_PRB:\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void test01 ( int prob, int nd )

/******************************************************************************/
/*
  Purpose:

    BARYCENTRIC_INTERP_1D_TEST01 tests LAGCHEBY1_INTERP_1D.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 September 2012

  Author:

    John Burkardt

  Parameters:

    Input, int PROB, the problem index.

    Input, int ND, the number of data points to use.
*/
{
  double a;
  double b;
  double int_error;
  int ni;
  double *xd;
  double *xi;
  double *yd;
  double *yi;

  printf ( "\n" );
  printf ( "BARYCENTRIC_INTERP_1D_TEST01:\n" );
  printf ( "  Interpolate data from TEST_INTERP_1D problem #%d\n", prob );
  printf ( "  Use Chebyshev Type 1 spacing for data points.\n" );
  printf ( "  Number of data points = %d\n", nd );
/*
  Define the data.
*/
  a =  0.0;
  b = +1.0;
  xd = r8vec_cheby1space_new ( nd, a, b );
  yd = p00_f ( prob, nd, xd );

  if ( nd < 10 )
  {
    r8vec2_print ( nd, xd, yd, "  Data array:" );
  }
/*
  #1:  Does the interpolant match the function at the interpolation points?
*/
  ni = nd;
  xi = r8vec_copy_new ( ni, xd );
  yi = lagcheby1_interp_1d ( nd, xd, yd, ni, xi );

  int_error = r8vec_norm_affine ( ni, yi, yd ) / ( double ) ( ni );

  printf ( "\n" );
  printf ( "  L2 interpolation error averaged per interpolant node = %g\n", int_error );

  free ( xd );
  free ( xi );
  free ( yd );
  free ( yi );

  return;
}
/******************************************************************************/

void test02 ( int prob, int nd )

/******************************************************************************/
/*
  Purpose:

    BARYCENTRIC_INTERP_1D_TEST02 tests LAGCHEBY2_INTERP_1D.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 September 2012

  Author:

    John Burkardt

  Parameters:

    Input, int PROB, the problem index.

    Input, int ND, the number of data points to use.
*/
{
  double a;
  double b;
  double int_error;
  int ni;
  double *xd;
  double *xi;
  double *yd;
  double *yi;

  printf ( "\n" );
  printf ( "BARYCENTRIC_INTERP_1D_TEST02:\n" );
  printf ( "  Interpolate data from TEST_INTERP_1D problem #%d\n", prob );
  printf ( "  Use Chebyshev Type 2 spacing for data points.\n" );
  printf ( "  Number of data points = %d\n", nd );
/*
  Define the data.
*/
  a =  0.0;
  b = +1.0;
  xd = r8vec_cheby2space_new ( nd, a, b );
  yd = p00_f ( prob, nd, xd );

  if ( nd < 10 )
  {
    r8vec2_print ( nd, xd, yd, "  Data array:" );
  }
/*
  #1:  Does the interpolant match the function at the interpolation points?
*/
  ni = nd;
  xi = r8vec_copy_new ( ni, xd );
  yi = lagcheby2_interp_1d ( nd, xd, yd, ni, xi );

  int_error = r8vec_norm_affine ( ni, yi, yd ) / ( double ) ( ni );

  printf ( "\n" );
  printf ( "  L2 interpolation error averaged per interpolant node = %g\n", int_error );

  free ( xd );
  free ( xi );
  free ( yd );
  free ( yi );

  return;
}
/******************************************************************************/

void test03 ( int prob, int nd )

/******************************************************************************/
/*
  Purpose:

    BARYCENTRIC_INTERP_1D_TEST03 tests LAGEVEN_INTERP_1D.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 September 2012

  Author:

    John Burkardt

  Parameters:

    Input, int PROB, the problem index.

    Input, int ND, the number of data points to use.
*/
{
  double a;
  double b;
  int i;
  double int_error;
  int ni;
  double *xd;
  double *xi;
  double *yd;
  double *yi;

  printf ( "\n" );
  printf ( "BARYCENTRIC_INTERP_1D_TEST03:\n" );
  printf ( "  Interpolate data from TEST_INTERP_1D problem #%d\n", prob );
  printf ( "  Use even spacing for data points.\n" );
  printf ( "  Number of data points = %d\n", nd );
/*
  Define the data.
*/
  a =  0.0;
  b = +1.0;
  xd = r8vec_midspace_new ( nd, a, b );
  yd = p00_f ( prob, nd, xd );

  if ( nd < 10 )
  {
    r8vec2_print ( nd, xd, yd, "  Data array:" );
  }
/*
  #1:  Does the interpolant match the function at the interpolation points?
*/
  ni = nd;
  xi = r8vec_copy_new ( ni, xd );
  yi = lageven_interp_1d ( nd, xd, yd, ni, xi );

  int_error = r8vec_norm_affine ( ni, yi, yd ) / ( double ) ( ni );

  printf ( "\n" );
  printf ( "  L2 interpolation error averaged per interpolant node = %g\n", int_error );

  free ( xd );
  free ( xi );
  free ( yd );
  free ( yi );

  return;
}
