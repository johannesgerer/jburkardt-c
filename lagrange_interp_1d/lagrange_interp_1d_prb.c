# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <string.h>

# include "lagrange_interp_1d.h"
# include "test_interp_1d.h"
# include "r8lib.h"

int main ( );
void test02 ( int prob, int nd );
void test03 ( int prob, int nd );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
//  Purpose:
//
//    LAGRANGE_INTERP_1D_TEST tests LAGRANGE_INTERP_1D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 October 2012
//
//  Author:
//
//    John Burkardt
*/
{
  int nd_test_num = 6;

  int j;
  int nd;
  int nd_test[6] = { 4, 8, 16, 32, 64, 256 };
  int prob;
  int prob_num;

  timestamp ( );
  printf ( "\n" );
  printf ( "LAGRANGE_INTERP_1D_TEST:\n" );
  printf ( "  C version\n" );
  printf ( "  Test the LAGRANGE_INTERP_1D library.\n" );
  printf ( "  The R8LIB library is needed.\n" );
  printf ( "  These tests need the TEST_INTERP_1D library.\n" );

  prob_num = p00_prob_num ( );

  for ( prob = 1; prob <= prob_num; prob++ )
  {
    for ( j = 0; j < nd_test_num; j++ )
    {
      nd = nd_test[j];
      test02 ( prob, nd );
    }
  }

  for ( prob = 1; prob <= prob_num; prob++ )
  {
    for ( j = 0; j < nd_test_num; j++ )
    {
      nd = nd_test[j];
      test03 ( prob, nd );
    }
  }
/*
//  Terminate.
*/
  printf ( "\n" );
  printf ( "LAGRANGE_INTERP_1D_TEST:\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void test02 ( int prob, int nd )

/******************************************************************************/
/*
//  Purpose:
//
//    TEST02 tests LAGRANGE_VALUE_1D with evenly spaced data.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 September 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int PROB, the problem index.
//
//    Input, int ND, the number of data points to use.
*/
{
  double a;
  double b;
  int i;
  double int_error;
  double ld;
  double li;
  int ni;
  double *xd;
  double *xi;
  double *yd;
  double *yi;
  double ymax;
  double ymin;

  printf ( "\n" );
  printf ( "TEST02:\n" );
  printf ( "  Interpolate data from TEST_INTERP_1D problem #%d\n", prob );
  printf ( "  Use even spacing for data points.\n" );
  printf ( "  Number of data points = %d\n", nd );

  a = 0.0;
  b = 1.0;
  
  xd = r8vec_linspace_new ( nd, a, b );

  yd = p00_f ( prob, nd, xd );

  if ( nd < 10 )
  {
    r8vec2_print ( nd, xd, yd, "  Data array:" );
  }
/*
//  #1:  Does interpolant match function at interpolation points?
*/
  ni = nd;
  xi = r8vec_copy_new ( ni, xd );
  yi = lagrange_value_1d ( nd, xd, yd, ni, xi );

  int_error = r8vec_norm_affine ( nd, yi, yd ) / ( double ) ( ni );

  printf ( "\n" );
  printf ( "  L2 interpolation error averaged per interpolant node = %g\n", int_error );

  free ( xi );
  free ( yi );
/*
//  #2: Compare estimated curve length to piecewise linear (minimal) curve length.
//  Assume data is sorted, and normalize X and Y dimensions by (XMAX-XMIN) and
//  (YMAX-YMIN).
*/
  ymin = r8vec_min ( nd, yd );
  ymax = r8vec_max ( nd, yd );

  ni = 501;
  xi = r8vec_linspace_new ( ni, a, b );
  yi = lagrange_value_1d ( nd, xd, yd, ni, xi );

  ld = 0.0;
  for ( i = 0; i < nd - 1; i++ )
  {
    ld = ld + sqrt ( pow ( ( xd[i+1] - xd[i] ) / ( b - a ), 2 )
                   + pow ( ( yd[i+1] - yd[i] ) / ( ymax - ymin ), 2 ) ); 
  }

  li = 0.0;
  for ( i = 0; i < ni - 1; i++ )
  {
    li = li + sqrt ( pow ( ( xi[i+1] - xi[i] ) / ( b - a ), 2 )
                   + pow ( ( yi[i+1] - yi[i] ) / ( ymax - ymin ), 2 ) );
  }

  printf ( "\n" );
  printf ( "  Normalized length of piecewise linear interpolant = %g\n", ld );
  printf ( "  Normalized length of polynomial interpolant       = %g\n", li );

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
//  Purpose:
//
//    TEST03 tests LAGRANGE_VALUE_1D with Chebyshev spaced data.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 September 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int PROB, the problem index.
//
//    Input, int ND, the number of data points to use.
*/
{
  double a;
  double b;
  int i;
  double int_error;
  double ld;
  double li;
  int ni;
  double *xd;
  double *xi;
  double *yd;
  double *yi;
  double ymax;
  double ymin;

  printf ( "\n" );
  printf ( "TEST02:\n" );
  printf ( "  Interpolate data from TEST_INTERP_1D problem #%d\n", prob );
  printf ( "  Use Chebyshev spacing for data points.\n" );
  printf ( "  Number of data points = %d\n", nd );

  a = 0.0;
  b = 1.0;
  
  xd = r8vec_chebyspace_new ( nd, a, b );

  yd = p00_f ( prob, nd, xd );

  if ( nd < 10 )
  {
    r8vec2_print ( nd, xd, yd, "  Data array:" );
  }
/*
//  #1:  Does interpolant match function at interpolation points?
*/
  ni = nd;
  xi = r8vec_copy_new ( ni, xd );
  yi = lagrange_value_1d ( nd, xd, yd, ni, xi );

  int_error = r8vec_norm_affine ( nd, yi, yd ) / ( double ) ( ni );

  printf ( "\n" );
  printf ( "  L2 interpolation error averaged per interpolant node = %g\n", int_error );

  free ( xi );
  free ( yi );
/*
//  #2: Compare estimated curve length to piecewise linear (minimal) curve length.
//  Assume data is sorted, and normalize X and Y dimensions by (XMAX-XMIN) and
//  (YMAX-YMIN).
*/
  ymin = r8vec_min ( nd, yd );
  ymax = r8vec_max ( nd, yd );

  ni = 501;
  xi = r8vec_linspace_new ( ni, a, b );
  yi = lagrange_value_1d ( nd, xd, yd, ni, xi );

  ld = 0.0;
  for ( i = 0; i < nd - 1; i++ )
  {
    ld = ld + sqrt ( pow ( ( xd[i+1] - xd[i] ) / ( b - a ), 2 )
                   + pow ( ( yd[i+1] - yd[i] ) / ( ymax - ymin ), 2 ) ); 
  }

  li = 0.0;
  for ( i = 0; i < ni - 1; i++ )
  {
    li = li + sqrt ( pow ( ( xi[i+1] - xi[i] ) / ( b - a ), 2 )
                   + pow ( ( yi[i+1] - yi[i] ) / ( ymax - ymin ), 2 ) );
  }

  printf ( "\n" );
  printf ( "  Normalized length of piecewise linear interpolant = %g\n", ld );
  printf ( "  Normalized length of polynomial interpolant       = %g\n", li );

  free ( xd );
  free ( xi );
  free ( yd );
  free ( yi );

  return;
}
