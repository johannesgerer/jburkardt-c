# include <stdlib.h>
# include <stdio.h>
# include <string.h>

# include "nearest_interp_1d.h"
# include "test_interp.h"
# include "r8lib.h"

int main ( );
void nearest_interp_1d_test01 ( int prob, int ni );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    NEAREST_INTERP_1D_TEST tests NEAREST_INTERP_1D.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 September 2012

  Author:

   John Burkardt
*/
{
  int ni;
  int prob;
  int prob_num;

  timestamp ( );
  printf ( "\n" );
  printf ( "NEAREST_INTERP_1D_TEST:\n" );
  printf ( "  C version\n" );
  printf ( "  Test the NEAREST_INTERP_1D library.\n" );
  printf ( "  The R8LIB library is needed.\n" );
  printf ( "  The test needs the TEST_INTERP library.\n" );

  prob_num = p00_prob_num ( );

  ni = 11;
  for ( prob = 1; prob <= prob_num; prob++ )
  {
    nearest_interp_1d_test01 ( prob, ni );
  }
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "NEAREST_INTERP_1D_TEST:\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void nearest_interp_1d_test01 ( int prob, int ni )

/******************************************************************************/
/*
  Purpose:

    NEAREST_INTERP_1D_TEST01 tests NEAREST_INTERP_1D.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 September 2012

  Author:

    John Burkardt

  Parameters:

    Input, int PROB, the index of the problem.

    Input, int NI, the number of interpolation points.
*/
{
  double *d;
  int j;
  int nd;
  char title[80];
  double *xd;
  double *xi;
  double xd_max;
  double xd_min;
  double *yd;
  double *yi;

  printf ( "\n" );
  printf ( "NEAREST_INTERP_1D_TEST01\n" );
  printf ( "  Sample the nearest neighbor interpolant for problem # %d\n", prob );

  nd = p00_data_num ( prob );

  d = p00_data ( prob, 2, nd );

  xd = ( double * ) malloc ( nd * sizeof ( double ) );
  yd = ( double * ) malloc ( nd * sizeof ( double ) );

  for ( j = 0; j < nd; j++ )
  {
    xd[j] = d[0+j*2];
    yd[j] = d[1+j*2];
  }

  xd_min = r8vec_min ( nd, xd );
  xd_max = r8vec_max ( nd, xd );

  xi = r8vec_linspace_new ( ni, xd_min, xd_max );
  yi = nearest_interp_1d ( nd, xd, yd, ni, xi );

  sprintf ( title, "X, Y for problem %d", prob );

  r8vec2_print ( ni, xi, yi, title );

  free ( d );
  free ( xd );
  free ( xi );
  free ( yd );
  free ( yi );

  return;
}

