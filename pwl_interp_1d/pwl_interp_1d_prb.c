# include <stdlib.h>
# include <stdio.h>
# include <string.h>

# include "pwl_interp_1d.h"
# include "test_interp.h"
# include "r8lib.h"

int main ( );
void pwl_interp_1d_test01 ( int prob );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    PWL_INTERP_1D_TEST tests PWL_INTERP_1D.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 September 2012

  Author:

    John Burkardt
*/
{
  int prob;
  int prob_num;

  timestamp ( );
  printf ( "\n" );
  printf ( "PWL_INTERP_1D_TEST:\n" );
  printf ( "  C version\n" );
  printf ( "  Test the PWL_INTERP_1D library.\n" );
  printf ( "  The R8LIB library is needed.\n" );
  printf ( "  The test needs the TEST_INTERP library.\n" );

  prob_num = p00_prob_num ( );
  for ( prob = 1; prob <= prob_num; prob++ )
  {
    pwl_interp_1d_test01 ( prob );
  }
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "PWL_INTERP_1D_TEST:\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void pwl_interp_1d_test01 ( int prob )

/******************************************************************************/
/*
  Purpose:

    PWL_INTERP_1D_TEST01 tests PWL_INTERP_1D.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 September 2012

  Author:

    John Burkardt

  Parameters:

    Input, int PROB, the problem index.
*/
{
  int i;
  double int_error;
  int nd;
  int ni;
  double *xd;
  double *xi;
  double *xy;
  double *yd;
  double *yi;

  printf ( "\n" );
  printf ( "PWL_INTERP_1D_TEST01:\n" );
  printf ( "  PWL_INTERP_1D evaluates the piecewise linear interpolant.\n" );
  printf ( "  Interpolate data from TEST_INTERP problem #%d\n", prob );

  nd = p00_data_num ( prob );
  printf ( "  Number of data points = %d\n", nd );

  xy = p00_data ( prob, 2, nd );
  
  r8mat_transpose_print ( 2, nd, xy, "  Data array:" );

  xd = ( double * ) malloc ( nd * sizeof ( double ) );
  yd = ( double * ) malloc ( nd * sizeof ( double ) );

  for ( i = 0; i < nd; i++ )
  {
    xd[i] = xy[0+2*i];
    yd[i] = xy[1+2*i];
  }
/*
  Does interpolant match function at interpolation points?
*/
  ni = nd;

  xi = ( double * ) malloc ( ni * sizeof ( double ) );
  for ( i = 0; i < ni; i++ )
  {
    xi[i] = xd[i];
  }

  yi = pwl_interp_1d ( nd, xd, yd, ni, xi );

  int_error = r8vec_diff_norm ( ni, yi, yd ) / ( double ) ( ni );

  printf ( "\n" );
  printf ( "  L2 interpolation error averaged per interpolant node = %g\n", int_error );

  free ( xd );
  free ( xi );
  free ( xy );
  free ( yd );
  free ( yi );

  return;
}
