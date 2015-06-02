# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>

# include "chebyshev_interp_1d.h"
# include "test_interp.h"
# include "r8lib.h"

int main ( );
void test01 ( int prob );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for CHEBYSHEV_INTERP_1D_PRB.

  Discussion:

    CHEBYSHEV_INTERP_1D_PRB tests the CHEBYSHEV_INTERP_1D library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    27 September 2012

  Author:

    John Burkardt
*/
{
  int prob;
  int prob_num;

  timestamp ( );
  printf ( "\n" );
  printf ( "CHEBYSHEV_INTERP_1D_PRB:\n" );
  printf ( "  C version\n" );
  printf ( "  Test the CHEBYSHEV_INTERP_1D library.\n" );
  printf ( "  The QR_SOLVE and R8LIB libaries are needed.\n" );
  printf ( "  The test needs the TEST_INTERP library.\n" );

  prob_num = p00_prob_num ( );
  for ( prob = 1; prob <= prob_num; prob++ )
  {
    test01 ( prob );
  }
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "CHEBYSHEV_INTERP_1D_PRB:\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void test01 ( int prob )

/******************************************************************************/
/*
  Purpose:

    TEST01 tests CHEBYSHEV_VALUE_1D.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    27 September 2012

  Author:

    John Burkardt
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
  printf ( "CHEBYSHEV_INTERP_1D_TEST01:\n" );
  printf ( "  Interpolate data from TEST_INTERP problem #%d\n", prob );

  nd = p00_data_num ( prob );
  printf ( "  Number of data points = %d\n", nd );

  xy = p00_data ( prob, 2, nd );
  
  r8mat_transpose_print ( 2, nd, xy, "  Data array:" );

  xd = ( double * ) malloc ( nd * sizeof ( double ) );
  yd = ( double * ) malloc ( nd * sizeof ( double ) );
  
  for ( i = 0; i < nd; i++ )
  {
    xd[i] = xy[0+i*2];
    yd[i] = xy[1+i*2];
  }
/*
  #1:  Does interpolant match function at interpolation points?
*/
  ni = nd;

  xi = ( double * ) malloc ( ni * sizeof ( double ) );

  for ( i = 0; i < ni; i++ )
  {
    xi[i] = xd[i];
  }
  yi = chebyshev_interp_1d ( nd, xd, yd, ni, xi );

  int_error = r8vec_norm_affine ( ni, yi, yd ) / ( double ) ( ni );

  printf ( "\n" );
  printf ( "  L2 interpolation error averaged per interpolant node = %g\n", int_error );

  free ( xd );
  free ( xi );
  free ( xy );
  free ( yd );
  free ( yi );

  return;
}
