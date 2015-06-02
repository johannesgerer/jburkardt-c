# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <string.h>

# include "pwl_interp_2d.h"
# include "test_interp_2d.h"
# include "r8lib.h"

int main ( );
void test01 ( int prob, int n );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for PWL_INTERP_2D_PRB.

  Discussion:

    PWL_INTERP_2D_PRB tests the PWL_INTERP_2D library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 October 2012

  Author:

    John Burkardt
*/
{
  int i;
  int n;
  int n_test[5] = { 2, 3, 4, 5, 9 };
  int n_test_num = 5;
  int prob;
  int prob_num;

  timestamp ( );
  printf ( "\n" );
  printf ( "PWL_INTERP_2D_PRB:\n" );
  printf ( "  C version\n" );
  printf ( "  Test the PWL_INTERP_2D library.\n" );
  printf ( "  The R8LIB library is needed.\n" );
  printf ( "  The test needs the TEST_INTERP_2D library.\n" );

  prob_num = f00_num ( );
/*
  Numerical tests.
*/
  for ( prob = 1; prob <= prob_num; prob++ )
  {
    for ( i = 0; i < n_test_num; i++ )
    {
      n = n_test[i];
      test01 ( prob, n );
    }
  }
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "PWL_INTERP_2D_PRB:\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void test01 ( int prob, int n )

/******************************************************************************/
/*
  Purpose:

    PWL_INTERP_2D_TEST01 tests PWL_INTERP_2D.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 October 2012

  Author:

    John Burkardt

  Parameters:

    Input, int PROB, the problem number.

    Input, int N, the grid size in each dimension.
*/
{
  double app_error;
  int i;
  int ij;
  double int_error;
  int j;
  int nd;
  int ni;
  int nxd;
  int nyd;
  double *xd;
  double *xd_1d;
  double *xi;
  double *xi_1d;
  double *yd;
  double *yd_1d;
  double *yi;
  double *yi_1d;
  double *zd;
  double *zdm;
  double *zi;

  nxd = n;
  nyd = n;

  printf ( "\n" );
  printf ( "PWL_INTERP_2D_TEST01:\n" );
  printf ( "  Interpolate data from TEST_INTERP_2D problem #%d\n", prob );
  printf ( "  Using polynomial interpolant of product degree %d x %d\n", nxd, nyd );

  nd = nxd * nyd;
  printf ( "  Number of data points = %d\n", nd );

  xd_1d = r8vec_linspace_new ( nxd, 0.0, 1.0 );
  yd_1d = r8vec_linspace_new ( nyd, 0.0, 1.0 );

  xd = ( double * ) malloc ( nxd * nyd * sizeof ( double ) );
  yd = ( double * ) malloc ( nxd * nyd * sizeof ( double ) );
  zd = ( double * ) malloc ( nxd * nyd * sizeof ( double ) );

  ij = 0;
  for ( j = 0; j < nyd; j++ )
  {
    for ( i = 0; i < nxd; i++ )
    {
      xd[ij] = xd_1d[i];
      yd[ij] = yd_1d[j];
      ij = ij + 1;
    }
  }

  f00_f0 ( prob, nd, xd, yd, zd );

  if ( nd <= 20 )
  {
    r8vec3_print ( nd, xd, yd, zd, "  X, Y, Z data:" );
  }
/*
  #1:  Does interpolant match function at data points?
*/
  ni = nd;
  xi = r8vec_copy_new ( ni, xd );
  yi = r8vec_copy_new ( ni, yd );

  zi = pwl_interp_2d ( nxd, nyd, xd_1d, yd_1d, zd, ni, xi, yi );

  if ( ni <= 20 )
  {
    r8vec3_print ( ni, xi, yi, zi, "  X, Y, Z interpolation:" );
  }

  int_error = r8vec_norm_affine ( ni, zi, zd ) / ( double ) ( ni );

  printf ( "\n" );
  printf ( "  RMS data interpolation error = %g\n", int_error );

  free ( xi );
  free ( yi );
  free ( zi );
/*
  #2:  Does interpolant approximate data at midpoints?
*/
  if ( 1 < nd )
  {
    xi_1d = ( double * ) malloc ( ( nxd - 1 ) * sizeof ( double ) );
    yi_1d = ( double * ) malloc ( ( nyd - 1 ) * sizeof ( double ) );

    for ( i = 0; i < nxd - 1; i++ )
    {
      xi_1d[i] = 0.5 * ( xd_1d[i] + xd_1d[i+1] );
    }
    for ( i = 0; i < nyd - 1; i++ )
    {
      yi_1d[i] = 0.5 * ( yd_1d[i] + yd_1d[i+1] );
    }

    ni = ( nxd - 1 ) * ( nyd - 1 );

    xi = ( double * ) malloc ( ni * sizeof ( double ) );
    yi = ( double * ) malloc ( ni * sizeof ( double ) );
    zdm = ( double * ) malloc ( ni * sizeof ( double ) );

    ij = 0;
    for ( j = 0; j < nyd - 1; j++ )
    {
      for ( i = 0; i < nxd - 1; i++ )
      {
        xi[ij] = xi_1d[i];
        yi[ij] = yi_1d[j];
        ij = ij + 1;
      }
    }

    f00_f0 ( prob, ni, xi, yi, zdm );

    zi = pwl_interp_2d ( nxd, nyd, xd_1d, yd_1d, zd, ni, xi, yi );

    app_error = r8vec_norm_affine ( ni, zi, zdm ) / ( double ) ( ni );

    printf ( "  RMS data approximation error = %g\n", app_error );

    free ( xi );
    free ( xi_1d );
    free ( yi );
    free ( yi_1d );
    free ( zdm );
    free ( zi );
  }

  free ( xd );
  free ( xd_1d );
  free ( yd );
  free ( yd_1d );
  free ( zd );

  return;
}
