# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <string.h>

# include "lagrange_interp_2d.h"
# include "test_interp_2d.h"
# include "r8lib.h"

int main ( );
void test01 ( int prob, int m );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    LAGRANGE_INTERP_2D_TEST tests LAGRANGE_INTERP_2D.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 September 2012

  Author:

    John Burkardt
*/
{
  int i;
  int m;
  int m_test[5] = { 1, 2, 3, 4, 8 };
  int m_test_num = 5;
  int prob;
  int prob_num;

  timestamp ( );
  printf ( "\n" );
  printf ( "LAGRANGE_INTERP_2D_TEST:\n" );
  printf ( "  C version\n" );
  printf ( "  Test the LAGRANGE_INTERP_2D library.\n" );
  printf ( "  The R8LIB library is needed.\n" );
  printf ( "  This test also needs the TEST_INTERP_2D library.\n" );

  prob_num = f00_num ( );
/*
  Numerical tests.
*/
  for ( prob = 1; prob <= prob_num; prob++ )
  {
    for ( i = 0; i < m_test_num; i++ )
    {
      m = m_test[i];
      test01 ( prob, m );
    }
  }
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "LAGRANGE_INTERP_2D_TEST:\n" );
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

    LAGRANGE_INTERP_2D_TEST01 tests LAGRANGE_INTERP_2D.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 September 2012

  Author:

    John Burkardt

  Parameters:

    Input, int PROB, the problem number.

    Input, int M, the polynomial degree in each dimension.
*/
{
  double app_error;
  int i;
  int ij;
  double int_error;
  int j;
  int mx;
  int my;
  int nd;
  int ni;
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

  mx = m;
  my = m;

  printf ( "\n" );
  printf ( "LAGRANGE_INTERP_2D_TEST01:\n" );
  printf ( "  Interpolate data from TEST_INTERP_2D problem #%d\n", prob );
  printf ( "  Using polynomial interpolant of product degree %d x %d\n", mx, my );

  nd = ( mx + 1 ) * ( my + 1 );
  printf ( "  Number of data points = %d\n", nd );

  xd_1d = r8vec_chebyspace_new ( mx + 1, 0.0, 1.0 );
  yd_1d = r8vec_chebyspace_new ( my + 1, 0.0, 1.0 );

  xd = ( double * ) malloc ( (mx+1)*(my+1) * sizeof ( double ) );
  yd = ( double * ) malloc ( (mx+1)*(my+1) * sizeof ( double ) );
  zd = ( double * ) malloc ( (mx+1)*(my+1) * sizeof ( double ) );

  ij = 0;
  for ( j = 0; j < my + 1; j++ )
  {
    for ( i = 0; i < mx + 1; i++ )
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

  xi = ( double * ) malloc ( ni * sizeof ( double ) );
  yi = ( double * ) malloc ( ni * sizeof ( double ) );

  for ( i = 0; i < ni; i++ )
  {
    xi[i] = xd[i];
    yi[i] = yd[i];
  }

  zi = lagrange_interp_2d ( mx, my, xd_1d, yd_1d, zd, ni, xi, yi );

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
    xi_1d = ( double * ) malloc ( mx * sizeof ( double ) );
    yi_1d = ( double * ) malloc ( my * sizeof ( double ) );

    for ( i = 0; i < mx; i++ )
    {
      xi_1d[i] = 0.5 * ( xd_1d[i] + xd_1d[i+1] );
    }
    for ( i = 0; i < my; i++ )
    {
      yi_1d[i] = 0.5 * ( yd_1d[i] + yd_1d[i+1] );
    }

    ni = mx * my;

    xi = ( double * ) malloc ( ni * sizeof ( double ) );
    yi = ( double * ) malloc ( ni * sizeof ( double ) );
    zdm = ( double * ) malloc ( ni * sizeof ( double ) );
    
    ij = 0;
    for ( j = 0; j < my; j++ )
    {
      for ( i = 0; i < mx; i++ )
      {
        xi[ij] = xi_1d[i];
        yi[ij] = yi_1d[j];
        ij = ij + 1;
      }
    }

    f00_f0 ( prob, ni, xi, yi, zdm );

    zi = lagrange_interp_2d ( mx, my, xd_1d, yd_1d, zd, ni, xi, yi );

    app_error = r8vec_norm_affine ( ni, zi, zdm ) / ( double ) ( ni );

    printf ( "\n" );
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
