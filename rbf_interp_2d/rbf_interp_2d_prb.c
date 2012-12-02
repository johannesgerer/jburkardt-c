# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>

# include "rbf_interp_2d.h"
# include "r8lib.h"
# include "test_interp_2d.h"

int main ( );
void test01 ( int prob, int g, 
  void phi ( int n, double r[], double r0, double v[] ), char *phi_name );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    RBF_INTERP_2D_TEST tests RBF_INTERP_2D.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 October 2012

  Author:

    John Burkardt
*/
{
  int g;
  int prob;
  int prob_num;

  timestamp ( );
  printf ( "\n" );
  printf ( "RBF_INTERP_2D_TEST:\n" );
  printf ( "  C version\n" );
  printf ( "  Test the RBF_INTERP_2D library.\n" );
  printf ( "  The R8LIB library is required.\n" );
  printf ( "  This test also needs the TEST_INTERP_2D library.\n" );

  prob_num = f00_num ( );
  g = 1;

  for ( prob = 1; prob <= prob_num; prob++ )
  {
    test01 ( prob, g, phi1, "phi1" );
    test01 ( prob, g, phi2, "phi2" );
    test01 ( prob, g, phi3, "phi3" );
    test01 ( prob, g, phi4, "phi4" );
  }
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "RBF_INTERP_2D_TEST:\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void test01 ( int prob, int g, 
  void phi ( int n, double r[], double r0, double v[] ), char *phi_name )

/******************************************************************************/
/*
  Purpose:

    RBF_INTERP_2D_TEST01 tests RBF_INTERP_2D.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 October 2012

  Author:

    John Burkardt

  Parameters:

    Input, int PROB, the index of the problem.

    Input, int G, the index of the grid.

    Input, void PHI ( int n, double r[], double r0, double v[] ), the 
    radial basis function.

    Input, char *PHI_NAME, the name of the radial basis function.
*/
{
  int debug = 0;
  double e;
  int i;
  double int_error;
  int m;
  int nd;
  int ni;
  double r0;
  double volume;
  double *w;
  double *xd;
  double xmax;
  double xmin;
  double *xyd;
  double *xyi;
  double *yd;
  double ymax;
  double ymin;
  double *zd;
  double *zi;

  printf ( "\n" );
  printf ( "RBF_INTERP_2D_TEST01:\n" );
  printf ( "  Interpolate data from TEST_INTERP_2D problem #%d\n", prob );
  printf ( "  using grid #%d\n", g );
  printf ( "  using radial basis function \"%s\".\n", phi_name );

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

  m = 2;
  xyd = ( double * ) malloc ( 2 * nd * sizeof ( double ) );

  for ( i = 0; i < nd; i++ )
  {
    xyd[0+i*2] = xd[i];
    xyd[1+i*2] = yd[i];
  }

  xmax = r8vec_max ( nd, xd );
  xmin = r8vec_min ( nd, xd );
  ymax = r8vec_max ( nd, yd );
  ymin = r8vec_min ( nd, yd );
  volume = ( xmax - xmin ) * ( ymax - ymin );

  e = 1.0 / ( double ) ( m );
  r0 = pow ( volume / nd, e );

  printf ( "  Setting R0 = %g\n", r0 );

  w = rbf_weight ( m, nd, xyd, r0, phi, zd );
/*
  #1:  Does interpolant match function at interpolation points?
*/
  ni = nd;
  xyi = r8mat_copy_new ( 2, ni, xyd );

  zi = rbf_interp ( m, nd, xyd, r0, phi, w, ni, xyi );

  int_error = r8vec_norm_affine ( ni, zi, zd ) / ( double ) ( ni );

  printf ( "\n" );
  printf ( "  L2 interpolation error averaged per interpolant node = %g\n", int_error );

  free ( w );
  free ( xd );
  free ( xyd );
  free ( xyi );
  free ( yd );
  free ( zd );
  free ( zi );

  return;
}
