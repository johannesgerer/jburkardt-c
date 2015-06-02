# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>

# include "rbf_interp_1d.h"
# include "test_interp.h"
# include "r8lib.h"

int main ( );
void test01 ( int prob, void phi ( int n, double r[], double r0, double v[] ), 
  char *phi_name, double r0 );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for RBF_INTERP_1D_PRB.

  Discussion:

    RBF_INTERP_1D_PRB tests the RBF_INTERP_1D library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 October 2012

  Author:

    John Burkardt
*/
{
  int i;
  int nd;
  int prob;
  int prob_num;
  double r0;
  double *xd;
  double xmax;
  double xmin;
  double *xy;

  timestamp ( );
  printf ( "\n" );
  printf ( "RBF_INTERP_1D_PRB:\n" );
  printf ( "  C version\n" );
  printf ( "  Test the RBF_INTERP_1D library.\n" );
  printf ( "  The R8LIB library is needed.\n" );
  printf ( "  The test requires the TEST_INTERP library.\n" );

  prob_num = p00_prob_num ( );

  for ( prob = 1; prob <= prob_num; prob++ )
  {
/*
  Determine an appropriate value of R0, the spacing parameter.
*/
    nd = p00_data_num ( prob );
    xy = p00_data ( prob, 2, nd );
    xd = ( double * ) malloc ( nd * sizeof ( double ) );
    for ( i = 0; i < nd; i++ )
    {
      xd[i] = xy[0+i*2];
    }
    xmax = r8vec_max ( nd, xd );
    xmin = r8vec_min ( nd, xd );
    r0 = ( xmax - xmin ) / ( double ) ( nd - 1 );
    free ( xd );
    free ( xy );
    test01 ( prob, phi1, "phi1", r0 );
    test01 ( prob, phi2, "phi2", r0 );
    test01 ( prob, phi3, "phi3", r0 );
    test01 ( prob, phi4, "phi4", r0 );
  }
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "RBF_INTERP_1D_PRB:\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void test01 ( int prob, void phi ( int n, double r[], double r0, double v[] ), 
  char *phi_name, double r0 )

/******************************************************************************/
/*
  Purpose:

    TEST01 tests RBF_INTERP_1D.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 October 2012

  Author:

    John Burkardt

  Parameters:

    Input, int PROB, the index of the problem.

    Input, double PHI ( int n, double r[], double r0, double v[] ), 
    the name of the radial basis function.

    Input, char *PHI_NAME, the name of the radial basis function.

    Input, double R0, the scale factor.  Typically, this might be
    a small multiple of the average distance between points.
*/
{
  int debug = 0;
  int i;
  double int_error;
  double ld;
  double li;
  int m;
  int nd;
  int ni;
  double *w;
  double *xd;
  double *xi;
  double xmax;
  double xmin;
  double *xy;
  double *yd;
  double *yi;
  double ymax;
  double ymin;

  printf ( "\n" );
  printf ( "TEST01:\n" );
  printf ( "  Interpolate data from TEST_INTERP problem #%d\n", prob );
  printf ( "  using radial basis function \"%s\".\n", phi_name );
  printf ( "  Scale factor R0 = %g\n", r0 );

  nd = p00_data_num ( prob );
  printf ( "  Number of data points = %d\n", nd );

  xy = p00_data ( prob, 2, nd );
  
  if ( debug )
  {
    r8mat_transpose_print ( 2, nd, xy, "  Data array:" );
  }

  xd = ( double * ) malloc ( nd * sizeof ( double ) );
  yd = ( double * ) malloc ( nd * sizeof ( double ) );
  for ( i = 0; i < nd; i++ )
  {
    xd[i] = xy[0+i*2];
    yd[i] = xy[1+i*2];
  }
  m = 1;
  w = rbf_weight ( m, nd, xd, r0, phi, yd );
/*
  #1:  Does interpolant match function at interpolation points?
*/
  ni = nd;
  xi = r8vec_copy_new ( ni, xd );
  yi = rbf_interp ( m, nd, xd, r0, phi, w, ni, xi );
  int_error = r8vec_norm_affine ( ni, yi, yd ) / ( double ) ( ni );

  printf ( "\n" );
  printf ( "  L2 interpolation error averaged per interpolant node = %g\n", int_error );

  free ( xi );
  free ( yi );
/*
  #2: Compare estimated curve length to piecewise linear (minimal) curve length.
  Assume data is sorted, and normalize X and Y dimensions by (XMAX-XMIN) and
  (YMAX-YMIN).
*/
  xmax = r8vec_max ( nd, xd );
  xmin = r8vec_min ( nd, xd );
  ymax = r8vec_max ( nd, yd );
  ymin = r8vec_min ( nd, yd );

  ni = 501;
  xi = r8vec_linspace_new ( ni, xmin, xmax );
  yi = rbf_interp ( m, nd, xd, r0, phi, w, ni, xi );

  ld = 0.0;
  for ( i = 0; i < nd - 1; i++ )
  {
    ld = ld + sqrt ( pow ( ( xd[i+1] - xd[i] ) / ( xmax - xmin ), 2 )
                   + pow ( ( yd[i+1] - yd[i] ) / ( ymax - ymin ), 2 ) );
  }
  li = 0.0;
  for ( i = 0; i < ni - 1; i++ )
  {
    li = li + sqrt ( pow ( ( xi[i+1] - xi[i] ) / ( xmax - xmin ), 2 )
                   + pow ( ( yi[i+1] - yi[i] ) / ( ymax - ymin ), 2 ) );
  }

  printf ( "  Normalized length of piecewise linear interpolant = %g\n", ld );
  printf ( "  Normalized length of polynomial interpolant       = %g\n", li );

  free ( w );
  free ( xd );
  free ( xi );
  free ( xy );
  free ( yd );
  free ( yi );

  return;
}
