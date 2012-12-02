# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "vandermonde_interp_1d.h"
# include "condition.h"
# include "qr_solve.h"
# include "test_interp.h"
# include "r8lib.h"

int main ( );
void test01 ( int prob );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    VANDERMONDE_INTERP_1D_TEST tests VANDERMONDE_INTERP_1D.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    06 October 2012

  Author:

    John Burkardt
*/
{
  int prob;
  int prob_num;

  timestamp ( );
  printf ( "\n" );
  printf ( "VANDERMONDE_INTERP_1D_TEST:\n" );
  printf ( "  C version\n" );
  printf ( "  Test the VANDERMONDE_INTERP_1D library.\n" );
  printf ( "  The QR_SOLVE library is needed.\n" );
  printf ( "  The R8LIB library is needed.\n" );
  printf ( "  This test needs the CONDITION library.\n" );
  printf ( "  This test needs the TEST_INTERP library.\n" );

  prob_num = p00_prob_num ( );
  for ( prob = 1; prob <= prob_num; prob++ )
  {
    test01 ( prob );
  }
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "VANDERMONDE_INTERP_1D_TEST:\n" );
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

    TEST01 tests VANDERMONDE_INTERP_1D_MATRIX.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    06 October 2012

  Author:

    John Burkardt
*/
{
  double *a;
  int debug = 0;
  double *c;
  double condition;
  int i;
  double int_error;
  double ld;
  double li;
  int m;
  int nd;
  int ni;
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
/*
  Choose the degree of the polynomial to be ND - 1.
*/
  m = nd - 1;
/*
  Compute Vandermonde matrix and get condition number.
*/
  a = vandermonde_interp_1d_matrix ( nd, xd );

  condition = condition_hager ( nd, a );

  printf ( "\n" );
  printf ( "  Condition of Vandermonde matrix is %g\n", condition );
/*
  Solve linear system.
*/
  c = qr_solve ( nd, nd, a, yd );
/*
  #1:  Does interpolant match function at interpolation points?
*/
  ni = nd;
  xi = r8vec_copy_new ( ni, xd );
  yi = r8poly_value ( m, c, ni, xi );

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
  xmin = r8vec_min ( nd, xd );
  xmax = r8vec_max ( nd, xd );
  ymin = r8vec_min ( nd, yd );
  ymax = r8vec_max ( nd, yd );

  ni = 501;
  xi = r8vec_linspace_new ( ni, xmin, xmax );
  yi = r8poly_value ( m, c, ni, xi );

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

  printf ( "\n" );
  printf ( "  Normalized length of piecewise linear interpolant = %g\n", ld );
  printf ( "  Normalized length of polynomial interpolant       = %g\n", li );

  free ( a );
  free ( c );
  free ( xd );
  free ( xi );
  free ( xy );
  free ( yd );
  free ( yi );

  return;
}
