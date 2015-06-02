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
void test02 ( int prob );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for VANDERMONDE_INTERP_1D_PRB.

  Discussion:

    VANDERMONDE_INTERP_1D_PRB tests the VANDERMONDE_INTERP_1D library.

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
  printf ( "VANDERMONDE_INTERP_1D_PRB:\n" );
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

  for ( prob = 1; prob <= prob_num; prob++ )
  {
    test02 ( prob );
  }
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "VANDERMONDE_INTERP_1D_PRB:\n" );
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
/******************************************************************************/

void test02 ( int prob )

/******************************************************************************/
/*
  Purpose:

    TEST02 tests VANDERMONDE_INTERP_1D_MATRIX.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 June 2012

  Author:

    John Burkardt

  Parameters:

    Input, int PROB, the problem index.
*/
{
  double *a;
  double *c;
  char command_filename[255];
  FILE *command_unit;
  char data_filename[255];
  FILE *data_unit;
  int i;
  char interp_filename[255];
  FILE *interp_unit;
  int j;
  int nd;
  int ni;
  char output_filename[255];
  char title[255];
  double *xd;
  double *xi;
  double *xy;
  double xmax;
  double xmin;
  double *yd;
  double *yi;

  printf ( "\n" );
  printf ( "TEST02:\n" );
  printf ( "  VANDERMONDE_INTERP_1D_MATRIX sets the Vandermonde linear system for the\n" );
  printf ( "  interpolating polynomial.\n" );
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
  Set up the Vandermonde matrix A.
*/
  a = vandermonde_interp_1d_matrix ( nd, xd );
/*
  Solve the linear system for the polynomial coefficients C.
*/
  c = qr_solve ( nd, nd, a, yd );
/*
  Create data file.
*/
  sprintf ( data_filename, "data%02d.txt", prob );
  data_unit = fopen ( data_filename, "wt" );
  for ( j = 0; j < nd; j++ )
  {
    fprintf ( data_unit, "  %14g  %14g\n", xd[j], yd[j] );
  }
  fclose ( data_unit );
  printf ( "\n" );
  printf ( "  Created graphics data file \"%s\".\n", data_filename );
/*
  Create interp file.
*/
  ni = 501;
  xmin = r8vec_min ( nd, xd );
  xmax = r8vec_max ( nd, xd );
  xi = r8vec_linspace_new ( ni, xmin, xmax );
  yi = r8poly_value ( nd - 1, c, ni, xi );

  sprintf ( interp_filename, "interp%02d.txt", prob );
  interp_unit = fopen ( interp_filename, "wt" );
  for ( j = 0; j < ni; j++ )
  {
    fprintf ( interp_unit, "  %g  %g\n", xi[j], yi[j] );
  }
  fclose ( interp_unit );
  printf ( "  Created graphics interp file \"%s\".\n", interp_filename );
/*
  Plot the data and the interpolant.
*/
  sprintf ( command_filename, "commands%02d.txt", prob );
  command_unit = fopen ( command_filename, "wt" );

  sprintf ( output_filename, "plot%02d.png", prob );

  fprintf ( command_unit, "# %s\n", command_filename );
  fprintf ( command_unit, "#\n" );
  fprintf ( command_unit, "# Usage:\n" );
  fprintf ( command_unit, "#  gnuplot < %s\n", command_filename );
  fprintf ( command_unit, "#\n" );
  fprintf ( command_unit, "set term png\n" );
  fprintf ( command_unit, "set output '%s'\n", output_filename );
  fprintf ( command_unit, "set xlabel '<---X--->'\n" );
  fprintf ( command_unit, "set ylabel '<---Y--->'\n" );
  fprintf ( command_unit, "set title 'Data versus Vandermonde polynomial interpolant'\n" );
  fprintf ( command_unit, "set grid\n" );
  fprintf ( command_unit, "set style data lines\n" );
  fprintf ( command_unit, "plot '%s' using 1:2 with points pt 7 ps 2 lc rgb 'blue',\\\n",
    data_filename );
  fprintf ( command_unit, "     '%s' using 1:2 lw 3 linecolor rgb 'red'\n", 
    interp_filename );

  fclose ( command_unit );
  printf ( "  Created graphics command file \"%s\".\n", command_filename );
/*
  Free memory.
*/
  free ( xd );
  free ( xi );
  free ( xy );
  free ( yd );
  free ( yi );

  return;
}
