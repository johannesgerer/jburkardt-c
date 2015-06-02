# include <stdlib.h>
# include <stdio.h>
# include <string.h>

# include "nearest_interp_1d.h"
# include "test_interp.h"
# include "r8lib.h"

int main ( );
void nearest_interp_1d_test01 ( int prob, int ni );
void nearest_interp_1d_test02 ( int prob );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for NEAREST_INTERP_1D_PRB.

  Discussion:

    NEAREST_INTERP_1D_PRB tests the NEAREST_INTERP_1D library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    28 May 2013

  Author:

   John Burkardt
*/
{
  int ni;
  int prob;
  int prob_num;

  timestamp ( );
  printf ( "\n" );
  printf ( "NEAREST_INTERP_1D_PRB:\n" );
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

  for ( prob = 1; prob <= prob_num; prob++ )
  {
    nearest_interp_1d_test02 ( prob );
  }
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "NEAREST_INTERP_1D_PRB:\n" );
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
/******************************************************************************/

void nearest_interp_1d_test02 ( int prob )

/******************************************************************************/
/*
  Purpose:

    NEAREST_INTERP_1D_TEST02 tests NEAREST_INTERP_1D.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    28 May 2013

  Author:

    John Burkardt

  Parameters:

    Input, int PROB, the index of the problem.
*/
{
  char command_filename[255];
  FILE *command_unit;
  char data_filename[255];
  FILE *data_unit;
  int i;
  double interp_error;
  char interp_filename[255];
  FILE *interp_unit;
  int j;
  int nd;
  int ni;
  char output_filename[255];
  char title[255];
  double *xd;
  double *xi;
  double xmax;
  double xmin;
  double *xy;
  double *yd;
  double *yi;

  printf ( "\n" );
  printf ( "NEAREST_INTERP_1D_TEST02:\n" );
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
  yi = nearest_interp_1d ( nd, xd, yd, ni, xi );

  interp_error = r8vec_diff_norm ( ni, yi, yd ) / ( double ) ( ni );

  printf ( "\n" );
  printf ( "  Node-averaged L2 interpolation error = %g\n", interp_error );
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
  free ( xi );
  free ( yi );

  ni = 501;
  xmin = r8vec_min ( nd, xd );
  xmax = r8vec_max ( nd, xd );
  xi = r8vec_linspace_new ( ni, xmin, xmax );
  yi = nearest_interp_1d ( nd, xd, yd, ni, xi );

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
  fprintf ( command_unit, "set title 'Data versus Nearest Neighbor Interpolant'\n" );
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
