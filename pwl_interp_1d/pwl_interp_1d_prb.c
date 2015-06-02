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

    MAIN is the main program for PWL_INTERP_1D_PRB.

  Discussion:

    PWL_INTERP_1D_PRB tests the PWL_INTERP_1D library.

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
  printf ( "PWL_INTERP_1D_PRB:\n" );
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
  printf ( "PWL_INTERP_1D_PRB:\n" );
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

    31 May 2012

  Author:

    John Burkardt

  Parameters:

    Input, int PROB, the problem index.
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
  double *xy;
  double xmax;
  double xmin;
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

  interp_error = r8vec_diff_norm ( ni, yi, yd ) / ( double ) ( ni );

  printf ( "\n" );
  printf ( "  L2 interpolation error averaged per interpolant node = %g\n", interp_error );
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
  yi = pwl_interp_1d ( nd, xd, yd, ni, xi );

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
  fprintf ( command_unit, "set title 'Data versus piecewise linear interpolant'\n" );
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
