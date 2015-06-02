# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <string.h>

# include "naca.h"

int main ( );
void test01 ( );
void test02 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for NACA_PRB.

  Discussion:

    NACA_PRB tests the NACA library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    22 May 2014

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "NACA_PRB:\n" );
  printf ( "  C version\n" );
  printf ( "  Test the NACA library.\n" );

  test01 ( );
  test02 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "NACA_PRB:\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void test01 ( )

/******************************************************************************/
/*
  Purpose:

    TEST01 tests NACA4_SYMMETRIC.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    22 May 2014

  Author:

    John Burkardt
*/
{
  double c;
  char command_filename[] = "symmetric_commands.txt";
  FILE *command_unit;
  char data_filename[] = "symmetric_data.txt";
  int i;
  int n = 51;
  double ratio;
  double t;
  double *x;
  double x_max;
  double x_min;
  double *xy;
  double *y;
  double y_max;
  double y_min;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  NACA4_SYMMETRIC evaluates y(x) for a NACA\n" );
  printf ( "  symmetric airfoil defined by a 4-digit code.\n" );

  c = 10.0;
  t = 0.15;
  x = r8vec_linspace_new ( n, 0.0, c );
  y = naca4_symmetric ( t, c, n, x );
/*
  Reorganize data into a single object.
*/
  xy = ( double * ) malloc ( 2 * 2 * n * sizeof ( double ) );

  for ( i = 0; i < n; i++ )
  {
    xy[0+i*2] = x[i];
    xy[1+i*2] = -y[i];
  }
  for ( i = 0; i < n; i++ )
  {
    xy[0+(n+i)*2] = x[n-1-i];
    xy[1+(n+i)*2] = y[n-1-i];
  }
/*
  Determine size ratio.
*/
  x_min = r8vec_min ( n, x );
  x_max = r8vec_max ( n, x );
  y_max = r8vec_max ( n, y );
  y_min = - y_max;
  ratio = ( y_max - y_min ) / ( x_max - x_min );
/*
  Save data to a file.
*/
  r8mat_write ( data_filename, 2, 2 * n, xy );
  printf ( "  Data saved in file '%s'\n", data_filename );
/*
  Create the command file.
*/
  command_unit = fopen ( command_filename, "wt" );
  fprintf ( command_unit, "set term png\n" );
  fprintf ( command_unit, "set grid\n" );
  fprintf ( command_unit, "set style data lines\n" );
  fprintf ( command_unit, "set size ratio %g\n", ratio );
  fprintf ( command_unit, "set timestamp\n" );
  fprintf ( command_unit, "unset key\n" );
  fprintf ( command_unit, "set output 'symmetric.png'\n" );
  fprintf ( command_unit, "set xlabel '<---X--->'\n" );
  fprintf ( command_unit, "set ylabel '<---Y--->'\n" );
  fprintf ( command_unit, "set title 'NACA Symmetric Airfoil'\n" );
  fprintf ( command_unit, "plot '%s' using 1:2 with lines lw 3\n", 
    data_filename );
  fprintf ( command_unit, "quit\n" );

  fclose ( command_unit );

  printf ( "  Created command file '%s'\n", command_filename );

  free ( x );
  free ( y );

  return;
}
/******************************************************************************/

void test02 ( )

/******************************************************************************/
/*
  Purpose:

    TEST02 tests NACA4_CAMBERED.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    21 May 2014

  Author:

    John Burkardt
*/
{
  double c;
  char command_filename[] = "cambered_commands.txt";
  FILE *command_unit;
  char data_filename[] = "cambered_data.txt";
  int i;
  double m;
  int n = 51;
  double p;
  double ratio;
  double t;
  double x_max;
  double x_min;
  double *xc;
  double *xl;
  double *xu;
  double *xy;
  double y_max;
  double y_min;
  double *yl;
  double *yu;

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  NACA4_CAMBERED evaluates (xu,yu) and (xl,yl) for a NACA\n" );
  printf ( "  cambered airfoil defined by a 4-digit code.\n" );

  m = 0.02;
  p = 0.4;
  t = 0.12;
  c = 10.0;

  xc = r8vec_linspace_new ( n, 0.0, c );

  xu = ( double * ) malloc ( n * sizeof ( double ) );
  xl = ( double * ) malloc ( n * sizeof ( double ) );
  yu = ( double * ) malloc ( n * sizeof ( double ) );
  yl = ( double * ) malloc ( n * sizeof ( double ) );

  naca4_cambered ( m, p, t, c, n, xc, xu, yu, xl, yl );
/*
  Reorganize data into a single object.
*/
  xy = ( double * ) malloc ( 2 * 2 * n * sizeof ( double ) );

  for ( i = 0; i < n; i++ )
  {
    xy[0+i*2] = xl[i];
    xy[1+i*2] = yl[i];
  }
  for ( i = 0; i < n; i++ )
  {
    xy[0+(n+i)*2] = xu[n-1-i];
    xy[1+(n+i)*2] = yu[n-1-i];
  }
/*
  Determine size ratio.
*/
  x_min = r8_min ( r8vec_min ( n, xl ), r8vec_min ( n, xu ) );
  x_max = r8_max ( r8vec_max ( n, xl ), r8vec_max ( n, xu ) );
  y_min = r8_min ( r8vec_min ( n, yl ), r8vec_min ( n, yu ) );
  y_max = r8_max ( r8vec_max ( n, yl ), r8vec_max ( n, yu ) );
  ratio = ( y_max - y_min ) / ( x_max - x_min );
/*
  Save data to a file.
*/
  r8mat_write ( data_filename, 2, 2 * n, xy );
  printf ( "  Data saved in file '%s'\n", data_filename );
/*
  Create the command file.
*/
  command_unit = fopen ( command_filename, "wt" );
  fprintf ( command_unit, "set term png\n" );
  fprintf ( command_unit, "set grid\n" );
  fprintf ( command_unit, "set style data lines\n" );
  fprintf ( command_unit, "set size ratio %g\n", ratio );
  fprintf ( command_unit, "set timestamp\n" );
  fprintf ( command_unit, "unset key\n" );
  fprintf ( command_unit, "set output 'cambered.png'\n" );
  fprintf ( command_unit, "set xlabel '<---X--->'\n" );
  fprintf ( command_unit, "set ylabel '<---Y--->'\n" );
  fprintf ( command_unit, "set title 'NACA Cambered Airfoil'\n" );
  fprintf ( command_unit, "plot '%s' using 1:2 with lines lw 3\n", 
    data_filename );
  fprintf ( command_unit, "quit\n" );

  fclose ( command_unit );

  printf ( "  Created command file '%s'\n", command_filename );

  free ( xc );
  free ( xl );
  free ( xu );
  free ( xy );
  free ( yl );
  free ( yu );

  return;
}

