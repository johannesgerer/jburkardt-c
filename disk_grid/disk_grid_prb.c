# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>

# include "disk_grid.h"

int main ( );
void disk_grid_test01 ( );
void disk_grid_test02 ( );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for DISK_GRID_PRB.

  Discussion:

    DISK_GRID_PRB tests the DISK_GRID library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 October 2013

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "DISK_GRID_PRB:\n" );
  printf ( "  C version\n" );
  printf ( "  Test the DISK_GRID library.\n" );

  disk_grid_test01 ( );
  disk_grid_test02 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "DISK_GRID_PRB:\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void disk_grid_test01 ( )

/******************************************************************************/
/*
  Purpose:

    DISK_GRID_TEST01 tests DISK_GRID.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 October 2013

  Author:

    John Burkardt
*/
{
  double c[2];
  double *cg;
  char boundary_filename[] = "disk_grid_test01_boundary.txt";
  FILE *boundary_unit;
  char command_filename[] = "disk_grid_test01_commands.txt";
  FILE *command_unit;
  char data_filename[] = "disk_grid_test01_data.txt";
  FILE *data_unit;
  char *filename = "disk_grid_test01.xy";
  int i;
  int n;
  int ng;
  const double pi = 3.141592653589793;
  double r;
  double t;

  printf ( "\n" );
  printf ( "TEST01:\n" );
  printf ( "  DISK_GRID can define a grid of points\n" );
  printf ( "  with N+1 points on a horizontal or vertical radius,\n" );
  printf ( "  based on any disk.\n" );

  n = 20;
  r = 2.0;
  c[0] = 1.0;
  c[1] = 5.0;

  printf ( "\n" );
  printf ( "  We use N = %d\n", n );
  printf ( "  Radius R = %g\n", r );
  printf ( "  Center C = (%g,%g)\n", c[0], c[1] );

  ng = disk_grid_count ( n, r, c );

  printf ( "\n" );
  printf ( "  Number of grid points will be %d\n", ng );

  cg = disk_grid ( n, r, c, ng );

  r82vec_print_part ( ng, cg, 20, "  Part of the grid point array:" );
/*
  Write the grid coordinates to a file.
*/
  r8mat_write ( filename, 2, ng, cg );

  printf ( "\n" );
  printf ( "  Data written to the file \"%s\".\n", filename );
/*
  Create graphics data files.
*/
  boundary_unit = fopen ( boundary_filename, "wt" );
  for ( i = 0; i <= 50; i++ )
  {
    t = 2.0 * pi * ( double ) ( i ) / 50.0;
    fprintf ( boundary_unit, "  %14.6g  %14.6g\n",
      c[0] + r * cos ( t ), c[1] + r * sin ( t ) );
  }
  fclose ( boundary_unit );
  printf ( "\n" );
  printf ( "  Created boundary file \"%s\".\n", boundary_filename );

  data_unit = fopen ( data_filename, "wt" );
  for ( i = 0; i < ng; i++ )
  {
    fprintf ( data_unit, "  %14.6g  %14.6g\n", cg[0+i*2], cg[1+i*2] );
  }
  fclose ( data_unit );
  printf ( "\n" );
  printf ( "  Created data file \"%s\"\n", data_filename );
/*
  Create graphics command file.
*/
  command_unit = fopen ( command_filename, "wt" );
  fprintf ( command_unit, "# %s\n", command_filename );
  fprintf ( command_unit, "#\n" );
  fprintf ( command_unit, "# Usage:\n" );
  fprintf ( command_unit, "#  gnuplot < %s\n", command_filename );
  fprintf ( command_unit, "#\n" );
  fprintf ( command_unit, "set term png\n" );
  fprintf ( command_unit, "set output 'disk_grid_test01.png'\n" );
  fprintf ( command_unit, "set xlabel '<--- X --->'\n" );
  fprintf ( command_unit, "set ylabel '<--- Y --->'\n" );
  fprintf ( command_unit, "set title 'Disk Grid'\n" );
  fprintf ( command_unit, "set grid\n" );
  fprintf ( command_unit, "set key off\n" );
  fprintf ( command_unit, "set size ratio -1\n" );
  fprintf ( command_unit, "set style data lines\n" );
  fprintf ( command_unit, "plot '%s' using 1:2 with points lt 3 pt 3,\\\n", 
    data_filename );
  fprintf ( command_unit, "    '%s' using 1:2 lw 3 linecolor rgb 'black'\n", 
    boundary_filename );
  fprintf ( command_unit, "quit\n" );
  fclose ( command_unit );

  printf ( "  Created command file \"%s\"\n", command_filename );

  free ( cg );

  return;
}
/******************************************************************************/

void disk_grid_test02 ( )

/******************************************************************************/
/*
  Purpose:

    DISK_GRID_TEST02 tests DISK_GRID_FIBONACCI.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 October 2013

  Author:

    John Burkardt
*/
{
  double c[2];
  char boundary_filename[] = "disk_grid_test02_boundary.txt";
  FILE *boundary_unit;
  char command_filename[] = "disk_grid_test02_commands.txt";
  FILE *command_unit;
  char data_filename[] = "disk_grid_test02_data.txt";
  FILE *data_unit;
  char *filename = "disk_grid_test02.xy";
  double *g;
  int i;
  int n;
  const double pi = 3.141592653589793;
  double r;
  double t;

  printf ( "\n" );
  printf ( "TEST02:\n" );
  printf ( "  DISK_GRID can define a grid of N points\n" );
  printf ( "  based on a Fibonacci spiral inside a disk.\n" );

  n = 1000;
  r = 2.0;
  c[0] = 1.0;
  c[1] = 5.0;

  printf ( "\n" );
  printf ( "  We use N = %d\n", n );
  printf ( "  Radius R = %g\n", r );
  printf ( "  Center C = (%g,%g)\n", c[0], c[1] );

  g = disk_grid_fibonacci ( n, r, c );

  r82vec_print_part ( n, g, 20, "  Part of the grid point array:" );
/*
  Write the grid coordinates to a file.
*/
  r8mat_write ( filename, 2, n, g );

  printf ( "\n" );
  printf ( "  Coordinate data written to the file \"%s\".\n", filename );
/*
  Create graphics data files.
*/
  boundary_unit = fopen ( boundary_filename, "wt" );
  for ( i = 0; i <= 50; i++ )
  {
    t = 2.0 * pi * ( double ) ( i ) / 50.0;
    fprintf ( boundary_unit, "  %14.6g  %14.6g\n",
      c[0] + r * cos ( t ), c[1] + r * sin ( t ) );
  }
  fclose ( boundary_unit );
  printf ( "\n" );
  printf ( "  Created boundary file \"%s\".\n", boundary_filename );

  data_unit = fopen ( data_filename, "wt" );
  for ( i = 0; i < n; i++ )
  {
    fprintf ( data_unit, "  %14.6g  %14.6g\n", g[0+i*2], g[1+i*2] );
  }
  fclose ( data_unit );
  printf ( "\n" );
  printf ( "  Created data file \"%s\"\n", data_filename );
/*
  Create graphics command file.
*/
  command_unit = fopen ( command_filename, "wt" );
  fprintf ( command_unit, "# %s\n", command_filename );
  fprintf ( command_unit, "#\n" );
  fprintf ( command_unit, "# Usage:\n" );
  fprintf ( command_unit, "#  gnuplot < %s\n", command_filename );
  fprintf ( command_unit, "#\n" );
  fprintf ( command_unit, "set term png\n" );
  fprintf ( command_unit, "set output 'disk_grid_test02.png'\n" );
  fprintf ( command_unit, "set xlabel '<--- X --->'\n" );
  fprintf ( command_unit, "set ylabel '<--- Y --->'\n" );
  fprintf ( command_unit, "set title 'Fibonacci Disk Grid'\n" );
  fprintf ( command_unit, "set grid\n" );
  fprintf ( command_unit, "set key off\n" );
  fprintf ( command_unit, "set size ratio -1\n" );
  fprintf ( command_unit, "set style data lines\n" );
  fprintf ( command_unit, "plot '%s' using 1:2 with points lt 3 pt 3,\\\n", 
    data_filename );
  fprintf ( command_unit, "    '%s' using 1:2 lw 3 linecolor rgb 'black'\n", 
    boundary_filename );
  fprintf ( command_unit, "quit\n" );
  fclose ( command_unit );

  printf ( "  Created command file \"%s\"\n", command_filename );

  free ( g );

  return;
}
