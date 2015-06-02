# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <string.h>

# include "wedge_grid.h"

int main ( );
void test01 ( );
void test02 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for WEDGE_GRID_PRB.

  Discussion:

    WEDGE_GRID_PRB tests the WEDGE_GRID library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    22 August 2014

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "WEDGE_GRID_PRB:\n" );
  printf ( "  C version\n" );
  printf ( "  Test the WEDGE_GRID library.\n" );

  test01 ( );
  test02 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "WEDGE_GRID_PRB:\n" );
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

    TEST01 tests WEDGE_GRID.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    22 August 2014

  Author:

    John Burkardt
*/
{
  double *g;
  int j;
  int n = 5;
  int ng;
  char output_filename[255];
  FILE *output_unit;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  WEDGE_GRID can define a grid of points\n" );
  printf ( "  with N+1 points on a side\n" );
  printf ( "  over the interior of the unit wedge in 3D.\n" );

  printf ( "\n" );
  printf ( "  Grid order N = %d\n", n );

  ng = wedge_grid_size ( n );

  printf ( "  Grid count NG = %d\n", ng );

  g = wedge_grid ( n, ng );

  printf ( "\n" );
  printf ( "     J      X                Y               Z\n" );
  printf ( "\n" );
  for ( j = 0; j < ng; j++ )
  {
    printf ( "  %4d  %14.6g  %14.6g  %14.6g\n", j, g[0+j*3], g[1+j*3], g[2+j*3] );
  }

  strcpy ( output_filename, "wedge_grid.xy" );

  output_unit = fopen ( output_filename, "wt" );
  for ( j = 0; j < ng; j++ )
  {
    fprintf ( output_unit, "%g  %g  %g\n", g[0+j*3], g[1+j*3], g[2+j*3] );
  }
  fclose ( output_unit );

  printf ( "\n" );
  printf ( "  Data written to '%s'\n", output_filename );

  free ( g );

  return;
}
/******************************************************************************/

void test02 ( )

/******************************************************************************/
/*
  Purpose:

    TEST02 tests WEDGE_GRID_PLOT.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    22 August 2014

  Author:

    John Burkardt
*/
{
  double *g;
  char header[255];
  int n = 5;
  int ng;

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  WEDGE_GRID_PLOT can create GNUPLOT data files\n" );
  printf ( "  for displaying a wedge grid.\n" );

  printf ( "\n" );
  printf ( "  Grid order N = %d\n", n );

  ng = wedge_grid_size ( n );

  printf ( "  Grid count NG = %d\n", ng );

  g = wedge_grid ( n, ng );

  strcpy ( header, "wedge" );

  wedge_grid_plot ( n, ng, g, header );

  free ( g );

  return;
}

