# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <string.h>

# include "ellipse_grid.h"

int main ( void );
void ellipse_grid_test01 ( void );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for ELLIPSE_GRID_PRB.

  Discussion:

    ELLIPSE_GRID_PRB tests the ELLIPSE_GRID library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 November 2011

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "ELLIPSE_GRID_PRB:\n" );
  printf ( "  C version\n" );
  printf ( "  Test the ELLIPSE_GRID library.\n" );

  ellipse_grid_test01 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "ELLIPSE_GRID_PRB:\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void ellipse_grid_test01 ( void )

/******************************************************************************/
/*
  Purpose:

    ELLIPSE_GRID_TEST01 tests ELLIPSE_GRID.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 November 2011

  Author:

    John Burkardt
*/
{
  double c[2];
  char *filename = "ellipse_grid_test01.xy";
  int n;
  int ng;
  double r[2];
  double *xy;

  printf ( "\n" );
  printf ( "TEST01:\n" );
  printf ( "  ELLIPSE_GRID can define a grid of points\n" );
  printf ( "  with N+1 points on the minor half axis,\n" );
  printf ( "  based on any ellipse.\n" );

  n = 8;
  r[0] = 2.0;
  r[1] = 1.0;
  c[0] = 1.0;
  c[1] = 2.0;

  printf ( "\n" );
  printf ( "  We use N = %d\n", n );
  printf ( "  Radius R = (%g,%g)\n", r[0], r[1] );
  printf ( "  Center C = (%g,%g)\n", c[0], c[1] );

  ng = ellipse_grid_count ( n, r, c );

  printf ( "\n" );
  printf ( "  Number of grid points will be %d\n", ng );

  xy = ellipse_grid ( n, r, c, ng );

  r82vec_print_part ( ng, xy, 20, "  Part of the grid point array:" );

  r8mat_write ( filename, 2, ng, xy );

  printf ( "\n" );
  printf ( "  Data written to the file \"%s\".\n", filename );

  free ( xy );

  return;
}
