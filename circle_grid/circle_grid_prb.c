# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>

# include "circle_grid.h"

int main ( void );
void circle_grid_test01 ( void );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    CIRCLE_GRID_TEST tests CIRCLE_GRID.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 November 2011

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "CIRCLE_GRID_TEST:\n" );
  printf ( "  C version\n" );
  printf ( "  Test the CIRCLE_GRID library.\n" );

  circle_grid_test01 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "CIRCLE_GRID_TEST:\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void circle_grid_test01 ( )

/******************************************************************************/
/*
  Purpose:

    CIRCLE_GRID_TEST01 tests CIRCLE_GRID.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 November 2011

  Author:

    John Burkardt
*/
{
  double c[2];
  double *cg;
  char *filename = "circle_grid_test01.xy";
  int n;
  int ng;
  double r;

  printf ( "\n" );
  printf ( "TEST01:\n" );
  printf ( "  CIRCLE_GRID can define a grid of points\n" );
  printf ( "  with N+1 points on a horizontal or vertical radius,\n" );
  printf ( "  based on any circle.\n" );

  n = 20;
  r = 2.0;
  c[0] = 1.0;
  c[1] = 5.0;

  printf ( "\n" );
  printf ( "  We use N = %d\n", n );
  printf ( "  Radius R = %g\n", r );
  printf ( "  Center C = (%g,%g)\n", c[0], c[1] );

  ng = circle_grid_count ( n, r, c );

  printf ( "\n" );
  printf ( "  Number of grid points will be %d\n", ng );

  cg = circle_grid ( n, r, c, ng );

  r82vec_print_part ( ng, cg, 20, "  Part of the grid point array:" );

  r8mat_write ( filename, 2, ng, cg );

  printf ( "\n" );
  printf ( "  Data written to the file \"%s\".\n", filename );

  free ( cg );

  return;
}
