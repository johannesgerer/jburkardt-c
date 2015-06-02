# include <stdlib.h>
# include <stdio.h>

# include "ball_grid.h"

int main ( void );
void ball_grid_test01 ( void );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for BALL_GRID_PRB.

  Discussion:

    BALL_GRID_PRB tests the BALL_GRID library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 November 2011

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "BALL_GRID_PRB:\n" );
  printf ( "  C version\n" );
  printf ( "  Test the BALL_GRID library.\n" );

  ball_grid_test01 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "BALL_GRID_PRB:\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void ball_grid_test01 ( void )

/******************************************************************************/
/*
  Purpose:

    BALL_GRID_TEST01 tests BALL_GRID.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    10 November 2011

  Author:

    John Burkardt
*/
{
  double *bg;
  double c[3];
  char *filename = "ball_grid_test01.xyz";
  int n;
  int ng;
  double r;

  printf ( "\n" );
  printf ( "TEST01:\n" );
  printf ( "  BALL_GRID can define a grid of points\n" );
  printf ( "  with N+1 points on a horizontal or vertical radius,\n" );
  printf ( "  based on any ball.\n" );

  n = 10;
  r = 2.0;
  c[0] = 1.0;
  c[1] = 5.0;
  c[2] = 2.0;

  printf ( "\n" );
  printf ( "  We use N = %d\n", n );
  printf ( "  Radius R = %g\n", r );
  printf ( "  Center C = (%g,%g,%g)\n", c[0], c[1], c[2] );

  ng = ball_grid_count ( n, r, c );

  printf ( "\n" );
  printf ( "  Number of grid points will be %d\n", ng );

  bg = ball_grid ( n, r, c, ng );

  r83vec_print_part ( ng, bg, 20, "  Part of the grid point array:" );

  r8mat_write ( filename, 3, ng, bg );

  printf ( "\n" );
  printf ( "  Data written to the file \"%s\".\n", filename );

  return;
}
