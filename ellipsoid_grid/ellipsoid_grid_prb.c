# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <string.h>

# include "ellipsoid_grid.h"

int main ( void );
void ellipsoid_grid_test01 ( void );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    ELLIPSOID_GRID_TEST tests ELLIPSOID_GRID.

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
  printf ( "ELLIPSOID_GRID_TEST:\n" );
  printf ( "  C version\n" );
  printf ( "  Test the ELLIPSOID_GRID library.\n" );

  ellipsoid_grid_test01 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "ELLIPSOID_GRID_TEST:\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void ellipsoid_grid_test01 ( void )

/******************************************************************************/
/*
  Purpose:

    ELLIPSOID_GRID_TEST01 tests ELLIPSOID_GRID.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 November 2011

  Author:

    John Burkardt
*/
{
  double c[3];
  char *filename = "ellipsoid_grid_test01.xyz";
  int n;
  int ng;
  double r[3];
  double *xyz;

  printf ( "\n" );
  printf ( "TEST01:\n" );
  printf ( "  ELLIPSOID_GRID can define a grid of points\n" );
  printf ( "  with N+1 points on the minor half axis,\n" );
  printf ( "  based on any ellipsoid.\n" );

  n = 4;
  r[0] = 2.0;
  r[1] = 1.0;
  r[2] = 1.5;
  c[0] = 1.0;
  c[1] = 2.0;
  c[2] = 1.5;

  printf ( "\n" );
  printf ( "  We use N = %d\n", n );
  printf ( "  Radius R = (%g,%g,%g)\n", r[0], r[1], r[2] );
  printf ( "  Center C = (%g,%g,%g)\n", c[0], c[1], c[2] );

  ng = ellipsoid_grid_count ( n, r, c );

  printf ( "\n" );
  printf ( "  Number of grid points will be %d\n", ng );

  xyz = ellipsoid_grid ( n, r, c, ng );

  r83vec_print_part ( ng, xyz, 20, "  Part of the grid point array:" );

  r8mat_write ( filename, 3, ng, xyz );

  printf ( "\n" );
  printf ( "  Data written to the file \"%s\".\n", filename );

  free ( xyz );

  return;
}
