# include <math.h>
# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <time.h>

# include "sphere_fibonacci_grid.h"

int main ( );
void sphere_fibonacci_grid_points_test ( );
void sphere_fibonacci_grid_display_test ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    SPHERE_FIBONACCI_GRID_TEST tests the SPHERE_FIBONACCI_GRID library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 May 2015

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "SPHERE_FIBONACCI_GRID_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the SPHERE_FIBONACCI_GRID library.\n" );

  sphere_fibonacci_grid_points_test ( );
  sphere_fibonacci_grid_display_test ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "SPHERE_FIBONACCI_GRID_PRB\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void sphere_fibonacci_grid_points_test ( )

/******************************************************************************/
/*
  Purpose:

    SPHERE_FIBONACCI_GRID_POINTS_TEST tests SPHERE_FIBONACCI_GRID_POINTS.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 May 2015

  Author:

    John Burkardt
*/
{
  char filename[255];
  int ng;
  double *xg;

  printf ( "\n" );
  printf ( "SPHERE_FIBONACCI_GRID_POINTS_TEST\n" );
  printf ( "  SPHERE_FIBONACCI_GRID_POINTS computes points on a sphere\n" );
  printf ( "  that lie on a Fibonacci spiral.\n" );

  ng = 1000;
  printf ( "\n" );
  printf ( "  Number of points NG = %d\n", ng );

  xg = sphere_fibonacci_grid_points ( ng );

  r8mat_transpose_print_some ( 3, ng, xg, 1, 1, 3, 10, 
    "  Part of the grid array:" );
/*
  Write the nodes to a file.
*/
  strcpy ( filename, "sphere_fibonacci_grid_n1000.xyz" );

  r8mat_write ( filename, 3, ng, xg );

  free ( xg );

  return;
}
/******************************************************************************/

void sphere_fibonacci_grid_display_test ( )

/******************************************************************************/
/*
  Purpose:

    SPHERE_FIBONACCI_GRID_DISPLAY_TEST tests SPHERE_FIBONACCI_GRID_DISPLAY.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 May 2015

  Author:

    John Burkardt
*/
{
  int ng;
  char prefix[255];
  double *xg;

  printf ( "\n" );
  printf ( "SPHERE_FIBONACCI_GRID_DISPLAY_TEST\n" );
  printf ( "  SPHERE_FIBONACCI_GRID_DISPLAY displays points\n" );
  printf ( "  on a sphere that lie on a Fibonacci spiral.\n" );

  ng = 1000;
  printf ( "\n" );
  printf ( "  Number of points NG = %d\n", ng );

  xg = sphere_fibonacci_grid_points ( ng );
/*
  Display the nodes on a sphere.
*/
  strcpy ( prefix, "sphere_fibonacci_grid_n1000" );

  sphere_fibonacci_grid_display ( ng, xg, prefix );

  free ( xg );

  return;
}
