# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <string.h>

# include "circle_arc_grid.h"

int main ( void );
void test01 ( void );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for CIRCLE_ARC_GRID_PRB.

  Discussion:

    CIRCLE_ARC_GRID_PRB tests the CIRCLE_ARC_GRID library.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    13 November 2011

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "CIRCLE_ARC_GRID_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the CIRCLE_ARC_GRID library.\n" );

  test01 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "CIRCLE_ARC_GRID_PRB\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void test01 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST01 demonstrates the use of CIRCLE_ARC_GRID.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    13 November 2011

  Author:

    John Burkardt
*/
{
  double a[2];
  double c[2];
  char *filename = "arc.txt";
  int n;
  double r;
  double *xy;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  Compute points along a 90 degree arc\n" );

  r = 2.0;
  c[0] = 5.0;
  c[1] = 5.0;
  a[0] = 0.0;
  a[1] = 90.0;
  n = 10;
/*
  Echo the input.
*/
  printf ( "\n" );
  printf ( "  Radius =           %g\n", r );
  printf ( "  Center =           %g  %g\n", c[0], c[1] );
  printf ( "  Angle 1 =          %g\n", a[0] );
  printf ( "  Angle 2 =          %g\n", a[1] );
  printf ( "  Number of points = %d\n", n );
/*
  Compute the data.
*/
  xy = circle_arc_grid ( r, c, a, n );
/*
  Print a little of the data.
*/
  r82vec_print_part ( n, xy, 5, "  A few of the points:" );
/*
  Write the data.
*/
  r8mat_write ( filename, 2, n, xy );
  printf ( "\n" );
  printf ( "  Data written to \"%s\".\n", filename );
/*
  Free memory.
*/
  free ( xy );

  return;
}
