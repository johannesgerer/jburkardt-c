# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "triangle_grid.h"

int main ( void );

void test01 ( void );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for TRIANGLE_GRID_PRB.

  Discussion:

    TRIANGLE_GRID_PRB tests the TRIANGLE_GRID library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 September 2010

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "TRIANGLE_GRID_PRB:\n" );
  printf ( "  C version\n" );
  printf ( "  Test the TRIANGLE_GRID library.\n" );

  test01 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "TRIANGLE_GRID_PRB:\n" );
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

    TEST01 tests TRIANGLE_GRID.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 September 2010

  Author:

    John Burkardt
*/
{
  int n = 10;
  int ng = ((n+1)*(n+2))/2;

  char *filename = "triangle_grid_test01.xy";
  int j;
  FILE *output;
  double t[2*3] = {
    0.0, 0.0, 
    1.0, 0.0, 
    0.5, 0.86602540378443860 };
  double *tg;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  TRIANGLE_GRID can define a triangular grid of points\n" );
  printf ( "  with N+1 points on a side, based on any triangle.\n" );

  printf ( "\n" );
  printf ( "  Defining triangle:\n" );
  printf ( "     J      X      Y\n" );
  printf ( "\n" );
  for ( j = 0; j < 3; j++ )
  {
    printf ( "  %4d  %12f  %12f\n", j, t[0+j*2], t[1+j*2] );
  }
  tg = triangle_grid ( n, t );

  printf ( "\n" );
  printf ( "     J      X      Y\n" );
  printf ( "\n" );
  for ( j = 0; j < ng; j++ )
  {
    printf ( "  %4d  %12f  %12f\n", j, tg[0+j*2], tg[1+j*2] );
  }

  output = fopen ( filename, "wt" );
  for ( j = 0; j < ng; j++ )
  {
    fprintf ( output, "  %12f  %12f\n", tg[0+j*2], tg[1+j*2] );
  }
  fclose ( output );

  printf ( "\n" );
  printf ( "  Data written to \"%s\"\n"; );

  free ( tg );

  return;
}


