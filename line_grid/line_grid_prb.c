# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>

# include "line_grid.h"

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );

/******************************************************************************/

int main ( ) 

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for LINE_GRID_PRB.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    02 September 2014

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "LINE_GRID_PRB:\n" );
  printf ( "  C version\n" );
  printf ( "  Test the LINE_GRID library.\n" );

  test01 ( );
  test02 ( );
  test03 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "LINE_GRID_PRB:\n" );
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

    TEST01 tests LINE_GRID using simple parameters.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    31 August 2014

  Author:

    John Burkardt
*/
{
  double a = -1.0;
  double b = +1.0;
  int c = 1;
  int n = 11;
  double *x;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  Create a grid using LINE_GRID.\n" );
  printf ( "  Use simple parameters.\n" );
  printf ( "  Number of grid points N = %d\n", n );
  printf ( "\n" );
  printf ( "     N      C      A         B\n" );
  printf ( "\n" );
  printf ( "  %4d  %4d  %8.4f  %8.4f\n", n, c, a, b );

  x = line_grid ( n, a, b, c );
  r8vec_print ( n, x, "  Grid points:" );
  free ( x );

  return;
}
/******************************************************************************/

void test02 ( )

/******************************************************************************/
/*
  Purpose:

    TEST02 changes the number of points.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    31 August 2014

  Author:

    John Burkardt
*/
{
  double a = 0.0;
  double b = 1.0;
  int c = 2;
  int n;
  int test;
  double *x;

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  Create a grid using LINE_GRID.\n" );
  printf ( "  Try an increasing number of points.\n" );

  n = 4;

  for ( test = 1; test <= 3; test++ )
  {
    n = 2 * n + 1;

    printf ( "\n" );
    printf ( "     N      C      A         B\n" );
    printf ( "\n" );
    printf ( "  %4d  %4d  %8.4f  %8.4f\n", n, c, a, b );

    x = line_grid ( n, a, b, c );
    r8vec_print ( n, x, "  Grid points:" );
    free ( x );
  }

  return;
}
/******************************************************************************/

void test03 ( )

/******************************************************************************/
/*
  Purpose:

    TEST03 tries all the centering options.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    31 August 2014

  Author:

    John Burkardt
*/
{
  double a = 0.0;
  double b = 100.0;
  int c;
  int n;
  double *x;

  printf ( "\n" );
  printf ( "TEST03\n" );
  printf ( "  Create a grid using LINE_GRID.\n" );
  printf ( "  Try the different centering options.\n" );

  for ( c = 1; c <= 5; c++ )
  {
    printf ( "\n" );
    printf ( "     N      C      A         B\n" );
    printf ( "\n" );
    printf ( "  %4d  %4d  %8.4f  %8.4f\n", n, c, a, b );

    x = line_grid ( n, a, b, c );
    r8vec_print ( n, x, "  Grid points:" );
    free ( x );
  }

  return;
}
