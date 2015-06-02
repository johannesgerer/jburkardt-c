# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>

# include "cube_grid.h"

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );

/******************************************************************************/

int main ( ) 

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for CUBE_GRID_PRB.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    31 August 2014

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "CUBE_GRID_PRB:\n" );
  printf ( "  C version\n" );
  printf ( "  Test the CUBE_GRID library.\n" );

  test01 ( );
  test02 ( );
  test03 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "CUBE_GRID_PRB:\n" );
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

    TEST01 tests CUBE_GRID using the same parameters for all dimensions.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    31 August 2014

  Author:

    John Burkardt
*/
{
  double a[3] = { -1.0, -1.0, -1.0 };
  double b[3] = { +1.0, +1.0, +1.0 };
  int c[3] = { 1, 1, 1 };
  int i;
  int n;
  int ns[3] = { 3, 3, 3 };
  double *x;

  n = ns[0] * ns[1] * ns[2];

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  Create a grid using CUBE_GRID.\n" );
  printf ( "  Use the same parameters in every dimension.\n" );
  printf ( "  Number of grid points N = %d\n", n );
  printf ( "\n" );
  printf ( "     I    NS     C      A         B\n" );
  printf ( "\n" );
  for ( i = 0; i < 3; i++ )
  {
    printf ( "  %4d  %4d  %4d  %8.4f  %8.4f\n",
      i, ns[i], c[i], a[i], b[i] );
  }

  x = cube_grid ( n, ns, a, b, c );
  r8mat_transpose_print ( 3, n, x, "  Grid points:" );
  free ( x );

  return;
}
/******************************************************************************/

void test02 ( )

/******************************************************************************/
/*
  Purpose:

    TEST02 uses a different number of points in each coordinate.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    31 August 2014

  Author:

    John Burkardt
*/
{
  double a[3] = { 0.0, 0.0, 0.0 };
  double b[3] = { 1.0, 1.0, 1.0  };
  int c[3] = { 2, 2, 2 };
  int i;
  int n;
  int ns[3] = { 4, 2, 3 };
  double *x;

  n = ns[0] * ns[1] * ns[2];

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  Create a grid using CUBE_GRID.\n" );
  printf ( "  se a different number of points in each dimension..\n" );
  printf ( "  Number of grid points N = %d\n", n );
  printf ( "\n" );
  printf ( "     I    NS     C      A         B\n" );
  printf ( "\n" );
  for ( i = 0; i < 3; i++ )
  {
    printf ( "  %4d  %4d  %4d  %8.4f  %8.4f\n",
      i, ns[i], c[i], a[i], b[i] );
  }

  x = cube_grid ( n, ns, a, b, c );
  r8mat_transpose_print ( 3, n, x, "  Grid points:" );
  free ( x );

  return;
}
/******************************************************************************/

void test03 ( )

/******************************************************************************/
/*
  Purpose:

    TEST03 uses a cube with different sizes in each dimension.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    31 August 2014

  Author:

    John Burkardt
*/
{
  double a[3] = {   0.0, -2.0, 50.0 };
  double b[3] = { +10.0, +2.0, 51.0 };
  int c[3] = { 3, 4, 5 };
  int i;
  int n;
  int ns[3] = { 3, 3, 3 };
  double *x;

  n = ns[0] * ns[1] * ns[2];

  printf ( "\n" );
  printf ( "TEST03\n" );
  printf ( "  Create a grid using CUBE_GRID.\n" );
  printf ( "  Use a different physical size in each dimension.\n" );
  printf ( "  Number of grid points N = %d\n", n );
  printf ( "\n" );
  printf ( "     I    NS     C      A         B\n" );
  printf ( "\n" );
  for ( i = 0; i < 3; i++ )
  {
    printf ( "  %4d  %4d  %4d  %8.4f  %8.4f\n",
      i, ns[i], c[i], a[i], b[i] );
  }

  x = cube_grid ( n, ns, a, b, c );
  r8mat_transpose_print ( 3, n, x, "  Grid points:" );
  free ( x );

  return;
}
