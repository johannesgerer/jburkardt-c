# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>

# include "square_grid.h"

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );

/******************************************************************************/

int main ( ) 

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for SQUARE_GRID_PRB.

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
  printf ( "SQUARE_GRID_PRB:\n" );
  printf ( "  C version\n" );
  printf ( "  Test the SQUARE_GRID library.\n" );

  test01 ( );
  test02 ( );
  test03 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "SQUARE_GRID_PRB:\n" );
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

    TEST01 tests SQUARE_GRID using the same parameters for all dimensions.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    31 August 2014

  Author:

    John Burkardt
*/
{
  double a[2] = { -1.0, -1.0 };
  double b[2] = { +1.0, +1.0 };
  int c[2] = { 1, 1 };
  int i;
  int n;
  int ns[2] = { 3, 3 };
  double *x;

  n = ns[0] * ns[1];

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  Create a grid using SQUARE_GRID.\n" );
  printf ( "  Use the same parameters in every dimension.\n" );
  printf ( "  Number of grid points N = %d\n", n );
  printf ( "\n" );
  printf ( "     I    NS     C      A         B\n" );
  printf ( "\n" );
  for ( i = 0; i < 2; i++ )
  {
    printf ( "  %4d  %4d  %4d  %8.4f  %8.4f\n",
      i, ns[i], c[i], a[i], b[i] );
  }

  x = square_grid ( n, ns, a, b, c );
  r8mat_transpose_print ( 2, n, x, "  Grid points:" );
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
  double a[2] = { 0.0, 0.0 };
  double b[2] = { 1.0, 1.0 };
  int c[2] = { 2, 2 };
  int i;
  int n;
  int ns[2] = { 4, 2 };
  double *x;

  n = ns[0] * ns[1];

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  Create a grid using SQUARE_GRID.\n" );
  printf ( "  se a different number of points in each dimension..\n" );
  printf ( "  Number of grid points N = %d\n", n );
  printf ( "\n" );
  printf ( "     I    NS     C      A         B\n" );
  printf ( "\n" );
  for ( i = 0; i < 2; i++ )
  {
    printf ( "  %4d  %4d  %4d  %8.4f  %8.4f\n",
      i, ns[i], c[i], a[i], b[i] );
  }

  x = square_grid ( n, ns, a, b, c );
  r8mat_transpose_print ( 2, n, x, "  Grid points:" );
  free ( x );

  return;
}
/******************************************************************************/

void test03 ( )

/******************************************************************************/
/*
  Purpose:

    TEST03 uses a square with different sizes in each dimension.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    31 August 2014

  Author:

    John Burkardt
*/
{
  double a[2] = {   0.0, -2.0 };
  double b[2] = { +10.0, +2.0 };
  int c[2] = { 3, 4 };
  int i;
  int n;
  int ns[2] = { 3, 3 };
  double *x;

  n = ns[0] * ns[1];

  printf ( "\n" );
  printf ( "TEST03\n" );
  printf ( "  Create a grid using SQUARE_GRID.\n" );
  printf ( "  Use a different physical size in each dimension.\n" );
  printf ( "  Number of grid points N = %d\n", n );
  printf ( "\n" );
  printf ( "     I    NS     C      A         B\n" );
  printf ( "\n" );
  for ( i = 0; i < 2; i++ )
  {
    printf ( "  %4d  %4d  %4d  %8.4f  %8.4f\n",
      i, ns[i], c[i], a[i], b[i] );
  }

  x = square_grid ( n, ns, a, b, c );
  r8mat_transpose_print ( 2, n, x, "  Grid points:" );
  free ( x );

  return;
}
