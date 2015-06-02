# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>

# include "hypercube_grid.h"

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );

/******************************************************************************/

int main ( ) 

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for HYPERCUBE_GRID_PRB.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    29 August 2014

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "HYPERCUBE_GRID_PRB:\n" );
  printf ( "  C version\n" );
  printf ( "  Test the HYPERCUBE_GRID library.\n" );

  test01 ( );
  test02 ( );
  test03 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "HYPERCUBE_GRID_PRB:\n" );
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

    TEST01 tests HYPERCUBE_GRID on a two dimensional example.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    28 August 2014

  Author:

    John Burkardt
*/
{
# define M 2

  double a[M] = { 0.0, 0.0 };
  double b[M] = { 1.0, 10.0 };
  int c[M] = { 2, 4 };
  int i;
  int m = M;
  int n;
  int ns[M] = { 4, 5 };
  double *x;

  n = i4vec_product ( m, ns );

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  Create a grid using HYPERCUBE_GRID.\n" );
  printf ( "  Spatial dimension M = %d\n", m );
  printf ( "  Number of grid points N = %d\n", n );
  printf ( "\n" );
  printf ( "     I    NS     C      A         B\n" );
  printf ( "\n" );
  for ( i = 0; i < m; i++ )
  {
    printf ( "  %4d  %4d  %4d  %8.4f  %8.4f\n",
      i, ns[i], c[i], a[i], b[i] );
  }

  x = hypercube_grid ( m, n, ns, a, b, c );
  r8mat_transpose_print ( m, n, x, "  Grid points:" );
  free ( x );

  return;
# undef M
}
/******************************************************************************/

void test02 ( )

/******************************************************************************/
/*
  Purpose:

    TEST02 tests HYPERCUBE_GRID on a five dimensional example.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    28 August 2014

  Author:

    John Burkardt
*/
{
# define M 5

  double a[M] = { 0.0, 0.0, 0.0, 0.0, 0.0 };
  double b[M] = { 1.0, 1.0, 1.0, 1.0, 1.0 };
  int c[M] = { 1, 2, 3, 4, 5 };
  int i;
  int m = M;
  int n;
  int ns[M] = { 2, 2, 2, 2, 2 };
  double *x;

  n = i4vec_product ( m, ns );

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  Create a grid using HYPERCUBE_GRID.\n" );
  printf ( "  Use a two point grid in each dimension.\n" );
  printf ( "  Use a different centering option in each dimension.\n" );
  printf ( "  Spatial dimension M = %d\n", m );
  printf ( "  Number of grid points N = %d\n", n );
  printf ( "\n" );
  printf ( "     I    NS     C      A         B\n" );
  printf ( "\n" );
  for ( i = 0; i < m; i++ )
  {
    printf ( "  %4d  %4d  %4d  %8.4f  %8.4f\n",
      i, ns[i], c[i], a[i], b[i] );
  }

  x = hypercube_grid ( m, n, ns, a, b, c );
  r8mat_transpose_print ( m, n, x, "  Grid points:" );
  free ( x );

  return;
# undef M
}
/******************************************************************************/

void test03 ( )

/******************************************************************************/
/*
  Purpose:

    TEST03 tests HYPERCUBE_GRID on a three dimensional example.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    28 August 2014

  Author:

    John Burkardt
*/
{
# define M 3

  double a[M] = { -1.0, -1.0, -1.0 };
  double b[M] = { +1.0, +1.0, +1.0 };
  int c[M] = { 1, 1, 1 };
  int i;
  int m = M;
  int n;
  int ns[M] = { 3, 3, 3 };
  double *x;

  n = i4vec_product ( m, ns );

  printf ( "\n" );
  printf ( "TEST03\n" );
  printf ( "  Create a grid using HYPERCUBE_GRID.\n" );
  printf ( "  Use the same parameters in every dimension.\n" );
  printf ( "  Spatial dimension M = %d\n", m );
  printf ( "  Number of grid points N = %d\n", n );
  printf ( "\n" );
  printf ( "     I    NS     C      A         B\n" );
  printf ( "\n" );
  for ( i = 0; i < m; i++ )
  {
    printf ( "  %4d  %4d  %4d  %8.4f  %8.4f\n",
      i, ns[i], c[i], a[i], b[i] );
  }

  x = hypercube_grid ( m, n, ns, a, b, c );
  r8mat_transpose_print ( m, n, x, "  Grid points:" );
  free ( x );

  return;
# undef M
}
