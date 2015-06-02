# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <string.h>

# include "simplex_grid.h"

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );
void test04 ( );
void test05 ( );
void test06 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    SIMPLEX_GRID_TEST tests the SIMPLEX_GRID library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 July 2014

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "SIMPLEX_GRID_TEST:\n" );
  printf ( "  C version\n" );
  printf ( "  Test the SIMPLEX_GRID library.\n" );

  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
  test05 ( );
  test06 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "SIMPLEX_GRID_TEST:\n" );
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

    TEST01 tests SIMPLEX_GRID_SIZE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 July 2014

  Author:

    John Burkardt
*/
{
  int m;
  int n;
  int ng;

  printf ( "\n" );
  printf ( "TEST01:\n" );
  printf ( "  SIMPLEX_GRID_SIZE counts the points in a regular grid\n" );
  printf ( "  with N+1 points on a side, in an M-dimensional simplex.\n" );
  printf ( "\n" );
  printf ( "        M: 0     1     2     3     4     5\n" );
  printf ( "    N:\n" );
  for ( n = 0; n <= 10; n++ )
  {
    printf ( "  %3d:", n );
    for ( m = 0; m <= 5; m++ )
    {
      ng = simplex_grid_size ( m, n );
      printf ( "%6d", ng );
    }
    printf ( "\n" );
  }
  return;
}
/******************************************************************************/

void test02 ( )

/******************************************************************************/
/*
  Purpose:

    TEST02 tests SIMPLEX_GRID_INDEX_NEXT.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 July 2014

  Author:

    John Burkardt
*/
{
  int *g;
  int i;
  int j;
  int m = 3;
  int n;

  printf ( "\n" );
  printf ( "TEST02:\n" );
  printf ( "  SIMPLEX_GRID_INDEX_NEXT lists, one by one, the indices\n" );
  printf ( "  of a simplex grid that uses N+1 points on a side,\n" );
  printf ( "  in an M-dimensional simplex.\n" );
  printf ( "\n" );
  printf ( "   #:  1  2  3  (*)\n" );
  printf ( "\n" );

  n = 3;

  j = 0;
  g = ( int * ) malloc ( ( m + 1 ) * sizeof ( int ) );
  for ( i = 0; i < m; i++ )
  {
    g[i] = 0;
  }
  g[m] = n;
  
  while ( 1 )
  {
    printf ( "  %2d:", j );
    for ( i = 0; i < m; i++ )
    {
      printf ( "%3d", g[i] );
    }
    printf ( " (%3d)\n", g[m] );

    if ( g[0] == n )
    {
      break;
    }

    simplex_grid_index_next ( m, n, g );

    j = j + 1;
  }

  free ( g );

  return;
}
/******************************************************************************/

void test03 ( )

/******************************************************************************/
/*
  Purpose:

    TEST03 tests SIMPLEX_GRID_INDEX_SAMPLE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 July 2014

  Author:

    John Burkardt
*/
{
  int *g;
  int i;
  int j;
  int m = 3;
  int n;
  int seed;

  printf ( "\n" );
  printf ( "TEST03:\n" );
  printf ( "  SIMPLEX_GRID_INDEX_SAMPLE returns a randomly selected\n" );
  printf ( "  index of a simplex grid that uses N+1 points on a side,\n" );
  printf ( "  in an M-dimensional simplex.\n" );
  printf ( "\n" );
  printf ( "   #:  1  2  3  (*)\n" );
  printf ( "\n" );

  n = 3;
  seed = 123456789;

  for ( j = 1; j <= 20; j++ )
  {
    g = simplex_grid_index_sample ( m, n, &seed );

    printf ( "  %2d:", j );
    for ( i = 0; i < m; i++ )
    {
      printf ( "%3d", g[i] );
    }
    printf ( " (%3d)\n", g[m] );

    free ( g );
  }

  return;
}
/******************************************************************************/

void test04 ( )

/******************************************************************************/
/*
  Purpose:

    TEST04 tests SIMPLEX_GRID_INDEX_TO_POINT.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 July 2014

  Author:

    John Burkardt
*/
{
  int *g;
  int i;
  int j;
  int m = 2;
  int n;
  int seed;
  double v[2*3] = {
    20.0,  0.0, 
    30.0, 40.0, 
    10.0, 20.0 };
  double *x;

  printf ( "\n" );
  printf ( "TEST04:\n" );
  printf ( "  SIMPLEX_GRID_INDEX_TO_POINT returns the physical point\n" );
  printf ( "  corresponding to a grid index of a simplex grid that\n" );
  printf ( "  that uses N+1 points on a side,\n" );
  printf ( "  in an M-dimensional simplex.\n" );

  n = 5;

  r8mat_transpose_print ( m, m + 1, v, "  Simplex vertices:" );

  printf ( "\n" );
  printf ( "  Choosing random simplex indices to convert:\n" );
  printf ( "   #:  1  2  3     X        Y\n" );
  printf ( "\n" );

  seed = 123456789;

  for ( j = 1; j <= 10; j++ )
  {
    g = simplex_grid_index_sample ( m, n, &seed );
    x = simplex_grid_index_to_point ( m, n, 1, g, v );

    printf ( "  %2d:", j );
    for ( i = 0; i <= m; i++ )
    {
      printf ( "%3d", g[i] );
    }
    printf ( "  " );
    for ( i = 0; i < m; i++ )
    {
      printf ( "  %8.4f", x[i] );
    }
    printf ( "\n" );

    free ( g );
    free ( x );
  }

  return;
}
/******************************************************************************/

void test05 ( )

/******************************************************************************/
/*
  Purpose:

    TEST05 tests SIMPLEX_GRID_INDEX_ALL.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 July 2014

  Author:

    John Burkardt
*/
{
  int *grid;
  int m;
  int n;
  int ng;

  printf ( "\n" );
  printf ( "TEST05:\n" );
  printf ( "  SIMPLEX_GRID_INDEX_ALL returns all the indices\n" );
  printf ( "  of a simplex grid that uses N+1 points on a side,\n" );
  printf ( "  in an M-dimensional simplex.\n" );

  m = 3;
  n = 3;
  ng = simplex_grid_size ( m, n );

  grid = simplex_grid_index_all ( m, n, ng );

  i4mat_transpose_print ( m + 1, ng, grid, 
    "  Transposed Simplex Grid Index Matrix:" );

  free ( grid );

  return;
}
/******************************************************************************/

void test06 ( )

/******************************************************************************/
/*
  Purpose:

    TEST06 tests SIMPLEX_GRID_INDEX_TO_POINT.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 July 2014

  Author:

    John Burkardt
*/
{
  int *grid;
  int m = 2;
  int n;
  int ng;
  double v[2*3] = {
    20.0,  0.0, 
    30.0, 40.0, 
    10.0, 20.0 };
  double *x;

  printf ( "\n" );
  printf ( "TEST06:\n" );
  printf ( "  SIMPLEX_GRID_INDEX_TO_POINT returns the physical point\n" );
  printf ( "  corresponding to a grid index of a simplex grid that\n" );
  printf ( "  that uses N+1 points on a side,\n" );
  printf ( "  in an M-dimensional simplex.\n" );

  n = 5;
  ng = simplex_grid_size ( m, n );

  r8mat_transpose_print ( m, m + 1, v, "  Simplex vertices:" );

  grid = simplex_grid_index_all ( m, n, ng );

  x = simplex_grid_index_to_point ( m, n, ng, grid, v );

  r8mat_transpose_print ( m, ng, x, "  Grid Point Coordinates:" );

  free ( grid );
  free ( x );

  return;
}

