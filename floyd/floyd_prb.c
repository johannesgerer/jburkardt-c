# include <stdlib.h>
# include <stdio.h>

# include "floyd.h"

int main ( main );

void test01 ( void );
void test02 ( void );
double test03 ( int n );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for FLOYD_PRB.

  Discussion:

    FLOYD_PRB calls a set of problems for FLOYD.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 July 2011

  Author:

    John Burkardt
*/
{
  int n;
  double ratio;
  double wtime;

  timestamp ( );

  printf ( "\n" );
  printf ( "FLOYD_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the FLOYD library.\n" );

  test01 ( );
  test02 ( );

  printf ( "\n" );
  printf ( "FLOYD_TEST03\n" );
  printf ( "  Test I4MAT_FLOYD on the MOD(I,J) matrix.\n" );
  printf ( "  The work is roughly N^3.\n" );
  printf ( "\n" );
  printf ( "         N   Time (seconds)  Time/N^3\n" );
  printf ( "\n" );

  n = 1;
  while ( n <= 2048 )
  {
    wtime = test03 ( n );
    ratio = 1000000.0 * wtime / ( double ) ( n ) / ( double ) ( n ) / ( double ) ( n );
    printf ( "  %8d  %14f  %14f\n", n, wtime, ratio );
    n = n * 2;
  }
//
//  Terminate.
//
  printf ( "\n" );
  printf ( "FLOYD_PRB\n" );
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

    TEST01 tests I4MAT_FLOYD.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 July 2011

  Author:

    John Burkardt
*/
{
# define N 6

  int a[N*N] = {
     0, -1, -1, -1, -1, -1,
     2,  0, -1, -1, -1,  5,
     5,  7,  0, -1,  2, -1,
    -1,  1,  4,  0, -1,  2,
    -1, -1, -1,  3,  0,  4,
    -1,  8, -1, -1,  3,  0  };
  int huge;
  int i;
  int j;
  int n = N;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  I4MAT_FLOYO uses Floyd's algorithm to find the\n" );
  printf ( "  shortest distance between all pairs of nodes\n" );
  printf ( "  in a directed graph, starting from the initial array\n" );
  printf ( "  of direct node-to-node distances.\n" );

  printf ( "\n" );
  printf ( "  In the initial direct distance array, if\n" );
  printf ( "    A(I,J) = -1,\n" );
  printf ( "  this indicates there is NO directed link from\n" );
  printf ( "  node I to node J.  In that case, the value of\n" );
  printf ( "  of A(I,J) is essentially \"infinity\".\n" );

  i4mat_print ( n, n, a, "  Initial direct distance array:" );

  huge = i4_huge ( ) / 2;

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      if ( a[i+j*n] == - 1 )
      {
        a[i+j*n] = huge;
      }
    }
  }

  i4mat_floyd ( n, a );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      if ( a[i+j*n] == huge )
      {
        a[i+j*n] = - 1;
      }
    }
  }

  printf ( "\n" );
  printf ( "  In the final shortest distance array, if\n" );
  printf ( "    A(I,J) = -1,\n" );
  printf ( "  this indicates there is NO directed path from\n" );
  printf ( "  node I to node J.\n" );

  i4mat_print ( n, n, a, "  Final shortest distance array:" );

  return;
# undef N
}
/******************************************************************************/

void test02 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST02 tests R8MAT_FLOYD.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 July 2011

  Author:

    John Burkardt
*/
{
# define N 6

  double a[N*N] = {
     0.0, -1.0, -1.0, -1.0, -1.0, -1.0,
     2.0,  0.0, -1.0, -1.0, -1.0,  5.0,
     5.0,  7.0,  0.0, -1.0,  2.0, -1.0,
    -1.0,  1.0,  4.0,  0.0, -1.0,  2.0,
    -1.0, -1.0, -1.0,  3.0,  0.0,  4.0,
    -1.0,  8.0, -1.0, -1.0,  3.0,  0.0  };
  double huge;
  int i;
  int j;
  int n = N;

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  R8MAT_FLOYO uses Floyd's algorithm to find the\n" );
  printf ( "  shortest distance between all pairs of nodes\n" );
  printf ( "  in a directed graph, starting from the initial array\n" );
  printf ( "  of direct node-to-node distances.\n" );

  printf ( "\n" );
  printf ( "  In the initial direct distance array, if\n" );
  printf ( "    A(I,J) = -1,\n" );
  printf ( "  this indicates there is NO directed link from\n" );
  printf ( "  node I to node J.  In that case, the value of\n" );
  printf ( "  of A(I,J) is essentially \"infinity\".\n" );

  r8mat_print ( n, n, a, "  Initial direct distance array:" );

  huge = r8_huge ( );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      if ( a[i+j*n] == - 1.0 )
      {
        a[i+j*n] = huge;
      }
    }
  }

  r8mat_floyd ( n, a );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      if ( a[i+j*n] == huge )
      {
        a[i+j*n] = - 1.0;
      }
    }
  }

  printf ( "\n" );
  printf ( "  In the final shortest distance array, if\n" );
  printf ( "    A(I,J) = -1,\n" );
  printf ( "  this indicates there is NO directed path from\n" );
  printf ( "  node I to node J.\n" );

  r8mat_print ( n, n, a, "  Final shortest distance array:" );

  return;
# undef N
}
/******************************************************************************/

double test03 ( int n )

/******************************************************************************/
/*
  Purpose:

    TEST03 tests I4MAT_FLOYD.

  Discussion:

    The matrix size is input by the user.

    The matrix A has the property that

      A(I,J) = 1 if I is divisible by J.

  Example:

    N = 6

    1 0 0 0 0 0
    1 1 0 0 0 0
    1 0 1 0 0 0
    1 1 0 1 0 0
    1 0 0 0 1 0
    1 1 1 0 0 1

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 July 2011

  Author:

    John Burkardt

  Parameters:

    Input, int N, the size of the matrix.

    Output, double TEST03, the CPU time required by I4MAT_FLOYD.
*/
{
  int *a;
  int huge;
  int i;
  int j;
  double time1;
  double time2;
  double wtime;

  a = ( int * ) malloc ( n * n * sizeof ( int ) );

  huge = i4_huge ( ) / 2;

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      if ( ( i + 1 ) % ( j + 1 ) == 0 )
      {
        a[i+j*n] = 1;
      }
      else
      {
        a[i+j*n] = huge;
      }
    }
  }

  time1 = cpu_time ( );

  i4mat_floyd ( n, a );

  time2 = cpu_time ( );

  wtime = time2 - time1;

  free ( a );

  return wtime;
}
