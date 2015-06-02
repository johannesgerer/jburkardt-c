# include <stdlib.h>
# include <stdio.h>

# include "latin_random.h"

int main ( );
void test01 ( int *seed );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for LATIN_RANDOM_PRB.

  Discussion:

    LATIN_RANDOM_PRB tests the LATIN_RANDOM library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 November 2014

  Author:

    John Burkardt
*/
{
  int seed;
  int test;

  timestamp ( );
  printf ( "\n" );
  printf ( "LATIN_RANDOM_PRB:\n" );
  printf ( "  C version\n" );
  printf ( "  Test the LATIN_RANDOM library.\n" );

  seed = 123456789;

  for ( test = 0; test < 3; test++ )
  {
    test01 ( &seed );
  }
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "LATIN_RANDOM_PRB:\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void test01 ( int *seed )

/******************************************************************************/
/*
  Purpose:

    TEST01 tests LATIN_RANDOM_NEW.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 November 2014

  Author:

    John Burkardt

  Parameters:

    Input/output, int *SEED, a seed for the random number generator.
*/
{
  int m = 2;
  int i;
  int j;
  int k;
  int kk;
  int n = 10;
  double *x;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  LATIN_RANDOM chooses a Latin Square cell arrangement,\n" );
  printf ( "  and then chooses a random point from each cell.\n" );
  printf ( "\n" );
  printf ( "  Spatial dimension = %d\n", m );
  printf ( "  Number of points =  %d\n", n );
  printf ( "  Initial seed for UNIFORM = %d\n", seed );

  x = latin_random_new ( m, n, seed );

  r8mat_transpose_print ( m, n, x, "  Latin Random Square:" );

  free ( x );

  return;
}
