# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <string.h>

# include "pyramid_monte_carlo.h"

int main ( );
void test01 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for PYRAMID_MONTE_CARLO_PRB.

  Discussion:

    PYRAMID_MONTE_CARLO_PRB tests the PYRAMID_MONTE_CARLO library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 April 2014

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "PYRAMID_MONTE_CARLO_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the PYRAMID_MONTE_CARLO library.\n" );

  test01 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "PYRAMID_MONTE_CARLO_PRB\n" );
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

    TEST01 compares exact and estimated monomial integrals.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 April 2014

  Author:

    John Burkardt
*/
{
# define M 3
# define TEST_NUM 10

  int e[M];
  int e_test[M*TEST_NUM] = {
    0, 0, 0, 
    0, 0, 1, 
    2, 0, 0, 
    0, 2, 0, 
    0, 0, 2, 
    2, 0, 1, 
    0, 2, 1, 
    0, 0, 3, 
    2, 2, 0, 
    2, 0, 2 };
  double error;
  double exact;
  int i;
  int j;
  int m = M;
  int n;
  double result;
  int seed;
  int test_num = TEST_NUM;
  double *value;
  double *x;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  Use PYRAMID01_SAMPLE to estimate integrals\n" );
  printf ( "  over the interior of the unit pyramid in 3D.\n" );

  seed = 123456789;

  printf ( "\n" );
  printf ( "         N" );
  printf ( "        1" );
  printf ( "               Z" );
  printf ( "             X^2" );
  printf ( "             Y^2" );
  printf ( "             Z^2" );
  printf ( "            X^2Z" );
  printf ( "            Y^2Z" );
  printf ( "             Z^3" );
  printf ( "          X^2Y^2" );
  printf ( "          X^2Z^2\n" );
  printf ( "\n" );

  n = 1;

  while ( n <= 65536 )
  {
    printf ( "  %8d", n );

    x = pyramid01_sample ( n, &seed );

    for ( j = 0; j < test_num; j++ )
    {
      for ( i = 0; i < m; i++ )
      {
        e[i] = e_test[i+j*m];
      }
      value = monomial_value ( m, n, e, x );

      result = pyramid01_volume ( ) * r8vec_sum ( n, value ) / ( double ) ( n );
      printf ( "  %14.6g", result );
    }
    printf ( "\n" );

    free ( value );
    free ( x );

    n = 2 * n;
  }

  printf ( "\n" );
  printf ( "     Exact" );

  for ( j = 0; j < 10; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      e[i] = e_test[i+j*m];
    }
    result = pyramid01_integral ( e );
    printf ( "  %14.6g", result );
  }
  printf ( "\n" );

  return;
}
