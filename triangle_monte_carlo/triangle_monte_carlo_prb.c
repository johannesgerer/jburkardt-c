# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "triangle_monte_carlo.h"

int main ( );
void test01 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for TRIANGLE_MONTE_CARLO_PRB.

  Discussion:

    TRIANGLE_MONTE_CARLO_PRB tests the TRIANGLE_MONTE_CARLO library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 January 2014

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "TRIANGLE_MONTE_CARLO_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the TRIANGLE_MONTE_CARLO library.\n" );

  test01 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "TRIANGLE_MONTE_CARLO_PRB\n" );
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

    TEST01 uses TRIANGLE_SAMPLE_01 with an increasing number of points.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 January 2014

  Author:

    John Burkardt
*/
{
  int e[2];
  int e_test[2*7] = {
    0, 0, 
    1, 0, 
    0, 1, 
    2, 0, 
    1, 1, 
    0, 2, 
    3, 0 };
  double error;
  double exact;
  int i;
  int j;
  int m = 2;
  int n;
  double result;
  int seed;
  double *value;
  double *x;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  Use TRIANGLE01_SAMPLE for a Monte Carlo estimate of an\n" );
  printf ( "  integral over the interior of the unit triangle in 2D.\n" );

  seed = 123456789;

  printf ( "\n" );
  printf ( "         N        1               X               Y " );
  printf ( "             X^2               XY             Y^2             X^3\n" );
  printf ( "\n" );

  n = 1;

  while ( n <= 65536 )
  {
    x = triangle01_sample ( n, &seed );
    printf ( "  %8d", n );

    for ( j = 0; j < 7; j++ )
    {
      for ( i = 0; i < m; i++ )
      {
        e[i] = e_test[i+j*m];
      }
      value = monomial_value ( m, n, e, x );

      result = triangle01_area ( ) * r8vec_sum ( n, value ) / ( double ) ( n );
      printf ( "  %14.6g", result );

      free ( value );
    }

    printf ( "\n" );

    free ( x );

    n = 2 * n;
  }

  printf ( "\n" );
  printf ( "     Exact" );
  for ( j = 0; j < 7; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      e[i] = e_test[i+j*m];
    }
    result = triangle01_monomial_integral ( e );
    printf ( "  %14.6g", result );
  }

  printf ( "\n" );

  return;
}
