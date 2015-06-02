# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "cube_integrals.h"

int main ( );
void test01 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for CUBE_INTEGRALS_PRB.

  Discussion:

    CUBE_INTEGRALS_PRB tests the CUBE_INTEGRALS library.
    
  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 January 2014

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "CUBE_INTEGRALS_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the CUBE_INTEGRALS library.\n" );

  test01 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "CUBE_INTEGRALS_PRB\n" );
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

    TEST01 estimates integrals over the unit cube in 3D.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 January 2014

  Author:

    John Burkardt
*/
{
  int *e;
  double error;
  double exact;
  int i;
  int j;
  int m = 3;
  int n = 4192;
  double result;
  int seed;
  int test;
  int test_num = 20;
  double *value;
  double *x;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  Compare exact and estimated integrals\n" );
  printf ( "  over the interior of the unit cube in 3D.\n" );
/*
  Get sample points.
*/
  seed = 123456789;
  x = cube01_sample ( n, &seed );
  printf ( "\n" );
  printf ( "  Number of sample points is %d\n", n );
/*
  Randomly choose exponents.
*/
  printf ( "\n" );
  printf ( "  Ex  Ey  Ez     MC-Estimate           Exact      Error\n" );
  printf ( "\n" );

  for ( test = 1; test <= test_num; test++ )
  {
    e = i4vec_uniform_ab_new ( m, 0, 7, &seed );

    value = monomial_value ( m, n, e, x );

    result = cube01_volume ( ) * r8vec_sum ( n, value ) / ( double ) ( n );
    exact = cube01_monomial_integral ( e );
    error = fabs ( result - exact );

    printf ( "  %2d  %2d  %2d  %14.6g  %14.6g  %10.2e\n",
      e[0], e[1], e[2], result, exact, error );

    free ( value );
  }

  free ( x );

  return;
}
