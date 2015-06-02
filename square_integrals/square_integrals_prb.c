# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "square_integrals.h"

int main ( );
void test01 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for SQUARE_INTEGRALS_PRB.

  Discussion:

    SQUARE_INTEGRALS_PRB tests the SQUARE_INTEGRALS library.
    
  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    18 January 2014

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "SQUARE_INTEGRALS_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the SQUARE_INTEGRALS library.\n" );

  test01 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "SQUARE_INTEGRALS_PRB\n" );
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

    TEST01 compares exact and estimated integrals over the unit square in 2D.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    18 January 2014

  Author:

    John Burkardt
*/
{
  int *e;
  double error;
  double exact;
  int i;
  int j;
  int m = 2;
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
  printf ( "  over the interior of the unit square in 2D.\n" );
/*
  Get sample points.
*/
  seed = 123456789;
  x = square01_sample ( n, &seed );
  printf ( "\n" );
  printf ( "  Number of sample points is %d\n", n );
/*
  Randomly choose exponents.
*/
  printf ( "\n" );
  printf ( "  Ex  Ey     MC-Estimate           Exact      Error\n" );
  printf ( "\n" );

  for ( test = 1; test <= test_num; test++ )
  {
    e = i4vec_uniform_ab_new ( m, 0, 7, &seed );

    value = monomial_value ( m, n, e, x );

    result = square01_area ( ) * r8vec_sum ( n, value ) / ( double ) ( n );
    exact = square01_monomial_integral ( e );
    error = fabs ( result - exact );

    printf ( "  %2d  %2d  %14.6g  %14.6g  %10.2e\n",
      e[0], e[1], result, exact, error );

    free ( value );
  }

  free ( x );

  return;
}
