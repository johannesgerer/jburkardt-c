# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "line_integrals.h"

int main ( );
void test01 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for LINE_INTEGRALS_PRB.

  Discussion:

    LINE_INTEGRALS_PRB tests the LINE_INTEGRALS library.
    
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
  printf ( "LINE_INTEGRALS_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the LINE_INTEGRALS library.\n" );

  test01 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "LINE_INTEGRALS_PRB\n" );
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

    18 January 2014

  Author:

    John Burkardt
*/
{
  int e;
  double error;
  double exact;
  int m = 1;
  int n = 4192;
  double result;
  int seed;
  int test;
  int test_num = 11;
  double *value;
  double *x;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  Compare exact and estimated integrals \n" );
  printf ( "  over the length of the unit line in 1D.\n" );
/*
  Get sample points.
*/
  seed = 123456789;
  x = line01_sample ( n, &seed );

  printf ( "\n" );
  printf ( "  Number of sample points used is %d\n", n );
  printf ( "\n" );
  printf ( "   E     MC-Estimate      Exact           Error\n" );
  printf ( "\n" );

  for ( test = 1; test <= test_num; test++ )
  {
    e = test - 1;

    value = monomial_value_1d ( n, e, x );

    result = line01_length ( ) * r8vec_sum ( n, value ) / ( double ) ( n );
    exact = line01_monomial_integral ( e );
    error = fabs ( result - exact );

    printf ( "  %2d  %14.6g  %14.6g  %10.2e\n",
      e, result, exact, error );

    free ( value );
  }

  free ( x );

  return;
}
