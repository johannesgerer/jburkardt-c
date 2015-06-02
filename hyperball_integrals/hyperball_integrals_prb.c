# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>

# include "hyperball_integrals.h"

int main ( );
void test01 ( );
void test02 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for HYPERBALL_INTEGRALS_PRB.

  Discussion:

    HYPERBALL_INTEGRALS_PRB tests the HYPERBALL_INTEGRALS library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 January 2014

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "HYPERBALL_INTEGRALS_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the HYPERBALL_INTEGRALS library.\n" );

  test01 ( );
  test02 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "HYPERBALL_INTEGRALS_PRB\n" );
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

    TEST01 uses HYPERBALL01_SAMPLE to compare exact and estimated integrals in 3D.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 January 2014

  Author:

    John Burkardt
*/
{
  int *e;
  double error;
  double exact;
  int i;
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
  printf ( "  Use the Monte Carlo method to estimate integrals over\n" );
  printf ( "  the interior of the unit hyperball in M dimensions.\n" );
  printf ( "\n" );
  printf ( "  Spatial dimension M = %d", m );
/*
  Get sample points.
*/
  seed = 123456789;
  x = hyperball01_sample ( m, n, &seed );
  printf ( "\n" );
  printf ( "  Number of sample points used is %d\n", n );
/*
  Randomly choose exponents between 0 and 8.
*/
  printf ( "\n" );
  printf ( "  If any exponent is odd, the integral is zero.\n" );
  printf ( "  We will restrict this test to randomly chosen even exponents.\n" );
  printf ( "\n" );
  printf ( "  Ex  Ey  Ez     MC-Estimate           Exact      Error\n" );
  printf ( "\n" );

  for ( test = 1; test <= test_num; test++ )
  {
    e = i4vec_uniform_ab_new ( m, 0, 4, &seed );

    for ( i = 0; i < m; i++ )
    {
      e[i] = e[i] * 2;
    }
    value = monomial_value ( m, n, e, x );

    result = hyperball01_volume ( m ) * r8vec_sum ( n, value )
      / ( double ) ( n );
    exact = hyperball01_monomial_integral ( m, e );
    error = fabs ( result - exact );

    for ( i = 0; i < m; i++ )
    {
      printf ( "  %2d", e[i] );
    }
    printf ( "  %14.6g  %14.6g  %10.2e\n", result, exact, error );

    free ( e );
    free ( value );
  }

  free ( x );

  return;
}
/******************************************************************************/

void test02 ( )

/******************************************************************************/
/*
  Purpose:

    TEST02 uses HYPERBALL01_SAMPLE to compare exact and estimated integrals in 6D.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 January 2014

  Author:

    John Burkardt
*/
{
  int *e;
  double error;
  double exact;
  int i;
  int m = 6;
  int n = 4192;
  double result;
  int seed;
  int test;
  int test_num = 20;
  double *value;
  double *x;

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  Use the Monte Carlo method to estimate integrals over\n" );
  printf ( "  the interior of the unit hyperball in M dimensions.\n" );
  printf ( "\n" );
  printf ( "  Spatial dimension M = %d", m );
/*
  Get sample points.
*/
  seed = 123456789;
  x = hyperball01_sample ( m, n, &seed );
  printf ( "\n" );
  printf ( "  Number of sample points used is %d\n", n );
/*
  Randomly choose exponents between 0 and 6.
*/
  printf ( "\n" );
  printf ( "  If any exponent is odd, the integral is zero.\n" );
  printf ( "  We will restrict this test to randomly chosen even exponents.\n" );
  printf ( "\n" );
  printf ( "  E1  E2  E3  E4  E5  E6     MC-Estimate           Exact      Error\n" );
  printf ( "\n" );

  for ( test = 1; test <= test_num; test++ )
  {
    e = i4vec_uniform_ab_new ( m, 0, 3, &seed );

    for ( i = 0; i < m; i++ )
    {
      e[i] = e[i] * 2;
    }
    value = monomial_value ( m, n, e, x );

    result = hyperball01_volume ( m ) * r8vec_sum ( n, value )
      / ( double ) ( n );
    exact = hyperball01_monomial_integral ( m, e );
    error = fabs ( result - exact );

    for ( i = 0; i < m; i++ )
    {
      printf ( "  %2d", e[i] );
    }
    printf ( "  %14.6g  %14.6g  %10.2e\n", result, exact, error );

    free ( e );
    free ( value );
  }

  free ( x );

  return;
}
