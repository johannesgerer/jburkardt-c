# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "hypercube_integrals.h"

int main ( );
void test01 ( );
void test02 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for HYPERCUBE_INTEGRALS_PRB.

  Discussion:

    HYPERCUBE_INTEGRALS_PRB tests the HYPERCUBE_INTEGRALS library.
    
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
  printf ( "HYPERCUBE_INTEGRALS_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the HYPERCUBE_INTEGRALS library.\n" );

  test01 ( );
  test02 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "HYPERCUBE_INTEGRALS_PRB\n" );
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

    TEST01 estimates integrals over the unit hypercube in 3D.

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
  printf ( "  over the interior of the unit hypercube in 3D.\n" );
/*
  Get sample points.
*/
  seed = 123456789;
  x = hypercube01_sample ( m, n, &seed );
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
    e = i4vec_uniform_ab_new ( m, 0, 4, &seed );

    value = monomial_value ( m, n, e, x );

    result = hypercube01_volume ( m ) * r8vec_sum ( n, value ) / ( double ) ( n );
    exact = hypercube01_monomial_integral ( m, e );
    error = fabs ( result - exact );

    printf ( "  %2d  %2d  %2d  %14.6g  %14.6g  %10.2e\n",
      e[0], e[1], e[2], result, exact, error );

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

    TEST02 estimates integrals over the unit hypercube in 6D.

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
  printf ( "  Compare exact and estimated integrals\n" );
  printf ( "  over the interior of the unit hypercube in 6D.\n" );
/*
  Get sample points.
*/
  seed = 123456789;
  x = hypercube01_sample ( m, n, &seed );
  printf ( "\n" );
  printf ( "  Number of sample points is %d\n", n );
/*
  Randomly choose exponents.
*/
  printf ( "\n" );
  printf ( "  E1  E2  E3  E4  E5  E6     MC-Estimate           Exact      Error\n" );
  printf ( "\n" );

  for ( test = 1; test <= test_num; test++ )
  {
    e = i4vec_uniform_ab_new ( m, 0, 3, &seed );

    value = monomial_value ( m, n, e, x );

    result = hypercube01_volume ( m ) * r8vec_sum ( n, value ) / ( double ) ( n );
    exact = hypercube01_monomial_integral ( m, e );
    error = fabs ( result - exact );

    printf ( "  %2d  %2d  %2d  %2d  %2d  %2d  %14.6g  %14.6g  %10.2e\n",
      e[0], e[1], e[2], e[3], e[4], e[5], result, exact, error );

    free ( value );
  }

  free ( x );

  return;
}
