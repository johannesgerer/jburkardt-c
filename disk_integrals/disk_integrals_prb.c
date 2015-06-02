# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "disk_integrals.h"

int main ( );
void test01 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for DISK_INTEGRALS_PRB.

  Discussion:

    DISK_INTEGRALS_PRB tests the DISK_INTEGRALS library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 January 2014

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "DISK_INTEGRALS_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the DISK_INTEGRALS library.\n" );

  test01 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "DISK_INTEGRALS_PRB\n" );
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

    TEST01 uses DISK01_SAMPLE to compare estimated and exact integrals.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 January 2014

  Author:

    John Burkardt
*/
{
  int *e;
  double error;
  double exact;
  int i;
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
  printf ( "  Estimate monomial integrals using Monte Carlo\n" );
  printf ( "  over the interior of the unit disk in 2D.\n" );
/*
  Get sample points.
*/
  seed = 123456789;
  x = disk01_sample ( n, &seed );

  printf ( "\n" );
  printf ( "  Number of sample points used is %d\n", n );
/*
  Randomly choose X,Y exponents between 0 and 8.
*/
  printf ( "\n" );
  printf ( "  If any exponent is odd, the integral is zero.\n" );
  printf ( "  We will restrict this test to randomly chosen even exponents.\n" );
  printf ( "\n" );
  printf ( "  Ex  Ey     MC-Estimate           Exact      Error\n" );
  printf ( "\n" );

  for ( test = 1; test <= test_num; test++ )
  {
    e = i4vec_uniform_ab_new ( m, 0, 4, &seed );
    for ( i = 0; i < m; i++ )
    {
      e[i] = e[i] * 2;
    }
    value = monomial_value ( m, n, e, x );

    result = disk01_area ( ) * r8vec_sum ( n, value ) / ( double ) ( n );
    exact = disk01_monomial_integral ( e );
    error = fabs ( result - exact );

    printf ( "  %2d  %2d  %14.6g  %14.6g  %10.2e\n",
      e[0], e[1], result, exact, error );

    free ( e );
    free ( value );
  }

  free ( x );

  return;
}
