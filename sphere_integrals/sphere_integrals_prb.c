# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "sphere_integrals.h"

int main ( );
void test01 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for SPHERE_INTEGRALS_PRB.

  Discussion:

    SPHERE_INTEGRALS_PRB tests the SPHERE_INTEGRALS library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    06 January 2014

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "SPHERE_INTEGRALS_PRB:\n" );
  printf ( "  C version\n" );
  printf ( "  Test the SPHERE_INTEGRALS library.\n" );

  test01 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "SPHERE_INTEGRALS_PRB\n" );
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

    TEST01 uses SPHERE01_SAMPLE to estimate monomial integrands.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    06 January 2014

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
  int n;
  double result;
  int seed;
  int test;
  const int test_num = 20;
  double *value;
  double *x;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  Estimate monomial integrands using Monte Carlo\n" );
  printf ( "  over the surface of the unit sphere in 3D.\n" );
/*
  Get sample points.
*/
  n = 8192;
  seed = 123456789;
  x = sphere01_sample ( n, &seed );
  printf ( "\n" );
  printf ( "  Number of sample points used is %d\n", n );
/*
  Randomly choose X,Y,Z exponents between (0,0,0) and (9,9,9).
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

    result = sphere01_area ( ) * r8vec_sum ( n, value ) / ( double ) ( n );
    exact = sphere01_monomial_integral ( e );
    error = fabs ( result - exact );

    printf ( "  %2d  %2d  %2d  %14.6g  %14.6g  %10.2e\n", 
      e[0], e[1], e[2], result, exact, error );

    free ( e );
    free ( value );
  }

  free ( x );

  return;
}

