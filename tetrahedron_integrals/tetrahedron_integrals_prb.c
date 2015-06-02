# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "tetrahedron_integrals.h"

int main ( );
void test01 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for TETRAHEDRON_INTEGRALS_PRB.

  Discussion:

    TETRAHEDRON_INTEGRALS_PRB tests the TETRAHEDRON_INTEGRALS library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 January 2014

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "TETRAHEDRON_INTEGRALS_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the TETRAHEDRON_INTEGRALS library.\n" );
/*
  Try each sampler on the unit tetrahedron, integrating quadratics.
*/
  test01 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "TETRAHEDRON_INTEGRALS_PRB\n" );
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

    TEST01 uses TETRAHEDRON_SAMPLE_01 to compare exact and estimated integrals.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 January 2014

  Author:

    John Burkardt
*/
{
  int e[3];
  double error;
  double exact;
  int i;
  int j;
  int k;
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
  printf ( "  Estimate monomial integrals using Monte Carlo\n" );
  printf ( "  over the interior of the unit tetrahedron in 3D.\n" );
/*
  Get sample points.
*/
  seed = 123456789;
  x = tetrahedron01_sample ( n, &seed );

  printf ( "\n" );
  printf ( "  Number of sample points used is %d\n", n );
/*
  Run through the exponents.
*/
  printf ( "\n" );
  printf ( "  Ex  Ey  Ez     MC-Estimate      Exact           Error\n" );
  printf ( "\n" );

  for ( i = 0; i <= 3; i++ )
  {
    e[0] = i;
    for ( j = 0; j <= 3; j++ )
    {
      e[1] = j;
      for ( k = 0; k <= 3; k++ )
      {
        e[2] = k;

        value = monomial_value ( m, n, e, x );

        result = tetrahedron01_volume ( ) * r8vec_sum ( n, value ) 
          / ( double ) ( n );
        exact = tetrahedron01_monomial_integral ( e );
        error = fabs ( result - exact );

        printf ( "  %2d  %2d  %2d  %14.6g  %14.6g  %10.2e\n",
          e[0], e[1], e[2], result, exact, error );

        free ( value );
      }
    }
  }

  free ( x );

  return;
}

