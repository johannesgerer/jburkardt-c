# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <string.h>

# include "pyramid_integrals.h"

int main ( );
void test01 ( );
void test02 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for PYRAMID_INTEGRALS_PRB.

  Discussion:

    PYRAMID_INTEGRALS_PRB tests the PYRAMID_INTEGRALS library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 April 2014

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "PYRAMID_INTEGRALS_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the PYRAMID_INTEGRALS library.\n" );

  test01 ( );
  test02 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "PYRAMID_INTEGRALS_PRB\n" );
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
  int e_max = 6;
  int e1;
  int e2;
  int e3;
  int expon[3];
  double error;
  double exact;
  int m = 3;
  static int n = 500000;
  double q;
  int seed;
  double *value;
  double *x;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  Compare exact and estimated integrals \n" );
  printf ( "  over the unit pyramid in 3D.\n" );
/*
  Get sample points.
*/
  seed = 123456789;
  x = pyramid01_sample ( n, &seed );

  printf ( "\n" );
  printf ( "  Number of sample points used is %d\n", n );
  printf ( "\n" );
  printf ( "   E1  E2  E3     MC-Estimate      Exact           Error\n" );
  printf ( "\n" );
/*
  Check all monomials, with only even dependence on X or Y, 
  up to total degree E_MAX.
*/
  for ( e3 = 0; e3 <= e_max; e3++ )
  {
    expon[2] = e3;
    for ( e2 = 0; e2 <= e_max - e3; e2 = e2 + 2 )
    {
      expon[1] = e2;
      for ( e1 = 0; e1 <= e_max - e3 - e2; e1 = e1 + 2 )
      {
        expon[0] = e1;

        value = monomial_value ( m, n, expon, x );

        q = pyramid01_volume ( ) * r8vec_sum ( n, value ) / ( double ) ( n );
        exact = pyramid01_integral ( expon );
        error = fabs ( q - exact );

        printf ( "  %2d  %2d  %2d  %14.6g  %14.6g  %10.2e\n",
          expon[0], expon[1], expon[2], q, exact, error );

        free ( value );
      }
    }
  }

  free ( x );

  return;
}
/******************************************************************************/

void test02 ( )

/******************************************************************************/
/*
  Purpose:

    TEST02 examines the sample points in the pyramid

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 April 2014

  Author:

    John Burkardt
*/
{
  int m = 3;
  int n = 20;
  int seed;
  double *x;

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  Print sample points in the unit pyramid in 3D.\n" );
  seed = 123456789;
  x = pyramid01_sample ( n, &seed );
  r8mat_transpose_print ( 3, n, x, "  Unit pyramid points" );

  free ( x );

  return;
}
