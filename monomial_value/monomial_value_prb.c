# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <string.h>
# include <time.h>

# include "monomial_value.h"

int main ( );
void test01 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for MONOMIAL_VALUE_PRB.

  Discussion:

    MONOMIAL_VALUE_PRB tests the MONOMIAL_VALUE library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 August 2014

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "MONOMIAL_VALUE_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the MONOMIAL_VALUE library.\n" );

  test01 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "MONOMIAL_VALUE_PRB\n" );
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

    TEST01 tests MONOMIAL_VALUE on sets of data in various dimensions.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 August 2014

  Author:

    John Burkardt
*/
{
  int *e;
  int e_max;
  int e_min;
  int i;
  int j;
  int m;
  int n;
  int seed;
  double *v;
  double *x;
  double x_max;
  double x_min;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  Usine monomial_value to evaluate monomials in\n" );
  printf ( "  dimensions 1 through 3.\n" ); 

  e_min = -3;
  e_max = 6;
  n = 5;
  seed = 123456789;
  x_min = -2.0;
  x_max = +10.0;

  for ( m = 1; m <= 3; m++ )
  {
    printf ( "\n" );
    printf ( "  Spatial dimension M = %d\n", m );

    e = i4vec_uniform_ab_new ( m, e_min, e_max, &seed );
    i4vec_transpose_print ( m, e, "  Exponents:" );
    x = r8mat_uniform_ab_new ( m, n, x_min, x_max, &seed );
/*
  To make checking easier, make the X values integers.
*/
    r8mat_nint ( m, n, x );
    v = monomial_value ( m, n, e, x );

    printf ( "\n" );
    printf ( "   V(X)         " );
    for ( i = 0; i < m; i++ )
    {
      printf ( "      X(%d)", i );
    }
    printf ( "\n" );
    printf ( "\n" );
    for ( j = 0; j < n; j++ )
    {
      printf ( "%14.6g  ", v[j] );
      for ( i = 0; i < m; i++ )
      {
        printf ( "%10.4f", x[i+j*m] );
      }
      printf ( "\n" );
    }

    free ( e );
    free ( x );
    free ( v );
  }

  return;
}
