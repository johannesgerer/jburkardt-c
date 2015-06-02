# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "square_monte_carlo.h"

int main ( );
void test01 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for SQUARE_MONTE_CARLO_PRB.

  Discussion:

    SQUARE_MONTE_CARLO_PRB tests the SQUARE_MONTE_CARLO library.
    
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
  printf ( "SQUARE_MONTE_CARLO_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the SQUARE_MONTE_CARLO library.\n" );

  test01 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "SQUARE_MONTE_CARLO_PRB\n" );
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

    TEST01 estimates integrals over the unit square in 2D.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 January 2014

  Author:

    John Burkardt
*/
{
  int e[2];
  int e_test[2*7] = {
    0, 0, 
    2, 0, 
    0, 2, 
    4, 0, 
    2, 2, 
    0, 4,
    6, 0 };
  double error;
  double exact;
  int i;
  int j;
  int m = 2;
  int n;
  double result;
  int seed;
  double *value;
  double *x;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  Use SQUARE01_SAMPLE to estimate integrals\n" );
  printf ( "  over the interior of the unit square in 2D.\n" );

  seed = 123456789;

  printf ( "\n" );
  printf ( "         N        1              X^2             Y^2" );
  printf ( "             X^4           X^2Y^2          Y^4          X^6\n" );
  printf ( "\n" );

  n = 1;

  while ( n <= 65536 )
  {
    x = square01_sample ( n, &seed );
    printf ( "  %8d", n );
    for ( j = 0; j < 7; j++ )
    {
      for ( i = 0; i < 2; i++ )
      {
        e[i] = e_test[i+j*2];
      }

      value = monomial_value ( 2, n, e, x );

      result = square01_area ( ) * r8vec_sum ( n, value ) / ( double ) ( n );
      printf ( "  %14.10g", result );

      free ( value );
    }
    printf ( "\n" );

    free ( x );

    n = 2 * n;
  }

  printf ( "\n" );
  printf ( "     Exact" );
  for ( j = 0; j < 7; j++ )
  {
    for ( i = 0; i < 2; i++ )
    {
      e[i] = e_test[i+j*2];
    }
    exact = square01_monomial_integral ( e );
    printf ( "  %14.10g", exact );
  }
  printf ( "\n" );

  return;
}
