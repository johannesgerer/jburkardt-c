# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "circle_monte_carlo.h"

int main ( );
void test01 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for CIRCLE_MONTE_CARLO_PRB.

  Discussion:

    CIRCLE_MONTE_CARLO_PRB tests the CIRCLE_MONTE_CARLO library.
    
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
  printf ( "CIRCLE_MONTE_CARLO_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the CIRCLE_MONTE_CARLO library.\n" );

  test01 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "CIRCLE_MONTE_CARLO_PRB\n" );
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

    TEST01 uses CIRCLE01_SAMPLE with an increasing number of points.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 January 2014

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
  int n;
  double result;
  int seed;
  double *value;
  double *x;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  Use CIRCLE01_SAMPLE to estimate integrals\n" );
  printf ( "  along the circumference of the unit circle in 2D.\n" );

  seed = 123456789;

  printf ( "\n" );
  printf ( "         N        1              X^2             Y^2" );
  printf ( "             X^4           X^2Y^2          Y^4          X^6\n" );
  printf ( "\n" );

  n = 1;

  while ( n <= 65536 )
  {
    x = circle01_sample ( n, &seed );
    printf ( "  %8d", n );
    for ( j = 0; j < 7; j++ )
    {
      for ( i = 0; i < 2; i++ )
      {
        e[i] = e_test[i+j*2];
      }

      value = monomial_value ( 2, n, e, x );

      result = circle01_length ( ) * r8vec_sum ( n, value ) / ( double ) ( n );
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
    exact = circle01_monomial_integral ( e );
    printf ( "  %14.10g", exact );
  }
  printf ( "\n" );

  return;
}
