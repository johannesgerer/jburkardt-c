# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "sphere_monte_carlo.h"

int main ( );
void test01 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for SPHERE_MONTE_CARLO_PRB.

  Discussion:

    SPHERE_MONTE_CARLO_PRB tests the SPHERE_MONTE_CARLO library.
    
  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 January 2014

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "SPHERE_MONTE_CARLO_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the SPHERE_MONTE_CARLO library.\n" );

  test01 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "SPHERE_MONTE_CARLO_PRB\n" );
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

    TEST01 uses SPHERE01_SAMPLE with an increasing number of points.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 January 2014

  Author:

    John Burkardt
*/
{
  int e[3];
  int e_test[3*7] = {
    0, 0, 0, 
    2, 0, 0, 
    0, 2, 0, 
    0, 0, 2, 
    4, 0, 0, 
    2, 2, 0, 
    0, 0, 4 };
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
  printf ( "  Use SPHERE01_SAMPLE to estimate integrals on the unit sphere surface.\n" );

  seed = 123456789;

  printf ( "\n" );
  printf ( "         N        1              X^2             Y^2" );
  printf ( "             Z^2             X^4            X^2Y^2          Z^4\n" );
  printf ( "\n" );

  n = 1;

  while ( n <= 65536 )
  {
    x = sphere01_sample ( n, &seed );
    printf ( "  %8d", n );
    for ( j = 0; j < 7; j++ )
    {
      for ( i = 0; i < 3; i++ )
      {
        e[i] = e_test[i+j*3];
      }

      value = monomial_value ( 3, n, e, x );

      result = sphere01_area ( ) * r8vec_sum ( n, value ) / ( double ) ( n );
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
    for ( i = 0; i < 3; i++ )
    {
      e[i] = e_test[i+j*3];
    }
    exact = sphere01_monomial_integral ( e );
    printf ( "  %14.10g", exact );
  }
  printf ( "\n" );

  return;
}
