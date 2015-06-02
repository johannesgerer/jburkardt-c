# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "cube_monte_carlo.h"

int main ( );
void test01 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for CUBE_MONTE_CARLO_PRB.

  Discussion:

    CUBE_MONTE_CARLO_PRB tests the CUBE_MONTE_CARLO library.
    
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
  printf ( "CUBE_MONTE_CARLO_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the CUBE_MONTE_CARLO library.\n" );

  test01 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "CUBE_MONTE_CARLO_PRB\n" );
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

    TEST01 estimates integrals over the unit cube in 3D.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 January 2014

  Author:

    John Burkardt
*/
{
  int e[3];
  int e_test[3*10] = {
    0, 0, 0, 
    1, 0, 0, 
    0, 1, 0, 
    0, 0, 1, 
    2, 0, 0, 
    1, 1, 0, 
    1, 0, 1, 
    0, 2, 0, 
    0, 1, 1, 
    0, 0, 2 };
  double error;
  double exact;
  int i;
  int j;
  int m = 3;
  int n;
  double result;
  int seed;
  double *value;
  double *x;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  Use CUBE01_SAMPLE to estimate integrals\n" );
  printf ( "  over the interior of the unit cube in 3D.\n" );

  seed = 123456789;

  printf ( "\n" );
  printf ( "         N" );
  printf ( "        1" );
  printf ( "               X" );
  printf ( "               Y " );
  printf ( "              Z" );
  printf ( "               X^2" );
  printf ( "              XY" );
  printf ( "             XZ" );
  printf ( "              Y^2" );
  printf ( "             YZ" );
  printf ( "               Z^2\n" );
  printf ( "\n" );

  n = 1;

  while ( n <= 65536 )
  {
    x = cube01_sample ( n, &seed );
    printf ( "  %8d", n );
    for ( j = 0; j < 10; j++ )
    {
      for ( i = 0; i < m; i++ )
      {
        e[i] = e_test[i+j*m];
      }

      value = monomial_value ( m, n, e, x );

      result = cube01_volume ( ) * r8vec_sum ( n, value ) / ( double ) ( n );
      printf ( "  %14.10g", result );
      free ( value );
    }
    printf ( "\n" );

    free ( x );

    n = 2 * n;
  }

  printf ( "\n" );
  printf ( "     Exact" );
  for ( j = 0; j < 10; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      e[i] = e_test[i+j*m];
    }
    exact = cube01_monomial_integral ( e );
    printf ( "  %14.10g", exact );
  }
  printf ( "\n" );

  return;
}
