# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>

# include "wedge_monte_carlo.h"

int main ( );
void test01 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for WEDGE_MONTE_CARLO_PRB.

  Discussion:

    WEDGE_MONTE_CARLO_PRB tests the WEDGE_MONTE_CARLO library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    18 August 2014

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "WEDGE_MONTE_CARLO_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the WEDGE_MONTE_CARLO library.\n" );

  test01 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "WEDGE_MONTE_CARLO_PRB\n" );
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

    TEST01 uses WEDGE01_SAMPLE with an increasing number of points.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    18 August 2014

  Author:

    John Burkardt
*/
{
  int e[3];
  int e_test[3*8] = {
    0, 0, 0, 
    1, 0, 0, 
    0, 1, 0, 
    0, 0, 1, 
    2, 0, 0, 
    1, 1, 0, 
    0, 0, 2, 
    3, 0, 0 };
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
  printf ( "  Use WEDGE01_SAMPLE for a Monte Carlo estimate of an\n" );
  printf ( "  integral over the interior of the unit wedge in 3D.\n" );

  seed = 123456789;

  printf ( "\n" );
  printf ( "         N        1               X               Y " );
  printf ( "              Z                X^2            XY              Z^2    " );
  printf ( "        X^3\n" );
  printf ( "\n" );

  n = 1;

  while ( n <= 65536 )
  {
    x = wedge01_sample ( n, &seed );

    printf ( "  %8d", n );

    for ( j = 0; j < 8; j++ )
    {
      for ( i = 0; i < m; i++ )
      {
        e[i] = e_test[i+j*m];
      }
 
      value = monomial_value ( m, n, e, x );

      result = wedge01_volume ( ) * r8vec_sum ( n, value ) / ( double ) ( n );
      printf ( "  %14.6g", result );
      free ( value );
    }

    printf ( "\n" );

    free ( x );

    n = 2 * n;
  }

  printf ( "     Exact" );

  for ( j = 0; j < 8; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      e[i] = e_test[i+j*m];
    }
    result = wedge01_integral ( e );
    printf ( "  %14.6g", result );
  }
  printf ( "\n" );

  return;
}
