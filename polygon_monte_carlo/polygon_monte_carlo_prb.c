# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>

# include "polygon_monte_carlo.h"

int main ( );
void test01 ( int nv, double v[] );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for POLYGON_MONTE_CARLO_PRB.

  Discussion:

    POLYGON_MONTE_CARLO_PRB tests the POLYGON_MONTE_CARLO library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 May 2014

  Author:

    John Burkardt
*/
{
  int nv1 = 4;
  double v1[2*4] = {
    -1.0, -1.0, 
     1.0, -1.0, 
     1.0,  1.0, 
    -1.0,  1.0 };

  timestamp ( );
  printf ( "\n" );
  printf ( "POLYGON_MONTE_CARLO_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the POLYGON_MONTE_CARLO library.\n" );

  test01 ( nv1, v1 );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "POLYGON_MONTE_CARLO_PRB\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void test01 ( int nv, double v[] )

/******************************************************************************/
/*
  Purpose:

    TEST01 estimates integrals over a polygon in 2D.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 May 2014

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
  int j;
  int n;
  double result;
  int seed;
  double *value;
  double *x;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  Use POLYGON_SAMPLE to estimate integrals\n" );
  printf ( "  over the interior of a polygon in 2D.\n" );

  seed = 123456789;

  printf ( "\n" );
  printf ( "         N" );
  printf ( "        1" );
  printf ( "              X^2 " );
  printf ( "             Y^2" );
  printf ( "             X^4" );
  printf ( "           X^2Y^2" );
  printf ( "             Y^4" );
  printf ( "           X^6\n" );
  printf ( "\n" );

  n = 1;

  while ( n <= 65536 )
  {
    x = polygon_sample ( nv, v, n, &seed );

    printf ( "  %8d", n );

    for ( j = 0; j < 7; j++ )
    {
      e[0] = e_test[0+j*2];
      e[1] = e_test[1+j*2];

      value = monomial_value ( 2, n, e, x );

      result = polygon_area ( nv, v ) * r8vec_sum ( n, value ) / ( double ) ( n );
      printf ( "  %14.6g", result );
    }

    printf ( "\n" );

    free ( value );
    free ( x );

    n = 2 * n;

  }

  printf ( "\n" );
  printf ( "     Exact" );
  for ( j = 0; j < 7; j++ )
  {
    e[0] = e_test[0+j*2];
    e[1] = e_test[1+j*2];

    result = polygon_monomial_integral ( nv, v, e );
    printf ( "  %14.6g", result );
  }

  printf ( "\n" );

  return;
}
