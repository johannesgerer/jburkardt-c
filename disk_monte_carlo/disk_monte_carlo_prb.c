# include <stdlib.h>
# include <stdio.h>

# include "disk_monte_carlo.h"

int main ( );
void test01 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for DISK_MONTE_CARLO_PRB.

  Discussion:

    DISK_MONTE_CARLO_PRB tests the DISK_MONTE_CARLO library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    03 January 2014

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "DISK_MONTE_CARLO_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the DISK_MONTE_CARLO library.\n" );

  test01 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "DISK_MONTE_CARLO_PRB\n" );
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

    TEST01 uses DISK01_SAMPLE with an increasing number of points.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    03 January 2014

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
  double result[7];
  int seed;
  double *value;
  double *x;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  Use DISK01_SAMPLE to estimate integrals in the unit disk.\n" );

  seed = 123456789;

  printf ( "\n" );
  printf ( "         N        1              X^2             Y^2" );
  printf ( "             X^4             X^2Y^2           Y^4             X^6\n" );
  printf ( "\n" );

  n = 1;

  while ( n <= 65536 )
  {
    x = disk01_sample ( n, &seed );

    printf ( "  %8d", n );
    for ( j = 0; j < 7; j++ )
    {
      for ( i = 0; i < 2; i++ )
      {
        e[i] = e_test[i+j*2];
      }
      value = monomial_value ( 2, n, e, x );

      result[j] = disk01_area ( ) * r8vec_sum ( n, value ) / ( double ) ( n );
      printf ( "  %14.6g", result[j] );

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
    result[j] = disk01_monomial_integral ( e );
    printf ( "  %14.6g", result[j] );
  }
  printf ( "\n" );

  return;
}
