# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "hypersphere_monte_carlo.h"

int main ( );
void test01 ( );
void test02 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for HYPERSPHERE_MONTE_CARLO_PRB.

  Discussion:

    HYPERSPHERE_MONTE_CARLO_PRB tests the HYPERSPHERE_MONTE_CARLO library.
    
  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    04 January 2014

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "HYPERSPHERE_MONTE_CARLO_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the HYPERSPHERE_MONTE_CARLO library.\n" );

  test01 ( );
  test02 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "HYPERSPHERE_MONTE_CARLO_PRB\n" );
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

    TEST01 uses HYPERSPHERE01_SAMPLE to estimate integrals in 3D.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    04 January 2014

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
  printf ( "  Use HYPERSPHERE01_SAMPLE to estimate integrals \n" );
  printf ( "  on the surface of the unit hypersphere in M dimensions.\n" );
  printf ( "\n" );
  printf ( "  The spatial dimension M = %d\n", m );

  seed = 123456789;

  printf ( "\n" );
  printf ( "         N        1              X^2             Y^2" );
  printf ( "             Z^2             X^4            X^2Y^2          Z^4\n" );
  printf ( "\n" );

  n = 1;

  while ( n <= 65536 )
  {
    x = hypersphere01_sample ( m, n, &seed );
    printf ( "  %8d", n );
    for ( j = 0; j < 7; j++ )
    {
      for ( i = 0; i < m; i++ )
      {
        e[i] = e_test[i+j*m];
      }

      value = monomial_value ( m, n, e, x );

      result = hypersphere01_area ( m ) * r8vec_sum ( n, value ) 
        / ( double ) ( n );
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
    for ( i = 0; i < m; i++ )
    {
      e[i] = e_test[i+j*m];
    }
    result = hypersphere01_monomial_integral ( m, e );
    printf ( "  %14.10g", result );
  }
  printf ( "\n" );

  return;
}
/******************************************************************************/

void test02 ( )

/******************************************************************************/
/*
  Purpose:

    TEST02 uses HYPERSPHERE01_SAMPLE to estimate integrals in 6D.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 January 2014

  Author:

    John Burkardt
*/
{
  int e[6];
  int e_test[6*7] = {
    0, 0, 0, 0, 0, 0, 
    1, 0, 0, 0, 0, 0, 
    0, 2, 0, 0, 0, 0, 
    0, 2, 2, 0, 0, 0, 
    0, 0, 0, 4, 0, 0, 
    2, 0, 0, 0, 2, 2, 
    0, 0, 0, 0, 0, 6 };
  int i;
  int j;
  int m = 6;
  int n;
  double result;
  int seed;
  double *value;
  double *x;

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  Use HYPERSPHERE01_SAMPLE to estimate integrals \n" );
  printf ( "  on the surface of the unit hypersphere in M dimensions.\n" );
  printf ( "\n" );
  printf ( "  The spatial dimension M = %d\n", m );

  seed = 123456789;

  printf ( "\n" );
  printf ( "         N" );
  printf ( "        1      " );
  printf ( "        U      " );
  printf ( "         V^2   " );
  printf ( "         V^2W^2" );
  printf ( "         X^4   " );
  printf ( "         Y^2Z^2" );
  printf ( "         Z^6\n" );
  printf ( "\n" );

  n = 1;

  while ( n <= 65536 )
  {
    x = hypersphere01_sample ( m, n, &seed );
    printf ( "  %8d", n );
    for ( j = 0; j < 7; j++ )
    {
      for ( i = 0; i < m; i++ )
      {
        e[i] = e_test[i+j*m];
      }

      value = monomial_value ( m, n, e, x );

      result = hypersphere01_area ( m ) * r8vec_sum ( n, value ) 
        / ( double ) ( n );
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
    for ( i = 0; i < m; i++ )
    {
      e[i] = e_test[i+j*m];
    }
    result = hypersphere01_monomial_integral ( m, e );
    printf ( "  %14.10g", result );
  }
  printf ( "\n" );

  return;
}
