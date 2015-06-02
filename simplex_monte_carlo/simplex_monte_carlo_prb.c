# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "simplex_monte_carlo.h"

int main ( );
void test01 ( );
void test02 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for SIMPLEX_MONTE_CARLO_PRB.

  Discussion:

    SIMPLEX_MONTE_CARLO_PRB tests the SIMPLEX_MONTE_CARLO library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 January 2014

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "SIMPLEX_MONTE_CARLO_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the SIMPLEX_MONTE_CARLO library.\n" );

  test01 ( );
  test02 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "SIMPLEX_MONTE_CARLO_PRB\n" );
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

    TEST01 uses SIMPLEX_SAMPLE_01 to estimate integrals in 3D.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 January 2014

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
  const int m = 3;
  int n;
  double result;
  int seed;
  double *value;
  double *x;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  Use SIMPLEX01_SAMPLE for a Monte Carlo estimate of an\n" );
  printf ( "  integral over the interior of the unit simplex in 3D.\n" );

  seed = 123456789;

  printf ( "\n" );
  printf ( "         N        1               X               Y " );
  printf ( "              Z               X^2              XY             XZ" );
  printf ( "              Y^2             YZ               Z^2\n" );
  printf ( "\n" );

  n = 1;

  while ( n <= 65536 )
  {
    x = simplex01_sample ( m, n, &seed );

    printf ( "  %8d", n );
    for ( j = 0; j < 10; j++ )
    {
      for ( i = 0; i < m; i++ )
      {
        e[i] = e_test[i+j*m];
      }
      value = monomial_value ( m, n, e, x );

      result = simplex01_volume ( m ) * r8vec_sum ( n, value ) 
        / ( double ) ( n );

      printf ( "  %14.6g", result );

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

    result = simplex01_monomial_integral ( m, e );
    printf ( "  %14.6g", result );
  }

  printf ( "\n" );

  return;
}
/******************************************************************************/

void test02 ( )

/******************************************************************************/
/*
  Purpose:

    TEST02 uses SIMPLEX_SAMPLE_01 to estimate integrals in 6D.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 January 2014

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
  double error;
  double exact;
  int i;
  int j;
  const int m = 6;
  int n;
  double result;
  int seed;
  double *value;
  double *x;

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  Use SIMPLEX01_SAMPLE for a Monte Carlo estimate of an\n" );
  printf ( "  integral over the interior of the unit simplex in 6D.\n" );

  seed = 123456789;

  printf ( "\n" );
  printf ( "         N" );
  printf ( "        1      " );
  printf ( "        U      " );
  printf ( "         V^2    " );
  printf ( "         V^2W^2 " );
  printf ( "         X^4    " );
  printf ( "         Y^2Z^2 " );
  printf ( "         Z^6\n" );
  printf ( "\n" );

  n = 1;

  while ( n <= 65536 )
  {
    x = simplex01_sample ( m, n, &seed );

    printf ( "  %8d", n );
    for ( j = 0; j < 7; j++ )
    {
      for ( i = 0; i < m; i++ )
      {
        e[i] = e_test[i+j*m];
      }
      value = monomial_value ( m, n, e, x );

      result = simplex01_volume ( m ) * r8vec_sum ( n, value ) 
        / ( double ) ( n );

      printf ( "  %14.6g", result );

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

    result = simplex01_monomial_integral ( m, e );
    printf ( "  %14.6g", result );
  }

  printf ( "\n" );

  return;
}

