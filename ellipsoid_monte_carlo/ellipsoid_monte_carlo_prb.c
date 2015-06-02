# include <stdlib.h>
# include <stdio.h>

# include "ellipsoid_monte_carlo.h"

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for ELLIPSOID_MONTE_CARLO_PRB.

  Discussion:

    ELLIPSOID_MONTE_CARLO_PRB tests the ELLIPSOID_MONTE_CARLO library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 August 2014

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "ELLIPSOID_MONTE_CARLO_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the ELLIPSOID_MONTE_CARLO library.\n" );

  test01 ( );
  test02 ( );
  test03 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "ELLIPSOID_MONTE_CARLO_PRB\n" );
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

    TEST01 uses ELLIPSOID_SAMPLE on a 2D ellipse centered at (0,0).

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 August 2014

  Author:

    John Burkardt
*/
{
# define M 2

  double a[M*M] = {
    9.0, 1.0, 
    1.0, 4.0 };
  int e[M];
  int e_test[M*7] = {
    0, 0, 
    1, 0, 
    0, 1, 
    2, 0, 
    1, 1, 
    0, 2, 
    3, 0 };
  int i;
  int j;
  int m = M;
  int n;
  double r = 2.0;
  double result;
  int seed;
  double v[M] = { 0.0, 0.0 };
  double *value;
  double volume;
  double *x;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  Use ELLIPSOID_SAMPLE to estimate integrals\n" );
  printf ( "  in a 2D ellipse x' * A * x <= r^2.\n" );

  printf ( "\n" );
  r8_print ( r, "  Ellipsoid radius R:" );
  r8vec_print ( m, v, "  Ellipsoid center V:" );
  r8mat_print ( m, m, a, "  Ellipsoid matrix A:" );

  volume = ellipsoid_volume ( m, a, v, r );
  printf ( "\n" );
  r8_print ( volume, "  Ellipsoid volume:" );

  seed = 123456789;

  printf ( "\n" );
  printf ( "         N        1              X               Y  " );
  printf ( "             X^2               XY             Y^2             X^3\n" );
  printf ( "\n" );

  n = 1;

  while ( n <= 65536 )
  {
    x = ellipsoid_sample ( m, n, a, v, r, &seed );

    printf ( "  %8d", n );
    for ( j = 0; j < 7; j++ )
    {
      for ( i = 0; i < m; i++ )
      {
        e[i] = e_test[i+j*m];
      }
      value = monomial_value ( m, n, e, x );

      result = volume * r8vec_sum ( n, value ) / ( double ) ( n );
      printf ( "  %14.6g", result );
      free ( value );
    }

    printf ( "\n" );

    free ( x );

    n = 2 * n;
  }

  return;
# undef M
}
/******************************************************************************/

void test02 ( )

/******************************************************************************/
/*
  Purpose:

    TEST02 uses ELLIPSOID_SAMPLE on a 2D ellipse centered at (2,3).

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 August 2014

  Author:

    John Burkardt
*/
{
# define M 2

  double a[M*M] = {
    9.0, 1.0, 
    1.0, 4.0 };
  int e[M];
  int e_test[M*7] = {
    0, 0, 
    1, 0, 
    0, 1, 
    2, 0, 
    1, 1, 
    0, 2, 
    3, 0 };
  int i;
  int j;
  int m = M;
  int n;
  double r = 0.5;
  double result;
  int seed;
  double v[M] = { 2.0, 3.0 };
  double *value;
  double volume;
  double *x;

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  Use ELLIPSOID_SAMPLE to estimate integrals\n" );
  printf ( "  in a 2D ellipse (x-v)' * A * (x-v) <= r^2.\n" );

  printf ( "\n" );
  r8_print ( r, "  Ellipsoid radius R:" );
  r8vec_print ( m, v, "  Ellipsoid center V:" );
  r8mat_print ( m, m, a, "  Ellipsoid matrix A:" );

  volume = ellipsoid_volume ( m, a, v, r );
  printf ( "\n" );
  r8_print ( volume, "  Ellipsoid volume:" );

  seed = 123456789;

  printf ( "\n" );
  printf ( "         N        1              X               Y  " );
  printf ( "             X^2               XY             Y^2             X^3\n" );
  printf ( "\n" );

  n = 1;

  while ( n <= 65536 )
  {
    x = ellipsoid_sample ( m, n, a, v, r, &seed );

    printf ( "  %8d", n );
    for ( j = 0; j < 7; j++ )
    {
      for ( i = 0; i < m; i++ )
      {
        e[i] = e_test[i+j*m];
      }
      value = monomial_value ( m, n, e, x );

      result = volume * r8vec_sum ( n, value ) / ( double ) ( n );
      printf ( "  %14.6g", result );
      free ( value );
    }

    printf ( "\n" );

    free ( x );

    n = 2 * n;
  }

  return;
# undef M
}
/******************************************************************************/

void test03 ( )

/******************************************************************************/
/*
  Purpose:

    TEST03 uses ELLIPSOID_SAMPLE on a 3D ellipse centered at (1,2,3).

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 August 2014

  Author:

    John Burkardt
*/
{
# define M 3

  double a[M*M] = {
    9.0, 6.0, 3.0, 
    6.0, 5.0, 4.0, 
    3.0, 4.0, 9.0 };
  int e[M];
  int e_test[M*7] = {
    0, 0, 0, 
    1, 0, 0, 
    0, 1, 0, 
    0, 0, 1, 
    2, 0, 0, 
    0, 2, 2, 
    0, 0, 3 };
  int i;
  int j;
  int m = M;
  int n;
  double r = 0.5;
  double result;
  int seed;
  double v[M] = { 1.0, 2.0, 3.0 };
  double *value;
  double volume;
  double *x;

  printf ( "\n" );
  printf ( "TEST03\n" );
  printf ( "  Use ELLIPSOID_SAMPLE to estimate integrals\n" );
  printf ( "  in a 3D ellipse (x-v)' * A * (x-v) <= r^2.\n" );

  printf ( "\n" );
  r8_print ( r, "  Ellipsoid radius R:" );
  r8vec_print ( m, v, "  Ellipsoid center V:" );
  r8mat_print ( m, m, a, "  Ellipsoid matrix A:" );

  volume = ellipsoid_volume ( m, a, v, r );
  printf ( "\n" );
  r8_print ( volume, "  Ellipsoid volume:" );

  seed = 123456789;

  printf ( "\n" );
  printf ( "         N        1              X               Y  " );
  printf ( "              Z                X^2            YZ              Z^3\n" );
  printf ( "\n" );

  n = 1;

  while ( n <= 65536 )
  {
    x = ellipsoid_sample ( m, n, a, v, r, &seed );

    printf ( "  %8d", n );
    for ( j = 0; j < 7; j++ )
    {
      for ( i = 0; i < m; i++ )
      {
        e[i] = e_test[i+j*m];
      }
      value = monomial_value ( m, n, e, x );

      result = volume * r8vec_sum ( n, value ) / ( double ) ( n );
      printf ( "  %14.6g", result );
      free ( value );
    }

    printf ( "\n" );

    free ( x );

    n = 2 * n;
  }

  return;
# undef M
}
