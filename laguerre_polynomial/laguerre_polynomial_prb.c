# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <string.h>

# include "laguerre_polynomial.h"

int main ( );
void laguerre_polynomial_test01 ( );
void laguerre_polynomial_test02 ( );
void laguerre_polynomial_test03 ( );
void laguerre_polynomial_test04 ( );
void laguerre_polynomial_test05 ( );
void laguerre_polynomial_test06 ( );
void laguerre_polynomial_test07 ( int p, double b );
void laguerre_polynomial_test08 ( int p, int e );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for LAGUERRE_POLYNOMIAL_PRB.

  Discussion:

    LAGUERRE_POLYNOMIAL_PRB tests the LAGUERRE_POLYNOMIAL library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 August 2013

  Author:

    John Burkardt
*/
{
  double b;
  int e;
  int p;

  timestamp ( );
  printf ( "\n" );
  printf ( "LAGUERRE_POLYNOMIAL_PRB:\n" );
  printf ( "  C version.\n" );
  printf ( "  Test the LAGUERRE_POLYNOMIAL library.\n" );

  laguerre_polynomial_test01 ( );
  laguerre_polynomial_test02 ( );
  laguerre_polynomial_test03 ( );
  laguerre_polynomial_test04 ( );
  laguerre_polynomial_test05 ( );
  laguerre_polynomial_test06 ( );

  p = 5;
  b = 0.0;
  laguerre_polynomial_test07 ( p, b );

  p = 5;
  b = 1.0;
  laguerre_polynomial_test07 ( p, b );

  p = 5;
  e = 0;
  laguerre_polynomial_test08 ( p, e );

  p = 5;
  e = 1;
  laguerre_polynomial_test08 ( p, e );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "LAGUERRE_POLYNOMIAL_PRB:\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void laguerre_polynomial_test01 ( )

/******************************************************************************/
/*
  Purpose:

    LAGUERRE_POLYNOMIAL_TEST01 tests L_POLYNOMIAL.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 August 2013

  Author:

    John Burkardt
*/
{
  int n_data;
  double e;
  double fx1;
  double fx2;
  double *fx2_vec;
  int n;
  double x;
  double x_vec[1];

  printf ( "\n" );
  printf ( "LAGUERRE_POLYNOMIAL_TEST01:\n" );
  printf ( "  L_POLYNOMIAL_VALUES stores values of\n" );
  printf ( "  the Laguerre polynomials.\n" );
  printf ( "  L_POLYNOMIAL evaluates the polynomial.\n" );
  printf ( "\n" );
  printf ( "                        Tabulated                 Computed\n" );
  printf ( "     N        X           L(N,X)                    L(N,X)                     Error\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    l_polynomial_values ( &n_data, &n, &x, &fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    x_vec[0] = x;
    fx2_vec = l_polynomial ( 1, n, x_vec );
    fx2 = fx2_vec[n];
    free ( fx2_vec );

    e = fx1 - fx2;

    printf ( "  %4d  %12g  %24.16g  %24.16g  %8g\n", n, x, fx1, fx2, e );
  }
  return;
}
/******************************************************************************/

void laguerre_polynomial_test02 ( )

/******************************************************************************/
/*
  Purpose:

    LAGUERRE_POLYNOMIAL_TEST02 tests L_POLYNOMIAL_COEFFICIENTS.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 August 2013

  Author:

    John Burkardt
*/
{
# define N 10

  double *c;
  int i;
  int j;

  printf ( "\n" );
  printf ( "LAGUERRE_POLYNOMIAL_TEST02\n" );
  printf ( "  L_POLYNOMIAL_COEFFICIENTS determines Laguerre polynomial coefficients.\n" );

  c = l_polynomial_coefficients ( N );

  for ( i = 0; i <= N; i++ )
  {
    printf ( "\n" );
    printf ( "  L(%d,x) =\n", i );
    printf ( "\n" );
    for ( j = i; 0 <= j; j-- )
    {
      if ( c[i+j*(N+1)] == 0.0 )
      {
      }
      else if ( j == 0 )
      {
        printf ( "%14g\n", c[i+j*(N+1)] );
      }
      else if ( j == 1 )
      {
        printf ( "%14g * x\n", c[i+j*(N+1)] );
      }
      else
      {
        printf ( "%14g * x^%d\n", c[i+j*(N+1)] , j );
      }
    }
  }

  free ( c );

  return;
# undef N
}
/******************************************************************************/

void laguerre_polynomial_test03 ( )

/******************************************************************************/
/*
  Purpose:

    LAGUERRE_POLYNOMIAL_TEST03 tests L_POLYNOMIAL_ZEROS.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 August 2013

  Author:

    John Burkardt
*/
{
  int degree;
  double *hz;
  char title[80];
  double *z;

  printf ( "\n" );
  printf ( "LAGUERRE_POLYNOMIAL_TEST03:\n" );
  printf ( "  L_POLYNOMIAL_ZEROS computes the zeros of L(n,x)\n" );
  printf ( "  Check by calling L_POLYNOMIAL there.\n" );

  for ( degree = 1; degree <= 5; degree++ )
  {
    z = l_polynomial_zeros ( degree );
    sprintf ( title, "  Computed zeros for L(%d,z):", degree );
    r8vec_print ( degree, z, title );

    hz = l_polynomial ( degree, degree, z );
    sprintf ( title, "  Evaluate L(%d,z):", degree );
    r8vec_print ( degree, hz+degree*degree, title );

    free ( hz );
    free ( z );
  }
  return;
}
/******************************************************************************/

void laguerre_polynomial_test04 ( )

/******************************************************************************/
/*
  Purpose:

    LAGUERRE_POLYNOMIAL_TEST04 tests L_QUADRATURE_RULE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 August 2013

  Author:

    John Burkardt
*/
{
  int e;
  double *f;
  int i;
  int n;
  double q;
  double q_exact;
  double *w;
  double *x;

  printf ( "\n" );
  printf ( "LAGUERRE_POLYNOMIAL_TEST04:\n" );
  printf ( "  L_QUADRATURE_RULE computes the quadrature rule\n" );
  printf ( "  associated with L(n,x)\n" );

  n = 7;
  x = ( double * ) malloc ( n * sizeof ( double ) );
  w = ( double * ) malloc ( n * sizeof ( double ) );

  l_quadrature_rule ( n, x, w );

  r8vec2_print ( n, x, w, "      X            W" );

  printf ( "\n" );
  printf ( "  Use the quadrature rule to estimate:\n" );
  printf ( "\n" );
  printf ( "    Q = Integral ( 0 <= X < +00 ) X^E exp(-X) dx\n" );
  printf ( "\n" );
  printf ( "   E       Q_Estimate      Q_Exact\n" );
  printf ( "\n" );

  f = ( double * ) malloc ( n * sizeof ( double ) );

  for ( e = 0; e <= 2 * n - 1; e++ )
  {
    if ( e == 0 )
    {
      for ( i = 0; i < n; i++ )
      {
        f[i] = 1.0;
      }
    }
    else
    {
      for ( i = 0; i < n; i++ )
      {
        f[i] = pow ( x[i], e );
      }
    }
    q = r8vec_dot_product ( n, w, f );
    q_exact = l_integral ( e );
    printf ( "  %2d  %14g  %14g\n", e, q, q_exact );
  }

  free ( f );
  free ( w );
  free ( x );

  return;
}
/******************************************************************************/

void laguerre_polynomial_test05 ( )

/******************************************************************************/
/*
  Purpose:

    LAGUERRE_POLYNOMIAL_TEST05 tests LM_POLYNOMIAL.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 August 2013

  Author:

    John Burkardt
*/
{
  int n_data;
  double e;
  double fx1;
  double fx2;
  double *fx2_vec;
  int m;
  int n;
  double x;
  double x_vec[1];

  printf ( "\n" );
  printf ( "LAGUERRE_POLYNOMIAL_TEST05:\n" );
  printf ( "  LM_POLYNOMIAL_VALUES stores values of\n" );
  printf ( "  the Laguerre polynomials.\n" );
  printf ( "  LM_POLYNOMIAL evaluates the polynomial.\n" );
  printf ( "\n" );
  printf ( "                              Tabulated                 Computed\n" );
  printf ( "     N     M        X         Lm(N,M,X)                  Lm(N,M,X)                     Error\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    lm_polynomial_values ( &n_data, &n, &m, &x, &fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    x_vec[0] = x;
    fx2_vec = lm_polynomial ( 1, n, m, x_vec );
    fx2 = fx2_vec[n];
    free ( fx2_vec );

    e = fx1 - fx2;

    printf ( "  %4d  %4d  %12g  %24.16g  %24.16g  %8g\n",
      n, m, x, fx1, fx2, e );
  }
  return;
}
/******************************************************************************/

void laguerre_polynomial_test06 ( )

/******************************************************************************/
/*
  Purpose:

    LAGUERRE_POLYNOMIAL_TEST06 tests LM_POLYNOMIAL_COEFFICIENTS.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 August 2013

  Author:

    John Burkardt
*/
{
# define N 5

  double *c;
  int i;
  int j;
  int m;

  printf ( "\n" );
  printf ( "LAGUERRE_POLYNOMIAL_TEST06\n" );
  printf ( "  LM_POLYNOMIAL_COEFFICIENTS determines Laguerre polynomial coefficients.\n" );

  for ( m = 0; m <= 4; m++ )
  {
    c = lm_polynomial_coefficients ( N, m );

    for ( i = 0; i <= N; i++ )
    {
      printf ( "\n" );
      printf ( "  Lm(%d,%d,x) =\n", i, m );
      printf ( "\n" );
      for ( j = i; 0 <= j; j-- )
      {
        if ( c[i+j*(N+1)] == 0.0 )
        {
        }
        else if ( j == 0 )
        {
          printf ( "%14g\n", c[i+j*(N+1)] );
        }
        else if ( j == 1 )
        {
          printf ( "%14g * x\n", c[i+j*(N+1)] );
        }
        else
        {
          printf ( "%14g * x^%d\n", c[i+j*(N+1)], j );
        }
      }
    }
    free ( c );
  }

  return;
# undef N
}
/******************************************************************************/

void laguerre_polynomial_test07 ( int p, double b )

/******************************************************************************/
/*
  Purpose:

    LAGUERRE_POLYNOMIAL_TEST07 tests L_EXPONENTIAL_PRODUCT.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 August 2013

  Author:

    John Burkardt

  Parameters:

    Input, int P, the maximum degree of the polynomial
    factors.

    Input, double B, the coefficient of X in the exponential factor.
*/
{
  double *table;

  printf ( "\n" );
  printf ( "LAGUERREE_POLYNOMIAL_TEST07\n" );
  printf ( "  Compute an exponential product table for L(n,x):\n" );
  printf ( "\n" );
  printf ( "  Tij = integral ( 0 <= x < +oo ) exp(b*x) Ln(i,x) Ln(j,x) exp(-x) dx\n" );
  printf ( "\n" );
  printf ( "  Maximum degree P = %d\n", p );
  printf ( "  Exponential argument coefficient B = %g\n", b );

  table = l_exponential_product ( p, b );

  r8mat_print ( p + 1, p + 1, table, "  Exponential product table:" );

  free ( table );

  return;
}
/******************************************************************************/

void laguerre_polynomial_test08 ( int p, int e )

/******************************************************************************/
/*
  Purpose:

    LAGUERRE_POLYNOMIAL_TEST08 tests L_POWER_PRODUCT.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 August 2013

  Author:

    John Burkardt

  Parameters:

    Input, int P, the maximum degree of the polynomial
    factors.

    Input, int E, the exponent of X.
*/
{
  double *table;

  printf ( "\n" );
  printf ( "LAGUERRE_POLYNOMIAL_TEST08\n" );
  printf ( "  Compute a power product table for L(n,x).\n" );
  printf ( "\n" );
  printf ( "  Tij = integral ( 0 <= x < +oo ) x^e L(i,x) L(j,x) exp(-x) dx\n" );

  printf ( "\n" );
  printf ( "  Maximum degree P = %d\n", p );
  printf ( "  Exponent of X, E = %d\n", e );

  table = l_power_product ( p, e );

  r8mat_print ( p + 1, p + 1, table, "  Power product table:" );

  free ( table );

  return;
}
