# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <string.h>

# include "hermite_polynomial.h"

int main ( );
void hermite_polynomial_test01 ( );
void hermite_polynomial_test02 ( );
void hermite_polynomial_test03 ( );
void hermite_polynomial_test04 ( );
void hermite_polynomial_test05 ( );
void hermite_polynomial_test06 ( );
void hermite_polynomial_test07 ( );
void hermite_polynomial_test08 ( int p, double b );
void hermite_polynomial_test09 ( int p, int e );
void hermite_polynomial_test10 ( int p, double b );
void hermite_polynomial_test11 ( int p, int e );
void hermite_polynomial_test12 ( int p, double b );
void hermite_polynomial_test13 ( int p, int e );
void hermite_polynomial_test14 ( );
void hermite_polynomial_test15 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for HERMITE_POLYNOMIAL_PRB.

  Discussion:

    HERMITE_POLYNOMIAL_PRB tests the HERMITE_POLYNOMIAL library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    10 August 2013

  Author:

    John Burkardt
*/
{
  double b;
  int e;
  int p;

  timestamp ( );
  printf ( "\n" );
  printf ( "HERMITE_POLYNOMIAL_PRB:\n" );
  printf ( "  C version.\n" );
  printf ( "  Test the HERMITE_POLYNOMIAL library.\n" );

  hermite_polynomial_test01 ( );
  hermite_polynomial_test02 ( );
  hermite_polynomial_test03 ( );
  hermite_polynomial_test04 ( );
  hermite_polynomial_test05 ( );
  hermite_polynomial_test06 ( );
  hermite_polynomial_test07 ( );

  p = 5;
  b = 0.0;
  hermite_polynomial_test08 ( p, b );

  p = 5;
  b = 1.0;
  hermite_polynomial_test08 ( p, b );

  p = 5;
  e = 0;
  hermite_polynomial_test09 ( p, e );

  p = 5;
  e = 1;
  hermite_polynomial_test09 ( p, e );

  p = 5;
  b = 0.0;
  hermite_polynomial_test10 ( p, b );

  p = 5;
  b = 1.0;
  hermite_polynomial_test10 ( p, b );

  p = 5;
  e = 0;
  hermite_polynomial_test11 ( p, e );

  p = 5;
  e = 1;
  hermite_polynomial_test11 ( p, e );

  p = 5;
  b = 0.0;
  hermite_polynomial_test12 ( p, b );

  p = 5;
  b = 1.0;
  hermite_polynomial_test12 ( p, b );

  p = 5;
  e = 0;
  hermite_polynomial_test13 ( p, e );

  p = 5;
  e = 1;
  hermite_polynomial_test13 ( p, e );

  hermite_polynomial_test14 ( );

  hermite_polynomial_test15 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "HERMITE_POLYNOMIAL_PRB:\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void hermite_polynomial_test01 ( )

/******************************************************************************/
/*
  Purpose:

    HERMITE_POLYNOMIAL_TEST01 tests H_POLYNOMIAL_VALUE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    10 August 2013

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
  printf ( "HERMITE_POLYNOMIAL_TEST01:\n" );
  printf ( "  H_POLYNOMIAL_VALUES stores values of\n" );
  printf ( "  the physicist's Hermite polynomials.\n" );
  printf ( "  H_POLYNOMIAL_VALUE evaluates the polynomial.\n" );
  printf ( "\n" );
  printf ( "                        Tabulated                 Computed\n" );
  printf ( "     N        X           H(N,X)                    H(N,X)                     Error\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    h_polynomial_values ( &n_data, &n, &x, &fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    x_vec[0] = x;
    fx2_vec = h_polynomial_value ( 1, n, x_vec );
    fx2 = fx2_vec[n];
    free ( fx2_vec );

    e = fx1 - fx2;

    printf ( "  %2d  %12g  %24.16g  %24.16g  %8g\n", n, x, fx1, fx2, e );
  }
  return;
}
/******************************************************************************/

void hermite_polynomial_test02 ( )

/******************************************************************************/
/*
  Purpose:

    HERMITE_POLYNOMIAL_TEST02 tests HE_POLYNOMIAL_VALUE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    10 August 2013

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
  printf ( "HERMITE_POLYNOMIAL_TEST02:\n" );
  printf ( "  HE_POLYNOMIAL_VALUES stores values of\n" );
  printf ( "  the probabilist's Hermite polynomials.\n" );
  printf ( "  HE_POLYNOMIAL_VALUE evaluates the polynomial.\n" );
  printf ( "\n" );
  printf ( "                        Tabulated                 Computed\n" );
  printf ( "     N        X          He(N,X)                   He(N,X)                     Error\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    he_polynomial_values ( &n_data, &n, &x, &fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    x_vec[0] = x;
    fx2_vec = he_polynomial_value ( 1, n, x_vec );
    fx2 = fx2_vec[n];
    free ( fx2_vec );

    e = fx1 - fx2;

    printf ( "  %4d  %12g  %24.16g  %24.16g  %8g\n", n, x, fx1, fx2, e );
  }
  return;
}
/******************************************************************************/

void hermite_polynomial_test03 ( )

/******************************************************************************/
/*
  Purpose:

    HERMITE_POLYNOMIAL_TEST03 tests HE_POLYNOMIAL.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    10 August 2013

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
  printf ( "HERMITE_POLYNOMIAL_TEST03:\n" );
  printf ( "  HF_FUNCTION_VALUES stores values of\n" );
  printf ( "  the Hermite function Hf(n,x).\n" );
  printf ( "  HF_FUNCTION_VALUE evaluates the function.\n" );
  printf ( "\n" );
  printf ( "                        Tabulated                 Computed\n" );
  printf ( "     N        X          Hf(N,X)                   Hf(N,X)                   Error\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    hf_function_values ( &n_data, &n, &x, &fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    x_vec[0] = x;
    fx2_vec = hf_function_value ( 1, n, x_vec );
    fx2 = fx2_vec[n];
    free ( fx2_vec );

    e = fx1 - fx2;

    printf ( "  %4d  %12g  %24.16g  %24.16g  %8g\n", n, x, fx1, fx2, e );
  }
  return;
}
/******************************************************************************/

void hermite_polynomial_test04 ( )

/******************************************************************************/
/*
  Purpose:

    HERMITE_POLYNOMIAL_TEST04 tests H_POLYNOMIAL_ZEROS.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    10 August 2013

  Author:

    John Burkardt
*/
{
  int degree;
  double *hz;
  char title[80];
  double *z;

  printf ( "\n" );
  printf ( "HERMITE_POLYNOMIAL_TEST04:\n" );
  printf ( "  H_POLYNOMIAL_ZEROS computes the zeros of H(n,x)\n" );
  printf ( "  Check by calling H_POLYNOMIAL there.\n" );

  for ( degree = 1; degree <= 5; degree++ )
  {
    z = h_polynomial_zeros ( degree );
    sprintf ( title, "  Computed zeros for H(%d,z):", degree );
    r8vec_print ( degree, z, title );

    hz = h_polynomial_value ( degree, degree, z );
    sprintf ( title, "  Evaluate H(%d,z):", degree );
    r8vec_print ( degree, hz+degree*degree, title );

    free ( hz );
    free ( z );
  }
  return;
}
/******************************************************************************/

void hermite_polynomial_test05 ( )

/******************************************************************************/
/*
  Purpose:

    HERMITE_POLYNOMIAL_TEST05 tests HE_POLYNOMIAL_ZEROS.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    10 August 2013

  Author:

    John Burkardt
*/
{
  int degree;
  double *hz;
  char title[80];
  double *z;

  printf ( "\n" );
  printf ( "HERMITE_POLYNOMIAL_TEST05:\n" );
  printf ( "  HE_POLYNOMIAL_ZEROS computes the zeros of He(n,x)\n" );
  printf ( "  Check by calling HE_POLYNOMIAL there.\n" );

  for ( degree = 1; degree <= 5; degree++ )
  {
    z = he_polynomial_zeros ( degree );
    sprintf ( title, "  Computed zeros for He(%d,z):", degree );
    r8vec_print ( degree, z, title );

    hz = he_polynomial_value ( degree, degree, z );
    sprintf ( title, "  Evaluate He(%d,z):", degree );
    r8vec_print ( degree, hz+degree*degree, title );

    free ( hz );
    free ( z );
  }
  return;
}
/******************************************************************************/

void hermite_polynomial_test06 ( )

/******************************************************************************/
/*
  Purpose:

    HERMITE_POLYNOMIAL_TEST06 tests H_QUADRATURE_RULE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    10 August 2013

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
  printf ( "HERMITE_POLYNOMIAL_TEST06:\n" );
  printf ( "  H_QUADRATURE_RULE computes the quadrature rule\n" );
  printf ( "  associated with H(n,x)\n" );

  n = 7;
  x = ( double * ) malloc ( n * sizeof ( double ) );
  w = ( double * ) malloc ( n * sizeof ( double ) );

  h_quadrature_rule ( n, x, w );

  r8vec2_print ( n, x, w, "      X            W" );

  printf ( "\n" );
  printf ( "  Use the quadrature rule to estimate:\n" );
  printf ( "\n" );
  printf ( "    Q = Integral ( -oo < X < +00 ) X^E exp(-X^2) dx\n" );
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
    q_exact = h_integral ( e );
    printf ( "  %2d  %14g  %14g\n", e, q, q_exact );
  }

  free ( f );
  free ( w );
  free ( x );

  return;
}
/******************************************************************************/

void hermite_polynomial_test07 ( )

/******************************************************************************/
/*
  Purpose:

    HERMITE_POLYNOMIAL_TEST07 tests HE_QUADRATURE_RULE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    10 August 2013

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
  printf ( "HERMITE_POLYNOMIAL_TEST07:\n" );
  printf ( "  HE_QUADRATURE_RULE computes the quadrature rule\n" );
  printf ( "  associated with He(n,x)\n" );

  n = 7;
  x = ( double * ) malloc ( n * sizeof ( double ) );
  w = ( double * ) malloc ( n * sizeof ( double ) );

  he_quadrature_rule ( n, x, w );

  r8vec2_print ( n, x, w, "      X            W" );

  printf ( "\n" );
  printf ( "  Use the quadrature rule to estimate:\n" );
  printf ( "\n" );
  printf ( "    Q = Integral ( -oo < X < +00 ) X^E exp(-X^2) dx\n" );
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
    q_exact = he_integral ( e );
    printf ( "  %2d  %14g  %14g\n", e, q, q_exact );
  }
  free ( f );
  free ( w );
  free ( x );

  return;
}
/******************************************************************************/

void hermite_polynomial_test08 ( int p, double b )

/******************************************************************************/
/*
  Purpose:

    HERMITE_POLYNOMIAL_TEST08 tests HN_EXPONENTIAL_PRODUCT.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    10 August 2013

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
  printf ( "HERMITE_POLYNOMIAL_TEST08\n" );
  printf ( "  Compute a normalized physicist''s Hermite exponential product table.\n" );
  printf ( "\n" );
  printf ( "  Tij = integral ( -oo < X < +oo ) exp(B*X) Hn(I,X) Hn(J,X) exp(-X*X) dx\n" );
  printf ( "\n" );
  printf ( "  where Hn(I,X) = normalized physicist''s Hermite polynomial of degree I.\n" );

  printf ( "\n" );
  printf ( "  Maximum degree P = %d\n", p );
  printf ( "  Exponential argument coefficient B = %g\n", b );

  table = hn_exponential_product ( p, b );

  r8mat_print ( p + 1, p + 1, table, "  Exponential product table:" );

  free ( table );

  return;
}
/******************************************************************************/

void hermite_polynomial_test09 ( int p, int e )

/******************************************************************************/
/*
  Purpose:

    HERMITE_POLYNOMIAL_TEST09 tests HN_POWER_PRODUCT.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    10 August 2013

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
  printf ( "HERMITE_POLYNOMIAL_TEST09\n" );
  printf ( "  Compute a normalized physicist''s Hermite power product table.\n" );
  printf ( "\n" );
  printf ( "  Tij = integral ( -oo < X < +oo ) X^E Hn(I,X) Hn(J,X) exp(-X*X) dx\n" );
  printf ( "\n" );
  printf ( "  where Hn(I,X) = normalized physicist''s Hermite polynomial of degree I.\n" );

  printf ( "\n" );
  printf ( "  Maximum degree P = %d\n", p );
  printf ( "  Exponent of X, E = %d\n", e );

  table = hn_power_product ( p, e );

  r8mat_print ( p + 1, p + 1, table, "  Power product table:" );

  free ( table );

  return;
}
/******************************************************************************/

void hermite_polynomial_test10 ( int p, double b )

/******************************************************************************/
/*
  Purpose:

    HERMITE_POLYNOMIAL_TEST10 tests HEN_EXPONENTIAL_PRODUCT.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    10 August 2013

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
  printf ( "HERMITE_POLYNOMIAL_TEST10\n" );
  printf ( "  Compute a normalized probabilist''s Hermite exponential product table.\n" );
  printf ( "\n" );
  printf ( "  Tij = integral ( -oo < X < +oo ) exp(B*X) Hen(I,X) Hen(J,X) exp(-0.5*X*X) dx\n" );
  printf ( "\n" );
  printf ( "  where Hen(I,X) = normalized probabilist''s Hermite polynomial of degree I.\n" );

  printf ( "\n" );
  printf ( "  Maximum degree P = %d\n", p );
  printf ( "  Exponential argument coefficient B = %g\n", b );

  table = hen_exponential_product ( p, b );

  r8mat_print ( p + 1, p + 1, table, "  Exponential product table:" );

  free ( table );

  return;
}
/******************************************************************************/

void hermite_polynomial_test11 ( int p, int e )

/******************************************************************************/
/*
  Purpose:

    HERMITE_POLYNOMIAL_TEST11 tests HEN_POWER_PRODUCT.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    10 August 2013

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
  printf ( "HERMITE_POLYNOMIAL_TEST11\n" );
  printf ( "  Compute a normalized probabilist''s Hermite power product table.\n" );
  printf ( "\n" );
  printf ( "  Tij = integral ( -oo < X < +oo ) X^E Hen(I,X) Hen(J,X) exp(-X*X) dx\n" );
  printf ( "\n" );
  printf ( "  where Hen(I,X) = normalized probabilist''s Hermite polynomial of degree I.\n" );

  printf ( "\n" );
  printf ( "  Maximum degree P = %d\n", p );
  printf ( "  Exponent of X, E = %d\n", e );

  table = hen_power_product ( p, e );

  r8mat_print ( p + 1, p + 1, table, "  Power product table:" );

  free ( table );

  return;
}
/******************************************************************************/

void hermite_polynomial_test12 ( int p, double b )

/******************************************************************************/
/*
  Purpose:

    HERMITE_POLYNOMIAL_TEST12 tests HF_EXPONENTIAL_PRODUCT.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    10 August 2013

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
  printf ( "HERMITE_POLYNOMIAL_TEST12\n" );
  printf ( "  Compute a Hermite function exponential product table.\n" );
  printf ( "\n" );
  printf ( "  Tij = integral ( -oo < X < +oo ) exp(B*X) Hf(I,X) Hf(J,X) dx\n" );
  printf ( "\n" );
  printf ( "  where Hf(I,X) = Hermite function of \"degree\" I.\n" );

  printf ( "\n" );
  printf ( "  Maximum degree P = %d\n", p );
  printf ( "  Exponential argument coefficient B = %g\n", b );

  table = hf_exponential_product ( p, b );

  r8mat_print ( p + 1, p + 1, table, "  Exponential product table:" );

  free ( table );

  return;
}
/******************************************************************************/

void hermite_polynomial_test13 ( int p, int e )

/******************************************************************************/
/*
  Purpose:

    HERMITE_POLYNOMIAL_TEST13 tests HF_POWER_PRODUCT.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    10 August 2013

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
  printf ( "HERMITE_POLYNOMIAL_TEST13\n" );
  printf ( "  Compute a Hermite function product table.\n" );
  printf ( "\n" );
  printf ( "  Tij = integral ( -oo < X < +oo ) X^E Hf(I,X) Hf(J,X) exp(-X*X) dx\n" );
  printf ( "\n" );
  printf ( "  where Hf(I,X) = Hermite function of \"degree\" I.\n" );

  printf ( "\n" );
  printf ( "  Maximum degree P = %d\n", p );
  printf ( "  Exponent of X, E = %d\n", e );

  table = hf_power_product ( p, e );

  r8mat_print ( p + 1, p + 1, table, "  Power product table:" );

  free ( table );

  return;
}
/******************************************************************************/

void hermite_polynomial_test14 ( )

/******************************************************************************/
/*
  Purpose:

    HERMITE_POLYNOMIAL_TEST14 tests H_POLYNOMIAL_COEFFICIENTS.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    10 August 2013

  Author:

    John Burkardt
*/
{
  double *c;
  int i;
  int j;
  int n = 10;

  printf ( "\n" );
  printf ( "HERMITE_POLYNOMIAL_TEST14\n" );
  printf ( "  H_POLYNOMIAL_COEFFICIENTS determines physicist's Hermite polynomial coefficients.\n" );

  c = h_polynomial_coefficients ( n );

  for ( i = 0; i <= n; i++ )
  {
    printf ( "\n" );
    printf ( "  H(%d,x) =\n", i );
    printf ( "\n" );
    for ( j = i; 0 <= j; j-- )
    {
      if ( c[i+j*(n+1)] == 0.0 )
      {
      }
      else if ( j == 0 )
      {
        printf ( "%g\n", c[i+j*(n+1)] );
      }
      else if ( j == 1 )
      {
        printf ( "%g * x\n", c[i+j*(n+1)] );
      }
      else
      {
        printf ( "%g * x^%d\n", c[i+j*(n+1)], j );
      }
    }
  }
  free ( c );

  return;
}
/******************************************************************************/

void hermite_polynomial_test15 ( )

/******************************************************************************/
/*
  Purpose:

    HERMITE_POLYNOMIAL_TEST15 tests HE_POLYNOMIAL_COEFFICIENTS.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    10 August 2013

  Author:

    John Burkardt
*/
{
  double *c;
  int i;
  int j;
  int n = 10;

  printf ( "\n" );
  printf ( "HERMITE_POLYNOMIAL_TEST15\n" );
  printf ( "  HE_POLYNOMIAL_COEFFICIENTS determines probabilist's Hermite polynomial coefficients.\n" );

  c = he_polynomial_coefficients ( n );

  for ( i = 0; i <= n; i++ )
  {
    printf ( "\n" );
    printf ( "  He(%d) =\n", i );
    printf ( "\n" );
    for ( j = i; 0 <= j; j-- )
    {
      if ( c[i+j*(n+1)] == 0.0 )
      {
      }
      else if ( j == 0 )
      {
        printf ( "%g\n", c[i+j*(n+1)] );
      }
      else if ( j == 1 )
      {
        printf ( "%g * x\n", c[i+j*(n+1)] );
      }
      else
      {
        printf ( "%g * x^%d\n", c[i+j*(n+1)], j );
      }
    }
  }
  free ( c );

  return;
}

