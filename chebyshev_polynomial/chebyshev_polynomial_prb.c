# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "chebyshev_polynomial.h"

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );
void test04 ( );
void test05 ( );
void test06 ( );
void test07 ( );
void test08 ( );
void test09 ( );
void test10 ( );
void test11 ( );
void test12 ( );
void test13 ( );
void test14 ( );
void test15 ( );
void test16 ( );
void test17 ( );

/******************************************************************************/

int main ( )

/****************************************************************************/
/*
  Purpose:

    MAIN is the main program for CHEBYSHEV_POLYNOMIAL_PRB.

  Discussion:

    CHEBYSHEV_POLYNOMIAL_PRB tests the CHEBYSHEV_POLYNOMIAL library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 August 2013

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "CHEBYSHEV_POLYNOMIAL_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the CHEBYSHEV_POLYNOMIAL library.\n" );

  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
  test05 ( );
  test06 ( );
  test07 ( );
  test08 ( );
  test09 ( );

  test10 ( );
  test11 ( );
  test12 ( );
  test13 ( );
  test14 ( );
  test15 ( );
  test16 ( );
  test17 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "CHEBYSHEV_POLYNOMIAL_PRB\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void test01 ( )

/****************************************************************************/
/*
  Purpose:

    TEST01 tests T_PROJECT_COEFFICIENTS_DATA.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 August 2013

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double *c;
  double *d;
  double *d2;
  int i;
  int m;
  int n;
  int seed;
  double *x;

  printf ( "\n" );
  printf ( "CHEBYSHEV_POLYNOMIAL_TEST01:\n" );
  printf ( "  T_PROJECT_COEFFICIENTS_DATA estimates the Chebyshev polynomial\n" );
  printf ( "  coefficients for a function given as data (x,fx).\n" );
  printf ( "\n" );
  printf ( "  Here, we use fx = f(x) = x^2 for the data.\n" );
  printf ( "\n" );
  printf ( "  Since T(0,x) = 1 and T(2,x) = 2*x^2 - 1, the correct expansion is\n" );
  printf ( "  f(x) = 1/2 T(0,x) + 0 T(1,x) + 1/2 T(2,x) + 0 * all other polys.\n" );
/*
  Data in [0,1];
*/
  a = 0.0;
  b = 1.0;
  m = 20;
  seed = 123456789;
  x = r8vec_uniform_01_new ( m, &seed );
  d = ( double * ) malloc ( m * sizeof ( double ) );
  for ( i = 0; i < m; i++ )
  {
    d[i] = x[i] * x[i];
  }

  r8vec2_print ( m, x, d, "  Data ( X, D ):" );

  n = 4;
  c = t_project_coefficients_data ( a, b, m, n, x, d );

  r8vec_print ( n, c, "  Coefficients of Chebyshev expansion of degree 4." );
/*
  Compare Chebyshev expansion and original function.
*/
  d2 = t_project_value ( m, n, x, c );

  printf ( "\n" );
  printf ( "   I      X(I)     Data(I)      Chebyshev(X(I))\n" );
  printf ( "\n" );
  for ( i = 0; i < m; i++ )
  {
    printf ( "  %2d  %12g  %12g  %12g\n", i, x[i], d[i], d2[i] );
  }
  free ( c );
  free ( d );
  free ( d2 );
  free ( x );

  return;
}
/******************************************************************************/

void test02 ( )

/****************************************************************************/
/*
  Purpose:

    TEST02 tests T_POLYNOMIAL_COEFFICIENTS.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 August 2013

  Author:

    John Burkardt
*/
{
  double *c;
  int i;
  int j;
  int n;

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  T_POLYNOMIAL_COEFFICIENTS determines the polynomial coefficients \n" );
  printf ( "  ot T(n,x).\n" );

  n = 5;

  c = t_polynomial_coefficients ( n );

  for ( i = 0; i <= n; i++ )
  {
    printf ( "\n" );
    printf ( "  T(%d,x)\n", i );
    printf ( "\n" );
    for ( j = i; 0 <= j; j-- )
    {
      if ( c[i+j*(n+1)] != 0.0 )
      {
        if ( j == 0 )
        {
          printf ( "%14.6g\n", c[i+j*(n+1)] );
        }
        else if ( j == 1 )
        {
          printf ( "%14.6g * x\n", c[i+j*(n+1)] );
        }
        else
        {
          printf ( "%14.6g * x^%d\n", c[i+j*(n+1)] );
        }
      }
    }
  }
  return;
}
/******************************************************************************/

void test03 ( )

/****************************************************************************/
/*
  Purpose:

    TEST03 tests T_POLYNOMIAL.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 August 2013

  Author:

    John Burkardt
*/
{
  double fx;
  double *fx2;
  int n;
  int n_data;
  double x;
  double x_vec[1];

  printf ( "\n" );
  printf ( "TEST03:\n" );
  printf ( "  T_POLYNOMIAL evaluates the Chebyshev polynomial T(n,x).\n" );
  printf ( "\n" );
  printf ( "                   Tabulated      Computed\n" );
  printf ( "     N      X        T(n,x)        T(n,x)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    t_polynomial_values ( &n_data, &n, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    x_vec[0] = x;
    fx2 = t_polynomial ( 1, n, x_vec );

    printf ( "  %8d  %14g  %14g  %14g\n", n, x, fx, fx2[n] );

    free ( fx2 );

  }

  return;
}
/******************************************************************************/

void test04 ( )

/****************************************************************************/
/*
  Purpose:

    TEST04 tests T_POLYNOMIAL_ZEROS.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 August 2013

  Author:

    John Burkardt
*/
{
  double *fx;
  int i;
  int n;
  int n_max = 5;
  double *z;

  printf ( "\n" );
  printf ( "TEST04:\n" );
  printf ( "  T_POLYNOMIAL_ZEROS returns zeroes of T(n,x).\n" );
  printf ( "\n" );
  printf ( "       N      X        T(n,x)\n" );
  printf ( "\n" );

  for ( n = 1; n <= n_max; n++ )
  {
    z = t_polynomial_zeros ( n );
    fx = t_polynomial ( n, n, z );
    for ( i = 0; i < n; i++ )
    {
      printf ( "  %8d  %8g  %14g\n", n, z[i], fx[i+n*n] );
    }
    printf ( "\n" );
    free ( fx );
    free ( z );
  }

  return;
}
/******************************************************************************/

void test05 ( )

/****************************************************************************/
/*
  Purpose:

    TEST05 tests T_QUADRATURE_RULE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 August 2013

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
  printf ( "TEST05:\n" );
  printf ( "  T_QUADRATURE_RULE computes the quadrature rule\n" );
  printf ( "  associated with T(n,x);\n" );

  n = 7;
  x = ( double * ) malloc ( n * sizeof ( double ) );
  w = ( double * ) malloc ( n * sizeof ( double ) );

  t_quadrature_rule ( n, x, w );

  r8vec2_print ( n, x, w, "    N      X            W" );

  printf ( "\n" );
  printf ( "  Use the quadrature rule to estimate:\n" );
  printf ( "\n" );
  printf ( "    Q = Integral ( -1 <= X <= +1 ) X^E / sqrt ( 1-x^2) dx\n" );
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
    q_exact = t_integral ( e );
    printf ( "  %2d  %14g  %14g\n", e, q, q_exact );
  }

  free ( f );
  free ( w );
  free ( x );

  return;
}
/******************************************************************************/

void test06 ( )

/****************************************************************************/
/*
  Purpose:

    CHEBYSHEV_POLYNOMIAL_TEST06 tests the projection of T(i,x) and T(j,x).

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 August 2013

  Author:

    John Burkardt
*/
{
  double *c;
  int i;
  int j;
  int k;
  int n;
  double *phi;
  char title[80];
  double *w;
  double *x;

  printf ( "\n" );
  printf ( "TEST06:\n" );
  printf ( "  As a sanity check, make sure that the projection of:\n" );
  printf ( "  T(i,x) onto T(j,x) is:\n" );
  printf ( "  0 if i is not equal to j;\n" );
  printf ( "  pi if i = j = 0;\n" );
  printf ( "  pi/2 if i = j =/= 0.\n" );

  n = 3;

  x = ( double * ) malloc ( (n+1) * sizeof ( double ) );
  w = ( double * ) malloc ( (n+1) * sizeof ( double ) );

  t_quadrature_rule ( n + 1, x, w );

  c = ( double * ) malloc ( (n+1) * sizeof ( double ) );

  phi = t_polynomial ( n + 1, n, x );

  for ( j = 0; j <= n; j++ )
  {

    for ( i = 0; i <= n; i++ )
    {
      c[i] = 0.0;
      for ( k = 0; k <= n; k++ )
      {
        c[i] = c[i] + w[k] * phi[k+i*(n+1)] * phi[k+j*(n+1)];
      }
    }

    sprintf ( title, "  Chebyshev expansion coefficients for T(%d,x)", j );
    r8vec_print ( n + 1, c, title );
  }

  free ( c );
  free ( phi );
  free ( w );
  free ( x );

  return;
}
/******************************************************************************/

void test07 ( )

/****************************************************************************/
/*
  Purpose:

    CHEBYSHEV_POLYNOMIAL_TEST07 tests T_PROJECT_COEFFICIENTS.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 August 2013

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double *c;
  int n;

  printf ( "\n" );
  printf ( "CHEBYSHEV_POLYNOMIAL_TEST07:\n" );
  printf ( "  T_PROJECT_COEFFICIENTS computes the Chebyshev coefficients\n" );
  printf ( "  of a function defined over [-1,+1].\n" );
  printf ( "  T_PROJECT_COEFFICIENTS_AB works in [A,B].\n" );

  n = 3;
  c = ( double * ) malloc ( ( n + 1 ) * sizeof ( double ) );
  c = t_project_coefficients ( n, exp );
  r8vec_print ( n + 1, c, "  Chebyshev coefficients for exp(x) in [-1,+1]" );
  free ( c );

  n = 5;
  c = ( double * ) malloc ( ( n + 1 ) * sizeof ( double ) );
  c = t_project_coefficients ( n, exp );
  r8vec_print ( n + 1, c, "  Chebyshev coefficients for exp(x) in [-1,+1]" );
  free ( c );

  n = 5;
  c = ( double * ) malloc ( (n+1) * sizeof ( double ) );
  c = t_project_coefficients ( n, sin );
  r8vec_print ( n + 1, c, "  Chebyshev coefficients for sin(x) in [-1,+1]" );
  free ( c );
/*
  Repeat calculation with T_PROJECT_COEFFICIENTS_AB.
*/
  n = 5;
  c = ( double * ) malloc ( (n+1) * sizeof ( double ) );
  a = -1.0;
  b = +1.0;
  c = t_project_coefficients_ab ( n, sin, a, b );
  r8vec_print ( n + 1, c, "  Chebyshev coefficients for sin(x) in [-1,+1]" );
  free ( c );
/*
  Now try a different interval.
*/
  n = 5;
  c = ( double * ) malloc ( (n+1) * sizeof ( double ) );
  a = 0.0;
  b = 1.0;
  c = t_project_coefficients_ab ( n, sqrt, a, b );
  r8vec_print ( n + 1, c, "  Chebyshev coefficients for sqrt(x) in [0,+1]" );
  free ( c );

  return;
}
/******************************************************************************/

void test08 ( )

/****************************************************************************/
/*
  Purpose:

    TEST08 tests T_PROJECT_COEFFICIENTS_DATA.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 August 2013

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double *c;
  double *d;
  int i;
  int m;
  int n;
  int seed;
  double *x;

  printf ( "\n" );
  printf ( "TEST08:\n" );
  printf ( "  T_PROJECT_COEFFICIENTS_DATA computes the Chebyshev\n" );
  printf ( "  coefficients of a function defined by data.\n" );
  printf ( "\n" );
  printf ( "  We are looking for an approximation that is good in [-1,+1].\n" );
  printf ( "\n" );
  printf ( "  Begin by using equally spaced points in [-1,+1].\n" );

  a = -1.0;
  b = +1.0;
  m = 10;
  x = r8vec_linspace_new ( m, a, b );
  d = ( double * ) malloc ( m * sizeof ( double ) );
  for ( i = 0; i < m; i++ )
  {
    d[i] = exp ( x[i] );
  }
  n = 3;
  c = t_project_coefficients_data ( a, b, m, n, x, d );
  r8vec_print ( n + 1, c, "  Chebyshev coefficients for exp(x) on [-1,+1]" );
  free ( c );
  free ( d );
  free ( x );

  a = -1.0;
  b = +1.0;
  m = 10;
  x = r8vec_linspace_new ( m, a, b );
  d = ( double * ) malloc ( m * sizeof ( double ) );
  for ( i = 0; i < m; i++ )
  {
    d[i] = exp ( x[i] );
  }
  n = 5;
  c = t_project_coefficients_data ( a, b, m, n, x, d );
  r8vec_print ( n + 1, c, "  Chebyshev coefficients for exp(x) on [-1,+1]" );
  free ( c );
  free ( d );
  free ( x );

  a = -1.0;
  b = +1.0;
  m = 10;
  x = r8vec_linspace_new ( m, a, b );
  d = ( double * ) malloc ( m * sizeof ( double ) );
  for ( i = 0; i < m; i++ )
  {
    d[i] = sin ( x[i] );
  }
  n = 5;
  c = t_project_coefficients_data ( a, b, m, n, x, d );
  r8vec_print ( n + 1, c, "  Chebyshev coefficients for sin(x) on [-1,+1]" );
  free ( c );
  free ( d );
  free ( x );

  printf ( "\n" );
  printf ( "  Now sample equally spaced points in [0,+1].\n" );
  printf ( "  The approximation still applies to the interval [-1,+1].\n" );

  a = 0.0;
  b = +1.0;
  m = 10;
  x = r8vec_linspace_new ( m, a, b );
  d = ( double * ) malloc ( m * sizeof ( double ) );
  for ( i = 0; i < m; i++ )
  {
    d[i] = sin ( x[i] );
  }
  n = 5;
  c = t_project_coefficients_data ( a, b, m, n, x, d );
  r8vec_print ( n + 1, c, "  Chebyshev coefficients for sin(x) on [0,+1]" );
  free ( c );
  free ( d );
  free ( x );

  a = 0.0;
  b = +1.0;
  m = 10;
  x = r8vec_linspace_new ( m, a, b );
  d = ( double * ) malloc ( m * sizeof ( double ) );
  for ( i = 0; i < m; i++ )
  {
    d[i] = sqrt ( x[i] );
  }
  n = 5;
  c = t_project_coefficients_data ( a, b, m, n, x, d );
  r8vec_print ( n + 1, c, "  Chebyshev coefficients for sqrt(x) on [0,+1]" );
  free ( c );
  free ( d );
  free ( x );

  printf ( "\n" );
  printf ( "  Now random points in [-1,+1].\n" );

  a = -1.0;
  b = +1.0;
  m = 10;
  seed = 123456789;
  x = r8vec_uniform_ab_new ( m, a, b, &seed );
  d = ( double * ) malloc ( m * sizeof ( double ) );
  for ( i = 0; i < m; i++ )
  {
    d[i] = sin ( x[i] );
  }
  n = 5;
  c = t_project_coefficients_data ( a, b, m, n, x, d );
  r8vec_print ( n + 1, c, "  Chebyshev coefficients for sin(x) on [-1,+1]" );
  free ( c );
  free ( d );
  free ( x );

  return;
}
/******************************************************************************/

void test09 ( )

/****************************************************************************/
/*
  Purpose:

    TEST09 compares a function and projection over [-1,+1].

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 August 2013

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double *c;
  int i;
  int m;
  int n;
  double r;
  double *v;
  double *x;

  printf ( "\n" );
  printf ( "TEST09:\n" );
  printf ( "  T_PROJECT_COEFFICIENTS computes the Chebyshev interpolant C(F)(n,x)\n" );
  printf ( "  of a function F(x) defined over [-1,+1].\n" );
  printf ( "  T_PROJECT_VALUE evaluates that projection.\n" );

  printf ( "\n" );
  printf ( "  Compute projections of order N to exp(x) over [-1,+1],\n" );
  printf ( "\n" );
  printf ( "   N   Max||F(x)-C(F)(n,x)||\n" );
  printf ( "\n" );

  a = -1.0;
  b = +1.0;

  for ( n = 0; n <= 10; n++ )
  {
    c = t_project_coefficients ( n, exp );
    m = 101;
    x = r8vec_linspace_new ( m, a, b );
    v = t_project_value ( m, n, x, c );
    r = 0.0;
    for ( i = 0; i < m; i++ )
    {
      r = r8_max ( r, r8_abs ( v[i] - exp ( x[i] ) ) );
    }
    printf ( "  %2d  %14g\n", n, r );
    free ( c );
    free ( v );
    free ( x );
  }

  return;
}
/******************************************************************************/

void test10 ( )

/****************************************************************************/
/*
  Purpose:

    TEST10 compares a function and projection over [A,B].

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 August 2013

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double *c;
  int i;
  int m;
  int n;
  double r;
  double *v;
  double *x;

  printf ( "\n" );
  printf ( "TEST10:\n" );
  printf ( "  T_PROJECT_COEFFICIENTS_AB computes the Chebyshev interpolant C(F)(n,x)\n" );
  printf ( "  of a function F(x) defined over [A,B].\n" );
  printf ( "  T_PROJECT_VALUE_AB evaluates that projection.\n" );

  a = 0.0;
  b = 1.5;

  printf ( "\n" );
  printf ( "  Compute projections of order N to exp(x) over [%g,%g]\n", a, b );
  printf ( "\n" );
  printf ( "   N   Max||F(x)-C(F)(n,x)||\n" );
  printf ( "\n" );

  for ( n = 0; n <= 10; n++ )
  {
    c = t_project_coefficients_ab ( n, exp, a, b );
    m = 101;
    x = r8vec_linspace_new ( m, a, b );
    v = t_project_value_ab ( m, n, x, c, a, b );
    r = 0.0;
    for ( i = 0; i < m; i++ )
    {
      r = r8_max ( r, r8_abs ( v[i] - exp ( x[i] ) ) );
    }
    printf ( "  %2d  %14g\n", n, r );
    free ( c );
    free ( v );
    free ( x );
  }

  return;
}
/******************************************************************************/

void test11 ( )

/****************************************************************************/
/*
  Purpose:

    TEST11 tests U_POLYNOMIAL_COEFFICIENTS.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 August 2013

  Author:

    John Burkardt
*/
{
  double *c;
  int i;
  int j;
  int n;

  printf ( "\n" );
  printf ( "TEST11\n" );
  printf ( "  U_POLYNOMIAL_COEFFICIENTS determines the polynomial coefficients \n" );
  printf ( "  of U(n,x).\n" );

  n = 5;

  c = u_polynomial_coefficients ( n );

  for ( i = 0; i <= n; i++ )
  {
    printf ( "\n" );
    printf ( "  U(%d,x)\n", i );
    printf ( "\n" );
    for ( j = i; 0 <= j; j-- )
    {
      if ( c[i+j*(n+1)] != 0.0 )
      {
        if ( j == 0 )
        {
          printf ( "%14g\n", c[i+j*(n+1)] );
        }
        else if ( j == 1 )
        {
          printf ( "%14g * x\n", c[i+j*(n+1)] );
        }
        else
        {
          printf ( "%14g * x^%d\n", c[i+j*(n+1)], j );
        }
      }
    }
  }
  return;
}
/******************************************************************************/

void test12 ( )

/****************************************************************************/
/*
  Purpose:

    TEST12 tests U_POLYNOMIAL.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 August 2013

  Author:

    John Burkardt
*/
{
  double fx;
  double *fx2;
  int n;
  int n_data;
  double x;
  double x_vec[1];

  printf ( "\n" );
  printf ( "TEST12:\n" );
  printf ( "  U_POLYNOMIAL evaluates the Chebyshev polynomial U(n,x).\n" );
  printf ( "\n" );
  printf ( "                   Tabulated      Computed\n" );
  printf ( "     N      X        U(n,x)        U(n,x)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    u_polynomial_values ( &n_data, &n, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    x_vec[0] = x;
    fx2 = u_polynomial ( 1, n, x_vec );

    printf ( "  %8d  %14g  %14g  %14g\n", n, x, fx, fx2[n] );

    free ( fx2 );

  }

  return;
}
/******************************************************************************/

void test13 ( )

/****************************************************************************/
/*
  Purpose:

    TEST13 tests U_POLYNOMIAL_ZEROS.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 August 2013

  Author:

    John Burkardt
*/
{
  double *fx;
  int i;
  int n;
  int n_max = 5;
  double *z;

  printf ( "\n" );
  printf ( "TEST13:\n" );
  printf ( "  U_POLYNOMIAL_ZEROS returns zeroes of U(n,x).\n" );
  printf ( "\n" );
  printf ( "       N      X        U(n,x)\n" );
  printf ( "\n" );

  for ( n = 1; n <= n_max; n++ )
  {
    z = u_polynomial_zeros ( n );
    fx = u_polynomial ( n, n, z );
    for ( i = 0; i < n; i++ )
    {
      printf ( "  %8d  %8g  %14g\n", n, z[i], fx[i+n*n] );
    }
    printf ( "\n" );
    free ( fx );
    free ( z );
  }

  return;
}
/******************************************************************************/

void test14 ( )

/****************************************************************************/
/*
  Purpose:

    TEST14 tests U_QUADRATURE_RULE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 August 2013

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
  printf ( "TEST14:\n" );
  printf ( "  U_QUADRATURE_RULE computes the quadrature rule\n" );
  printf ( "  associated with U(n,x);\n" );

  n = 7;
  x = ( double * ) malloc ( n * sizeof ( double ) );
  w = ( double * ) malloc ( n * sizeof ( double ) );

  u_quadrature_rule ( n, x, w );

  r8vec2_print ( n, x, w, "    N      X            W" );

  printf ( "\n" );
  printf ( "  Use the quadrature rule to estimate:\n" );
  printf ( "\n" );
  printf ( "    Q = Integral ( -1 <= X <= +1 ) X^E * sqrt ( 1-x^2) dx\n" );
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
    q_exact = u_integral ( e );
    printf ( " %2d  %14g  %14g\n", e, q, q_exact );
  }

  free ( f );
  free ( w );
  free ( x );

  return;
}
/******************************************************************************/

void test15 ( )

/****************************************************************************/
/*
  Purpose:

    TEST15 tests V_POLYNOMIAL.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 August 2013

  Author:

    John Burkardt
*/
{
  double fx;
  double *fx2;
  int n;
  int n_data;
  double x;
  double x_vec[1];

  printf ( "\n" );
  printf ( "TEST15:\n" );
  printf ( "  V_POLYNOMIAL evaluates the Chebyshev polynomial V(n,x).\n" );
  printf ( "\n" );
  printf ( "                   Tabulated      Computed\n" );
  printf ( "     N      X        V(n,x)        V(n,x)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    v_polynomial_values ( &n_data, &n, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    x_vec[0] = x;
    fx2 = v_polynomial ( 1, n, x_vec );

    printf ( "  %8d  %8g  %14g  %14g\n", n, x, fx, fx2[n] );

    free ( fx2 );

  }

  return;
}
/******************************************************************************/

void test16 ( )

/****************************************************************************/
/*
  Purpose:

    TEST16 tests W_POLYNOMIAL.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 August 2013

  Author:

    John Burkardt
*/
{
  double fx;
  double *fx2;
  int n;
  int n_data;
  double x;
  double x_vec[1];

  printf ( "\n" );
  printf ( "TEST16:\n" );
  printf ( "  W_POLYNOMIAL evaluates the Chebyshev polynomial W(n,x).\n" );
  printf ( "\n" );
  printf ( "                   Tabulated      Computed\n" );
  printf ( "     N      X        W(n,x)        W(n,x)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    w_polynomial_values ( &n_data, &n, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    x_vec[0] = x;
    fx2 = w_polynomial ( 1, n, x_vec );

    printf ( "  %8d  %8g  %14g  %14g\n", n, x, fx, fx2[n] );

    free ( fx2 );

  }

  return;
}
/******************************************************************************/

void test17 ( )

/****************************************************************************/
/*
  Purpose:

    TEST17 tests T_TRIPLE_PRODUCT_INTEGRAL.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 August 2013

  Author:

    John Burkardt
*/
{
  double fx1;
  double fx2;
  int i;
  int j;
  int k;
  int l;
  int n;
  int seed;
  int test;
  int test_num = 20;
  double ti;
  double tj;
  double tk;
  double *w;
  double *x;

  printf ( "\n" );
  printf ( "TEST17:\n" );
  printf ( "  T_TRIPLE_PRODUCT_INTEGRAL computes the triple integral\n" );
  printf ( "  Tijk = integral ( -1 <= x <= 1 ) T(i,x) T(j,x) T(k,x) / sqrt ( 1-x^2) dx\n" );
  printf ( "\n" );
  printf ( "   I   J   K     Tijk           Tijk\n" );
  printf ( "                 computed       exact\n" );
  printf ( "\n" );

  n = 15;
  x = ( double * ) malloc ( n * sizeof ( double ) );
  w = ( double * ) malloc ( n * sizeof ( double ) );

  t_quadrature_rule ( n, x, w );

  seed = 123456789;

  for ( test = 1; test <= test_num; test++ )
  {
    i = i4_uniform_ab ( 2, 6, &seed );
    j = i4_uniform_ab ( 1, 3, &seed );
    k = i4_uniform_ab ( 0, 4, &seed );
    fx1 = t_triple_product_integral ( i, j, k );
    fx2 = 0.0;
    for ( l = 0; l < n; l++ )
    {
      ti = t_polynomial_value ( i, x[l] );
      tj = t_polynomial_value ( j, x[l] );
      tk = t_polynomial_value ( k, x[l] );
      fx2 = fx2 + w[l] * ti * tj * tk;
    }
    printf ( "  %2d  %2d  %2d  %14g  %14g\n", i, j, k, fx1, fx2 );
  }

  free ( x );
  free ( w );

  return;
}
