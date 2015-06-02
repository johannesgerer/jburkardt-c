# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "hermite.h"

int main ( void );
void test01 ( void );
void test02 ( void );
void test03 ( void );
void test04 ( void );
void test05 ( void );
void test06 ( void );
void test07 ( void );
void test08 ( void );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for HERMITE_PRB.

  Discussion:

    HERMITE_PRB tests the HERMITE library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 June 2011

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "HERMITE_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the HERMITE library.\n" );

  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
  test05 ( );
  test06 ( );
  test07 ( );
  test08 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "HERMITE_PRB\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void test01 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST01 uses f(x) = 1 + 2x + 3x^2 at x = 0, 1, 2.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 May 2011

  Author:

    John Burkardt
*/
{
# define N 3

  int n = N;
  double x[N] = { 0.0, 1.0,  2.0 };
  double y[N] = { 1.0, 6.0, 17.0 };
  double yp[N] = { 2.0, 8.0, 14.0 };

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  HERMITE computes the Hermite interpolant to data.\n" );
  printf ( "  Here, f(x) = 1 + 2x + 3x^2.\n" );

  hermite_demo ( n, x, y, yp );

  return;
# undef N
}
/******************************************************************************/

void test02 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST02 uses f(x) = 6 + 5x + 4x^2 + 3x^3 + 2x^4 + x^5 at x = 0, 1, 2.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 May 2011

  Author:

    John Burkardt
*/
{
# define N 3

  int i;
  int n = N;
  double *x;
  double *y;
  double *yp;

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  HERMITE computes the Hermite interpolant to data.\n" );
  printf ( "  Here, f(x) = 6 + 5x + 4x^2 + 3x^3 + 2x^4 + x^5.\n" );

  x = ( double * ) malloc ( n * sizeof ( double ) );
  y = ( double * ) malloc ( n * sizeof ( double ) );
  yp = ( double * ) malloc ( n * sizeof ( double ) );

  for ( i = 0; i < n; i++ )
  {
    x[i] = ( double ) ( i );
    
    y[i] = 6.0 + x[i] * ( 
           5.0 + x[i] * ( 
           4.0 + x[i] * ( 
           3.0 + x[i] * ( 
           2.0 + x[i] ) ) ) );

    yp[i] = 5.0 + x[i] * ( 
            8.0 + x[i] * ( 
            9.0 + x[i] * ( 
            8.0 + x[i] *   
            5.0 ) ) );
  }
  hermite_demo ( n, x, y, yp );

  free ( x );
  free ( y );
  free ( yp );

  return;
# undef N
}
/******************************************************************************/

void test03 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST03 uses f(x) = r1 + r2x + r3x^2 + r4x^3 + r5x^4 + r6x^5 at x = r7 r8 r9

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 May 2011

  Author:

    John Burkardt
*/
{
# define N 3

  double *c;
  int i;
  int n = N;
  int seed;
  double *x;
  double *y;
  double *yp;

  printf ( "\n" );
  printf ( "TEST03\n" );
  printf ( "  HERMITE computes the Hermite interpolant to data.\n" );
  printf ( "  Here, f(x) is a fifth order polynomial with random\n" );
  printf ( "  coefficients, and the abscissas are random.\n" );

  c = ( double * ) malloc ( 2 * n * sizeof ( double ) );
  x = ( double * ) malloc ( n * sizeof ( double ) );
  y = ( double * ) malloc ( n * sizeof ( double ) );
  yp = ( double * ) malloc ( n * sizeof ( double ) );

  seed = 123456789;

  r8vec_uniform_01 ( n, &seed, x );
  r8vec_print ( n, x, "  Random abscissas" );

  r8vec_uniform_01 ( 2 * n, &seed, c );
  r8vec_print ( 2 * n, c, "  Random polynomial coefficients." );

  for ( i = 0; i < n; i++ )
  {
    y[i] = c[0] + x[i] * ( 
           c[1] + x[i] * ( 
           c[2] + x[i] * ( 
           c[3] + x[i] * ( 
           c[4] + x[i] * ( 
           c[5] ) ) ) ) );

    yp[i] = c[1]       + x[i] * ( 
            c[2] * 2.0 + x[i] * ( 
            c[3] * 3.0 + x[i] * ( 
            c[4] * 4.0 + x[i] *   
            c[5] * 5.0 ) ) );
  }

  hermite_demo ( n, x, y, yp );

  free ( c );
  free ( x );
  free ( y );
  free ( yp );

  return;
# undef N
}
/******************************************************************************/

void test04 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST04 interpolates the Runge function with equally spaced data.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 October 2011

  Author:

    John Burkardt
*/
{
  int i;
  double max_dif;
  int n;
  int nd;
  int ndp;
  int ns;
  double *x;
  double *xd;
  double *xdp;
  double *xs;
  double xhi;
  double xlo;
  double xt;
  double *y;
  double *yd;
  double *ydp;
  double *yp;
  double *ys;
  double yt;

  printf ( "\n" );
  printf ( "TEST04\n" );
  printf ( "  HERMITE computes the Hermite interpolant to data.\n" );
  printf ( "  Here, f(x) is the Runge function\n" );
  printf ( "  and the data is evaluated at equally spaced points.\n" );
  printf ( "  As N increases, the maximum error grows.\n" );
  printf ( "\n" );
  printf ( "     N     Max | F(X) - H(F(X)) |\n" );
  printf ( "\n" );

  for ( n = 3; n <= 15; n = n + 2 )
  {
    y = ( double * ) malloc ( n * sizeof ( double ) );
    yp = ( double * ) malloc ( n * sizeof ( double ) );

    nd = 2 * n;

    xd = ( double * ) malloc ( nd * sizeof ( double ) );
    yd = ( double * ) malloc ( nd * sizeof ( double ) );

    ndp = 2 * n - 1;

    xdp = ( double * ) malloc ( ndp * sizeof ( double ) );
    ydp = ( double * ) malloc ( ndp * sizeof ( double ) );

    ns = 10 * ( n - 1 ) + 1;

    xlo = -5.0;
    xhi = +5.0;
    x = r8vec_linspace_new ( n, xlo, xhi );

    for ( i = 0; i < n; i++ )
    {
      y[i] = 1.0 / ( 1.0 + x[i] * x[i] );
      yp[i] = - 2.0 * x[i] / ( 1.0 + x[i] * x[i] ) / ( 1.0 + x[i] * x[i] );
    }

    hermite_interpolant ( n, x, y, yp, xd, yd, xdp, ydp );
/*
  Compare exact and interpolant at sample points.
*/
    xs = r8vec_linspace_new ( ns, xlo, xhi );

    ys = dif_vals ( nd, xd, yd, ns, xs );

    max_dif = 0.0;
    for ( i = 0; i < ns; i++ )
    {
      xt = xs[i];
      yt = 1.0 / ( 1.0 + xt * xt );
      max_dif = r8_max ( max_dif, r8_abs ( ys[i] - yt ) );
    }

    printf ( "  %4d  %14g\n", n, max_dif );

    free ( x );
    free ( xd );
    free ( xdp );
    free ( xs );
    free ( y );
    free ( yd );
    free ( ydp );
    free ( yp );
    free ( ys );
  }

  return;
}
/******************************************************************************/

void test05 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST05 interpolates the Runge function with Chebyshev spaced data.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 October 2011

  Author:

    John Burkardt
*/
{
  int i;
  double max_dif;
  int n;
  int nd;
  int ndp;
  int ns;
  double *x;
  double *xd;
  double *xdp;
  double *xs;
  double xhi;
  double xlo;
  double xt;
  double *y;
  double *yd;
  double *ydp;
  double *yp;
  double *ys;
  double yt;

  printf ( "\n" );
  printf ( "TEST05\n" );
  printf ( "  HERMITE computes the Hermite interpolant to data.\n" );
  printf ( "  Here, f(x) is the Runge function\n" );
  printf ( "  and the data is evaluated at Chebyshev spaced points.\n" );
  printf ( "  As N increases, the maximum error decreases.\n" );
  printf ( "\n" );
  printf ( "     N     Max | F(X) - H(F(X)) |\n" );
  printf ( "\n" );

  for ( n = 3; n <= 15; n = n + 2 )
  {
    y = ( double * ) malloc ( n * sizeof ( double ) );
    yp = ( double * ) malloc ( n * sizeof ( double ) );

    nd = 2 * n;

    xd = ( double * ) malloc ( nd * sizeof ( double ) );
    yd = ( double * ) malloc ( nd * sizeof ( double ) );

    ndp = 2 * n - 1;

    xdp = ( double * ) malloc ( ndp * sizeof ( double ) );
    ydp = ( double * ) malloc ( ndp * sizeof ( double ) );


    ns = 10 * ( n - 1 ) + 1;

    xlo = -5.0;
    xhi = +5.0;
    x = r8vec_chebyshev_new ( n, xlo, xhi );

    for ( i = 0; i < n; i++ )
    {
      y[i] = 1.0 / ( 1.0 + x[i] * x[i] );
      yp[i] = - 2.0 * x[i] / ( 1.0 + x[i] * x[i] ) / ( 1.0 + x[i] * x[i] );
    }

    hermite_interpolant ( n, x, y, yp, xd, yd, xdp, ydp );
/*
  Compare exact and interpolant at sample points.
*/
    xs = r8vec_linspace_new ( ns, xlo, xhi );

    ys = dif_vals ( nd, xd, yd, ns, xs );

    max_dif = 0.0;
    for ( i = 0; i < ns; i++ )
    {
      xt = xs[i];
      yt = 1.0 / ( 1.0 + xt * xt );
      max_dif = r8_max ( max_dif, r8_abs ( ys[i] - yt ) );
    }

    printf ( "  %4d  %14g\n", n, max_dif );

    free ( x );
    free ( xd );
    free ( xdp );
    free ( xs );
    free ( y );
    free ( yd );
    free ( ydp );
    free ( yp );
    free ( ys );
  }

  return;
}
/******************************************************************************/

void test06 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST06 tests HERMITE_BASIS_0 and HERMITE_BASIS_1.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    28 May 2011

  Author:

    John Burkardt
*/
{
# define ND 2

  double f01;
  double f02;
  double f11;
  double f12;
  int i;
  int j;
  int nd = ND;
  double xd[ND] = { 0.0, 10.0 };
  double xv;
  double yd[ND];
  double yh;
  double ypd[ND];
  double yv;

  printf ( "\n" );
  printf ( "TEST06:\n" );
  printf ( "  HERMITE_BASIS_0 and HERMITE_BASIS_1 evaluate the\n" );
  printf ( "  Hermite global polynomial basis functions\n" );
  printf ( "  of type 0: associated with function values, and\n" );
  printf ( "  of type 1: associated with derivative values.\n" );
/*
  Let y = x^3 + x^2 + x + 1,
  and compute the Hermite global polynomial interpolant based on two 
  abscissas:
*/
  for ( j = 0; j < nd; j++ )
  {
    yd[j] = pow ( xd[j], 3 ) + pow ( xd[j], 2 ) + xd[j] + 1.0;
    ypd[j] = 3.0 * pow ( xd[j], 2 ) + 2.0 * xd[j] + 1.0;
  }

  printf ( "\n" );
  printf ( "  Interpolate y = x^3 + x^2 + x + 1.\n" );
  printf ( "\n" );
  printf ( "     XD         Y(XD)      Y'(XD)\n" );
  printf ( "\n" );
  for ( j = 0; j < nd; j++ )
  {
    printf ( "  %10.4g  %10.4g  %10.4g\n", xd[j], yd[j], ypd[j] );
  }

  printf ( "\n" );
  printf ( "     XV         Y(XV)      H(XV)\n" );
  printf ( "\n" );

  for ( j = 0; j <= 10; j++ )
  {
    xv = ( double ) ( j );

    yv = pow ( xv, 3 ) + pow ( xv, 2 ) + xv + 1.0;

    f01 = hermite_basis_0 ( 2, xd, 0, xv );
    f11 = hermite_basis_1 ( 2, xd, 0, xv );
    f02 = hermite_basis_0 ( 2, xd, 1, xv );
    f12 = hermite_basis_1 ( 2, xd, 1, xv );

    yh = yd[0] * f01 + ypd[0] * f11 + yd[1] * f02 + ypd[1] * f12;

    printf ( "  %10.4g  %10.4g  %10.4g\n", xv, yv, yh );
  }
  return;
# undef ND
}
/******************************************************************************/

void test07 ( )

/******************************************************************************/
/*
  Purpose:

    TEST07 tests HERMITE_INTERPOLANT_RULE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 June 2011

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  int i;
  int k;
  int n;
  double q;
  double *w;
  double *x;

  printf ( "\n" );
  printf ( "TEST07:\n" );
  printf ( "  HERMITE_INTERPOLANT_RULE\n" );
  printf ( "  is given a set of N abscissas for a Hermite interpolant\n" );
  printf ( "  and returns N pairs of quadrature weights\n" );
  printf ( "  for function and derivative values at the abscissas.\n" );

  n = 3;
  a = 0.0;
  b = 10.0;
  x = r8vec_linspace_new ( n, a, b );
  w = hermite_interpolant_rule ( n, a, b, x );

  printf ( "\n" );
  printf ( "     I       X               W(F(X))        W(F'(X))\n" );
  printf ( "\n" );
  k = 0;
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %4d  %14g  %14g  %14g\n", i, x[i], w[k], w[k+1] );
    k = k + 2;
  }

  printf ( "\n" );
  printf ( "  Use the quadrature rule over interval %g to %g\n", a, b );
  printf ( "\n" );

  q = 0.0;
  k = 0;
  for ( i = 0; i < n; i++ )
  {
    q = q + w[k] * 1 + w[k+1] * 0.0;
    k = k + 2;
  }
  printf ( "  Estimate integral of 1 = %14g\n", q );

  q = 0.0;
  k = 0;
  for ( i = 0; i < n; i++ )
  {
    q = q + w[k] * x[i] + w[k+1] * 1.0;
    k = k + 2;
  }
  printf ( "  Estimate integral of X = %14g\n", q );

  q = 0.0;
  k = 0;
  for ( i = 0; i < n; i++ )
  {
    q = q + w[k] * x[i] * x[i] + w[k+1] * 2.0 * x[i];
    k = k + 2;
  }
  printf ( "  Estimate integral of X^2 = %14g\n", q );

  q = 0.0;
  k = 0;
  for ( i = 0; i < n; i++ )
  {
    q = q + w[k] * sin ( x[i] ) + w[k+1] * cos ( x[i] );
    k = k + 2;
  }
  printf ( "  Estimate integral of SIN(X) = %14g\n", q );

  free ( w );
  free ( x );

  n = 3;
  a = 0.0;
  b = 1.0;
  x = r8vec_linspace_new ( n, a, b );
  w = hermite_interpolant_rule ( n, a, b, x );

  printf ( "\n" );
  printf ( "     I       X               W(F(X))        W(F'(X))\n" );
  printf ( "\n" );
  k = 0;
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %4d  %14g  %14g  %14g\n", i, x[i], w[k], w[k+1] );
    k = k + 2;
  }

  printf ( "\n" );
  printf ( "  Use the quadrature rule over interval %g to %g\n", a, b );
  printf ( "\n" );

  q = 0.0;
  k = 0;
  for ( i = 0; i < n; i++ )
  {
    q = q + w[k] * 1 + w[k+1] * 0.0;
    k = k + 2;
  }
  printf ( "  Estimate integral of 1 = %14g\n", q );

  q = 0.0;
  k = 0;
  for ( i = 0; i < n; i++ )
  {
    q = q + w[k] * x[i] + w[k+1] * 1.0;
    k = k + 2;
  }
  printf ( "  Estimate integral of X = %14g\n", q );

  q = 0.0;
  k = 0;
  for ( i = 0; i < n; i++ )
  {
    q = q + w[k] * x[i] * x[i] + w[k+1] * 2.0 * x[i];
    k = k + 2;
  }
  printf ( "  Estimate integral of X^2 = %14g\n", q );

  q = 0.0;
  k = 0;
  for ( i = 0; i < n; i++ )
  {
    q = q + w[k] * sin ( x[i] ) + w[k+1] * cos ( x[i] );
    k = k + 2;
  }
  printf ( "  Estimate integral of SIN(X) = %14g\n", q );

  free ( w );
  free ( x );

  n = 11;
  a = 0.0;
  b = 10.0;
  x = r8vec_linspace_new ( n, a, b );
  w = hermite_interpolant_rule ( n, a, b, x );

  printf ( "\n" );
  printf ( "     I       X               W(F(X))        W(F'(X))\n" );
  printf ( "\n" );
  k = 0;
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %4d  %14g  %14g  %14g\n", i, x[i], w[k], w[k+1] );
    k = k + 2;
  }

  printf ( "\n" );
  printf ( "  Use the quadrature rule over interval %g to %g\n", a, b );
  printf ( "\n" );

  q = 0.0;
  k = 0;
  for ( i = 0; i < n; i++ )
  {
    q = q + w[k] * 1 + w[k+1] * 0.0;
    k = k + 2;
  }
  printf ( "  Estimate integral of 1 = %14g\n", q );

  q = 0.0;
  k = 0;
  for ( i = 0; i < n; i++ )
  {
    q = q + w[k] * x[i] + w[k+1] * 1.0;
    k = k + 2;
  }
  printf ( "  Estimate integral of X = %14g\n", q );

  q = 0.0;
  k = 0;
  for ( i = 0; i < n; i++ )
  {
    q = q + w[k] * x[i] * x[i] + w[k+1] * 2.0 * x[i];
    k = k + 2;
  }
  printf ( "  Estimate integral of X^2 = %14g\n", q );

  q = 0.0;
  k = 0;
  for ( i = 0; i < n; i++ )
  {
    q = q + w[k] * sin ( x[i] ) + w[k+1] * cos ( x[i] );
    k = k + 2;
  }
  printf ( "  Estimate integral of SIN(X) = %14g\n", q );

  free ( w );
  free ( x );

  n = 11;
  a = 0.0;
  b = 1.0;
  x = r8vec_linspace_new ( n, a, b );
  w = hermite_interpolant_rule ( n, a, b, x );

  printf ( "\n" );
  printf ( "     I       X               W(F(X))        W(F'(X))\n" );
  printf ( "\n" );
  k = 0;
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %4d  %14g  %14g  %14g\n", i, x[i], w[k], w[k+1] );
    k = k + 2;
  }

  printf ( "\n" );
  printf ( "  Use the quadrature rule over interval %g to %g\n", a, b );
  printf ( "\n" );

  q = 0.0;
  k = 0;
  for ( i = 0; i < n; i++ )
  {
    q = q + w[k] * 1 + w[k+1] * 0.0;
    k = k + 2;
  }
  printf ( "  Estimate integral of 1 = %14g\n", q );

  q = 0.0;
  k = 0;
  for ( i = 0; i < n; i++ )
  {
    q = q + w[k] * x[i] + w[k+1] * 1.0;
    k = k + 2;
  }
  printf ( "  Estimate integral of X = %14g\n", q );

  q = 0.0;
  k = 0;
  for ( i = 0; i < n; i++ )
  {
    q = q + w[k] * x[i] * x[i] + w[k+1] * 2.0 * x[i];
    k = k + 2;
  }
  printf ( "  Estimate integral of X^2 = %14g\n", q );

  q = 0.0;
  k = 0;
  for ( i = 0; i < n; i++ )
  {
    q = q + w[k] * sin ( x[i] ) + w[k+1] * cos ( x[i] );
    k = k + 2;
  }
  printf ( "  Estimate integral of SIN(X) = %14g\n", q );

  free ( w );
  free ( x );

  n = 11;
  a = 0.0;
  b = 1.0;
  x = r8vec_chebyshev_new ( n, a, b );
  w = hermite_interpolant_rule ( n, a, b, x );

  printf ( "\n" );
  printf ( "     I       X               W(F(X))        W(F'(X))\n" );
  printf ( "\n" );
  k = 0;
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %4d  %14g  %14g  %14g\n", i, x[i], w[k], w[k+1] );
    k = k + 2;
  }

  printf ( "\n" );
  printf ( "  Use the quadrature rule over interval %g to %g\n", a, b );
  printf ( "\n" );

  q = 0.0;
  k = 0;
  for ( i = 0; i < n; i++ )
  {
    q = q + w[k] * 1 + w[k+1] * 0.0;
    k = k + 2;
  }
  printf ( "  Estimate integral of 1 = %14g\n", q );

  q = 0.0;
  k = 0;
  for ( i = 0; i < n; i++ )
  {
    q = q + w[k] * x[i] + w[k+1] * 1.0;
    k = k + 2;
  }
  printf ( "  Estimate integral of X = %14g\n", q );

  q = 0.0;
  k = 0;
  for ( i = 0; i < n; i++ )
  {
    q = q + w[k] * x[i] * x[i] + w[k+1] * 2.0 * x[i];
    k = k + 2;
  }
  printf ( "  Estimate integral of X^2 = %14g\n", q );

  q = 0.0;
  k = 0;
  for ( i = 0; i < n; i++ )
  {
    q = q + w[k] * sin ( x[i] ) + w[k+1] * cos ( x[i] );
    k = k + 2;
  }
  printf ( "  Estimate integral of SIN(X) = %14g\n", q );

  free ( w );
  free ( x );

  return;
}
/******************************************************************************/

void test08 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST08 tabulates the interpolant and its derivative. 

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 November 2011

  Author:

    John Burkardt
*/
{
  int i;
  int n;
  int nd;
  int ndp;
  int ns;
  double *x;
  double *xd;
  double *xdp;
  double *xs;
  double *y;
  double *yd;
  double *ydp;
  double *yp;
  double *ys;
  double *ysp;

  printf ( "\n" );
  printf ( "TEST08\n" );
  printf ( "  HERMITE_INTERPOLANT sets up the Hermite interpolant.\n" );
  printf ( "  HERMITE_INTERPOLANT_VALUE evaluates it.\n" );
  printf ( "  Consider data for y=sin(x) at x=0,1,2,3,4.\n" );

  n = 5;
  y = ( double * ) malloc ( n * sizeof ( double ) );
  yp = ( double * ) malloc ( n * sizeof ( double ) );

  nd = 2 * n;
  xd = ( double * ) malloc ( nd * sizeof ( double ) );
  yd = ( double * ) malloc ( nd * sizeof ( double ) );

  ndp = 2 * n - 1;
  xdp = ( double * ) malloc ( ndp * sizeof ( double ) );
  ydp = ( double * ) malloc ( ndp * sizeof ( double ) );

  x = r8vec_linspace_new ( n, 0.0, 4.0 );
  for ( i = 0; i < n; i++ )
  {
    y[i] = sin ( x[i] );
    yp[i] = cos ( x[i] );
  }

  hermite_interpolant ( n, x, y, yp, xd, yd, xdp, ydp );
/*
  Now sample the interpolant at NS points, which include data values.
*/
  ns = 4 * ( n - 1 ) + 1;
  ys = ( double * ) malloc ( ns * sizeof ( double ) );
  ysp = ( double * ) malloc ( ns * sizeof ( double ) );

  xs = r8vec_linspace_new ( ns, 0.0, 4.0 );

  hermite_interpolant_value ( nd, xd, yd, xdp, ydp, ns, xs, ys, ysp );

  printf ( "\n" );
  printf ( "  In the following table, there should be perfect\n" );
  printf ( "  agreement between F and H, and F' and H'\n" );
  printf ( "  at the data points X = 0, 1, 2, 3, and 4.\n" );
  printf ( "\n" );
  printf ( "  In between, H and H' approximate F and F'.\n" );
  printf ( "\n" );
  printf ( "     I       X(I)          F(X(I))         H(X(I)) " );
  printf ( "        F'(X(I))        H'(X(I))\n" );
  printf ( "\n" );
  for ( i = 0; i < ns; i++ )
  {
    printf ( "  %4d  %14.6g  %14.6g  %14.6g  %14.6g  %14.6g\n",
      i, xs[i], sin ( xs[i] ), ys[i], cos ( xs[i] ), ysp[i] );
  }

  free ( x );
  free ( xd );
  free ( xdp );
  free ( xs );
  free ( y );
  free ( yd );
  free ( ydp );
  free ( yp );
  free ( ys );
  free ( ysp );

  return;
}
