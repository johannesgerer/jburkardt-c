# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "hermite_cubic.h"

int main ( void );
void test01 ( void );
void test02 ( void );
void test03 ( void );
void test04 ( void );
void test05 ( void );
void test06 ( void );
void test07 ( void );
void test08 ( void );
void test09 ( void );
void test10 ( void );
void test11 ( void );
void test12 ( void );
void test13 ( void );
void test14 ( void );
void test15 ( void );
double cubic_antiderivative ( double x );
double cubic_integrate ( double a, double b );
void cubic_value ( double x, double *f, double *d, double *s, double *t );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for HERMITE_CUBIC_PRB.

  Discussion:

    HERMITE_CUBIC_PRB calls the HERMITE_CUBIC tests.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 February 2012

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "HERMITE_CUBIC_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the HERMITE_CUBIC library.\n" );

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
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "HERMITE_CUBIC_PRB\n" );
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

    TEST01 tests HERMITE_CUBIC_VALUE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 February 2012

  Author:

    John Burkardt
*/
{
  double d[1];
  double d1;
  double d2;
  double f[1];
  double f1;
  double f2;
  int i;
  int j;
  int n = 1;
  double s[1];
  double t[1];
  double x[1];
  int x_interval;
  double x1;
  double x2;

  printf ( "\n" );
  printf ( "TEST01:\n" );
  printf ( "  HERMITE_CUBIC_VALUE evaluates a Hermite cubic polynomial.\n" );
  printf ( "  Try out four sets of data:\n" );
  printf ( "  (F1,D1,F2,D2) = (1,0,0,0), (0,1,0,0), (0,0,1,0), (0,0,0,1)\n" );
  printf ( "  on [0,1] and [1.0,-2.0] (interval reversed)\n" );

  for ( x_interval = 1; x_interval <= 2; x_interval++ )
  {
    if ( x_interval == 1 )
    {
      x1 = 0.0;
      x2 = 1.0;
    }
    else
    {
      x1 = 1.0;
      x2 = -2.0;
    }

    for ( i = 1; i <= 4; i++ )
    {
      f1 = 0.0;
      d1 = 0.0;
      f2 = 0.0;
      d2 = 0.0;

      if ( i == 1 )
      {
        f1 = 1.0;
      }
      else if ( i == 2 )
      {
        d1 = 1.0;
      }
      else if ( i == 3 )
      {
        f2 = 1.0;
      }
      else if ( i == 4 )
      {
        d2 = 1.0;
      }

      printf ( "\n" );
      printf ( "    J      X           F           D\n" );
      printf ( "\n" );

      for ( j = -3; j <= 12; j++ )
      {
        x[0] = ( ( double ) ( 10 - j ) * x1
               + ( double ) (      j ) * x2 )
               / ( double ) ( 10     );

        hermite_cubic_value ( x1, f1, d1, x2, f2, d2, n, x, f, d, s, t );

        if ( j == 0 )
        {
          printf ( "*Data  %10g  %10g  %10g\n", x1, f1, d1 );
        }
        printf ( "  %3d  %10g  %10g  %10g\n", j, x[0], f[0], d[0] );
        if ( j == 10 )
        {
          printf ( "*Data  %10g  %10g  %10g\n", x2, f2, d2 );
        }
      }
    }
  }
  return;
}
/******************************************************************************/

void test02 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST02 tests HERMITE_CUBIC_VALUE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 February 2012

  Author:

    John Burkardt
*/
{
  double d[1];
  double dc;
  double d1;
  double d2;
  double f[1];
  double fc;
  double f1;
  double f2;
  int j;
  int n = 1;
  double s[1];
  double s1;
  double s2;
  double sc;
  double t[1];
  double t1;
  double t2;
  double tc;
  double x[1];
  int x_interval;
  double x1;
  double x2;

  printf ( "\n" );
  printf ( "TEST02:\n" );
  printf ( "  HERMITE_CUBIC_VALUE evaluates a Hermite cubic polynomial.\n" );
  printf ( "  Try out data from a cubic function:\n" );
  printf ( "  on [0,10] and [-1.0,1.0] and [0.5,0.75]\n" );

  for ( x_interval = 1; x_interval <= 3; x_interval++ )
  {
    if ( x_interval == 1 )
    {
      x1 = 0.0;
      x2 = 10.0;
    }
    else if ( x_interval == 2 )
    {
      x1 = -1.0;
      x2 = +1.0;
    }
    else if ( x_interval == 3 )
    {
      x1 = 0.5;
      x2 = 0.75;
    }

    cubic_value ( x1, &f1, &d1, &s1, &t1 );
    cubic_value ( x2, &f2, &d2, &s2, &t2 );

    printf ( "\n" );
    printf ( "    J      X           F           D           S           T\n" );
    printf ( "\n" );

    for ( j = -3; j <= 12; j++ )
    {
      x[0] = ( ( double ) ( 10 - j ) * x1
             + ( double ) (      j ) * x2 )
             / ( double ) ( 10 );

      hermite_cubic_value ( x1, f1, d1, x2, f2, d2, n, x, f, d, s, t );
      cubic_value ( x[0], &fc, &dc, &sc, &tc );

      if ( j == 0 )
      {
        printf ( "*Data  %10g  %10g  %10g\n", x1, f1, d1 );
      }
      printf ( "Exact  %10g  %10g  %10g  %10g  %10g\n",
        x[0], fc, dc, sc, tc );
      printf ( "  %3d  %10g  %10g  %10g  %10g  %10g\n",
        j, x[0], f[0], d[0], s[0], t[0] );
      if ( j == 10 )
      {
        printf ( "*Data  %10g  %10g  %10g\n", x2, f2, d2 );
      }
    }
  }
  return;
}
/******************************************************************************/

void test03 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST03 tests HERMITE_CUBIC_INTEGRATE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 February 2012

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double d1;
  double d2;
  double f1;
  double f2;
  int j;
  double q_computed;
  double q_exact;
  double s1;
  double s2;
  double t1;
  double t2;
  int x_interval;
  double x1;
  double x2;

  printf ( "\n" );
  printf ( "TEST03:\n" );
  printf ( "  HERMITE_CUBIC_INTEGRATE integrates a Hermite cubic\n" );
  printf ( "  polynomial from A to B.\n" );

  for ( x_interval = 1; x_interval <= 3; x_interval++ )
  {
    if ( x_interval == 1 )
    {
      x1 = 0.0;
      x2 = 10.0;
    }
    else if ( x_interval == 2 )
    {
      x1 = -1.0;
      x2 = +1.0;
    }
    else if ( x_interval == 3 )
    {
      x1 = 0.5;
      x2 = 0.75;
    }

    cubic_value ( x1, &f1, &d1, &s1, &t1 );
    cubic_value ( x2, &f2, &d2, &s2, &t2 );

    printf ( "\n" );
    printf ( "                                     Exact           Computed\n" );
    printf ( "    J          A           B         Integral        Integral\n" );
    printf ( "\n" );

    a = x1 - 1.0;

    for ( j = - 3; j <= 12; j++ )
    {
      b = ( ( double ) ( 10 - j ) * x1
          + ( double ) (      j ) * x2 )
          / ( double ) ( 10     );

      q_exact = cubic_integrate ( a, b );

      q_computed = hermite_cubic_integrate ( x1, f1, d1, x2, f2, d2, a, b );

      printf ( "  %3d    %10g  %10g  %14g  %14g\n",
        j, a, b, q_exact, q_computed );
    }
  }
  return;
}
/******************************************************************************/

void test04 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST04 tests HERMITE_CUBIC_SPLINE_VALUE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 February 2012

  Author:

    John Burkardt
*/
{
  double *d;
  double *dn;
  double *f;
  double *fn;
  int i;
  int n = 51;
  int nn = 11;
  double *s;
  double *t;
  double u;
  double v;
  double x1;
  double x2;
  double *x;
  double *xn;

  printf ( "\n" );
  printf ( "TEST04:\n" );
  printf ( "  HERMITE_CUBIC_SPLINE_VALUE evaluates a Hermite cubic spline.\n" );

  x1 = 0.0;
  x2 = 10.0;

  xn = r8vec_even_new ( nn, x1, x2 );
  fn = ( double * ) malloc ( nn * sizeof ( double ) );
  dn = ( double * ) malloc ( nn * sizeof ( double ) );

  for ( i = 0; i < nn; i++ )
  {
    fn[i] = sin ( xn[i] );
    dn[i] = cos ( xn[i] );
  }

  x = r8vec_even_new ( n, x1, x2 );
  f = ( double * ) malloc ( n * sizeof ( double ) );
  d = ( double * ) malloc ( n * sizeof ( double ) );
  s = ( double * ) malloc ( n * sizeof ( double ) );
  t = ( double * ) malloc ( n * sizeof ( double ) );

  hermite_cubic_spline_value ( nn, xn, fn, dn, n, x, f, d, s, t );

  printf ( "\n" );
  printf ( "     I      X       F computed     F exact      Error\n" );
  printf ( "\n" );

  for ( i = 0; i < n; i++ )
  {
    u = sin ( x[i] );
    v = r8_abs ( f[i] - u );
    printf ( "  %4d    %10g  %10g  %10g  %14g\n", i, x[i], f[i], u, v );
  }

  printf ( "\n" );
  printf ( "     I      X       D computed     D exact      Error\n" );
  printf ( "\n" );

  for ( i = 0; i < n; i++ )
  {
    u = cos ( x[i] );
    v = r8_abs ( d[i] - u );
    printf ( "  %4d  %10g  %10g  %10g  %14g\n", i, x[i], d[i], u, v );
  }

  free ( d );
  free ( dn );
  free ( f );
  free ( fn );
  free ( s );
  free ( t );
  free ( x );
  free ( xn );

  return;
}
/******************************************************************************/

void test05 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST05 tests HERMITE_CUBIC_TO_POWER_CUBIC

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 February 2012

  Author:

    John Burkardt
*/
{
  double c0;
  double c1;
  double c2;
  double c3;
  double d[1];
  double d1;
  double d1r;
  double d2;
  double d2r;
  double f[1];
  double f1;
  double f1r;
  double f2;
  double f2r;
  double fp;
  int j;
  int n = 1;
  double s[1];
  double s1;
  double s2;
  double t[1];
  double t1;
  double t2;
  double x[1];
  double x1;
  double x2;

  printf ( "\n" );
  printf ( "TEST05:\n" );
  printf ( "  HERMITE_CUBIC_TO_POWER_CUBIC converts the Hermite data\n" );
  printf ( "  to the coefficients of the power form of the polynomial\n" );
  printf ( "  POWER_CUBIC_TO_HERMITE_CUBIC converts the power form\n" );
  printf ( "  to Hermite form\n" );

  x1 = -1.0;
  x2 = +1.0;

  cubic_value ( x1, &f1, &d1, &s1, &t1 );
  cubic_value ( x2, &f2, &d2, &s2, &t2 );

  printf ( "\n" );
  printf ( "  Hermite data:\n" );
  printf ( "\n" );
  printf ( "  X1, F1, D1:  %10g  %10g  %10g\n", x1, f1, d1 );
  printf ( "  X2, F2, D2:  %10g  %10g  %10g\n", x2, f2, d2 );

  hermite_cubic_to_power_cubic ( x1, f1, d1, x2, f2, d2, &c0, &c1, &c2, &c3 );

  printf ( "\n" );
  printf ( "  Power form:\n" );
  printf ( "  p(x) = %g + %g * x + %g * x^2 + %g * x^3\n", c0, c1, c2, c3 );
  printf ( "\n" );
  printf ( "      X       F (Hermite)  F (power)\n" );
  printf ( "\n" );

  for ( j = -3; j <= 12; j++ )
  {
    x[0] = ( ( double ) ( 10 - j ) * x1
           + ( double ) (      j ) * x2 )
           / ( double ) ( 10     );

    hermite_cubic_value ( x1, f1, d1, x2, f2, d2, n, x, f+0, d+0, s+0, t+0 );

    fp = c0 + x[0] * ( c1 + x[0] * ( c2 + x[0] * c3 ) );

    printf ( "    %10g  %10g  %10g\n", x[0], f[0], fp );
  }

  power_cubic_to_hermite_cubic ( c0, c1, c2, c3, x1, x2, &f1r, &d1r,
    &f2r, &d2r );

  printf ( "\n" );
  printf ( "  Use POWER_CUBIC_TO_HERMITE_CUBIC to recover the\n" );
  printf ( "  original Hermite data:\n" );
  printf ( "\n" );
  printf ( "         Original   Recovered\n" );
  printf ( "\n" );
  printf ( "  F1:    %10g  %10g\n", f1, f1r );
  printf ( "  D1:    %10g  %10g\n", d1, d1r );
  printf ( "  F2:    %10g  %10g\n", f2, f2r );
  printf ( "  D2:    %10g  %10g\n", d2, d2r );
  return;
}
/******************************************************************************/

void test06 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST06 tests HERMITE_CUBIC_INTEGRATE using vectors.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 February 2012

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double b_hi;
  double b_lo;
  double d1;
  double d2;
  double f1;
  double f2;
  int i;
  double q_computed;
  double q_exact;
  double s1;
  double s2;
  double t1;
  double t2;
  double x1;
  double x2;

  printf ( "\n" );
  printf ( "TEST06:\n" );
  printf ( "  HERMITE_CUBIC_INTEGRATE integrates a Hermite cubic\n" );
  printf ( "  polynomial from A to B.\n" );
  printf ( "  Use A, B vectors for the calculation.\n" );

  x1 = 0.0;
  x2 = 10.0;

  cubic_value ( x1, &f1, &d1, &s1, &t1 );
  cubic_value ( x2, &f2, &d2, &s2, &t2 );

  printf ( "\n" );
  printf ( "                                 Exact       Computed\n" );
  printf ( "    J      A           B         Integral    Integral\n" );
  printf ( "\n" );

  for ( i = -3; i <= 12; i++ )
  {
    a = x1 - 1.0;
    b = ( ( double ) ( 10 - i ) * x1
        + ( double ) (      i ) * x2 )
        / ( double ) ( 10     );

    q_exact = cubic_integrate ( a, b );

    q_computed = hermite_cubic_integrate ( x1, f1, d1, x2, f2, d2, a, b );

    printf ( "  %3d  %10g  %10g  %14g  %14g\n", i, a, b, q_exact, q_computed );
  }
  return;
}
/******************************************************************************/

void test07 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST07 tests HERMITE_CUBIC_INTEGRAL.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 February 2012

  Author:

    John Burkardt
*/
{
  double d1;
  double d2;
  double f1;
  double f2;
  double q_computed;
  double q_exact;
  double s1;
  double s2;
  double t1;
  double t2;
  int x_interval;
  double x1;
  double x2;

  printf ( "\n" );
  printf ( "TEST07:\n" );
  printf ( "  HERMITE_CUBIC_INTEGRAL integrates a Hermite cubic\n" );
  printf ( "  polynomial over the definition interval [X1,X2].\n" );
  printf ( "\n" );
  printf ( "                            Exact       Computed\n" );
  printf ( "     X1          X2         Integral    Integral\n" );
  printf ( "\n" );

  for ( x_interval = 1; x_interval <= 3; x_interval++ )
  {
    if ( x_interval == 1 )
    {
      x1 = 0.0;
      x2 = 10.0;
    }
    else if ( x_interval == 2 )
    {
      x1 = -1.0;
      x2 = +1.0;
    }
    else if ( x_interval == 3 )
    {
      x1 = 0.5;
      x2 = 0.75;
    }

    cubic_value ( x1, &f1, &d1, &s1, &t1 );
    cubic_value ( x2, &f2, &d2, &s2, &t2 );

    q_exact = cubic_integrate ( x1, x2 );

    q_computed = hermite_cubic_integral ( x1, f1, d1, x2, f2, d2 );

    printf ( "    %10g  %10g  %14g  %14g\n", x1, x2, q_exact, q_computed );
  }
  return;
}
/******************************************************************************/

void test08 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST08 tests HERMITE_CUBIC_SPLINE_INTEGRAL.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 February 2012

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double *dn;
  double *fn;
  int i;
  int nn = 11;
  double pi = 3.141592653589793;
  double q_computed;
  double q_exact;
  int test;
  double *xn;

  printf ( "\n" );
  printf ( "TEST08:\n" );
  printf ( "  HERMITE_CUBIC_SPLINE_INTEGRAL integrates a Hermite\n" );
  printf ( "  cubic spline over the definition interval [X1,XNN].\n" );
  printf ( "\n" );
  printf ( "                            Exact       Computed\n" );
  printf ( "     X1          XNN        Integral    Integral\n" );
  printf ( "\n" );

  fn = ( double * ) malloc ( nn * sizeof ( double ) );
  dn = ( double * ) malloc ( nn * sizeof ( double ) );

  for ( test = 1; test <= 3; test++ )
  {
    if ( test == 1 )
    {
      a = 0.0;
      b = 1.0;

      xn = r8vec_even_new ( nn, a, b );
      for ( i = 0; i < nn; i++ )
      {
        fn[i] = xn[i] * ( 4.0 * xn[i] - 1.0 ) * ( xn[i] - 1.0 );
        dn[i] = 1.0 + xn[i] * ( - 10.0 + xn[i] * 12.0 );
      }
      q_exact =
        ( xn[nn-1] * xn[nn-1] * ( 0.5 + xn[nn-1] * ( - ( 5.0 / 3.0 ) + xn[nn-1] ) ) )
      - ( xn[0]  * xn[0]  * ( 0.5 + xn[0]  * ( - ( 5.0 / 3.0 ) + xn[0]  ) ) );
    }
/*
  Use variable spacing.
*/
    else if ( test == 2 )
    {
      a = 0.0;
      b = 1.0;

      xn = r8vec_even_new ( nn, a, b );
      for ( i = 0; i < nn; i++ )
      {
        xn[i] = sqrt ( xn[i] );
        fn[i] = xn[i] * ( 4.0 * xn[i] - 1.0 ) * ( xn[i] - 1.0 );
        dn[i] = 1.0 + xn[i] * ( - 10.0 + xn[i] * 12.0 );
      }
      q_exact =
        ( xn[nn-1] * xn[nn-1] * ( 0.5 + xn[nn-1] * ( - ( 5.0 / 3.0 ) + xn[nn-1] ) ) )
      - ( xn[0]  * xn[0]  * ( 0.5 + xn[0]  * ( - ( 5.0 / 3.0 ) + xn[0]  ) ) );
    }
/*
  Try a non-cubic.
*/
    else if ( test == 3 )
    {
      a = 0.0;
      b = pi;

      xn = r8vec_even_new ( nn, a, b );
      for ( i = 0; i < nn; i++ )
      {
        fn[i] = sin ( xn[i] );
        dn[i] = cos ( xn[i] );
      }
      q_exact = - cos ( xn[nn-1] ) + cos ( xn[0] );
    }

    q_computed = hermite_cubic_spline_integral ( nn, xn, fn, dn );

    printf ( "  %10g  %10g  %14g  %14g\n", xn[0], xn[nn-1], q_exact, q_computed );

    free ( xn );
  }
  free ( dn );
  free ( fn );

  return;
}
/******************************************************************************/

void test09 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST09 tests HERMITE_CUBIC_SPLINE_INTEGRATE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 February 2012

  Author:

    John Burkardt
*/
{
  double *a;
  double *b;
  double *dn;
  double *fn;
  int i;
  int n = 25;
  int nn = 11;
  double *q;
  double q_exact;
  double *sn;
  double *tn;
  double x1;
  double x2;
  double *xn;

  printf ( "\n" );
  printf ( "TEST09:\n" );
  printf ( "  HERMITE_CUBIC_SPLINE_INTEGRATE integrates a Hermite\n" );
  printf ( "  cubic spline from A to B.\n" );
/*
  Define the cubic spline.
*/
  x1 = 0.0;
  x2 = 10.0;

  xn = r8vec_even_new ( nn, x1, x2 );
  fn = ( double * ) malloc ( nn * sizeof ( double ) );
  dn = ( double * ) malloc ( nn * sizeof ( double ) );
  sn = ( double * ) malloc ( nn * sizeof ( double ) );
  tn = ( double * ) malloc ( nn * sizeof ( double ) );

  for ( i = 0; i < nn; i++ )
  {
    cubic_value ( xn[i], fn+i, dn+i, sn+i, tn+i );
  }

  a = ( double * ) malloc ( n * sizeof ( double ) );
  for ( i = 0; i < n; i++ )
  {
    a[i] = 2.5;
  }
  b = r8vec_even_new ( n, x1 - 1.0, x2 + 1.0 );

  q = hermite_cubic_spline_integrate ( nn, xn, fn, dn, n, a, b );

  printf ( "\n" );
  printf ( "                                 Exact       Computed\n" );
  printf ( "    I      A           B         Integral    Integral\n" );
  printf ( "\n" );

  for  ( i = 0; i < n; i++ )
  {
    q_exact = cubic_integrate ( a[i], b[i] );

    printf ( "  %3d    %10g  %10g  %10g  %10g\n", i, a[i], b[i], q_exact, q[i] );
  }

  free ( a );
  free ( b );
  free ( dn );
  free ( fn );
  free ( q );
  free ( sn );
  free ( tn );
  free ( xn );

  return;
}
/******************************************************************************/

void test10 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST10 tests HERMITE_CUBIC_SPLINE_INTEGRAL.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 February 2012

  Author:

    John Burkardt
*/
{
  char comment[100];
  double *dn;
  double *fn;
  int i;
  int nn = 11;
  double pi = 3.141592653589793;
  double q_computed;
  double q_exact;
  int seed;
  int test;
  double *xn;

  seed = 123456789;

  printf ( "\n" );
  printf ( "TEST10:\n" );
  printf ( "  HERMITE_CUBIC_SPLINE_INTEGRAL integrates a Hermite\n" );
  printf ( "  cubic spline over the definition interval [X1,XNN].\n" );
  printf ( "\n" );
  printf ( "  If the subintervals are equally spaced, the derivative\n" );
  printf ( "  information has no effect on the result, except for\n" );
  printf ( "  the first and last values, DN(1) and DN(NN).\n" );
  printf ( "\n" );
  printf ( "                            Exact       Computed\n" );
  printf ( "     X1          XNN        Integral    Integral  Comment\n" );
  printf ( "\n" );

  fn = ( double * ) malloc ( nn * sizeof ( double ) );
  dn = ( double * ) malloc ( nn * sizeof ( double ) );

  for ( test = 1; test <= 5; test++ )
  {
/*
  Equal spacing.
*/
    if ( test == 1 )
    {
      xn = r8vec_even_new ( nn, 0.0, pi );
      for ( i = 0; i < nn; i++ )
      {
        fn[i] = sin ( xn[i] );
        dn[i] = cos ( xn[i] );
      }
      q_exact = - cos ( xn[nn-1] ) + cos ( xn[0] );
      strcpy ( comment, "Equal spacing, correct DN" );
    }
/*
  Equal spacing, reset DN(2:NN-1) to random numbers.
*/
    else if ( test == 2 )
    {
      xn = r8vec_even_new ( nn, 0.0, pi );
      for ( i = 0; i < nn; i++ )
      {
        fn[i] = sin ( xn[i] );
        if ( i == 0 || i == nn - 1 )
        {
          dn[i] = cos ( xn[i] );
        }
        else
        {
          dn[i] = 1000.0 * r8_uniform_01 ( &seed );
        }
      }
      q_exact = - cos ( xn[nn-1] ) + cos ( xn[0] );
      strcpy ( comment, "Equal spacing, DN(2:N-1) random" );
    }
/*
  Equal spacing, now reset all of DN to random numbers.
*/
    else if ( test == 3 )
    {
      xn = r8vec_even_new ( nn, 0.0, pi );
      for ( i = 0; i < nn; i++ )
      {
        fn[i] = sin ( xn[i] );
        dn[i] = 1000.0 * r8_uniform_01 ( &seed );
      }
      q_exact = - cos ( xn[nn-1] ) + cos ( xn[0] );
      strcpy ( comment, "Equal spacing, DN(1:N) random" );
    }
/*
  Variable spacing, correct data.
*/
    else if ( test == 4 )
    {
      xn = r8vec_even_new ( nn, 0.0, pi * pi );
      for ( i = 0; i < nn; i++ )
      {
        xn[i] = sqrt ( xn[i] );
        fn[i] = sin ( xn[i] );
        dn[i] = cos ( xn[i] );
      }
      q_exact = - cos ( xn[nn-1] ) + cos ( xn[0] );
      strcpy ( comment, "Variable spacing, correct DN" );
    }
/*
  Variable spacing, change one entry in DN.
*/
    else if ( test == 5 )
    {
      xn = r8vec_even_new ( nn, 0.0, pi * pi );
      for ( i = 0; i < nn; i++ )
      {
        xn[i] = sqrt ( xn[i] );
        fn[i] = sin ( xn[i] );
        dn[i] = cos ( xn[i] );
      }
      dn[ ( nn - 1 ) / 2 ] = 1000.0 * r8_uniform_01 ( &seed );
      q_exact = - cos ( xn[nn-1] ) + cos ( xn[0] );
      strcpy ( comment, "Variable spacing, a single internal DN randomized." );
    }

    q_computed = hermite_cubic_spline_integral ( nn, xn, fn, dn );

    printf ( "  %10g  %10g  %14g  %14g  %s\n",
      xn[0], xn[nn-1], q_exact, q_computed, comment );

    free ( xn );
  }

  free ( dn );
  free ( fn );

  return;
}
/******************************************************************************/

void test11 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST11 tests HERMITE_CUBIC_LAGRANGE_VALUE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 February 2012

  Author:

    John Burkardt
*/
{
  double *d;
  double *f;
  int j;
  int n = 11;
  double *s;
  double *t;
  double *x;
  double x1;
  double x2;

  printf ( "\n" );
  printf ( "TEST11:\n" );
  printf ( "  HERMITE_CUBIC_LAGRANGE_VALUE evaluates the four\n" );
  printf ( "  Lagrange basis functions associated with F1, D1,\n" );
  printf ( "  F2 and D2 such that\n" );
  printf ( "\n" );
  printf ( "  P(X) = F1 * LF1(X) + D1 * LD1(X)\n" );
  printf ( "       + F2 * LF2(X) + D2 * LD2(X).\n" );
  printf ( "\n" );
  printf ( "  The first, second and third derivatives of these four\n" );
  printf ( "  Lagrange basis functions are also computed.\n" );

  x1 = 1.0;
  x2 = 2.0;
  x = r8vec_even_new ( n, 0.0, 2.5 );

  f = ( double * ) malloc ( 4 * n * sizeof ( double ) );
  d = ( double * ) malloc ( 4 * n * sizeof ( double ) );
  s = ( double * ) malloc ( 4 * n * sizeof ( double ) );
  t = ( double * ) malloc ( 4 * n * sizeof ( double ) );

  hermite_cubic_lagrange_value ( x1, x2, n, x, f, d, s, t );

  printf ( "\n" );
  printf ( "  The Lagrange basis functions:\n" );
  printf ( "\n" );
  printf ( "     I        X           LF1         LD1         LF2         LD2\n" );
  printf ( "\n" );
  for ( j = 0; j < n; j++ )
  {
    printf ( "  %4d  %10g  %10g  %10g  %10g  %10g\n", 
      j, x[j], f[0+j*4], f[1+j*4], f[2+j*4], f[3+j*4] );
  }

  printf ( "\n" );
  printf ( "  The derivative of the Lagrange basis functions:\n" );
  printf ( "\n" );
  printf ( "     I        X           LF1         LD1         LF2         LD2\n" );
  printf ( "\n" );
  for ( j = 0; j < n; j++ )
  {
    printf ( "  %4d  %10g  %10g  %10g  %10g  %10g\n",
      j, x[j], d[0+j*4], d[1+j*4], d[2+j*4], d[3+j*4] );
  }

  free ( d );
  free ( f );
  free ( s );
  free ( t );
  free ( x );

  return;
}
/******************************************************************************/

void test12 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST12 tests HERMITE_CUBIC_LAGRANGE_INTEGRAL.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 February 2012

  Author:

    John Burkardt
*/
{
  int i;
  double *q;
  double x1;
  double x2;

  printf ( "\n" );
  printf ( "TEST12:\n" );
  printf ( "  HERMITE_CUBIC_LAGRANGE_INTEGRAL returns the integrals\n" );
  printf ( "  of the four Lagrange basis functions associated\n" );
  printf ( "  with F1, D1, F2 and D2 such that\n" );
  printf ( "\n" );
  printf ( "  P(X) = F1 * LF1(X) + D1 * LD1(X)\n" );
  printf ( "       + F2 * LF2(X) + D2 * LD2(X).\n" );
  printf ( "\n" );
  printf ( "  The Lagrange basis function integrals:\n" );
  printf ( "\n" );
  printf ( "        X1          X2          LF1         LD1         LF2         LD2\n" );
  printf ( "\n" );

  x2 = 1.0;

  for ( i = -6; i <= 2; i++ )
  {
    x1 = ( double ) ( i );
    q = hermite_cubic_lagrange_integral ( x1, x2 );
    printf ( "  %10g  %10g  %10g  %10g  %10g  %10g\n",
      x1, x2, q[0], q[1], q[2], q[3] );
    free ( q );
  }
  return;
}
/******************************************************************************/

void test13 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST13 tests HERMITE_CUBIC_LAGRANGE_INTEGRATE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 February 2012

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double d1;
  double d2;
  double f1;
  double f2;
  int j;
  double p[4];
  double *q;
  double x1;
  double x2;

  printf ( "\n" );
  printf ( "TEST13:\n" );
  printf ( "  HERMITE_CUBIC_LAGRANGE_INTEGRATE integrates a Hermite cubic\n" );
  printf ( "  Lagrange polynomial from A to B.\n" );
  printf ( "\n" );
  printf ( "  Compute each result TWICE:\n" );
  printf ( "  First row computed using HERMITE_CUBIC_INTEGRATE.\n" );
  printf ( "  Second row computed using HERMITE_CUBIC_LAGRANGE_INTEGRATE.\n" );

  x1 = 0.0;
  x2 = 10.0;

  printf ( "\n" );
  printf ( "        A           B           LF1         LD1         LF2         LD2\n" );
  printf ( "\n" );

  a = x1 - 1.0;

  for ( j = -3; j <= 12; j++ )
  {
    b = ( ( double ) ( 10 - j ) * x1
        + ( double ) (      j ) * x2 )
        / ( double ) ( 10     );

    f1 = 1.0;
    d1 = 0.0;
    f2 = 0.0;
    d2 = 0.0;
    p[0] = hermite_cubic_integrate ( x1, f1, d1, x2, f2, d2, a, b );

    f1 = 0.0;
    d1 = 1.0;
    f2 = 0.0;
    d2 = 0.0;
    p[1] = hermite_cubic_integrate ( x1, f1, d1, x2, f2, d2, a, b );

    f1 = 0.0;
    d1 = 0.0;
    f2 = 1.0;
    d2 = 0.0;
    p[2] = hermite_cubic_integrate ( x1, f1, d1, x2, f2, d2, a, b );

    f1 = 0.0;
    d1 = 0.0;
    f2 = 0.0;
    d2 = 1.0;
    p[3] = hermite_cubic_integrate ( x1, f1, d1, x2, f2, d2, a, b );

    q = hermite_cubic_lagrange_integrate ( x1, x2, a, b );

    printf ( "  %10g  %10g  %10g  %10g  %10g  %10g\n",
      a, b, p[0], p[1], p[2], p[3] );
    printf ( "                          %10g  %10g  %10g  %10g\n",
      q[0], q[1], q[2], q[3] );

    free ( q );
  }

  return;
}
/******************************************************************************/

void test14 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST14 tests HERMITE_CUBIC_SPLINE_QUAD_RULE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 February 2012

  Author:

    John Burkardt
*/
{
# define N 11

  double dn[N];
  double fn[N];
  int i;
  int j;
  int k;
  int l;
  int n = N;
  double q;
  double *r;
  int seed;
  double *w;
  double x[N];

  printf ( "\n" );
  printf ( "TEST14:\n" );
  printf ( "  HERMITE_CUBIC_SPLINE_QUAD_RULE returns a quadrature rule\n" );
  printf ( "  for Hermite cubic splines.\n" );

  seed = 123456789;

  for ( k = 1; k <= 2; k++ )
  {
    printf ( "\n" );
    if ( k == 1 )
    {
      printf ( "  Case 1: Random spacing\n" );
      r = r8vec_uniform_01_new ( n, &seed );

      x[0] = r[0];
      for ( i = 1; i < n; i++ )
      {
        x[i] = x[i-1] + r[i];
      }
      free ( r );
    }
    else if ( k == 2 )
    {
      printf ( "  Case 2: Uniform spacing\n" );
      printf ( "  F(2:N-1) have equal weight.\n" );
      printf ( "  D(2:N-1) have zero weight.\n" );
      for ( i = 0; i < n; i++ )
      {
        x[i] = ( double ) ( 10 + i ) / 20.0;
      }
    }

    w = hermite_cubic_spline_quad_rule ( n, x );

    printf ( "\n" );
    printf ( "   I   J        X         W                Q\n" );
    printf ( "\n" );

    for ( i = 0; i <= 1; i++ )
    {
      for ( j = 0; j < n; j++ )
      {
        for ( l = 0; l < n; l++ )
        {
          fn[l] = 0.0;
          dn[l] = 0.0;
        }
        if ( i == 0 )
        {
          fn[j] = 1.0;
        }
        else
        {
          dn[j] = 1.0;
        }

        q = hermite_cubic_spline_integral ( n, x, fn, dn );

        printf ( "  %2d  %2d    %10g  %14g  %14g\n",
          i, j, x[j], q, w[i+j*2]  );
      }
    }
    free ( w );
  }
  return;
# undef N
}
/******************************************************************************/

void test15 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST15 tests HERMITE_CUBIC_SPLINE_QUAD_RULE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 February 2012

  Author:

    John Burkardt
*/
{
# define N 11

  double dn;
  double fn;
  int j;
  int n = N;
  double q;
  double q_exact;
  double *r;
  double s;
  int seed;
  double t;
  double *w;
  double x[N];

  printf ( "\n" );
  printf ( "TEST15:\n" );
  printf ( "  HERMITE_CUBIC_SPLINE_QUAD_RULE returns a quadrature rule\n" );
  printf ( "  for Hermite cubic splines.\n" );

  seed = 123456789;

  r = r8vec_uniform_01_new ( n, &seed );

  x[0] = r[0];
  for ( j = 1; j < n; j++ )
  {
    x[j] = x[j-1] + r[j];
  }
  free ( r );

  printf ( "\n" );
  printf ( "  Random spacing\n" );
  printf ( "  Number of points N = %d\n", n );
  printf ( "  Interval = [%g, %g]\n", x[0], x[n-1] );

  w = hermite_cubic_spline_quad_rule ( n, x );

  q = 0.0;

  for ( j = 0; j < n; j++ )
  {
    cubic_value ( x[j], &fn, &dn, &s, &t );
    q = q + w[0+j*2] * fn + w[1+j*2] * dn;
  }

  q_exact = cubic_integrate ( x[0], x[n-1] );

  printf ( "\n" );
  printf ( "  Q         = %g\n", q );
  printf ( "  Q (exact) = %g\n", q_exact );

  free ( w );

  return;
# undef N
}
/******************************************************************************/

double cubic_antiderivative ( double x )

/******************************************************************************/
/*
  Purpose:

    CUBIC_ANTIDERIVATIVE evaluates the antiderivative function of a cubic.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 February 2012

  Author:

    John Burkardt

  Parameters:

    Input, double X, the argument.

    Output, double CUBIC_ANTIDERIVATIVE, the value.
*/
{
  double value;

  value = x * x * ( 5.0 + x * ( - 7.0 / 3.0 + x * 1.0 / 4.0 ) );

  return value;
}
/******************************************************************************/

double cubic_integrate ( double a, double b )

/******************************************************************************/
/*
  Purpose:

    CUBIC_INTEGRATE integrates the cubic from A to B.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 February 2012

  Author:

    John Burkardt

  Parameters:

    Input, double A, B, the integration interval.

    Output, double Q, the integral from A to B.
*/
{
  double q;

  q = cubic_antiderivative ( b ) - cubic_antiderivative ( a );

  return q;
}
/******************************************************************************/

void cubic_value ( double x, double *f, double *d, double *s, double *t )

/******************************************************************************/
/*
  Purpose:

    CUBIC_VALUE evaluates a cubic function.

  Discussion:

    f(x) =   x^3 -  7 x^2 + 10 x
    d(x) = 3 x^2 - 14 x   + 10
    s(x) = 6 x   - 14
    t(x) = 6

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 February 2012

  Author:

    John Burkardt

  Parameters:

    Input, double X, the argument.

    Output, double F, D, S, T, the value and first three
    derivatives of the cubic function.
*/
{
  *f = 0.0 + x * ( 10.0 + x * ( -  7.0 + x * 1.0 ) );
  *d =             10.0 + x * ( - 14.0 + x * 3.0 );
  *s =                          - 14.0 + x * 6.0;
  *t =                                       6.0;

  return;
}
