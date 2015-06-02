# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "piecewise_linear_product_integral.h"

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );
void test04 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for PIECEWISE_LINEAR_PRODUCT_INTEGRAL_PRB.

  Discussion:

    PIECEWISE_LINEAR_PRODUCT_INTEGRAL_PRB tests the 
    PIECEWISE_LINEAR_PRODUCT_INTEGRAL library.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    04 July 2013

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "PIECEWISE_LINEAR_PRODUCT_INTEGRAL_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the PIECEWISE_LINEAR_PRODUCT_INTEGRAL_INTEGRAL library.\n" );
 
  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "PIECEWISE_LINEAR_PRODUCT_INTEGRAL_PRB\n" );
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

    TEST01 tests PIECEWISE_LINEAR_PRODUCT_INTEGRAL.

  Discussion:

    For the first test, we use the same single "piece" for both F and G.
    Hence, we are actually integrating X^2.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    04 July 2013

  Author:

    John Burkardt
*/
{
# define F_NUM 2
# define G_NUM 2

  double a;
  double b;
  double exact;
  int f_num = F_NUM;
  double f_v[F_NUM] = { 0.0, 5.0 };
  double f_x[F_NUM] = { 0.0, 5.0 };
  int g_num = G_NUM;
  double g_v[G_NUM] = { 0.0, 5.0 };
  double g_x[G_NUM] = { 0.0, 5.0 };
  int i;
  double integral;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  Test PIECEWISE_LINEAR_PRODUCT_INTEGRAL on a very simple problem.\n" );
  printf ( "  F and G are both defined over a single common\n" );
  printf ( "  interval, so that F(X) = G(X) = X.\n" );
  printf ( "\n" );
  printf ( "           A           B      Integral        Exact\n" );
  printf ( "\n" );

  a = 1.0;
  for ( i = 1; i <= 5; i++ )
  {
    b = ( double ) ( i );
    integral = piecewise_linear_product_integral ( a, b, f_num, f_x, f_v, g_num, 
      g_x, g_v );
    exact = ( b * b * b - a * a * a ) / 3.0;
    printf ( "  %10g  %10g  %14g  %14g\n", a, b, integral, exact );
  }

  return;
# undef F_NUM
# undef G_NUM
}
/******************************************************************************/

void test02 ( )

/******************************************************************************/
/*
  Purpose:

    TEST02 tests PIECEWISE_LINEAR_PRODUCT_INTEGRAL.

  Discussion:

    For this test, we use multiple "pieces" for both F and G,
    but we define the values so that we are still actually integrating X^2.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    04 July 2013

  Author:

    John Burkardt
*/
{
# define F_NUM 3
# define G_NUM 4

  double a;
  double b;
  double exact;
  int f_num = F_NUM;
  double f_v[F_NUM] = { 0.0, 2.0, 5.0 };
  double f_x[F_NUM] = { 0.0, 2.0, 5.0 };
  int g_num = G_NUM;
  double g_v[G_NUM] = { 0.0, 1.5, 3.0, 5.0 };
  double g_x[G_NUM] = { 0.0, 1.5, 3.0, 5.0 };
  int i;
  double integral;

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  Test PIECEWISE_LINEAR_PRODUCT_INTEGRAL on a simple problem.\n" );
  printf ( "  F and G are both defined over separate, multiple\n" );
  printf ( "  intervals, but still true that F(X) = G(X) = X.\n" );
  printf ( "\n" );
  printf ( "           A           B      Integral        Exact\n" );
  printf ( "\n" );

  a = 1.0;
  for ( i = 1; i <= 5; i++ )
  {
    b = ( double ) ( i );
    integral = piecewise_linear_product_integral ( a, b, f_num, f_x, f_v, g_num, 
      g_x, g_v );
    exact = ( b * b * b - a * a * a ) / 3.0;
    printf ( "  %10g  %10g  %14g  %14g\n", a, b, integral, exact );
  }

  return;
# undef F_NUM
# undef G_NUM
}
/******************************************************************************/

void test03 ( )

/******************************************************************************/
/*
  Purpose:

    TEST03 tests PIECEWISE_LINEAR_PRODUCT_INTEGRAL.

  Discussion:

    For this test, F(X) and G(X) are piecewise linear interpolants to
    SIN(X) and 2 * COS(X), so we know the exact value of the integral
    of the product of the original functions, but this is only an estimate 
    of the exact value of the integral of the product of the interpolants.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    28 April 2009

  Author:

    John Burkardt
*/
{
# define F_NUM 11
# define G_NUM 31

  double a;
  double b;
  double exact;
  int f_num = F_NUM;
  double f_v[F_NUM];
  double f_x[F_NUM];
  int g_num = G_NUM;
  double g_v[G_NUM];
  double g_x[G_NUM];
  int i;
  double integral;
  double pi = 3.141592653589793;
  double quad;
  int quad_num;

  printf ( "\n" );
  printf ( "TEST03\n" );
  printf ( "  Test PIECEWISE_LINEAR_PRODUCT_INTEGRAL on a simple problem.\n" );
  printf ( "  F and G are defined over separate, multiple\n" );
  printf ( "  intervals.\n" );
  printf ( "\n" );
  printf ( "  F(X) interpolates SIN(X),\n" );
  printf ( "  G(X) interpolates 2*COS(X).\n" );
  printf ( "\n" );
  printf ( "  We compare:\n" );
  printf ( "\n" );
  printf ( "  INTEGRAL, our value for the integral,\n" );
  printf ( "  QUAD, a quadrature estimate for the integral, and\n" );
  printf ( "  CLOSE, the value of the integral of 2*COS(X)*SIN(X)\n" );
  printf ( "\n" );
  printf ( "           A           B      Integral        Quad            Close\n" );
  printf ( "\n" );

  for ( i = 0; i < f_num; i++ )
  {
    f_x[i] = ( ( f_num - i - 1 ) * 0.0
             + (         i     ) * pi )
             / ( f_num     - 1 );
    f_v[i] = sin ( f_x[i] );
  }

  for ( i = 0; i < g_num; i++ )
  {
    g_x[i] = ( ( g_num - i - 1 ) * 0.0
             + (         i     ) * pi )
             / ( g_num     - 1 );
    g_v[i] = 2.0 * cos ( g_x[i] );
  }

  a = 0.0;
  for ( i = 0; i <= 6; i++ )
  {
    b = ( double ) ( i ) * pi / 6.0;
    integral = piecewise_linear_product_integral ( a, b, f_num, f_x, f_v, 
      g_num, g_x, g_v );
    exact = - ( cos ( 2.0 * b ) - cos ( 2.0 * a ) ) / 2.0;
    quad_num = 2000;
    quad = piecewise_linear_product_quad ( a, b, f_num, f_x, f_v, g_num, 
      g_x, g_v, quad_num );
    printf ( "  %10g  %10g  %14g  %14g  %14g\n", a, b, integral, quad, exact );
  }

  return;
# undef F_NUM
# undef G_NUM
}
/******************************************************************************/

void test04 ( )

/******************************************************************************/
/*
  Purpose:

    TEST04 tests PIECEWISE_LINEAR_PRODUCT_INTEGRAL.

  Discussion:

    For this test, we compute the integrals of a hat function with itself,
    and a hat function with its neighbor.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    04 July 2013

  Author:

    John Burkardt
*/
{
# define F_NUM 3
# define G_NUM 3

  double a;
  double b;
  double exact;
  int f_num = F_NUM;
  double f_v[F_NUM] = { 0.0, 1.0, 0.0 };
  double f_x[F_NUM] = { 0.0, 1.0, 2.0 };
  int g_num = G_NUM;
  double g_v[G_NUM] = { 1.0, 0.0, 0.0 };
  double g_x[G_NUM] = { 0.0, 1.0, 2.0 };
  int i;
  double integral;

  printf ( "\n" );
  printf ( "TEST04\n" );
  printf ( "  Test PIECEWISE_LINEAR_PRODUCT_INTEGRAL.\n" );
  printf ( "  The nodes are at 0, 1, and 2.\n" );
  printf ( "  F(X) = ( 0, 1, 0 ).\n" );
  printf ( "  G(X) = ( 1, 0, 0 ).\n" );
  printf ( "\n" );

  a = 0.0;
  b = 2.0;

  integral = piecewise_linear_product_integral ( a, b, f_num, f_x, f_v, f_num, 
    f_x, f_v );

  printf ( "  Integral F(X) * F(X) dx = %g\n", integral );

  integral = piecewise_linear_product_integral ( a, b, f_num, f_x, f_v, g_num, 
    g_x, g_v );

  printf ( "  Integral F(X) * G(X) dx = %g\n", integral );

  integral = piecewise_linear_product_integral ( a, b, g_num, g_x, g_v, g_num, 
    g_x, g_v );

  printf ( "  Integral G(X) * G(X) dx = %g\n", integral );

  return;
# undef F_NUM
# undef G_NUM
}
