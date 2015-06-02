# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "cordic.h"

int main ( );

void test001 ( );
void test002 ( );
void test003 ( );
void test004 ( );
void test005 ( );
void test006 ( );
void test007 ( );
void test008 ( );
void test009 ( );
void test010 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for CORDIC_PRB.

  Discussion:

    CORDIC_PRB tests the CORDIC library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 January 2012

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "CORDIC_PRB:\n" );
  printf ( "  C version,\n" );
  printf ( "  Test the CORDIC library.\n" );

  test001 ( );
  test002 ( );
  test003 ( );
  test004 ( );
  test005 ( );
  test006 ( );
  test007 ( );
  test008 ( );
  test009 ( );
  test010 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "CORDIC_PRB:\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void test001 ( )

/******************************************************************************/
/*
  Purpose:

    TEST001 demonstrates the use of COSSIN_CORDIC.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 January 2012

  Author:

    John Burkardt
*/
{
  double a;
  double c1;
  double c2;
  double d;
  int n;
  int n_data;
  double s2;

  printf ( "\n" );
  printf ( "TEST001:\n" );
  printf ( "  COSSIN_CORDIC computes the cosine and sine\n" );
  printf ( "  using the CORDIC algorithm.\n" );
  printf ( "\n" );
  printf ( "          A        N      Cos(A)           Cos(A)           Difference\n" );
  printf ( "                          Tabulated        CORDIC\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    cos_values ( &n_data, &a, &c1 );

    if ( n_data == 0 )
    {
      break;
    }

    printf ( "\n" );
    for ( n = 0; n <= 25; n = n + 5 )
    {
      cossin_cordic ( a, n, &c2, &s2 );

      d = c1 - c2;

      printf ( "  %12.4g  %4d  %16.8g  %16.8g  %12.4g\n", a, n, c1, c2, d );
    }
  }
  return;
}
/******************************************************************************/

void test002 ( )

/******************************************************************************/
/*
  Purpose:

    TEST002 demonstrates the use of COSSIN_CORDIC.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 January 2012

  Author:

    John Burkardt
*/
{
  double a;
  double c2;
  double d;
  int n;
  int n_data;
  double s1;
  double s2;

  printf ( "\n" );
  printf ( "TEST002:\n" );
  printf ( "  COSSIN_CORDIC computes the cosine and sine\n" );
  printf ( "  using the CORDIC algorithm.\n" );
  printf ( "\n" );
  printf ( "          A        N      Sin(A)           Sin(A)           Difference\n" );
  printf ( "                          Tabulated        CORDIC\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    sin_values ( &n_data, &a, &s1 );

    if ( n_data == 0 )
    {
      break;
    }

    printf ( "\n" );
    for ( n = 0; n <= 25; n = n + 5 )
    {
      cossin_cordic ( a, n, &c2, &s2 );

      d = s1 - s2;

      printf ( "  %12.4g  %4d  %16.8g  %16.8g  %12.4g\n", a, n, s1, s2, d );
    }
  }
  return;
}
/******************************************************************************/

void test003 ( )

/******************************************************************************/
/*
  Purpose:

    TEST003 demonstrates the use of ARCTAN_CORDIC.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 January 2012

  Author:

    John Burkardt
*/
{
  double a1;
  double a2;
  double d;
  int n;
  int n_data;
  double r;
  double s;
  int seed;
  double x;
  double y;
  double z;

  printf ( "\n" );
  printf ( "TEST003:\n" );
  printf ( "  ARCTAN_CORDIC computes the arctangent of Y/X\n" );
  printf ( "  using the CORDIC algorithm.\n" );
  printf ( "\n" );
  printf ( "      X      Y    N       ArcTan(Y/X) ArcTan(Y/X)      Difference\n" );
  printf ( "                           Tabulated   CORDIC\n" );
  printf ( "\n" );

  seed = 123456789;
  n_data = 0;

  for ( ; ; )
  {
    arctan_values ( &n_data, &z, &a1 );

    if ( n_data == 0 )
    {
      break;
    }

    r = r8_uniform_01 ( &seed );

    x = r;
    y = r * z;

    s = r8_uniform_01 ( &seed );

    if ( s < 0.5 )
    {
      x = -x;
      y = -y;
    }

    printf ( "\n" );
    for ( n = 0; n <= 25; n = n + 5 )
    {
      a2 = arctan_cordic ( x, y, n );

      d = a1 - a2;

      printf ( "  %12.4g  %12.4g  %4d  %16.8g  %16.8g  %12.4g\n", x, y, n, a1, a2, d );
    }
  }
  return;
}
/******************************************************************************/

void test004 ( )

/******************************************************************************/
/*
  Purpose:

    TEST004 demonstrates the use of ARCCOS_CORDIC.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 January 2012

  Author:

    John Burkardt
*/
{
  double a1;
  double a2;
  double d;
  int n;
  int n_data;
  double t;

  printf ( "\n" );
  printf ( "TEST004:\n" );
  printf ( "  ARCCOS_CORDIC computes the arccosine of T\n" );
  printf ( "  using the CORDIC algorithm.\n" );
  printf ( "\n" );
  printf ( "      T    N        ArcCos(T)  ArcCos(T)      Difference\n" );
  printf ( "                   Tabulated   CORDIC\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    arccos_values ( &n_data, &t, &a1 );

    if ( n_data == 0 )
    {
      break;
    }

    printf ( "\n" );
    for ( n = 0; n <= 25; n = n + 5 )
    {
      a2 = arccos_cordic ( t, n );

      d = a1 - a2;

      printf ( "  %12.4g  %4d  %16.8g  %16.8g  %12.4g\n", t, n, a1, a2, d );
    }
  }
  return;
}
/******************************************************************************/

void test005 ( )

/******************************************************************************/
/*
  Purpose:

    TEST005 demonstrates the use of ARCSIN_CORDIC.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 January 2012

  Author:

    John Burkardt
*/
{
  double a1;
  double a2;
  double d;
  int n;
  int n_data;
  double t;

  printf ( "\n" );
  printf ( "TEST005:\n" );
  printf ( "  ARCSIN_CORDIC computes the arcsine of T\n" );
  printf ( "  using the CORDIC algorithm.\n" );
  printf ( "\n" );
  printf ( "      T    N        ArcSin(T)  ArcSin(T)      Difference\n" );
  printf ( "                   Tabulated   CORDIC\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    arcsin_values ( &n_data, &t, &a1 );

    if ( n_data == 0 )
    {
      break;
    }

    printf ( "\n" );
    for ( n = 0; n <= 25; n = n + 5 )
    {
      a2 = arcsin_cordic ( t, n );

      d = a1 - a2;

      printf ( "  %12.4g  %4d  %16.8g  %16.8g  %12.4g\n", t, n, a1, a2, d );
    }
  }
  return;
}
/******************************************************************************/

void test006 ( )

/******************************************************************************/
/*
  Purpose:

    TEST006 demonstrates the use of TAN_CORDIC.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 January 2012

  Author:

    John Burkardt
*/
{
  double d;
  int n;
  int n_data;
  double t1;
  double t2;
  double theta;

  printf ( "\n" );
  printf ( "TEST006:\n" );
  printf ( "  TAN_CORDIC computes the tangent of THETA\n" );
  printf ( "  using the CORDIC algorithm.\n" );
  printf ( "\n" );
  printf ( "  THETA    N         Tan(THETA)  Tan(THETA)      Difference\n" );
  printf ( "                     Tabulated   CORDIC\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    tan_values ( &n_data, &theta, &t1 );

    if ( n_data == 0 )
    {
      break;
    }

    printf ( "\n" );
    for ( n = 0; n <= 25; n = n + 5 )
    {
      t2 = tan_cordic ( theta, n );

      d = t1 - t2;

      printf ( "  %12.4g  %4d  %16.8g  %16.8g  %12.4g\n", theta, n, t1, t2, d );
    }
  }
  return;
}
/******************************************************************************/

void test007 ( )

/******************************************************************************/
/*
  Purpose:

    TEST007 demonstrates the use of EXP_CORDIC.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 January 2012

  Author:

    John Burkardt
*/
{
  double d;
  double fx1;
  double fx2;
  int n;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST007:\n" );
  printf ( "  EXP_CORDIC computes the exponential function\n" );
  printf ( "  using the CORDIC algorithm.\n" );
  printf ( "\n" );
  printf ( "    X      N           Exp(X)      Exp(X)        Difference\n" );
  printf ( "                     Tabulated   CORDIC\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    exp_values ( &n_data, &x, &fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    printf ( "\n" );
    for ( n = 0; n <= 25; n = n + 5 )
    {
      fx2 = exp_cordic ( x, n );

      d = fx1 - fx2;

      printf ( "  %12.4g  %4d  %16.8g  %16.8g  %12.4g\n", x, n, fx1, fx2, d );
    }
  }
  return;
}
/******************************************************************************/

void test008 ( )

/******************************************************************************/
/*
  Purpose:

    TEST008 demonstrates the use of LN_CORDIC.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    21 January 2012

  Author:

    John Burkardt
*/
{
  double d;
  double fx1;
  double fx2;
  int n;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST008:\n" );
  printf ( "  LN_CORDIC computes the natural logarithm\n" );
  printf ( "  using the CORDIC algorithm.\n" );
  printf ( "\n" );
  printf ( "    X      N            Ln(X)       Ln(X)        Difference\n" );
  printf ( "                     Tabulated   CORDIC\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    ln_values ( &n_data, &x, &fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    printf ( "\n" );
    for ( n = 0; n <= 25; n = n + 5 )
    {
      fx2 = ln_cordic ( x, n );

      d = fx1 - fx2;

      printf ( "  %12.4g  %4d  %16.8g  %16.8g  %12.4g\n", x, n, fx1, fx2, d );
    }
  }
  return;
}
/******************************************************************************/

void test009 ( )

/******************************************************************************/
/*
  Purpose:

    TEST009 demonstrates the use of SQRT_CORDIC.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    21 January 2012

  Author:

    John Burkardt
*/
{
  double d;
  double fx1;
  double fx2;
  int n;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST009:\n" );
  printf ( "  SQRT_CORDIC computes the square root\n" );
  printf ( "  using the CORDIC algorithm.\n" );
  printf ( "\n" );
  printf ( "    X      N          Sqrt(X)     Sqrt(X)        Difference\n" );
  printf ( "                     Tabulated   CORDIC\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    sqrt_values ( &n_data, &x, &fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    printf ( "\n" );
    for ( n = 0; n <= 25; n = n + 5 )
    {
      fx2 = sqrt_cordic ( x, n );

      d = fx1 - fx2;

      printf ( "  %12.4g  %4d  %16.8g  %16.8g  %12.4g\n", x, n, fx1, fx2, d );
    }
  }
  return;
}
/******************************************************************************/

void test010 ( )

/******************************************************************************/
/*
  Purpose:

    TEST010 demonstrates the use of CBRT_CORDIC.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 January 2012

  Author:

    John Burkardt
*/
{
  double d;
  double fx1;
  double fx2;
  int n;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST010:\n" );
  printf ( "  CBRT_CORDIC computes the cube root\n" );
  printf ( "  using the CORDIC algorithm.\n" );
  printf ( "\n" );
  printf ( "    X      N          Cbrt(X)     Cbrt(X)        Difference\n" );
  printf ( "                     Tabulated   CORDIC\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    cbrt_values ( &n_data, &x, &fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    printf ( "\n" );
    for ( n = 0; n <= 25; n = n + 5 )
    {
      fx2 = cbrt_cordic ( x, n );

      d = fx1 - fx2;

      printf ( "  %12.4g  %4d  %16.8g  %16.8g  %12.4g\n", x, n, fx1, fx2, d );
    }
  }
  return;
}
