# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <complex.h>

# include "fn.h"

int main ( void );

void float_sin_test ( void );
void double_sin_test ( void );
void r8_sin_test ( void );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for FN_PRB.

  Discussion:

    FN_PRB calls the FN tests.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    10 March 2010

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "FN_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the FN library.\n" );

  float_sin_test ( );
  double_sin_test ( );

  r8_sin_test ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "FN_PRB\n" );
  printf ( "  Normal end of execution.\n" );

  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void float_sin_test ( void )

/******************************************************************************/
/*
  Purpose:

    FLOAT_SIN_TEST tests the system SIN and COS functions with float arguments.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    10 March 2010

  Author:

    Original FORTRAN77 version by Wayne Cody.
    C version by John Burkardt
*/
{
  float a;
  float ait;
  float albeta;
  float b;
  float beta;
  float betap;
  float c;
  float del;
  float eps;
  float epsneg;
  float expon;
  int i;
  long int ibeta;
  long int iexp;
  long int irnd;
  long int it;
  int i1;
  int j;
  int k1;
  int k2;
  int k3;
  long int machep;
  long int maxexp;
  long int minexp;
  int n;
  long int negep;
  long int ngrd;
  float r6;
  float r7;
  float w;
  float x;
  float x1;
  float xl;
  float xmax;
  float xmin;
  float y;
  float z;
  float zz;

  printf ( "\n" );
  printf ( "REAL_SIN_TEST\n" );
  printf ( "  Test SIN and COS on float arguments.\n" );

  r4_machar ( &ibeta, &it, &irnd, &ngrd, &machep, &negep, &iexp, &minexp,
    &maxexp, &eps, &epsneg, &xmin, &xmax );

  beta = ( float ) ( ibeta );
  albeta = log ( beta );
  ait = ( float ) ( it );
  a = 0.0;
  b = 1.570796327;
  c = 1.570796327;
  n = 2000;
  i1 = 0;
/*
  Random argument accuracy tests.
*/
  for ( j = 1; j <= 3; j++ )
  {
    k1 = 0;
    k2 = 0;
    k3 = 0;
    x1 = 0.0;
    r6 = 0.0;
    r7 = 0.0;
    del = ( b - a ) / ( float ) ( n );
    xl = a;

    for ( i = 1; i < n; i++ )
    {
      x = del * r4_ren ( ) + xl;
      y = x / 3.0;
      y = ( x + y ) - x;
      x = 3.0 * y;

      if ( j < 3 )
      {
        z = sin ( x );
        zz = sin ( y );
        w = 1.0;

        if ( z != 0.0 )
        {
          w = ( z - zz * ( 3.0 - 4.0 * zz * zz ) ) / z;
        }
      }
      else
      {
        z = cos ( x );
        zz = cos ( y );
        w = 1.0;

        if ( z != 0.0 )
        {
          w = ( z + zz * ( 3.0 - 4.0 * zz * zz ) ) / z;
        }

      }

      if ( 0.0 < w )
      {
        k1 = k1 + 1;
      }
      else if ( w == 0.0 )
      {
        k2 = k2 + 1;
      }
      else if ( w < 0.0 )
      {
        k3 = k3 + 1;
      }

      w = r4_abs ( w );

      if ( r6 < w )
      {
        r6 = w;
        x1 = x;
      }

      r7 = r7 + w * w;
      xl = xl + del;

    }

    r7 = sqrt ( r7 / ( float ) ( n ) );

    if ( j < 3 )
    {
      printf ( "\n" );
      printf ( "  Test of sin(x) vs 3*sin(x/3)-4*sin(x/3)^3\n" );
      printf ( "\n" );
      printf ( "  %d random arguments were tested from the interval\n", n );
      printf ( "  ( %g,%g )\n", a, b );
      printf ( "\n" );
      printf ( "  sin(x) was larger %d times,\n", k1 );
      printf ( "  agreed %d times, and \n", k2 );
      printf ( "  was smaller %d times.\n", k3 );
    }
    else
    {
      printf ( "\n" );
      printf ( "  Test of cos(x) vs 4*cos(x/3)^3-3*cos(x/3)\n" );
      printf ( "\n" );
      printf ( "  %d random arguments were tested from the interval\n", n );
      printf ( "  ( %d,%d )\n", a, b );
      printf ( "\n" );
      printf ( "  cos(x) was larger %d times,\n", k1 );
      printf ( "  agreed %d times, and \n", k2 );
      printf ( "  was smaller %d times.\n", k3 );
    }

    printf ( "\n" );
    printf ( "  There are %d base %d significant digits in a floating-point number.\n", 
      it, ibeta );

    w = - 999.0;
    if ( r6 != 0.0 )
    {
      w = log ( r8_abs ( r6 ) ) / albeta;
    }
    printf ( "  The maximum relative error of %d = %d ^ %g\n", r6, ibeta, w );
    printf ( "  occurred for X = %g\n", x1 );
    w = r4_max ( ait + w, 0.0 );
    printf ( "  The estimated loss of base %d significant digits is %g\n", ibeta, w );

    w = -999.0;
    if ( r7 != 0.0 )
    {
      w = log ( r4_abs ( r7 ) ) / albeta;
    }
    printf ( "  The root mean square relative error was %g = %d ^ %g\n", r7, ibeta, w );
    w = r4_max ( ait + w, 0.0 );
    printf ( "  The estimated loss of base %d significant digits is %g\n", ibeta, w );

    a = 18.84955592;
    if ( j == 2 )
    {
      a = b + c;
    }
    b = a + c;
  }
/*
  Special tests.
*/
  printf ( "\n" );
  printf ( "  Special tests:\n" );
  printf ( "\n" );

  c = 1.0 / pow ( beta, ( int ) ( it ) / 2 );
  z = ( sin ( a + c ) - sin ( a - c ) ) / ( c + c );

  printf ( "  if %g is not almost 1.0,\n", z );
  printf ( "  sin has the wrong period.\n" );

  printf ( "  The identity   sin(-x) = -sin(x)   will be tested.\n" );
  printf ( "        x         f(x) + f(-x)\n" );
  printf ( "\n" );
  for ( i = 1; i <= 5; i++ )
  {
    x = r4_ren ( ) * a;
    z = sin ( x ) + sin ( - x );
    printf ( "  %15g  %15g \n", x, z );
  }

  printf ( "\n" );
  printf ( "  The identity sin(x) = x , x small, will be tested.\n" );
  printf ( "        x         x - f(x)\n" );
  printf ( "\n" );
  betap = pow ( beta, ( int ) it );
  x = r4_ren ( ) / betap;

  for ( i = 1; i <= 5; i++ )
  {
    z = x - sin ( x );
    printf ( "  %15g  %15g \n", x, z );
    x = x / beta;
  }

  printf ( "\n" );
  printf ( "  The identity   cos(-x) = cos(x)   will be tested.\n" );
  printf ( "         x         f(x) - f(-x)\n" );
  printf ( "\n" );
  for ( i = 1; i <= 5; i++ )
  {
    x = r4_ren ( ) * a;
    z = cos ( x ) - cos ( - x );
    printf ( "  %15g  %15g \n", x, z );
  }

  printf ( "\n" );
  printf ( "  Test of underflow for very small argument.\n" );
  printf ( "\n" );
  expon = ( float ) ( minexp ) * 0.75;
  x = pow ( beta, expon );
  y = sin ( x );
  printf ( "  sin ( %g ) = %g\n", x, y );
  printf ( "\n" );
  printf ( "  The following three lines illustrate the loss in \n" );
  printf ( "  significance for large arguments.  The arguments,\n" );
  printf ( "  are consecutive.\n" );
  printf ( "\n" );
  z = sqrt ( betap );
  x = z * ( 1.0 - epsneg );
  y = sin ( x );
  printf ( "  sin ( %20.16g ) = %20.16g\n", x, y );
  y = sin ( z );
  printf ( "  sin ( %20.16g ) = %20.16g\n", x, y );
  x = z * ( 1.0 + eps );
  y = sin ( x );
  printf ( "  sin ( %20.16g ) = %20.16g\n", x, y );
/*
  Test of error returns.
*/
  printf ( "\n" );
  printf ( "  Test of error returns:\n" );
  printf ( "\n" );
  x = betap;
  printf ( "  SIN will be called with the argument %g\n", x );
  printf ( "  This should trigger an error message.\n" );
  y = sin ( x );
  printf ( "\n" );
  printf ( "  SIN returned the value %g\n", y );
  printf ( "\n" );
  printf ( "  This concludes the tests.\n" );

  return;
}
/******************************************************************************/

void double_sin_test ( void )

/******************************************************************************/
/*
  Purpose:

    DOUBLE_SIN_TEST tests the system SIN and COS functions with double arguments.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    10 March 2010

  Author:

    Original FORTRAN77 version by Wayne Cody.
    C version by John Burkardt
*/
{
  double a;
  double ait;
  double albeta;
  double b;
  double beta;
  double betap;
  double c;
  double del;
  double eps;
  double epsneg;
  double expon;
  int i;
  long int ibeta;
  long int iexp;
  long int irnd;
  long int it;
  int i1;
  int j;
  int k1;
  int k2;
  int k3;
  long int machep;
  long int maxexp;
  long int minexp;
  int n;
  long int negep;
  long int ngrd;
  double r6;
  double r7;
  double w;
  double x;
  double x1;
  double xl;
  double xmax;
  double xmin;
  double y;
  double z;
  double zz;

  printf ( "\n" );
  printf ( "DOUBLE_SIN_TEST\n" );
  printf ( "  Test SIN and COS on double arguments.\n" );

  r8_machar ( &ibeta, &it, &irnd, &ngrd, &machep, &negep, &iexp, &minexp,
    &maxexp, &eps, &epsneg, &xmin, &xmax );

  beta = ( double ) ( ibeta );
  albeta = log ( beta );
  ait = ( double ) ( it );
  a = 0.0;
  b = 1.570796327;
  c = 1.570796327;
  n = 2000;
  i1 = 0;
/*
  Random argument accuracy tests.
*/
  for ( j = 1; j <= 3; j++ )
  {
    k1 = 0;
    k2 = 0;
    k3 = 0;
    x1 = 0.0;
    r6 = 0.0;
    r7 = 0.0;
    del = ( b - a ) / ( double ) ( n );
    xl = a;

    for ( i = 1; i < n; i++ )
    {
      x = del * r8_ren ( ) + xl;
      y = x / 3.0;
      y = ( x + y ) - x;
      x = 3.0 * y;

      if ( j < 3 )
      {
        z = sin ( x );
        zz = sin ( y );
        w = 1.0;

        if ( z != 0.0 )
        {
          w = ( z - zz * ( 3.0 - 4.0 * zz * zz ) ) / z;
        }
      }
      else
      {
        z = cos ( x );
        zz = cos ( y );
        w = 1.0;

        if ( z != 0.0 )
        {
          w = ( z + zz * ( 3.0 - 4.0 * zz * zz ) ) / z;
        }

      }

      if ( 0.0 < w )
      {
        k1 = k1 + 1;
      }
      else if ( w == 0.0 )
      {
        k2 = k2 + 1;
      }
      else if ( w < 0.0 )
      {
        k3 = k3 + 1;
      }

      w = r8_abs ( w );

      if ( r6 < w )
      {
        r6 = w;
        x1 = x;
      }

      r7 = r7 + w * w;
      xl = xl + del;

    }

    r7 = sqrt ( r7 / ( double ) ( n ) );

    if ( j < 3 )
    {
      printf ( "\n" );
      printf ( "  Test of sin(x) vs 3*sin(x/3)-4*sin(x/3)^3\n" );
      printf ( "\n" );
      printf ( "  %d random arguments were tested from the interval\n", n );
      printf ( "  ( %g,%g )\n", a, b );
      printf ( "\n" );
      printf ( "  sin(x) was larger %d times,\n", k1 );
      printf ( "  agreed %d times, and \n", k2 );
      printf ( "  was smaller %d times.\n", k3 );
    }
    else
    {
      printf ( "\n" );
      printf ( "  Test of cos(x) vs 4*cos(x/3)^3-3*cos(x/3)\n" );
      printf ( "\n" );
      printf ( "  %d random arguments were tested from the interval\n", n );
      printf ( "  ( %g,%g )\n", a, b );
      printf ( "\n" );
      printf ( "  cos(x) was larger %d times,\n", k1 );
      printf ( "  agreed %d times, and \n", k2 );
      printf ( "  was smaller %d times.\n", k3 );
    }

    printf ( "\n" );
    printf ( "  There are %d base %d significant digits in a doubleing-point number.\n", 
      it, ibeta );

    w = - 999.0;
    if ( r6 != 0.0 )
    {
      w = log ( r8_abs ( r6 ) ) / albeta;
    }
    printf ( "  The maximum relative error of %g = %d ^ %g\n", r6, ibeta, w );
    printf ( "  occurred for X = %g\n", x1 );
    w = r8_max ( ait + w, 0.0 );
    printf ( "  The estimated loss of base %d significant digits is %g\n", ibeta, w );

    w = -999.0;
    if ( r7 != 0.0 )
    {
      w = log ( r8_abs ( r7 ) ) / albeta;
    }
    printf ( "  The root mean square relative error was %g = %d ^ %g\n", r7, ibeta, w );
    w = r8_max ( ait + w, 0.0 );
    printf ( "  The estimated loss of base %d significant digits is %g\n", ibeta, w );

    a = 18.84955592;
    if ( j == 2 )
    {
      a = b + c;
    }
    b = a + c;
  }
/*
  Special tests.
*/
  printf ( "\n" );
  printf ( "  Special tests:\n" );
  printf ( "\n" );

  c = 1.0 / pow ( beta, ( int ) ( it ) / 2 );
  z = ( sin ( a + c ) - sin ( a - c ) ) / ( c + c );

  printf ( "  if %g is not almost 1.0,\n", z );
  printf ( "  sin has the wrong period.\n" );

  printf ( "  The identity   sin(-x) = -sin(x)   will be tested.\n" );
  printf ( "        x         f(x) + f(-x)\n" );
  printf ( "\n" );
  for ( i = 1; i <= 5; i++ )
  {
    x = r8_ren ( ) * a;
    z = sin ( x ) + sin ( - x );
    printf ( "  %15g  %15g \n", x, z );
  }

  printf ( "\n" );
  printf ( "  The identity sin(x) = x , x small, will be tested.\n" );
  printf ( "        x         x - f(x)\n" );
  printf ( "\n" );
  betap = pow ( beta, ( int ) it );
  x = r8_ren ( ) / betap;

  for ( i = 1; i <= 5; i++ )
  {
    z = x - sin ( x );
    printf ( "  %15g  %15g \n", x, z );
    x = x / beta;
  }

  printf ( "\n" );
  printf ( "  The identity   cos(-x) = cos(x)   will be tested.\n" );
  printf ( "         x         f(x) - f(-x)\n" );
  printf ( "\n" );
  for ( i = 1; i <= 5; i++ )
  {
    x = r8_ren ( ) * a;
    z = cos ( x ) - cos ( - x );
    printf ( "  %15g  %15g \n", x, z );
  }

  printf ( "\n" );
  printf ( "  Test of underflow for very small argument.\n" );
  printf ( "\n" );
  expon = ( double ) ( minexp ) * 0.75;
  x = pow ( beta, expon );
  y = sin ( x );
  printf ( "  sin ( %g ) = %g\n", x, y );
  printf ( "\n" );
  printf ( "  The following three lines illustrate the loss in \n" );
  printf ( "  significance for large arguments.  The arguments,\n" );
  printf ( "  are consecutive.\n" );
  printf ( "\n" );
  z = sqrt ( betap );
  x = z * ( 1.0 - epsneg );
  y = sin ( x );
  printf ( "  sin ( %20.16g ) = %20.16g\n", x, y );
  y = sin ( z );
  printf ( "  sin ( %20.16g ) = %20.16g\n", x, y );
  x = z * ( 1.0 + eps );
  y = sin ( x );
  printf ( "  sin ( %20.16g ) = %20.16g\n", x, y );
/*
  Test of error returns.
*/
  printf ( "\n" );
  printf ( "  Test of error returns:\n" );
  printf ( "\n" );
  x = betap;
  printf ( "  SIN will be called with the argument %g\n", x );
  printf ( "  This should trigger an error message.\n" );
  y = sin ( x );
  printf ( "\n" );
  printf ( "  SIN returned the value %g\n", y );
  printf ( "\n" );
  printf ( "  This concludes the tests.\n" );

  return;
}
/******************************************************************************/

void r8_sin_test ( void )

/******************************************************************************/
/*
  Purpose:

    R8_SIN_TEST tests the R8_SIN and R8_COS functions.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    10 March 2010

  Author:

    Original FORTRAN77 version by Wayne Cody.
    C version by John Burkardt
*/
{
  double a;
  double ait;
  double albeta;
  double b;
  double beta;
  double betap;
  double c;
  double del;
  double eps;
  double epsneg;
  double expon;
  int i;
  long int ibeta;
  long int iexp;
  long int irnd;
  long int it;
  int i1;
  int j;
  int k1;
  int k2;
  int k3;
  long int machep;
  long int maxexp;
  long int minexp;
  int n;
  long int negep;
  long int ngrd;
  double r6;
  double r7;
  double w;
  double x;
  double x1;
  double xl;
  double xmax;
  double xmin;
  double y;
  double z;
  double zz;

  printf ( "\n" );
  printf ( "R8_SIN_TEST\n" );
  printf ( "  Test R8_SIN and R8_COS on double arguments.\n" );

  r8_machar ( &ibeta, &it, &irnd, &ngrd, &machep, &negep, &iexp, &minexp,
    &maxexp, &eps, &epsneg, &xmin, &xmax );

  beta = ( double ) ( ibeta );
  albeta = log ( beta );
  ait = ( double ) ( it );
  a = 0.0;
  b = 1.570796327;
  c = 1.570796327;
  n = 2000;
  i1 = 0;
/*
  Random argument accuracy tests.
*/
  for ( j = 1; j <= 3; j++ )
  {
    k1 = 0;
    k2 = 0;
    k3 = 0;
    x1 = 0.0;
    r6 = 0.0;
    r7 = 0.0;
    del = ( b - a ) / ( double ) ( n );
    xl = a;

    for ( i = 1; i < n; i++ )
    {
      x = del * r8_ren ( ) + xl;
      y = x / 3.0;
      y = ( x + y ) - x;
      x = 3.0 * y;

      if ( j < 3 )
      {
        z = r8_sin ( x );
        zz = r8_sin ( y );
        w = 1.0;

        if ( z != 0.0 )
        {
          w = ( z - zz * ( 3.0 - 4.0 * zz * zz ) ) / z;
        }
      }
      else
      {
        z = r8_cos ( x );
        zz = r8_cos ( y );
        w = 1.0;

        if ( z != 0.0 )
        {
          w = ( z + zz * ( 3.0 - 4.0 * zz * zz ) ) / z;
        }

      }

      if ( 0.0 < w )
      {
        k1 = k1 + 1;
      }
      else if ( w == 0.0 )
      {
        k2 = k2 + 1;
      }
      else if ( w < 0.0 )
      {
        k3 = k3 + 1;
      }

      w = r8_abs ( w );

      if ( r6 < w )
      {
        r6 = w;
        x1 = x;
      }

      r7 = r7 + w * w;
      xl = xl + del;

    }

    r7 = sqrt ( r7 / ( double ) ( n ) );

    if ( j < 3 )
    {
      printf ( "\n" );
      printf ( "  Test of sin(x) vs 3*sin(x/3)-4*sin(x/3)^3\n" );
      printf ( "\n" );
      printf ( "  %d random arguments were tested from the interval\n", n );
      printf ( "  ( %g, %g )\n", a, b );
      printf ( "\n" );
      printf ( "  sin(x) was larger %d times,\n", k1 );
      printf ( "  agreed %d times, and \n", k2 );
      printf ( "  was smaller %d times.\n", k3 );
    }
    else
    {
      printf ( "\n" );
      printf ( "  Test of cos(x) vs 4*cos(x/3)^3-3*cos(x/3)\n" );
      printf ( "\n" );
      printf ( "  %d random arguments were tested from the interval\n", n );
      printf ( "  ( %g,%g )\n", a, b );
      printf ( "\n" );
      printf ( "  cos(x) was larger %d times,\n", k1 );
      printf ( "  agreed %d times, and \n", k2 );
      printf ( "  was smaller %d times.\n", k3 );
    }

    printf ( "\n" );
    printf ( "  There are %d base %d significant digits in a double.\n", it, ibeta );

    w = - 999.0;
    if ( r6 != 0.0 )
    {
      w = log ( r8_abs ( r6 ) ) / albeta;
    }
    printf ( "  The maximum relative error of %g = %d ^ %g\n", r6, ibeta, w );
    printf ( "  occurred for X = %g\n", x1 );
    w = r8_max ( ait + w, 0.0 );
    printf ( "  The estimated loss of base %d significant digits is %g\n", ibeta, w );

    w = -999.0;
    if ( r7 != 0.0 )
    {
      w = log ( r8_abs ( r7 ) ) / albeta;
    }
    printf ( "  The root mean square relative error was %g = %d ^ %g\n", r7, ibeta, w );
    w = r8_max ( ait + w, 0.0 );
    printf ( "  The estimated loss of base %d significant digits is %g\n", ibeta, w );

    a = 18.84955592;
    if ( j == 2 )
    {
      a = b + c;
    }
    b = a + c;
  }
/*
  Special tests.
*/
  printf ( "\n" );
  printf ( "  Special tests:\n" );
  printf ( "\n" );

  c = 1.0 / pow ( beta, ( int ) ( it ) / 2 );
  z = ( r8_sin ( a + c ) - r8_sin ( a - c ) ) / ( c + c );

  printf ( "  if %g is not almost 1.0,\n", z );
  printf ( "  sin has the wrong period.\n" );

  printf ( "  The identity   sin(-x) = -sin(x)   will be tested.\n" );
  printf ( "        x         f(x) + f(-x)\n" );
  printf ( "\n" );
  for ( i = 1; i <= 5; i++ )
  {
    x = r8_ren ( ) * a;
    z = r8_sin ( x ) + r8_sin ( - x );
    printf ( "  %15g  %15g \n", x, z );
  }

  printf ( "\n" );
  printf ( "  The identity sin(x) = x , x small, will be tested.\n" );
  printf ( "        x         x - f(x)\n" );
  printf ( "\n" );
  betap = pow ( beta, ( int ) it );
  x = r8_ren ( ) / betap;

  for ( i = 1; i <= 5; i++ )
  {
    z = x - r8_sin ( x );
    printf ( "  %15g  %15g \n", x, z );
    x = x / beta;
  }

  printf ( "\n" );
  printf ( "  The identity   cos(-x) = cos(x)   will be tested.\n" );
  printf ( "         x         f(x) - f(-x)\n" );
  printf ( "\n" );
  for ( i = 1; i <= 5; i++ )
  {
    x = r8_ren ( ) * a;
    z = r8_cos ( x ) - r8_cos ( - x );
    printf ( "  %15g  %15g \n", x, z );
  }

  printf ( "\n" );
  printf ( "  Test of underflow for very small argument.\n" );
  printf ( "\n" );
  expon = ( double ) ( minexp ) * 0.75;
  x = pow ( beta, expon );
  y = r8_sin ( x );
  printf ( "  r8_sin ( %g ) = %g\n", x, y );
  printf ( "\n" );
  printf ( "  The following three lines illustrate the loss in \n" );
  printf ( "  significance for large arguments.  The arguments,\n" );
  printf ( "  are consecutive.\n" );
  printf ( "\n" );
  z = sqrt ( betap );
  x = z * ( 1.0 - epsneg );
  y = r8_sin ( x );
  printf ( "  sin ( %20.16g ) = %20.16g\n", x, y );
  y = r8_sin ( z );
  printf ( "  sin ( %20.16g ) = %20.16g\n", x, y );
  x = z * ( 1.0 + eps );
  y = r8_sin ( x );
  printf ( "  sin ( %20.16g ) = %20.16g\n", x, y );
/*
  Test of error returns.
*/
  printf ( "\n" );
  printf ( "  Test of error returns:\n" );
  printf ( "\n" );
  x = betap;
  printf ( "  SIN will be called with the argument %g\n", x );
  printf ( "  This should trigger an error message.\n" );
  y = r8_sin ( x );
  printf ( "\n" );
  printf ( "  SIN returned the value %g\n", y );
  printf ( "\n" );
  printf ( "  This concludes the tests.\n" );

  return;
}
