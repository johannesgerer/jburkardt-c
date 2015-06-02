# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "bisection_rc.h"

int main ( );
void test01 ( );
double f01 ( double x );
void test02 ( );
void test03 ( );
void test04 ( );
void test05 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for BISECTION_RC_PRB.

  Discussion:

    BISECTION_RC_PRB tests the BISECTION_RC library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 January 2015

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "BISECTION_RC_PRB:\n" );
  printf ( "  C version.\n" );
  printf ( "  Test the BISECTION_RC library.\n" );

  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
  test05 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "BISECTION_RC_PRB:\n" );
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

    TEST01 tests BISECTION_RC, evaluating the function in a separate routine.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 January 2013

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double dx;
  double dx_tol;
  double fx;
  double fx_tol;
  int it;
  int it_max;
  int job;
  double x;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  Demonstrate BISECTION_RC on a simple example.\n" );
  printf ( "  The function is evaluated in a separate routine.\n" );

  fx_tol = 1.0E-08;
  dx_tol = 1.0E-06;
  it = 0;
  it_max = 30;

  a = 0.0;
  b = 1.0;
  fx = 0.0;
  job = 0;

  printf ( "\n" );
  printf ( "     I      X               FX              DX\n" );
  printf ( "\n" );

  for ( ; ; )
  {
    x = bisection_rc ( &a, &b, fx, &job );

    if ( job < 0 )
    {
      printf ( "\n" );
      printf ( "  Error return.\n" );
      break;
    }

    it = it + 1;

    fx = f01 ( x );

    if ( it <= 2 )
    {
      dx = fabs ( b - a );
    }
    else
    {
      dx = 0.5 * fabs ( b - a );
    }

    printf ( "  %4d  %14.6g  %14.6g  %14.6g\n", it, x, fx, dx );

    if ( fabs ( fx ) <= fx_tol )
    {
      printf ( "\n" );
      printf ( "  Function is small.\n" );
      break;
    }

    if ( dx <= dx_tol )
    {
      printf ( "\n" );
      printf ( "  Interval is tiny.\n" );
      break;
    }

    if ( it_max <= it )
    {
      printf ( "\n" );
      printf ( "  Reached iteration limit.\n" );
      break;
    }

  }

  printf ( "\n" );
  printf ( "  A = %14.6g, F(A) = %14.6g\n", a, f01 ( a ) );
  printf ( "  X = %14.6g, F(X) = %14.6g\n", x, f01 ( x ) );
  printf ( "  B = %14.6g, F(B) = %14.6g\n", b, f01 ( b ) );

  return;
}
/******************************************************************************/

double f01 ( double x )

/******************************************************************************/
/*
  Purpose:

    F01 evaluates the function f(x) = cos ( x ) - x which is zero around 0.74

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 January 2013

  Author:

    John Burkardt

  Parameters:

    Input, double X, the evaluation point.

    Output, double F01, the function value.
*/
{
  double value;

  value = cos ( x ) - x;

  return value;
}
/******************************************************************************/

void test02 ( )

/******************************************************************************/
/*
  Purpose:

    TEST02 tests BISECTION_RC, evaluating the function within the routine.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 January 2013

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double dx;
  double dx_tol;
  double fx;
  double fx_tol;
  int it;
  int it_max;
  int job;
  double x;

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  Demonstrate BISECTION_RC on a simple example.\n" );
  printf ( "  The function is evaluated within this routine.\n" );

  fx_tol = 1.0E-09;
  dx_tol = 1.0E-09;
  it = 0;
  it_max = 30;

  a = 0.0;
  b = 1.0;
  fx = 0.0;
  job = 0;

  printf ( "\n" );
  printf ( "     I      X               FX              DX\n" );
  printf ( "\n" );

  for ( ; ; )
  {
    x = bisection_rc ( &a, &b, fx, &job );

    if ( job < 0 )
    {
      printf ( "\n" );
      printf ( "  Error return.\n" );
      break;
    }

    it = it + 1;

    fx = cos ( 100.0 * x ) - 4.0 * erf ( 30.0 * x - 10.0 );

    if ( it <= 2 )
    {
      dx = fabs ( b - a );
    }
    else
    {
      dx = 0.5 * fabs ( b - a );
    }

    printf ( "  %4d  %14.6g  %14.6g  %14.6g\n", it, x, fx, dx );

    if ( fabs ( fx ) <= fx_tol )
    {
      printf ( "\n" );
      printf ( "  Function is small.\n" );
      break;
    }

    if ( dx <= dx_tol )
    {
      printf ( "\n" );
      printf ( "  Interval is tiny.\n" );
      break;
    }

    if ( it_max <= it )
    {
      printf ( "\n" );
      printf ( "  Reached iteration limit.\n" );
      break;
    }

  }

  printf ( "\n" );
  fx = cos ( 100.0 * a ) - 4.0 * erf ( 30.0 * a - 10.0 );
  printf ( "  A = %14.6g, F(A) = %14.6g\n", a, fx );
  fx = cos ( 100.0 * x ) - 4.0 * erf ( 30.0 * x - 10.0 );
  printf ( "  X = %14.6g, F(X) = %14.6g\n", x, fx );
  fx = cos ( 100.0 * b ) - 4.0 * erf ( 30.0 * b - 10.0 );
  printf ( "  A = %14.6g, F(A) = %14.6g\n", b, fx );

  return;
}
/******************************************************************************/

void test03 ( )

/******************************************************************************/
/*
  Purpose:

    TEST03 tests BISECTION_RC, to invert the cardioid CDF.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 January 2013

  Author:

    John Burkardt
*/
{
  double a;
  double alpha = 0.0;
  double b;
  double beta = 0.25;
  double cdf;
  double dx;
  double dx_tol;
  double fx;
  double fx_tol;
  int it;
  int it_max;
  int job;
  const double r8_pi = 3.141592653589793;
  double x;

  printf ( "\n" );
  printf ( "TEST03\n" );
  printf ( "  Demonstrate BISECTION_RC on a probability example.\n" );
  printf ( "\n" );
  printf ( "  The cardioid probability density function has a\n" );
  printf ( "  cumulative density function of the form:\n" );
  printf ( "    CDF(X) = ( pi + x - alpha + 2 beta * sin ( x - alpha ) ) / ( 2 * pi )\n" );
  printf ( "  where alpha and beta are parameters, and x is a value\n" );
  printf ( "  in the range -pi <= x <= +pi.\n" );
  printf ( "\n" );
  printf ( "  CDF(X) is the probability that a random sample will have\n" );
  printf ( "  a value less than or equal to X.\n" );
  printf ( "\n" );
  printf ( "  As X moves from -pi to +pi,\n" );
  printf ( "  the CDF rises from 0 (no probability)\n" );
  printf ( "  to 1 (certain probability).\n" );
  printf ( "\n" );
  printf ( "  Assuming that:\n" );
  printf ( "  * ALPHA = %g\n", alpha ); 
  printf ( "  * BETA =  %g\n", beta );
  printf ( "  determine the value X where the Cardioid CDF is exactly 0.75.\n" );

  fx_tol = 1.0E-05;
  dx_tol = 1.0E-08;
  it = 0;
  it_max = 30;

  job = 0;
  a = - r8_pi;
  b = + r8_pi;

  fx = 0.0;

  printf ( "\n" );
  printf ( "     I      X               FX              DX\n" );
  printf ( "\n" );

  for ( ; ; )
  {
    x = bisection_rc ( &a, &b, fx, &job );

    if ( job < 0 )
    {
      printf ( "\n" );
      printf ( "  Error return.\n" );
      break;
    }

    it = it + 1;

    cdf = ( r8_pi + x - alpha + 2.0 * beta * sin ( x - alpha ) ) / ( 2.0 * r8_pi );
    fx = cdf - 0.75;

    if ( it <= 2 )
    {
      dx = fabs ( b - a );
    }
    else
    {
      dx = 0.5 * fabs ( b - a );
    }

    printf ( "  %4d  %14.6g  %14.6g  %14.6g\n", it, x, fx, dx );

    if ( fabs ( fx ) <= fx_tol )
    {
      printf ( "\n" );
      printf ( "  Function is small.\n" );
      break;
    }

    if ( dx <= dx_tol )
    {
      printf ( "\n" );
      printf ( "  Interval is tiny.\n" );
      break;
    }

    if ( it_max <= it )
    {
      printf ( "\n" );
      printf ( "  Reached iteration limit.\n" );
      break;
    }

  }

  printf ( "\n" );
  cdf = ( r8_pi + a - alpha + 2.0 * beta * sin ( a - alpha ) ) / ( 2.0 * r8_pi );
  fx = cdf - 0.75;
  printf ( "  A = %14.6g, F(A) = %14.6g\n", a, fx );
  cdf = ( r8_pi + x - alpha + 2.0 * beta * sin ( x - alpha ) ) / ( 2.0 * r8_pi );
  fx = cdf - 0.75;
  printf ( "  X = %14.6g, F(X) = %14.6g\n", x, fx );
  cdf = ( r8_pi + b - alpha + 2.0 * beta * sin ( b - alpha ) ) / ( 2.0 * r8_pi );
  fx = cdf - 0.75;
  printf ( "  B = %14.6g, F(B) = %14.6g\n", b, fx );

  printf ( "\n" );
  printf ( "  Look at the actual cardioid CDF value now:\n" );
  printf ( "\n" );
  cdf = ( r8_pi + x - alpha + 2.0 * beta * sin ( x - alpha ) ) / ( 2.0 * r8_pi );
  printf ( "  Cardioid(%g) = %g\n", x, cdf );
  return;
}
/******************************************************************************/

void test04 ( )

/******************************************************************************/
/*
  Purpose:

    TEST04 tests BISECTION_RC for the pipe freezing problem.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 January 2015

  Author:

    John Burkardt

  Reference:

    Cleve Moler,
    Numerical Computing with MATLAB,
    SIAM, 2004,
    ISBN13: 978-0-898716-60-3,
    LC: QA297.M625,
    ebook: http://www.mathworks.com/moler/chapters.html
*/
{
  double a;
  double alpha;
  double b;
  double dx;
  double dx_tol;
  double fx;
  double fx_tol;
  int it;
  int it_max;
  int job;
  double t;
  double tc;
  double ti;
  double x;

  printf ( "\n" );
  printf ( "BISECTION_RC_TEST04\n" );
  printf ( "  The freezing pipe problem.\n" );
  printf ( "\n" );
  printf ( "  At the beginning of a cold spell, the soil is at a uniform\n" );
  printf ( "  temperature of Ti.  The cold spell applies a uniform air\n" );
  printf ( "  temperature of Tc, which begins to cool the soil.\n" );
  printf ( "  As a function of depth x and time t, the soil temperature\n" );
  printf ( "  will now cool down as:\n" );
  printf ( "    ( T(x,t) - Tc ) / ( Ti - Tc ) = erf ( 0.5 * x / sqrt ( alpha * t ) ).\n" );
  printf ( "  where:\n" );
  printf ( "    Ti =  20 degrees centigrade,\n" );
  printf ( "    Tc = -15 degrees centigrade,\n" );
  printf ( "    alpha = 0.000000138 meter^2 / second, thermal conductivity;\n" );
  printf ( "    and erf() is the error function.\n" );
  printf ( "  Water freezes at 0 degrees centigrade.\n" );
  printf ( "\n" );
  printf ( "  What depth x in meters must a water pipe be buried so that it will\n" );
  printf ( "  not freeze even if this cold snap lasts for 60 days?\n" );
/*
  Problem parameters.
*/
  ti = 20.0;
  tc = -15.0;
  t = 60.0 * 24.0 * 60.0 * 60.0;
  alpha = 0.000000138;
/*
  Iteration parameters.
*/
  fx_tol = 1.0E-09;
  dx_tol = 1.0E-09;
  it = 0;
  it_max = 30;
  job = 0;
  fx = 0.0;
/*
  Initial guess for interval.
*/
  a = 0.0;
  b = 1000.0;

  printf ( "\n" );
  printf ( "     I      X               FX              DX\n" );
  printf ( "\n" );

  for ( ; ; )
  {
    x = bisection_rc ( &a, &b, fx, &job );

    if ( job < 0 )
    {
      printf ( "\n" );
      printf ( "  Error return.\n" );
      break;
    }

    it = it + 1;

    fx = tc + ( ti - tc ) * erf ( 0.5 * x / sqrt ( alpha * t ) );

    if ( it <= 2 )
    {
      dx = fabs ( b - a );
    }
    else
    {
      dx = 0.5 * fabs ( b - a );
    }

    printf ( "  %4d  %14.6g  %14.6g  %14.6g\n", it, x, fx, dx );

    if ( fabs ( fx ) <= fx_tol )
    {
      printf ( "\n" );
      printf ( "  Function is small.\n" );
      break;
    }

    if ( dx <= dx_tol )
    {
      printf ( "\n" );
      printf ( "  Interval is tiny.\n" );
      break;
    }

    if ( it_max <= it )
    {
      printf ( "\n" );
      printf ( "  Reached iteration limit.\n" );
      break;
    }

  }

  printf ( "\n" );
  fx = tc + ( ti - tc ) * erf ( 0.5 * a / sqrt ( alpha * t ) );
  printf ( "  A = %g, F(A) = %g\n", a, fx );
  fx = tc + ( ti - tc ) * erf ( 0.5 * x / sqrt ( alpha * t ) );
  printf ( "  X = %g, F(X) = %g\n", x, fx );
  fx = tc + ( ti - tc ) * erf ( 0.5 * b / sqrt ( alpha * t ) );
  printf ( "  B = %g, F(B) = %g\n", b, fx );

  return;
}
/******************************************************************************/

void test05 ( )

/******************************************************************************/
/*
  Purpose:

    TEST05 tests BISECTION_RC for Kepler's problem.

  Discussion:

    Kepler's equation has the form:

      X = M + E * sin ( X )

    X represents the eccentric anomaly of a planet, the angle between the
    perihelion (the point on the orbit nearest to the sun) through the sun 
    to the center of the ellipse, and the line from the center of the ellipse
    to the planet.

    There are two parameters, E and M:

    * E is the eccentricity of the orbit, which should be between 0 and 1.0;

    * M is the angle from the perihelion made by a fictitious planet traveling
      on a circular orbit centered at the sun, and traveling at a constant
      angular velocity equal to the average angular velocity of the true
      planet.  M is usually between 0 and 180 degrees, but can have any value.

    For convenience, X and M are measured in degrees.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 January 2015

  Author:

    John Burkardt

  Reference:

    Cleve Moler,
    Numerical Computing with MATLAB,
    SIAM, 2004,
    ISBN13: 978-0-898716-60-3,
    LC: QA297.M625,
    ebook: http://www.mathworks.com/moler/chapters.html
*/
{
  double ad;
  double ar;
  double bd;
  double br;
  double dx;
  double dx_tol;
  double e;
  double fx;
  double fx_tol;
  int it;
  int it_max;
  int job;
  double md;
  double mr;
  const double r8_pi = 3.141592653589793;
  double xd;
  double xr;

  printf ( "\n" );
  printf ( "TEST05\n" );
  printf ( "  The Kepler equation.\n" );
  printf ( "\n" );
  printf ( "  Kepler's equation has the form'\n" );
  printf ( "\n" );
  printf ( "    X = M + E * sin ( X )\n" );
  printf ( "\n" );
  printf ( "  X represents the eccentric anomaly of a planet, the angle between the\n" );
  printf ( "  perihelion (the point on the orbit nearest to the sun) through the sun\n" );
  printf ( "  to the center of the ellipse, and the line from the center of the ellipse\n" );
  printf ( "  to the planet.\n" );
  printf ( "\n" );
  printf ( "  There are two parameters, E and M:\n" );
  printf ( "\n" );
  printf ( "  * E is the eccentricity of the orbit, which should be between 0 and 1.0;\n" );
  printf ( "\n" );
  printf ( "  * M is the angle from the perihelion made by a fictitious planet traveling\n" );
  printf ( "    on a circular orbit centered at the sun, and traveling at a constant\n" );
  printf ( "    angular velocity equal to the average angular velocity of the true\n" );
  printf ( "    planet.  M is usually between 0 and 180 degrees, but can have any value.\n" );
  printf ( "\n" );
  printf ( "  For convenience, X and M are measured in degrees.\n" );
/*
  Problem parameters.
*/
  md = 24.851090;
  mr = md * r8_pi / 180.0;
  e = 0.1;

  printf ( "\n" );
  printf ( "  Given eccentricity E = %g\n", e );
  printf ( "  Given angle M = %g (degrees)\n", md );
  printf ( "                = %g (radians)\n", mr );
  printf ( "\n" );
  printf ( "  Given E and M, find corresponding X.\n"  );
/*
  Iteration parameters.
*/
  fx_tol = 1.0E-09;
  dx_tol = 1.0E-09;
  it = 0;;
  it_max = 30;
  job = 0;
  fx = 0.0;
/*
  Initial guess for interval.
*/
  ad = 0.0;
  bd = 180.0;

  ar = ad * r8_pi / 180.0;
  br = bd * r8_pi / 180.0;

  printf ( "\n" );
  printf ( "     I      X               FX              DX\n" );
  printf ( "\n" );

  for ( ; ; )
  {
    xr = bisection_rc ( &ar, &br, fx, &job );

    if ( job < 0 )
    {
      printf ( "\n" );
      printf ( "  Error return.\n" );
      break;
    }

    it = it + 1;

    fx = xr - mr - e * sin ( xr );

    if ( it <= 2 )
    {
      dx = fabs ( br - ar );
    }
    else
    {
      dx = 0.5 * fabs ( br - ar );
    }

    printf ( "  %4d  %14.6g  %14.6g  %14.6g\n", it, xr, fx, dx );

    if ( fabs ( fx ) <= fx_tol )
    {
      printf ( "\n" );
      printf ( "  Function is small.\n" );
      break;
    }

    if ( dx <= dx_tol )
    {
      printf ( "\n" );
      printf ( "  Interval is tiny.\n" );
      break;
    }

    if ( it_max <= it )
    {
      printf ( "\n" );
      printf ( "  Reached iteration limit.\n" );
      break;
    }

  }

  printf ( "\n" );
  printf ( "  In Radians:\n" );
  printf ( "\n" );
  fx = ar - mr - e * sin ( ar );
  printf ( "  A = %g, F(A) = %g\n", ar, fx );
  fx = xr - mr - e * sin ( xr );
  printf ( "  X = %g, F(X) = %g\n", xr, fx );
  fx = br - mr - e * sin ( br );
  printf ( "  B = %g, F(B) = %g\n", br, fx );

  ad = ar * 180.0 / r8_pi;
  xd = xr * 180.0 / r8_pi;
  bd = br * 180.0 / r8_pi;

  printf ( "\n" );
  printf ( "  In Degrees:\n" );
  printf ( "\n" );
  fx = ( ad - md ) * r8_pi / 180.0 - e * sin ( ad * r8_pi / 180.0 );
  printf ( "  A = %g, F(A) = %g\n", ad, fx );
  fx = ( xd - md ) * r8_pi / 180.0 - e * sin ( xd * r8_pi / 180.0 );
  printf ( "  X = %g, F(X) = %g\n", xd, fx );
  fx = ( bd - md ) * r8_pi / 180.0 - e * sin ( bd * r8_pi / 180.0);
  printf ( "  B = %g, F(B) = %g\n", bd, fx );

  return;
}

