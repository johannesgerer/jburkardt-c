# include <math.h>
# include <stdlib.h>
# include <stdio.h>
# include <string.h>
# include <time.h>

# include "cpv.h"

int main ( );
void cpv_test01 ( );
double f01 ( double t );
void cpv_test02 ( );
double f02 ( double t );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    CPV_PRB tests the CPV library.

  Location:

    http://people.sc.fsu.edu/~jburkardt/c_src/cauchy_principal_value/cpv_prb.c

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 April 2015

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "CPV_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the CPV library.\n" );

  cpv_test01 ( );
  cpv_test02 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "CPV_PRB\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void cpv_test01 ( )

/******************************************************************************/
/*
  Purpose:

    CPV_TEST01 seeks the CPV of Integral ( -1 <= t <= 1 ) exp(t) / t dt

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 April 2015

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double exact;
  int n;
  double value;

  printf ( "\n" );
  printf ( "CPV_TEST01:\n" );
  printf ( "  CPV of Integral ( -1 <= t <= 1 ) exp(t) / t dt\n" );

  printf ( "\n" );
  printf ( "   N           Estimate             Error\n" );
  printf ( "\n" );

  exact = 2.11450175075;
  a = -1.0;
  b = +1.0;
  for ( n = 2; n <= 8; n = n + 2 )
  {
    value = cpv ( f01, a, b, n );
    printf ( "  %2d  %24.16g  %14.6g\n", n, value, fabs ( value - exact ) );
  }

  return;
}
/******************************************************************************/

double f01 ( double t )

/******************************************************************************/
/*
  Purpose:

    F01 evaluates the integrand of Integral ( -1 <= t <= 1 ) exp(t) / t dt

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 April 2015

  Author:

    John Burkardt

  Parameters:

    Input, real T, the argument.

    Output, real VALUE, the value of the integrand.
*/
{
  double value;

  value = exp ( t );

  return value;
}
/******************************************************************************/

void cpv_test02 ( )

/******************************************************************************/
/*
  Purpose:

    CPV_TEST02 is another test.

 Discussion:

    We seek
      CPV ( Integral ( 1-delta <= t <= 1+delta ) 1/(1-t)^3 dt )
    which we must rewrite as
      CPV ( Integral ( 1-delta <= t <= 1+delta ) 1/(1+t+t^2) 1/(1-t) dt )
    so that our "integrand" is 1/(1+t+t^2).

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 April 2015

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double delta;
  double exact;
  int k;
  int n;
  double r1;
  double r2;
  double r3;
  double r4;
  double value;

  printf ( "\n" );
  printf ( "CPV_TEST02:\n" );
  printf ( "  Compute CPV ( Integral ( 1-delta <= t <= 1+delta ) 1/(1-t)^3 dt )\n" );
  printf ( "  Try this for delta = 1, 1/2, 1/4.\n" );
  printf ( "\n" );
  printf ( "   N          Estimate                  Exact                  Error\n" );
  delta = 1.0;
  for ( k = 1; k <= 3; k++ )
  {
    printf ( "\n" );
    r1 = pow (   delta + 1.5, 2 ) + 0.75;
    r2 = pow ( - delta + 1.5, 2 ) + 0.75;
    r3 = atan ( sqrt ( 0.75 ) / (   delta + 1.5 ) );
    r4 = atan ( sqrt ( 0.75 ) / ( - delta + 1.5 ) );
    exact = - log ( r1 / r2 ) / 6.0 + ( r3 - r4 ) / sqrt ( 3.0 );
    for ( n = 2; n <= 8; n = n + 2 )
    {
      a = 1.0 - delta;
      b = 1.0 + delta;
      value = cpv ( f02, a, b, n );
      printf ( "  %2d  %24.16g  %24.16g  %14.6g\n",
        n, value, exact, fabs ( exact - value ) );
    }
    delta = delta / 2.0;
  }

  return;
}
/******************************************************************************/

double f02 ( double t )

/******************************************************************************/
/*
  Purpose:

    F02: integrand of Integral ( 1-delta <= t <= 1+delta ) 1/(1-t^3) dt

  Discussion:

    1/(1-t^3) = 1/(1+t+t^2) * 1/(1-t)

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 April 2015

  Author:

    John Burkardt

  Parameters:

    Input, double T, the evaluation point.

    Output, double F02, the value of the integrand at T.
*/
{
  double value;

  value = 1.0 / ( 1.0 + t + t * t );

  return value;
}
