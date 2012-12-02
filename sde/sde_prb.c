# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <string.h>

# include "sde.h"

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

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for SDE_PRB.

  Discussion:

    SDE_PRB demonstrates the use of SDE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    27 September 2012

  Author:

    John Burkardt
*/
{
  timestamp ( );

  printf ( "\n" );
  printf ( "SDE_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the SDE library.\n" );

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
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "SDE_PRB\n" );
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

    TEST01 tests BPATH.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    23 September 2012

  Author:

    John Burkardt
*/
{
  int n = 500;
  int seed;
  double *w;

  printf ( "\n" );
  printf ( "TEST01:\n" );
  printf ( "  BPATH generates a sample Brownian motion path.\n" );

  seed = 123456789;

  w = bpath ( &seed, n );

  bpath_gnuplot ( n, w );

  free ( w );

  return;
}
/******************************************************************************/

void test02 ( )

/******************************************************************************/
/*
  Purpose:

    TEST02 tests BPATH_AVERAGE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    23 September 2012

  Author:

    John Burkardt
*/
{
  double error;
  int m = 1000;
  int n = 500;
  int seed;
  double *u;
  double *umean;

  printf ( "\n" );
  printf ( "TEST02:\n" );
  printf ( "  BPATH_AVERAGE generates many Brownian paths\n" );
  printf ( "  and averages them.\n" );

  seed = 123456789;
  u = ( double * ) malloc ( m * ( n + 1 ) * sizeof ( double ) );
  umean = ( double * ) malloc ( ( n + 1 ) * sizeof ( double ) );

  bpath_average ( &seed, m, n, u, umean, &error );

  bpath_average_gnuplot ( m, n, u, umean );

  free ( u );
  free ( umean );

  return;
}
/******************************************************************************/

void test03 ( )

/******************************************************************************/
/*
  Purpose:

    TEST03 tests CHAIN.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    23 September 2012

  Author:

    John Burkardt
*/
{
  double diff;
  int n = 200;
  int seed;
  double *vem;
  double *xem;

  printf ( "\n" );
  printf ( "TEST03:\n" );
  printf ( "  CHAIN solves a stochastic differential equation for\n" );
  printf ( "  a function of a stochastic variable X.\n" );
  printf ( "  We can solve for X(t), and then evaluate V(X(t)).\n" );
  printf ( "  Or, we can apply the stochastic chain rule to derive an\n" );
  printf ( "  an SDE for V, and solve that.\n" );

  seed = 123456789;
  xem = ( double * ) malloc ( ( n + 1 ) * sizeof ( double ) );
  vem = ( double * ) malloc ( ( n + 1 ) * sizeof ( double ) );
  chain ( &seed, n, xem, vem, &diff );

  printf ( "\n" );
  printf ( "  Maximum | Sqrt(X) - V | = %g\n", diff );

  chain_gnuplot ( n, xem, vem );

  free ( vem );
  free ( xem );

  return;
}
/******************************************************************************/

void test04 ( )

/******************************************************************************/
/*
  Purpose:

    TEST04 tests EM.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    25 September 2012

  Author:

    John Burkardt
*/
{
  double diff;
  int n = 256;
  int seed;
  double *t;
  double *t2;
  double *xem;
  double *xtrue;

  printf ( "\n" );
  printf ( "TEST04:\n" );
  printf ( "  EM solves a stochastic differential equation\n" );
  printf ( "  using the Euler-Maruyama method.\n" );

  seed = 123456789;

  t = ( double * ) malloc ( ( n + 1 ) * sizeof ( double ) );
  t2 = ( double * ) malloc ( ( n / 4 + 1 ) * sizeof ( double ) );
  xem = ( double * ) malloc ( ( n / 4 + 1 ) * sizeof ( double ) );
  xtrue = ( double * ) malloc ( ( n + 1 ) * sizeof ( double ) );

  em ( &seed, n, t, xtrue, t2, xem, &diff );

  printf ( "\n" );
  printf ( "  | Exact X(T) - EM X(T) | = %14.6g\n", diff );

  em_gnuplot ( n, t, xtrue, t2, xem );

  free ( t );
  free ( t2 );
  free ( xem );
  free ( xtrue );

  return;
}
/******************************************************************************/

void test05 ( )

/******************************************************************************/
/*
  Purpose:

    TEST05 tests EMSTRONG.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    25 September 2012

  Author:

    John Burkardt
*/
{
  double *dtvals;
  int m = 100;
  int n = 512;
  int p_max = 6;
  int seed;
  double *xerr;

  dtvals = ( double * ) malloc ( p_max * sizeof ( double ) );
  xerr = ( double * ) malloc ( p_max * sizeof ( double ) );

  printf ( "\n" );
  printf ( "TEST05:\n" );
  printf ( "  EMSTRONG investigates the strong convergence\n" );
  printf ( "  of the Euler-Maruyama method.\n" );

  seed = 123456789;

  emstrong ( &seed, m, n, p_max, dtvals, xerr );

  emstrong_gnuplot ( p_max, dtvals, xerr );

  free ( dtvals );
  free ( xerr );

  return;
}
/******************************************************************************/

void test06 ( )

/******************************************************************************/
/*
  Purpose:

    TEST06 tests EMWEAK.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    25 September 2012

  Author:

    John Burkardt
*/
{
  double *dtvals;
  int m = 50000;
  int method;
  int p_max = 5;
  int seed;
  double *xerr;

  dtvals = ( double * ) malloc ( p_max * sizeof ( double ) );
  xerr = ( double * ) malloc ( p_max * sizeof ( double ) );

  printf ( "\n" );
  printf ( "TEST06:\n" );
  printf ( "  EMWEAK investigates the weak convergence\n" );
  printf ( "  of the Euler-Maruyama method.\n" );

  seed = 123456789;
  method = 0;

  emweak ( &seed, method, m, p_max, dtvals, xerr );

  emweak_gnuplot ( p_max, dtvals, xerr, method );

  seed = 123456789;
  method = 1;

  emweak ( &seed, method, m, p_max, dtvals, xerr );

  emweak_gnuplot ( p_max, dtvals, xerr, method );

  free ( dtvals );
  free ( xerr );

  return;
}
/******************************************************************************/

void test07 ( )

/******************************************************************************/
/*
  Purpose:

    TEST07 tests MILSTRONG.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    27 September 2012

  Author:

    John Burkardt
*/
{
  double *dtvals;
  int p_max = 4;
  int seed;
  double *xerr;

  dtvals = ( double * ) malloc ( p_max * sizeof ( double ) );
  xerr = ( double * ) malloc ( p_max * sizeof ( double ) );

  printf ( "\n" );
  printf ( "TEST07:\n" );
  printf ( "  MILSTRONG investigates the strong convergence\n" );
  printf ( "  of the Milstein method.\n" );

  seed = 123456789;

  milstrong ( &seed, p_max, dtvals, xerr );

  milstrong_gnuplot ( p_max, dtvals, xerr );

  free ( dtvals );
  free ( xerr );

  return;
}
/******************************************************************************/

void test08 ( )

/******************************************************************************/
/*
  Purpose:

    TEST08 tests STAB_ASYMPTOTIC.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    27 September 2012

  Author:

    John Burkardt
*/
{
  int n = 1000;
  int p_max = 3;
  int seed;

  printf ( "\n" );
  printf ( "TEST08:\n" );
  printf ( "  STAB_ASYMPTOTIC investigates the asymptotic\n" );
  printf ( "  stability of the Euler-Maruyama method.\n" );
  printf ( "\n" );
  printf ( "  For technical reasons, the plotting is done\n" );
  printf ( "  in the same routine as the computations.\n" );

  seed = 123456789;

  stab_asymptotic ( &seed, n, p_max );

  return;
}

/******************************************************************************/

void test09 ( )

/******************************************************************************/
/*
  Purpose:

    TEST09 tests STAB_MEANSQUARE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    25 September 2012

  Author:

    John Burkardt
*/
{
  int seed;

  printf ( "\n" );
  printf ( "TEST09:\n" );
  printf ( "  STAB_MEANSQUARE investigates the mean square\n" );
  printf ( "  stability of the Euler-Maruyama method.\n" );
  printf ( "\n" );
  printf ( "  For technical reasons, the plotting is done\n" );
  printf ( "  in the same routine as the computations.\n" );

  seed = 123456789;

  stab_meansquare ( &seed );

  return;
}
/******************************************************************************/

void test10 ( )

/******************************************************************************/
/*
  Purpose:

    TEST10 tests STOCHASTIC_INTEGRAL_ITO.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    24 September 2012

  Author:

    John Burkardt
*/
{
  double error;
  double estimate;
  double exact;
  int i;
  int n;
  int seed;

  printf ( "\n" );
  printf ( "TEST10:\n" );
  printf ( "  Estimate the Ito integral of W(t) dW over [0,1].\n" );
  printf ( "\n" );
  printf ( "                                                 Abs          Rel\n" );
  printf ( "         N        Exact        Estimate          Error        Error\n" );
  printf ( "\n" );

  n = 100;
  seed = 123456789;

  for ( i = 1; i <= 7; i++ )
  {
    stochastic_integral_ito ( n, &seed, &estimate, &exact, &error );

    printf ( "  %8d  %16.8g  %16.8g  %10.2g  %10.2g\n",
      n, exact, estimate, error, error / exact );

    n = n * 4;
  }
  return;
}
/******************************************************************************/

void test11 ( )

/******************************************************************************/
/*
  Purpose:

    TEST11 tests STOCHASTIC_INTEGRAL_STRAT.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    24 September 2012

  Author:

    John Burkardt
*/
{
  double error;
  double estimate;
  double exact;
  int i;
  int n;
  int seed;

  printf ( "\n" );
  printf ( "TEST11:\n" );
  printf ( "  Estimate the Stratonovich integral of W(t) dW over [0,1].\n" );
  printf ( "\n" );
  printf ( "                                                 Abs          Rel\n" );
  printf ( "         N        Exact        Estimate          Error        Error\n" );
  printf ( "\n" );

  n = 100;
  seed = 123456789;

  for ( i = 1; i <= 7; i++ )
  {
    stochastic_integral_strat ( n, &seed, &estimate, &exact, &error );

    printf ( "  %8d  %16.8g  %16.8g  %10.2g  %10.2g\n",
      n, exact, estimate, error, error / exact );

    n = n * 4;
  }
  return;
}
