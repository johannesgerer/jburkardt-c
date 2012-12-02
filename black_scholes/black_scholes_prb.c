# include <stdlib.h>
# include <stdio.h>

# include "black_scholes.h"

int main ( void );
void asset_path_test ( void );
void binomial_test ( void );
void bsf_test ( void );
void forward_test ( void );
void mc_test ( void );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    BLACK_SCHOLES_PRB tests BLACK_SCHOLES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    18 February 2012

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "BLACK_SCHOLES_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the BLACK_SCHOLES library.\n" );

  asset_path_test ( );
  binomial_test ( );
  bsf_test ( );
  forward_test ( );
  mc_test ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "BLACK_SCHOLES_PRB\n" );
  printf ( "  Normal end of execution.\n" );

  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void asset_path_test ( )

/******************************************************************************/
/*
  Purpose:

    ASSET_PATH_TEST tests ASSET_PATH.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    18 February 2012

  Author:

    John Burkardt
*/
{
  double mu;
  int n = 100;
  char output_filename[100] = "asset_path.txt";;
  double *s;
  double s0;
  int seed;
  double sigma;
  double t1;

  printf ( "\n" );
  printf ( "ASSET_PATH_TEST:\n" );
  printf ( "  Demonstrate the simulated of an asset price path.\n" );

  s0 = 2.0;
  mu = 0.1;
  sigma = 0.3;
  t1 = 1.0;
  seed = 123456789;

  printf ( "\n" );
  printf ( "  The asset price at time 0      S0    = %g\n", s0 );
  printf ( "  The asset expected growth rate MU    = %g\n", mu );
  printf ( "  The asset volatility           SIGMA = %g\n", sigma );
  printf ( "  The expiry date                T1    = %g\n", t1 );
  printf ( "  The number of time steps       N     = %d\n", n );
  printf ( "  The random number seed was     SEED  = %d\n", seed );

  s = asset_path ( s0, mu, sigma, t1, n, &seed );

  r8vec_print_part ( n + 1, s, 10, "  Partial results:" );

  r8vec_write ( output_filename, n + 1, s );

  printf ( "\n" );
  printf ( "  Full results written to \"%s\".\n", output_filename );

  free ( s );

  return;
}
/******************************************************************************/

void binomial_test ( )

/******************************************************************************/
/*
  Purpose:

    BINOMIAL_TEST tests BINOMIAL.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    18 February 2012

  Author:

    John Burkardt
*/
{
  double c;
  double e;
  int m;
  double r;
  double s0;
  double sigma;
  double t1;

  printf ( "\n" );
  printf ( "BINOMIAL_TEST:\n" );
  printf ( "  A demonstration of the binomial method\n" );
  printf ( "  for option valuation.\n" );

  s0 = 2.0;
  e = 1.0;
  r = 0.05;
  sigma = 0.25;
  t1 = 3.0;
  m = 256;

  printf ( "\n" );
  printf ( "  The asset price at time 0 S0    = %g\n", s0 );
  printf ( "  The exercise price        E     = %g\n", e );
  printf ( "  The interest rate         R     = %g\n", r );
  printf ( "  The asset volatility      SIGMA = %g\n", sigma );
  printf ( "  The expiry date           T1    = %g\n", t1 );
  printf ( "  The number of intervals   M     = %d\n", m );

  c = binomial ( s0, e, r, sigma, t1, m );

  printf ( "\n" );
  printf ( "  The option value is %g\n", c );

  return;
}
/******************************************************************************/

void bsf_test ( )

/******************************************************************************/
/*
  Purpose:

    BSF_TEST tests BSF.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    18 February 2012

  Author:

    John Burkardt
*/
{
  double c;
  double e;
  double r;
  double s0;
  double sigma;
  double t0;
  double t1;

  printf ( "\n" );
  printf ( "BSF_TEST:\n" );
  printf ( "  A demonstration of the Black-Scholes formula\n" );
  printf ( "  for option valuation.\n" );

  s0 = 2.0;
  t0 = 0.0;
  e = 1.0;
  r = 0.05;
  sigma = 0.25;
  t1 = 3.0;

  printf ( "\n" );
  printf ( "  The asset price at time T0 S0    = %g\n", s0 );
  printf ( "  The time                   T0    = %g\n", t0 );
  printf ( "  The exercise price         E     = %g\n", e );
  printf ( "  The interest rate          R     = %g\n", r );
  printf ( "  The asset volatility       SIGMA = %g\n", sigma );
  printf ( "  The expiry date            T1    = %g\n", t1 );

  c = bsf ( s0, t0, e, r, sigma, t1 );

  printf ( "\n" );
  printf ( "  The option value C = %g\n\n", c );

  return;
}
/******************************************************************************/

void forward_test ( )

/******************************************************************************/
/*
  Purpose:

    FORWARD_TEST tests FORWARD.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    18 February 2012

  Author:

    John Burkardt
*/
{
  double e;
  int i;
  int nt;
  int nx;
  double r;
  double s;
  double sigma;
  double smax;
  double smin;
  double t1;
  double *u;

  printf ( "\n" );
  printf ( "FORWARD_TEST:\n" );
  printf ( "  A demonstration of the forward difference method\n" );
  printf ( "  for option valuation.\n" );

  e = 4.0;
  r = 0.03;
  sigma = 0.50;
  t1 = 1.0;
  nx = 11;
  nt = 29;
  smax = 10.0;

  printf ( "\n" );
  printf ( "  The exercise price        E =     %g\n", e );
  printf ( "  The interest rate         R =     %g\n", r );
  printf ( "  The asset volatility      SIGMA = %g\n", sigma );
  printf ( "  The expiry date           T1 =    %g\n", t1 );
  printf ( "  The number of space steps NX =    %d\n", nx );
  printf ( "  The number of time steps  NT =    %d\n", nt );
  printf ( "  The value of              SMAX =  %g\n", smax );

  u = forward ( e, r, sigma, t1, nx, nt, smax );

  printf ( "\n" );
  printf ( "         Initial          Option\n" );
  printf ( "           Value           Value\n" ); 
  printf ( "\n" );

  smin = 0.0;
  for ( i = 0; i < nx - 1; i++ )
  {
    s = ( ( nx - i - 2 ) * smin +  ( i + 1 ) * smax ) / ( double ) ( nx - 1 );
    printf ( "  %14g  %14g\n", s, u[i+nt*(nx-1)] );
  }

  free ( u );

  return;
}
/******************************************************************************/

void mc_test ( )

/******************************************************************************/
/*
  Purpose:

    MC_TEST tests MC.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    18 February 2012

  Author:

    John Burkardt
*/
{
  double *conf;
  double e;
  int m;
  double r;
  double s0;
  int seed;
  double sigma;
  double t1;

  printf ( "\n" );
  printf ( "MC_TEST:\n" );
  printf ( "  A demonstration of the Monte Carlo method\n" );
  printf ( "  for option valuation.\n" );

  s0 = 2.0;
  e = 1.0;
  r = 0.05;
  sigma = 0.25;
  t1 = 3.0;
  m = 1000000;
  seed = 123456789;

  printf ( "\n" );
  printf ( "  The asset price at time 0, S0    = %g\n", s0 );
  printf ( "  The exercise price         E     = %g\n", e );
  printf ( "  The interest rate          R     = %g\n", r );
  printf ( "  The asset volatility       SIGMA = %g\n", sigma );
  printf ( "  The expiry date            T1    = %g\n", t1 );
  printf ( "  The number of simulations  M     = %d\n", m );

  conf = mc ( s0, e, r, sigma, t1, m, &seed );

  printf ( "\n" );
  printf ( "  The confidence interval is [%g,%g].\n", conf[0], conf[1] );

  return;
}
