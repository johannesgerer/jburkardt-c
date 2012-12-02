# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>

# include "pce_ode_hermite.h"

int main ( );
void pce_ode_hermite_test01 ( );
void pce_ode_hermite_test02 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for PCE_ODE_HERMITE_TEST.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 September 2012

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "PCE_ODE_HERMITE_TEST:\n" );
  printf ( "  C version\n" );
  printf ( "  Test PCE_ODE_HERMITE.\n" );

  pce_ode_hermite_test01 ( );
  pce_ode_hermite_test02 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "PCE_ODE_HERMITE_TEST:\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void pce_ode_hermite_test01 ( )

/******************************************************************************/
/*
  Purpose:

    PCE_ODE_HERMITE_TEST01 runs a test problem with PCE_ODE_HERMITE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 September 2012

  Author:

    John Burkardt
*/
{
  double alpha_mu;
  double alpha_sigma;
  int i;
  int np = 4;
  int nt = 200;
  double *t;
  double tf;
  double ti;
  double *u;
  double *uex;
  double ui;

  t = ( double * ) malloc ( ( nt + 1 ) * sizeof ( double ) );
  u = ( double * ) malloc ( ( nt + 1 ) * ( np + 1 ) * sizeof ( double ) );
  uex = ( double * ) malloc ( ( nt + 1 ) * sizeof ( double ) );

  printf ( "\n" );
  printf ( "PCE_ODE_HERMITE_TEST01:\n" );
  printf ( "  Call PCE_ODE_HERMITE to compute a polynomial chaos expansion\n" );
  printf ( "  for the ODE:\n" );
  printf ( "\n" );
  printf ( "    u' = - alpha * u,\n" );
  printf ( "    u(0) = 1.\n" );

  ti = 0.0;
  tf = 2.0;
  ui = 1.0;
  alpha_mu = 0.0;
  alpha_sigma = 1.0;

  printf ( "\n" );
  printf ( "  Initial time         TI = %g\n", ti );
  printf ( "  Final time           TF = %g\n", tf );
  printf ( "  Number of time steps NT = %d\n", nt );
  printf ( "  Initial condition    UI = %g\n", ui );
  printf ( "  Expansion degree     NP = %d\n", np );
  printf ( "  E(ALPHA)       ALPHA_MU = %g\n", alpha_mu );
  printf ( "  STD(ALPHA)  ALPHA_SIGMA = %g\n", alpha_sigma );

  pce_ode_hermite ( ti, tf, nt, ui, np, alpha_mu, alpha_sigma, t, u );
/*
  Evaluate the exact expected value function.
*/
  for ( i = 0; i <= nt; i++ )
  {
    uex[i] = ui * exp ( t[i] * t[i] / 2.0 );
  }
/*
  Compare the first computed component against the exact expected value.
*/
  printf ( "\n" );
  printf ( "     i    T(i)        E(U(T(i)))      U(T(i),0)       Error\n" );
  printf ( "\n" );
  for ( i = 0; i <= nt; i = i + 10 )
  {
    printf ( "  %4d  %6g  %14g  %14g  %14g\n", 
      i, t[i], uex[i], u[i+0*(nt+1)], r8_abs ( uex[i] - u[i+0*(nt+1)] ) );
  }

  free ( t );
  free ( u );
  free ( uex );

  return;
}
/******************************************************************************/

void pce_ode_hermite_test02 ( )

/******************************************************************************/
/*
  Purpose:

    PCE_ODE_HERMITE_TEST02 looks at convergence behavior for a fixed time.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 September 2012

  Author:

    John Burkardt
*/
{
  double alpha_mu;
  double alpha_sigma;
  double ep[6];
  int i;
  int np;
  int nt = 2000;
  double *t;
  double tf;
  double ti;
  double *u;
  double uexf;
  double ui;

  t = ( double * ) malloc ( ( nt + 1 ) * sizeof ( double ) );

  printf ( "\n" );
  printf ( "PCE_ODE_HERMITE_TEST02:\n" );
  printf ( "  Examine convergence behavior as the PCE degree increases:\n" );
  printf ( "\n" );
  printf ( "    u' = - alpha * u,\n" );
  printf ( "    u(0) = 1.\n" );

  ti = 0.0;
  tf = 2.0;
  ui = 1.0;
  alpha_mu = 0.0;
  alpha_sigma = 1.0;

  printf ( "\n" );
  printf ( "  Initial time         TI = %g\n", ti );
  printf ( "  Final time           TF = %g\n", tf );
  printf ( "  Number of time steps NT = %d\n", nt );
  printf ( "  Initial condition    UI = %g\n", ui );
  printf ( "  E(ALPHA)       ALPHA_MU = %g\n", alpha_mu );
  printf ( "  STD(ALPHA)  ALPHA_SIGMA = %g\n", alpha_sigma );

  uexf = ui * exp ( tf * tf / 2.0 );

  for ( np = 0; np <= 5; np++ )
  {
    u = ( double * ) malloc ( (nt+1)*(np+1) * sizeof ( double ) );

    pce_ode_hermite ( ti, tf, nt, ui, np, alpha_mu, alpha_sigma, t, u );

    ep[np] = r8_abs ( uexf - u[nt+0*(nt+1)] );

    free ( u );
  }
/*
  Plot error in expected value as a function of the PCE degree.
*/
  printf ( "\n" );
  printf ( "    NP     Error(NP)     Log(Error(NP))\n" );
  printf ( "\n" );
  for ( np = 0; np <= 5; np++ )
  {
    printf ( "  %4d  %14g  %14g\n", np, ep[np], log ( ep[np] ) );
  }

  free ( t );

  return;
}
