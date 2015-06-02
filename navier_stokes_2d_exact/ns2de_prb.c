# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <string.h>
# include <time.h>

# include "ns2de.h"

int main ( );
void uvp_taylor_test ( );
void uvp_taylor_test2 ( );
void rhs_taylor_test ( );
void resid_taylor_test ( );
void gnuplot_taylor_test ( );
void parameter_taylor_test ( );
void uvp_spiral_test ( );
void uvp_spiral_test2 ( );
void rhs_spiral_test ( );
void resid_spiral_test ( );
void gnuplot_spiral_test ( );
void parameter_spiral_test ( );
void uvp_lucas_test ( );
void uvp_lucas_test2 ( );
void rhs_lucas_test ( );
void resid_lucas_test ( );
void gnuplot_lucas_test ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    NS2DE_TEST tests the NS2DE library.

  Location:

    http://people.sc.fsu.edu/~jburkardt/c_src/navier_stokes_2d_exact/ns2de_prb.c

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    06 March 2015

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "NS2DE_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the NS2DE library.\n" );
/*
  Taylor Vortex Flow.
*/
  uvp_taylor_test ( );
  uvp_taylor_test2 ( );
  rhs_taylor_test ( );
  resid_taylor_test ( );
  gnuplot_taylor_test ( );
  parameter_taylor_test ( );
/*
  Spiral Flow.
*/
  uvp_spiral_test ( );
  uvp_spiral_test2 ( );
  rhs_spiral_test ( );
  resid_spiral_test ( );
  gnuplot_spiral_test ( );
  parameter_spiral_test ( );
/*
  Lucas Bystricky Flow.
*/
  uvp_lucas_test ( );
  uvp_lucas_test2 ( );
  rhs_lucas_test ( );
  resid_lucas_test ( );
  gnuplot_lucas_test ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "NS2DE_TEST\n" );
  printf ( "  Normal end of execution.\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void uvp_taylor_test ( )

/******************************************************************************/
/*
  Purpose:

    UVP_TAYLOR_TEST samples the solution at the initial time.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 January 2015

  Author:

    John Burkardt
*/
{
  int n;
  double nu;
  double *p;
  double rho;
  int seed;
  double t;
  double *u;
  double *v;
  double *x;
  double r8_hi;
  double r8_lo;
  double *y;

  nu = 1.0;
  rho = 1.0;

  printf ( "\n" );
  printf ( "UVP_TAYLOR_TEST\n" );
  printf ( "  Taylor Vortex Flow:\n" );
  printf ( "  Estimate the range of velocity and pressure\n" );
  printf ( "  at the initial time T = 0, using a region that is\n" );
  printf ( "  the square centered at (1.5,1.5) with 'radius' 1.0,\n" );
  printf ( "  Kinematic viscosity NU = %g\n", nu );
  printf ( "  Fluid density RHO = %g\n", rho );

  n = 1000;

  p = ( double * ) malloc ( n * sizeof ( double ) );
  u = ( double * ) malloc ( n * sizeof ( double ) );
  v = ( double * ) malloc ( n * sizeof ( double ) );

  r8_lo = 0.5;
  r8_hi = 2.5;
  seed = 123456789;

  x = r8vec_uniform_ab_new ( n, r8_lo, r8_hi, &seed );
  y = r8vec_uniform_ab_new ( n, r8_lo, r8_hi, &seed );
  t = 0.0;

  uvp_taylor ( nu, rho, n, x, y, t, u, v, p );

  printf ( "\n" );
  printf ( "           Minimum       Maximum\n" );
  printf ( "\n" );
  printf ( "  U:  %14.6g  %14.6g\n", r8vec_min ( n, u ), r8vec_max ( n, u ) );
  printf ( "  V:  %14.6g  %14.6g\n", r8vec_min ( n, v ), r8vec_max ( n, v ) );
  printf ( "  P:  %14.6g  %14.6g\n", r8vec_min ( n, p ), r8vec_max ( n, p ) );

  free ( p );
  free ( u );
  free ( v );
  free ( x );
  free ( y );

  return;
}
/******************************************************************************/

void uvp_taylor_test2 ( )

/******************************************************************************/
/*
  Purpose:

    UVP_TAYLOR_TEST2 samples the solution on the boundary at the initial time.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    06 March 2015

  Author:

    John Burkardt
*/
{
  int i;
  int n;
  double nu;
  double *p;
  double r8_hi;
  double r8_lo;
  double rho;
  int seed;
  double t;
  double *u;
  double *v;
  double *x;
  double *y;

  nu = 1.0;
  rho = 1.0;
  t = 0.0;

  printf ( "\n" );
  printf ( "UVP_TAYLOR_TEST2\n" );
  printf ( "  Taylor Vortex Flow:\n" );
  printf ( "  Estimate the range of velocity and pressure\n" );
  printf ( "  on the boundary,\n" );
  printf ( "  at the initial time T = 0, using a region that is\n" );
  printf ( "  the square centered at (1.5,1.5) with 'radius' 1.0,\n" );
  printf ( "  Kinematic viscosity NU = %g\n", nu );
  printf ( "  Fluid density RHO = %g\n", rho );

  n = 400;

  p = ( double * ) malloc ( n * sizeof ( double ) );
  u = ( double * ) malloc ( n * sizeof ( double ) );
  v = ( double * ) malloc ( n * sizeof ( double ) );
  x = ( double * ) malloc ( n * sizeof ( double ) );
  y = ( double * ) malloc ( n * sizeof ( double ) );

  r8_lo = 0.5;
  r8_hi = 2.5;

  r8vec_linspace ( 100, r8_lo, r8_hi, x );
  for ( i = 0; i < 100; i++ )
  {
    y[i] = r8_lo;
  }

  for ( i = 100; i < 200; i++ )
  {
    x[i] = r8_hi;
  }
  r8vec_linspace ( 100, r8_lo, r8_hi, y+100 );

  r8vec_linspace ( 100, r8_hi, r8_lo, x+200 );
  for ( i = 200; i < 300; i++ )
  {
    y[i] = r8_hi;
  }

  for ( i = 300; i < 400; i++ )
  {
    x[i] = r8_lo;
  }
  r8vec_linspace ( 100, r8_lo, r8_hi, y+300 );

  uvp_taylor ( nu, rho, n, x, y, t, u, v, p );

  printf ( "\n" );
  printf ( "           Minimum       Maximum\n" );
  printf ( "\n" );
  printf ( "  U:  %14.6g  %14.6g\n", r8vec_min ( n, u ), r8vec_max ( n, u ) );
  printf ( "  V:  %14.6g  %14.6g\n", r8vec_min ( n, v ), r8vec_max ( n, v ) );
  printf ( "  P:  %14.6g  %14.6g\n", r8vec_min ( n, p ), r8vec_max ( n, p ) );

  free ( p );
  free ( u );
  free ( v );
  free ( x );
  free ( y );

  return;
}
/******************************************************************************/

void rhs_taylor_test ( )

/******************************************************************************/
/*
  Purpose:

    RHS_TAYLOR_TEST samples the right hand sides at the initial time.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 January 2015

  Author:

    John Burkardt
*/
{
  double *f;
  double *g;
  double *h;
  int n;
  double nu;
  double rho;
  int seed;
  double t;
  double *x;
  double r8_hi;
  double r8_lo;
  double *y;

  nu = 1.0;
  rho = 1.0;

  printf ( "\n" );
  printf ( "RHS_TAYLOR_TEST\n" );
  printf ( "  Taylor Vortex Flow:\n" );
  printf ( "  Sample the Navier-Stokes right hand sides\n" );
  printf ( "  at the initial time T = 0, using a region that is\n" );
  printf ( "  the square centered at (1.5,1.5) with 'radius' 1.0,\n" );
  printf ( "  Kinematic viscosity NU = %g\n", nu );
  printf ( "  Fluid density RHO = %g\n", rho );

  n = 1000;

  f = ( double * ) malloc ( n * sizeof ( double ) );
  g = ( double * ) malloc ( n * sizeof ( double ) );
  h = ( double * ) malloc ( n * sizeof ( double ) );

  r8_lo = 0.5;
  r8_hi = +2.5;
  seed = 123456789;

  x = r8vec_uniform_ab_new ( n, r8_lo, r8_hi, &seed );
  y = r8vec_uniform_ab_new ( n, r8_lo, r8_hi, &seed );
  t = 0.0;

  rhs_taylor ( nu, rho, n, x, y, t, f, g, h );

  printf ( "\n" );
  printf ( "           Minimum       Maximum\n" );
  printf ( "\n" );
  printf ( "  F:  %14.6g  %14.6g\n", r8vec_min ( n, f ), r8vec_max ( n, f ) );
  printf ( "  G:  %14.6g  %14.6g\n", r8vec_min ( n, g ), r8vec_max ( n, g ) );
  printf ( "  H:  %14.6g  %14.6g\n", r8vec_min ( n, h ), r8vec_max ( n, h ) );

  free ( f );
  free ( g );
  free ( h );
  free ( x );
  free ( y );

  return;
}
/******************************************************************************/

void resid_taylor_test ( )

/******************************************************************************/
/*
  Purpose:

    RESID_TAYLOR_TEST samples the residual at the initial time.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 January 2015

  Author:

    John Burkardt
*/
{
  int n;
  double nu;
  double *pr;
  double rho;
  int seed;
  double t;
  double *ur;
  double *vr;
  double *x;
  double r8_hi;
  double r8_lo;
  double *y;

  nu = 1.0;
  rho = 1.0;

  printf ( "\n" );
  printf ( "RESID_TAYLOR_TEST\n" );
  printf ( "  Taylor Vortex Flow:\n" );
  printf ( "  Sample the Navier-Stokes residuals\n" );
  printf ( "  at the initial time T = 0, using a region that is\n" );
  printf ( "  the square centered at (1.5,1.5) with 'radius' 1.0,\n" );
  printf ( "  Kinematic viscosity NU = %g\n", nu );
  printf ( "  Fluid density RHO = %g\n", rho );

  n = 1000;

  pr = ( double * ) malloc ( n * sizeof ( double ) );
  ur = ( double * ) malloc ( n * sizeof ( double ) );
  vr = ( double * ) malloc ( n * sizeof ( double ) );

  r8_lo = 0.5;
  r8_hi = +2.5;
  seed = 123456789;

  x = r8vec_uniform_ab_new ( n, r8_lo, r8_hi, &seed );
  y = r8vec_uniform_ab_new ( n, r8_lo, r8_hi, &seed );
  t = 0.0;

  resid_taylor ( nu, rho, n, x, y, t, ur, vr, pr );

  printf ( "\n" );
  printf ( "           Minimum       Maximum\n" );
  printf ( "\n" );
  printf ( "  Ur:  %14.6g  %14.6g\n", r8vec_amin ( n, ur ), r8vec_amax ( n, ur ) );
  printf ( "  Vr:  %14.6g  %14.6g\n", r8vec_amin ( n, vr ), r8vec_amax ( n, vr ) );
  printf ( "  Pr:  %14.6g  %14.6g\n", r8vec_amin ( n, pr ), r8vec_amax ( n, pr ) );

  free ( pr );
  free ( ur );
  free ( vr );
  free ( x );
  free ( y );

  return;
}
/******************************************************************************/

void gnuplot_taylor_test ( )

/******************************************************************************/
/*
  Purpose:

    GNUPLOT_TAYLOR_TEST generates a field on a regular grid and plots it.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    21 January 2015

  Author:

    John Burkardt
*/
{
  char header[255];
  int n;
  double nu;
  double *p;
  double rho;
  double s;
  int seed;
  double t;
  double *u;
  double *v;
  double *x;
  double x_hi;
  double x_lo;
  int x_num = 21;
  double *y;
  double y_hi;
  double y_lo;
  int y_num = 21;

  printf ( "\n" );
  printf ( "GNUPLOT_TAYLOR_TEST:\n" );
  printf ( "  Taylor Vortex Flow:\n" );
  printf ( "  Generate a velocity field on a regular grid.\n" );
  printf ( "  Store in GNUPLOT data and command files.\n" );

  x_lo = 0.5;
  x_hi = 2.5;

  y_lo = 0.5;
  y_hi = 2.5;

  x = ( double * ) malloc ( x_num * y_num * sizeof ( double ) );
  y = ( double * ) malloc ( x_num * y_num * sizeof ( double ) );

  grid_2d ( x_num, x_lo, x_hi, y_num, y_lo, y_hi, x, y );

  nu = 1.0;
  rho = 1.0;
  n = x_num * y_num;
  t = 0.0;

  u = ( double * ) malloc ( x_num * y_num * sizeof ( double ) );
  v = ( double * ) malloc ( x_num * y_num * sizeof ( double ) );
  p = ( double * ) malloc ( x_num * y_num * sizeof ( double ) );

  uvp_taylor ( nu, rho, n, x, y, t, u, v, p );

  strcpy ( header, "taylor" );
  s = 0.10;
  ns2de_gnuplot ( header, n, x, y, u, v, s );

  free ( p );
  free ( u );
  free ( v );
  free ( x );
  free ( y );

  return;
}
/******************************************************************************/

void parameter_taylor_test ( )

/******************************************************************************/
/*
  Purpose:

    PARAMETER_TAYLOR_TEST monitors solution norms for various values of NU, RHO.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    22 January 2015

  Author:

    John Burkardt
*/
{
  int i;
  int j;
  int k;
  int n;
  double nu;
  double *p;
  double p_norm;
  double rho;
  int seed;
  double t;
  double *u;
  double u_norm;
  double *v;
  double v_norm;
  double *x;
  double r8_hi;
  double r8_lo;
  double *y;

  printf ( "\n" );
  printf ( "PARAMETER_TAYLOR_TEST\n" );
  printf ( "  Taylor Vortex Flow:\n" );
  printf ( "  Monitor solution norms over time for various\n" );
  printf ( "  values of NU, RHO.\n" );

  n = 1000;

  u = ( double * ) malloc ( n * sizeof ( double ) );
  v = ( double * ) malloc ( n * sizeof ( double ) );
  p = ( double * ) malloc ( n * sizeof ( double ) );

  r8_lo = 0.5;
  r8_hi = +2.5;
  seed = 123456789;

  x = r8vec_uniform_ab_new ( n, r8_lo, r8_hi, &seed );
  y = r8vec_uniform_ab_new ( n, r8_lo, r8_hi, &seed );
/*
  Vary RHO.
*/
  printf ( "\n" );
  printf ( "  RHO affects the pressure scaling.\n" );
  printf ( "\n" );
  printf ( "     RHO         NU           T     ||U||       ||V||       ||P||\n" );
  printf ( "\n" );

  nu = 1.0;
  rho = 1.0;

  for ( j = 1; j <= 3; j++ )
  {
    for ( k = 0; k <= 5; k++ )
    {
      t = ( double ) ( k ) / 5.0;

      uvp_taylor ( nu, rho, n, x, y, t, u, v, p );

      u_norm = r8vec_norm_l2 ( n, u ) / ( double ) ( n );
      v_norm = r8vec_norm_l2 ( n, v ) / ( double ) ( n );
      p_norm = r8vec_norm_l2 ( n, p ) / ( double ) ( n );

      printf ( "  %10.4g  %10.4g  %8.4f  %10.4g  %10.4g  %10.4g\n",
        rho, nu, t, u_norm, v_norm, p_norm );
    }
    printf ( "\n" );
    rho = rho / 100.0;
  }
/*
  Vary NU.
*/
  printf ( "\n" );
  printf ( "  NU affects the time scaling.\n" );
  printf ( "\n" );
  printf ( "     RHO         NU           T     ||U||       ||V||       ||P||\n" );
  printf ( "\n" );

  nu = 1.0;
  rho = 1.0;

  for ( i = 1; i <= 4; i++ )
  {
    for ( k = 0; k <= 5; k++ )
    {
      t = ( double ) ( k ) / 5.0;

      uvp_taylor ( nu, rho, n, x, y, t, u, v, p );

      u_norm = r8vec_norm_l2 ( n, u ) / ( double ) ( n );
      v_norm = r8vec_norm_l2 ( n, v ) / ( double ) ( n );
      p_norm = r8vec_norm_l2 ( n, p ) / ( double ) ( n );

      printf ( "  %10.4g  %10.4g  %8.4f  %10.4g  %10.4g  %10.4g\n",
        rho, nu, t, u_norm, v_norm, p_norm );
    }

    printf ( "\n" );

    nu = nu / 10.0;
  }

  free ( p );
  free ( u );
  free ( v );
  free ( x );
  free ( y );

  return;
}
/******************************************************************************/

void uvp_spiral_test ( )

/******************************************************************************/
/*
  Purpose:

    UVP_SPIRAL_TEST samples the spiral flow solution at the initial time.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 January 2015

  Author:

    John Burkardt
*/
{
  int n;
  double nu;
  double *p;
  double rho;
  int seed;
  double t;
  double *u;
  double *v;
  double *x;
  double r8_hi;
  double r8_lo;
  double *y;

  nu = 1.0;
  rho = 1.0;

  printf ( "\n" );
  printf ( "UVP_SPIRAL_TEST\n" );
  printf ( "  Spiral flow:\n" );
  printf ( "  Estimate the range of velocity and pressure\n" );
  printf ( "  at the initial time T = 0, using a region that is\n" );
  printf ( "  the unit square.\n" );
  printf ( "  Kinematic viscosity NU = %g\n", nu );
  printf ( "  Fluid density RHO = %g\n", rho );

  n = 1000;

  p = ( double * ) malloc ( n * sizeof ( double ) );
  u = ( double * ) malloc ( n * sizeof ( double ) );
  v = ( double * ) malloc ( n * sizeof ( double ) );

  r8_lo = 0.0;
  r8_hi = 1.0;
  seed = 123456789;

  x = r8vec_uniform_ab_new ( n, r8_lo, r8_hi, &seed );
  y = r8vec_uniform_ab_new ( n, r8_lo, r8_hi, &seed );
  t = 0.0;

  uvp_spiral ( nu, rho, n, x, y, t, u, v, p );

  printf ( "\n" );
  printf ( "           Minimum       Maximum\n" );
  printf ( "\n" );
  printf ( "  U:  %14.6g  %14.6g\n", r8vec_min ( n, u ), r8vec_max ( n, u ) );
  printf ( "  V:  %14.6g  %14.6g\n", r8vec_min ( n, v ), r8vec_max ( n, v ) );
  printf ( "  P:  %14.6g  %14.6g\n", r8vec_min ( n, p ), r8vec_max ( n, p ) );

  free ( p );
  free ( u );
  free ( v );
  free ( x );
  free ( y );

  return;
}
/******************************************************************************/

void uvp_spiral_test2 ( )

/******************************************************************************/
/*
  Purpose:

    UVP_SPIRAL_TEST2 samples the solution on the boundary at the initial time.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    06 March 2015

  Author:

    John Burkardt
*/
{
  int i;
  int n;
  double nu;
  double *p;
  double rho;
  int seed;
  double t;
  double *u;
  double *v;
  double *x;
  double r8_hi;
  double r8_lo;
  double *y;

  nu = 1.0;
  rho = 1.0;
  t = 0.0;

  printf ( "\n" );
  printf ( "UVP_SPIRAL_TEST2\n" );
  printf ( "  Spiral Flow:\n" );
  printf ( "  Estimate the range of velocity and pressure\n" );
  printf ( "  on the boundary,\n" );
  printf ( "  at the initial time T = 0, using the unit square.\n" );
  printf ( "  Kinematic viscosity NU = %g\n", nu );
  printf ( "  Fluid density RHO = %g\n", rho );

  n = 400;

  p = ( double * ) malloc ( n * sizeof ( double ) );
  u = ( double * ) malloc ( n * sizeof ( double ) );
  v = ( double * ) malloc ( n * sizeof ( double ) );
  x = ( double * ) malloc ( n * sizeof ( double ) );
  y = ( double * ) malloc ( n * sizeof ( double ) );

  r8_lo = 0.0;
  r8_hi = 1.0;

  r8vec_linspace ( 100, r8_lo, r8_hi, x );
  for ( i = 0; i < 100; i++ )
  {
    y[i] = r8_lo;
  }

  for ( i = 100; i < 200; i++ )
  {
    x[i] = r8_hi;
  }
  r8vec_linspace ( 100, r8_lo, r8_hi, y+100 );

  r8vec_linspace ( 100, r8_hi, r8_lo, x+200 );
  for ( i = 200; i < 300; i++ )
  {
    y[i] = r8_hi;
  }

  for ( i = 300; i < 400; i++ )
  {
    x[i] = r8_lo;
  }
  r8vec_linspace ( 100, r8_lo, r8_hi, y+300 );

  uvp_spiral ( nu, rho, n, x, y, t, u, v, p );

  printf ( "\n" );
  printf ( "           Minimum       Maximum\n" );
  printf ( "\n" );
  printf ( "  U:  %14.6g  %14.6g\n", r8vec_min ( n, u ), r8vec_max ( n, u ) );
  printf ( "  V:  %14.6g  %14.6g\n", r8vec_min ( n, v ), r8vec_max ( n, v ) );
  printf ( "  P:  %14.6g  %14.6g\n", r8vec_min ( n, p ), r8vec_max ( n, p ) );

  free ( p );
  free ( u );
  free ( v );
  free ( x );
  free ( y );

  return;
}
/******************************************************************************/

void rhs_spiral_test ( )

/******************************************************************************/
/*
  Purpose:

    RHS_SPIRAL_TEST samples the right hand sides at the initial time.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 January 2015

  Author:

    John Burkardt
*/
{
  double *f;
  double *g;
  double *h;
  int n;
  double nu;
  double rho;
  int seed;
  double t;
  double *x;
  double r8_hi;
  double r8_lo;
  double *y;

  nu = 1.0;
  rho = 1.0;

  printf ( "\n" );
  printf ( "RHS_SPIRAL_TEST\n" );
  printf ( "  Spiral Flow:\n" );
  printf ( "  Sample the Navier-Stokes right hand sides\n" );
  printf ( "  at the initial time T = 0, using the unit square.\n" );
  printf ( "  Kinematic viscosity NU = %g\n", nu );
  printf ( "  Fluid density RHO = %g\n", rho );

  n = 1000;

  f = ( double * ) malloc ( n * sizeof ( double ) );
  g = ( double * ) malloc ( n * sizeof ( double ) );
  h = ( double * ) malloc ( n * sizeof ( double ) );

  r8_lo = 0.0;
  r8_hi = 1.0;
  seed = 123456789;

  x = r8vec_uniform_ab_new ( n, r8_lo, r8_hi, &seed );
  y = r8vec_uniform_ab_new ( n, r8_lo, r8_hi, &seed );
  t = 0.0;

  rhs_spiral ( nu, rho, n, x, y, t, f, g, h );

  printf ( "\n" );
  printf ( "           Minimum       Maximum\n" );
  printf ( "\n" );
  printf ( "  F:  %14.6g  %14.6g\n", r8vec_min ( n, f ), r8vec_max ( n, f ) );
  printf ( "  G:  %14.6g  %14.6g\n", r8vec_min ( n, g ), r8vec_max ( n, g ) );
  printf ( "  H:  %14.6g  %14.6g\n", r8vec_min ( n, h ), r8vec_max ( n, h ) );

  free ( f );
  free ( g );
  free ( h );
  free ( x );
  free ( y );

  return;
}
/******************************************************************************/

void resid_spiral_test ( )

/******************************************************************************/
/*
  Purpose:

    SPIRAL_TAYLOR_TEST samples the residual at the initial time.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 January 2015

  Author:

    John Burkardt
*/
{
  int n;
  double nu;
  double *pr;
  double rho;
  int seed;
  double t;
  double *ur;
  double *vr;
  double *x;
  double r8_hi;
  double r8_lo;
  double *y;

  nu = 1.0;
  rho = 1.0;

  printf ( "\n" );
  printf ( "RESID_SPIRAL_TEST\n" );
  printf ( "  Spiral Flow:\n" );
  printf ( "  Sample the Navier-Stokes residuals\n" );
  printf ( "  at the initial time T = 0, over the unit square.\n" );
  printf ( "  Kinematic viscosity NU = %g\n", nu );
  printf ( "  Fluid density RHO = %g\n", rho );

  n = 1000;

  pr = ( double * ) malloc ( n * sizeof ( double ) );
  ur = ( double * ) malloc ( n * sizeof ( double ) );
  vr = ( double * ) malloc ( n * sizeof ( double ) );

  r8_lo = 0.0;
  r8_hi = 1.0;
  seed = 123456789;

  x = r8vec_uniform_ab_new ( n, r8_lo, r8_hi, &seed );
  y = r8vec_uniform_ab_new ( n, r8_lo, r8_hi, &seed );
  t = 0.0;

  resid_spiral ( nu, rho, n, x, y, t, ur, vr, pr );

  printf ( "\n" );
  printf ( "           Minimum       Maximum\n" );
  printf ( "\n" );
  printf ( "  Ur:  %14.6g  %14.6g\n", r8vec_amin ( n, ur ), r8vec_amax ( n, ur ) );
  printf ( "  Vr:  %14.6g  %14.6g\n", r8vec_amin ( n, vr ), r8vec_amax ( n, vr ) );
  printf ( "  Pr:  %14.6g  %14.6g\n", r8vec_amin ( n, pr ), r8vec_amax ( n, pr ) );

  free ( pr );
  free ( ur );
  free ( vr );
  free ( x );
  free ( y );

  return;
}
/******************************************************************************/

void gnuplot_spiral_test ( )

/******************************************************************************/
/*
  Purpose:

    GNUPLOT_SPIRAL_TEST generates a field on a regular grid and plots it.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 January 2015

  Author:

    John Burkardt
*/
{
  char header[255];
  int n;
  double nu;
  double *p;
  double rho;
  double s;
  int seed;
  double t;
  double *u;
  double *v;
  double *x;
  double x_hi;
  double x_lo;
  int x_num = 21;
  double *y;
  double y_hi;
  double y_lo;
  int y_num = 21;

  printf ( "\n" );
  printf ( "GNUPLOT_SPIRAL_TEST:\n" );
  printf ( "  Spiral Flow:\n" );
  printf ( "  Generate a velocity field on a regular grid.\n" );
  printf ( "  Store in GNUPLOT data and command files.\n" );

  x_lo = 0.0;
  x_hi = 1.0;

  y_lo = 0.0;
  y_hi = 1.0;

  x = ( double * ) malloc ( x_num * y_num * sizeof ( double ) );
  y = ( double * ) malloc ( x_num * y_num * sizeof ( double ) );

  grid_2d ( x_num, x_lo, x_hi, y_num, y_lo, y_hi, x, y );

  nu = 1.0;
  rho = 1.0;
  n = x_num * y_num;
  t = 0.0;

  u = ( double * ) malloc ( x_num * y_num * sizeof ( double ) );
  v = ( double * ) malloc ( x_num * y_num * sizeof ( double ) );
  p = ( double * ) malloc ( x_num * y_num * sizeof ( double ) );

  uvp_spiral ( nu, rho, n, x, y, t, u, v, p );

  strcpy ( header, "spiral" );
  s = 5.00;
  ns2de_gnuplot ( header, n, x, y, u, v, s );

  free ( p );
  free ( u );
  free ( v );
  free ( x );
  free ( y );

  return;
}
/******************************************************************************/

void parameter_spiral_test ( )

/******************************************************************************/
/*
  Purpose:

    PARAMETER_SPIRAL_TEST monitors solution norms for various values of NU, RHO.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 January 2015

  Author:

    John Burkardt
*/
{
  int i;
  int j;
  int k;
  int n;
  double nu;
  double *p;
  double p_norm;
  double rho;
  int seed;
  double t;
  double *u;
  double u_norm;
  double *v;
  double v_norm;
  double *x;
  double r8_hi;
  double r8_lo;
  double *y;

  printf ( "\n" );
  printf ( "PARAMETER_SPIRAL_TEST\n" );
  printf ( "  Spiral Flow:\n" );
  printf ( "  Monitor solution norms over time for various\n" );
  printf ( "  values of NU, RHO.\n" );

  n = 1000;

  u = ( double * ) malloc ( n * sizeof ( double ) );
  v = ( double * ) malloc ( n * sizeof ( double ) );
  p = ( double * ) malloc ( n * sizeof ( double ) );

  r8_lo = 0.0;
  r8_hi = 1.0;
  seed = 123456789;

  x = r8vec_uniform_ab_new ( n, r8_lo, r8_hi, &seed );
  y = r8vec_uniform_ab_new ( n, r8_lo, r8_hi, &seed );
/*
  Vary RHO.
*/
  printf ( "\n" );
  printf ( "  RHO affects the pressure scaling.\n" );
  printf ( "\n" );
  printf ( "     RHO         NU           T     ||U||       ||V||       ||P||\n" );
  printf ( "\n" );

  nu = 1.0;
  rho = 1.0;

  for ( j = 1; j <= 3; j++ )
  {
    for ( k = 0; k <= 5; k++ )
    {
      t = ( double ) ( k ) / 5.0;

      uvp_spiral ( nu, rho, n, x, y, t, u, v, p );

      u_norm = r8vec_norm_l2 ( n, u ) / ( double ) ( n );
      v_norm = r8vec_norm_l2 ( n, v ) / ( double ) ( n );
      p_norm = r8vec_norm_l2 ( n, p ) / ( double ) ( n );

      printf ( "  %10.4g  %10.4g  %8.4f  %10.4g  %10.4g  %10.4g\n",
        rho, nu, t, u_norm, v_norm, p_norm );
    }
    printf ( "\n" );
    rho = rho / 100.0;
  }
/*
  Vary NU.
*/
  printf ( "\n" );
  printf ( "  NU affects the time scaling.\n" );
  printf ( "\n" );
  printf ( "     RHO         NU           T     ||U||       ||V||       ||P||\n" );
  printf ( "\n" );

  nu = 1.0;
  rho = 1.0;

  for ( i = 1; i <= 4; i++ )
  {
    for ( k = 0; k <= 5; k++ )
    {
      t = ( double ) ( k ) / 5.0;

      uvp_spiral ( nu, rho, n, x, y, t, u, v, p );

      u_norm = r8vec_norm_l2 ( n, u ) / ( double ) ( n );
      v_norm = r8vec_norm_l2 ( n, v ) / ( double ) ( n );
      p_norm = r8vec_norm_l2 ( n, p ) / ( double ) ( n );

      printf ( "  %10.4g  %10.4g  %8.4f  %10.4g  %10.4g  %10.4g\n",
        rho, nu, t, u_norm, v_norm, p_norm );
    }

    printf ( "\n" );

    nu = nu / 10.0;
  }

  free ( p );
  free ( u );
  free ( v );
  free ( x );
  free ( y );

  return;
}
/******************************************************************************/

void uvp_lucas_test ( )

/******************************************************************************/
/*
  Purpose:

    UVP_LUCAS_TEST samples the lucas flow solution at the initial time.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    06 March 2015

  Author:

    John Burkardt
*/
{
  int n;
  double nu;
  double *p;
  double rho;
  int seed;
  double t;
  double *u;
  double *v;
  double *x;
  double r8_hi;
  double r8_lo;
  double *y;

  nu = 1.0;
  rho = 1.0;

  printf ( "\n" );
  printf ( "UVP_LUCAS_TEST\n" );
  printf ( "  Lucas Bystricky Flow:\n" );
  printf ( "  Estimate the range of velocity and pressure\n" );
  printf ( "  at the initial time T = 0, using a region that is\n" );
  printf ( "  the unit square.\n" );
  printf ( "  Kinematic viscosity NU = %g\n", nu );
  printf ( "  Fluid density RHO = %g\n", rho );

  n = 1000;

  p = ( double * ) malloc ( n * sizeof ( double ) );
  u = ( double * ) malloc ( n * sizeof ( double ) );
  v = ( double * ) malloc ( n * sizeof ( double ) );

  r8_lo = 0.0;
  r8_hi = 1.0;
  seed = 123456789;

  x = r8vec_uniform_ab_new ( n, r8_lo, r8_hi, &seed );
  y = r8vec_uniform_ab_new ( n, r8_lo, r8_hi, &seed );
  t = 0.0;

  uvp_lucas ( nu, rho, n, x, y, t, u, v, p );

  printf ( "\n" );
  printf ( "           Minimum       Maximum\n" );
  printf ( "\n" );
  printf ( "  U:  %14.6g  %14.6g\n", r8vec_min ( n, u ), r8vec_max ( n, u ) );
  printf ( "  V:  %14.6g  %14.6g\n", r8vec_min ( n, v ), r8vec_max ( n, v ) );
  printf ( "  P:  %14.6g  %14.6g\n", r8vec_min ( n, p ), r8vec_max ( n, p ) );

  free ( p );
  free ( u );
  free ( v );
  free ( x );
  free ( y );

  return;
}
/******************************************************************************/

void uvp_lucas_test2 ( )

/******************************************************************************/
/*
  Purpose:

    UVP_LUCAS_TEST2 samples the solution on the boundary at the initial time.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    06 March 2015

  Author:

    John Burkardt
*/
{
  int i;
  int n;
  double nu;
  double *p;
  double rho;
  int seed;
  double t;
  double *u;
  double *v;
  double *x;
  double r8_hi;
  double r8_lo;
  double *y;

  nu = 1.0;
  rho = 1.0;
  t = 0.0;

  printf ( "\n" );
  printf ( "UVP_LUCAS_TEST2\n" );
  printf ( "  Lucas Bystricky Flow:\n" );
  printf ( "  Estimate the range of velocity and pressure\n" );
  printf ( "  on the boundary,\n" );
  printf ( "  at the initial time T = 0, using the unit square.\n" );
  printf ( "  Kinematic viscosity NU = %g\n", nu );
  printf ( "  Fluid density RHO = %g\n", rho );

  n = 400;

  p = ( double * ) malloc ( n * sizeof ( double ) );
  u = ( double * ) malloc ( n * sizeof ( double ) );
  v = ( double * ) malloc ( n * sizeof ( double ) );
  x = ( double * ) malloc ( n * sizeof ( double ) );
  y = ( double * ) malloc ( n * sizeof ( double ) );

  r8_lo = 0.0;
  r8_hi = 1.0;

  r8vec_linspace ( 100, r8_lo, r8_hi, x );
  for ( i = 0; i < 100; i++ )
  {
    y[i] = r8_lo;
  }

  for ( i = 100; i < 200; i++ )
  {
    x[i] = r8_hi;
  }
  r8vec_linspace ( 100, r8_lo, r8_hi, y+100 );

  r8vec_linspace ( 100, r8_hi, r8_lo, x+200 );
  for ( i = 200; i < 300; i++ )
  {
    y[i] = r8_hi;
  }

  for ( i = 300; i < 400; i++ )
  {
    x[i] = r8_lo;
  }
  r8vec_linspace ( 100, r8_lo, r8_hi, y+300 );

  uvp_lucas ( nu, rho, n, x, y, t, u, v, p );

  printf ( "\n" );
  printf ( "           Minimum       Maximum\n" );
  printf ( "\n" );
  printf ( "  U:  %14.6g  %14.6g\n", r8vec_min ( n, u ), r8vec_max ( n, u ) );
  printf ( "  V:  %14.6g  %14.6g\n", r8vec_min ( n, v ), r8vec_max ( n, v ) );
  printf ( "  P:  %14.6g  %14.6g\n", r8vec_min ( n, p ), r8vec_max ( n, p ) );

  free ( p );
  free ( u );
  free ( v );
  free ( x );
  free ( y );

  return;
}
/******************************************************************************/

void rhs_lucas_test ( )

/******************************************************************************/
/*
  Purpose:

    RHS_LUCAS_TEST samples the right hand sides at the initial time.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    06 March 2015

  Author:

    John Burkardt
*/
{
  double *f;
  double *g;
  double *h;
  int n;
  double nu;
  double rho;
  int seed;
  double t;
  double *x;
  double r8_hi;
  double r8_lo;
  double *y;

  nu = 1.0;
  rho = 1.0;

  printf ( "\n" );
  printf ( "RHS_LUCAS_TEST\n" );
  printf ( "  Lucas Bystricky Flow:\n" );
  printf ( "  Sample the Navier-Stokes right hand sides\n" );
  printf ( "  at the initial time T = 0, using the unit square.\n" );
  printf ( "  Kinematic viscosity NU = %g\n", nu );
  printf ( "  Fluid density RHO = %g\n", rho );

  n = 1000;

  f = ( double * ) malloc ( n * sizeof ( double ) );
  g = ( double * ) malloc ( n * sizeof ( double ) );
  h = ( double * ) malloc ( n * sizeof ( double ) );

  r8_lo = 0.0;
  r8_hi = 1.0;
  seed = 123456789;

  x = r8vec_uniform_ab_new ( n, r8_lo, r8_hi, &seed );
  y = r8vec_uniform_ab_new ( n, r8_lo, r8_hi, &seed );
  t = 0.0;

  rhs_lucas ( nu, rho, n, x, y, t, f, g, h );

  printf ( "\n" );
  printf ( "           Minimum       Maximum\n" );
  printf ( "\n" );
  printf ( "  F:  %14.6g  %14.6g\n", r8vec_min ( n, f ), r8vec_max ( n, f ) );
  printf ( "  G:  %14.6g  %14.6g\n", r8vec_min ( n, g ), r8vec_max ( n, g ) );
  printf ( "  H:  %14.6g  %14.6g\n", r8vec_min ( n, h ), r8vec_max ( n, h ) );

  free ( f );
  free ( g );
  free ( h );
  free ( x );
  free ( y );

  return;
}
/******************************************************************************/

void resid_lucas_test ( )

/******************************************************************************/
/*
  Purpose:

    LUCAS_TAYLOR_TEST samples the residual at the initial time.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    06 March 2015

  Author:

    John Burkardt
*/
{
  int n;
  double nu;
  double *pr;
  double rho;
  int seed;
  double t;
  double *ur;
  double *vr;
  double *x;
  double r8_hi;
  double r8_lo;
  double *y;

  nu = 1.0;
  rho = 1.0;

  printf ( "\n" );
  printf ( "RESID_LUCAS_TEST\n" );
  printf ( "  Lucas Bystricky Flow:\n" );
  printf ( "  Sample the Navier-Stokes residuals\n" );
  printf ( "  at the initial time T = 0, over the unit square.\n" );
  printf ( "  Kinematic viscosity NU = %g\n", nu );
  printf ( "  Fluid density RHO = %g\n", rho );

  n = 1000;

  pr = ( double * ) malloc ( n * sizeof ( double ) );
  ur = ( double * ) malloc ( n * sizeof ( double ) );
  vr = ( double * ) malloc ( n * sizeof ( double ) );

  r8_lo = 0.0;
  r8_hi = 1.0;
  seed = 123456789;

  x = r8vec_uniform_ab_new ( n, r8_lo, r8_hi, &seed );
  y = r8vec_uniform_ab_new ( n, r8_lo, r8_hi, &seed );
  t = 0.0;

  resid_lucas ( nu, rho, n, x, y, t, ur, vr, pr );

  printf ( "\n" );
  printf ( "           Minimum       Maximum\n" );
  printf ( "\n" );
  printf ( "  Ur:  %14.6g  %14.6g\n", r8vec_amin ( n, ur ), r8vec_amax ( n, ur ) );
  printf ( "  Vr:  %14.6g  %14.6g\n", r8vec_amin ( n, vr ), r8vec_amax ( n, vr ) );
  printf ( "  Pr:  %14.6g  %14.6g\n", r8vec_amin ( n, pr ), r8vec_amax ( n, pr ) );

  free ( pr );
  free ( ur );
  free ( vr );
  free ( x );
  free ( y );

  return;
}
/******************************************************************************/

void gnuplot_lucas_test ( )

/******************************************************************************/
/*
  Purpose:

    GNUPLOT_LUCAS_TEST generates a field on a regular grid and plots it.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    06 March 2015

  Author:

    John Burkardt
*/
{
  char header[255];
  int n;
  double nu;
  double *p;
  double rho;
  double s;
  int seed;
  double t;
  double *u;
  double *v;
  double *x;
  double x_hi;
  double x_lo;
  int x_num = 21;
  double *y;
  double y_hi;
  double y_lo;
  int y_num = 21;

  printf ( "\n" );
  printf ( "GNUPLOT_LUCAS_TEST:\n" );
  printf ( "  Lucas Bystricky Flow:\n" );
  printf ( "  Generate a velocity field on a regular grid.\n" );
  printf ( "  Store in GNUPLOT data and command files.\n" );

  x_lo = 0.0;
  x_hi = 1.0;

  y_lo = 0.0;
  y_hi = 1.0;

  x = ( double * ) malloc ( x_num * y_num * sizeof ( double ) );
  y = ( double * ) malloc ( x_num * y_num * sizeof ( double ) );

  grid_2d ( x_num, x_lo, x_hi, y_num, y_lo, y_hi, x, y );

  nu = 1.0;
  rho = 1.0;
  n = x_num * y_num;
  t = 0.0;

  u = ( double * ) malloc ( x_num * y_num * sizeof ( double ) );
  v = ( double * ) malloc ( x_num * y_num * sizeof ( double ) );
  p = ( double * ) malloc ( x_num * y_num * sizeof ( double ) );

  uvp_lucas ( nu, rho, n, x, y, t, u, v, p );

  strcpy ( header, "lucas" );
  s = 0.25;
  ns2de_gnuplot ( header, n, x, y, u, v, s );

  free ( p );
  free ( u );
  free ( v );
  free ( x );
  free ( y );

  return;
}
