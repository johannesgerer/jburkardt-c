# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>

# include "ns3de.h"

int main ( );
void test01 ( );
void test02 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    NS3DE_PRB tests the NS3DE library.

  Location:

    http://people.sc.fsu.edu/~jburkardt/c_src/navier_stokes_3d_exact/ns3de_prb.c

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
  printf ( "NS3DE_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the NS3DE library.\n" );

  test01 ( );
  test02 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "NS3DE_PRB\n" );
  printf ( "  Normal end of execution.\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void test01 ( )

/******************************************************************************/
/*
  Purpose:

    TEST01 samples the solution at the initial time.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 January 2015

  Author:

    John Burkardt
*/
{
  double a;
  double d;
  int n;
  double *p;
  const double r8_pi = 3.141592653589793;
  int seed;
  double t;
  double *u;
  double *v;
  double *w;
  double *x;
  double xyz_hi;
  double xyz_lo;
  double *y;
  double *z;

  a = r8_pi / 4.0;
  d = r8_pi / 2.0;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  Estimate the range of velocity and pressure\n" );
  printf ( "  at the initial time T = 0, in a region that is the\n" );
  printf ( "  cube centered at (0,0,0) with 'radius' 1.0.\n" );
  printf ( "  Parameter A = %g\n", a );
  printf ( "  Parameter D = %g\n", d );

  n = 1000;

  p = ( double * ) malloc ( n * sizeof ( double ) );
  u = ( double * ) malloc ( n * sizeof ( double ) );
  v = ( double * ) malloc ( n * sizeof ( double ) );
  w = ( double * ) malloc ( n * sizeof ( double ) );

  xyz_lo = -1.0;
  xyz_hi = +1.0;
  seed = 123456789;

  x = r8vec_uniform_ab_new ( n, xyz_lo, xyz_hi, &seed );
  y = r8vec_uniform_ab_new ( n, xyz_lo, xyz_hi, &seed );
  z = r8vec_uniform_ab_new ( n, xyz_lo, xyz_hi, &seed );
  t = 0.0;

  uvwp_ethier ( a, d, n, x, y, z, t, u, v, w, p );

  printf ( "\n" );
  printf ( "           Minimum       Maximum\n" );
  printf ( "\n" );
  printf ( "  U:    %14.6g  %14.6g\n", r8vec_min ( n, u ), r8vec_max ( n, u ) );
  printf ( "  V:    %14.6g  %14.6g\n", r8vec_min ( n, v ), r8vec_max ( n, v ) );
  printf ( "  W:    %14.6g  %14.6g\n", r8vec_min ( n, w ), r8vec_max ( n, w ) );
  printf ( "  P:    %14.6g  %14.6g\n", r8vec_min ( n, p ), r8vec_max ( n, p ) );

  free ( p );
  free ( u );
  free ( v );
  free ( w );
  free ( x );
  free ( y );
  free ( z );

  return;
}
/******************************************************************************/

void test02 ( )

/******************************************************************************/
/*
  Purpose:

    TEST02 samples the residual at the initial time.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 January 2015

  Author:

    John Burkardt
*/
{
  double a;
  double d;
  int n;
  double *pr;
  const double r8_pi = 3.141592653589793;
  int seed;
  double t;
  double *ur;
  double *vr;
  double *wr;
  double *x;
  double xyz_hi;
  double xyz_lo;
  double *y;
  double *z;

  a = r8_pi / 4.0;
  d = r8_pi / 2.0;

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  Sample the Navier-Stokes residuals\n" );
  printf ( "  at the initial time T = 0, using a region that is\n" );
  printf ( "  the cube centered at (0,0,0) with 'radius' 1.0,\n" );
  printf ( "  Parameter A = %g\n", a );
  printf ( "  Parameter D = %g\n", d );

  n = 1000;

  pr = ( double * ) malloc ( n * sizeof ( double ) );
  ur = ( double * ) malloc ( n * sizeof ( double ) );
  vr = ( double * ) malloc ( n * sizeof ( double ) );
  wr = ( double * ) malloc ( n * sizeof ( double ) );

  xyz_lo = -1.0;
  xyz_hi = +1.0;
  seed = 123456789;

  x = r8vec_uniform_ab_new ( n, xyz_lo, xyz_hi, &seed );
  y = r8vec_uniform_ab_new ( n, xyz_lo, xyz_hi, &seed );
  z = r8vec_uniform_ab_new ( n, xyz_lo, xyz_hi, &seed );
  t = 0.0;

  resid_ethier ( a, d, n, x, y, z, t, ur, vr, wr, pr );

  printf ( "\n" );
  printf ( "           Minimum       Maximum\n" );
  printf ( "\n" );
  printf ( "  Ur:    %14.6g  %14.6g\n", r8vec_amin ( n, ur ), r8vec_amax ( n, ur ) );
  printf ( "  Vr:    %14.6g  %14.6g\n", r8vec_amin ( n, vr ), r8vec_amax ( n, vr ) );
  printf ( "  Wr:    %14.6g  %14.6g\n", r8vec_amin ( n, wr ), r8vec_amax ( n, wr ) );
  printf ( "  Pr:    %14.6g  %14.6g\n", r8vec_amin ( n, pr ), r8vec_amax ( n, pr ) );

  free ( pr );
  free ( ur );
  free ( vr );
  free ( wr );
  free ( x );
  free ( y );
  free ( z );

  return;
}

