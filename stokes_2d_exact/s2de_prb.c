# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <string.h>
# include <time.h>

# include "s2de.h"

int main ( );
void uvp_stokes1_test ( );
void resid_stokes1_test ( );
void gnuplot_stokes1_test ( );
void uvp_stokes2_test ( );
void resid_stokes2_test ( );
void gnuplot_stokes2_test ( );
void uvp_stokes3_test ( );
void resid_stokes3_test ( );
void gnuplot_stokes3_test ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    S2DE_PRB tests the S2DE library.

  Location:

    http://people.sc.fsu.edu/~jburkardt/c_src/stokes_2d_exact/s2de_prb.c

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    24 January 2015

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "S2DE_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the S2DE library.\n" );

  uvp_stokes1_test ( );
  resid_stokes1_test ( );
  gnuplot_stokes1_test ( );

  uvp_stokes2_test ( );
  resid_stokes2_test ( );
  gnuplot_stokes2_test ( );

  uvp_stokes3_test ( );
  resid_stokes3_test ( );
  gnuplot_stokes3_test ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "S2DE_PRB\n" );
  printf ( "  Normal end of execution.\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void uvp_stokes1_test ( )

/******************************************************************************/
/*
  Purpose:

    UVP_STOKES1_TEST samples the solution #1.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    24 January 2015

  Author:

    John Burkardt
*/
{
  int n = 1000;
  double *p;
  int seed;
  double *u;
  double *v;
  double *x;
  double xy_hi;
  double xy_lo;
  double *y;

  printf ( "\n" );
  printf ( "UVP_STOKES1_TEST\n" );
  printf ( "  Exact Stokes solution #1:\n" );
  printf ( "  Estimate the range of velocity and pressure\n" );
  printf ( "  using a region that is the unit square.\n" );

  xy_lo = 0.0;
  xy_hi = 1.0;
  seed = 123456789;

  x = r8vec_uniform_ab_new ( n, xy_lo, xy_hi, &seed );
  y = r8vec_uniform_ab_new ( n, xy_lo, xy_hi, &seed );

  u = ( double * ) malloc ( n * sizeof ( double ) );
  v = ( double * ) malloc ( n * sizeof ( double ) );
  p = ( double * ) malloc ( n * sizeof ( double ) );

  uvp_stokes1 ( n, x, y, u, v, p );

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

void resid_stokes1_test ( )

/******************************************************************************/
/*
  Purpose:

    RESID_STOKES1_TEST samples the residual for solution #1.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    24 January 2015

  Author:

    John Burkardt
*/
{
  int n = 1000;
  double *pr;
  int seed;
  double *ur;
  double *vr;
  double *x;
  double xy_hi;
  double xy_lo;
  double *y;

  printf ( "\n" );
  printf ( "RESID_STOKES1_TEST\n" );
  printf ( "  Exact Stokes solution #1:\n" );
  printf ( "  Sample the Stokes residuals\n" );
  printf ( "  using a region that is the unit square.\n" );

  xy_lo = 0.0;
  xy_hi = 1.0;
  seed = 123456789;

  x = r8vec_uniform_ab_new ( n, xy_lo, xy_hi, &seed );
  y = r8vec_uniform_ab_new ( n, xy_lo, xy_hi, &seed );

  ur = ( double * ) malloc ( n * sizeof ( double ) );
  vr = ( double * ) malloc ( n * sizeof ( double ) );
  pr = ( double * ) malloc ( n * sizeof ( double ) );

  resid_stokes1 ( n, x, y, ur, vr, pr );

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

void gnuplot_stokes1_test ( )

/******************************************************************************/
/*
  Purpose:

    GNUPLOT_STOKES1_TEST plots solution #1 on a regular grid.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    24 January 2015

  Author:

    John Burkardt
*/
{
  char header[255];
  int n;
  double *p;
  double s;
  int seed;
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
  printf ( "GNUPLOT_STOKES1_TEST:\n" );
  printf ( "  Exact Stokes solution #1:\n" );
  printf ( "  Generate a Stokes velocity field on a regular grid.\n" );
  printf ( "  Store in GNUPLOT data and command files.\n" );
   
  x_lo = 0.0;
  x_hi = 1.0;

  y_lo = 0.0;
  y_hi = 1.0;

  x = ( double * ) malloc ( x_num * y_num * sizeof ( double ) );
  y = ( double * ) malloc ( x_num * y_num * sizeof ( double ) );

  grid_2d ( x_num, x_lo, x_hi, y_num, y_lo, y_hi, x, y );

  n = x_num * y_num;

  u = ( double * ) malloc ( n * sizeof ( double ) );
  v = ( double * ) malloc ( n * sizeof ( double ) );
  p = ( double * ) malloc ( n * sizeof ( double ) );

  uvp_stokes1 ( n, x, y, u, v, p );

  strcpy ( header, "stokes1" );
  s = 4.0;
  stokes_gnuplot ( header, n, x, y, u, v, s );

  free ( p );
  free ( u );
  free ( v );
  free ( x );
  free ( y );

  return;
}
/******************************************************************************/

void uvp_stokes2_test ( )

/******************************************************************************/
/*
  Purpose:

    UVP_STOKES2_TEST samples the solution #2.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    24 January 2015

  Author:

    John Burkardt
*/
{
  int n = 1000;
  double *p;
  int seed;
  double *u;
  double *v;
  double *x;
  double xy_hi;
  double xy_lo;
  double *y;

  printf ( "\n" );
  printf ( "UVP_STOKES2_TEST\n" );
  printf ( "  Exact Stokes solution #2:\n" );
  printf ( "  Estimate the range of velocity and pressure\n" );
  printf ( "  using a region that is the unit square.\n" );

  xy_lo = 0.0;
  xy_hi = 1.0;
  seed = 123456789;

  x = r8vec_uniform_ab_new ( n, xy_lo, xy_hi, &seed );
  y = r8vec_uniform_ab_new ( n, xy_lo, xy_hi, &seed );

  u = ( double * ) malloc ( n * sizeof ( double ) );
  v = ( double * ) malloc ( n * sizeof ( double ) );
  p = ( double * ) malloc ( n * sizeof ( double ) );

  uvp_stokes2 ( n, x, y, u, v, p );

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

void resid_stokes2_test ( )

/******************************************************************************/
/*
  Purpose:

    RESID_STOKES2_TEST samples the residual for solution #1.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    24 January 2015

  Author:

    John Burkardt
*/
{
  int n = 1000;
  double *pr;
  int seed;
  double *ur;
  double *vr;
  double *x;
  double xy_hi;
  double xy_lo;
  double *y;

  printf ( "\n" );
  printf ( "RESID_STOKES2_TEST\n" );
  printf ( "  Exact Stokes solution #2:\n" );
  printf ( "  Sample the Stokes residuals\n" );
  printf ( "  using a region that is the unit square.\n" );

  xy_lo = 0.0;
  xy_hi = 1.0;
  seed = 123456789;

  x = r8vec_uniform_ab_new ( n, xy_lo, xy_hi, &seed );
  y = r8vec_uniform_ab_new ( n, xy_lo, xy_hi, &seed );

  ur = ( double * ) malloc ( n * sizeof ( double ) );
  vr = ( double * ) malloc ( n * sizeof ( double ) );
  pr = ( double * ) malloc ( n * sizeof ( double ) );

  resid_stokes2 ( n, x, y, ur, vr, pr );

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

void gnuplot_stokes2_test ( )

/******************************************************************************/
/*
  Purpose:

    GNUPLOT_STOKES2_TEST plots solution #2 on a regular grid.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    24 January 2015

  Author:

    John Burkardt
*/
{

  char header[255];
  int n;
  double *p;
  double s;
  int seed;
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
  printf ( "GNUPLOT_STOKES2_TEST:\n" );
  printf ( "  Exact Stokes solution #2:\n" );
  printf ( "  Generate a Stokes velocity field on a regular grid.\n" );
  printf ( "  Store in GNUPLOT data and command files.\n" );
   
  x_lo = 0.0;
  x_hi = 1.0;

  y_lo = 0.0;
  y_hi = 1.0;

  x = ( double * ) malloc ( x_num * y_num * sizeof ( double ) );
  y = ( double * ) malloc ( x_num * y_num * sizeof ( double ) );

  grid_2d ( x_num, x_lo, x_hi, y_num, y_lo, y_hi, x, y );

  n = x_num * y_num;

  u = ( double * ) malloc ( n * sizeof ( double ) );
  v = ( double * ) malloc ( n * sizeof ( double ) );
  p = ( double * ) malloc ( n * sizeof ( double ) );

  uvp_stokes2 ( n, x, y, u, v, p );

  strcpy ( header, "stokes2" );
  s = 0.05;
  stokes_gnuplot ( header, n, x, y, u, v, s );

  free ( p );
  free ( u );
  free ( v );
  free ( x );
  free ( y );

  return;
}
/******************************************************************************/

void uvp_stokes3_test ( )

/******************************************************************************/
/*
  Purpose:

    UVP_STOKES3_TEST samples the solution #3.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 February 2015

  Author:

    John Burkardt
*/
{
  int n = 1000;
  double *p;
  int seed;
  double *u;
  double *v;
  double *x;
  double xy_hi;
  double xy_lo;
  double *y;

  printf ( "\n" );
  printf ( "UVP_STOKES3_TEST\n" );
  printf ( "  Exact Stokes solution #3:\n" );
  printf ( "  Estimate the range of velocity and pressure\n" );
  printf ( "  using a region that is [-1,+1]x[-1,+1].\n" );

  xy_lo = -1.0;
  xy_hi = +1.0;
  seed = 123456789;

  x = r8vec_uniform_ab_new ( n, xy_lo, xy_hi, &seed );
  y = r8vec_uniform_ab_new ( n, xy_lo, xy_hi, &seed );

  u = ( double * ) malloc ( n * sizeof ( double ) );
  v = ( double * ) malloc ( n * sizeof ( double ) );
  p = ( double * ) malloc ( n * sizeof ( double ) );

  uvp_stokes3 ( n, x, y, u, v, p );

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

void resid_stokes3_test ( )

/******************************************************************************/
/*
  Purpose:

    RESID_STOKES3_TEST samples the residual for solution #3.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 February 2015

  Author:

    John Burkardt
*/
{
  int n = 1000;
  double *pr;
  int seed;
  double *ur;
  double *vr;
  double *x;
  double xy_hi;
  double xy_lo;
  double *y;

  printf ( "\n" );
  printf ( "RESID_STOKES3_TEST\n" );
  printf ( "  Exact Stokes solution #3:\n" );
  printf ( "  Sample the Stokes residuals\n" );
  printf ( "  using a region that is [-1,+1]x[-1,+1].\n" );

  xy_lo = -1.0;
  xy_hi = +1.0;
  seed = 123456789;

  x = r8vec_uniform_ab_new ( n, xy_lo, xy_hi, &seed );
  y = r8vec_uniform_ab_new ( n, xy_lo, xy_hi, &seed );

  ur = ( double * ) malloc ( n * sizeof ( double ) );
  vr = ( double * ) malloc ( n * sizeof ( double ) );
  pr = ( double * ) malloc ( n * sizeof ( double ) );

  resid_stokes3 ( n, x, y, ur, vr, pr );

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

void gnuplot_stokes3_test ( )

/******************************************************************************/
/*
  Purpose:

    GNUPLOT_STOKES3_TEST plots solution #3 on a regular grid.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 February 2015

  Author:

    John Burkardt
*/
{
  char header[255];
  int n;
  double *p;
  double s;
  int seed;
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
  printf ( "GNUPLOT_STOKES3_TEST:\n" );
  printf ( "  Exact Stokes solution #3:\n" );
  printf ( "  Generate a Stokes velocity field on [-1,+1]x[-1,+1].\n" );
  printf ( "  Store in GNUPLOT data and command files.\n" );
   
  x_lo = -1.0;
  x_hi = +1.0;

  y_lo = -1.0;
  y_hi = +1.0;

  x = ( double * ) malloc ( x_num * y_num * sizeof ( double ) );
  y = ( double * ) malloc ( x_num * y_num * sizeof ( double ) );

  grid_2d ( x_num, x_lo, x_hi, y_num, y_lo, y_hi, x, y );

  n = x_num * y_num;

  u = ( double * ) malloc ( n * sizeof ( double ) );
  v = ( double * ) malloc ( n * sizeof ( double ) );
  p = ( double * ) malloc ( n * sizeof ( double ) );

  uvp_stokes3 ( n, x, y, u, v, p );

  strcpy ( header, "stokes3" );
  s = 0.05;
  stokes_gnuplot ( header, n, x, y, u, v, s );

  free ( p );
  free ( u );
  free ( v );
  free ( x );
  free ( y );

  return;
}