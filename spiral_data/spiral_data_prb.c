# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <string.h>
# include <time.h>

# include "spiral_data.h"

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for SPIRAL_DATA_PRB.

  Discussion:

    SPIRAL_DATA_PRB tests the SPIRAL_DATA library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 January 2015

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "SPIRAL_DATA_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the SPIRAL_DATA library.\n" );

  test01 ( );
  test02 ( );
  test03 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "SPIRAL_DATA_PRB\n" );
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

    TEST01 generates a field and estimates its range.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 January 2015

  Author:

    John Burkardt
*/
{
  int n;
  double c;
  int seed;
  double *u;
  double *v;
  double *x;
  double xy_hi;
  double xy_lo;
  double *y;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  Sample a spiral velocity field and estimate\n" );
  printf ( "  the range of the solution values.\n" );

  n = 1000;

  xy_lo = +0.0;
  xy_hi = +1.0;
  seed = 123456789;

  x = r8vec_uniform_ab_new ( n, xy_lo, xy_hi, &seed );
  y = r8vec_uniform_ab_new ( n, xy_lo, xy_hi, &seed );
  c = 1.0;

  u = ( double * ) malloc ( n * sizeof ( double ) );
  v = ( double * ) malloc ( n * sizeof ( double ) );
  uv_spiral ( n, x, y, c, u, v );

  printf ( "\n" );
  printf ( "           Minimum       Maximum\n" );
  printf ( "\n" );
  printf ( "  U:    %14.6g  %14.6g\n", r8vec_min ( n, u ), r8vec_max ( n, u ) );
  printf ( "  V:    %14.6g  %14.6g\n", r8vec_min ( n, v ), r8vec_max ( n, v ) );

  free ( u );
  free ( v );
  free ( x );
  free ( y );

  return;
}
/******************************************************************************/

void test02 ( )

/******************************************************************************/
/*
  Purpose:

    TEST02 generates a field and samples its residuals.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 January 2015

  Author:

    John Burkardt
*/
{
  int n;
  double c;
  double *pr;
  int seed;
  double *x;
  double xy_hi;
  double xy_lo;
  double *y;

  printf ( "\n" );
  printf ( "TEST02:\n" );
  printf ( "  Sample a spiral velocity field and estimate the\n" );
  printf ( "  range of residuals in the continuity equation.\n" );

  n = 1000;

  xy_lo = +0.0;
  xy_hi = +1.0;
  seed = 123456789;

  x = r8vec_uniform_ab_new ( n, xy_lo, xy_hi, &seed );
  y = r8vec_uniform_ab_new ( n, xy_lo, xy_hi, &seed );
  c = 1.0;

  pr = ( double * ) malloc ( n * sizeof ( double ) );
  resid_spiral ( n, x, y, c, pr );

  printf ( "\n" );
  printf ( "           Minimum       Maximum\n" );
  printf ( "\n" );
  printf ( "  Pr:    %14.6g  %14.6g\n", r8vec_amin ( n, pr ), r8vec_amax ( n, pr ) );

  free ( pr );
  free ( x );
  free ( y );

  return;
}
/******************************************************************************/

void test03 ( )

/******************************************************************************/
/*
  Purpose:

    TEST03 generates a field on a regular grid and plots it.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 January 2015

  Author:

    John Burkardt
*/
{
  double c;
  char header[255];
  int n;
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
  printf ( "TEST03:\n" );
  printf ( "  Generate a spiral velocity field on a regular grid.\n" );
  printf ( "  Store in GNUPLOT data and command files.\n" );

  x_lo = 0.0;
  x_hi = 1.0;

  y_lo = 0.0;
  y_hi = 1.0;

  x = ( double * ) malloc ( x_num * y_num * sizeof ( double ) );
  y = ( double * ) malloc ( x_num * y_num * sizeof ( double ) );

  grid_2d ( x_num, x_lo, x_hi, y_num, y_lo, y_hi, x, y );

  n = x_num * y_num;
  c = 1.0;

  u = ( double * ) malloc ( x_num * y_num * sizeof ( double ) );
  v = ( double * ) malloc ( x_num * y_num * sizeof ( double ) );

  uv_spiral ( n, x, y, c, u, v );

  strcpy ( header, "spiral" );
  s = 0.05;
  spiral_gnuplot ( header, n, x, y, u, v, s );

  free ( u );
  free ( v );
  free ( x );
  free ( y );

  return;
}