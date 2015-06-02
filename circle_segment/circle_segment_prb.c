# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>

# include "circle_segment.h"

int main ( );
void test01 ( );
void test05 ( );
void test06 ( );
void test07 ( );
void test08 ( );
void test09 ( );
void test11 ( );
void test13 ( );
void test14 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for CIRCLE_SEGMENT_PRB.

  Discussion:

    CIRCLE_SEGMENT_PRB tests the CIRCLE_SEGMENT library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 July 2013

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "CIRCLE_SEGMENT_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the CIRCLE_SEGMENT library.\n" );

  test01 ( );
  test05 ( );
  test06 ( );
  test07 ( );
  test08 ( );
  test09 ( );
  test11 ( );
  test13 ( );
  test14 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "CIRCLE_SEGMENT_PRB\n" );
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

    TEST01 tests CIRCLE_SEGMENT_AREA_FROM_HEIGHT.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    06 July 2013

  Author:

    John Burkardt
*/
{
  double area;
  double h;
  int i;
  double r;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  CIRCLE_SEGMENT_AREA_FROM_HEIGHT computes the area of a circle segment.\n" );

  printf ( "\n" );
  printf ( "          R               H               Area\n" );
  printf ( "\n" );
  r = 1.0;
  h = 1.0;
  for ( i = 0; i <= 10; i++ )
  {
    area = circle_segment_area_from_height ( r, h );
    printf ( "  %14.6f  %14.6f  %14.6f\n", r, h, area );
    h = h / 2.0;
  }

  return;
}
/******************************************************************************/

void test05 ( )

/******************************************************************************/
/*
  Purpose:

    TEST05 tests the AREA and HEIGHT functions.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    06 July 2013

  Author:

    John Burkardt
*/
{
  double a;
  double a2;
  double h;
  double h2;
  const double pi = 3.141592653589793;
  double r;
  int seed;
  int test;

  printf ( "\n" );
  printf ( "CIRCLE_SEGMENT_TEST05\n" );
  printf ( "  For circle segment with a given radius R,\n" );
  printf ( "  CIRCLE_SEGMENT_AREA_FROM_HEIGHT computes the area A, given the height.\n" );
  printf ( "  CIRCLE_SEGMENT_HEIGHT_FROM_AREA computes height H, given the area.\n" );
  printf ( "  Check that these functions are inverses of each other\n" );
  printf ( "  using random values of R, A, and H.\n" );

  printf ( "\n" );
  printf ( "        R             H      =>     A    =>       H2\n" );
  printf ( "\n" );

  seed = 123456789;

  for ( test = 1; test <= 5; test++ )
  {
    r = 5.0 * r8_uniform_01 ( &seed );
    h = 2.0 * r * r8_uniform_01 ( &seed );
    a = circle_segment_area_from_height ( r, h );
    h2 = circle_segment_height_from_area ( r, a );
    printf ( "  %12.6f  %12.6f  %12.6f  %12.6f\n", r, h, a, h2 );
  }

  printf ( "\n" );
  printf ( "        R             A      =>     H    =>       A2\n" );
  printf ( "\n" );

  for ( test = 1; test <= 5; test++ )
  {
    r = 5.0 * r8_uniform_01 ( &seed );
    a = pi * r * r * r8_uniform_01 ( &seed );
    h = circle_segment_height_from_area ( r, a );
    a2 = circle_segment_area_from_height ( r, h );
    printf ( "  %12.6f  %12.6f  %12.6f  %12.6f\n", r, a, h, a2 );
  }

  return;
}
/******************************************************************************/

void test06 ( )

/******************************************************************************/
/*
  Purpose:

    TEST06 samples using CIRCLE_SEGMENT_SAMPLE_FROM_HEIGHT.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    06 July 2013

  Author:

    John Burkardt
*/
{
  double *an;
  int an_num = 51;
  char boundary_filename[] = "sample00_boundary.txt";
  FILE *boundary_unit;
  double *boundary_x;
  double *boundary_y;
  char command_filename[] = "sample00_commands.txt";
  FILE *command_unit;
  char data_filename[] = "sample00_data.txt";
  FILE *data_unit;
  int data_num = 100;
  double *data_x;
  double *data_y;
  char graphics_filename[] = "sample00.png";
  double h;
  int i;
  const double pi = 3.141592653589793;
  double r;
  int seed;
  int test;
  double theta;
  double thetah;

  seed = 123456789;

  printf ( "\n" );
  printf ( "CIRCLE_SEGMENT_TEST06\n" );
  printf ( "  CIRCLE_SEGMENT_SAMPLE_FROM_HEIGHT samples a circle segment.\n" );
  printf ( "\n" );
  printf ( "  Plot %d points from several segments.\n", data_num );
  printf ( "\n" );

  r = 1.0;
  theta = pi;

  for ( test = 1; test <= 4; test++ )
  {
    h = circle_segment_height_from_angle ( r, theta );

    thetah = theta / 2.0;
/*
  Create boundary.
*/
    an = r8vec_linspace_new ( an_num, -thetah, +thetah );
    for ( i = 0; i < an_num; i++ )
    {
      an[i] = an[i] + 0.5 * pi;
    }

    boundary_x = ( double * ) malloc ( ( an_num + 1 ) * sizeof ( double ) );
    boundary_y = ( double * ) malloc ( ( an_num + 1 ) * sizeof ( double ) );

    for ( i = 0; i < an_num; i++ )
    {
      boundary_x[i] = r * cos ( an[i] );
      boundary_y[i] = r * sin ( an[i] );
    }
    boundary_x[an_num] = boundary_x[0];
    boundary_y[an_num] = boundary_y[0];

    filename_inc ( boundary_filename );
    boundary_unit = fopen ( boundary_filename, "wt" );
    for ( i = 0; i <= an_num; i++ )
    {
      fprintf ( boundary_unit, "  %14.6g  %14.6g\n", boundary_x[i], boundary_y[i] );
    }
    fclose ( boundary_unit );
    printf ( "\n" );
    printf ( "  Created boundary file \"%s\".\n", boundary_filename );
/*
  Create data.
*/
    data_x = ( double * ) malloc ( ( data_num + 1 ) * sizeof ( double ) );
    data_y = ( double * ) malloc ( ( data_num + 1 ) * sizeof ( double ) );

    circle_segment_sample_from_height ( r, h, data_num, &seed, data_x, data_y );

    filename_inc ( data_filename );
    data_unit = fopen ( data_filename, "wt" );
    for ( i = 0; i < data_num; i++ )
    {
      fprintf ( data_unit, "  %14.6g  %14.6g\n", data_x[i], data_y[i] );
    }
    fclose ( data_unit );
    printf ( "\n" );
    printf ( "  Created data file \"%s\".\n", data_filename );
/*
  Create commands.
*/
    filename_inc ( command_filename );
    command_unit = fopen ( command_filename, "wt" );
    fprintf ( command_unit, "# %s\n", command_filename );
    fprintf ( command_unit, "#\n" );
    fprintf ( command_unit, "# Usage:\n" );
    fprintf ( command_unit, "#  gnuplot < %s\n", command_filename );
    fprintf ( command_unit, "#\n" );
    fprintf ( command_unit, "set term png\n" );
    filename_inc ( graphics_filename );
    fprintf ( command_unit, "set output '%s'\n", graphics_filename );
    fprintf ( command_unit, "set xlabel '<--- X --->'\n" );
    fprintf ( command_unit, "set ylabel '<--- Y --->'\n" );
    fprintf ( command_unit, "set title 'Circle Segment Sample'\n" );
    fprintf ( command_unit, "set grid\n" );
    fprintf ( command_unit, "set key off\n" );
    fprintf ( command_unit, "set size ratio -1\n" );
    fprintf ( command_unit, "set style data lines\n" );
    fprintf ( command_unit, "plot '%s' using 1:2 with points lt 3 pt 3,\\\n", data_filename );
    fprintf ( command_unit, "    '%s' using 1:2 lw 3 linecolor rgb 'black'\n", boundary_filename );
    fprintf ( command_unit, "quit\n" );
    fclose ( command_unit );

    printf ( "  Created command file \"%s\".\n", command_filename );

    theta = theta / 2.0;

    free ( an );
    free ( boundary_x );
    free ( boundary_y );
    free ( data_x );
    free ( data_y );
  }
 
  return;
}
/******************************************************************************/

void test07 ( )

/******************************************************************************/
/*
  Purpose:

    TEST07 tests the ANGLE and HEIGHT functions.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    06 July 2013

  Author:

    John Burkardt
*/
{
  double h;
  double h2;
  const double pi = 3.141592653589793;
  double r;
  int seed;
  double t;
  double t2;
  int test;

  printf ( "\n" );
  printf ( "TEST07\n" );
  printf ( "  For circle segment with a given radius R,\n" );
  printf ( "  CIRCLE_SEGMENT_ANGLE_FROM_HEIGHT computes the angle THETA, given the height.\n" );
  printf ( "  CIRCLE_SEGMENT_HEIGHT_FROM_ANGLE computes height H, given the angle.\n" );
  printf ( "  Check that these functions are inverses of each other\n" );
  printf ( "  using random values of R, T, and H.\n" );
  printf ( "\n" );
  printf ( "        R             H      =>     T    =>       H2\n" );
  printf ( "\n" );

  seed = 123456789;

  for ( test = 1; test <= 5; test++ )
  {
    r = 5.0 * r8_uniform_01 ( &seed );
    h = 2.0 * r * r8_uniform_01 ( &seed );
    t = circle_segment_angle_from_height ( r, h );
    h2 = circle_segment_height_from_angle ( r, t );
    printf ( "  %12.6f  %12.6f  %12.6f  %12.6f\n", r, h, t, h2 );
  }

  printf ( "\n" );
  printf ( "        R             T      =>     H    =>       T2\n" );
  printf ( "\n" );
  for ( test = 1; test <= 5; test++ )
  {
    r = 5.0 * r8_uniform_01 ( &seed );
    t = 2.0 * pi * r8_uniform_01 ( &seed );
    h = circle_segment_height_from_angle ( r, t );
    t2 = circle_segment_angle_from_height ( r, h );
    printf ( "  %12.6f  %12.6f  %12.6f  %12.6f\n", r, t, h, t2 );
  }

  return;
}
/******************************************************************************/

void test08 ( )

/******************************************************************************/
/*
  Purpose:

    TEST08 tests CIRCLE_SEGMENT_CONTAINS_POINT.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    06 July 2013

  Author:

    John Burkardt
*/
{
  double area;
  double area_est;
  double c[2];
  int i;
  int *inout;
  int j;
  int n = 1000;
  double omega1;
  double omega2;
  const double pi = 3.141592653589793;
  double r;
  int seed;
  int test;
  double theta;
  double *xy;

  printf ( "\n" );
  printf ( "TEST08\n" );
  printf ( "  CIRCLE_SEGMENT_CONTAINS_POINT reports whether\n" );
  printf ( "  a circle segment contains a point.\n" );
  printf ( "\n" );
  printf ( "  Pick a circle segment at random.\n" );
  printf ( "  Compute %d sample points in the surrounding box.\n", n );
  printf ( "  Compare the area of the segment to the percentage of points\n" );
  printf ( "  contained in the circle segment.\n" );
  printf ( "\n" );
  printf ( "       N       Omega1          Omega2           Area         Estimate\n" );
  printf ( "\n" );

  r = 1.0;
  c[0] = 0.0;
  c[1] = 0.0;
  seed = 123456789;
  inout = ( int * ) malloc ( n * sizeof ( int ) );

  for ( test = 1; test <= 5; test++ )
  {
    omega1 = 2.0 * pi * r8_uniform_01 ( &seed );
    omega2 = 2.0 * pi * r8_uniform_01 ( &seed );
  
    if ( omega2 < omega1 )
    {
      omega2 = omega2 + 2.0 * pi;
    }

    xy = r8mat_uniform_01_new ( 2, n, &seed );
    for ( j = 0; j < n; j++ )
    {
      for ( i = 0; i < 2; i++ )
      {
        xy[i+j*2] = 2.0 * xy[i+j*2] - 1.0;
      }
    }

    for ( j = 0; j < n; j++ )
    {
      inout[j] = circle_segment_contains_point ( r, c, omega1, omega2, xy + j * 2 );
    }

    theta = circle_segment_angle_from_chord_angles ( omega1, omega2 );
    area = circle_segment_area_from_angle ( r, theta );
    area_est = 4.0 * ( double ) ( i4vec_sum ( n, inout ) ) / ( double ) ( n );

    printf ( "  %6d  %14.6g  %14.6g  %14.6g  %14.6g\n",
      n, omega1, omega2, area, area_est );

    free ( xy );
  }

  free ( inout );

  return;
}
/******************************************************************************/

void test09 ( )

/******************************************************************************/
/*
  Purpose:

    CIRCLE_SEGMENT_TEST09 looks at the area and centroid calculations.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    06 July 2013

  Author:

    John Burkardt
*/
{
  double a1;
  double a2;
  double a3;
  double c[2];
  double *d1;
  double *d2;
  double *d3;
  double h;
  int n;
  double omega1;
  double omega2;
  double p1[2];
  double p2[2];
  const double pi = 3.141592653589793;
  double r;
  int seed;
  double theta;

  printf ( "\n" );
  printf ( "CIRCLE_SEGMENT_TEST09\n" );
  printf ( "  CIRCLE_SEGMENT_AREA_FROM_CHORD and\n" );
  printf ( "  CIRCLE_SEGMENT_CENTROID_FROM_CHORD evaluate the area\n" );
  printf ( "  and centroid of a circle segment, given R, C and P1:P2.\n" );
  printf ( "\n" );
  printf ( "  CIRCLE_SEGMENT_AREA_FROM_SAMPLE and\n" );
  printf ( "  CIRCLE_SEGMENT_CENTROID_FROM_SAMPLE give us Monte Carlo estimates.\n" );
  printf ( "\n" );
  printf ( "  GQCIRCSEGM can estimate these values by quadrature.\n" );
  printf ( "\n" );
  printf ( "  Start easy, with R = 1, C = (0,0), and Theta centered.\n" );

  seed = 123456789;
  r = 1.0;
  c[0] = 0.0;
  c[1] = 0.0;
  theta = pi / 4.0;
  h = circle_segment_height_from_angle ( r, theta );
  omega1 = - theta / 2.0;
  omega2 = + theta / 2.0;
  p1[0] = c[0] + r * cos ( omega1 );
  p1[1] = c[1] + r * sin ( omega1 );
  p2[0] = c[0] + r * cos ( omega2 );
  p2[1] = c[1] + r * sin ( omega2 );

  a1 = circle_segment_area_from_chord ( r, c, p1, p2 );
  d1 = circle_segment_centroid_from_chord ( r, c, p1, p2 );

  printf ( "\n" );
  printf ( "         Area          CentroidX    CentroidY\n" );
  printf ( "\n" );
  printf ( "  %14.6g  %14.6g  %14.6g\n", a1, d1[0], d1[1] );
/*
  This only works because the centroid of the height-based circle segment 
  is easily transformed to the centroid of the chord based circle segment.
*/
  a2 = circle_segment_area_from_height ( r, h );
  d2 = circle_segment_centroid_from_height ( r, h );
  printf ( "  %14.6g  %14.6g  %14.6g\n", a2, d2[1], -d2[0] );

  n = 10000;
  a3 = circle_segment_area_from_sample ( r, c, p1, p2, n, &seed );
  d3 = circle_segment_centroid_from_sample ( r, c, p1, p2, n, &seed );
  printf ( "  %14.6g  %14.6g  %14.6g\n", a3, d3[0], d3[1] );

  free ( d1 );
  free ( d2 );
  free ( d3 );

  return;
}
/******************************************************************************/

void test11 ( )

/******************************************************************************/
/*
  Purpose:

    TEST11 demonstrates CIRCLE_SEGMENT_ROTATION_FROM_CHORD.

  Discussion:

    We make a table of all pairs of angles that are multiples of pi/12.

    For each pair, we compute the rotation, that is, the angle of the
    central radius of the circle segment.  We print out the result in
    terms of multiples of pi/12.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 July 2013

  Author:

    John Burkardt
*/
{
  double alpha;
  double c[2];
  int i;
  int j;
  double p1[2];
  double p2[2];
  double pi = 3.141592653589793;
  double r;
  double rho1;
  double rho2;
  double t;

  printf ( "\n" );
  printf ( "TEST11:\n" );
  printf ( "  CIRCLE_SEGMENT_ROTATION_FROM_CHORD is given the endpoints\n" );
  printf ( "  of a chord, and is asked to determine the angle of the\n" );
  printf ( "  central radius vector.\n" );
  printf ( "\n" );
  printf ( "  We make a table of all pairs of angles that are multiples\n" );
  printf ( "  of pi/12, determine the corresponding chord endpoints, and\n" );
  printf ( "  compute the rotation angle, also printed as a multiple of pi/12.\n" );

  r = 2.0;
  c[0] = 3.0;
  c[1] = 5.0;
  printf ( "\n" );
  printf ( "     0.0   1.0   2.0   3.0   4.0   5.0   6.0   7.0" );
  printf ( "   8.0   9.0  10.0  11.0  12.0\n" );
  printf ( "\n" );
  for ( i = 0; i <= 12; i++ )
  {
    rho1 = ( double ) ( i ) * pi / 6.0;
    p1[0] = c[0] + r * cos ( rho1 );
    p1[1] = c[1] + r * sin ( rho1 );
    printf ( "%2d", i );
    for ( j = 0; j <= 12; j++ )
    {
      rho2 = ( double ) ( j ) * pi / 6.0;
      p2[0] = c[0] + r * cos ( rho2 );
      p2[1] = c[1] + r * sin ( rho2 );
      alpha = circle_segment_rotation_from_chord ( r, c, p1, p2 );
      t = 6.0 * alpha / pi;
      printf ( "  %4.1f", t );
    }
    printf ( "\n" );
  }
  return;
}
/******************************************************************************/

void test13 ( )

/******************************************************************************/
/*
  Purpose:

    TEST13 demonstrates GAUSS for quadrature rules.

  Discussion:

    Some recursion coefficients ALPHA and BETA are listed in Kautsky
    and Elhay.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 July 2013

  Author:

    John Burkardt

  Reference

    Jaroslav Kautsky, Sylvan Elhay,
    Calculation of the Weights of Interpolatory Quadratures,
    Numerische Mathematik,
    Volume 40, Number 3, October 1982, pages 407-422.
*/
{
  double *alpha;
  double *beta;
  int i;
  int n;
  double pi = 3.141592653589793;
  double *w;
  double *x;

  printf ( "\n" );
  printf ( "TEST13\n" );
  printf ( "  GAUSS computes the points and weights for a\n" );
  printf ( "  Gauss quadrature rule, given the ALPHA and BETA\n" );
  printf ( "  recursion coefficients.\n" );
/*
  Legendre rule.
*/
  n = 10;

  alpha = ( double * ) malloc ( n * sizeof ( double ) );
  beta = ( double * ) malloc ( n * sizeof ( double ) );

  for ( i = 0; i < n; i++ )
  {
    alpha[i] = 0.0;
    if ( i == 0 )
    {
      beta[i] = 2.0;
    }
    else
    {
      beta[i] = 1.0 / ( 4.0 - 1.0 / ( double ) ( i * i ) );
    }
  }

  x = ( double * ) malloc ( n * sizeof ( double ) );
  w = ( double * ) malloc ( n * sizeof ( double ) );

  gauss ( n, alpha, beta, x, w );

  printf ( "\n" );
  printf ( "  LEGENDRE RULE\n" );
  printf ( "  Point   Weight\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %14.6g  %14.6g\n", x[i], w[i] );
  }

  free ( alpha );
  free ( beta );
  free ( w );
  free ( x );
/*
  Hermite rule.
*/
  n = 10;

  alpha = ( double * ) malloc ( n * sizeof ( double ) );
  beta = ( double * ) malloc ( n * sizeof ( double ) );

  for ( i = 0; i < n; i++ )
  {
    alpha[i] = 0.0;
    if ( i == 0 )
    {
      beta[i] = sqrt ( pi );
    }
    else
    {
      beta[i] = ( double ) ( i ) / 2.0;
    }
  }

  x = ( double * ) malloc ( n * sizeof ( double ) );
  w = ( double * ) malloc ( n * sizeof ( double ) );

  gauss ( n, alpha, beta, x, w );

  printf ( "\n" );
  printf ( "  HERMITE RULE\n" );
  printf ( "  Point   Weight\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %14.6g  %14.6g\n", x[i], w[i] );
  }

  free ( alpha );
  free ( beta );
  free ( w );
  free ( x );
/*
  Laguerre rule.
*/
  n = 10;

  alpha = ( double * ) malloc ( n * sizeof ( double ) );
  beta = ( double * ) malloc ( n * sizeof ( double ) );

  for ( i = 0; i < n; i++ )
  {
    alpha[i] = 2.0 * ( double ) ( i + 1 ) - 1.0;
    if ( i == 0 )
    {
      beta[i] = 1.0;
    }
    else
    {
      beta[i] = ( double ) ( i * i );
    }
  }

  x = ( double * ) malloc ( n * sizeof ( double ) );
  w = ( double * ) malloc ( n * sizeof ( double ) );

  gauss ( n, alpha, beta, x, w );

  printf ( "\n" );
  printf ( "  LAGUERRE RULE\n" );
  printf ( "  Point   Weight\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %14.6g  %14.6g\n", x[i], w[i] );
  }

  free ( alpha );
  free ( beta );
  free ( w );
  free ( x );

  return;
}
/******************************************************************************/

void test14 ( )

/******************************************************************************/
/*
  Purpose:

    TEST14 demonstrates R_JACOBI.

  Discussion:

    R_JACOBI returns recursion coefficients ALPHA and BETA for rules
    using a Jacobi type weight w(x) = (1-x)^A * (1+x)^B.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 July 2013

  Author:

    John Burkardt

  Reference

    Walter Gautschi,
    Orthogonal Polynomials: Computation and Approximation,
    Oxford, 2004,
    ISBN: 0-19-850672-4,
    LC: QA404.5 G3555.
*/
{
  double a;
  double *alpha;
  double b;
  double *beta;
  int i;
  int n;

  printf ( "\n" );
  printf ( "TEST14\n" );
  printf ( "  R_JACOBI computes recursion coefficients ALPHA and BETA\n" );
  printf ( "  Gauss quadrature rule, given the ALPHA and BETA\n" );
  printf ( "  recursion coefficients.\n" );
/*
  Legendre rule.
*/
  n = 10;

  a = 0.0;
  b = 0.0;
  alpha = ( double * ) malloc ( n * sizeof ( double ) );
  beta = ( double * ) malloc ( n * sizeof ( double ) );

  r_jacobi ( n, a, b, alpha, beta );

  printf ( "\n" );
  printf ( "  Legendre weight\n" );
  printf ( "  A = %g,  B = %g\n", a, b );
  printf ( "  Alpha          Beta\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %14.6g  %14.6g\n", alpha[i], beta[i] );
  }
  free ( alpha );
  free ( beta );
/*
  Chebyshev Type 1 rule.
*/
  n = 10;

  a = -0.5;
  b = -0.5;
  alpha = ( double * ) malloc ( n * sizeof ( double ) );
  beta = ( double * ) malloc ( n * sizeof ( double ) );

  r_jacobi ( n, a, b, alpha, beta );

  printf ( "\n" );
  printf ( "  Chebyshev Type 1 weight\n" );
  printf ( "  A = %g,  B = %g\n", a, b );
  printf ( "  Alpha          Beta\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %14.6g  %14.6g\n", alpha[i], beta[i] );
  }
  free ( alpha );
  free ( beta );
/*
  Chebyshev Type 2 rule.
*/
  n = 10;

  a = +0.5;
  b = +0.5;
  alpha = ( double * ) malloc ( n * sizeof ( double ) );
  beta = ( double * ) malloc ( n * sizeof ( double ) );

  r_jacobi ( n, a, b, alpha, beta );

  printf ( "\n" );
  printf ( "  Chebyshev Type 2 weight\n" );
  printf ( "  A = %g,  B = %g\n", a, b );
  printf ( "  Alpha          Beta\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %14.6g  %14.6g\n", alpha[i], beta[i] );
  }
  free ( alpha );
  free ( beta );
/*
  General Jacobi rule.
*/
  n = 10;

  a = +0.5;
  b = +1.5;
  alpha = ( double * ) malloc ( n * sizeof ( double ) );
  beta = ( double * ) malloc ( n * sizeof ( double ) );

  r_jacobi ( n, a, b, alpha, beta );

  printf ( "\n" );
  printf ( "  General Jacobi weight\n" );
  printf ( "  A = %g,  B = %g\n", a, b );
  printf ( "  Alpha          Beta\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %14.6g  %14.6g\n", alpha[i], beta[i] );
  }
  free ( alpha );
  free ( beta );

  return;
}
