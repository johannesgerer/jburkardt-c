# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <string.h>
# include <complex.h>

# include "subpak.h"

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
void test12 ( );
void test13 ( );
void test14 ( );
void test15 ( );
void test16 ( );
void test17 ( );
void test18 ( );
void test19 ( );

void test20 ( );
void test21 ( );
void test22 ( );
void test225 ( );
void test23 ( );
void test24 ( );
void test25 ( );
void test26 ( );
void test27 ( );
void test29 ( );

void test30 ( );
void test31 ( );
void test32 ( );
void test33 ( );
void test34 ( );
void test35 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for SUBPAK_PRB.

  Discussion:

    SUBPAK_PRB tests the SUBPAK library.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    03 November 2011

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "SUBPAK_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the SUBPAK routines.\n" );

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
  test12 ( );
  test13 ( );
  test14 ( );
  test15 ( );
  test16 ( );
  test17 ( );
  test18 ( );
  test19 ( );

  test20 ( );
  test21 ( );
  test22 ( );
  test225 ( );
  test23 ( );
  test24 ( );
  test25 ( );
  test26 ( );
  test27 ( );
  test29 ( );

  test30 ( );
  test31 ( );
  test32 ( );
  test33 ( );
  test34 ( );
  test35 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "SUBPAK_PRB\n" );
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

    TEST01 tests ANGLE_SHIFT.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    20 July 2010

  Author:

    John Burkardt
*/
{
  double alpha;
  double angle_hi;
  double angle_lo;
  double beta;
  double gamma;
  double pi = 3.141592653589793;
  int seed;
  int test;
  int test_num = 10;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  ANGLE_SHIFT shifts an angle by multiples of\n" );
  printf ( "  2 Pi until it lies between BETA and BETA+2Pi.\n" );
  printf ( "\n" );
  printf ( "     ALPHA      BETA     GAMMA   BETA+2Pi\n" );
  printf ( "\n" );

  angle_lo = -4.0 * pi;
  angle_hi = +4.0 * pi;

  seed = 123456789;

  for ( test = 1; test <= test_num; test++ )
  {
    alpha = r8_uniform ( angle_lo, angle_hi, &seed );

    beta = r8_uniform ( angle_lo, angle_hi, &seed );

    gamma = angle_shift ( alpha, beta );

    printf ( "  %8.2f  %8.2f  %8.2f  %8.2f\n", 
      alpha, beta, gamma, beta + 2.0 * pi );
  }

  return;
}
/******************************************************************************/

void test02 ( )

/******************************************************************************/
/*
  Purpose:

    TEST02 tests ANGLE_SHIFT_DEG.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    20 July 2010

  Author:

    John Burkardt
*/
{
  double alpha;
  double angle_hi;
  double angle_lo;
  double beta;
  double gamma;
  int seed;
  int test;
  int test_num = 10;

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  ANGLE_SHIFT_DEG shifts an angle by multiples of\n" );
  printf ( "  360 until it lies between BETA and BETA+360.\n" );
  printf ( "\n" );
  printf ( "     ALPHA      BETA     GAMMA   BETA+360\n" );
  printf ( "\n" );

  angle_lo = -720.0;
  angle_hi = +720.0;

  seed = 123456789;

  for ( test = 1; test <= test_num; test++ )
  {
    alpha = r8_uniform ( angle_lo, angle_hi, &seed );

    beta = r8_uniform ( angle_lo, angle_hi, &seed );

    gamma = angle_shift_deg ( alpha, beta );

    printf ( "  %8.2f  %8.2f  %8.2f  %8.2f\n", 
      alpha, beta, gamma, beta + 360.0 );
  }

  return;
}
/******************************************************************************/

void test03 ( )

/******************************************************************************/
/*
  Purpose:

    TEST03 tests ANGLE_TO_RGB.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    05 May 2006

  Author:

    John Burkardt
*/
{
  double angle;
  double angle_lo = 0.0;
  double angle_hi = 360.0;
  int i;
  double *rgb;
  int seed;
  int test;
  int test_num = 10;

  printf ( "\n" );
  printf ( "TEST03\n" );
  printf ( "  ANGLE_TO_RGB converts an angle into an RGB color.\n" );
  printf ( "\n" );
  printf ( "     ANGLE         R         G         B\n" );
  printf ( "\n" );

  seed = 123456789;

  for ( test = 0; test < test_num; test++ )
  {
    angle = r8_uniform ( angle_lo, angle_hi, &seed );

    rgb = angle_to_rgb ( angle );

    printf ( "  %8.2f  %8.2f  %8.2f  %8.2f\n", angle, rgb[0], rgb[1], rgb[2] );

    free ( rgb );
  }

  return;
}
/******************************************************************************/

void test04 ( )

/******************************************************************************/
/*
  Purpose:

    TEST04 tests AXIS_LIMITS.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    21 July 2010

  Author:

    John Burkardt
*/
{
  int ndivs;
  int nticks;
  double pxdiv;
  double pxmax;
  double pxmin;
  double xmax;
  double xmin;

  xmin = 67.3;
  xmax = 114.7;
  ndivs = 6;

  axis_limits ( xmin, xmax, ndivs, &pxmin, &pxmax, &pxdiv, &nticks );

  printf ( "\n" );
  printf ( "TEST04\n" );
  printf ( "  AXIS_LIMITS adjusts plot limits to \"nicer\" values.\n" );
  printf ( "\n" );
  printf ( "  Input XMIN =    %f\n", xmin );
  printf ( "  Input XMAX =    %f\n", xmax );
  printf ( "  Input NDIVS =   %d\n", ndivs );
  printf ( "\n" );
  printf ( "  Output PXMIN =  %f\n", pxmin );
  printf ( "  Output PXMAX =  %f\n", pxmax );
  printf ( "  Output PXDIV =  %f\n", pxdiv );
  printf ( "  Output NTICKS = %d\n", nticks );

  xmin = -26.0;
  xmax = +26.0;
  ndivs = 10;

  axis_limits ( xmin, xmax, ndivs, &pxmin, &pxmax, &pxdiv, &nticks );

  printf ( "\n" );
  printf ( "  Input XMIN =    %f\n", xmin );
  printf ( "  Input XMAX =    %f\n", xmax );
  printf ( "  Input NDIVS =   %d\n", ndivs );
  printf ( "\n" );
  printf ( "  Output PXMIN =  %f\n", pxmin );
  printf ( "  Output PXMAX =  %f\n", pxmax );
  printf ( "  Output PXDIV =  %f\n", pxdiv );
  printf ( "  Output NTICKS = %d\n", nticks );

  return;
}
/******************************************************************************/

void test05 ( )

/******************************************************************************/
/*
  Purpose:

    TEST05 tests AXIS_LIMITS.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    21 July 2010

  Author:

    John Burkardt
*/
{
# define TEST_NUM 5

  int i;
  int ndivs;
  int nticks;
  double pxdiv;
  double pxmax;
  double pxmin;
  double test_max[TEST_NUM] = {
    9.0, 4.125, 193.75, 2000.250, 12.0 };
  double test_min[TEST_NUM] = {
    1.0, 1.003, 101.25, 2000.125, -7.0 };
  double xmax;
  double xmin;

  printf ( "\n" );
  printf ( "TEST05\n" );
  printf ( "  AXIS_LIMITS computes \"nice\" limits for a graph\n" );
  printf ( "    that must include a given range.\n" );

  ndivs = 5;

  printf ( "\n" );
  printf ( "  All tests use NDIVS = %d\n", ndivs );
  printf ( "\n" );
  printf ( "          XMIN          XMAX         PXMIN" );
  printf ( "         PXMAX         PXDIV  NTICKS\n" );
  printf ( "\n" );

  for ( i = 0; i < TEST_NUM; i++ )
  {
    xmin = test_min[i];
    xmax = test_max[i];

    axis_limits ( xmin, xmax, ndivs, &pxmin, &pxmax, &pxdiv, &nticks );

    printf ( "  %12f  %12f  %12f  %12f  %12f  %6d\n",
      xmin, xmax, pxmin, pxmax, pxdiv, nticks );
  }

  return;
# undef TEST_NUM
}
/******************************************************************************/

void test06 ( )

/******************************************************************************/
/*
//  Purpose:
//
//    TEST06 tests BAR_CHECK, BAR_CODE, BAR_DIGIT_CODE_LEFT, BAR_DIGIT_CODE_RIGHT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 July 2010
//
//  Author:
//
//    John Burkardt
*/
{
  char *bar;
  int check;
  int digit[12];
  char *codel;
  char *coder;
  int i;

  printf ( "\n" );
  printf ( "TEST06\n" );
  printf ( "  BAR_CHECK checks digits for a barcode;\n" );
  printf ( "  BAR_CODE computes the barcode for a string of 11 digits;\n" );
  printf ( "  BAR_DIGIT_CODE_LEFT returns the left digit code.\n" );
  printf ( "  BAR_DIGIT_CODE_RIGHT returns the right digit code.\n" );

  for ( i = 0; i <= 10; i++ )
  {
    digit[i] = ( i % 10 );
  }
 
  check = bar_check ( digit );
 
  printf ( "\n" );
  printf ( "  The check digit is %d\n", check );

  digit[11] = check;
 
  printf ( "\n" );
  printf ( "  The left and right digit codes:\n" );
  printf ( "\n" );

  for( i = 0; i <= 9; i++ )
  {
    codel = bar_digit_code_left ( i );
    coder = bar_digit_code_right ( i );
    printf ( "  %2d  %s  %s\n", i, codel, coder );
    free ( codel );
    free ( coder );
  }
 
  bar = bar_code ( digit );
 
  printf ( "\n" );
  printf ( "  Bar code:\n" );
  printf ( "\n" );
  for ( i = 0; i < 9; i++ )
  {
    printf ( "%c", bar[i] );
  }
  printf ( "\n" );
  for ( i = 9; i < 12; i++ )
  {
    printf ( "%c", bar[i] );
  }
  printf ( "\n" );
  for ( i = 12; i < 19; i++ )
  {
    printf ( "%c", bar[i] );
  }
  printf ( "\n" );

  for ( i = 19; i < 26; i++ )
  {
    printf ( "%c", bar[i] );
  }
  printf ( "\n" );
  for ( i = 16; i < 33; i++ )
  {
    printf ( "%c", bar[i] );
  }
  printf ( "\n" );
  for ( i = 33; i < 40; i++ )
  {
    printf ( "%c", bar[i] );
  }
  printf ( "\n" );
  for ( i = 40; i < 47; i++ )
  {
    printf ( "%c", bar[i] );
  }
  printf ( "\n" );
  for ( i = 47; i < 54; i++ )
  {
    printf ( "%c", bar[i] );
  }
  printf ( "\n" );
  for ( i = 54; i < 59; i++ )
  {
    printf ( "%c", bar[i] );
  }
  printf ( "\n" );
  for ( i = 59; i < 66; i++ )
  {
    printf ( "%c", bar[i] );
  }
  printf ( "\n" );
  for ( i = 66; i < 73; i++ )
  {
    printf ( "%c", bar[i] );
  }
  printf ( "\n" );
  for ( i = 73; i < 80; i++ )
  {
    printf ( "%c", bar[i] );
  }
  printf ( "\n" );
  for ( i = 80; i < 87; i++ )
  {
    printf ( "%c", bar[i] );
  }
  printf ( "\n" );
  for ( i = 87; i < 94; i++ )
  {
    printf ( "%c", bar[i] );
  }
  printf ( "\n" );
  for ( i = 94; i < 101; i++ )
  {
    printf ( "%c", bar[i] );
  }
  printf ( "\n" );
  for ( i = 101; i < 104; i++ )
  {
    printf ( "%c", bar[i] );
  }
  printf ( "\n" );
  for ( i = 104; i < 113; i++ )
  {
    printf ( "%c", bar[i] );
  }
  printf ( "\n" );

  free ( bar );

  return;
}
/******************************************************************************/

void test07 ( )

/******************************************************************************/
/*
  Purpose:

    TEST07 tests BMI_ENGLISH.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    21 July 2010

  Author:

    John Burkardt
*/
{
# define TEST_NUM 10

  double b;
  double bmi;
  double c;
  double h;
  double h_ft;
  double h_in;
  int seed;
  int test;
  double w;

  printf ( "\n" );
  printf ( "TEST07\n" );
  printf ( "  BMI_ENGLISH computes the Body Mass Index\n" );
  printf ( "  given body measurements in English Units.\n" );
  printf ( "\n" );
  printf ( "      Weight               Height            BMI\n" );
  printf ( "       (LB)          (FT          IN)\n" );
  printf ( "\n" );

  seed = 123456789;

  for ( test = 0; test < TEST_NUM; test++ )
  {
    b = 100.0;
    c = 250.0;

    w = r8_uniform ( b, c, &seed );

    b = 4.0;
    c = 6.75;

    h = r8_uniform ( b, c, &seed );

    h_ft = ( int ) ( h );
    h_in = ( double ) ( ( int ) ( 12.0 * ( h - h_ft ) ) );
 
    bmi = bmi_english ( w, h_ft, h_in );
    printf ( "  %10f  %10f  %10f  %10f\n", w, h_ft, h_in, bmi );
  }

  return;
# undef TEST_NUM
}
/******************************************************************************/

void test08 ( )

/******************************************************************************/
/*
  Purpose:

    TEST08 tests FAC_DIV, FAC_GCD, FAC_LCM, FAC_MUL, FAC_TO_I4, I4_TO_FAC.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    21 July 2010

  Author:

    John Burkardt
*/
{
# define PRIME_NUM 5

  int bot;
  int i1;
  int i2;
  int *npower1;
  int *npower2;
  int npower3[PRIME_NUM];
  int top;

  printf ( "\n" );
  printf ( "TEST08\n" );
  printf ( "  For products of prime factors:\n" );
  printf ( "  FAC_DIV computes a quotient;\n" );
  printf ( "  FAC_MUL multiplies;\n" );
  printf ( "  FAC_LCM computes the LCM;\n" );
  printf ( "  FAC_GCD computes the GCD;\n" );
  printf ( "  I4_TO_FAC converts an integer;\n" );
  printf ( "  FAC_TO_I4 converts to an integer.\n" );
  printf ( "  FAC_TO_RAT converts to a ratio.\n" );

  i1 = 720;
  i2 = 42;

  npower1 = i4_to_fac ( i1, PRIME_NUM );

  printf ( "\n" );
  printf ( "  Representation of I1 = %d\n", i1 );
  printf ( "\n" );

  fac_print ( PRIME_NUM, npower1 );

  npower2 = i4_to_fac ( i2, PRIME_NUM );

  printf ( "\n" );
  printf ( "  Representation of I2 = %d\n", i2 );
  printf ( "\n" );

  fac_print ( PRIME_NUM, npower2 );

  fac_lcm ( PRIME_NUM, npower1, npower2, npower3 );

  printf ( "\n" );
  printf ( "  LCM of I1, I2:\n" );
  printf ( "\n" );

  fac_print ( PRIME_NUM, npower3 );

  fac_gcd ( PRIME_NUM, npower1, npower2, npower3 );

  printf ( "\n" );
  printf ( "  GCD of I1, I2:\n" );
  printf ( "\n" );

  fac_print ( PRIME_NUM, npower3 );

  fac_mul ( PRIME_NUM, npower1, npower2, npower3 );

  printf ( "\n" );
  printf ( "  Product of I1, I2:\n" );
  printf ( "\n" );

  fac_print ( PRIME_NUM, npower3 );

  fac_div ( PRIME_NUM, npower2, npower1, npower3 );

  printf ( "\n" );
  printf ( "  Quotient of I2 / I1:\n" );
  printf ( "\n" );

  fac_print ( PRIME_NUM, npower3 );

  fac_to_rat ( PRIME_NUM, npower3, &top, &bot );

  printf ( "\n" );
  printf ( "  Quotient as a rational: %d / %d\n", top, bot );

  free ( npower1 );
  free ( npower2 );

  return;
# undef PRIME_NUM
}
/******************************************************************************/

void test09 ( )

/******************************************************************************/
/*
  Purpose:

    TEST09 tests GAUSS_SUM

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    21 July 2010

  Author:

    John Burkardt
*/
{
# define DIM_NUM 2
# define N 3

  double amplitude[N] = { 10.0, 5.0, -3.0 };
  double center[DIM_NUM*N] = { 2.0, 3.0,  5.0, 8.0,  7.0, 5.0 };
  double gxy;
  int i;
  int j;
  double width[N] = { 1.0, 2.0, 4.0 };
  double x[DIM_NUM];

  printf ( "\n" );
  printf ( "TEST09\n" );
  printf ( "  GAUSS_SUM evaluates a function which is the sum of\n" );
  printf ( "  Gaussian functions.\n" );
  printf ( "\n" );
  printf ( "  Number of component Gaussians = %d\n", N );
  printf ( "\n" );
  printf ( "          Center    Amplitude  Width\n" );
  printf ( "        X       Y\n" );
  printf ( "\n" );
  for ( j = 0; j < N; j++ )
  {
    printf ( "  %2d  %6.2f  %6.2f  %6.2f  %6.2f\n", 
      j, center[0+j*DIM_NUM], center[1+j*DIM_NUM], amplitude[j], width[j] );
  }

  printf ( "\n" );
  printf ( "      X       Y        Gauss_Sum(X,Y)\n" );
  printf ( "\n" );

  for ( i = 0; i <= 10; i++ )
  {
    x[0] = ( double ) i;
    for ( j = 0; j <= 10; j++ )
    {
      x[1] = ( double ) j;
      gxy = gauss_sum ( DIM_NUM, N, amplitude, center, width, x );
      printf ( "  %6.2f  %6.2f  %14.6f\n", x[0], x[1], gxy );
    }
  }

  return;
# undef DIM_NUM
# undef N
}
/******************************************************************************/

void test10 ( )

/******************************************************************************/
/*
  Purpose:

    TEST10 tests GET_SEED.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    22 July 2010

  Author:

    John Burkardt
*/
{
# define TEST_NUM 10

  int seed;
  int seed_0;
  int seed_1;
  int seed_2;
  int seed_3;
  int test;
  double x;

  printf ( "\n" );
  printf ( "TEST10\n" );
  printf ( "  GET_SEED gets a seed for the random number\n" );
  printf ( "    generator.  These values are computed from\n" );
  printf ( "    the time and date.  Values computed nearby\n" );
  printf ( "    in time will be near to each other, and\n" );
  printf ( "    should be passed through a random number\n" );
  printf ( "    generator a few times before use.\n" );
  printf ( "\n" );
  printf ( "     I	     R(I)	 R2(I)        R3(I)\n" );
  printf ( "\n" );

  for ( test = 0; test < TEST_NUM; test++ )
  {
    seed = ( int ) get_seed ( );
    seed_0 = seed;
    x = r8_uniform_01 ( &seed );
    seed_1 = seed;
    x = r8_uniform_01 ( &seed );
    seed_2 = seed;
    x = r8_uniform_01 ( &seed );
    seed_3 = seed;

    printf ( "  %12d  %12d  %12d  %12d\n", seed_0, seed_1, seed_2, seed_3 );
  }

  return;
# undef TEST_NUM
}
/******************************************************************************/

void test11 ( )

/******************************************************************************/
/*
  Purpose:

    TEST11 tests GRID1.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    22 July 2010

  Author:

    John Burkardt
*/
{
# define DIM_NUM 5
# define NSTEP 11

  int i;
  double *x;
  double x1[DIM_NUM] = { 1.0,  0.0, 20.0, -5.0, 1.0 };
  double x2[DIM_NUM] = { 1.0, 10.0,  0.0,  5.0, 2.0 };

  printf ( "\n" );
  printf ( "TEST11\n" );
  printf ( "  GRID1 computes a 1D grid between\n" );
  printf ( "    two DIM_NUM dimensional points X1 and X2.\n" );
  printf ( "\n" );
  printf ( "  Here, we will use %d steps\n", NSTEP );
  printf ( "  going from: \n" );
  for ( i = 0; i < DIM_NUM; i++ )
  {
    printf ( "  %f", x1[i] );
  }
  printf ( "\n" );
  printf ( "  to:\n" );
  for ( i = 0; i < DIM_NUM; i++ )
  {
    printf ( "  %f", x2[i] );
  }
  printf ( "\n" );
  printf ( "\n" );

  x = grid1 ( DIM_NUM, NSTEP, x1, x2 );

  r8mat_transpose_print ( DIM_NUM, NSTEP, x, "  The grid matrix:" );

  free ( x );

  return;
# undef DIM_NUM
# undef NSTEP
}
/******************************************************************************/

void test12 ( )

/******************************************************************************/
/*
  Purpose:

    TEST12 tests GRID1N.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    23 July 2010

  Author:

    John Burkardt
*/
{
# define DIM_NUM 5
# define NSTEP 11

  int i;
  int j;
  double *x;
  double x1[DIM_NUM] = { 1.0,  0.0, 20.0, -5.0, 1.0 };
  double x2[DIM_NUM] = { 1.0, 10.0,  0.0,  5.0, 2.0 };

  printf ( "\n" );
  printf ( "TEST12\n" );
  printf ( "  GRID1N computes a 1D grid between\n" );
  printf ( "    two DIM_NUM dimensional points X1 and X2,\n" );
  printf ( "    one point at a time.\n" );
  printf ( "\n" );
  printf ( "  Here, we will use %d steps\n", NSTEP );
  printf ( "  going from \n" );
  for ( i = 0; i < DIM_NUM; i++ )
  {
    printf ( "  %12.4f", x1[i] );
  }
  printf ( "\n" );
  printf ( "  to\n" );
  for ( i = 0; i < DIM_NUM; i++ )
  {
    printf ( "  %12.4f", x2[i] );
  }
  printf ( "\n" );
  printf ( "\n" );
 
  for ( j = 1; j <= NSTEP; j++ )
  {
    x = grid1n ( j, DIM_NUM, NSTEP, x1, x2 );

    printf ( "  %6d  ", j );
    for ( i = 0; i < DIM_NUM; i++ )
    {
      printf ( "  %12.4f", x[i] );
    }
    printf ( "\n" );
    free ( x );
  }
 
  return;
# undef DIM_NUM
# undef NSTEP
}
/******************************************************************************/

void test13 ( )

/******************************************************************************/
/*
  Purpose:

    TEST13 tests GRID2.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    23 July 2010

  Author:

    John Burkardt
*/
{
# define DIM_NUM 5
# define NSTEP 20

  int i;
  int j1;
  int j2;
  double *x;
  double x1[DIM_NUM] = { 1.0,  0.0, 20.0, -5.0, 1.0 };
  double x2[DIM_NUM] = { 1.0, 10.0,  0.0,  5.0, 2.0 };

  j1 = 3;
  j2 = 13;

  printf ( "\n" );
  printf ( "TEST13\n" );
  printf ( "  GRID2 computes a 1 D grid between\n" );
  printf ( "    two DIM_NUM dimensional points X1 and X2,\n" );
  printf ( "    computing X1 and X2 at user specified times.\n" );
  printf ( "\n" );
  printf ( "  Here, we will use %d steps.\n", NSTEP );
  printf ( "  and on step %d we will compute\n", j1 );
  for ( i = 0; i < DIM_NUM; i++ )
  {
    printf ( "  %10.4f", x1[i] );
  }
  printf ( "\n" );
  printf ( "  and on step %d we will compute\n", j2 );
  for ( i = 0; i < DIM_NUM; i++ )
  {
    printf ( "  %10.4f", x2[i] );
  }
  printf ( "\n" );
  printf ( "\n" );
 
  x = grid2 ( j1, j2, DIM_NUM, NSTEP, x1, x2 );
 
  r8mat_print ( DIM_NUM, NSTEP, x, "  The grid matrix:" );

  free ( x );

  return;
# undef DIM_NUM
# undef NSTEP
}
/******************************************************************************/

void test14 ( )

/******************************************************************************/
/*
  Purpose:

    TEST14 tests GRID2N.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    23 July 2010

  Author:

    John Burkardt
*/
{
# define DIM_NUM 5

  int i;
  int j;
  int j1;
  int j2;
  double *x;
  double x1[DIM_NUM] = { 1.0,  0.0, 20.0, -5.0, 1.0 };
  double x2[DIM_NUM] = { 1.0, 10.0,  0.0,  5.0, 2.0 };

  j1 = 3;
  j2 = 13;

  printf ( "\n" );
  printf ( "TEST14\n" );
  printf ( "  GRID2N computes points from a 1D grid\n" );
  printf ( "    between two DIM_NUM dimensional points\n" );
  printf ( "    X1 and X2, one at a time, with X1 and X2\n" );
  printf ( "    having user specified J coordinates.\n" );
  printf ( "\n" );
  printf ( "  On step %d we will compute\n", j1 );
  for ( i = 0; i < DIM_NUM; i++ )
  {
    printf ( "  %10f", x1[i] );
  }
  printf ( "\n" );
  printf ( "  and on step %d we will compute\n", j2 );
  for ( i = 0; i < DIM_NUM; i++ )
  {
    printf ( "  %10.4f", x2[i] );
  }
  printf ( "\n" );
  printf ( "\n" );
 
  for ( j = 1; j <= 20; j++ )
  {
    x = grid2n ( j, j1, j2, DIM_NUM, x1, x2 );

    printf ( "  %6d", j );

    for ( i = 0; i < DIM_NUM; i++ )
    {
      printf ( "  %10.4f", x[i] );
    }
    printf ( "\n" );

    free ( x );
  }
 
  return;
# undef DIM_NUM
}
/******************************************************************************/

void test15 ( )

/******************************************************************************/
/*
  Purpose:

    TEST15 tests GRID3.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    25 July 2010

  Author:

    John Burkardt
*/
{
# define DIM_NUM 5

  int i;
  int j;
  int k;
  int nstep1 = 3;
  int nstep2 = 6;
  double *x;
  double x1[DIM_NUM] = { 1.0,  0.0, 20.0, -5.0, 1.0 };
  double x2[DIM_NUM] = { 1.0, 10.0,  0.0,  5.0, 2.0 };
  double x3[DIM_NUM] = { 1.0,  5.0,  0.0,  0.0, 3.0 };
 
  printf ( "\n" );
  printf ( "TEST15\n" );
  printf ( "  GRID3 computes a 2D grid in the plane\n" );
  printf ( "  containing the DIM_NUM-dimensional\n" );
  printf ( "  points X1, X2 and X3.\n" );
  printf ( "\n" );
  printf ( "  Here, we will use %d steps\n", nstep1 );
  printf ( "  going from \n" );
  for ( i = 0; i < DIM_NUM; i++ )
  {
    printf ( "  %10f", x1[i] );
  }
  printf ( "\n" );
  printf ( "  to\n" );
  for ( i = 0; i < DIM_NUM; i++ )
  {
    printf ( "  %10f", x2[i] );
  }
  printf ( "\n" );
  printf ( "  and %d steps going to \n", nstep2 );
  for ( i = 0; i < DIM_NUM; i++ )
  {
    printf ( "  %10f", x3[i] );
  }
  printf ( "\n" );
 
  x = grid3 ( DIM_NUM, nstep1, nstep2, x1, x2, x3 );
 
  for ( j = 1; j <= nstep1; j++ )
  {
    printf ( "\n" );
    for ( k = 1; k <= nstep2; k++ )
    { 
      printf ( "  %3d  %3d", j, k );
      for ( i = 0; i < DIM_NUM; i++ )
      {
        printf ( "  %10f", x[i+(j-1)*DIM_NUM+(j-1)*(k-1)*DIM_NUM*nstep1] );
      }
      printf ( "\n" );
    }
  }

  free ( x );
 
  return;
# undef DIM_NUM
}
/******************************************************************************/

void test16 ( )

/******************************************************************************/
/*
  Purpose:

    TEST16 tests GRID3N.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    25 July 2010

  Author:

    John Burkardt
*/
{
# define DIM_NUM 5

  int i;
  int j;
  int k;
  int nstep1 = 3;
  int nstep2 = 6;
  double *x;
  double x1[DIM_NUM] = { 1.0,  0.0, 20.0, -5.0, 1.0 };
  double x2[DIM_NUM] = { 1.0, 10.0,  0.0,  5.0, 2.0 };
  double x3[DIM_NUM] = { 1.0,  5.0,  0.0,  0.0, 3.0 };

  printf ( "\n" );
  printf ( "TEST16\n" );
  printf ( "  GRID3N computes a point from a 2D\n" );
  printf ( "  grid in the plane containing the \n" );
  printf ( "  DIM_NUM-dimensional points X1, X2 and X3.\n" );
  printf ( "\n" );
  printf ( "  We use %d steps from \n", nstep1 );
  for ( i = 0; i < DIM_NUM; i++ )
  {
    printf ( "  %10f", x1[i] );
  }
  printf ( "\n" );
  printf ( "  to\n" );
  for ( i = 0; i < DIM_NUM; i++ )
  {
    printf ( "  %10f", x2[i] );
  }
  printf ( "\n" );
  printf ( "  and %d steps going to \n", nstep2 );
  for ( i = 0; i < DIM_NUM; i++ )
  {
    printf ( "  %10f", x3[i] );
  }
  printf ( "\n" );
 
  for ( j = 1; j <= nstep1; j++ )
  {
    printf ( "\n" );
    for ( k = 1; k <= nstep2; k++ )
    {
      x = grid3n ( j, k, DIM_NUM, nstep1, nstep2, x1, x2, x3 );
      printf ( "  %3d  %3d", j, k );
      for ( i = 0; i < DIM_NUM; i++ )
      {
        printf ( "  %10f", x[i] );
      }
      printf ( "\n" );
      free ( x );
    }
  }
 
  return;
#  undef DIM_NUM
}
/******************************************************************************/

void test17 ( )

/******************************************************************************/
/*
//  Purpose:
//
//    TEST17 tests GRID4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 July 2010
//
//  Author:
//
//    John Burkardt
*/
{
# define DIM_NUM 5

  int i;
  int j;
  int j1;
  int j2;
  int k;
  int k1;
  int k2;
  int nstep1 = 6;
  int nstep2 = 10;
  double *x;
  double x1[DIM_NUM] = { 1.0,  0.0, 20.0, -5.0, 1.0 };
  double x2[DIM_NUM] = { 1.0, 10.0,  0.0,  5.0, 2.0 };
  double x3[DIM_NUM] = { 1.0,  5.0,  0.0,  0.0, 3.0 };

  j1 = 2;
  j2 = 5;
  k1 = 3;
  k2 = 9;

  printf ( "\n" );
  printf ( "TEST17\n" );
  printf ( "  GRID4 computes a 2D planar grid\n" );
  printf ( "  containing the DIM_NUM-dimensional\n" );
  printf ( "  points X1, X2 and X3.\n" );
  printf ( "\n" );
  printf ( "  We compute the points on the following steps:\n" );
  printf ( "\n" );
  printf ( "  X1 on step %d  %d\n", j1, k1 );
  printf ( "  X2 on step %d  %d\n", j2, k1 );
  printf ( "  X3 on step %d  %d\n", j1, k2 );
  printf ( "\n" );
  printf ( "  We use %d steps in the J direction\n", nstep1 );
  printf ( "  and %d steps in the K direction.\n", nstep2 );
  printf ( "\n" );
  printf ( "  The points X1, X2 and X3 are:\n" );
  printf ( "\n" );
  for ( i = 0; i < DIM_NUM; i++ )
  {
    printf ( "  %10f", x1[i] );
  }
  printf ( "\n" );
  printf ( "\n" );
  for ( i = 0; i < DIM_NUM; i++ )
  {
    printf ( "  %10f", x2[i] );
  }
  printf ( "\n" );
  printf ( "\n" );
  for ( i = 0; i < DIM_NUM; i++ )
  {
    printf ( "  %10f", x3[i] );
  }
  printf ( "\n" );
 
  x = grid4 ( j1, j2, k1, k2, DIM_NUM, nstep1, nstep2, x1, x2, x3 );
 
  for ( j = 1; j <= nstep1; j++ )
  {
    printf ( "\n" );
    for ( k = 1; k <= nstep2; k++ )
    { 
      printf ( "  %3d  %3d", j, k );
      for ( i = 0; i < DIM_NUM; i++ )
      {
        printf ( "  %10f", x[i+j*DIM_NUM+k*DIM_NUM*nstep1] );
      }
      printf ( "\n" ); 
    }
  }

  free ( x );

  return;
# undef DIM_NUM
}
/******************************************************************************/

void test18 ( )

/******************************************************************************/
/*
//  Purpose:
//
//    TEST18 tests GRID4N.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 July 2010
//
//  Author:
//
//    John Burkardt
*/
{
# define DIM_NUM 5

  int i;
  int j;
  int j1;
  int j2;
  int k;
  int k1;
  int k2;
  int nstep1 = 6;
  int nstep2 = 10;
  double *x;
  double x1[DIM_NUM] = { 1.0,  0.0, 20.0, -5.0, 1.0 };
  double x2[DIM_NUM] = { 1.0, 10.0,  0.0,  5.0, 2.0 };
  double x3[DIM_NUM] = { 1.0,  5.0,  0.0,  0.0, 3.0 };

  j1 = 2;
  j2 = 5;
  k1 = 3;
  k2 = 9;
  
  printf ( "\n" );
  printf ( "TEST18\n" );
  printf ( "  GRID4N computes, one at a time, points\n" );
  printf ( "  on a 2D grid in the plane containing\n" );
  printf ( "  the DIM_NUM-dimensional points X1, X2 and X3.\n" );
  printf ( "\n" );
  printf ( "  We compute the points on the following steps:\n" );
  printf ( "\n" );
  printf ( "  X1 on step %d  %d\n", j1, k1 );
  printf ( "  X2 on step %d  %d\n", j2, k1 );
  printf ( "  X3 on step %d  %d\n", j1, k2 );
  printf ( "\n" );
  printf ( "  We use %d steps in the J direction\n", nstep1 );
  printf ( "  and %d steps in the K direction.\n", nstep2 );
  printf ( "\n" );
  printf ( "  The points X1, X2 and X3 are:\n" );
  printf ( "\n" );
  for ( i = 0; i < DIM_NUM; i++ )
  {
    printf ( "  %10f", x1[i] );
  }
  printf ( "\n" );
  printf ( "\n" );
  for ( i = 0; i < DIM_NUM; i++ )
  {
    printf ( "  %10f", x2[i] );
  }
  printf ( "\n" );
  printf ( "\n" );
  for ( i = 0; i < DIM_NUM; i++ )
  {
    printf ( "  %10f", x3[i] );
  }
  printf ( "\n" );

  for ( j = 1; j <= nstep1; j++ )
  {
    printf ( "\n" );
    for ( k = 1; k <= nstep2; k++ )
    {
      x = grid4n ( j, j1, j2, k, k1, k2, DIM_NUM, nstep1, nstep2, 
        x1, x2, x3 );
      printf ( "  %3d  %3d", j, k );
      for ( i = 1; i <= DIM_NUM; i++ )
      {
        printf ( "  %10f", x[i-1] );
      }
      printf ( "\n" );
      free ( x );
    }
  }
 
  return;
# undef DIM_NUM
}
/******************************************************************************/

void test19 ( )

/******************************************************************************/
/*
  Purpose:

    TEST19 tests INDEX1_COL, INDEX1_ROW, and related functions.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    03 April 2010

  Author:

    John Burkardt
*/
{
  int i;
  int i_max;
  int i_min;
  int in[4];
  int in_max[4];
  int in_min[4];
  int index_min;
  int j;
  int j_max;
  int j_min;
  int k;
  int k_max;
  int k_min;
  int l;
  int l_max;
  int l_min;
  int n;
  int value;

  printf ( "\n" );
  printf ( "TEST19\n" );
  printf ( "  INDEX1_COL column indexes a 1D array,\n" );
  printf ( "  INDEX1_ROW row indexes a 1D array,\n" );
  printf ( "  and there are several more versions of these functions.\n" );

  printf ( "\n" );
  printf ( "  By COLS:\n" );
  printf ( "\n" );
  printf ( "  Imin     I  Imax  Xmin Index\n" );
  printf ( "\n" );

  i_min = 1;
  i = 3;
  i_max = 5;
  index_min = 0;

  value = index1_col ( i_min, i, i_max, index_min );
  printf ( "\n" );
  printf ( "  %4d  %4d %4d\n", i_min, i, i_max );
  printf ( "        INDEX1_COL  %4d  %4d\n", index_min, value );

  n = 1;
  in_min[0] = 1;
  in[0] = 3;
  in_max[0] = 5;
  index_min = 0;
  value = indexn_col ( n, in_min, in, in_max, index_min );
  printf ( "        INDEXN_COL  %4d  %4d\n", index_min, value );

  i_min = 1;
  i = 3;
  i_max = 5;
  j_min = 1;
  j = 2;
  j_max = 4;
  index_min = 0;
  value = index2_col ( i_min, i, i_max, j_min, j, j_max, index_min );
  printf ( "\n" );
  printf ( "  %4d  %4d %4d\n", i_min, i, i_max );
  printf ( "  %4d  %4d %4d\n", j_min, j, j_max );
  printf ( "        INDEX2_COL  %4d  %4d\n", index_min, value );

  n = 2;
  in_min[0] = 1;
  in[0] = 3;
  in_max[0] = 5;
  in_min[1] = 1;
  in[1] = 2;
  in_max[1] = 4;
  index_min = 0;
  value = indexn_col ( n, in_min, in, in_max, index_min );
  printf ( "        INDEXN_COL  %4d  %4d\n", index_min, value );

  i_min = 1;
  i = 3;
  i_max = 5;
  j_min = 1;
  j = 2;
  j_max = 4;
  k_min = 1;
  k = 1;
  k_max = 3;
  index_min = 0;
  value = index3_col ( i_min, i, i_max, j_min, j, j_max, k_min, k, k_max, 
   index_min );
  printf ( "\n" );
  printf ( "  %4d  %4d %4d\n", i_min, i, i_max );
  printf ( "  %4d  %4d %4d\n", j_min, j, j_max );
  printf ( "  %4d  %4d %4d\n", k_min, k, k_max );
  printf ( "        INDEX3_COL  %4d  %4d\n", index_min, value );

  n = 3;
  in_min[0] = 1;
  in[0] = 3;
  in_max[0] = 5;
  in_min[1] = 1;
  in[1] = 2;
  in_max[1] = 4;
  in_min[2] = 1;
  in[2] = 1;
  in_max[2] = 3;
  index_min = 0;
  value = indexn_col ( n, in_min, in, in_max, index_min );
  printf ( "        INDEXN_COL  %4d  %4d\n", index_min, value );

  i_min = 1;
  i = 3;
  i_max = 5;
  j_min = 1;
  j = 2;
  j_max = 4;
  k_min = 1;
  k = 1;
  k_max = 3;
  l_min = 1;
  l = 2;
  l_max = 2;
  index_min = 0;
  value = index4_col ( i_min, i, i_max, j_min, j, j_max, k_min, k, k_max, 
    l_min, l, l_max, index_min );
  printf ( "\n" );
  printf ( "  %4d  %4d %4d\n", i_min, i, i_max );
  printf ( "  %4d  %4d %4d\n", j_min, j, j_max );
  printf ( "  %4d  %4d %4d\n", k_min, k, k_max );
  printf ( "  %4d  %4d %4d\n", l_min, l, l_max );
  printf ( "        INDEX4_COL  %4d  %4d\n", index_min, value );

  n = 4;
  in_min[0] = 1;
  in[0] = 3;
  in_max[0] = 5;
  in_min[1] = 1;
  in[1] = 2;
  in_max[1] = 4;
  in_min[2] = 1;
  in[2] = 1;
  in_max[2] = 3;
  in_min[3] = 1;
  in[3] = 2;
  in_max[3] = 2;
  index_min = 0;
  value = indexn_col ( n, in_min, in, in_max, index_min );
  printf ( "        INDEXN_COL  %4d  %4d\n", index_min, value );

  printf ( "\n" );
  printf ( "  By ROWS:\n" );
  printf ( "\n" );
  printf ( "  Imin     I  Imax  Xmin Index\n" );
  printf ( "\n" );

  i_min = 1;
  i = 3;
  i_max = 5;
  index_min = 0;
  value = index1_row ( i_min, i, i_max, index_min );
  printf ( "\n" );
  printf ( "  %4d  %4d %4d\n", i_min, i, i_max );
  printf ( "        INDEX1_ROW  %4d  %4d\n", index_min, value );

  n = 1;
  in_min[0] = 1;
  in[0] = 3;
  in_max[0] = 5;
  index_min = 0;
  value = indexn_row ( n, in_min, in, in_max, index_min );
  printf ( "        INDEXN_ROW  %4d  %4d\n", index_min, value );

  i_min = 1;
  i = 3;
  i_max = 5;
  j_min = 1;
  j = 2;
  j_max = 4;
  index_min = 0;
  value = index2_row ( i_min, i, i_max, j_min, j, j_max, index_min );
  printf ( "\n" );
  printf ( "  %4d  %4d %4d\n", i_min, i, i_max );
  printf ( "  %4d  %4d %4d\n", j_min, j, j_max );
  printf ( "        INDEX2_ROW  %4d  %4d\n", index_min, value );

  n = 2;
  in_min[0] = 1;
  in[0] = 3;
  in_max[0] = 5;
  in_min[1] = 1;
  in[1] = 2;
  in_max[1] = 4;
  index_min = 0;
  value = indexn_row ( n, in_min, in, in_max, index_min );
  printf ( "        INDEXN_ROW  %4d  %4d\n", index_min, value );

  i_min = 1;
  i = 3;
  i_max = 5;
  j_min = 1;
  j = 2;
  j_max = 4;
  k_min = 1;
  k = 1;
  k_max = 3;
  index_min = 0;
  value = index3_row ( i_min, i, i_max, j_min, j, j_max, k_min, k, k_max, 
    index_min );
  printf ( "\n" );
  printf ( "  %4d  %4d %4d\n", i_min, i, i_max );
  printf ( "  %4d  %4d %4d\n", j_min, j, j_max );
  printf ( "  %4d  %4d %4d\n", k_min, k, k_max );
  printf ( "        INDEX3_ROW  %4d  %4d\n", index_min, value );

  n = 3;
  in_min[0] = 1;
  in[0] = 3;
  in_max[0] = 5;
  in_min[1] = 1;
  in[1] = 2;
  in_max[1] = 4;
  in_min[2] = 1;
  in[2] = 1;
  in_max[2] = 3;
  index_min = 0;
  value = indexn_row ( n, in_min, in, in_max, index_min );
  printf ( "        INDEXN_ROW  %4d  %4d\n", index_min, value );

  i_min = 1;
  i = 3;
  i_max = 5;
  j_min = 1;
  j = 2;
  j_max = 4;
  k_min = 1;
  k = 1;
  k_max = 3;
  l_min = 1;
  l = 2;
  l_max = 2;
  index_min = 0;
  value = index4_row ( i_min, i, i_max, j_min, j, j_max, k_min, k, k_max, 
    l_min, l, l_max, index_min );
  printf ( "\n" );
  printf ( "  %4d  %4d %4d\n", i_min, i, i_max );
  printf ( "  %4d  %4d %4d\n", j_min, j, j_max );
  printf ( "  %4d  %4d %4d\n", k_min, k, k_max );
  printf ( "  %4d  %4d %4d\n", l_min, l, l_max );
  printf ( "        INDEX4_ROW  %4d  %4d\n", index_min, value );

  n = 4;
  in_min[0] = 1;
  in[0] = 3;
  in_max[0] = 5;
  in_min[1] = 1;
  in[1] = 2;
  in_max[1] = 4;
  in_min[2] = 1;
  in[2] = 1;
  in_max[2] = 3;
  in_min[3] = 1;
  in[3] = 2;
  in_max[3] = 2;
  index_min = 0;
  value = indexn_row ( n, in_min, in, in_max, index_min );
  printf ( "        INDEXN_ROW  %4d  %4d\n", index_min, value );

  return;
}
/******************************************************************************/

void test20 ( )

/******************************************************************************/
/*
  Purpose:

    TEST20 tests ISBN_CHECK.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    25 July 2010

  Author:

    John Burkardt
*/
{
# define TEST_NUM 8

  int check;
  char isbn[14];
  int test;

  printf ( "\n" );
  printf ( "TEST20\n" );
  printf ( "  ISBN_CHECK checks ISBN's.\n" );
  printf ( "\n" );
  printf ( "  A correct ISBN has a checksum of 0.\n" );
  printf ( "\n" );
//
//  Sorry, but until I figure out a decent way to set up and
//  access an array of strings we'll have to do this the
//  bonehead way.
//
  for ( test = 0; test < TEST_NUM; test++ )
  {
    if ( test == 0 )
    {
      strcpy ( isbn, "0-8493-9640-9" );
    }
    else if ( test == 1 ) 
    {
      strcpy ( isbn, "0-201-54275-7" );
    }
    else if ( test == 2 ) 
    {
      strcpy ( isbn, "0-521-35796-9" );
    }
    else if ( test == 3 ) 
    {
      strcpy ( isbn, "0-07-034025-0" );
    }
    else if ( test == 4 ) 
    {
      strcpy ( isbn, "0-7493-9640-9" );
    }
    else if ( test == 5 ) 
    {
      strcpy ( isbn, "0-201-54275-X" );
    }
    else if ( test == 6 ) 
    {
      strcpy ( isbn, "0-521-X5796-9" );
    }
    else if ( test == 7 ) 
    {
      strcpy ( isbn, "0-37-034025-0" );
    }
    check = isbn_check ( isbn );

    printf ( "  %s  %d\n", isbn, check );
  }

  return;
# undef TEST_NUM
}
/******************************************************************************/

void test21 ( )

/******************************************************************************/
/*
  Purpose:

    TEST21 tests ISBN_FILL.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    24 July 2010

  Author:

    John Burkardt
*/
{
# define TEST_NUM 5

  int check;
  char isbn[20];
  char isbn_test[20];    
  int test;

  printf ( "\n" );
  printf ( "TEST21\n" );
  printf ( "  ISBN_FILL can fill in a single missing digit\n" );
  printf ( "  in an ISBN.\n" );
  printf ( "\n" );

  for ( test = 0; test < TEST_NUM; test++ )
  {
    if ( test == 0 )
    {
      strcpy ( isbn_test, "0-?493-9640-9" );
    }
    else if ( test == 1 )
    {
      strcpy ( isbn_test, "0-201-5427?-7" );
    }
    else if ( test == 2 )
    {
      strcpy ( isbn_test, "0-521-35796-?" );
    }
    else if ( test == 3 )
    {
      strcpy ( isbn_test, "?-07-034025-0" );
    }
    else if ( test == 4 )
    {
      strcpy ( isbn_test, "0-07-05?489-2" );
    }
    strcpy ( isbn, isbn_test );

    isbn_fill ( isbn );

    check = isbn_check ( isbn );

    printf ( "  %s  %s  %d\n", isbn_test, isbn, check );
  }

  return;
# undef TEST_NUM
}
/******************************************************************************/

void test22 ( )

/******************************************************************************/
/*
  Purpose:

    TEST22 tests LCM_12N.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    24 July 2010

  Author:

    John Burkardt
*/
{
  int n;

  printf ( "\n" );
  printf ( "TEST22\n" );
  printf ( "  LCM_12N computes the least common multiple of the\n" );
  printf ( "  integers 1 through N.\n" );
  printf ( "\n" );
  printf ( "  N     LCM_12N ( N )\n" );
  printf ( "\n" );
  for ( n = 1; n <= 12; n++ )
  {
    printf ( "  %3d  %8d\n", n, lcm_12n ( n ) );
  }

  return;
}
/******************************************************************************/

void test225 ( )

/******************************************************************************/
/*
  Purpose:

    TEST225 tests L4MAT_PRINT.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    03 November 2011

  Author:

    John Burkardt
*/
{
  int *a;
  int i;
  int j;
  int m = 20;
  int n = 50;

  a = ( int * ) malloc ( m * n * sizeof ( int ) );

  printf ( "\n" );
  printf ( "TEST225\n" );
  printf ( "  L4MAT_PRINT prints a logical matrix.\n" );

  for ( i = 0; i < m; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      a[i+j*m] = ( ( i + 1 ) % ( j + 1 ) == 0 );
    }
  }

  l4mat_print ( m, n, a, "  A(I,J) = I+1 is divisible by J+1" );

  free ( a );

  return;
}
/******************************************************************************/

void test23 ( )

/******************************************************************************/
/*
  Purpose:

    TEST23 tests LUHN_CHECK.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    24 July 2010

  Author:

    John Burkardt
*/
{
# define TEST_NUM 4

  int check_sum;
  int check_sum_test[TEST_NUM] = {
     6, 
    20, 
    40, 
    80 };
  int *digit;
  int digit_num;
  int digit_num_test[TEST_NUM] = {
     4, 
     4, 
     9, 
    15 };
  int digit_test[32] = {
    1, 1, 1, 1, 
    8, 7, 6, 3, 
    4, 4, 6, 6, 6, 7, 6, 5, 1, 
    3, 7, 7, 9, 5, 6, 5, 7, 0, 9, 4, 4, 7, 2, 6 };
  int i;
  int ilo;
  int test;

  printf ( "\n" );
  printf ( "TEST23\n" );
  printf ( "  LUHN_CHECK computes the Luhn checksum\n" );
  printf ( "  for a string of digits.\n" );
  printf ( "\n" );
  printf ( "  A correct string has a checksum divisible by 10.\n" );

  ilo = 0;
  for ( test = 0; test < TEST_NUM; test++ )
  {
    digit_num = digit_num_test[test];

    digit = ( int * ) malloc ( digit_num * sizeof ( int ) );
    for ( i = 0; i < digit_num; i++ )
    {
      digit[i] = digit_test[i+ilo];
    }
    ilo = ilo + digit_num;

    check_sum = luhn_check ( digit_num, digit );

    printf ( "\n" );
    printf ( "  Test number %d\n", test );
    printf ( "  Number of digits = %d\n", digit_num );
    printf ( "  Digits = " );
    for ( i = 0; i < digit_num; i++ )
    {
      printf ( "%d", digit[i] );
    }
    printf ( "\n" );
    printf ( "  Computed check sum = %d\n", check_sum );
    printf ( "  Correct check sum =  %d\n", check_sum_test[test] );

    free ( digit );
  }

  return;
# undef TEST_NUM
}
/******************************************************************************/

void test24 ( )

/******************************************************************************/
/*
  Purpose:

    TEST24 tests PERM_INVERSE;

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    24 July 2010

  Author:

    John Burkardt
*/
{
# define N 7

  int p[N] = { 4, 3, 5, 1, 7, 6, 2 };

  printf ( "\n" );
  printf ( "TEST24\n" );
  printf ( "  PERM_INVERSE inverts a permutation in place;\n" );
  printf ( "\n" );

  perm_print ( N, p, "  The original permutation:" );

  perm_inverse ( N, p );

  perm_print ( N, p, "  The inverted permutation:" );

  return;
# undef N
}
/******************************************************************************/

void test25 ( )

/******************************************************************************/
/*
  Purpose:

    TEST25 tests PRIME_GE.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    24 July 2010

  Author:

    John Burkardt
*/
{
  int n;
  int p;

  printf ( "\n" );
  printf ( "TEST25\n" );
  printf ( "  PRIME_GE returns the smallest prime number greater\n" );
  printf ( "  than or equal to N.\n" );
  printf ( "\n" );
  printf ( "  N   PRIME_GE\n" );
  printf ( "\n" );
  for ( n = 1; n <= 10; n++ )
  {
    p = prime_ge ( n );
    printf ( "  %6d  %6d\n", n, p );
  }

  return;
}
/******************************************************************************/

void test26 ( )

/******************************************************************************/
/*
  Purpose:

    TEST26 tests RANDOM_INITIALIZE.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    25 July 2010

  Author:

    John Burkardt
*/
{
  int i;
  double r1;
  double r2;
  double r3;
  unsigned long seed;
  unsigned long seed_in;
  unsigned long seed_out;

  printf ( "\n" );
  printf ( "TEST26\n" );
  printf ( "  RANDOM_INITIALIZE can make up a seed for the C\n" );
  printf ( "  random number generator RANDOM, or use a\n" );
  printf ( "  single SEED value from the user.\n" );
  printf ( "\n" );
  printf ( "  Calling RANDOM_INITIALIZE with a zero input value of SEED\n" );
  printf ( "  tells the routine to make up a seed.  And, at least for\n" );
  printf ( "  calls a few milliseconds apart, the output SEED should\n" );
  printf ( "  be different.\n" );
  printf ( "\n" );
  printf ( "  In any case, if RANDOM is restarted by calling\n" );
  printf ( "  RANDOM_INITIALIZE with a nonzero input SEED, then\n" );
  printf ( "  the random number sequence should repeat.\n" );
  printf ( "\n" );
  printf ( "  Call RANDOM_INITIALIZE 10 times, with a zero input SEED.\n" );
  printf ( "  Also, get the first three real random values.\n" );
  printf ( "\n" );
  printf ( "    SEED_IN         SEED_OUT     Random 1, 2, 3\n" );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    seed_in = 0;
    seed = seed_in;
    seed_out = random_initialize ( seed );

    r1 = ( double ) random ( ) / ( double ) ( RAND_MAX );
    r2 = ( double ) random ( ) / ( double ) ( RAND_MAX );
    r3 = ( double ) random ( ) / ( double ) ( RAND_MAX );

    printf ( "  %12d  %12d  %12f  %12f  %12f\n", seed_in,seed_out, r1, r2, r3 );
  }

  printf ( "\n" );
  printf ( "  Now call RANDOM_INITIALIZE with SEED = 5, 95, 5, 95.\n" );
  printf ( "  We promise the random numbers will repeat the second time.\n" );
  printf ( "\n" );
  printf ( "    SEED_IN         SEED_OUT     Random 1, 2, 3\n" );
  printf ( "\n" );

  seed_in = 5;

  for ( i = 1; i <= 4; i++ )
  {
    seed = seed_in;
    seed_out = random_initialize ( seed );

    r1 = ( double ) random ( ) / ( double ) ( RAND_MAX );
    r2 = ( double ) random ( ) / ( double ) ( RAND_MAX );
    r3 = ( double ) random ( ) / ( double ) ( RAND_MAX );

    printf ( "  %12d  %12d  %12f  %12f  %12f\n", seed_in,seed_out, r1, r2, r3 );

    seed_in = 100 - seed_in;
  }

  return;
}
/******************************************************************************/

void test27 ( )

/******************************************************************************/
/*
  Purpose:

    TEST27 tests RANDOM.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    24 July 2010

  Author:

    John Burkardt
*/
{
# define N_MAX 1000

  int i;
  int n;
  int seed = 123456789;
  double x[N_MAX];
  double x_max;
  double x_mean;
  double x_min;
  double x_var;

  printf ( "\n" );
  printf ( "TEST27\n" );
  printf ( "  RANDOM is an intrinsic C routine\n" );
  printf ( "  to computer uniform random numbers.\n" );
  printf ( "  Using initial random number seed = %d\n", seed );
//
//  Set the random number seed.
//
  srandom ( seed );
//
//  Test 1:
//  Simply call 5 times for 1 value, and print.
//
  printf ( "\n" );
  printf ( "  Test #1: Call 5 times, 1 value each time.\n" );
  printf ( "\n" );

  n = 1;
  for ( i = 1; i <= 5; i++ )
  {
    x[0] = ( double ) random ( ) / ( double ) ( RAND_MAX );
    printf ( "  %6d  %12f\n", i, x[0] );
  }
//
//  Test 2:
//  Restore the random number seed, and repeat.
//
  printf ( "\n" );
  printf ( "  Test #2: Restore the random number seed.\n" );
  printf ( "  Call 5 times, 1 value each time.\n" );
  printf ( "  The results should be identical.\n" );
  printf ( "\n" );

  seed = 123456789;
  srandom ( seed );

  n = 1;
  for ( i = 1; i <= 5; i++ )
  {
    x[0] = ( double ) random ( ) / ( double ) ( RAND_MAX );
    printf ( "  %6d  %12f\n", i, x[0] );
  }
//
//  Test 5:
//  Determine the minimum, maximum, mean and variance.
//
  n = N_MAX;
  for ( i = 0; i < n; i++ )
  {
    x[i] = ( double ) random ( ) / ( double ) ( RAND_MAX );
  }
  x_min = r8vec_min ( n, x );
  x_max = r8vec_max ( n, x );
  x_mean = r8vec_mean ( n, x );
  x_var = r8vec_variance ( n, x );

  printf ( "\n" );
  printf ( "  Test #5:\n" );
  printf ( "  Number of samples was %d\n", n );
  printf ( "  Minimum value was %f\n", x_min );
  printf ( "  Maximum value was %f\n", x_max );
  printf ( "  Average value was %f\n", x_mean );
  printf ( "  Variance was      %f\n", x_var );
  printf ( "  Expected average  %f\n", 0.5 );
  printf ( "  Expected variance %f\n", 1.0 / 12.0 );

  return;
# undef N_MAX
}
/******************************************************************************/

void test29 ( )

/******************************************************************************/
/*
  Purpose:

    TEST29 tests RAT_FACTOR.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    24 July 2010

  Author:

    John Burkardt
*/
{
# define FACTOR_MAX 10

  int factor[FACTOR_MAX];
  int factor_num;
  int i;
  int m;
  int mleft;
  int n;
  int nleft;
  int power[FACTOR_MAX];

  printf ( "\n" );
  printf ( "TEST29\n" );
  printf ( "  RAT_FACTOR factors a rational value.\n" );

  m = 13 * 7 * 9 * 2;
  n = 12;

  printf ( "\n" );
  printf ( "  Rational value is %d / %d\n", m, n );

  rat_factor ( m, n, FACTOR_MAX, &factor_num, factor, power, &mleft, &nleft );

  printf ( "\n" );
  printf ( "  Prime representation:\n" );
  printf ( "\n" );
  printf ( "  I, FACTOR(I), POWER(I)\n" );
  printf ( "\n" );

  if ( mleft != 1 || nleft != 1 )
  {
    printf ( "       0  %6d  %6d (UNFACTORED PORTION)\n", mleft, nleft );
    printf ( "\n" );
  }

  for ( i = 0; i < factor_num; i++ )
  {
    printf ( "  %6d  %6d  %6d\n", i + 1, factor[i], power[i] );
  }
 
  return;
# undef FACTOR_MAX
}
/******************************************************************************/

void test30 ( )

/******************************************************************************/
/*
  Purpose:

    TEST30 tests ROOTS_TO_R8POLY.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    23 July 2010

  Author:

    John Burkardt
*/
{
# define N 4

  double *c;
  double *x;

  printf ( "\n" );
  printf ( "TEST30\n" );
  printf ( "  ROOTS_TO_R8POLY computes the coefficients of\n" );
  printf ( "  a polynomial from its roots.\n" );
  printf ( "  R8POLY_PRINT prints a polynomial.\n" );

  x = r8vec_indicator_new ( N );

  r8vec_print ( N, x, "  Roots:" );

  c = roots_to_r8poly ( N, x );

  r8poly_print ( N, c, "  The polynomial" );

  free ( c );
  free ( x );

  return;
# undef N
}
/******************************************************************************/

void test31 ( )

/******************************************************************************/
/*
  Purpose:

    TEST31 tests SORT_HEAP_EXTERNAL.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    23 July 2010

  Author:

    John Burkardt
*/
{
# define N 20

  int a[N];
  int b;
  int c;
  int i;
  int indx;
  int isgn;
  int itemp;
  int j;
  int seed;

  printf ( "\n" );
  printf ( "TEST31\n" );
  printf ( "  SORT_HEAP_EXTERNAL sorts objects externally.\n" );

  indx = 0;
  i = 0;
  j = 0;
  isgn = 0;

  b = 1;
  c = N;
  seed = 123456789;

  for ( i = 0; i < N; i++ )
  {
    a[i] = i4_uniform ( b, c, &seed );
  }
 
  i4vec_print ( N, a, "  Unsorted array:" );
/*
  Call the sort routine over and over.
*/
  for ( ;; )
  {
    sort_heap_external ( N, &indx, &i, &j, isgn );
/*
  If the return value of INDX is negative, we're asked to compare
  array elements I and J;
*/
    if ( indx < 0 )
    {
      if ( a[i] <= a[j] )
      {
        isgn = -1;
      }
      else
      {
        isgn = 1;
      }

    }
/*
  ...and if the return value of INDX is positive, we're asked to switch
  array elements I and J;
*/
    else if ( 0 < indx )
    {
      i4_swap ( &a[i], &a[j] );
/*
  ...and if the return value of INDX is 0, we're done.
*/
    } 
    else
    {
      break;
    }

  }

  i4vec_print ( N, a, "  Sorted array:" );
 
  return;

# undef N
}
/******************************************************************************/

void test32 ( )

/******************************************************************************/
/*
  Purpose:

    TEST32 tests TVEC_EVEN, TVEC_EVEN2 and TVEC_EVEN3.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    23 July 2010

  Author:

    John Burkardt
*/
{
  int nt;
  double *t;

  printf ( "\n" );
  printf ( "TEST32\n" );
  printf ( "  For evenly spaced angles between 0 and 2*PI:\n" );
  printf ( "  TVEC_EVEN\n" );
  printf ( "  TVEC_EVEN2\n" );
  printf ( "  TVEC_EVEN3\n" );

  nt = 4;

  t = tvec_even ( nt );

  r8vec_print ( nt, t, "  TVEC_EVEN" );

  free ( t );

  nt = 4;

  t = tvec_even2 ( nt );

  r8vec_print ( nt, t, "  TVEC_EVEN2" );

  free ( t );

  nt = 4;

  t = tvec_even3 ( nt );

  r8vec_print ( nt, t, "  TVEC_EVEN3" );

  free ( t );

  return;
}
/******************************************************************************/

void test33 ( )

/******************************************************************************/
/*
  Purpose:

    TEST33 tests TVEC_EVEN_BRACKET, TVEC_EVEN_BRACKET2, TVEC_EVEN_BRACKET3.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    23 July 2010

  Author:

    John Burkardt
*/
{
  int i;
  int nt;
  double *t;
  double theta1;
  double theta2;

  printf ( "\n" );
  printf ( "TEST33\n" );
  printf ( "  For evenly spaced angles between THETA1 and THETA2:\n" );
  printf ( "  TVEC_EVEN_BRACKET\n" );
  printf ( "  TVEC_EVEN_BRACKET2.\n" );
  printf ( "  TVEC_EVEN_BRACKET3.\n" );

  nt = 4;
  theta1 = 30.0;
  theta2 = 90.0;

  t = tvec_even_bracket ( nt, theta1, theta2 );

  printf ( "\n" );
  printf ( "    NT = %d\n", nt );
  printf ( "    THETA1 = %f\n", theta1 );
  printf ( "    THETA2 = %f\n", theta2 );

  r8vec_print ( nt, t, "  TVEC_BRACKET" );

  nt = 5;

  t = tvec_even_bracket2 ( nt, theta1, theta2 );

  printf ( "\n" );
  printf ( "    NT = %d\n", nt );
  printf ( "    THETA1 = %f\n", theta1 );
  printf ( "    THETA2 = %f\n", theta2 );

  r8vec_print ( nt, t, "  TVEC_EVEN_BRACKET2" );

  nt = 3;

  t = tvec_even_bracket3 ( nt, theta1, theta2 );

  printf ( "\n" );
  printf ( "    NT = %d\n", nt );
  printf ( "    THETA1 = %f\n", theta1 );
  printf ( "    THETA2 = %f\n", theta2 );

  r8vec_print ( nt, t, "  TVEC_EVEN_BRACKET3" );

  return;
}
/******************************************************************************/

void test34 ( )

/******************************************************************************/
/*
  Purpose:

    TEST34 tests UPC_CHECK_DIGIT.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    21 July 2010

  Author:

    John Burkardt
*/
{
# define TEST_NUM 2

  int c;
  int l;
  int l_test[TEST_NUM] = { 72890, 12345 };
  int p;
  int p_test[TEST_NUM] = { 0, 0 };
  int r;
  int r_test[TEST_NUM] = { 11, 67890 };
  int test;

  printf ( "\n" );
  printf ( "TEST34\n" );
  printf ( "  UPC_CHECK_DIGIT determines the check digit for a UPC.\n" );
  printf ( "\n" );
  printf ( "  P-LLLLL-RRRRR-C\n" );
  printf ( "\n" );

  for ( test = 0; test < TEST_NUM; test++ )
  {
    p = p_test[test];
    l = l_test[test];
    r = r_test[test];

    c = upc_check_digit ( p, l, r );

    printf ( "  %1d-%05d-%05d-%1d\n", p, l, r, c );
  }
  return;
# undef TEST_NUM
}
/******************************************************************************/

void test35 ( )

/******************************************************************************/
/*
  Purpose:

    TEST35 calls VERSINE_PULSE.

  Modified:

    26 February 2010

  Author:

    John Burkardt

  Parameters:

    Input, double T, the current time.

    Input, double TA, the time at which the pulse begins.

    Input, double TB, the time at which the pulse finishes
    its first period.

    Input, double V1, the constant value.

    Input, double AMP, the amplitude of the pulse.  Actually,
    the pulse will vary from 0 to 2 * AMP.
*/
{
  double amp;
  int i;
  double t;
  double ta;
  double tb;
  double v;
  double v1;

  printf ( "\n" );
  printf ( "TEST35\n" );
  printf ( "  VERSINE_PULSE adds a versine pulse to a constant signal.\n" );
  printf ( "\n" );

  ta = 2.0;
  tb = 4.0;
  v1 = 1.0;
  amp = 3.0;

  for ( i = 0; i <= 100; i++ )
  {
    t = ( double ) i / 10.0;
    v = versine_pulse ( t, ta, tb, v1, amp );
    printf ( "  %d  %f  %f\n", i, t, v );
  }
  return;
}
