# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <complex.h>

# include "r4lib.h"

int main ( );

void test001 ( );
void test002 ( );
void test003 ( );
void test004 ( );
void test005 ( );
void test006 ( );
void test007 ( );
void test009 ( );
void test010 ( );
void test023 ( );
void test0235 ( );
void test026 ( );
void test027 ( );
void test028 ( );
void test12555 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for R4LIB_PRB.

  Discussion:

    R4LIB_PRB tests the R4LIB library.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    28 September 2014

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "R4LIB_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the R4LIB library.\n" );

  test001 ( );
  test002 ( );
  test003 ( );
  test004 ( );
  test005 ( );
  test006 ( );
  test007 ( );
  test009 ( );

  test010 ( );

  test023 ( );
  test0235 ( );
  test026 ( );
  test027 ( );
  test028 ( );

  test12555 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "R4LIB_PRB\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void test001 ( )

/******************************************************************************/
/*
  Purpose:

    TEST001 tests R4_ABS.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    26 June 2008

  Author:

    John Burkardt
*/
{
  float r4;
  float r4_absolute;
  float r4_hi = 5.0;
  float r4_lo = -3.0;
  int seed;
  int test;
  int test_num = 10;

  seed = 123456789;

  printf ( "\n" );
  printf ( "TEST001\n" );
  printf ( "  R4_ABS returns the absolute value of an R4.\n" );
  printf ( "\n" );

  for ( test = 1; test <= test_num; test++ )
  {
    r4 = r4_uniform_ab ( r4_lo, r4_hi, &seed );
    r4_absolute = r4_abs ( r4 );
    printf ( "  %10.6f  %10.6f\n", r4, r4_absolute );
  }

  return;
}
/******************************************************************************/

void test002 ( )

/******************************************************************************/
/*
  Purpose:

    TEST002 tests R4_ATAN.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    13 September 2006

  Author:

    John Burkardt
*/
{
  float x;
  float y;

  fprintf ( stdout, "\n" );
  fprintf ( stdout, "TEST002\n" );
  fprintf ( stdout, "  R4_ATAN computes the arc-tangent given Y and X;\n" );
  fprintf ( stdout, "  ATAN2 is the system version of this routine.\n" );
  fprintf ( stdout, "\n" );
  fprintf ( stdout, "           X           Y  ATAN2(Y,X)  ATAN4(Y,X)\n" );
  fprintf ( stdout, "\n" );

  x = 1.0;
  y = 0.0;
  fprintf ( stdout, "  %10f  %10f  %10f  %10f\n", 
    x, y, atan2 ( y, x ), r4_atan ( y, x ) );

  x = 1.0;
  y = 1.0;
  fprintf ( stdout, "  %10f  %10f  %10f  %10f\n", 
    x, y, atan2 ( y, x ), r4_atan ( y, x ) );

  x = 0.0;
  y = 1.0;
  fprintf ( stdout, "  %10f  %10f  %10f  %10f\n", 
    x, y, atan2 ( y, x ), r4_atan ( y, x ) );

  x = -1.0;
  y = 1.0;
  fprintf ( stdout, "  %10f  %10f  %10f  %10f\n", 
    x, y, atan2 ( y, x ), r4_atan ( y, x ) );

  x = -1.0;
  y = 0.0;
  fprintf ( stdout, "  %10f  %10f  %10f  %10f\n", 
    x, y, atan2 ( y, x ), r4_atan ( y, x ) );

  x = - 1.0;
  y = - 1.0;
  fprintf ( stdout, "  %10f  %10f  %10f  %10f\n", 
    x, y, atan2 ( y, x ), r4_atan ( y, x ) );

  x =   0.0;
  y = - 1.0;
  fprintf ( stdout, "  %10f  %10f  %10f  %10f\n", 
    x, y, atan2 ( y, x ), r4_atan ( y, x ) );

  x =   1.0;
  y = - 1.0;
  fprintf ( stdout, "  %10f  %10f  %10f  %10f\n", 
    x, y, atan2 ( y, x ), r4_atan ( y, x ) );

  return;
}
/******************************************************************************/

void test003 ( )

/******************************************************************************/
/*
  Purpose:

    TEST003 tests R4_CAS.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 August 2010

  Author:

    John Burkardt
*/
{
# define TEST_NUM 12

  int test;
  float x;

  printf ( "\n" );
  printf ( "TEST003\n" );
  printf ( "  R4_CAS evaluates the casine of a number.\n" );
  printf ( "\n" );
  printf ( "        X           R4_CAS ( X )\n" );
  printf ( "\n" );

  for ( test = 0; test <= TEST_NUM; test++ )
  {
    x = r4_pi ( ) * ( float ) ( test ) / ( float ) ( TEST_NUM );
    printf ( "  %14f  %14f\n", x, r4_cas ( x ) );
  }

  return;
# undef TEST_NUM
}
/******************************************************************************/

void test004 ( )

/******************************************************************************/
/*
  Purpose:

    TEST004 tests R4_CEILING.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 August 2010

  Author:

    John Burkardt
*/
{
  int i;
  float rval;
  float rval_rounded;

  printf ( "\n" );
  printf ( "TEST004\n" );
  printf ( "  R4_CEILING rounds a value up.\n" );
  printf ( "\n" );

  for ( i = -6; i <= 6; i++ )
  {
    rval = ( float ) ( i ) / 5.0;
    rval_rounded = r4_ceiling ( rval );
    printf ( "  %14f  %14f\n", rval, rval_rounded );
  }

  return;
}
/******************************************************************************/

void test005 ( )

/******************************************************************************/
/*
  Purpose:

    TEST005 tests R4_DIFF.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 August 2010

  Author:

    John Burkardt
*/
{
# define TEST_NUM 15

  int ndig = 3;
  int test;
  float x = 1.0;
  float y;
  float y_test[TEST_NUM] = {
    0.0625, 0.125, 0.25, 0.50,  0.874, 
    0.876,  0.90,  0.95, 0.99,  1.0, 
    1.01,   1.05,  1.10, 3.0,  10.0 };

  printf ( "\n" );
  printf ( "TEST005\n" );
  printf ( "  R4_DIFF computes a difference X-Y to a given\n" );
  printf ( "    number of binary places.\n" );
  printf ( "\n" );
  printf ( "  For this test, we use %d binary places.\n", ndig );
  printf ( "\n" );
  printf ( "       X       Y       X-Y     R4_DIFF(X,Y)\n" );
  printf ( "\n" );

  for ( test = 0; test < TEST_NUM; test++ )
  {
    y = y_test[test];
    printf ( "  %10f  %10f  %10f  %10f\n", x, y, x - y, r4_diff ( x, y, ndig ) );
  }
 
  return;
# undef TEST_NUM
}
/******************************************************************************/

void test006 ( )

/******************************************************************************/
/*
  Purpose:

    TEST006 tests R4_DIGIT.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 August 2010

  Author:

    John Burkardt
*/
{
# define MAXDIG 20

  int idigit;
  float x;

  x = r4_pi ( );

  printf ( "\n" );
  printf ( "TEST006\n" );
  printf ( "  R4_DIGIT extracts decimal digits.\n" );
  printf ( "\n" );
  printf ( "  Here, we get digits of %24.16f\n", x );
  printf ( "\n" );
 
  printf ( "  " );
  for ( idigit = -2; idigit <= MAXDIG; idigit++ )
  {
    printf ( " %2d", idigit );
  }
  printf ( "\n" );

  printf ( "  " );
  for ( idigit = -2; idigit <= MAXDIG; idigit++ )
  {
    printf ( "  %1d", r4_digit ( x, idigit ) );
  }
  printf ( "\n" );
 
  return;
# undef MAXDIG
}
/******************************************************************************/

void test007 ( )

/******************************************************************************/
/*
  Purpose:

    TEST007 tests R4_EPSILON

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    25 June 2010

  Author:

    John Burkardt
*/
{
  float r;
  float r2;
  float s;
  float t;

  printf ( "\n" );
  printf ( "TEST007\n" );
  printf ( "  R4_EPSILON produces the R4 roundoff unit.\n" );
  printf ( "\n" );

  r = r4_epsilon ( );
  printf ( "  R = R4_EPSILON()  = %14e\n", r );

  s = 1.0 + r;
  t = s - 1.0;
  printf ( "  ( 1 + R ) - 1     = %14e\n", t );

  r2 = r / 2;
  s = 1.0 + r2;
  t = s - 1.0;
  printf ( "  ( 1 + (R/2) ) - 1 = %14e\n", t );

  return;
}
/******************************************************************************/

void test009 ( )

/******************************************************************************/
/*
  Purpose:

    TEST009 tests R4_HUGE.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    25 June 2010

  Author:

    John Burkardt
*/
{
  printf ( "\n" );
  printf ( "TEST009\n" );
  printf ( "  R4_HUGE returns a large R4 value;\n" );
  printf ( "\n" );
  printf ( "  R4_HUGE = %g\n", r4_huge ( ) );

  return;
}
/******************************************************************************/

void test010 ( )

/******************************************************************************/
/*
  Purpose:

    TEST010 tests R4_LOG_2.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    28 September 2014

  Author:

    John Burkardt
*/
{
  int test;
  const int test_num = 18;
  float x;
  float x_test[18] = {
    0.0E+00,  1.0E+00,  2.0E+00,   3.0E+00,  9.0E+00, 
   10.0E+00, 11.0E+00, 99.0E+00, 101.0E+00, -1.0E+00, 
   -2.0E+00, -3.0E+00, -9.0E+00,   0.5E+00,  0.33E+00, 
    0.25E+00, 0.20E+00, 0.01E+00 };

  printf ( "\n" );
  printf ( "TEST010\n" );
  printf ( "  R4_LOG_2 computes the logarithm base 2.\n" );
  printf ( "\n" );
  printf ( "  X, R4_LOG_2\n" );
  printf ( "\n" );

  for ( test = 0; test < test_num; test++ )
  {
    x = x_test[test];
    printf ( "  %14.6g  %14.6g\n", x, r4_log_2 ( x ) );
  }

  return;
}
/******************************************************************************/

void test023 ( )

/******************************************************************************/
/*
  Purpose:

    TEST023 tests R4_SIGN and R4_SIGN3.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    28 September 2014

  Author:

    John Burkardt
*/
{
  float r4;
  float r4_test[5] = { -1.25E+00, -0.25E+00, 0.0E+00, +0.5E+00, +9.0E+00 };
  float s1;
  float s2;
  int test;
  const int test_num = 5;

  printf ( "\n" );
  printf ( "TEST023\n" );
  printf ( "  R4_SIGN returns the sign of an R4.\n" );
  printf ( "  R4_SIGN3 returns the three-way sign of an R4.\n" );
  printf ( "\n" );
  printf ( "      R4    R4_SIGN(R4)  R4_SIGN3(R4)\n" );
  printf ( "\n" );

  for ( test = 0; test < test_num; test++ )
  {
    r4 = r4_test[test];
    s1 = r4_sign ( r4 );
    s2 = r4_sign3 ( r4 );
    printf ( "  %8.4f  %8.0f  %8.0f\n", r4, s1, s2 );
  }

  return;
}
/******************************************************************************/

void test0235 ( )

/******************************************************************************/
/*
  Purpose:

    TEST0235 tests R4_SWAP.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    25 June 2010

  Author:

    John Burkardt
*/
{
  float x;
  float y;

  fprintf ( stdout, "\n" );
  fprintf ( stdout, "TEST0235\n" );
  fprintf ( stdout, "  R4_SWAP swaps two reals.\n" );

  x = 1.0;
  y = 3.14159;

  fprintf ( stdout, "\n" );
  fprintf ( stdout, "  Before swapping: \n" );
  fprintf ( stdout, "\n" );
  fprintf ( stdout, "    X = %f\n", x );
  fprintf ( stdout, "    Y = %f\n", y );

  r4_swap ( &x, &y );

  fprintf ( stdout, "\n" );
  fprintf ( stdout, "  After swapping: \n" );
  fprintf ( stdout, "\n" );
  fprintf ( stdout, "    X = %f\n", x );
  fprintf ( stdout, "    Y = %f\n", y );

  return;
}
/******************************************************************************/

void test026 ( )

/******************************************************************************/
/*
  Purpose:

    TEST026 tests R4_UNIFORM_AB.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    05 May 2006

  Author:

    John Burkardt
*/
{
  float a;
  float b;
  float c;
  int i;
  int seed;

  b = 10.0;
  c = 25.0;
  seed = 17;

  printf ( "\n" );
  printf ( "TEST026\n" );
  printf ( "  R4_UNIFORM_AB produces a random real in a given range.\n" );
  printf ( "\n" );
  printf ( "  Using range %f <= A <= %f.\n", b, c );
  printf ( "\n" );

  printf ( "\n" );
  printf ( "  I   A\n" );
  printf ( "\n" );
  for ( i = 0; i < 10; i++ )
  {
    a = r4_uniform_ab ( b, c, &seed );
    printf ( "%6d  %10f\n", i, a );
  }

  return;
}
/******************************************************************************/

void test027 ( )

/******************************************************************************/
/*
  Purpose:

    TEST027 tests R4_UNIFORM_01.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    05 May 2006

  Author:

    John Burkardt
*/
{
  int i;
  int seed;
  float x;

  printf ( "\n" );
  printf ( "TEST027\n" );
  printf ( "  R4_UNIFORM_01 produces a sequence of random values.\n" );

  seed = 123456789;

  printf ( "\n" );
  printf ( "  Using random seed %d.\n", seed );

  printf ( "\n" );
  printf ( "  SEED   R4_UNIFORM_01(SEED)\n" );
  printf ( "\n" );
  for ( i = 0; i < 10; i++ )
  {
    printf ( "%12d", seed );
    x = r4_uniform_01 ( &seed );
    printf ( "  %f\n", x );
  }

  printf ( "\n" );
  printf ( "  Verify that the sequence can be restarted.\n" );
  printf ( "  Set the seed back to its original value, and see that\n" );
  printf ( "  we generate the same sequence.\n" );

  seed = 123456789;
  printf ( "\n" );
  printf ( "  SEED   R4_UNIFORM_01(SEED)\n" );
  printf ( "\n" );

  for ( i = 0; i < 10; i++ )
  {
    printf ( "%12d", seed );
    x = r4_uniform_01 ( &seed );
    printf ( "  %f\n", x );
  }

  return;
}
/******************************************************************************/

void test028 ( )

/******************************************************************************/
/*
  Purpose:

    TEST028 tests R4_UNIFORM_01.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    05 May 2006

  Author:

    John Burkardt
*/
{
# define N 1000

  int i;
  float max;
  float mean;
  float min;
  int n;
  int seed = 123456789;
  float x[N];
  float variance;

  printf ( "\n" );
  printf ( "TEST028\n" );
  printf ( "  R4_UNIFORM_01 samples a uniform random distribution in [0,1].\n" );
  printf ( "  distributed random numbers.\n" );
  printf ( "  Using initial random number seed = %d\n", seed );

  for ( i = 0; i < N; i++ )
  {
    x[i] = r4_uniform_01 ( &seed );
  }

  printf ( "\n" );
  printf ( "  First few values:\n" );
  printf ( "\n" );
  for ( i = 0; i < 10; i++ )
  {
    printf ( "  %6d  %14f\n", i, x[i] );
  }
  min = r4vec_min ( N, x );
  max = r4vec_max ( N, x );
  mean = r4vec_mean ( N, x );
  variance = r4vec_variance ( N, x );

  printf ( "\n" );
  printf ( "  Number of samples was %d\n", N );
  printf ( "  Minimum value was %f\n", min );
  printf ( "  Maximum value was %f\n", max );
  printf ( "  Average value was %f\n", mean );
  printf ( "  Variance was      %f\n", variance );

  return;
# undef N
}
/******************************************************************************/

void test12555 ( )

/******************************************************************************/
/*
  Purpose:

    TEST12555 tests R4VEC_INDICATOR0_NEW.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    27 September 2014

  Author:

    John Burkardt
*/
{
  int n;
  float *v;

  printf ( "\n" );
  printf ( "TEST12555\n" );
  printf ( "  R4VEC_INDICATOR0_NEW returns an indicator vector.\n" );

  n = 10;
  v = r4vec_indicator0_new ( n );
  r4vec_print ( n, v, "  Indicator0 vector:" );
  free ( v );

  return;
# undef N
}
