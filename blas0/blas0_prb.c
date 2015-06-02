# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <complex.h>

# include "blas0.h"

int main ( );
void test01 ( );
void test015 ( );
void test02 ( );
void test03 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for BLAS0_PRB.

  Discussion:

    BLAS0_PRB tests the BLAS0 library.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    20 March 2014

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "BLAS0_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the BLAS0 library.\n" );

  test01 ( );
  test015 ( );
  test02 ( );
  test03 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "BLAS0_PRB\n" );
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

    TEST01 tests R4_ABS.

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
  printf ( "TEST01\n" );
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

void test015 ( )

/******************************************************************************/
/*
  Purpose:

    TEST015 tests R4_SIGN.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    29 June 2011

  Author:

    John Burkardt
*/
{
# define TEST_NUM 5

  int test;
  float x;
  float x_test[TEST_NUM] = { -1.25, -0.25, 0.0, +0.5, +9.0 };

  printf ( "\n" );
  printf ( "TEST015\n" );
  printf ( "  R4_SIGN returns the sign of a number.\n" );
  printf ( "\n" );

  for ( test = 0; test < TEST_NUM; test++ )
  {
    x = x_test[test];
    printf ( "  %8f  %8f\n", x, r4_sign ( x ) );
  }

  return;
# undef TEST_N
}
/******************************************************************************/

void test02 ( )

/******************************************************************************/
/*
  Purpose:

    TEST02 tests R8_ABS.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    18 August 2010

  Author:

    John Burkardt
*/
{
  double r8;
  double r8_absolute;
  double r8_hi = 5.0;
  double r8_lo = -3.0;
  int seed;
  int test;
  int test_num = 10;

  seed = 123456789;

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  R8_ABS returns the absolute value of an R8.\n" );
  printf ( "\n" );
  printf ( "      X         R8_ABS(X)\n" );
  printf ( "\n" );

  for ( test = 1; test <= test_num; test++ )
  {
    r8 = r8_uniform_ab ( r8_lo, r8_hi, &seed );
    r8_absolute = r8_abs ( r8 );
    printf ( "  %10.6f  %10.6f\n", r8, r8_absolute );
  }

  return;
}
/******************************************************************************/

void test03 ( )

/******************************************************************************/
/*
  Purpose:

    TEST03 tests R8_SIGN.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    29 June 2011

  Author:

    John Burkardt
*/
{
# define TEST_NUM 5

  int test;
  double x;
  double x_test[TEST_NUM] = { -1.25, -0.25, 0.0, +0.5, +9.0 };

  printf ( "\n" );
  printf ( "TEST03\n" );
  printf ( "  R8_SIGN returns the sign of a number.\n" );
  printf ( "\n" );

  for ( test = 0; test < TEST_NUM; test++ )
  {
    x = x_test[test];
    printf ( "  %8f  %8f\n", x, r8_sign ( x ) );
  }

  return;
# undef TEST_NUM
}
