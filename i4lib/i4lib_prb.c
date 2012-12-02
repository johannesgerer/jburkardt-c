# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "i4lib.h"

int main ( void );

void test01 ( void );
void test02 ( void );
void test03 ( void );
void test04 ( void );
void test05 ( void );
void test06 ( void );
void test07 ( void );
void test08 ( void );
void test09 ( void );

void test10 ( void );
void test11 ( void );
void test12 ( void );
void test13 ( void );
void test14 ( void );
void test15 ( void );
void test16 ( void );
void test17 ( void );
void test18 ( void );
void test19 ( void );

void test20 ( void );
void test21 ( void );
void test22 ( void );
void test23 ( void );
void test24 ( void );
void test243 ( void );
void test245 ( void );
void test25 ( void );
void test26 ( void );
void test27 ( void );
void test28 ( void );
void test29 ( void );

void test30 ( void );
void test31 ( void );
void test32 ( void );
void test33 ( void );
void test335 ( void );

void test50 ( void );

void test602 ( void );
void test605 ( void );

void test73 ( void );
void test76 ( void );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    I4LIB_PRB calls the I4LIB tests.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    25 June 2012

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "I4LIB_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the I4LIB library.\n" );

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
  test23 ( );
  test24 ( );
  test243 ( );
  test245 ( );
  test25 ( );
  test26 ( );
  test27 ( );
  test28 ( );
  test29 ( );

  test30 ( );
  test31 ( );
  test32 ( );
  test33 ( );
  test335 ( );

  test50 ( );

  test602 ( );
  test605 ( );

  test73 ( );
  test76 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "I4LIB_PRB\n" );
  printf ( "  Normal end of execution.\n" );

  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void test01 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST01 tests I4_BIT_HI1.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 July 2010

  Author:

    John Burkardt
*/
{
  int i;
  int j;
  int seed = 123456789;
  int test;
  int test_num = 10;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  I4_BIT_HI1 returns the location of the high 1 bit.\n" );
  printf ( "\n" );
  printf ( "       I  I4_BIT_HI1(I)\n" );
  printf ( "\n" );

  for ( test = 1; test <= test_num; test++ )
  {
    i = i4_uniform_ab ( 0, 100, &seed );
    j = i4_bit_hi1 ( i );
    printf ( "  %6d  %6d\n", i, j );
  }

  return;
}
/******************************************************************************/

void test02 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST02 tests I4_BIT_LO0.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    27 July 2010

  Author:

    John Burkardt
*/
{
  int i;
  int j;
  int seed = 123456789;
  int test;
  int test_num = 10;

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  I4_BIT_LO0 returns the location of the low 0 bit.\n" );
  printf ( "\n" );
  printf ( "       I  I4_BIT_LO0(I)\n" );
  printf ( "\n" );

  for ( test = 1; test <= test_num; test++ )
  {
    i = i4_uniform_ab ( 0, 100, &seed );
    j = i4_bit_lo0 ( i );
    printf ( "  %6d  %6d\n", i, j );
  }

  return;
}
/******************************************************************************/

void test03 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST03 tests I4_BIT_LO1.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    29 July 2010

  Author:

    John Burkardt
*/
{
  int i;
  int j;
  int seed = 123456789;
  int test;
  int test_num = 10;

  printf ( "\n" );
  printf ( "TEST03\n" );
  printf ( "  I4_BIT_LO1 returns the location of the low 1 bit.\n" );
  printf ( "\n" );
  printf ( "       I  I4_BIT_LO1(I)\n" );
  printf ( "\n" );

  for ( test = 1; test <= test_num; test++ )
  {
    i = i4_uniform_ab ( 0, 100, &seed );
    j = i4_bit_lo1 ( i );
    printf ( "  %6d  %6d\n", i, j );
  }

  return;
}
/******************************************************************************/

void test04 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST04 tests I4_BIT_REVERSE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    29 July 2010

  Author:

    John Burkardt
*/
{
  int i;
  int i_hi;
  int j;
  int k;

  printf ( "\n" );
  printf ( "TEST04\n" );
  printf ( "  I4_BIT_REVERSE bit reverses I with respect to 2^J\n" );
  printf ( "\n" );
  printf ( "         I         J  I4_BIT_REVERSE(I,J)\n" );
  printf ( "\n" );

  for ( j = 0; j <= 4; j++ )
  {
    i_hi = i4_power ( 2, j ) - 1;
    for ( i = 0; i <= i_hi; i++ )
    {
      k = i4_bit_reverse ( i, j );
      printf ( "  %8d  %8d  %8d\n", i, j, k );
    }
  }
  return;
}
/******************************************************************************/

void test05 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST05 tests I4_CHARACTERISTIC.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    29 July 2010

  Author:

    John Burkardt
*/
{
  int i;

  printf ( "TEST05\n" );
  printf ( "  I4_CHARACTERISTIC computes the characteristic\n" );
  printf ( "  of an integer Q, which is  \n" );
  printf ( "    Q if Q is prime;\n" );
  printf ( "    P, if Q = P**N for some prime P;\n" );
  printf ( "    0, if Q is negative, 0, 1, or the product of \n" );
  printf ( "      more than 1 distinct prime.\n" );
  printf ( "\n" );
  printf ( "  I, I4_CHARACTERISTIC\n" );
  printf ( "\n" );

  for ( i = 1; i <= 50; i++)
  {
    printf ( "  %2d  %4d\n", i, i4_characteristic ( i ) );
  }

  return;
}
/******************************************************************************/

void test06 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST06 tests I4_DIV_ROUNDED.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    29 July 2010

  Author:

    John Burkardt
*/
{
  int a;
  int a_hi =  100;
  int a_lo = -100;
  int b;
  int b_hi =  10;
  int b_lo = -10;
  double c0;
  int c1;
  int c2;
  int c3;
  int c4;
  int seed;
  int test;
  int test_num = 20;

  printf ( "\n" );
  printf ( "TEST06\n" );
  printf ( "  I4_DIV_ROUNDED performs rounded integer division.\n" );
  printf ( "\n" );
  printf ( "  C0 = ( double ) ( a ) / ( double ) ( b )\n" );
  printf ( "  C1 = I4_DIV_ROUNDED ( A, B )\n" );
  printf ( "  C2 = r8_nint ( ( double ) ( a ) / ( double ) ( b ) )\n" );
  printf ( "  C3 = A / B\n" );
  printf ( "  C4 = ( int ) ( ( double ) ( a ) / ( double ) ( b ) )\n" );
  printf ( "\n" );
  printf ( "  C1 and C2 should be equal;\n" );
  printf ( "  C3 and C4 should be equal.\n" );
  printf ( "\n" );
  printf ( "     A     B           C0         C1    C2      C3    C4\n" );
  printf ( "\n" );

  seed = 123456789;

  for ( test = 1; test <= test_num; test++ )
  {
    a = i4_uniform_ab ( a_lo, a_hi, &seed );
    b = i4_uniform_ab ( b_lo, b_hi, &seed );
    if ( b == 0 )
    {
      b = 7;
    }
    c0 = ( double ) ( a ) / ( double ) ( b );
    c1 = i4_div_rounded ( a, b );
    c2 = r8_nint ( ( double ) ( a ) / ( double ) ( b ) );
    c3 = a / b;
    c4 = ( int ) ( ( double ) ( a ) / ( double ) ( b ) );
    printf ( "  %4d  %4d  %10.4f  %4d  %4d  %4d  %4d\n",
      a, b, c0, c1, c2, c3, c4 );
  }
  return;
}
/******************************************************************************/

void test07 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST07 tests I4_DIVP.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    29 July 2010

  Author:

    John Burkardt
*/
{
  int a;
  int a_hi =  100;
  int a_lo = -100;
  int b;
  int b_hi =  10;
  int b_lo = -10;
  int c;
  int d;
  int seed;
  int test;
  int test_num = 20;

  printf ( "\n" );
  printf ( "TEST07\n" );
  printf ( "  I4_DIVP(A,B) returns the smallest multiplier of J\n" );
  printf ( "  that is less than I\n" );
  printf ( "\n" );
  printf ( "     A     B     C     D\n" );
  printf ( "\n" );

  seed = 123456789;

  for ( test = 1; test <= test_num; test++ )
  {
    a = i4_uniform_ab ( a_lo, a_hi, &seed );
    b = i4_uniform_ab ( b_lo, b_hi, &seed );
    if ( b == 0 )
    {
      b = 7;
    }
    c = i4_divp ( a, b );
    d = c * b;
    printf (  "  %4d  %4d  %4d  %4d\n", a, b, c, d );
  }

  return;
}
/******************************************************************************/

void test08 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST08 tests I4_GCD.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    29 July 2010

  Author:

    John Burkardt
*/
{
# define TEST_NUM 7

  int i;
  int i_test[TEST_NUM] = { 36, 49, 0, 12, 36, 1, 91 };
  int j;
  int j_test[TEST_NUM] = { 30, -7, 71, 12, 49, 42, 28 };
  int test;

  printf ( "\n" );
  printf ( "TEST08\n" );
  printf ( "  I4_GCD computes the greatest common factor,\n" );
  printf ( "\n" );
  printf ( "     I     J   I4_GCD\n" );
  printf ( "\n" );

  for ( test = 0; test < TEST_NUM; test++ )
  {
    i = i_test[test];
    j = j_test[test];
    printf ( "  %6d  %6d  %6d\n", i, j, i4_gcd ( i, j ) );
  }

  return;
# undef TEST_NUM
}
/******************************************************************************/

void test09 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST09 tests I4_HUGE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 January 2007

  Author:

    John Burkardt
*/
{
  printf ( "\n" );
  printf ( "TEST09\n" );
  printf ( "  I4_HUGE returns a huge integer.\n" );
  printf ( "\n" );
  printf ( "  I4_HUGE() = %d\n", i4_huge ( ) );

  return;
}
/******************************************************************************/

void test10 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST10 tests I4_HUGE_NORMALIZER.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 January 2007

  Author:

    John Burkardt
*/
{
  int i4;
  double r8;
  double value;

  i4 = i4_huge ( );
  r8 = i4_huge_normalizer ( );

  printf ( "\n" );
  printf ( "TEST10\n" );
  printf ( "  I4_HUGE_NORMALIZER returns 1/(I4_HUGE+1).\n" );
  printf ( "\n" );
  printf ( "  I4_HUGE() = %d\n", i4 );
  printf ( "  I4_HUGE_NORMALIZER() = %e\n", r8 );

  value = ( ( double ) ( i4 ) ) * r8;

  printf ( "\n" );
  printf ( "  I4_HUGE * I4_HUGE_NORMALIZER = %e\n", value );

  return;
}
/******************************************************************************/

void test11 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST11 tests I4_IS_PRIME.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    29 July 2010

  Author:

    John Burkardt
*/
{
  int i;

  printf ( "\n" );
  printf ( "TEST11\n" );
  printf ( "  I4_IS_PRIME reports whether an integer is prime.\n" );
  printf ( "\n" );
  printf ( "  I     I4_IS_PRIME(I)\n" );
  printf ( "\n" );

  for ( i = -2; i <= 25; i++ )
  {
    printf ( "   %6d  %6d\n", i, i4_is_prime ( i ) );
  }

  return;
}
/******************************************************************************/

void test12 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST12 tests I4_LCM.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    29 July 2010

  Author:

    John Burkardt
*/
{
# define TEST_NUM 7

  int i;
  int i_test[TEST_NUM] = { 36, 49,  0, 12, 36,  1, 91 };
  int j;
  int j_test[TEST_NUM] = { 30, -7, 71, 12, 49, 42, 28 };
  int test;

  printf ( "\n" );
  printf ( "TEST12\n" );
  printf ( "  I4_LCM computes the least common multiple.\n" );
  printf ( "\n" );
  printf ( "     I     J   I4_LCM\n" );
  printf ( "\n" );

  for ( test = 0; test < TEST_NUM; test++ )
  {
    i = i_test[test];
    j = j_test[test];
    printf ( "  %6d  %6d  %6d\n", i, j, i4_lcm ( i, j ) );
  }

  return;
# undef TEST_NUM
}
/******************************************************************************/

void test13 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST13 tests I4_LOG_10.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    29 July 2010

  Author:

    John Burkardt
*/
{
# define N 13

  int i;
  int x[N] = { 0, 1, 2, 3, 9, 10, 11, 99, 101, -1, -2, -3, -9 };

  printf ( "\n" );
  printf ( "TEST13\n" );
  printf ( "  I4_LOG_10: whole part of log base 10,\n" );
  printf ( "\n" );
  printf ( "  X, I4_LOG_10\n" );
  printf ( "\n" );

  for ( i = 0; i < N; i++ )
  {
    printf ( "  %6d  %6d\n", x[i], i4_log_10 ( x[i] ) );
  }

  return;
# undef N
}
/******************************************************************************/

void test14 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST14 tests I4_LOG_2.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    29 July 2010

  Author:

    John Burkardt
*/
{
# define TEST_NUM 17

  int test;
  int x;
  int x_test[TEST_NUM] = {
      0,    1,    2,    3,    9,
     10,   11,   99,  101,   -1,
     -2,   -3,   -9, 1000, 1023,
   1024, 1025 };

  printf ( "\n" );
  printf ( "TEST14\n" );
  printf ( "  I4_LOG_2: whole part of log base 2.\n" );
  printf ( "\n" );
  printf ( "       X     I_LOG_2\n" );
  printf ( "\n" );

  for ( test = 0; test < TEST_NUM; test++ )
  {
    x = x_test[test];
    printf ( "  %6d  %12d\n", x, i4_log_2 ( x ) );
  }

  return;
# undef TEST_NUM
}
/******************************************************************************/

void test15 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST15 tests I4_LOG_I4.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    29 August 2010

  Author:

    John Burkardt
*/
{
  int i4;
  int j4;

  printf ( "\n" );
  printf ( "TEST15\n" );
  printf ( "  I4_LOG_I4: whole part of log base B,\n" );
  printf ( "\n" );
  printf ( "        I4        J4 I4_LOG_I4\n" );
  printf ( "\n" );

  for ( j4 = 2; j4 <= 5; j4++ )
  {
    for ( i4 = 0; i4 <= 10; i4++ )
    {
      printf ( "  %8d  %8d  %8d\n", i4, j4, i4_log_i4 ( i4, j4 ) );
    }
    printf ( "\n" );
  }

  return;
}
/******************************************************************************/

void test16 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST16 tests I4_LOG_R8.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 July 2010

  Author:

    John Burkardt
*/
{
# define TEST_NUM 10

  double b;
  double b_test[TEST_NUM] = {
    2.0, 3.0,  4.0,  5.0,   6.0,
    7.0, 8.0, 16.0, 32.0, 256.0 };
  int test;
  int x;

  x = 16;

  printf ( "\n" );
  printf ( "TEST16\n" );
  printf ( "  I4_LOG_R8: whole part of log base B,\n" );
  printf ( "\n" );
  printf ( "  X  B  I4_LOG_R8\n" );
  printf ( "\n" );

  for ( test = 0; test < TEST_NUM; test++ )
  {
    b = b_test[test];

    printf ( "  %6d  %14f  %12d\n", x, b, i4_log_r8 ( x, b ) );
  }

  return;
# undef TEST_NUM
}
/******************************************************************************/

void test17 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST17 tests I4_MANT.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 July 2010

  Author:

    John Burkardt
*/
{
  int is;
  int j;
  int k;
  int l;
  double x;

  x = -314.159;

  printf ( "\n" );
  printf ( "TEST17\n" );
  printf ( "  I4_MANT decomposes an integer,\n" );
  printf ( "\n" );
  printf ( "  Number to be decomposed is X = %f\n", x );

  i4_mant ( x, &is, &j, &k, &l );

  printf ( "\n" );
  printf ( "  X = %d * ( %d / %d ) * 2 ^ %d\n", is, j, k, l );

  return;
}
/******************************************************************************/

void test18 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST18 tests I4_MODDIV;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 July 2010

  Author:

    John Burkardt
*/
{
# define TEST_NUM 4

  int ndivid[TEST_NUM] = { 50, -50, 50, -50 };
  int nmult;
  int nrem;
  int number[TEST_NUM] = { 107, 107, -107, -107 };
  int test;

  printf ( "\n" );
  printf ( "TEST18\n" );
  printf ( "  I4_MODDIV factors a number\n" );
  printf ( "  into a multiple and a remainder.\n" );
  printf ( "\n" );
  printf ( "    Number   Divisor  Multiple Remainder\n" );
  printf ( "\n" );

  for ( test = 0; test < TEST_NUM; test++ )
  {
    i4_moddiv ( number[test], ndivid[test], &nmult, &nrem );

    printf ( "  %10d  %10d  %10d  %10d\n",
      number[test], ndivid[test], nmult, nrem );
  }

  printf ( "\n" );
  printf ( "  Repeat using C percent\n" );
  printf ( "\n" );

  for ( test = 0; test < TEST_NUM; test++ )
  {
    nrem = ( number[test] % ndivid[test] );
    nmult = number[test] / ndivid[test];

    printf ( "  %10d  %10d  %10d  %10d\n",
      number[test], ndivid[test], nmult, nrem );
  }

  return;
# undef TEST_NUM
}
/******************************************************************************/

void test19 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST19 tests I4_MODP.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 July 2010

  Author:

    John Burkardt
*/
{
# define TEST_NUM 4

  int ndivid[TEST_NUM] = { 50, -50, 50, -50 };
  int nmult;
  int nrem;
  int number[TEST_NUM] = { 107, 107, -107, -107 };
  int test;

  printf ( "\n" );
  printf ( "TEST19\n" );
  printf ( "  I4_MODP factors a number\n" );
  printf ( "  into a multiple and a remainder.\n" );
  printf ( "\n" );
  printf ( "    Number   Divisor  Multiple Remainder\n" );
  printf ( "\n" );

  for ( test = 0; test < TEST_NUM; test++ )
  {
    nrem = i4_modp ( number[test], ndivid[test] );
    nmult = number[test] / ndivid[test];

    printf ( "  %10d  %10d  %10d  %10d\n",
      number[test], ndivid[test], nmult, nrem );
  }

  printf ( "\n" );
  printf ( "  Repeat using C percent operator:\n" );
  printf ( "\n" );

  for ( test = 0; test < TEST_NUM; test++ )
  {
    nrem = ( number[test] % ndivid[test] );
    nmult = number[test] / ndivid[test];

    printf ( "  %10d  %10d  %10d  %10d\n",
      number[test], ndivid[test], nmult, nrem );
  }

  return;
# undef TEST_NUM
}
/******************************************************************************/

void test20 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST20 tests I4_SIGN.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 July 2010

  Author:

    John Burkardt
*/
{
# define TEST_NUM 5

  int test;
  int x;
  int x_test[TEST_NUM] = { -10, -7, 0, 5, 9 };

  printf ( "\n" );
  printf ( "TEST20\n" );
  printf ( "  I4_SIGN returns the sign of a number.\n" );
  printf ( "\n" );

  for ( test = 0; test < TEST_NUM; test++ )
  {
    x = x_test[test];
    printf ( "  %6d  %6d\n", x, i4_sign ( x ) );
  }

  return;
# undef TEST_NUM
}
/******************************************************************************/

void test21 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST21 tests I4_SWAP.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 July 2010

  Author:

    John Burkardt
*/
{
  int i;
  int j;

  printf ( "\n" );
  printf ( "TEST21\n" );
  printf ( "  I4_SWAP swaps two integers.\n" );

  i = 1;
  j = 202;

  printf ( "\n" );
  printf ( "  Before swapping: \n" );
  printf ( "\n" );
  printf ( "    I = %d\n", i );
  printf ( "    J = %d\n", j );

  i4_swap ( &i, &j );

  printf ( "\n" );
  printf ( "  After swapping: \n" );
  printf ( "\n" );
  printf ( "    I = %d\n", i );
  printf ( "    J = %d\n", j );

  return;
}
/******************************************************************************/

void test22 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST22 tests I4_WALSH_1D;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 July 2010

  Author:

    John Burkardt
*/
{
  int i;
  int w0;
  int wm1;
  int wm2;
  int wm3;
  int wp1;
  int wp2;
  double x;

  printf ( "\n" );
  printf ( "TEST22\n" );
  printf ( "  I4_WALSH_1D evaluates 1D Walsh functions:\n" );
  printf ( "\n" );
  printf ( "X  W(+2) W(+1) W(0) W(-1) W(-2) W(-3)\n" );
  printf ( "\n" );

  for ( i = 0; i <= 32; i++ )
  {
    x = ( double ) i / 4.0;

    wp2 = i4_walsh_1d ( x,  2 );
    wp1 = i4_walsh_1d ( x,  1 );
    w0  = i4_walsh_1d ( x,  0 );
    wm1 = i4_walsh_1d ( x, -1 );
    wm2 = i4_walsh_1d ( x, -2 );
    wm3 = i4_walsh_1d ( x, -3 );

    printf ( "  %10.4f  %10d  %10d  %10d  %10d  %10d  %10d\n",
      x, wp2, wp1, w0, wm1, wm2, wm3 );
  }

  return;
}
/******************************************************************************/

void test23 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST23 tests I4_WRAP.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 July 2010

  Author:

    John Burkardt
*/
{
  int i;
  int ihi = 8;
  int ilo = 4;

  printf ( "\n" );
  printf ( "TEST23\n" );
  printf ( "  I4_WRAP forces an integer to lie within given limits.\n" );
  printf ( "\n" );
  printf ( "  ILO = %d\n", ilo );
  printf ( "  IHI = %d\n", ihi );
  printf ( "\n" );
  printf ( "     I  I4_WRAP(I)\n" );
  printf ( "\n" );

  for ( i = -10; i <= 20; i++ )
  {
    printf ( "  %6d  %6d\n", i, i4_wrap ( i, ilo, ihi )  );
  }

  return;
}
/******************************************************************************/

void test24 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST24 tests I4_XOR.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 January 2010

  Author:

    John Burkardt
*/
{
  int i;
  int i_lo = 0;
  int i_hi = 100;
  int j;
  int k;
  int l;
  int seed;
  int test;
  int test_num = 10;

  seed = 123456789;

  printf ( "\n" );
  printf ( "TEST24\n" );
  printf ( "  I4_XOR returns the bitwise exclusive OR of\n" );
  printf ( "  two I4's.\n" );
  printf ( "  The operator ^ should generally be used instead.\n" );
  printf ( "\n" );
  printf ( "       I       J  I4_XOR     I^J\n" );
  printf ( "\n" );

  for ( test = 1; test <= test_num; test++ )
  {
    i = i4_uniform_ab ( i_lo, i_hi, &seed );
    j = i4_uniform_ab ( i_lo, i_hi, &seed );
    k = i4_xor ( i, j );
    l = i ^ j;

    printf ( "  %6d  %6d  %6d  %6d\n", i, j, k, l );
  }

  return;
}
/******************************************************************************/

void test243 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST243 tests I4BLOCK_NEW and I4BLOCK_DELETE.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    02 March 2012

  Author:

    John Burkardt
*/
{
  int ***a;
  int i;
  int j;
  int k;
  int l;
  int m;
  int n;

  printf ( "\n" );
  printf ( "TEST243:\n" );
  printf ( "  I4BLOCK_NEW dynamically creates a 3D array.\n" );
  printf ( "  I4BLOCK_DELETE deletes it.\n" );
  printf ( "  Array entries can be addressed using the\n" );
  printf ( "  notation \"a[i][j][k]\".\n" );
//
//  These dimensions could be entered by the user; they could depend on
//  some other calculation; or they could be changed repeatedly during this
//  computation, as long as old memory is deleted by I4BLOCK_DELETE and new memory
//  requested by I4BLOCK_NEW.
//
  l = 2;
  m = 3;
  n = 2;
//
//  Allocate memory.
//
  printf ( "\n" );
  printf ( "  Allocating memory for array A of size %d by %d by %d.\n", l, m, n );

  a = i4block_new ( l, m, n );

  printf ( "\n" );
  printf ( "  Assigning values to A.\n" );
//
//  Store values in A.
//
  for ( i = 0; i < l; i++ )
  {
    for ( j = 0; j < m; j++ )
    {
      for ( k = 0; k < n; k++ )
      {
        a[i][j][k] = 100 * i + 10 * j + k;
      }
    }
  }
//
//  Print A.
//
  printf ( "\n" );
  printf ( "  Dynamically allocated matrix A:\n" );
  printf ( "\n" );
  for ( i = 0; i < l; i++ )
  {
    for ( j = 0; j < m; j++ )
    {
      for ( k = 0; k < n; k++ )
      {
        printf ( "  %8d", a[i][j][k] );
      }
      printf ( "\n" );
    }
    printf ( "\n" );
  }
//
//  Free memory.
//
  i4block_delete ( a, l, m, n );

  return;
}
/******************************************************************************/

void test245 ( )

/******************************************************************************/
/*
  Purpose:

    TEST245 tests I4BLOCK_PRINT.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    24 June 2012

  Author:

    John Burkardt
*/
{
  int l = 4;
  int m = 3;
  int n = 2;
  int x[4*3*2] = {
        1,  2,  3,   4,  1, 
        4,  9, 16,   1,  8, 
       27, 64,  2,   4,  6, 
        8,  2,  8,  18, 32, 
        2, 16, 54, 128 };

  printf ( "\n" );
  printf ( "TEST245\n" );
  printf ( "  I4BLOCK_PRINT prints an I4BLOCK.\n" );

  i4block_print ( l, m, n, x, "  The 3D array:" );

  return;
}
/******************************************************************************/

void test25 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST25 tests I4COL_FIND_ITEM.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    23 December 2010

  Author:

    John Burkardt
*/
{
# define M 5
# define N 4
# define TEST_NUM 3

  int a[M*N];
  int col;
  int i;
  int item;
  int item_test[TEST_NUM] = { 34, 12, 90 };
  int j;
  int row;
  int test;

  printf ( "\n" );
  printf ( "TEST25\n" );
  printf ( "  I4COL_FIND_ITEM finds the first occurrence of\n" );
  printf ( "  an item in an integer array of columns.\n" );

  for ( i = 0; i < M; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      a[i+j*M] = 10 * ( i + 1 ) + ( j + 1 );
    }
  }

  i4mat_print ( M, N, a, "  The matrix of columns:" );

  for ( test = 0; test < TEST_NUM; test++ )
  {
    item = item_test[test];

    i4col_find_item ( M, N, a, item, &row, &col );

    printf ( "  Item %d occurs in row %d and column %d\n", item, row, col );
  }

  return;
# undef M
# undef N
# undef TEST_NUM
}
/******************************************************************************/

void test26 ( )

/******************************************************************************/
/*
  Purpose:

    TEST26 tests I4COL_FIND_PAIR_WRAP.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 March 2012

  Author:

    John Burkardt
*/
{
# define M 5
# define N 4
# define TEST_NUM 5

  int a[M*N];
  int col;
  int i;
  int item1;
  int item1_test[TEST_NUM] = { 22, 32, 22, 54, 54 };
  int item2;
  int item2_test[TEST_NUM] = { 32, 22, 23, 14, 11 };
  int j;
  int row;
  int test;

  printf ( "\n" );
  printf ( "TEST26\n" );
  printf ( "  I4COL_FIND_PAIR_WRAP finds the first occurrence of\n" );
  printf ( "  a pair of item in an integer array of columns.\n" );
  printf ( "  Items in the array are ordered by column, and\n" );
  printf ( "  wraparound is allowed.\n" );

  for ( i = 0; i < M; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      a[i+j*M] = 10 * ( i + 1 ) + ( j + 1 );
    }
  }

  i4mat_print ( M, N, a, "  The matrix of columns:" );

  for ( test = 0; test < TEST_NUM; test++ )
  {
    item1 = item1_test[test];
    item2 = item2_test[test];

    i4col_find_pair_wrap ( M, N, a, item1, item2, &row, &col );

    printf ( "  Item %d followed by item %d occurs in row %d and column %d\n", 
      item1, item2, row, col );
  }

  return;
# undef M
# undef N
# undef TEST_NUM
}
/******************************************************************************/

void test27 ( )

/******************************************************************************/
/*
  Purpose:

    TEST27 tests I4COL_SORT_A and I4COL_SORT_D.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 March 2012

  Author:

    John Burkardt
*/
{
  int *a;
  int b = 1;
  int c = 10;
  int m = 5;
  int n = 4;
  int seed;

  printf ( "\n" );
  printf ( "TEST27\n" );
  printf ( "  I4COL_SORT_A ascending sorts an integer array\n" );
  printf ( "  as a table of columns.\n" );
  printf ( "  I4COL_SORT_D descending sorts an integer array\n" );
  printf ( "  as a table of columns.\n" );

  seed = 123456789;

  a = i4mat_uniform_new ( m, n, b, c, &seed );

  i4mat_print ( m, n, a, "  The original matrix:" );

  i4col_sort_a ( m, n, a );

  i4mat_print ( m, n, a, "  Ascending sorted:" );

  i4col_sort_d ( m, n, a );

  i4mat_print ( m, n, a, "  Descending sorted:" );

  free ( a );

  return;
}
/******************************************************************************/

void test28 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST28 tests I4COL_SORT2_A;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 March 2012

  Author:

    John Burkardt
*/
{
  int *a;
  int b = 0;
  int c = 20;
  int m = 6;
  int n = 4;
  int seed = 123456789;

  printf ( "\n" );
  printf ( "TEST28\n" );
  printf ( "  For a rectangular integer matrix:\n" );
  printf ( "  I4COL_SORT2_D sorts the elements of the columns.\n" );

  a = i4mat_uniform_new ( m, n, b, c, &seed );

  i4mat_print ( m, n, a, "  The matrix:" );

  i4col_sort2_a ( m, n, a );

  i4mat_print ( m, n, a, "  The element-sorted column matrix:" );

  free ( a );

  return;
}
/******************************************************************************/

void test29 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST29 tests I4COL_SORTED_SINGLETON_COUNT;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 March 2012

  Author:

    John Burkardt
*/
{
  int *a;
  int b;
  int c;
  int m = 3;
  int n = 10;
  int seed;
  int singleton_num;
  int test;
  int test_num = 2;

  printf ( "\n" );
  printf ( "TEST29\n" );
  printf ( "  I4COL_SORTED_SINGLETON_COUNT counts singletons\n" );
  printf ( "  in a sorted ICOL;\n" );

  seed = 123456789;

  for ( test = 1; test <= test_num; test++ )
  {
    b = 0;
    c = 3;

    a = i4mat_uniform_new ( m, n, b, c, &seed );

    i4col_sort_a ( m, n, a );

    i4mat_print ( m, n, a, "  Ascending sorted ICOL:" );

    singleton_num = i4col_sorted_singleton_count ( m, n, a );

    printf ( "\n" );
    printf ( "  Number of singletons = %d\n", singleton_num );

    free ( a );
  }

  return;
}
/******************************************************************************/

void test30 ( )

/******************************************************************************/
/*
  Purpose:

    TEST30 tests I4COL_SORTED_UNIQUE_COUNT;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 March 2012

  Author:

    John Burkardt
*/
{
  int *a;
  int b;
  int c;
  int m = 3;
  int n = 10;
  int seed;
  int unique_num;
  int test;
  int test_num = 2;

  printf ( "\n" );
  printf ( "TEST30\n" );
  printf ( "  I4COL_SORTED_UNIQUE_COUNT counts the unique entries\n" );
  printf ( "  of a sorted ICOL;\n" );

  seed = 123456789;

  for ( test = 1; test <= test_num; test++ )
  {
    b = 0;
    c = 3;

    a = i4mat_uniform_new ( m, n, b, c, &seed );

    i4col_sort_a ( m, n, a );

    i4mat_print ( m, n, a, "  Ascending sorted ICOL:" );

    unique_num = i4col_sorted_unique_count ( m, n, a );

    printf ( "\n" );
    printf ( "  Number of unique entries = %d\n", unique_num );

    free ( a );
  }

  return;
}
/******************************************************************************/

void test31 ( )

/******************************************************************************/
/*
  Purpose:

    TEST31 tests I4MAT_ELIM and I4MAT_RED.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 March 2012

  Author:

    John Burkardt
*/
{
# define M 5
# define N 5

  int a[M*N];
  int col[N];
  int factor;
  int i;
  int j;
  int k;
  int row[M];
  int test;
  int test_num = 3;

  printf ( "\n" );
  printf ( "TEST31\n" );
  printf ( "  I4MAT_ELIM does exact Gauss elimination.\n" );
  printf ( "  I4MAT_RED divides common factors in a matrix;\n" );

  for ( test = 1; test <= test_num; test++ )
  {
    if ( test == 1 )
    {
      k = 0;
      for ( i = 0; i < M; i++ )
      {
        for ( j = 0; j < N; j++ )
        {
          k = k + 1;
          a[i+j*M] = k;
        }
      }
    }
    else if ( test == 2 )
    {
      factor = 8 * 7 * 6 * 5 * 4 * 3 * 2;

      for ( i = 0; i < M; i++ )
      {
        for ( j = 0; j < N; j++ )
        {
          a[i+j*M] = factor / ( i + j + 1 );
        }
      }
    }
    else if ( test == 3 )
    {
      for ( i = 0; i < M; i++ )
      {
        for ( j = 0; j < N; j++ )
        {
          a[i+j*M] = ( i + 1 ) * ( j + 1 );
        }
      }
    }

    i4mat_print ( M, N, a, "  The original matrix:" );

    i4mat_red ( M, N, a, row, col );

    printf ( "\n" );
    printf ( "  The matrix, as returned by I4MAT_RED:\n" );
    printf ( "  (Factors are displayed in an extra row and column.)\n" );
    printf ( "\n" );
    for ( i = 0; i < M; i++ )
    {
      for ( j = 0; j < N; j++ )
      {
        printf ( "  %6d", a[i+j*M] );
      }
      printf ( "  %6d\n", row[i] );
    }
    for ( j = 0; j < N; j++ )
    {
      printf ( "  %6d", col[j] );
    }
    printf ( "\n" );

    i4mat_elim ( M, N, a );

    i4mat_print ( M, N, a, "  The matrix returned by I4MAT_ELIM:" );
  }

  return;
# undef M
# undef N
}
/******************************************************************************/

void test32 ( )

/******************************************************************************/
/*
  Purpose:

    TEST32 tests I4MAT_MAX_INDEX and I4MAT_MIN_INDEX.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 March 2012

  Author:

    John Burkardt
*/
{
  int *a;
  int b = 0;
  int c = 10;
  int i;
  int j;
  int m = 5;
  int n = 7;
  int seed;

  printf ( "\n" );
  printf ( "TEST32\n" );
  printf ( "  I4MAT_MAX_INDEX locates the maximum;\n" );
  printf ( "  I4MAT_MIN_INDEX locates the minimum;\n" );

  seed = 123456789;

  a = i4mat_uniform_new ( m, n, b, c, &seed );

  i4mat_print ( m, n, a, "  Random array:" );

  printf ( "\n" );
  i4mat_max_index ( m, n, a, &i, &j );
  printf ( "  Maximum I,J indices            %d  %d\n", i, j );
  i4mat_min_index ( m, n, a, &i, &j );
  printf ( "  Minimum I,J indices            %d  %d\n", i, j );

  free ( a );

  return;
}
/******************************************************************************/

void test33 ( )

/******************************************************************************/
/*
  Purpose:

    TEST33 tests I4MAT_L1_INVERSE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 March 2012

  Author:

    John Burkardt
*/
{
# define N 6
/*
  Each row of this definition is a COLUMN of the matrix.
*/
  int a[N*N] = {
     1,  2,  0,  5,  0, 75,
     0,  1,  0,  0,  0,  0,
     0,  0,  1,  3,  0,  0,
     0,  0,  0,  1,  0,  6,
     0,  0,  0,  0,  1,  4,
     0,  0,  0,  0,  0,  1 };
  int *b;
  int *c;

  printf ( "\n" );
  printf ( "TEST33\n" );
  printf ( "  I4MAT_L1_INVERSE inverts a unit lower triangular matrix.\n" );

  i4mat_print ( N, N, a, "  The original matrix:" );

  b = i4mat_l1_inverse ( N, a );

  i4mat_print ( N, N, b, "  The inverse matrix:" );

  c = i4mat_mm ( N, N, N, a, b );

  i4mat_print ( N, N, c, "  Product C = A * B:" );

  free ( b );
  free ( c );

  return;
# undef N
}
/******************************************************************************/

void test335 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST335 tests I4MAT_NEW and I4MAT_DELETE.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    02 March 2012

  Author:

    John Burkardt
*/
{
  int **a;
  int **b;
  int i;
  int j;
  int k;
  int m;
  int n;

  printf ( "\n" );
  printf ( "TEST335:\n" );
  printf ( "  I4MAT_NEW dynamically creates a 2D array.\n" );
  printf ( "  I4MAT_DELETE deletes it.\n" );
  printf ( "  Array entries can be addressed using the\n" );
  printf ( "  notation \"a[i][j]\".\n" );
/*
  These dimensions could be entered by the user; they could depend on
  some other calculation; or they could be changed repeatedly during this
  computation, as long as old memory is deleted by I4MAT_DELETE and new memory
  requested by I4MAT_NEW.
*/
  m = 4;
  n = 5;
/*
  Allocate memory.
*/
  printf ( "\n" );
  printf ( "  Allocating memory for array A of size %d by %d.\n", m, n );

  a = i4mat_new ( m, n );

  printf ( "\n" );
  printf ( "  Assigning values to A.\n" );
/*
  Store values in A.
*/
  for ( i = 0; i < m; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      a[i][j] = 10 * i + j;
    }
  }
/*
  Print A.
*/
  printf ( "\n" );
  printf ( "  Dynamically allocated matrix A:\n" );
  printf ( "\n" );
  for ( i = 0; i < m; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      printf ( "  %8d", a[i][j] );
    }
    printf ( "\n" );
  }
/*
  Create a new matrix B to store A' * A.
*/
  b = i4mat_new ( n, n );

  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      b[i][j] = 0;
      for ( k = 0; k < m; k++ )
      {
        b[i][j] = b[i][j] + a[k][i] * a[k][j];
      }
    }
  }
/*
  Print the matrix.
*/
  printf ( "\n" );
  printf ( "  Dynamically allocated matrix B = A' * A:\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      printf ( "  %8d\n", b[i][j] );
    }
    printf ( "\n" );
  }
/*
  Free memory.
*/
  i4mat_delete ( a, m, n );
  i4mat_delete ( b, n, n );

  return;
}
/******************************************************************************/

void test50 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST50 tests I4VEC_CUM_NEW and I4VEC_CUM0_NEW.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    22 December 2010

  Author:

    John Burkardt
*/
{
# define N 10

  int *a;
  int *a_cum;
  int *a_cum0;
  int b;
  int c;
  int n = N;
  int seed;

  printf ( "\n" );
  printf ( "TEST50\n" );
  printf ( "  For an integer vector:\n" );
  printf ( "  I4VEC_CUM_NEW:   cumulative sum;\n" );
  printf ( "  I4VEC_CUM0_NEW:  cumulative sum, zero based;\n" );

  seed = 123456789;

  b = -n;
  c = n;

  a = i4vec_uniform_new ( n, b, c, &seed );

  i4vec_print ( n, a, "  Input vector:" );

  a_cum = i4vec_cum_new ( n, a );

  i4vec_print ( n, a_cum, "  Cumulative sums:" );

  a_cum0 = i4vec_cum0_new ( n, a );

  i4vec_print ( n + 1, a_cum0, "  0-based Cumulative sums:" );

  free ( a );
  free ( a_cum );
  free ( a_cum0 );

  return;
# undef N
}
/******************************************************************************/

void test602 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST602 tests I4VEC_INDEXED_HEAP_D;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 August 2010

  Author:

    John Burkardt
*/
{
  int a[20] = {
    101, 102, 103, 104, 105, 106, 107, 108, 109, 110,
    111, 112, 113, 114, 115, 116, 117, 118, 119, 120 };
  int i;
  int indx[10] = {
    0, 10, 16, 4, 6, 12, 14, 2, 18, 8 };
  int m = 20;
  int n = 10;

  printf ( "\n" );
  printf ( "TEST602\n" );
  printf ( "  I4VEC_INDEXED_HEAP_D creates a descending heap\n" );
  printf ( "  from an indexed vector.\n" );
/*
  Print before.
*/
  i4vec_print ( m, a, "  The data vector:" );
  i4vec_print ( n, indx, "  The index vector:" );
  printf ( "\n" );
  printf ( "  A(INDX):\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %4d  %4d\n", i, a[indx[i]]  );
  }
/*
  Heap the data.
*/
  i4vec_indexed_heap_d ( n, a, indx );
/*
  Print afterwards.  Only INDX should change.
*/
  i4vec_print ( m, a, "  The data vector (should NOT change):" );
  i4vec_print ( n, indx, "  The index vector (may change):" );
  printf ( "\n" );
  printf ( "  A(INDX) is now a descending heap:\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %4d  %4d\n", i, a[indx[i]]  );
  }

  return;
}
/******************************************************************************/

void test605 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST605 tests I4VEC_INDEXED_HEAP_D_EXTRACT and related routines.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 August 2010

  Author:

    John Burkardt
*/
{
  int *a;
  int i;
  int indx[20];
  int indx_extract;
  int indx_insert;
  int indx_max;
  int m = 20;
  int n;
  int n_max = 20;

  printf ( "\n" );
  printf ( "TEST605\n" );
  printf ( "  For an indexed I4VEC,\n" );
  printf ( "  I4VEC_INDEXED_HEAP_D_INSERT inserts a value into the heap.\n" );
  printf ( "  I4VEC_INDEXED_HEAP_D_EXTRACT extracts the maximum value;\n" );
  printf ( "  I4VEC_INDEXED_HEAP_D_MAX reports the maximum value.\n" );
  printf ( "\n" );
  printf ( "  These 3 operations are enough to model a priority queue.\n" );
/*
  Set the data array.  To keep things easy, we will use the indicator vector.
*/
  a = i4vec_indicator_new ( m );
/*
  The index array will initially be a random subset of the numbers 1 to M,
  in random order.
*/
  n = 5;
  indx[0]  =  8;
  indx[1]  =  1;
  indx[2]  =  7;
  indx[3]  = 13;
  indx[4]  =  4;
  indx[5]  =  6;
  indx[6]  = 14;
  indx[7]  =  0;
  indx[8]  = 18;
  indx[9]  = 19;
  indx[10] =  2;

  i4vec_print ( m, a, "  The data vector:" );
  i4vec_print ( n, indx, "  The index vector:" );
  printf ( "\n" );
  printf ( "  A(INDX):\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %4d  %4d\n", i, a[indx[i]]  );
  }
/*
  Create a descending heap from the indexed array.
*/
  i4vec_indexed_heap_d ( n, a, indx );

  i4vec_print ( n, indx, "  The index vector after heaping:" );
  printf ( "\n" );
  printf ( "  A(INDX) after heaping:\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %4d  %4d\n", i, a[indx[i]]  );
  }
/*
  Insert five entries, and monitor the maximum.
*/
  for ( i = 0; i < 5; i++ )
  {
    indx_insert = indx[n];

    printf ( "\n" );
    printf ( "  Inserting value %d\n", a[indx_insert] );

    i4vec_indexed_heap_d_insert ( &n, a, indx, indx_insert );

    indx_max = i4vec_indexed_heap_d_max ( n, a, indx );

    printf ( "  Current maximum is %d\n", a[indx_max] );
  }
  i4vec_print ( m, a, "  The data vector after insertions:" );
  i4vec_print ( n, indx, "  The index vector after insertions:" );
  printf ( "\n" );
  printf ( "  A(INDX) after insertions:\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %4d  %4d\n", i, a[indx[i]]  );
  }
/*
  Extract the first 5 largest elements.
*/
  printf ( "\n" );
  printf ( "  Now extract the maximum several times.\n" );
  printf ( "\n" );

  for ( i = 0; i < 5; i++ )
  {
    indx_extract = i4vec_indexed_heap_d_extract ( &n, a, indx );
    printf ( "  Extracting maximum element A[%d] = %d\n",
      indx_extract, a[indx_extract] );
  }
  i4vec_print ( m, a, "  The data vector after extractions:" );
  i4vec_print ( n, indx, "  The index vector after extractions:" );
  printf ( "\n" );
  printf ( "  A(INDX) after extractions:\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %4d  %4d\n", i, a[indx[i]]  );
  }

  free ( a );

  return;
}
/******************************************************************************/

void test73 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST73 tests I4VEC_RUN_COUNT.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 January 2007

  Author:

    John Burkardt
*/
{
  int *a;
  int j;
  int n = 20;
  int run_count;
  int seed;
  int test;
  int test_num = 10;

  printf ( "\n" );
  printf ( "TEST73\n" );
  printf ( "  I4VEC_RUN_COUNT counts runs in an I4VEC\n" );
  printf ( "\n" );
  printf ( " Run Count        Sequence\n" );
  printf ( "\n" );

  seed = 123456789;

  for ( test = 1; test <= test_num; test++ )
  {
    a = i4vec_uniform_new ( n, 0, 1, &seed );

    run_count = i4vec_run_count ( n, a );

    printf ( "  %8d", run_count );
    printf ( "        " );
    for ( j = 0; j < n; j++ )
    {
      printf ( "%2d", a[j] );
    }
    printf ( "\n" );
    free ( a );
  }

  return;
}
/******************************************************************************/

void test76 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST76 tests I4VEC_SORT_HEAP_A;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    29 August 2006

  Author:

    John Burkardt
*/
{
# define N 30

  int *a;
  int b;
  int c;
  int i;
  int seed;

  printf ( "\n" );
  printf ( "TEST76\n" );
  printf ( "  I4VEC_SORT_HEAP_A sorts an integer array;\n" );

  b = 0;
  c = N;
  seed = 123456789;

  printf ( "\n" );
  printf ( "  Using random seed %d\n", seed );

  a = i4vec_uniform_new ( N, b, c, &seed );

  i4vec_print ( N, a, "  Unsorted array:" );

  i4vec_sort_heap_a ( N, a );

  i4vec_print ( N, a, "  Sorted array:" );

  free ( a );

  return;
# undef N
}

