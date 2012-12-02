# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>

int main ( int argc, char *argv[] );
void test_abs ( void );
void test_acos ( void );
void test_asin ( void );
void test_atan ( void );
void test_atan2 ( void );
void test_ceil ( void );
void test_cos ( void );
void test_cosh ( void );
void test_exp ( void );
void test_fabs ( void );
void test_floor ( void );
void test_fmod ( void );
void test_frexp ( void );
void test_ldexp ( void );
void test_log ( void );
void test_log10 ( void );
void test_modf ( void );
void test_pow ( void );
void test_sin ( void );
void test_sinh ( void );
void test_sqrt ( void );
void test_tan ( void );
void test_tanh ( void );
int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
int i4_uniform ( int a, int b, int *seed );
int r4_nint ( float x );
float r4_uniform ( float b, float c, int *seed );
float r4_uniform_01 ( int *seed );
void timestamp ( void );

/******************************************************************************/

int main ( int argc, char *argv[] )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for C_INTRINSICS.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 April 2011

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "C_INTRINSICS:\n" );
  printf ( "  Test the C intrinsic library.\n" );

  test_abs ( );
  test_acos ( );
  test_asin ( );
  test_atan ( );
  test_atan2 ( );
  test_ceil ( );
  test_cos ( );
  test_cosh ( );
  test_exp ( );
  test_fabs ( );
  test_floor ( );
  test_fmod ( );
  test_frexp ( );
  test_ldexp ( );
  test_log ( );
  test_log10 ( );
  test_modf ( );
  test_pow ( );
  test_sin ( );
  test_sinh ( );
  test_sqrt ( );
  test_tan ( );
  test_tanh ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "C_INTRINSICS:\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void test_abs ( void )

/******************************************************************************/
/*
  Purpose:

    TEST_ABS tests ABS.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 April 2011

  Author:

    John Burkardt
*/
{
  int arg;
  int result;
  int seed = 123456789;
  int test;
  int test_num = 10;

  printf ( "\n" );
  printf ( "TEST_ABS:\n" );
  printf ( "  Test ABS, which evaluates the absolute value of an int.\n" );
  printf ( "\n" );
  printf ( "         I         ABS(I)\n" );
  printf ( "\n" );

  for ( test = 1; test <= test_num; test++ )
  {
    arg = i4_uniform ( -10000, +10000, &seed );
    result = abs ( arg );

    printf ( "  %10d  %10d\n", arg, result );
  }
  return;
}
/******************************************************************************/

void test_acos ( void )

/******************************************************************************/
/*
  Purpose:

    TEST_ACOS tests ACOS.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 April 2011

  Author:

    John Burkardt
*/
{
  int seed = 123456789;
  int test;
  int test_num = 10;
  float x;
  float y;
  float z;

  printf ( "\n" );
  printf ( "TEST_ACOS:\n" );
  printf ( "  Test ACOS, which evaluates the arc-cosine function.\n" );
  printf ( "\n" );
  printf ( "         X          Y           Z\n" );
  printf ( "                 ACOS(X)      COS(Y)\n" );
  printf ( "\n" );

  for ( test = 1; test <= test_num; test++ )
  {
    x = r4_uniform ( -1.0, +1.0, &seed );
    y = acos ( x );
    z = cos ( y );
    printf ( "  %10.4f  %10.4f  %10.4f\n", x, y, z );
  }
  return;
}
/******************************************************************************/

void test_asin ( void )

/******************************************************************************/
/*
  Purpose:

    TEST_ASIN tests ASIN.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 April 2011

  Author:

    John Burkardt
*/
{
  int seed = 123456789;
  int test;
  int test_num = 10;
  float x;
  float y;
  float z;

  printf ( "\n" );
  printf ( "TEST_ASIN:\n" );
  printf ( "  Test ASIN, which evaluates the arc-sine function.\n" );
  printf ( "\n" );
  printf ( "         X          Y           Z\n" );
  printf ( "                 ASIN(X)      SIN(Y)\n" );
  printf ( "\n" );

  for ( test = 1; test <= test_num; test++ )
  {
    x = r4_uniform ( -1.0, +1.0, &seed );
    y = asin ( x );
    z = sin ( y );
    printf ( "  %10.4f  %10.4f  %10.4f\n", x, y, z );
  }
  return;
}
/******************************************************************************/

void test_atan ( void )

/******************************************************************************/
/*
  Purpose:

    TEST_ATAN tests ATAN.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 April 2011

  Author:

    John Burkardt
*/
{
  int seed = 123456789;
  int test;
  int test_num = 10;
  float x;
  float y;
  float z;

  printf ( "\n" );
  printf ( "TEST_ATAN:\n" );
  printf ( "  Test ATAN, which evaluates the arc-tangent function.\n" );
  printf ( "\n" );
  printf ( "         X          Y           Z\n" );
  printf ( "                 ATAN(X)      TAN(Y)\n" );
  printf ( "\n" );

  for ( test = 1; test <= test_num; test++ )
  {
    x = r4_uniform ( -1.0, +1.0, &seed );
    y = atan ( x );
    z = tan ( y );
    printf ( "  %10.4f  %10.4f  %10.4f\n", x, y, z );
  }
  return;
}
/******************************************************************************/

void test_atan2 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST_ATAN2 tests ATAN2.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 April 2011

  Author:

    John Burkardt
*/
{
  int seed = 123456789;
  int test;
  int test_num = 10;
  float x1;
  float x2;
  float y;
  float z;

  printf ( "\n" );
  printf ( "TEST_ATAN2:\n" );
  printf ( "  Test ATAN2, which evaluates the arc-tangent function.\n" );
  printf ( "\n" );
  printf ( "         X1        X2          X1/X2        Y           Z\n" );
  printf ( "                                         ATAN(X1,X2)  TAN(Y)\n" );
  printf ( "\n" );

  for ( test = 1; test <= test_num; test++ )
  {
    x1 = r4_uniform ( -1.0, +1.0, &seed );
    x2 = r4_uniform ( -1.0, +1.0, &seed );
    y = atan2 ( x1, x2 );
    z = tan ( y );
    printf ( "  %10.4f  %10.4f  %10.4f  %10.4f  %10.4f\n", x1, x2, x1 / x2, y, z );
  }
  return;
}
/******************************************************************************/

void test_ceil ( void )

/******************************************************************************/
/*
  Purpose:

    TEST_CEIL tests CEIL.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 April 2011

  Author:

    John Burkardt
*/
{
  int seed = 123456789;
  int test;
  int test_num = 10;
  float x1;
  float x2;

  printf ( "\n" );
  printf ( "TEST_CEIL:\n" );
  printf ( "  Test CEIL, which evaluates the ceiling function.\n" );
  printf ( "\n" );
  printf ( "         X           Y\n" );
  printf ( "                    CEIL(X)\n" );
  printf ( "\n" );

  for ( test = 1; test <= test_num; test++ )
  {
    x1 = r4_uniform ( -20.0, +20.0, &seed );
    x2 = ceil ( x1 );

    printf ( "  %10.4f  %10.4f\n", x1, x2 );
  }
  return;
}
/******************************************************************************/

void test_cos ( void )

/******************************************************************************/
/*
  Purpose:

    TEST_COS tests COS.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 April 2011

  Author:

    John Burkardt
*/
{
  float pi = 3.141592653589793;
  int seed = 123456789;
  int test;
  int test_num = 10;
  float x;
  float y;

  printf ( "\n" );
  printf ( "TEST_COS:\n" );
  printf ( "  Test COS, which evaluates the cosine function.\n" );
  printf ( "\n" );
  printf ( "         X          Y\n" );
  printf ( "                 COS(X)\n" );
  printf ( "\n" );

  for ( test = 1; test <= test_num; test++ )
  {
    x = r4_uniform ( -pi, +pi, &seed );
    y = cos ( x );
    printf ( "  %10.4f  %10.4f\n", x, y );
  }
  return;
}
/******************************************************************************/

void test_cosh ( void )

/******************************************************************************/
/*
  Purpose:

    TEST_COSH tests COSH.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 April 2011

  Author:

    John Burkardt
*/
{
  float pi = 3.141592653589793;
  int seed = 123456789;
  int test;
  int test_num = 10;
  float x;
  float y;

  printf ( "\n" );
  printf ( "TEST_COSH:\n" );
  printf ( "  Test COSH, which evaluates the hyperbolic cosine function.\n" );
  printf ( "\n" );
  printf ( "         X          Y\n" );
  printf ( "                 COSH(X)\n" );
  printf ( "\n" );

  for ( test = 1; test <= test_num; test++ )
  {
    x = r4_uniform ( -5.0, +5.0, &seed );
    y = cosh ( x );
    printf ( "  %10.4f  %10.4f\n", x, y );
  }
  return;
}
/******************************************************************************/

void test_exp ( void )

/******************************************************************************/
/*
  Purpose:

    TEST_EXP tests EXP.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 April 2011

  Author:

    John Burkardt
*/
{
  int seed = 123456789;
  int test;
  int test_num = 10;
  float x;
  float y;
  float z;

  printf ( "\n" );
  printf ( "TEST_EXP:\n" );
  printf ( "  Test EXP, which evaluates the exponential function.\n" );
  printf ( "\n" );
  printf ( "         X          Y           Z\n" );
  printf ( "                  EXP(X)      LOG(Y)\n" );
  printf ( "\n" );

  for ( test = 1; test <= test_num; test++ )
  {
    x = r4_uniform ( -5.0, +10.0, &seed );
    y = exp ( x );
    z = log ( y );
    printf ( "  %10.4f  %10.4e  %10.4f\n", x, y, z );
  }
  return;
}
/******************************************************************************/

void test_fabs ( void )

/******************************************************************************/
/*
  Purpose:

    TEST_FABS tests FABS.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 April 2011

  Author:

    John Burkardt
*/
{
  int seed = 123456789;
  int test;
  int test_num = 10;
  float x;
  float y;

  printf ( "\n" );
  printf ( "TEST_FABS:\n" );
  printf ( "  Test FABS, which evaluates the absolute value of a real quantity.\n" );
  printf ( "\n" );
  printf ( "         X           Y\n" );
  printf ( "                   FABS(X)\n" );
  printf ( "\n" );

  for ( test = 1; test <= test_num; test++ )
  {
    x = r4_uniform ( -10000.0, +10000.0, &seed );
    y = fabs ( x );

    printf ( "  %10.4f  %10.4f\n", x, y );
  }
  return;
}
/******************************************************************************/

void test_floor ( void )

/******************************************************************************/
/*
  Purpose:

    TEST_FLOOR tests FLOOR.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 April 2011

  Author:

    John Burkardt
*/
{
  int seed = 123456789;
  int test;
  int test_num = 10;
  float x1;
  float x2;

  printf ( "\n" );
  printf ( "TEST_FLOOR:\n" );
  printf ( "  Test FLOOR, which evaluates the floor function.\n" );
  printf ( "\n" );
  printf ( "         X           Y\n" );
  printf ( "                    FLOOR(X)\n" );
  printf ( "\n" );

  for ( test = 1; test <= test_num; test++ )
  {
    x1 = r4_uniform ( -20.0, +20.0, &seed );
    x2 = floor ( x1 );

    printf ( "  %10.4f  %10.4f\n", x1, x2 );
  }
  return;
}
/******************************************************************************/

void test_fmod ( void )

/******************************************************************************/
/*
  Purpose:

    TEST_FMOD tests FMOD.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 April 2011

  Author:

    John Burkardt
*/
{
  int seed = 123456789;
  int test;
  int test_num = 10;
  float x1;
  float x2;
  float y;
  float z;
  float z1;
  double z2;

  printf ( "\n" );
  printf ( "TEST_FMOD:\n" );
  printf ( "  Test FMOD, which returns the remainder of X1 / X2.\n" );
  printf ( "\n" );
  printf ( "         X1        X2          X1/X2        Y           Z\n" );
  printf ( "                                         FMOD(X1,X2)  X2*MODF(X1/X2,*)\n" );
  printf ( "\n" );

  for ( test = 1; test <= test_num; test++ )
  {
    x1 = r4_uniform ( -5.0, +5.0, &seed );
    x2 = r4_uniform ( -5.0, +5.0, &seed );
    y = fmod ( x1, x2 );
    z1 = modf ( x1 / x2, &z2 ); 
    z = z1 * x2;
    printf ( "  %10.4f  %10.4f  %10.4f  %10.4f  %10.4f\n", x1, x2, x1 / x2, y, z );
  }
  return;
}
/******************************************************************************/

void test_frexp ( void )

/******************************************************************************/
/*
  Purpose:

    TEST_FREXP tests FREXP.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 April 2011

  Author:

    John Burkardt
*/
{
  int n;
  int seed = 123456789;
  int test;
  int test_num = 10;
  float x;
  float y;
  float z;

  printf ( "\n" );
  printf ( "TEST_FREXP:\n" );
  printf ( "  Test FREXP, which splits X into a normalized fraction\n" );
  printf ( "  and a power of 2.\n" );
  printf ( "\n" );
  printf ( "        X           Y              N         Z\n" );
  printf ( "                                           Y*2^N\n" );
  printf ( "\n" );

  for ( test = 1; test <= test_num; test++ )
  {
    x = r4_uniform ( -20.0, +20.0, &seed );
    y = frexp ( x, &n );
    z = y * pow ( 2.0, n );
    printf ( "  %10.4f  %10.4f  %10d  %10.4f\n", x, y, n, z );
  }
  return;
}
/******************************************************************************/

void test_ldexp ( void )

/******************************************************************************/
/*
  Purpose:

    TEST_LDEXP tests LDEXP.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 April 2011

  Author:

    John Burkardt
*/
{
  int n;
  int result;
  int seed = 123456789;
  int test;
  int test_num = 10;
  float x;
  float y;

  printf ( "\n" );
  printf ( "TEST_LDEXP:\n" );
  printf ( "  Test LDEXP, which evaluates X*2^N.\n" );
  printf ( "\n" );
  printf ( "         X           N           Y\n" );
  printf ( "                             LDEXP(X,N)\n" );
  printf ( "\n" );

  for ( test = 1; test <= test_num; test++ )
  {
    x = r4_uniform ( -1.0, +1.0, &seed );
    n = i4_uniform ( -10, +10, &seed );
    y = ldexp ( x, n );

    printf ( "  %10.4f  %10d  %10.4f\n", x, n, y );
  }
  return;
}
/******************************************************************************/

void test_log ( void )

/******************************************************************************/
/*
  Purpose:

    TEST_LOG tests LOG.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 April 2011

  Author:

    John Burkardt
*/
{
  int seed = 123456789;
  int test;
  int test_num = 10;
  float x;
  float y;
  float z;

  printf ( "\n" );
  printf ( "TEST_LOG:\n" );
  printf ( "  Test LOG, which evaluates the logarithm function.\n" );
  printf ( "\n" );
  printf ( "         X          Y           Z\n" );
  printf ( "                  LOG(X)      EXP(Y)\n" );
  printf ( "\n" );

  for ( test = 1; test <= test_num; test++ )
  {
    x = r4_uniform ( 0, +10000.0, &seed );
    y = log ( x );
    z = exp ( y );
    printf ( "  %10.4e  %10.4f  %10.4e\n", x, y, z );
  }
  return;
}
/******************************************************************************/

void test_log10 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST_LOG10 tests LOG10.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 April 2011

  Author:

    John Burkardt
*/
{
  int seed = 123456789;
  int test;
  int test_num = 10;
  float x;
  float y;
  float z;

  printf ( "\n" );
  printf ( "TEST_LOG10:\n" );
  printf ( "  Test LOG10, which evaluates the logarithm base 10 function.\n" );
  printf ( "\n" );
  printf ( "         X          Y           Z\n" );
  printf ( "                 LOG10(X)    POW(10,Y)\n" );
  printf ( "\n" );

  for ( test = 1; test <= test_num; test++ )
  {
    x = r4_uniform ( 0, +10000.0, &seed );
    y = log10 ( x );
    z = pow ( 10.0, y );
    printf ( "  %10.4e  %10.4f  %10.4e\n", x, y, z );
  }
  return;
}
/******************************************************************************/

void test_modf ( void )

/******************************************************************************/
/*
  Purpose:

    TEST_MODF tests MODF.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 April 2011

  Author:

    John Burkardt
*/
{
  int seed = 123456789;
  int test;
  int test_num = 10;
  float x;
  float y1;
  double y2;
  float z;

  printf ( "\n" );
  printf ( "TEST_MODF:\n" );
  printf ( "  Test MODF, which splits X into integer and fractional\n" );
  printf ( "  parts.\n" );
  printf ( "\n" );
  printf ( "        X           Y1          Y2         Z\n" );
  printf ( "                                         Y1+Y2\n" );
  printf ( "\n" );

  for ( test = 1; test <= test_num; test++ )
  {
    x = r4_uniform ( -20.0, +20.0, &seed );
    y1 = modf ( x, &y2 );
    z = y1 + y2;
    printf ( "  %10.4f  %10.4f  %10.4f  %10.4f\n", x, y1, y2, z );
  }
  return;
}
/******************************************************************************/

void test_pow ( void )

/******************************************************************************/
/*
  Purpose:

    TEST_POW tests POW.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 April 2011

  Author:

    John Burkardt
*/
{
  int seed = 123456789;
  int test;
  int test_num = 10;
  float x1;
  float x2;
  float y;
  float z;

  printf ( "\n" );
  printf ( "TEST_POW:\n" );
  printf ( "  Test POW, which evaluates the power function X1^X2.\n" );
  printf ( "\n" );
  printf ( "         X1        X2          Y\n" );
  printf ( "                           POW(X1,X2)\n" );
  printf ( "\n" );

  for ( test = 1; test <= test_num; test++ )
  {
    x1 = r4_uniform ( 0.0, +10.0, &seed );
    x2 = r4_uniform ( -2.0, +10.0, &seed );
    y = pow ( x1, x2 );
    printf ( "  %10.4f  %10.4f  %10.4e\n", x1, x2, y );
  }
  return;
}
/******************************************************************************/

void test_sin ( void )

/******************************************************************************/
/*
  Purpose:

    TEST_SIN tests SIN.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 April 2011

  Author:

    John Burkardt
*/
{
  float pi = 3.141592653589793;
  int seed = 123456789;
  int test;
  int test_num = 10;
  float x;
  float y;

  printf ( "\n" );
  printf ( "TEST_SIN:\n" );
  printf ( "  Test SIN, which evaluates the sine function.\n" );
  printf ( "\n" );
  printf ( "         X          Y\n" );
  printf ( "                 SIN(X)\n" );
  printf ( "\n" );

  for ( test = 1; test <= test_num; test++ )
  {
    x = r4_uniform ( -pi, +pi, &seed );
    y = sin ( x );
    printf ( "  %10.4f  %10.4f\n", x, y );
  }
  return;
}
/******************************************************************************/

void test_sinh ( void )

/******************************************************************************/
/*
  Purpose:

    TEST_SINH tests SINH.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 April 2011

  Author:

    John Burkardt
*/
{
  float pi = 3.141592653589793;
  int seed = 123456789;
  int test;
  int test_num = 10;
  float x;
  float y;

  printf ( "\n" );
  printf ( "TEST_SINH:\n" );
  printf ( "  Test SINH, which evaluates the hyperbolic sine function.\n" );
  printf ( "\n" );
  printf ( "         X          Y\n" );
  printf ( "                 SINH(X)\n" );
  printf ( "\n" );

  for ( test = 1; test <= test_num; test++ )
  {
    x = r4_uniform ( -5.0, +5.0, &seed );
    y = sinh ( x );
    printf ( "  %10.4f  %10.4f\n", x, y );
  }
  return;
}
/******************************************************************************/

void test_sqrt ( void )

/******************************************************************************/
/*
  Purpose:

    TEST_SQRT tests SQRT.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 April 2011

  Author:

    John Burkardt
*/
{
  int seed = 123456789;
  int test;
  int test_num = 10;
  float x;
  float y;
  float z;

  printf ( "\n" );
  printf ( "TEST_SQRT:\n" );
  printf ( "  Test SQRT, which evaluates the square root function.\n" );
  printf ( "\n" );
  printf ( "         X          Y           Z\n" );
  printf ( "                 SQRT(X)       Y*Y\n" );
  printf ( "\n" );

  for ( test = 1; test <= test_num; test++ )
  {
    x = r4_uniform ( 0.0, +100.0, &seed );
    y = sqrt ( x );
    z = y * y;
    printf ( "  %10.4f  %10.4f  %10.4f\n", x, y, z );
  }
  return;
}
/******************************************************************************/

void test_tan ( void )

/******************************************************************************/
/*
  Purpose:

    TEST_TAN tests TAN.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 April 2011

  Author:

    John Burkardt
*/
{
  float pi = 3.141592653589793;
  int seed = 123456789;
  int test;
  int test_num = 10;
  float x;
  float y;

  printf ( "\n" );
  printf ( "TEST_TAN:\n" );
  printf ( "  Test TAN, which evaluates the tangent function.\n" );
  printf ( "\n" );
  printf ( "         X          Y\n" );
  printf ( "                 TAN(X)\n" );
  printf ( "\n" );

  for ( test = 1; test <= test_num; test++ )
  {
    x = r4_uniform ( -pi/2.0, +pi/2.0, &seed );
    y = tan ( x );
    printf ( "  %10.4f  %10.4f\n", x, y );
  }
  return;
}
/******************************************************************************/

void test_tanh ( void )

/******************************************************************************/
/*
  Purpose:

    TEST_TANH tests TANH.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 April 2011

  Author:

    John Burkardt
*/
{
  float pi = 3.141592653589793;
  int seed = 123456789;
  int test;
  int test_num = 10;
  float x;
  float y;

  printf ( "\n" );
  printf ( "TEST_TANH:\n" );
  printf ( "  Test TANH, which evaluates the hyperbolic tangent function.\n" );
  printf ( "\n" );
  printf ( "         X          Y\n" );
  printf ( "                 TANH(X)\n" );
  printf ( "\n" );

  for ( test = 1; test <= test_num; test++ )
  {
    x = r4_uniform ( -5.0, +5.0, &seed );
    y = tanh ( x );
    printf ( "  %10.4f  %10.4f\n", x, y );
  }
  return;
}
/******************************************************************************/

int i4_max ( int i1, int i2 )

/******************************************************************************/
/*
  Purpose:

    I4_MAX returns the maximum of two I4's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    29 August 2006

  Author:

    John Burkardt

  Parameters:

    Input, int I1, I2, are two integers to be compared.

    Output, int I4_MAX, the larger of I1 and I2.
*/
{
  int value;

  if ( i2 < i1 )
  {
    value = i1;
  }
  else
  {
    value = i2;
  }
  return value;
}
/******************************************************************************/

int i4_min ( int i1, int i2 )

/******************************************************************************/
/*
  Purpose:

    I4_MIN returns the smaller of two I4's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    29 August 2006

  Author:

    John Burkardt

  Parameters:

    Input, int I1, I2, two integers to be compared.

    Output, int I4_MIN, the smaller of I1 and I2.
*/
{
  int value;

  if ( i1 < i2 )
  {
    value = i1;
  }
  else
  {
    value = i2;
  }
  return value;
}
/******************************************************************************/

int i4_uniform ( int a, int b, int *seed )

/******************************************************************************/
/*
  Purpose:

    I4_UNIFORM returns a scaled pseudorandom I4.

  Discussion:

    The pseudorandom number should be uniformly distributed
    between A and B.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 November 2006

  Author:

    John Burkardt

  Reference:

    Paul Bratley, Bennett Fox, Linus Schrage,
    A Guide to Simulation,
    Springer Verlag, pages 201-202, 1983.

    Pierre L'Ecuyer,
    Random Number Generation,
    in Handbook of Simulation,
    edited by Jerry Banks,
    Wiley Interscience, page 95, 1998.

    Bennett Fox,
    Algorithm 647:
    Implementation and Relative Efficiency of Quasirandom
    Sequence Generators,
    ACM Transactions on Mathematical Software,
    Volume 12, Number 4, pages 362-376, 1986.

    Peter Lewis, Allen Goodman, James Miller
    A Pseudo-Random Number Generator for the System/360,
    IBM Systems Journal,
    Volume 8, pages 136-143, 1969.

  Parameters:

    Input, int A, B, the limits of the interval.

    Input/output, int *SEED, the "seed" value, which should NOT be 0.
    On output, SEED has been updated.

    Output, int I4_UNIFORM, a number between A and B.
*/
{
  int k;
  float r;
  int value;

  if ( *seed == 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "I4_UNIFORM - Fatal error!\n" );
    fprintf ( stderr, "  Input value of SEED = 0.\n" );
    exit ( 1 );
  }

  k = *seed / 127773;

  *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

  if ( *seed < 0 )
  {
    *seed = *seed + 2147483647;
  }

  r = ( float ) ( *seed ) * 4.656612875E-10;
/*
  Scale R to lie between A-0.5 and B+0.5.
*/
  r = ( 1.0 - r ) * ( ( float ) ( i4_min ( a, b ) ) - 0.5 )
    +         r   * ( ( float ) ( i4_max ( a, b ) ) + 0.5 );
/*
  Use rounding to convert R to an integer between A and B.
*/
  value = r4_nint ( r );

  value = i4_max ( value, i4_min ( a, b ) );
  value = i4_min ( value, i4_max ( a, b ) );

  return value;
}
/******************************************************************************/

int r4_nint ( float x )

/******************************************************************************/
/*
  Purpose:

    R4_NINT returns the nearest integer to an R4.

  Example:

        X         R4_NINT

      1.3         1
      1.4         1
      1.5         1 or 2
      1.6         2
      0.0         0
     -0.7        -1
     -1.1        -1
     -1.6        -2

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 May 2006

  Author:

    John Burkardt

  Parameters:

    Input, float X, the value.

    Output, int R4_NINT, the nearest integer to X.
*/
{
  int s;
  int value;

  if ( x < 0.0 )
  {
    s = - 1;
  }
  else
  {
    s = + 1;
  }
  value = s * ( int ) ( fabs ( x ) + 0.5 );

  return value;
}
/******************************************************************************/

float r4_uniform ( float b, float c, int *seed )

/******************************************************************************/
/*
  Purpose:

    R4_UNIFORM returns a scaled pseudorandom R4.

  Discussion:

    The pseudorandom number should be uniformly distributed
    between A and B.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    21 November 2004

  Author:

    John Burkardt

  Reference:

    Paul Bratley, Bennett Fox, Linus Schrage,
    A Guide to Simulation,
    Second Edition,
    Springer, 1987,
    ISBN: 0387964673,
    LC: QA76.9.C65.B73.

    Bennett Fox,
    Algorithm 647:
    Implementation and Relative Efficiency of Quasirandom
    Sequence Generators,
    ACM Transactions on Mathematical Software,
    Volume 12, Number 4, December 1986, pages 362-376.

    Pierre L'Ecuyer,
    Random Number Generation,
    in Handbook of Simulation,
    edited by Jerry Banks,
    Wiley, 1998,
    ISBN: 0471134031,
    LC: T57.62.H37.

    Peter Lewis, Allen Goodman, James Miller,
    A Pseudo-Random Number Generator for the System/360,
    IBM Systems Journal,
    Volume 8, Number 2, 1969, pages 136-143.

  Parameters:

    Input, float B, C, the limits of the interval.

    Input/output, int *SEED, the "seed" value, which should NOT be 0.
    On output, SEED has been updated.

    Output, float R4_UNIFORM, a number strictly between A and B.
*/
{
  float value;

  if ( *seed == 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R4_UNIFORM - Fatal error!\n" );
    fprintf ( stderr, "  Input value of SEED = 0.\n" );
    exit ( 1 );
  }

  value = b + ( c - b ) * r4_uniform_01 ( seed );

  return value;
}
/******************************************************************************/

float r4_uniform_01 ( int *seed )

/******************************************************************************/
/*
  Purpose:

    R4_UNIFORM_01 returns a real pseudorandom R4.

  Discussion:

    This routine implements the recursion

      seed = 16807 * seed mod ( 2**31 - 1 )
      r4_uniform_01 = seed / ( 2**31 - 1 )

    The integer arithmetic never requires more than 32 bits,
    including a sign bit.

    If the initial seed is 12345, then the first three computations are

      Input     Output      R4_UNIFORM_01
      SEED      SEED

         12345   207482415  0.096616
     207482415  1790989824  0.833995
    1790989824  2035175616  0.947702

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 November 2004

  Author:

    John Burkardt

  Reference:

    Paul Bratley, Bennett Fox, Linus Schrage,
    A Guide to Simulation,
    Springer Verlag, pages 201-202, 1983.

    Pierre L'Ecuyer,
    Random Number Generation,
    in Handbook of Simulation
    edited by Jerry Banks,
    Wiley Interscience, page 95, 1998.

    Bennett Fox,
    Algorithm 647:
    Implementation and Relative Efficiency of Quasirandom
    Sequence Generators,
    ACM Transactions on Mathematical Software,
    Volume 12, Number 4, pages 362-376, 1986.

    Peter Lewis, Allen Goodman, James Miller,
    A Pseudo-Random Number Generator for the System/360,
    IBM Systems Journal,
    Volume 8, pages 136-143, 1969.

  Parameters:

    Input/output, int *SEED, the "seed" value.  Normally, this
    value should not be 0.  On output, SEED has been updated.

    Output, float R4_UNIFORM_01, a new pseudorandom variate, strictly between
    0 and 1.
*/
{
  int k;
  float r;

  k = *seed / 127773;

  *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

  if ( *seed < 0 )
  {
    *seed = *seed + 2147483647;
  }
/*
  Although SEED can be represented exactly as a 32 bit integer,
  it generally cannot be represented exactly as a 32 bit real number!
*/
  r = ( float ) ( *seed ) * 4.656612875E-10;

  return r;
}
/******************************************************************************/

void timestamp ( void )

/******************************************************************************/
/*
  Purpose:

    TIMESTAMP prints the current YMDHMS date as a time stamp.

  Example:

    31 May 2001 09:45:54 AM

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    24 September 2003

  Author:

    John Burkardt

  Parameters:

    None
*/
{
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct tm *tm;
  size_t len;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  len = strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  printf ( "%s\n", time_buffer );

  return;
# undef TIME_SIZE
}
