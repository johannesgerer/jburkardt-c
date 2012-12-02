# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <complex.h>
# include <time.h>

# include "c8lib.h"

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
void test25 ( void );
void test26 ( void );
void test27 ( void );
void test28 ( void );
void test102 ( void );
void test103 ( void );
void test104 ( void );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    C8LIB_PRB calls the C8LIB tests.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2010

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "C8LIB_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the C8LIB library.\n" );

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
  test25 ( );
  test26 ( );
  test27 ( );
  test28 ( );

  test102 ( );
  test103 ( );
  test104 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "C8LIB_PRB\n" );
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

    TEST01 tests C8_ABS.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2010

  Author:

    John Burkardt
*/
{
  double complex c1;
  int i;
  double r2;
  int seed;

  seed = 123456789;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  C8_ABS computes the absolute value of a C8.\n" );
  printf ( "\n" );
  printf ( "       C1=C8_UNIFORM_01          R2=C8_ABS(C1)\n" );
  printf ( "     ---------------------     ---------------------\n" );
  printf ( "\n" );
  for ( i = 0; i < 10; i++ )
  {
    c1 = c8_uniform_01 ( &seed );
    r2 = c8_abs ( c1 );

    printf ( "  %12.6f%12.6f  %12.6f\n", 
      creal ( c1 ), cimag ( c1 ), r2 );
  }

  return;
}
/******************************************************************************/

void test02 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST02 tests C8_ACOS and C8_COS.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2010

  Author:

    John Burkardt
*/
{
  double complex c1;
  double complex c2;
  double complex c3;
  int i;
  int seed;

  seed = 123456789;

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  C8_ACOS computes the inverse cosine;\n" );
  printf ( "  C8_COS computes the cosine;\n" );
  printf ( "\n" );
  printf ( "       C1=C8_UNIFORM_01          C2 = C8_ACOS(C1)           C3 = C8_COS(C2)\n" );
  printf ( "     ---------------------     ---------------------     ---------------------\n" );
  printf ( "\n" );

  for ( i = 0; i < 10; i++ )
  {
    c1 = c8_uniform_01 ( &seed );
    c2 = c8_acos ( c1 );
    c3 = c8_cos ( c2 );
    printf ( "  %12.6f%12.6f  %12.6f%12.6f  %12.6f%12.6f\n", 
      creal ( c1 ), cimag ( c1 ), creal ( c2 ), cimag ( c2 ), 
      creal ( c3 ), cimag ( c3 ) );
  }

  return;
}
/******************************************************************************/

void test03 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST03 tests C8_ACOSH and C8_COSH.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2010

  Author:

    John Burkardt
*/
{
  double complex c1;
  double complex c2;
  double complex c3;
  int i;
  int seed;

  seed = 123456789;

  printf ( "\n" );
  printf ( "TEST03\n" );
  printf ( "  C8_ACOSH computes the inverse hyperbolic cosine;\n" );
  printf ( "  C8_COSH computes the hyperbolic cosine;\n" );
  printf ( "\n" );
  printf ( "       C1=C8_UNIFORM_01          C2 = C8_ACOSH(C1)          C3 = C8_COSH(C2)\n" );
  printf ( "     ---------------------     ---------------------     ---------------------\n" );
  printf ( "\n" );

  for ( i = 0; i < 10; i++ )
  {
    c1 = c8_uniform_01 ( &seed );
    c2 = c8_acosh ( c1 );
    c3 = c8_cosh ( c2 );

    printf ( "  %12.6f%12.6f  %12.6f%12.6f  %12.6f%12.6f\n", 
      creal ( c1 ), cimag ( c1 ), creal ( c2 ), cimag ( c2 ), 
      creal ( c3 ), cimag ( c3 ) );
  }

  return;
}
/******************************************************************************/

void test04 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST04 tests C8_ADD.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2010

  Author:

    John Burkardt
*/
{
  double complex c1;
  double complex c2;
  double complex c3;
  int i;
  int seed;

  seed = 123456789;

  printf ( "\n" );
  printf ( "TEST04\n" );
  printf ( "  C8_ADD computes C3 = C1 + C2.\n" );
  printf ( "\n" );
  printf ( "       C1=C8_UNIFORM_01          C2=C8_UNIFORM_01          C3 = C8_ADD(C1,C2)\n" );
  printf ( "     ---------------------     ---------------------     ---------------------\n" );
  printf ( "\n" );

  for ( i = 0; i < 10; i++ )
  {
    c1 = c8_uniform_01 ( &seed );
    c2 = c8_uniform_01 ( &seed );
    c3 = c8_add ( c1, c2 );

    printf ( "  %12.6f%12.6f  %12.6f%12.6f  %12.6f%12.6f\n", 
      creal ( c1 ), cimag ( c1 ), creal ( c2 ), cimag ( c2 ), 
      creal ( c3 ), cimag ( c3 ) );
  }

  return;
}
/******************************************************************************/

void test05 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST05 tests C8_ARG.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2010

  Author:

    John Burkardt
*/
{
  double complex c1;
  double r2;
  int i;
  int seed;

  seed = 123456789;

  printf ( "\n" );
  printf ( "TEST05\n" );
  printf ( "  C8_ARG computes the argument of a C8 value.\n" );
  printf ( "\n" );
  printf ( "       C1=C8_UNIFORM_01          R2=C8_ARG(C1)\n" );
  printf ( "     ---------------------     ---------------------\n" );
  printf ( "\n" );
  for ( i = 0; i < 10; i++ )
  {
    c1 = c8_uniform_01 ( &seed );
    r2 = c8_arg ( c1 );

    printf ( "  %12.6f%12.6f  %12.6f\n", 
      creal ( c1 ), cimag ( c1 ), r2 );
  }

  return;
}
/******************************************************************************/

void test06 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST06 tests C8_ASIN and C8_SIN.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2010

  Author:

    John Burkardt
*/
{
  double complex c1;
  double complex c2;
  double complex c3;
  int i;
  int seed;

  seed = 123456789;

  printf ( "\n" );
  printf ( "TEST06\n" );
  printf ( "  C8_ASIN computes the inverse sine;\n" );
  printf ( "  C8_SIN computes the sine;\n" );
  printf ( "\n" );
  printf ( "       C1=C8_UNIFORM_01          C2 = C8_ASIN(C1)           C3 = C8_SIN(C2)\n" );
  printf ( "     ---------------------     ---------------------     ---------------------\n" );
  printf ( "\n" );

  for ( i = 0; i < 10; i++ )
  {
    c1 = c8_uniform_01 ( &seed );
    c2 = c8_asin ( c1 );
    c3 = c8_sin ( c2 );
    printf ( "  %12.6f%12.6f  %12.6f%12.6f  %12.6f%12.6f\n", 
      creal ( c1 ), cimag ( c1 ), creal ( c2 ), cimag ( c2 ), 
      creal ( c3 ), cimag ( c3 ) );
  }

  return;
}
/******************************************************************************/

void test07 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST07 tests C8_ASINH and C8_SINH.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2010

  Author:

    John Burkardt
*/
{
  double complex c1;
  double complex c2;
  double complex c3;
  int i;
  int seed;

  seed = 123456789;

  printf ( "\n" );
  printf ( "TEST07\n" );
  printf ( "  C8_ASINH computes the inverse hyperbolic sine;\n" );
  printf ( "  C8_SINH computes the hyperbolic sine;\n" );
  printf ( "\n" );
  printf ( "       C1=C8_UNIFORM_01          C2 = C8_ASINH(C1)          C3 = C8_SINH(C2)\n" );
  printf ( "     ---------------------     ---------------------     ---------------------\n" );
  printf ( "\n" );

  for ( i = 0; i < 10; i++ )
  {
    c1 = c8_uniform_01 ( &seed );
    c2 = c8_asinh ( c1 );
    c3 = c8_sinh ( c2 );
    printf ( "  %12.6f%12.6f  %12.6f%12.6f  %12.6f%12.6f\n", 
      creal ( c1 ), cimag ( c1 ), creal ( c2 ), cimag ( c2 ), 
      creal ( c3 ), cimag ( c3 ) );
  }

  return;
}
/******************************************************************************/

void test08 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST08 tests C8_ATAN and C8_TAN.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2010

  Author:

    John Burkardt
*/
{
  double complex c1;
  double complex c2;
  double complex c3;
  int i;
  int seed;

  seed = 123456789;

  printf ( "\n" );
  printf ( "TEST08\n" );
  printf ( "  C8_ATAN computes the inverse tangent;\n" );
  printf ( "  C8_TAN computes the tangent;\n" );
  printf ( "\n" );
  printf ( "       C1=C8_UNIFORM_01          C2 = C8_ATAN(C1)          C3 = C8_TAN(C2)\n" );
  printf ( "     ---------------------     ---------------------     ---------------------\n" );
  printf ( "\n" );

  for ( i = 0; i < 10; i++ )
  {
    c1 = c8_uniform_01 ( &seed );
    c2 = c8_atan ( c1 );
    c3 = c8_tan ( c2 );
    printf ( "  %12.6f%12.6f  %12.6f%12.6f  %12.6f%12.6f\n", 
      creal ( c1 ), cimag ( c1 ), creal ( c2 ), cimag ( c2 ), 
      creal ( c3 ), cimag ( c3 ) );
  }

  return;
}
/******************************************************************************/

void test09 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST09 tests C8_ATANH and C8_TANH.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2010

  Author:

    John Burkardt
*/
{
  double complex c1;
  double complex c2;
  double complex c3;
  int i;
  int seed;

  seed = 123456789;

  printf ( "\n" );
  printf ( "TEST09\n" );
  printf ( "  C8_ATANH computes the inverse hyperbolic tangent;\n" );
  printf ( "  C8_TANH computes the hyperbolic tangent;\n" );
  printf ( "\n" );
  printf ( "       C1=C8_UNIFORM_01          C2 = C8_ATANH(C1)         C3 = C8_TANH(C2)\n" );
  printf ( "     ---------------------     ---------------------     ---------------------\n" );
  printf ( "\n" );

  for ( i = 0; i < 10; i++ )
  {
    c1 = c8_uniform_01 ( &seed );
    c2 = c8_atanh ( c1 );
    c3 = c8_tanh ( c2 );
    printf ( "  %12.6f%12.6f  %12.6f%12.6f  %12.6f%12.6f\n", 
      creal ( c1 ), cimag ( c1 ), creal ( c2 ), cimag ( c2 ), 
      creal ( c3 ), cimag ( c3 ) );
  }

  return;
}
/******************************************************************************/

void test10 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST10 tests C8_CONJ.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2010

  Author:

    John Burkardt
*/
{
  double complex c1;
  double complex c2;
  int i;
  int seed;

  seed = 123456789;

  printf ( "\n" );
  printf ( "TEST10\n" );
  printf ( "  C8_CONJ computes C2 = conj ( C1 ).\n" );
  printf ( "\n" );
  printf ( "       C1=C8_UNIFORM_01          C2=C8_CONJ(C1)\n" );
  printf ( "     ---------------------     ---------------------\n" );
  printf ( "\n" );

  for ( i = 0; i < 10; i++ )
  {
    c1 = c8_uniform_01 ( &seed );
    c2 = c8_conj ( c1 );

    printf ( "  %12.6f%12.6f  %12.6f%12.6f\n", 
      creal ( c1 ), cimag ( c1 ), creal ( c2 ), cimag ( c2 ) );
  }

  return;
}
/******************************************************************************/

void test11 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST11 tests C8_COS.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2010

  Author:

    John Burkardt
*/
{
  double complex c1;
  double complex c2;
  int i;
  int seed;

  seed = 123456789;

  printf ( "\n" );
  printf ( "TEST11\n" );
  printf ( "  C8_COS computes the cosine of a C8;\n" );
  printf ( "\n" );
  printf ( "       C1=C8_UNIFORM_01          C2 = C8_COS(C1)\n" );
  printf ( "     ---------------------     ---------------------\n" );
  printf ( "\n" );

  for ( i = 0; i < 10; i++ )
  {
    c1 = c8_uniform_01 ( &seed );
    c2 = c8_cos ( c1 );

    printf ( "  %12.6f%12.6f  %12.6f%12.6f\n", 
      creal ( c1 ), cimag ( c1 ), creal ( c2 ), cimag ( c2 ) );
  }

  return;
}
/******************************************************************************/

void test12 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST12 tests C8_COSH.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2010

  Author:

    John Burkardt
*/
{
  double complex c1;
  double complex c2;
  int i;
  int seed;

  seed = 123456789;

  printf ( "\n" );
  printf ( "TEST12\n" );
  printf ( "  C8_COSH computes the hyperbolic cosine of a C8;\n" );
  printf ( "\n" );
  printf ( "       C1=C8_UNIFORM_01          C2 = C8_COSH(C1)\n" );
  printf ( "     ---------------------     ---------------------\n" );
  printf ( "\n" );

  for ( i = 0; i < 10; i++ )
  {
    c1 = c8_uniform_01 ( &seed );
    c2 = c8_cosh ( c1 );

    printf ( "  %12.6f%12.6f  %12.6f%12.6f\n", 
      creal ( c1 ), cimag ( c1 ), creal ( c2 ), cimag ( c2 ) );
  }

  return;
}
/******************************************************************************/

void test13 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST13 tests C8_CUBE_ROOT.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2010

  Author:

    John Burkardt
*/
{
  double complex c1;
  double complex c2;
  double complex c3;
  double complex c4;
  int i;
  int seed;

  seed = 123456789;

  printf ( "\n" );
  printf ( "TEST13\n" );
  printf ( "  C8_CUBE_ROOT computes C2 = cube root ( C1 ).\n" );
  printf ( "  Check by C3 = C2 * C2.\n" );
  printf ( "\n" );
  printf ( "        C1=C8_UNIFORM_01          C2=C8_CUBE_ROOT(C1)      C3=C2*C2*C2\n" );
  printf ( "     ---------------------     ---------------------     ---------------------\n" );
  printf ( "\n" );

  for ( i = 0; i < 10; i++ )
  {
    c1 = c8_uniform_01 ( &seed );
    c2 = c8_cube_root ( c1 );
    c3 = c8_mul ( c2, c2 );
    c4 = c8_mul ( c3, c2 );

    printf ( "  %12.6f%12.6f  %12.6f%12.6f  %12.6f%12.6f\n", 
      creal ( c1 ), cimag ( c1 ), creal ( c2 ), cimag ( c2 ), 
      creal ( c4 ), cimag ( c4 ) );
  }

  return;
}
/******************************************************************************/

void test14 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST14 tests C8_DIV.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2010

  Author:

    John Burkardt
*/
{
  double complex c1;
  double complex c2;
  double complex c3;
  int i;
  int seed;

  seed = 123456789;

  printf ( "\n" );
  printf ( "TEST14\n" );
  printf ( "  C8_DIV computes C3 = C1 / C2.\n" );
  printf ( "\n" );
  printf ( "        C1=C8_UNIFORM_01          C2=C8_UNIFORM_01          C3=C8_DIV(C1,C2)\n" );
  printf ( "     ---------------------     ---------------------     ---------------------\n" );
  printf ( "\n" );

  for ( i = 0; i < 10; i++ )
  {
    c1 = c8_uniform_01 ( &seed );
    c2 = c8_uniform_01 ( &seed );
    c3 = c8_div ( c1, c2 );

    printf ( "  %12.6f%12.6f  %12.6f%12.6f  %12.6f%12.6f\n", 
      creal ( c1 ), cimag ( c1 ), creal ( c2 ), cimag ( c2 ), 
      creal ( c3 ), cimag ( c3 ) );
  }

  return;
}
/******************************************************************************/

void test15 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST15 tests C8_EXP.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2010

  Author:

    John Burkardt
*/
{
  double complex c1;
  double complex c2;
  double complex c3;
  int i;
  int seed;

  seed = 123456789;

  printf ( "\n" );
  printf ( "TEST15\n" );
  printf ( "  C8_EXP computes C2 = e ^ C1.\n" );
  printf ( "\n" );
  printf ( "        C1=C8_UNIFORM_01          C2=C8_EXP(C1)\n" );
  printf ( "     ---------------------     ---------------------\n" );
  printf ( "\n" );

  for ( i = 0; i < 10; i++ )
  {
    c1 = c8_uniform_01 ( &seed );
    c2 = c8_exp ( c1 );

    printf ( "  %12.6f%12.6f  %12.6f%12.6f\n", 
      creal ( c1 ), cimag ( c1 ), creal ( c2 ), cimag ( c2 ) );
  }
  return;
}
/******************************************************************************/

void test16 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST16 tests C8_INV.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2010

  Author:

    John Burkardt
*/
{
  double complex c1;
  double complex c2;
  double complex c3;
  int i;
  int seed;

  seed = 123456789;

  printf ( "\n" );
  printf ( "TEST16\n" );
  printf ( "  C8_INV computes C2 = 1 / C1.\n" );
  printf ( "  Check by C3 = 1 / C2.\n" );
  printf ( "\n" );
  printf ( "        C1=C8_UNIFORM_01          C2=C8_INV(C1)             C3=C8_INV(C2)\n" );
  printf ( "     ---------------------     ---------------------     ---------------------\n" );
  printf ( "\n" );

  for ( i = 0; i < 10; i++ )
  {
    c1 = c8_uniform_01 ( &seed );
    c2 = c8_inv ( c1 );
    c3 = c8_inv ( c2 );

    printf ( "  %12.6f%12.6f  %12.6f%12.6f  %12.6f%12.6f\n", 
      creal ( c1 ), cimag ( c1 ), creal ( c2 ), cimag ( c2 ), 
      creal ( c3 ), cimag ( c3 ) );
  }

  return;
}
/******************************************************************************/

void test17 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST17 tests C8_LOG.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2010

  Author:

    John Burkardt
*/
{
  double complex c1;
  double complex c2;
  double complex c3;
  int i;
  int seed;

  seed = 123456789;

  printf ( "\n" );
  printf ( "TEST17\n" );
  printf ( "  C8_LOG computes log ( Z ).\n" );
  printf ( "\n" );
  printf ( "        C1=C8_UNIFORM_01          C2=C8_LOG(C1)             C3=C8_EXP(C2)\n" );
  printf ( "     ---------------------     ---------------------     ---------------------\n" );
  printf ( "\n" );
  for ( i = 0; i < 10; i++ )
  {
    c1 = c8_uniform_01 ( &seed );
    c2 = c8_log ( c1 );
    c3 = c8_exp ( c2 );

    printf ( "  %12.6f%12.6f  %12.6f%12.6f  %12.6f%12.6f\n", 
      creal ( c1 ), cimag ( c1 ), creal ( c2 ), cimag ( c2 ), 
      creal ( c3 ), cimag ( c3 ) );
  }

  return;
}
/******************************************************************************/

void test18 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST18 tests C8_MAG.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2010

  Author:

    John Burkardt
*/
{
  double complex c1;
  double r2;
  int i;
  int seed;

  seed = 123456789;

  printf ( "\n" );
  printf ( "TEST18\n" );
  printf ( "  C8_MAG computes the magnitude of a C8.\n" );
  printf ( "\n" );
  printf ( "       C1=C8_UNIFORM_01          R2=C8_MAG(C1)\n" );
  printf ( "     ---------------------     ---------------------\n" );
  printf ( "\n" );
  for ( i = 0; i < 10; i++ )
  {
    c1 = c8_uniform_01 ( &seed );
    r2 = c8_mag ( c1 );

    printf ( "  %12.6f%12.6f  %12.6f\n", 
      creal ( c1 ), cimag ( c1 ), r2 );
  }

  return;
}
/******************************************************************************/

void test19 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST19 tests C8_MUL.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2010

  Author:

    John Burkardt
*/
{
  double complex c1;
  double complex c2;
  double complex c3;
  int i;
  int seed;

  seed = 123456789;

  printf ( "\n" );
  printf ( "TEST19\n" );
  printf ( "  C8_MUL computes C3 = C1 * C2.\n" );
  printf ( "\n" );
  printf ( "        C1=C8_UNIFORM_01          C2=C8_UNIFORM_01          C3=C8_MUL(C1,C2)\n" );
  printf ( "     ---------------------     ---------------------     ---------------------\n" );
  printf ( "\n" );

  for ( i = 0; i < 10; i++ )
  {
    c1 = c8_uniform_01 ( &seed );
    c2 = c8_uniform_01 ( &seed );
    c3 = c8_mul ( c1, c2 );

    printf ( "  %12.6f%12.6f  %12.6f%12.6f  %12.6f%12.6f\n", 
      creal ( c1 ), cimag ( c1 ), creal ( c2 ), cimag ( c2 ), 
      creal ( c3 ), cimag ( c3 ) );
  }

  return;
}
/******************************************************************************/

void test20 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST20 tests C8_NORMAL_01.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2010

  Author:

    John Burkardt
*/
{
  double complex c;
  int i;
  int seed;

  seed = 123456789;

  printf ( "\n" );
  printf ( "TEST20\n" );
  printf ( "  C8_NORMAL_01 generates unit pseudonormal C8's\n" );
  printf ( "\n" );
  printf ( "       C1=C8_NORMAL_01(SEED)\n" );
  printf ( "     ---------------------\n" );
  printf ( "\n" );

  for ( i = 0; i < 10; i++ )
  {
    c = c8_normal_01 ( &seed );

    printf ( "  %12.6f%12.6f\n", creal ( c ), cimag ( c ) );
  }

  return;
}
/******************************************************************************/

void test21 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST21 tests C8_SIN.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2010

  Author:

    John Burkardt
*/
{
  double complex c1;
  double complex c2;
  int i;
  int seed;

  seed = 123456789;

  printf ( "\n" );
  printf ( "TEST21\n" );
  printf ( "  C8_SIN computes the sine of a C8;\n" );
  printf ( "\n" );
  printf ( "       C1=C8_UNIFORM_01          C2 = C8_SIN(C1)\n" );
  printf ( "     ---------------------     ---------------------\n" );
  printf ( "\n" );

  for ( i = 0; i < 10; i++ )
  {
    c1 = c8_uniform_01 ( &seed );
    c2 = c8_sin ( c1 );

    printf ( "  %12.6f%12.6f  %12.6f%12.6f\n", 
      creal ( c1 ), cimag ( c1 ), creal ( c2 ), cimag ( c2 ) );
  }

  return;
}
/******************************************************************************/

void test22 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST22 tests C8_SINH.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2010

  Author:

    John Burkardt
*/
{
  double complex c1;
  double complex c2;
  int i;
  int seed;

  seed = 123456789;

  printf ( "\n" );
  printf ( "TEST22\n" );
  printf ( "  C8_SINH computes the hyperbolic sine of a C8;\n" );
  printf ( "\n" );
  printf ( "       C1=C8_UNIFORM_01          C2 = C8_SINH(C1)\n" );
  printf ( "     ---------------------     ---------------------\n" );
  printf ( "\n" );

  for ( i = 0; i < 10; i++ )
  {
    c1 = c8_uniform_01 ( &seed );
    c2 = c8_sinh ( c1 );

    printf ( "  %12.6f%12.6f  %12.6f%12.6f\n", 
      creal ( c1 ), cimag ( c1 ), creal ( c2 ), cimag ( c2 ) );
  }

  return;
}
/******************************************************************************/

void test23 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST23 tests C8_SQRT.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2010

  Author:

    John Burkardt
*/
{
  double complex c1;
  double complex c2;
  double complex c3;
  int i;
  int seed;

  seed = 123456789;

  printf ( "\n" );
  printf ( "TEST23\n" );
  printf ( "  C8_SQRT computes C2 = sqrt ( C1 ).\n" );
  printf ( "  Check by C3 = C2 * C2.\n" );
  printf ( "\n" );
  printf ( "        C1=C8_UNIFORM_01          C2=C8_SQRT(C1)            C3=C2*C2\n" );
  printf ( "     ---------------------     ---------------------     ---------------------\n" );
  printf ( "\n" );

  for ( i = 0; i < 10; i++ )
  {
    c1 = c8_uniform_01 ( &seed );
    c2 = c8_sqrt ( c1 );
    c3 = c8_mul ( c2, c2 );

    printf ( "  %12.6f%12.6f  %12.6f%12.6f  %12.6f%12.6f\n", 
      creal ( c1 ), cimag ( c1 ), creal ( c2 ), cimag ( c2 ), 
      creal ( c3 ), cimag ( c3 ) );
  }

  return;
}
/******************************************************************************/

void test24 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST24 tests C8_SUB.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2010

  Author:

    John Burkardt
*/
{
  double complex c1;
  double complex c2;
  double complex c3;
  int i;
  int seed;

  seed = 123456789;

  printf ( "\n" );
  printf ( "TEST24\n" );
  printf ( "  C8_SUB computes C3 = C1 - C2.\n" );
  printf ( "\n" );
  printf ( "       C1=C8_UNIFORM_01          C2=C8_UNIFORM_01          C3 = C8_SUB(C1,C2)\n" );
  printf ( "     ---------------------     ---------------------     ---------------------\n" );
  printf ( "\n" );

  for ( i = 0; i < 10; i++ )
  {
    c1 = c8_uniform_01 ( &seed );
    c2 = c8_uniform_01 ( &seed );
    c3 = c8_sub ( c1, c2 );

    printf ( "  %12.6f%12.6f  %12.6f%12.6f  %12.6f%12.6f\n", 
      creal ( c1 ), cimag ( c1 ), creal ( c2 ), cimag ( c2 ), 
      creal ( c3 ), cimag ( c3 ) );
  }

  return;
}
/******************************************************************************/

void test25 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST25 tests C8_TAN.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2010

  Author:

    John Burkardt
*/
{
  double complex c1;
  double complex c2;
  int i;
  int seed;

  seed = 123456789;

  printf ( "\n" );
  printf ( "TEST25\n" );
  printf ( "  C8_TAN computes the tangent of a C8;\n" );
  printf ( "\n" );
  printf ( "       C1=C8_UNIFORM_01          C2 = C8_TAN(C1)\n" );
  printf ( "     ---------------------     ---------------------\n" );
  printf ( "\n" );

  for ( i = 0; i < 10; i++ )
  {
    c1 = c8_uniform_01 ( &seed );
    c2 = c8_tan ( c1 );

    printf ( "  %12.6f%12.6f  %12.6f%12.6f\n", 
      creal ( c1 ), cimag ( c1 ), creal ( c2 ), cimag ( c2 ) );
  }

  return;
}
/******************************************************************************/

void test26 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST26 tests C8_TO_CARTESIAN and CARTESIAN_TO_C8.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2010

  Author:

    John Burkardt
*/
{
  double complex c1;
  double complex c3;
  int i;
  int seed;
  double x2;
  double y2;

  seed = 123456789;

  printf ( "\n" );
  printf ( "TEST26\n" );
  printf ( "  C8_TO_CARTESIAN computes C8 -> ( X, Y ).\n" );
  printf ( "  CARTESIAN_TO_C8 computes ( X, Y ) -> C8.\n" );
  printf ( "\n" );
  printf ( "        C1=C8_UNIFORM_01    (X2,Y2)=C8_TO_CARTESIAN(C1)     C3=CARTESIAN_TO_C8(X2,Y2)\n" );
  printf ( "     ---------------------     ---------------------     ---------------------\n" );
  printf ( "\n" );

  for ( i = 0; i < 10; i++ )
  {
    c1 = c8_uniform_01 ( &seed );
    c8_to_cartesian ( c1, &x2, &y2 );
    c3 = cartesian_to_c8 ( x2, y2 );

    printf ( "  %12.6f%12.6f  %12.6f%12.6f  %12.6f%12.6f\n", 
      creal ( c1 ), cimag ( c1 ), x2, y2, creal ( c3 ), cimag ( c3 ) );
  }

  return;
}
/******************************************************************************/

void test27 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST27 tests C8_TO_POLAR and POLAR_TO_C8.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2010

  Author:

    John Burkardt
*/
{
  double complex c1;
  double complex c3;
  int i;
  double r2;
  int seed;
  double t2;

  seed = 123456789;

  printf ( "\n" );
  printf ( "TEST27\n" );
  printf ( "  C8_TO_POLAR computes C8 -> ( R, T ).\n" );
  printf ( "  POLAR_TO_C8 computes ( R, T ) -> C8.\n" );
  printf ( "\n" );
  printf ( "        C1=C8_UNIFORM_01     (X2,Y2)=C8_TO_POLAR(C1)     C3=POLAR_TO_C8(X2,Y2)\n" );
  printf ( "     ---------------------     ---------------------     ---------------------\n" );

  printf ( "\n" );

  for ( i = 0; i < 10; i++ )
  {
    c1 = c8_uniform_01 ( &seed );
    c8_to_polar ( c1, &r2, &t2 );
    c3 = polar_to_c8 ( r2, t2 );

    printf ( "  %12.6f%12.6f  %12.6f%12.6f  %12.6f%12.6f\n", 
      creal ( c1 ), cimag ( c1 ), r2, t2, creal ( c3 ), cimag ( c3 ) );
  }

  return;
}
/******************************************************************************/

void test28 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST28 tests C8_UNIFORM_01.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2010

  Author:

    John Burkardt
*/
{
  double complex c;
  int i;
  int seed;

  seed = 123456789;

  printf ( "\n" );
  printf ( "TEST28\n" );
  printf ( "  C8_UNIFORM_01 returns a uniformly random \"unit\" C8\n" );
  printf ( "\n" );
  printf ( "       C1=C8_UNIFORM_01(SEED)\n" );
  printf ( "     ---------------------\n" );
  printf ( "\n" );

  for ( i = 0; i < 10; i++ )
  {
    c = c8_uniform_01 ( &seed );

    printf ( "  %12.6f%12.6f\n", creal ( c ), cimag ( c ) );
  }

  return;
}
/******************************************************************************/

void test102 ( )

/******************************************************************************/
/*
  Purpose:

    TEST102 tests R8POLY2_ROOT.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 June 2012

  Author:

    John Burkardt
*/
{
# define TEST_NUM 3

  double a;
  double a_test[TEST_NUM] = { 2.0, 1.0, 1.0 };
  double b;
  double b_test[TEST_NUM] = { -2.0, -20.0, -2.0 };
  double c;
  double c_test[TEST_NUM] = { -24.0, 100.0, 10.0 };
  double complex r1;
  double complex r2;
  int test;

  printf ( "\n" );
  printf ( "TEST102\n" );
  printf ( "  R8POLY2_ROOT finds quadratic equation roots.\n" );
  printf ( "\n" );
  printf ( "         A         B         C     R1         R2\n" );
  printf ( "\n" );

  for ( test = 0; test < TEST_NUM; test++ )
  {
    a = a_test[test];
    b = b_test[test];
    c = c_test[test];

    r8poly2_root ( a, b, c, &r1, &r2 );

    printf ( "  %8g  %8g  (%8g,%8g)  (%g,%g)\n", b, c, creal ( r1 ), cimag ( r1 ),
      creal ( r2 ), cimag ( r2 ) );
  }

  return;
# undef TEST_NUM
}
/******************************************************************************/

void test103 ( )

/******************************************************************************/
/*
  Purpose:

    TEST103 tests R8POLY3_ROOT.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 June 2012

  Author:

    John Burkardt
*/
{
# define TEST_NUM 4

  double a;
  double a_test[TEST_NUM] = { 1.0, 9.0, 1.0, 1.0 };
  double b;
  double b_test[TEST_NUM] = { -6.0, -36.0, -5.0, -8.0 };
  double c;
  double c_test[TEST_NUM] = { 11.0, 54.0, 8.0, 25.0 };
  double d;
  double r8_test[TEST_NUM] = { -6.0, -27.0, -4.0, -26.0 };
  double complex r1;
  double complex r2;
  double complex r3;
  int test;
/*
  1: Three distinct real roots, 1, 2, 3.
  2: One repeated real root, 1.5, 1.5, 1.5.
  3: Two real roots, one repeated, 1, 2, 2.
  4: One real root, a complex conjugate pair, 2, 3+2I, 3-2I.
*/
  printf ( "\n" );
  printf ( "TEST103\n" );
  printf ( "  R8POLY3_ROOT finds roots of cubic equations.\n" );
  printf ( "\n" );

  for ( test = 0; test < TEST_NUM; test++ )
  {
    a = a_test[test];
    b = b_test[test];
    c = c_test[test];
    d = r8_test[test];

    printf ( "\n" );
    printf ( "  Polynomial coefficients A, B, C, D:\n" );
    printf ( "\n" );

    printf ( "  A = %g\n", a );
    printf ( "  B = %g\n", b );
    printf ( "  C = %g\n", c );
    printf ( "  D = %g\n", d );

    r8poly3_root ( a, b, c, d, &r1, &r2, &r3 );

    printf ( "\n" );
    printf ( "  Roots:\n" );
    printf ( "\n" );
    printf ( "  (%g,%g)\n", creal ( r1 ), cimag ( r1 ) );
    printf ( "  (%g,%g)\n", creal ( r2 ), cimag ( r2 ) );
    printf ( "  (%g,%g)\n", creal ( r3 ), cimag ( r3 ) );
  }

  return;
# undef TEST_NUM
}
/******************************************************************************/

void test104 ( )

/******************************************************************************/
/*
  Purpose:

    TEST104 tests R8POLY4_ROOT.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 June 2012

  Author:

    John Burkardt
*/
{
# define TEST_NUM 7

  double a;
  double a_test[TEST_NUM] = {
    1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0 };
  double b;
  double b_test[TEST_NUM] = {
    -10.0, -5.0, -22.0, -16.0, -20.0,
    2.0, 0.0 };
  double c;
  double c_test[TEST_NUM] = {
    35.0, 1.0, 141.0, 72.0, 150.0,
    1.0, 13.0 };
  double d;
  double r8_test[TEST_NUM] = {
    -50.0, 21.0, -220.0, -128.0, -500.0,
    8.0, 0.0 };
  double e;
  double e_test[TEST_NUM] = {
    24.0, -18.0, +100.0, 80.0, 625.0,
    -12.0, 36.0 };
  double complex r1;
  double complex r2;
  double complex r3;
  double complex r4;
  int test;
/*
  1: Four distinct real roots, 1, 2, 3, 4.
  2: Three distinct real roots, 1, -2, 3, 3
  3: Two distinct real roots, 1, 1, 10, 10.
  4: Two distinct real roots, 2, 2, 2, 10
  5: One real root, 5, 5, 5, 5
  6: Two distinct real roots, one complex conjugate pair.
  7: Two distinct complex conjugate pairs.
*/
  printf ( "\n" );
  printf ( "TEST104\n" );
  printf ( "  R8POLY4_ROOT finds roots of quartic equations.\n" );
  printf ( "\n" );

  for ( test = 0; test < TEST_NUM; test++ )
  {
    a = a_test[test];
    b = b_test[test];
    c = c_test[test];
    d = r8_test[test];
    e = e_test[test];

    printf ( "\n" );
    printf ( "  A = %g\n", a );
    printf ( "  B = %g\n", b );
    printf ( "  C = %g\n", c );
    printf ( "  D = %g\n", d );
    printf ( "  E = %g\n", e );

    r8poly4_root ( a, b, c, d, e, &r1, &r2, &r3, &r4 );

    printf ( "\n" );
    printf ( "  Roots:\n" );
    printf ( "\n" );
    printf ( "  (%g,%g)\n", creal ( r1 ), cimag ( r1 ) );
    printf ( "  (%g,%g)\n", creal ( r2 ), cimag ( r2 ) );
    printf ( "  (%g,%g)\n", creal ( r3 ), cimag ( r3 ) );
    printf ( "  (%g,%g)\n", creal ( r4 ), cimag ( r4 ) );
  }

  return;
# undef TEST_NUM
}
