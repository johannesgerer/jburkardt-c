# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <complex.h>

# include "blas1_c.h"

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
_Complex float c4_uniform_01 ( int *seed );
void timestamp ( void );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for BLAS1_C_PRB.

  Discussion:

    BLAS1_C_PRB tests the BLAS1 single precision complex routines.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 March 2007

  Author:

    John Burkardt
*/
{
  timestamp ( );

  printf ( "\n" );
  printf ( "BLAS1_C_PRB:\n" );
  printf ( "  C version\n" );
  printf ( "  Test the BLAS1_C library.\n" );

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
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "BLAS1_C_PRB:\n" );
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

    TEST01 tests CABS1.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 March 2007

  Author:

    John Burkardt
*/
{
  _Complex float c;
  float c_norm;
  int i;
  int seed = 123456789;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  CABS1 returns the L1 norm of a complex number.\n" );
  printf ( "\n" );
  printf ( "      Real      Imaginary\n" );
  printf ( "      Part      Part           CABS1(Z)\n" );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    c = c4_uniform_01 ( &seed );

    c = 5.0 * c;

    c_norm = cabs1 ( c );

    printf ( "  %10f  %10f    %10f\n", crealf ( c ), cimagf ( c ), c_norm );
  }

  return;
}
/******************************************************************************/

void test02 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST02 tests CABS2.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 March 2007

  Author:

    John Burkardt
*/
{
  _Complex float c;
  float c_norm;
  int i;
  int seed = 123456789;

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  CABS2 returns the L2 norm of a complex number.\n" );
  printf ( "\n" );
  printf ( "      Real      Imaginary\n" );
  printf ( "      Part      Part           CABS2(Z)\n" );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    c = c4_uniform_01 ( &seed );

    c = 5.0 * c;

    c_norm = cabs2 ( c );

    printf ( "  %10f  %10f    %10f\n", crealf ( c ), cimagf ( c ), c_norm );
  }

  return;
}
/******************************************************************************/

void test03 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST03 tests CAXPY.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 March 2007

  Author:

    John Burkardt
*/
{
# define N 5

  int i;
  _Complex float s;
  _Complex float x[N] = {
     2.0 - 1.0 * _Complex_I,
    -4.0 - 2.0 * _Complex_I,
     3.0 + 1.0 * _Complex_I,
     2.0 + 2.0 * _Complex_I,
    -1.0 - 1.0 * _Complex_I };
  _Complex float y[N] = {
    -1.0 + 0.0 * _Complex_I,
     0.0 - 3.0 * _Complex_I,
     4.0 + 0.0 * _Complex_I,
    -3.0 + 4.0 * _Complex_I,
    -2.0 + 0.0 * _Complex_I };

  printf ( "\n" );
  printf ( "TEST03\n" );
  printf ( "  CAXPY adds a multiple of one complex vector to another.\n" );

  printf ( "\n" );
  printf ( "  X =\n" );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    printf ( "  %6d  %6f  %6f\n", i, crealf ( x[i] ), cimagf ( x[i] ) );
  }

  printf ( "\n" );
  printf ( "  Y =\n" );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    printf ( "  %6d  %6f  %6f\n", i, crealf ( y[i] ), cimagf ( y[i] ) );
  }

  s = 0.50 - 1.00 * _Complex_I;

  printf ( "\n" );
  printf ( "  The scalar multiplier is: %f  %f\n", crealf ( s ), cimagf ( s ) );

  caxpy ( N, s, x, 1, y, 1 );

  printf ( "\n" );
  printf ( "  A * X + Y =\n" );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    printf ( "  %6d  %6f  %6f\n", i, crealf ( y[i] ), cimagf ( y[i] ) );
  }

  return;
# undef N
}
/******************************************************************************/

void test04 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST04 tests CCOPY.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 March 2007

  Author:

    John Burkardt
*/
{
# define N1 5
# define N2 5
# define N 10

  _Complex float a[N1*N2];
  int i;
  int j;
  _Complex float x[N];
  _Complex float y[N];

  printf ( "\n" );
  printf ( "TEST04\n" );
  printf ( "  CCOPY copies one complex vector into another.\n" );

  for ( i = 0; i < N; i++ )
  {
    x[i] = 10 * ( i + 1 ) + ( i + 1 ) * _Complex_I;
  }

  for ( i = 0; i < N; i++ )
  {
    y[i] = 20 * ( i + 1 ) + 2 * ( i + 1 ) * _Complex_I;
  }

  for ( i = 0; i < N1; i++ )
  {
    for ( j = 0; j < N2; j++ )
    {
      a[i+j*N1] = 10 * ( i + 1 ) + ( j + 1 ) * _Complex_I;
    }
  }

  printf ( "\n" );
  printf ( "  X =\n" );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    printf ( "  %6d  %6f  %6f\n", i, crealf ( x[i] ), cimagf ( x[i] ) );
  }
  printf ( "\n" );
  printf ( "  Y =\n" );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    printf ( "  %6d  %6f  %6f\n", i, crealf ( y[i] ), cimagf ( y[i] ) );
  }
  printf ( "\n" );
  printf ( "  A =\n" );
  printf ( "\n" );
  for ( i = 0; i < N1; i++ )
    {
    for ( j = 0; j < N2; j++ )
    {
      printf ( "  %5f  %5f\n", crealf ( a[i+j*N1] ), cimagf ( a[i+j*N1] ) );
    }
    printf ( "\n" );
  }

  ccopy ( 5, x, 1, y, 1 );
  printf ( "\n" );
  printf ( "  CCOPY ( 5, X, 1, Y, 1 )\n" );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    printf ( "  %6d  %6f  %6f\n", i, crealf ( y[i] ), cimagf ( y[i] ) );
  }

  for ( i = 0; i < N; i++ )
  {
    y[i] = 20 * ( i + 1 ) + 2 * ( i + 1 ) * _Complex_I;
  }

  ccopy ( 3, x, 2, y, 3 );

  printf ( "\n" );
  printf ( "  CCOPY ( 3, X, 2, Y, 3 )\n" );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    printf ( "  %6d  %6f  %6f\n", i, crealf ( y[i] ), cimagf ( y[i] ) );
  }

  ccopy ( 5, x, 1, a, 1 );

  printf ( "\n" );
  printf ( "  CCOPY ( 5, X, 1, A, 1 )\n" );
  printf ( "\n" );
  printf ( "\n" );
  printf ( "  A =\n" );
  printf ( "\n" );
  for ( i = 0; i < N1; i++ )
    {
    for ( j = 0; j < N2; j++ )
    {
      printf ( "  %5f  %5f\n", crealf ( a[i+j*N1] ), cimagf ( a[i+j*N1] ) );
    }
    printf ( "\n" );
  }

  for ( i = 0; i < N1; i++ )
  {
    for ( j = 0; j < N2; j++ )
    {
      a[i+j*N1] = 10 * ( i + 1 ) + ( j + 1 ) * _Complex_I;
    }
  }

  ccopy ( 5, x, 2, a, 5 );

  printf ( "\n" );
  printf ( "  CCOPY ( 5, X, 2, A, 5 )\n" );
  printf ( "\n" );
  printf ( "  A =\n" );
  printf ( "\n" );
  for ( i = 0; i < N1; i++ )
    {
    for ( j = 0; j < N2; j++ )
    {
      printf ( "  %5f  %5f\n", crealf ( a[i+j*N1] ), cimagf ( a[i+j*N1] ) );
    }
    printf ( "\n" );
  }

  return;
# undef N
# undef N1
# undef N2
}
/******************************************************************************/

void test05 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST05 tests CDOTC.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 March 2007

  Author:

    John Burkardt
*/
{
# define N 5

  int i;
  _Complex float x_norm;
  _Complex float xy_dot;
  _Complex float x[N] = {
     2.0 - 1.0 * _Complex_I,
    -4.0 - 2.0 * _Complex_I,
     3.0 + 1.0 * _Complex_I,
     2.0 + 2.0 * _Complex_I,
    -1.0 - 1.0 * _Complex_I };
  _Complex float y[N] = {
    -1.0 + 0.0 * _Complex_I,
     0.0 - 3.0 * _Complex_I,
     4.0 + 0.0 * _Complex_I,
    -3.0 + 4.0 * _Complex_I,
    -2.0 + 0.0 * _Complex_I };

  printf ( "\n" );
  printf ( "TEST05\n" );
  printf ( "  CDOTC computes the conjugated dot product of\n" );
  printf ( "  two complex vectors.\n" );

  printf ( "\n" );
  printf ( "  X =\n" );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    printf ( "  %6d  %6f  %6f\n", i, crealf ( x[i] ), cimagf ( x[i] ) );
  }

  x_norm = cdotc ( N, x, 1, x, 1 );

  printf ( "\n" );
  printf ( "  The square of the norm of X, computed as\n" );
  printf ( "  CDOTC(X,X) = (%f,  %f)\n", crealf ( x_norm ), cimagf ( x_norm ) );

  xy_dot = cdotc ( N, x, 1, y, 1 );

  printf ( "\n" );
  printf ( "  Y =\n" );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    printf ( "  %6d  %6f  %6f\n", i, crealf ( y[i] ), cimagf ( y[i] ) );
  }
  printf ( "\n" );
  printf ( "  The dot product X.Y* is (%f,  %f)\n",
    crealf ( xy_dot ), cimagf ( xy_dot ) );

  return;
# undef N
}
/******************************************************************************/

void test06 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST06 tests CDOTU.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 March 2007

  Author:

    John Burkardt
*/
{
# define N 5

  int i;
  _Complex float x_norm;
  _Complex float xy_dot;
  _Complex float x[N] = {
    2.0 - 1.0 * _Complex_I,
   -4.0 - 2.0 * _Complex_I,
    3.0 + 1.0 * _Complex_I,
    2.0 + 2.0 * _Complex_I,
   -1.0 - 1.0 * _Complex_I };
  _Complex float y[N] = {
   -1.0 + 0.0 * _Complex_I,
    0.0 - 3.0 * _Complex_I,
    4.0 + 0.0 * _Complex_I,
   -3.0 + 4.0 * _Complex_I,
   -2.0 + 0.0 * _Complex_I };

  printf ( "\n" );
  printf ( "TEST06\n" );
  printf ( "  CDOTU computes the unconjugated dot product of\n" );
  printf ( "  two complex vectors.\n" );

  printf ( "\n" );
  printf ( "  X =\n" );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    printf ( "  %6d  %6f  %6f\n", i, crealf ( x[i] ), cimagf ( x[i] ) );
  }

  x_norm = cdotu ( N, x, 1, x, 1 );

  printf ( "\n" );
  printf ( "  The unconjugated dot product ( X dot X )\n" );
  printf ( "  (which is NOT the square of the norm of X!):\n" );
  printf ( "  CDOTU(X,X) = (%f,  %f)\n", crealf ( x_norm ), cimagf ( x_norm ) );

  xy_dot = cdotu ( N, x, 1, y, 1 );

  printf ( "\n" );
  printf ( "  Y =\n" );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    printf ( "  %6d  %6f  %6f\n", i, crealf ( y[i] ), cimagf ( y[i] ) );
  }

  printf ( "\n" );
  printf ( "  The dot product ( X dot Y ) is (%f,  %f)\n",
    crealf ( xy_dot ), cimagf ( xy_dot ) );

  return;
# undef N
}
/******************************************************************************/

void test07 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST07 tests CMACH.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 March 2007

  Author:

    John Burkardt
*/
{
  printf ( "\n" );
  printf ( "TEST07\n" );
  printf ( "  CMACH computes several machine-dependent\n" );
  printf ( "  complex arithmetic parameters.\n" );

  printf ( "\n" );
  printf ( "  CMACH(1)  = machine epsilon = %e\n", cmach ( 1 ) );
  printf ( "  CMACH(2)  = a tiny value    = %e\n", cmach ( 2 ) );
  printf ( "  CMACH(3)  = a huge value    = %e\n", cmach ( 3 ) );

  return;
}
/******************************************************************************/

void test08 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST08 tests CROTG.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 March 2007

  Author:

    John Burkardt
*/
{
  _Complex float a;
  _Complex float b;
  float c;
  _Complex float r;
  _Complex float s;
  _Complex float sa;
  _Complex float sb;
  int seed;
  _Complex float t;
  int test;
  int test_num = 5;

  printf ( "\n" );
  printf ( "TEST08\n" );
  printf ( "  CROTG generates a complex Givens rotation\n" );
  printf ( "    (  C  S ) * ( A ) = ( R )\n" );
  printf ( "    ( -S  C )   ( B )   ( 0 )\n" );
  printf ( "\n" );

  seed = 123456789;

  for ( test = 1; test <= test_num; test++ )
  {
    a = c4_uniform_01 ( &seed );
    b = c4_uniform_01 ( &seed );

    sa = a;
    sb = b;

    crotg ( &sa, sb, &c, &s );

    r = sa;

    printf ( "\n" );
    printf ( "  A =  ( %f. %f )\n", crealf ( a ), cimagf ( a ) );
    printf ( "  B =  ( %f. %f )\n", crealf ( b ), cimagf ( b ) );
    printf ( "  C =    %f\n", c );
    printf ( "  S =  ( %f. %f )\n", crealf ( s ), cimagf ( s ) );
    printf ( "  R =  ( %f. %f )\n", crealf ( r ), cimagf ( r ) );
    t = c * a + s * b;
    printf ( "         C *A+S*B = ( %f. %f )\n", crealf ( t ), cimagf ( t ) );
    t = - ( ~s ) * a + c * b;
    printf ( "  -conjg(S)*A+C*B = ( %f. %f )\n", crealf ( t ), cimagf ( t ) );
  }

  return;
}
/******************************************************************************/

void test09 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST09 tests CSCAL.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 March 2007

  Author:

    John Burkardt
*/
{
# define N 6

  _Complex float da;
  int i;
  _Complex float x[N];

  for ( i = 0; i < N; i++ )
  {
    x[i] = 10 * ( i + 1 ) + ( i + 1 ) * _Complex_I;
  }

  printf ( "\n" );
  printf ( "TEST09\n" );
  printf ( "  CSCAL multiplies a complex scalar times a vector.\n" );

  printf ( "\n" );
  printf ( "  X =\n" );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    printf ( "  %6d  %6f  %6f\n", i, crealf ( x[i] ), cimagf ( x[i] ) );
  }

  da = 5.0;
  cscal ( N, da, x, 1 );
  printf ( "\n" );
  printf ( "  CSCAL ( N, (%f, %f), X, 1 )\n", crealf ( da ), cimagf ( da ) );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    printf ( "  %6d  %6f  %6f\n", i, crealf ( x[i] ), cimagf ( x[i] ) );
  }

  for ( i = 0; i < N; i++ )
  {
    x[i] = 10 * ( i + 1 ) + ( i + 1 ) * _Complex_I;
  }

  da = -2.0 + 1.0 * _Complex_I;
  cscal ( 3, da, x, 2 );
  printf ( "\n" );
  printf ( "  CSCAL ( 3, (%f, %f), X, 2 )\n", crealf ( da ), cimagf ( da ) );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    printf ( "  %6d  %6f  %6f\n", i, crealf ( x[i] ), cimagf ( x[i] ) );
  }

  return;
# undef N
}
/******************************************************************************/

void test10 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST10 tests CSIGN1.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 March 2007

  Author:

    John Burkardt
*/
{
  _Complex float c1;
  _Complex float c2;
  _Complex float c3;
  int i;
  int seed = 123456789;

  printf ( "\n" );
  printf ( "TEST10\n" );
  printf ( "  CSIGN1 ( C1, C2 ) transfers the sign of complex C2\n" );
  printf ( "  to the CABS1 magnitude of C1.\n" );
  printf ( "\n" );
  printf ( "           C1                    C2                    C3\n" );
  printf (
    "  --------------------  --------------------  --------------------\n" );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    c1 = 5.0 * c4_uniform_01 ( &seed );
    c2 = 5.0 * c4_uniform_01 ( &seed );
    c3 = csign1 ( c1, c2 );

    printf ( "  (%9f  %9f)  (%9f  %9f)   (%9f  %9f)\n",
      crealf ( c1 ), cimagf ( c1 ),
      crealf ( c2 ), cimagf ( c2 ),
      crealf ( c3 ), cimagf ( c3 ) );
  }

  return;
}
/******************************************************************************/

void test11 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST11 tests CSIGN2.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 March 2007

  Author:

    John Burkardt
*/
{
  _Complex float c1;
  _Complex float c2;
  _Complex float c3;
  int i;
  int seed = 123456789;

  printf ( "\n" );
  printf ( "TEST11\n" );
  printf ( "  CSIGN2 ( C1, C2 ) transfers the sign of complex C2\n" );
  printf ( "  to the CABS2 magnitude of C1.\n" );
  printf ( "\n" );
  printf (
    "           C1                    C2                    C3\n" );
  printf (
    "  --------------------  --------------------  --------------------\n" );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    c1 = 5.0 * c4_uniform_01 ( &seed );
    c2 = 5.0 * c4_uniform_01 ( &seed );
    c3 = csign2 ( c1, c2 );

    printf ( "  (%9f  %9f)  (%9f  %9f)   (%9f  %9f)\n",
      crealf ( c1 ), cimagf ( c1 ),
      crealf ( c2 ), cimagf ( c2 ),
      crealf ( c3 ), cimagf ( c3 ) );
  }

  return;
}
/******************************************************************************/

void test12 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST12 tests CSROT.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 March 2007

  Author:

    John Burkardt
*/
{
# define N 6

  float c;
  int i;
  float s;
  _Complex float x[N];
  _Complex float y[N];

  for ( i = 0; i < N; i++ )
  {
    x[i] = 10 * ( i + 1 ) + ( i + 1 ) * _Complex_I;
  }

  for ( i = 0; i < N; i++ )
  {
    y[i] = 20 * ( i + 1 ) + 2 * ( i + 1 ) * _Complex_I;
  }

  printf ( "\n" );
  printf ( "TEST12\n" );
  printf ( "  CSROT carries out a Givens rotation\n" );
  printf ( "  on a complex vector.\n" );
  printf ( "\n" );
  printf ( "  X and Y\n" );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    printf ( "  %6d  (%9f  %9f)  (%9f  %9f)\n",
    i, crealf ( x[i] ), cimagf ( x[i] ), crealf ( y[i] ), cimagf ( y[i] ) );
  }

  c = 0.5;
  s = sqrt ( 1.0 - c * c );
  csrot ( N, x, 1, y, 1, c, s );
  printf ( "\n" );
  printf ( "  CSROT ( N, X, 1, Y, 1, %f, %f )\n", c, s );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    printf ( "  %6d  (%9f  %9f)  (%9f  %9f)\n",
    i, crealf ( x[i] ), cimagf ( x[i] ), crealf ( y[i] ), cimagf ( y[i] ) );
  }

  return;
# undef N
}
/******************************************************************************/

void test13 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST13 tests CSSCAL.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 March 2007

  Author:

    John Burkardt
*/
{
# define N 6

  float da;
  int i;
  _Complex float x[N];

  for ( i = 0; i < N; i++ )
  {
    x[i] = 10 * ( i + 1 ) + ( i + 1 ) * _Complex_I;
  }

  printf ( "\n" );
  printf ( "TEST13\n" );
  printf ( "  CSSCAL multiplies a real scalar times a complex vector.\n" );

  printf ( "\n" );
  printf ( "  X =\n" );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    printf ( "  %6d  %6f  %6f\n", i, crealf ( x[i] ), cimagf ( x[i] ) );
  }

  da = 5.0;
  csscal ( N, da, x, 1 );
  printf ( "\n" );
  printf ( "  CSSCAL ( N, %f, X, 1 )\n", da );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    printf ( "  %6d  %6f  %6f\n", i, crealf ( x[i] ), cimagf ( x[i] ) );
  }

  for ( i = 0; i < N; i++ )
  {
    x[i] = 10 * ( i + 1 ) + ( i + 1 ) * _Complex_I;
  }

  da = -2.0;
  csscal ( 3, da, x, 2 );
  printf ( "\n" );
  printf ( "  CSSCAL ( 3, %f, X, 2 )\n", da );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    printf ( "  %6d  %6f  %6f\n", i, crealf ( x[i] ), cimagf ( x[i] ) );
  }

  return;
# undef N
}
/******************************************************************************/

void test14 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST14 tests CSWAP.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 March 2007

  Author:

    John Burkardt
*/
{
# define N 5

  int i;
  _Complex float x[N];
  _Complex float y[N];

  for ( i = 0; i < N; i++ )
  {
    x[i] = 10 * ( i + 1 ) + ( i + 1 ) * _Complex_I;
  }

  for ( i = 0; i < N; i++ )
  {
    y[i] = 20 * ( i + 1 ) + 2 * ( i + 1 ) * _Complex_I;
  }

  printf ( "\n" );
  printf ( "TEST14\n" );
  printf ( "  CSWAP swaps two complex vectors.\n" );
  printf ( "\n" );
  printf ( "  X and Y\n" );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    printf ( "  %6d  %6f  %6f    %6f  %6f\n",
    i, crealf ( x[i] ), cimagf ( x[i] ), crealf ( y[i] ), cimagf ( y[i] ) );
  }

  cswap ( N, x, 1, y, 1 );
  printf ( "\n" );
  printf ( "  CSWAP ( N, X, 1, Y, 1 )\n" );
  printf ( "\n" );
  printf ( "  X and Y\n" );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    printf ( "  %6d  %6f  %6f    %6f  %6f\n",
    i, crealf ( x[i] ), cimagf ( x[i] ), crealf ( y[i] ), cimagf ( y[i] ) );
  }

  for ( i = 0; i < N; i++ )
  {
    x[i] = 10 * ( i + 1 ) + ( i + 1 ) * _Complex_I;
  }

  for ( i = 0; i < N; i++ )
  {
    y[i] = 20 * ( i + 1 ) + 2 * ( i + 1 ) * _Complex_I;
  }

  cswap ( 3, x, 2, y, 1 );
  printf ( "\n" );
  printf ( "  CSWAP ( 3, X, 2, Y, 1 )\n" );
  printf ( "\n" );
  printf ( "  X and Y\n" );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    printf ( "  %6d  %6f  %6f    %6f  %6f\n",
    i, crealf ( x[i] ), cimagf ( x[i] ), crealf ( y[i] ), cimagf ( y[i] ) );
  }

  return;
}
/******************************************************************************/

void test15 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST15 tests ICAMAX.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 March 2007

  Author:

    John Burkardt
*/
{
# define N 5

  int i;
  int incx;
  _Complex float x[N] = {
     2.0 - 1.0 * _Complex_I,
    -4.0 - 2.0 * _Complex_I,
     3.0 + 1.0 * _Complex_I,
     2.0 + 2.0 * _Complex_I,
    -1.0 - 1.0 * _Complex_I };

  printf ( "\n" );
  printf ( "TEST15\n" );
  printf ( "  ICAMAX returns the index of the entry of\n" );
  printf ( "  maximum magnitude in a complex vector.\n" );

  printf ( "\n" );
  printf ( "  The entries and CABS1 magnitudes:\n" );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    printf ( "  %6d  %6f  %6f    %6f\n",
    i, crealf ( x[i] ), cimagf ( x[i] ), cabs1 ( x[i] ) );
  }

  incx = 1;

  i = icamax ( N, x, incx );

  printf ( "\n" );
  printf ( "  The index of maximum magnitude = %d\n", i );
  printf ( "\n" );
  printf ( "  Note that this is a 1-based index.\n" );
  printf ( "  Note that the L1 norm is used.\n" );

  return;
# undef N
}
/******************************************************************************/

void test16 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST16 tests SCASUM.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 March 2007

  Author:

    John Burkardt
*/
{
# define MA 5
# define NA 4
# define NX 8

  _Complex float a[MA*NA] = {
    -3.0 + 4.0 * _Complex_I,
     2.0 + 0.0 * _Complex_I,
     3.0 - 4.0 * _Complex_I,
     2.0 + 0.0 * _Complex_I,
     2.0 - 1.0 * _Complex_I,
    -1.0 + 1.0 * _Complex_I,
     0.0 + 5.0 * _Complex_I,
    -4.0 - 2.0 * _Complex_I,
    -4.0 + 1.0 * _Complex_I,
    -4.0 - 3.0 * _Complex_I,
     0.0 - 2.0 * _Complex_I,
     1.0 + 3.0 * _Complex_I,
    -3.0 + 3.0 * _Complex_I,
    -3.0 + 3.0 * _Complex_I,
    -1.0 - 2.0 * _Complex_I,
    -1.0 + 2.0 * _Complex_I,
     2.0 - 4.0 * _Complex_I,
     0.0 - 1.0 * _Complex_I,
     0.0 - 1.0 * _Complex_I,
    -2.0 + 4.0 * _Complex_I };
  int i;
  int j;
  _Complex float x[NX] = {
    2.0 - 1.0 * _Complex_I,
   -4.0 - 2.0 * _Complex_I,
    3.0 + 1.0 * _Complex_I,
    2.0 + 2.0 * _Complex_I,
   -1.0 - 1.0 * _Complex_I,
   -1.0 + 0.0 * _Complex_I,
    0.0 - 3.0 * _Complex_I,
    4.0 + 0.0 * _Complex_I };

  printf ( "\n" );
  printf ( "TEST16\n" );
  printf ( "  SCASUM adds the absolute values of elements\n" );
  printf ( "  of a complex vector.\n" );
  printf ( "\n" );
  printf ( "  X =\n" );
  printf ( "\n" );
  for ( i = 0; i < NX; i++ )
  {
    printf ( "  %6d  (%6.1f  %6.1f )\n", i, crealf ( x[i] ), cimagf ( x[i] ) );
  }

  printf ( "\n" );
  printf ( "  SCASUM ( NX,   X, 1    ) = %f\n",
    scasum ( NX,   x, 1 ) );

  printf ( "  SCASUM ( NX/2, X, 2    ) = %f\n",
    scasum ( NX/2, x, 2 ) );

  printf ( "  SCASUM ( 2,    X, NX/2 ) = %f\n",
    scasum ( 2,    x, NX/2 ) );

  printf ( "\n" );
  printf ( "  Demonstrate with a matrix A:\n" );
  printf ( "\n" );
  for ( i = 0; i < MA; i++ )
  {
    for ( j = 0; j < NA; j++ )
    {
      printf ( "  (%6.1f  %6.1f )",
        crealf ( a[i+j*MA] ), cimagf ( a[i+j*MA] ) );
    }
    printf ( "\n" );
  }

  printf ( "\n" );
  printf ( "  SCASUM ( MA, A[1,2], 1 )   = %f\n",
    scasum ( MA, a+0+1*MA, 1 ) );

  printf ( "  SCASUM ( NA, A[2,1], MA ) = %f\n",
    scasum ( NA, a+1+0*MA, MA ) );

  return;
# undef MA
# undef NA
# undef NX
}
/******************************************************************************/

void test17 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST17 tests SCNRM2.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 March 2007

  Author:

    John Burkardt
*/
{
# define N 5

  int i;
  int incx;
  float norm;
  _Complex float x[N] = {
    2.0 - 1.0 * _Complex_I,
   -4.0 - 2.0 * _Complex_I,
    3.0 + 1.0 * _Complex_I,
    2.0 + 2.0 * _Complex_I,
   -1.0 - 1.0 * _Complex_I };

  printf ( "\n" );
  printf ( "TEST17\n" );
  printf ( "  SCNRM2 returns the Euclidean norm of a complex vector.\n" );

  printf ( "\n" );
  printf ( "  The vector X:\n" );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    printf ( "  %6d  %6f  %6f\n",
    i, crealf ( x[i] ), cimagf ( x[i] ) );
  }

  incx = 1;
  norm = scnrm2 ( N, x, incx );

  printf ( "\n" );
  printf ( "  The L2 norm of X is %f\n", norm );

  return;
# undef N
}
/******************************************************************************/

_Complex float c4_uniform_01 ( int *seed )

/******************************************************************************/
/*
  Purpose:

    C4_UNIFORM_01 returns a unit complex pseudorandom number.

  Discussion:

    The angle should be uniformly distributed between 0 and 2 * PI,
    the square root of the radius uniformly distributed between 0 and 1.

    This results in a uniform distribution of values in the unit circle.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 March 2007

  Author:

    John Burkardt

  Parameters:

    Input/output, int *SEED, the seed value, which should NOT be 0.
    On output, SEED has been updated.

    Output, complex C4_UNIFORM_01, a pseudorandom complex value.
*/
{
  float r;
  int k;
  float pi = 3.1415926;
  float theta;
  _Complex float value;

  k = *seed / 127773;

  *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

  if ( *seed < 0 )
  {
    *seed = *seed + 2147483647;
  }

  r = sqrt ( ( float ) ( ( double ) ( *seed ) * 4.656612875E-10 ) );

  k = *seed / 127773;

  *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

  if ( *seed < 0 )
  {
    *seed = *seed + 2147483647;
  }

  theta = 2.0 * pi * ( float )
    ( ( double ) ( *seed ) * 4.656612875E-10 );

  value = r * cos ( theta ) + ( r * sin ( theta ) ) * _Complex_I;

  return value;
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
