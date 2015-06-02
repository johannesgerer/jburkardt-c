# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <complex.h>

# include "blas1_z.h"

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
double complex c8_uniform_01 ( int *seed );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for BLAS1_Z_PRB.

  Discussion:

    BLAS1_Z_PRB tests the BLAS1_Z library.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    31 March 2007

  Author:

    John Burkardt
*/
{
  timestamp ( );

  printf ( "\n" );
  printf ( "BLAS1_Z_PRB:\n" );
  printf ( "  C version\n" );
  printf ( "  Test the BLAS1_Z library.\n" );
 
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
  printf ( "BLAS1_Z_PRB:\n" );
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

    TEST01 tests ZABS1.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    31 March 2007

  Author:

    John Burkardt
*/
{
  double complex c;
  double c_norm;
  int i;
  int seed = 123456789;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  ZABS1 returns the L1 norm of a complex number.\n" );
  printf ( "\n" );
  printf ( "      Real      Imaginary\n" );
  printf ( "      Part      Part           ZABS1(Z)\n" );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    c = c8_uniform_01 ( &seed );
    c = 5.0 * c;

    c_norm = zabs1 ( c );

    printf ( "  %10f  %10f    %10f\n", creal ( c ), cimag ( c ), c_norm );
  }

  return;
}
/******************************************************************************/

void test02 ( )

/******************************************************************************/
/*
  Purpose:

    TEST02 tests ZABS2.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    31 March 2007

  Author:

    John Burkardt
*/
{
  double complex c;
  double c_norm;
  int i;
  int seed = 123456789;

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  ZABS2 returns the L2 norm of a complex number.\n" );
  printf ( "\n" );
  printf ( "      Real      Imaginary\n" );
  printf ( "      Part      Part           ZABS2(Z\n" );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    c = c8_uniform_01 ( &seed );
    c = 5.0 * c;

    c_norm = zabs2 ( c );

    printf ( "  %10f  %10f    %10f\n", creal ( c ), cimag ( c ), c_norm );
  }

  return;
}
/******************************************************************************/

void test03 ( )

/******************************************************************************/
/*
  Purpose:

    TEST03 tests ZAXPY.

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
  double complex s;
  double complex x[N] = {
     2.0 - 1.0 * I, 
    -4.0 - 2.0 * I, 
     3.0 + 1.0 * I, 
     2.0 + 2.0 * I, 
    -1.0 - 1.0 * I };
  double complex y[N] = {
    -1.0 + 0.0 * I, 
     0.0 - 3.0 * I,
     4.0 + 0.0 * I, 
    -3.0 + 4.0 * I, 
    -2.0 + 0.0 * I };

  printf ( "\n" );
  printf ( "TEST03\n" );
  printf ( "  ZAXPY adds a multiple of one complex vector to another.\n" );
 
  printf ( "\n" );
  printf ( "  X =\n" );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    printf ( "  %6d  %6f  %6f\n", i, creal ( x[i] ), cimag ( x[i] ) );
  }

  printf ( "\n" );
  printf ( "  Y =\n" );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    printf ( "  %6d  %6f  %6f\n", i, creal ( y[i] ), cimag ( y[i] ) );
  }

  s = 0.50 - 1.00 * I;

  printf ( "\n" );
  printf ( "  The scalar multiplier is: %f  %f\n", creal ( s ), cimag ( s ) );

  zaxpy ( N, s, x, 1, y, 1 );

  printf ( "\n" );
  printf ( "  A * X + Y =\n" );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    printf ( "  %6d  %6f  %6f\n", i, creal ( y[i] ), cimag ( y[i] ) );
  }

  return;
# undef N
}
/******************************************************************************/

void test04 ( )

/******************************************************************************/
/*
  Purpose:

    TEST04 tests ZCOPY.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    31 March 2007

  Author:

    John Burkardt
*/
{
# define N1 5
# define N2 5
# define N 10

  double complex a[N1*N2];
  int i;
  int j;
  double complex x[N];
  double complex y[N];

  printf ( "\n" );
  printf ( "TEST04\n" );
  printf ( "  ZCOPY copies one complex vector into another.\n" );
 
  for ( i = 0; i < N; i++ )
  {
    x[i] = 10 * ( i + 1 ) + ( i + 1 ) * I;
  }

  for ( i = 0; i < N; i++ )
  {
    y[i] = 20 * ( i + 1 ) + 2 * ( i + 1 ) * I;
  }

  for ( i = 0; i < N1; i++ )
  {
    for ( j = 0; j < N2; j++ )
    {
      a[i+j*N1] = 10 * ( i + 1 ) + ( j + 1 ) * I;
    }
  }

  printf ( "\n" );
  printf ( "  X =\n" );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    printf ( "  %6d  %6f  %6f\n", i, creal ( x[i] ), cimag ( x[i] ) );
  }
  printf ( "\n" );
  printf ( "  Y =\n" );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    printf ( "  %6d  %6f  %6f\n", i, creal ( y[i] ), cimag ( y[i] ) );
  }
  printf ( "\n" );
  printf ( "  A =\n" );
  printf ( "\n" );
  for ( i = 0; i < N1; i++ )
    {
    for ( j = 0; j < N2; j++ )
    {
      printf ( "  %5f  %5f\n", creal ( a[i+j*N1] ), cimag ( a[i+j*N1] ) );
    }
    printf ( "\n" );
  }

  zcopy ( 5, x, 1, y, 1 );
  printf ( "\n" );
  printf ( "  ZCOPY ( 5, X, 1, Y, 1 )\n" );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    printf ( "  %6d  %6f  %6f\n", i, creal ( y[i] ), cimag ( y[i] ) );
  }

  for ( i = 0; i < N; i++ )
  {
    y[i] = 20 * ( i + 1 ) + 2 * ( i + 1 ) * I;
  }

  zcopy ( 3, x, 2, y, 3 );

  printf ( "\n" );
  printf ( "  ZCOPY ( 3, X, 2, Y, 3 )\n" );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    printf ( "  %6d  %6f  %6f\n", i, creal ( y[i] ), cimag ( y[i] ) );
  }

  zcopy ( 5, x, 1, a, 1 );

  printf ( "\n" );
  printf ( "  ZCOPY ( 5, X, 1, A, 1 )\n" );
  printf ( "\n" );
  printf ( "\n" );
  printf ( "  A =\n" );
  printf ( "\n" );
  for ( i = 0; i < N1; i++ )
    {
    for ( j = 0; j < N2; j++ )
    {
      printf ( "  %5f  %5f\n", creal ( a[i+j*N1] ), cimag ( a[i+j*N1] ) );
    }
    printf ( "\n" );
  }

  for ( i = 0; i < N1; i++ )
  {
    for ( j = 0; j < N2; j++ )
    {
      a[i+j*N1] = 10 * ( i + 1 ) + ( j + 1 ) * I;
    }
  }

  zcopy ( 5, x, 2, a, 5 );

  printf ( "\n" );
  printf ( "  ZCOPY ( 5, X, 2, A, 5 )\n" );
  printf ( "\n" );
  printf ( "  A =\n" );
  printf ( "\n" );
  for ( i = 0; i < N1; i++ )
    {
    for ( j = 0; j < N2; j++ )
    {
      printf ( "  %5f  %5f\n", creal ( a[i+j*N1] ), cimag ( a[i+j*N1] ) );
    }
    printf ( "\n" );
  }

  return;
# undef N
# undef N1
# undef N2
}
/******************************************************************************/

void test05 ( )

/******************************************************************************/
/*
  Purpose:

    TEST05 tests ZDOTC.

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
  double complex x_norm;
  double complex xy_dot;
  double complex x[N] = {
     2.0 - 1.0 * I, 
    -4.0 - 2.0 * I, 
     3.0 + 1.0 * I, 
     2.0 + 2.0 * I, 
    -1.0 - 1.0 * I };
  double complex y[N] = {
    -1.0 + 0.0 * I, 
     0.0 - 3.0 * I, 
     4.0 + 0.0 * I, 
    -3.0 + 4.0 * I, 
    -2.0 + 0.0 * I };

  printf ( "\n" );
  printf ( "TEST05\n" );
  printf ( "  ZDOTC computes the conjugated dot product of\n" );
  printf ( "  two complex vectors.\n" );
 
  printf ( "\n" );
  printf ( "  X =\n" );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    printf ( "  %6d  %6f  %6f\n", i, creal ( x[i] ), cimag ( x[i] ) );
  }

  x_norm = zdotc ( N, x, 1, x, 1 );

  printf ( "\n" );
  printf ( "  The square of the norm of X, computed as\n" );
  printf ( "  ZDOTC(X,X) = (%f,  %f)\n", creal ( x_norm ), cimag ( x_norm ) );

  xy_dot = zdotc ( N, x, 1, y, 1 );

  printf ( "\n" );
  printf ( "  Y =\n" );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    printf ( "  %6d  %6f  %6f\n", i, creal ( y[i] ), cimag ( y[i] ) );
  }
  printf ( "\n" );
  printf ( "  The dot product X.Y* is (%f,  %f)\n", 
    creal ( xy_dot ), cimag ( xy_dot ) );

  return;
# undef N
}
/******************************************************************************/

void test06 ( )

/******************************************************************************/
/*
  Purpose:

    TEST06 tests ZDOTU.

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
  double complex x_norm;
  double complex xy_dot;
  double complex x[N] = {
    2.0 - 1.0 * I, 
   -4.0 - 2.0 * I, 
    3.0 + 1.0 * I, 
    2.0 + 2.0 * I, 
   -1.0 - 1.0 * I };
  double complex y[N] = {
   -1.0 + 0.0 * I, 
    0.0 - 3.0 * I, 
    4.0 + 0.0 * I, 
   -3.0 + 4.0 * I, 
   -2.0 + 0.0 * I };

  printf ( "\n" );
  printf ( "TEST06\n" );
  printf ( "  ZDOTU computes the unconjugated dot product of\n" );
  printf ( "  two complex vectors.\n" );
 
  printf ( "\n" );
  printf ( "  X =\n" );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    printf ( "  %6d  %6f  %6f\n", i, creal ( x[i] ), cimag ( x[i] ) );
  }

  x_norm = zdotu ( N, x, 1, x, 1 );

  printf ( "\n" );
  printf ( "  The unconjugated dot product ( X dot X )\n" );
  printf ( "  (which is NOT the square of the norm of X!):\n" );
  printf ( "  ZDOTU(X,X) = (%f,  %f)\n", creal ( x_norm ), cimag ( x_norm ) );

  xy_dot = zdotu ( N, x, 1, y, 1 );

  printf ( "\n" );
  printf ( "  Y =\n" );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    printf ( "  %6d  %6f  %6f\n", i, creal ( y[i] ), cimag ( y[i] ) );
  }

  printf ( "\n" );
  printf ( "  The dot product ( X dot Y ) is (%f,  %f)\n",
    creal ( xy_dot ), cimag ( xy_dot ) );

  return;
# undef N
}
/******************************************************************************/

void test07 ( )

/******************************************************************************/
/*
  Purpose:

    TEST07 tests ZMACH.

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
  printf ( "  ZMACH computes several machine-dependent\n" );
  printf ( "  complex arithmetic parameters.\n" );

  printf ( "\n" );
  printf ( "  ZMACH(1)  = machine epsilon = %e\n", zmach ( 1 ) );
  printf ( "  ZMACH(2)  = a tiny value    = %e\n", zmach ( 2 ) );
  printf ( "  ZMACH(3)  = a huge value    = %e\n", zmach ( 3 ) );

  return;
}
/******************************************************************************/

void test08 ( )

/******************************************************************************/
/*
  Purpose:

    TEST08 tests ZSCAL.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    31 March 2007

  Author:

    John Burkardt
*/
{
# define N 6

  double complex da;
  int i;
  double complex x[N];

  for ( i = 0; i < N; i++ )
  {
    x[i] = 10 * ( i + 1 ) + ( i + 1 ) * I;
  }

  printf ( "\n" );
  printf ( "TEST08\n" );
  printf ( "  ZSCAL multiplies a complex scalar times a vector.\n" );

  printf ( "\n" );
  printf ( "  X =\n" );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    printf ( "  %6d  %6f  %6f\n", i, creal ( x[i] ), cimag ( x[i] ) );
  }

  da = 5.0;
  zscal ( N, da, x, 1 );
  printf ( "\n" );
  printf ( "  ZSCAL ( N, (%f, %f), X, 1 )\n", creal ( da ), cimag ( da ) );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    printf ( "  %6d  %6f  %6f\n", i, creal ( x[i] ), cimag ( x[i] ) );
  }

  for ( i = 0; i < N; i++ )
  {
    x[i] = 10 * ( i + 1 ) + ( i + 1 ) * I;
  }

  da = -2.0 + 1.0 * I;
  zscal ( 3, da, x, 2 );
  printf ( "\n" );
  printf ( "  ZSCAL ( 3, (%f, %f), X, 2 )\n", creal ( da ), cimag ( da ) );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    printf ( "  %6d  %6f  %6f\n", i, creal ( x[i] ), cimag ( x[i] ) );
  }

  return;
# undef N
}
/******************************************************************************/

void test09 ( )

/******************************************************************************/
/*
  Purpose:

    TEST09 tests ZSIGN1.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    31 March 2007

  Author:

    John Burkardt
*/
{
  double complex c1;
  double complex c2;
  double complex c3;
  int i;
  int seed = 123456789;

  printf ( "\n" );
  printf ( "TEST09\n" );
  printf ( "  ZSIGN1 ( C1, C2 ) transfers the sign of complex C2\n" );
  printf ( "  to the ZABS1 magnitude of C1.\n" );
  printf ( "\n" );
  printf ( "\n" );
  printf ( "           C1                    C2                    C3\n" );
  printf ( 
    "  --------------------  --------------------  --------------------\n" );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    c1 = 5.0 * c8_uniform_01 ( &seed );
    c2 = 5.0 * c8_uniform_01 ( &seed );
    c3 = zsign1 ( c1, c2 );

    printf ( "  (%9f  %9f)  (%9f  %9f)   (%9f  %9f)\n", 
      creal ( c1 ), cimag ( c1 ), 
      creal ( c2 ), cimag ( c2 ), 
      creal ( c3 ), cimag ( c3 ) );
  }

  return;
}
/******************************************************************************/

void test10 ( )

/******************************************************************************/
/*
  Purpose:

    TEST10 tests ZSIGN2.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    31 March 2007

  Author:

    John Burkardt
*/
{
  double complex c1;
  double complex c2;
  double complex c3;
  int i;
  int seed = 123456789;

  printf ( "\n" );
  printf ( "TEST10\n" );
  printf ( "  ZSIGN2 ( C1, C2 ) transfers the sign of complex C2\n" );
  printf ( "  to the ZABS2 magnitude of C1.\n" );
  printf ( 
    "           C1                    C2                    C3\n" );
  printf ( 
    "  --------------------  --------------------  --------------------\n" );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    c1 = 5.0 * c8_uniform_01 ( &seed );
    c2 = 5.0 * c8_uniform_01 ( &seed );
    c3 = zsign2 ( c1, c2 );

    printf ( "  (%9f  %9f)  (%9f  %9f)   (%9f  %9f)\n", 
      creal ( c1 ), cimag ( c1 ), 
      creal ( c2 ), cimag ( c2 ), 
      creal ( c3 ), cimag ( c3 ) );
  }

  return;
}
/******************************************************************************/

void test11 ( )

/******************************************************************************/
/*
  Purpose:

    TEST11 tests ZDROT.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    31 March 2007

  Author:

    John Burkardt
*/
{
# define N 6

  double c;
  int i;
  double s;
  double complex x[N];
  double complex y[N];

  for ( i = 0; i < N; i++ )
  {
    x[i] = 10 * ( i + 1 ) + ( i + 1 ) * I;
  }

  for ( i = 0; i < N; i++ )
  {
    y[i] = 20 * ( i + 1 ) + 2 * ( i + 1 ) * I;
  }

  printf ( "\n" );
  printf ( "TEST11\n" );
  printf ( "  ZDROT carries out a Givens rotation\n" );
  printf ( "  on a complex vector.\n" );
  printf ( "\n" );
  printf ( "  X and Y\n" );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    printf ( "  %6d  (%9f  %9f)  (%9f  %9f)\n", 
    i, creal ( x[i] ), cimag ( x[i] ), creal ( y[i] ), cimag ( y[i] ) );
  }

  c = 0.5;
  s = sqrt ( 1.0 - c * c );
  zdrot ( N, x, 1, y, 1, c, s );
  printf ( "\n" );
  printf ( "  ZDROT ( N, X, 1, Y, 1, %f, %f )\n", c, s );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    printf ( "  %6d  (%9f  %9f)  (%9f  %9f)\n", 
    i, creal ( x[i] ), cimag ( x[i] ), creal ( y[i] ), cimag ( y[i] ) );
  }

  return;
# undef N
}
/******************************************************************************/

void test12 ( )

/******************************************************************************/
/*
  Purpose:

    TEST12 tests ZDSCAL.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    31 March 2007

  Author:

    John Burkardt
*/
{
# define N 6

  double da;
  int i;
  double complex x[N];

  for ( i = 0; i < N; i++ )
  {
    x[i] = 10 * ( i + 1 ) + ( i + 1 ) * I;
  }

  printf ( "\n" );
  printf ( "TEST12\n" );
  printf ( "  ZDSCAL multiplies a real scalar times a complex vector.\n" );

  printf ( "\n" );
  printf ( "  X =\n" );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    printf ( "  %6d  %6f  %6f\n", i, creal ( x[i] ), cimag ( x[i] ) );
  }

  da = 5.0;
  zdscal ( N, da, x, 1 );
  printf ( "\n" );
  printf ( "  ZDSCAL ( N, %f, X, 1 )\n", da );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    printf ( "  %6d  %6f  %6f\n", i, creal ( x[i] ), cimag ( x[i] ) );
  }

  for ( i = 0; i < N; i++ )
  {
    x[i] = 10 * ( i + 1 ) + ( i + 1 ) * I;
  }

  da = -2.0;
  zdscal ( 3, da, x, 2 );
  printf ( "\n" );
  printf ( "  ZDSCAL ( 3, %f, X, 2 )\n", da );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    printf ( "  %6d  %6f  %6f\n", i, creal ( x[i] ), cimag ( x[i] ) );
  }

  return;
# undef N
}
/******************************************************************************/

void test13 ( )

/******************************************************************************/
/*
  Purpose:

    TEST13 tests ZROTG.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    31 March 2007

  Author:

    John Burkardt
*/
{
  double complex a;
  double complex b;
  double c;
  double complex r;
  double complex s;
  double complex sa;
  double complex sb;
  double complex t;
  int seed;
  int test;
  int test_num = 5;

  printf ( "\n" );
  printf ( "TEST13\n" );
  printf ( "  ZROTG generates a complex Givens rotation\n" );
  printf ( "    (  C  S ) * ( A ) = ( R )\n" );
  printf ( "    ( -S  C )   ( B )   ( 0 )\n" );
  printf ( "\n" );

  seed = 123456789;

  for ( test = 1; test <= test_num; test++ )
  {
    a = c8_uniform_01 ( &seed );
    b = c8_uniform_01 ( &seed );

    sa = a;
    sb = b;

    zrotg ( &sa, sb, &c, &s );

    r = sa;

    printf ( "\n" );
    printf ( "  A =  ( %f. %f )\n", creal ( a ), cimag ( a ) );
    printf ( "  B =  ( %f. %f )\n", creal ( b ), cimag ( b ) );
    printf ( "  C =    %f\n", c );
    printf ( "  S =  ( %f. %f )\n", creal ( s ), cimag ( s ) );
    printf ( "  R =  ( %f. %f )\n", creal ( r ), cimag ( r ) );
    t = c * a + s * b;
    printf ( "         C *A+S*B = ( %f. %f )\n", creal ( t ), cimag ( t ) );
    t = - ( ~s ) * a + c * b;
    printf ( "  -conjg(S)*A+C*B = ( %f. %f )\n", creal ( t ), cimag ( t ) );
  }

  return;
}
/******************************************************************************/

void test14 ( )

/******************************************************************************/
/*
  Purpose:

    TEST14 tests ZSWAP.

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
  double complex x[N];
  double complex y[N];

  for ( i = 0; i < N; i++ )
  {
    x[i] = 10 * ( i + 1 ) + ( i + 1 ) * I;
  }

  for ( i = 0; i < N; i++ )
  {
    y[i] = 20 * ( i + 1 ) + 2 * ( i + 1 ) * I;
  }

  printf ( "\n" );
  printf ( "TEST14\n" );
  printf ( "  ZSWAP swaps two complex vectors.\n" );
  printf ( "\n" );
  printf ( "  X and Y\n" );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    printf ( "  %6d  %6f  %6f    %6f  %6f\n", 
    i, creal ( x[i] ), cimag ( x[i] ), creal ( y[i] ), cimag ( y[i] ) );
  }

  zswap ( N, x, 1, y, 1 );
  printf ( "\n" );
  printf ( "  ZSWAP ( N, X, 1, Y, 1 )\n" );
  printf ( "\n" );
  printf ( "  X and Y\n" );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    printf ( "  %6d  %6f  %6f    %6f  %6f\n", 
    i, creal ( x[i] ), cimag ( x[i] ), creal ( y[i] ), cimag ( y[i] ) );
  }

  for ( i = 0; i < N; i++ )
  {
    x[i] = 10 * ( i + 1 ) + ( i + 1 ) * I;
  }

  for ( i = 0; i < N; i++ )
  {
    y[i] = 20 * ( i + 1 ) + 2 * ( i + 1 ) * I;
  }

  zswap ( 3, x, 2, y, 1 );
  printf ( "\n" );
  printf ( "  ZSWAP ( 3, X, 2, Y, 1 )\n" );
  printf ( "\n" );
  printf ( "  X and Y\n" );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    printf ( "  %6d  %6f  %6f    %6f  %6f\n", 
    i, creal ( x[i] ), cimag ( x[i] ), creal ( y[i] ), cimag ( y[i] ) );
  }

  return;
}
/******************************************************************************/

void test15 ( )

/******************************************************************************/
/*
  Purpose:

    TEST15 tests IZAMAX.

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
  int incx;
  double complex x[N] = {
     2.0 - 1.0 * I, 
    -4.0 - 2.0 * I, 
     3.0 + 1.0 * I, 
     2.0 + 2.0 * I, 
    -1.0 - 1.0 * I };

  printf ( "\n" );
  printf ( "TEST15\n" );
  printf ( "  IZAMAX returns the index of the entry of\n" );
  printf ( "  maximum magnitude in a complex vector.\n" );
 
  printf ( "\n" );
  printf ( "  The entries and ZABS1 magnitudes:\n" );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    printf ( "  %6d  %6f  %6f    %6f\n", 
    i, creal ( x[i] ), cimag ( x[i] ), zabs1 ( x[i] ) );
  }

  incx = 1;

  i = izamax ( N, x, incx );

  printf ( "\n" );
  printf ( "  The index of maximum magnitude = %d\n", i );
  printf ( "\n" );
  printf ( "  Note that this is a 1-based index.\n" );
  printf ( "  Note that the L1 norm is used.\n" );

  return;
# undef N
}
/******************************************************************************/

void test16 ( )

/******************************************************************************/
/*
  Purpose:

    TEST16 tests DZASUM.

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

  double complex a[MA*NA] = {
    -3.0 + 4.0 * I, 
     2.0 + 0.0 * I, 
     3.0 - 4.0 * I, 
     2.0 + 0.0 * I, 
     2.0 - 1.0 * I, 
    -1.0 + 1.0 * I, 
     0.0 + 5.0 * I, 
    -4.0 - 2.0 * I, 
    -4.0 + 1.0 * I, 
    -4.0 - 3.0 * I, 
     0.0 - 2.0 * I, 
     1.0 + 3.0 * I, 
    -3.0 + 3.0 * I, 
    -3.0 + 3.0 * I, 
    -1.0 - 2.0 * I, 
    -1.0 + 2.0 * I, 
     2.0 - 4.0 * I, 
     0.0 - 1.0 * I, 
     0.0 - 1.0 * I, 
    -2.0 + 4.0 * I };
  int i;
  int j;
  double complex x[NX] = {
    2.0 - 1.0 * I, 
   -4.0 - 2.0 * I, 
    3.0 + 1.0 * I, 
    2.0 + 2.0 * I, 
   -1.0 - 1.0 * I, 
   -1.0 + 0.0 * I, 
    0.0 - 3.0 * I, 
    4.0 + 0.0 * I };

  printf ( "\n" );
  printf ( "TEST16\n" );
  printf ( "  DZASUM adds the absolute values of elements\n" );
  printf ( "  of a complex vector.\n" );
  printf ( "\n" );
  printf ( "  X =\n" );
  printf ( "\n" );
  for ( i = 0; i < NX; i++ )
  {
    printf ( "  %6d  (%6.1f  %6.1f )\n", i, creal ( x[i] ), cimag ( x[i] ) );
  }

  printf ( "\n" );
  printf ( "  DZASUM ( NX,   X, 1    ) = %f\n",
    dzasum ( NX,   x, 1 ) );
  printf ( "  DZASUM ( NX/2, X, 2    ) = %f\n",
    dzasum ( NX/2, x, 2 ) );
  printf ( "  DZASUM ( 2,    X, NX/2 ) = %f\n",
    dzasum ( 2,    x, NX/2 ) );

  printf ( "\n" );
  printf ( "  Demonstrate with a matrix A:\n" );
  printf ( "\n" );
  for ( i = 0; i < MA; i++ )
  {
    for ( j = 0; j < NA; j++ )
    {
      printf ( "  (%6.1f  %6.1f )", 
        creal ( a[i+j*MA] ), cimag ( a[i+j*MA] ) );
    }
    printf ( "\n" );
  }

  printf ( "\n" );
  printf ( "  DZASUM ( MA, A[1,2], 1 )   = %f\n",
    dzasum ( MA, a+0+1*MA, 1 ) );
  printf ( "  DZASUM ( NA, A[2,1], MA ) = %f\n",
    dzasum ( NA, a+1+0*MA, MA ) );

  return;
# undef MA
# undef NA
# undef NX
}
/******************************************************************************/

void test17 ( )

/******************************************************************************/
/*
  Purpose:

    TEST17 tests DZNRM2.

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
  int incx;
  double norm;
  double complex x[N] = {
    2.0 - 1.0 * I, 
   -4.0 - 2.0 * I, 
    3.0 + 1.0 * I, 
    2.0 + 2.0 * I, 
   -1.0 - 1.0 * I };

  printf ( "\n" );
  printf ( "TEST17\n" );
  printf ( "  DZNRM2 returns the Euclidean norm of a complex vector.\n" );
 
  printf ( "\n" );
  printf ( "  The vector X:\n" );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    printf ( "  %6d  %6f  %6f\n", 
    i, creal ( x[i] ), cimag ( x[i] ) );
  }

  incx = 1;
  norm = dznrm2 ( N, x, incx );

  printf ( "\n" ); 
  printf ( "  The L2 norm of X is %f\n", norm );

  return;
# undef N
}

