# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <complex.h>

# include "blas0.h"
# include "blas1_d.h"

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

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for BLAS1_D_PRB.

  Discussion:

    BLAS1_D_PRB tests the BLAS1_D library.

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
  printf ( "BLAS1_D_PRB:\n" );
  printf ( "  C version\n" );
  printf ( "  Test the BLAS1_D library.\n" );

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
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "BLAS1_D_PRB:\n" );
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

    TEST01 tests DASUM.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 March 2007

  Author:

    John Burkardt
*/
{
# define LDA 6
# define MA 5
# define NA 4
# define NX 10

  double a[LDA*NA];
  int i;
  int j;
  double x[NX];

  for ( i = 0; i < NX; i++ )
  {
    x[i] = pow ( -1.0, i + 1 ) * ( double ) ( 2 * ( i + 1 ) );
  }

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  DASUM adds the absolute values of elements of a vector.\n" );
  printf ( "\n" );
  printf ( "  X = \n" );
  printf ( "\n" );
  for ( i = 0; i < NX; i++ )
  {
    printf ( "  %6d  %14f\n", i+1, x[i] );
  }

  printf ( "\n" );
  printf ( "  DASUM ( NX,   X, 1 ) =    %f\n", dasum ( NX,   x, 1 ) );
  printf ( "  DASUM ( NX/2, X, 2 ) =    %f\n", dasum ( NX/2, x, 2 ) );
  printf ( "  DASUM ( 2,    X, NX/2 ) = %f\n", dasum ( 2,    x, NX/2 ) );

  for ( i = 0; i < MA; i++ )
  {
    for ( j = 0; j < NA; j++ )
    {
      a[i+j*LDA] = pow ( -1.0, i + 1 + j + 1)
        * ( double ) ( 10 * ( i + 1 ) + j + 1 );
    }
  }

  printf ( "\n" );
  printf ( "  Demonstrate with a matrix A:\n" );
  printf ( "\n" );
  for ( i = 0; i < MA; i++ )
  {
    for ( j = 0; j < NA; j++ )
    {
      printf ( "  %14f", a[i+j*LDA] );
    }
    printf ( "\n" );
  }

  printf ( "\n" );
  printf ( "  DASUM(MA,A(1,2),1) =   %f\n", dasum ( MA, a+0+1*LDA, 1 ) );
  printf ( "  DASUM(NA,A(2,1),LDA) = %f\n", dasum ( NA, a+1+0*LDA, LDA ) );

  return;
# undef LDA
# undef MA
# undef NA
# undef NX
}
/******************************************************************************/

void test02 ( )

/******************************************************************************/
/*
  Purpose:

    TEST02 tests DAXPY.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 March 2007

  Author:

    John Burkardt
*/
{
# define N 6

  double da;
  int i;
  double x[N];
  double y[N];

  for ( i = 0; i < N; i++ )
  {
    x[i] = ( double ) ( i + 1 );
  }

  for ( i = 0; i < N; i++ )
  {
    y[i] = ( double ) ( 100 * ( i + 1 ) );
  }

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  DAXPY adds a multiple of vector X to vector Y.\n" );
  printf ( "\n" );
  printf ( "  X =\n" );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    printf ( "  %6d  %14f\n", i+1, x[i] );
  }
  printf ( "\n" );
  printf ( "  Y =\n" );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    printf ( "  %6d  %14f\n", i+1, y[i] );
  }

  da = 1.0;
  daxpy ( N, da, x, 1, y, 1 );
  printf ( "\n" );
  printf ( "  DAXPY ( N, %f, X, 1, Y, 1 )\n", da );
  printf ( "\n" );
  printf ( "  Y =\n" );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    printf ( "  %6d  %14f\n", i+1, y[i] );
  }

  for ( i = 0; i < N; i++ )
  {
    y[i] = ( double ) ( 100 * ( i + 1 ) );
  }

  da = -2.0;
  daxpy ( N, da, x, 1, y, 1 );
  printf ( "\n" );
  printf ( "  DAXPY ( N, %14f, X, 1, Y, 1 )\n", da );
  printf ( "\n" );
  printf ( "  Y =\n" );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    printf ( "  %6d  %14f\n", i+1, y[i] );
  }

  for ( i = 0; i < N; i++ )
  {
    y[i] = ( double ) ( 100 * ( i + 1 ) );
  }

  da = 3.0;
  daxpy ( 3, da, x, 2, y, 1 );
  printf ( "\n" );
  printf ( "  DAXPY ( 3, %f, X, 2, Y, 1 )\n", da );
  printf ( "\n" );
  printf ( "  Y =\n" );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    printf ( "  %6d  %14f\n", i+1, y[i] );
  }

  for ( i = 0; i < N; i++ )
  {
    y[i] = ( double ) ( 100 * ( i + 1 ) );
  }

  da = -4.0;
  daxpy ( 3, da, x, 1, y, 2 );
  printf ( "\n" );
  printf ( "  DAXPY ( 3, %f, X, 1, Y, 2 )\n", da );
  printf ( "\n" );
  printf ( "  Y =\n" );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    printf ( "  %6d  %14f\n", i+1, y[i] );
  }

  return;
# undef N
}
/******************************************************************************/

void test03 ( )

/******************************************************************************/
/*
  Purpose:

    TEST03 demonstrates DCOPY.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 March 2007

  Author:

    John Burkardt
*/
{
  double a[5*5];
  int i;
  int j;
  double x[10];
  double y[10];

  printf ( "\n" );
  printf ( "TEST03\n" );
  printf ( "  DCOPY copies one vector into another.\n" );
  printf ( "\n" );

  for ( i = 0; i < 10; i++ )
  {
    x[i] = ( double ) ( i + 1 );
  }

  for ( i = 0; i < 10; i++ )
  {
    y[i] = ( double ) ( 10 * ( i + 1 ) );
  }

  for ( i = 0; i < 5; i++ )
  {
    for ( j = 0; j < 5; j++ )
    {
      a[i+j*5] = ( double ) ( 10 * ( i + 1 ) + j + 1 );
    }
  }

  printf ( "\n" );
  printf ( "  X =\n" );
  printf ( "\n" );
  for ( i = 0; i < 10; i++ )
  {
    printf ( "  %6d  %14f\n", i+1, x[i] );
  }
  printf ( "\n" );
  printf ( "  Y =\n" );
  printf ( "\n" );
  for ( i = 0; i < 10; i++ )
  {
    printf ( "  %6d  %14f\n", i+1, y[i] );
  }
  printf ( "\n" );
  printf ( "  A =\n" );
  printf ( "\n" );
  for ( i = 0; i < 5; i++ )
  {
    for ( j = 0; j < 5; j++ )
    {
      printf ( "  %14f", a[i+j*5] );
    }
      printf ( "\n" );
  }

  dcopy ( 5, x, 1, y, 1 );
  printf ( "\n" );
  printf ( "  DCOPY ( 5, X, 1, Y, 1 )\n" );
  printf ( "\n" );
  for ( i = 0; i < 10; i++ )
  {
    printf ( "  %6d  %14f\n", i+1, y[i] );
  }

  for ( i = 0; i < 10; i++ )
  {
    y[i] = ( double ) ( 10 * ( i + 1 ) );
  }

  dcopy ( 3, x, 2, y, 3 );
  printf ( "\n" );
  printf ( "  DCOPY ( 3, X, 2, Y, 3 )\n" );
  printf ( "\n" );
  for ( i = 0; i < 10; i++ )
  {
    printf ( "  %6d  %14f\n", i+1, y[i] );
  }

  dcopy ( 5, x, 1, a, 1 );
  printf ( "\n" );
  printf ( "  DCOPY ( 5, X, 1, A, 1 )\n" );
  printf ( "\n" );
  printf ( "  A =\n" );
  printf ( "\n" );
  for ( i = 0; i < 5; i++ )
  {
    for ( j = 0; j < 5; j++ )
    {
      printf ( "  %14f", a[i+j*5] );
    }
      printf ( "\n" );
  }

  for ( i = 0; i < 5; i++ )
  {
    for ( j = 0; j < 5; j++ )
    {
      a[i+j*5] = ( double ) ( 10 * ( i + 1 ) + j + 1 );
    }
  }

  dcopy ( 5, x, 2, a, 5 );
  printf ( "\n" );
  printf ( "  DCOPY ( 5, X, 2, A, 5 )\n" );
  printf ( "\n" );
  printf ( "  A =\n" );
  printf ( "\n" );
  for ( i = 0; i < 5; i++ )
  {
    for ( j = 0; j < 5; j++ )
    {
      printf ( "  %14f", a[i+j*5] );
    }
      printf ( "\n" );
  }

  return;
}
/******************************************************************************/

void test04 ( )

/******************************************************************************/
/*
  Purpose:

    TEST04 demonstrates DDOT.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 March 2007

  Author:

    John Burkardt
*/
{
# define N 5
# define LDA 10
# define LDB 7
# define LDC 6

  double a[LDA*LDA];
  double b[LDB*LDB];
  double c[LDC*LDC];
  int i;
  int j;
  double sum1;
  double x[N];
  double y[N];

  printf ( "\n" );
  printf ( "TEST04\n" );
  printf ( "  DDOT computes the dot product of vectors.\n" );
  printf ( "\n" );

  for ( i = 0; i < N; i++ )
  {
    x[i] = ( double ) ( i + 1 );
  }

  for ( i = 0; i < N; i++ )
  {
    y[i] = - ( double ) ( i + 1 );
  }

  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      a[i+j*LDA] = ( double ) ( i + 1 + j + 1 );
    }
  }

  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      b[i+j*LDB] = ( double ) ( ( i + 1 ) - ( j + 1 ) );
    }
  }

  sum1 = ddot ( N, x, 1, y, 1 );

  printf ( "\n" );
  printf ( "  Dot product of X and Y is %f\n", sum1 );
/*
  To multiply a ROW of a matrix A times a vector X, we need to
  specify the increment between successive entries of the row of A:
*/
  sum1 = ddot ( N, a+1+0*LDA, LDA, x, 1 );

  printf ( "\n" );
  printf ( "  Product of row 2 of A and X is %f\n", sum1 );
/*
  Product of a column of A and a vector is simpler:
*/
  sum1 = ddot ( N, a+0+1*LDA, 1, x, 1 );

  printf ( "\n" );
  printf ( "  Product of column 2 of A and X is %f\n", sum1 );
/*
  Here's how matrix multiplication, c = a*b, could be done
  with DDOT:
*/
  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      c[i+j*LDC] = ddot ( N, a+i, LDA, b+0+j*LDB, 1 );
    }
  }

  printf ( "\n" );
  printf ( "  Matrix product computed with DDOT:\n" );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      printf ( "  %14f", c[i+j*LDC] );
    }
    printf ( "\n" );
  }

  return;
# undef N
# undef LDA
# undef LDB
# undef LDC
}
/******************************************************************************/

void test05 ( )

/******************************************************************************/
/*
  Purpose:

    TEST05 demonstrates DMACH.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 March 2007

  Author:

    John Burkardt
*/
{
  int job;

  printf ( "\n" );
  printf ( "TEST05\n" );
  printf ( "  DMACH returns some approximate machine numbers.\n" );
  printf ( "\n" );
  job = 1;
  printf ( "  DMACH(1) = EPS =  %e\n", dmach ( job ) );
  job = 2;
  printf ( "  DMACH(2) = TINY = %e\n", dmach ( job ) );
  job = 3;
  printf ( "  DMACH(3) = HUGE = %e\n", dmach ( job ) );

  return;
}
/******************************************************************************/

void test06 ( )

/******************************************************************************/
/*
  Purpose:

    TEST06 demonstrates DNRM2.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 March 2007

  Author:

    John Burkardt
*/
{
# define N 5
# define LDA 10
/*
  These parameters illustrate the fact that matrices are typically
  dimensioned with more space than the user requires.
*/
  double a[LDA*LDA];
  int i;
  int incx;
  int j;
  double sum1;
  double x[N];

  printf ( "\n" );
  printf ( "TEST06\n" );
  printf ( "  DNRM2 computes the Euclidean norm of a vector.\n" );
  printf ( "\n" );
/*
  Compute the euclidean norm of a vector:
*/
  for ( i = 0; i < N; i++ )
  {
    x[i] = ( double ) ( i + 1 );
  }

  printf ( "\n" );
  printf ( "  X =\n" );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    printf ( "  %6d  %14f\n", i+1, x[i] );
  }
  printf ( "\n" );
  printf ( "  The 2-norm of X is %f\n", dnrm2 ( N, x, 1 ) );
/*
  Compute the euclidean norm of a row or column of a matrix:
*/
  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      a[i+j*LDA] = ( double ) ( i + 1 + j + 1 );
    }
  }

  printf ( "\n" );
  printf ( "  The 2-norm of row 2 of A is %f\n",
    dnrm2 ( N, a+1+0*LDA, LDA ) );

  printf ( "\n" );
  printf ( "  The 2-norm of column 2 of A is %f\n",
    dnrm2 ( N, a+0+1*LDA, 1 ) );

  return;
# undef N
# undef LDA
}
/******************************************************************************/

void test07 ( )

/******************************************************************************/
/*
  Purpose:

    TEST07 tests DROT.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 March 2007

  Author:

    John Burkardt
*/
{
# define N 6

  double c;
  int i;
  double s;
  double x[N];
  double y[N];

  for ( i = 0; i < N; i++ )
  {
    x[i] = ( double ) ( i + 1 );
  }

  for ( i = 0; i < N; i++ )
  {
    y[i] = ( double ) ( ( i + 1 ) * ( i + 1 ) - 12 );
  }

  printf ( "\n" );
  printf ( "TEST07\n" );
  printf ( "  DROT carries out a Givens rotation.\n" );
  printf ( "\n" );
  printf ( "  X and Y\n" );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    printf ( "  %6d  %14f  %14f\n", i+1, x[i], y[i] );
  }

  c = 0.5;
  s = sqrt ( 1.0 - c * c );
  drot ( N, x, 1, y, 1, c, s );
  printf ( "\n" );
  printf ( "  DROT ( N, X, 1, Y, 1, %f, %f )\n", c, s );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    printf ( "  %6d  %14f  %14f\n", i+1, x[i], y[i] );
  }

  for ( i = 0; i < N; i++ )
  {
    x[i] = ( double ) ( i + 1 );
  }

  for ( i = 0; i < N; i++ )
  {
    y[i] = ( double ) ( ( i + 1 ) * ( i + 1 ) - 12 );
  }

  c = x[0] / sqrt ( x[0] * x[0] + y[0] * y[0] );
  s = y[0] / sqrt ( x[0] * x[0] + y[0] * y[0] );
  drot ( N, x, 1, y, 1, c, s );
  printf ( "\n" );
  printf ( "  DROT ( N, X, 1, Y, 1, %f, %f )\n", c, s );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    printf ( "  %6d  %14f  %14f\n", i+1, x[i], y[i] );
  }

  return;
# undef N
}
/******************************************************************************/

void test08 ( )

/******************************************************************************/
/*
  Purpose:

    TEST08 tests DROTG.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 March 2007

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double c;
  double r;
  double s;
  double sa;
  double sb;
  int seed;
  int test;
  int test_num = 5;
  double z;

  printf ( "\n" );
  printf ( "TEST08\n" );
  printf ( "  DROTG generates a real Givens rotation\n" );
  printf ( "    (  C  S ) * ( A ) = ( R )\n" );
  printf ( "    ( -S  C )   ( B )   ( 0 )\n" );
  printf ( "\n" );

  seed = 123456789;

  for ( test = 1; test <= test_num; test++ )
  {
    a = r8_uniform_01 ( &seed );
    b = r8_uniform_01 ( &seed );

    sa = a;
    sb = b;

    drotg ( &sa, &sb, &c, &s );

    r = sa;
    z = sb;

    printf ( "\n" );
    printf ( "  A =  %f,  B =  %f\n", a, b );
    printf ( "  C =  %f,  S =  %f\n" );
    printf ( "  R =  %f,  Z =  %f\n" );
    printf ( "   C*A+S*B = %f\n",  c * a + s * b );
    printf ( "  -S*A+C*B = %f\n", -s * a + c * b );
  }

  return;
}
/******************************************************************************/

void test09 ( )

/******************************************************************************/
/*
  Purpose:

    TEST09 tests DSCAL.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 March 2007

  Author:

    John Burkardt
*/
{
# define N 6

  double da;
  int i;
  double x[N];

  for ( i = 0; i < N; i++ )
  {
    x[i] = ( double ) ( i + 1 );
  }

  printf ( "\n" );
  printf ( "TEST09\n" );
  printf ( "  DSCAL multiplies a vector by a scalar.\n" );
  printf ( "\n" );
  printf ( "  X =\n" );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    printf ( "  %6d  %14f\n", i+1, x[i] );
  }

  da = 5.0;
  dscal ( N, da, x, 1 );
  printf ( "\n" );
  printf ( "  DSCAL ( N, %f, X, 1 )\n", da );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    printf ( "  %6d  %14f\n", i+1, x[i] );
  }

  for ( i = 0; i < N; i++ )
  {
    x[i] = ( double ) ( i + 1 );
  }

  da = -2.0;
  dscal ( 3, da, x, 2 );
  printf ( "\n" );
  printf ( "  DSCAL ( 3, %f, X, 2 )\n", da );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    printf ( "  %6d  %14f\n", i+1, x[i] );
  }

  return;
# undef N
}
/******************************************************************************/

void test10 ( )

/******************************************************************************/
/*
  Purpose:

    TEST10 tests DSWAP.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 March 2007

  Author:

    John Burkardt
*/
{
# define N 6

  int i;
  double x[N];
  double y[N];

  for ( i = 0; i < N; i++ )
  {
    x[i] = ( double ) ( i + 1 );
  }

  for ( i = 0; i < N; i++ )
  {
    y[i] = ( double ) ( 100 * ( i + 1 ) );
  }

  printf ( "\n" );
  printf ( "TEST10\n" );
  printf ( "  DSWAP swaps two vectors.\n" );
  printf ( "\n" );
  printf ( "  X and Y:\n" );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    printf ( "  %6d  %14f  %14f\n", i+1, x[i], y[i] );
  }

  dswap ( N, x, 1, y, 1 );
  printf ( "\n" );
  printf ( "  DSWAP ( N, X, 1, Y, 1 )\n" );
  printf ( "\n" );
  printf ( "  X and Y:\n" );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    printf ( "  %6d  %14f  %14f\n", i+1, x[i], y[i] );
  }

  for ( i = 0; i < N; i++ )
  {
    x[i] = ( double ) ( i + 1 );
  }

  for ( i = 0; i < N; i++ )
  {
    y[i] = ( double ) ( 100 * ( i + 1 ) );
  }

  dswap ( 3, x, 2, y, 1 );
  printf ( "\n" );
  printf ( "  DSWAP ( 3, X, 2, Y, 1 )\n" );

  printf ( "\n" );
  printf ( "  X and Y:\n" );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    printf ( "  %6d  %14f  %14f\n", i+1, x[i], y[i] );
  }

  return;
# undef N
}
/******************************************************************************/

void test11 ( )

/******************************************************************************/
/*
  Purpose:

    TEST11 demonstrates IDAMAX.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 March 2007

  Author:

    John Burkardt
*/
{
# define N 11

  int i;
  int i1;
  int incx;
  double x[N];

  printf ( "\n" );
  printf ( "TEST11\n" );
  printf ( "  IDAMAX returns the index of maximum magnitude;\n" );

  for ( i = 1; i <= N; i++ )
  {
    x[i-1] = ( double ) ( ( 7 * i ) % 11 ) - ( double ) ( N / 2 );
  }

  printf ( "\n" );
  printf ( "  The vector X:\n" );
  printf ( "\n" );

  for ( i = 0; i < N; i++ )
  {
    printf ( "  %6d  %14f\n", i+1, x[i] );
  }

  incx = 1;

  i1 = idamax ( N, x, incx );

  printf ( "\n" );
  printf ( "  The index of maximum magnitude = %d\n", i1 );

  return;
# undef N
}
/******************************************************************************/

void test12 ( )

/******************************************************************************/
/*
  Purpose:

    TEST12 tests IDAMAX, DAXPY and DSCAL.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 March 2007

  Author:

    John Burkardt
*/
{
# define N 10
# define LDA N

  double a[LDA*N];
  double b[N];
  int i;
  int info;
  int ipvt[N];
  int j;
  int k;
  int l;
  double t;

  printf ( "\n" );
  printf ( "TEST12\n" );
  printf ( "  Use IDAMAX, DAXPY and DSCAL\n" );
  printf ( "  in a Gauss elimination routine.\n" );
/*
  Set the matrix.
*/
  for ( i = 1; i <= N; i++ )
  {
    for ( j = 1; j <= N; j++ )
    {
      if ( i == j )
      {
        a[i-1+(j-1)*LDA] = 2.0;
      }
      else if ( i == j + 1 )
      {
        a[i-1+(j-1)*LDA] = -1.0;
      }
      else if ( i == j - 1 )
      {
        a[i-1+(j-1)*LDA] = -1.0;
      }
      else
      {
        a[i-1+(j-1)*LDA] = 0.0;
      }
    }
  }
/*
  Set the right hand side.
*/
  for ( i = 1; i <= N-1; i++ )
  {
    b[i-1] = 0.0;
  }
  b[N-1] = ( double ) ( N + 1 );

  info = 0;

  for ( k = 1; k <= N-1; k++ )
  {
    l = idamax ( N-k+1, a+(k-1)+(k-1)*LDA, 1 ) + k - 1;
    ipvt[k-1] = l;

    if ( a[l-1+(k-1)*LDA] == 0.0 )
    {
      info = k;
    }
    else
    {
      if ( l != k )
      {
        t = a[l-1+(k-1)*LDA];
        a[l-1+(k-1)*LDA] = a[k-1+(k-1)*LDA];
        a[k-1+(k-1)*LDA] = t;
      }

      t = -1.0 / a[k-1+(k-1)*LDA];
      dscal ( N-k, t, a+k+(k-1)*LDA, 1 );

      for ( j = k+1; j <= N; j++ )
      {
        t = a[l-1+(j-1)*LDA];
        if ( l != k )
        {
          a[l-1+(j-1)*LDA] = a[k-1+(j-1)*LDA];
          a[k-1+(j-1)*LDA] = t;
        }
        daxpy ( N-k, t, a+k+(k-1)*LDA, 1, a+k+(j-1)*LDA, 1 );
      }
    }
  }

  ipvt[N-1] = N;
  if ( a[N-1+(N-1)*LDA] == 0.0 )
  {
    info = N;
  }

  if ( info != 0 )
  {
    printf ( "\n" );
    printf ( "  The matrix is singular.\n" );
    return;
  }

  for ( k = 1; k <= N-1; k++ )
  {
    l = ipvt[k-1];
    t = b[l-1];
    if ( l != k )
    {
      b[l-1] = b[k-1];
      b[k-1] = t;
    }
    daxpy ( N-k, t, a+k+(k-1)*LDA, 1, b+k, 1 );
  }

  for ( k = N; 1 <= k; k-- )
  {
    b[k-1] = b[k-1] / a[k-1+(k-1)*LDA];
    t = -b[k-1];
    daxpy ( k-1, t, a+0+(k-1)*LDA, 1, b, 1 );
  }

  printf ( "\n" );
  printf ( "  First five entries of solution:\n" );
  printf ( "\n" );
  for ( i = 1; i <= 5; i++ )
  {
    printf ( "  %14f", b[i-1] );
  }
  printf ( "\n" );

  return;
# undef LDA
# undef N
}

