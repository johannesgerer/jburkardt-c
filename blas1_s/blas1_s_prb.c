# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>

# include "blas1_s.h"

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
float r4_uniform_01 ( int *seed );
void timestamp ( void );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for BLAS1_S_PRB.

  Discussion:

    BLAS1_S_PRB tests the BLAS1 single precision real routines.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    29 March 2007

  Author:

    John Burkardt
*/
{
  timestamp ( );

  printf ( "\n" );
  printf ( "BLAS1_S_PRB:\n" );
  printf ( "  C version\n" );
  printf ( "  Test the BLAS1_S library.\n" );

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
  printf ( "BLAS1_S_PRB:\n" );
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

    TEST01 demonstrates ISAMAX.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    29 March 2007

  Author:

    John Burkardt
*/
{
# define N 11

  int i;
  int i1;
  int incx;
  float x[N];

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  ISAMAX returns the index of maximum magnitude;\n" );

  for ( i = 1; i <= N; i++ )
  {
    x[i-1] = ( float ) ( ( 7 * i ) % 11 ) - ( float ) ( N / 2 );
  }

  printf ( "\n" );
  printf ( "  The vector X:\n" );
  printf ( "\n" );

  for ( i = 1; i <= N; i++ )
  {
    printf ( "  %6d  %8f\n", i, x[i-1] );
  }

  incx = 1;

  i1 = isamax ( N, x, incx );

  printf ( "\n" );
  printf ( "  The index of maximum magnitude = %d\n", i1 );

  return;
# undef N
}
/******************************************************************************/

void test02 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST02 tests ISAMAX, SAXPY and SSCAL.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    29 March 2007

  Author:

    John Burkardt
*/
{
# define N 10
# define LDA N

  float a[LDA*N];
  float b[N];
  int i;
  int info;
  int ipvt[N];
  int j;
  int k;
  int l;
  float t;

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  Use ISAMAX, SAXPY and SSCAL\n" );
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
  b[N-1] = ( float ) ( N + 1 );

  info = 0;

  for ( k = 1; k <= N-1; k++ )
  {
    l = isamax ( N-k+1, a+(k-1)+(k-1)*LDA, 1 ) + k - 1;
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
      sscal ( N-k, t, a+k+(k-1)*LDA, 1 );

      for ( j = k+1; j <= N; j++ )
      {
        t = a[l-1+(j-1)*LDA];
        if ( l != k )
        {
          a[l-1+(j-1)*LDA] = a[k-1+(j-1)*LDA];
          a[k-1+(j-1)*LDA] = t;
        }
        saxpy ( N-k, t, a+k+(k-1)*LDA, 1, a+k+(j-1)*LDA, 1 );
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
    saxpy ( N-k, t, a+k+(k-1)*LDA, 1, b+k, 1 );
  }

  for ( k = N; 1 <= k; k-- )
  {
    b[k-1] = b[k-1] / a[k-1+(k-1)*LDA];
    t = -b[k-1];
    saxpy ( k-1, t, a+0+(k-1)*LDA, 1, b, 1 );
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
/******************************************************************************/

void test03 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST03 tests SASUM.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    29 March 2007

  Author:

    John Burkardt
*/
{
# define LDA 6
# define MA 5
# define NA 4
# define NX 10

  float a[LDA*NA];
  int i;
  int j;
  float x[NX];

  for ( i = 0; i < NX; i++ )
  {
    x[i] = pow ( -1.0, i + 1 ) * ( float ) ( 2 * ( i + 1 ) );
  }

  printf ( "\n" );
  printf ( "TEST03\n" );
  printf ( "  SASUM adds the absolute values of elements of a vector.\n" );
  printf ( "\n" );
  printf ( "  X = \n" );
  printf ( "\n" );
  for ( i = 0; i < NX; i++ )
  {
    printf ( "  %6d  %14f\n", i + 1, x[i] );
  }

  printf ( "\n" );
  printf ( "  SASUM ( NX,   X, 1 ) =    %f\n", sasum ( NX,   x, 1 ) );
  printf ( "  SASUM ( NX/2, X, 2 ) =    %f\n", sasum ( NX/2, x, 2 ) );
  printf ( "  SASUM ( 2,    X, NX/2 ) = %f\n", sasum ( 2,    x, NX/2 ) );

  for ( i = 0; i < MA; i++ )
  {
    for ( j = 0; j < NA; j++ )
    {
      a[i+j*LDA] = pow ( -1.0, i + 1 + j + 1)
        * ( float ) ( 10 * ( i + 1 ) + j + 1 );
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
  printf ( "  SASUM(MA,A(1,2),1) =   %f\n", sasum ( MA, a+0+1*LDA, 1 ) );
  printf ( "  SASUM(NA,A(2,1),LDA) = %f\n", sasum ( NA, a+1+0*LDA, LDA ) );

  return;
# undef LDA
# undef MA
# undef NA
# undef NX
}
/******************************************************************************/

void test04 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST04 tests SAXPY.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    29 March 2007

  Author:

    John Burkardt
*/
{
# define N 6

  float da;
  int i;
  float x[N];
  float y[N];

  for ( i = 0; i < N; i++ )
  {
    x[i] = ( float ) ( i + 1 );
  }

  for ( i = 0; i < N; i++ )
  {
    y[i] = ( float ) ( 100 * ( i + 1 ) );
  }

  printf ( "\n" );
  printf ( "TEST04\n" );
  printf ( "  SAXPY adds a multiple of vector X to vector Y.\n" );
  printf ( "\n" );
  printf ( "  X =\n" );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    printf ( "  %6d  %14f\n", i + 1, x[i] );
  }
  printf ( "\n" );
  printf ( "  Y =\n" );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    printf ( "  %6d  %14f\n", i + 1, y[i] );
  }

  da = 1.0;
  saxpy ( N, da, x, 1, y, 1 );
  printf ( "\n" );
  printf ( "  SAXPY ( N, %f, X, 1, Y, 1 )\n", da );
  printf ( "\n" );
  printf ( "  Y =\n" );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    printf ( "  %6d  %14f\n", i + 1, y[i] );
  }

  for ( i = 0; i < N; i++ )
  {
    y[i] = ( float ) ( 100 * ( i + 1 ) );
  }

  da = -2.0;
  saxpy ( N, da, x, 1, y, 1 );
  printf ( "\n" );
  printf ( "  SAXPY ( N, %f, X, 1, Y, 1 )\n", da );
  printf ( "\n" );
  printf ( "  Y =\n" );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    printf ( "  %6d  %14f\n", i + 1, y[i] );
  }

  for ( i = 0; i < N; i++ )
  {
    y[i] = ( float ) ( 100 * ( i + 1 ) );
  }

  da = 3.0;
  saxpy ( 3, da, x, 2, y, 1 );
  printf ( "\n" );
  printf ( "  SAXPY ( 3, %14f, X, 2, Y, 1 )\n", da );
  printf ( "\n" );
  printf ( "  Y =\n" );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    printf ( "  %6d  %14f\n", i + 1, y[i] );
  }

  for ( i = 0; i < N; i++ )
  {
    y[i] = ( float ) ( 100 * ( i + 1 ) );
  }

  da = -4.0;
  saxpy ( 3, da, x, 1, y, 2 );
  printf ( "\n" );
  printf ( "  SAXPY ( 3, %14f, X, 1, Y, 2 )\n", da );
  printf ( "\n" );
  printf ( "  Y =\n" );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    printf ( "  %6d  %14f\n", i + 1, y[i] );
  }

  return;
# undef N
}
/******************************************************************************/

void test05 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST05 demonstrates SCOPY.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    29 March 2007

  Author:

    John Burkardt
*/
{
  float a[5*5];
  int i;
  int j;
  float x[10];
  float y[10];

  printf ( "\n" );
  printf ( "TEST05\n" );
  printf ( "  SCOPY copies one vector into another.\n" );
  printf ( "\n" );

  for ( i = 0; i < 10; i++ )
  {
    x[i] = ( float ) ( i + 1 );
  }

  for ( i = 0; i < 10; i++ )
  {
    y[i] = ( float ) ( 10 * ( i + 1 ) );
  }

  for ( i = 0; i < 5; i++ )
  {
    for ( j = 0; j < 5; j++ )
    {
      a[i+j*5] = ( float ) ( 10 * ( i + 1 ) + j + 1 );
    }
  }

  printf ( "\n" );
  printf ( "  X =\n" );
  printf ( "\n" );
  for ( i = 0; i < 10; i++ )
  {
    printf ( "  %6d  %14f\n", i + 1, x[i] );
  }
  printf ( "\n" );
  printf ( "  Y =\n" );
  printf ( "\n" );
  for ( i = 0; i < 10; i++ )
  {
    printf ( "  %6d  %14f\n", i + 1, y[i] );
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

  scopy ( 5, x, 1, y, 1 );
  printf ( "\n" );
  printf ( "  SCOPY ( 5, X, 1, Y, 1 )\n" );
  printf ( "\n" );
  for ( i = 0; i < 10; i++ )
  {
    printf ( "  %6d  %14f\n", i + 1, y[i] );
  }

  for ( i = 0; i < 10; i++ )
  {
    y[i] = ( float ) ( 10 * ( i + 1 ) );
  }

  scopy ( 3, x, 2, y, 3 );
  printf ( "\n" );
  printf ( "  SCOPY ( 3, X, 2, Y, 3 )\n" );
  printf ( "\n" );
  for ( i = 0; i < 10; i++ )
  {
    printf ( "  %6d  %14f\n", i + 1, y[i] );
  }

  scopy ( 5, x, 1, a, 1 );
  printf ( "\n" );
  printf ( "  SCOPY ( 5, X, 1, A, 1 )\n" );
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
      a[i+j*5] = ( float ) ( 10 * ( i + 1 ) + j + 1 );
    }
  }

  scopy ( 5, x, 2, a, 5 );
  printf ( "\n" );
  printf ( "  SCOPY ( 5, X, 2, A, 5 )\n" );
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

void test06 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST06 demonstrates SDOT.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    29 March 2007

  Author:

    John Burkardt
*/
{
# define N 5
# define LDA 10
# define LDB 7
# define LDC 6

  float a[LDA*LDA];
  float b[LDB*LDB];
  float c[LDC*LDC];
  int i;
  int j;
  float sum1;
  float x[N];
  float y[N];

  printf ( "\n" );
  printf ( "TEST06\n" );
  printf ( "  SDOT computes the dot product of vectors.\n" );
  printf ( "\n" );

  for ( i = 0; i < N; i++ )
  {
    x[i] = ( float ) ( i + 1 );
  }

  for ( i = 0; i < N; i++ )
  {
    y[i] = - ( float ) ( i + 1 );
  }

  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      a[i+j*LDA] = ( float ) ( i + 1 + j + 1 );
    }
  }

  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      b[i+j*LDB] = ( float ) ( ( i + 1 ) - ( j + 1 ) );
    }
  }

  sum1 = sdot ( N, x, 1, y, 1 );

  printf ( "\n" );
  printf ( "  Dot product of X and Y is %f\n", sum1 );
/*
  To multiply a ROW of a matrix A times a vector X, we need to
  specify the increment between successive entries of the row of A:
*/
  sum1 = sdot ( N, a+1+0*LDA, LDA, x, 1 );

  printf ( "\n" );
  printf ( "  Product of row 2 of A and X is %14f\n", sum1 );
/*
  Product of a column of A and a vector is simpler:
*/
  sum1 = sdot ( N, a+0+1*LDA, 1, x, 1 );

  printf ( "\n" );
  printf ( "  Product of column 2 of A and X is %14f\n", sum1 );
/*
  Here's how matrix multiplication, c = a*b, could be done
  with SDOT:
*/
  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      c[i+j*LDC] = sdot ( N, a+i, LDA, b+0+j*LDB, 1 );
    }
  }

  printf ( "\n" );
  printf ( "  Matrix product computed with SDOT:\n" );
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

void test07 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST07 demonstrates SMACH.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    29 March 2007

  Author:

    John Burkardt
*/
{
  int job;

  printf ( "\n" );
  printf ( "TEST07\n" );
  printf ( "  SMACH returns some approximate machine numbers.\n" );
  printf ( "\n" );
  job = 1;
  printf ( "  SMACH(1) = EPS =  %e\n", smach ( job ) );
  job = 2;
  printf ( "  SMACH(2) = TINY = %e\n", smach ( job ) );
  job = 3;
  printf ( "  SMACH(3) = HUGE = %e\n", smach ( job ) );

  return;
}
/******************************************************************************/

void test08 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST08 demonstrates SNRM2.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    29 March 2007

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
  float a[LDA*LDA];
  int i;
  int incx;
  int j;
  float sum1;
  float x[N];

  printf ( "\n" );
  printf ( "TEST08\n" );
  printf ( "  SNRM2 computes the Euclidean norm of a vector.\n" );
  printf ( "\n" );
/*
  Compute the euclidean norm of a vector:
*/
  for ( i = 0; i < N; i++ )
  {
    x[i] = ( float ) ( i + 1 );
  }

  printf ( "\n" );
  printf ( "  X =\n" );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    printf ( "  %6d  %14d\n", i + 1, x[i] );
  }
  printf ( "\n" );
  printf ( "  The 2-norm of X is %f\n", snrm2 ( N, x, 1 ) );
/*
  Compute the euclidean norm of a row or column of a matrix:
*/
  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      a[i+j*LDA] = ( float ) ( i + 1 + j + 1 );
    }
  }

  printf ( "\n" );
  printf ( "  The 2-norm of row 2 of A is %f\n", snrm2 ( N, a+1+0*LDA, LDA ) );

  printf ( "\n" );
  printf ( "  The 2-norm of column 2 of A is %f\n" ,
       snrm2 ( N, a+0+1*LDA, 1 ) );

  return;
# undef N
# undef LDA
}
/******************************************************************************/

void test09 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST09 tests SROT.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    29 March 2007

  Author:

    John Burkardt
*/
{
# define N 6

  float c;
  int i;
  float s;
  float x[N];
  float y[N];

  for ( i = 0; i < N; i++ )
  {
    x[i] = ( float ) ( i + 1 );
  }

  for ( i = 0; i < N; i++ )
  {
    y[i] = ( float ) ( ( i + 1 ) * ( i + 1 ) - 12 );
  }

  printf ( "\n" );
  printf ( "TEST09\n" );
  printf ( "  SROT carries out a Givens rotation.\n" );
  printf ( "\n" );
  printf ( "  X and Y\n" );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    printf ( "  %6d  %14f  %14f\n", i + 1, x[i], y[i] );
  }

  c = 0.5;
  s = sqrt ( 1.0 - c * c );
  srot ( N, x, 1, y, 1, c, s );
  printf ( "\n" );
  printf ( "  SROT ( N, X, 1, Y, 1, %f, %f )\n", c, s );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    printf ( "  %6d  %14f  %14f\n", i + 1, x[i], y[i] );
  }

  for ( i = 0; i < N; i++ )
  {
    x[i] = ( float ) ( i + 1 );
  }

  for ( i = 0; i < N; i++ )
  {
    y[i] = ( float ) ( ( i + 1 ) * ( i + 1 ) - 12 );
  }

  c = x[0] / sqrt ( x[0] * x[0] + y[0] * y[0] );
  s = y[0] / sqrt ( x[0] * x[0] + y[0] * y[0] );
  srot ( N, x, 1, y, 1, c, s );
  printf ( "\n" );
  printf ( "  SROT ( N, X, 1, Y, 1, %f, %f )\n", c, s );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    printf ( "  %6d  %14f  %14f\n", i + 1, x[i], y[i] );
  }

  return;
# undef N
}
/******************************************************************************/

void test10 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST10 tests SROTG.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    29 March 2007

  Author:

    John Burkardt
*/
{
  float a;
  float b;
  float c;
  float r;
  float s;
  float sa;
  float sb;
  int seed;
  int test;
  int test_num = 5;
  float z;

  printf ( "\n" );
  printf ( "TEST10\n" );
  printf ( "  SROTG generates a real Givens rotation\n" );
  printf ( "    (  C  S ) * ( A ) = ( R )\n" );
  printf ( "    ( -S  C )   ( B )   ( 0 )\n" );
  printf ( "\n" );

  seed = 123456789;

  for ( test = 1; test <= test_num; test++ )
  {
    a = r4_uniform_01 ( &seed );
    b = r4_uniform_01 ( &seed );

    sa = a;
    sb = b;

    srotg ( &sa, &sb, &c, &s );

    r = sa;
    z = sb;

    printf ( "\n" );
    printf ( "  A =  %14f,  B =  %14f\n", a, b );
    printf ( "  C =  %14f   S =  %14f\n", c, s );
    printf ( "  R =  %14f   Z =  %14f\n", r, z );
    printf ( "   C*A+S*B = %f\n",  c * a + s * b );
    printf ( "  -S*A+C*B = %f\n", -s * a + c * b );
  }

  return;
}
/******************************************************************************/

void test11 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST11 tests SSCAL.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    29 March 2007

  Author:

    John Burkardt
*/
{
# define N 6

  float da;
  int i;
  float x[N];

  for ( i = 0; i < N; i++ )
  {
    x[i] = ( float ) ( i + 1 );
  }

  printf ( "\n" );
  printf ( "TEST11\n" );
  printf ( "  SSCAL multiplies a vector by a scalar.\n" );
  printf ( "\n" );
  printf ( "  X =\n" );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    printf ( "  %6d  %14f\n", i + 1, x[i] );
  }

  da = 5.0;
  sscal ( N, da, x, 1 );
  printf ( "\n" );
  printf ( "  SSCAL ( N, %f, X, 1 )\n", da );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    printf ( "  %6d  %14f\n", i + 1, x[i] );
  }

  for ( i = 0; i < N; i++ )
  {
    x[i] = ( float ) ( i + 1 );
  }

  da = -2.0;
  sscal ( 3, da, x, 2 );
  printf ( "\n" );
  printf ( "  SSCAL ( 3, %f, X, 2 )\n", da );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    printf ( "  %6d  %14f\n", i + 1, x[i] );
  }

  return;
# undef N
}
/******************************************************************************/

void test12 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST12 tests SSWAP.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    29 March 2007

  Author:

    John Burkardt
*/
{
# define N 6

  int i;
  float x[N];
  float y[N];

  for ( i = 0; i < N; i++ )
  {
    x[i] = ( float ) ( i + 1 );
  }

  for ( i = 0; i < N; i++ )
  {
    y[i] = ( float ) ( 100 * ( i + 1 ) );
  }

  printf ( "\n" );
  printf ( "TEST12\n" );
  printf ( "  SSWAP swaps two vectors.\n" );
  printf ( "\n" );
  printf ( "  X and Y:\n" );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    printf ( "  %6d  %14f  %14f\n", i + 1, x[i], y[i] );
  }

  sswap ( N, x, 1, y, 1 );
  printf ( "\n" );
  printf ( "  SSWAP ( N, X, 1, Y, 1 )\n" );
  printf ( "\n" );
  printf ( "  X and Y:\n" );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    printf ( "  %6d  %14f  %14f\n", i + 1, x[i], y[i] );
  }

  for ( i = 0; i < N; i++ )
  {
    x[i] = ( float ) ( i + 1 );
  }

  for ( i = 0; i < N; i++ )
  {
    y[i] = ( float ) ( 100 * ( i + 1 ) );
  }

  sswap ( 3, x, 2, y, 1 );
  printf ( "\n" );
  printf ( "  SSWAP ( 3, X, 2, Y, 1 )\n" );

  printf ( "\n" );
  printf ( "  X and Y:\n" );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    printf ( "  %6d  %14f  %14f\n", i + 1, x[i], y[i] );
  }

  return;
# undef N
}
/******************************************************************************/

float r4_uniform_01 ( int *seed )

/******************************************************************************/
/*
  Purpose:

    R4_UNIFORM_01 returns a real pseudorandom number.

  Discussion:

    This routine implements the recursion

      seed = 16807 * seed mod ( 2**31 - 1 )
      r4_uniform_01 = seed / ( 2**31 - 1 )

    The integer arithmetic never requires more than 32 bits,
    including a sign bit.

    If the initial seed is 12345, then the first three computations are

      Input     Output      R_UNIFORM_01
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

    Philip Lewis, Allen Goodman, James Miller,
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
