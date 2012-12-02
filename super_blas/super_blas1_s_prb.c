# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>

# include "super_blas.h"

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
void test11 ( void );
float r4_uniform_01 ( int *seed );
float smach ( int job );
void timestamp ( void );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for SUPER_BLAS1_S_PRB.

  Discussion:

    SUPER_BLAS1_S_PRB tests the BLAS1 single precision real routines.

  Modified:

    29 March 2007

  Author:

    John Burkardt
*/
{
  timestamp ( );

  printf ( "\n" );
  printf ( "SUPER_BLAS1_S_PRB:\n" );
  printf ( "  C version\n" );
  printf ( "  Single precision real arithmetic.\n" );
  printf ( "  Test routines in the SUPER_BLAS library,\n" );
  printf ( "  the Level 1 Basic Linear Algebra Subprograms.\n" );
 
  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
  test05 ( );
  test06 ( );
  test07 ( );
  test08 ( );
  test09 ( );

  test11 ( );

  printf ( "\n" );
  printf ( "SUPER_BLAS1_S_PRB:\n" );
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

  Modified:

    29 March 2007

  Author:

    John Burkardt
*/
{
  int i;
  int i1;
  int incx;
  int n = 11;
  float *x;

  x = malloc ( n * sizeof ( float ) );

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  ISAMAX returns the index of maximum magnitude;\n" );
 
  for ( i = 1; i <= n; i++ )
  {
    x[i-1] = ( float ) ( ( 7 * i ) % 11 ) - ( float ) ( n / 2 );
  }

  printf ( "\n" );
  printf ( "  The vector X:\n" );
  printf ( "\n" );

  for ( i = 1; i <= n; i++ )
  {
    printf ( "  %6d  %8f\n", i, x[i-1] );
  }

  incx = 1;

  i1 = isamax_ ( &n, x, &incx );

  printf ( "\n" );
  printf ( "  The index of maximum magnitude = %d\n", i1 );

  free ( x );

  return;
}
/******************************************************************************/

void test02 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST02 tests ISAMAX, SAXPY and SSCAL.

  Modified:

    29 March 2007

  Author:

    John Burkardt
*/
{
  float *a;
  float *b;
  int i;
  int inc1;
  int inc2;
  int info;
  int *ipvt;
  int j;
  int k;
  int l;
  int lda;
  int n = 10;
  int ncopy;
  float t;

  lda = n;

  a = malloc ( lda * n * sizeof ( float ) );
  b = malloc ( n * sizeof ( float ) );
  ipvt = malloc ( n * sizeof ( int ) );

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  Use ISAMAX, SAXPY and SSCAL\n" );
  printf ( "  in a Gauss elimination routine.\n" );
/*
  Set the matrix.
*/
  for ( i = 1; i <= n; i++ )
  {
    for ( j = 1; j <= n; j++ )
    {
      if ( i == j )
      {
        a[i-1+(j-1)*lda] = 2.0;
      }
      else if ( i == j + 1 )
      {
        a[i-1+(j-1)*lda] = -1.0;
      }
      else if ( i == j - 1 )
      {
        a[i-1+(j-1)*lda] = -1.0;
      }
      else
      {
        a[i-1+(j-1)*lda] = 0.0;
      }
    }
  }
/*
  Set the right hand side.
*/
  for ( i = 1; i <= n-1; i++ )
  {
    b[i-1] = 0.0;
  }
  b[n-1] = ( float ) ( n + 1 );

  info = 0;

  for ( k = 1; k <= n-1; k++ )
  {
    ncopy = n - k + 1;
    inc1 = 1;

    l = isamax_ ( &ncopy, a+(k-1)+(k-1)*lda, &inc1 ) + k - 1;

    ipvt[k-1] = l;

    if ( a[l-1+(k-1)*lda] == 0.0 )
    {
      info = k;
    }
    else
    {
      if ( l != k )
      {
        t = a[l-1+(k-1)*lda];
        a[l-1+(k-1)*lda] = a[k-1+(k-1)*lda];
        a[k-1+(k-1)*lda] = t;
      }

      t = -1.0 / a[k-1+(k-1)*lda];

      ncopy = n - k;
      inc1 = 1;

      sscal_ ( &ncopy, &t, a+k+(k-1)*lda, &inc1 );

      for ( j = k+1; j <= n; j++ )
      {
        t = a[l-1+(j-1)*lda];
        if ( l != k )
        {
          a[l-1+(j-1)*lda] = a[k-1+(j-1)*lda];
          a[k-1+(j-1)*lda] = t;
        }

        ncopy = n - k;
        inc1 = 1;
        inc2 = 1;

        saxpy_ ( &ncopy, &t, a+k+(k-1)*lda, &inc1, a+k+(j-1)*lda, &inc2 );
      }
    }
  }

  ipvt[n-1] = n;
  if ( a[n-1+(n-1)*lda] == 0.0 )
  {
    info = n;
  }

  if ( info != 0 )
  {
    printf ( "\n" );
    printf ( "  The matrix is singular.\n" );
    return;
  }

  for ( k = 1; k <= n-1; k++ )
  {
    l = ipvt[k-1];
    t = b[l-1];
    if ( l != k )
    {
      b[l-1] = b[k-1];
      b[k-1] = t;
    }

    ncopy = n - k;
    inc1 = 1;
    inc2 = 1;

    saxpy_ ( &ncopy, &t, a+k+(k-1)*lda, &inc1, b+k, &inc2 );
  }

  for ( k = n; 1 <= k; k-- )
  {
    b[k-1] = b[k-1] / a[k-1+(k-1)*lda];
    t = -b[k-1];

    ncopy = k - 1;
    inc1 = 1;
    inc2 = 1;

    saxpy_ ( &ncopy, &t, a+0+(k-1)*lda, &inc1, b, &inc2 );
  }

  printf ( "\n" );
  printf ( "  First five entries of solution:\n" );
  printf ( "\n" );
  for ( i = 1; i <= 5; i++ )
  {
    printf ( "  %14f", b[i-1] );
  }
  printf ( "\n" );

  free ( a );
  free ( b );
  free ( ipvt );

  return;
}
/******************************************************************************/

void test03 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST03 tests SASUM.

  Modified:

    29 March 2007

  Author:

    John Burkardt
*/
{
  float *a;
  int i;
  int inc;
  int j;
  int lda = 6;
  int ma = 5;
  int na = 4;
  int nx = 10;
  int ncopy;
  float *x;

  a = malloc ( lda * na * sizeof ( float ) );
  x = malloc ( nx * sizeof ( float ) );

  for ( i = 0; i < nx; i++ )
  {
    x[i] = pow ( -1.0, i + 1 ) * ( float ) ( 2 * ( i + 1 ) );
  }

  printf ( "\n" );
  printf ( "TEST03\n" );
  printf ( "  SASUM adds the absolute values of elements of a vector.\n" );
  printf ( "\n" );
  printf ( "  X = \n" );
  printf ( "\n" );
  for ( i = 0; i < nx; i++ )
  {
    printf ( "  %6d  %14f\n", i + 1, x[i] );
  }

  printf ( "\n" );

  ncopy = nx;
  inc = 1;
  printf ( "  SASUM ( NX,   X, 1 ) =    %f\n", sasum_ ( &ncopy, x, &inc ) );

  ncopy = nx / 2;
  inc = 2;
  printf ( "  SASUM ( NX/2, X, 2 ) =    %f\n", sasum_ ( &ncopy, x, &inc ) );

  ncopy = 2;
  inc = nx / 2;
  printf ( "  SASUM ( 2,    X, NX/2 ) = %f\n", sasum_ ( &ncopy, x, &inc ) );

  for ( i = 0; i < ma; i++ )
  {
    for ( j = 0; j < na; j++ )
    {
      a[i+j*lda] = pow ( -1.0, i + 1 + j + 1) 
        * ( float ) ( 10 * ( i + 1 ) + j + 1 );
    }
  }

  printf ( "\n" );
  printf ( "  Demonstrate with a matrix A:\n" );
  printf ( "\n" );
  for ( i = 0; i < ma; i++ )
  {
    for ( j = 0; j < na; j++ )
    {
      printf ( "  %14f", a[i+j*lda] );
    }
    printf ( "\n" );
  }

  printf ( "\n" );

  ncopy = ma;
  inc = 1;
  printf ( "  SASUM(MA,A(1,2),1) =   %f\n", sasum_ ( &ncopy, a+0+1*lda, &inc ) );

  ncopy = na;
  inc = lda;
  printf ( "  SASUM(NA,A(2,1),LDA) = %f\n", sasum_ ( &ncopy, a+1+0*lda, &inc ) );

  free ( a );
  free ( x );

  return;
}
/******************************************************************************/

void test04 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST04 tests SAXPY.

  Modified:

    29 March 2007

  Author:

    John Burkardt
*/
{
  float da;
  int i;
  int inc1;
  int inc2;
  int n = 6;
  int ncopy;
  float *x;
  float *y;

  x = malloc ( n * sizeof ( float ) );
  y = malloc ( n * sizeof ( float ) );

  for ( i = 0; i < n; i++ )
  {
    x[i] = ( float ) ( i + 1 );
  }

  for ( i = 0; i < n; i++ )
  {
    y[i] = ( float ) ( 100 * ( i + 1 ) );
  }

  printf ( "\n" );
  printf ( "TEST04\n" );
  printf ( "  SAXPY adds a multiple of vector X to vector Y.\n" );
  printf ( "\n" );
  printf ( "  X =\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %6d  %14f\n", i + 1, x[i] );
  }
  printf ( "\n" );
  printf ( "  Y =\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %6d  %14f\n", i + 1, y[i] );
  }

  ncopy = n;
  da = 1.0;
  inc1 = 1;
  inc2 = 1;
  saxpy_ ( &ncopy, &da, x, &inc1, y, &inc2 );
  printf ( "\n" );
  printf ( "  SAXPY ( N, %f, X, 1, Y, 1 )\n", da );
  printf ( "\n" );
  printf ( "  Y =\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %6d  %14f\n", i + 1, y[i] );
  }

  for ( i = 0; i < n; i++ )
  {
    y[i] = ( float ) ( 100 * ( i + 1 ) );
  }

  ncopy = n;
  da = -2.0;
  inc1 = 1;
  inc2 = 1;

  saxpy_ ( &ncopy, &da, x, &inc1, y, &inc2 );
  printf ( "\n" );
  printf ( "  SAXPY ( N, %f, X, 1, Y, 1 )\n", da );
  printf ( "\n" );
  printf ( "  Y =\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %6d  %14f\n", i + 1, y[i] );
  }

  for ( i = 0; i < n; i++ )
  {
    y[i] = ( float ) ( 100 * ( i + 1 ) );
  }

  ncopy = 3;
  da = 3.0;
  inc1 = 2;
  inc2 = 1;

  saxpy_ ( &ncopy, &da, x, &inc1, y, &inc2 );

  printf ( "\n" );
  printf ( "  SAXPY ( 3, %14f, X, 2, Y, 1 )\n", da );
  printf ( "\n" );
  printf ( "  Y =\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %6d  %14f\n", i + 1, y[i] );
  }

  for ( i = 0; i < n; i++ )
  {
    y[i] = ( float ) ( 100 * ( i + 1 ) );
  }

  ncopy = 3;
  da = -4.0;
  inc1 = 1;
  inc2 = 2;

  saxpy_ ( &ncopy, &da, x, &inc1, y, &inc2 );

  printf ( "\n" );
  printf ( "  SAXPY ( 3, %14f, X, 1, Y, 2 )\n", da );
  printf ( "\n" );
  printf ( "  Y =\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %6d  %14f\n", i + 1, y[i] );
  }

  free ( x );
  free ( y );

  return;
}
/******************************************************************************/

void test05 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST05 demonstrates SCOPY.

  Modified:

    29 March 2007

  Author:

    John Burkardt
*/
{
  float a[5*5];
  int i;
  int inc1;
  int inc2;
  int j;
  int ncopy;
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

  ncopy = 5;
  inc1 = 1;
  inc2 = 1;

  scopy_ ( &ncopy, x, &inc1, y, &inc2 );

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

  ncopy = 3;
  inc1 = 2;
  inc2 = 3;

  scopy_ ( &ncopy, x, &inc1, y, &inc2 );

  printf ( "\n" );
  printf ( "  SCOPY ( 3, X, 2, Y, 3 )\n" );
  printf ( "\n" );
  for ( i = 0; i < 10; i++ )
  {
    printf ( "  %6d  %14f\n", i + 1, y[i] );
  }

  ncopy = 5;
  inc1 = 1;
  inc2 = 1;

  scopy_ ( &ncopy, x, &inc1, a, &inc2 );

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

  ncopy = 5;
  inc1 = 2;
  inc2 = 5;

  scopy_ ( &ncopy, x, &inc1, a, &inc2 );

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

  Modified:

    29 March 2007

  Author:

    John Burkardt
*/
{
  float *a;
  float *b;
  float *c;
  int i;
  int inc1;
  int inc2;
  int j;
  int lda = 10;
  int ldb = 7;
  int ldc = 6;
  int n = 5;
  int ncopy;
  float sum1;
  float *x;
  float *y;

  a = malloc ( lda * lda * sizeof ( float ) );
  b = malloc ( ldb * ldb * sizeof ( float ) );
  c = malloc ( ldc * ldc * sizeof ( float ) );
  x = malloc ( n * sizeof ( float ) );
  y = malloc ( n * sizeof ( float ) );

  printf ( "\n" );
  printf ( "TEST06\n" );
  printf ( "  SDOT computes the dot product of vectors.\n" );
  printf ( "\n" );
 
  for ( i = 0; i < n; i++ )
  {
    x[i] = ( float ) ( i + 1 );
  }

  for ( i = 0; i < n; i++ )
  {
    y[i] = - ( float ) ( i + 1 );
  }

  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      a[i+j*lda] = ( float ) ( i + 1 + j + 1 );
    }
  }

  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      b[i+j*ldb] = ( float ) ( ( i + 1 ) - ( j + 1 ) );
    }
  }

  ncopy = n;
  inc1 = 1;
  inc2 = 1;

  sum1 = sdot_ ( &ncopy, x, &inc1, y, &inc2 );

  printf ( "\n" );
  printf ( "  Dot product of X and Y is %f\n", sum1 );
/*
  To multiply a ROW of a matrix A times a vector X, we need to
  specify the increment between successive entries of the row of A:
*/
  ncopy = n;
  inc1 = lda;
  inc2 = 1;

  sum1 = sdot_ ( &ncopy, a+1+0*lda, &inc1, x, &inc2 );

  printf ( "\n" );
  printf ( "  Product of row 2 of A and X is %14f\n", sum1 );
/*
  Product of a column of A and a vector is simpler:
*/
  ncopy = n;
  inc1 = 1;
  inc2 = 1;

  sum1 = sdot_ ( &ncopy, a+0+1*lda, &inc1, x, &inc2 );

  printf ( "\n" );
  printf ( "  Product of column 2 of A and X is %14f\n", sum1 );
/*
  Here's how matrix multiplication, c = a*b, could be done
  with SDOT:
*/
  ncopy = n;
  inc1 = lda;
  inc2 = 1;

  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      c[i+j*ldc] = sdot_ ( &ncopy, a+i, &inc1, b+0+j*ldb, &inc2 );
    }
  }
 
  printf ( "\n" );
  printf ( "  Matrix product computed with SDOT:\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      printf ( "  %14f", c[i+j*ldc] );
    }
    printf ( "\n" );
  }

  free ( a );
  free ( b );
  free ( c );
  free ( x );
  free ( y );
 
  return;
}
/******************************************************************************/

void test07 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST07 demonstrates SMACH.

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

  Modified:

    29 March 2007

  Author:

    John Burkardt
*/
{
/*
  These parameters illustrate the fact that matrices are typically
  dimensioned with more space than the user requires.
*/
  float *a;
  int i;
  int inc;
  int j;
  int lda = 10;
  int n = 5;
  int ncopy;
  float sum1;
  float *x;

  a = malloc ( lda * lda * sizeof ( float ) );
  x = malloc ( n * sizeof ( float ) );

  printf ( "\n" );
  printf ( "TEST08\n" );
  printf ( "  SNRM2 computes the Euclidean norm of a vector.\n" );
  printf ( "\n" );
/*
  Compute the euclidean norm of a vector:
*/
  for ( i = 0; i < n; i++ )
  {
    x[i] = ( float ) ( i + 1 );
  }

  printf ( "\n" );
  printf ( "  X =\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %6d  %14d\n", i + 1, x[i] );
  }
  printf ( "\n" );
  ncopy = n;
  inc = 1;
  printf ( "  The 2-norm of X is %f\n", snrm2_ ( &ncopy, x, &inc ) );
/*
  Compute the euclidean norm of a row or column of a matrix:
*/
  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      a[i+j*lda] = ( float ) ( i + 1 + j + 1 );
    }
  }

  printf ( "\n" );
  ncopy = n;
  inc = lda;
  printf ( "  The 2-norm of row 2 of A is %f\n", 
    snrm2_ ( &ncopy, a+1+0*lda, &inc ) );

  printf ( "\n" );
  ncopy = n;
  inc = 1;
  printf ( "  The 2-norm of column 2 of A is %f\n" ,
       snrm2_ ( &ncopy, a+0+1*lda, &inc ) );

  free ( a );
  free ( x );
 
  return;
}
/******************************************************************************/

void test09 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST09 tests SROT.

  Modified:

    29 March 2007

  Author:

    John Burkardt
*/
{
  float c;
  int i;
  int inc1;
  int inc2;
  int n = 6;
  int ncopy;
  float s;
  float *x;
  float *y;

  x = malloc ( n * sizeof ( float ) );
  y = malloc ( n * sizeof ( float ) );

  for ( i = 0; i < n; i++ )
  {
    x[i] = ( float ) ( i + 1 );
  }

  for ( i = 0; i < n; i++ )
  {
    y[i] = ( float ) ( ( i + 1 ) * ( i + 1 ) - 12 );
  }

  printf ( "\n" );
  printf ( "TEST09\n" );
  printf ( "  SROT carries out a Givens rotation.\n" );
  printf ( "\n" );
  printf ( "  X and Y\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %6d  %14f  %14f\n", i + 1, x[i], y[i] );
  }

  ncopy = n;
  inc1 = 1;
  inc2 = 1;
  c = 0.5;
  s = sqrt ( 1.0 - c * c );

  srot_ ( &ncopy, x, &inc1, y, &inc2, &c, &s );

  printf ( "\n" );
  printf ( "  SROT ( N, X, 1, Y, 1, %f, %f )\n", c, s );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %6d  %14f  %14f\n", i + 1, x[i], y[i] );
  }

  for ( i = 0; i < n; i++ )
  {
    x[i] = ( float ) ( i + 1 );
  }

  for ( i = 0; i < n; i++ )
  {
    y[i] = ( float ) ( ( i + 1 ) * ( i + 1 ) - 12 );
  }

  c = x[0] / sqrt ( x[0] * x[0] + y[0] * y[0] );
  s = y[0] / sqrt ( x[0] * x[0] + y[0] * y[0] );

  ncopy = n;
  inc1 = 1;
  inc2 = 1;

  srot_ ( &ncopy, x, &inc1, y, &inc2, &c, &s );

  printf ( "\n" );
  printf ( "  SROT ( N, X, 1, Y, 1, %f, %f )\n", c, s );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %6d  %14f  %14f\n", i + 1, x[i], y[i] );
  }

  free ( x );
  free ( y );

  return;
}
/******************************************************************************/

void test11 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST11 tests SSCAL.

  Modified:

    29 March 2007

  Author:

    John Burkardt
*/
{
  float da;
  int i;
  int inc;
  int n = 6;
  int ncopy;
  float *x;

  x = malloc ( n * sizeof ( float ) );

  for ( i = 0; i < n; i++ )
  {
    x[i] = ( float ) ( i + 1 );
  }

  printf ( "\n" );
  printf ( "TEST11\n" );
  printf ( "  SSCAL multiplies a vector by a scalar.\n" );
  printf ( "\n" );
  printf ( "  X =\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %6d  %14f\n", i + 1, x[i] );
  }

  ncopy = n;
  da = 5.0;
  inc = 1;

  sscal_ ( &ncopy, &da, x, &inc );

  printf ( "\n" );
  printf ( "  SSCAL ( N, %f, X, 1 )\n", da );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %6d  %14f\n", i + 1, x[i] );
  }

  for ( i = 0; i < n; i++ )
  {
    x[i] = ( float ) ( i + 1 );
  }

  ncopy = 3;
  da = -2.0;
  inc = 2;

  sscal_ ( &ncopy, &da, x, &inc );

  printf ( "\n" );
  printf ( "  SSCAL ( 3, %f, X, 2 )\n", da );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %6d  %14f\n", i + 1, x[i] );
  }

  free ( x );

  return;
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

float smach ( int job )

/******************************************************************************/
/*
  Purpose:

    SMACH computes machine parameters of single precision real arithmetic.

  Discussion:

    This routine is for testing only.  It is not required by LINPACK.

    If there is trouble with the automatic computation of these quantities,
    they can be set by direct assignment statements.

    We assume the computer has

      B = base of arithmetic;
      T = number of base B digits;
      L = smallest possible exponent;
      U = largest possible exponent.

    then

      EPS = B**(1-T)
      TINY = 100.0 * B**(-L+T)
      HUGE = 0.01 * B**(U-T)

  Modified:

    29 March 2007

  Author:

    C version by John Burkardt

  Reference:

    Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart,
    LINPACK User's Guide,
    SIAM, 1979.

    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
    Basic Linear Algebra Subprograms for Fortran Usage,
    Algorithm 539, 
    ACM Transactions on Mathematical Software, 
    Volume 5, Number 3, September 1979, pages 308-323.

  Parameters:

    Input, int JOB:
    1: requests EPS;
    2: requests TINY;
    3: requests HUGE.

    Output, float SMACH, the requested value.
*/
{
  float eps;
  float huge;
  float s;
  float tiny;
  float value;

  eps = 1.0;
  for ( ; ; )
  {
    value = 1.0 + ( eps / 2.0 );
    if ( value <= 1.0 )
    {
      break;
    }
    eps = eps / 2.0;
  }

  s = 1.0;

  for ( ; ; )
  {
    tiny = s;
    s = s / 16.0;

    if ( s * 1.0 == 0.0 )
    {
      break;
    }

  }

  tiny = ( tiny / eps ) * 100.0;
  huge = 1.0 / tiny;

  if ( job == 1 )
  {
    value = eps;
  }
  else if ( job == 2 )
  {
    value = tiny;
  }
  else if ( job == 3 )
  {
    value = huge;
  }
  else
  {
    printf ( "\n" );
    printf ( "SMACH - Fatal error!\n" );
    printf ( "  Illegal input value of JOB = %d\n", job );
    exit ( 1 );
  }

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
