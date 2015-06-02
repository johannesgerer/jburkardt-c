# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <complex.h>

# include "blas0.h"
# include "blas2.h"

int main ( );
void test01 ( );
void test05 ( );
void test06 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for BLAS2_PRB.

  Discussion:

    BLAS2_PRB tests the BLAS library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    04 April 2014

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "BLAS2_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the BLAS library.\n" );

  test01 ( );
  test05 ( );
  test06 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "BLAS2_PRB\n" );
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

    TEST01 tests DGEMV.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 April 2014

  Author:

    John Burkardt
*/
{
  double *a;
  double alpha;
  double beta;
  int i;
  int incx;
  int incy;
  int j;
  int lda;
  int m;
  int n;
  char trans;
  double *x;
  double *y;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  For a general matrix A,\n" );
  printf ( "  DGEMV computes y := alpha * A * x + beta * y\n" );
  printf ( "  or             y := alpha * A'' * x + beta * y.\n" );
/*
  y = alpha * A * x + beta * y
*/
  trans = 'N';
  m = 5;
  n = 4;
  alpha = 2.0;
  lda = m;
  a = r8mat_test ( trans, lda, m, n );
  x = ( double * ) malloc ( n * sizeof ( double ) );
  for ( i = 0; i < n; i++ )
  {
    x[i] = ( double ) ( i + 1 );
  }
  incx = 1;
  beta = 3.0;
  y = ( double * ) malloc ( m * sizeof ( double ) );
  for ( i = 0; i < m; i++ )
  {
    y[i] = ( double ) ( 10 * ( i + 1 ) );
  }
  incy = 1;

  r8mat_print ( m, n, a, "  Matrix A:" );
  r8vec_print ( n, x, "  Vector X:" );
  r8vec_print ( m, y, "  Vector Y:" );

  dgemv ( trans, m, n, alpha, a, lda, x, incx, beta, y, incy );

  r8vec_print ( m, y, "  Result Y = alpha * A  * x + beta * y" );

  free ( a );
  free ( x );
  free ( y );
/*
  y = alpha * A' * x + beta * y
*/
  trans = 'T';
  m = 5;
  n = 4;
  alpha = 2.0;
  lda = m;
  a = r8mat_test ( trans, lda, n, m );
  x = ( double * ) malloc ( m * sizeof ( double ) );
  for ( i = 0; i < m; i++ )
  {
    x[i] = ( double ) ( i + 1 );
  }
  incx = 1;
  beta = 3.0;
  y = ( double * ) malloc ( n * sizeof ( double ) );
  for ( i = 0; i < n; i++ )
  {
    y[i] = ( double ) ( 10 * ( i + 1 ) );
  }
  incy = 1;

  r8mat_print ( m, n, a, "  Matrix A:" );
  r8vec_print ( m, x, "  Vector X:" );
  r8vec_print ( n, y, "  Vector Y:" );

  dgemv ( trans, m, n, alpha, a, lda, x, incx, beta, y, incy );

  r8vec_print ( n, y, "  Result Y = alpha * A  * x + beta * y" );

  free ( a );
  free ( x );
  free ( y );

  return;
}
/******************************************************************************/

void test05 ( )

/******************************************************************************/
/*
  Purpose:

    TEST05 tests DGER.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    04 April 2014

  Author:

    John Burkardt
*/
{
  double *a;
  double alpha;
  int i;
  int incx;
  int incy;
  int lda;
  int m;
  int n;
  char trans;
  double *x;
  double *y;

  printf ( "\n" );
  printf ( "TEST05\n" );
  printf ( "  For a general matrix A,\n" );
  printf ( "  DGER computes A := A + alpha * x * y'\n" );

  m = 5;
  n = 4;
  alpha = 2.0;
  trans = 'N';
  lda = m;
  a = r8mat_test ( trans, lda, m, n );

  x = ( double * ) malloc ( m * sizeof ( double ) );
  for ( i = 0; i < m; i++ )
  {
    x[i] = ( double ) ( i + 1 );
  }
  incx = 1;

  y = ( double * ) malloc ( n * sizeof ( double ) );
  for ( i = 0; i < n; i++ )
  {
    y[i] = ( double ) ( 10 * ( i + 1 ) );
  }
  incy = 1;

  r8mat_print ( m, n, a, "  Matrix A:" );
  r8vec_print ( m, x, "  Vector X:" );
  r8vec_print ( n, y, "  Vector Y:" );

  dger ( m, n, alpha, x, incx, y, incy, a, lda );

  r8mat_print ( m, n, a, "  Result A = A + alpha * x * y" );

  free ( a );
  free ( x );
  free ( y );

  return;
}
/******************************************************************************/

void test06 ( )

/******************************************************************************/
/*
  Purpose:

    TEST06 tests DTRMV.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    04 April 2014

  Author:

    John Burkardt
*/
{
  double *a;
  char diag;
  int i;
  int incx;
  int j;
  int lda = 5;
  int m = 5;
  int n = 5;
  int test;
  char trans;
  char uplo;
  double *x;

  a = ( double * ) malloc ( lda * n * sizeof ( double ) );
  x = ( double * ) malloc ( n * sizeof ( double ) );

  printf ( "\n" );
  printf ( "TEST06\n" );
  printf ( "  For a triangular matrix A,\n" );
  printf ( "  DTRMV computes y := A * x or y := A' * x\n" );

  for ( test = 1; test <= 2; test++ )
  {
    uplo = 'U';

    if ( test == 1 )
    {
      trans = 'N';
    }
    else
    {
      trans = 'T';
    }

    diag = 'N';

    for ( j = 0; j < n; j++ )
    {
      for ( i = 0; i <= j; i++ )
      {
        a[i+j*m] = ( double ) ( i + j + 2 );
      }
      for ( i = j + 1; i < m; i++ )
      {
        a[i+j*m] = 0.0;
      }
    }

    incx = 1;
    for ( i = 0; i < n; i++ )
    {
      x[i] = ( double ) ( i + 1 );
    }

    dtrmv ( uplo, trans, diag, n, a, lda, x, incx );

    if ( trans == 'N' )
    {
      r8vec_print ( n, x, "  Result y = A * x" );
    }
    else
    {
      r8vec_print ( n, x, "  Result y = A' * x" );
    }
  }

  free ( a );
  free ( x );

  return;
}