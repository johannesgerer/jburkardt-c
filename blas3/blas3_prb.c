# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <complex.h>

# include "blas0.h"
# include "blas3.h"

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for BLAS3_PRB.

  Discussion:

    BLAS3_PRB tests the BLAS3 library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 April 2014

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "BLAS3_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the BLAS3 library.\n" );

  test01 ( );
  test02 ( );
  test03 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "BLAS3_PRB\n" );
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

    TEST01 tests DGEMM.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    10 February 2014

  Author:

    John Burkardt
*/
{
  double *a;
  double alpha;
  double *b;
  double beta;
  double *c;
  int k;
  int lda;
  int ldb;
  int ldc;
  int m;
  int n;
  char transa;
  char transb;
  char transc;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  DGEMM carries out matrix multiplications\n" );
  printf ( "  for double precision real matrices.\n" );
  printf ( "\n" );
  printf ( "  1: C = alpha * A  * B  + beta * C;\n" );
  printf ( "  2: C = alpha * A' * B  + beta * C;\n" );
  printf ( "  3: C = alpha * A  * B' + beta * C;\n" );
  printf ( "  4: C = alpha * A' * B' + beta * C;\n" );
  printf ( "\n" );
  printf ( "  We carry out all four calculations, but in each case,\n" );
  printf ( "  we choose our input matrices so that we get the same result.\n" );
/*
  C = alpha * A * B + beta * C.
*/
  transa = 'N';
  transb = 'N';
  transc = 'N';
  m = 4;
  n = 5;
  k = 3;
  alpha = 2.0;
  lda = m;
  a = r8mat_test ( transa, lda, m, k );
  ldb = k;
  b = r8mat_test ( transb, ldb, k, n );
  beta = 3.0;
  ldc = m;
  c = r8mat_test ( transc, ldc, m, n );

  dgemm ( transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc );

  r8mat_print ( m, n, c, "  C = alpha * A * B + beta * C:" );

  free ( a );
  free ( b );
  free ( c );
/*
  C = alpha * A' * B + beta * C.
*/
  transa = 'T';
  transb = 'N';
  transc = 'N';
  m = 4;
  n = 5;
  k = 3;
  alpha = 2.0;
  lda = k;
  a = r8mat_test ( transa, lda, m, k );
  ldb = k;
  b = r8mat_test ( transb, ldb, k, n );
  beta = 3.0;
  ldc = m;
  c = r8mat_test ( transc, ldc, m, n );

  dgemm ( transa, transb,  m, n, k, alpha, a, lda, b, ldb, beta, c, ldc );

  r8mat_print ( m, n, c, "  C = alpha * A' * B + beta * C:" );

  free ( a );
  free ( b );
  free ( c );
/*
  C = alpha * A * B' + beta * C.
*/
  transa = 'N';
  transb = 'T';
  transc = 'N';
  m = 4;
  n = 5;
  k = 3;
  alpha = 2.0;
  lda = m;
  a = r8mat_test ( transa, lda, m, k );
  ldb = n;
  b = r8mat_test ( transb, ldb, k, n );
  beta = 3.0;
  ldc = m;
  c = r8mat_test ( transc, ldc, m, n );

  dgemm ( transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc );

  r8mat_print ( m, n, c, "  C = alpha * A * B' + beta * C:" );

  free ( a );
  free ( b );
  free ( c );
/*
  C = alpha * A' * B' + beta * C.
*/
  transa = 'T';
  transb = 'T';
  transc = 'N';
  m = 4;
  n = 5;
  k = 3;
  alpha = 2.0;
  lda = k;
  a = r8mat_test ( transa, lda, m, k );
  ldb = n;
  b = r8mat_test ( transb, ldb, k, n );
  beta = 3.0;
  ldc = m;
  c = r8mat_test ( transc, ldc, m, n );

  dgemm ( transa, transb,  m, n, k, alpha, a, lda, b, ldb, beta, c, ldc );

  r8mat_print ( m, n, c, "  C = alpha * A' * B' + beta * C:" );

  free ( a );
  free ( b );
  free ( c );

  return;
}
/******************************************************************************/

void test02 ( )

/******************************************************************************/
/*
  Purpose:

    TEST02 tests DTRMM.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 April 2014

  Author:

    John Burkardt
*/
{
  double *a;
  double alpha;
  double *b;
  char diag;
  int i;
  int j;
  int lda;
  int ldb;
  int m;
  int n;
  char side;
  char transa;
  char transb;
  char uplo;

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  DTRMM multiplies a triangular matrix A and a\n" );
  printf ( "  rectangular matrix B\n" );
  printf ( "\n" );
  printf ( "  1: B = alpha * A  * B;\n" );
  printf ( "  2: B = alpha * A' * B;\n" );
/*
  B = alpha * A * B.
*/
  side = 'L';
  uplo = 'U';
  transa = 'N';
  diag = 'N';
  m = 4;
  n = 5;
  alpha = 2.0;
  lda = m;
  ldb = m;

  a = ( double * ) malloc ( lda * m * sizeof ( double ) );
  for ( j = 0; j < m; j++ )
  {
    for ( i = 0; i <= j; i++ )
    {
      a[i+j*lda] = ( double ) ( i + j + 2 );
    }
    for ( i = j + 1; i < m; i++ )
    {
      a[i+j*lda] = 0.0;
    }
  }

  transb = 'N';
  b = r8mat_test ( transb, ldb, m, n );

  dtrmm ( side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb );

  r8mat_print ( m, n, b, "  B = alpha * A * B:" );

  free ( a );
  free ( b );
/*
  B = alpha * A' * B.
*/
  side = 'L';
  uplo = 'U';
  transa = 'T';
  diag = 'N';
  m = 4;
  n = 5;
  alpha = 2.0;
  lda = m;
  ldb = m;

  a = ( double * ) malloc ( lda * m * sizeof ( double ) );
  for ( j = 0; j < m; j++ )
  {
    for ( i = 0; i <= j; i++ )
    {
      a[i+j*lda] = ( double ) ( i + j + 2 );
    }
    for ( i = j + 1; i < m; i++ )
    {
      a[i+j*lda] = 0.0;
    }
  }

  transb = 'N';
  b = r8mat_test ( transb, ldb, m, n );

  dtrmm ( side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb );

  r8mat_print ( m, n, b, "  B = alpha * A' * B:" );

  free ( a );
  free ( b );

  return;
}
/******************************************************************************/

void test03 ( )

/******************************************************************************/
/*
  Purpose:

    TEST03 tests DTRSM.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 April 2014

  Author:

    John Burkardt
*/
{
  double *a;
  double alpha;
  double *b;
  char diag;
  int i;
  int j;
  int lda;
  int ldb;
  int m;
  int n;
  char side;
  char transa;
  char transb;
  char uplo;

  printf ( "\n" );
  printf ( "TEST03\n" );
  printf ( "  DTRSM solves a linear system involving a triangular\n" );
  printf ( "  matrix A and a rectangular matrix B.\n" );
  printf ( "\n" );
  printf ( "  1: Solve A  * X  = alpha * B;\n" );
  printf ( "  2: Solve A' * X  = alpha * B;\n" );
  printf ( "  3: Solve X  * A  = alpha * B;\n" );
  printf ( "  4: Solve X  * A' = alpha * B;\n" );
/*
  Solve A * X = alpha * B.
*/
  side = 'L';
  uplo = 'U';
  transa = 'N';
  diag = 'N';
  m = 4;
  n = 5;
  alpha = 2.0;
  lda = m;
  ldb = m;

  a = ( double * ) malloc ( lda * m * sizeof ( double ) );

  for ( j = 0; j < m; j++ )
  {
    for ( i = 0; i <= j; i++ )
    {
      a[i+j*lda] = ( double ) ( i + j + 2 );
    }
    for ( i = j + 1; i < m; i++ )
    {
      a[i+j*lda] = 0.0;
    }
  }

  transb = 'N';
  b = r8mat_test ( transb, ldb, m, n );

  dtrsm ( side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb );

  r8mat_print ( m, n, b, "  X = inv ( A ) * alpha * B:" );

  free ( a );
  free ( b );
/*
  Solve A' * X = alpha * B.
*/
  side = 'L';
  uplo = 'U';
  transa = 'T';
  diag = 'N';
  m = 4;
  n = 5;
  alpha = 2.0;
  lda = m;
  ldb = m;

  a = ( double * ) malloc ( lda * m * sizeof ( double ) );
  for ( j = 0; j < m; j++ )
  {
    for ( i = 0; i <= j; i++ )
    {
      a[i+j*lda] = ( double ) ( i + j + 2 );
    }
    for ( i = j + 1; i < m; i++ )
    {
      a[i+j*lda] = 0.0;
    }
  }

  transb = 'N';
  b = r8mat_test ( transb, ldb, m, n );

  dtrsm ( side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb );

  r8mat_print ( m, n, b, "  X = inv ( A' ) * alpha * B:" );

  free ( a );
  free ( b );
/*
  Solve X * A = alpha * B.
*/
  side = 'R';
  uplo = 'U';
  transa = 'N';
  diag = 'N';
  m = 4;
  n = 5;
  alpha = 2.0;
  lda = n;
  ldb = m;

  a = ( double * ) malloc ( lda * n * sizeof ( double ) );
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i <= j; i++ )
    {
      a[i+j*lda] = ( double ) ( i + j + 2 );
    }
    for ( i = j + 1; i < n; i++ )
    {
      a[i+j*lda] = 0.0;
    }
  }

  transb = 'N';
  b = r8mat_test ( transb, ldb, m, n );

  dtrsm ( side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb );

  r8mat_print ( m, n, b, "  X = alpha * B * inv ( A ):" );

  free ( a );
  free ( b );
/*
  Solve X * A'' = alpha * B.
*/
  side = 'R';
  uplo = 'U';
  transa = 'T';
  diag = 'N';
  m = 4;
  n = 5;
  alpha = 2.0;
  lda = n;
  ldb = m;

  a = ( double * ) malloc ( lda * n * sizeof ( double ) );
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i <= j; i++ )
    {
      a[i+j*lda] = ( double ) ( i + j + 2 );
    }
    for ( i = j + 1; i < n; i++ )
    {
      a[i+j*lda] = 0.0;
    }
  }

  transb = 'N';
  b = r8mat_test ( transb, ldb, m, n );

  dtrsm ( side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb );

  r8mat_print ( m, n, b, "  X = alpha * B * inv ( A' ):" );

  free ( a );
  free ( b );

  return;
}

