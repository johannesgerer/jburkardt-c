# include <stdlib.h>
# include <stdio.h>
# include <time.h>

# include "eispack.h"

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );
void test04 ( );
void test05 ( );
void test06 ( );
void test065 ( );
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

/******************************************************************************/

int main ()

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for EISPACK_PRB1.

  Discussion:

    EISPACK_PRB1 calls the EISPACK sample programs.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 November 2012

  Author:

    John Burkardt
*/
{
  timestamp ( );

  printf ( "\n" );
  printf ( "EISPACK_PRB\n" );
  printf ( "  C version.\n" );
  printf ( "  Test the EISPACK library.\n" );
/*
  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
  test05 ( );
*/
  test06 ( );
  test065 ( );
  test07 ( );
/*
  test08 ( );
  test09 ( );
  test10 ( );

  test11 ( );
  test12 ( );
  test13 ( );
  test14 ( );
  test15 ( );
  test16 ( );
*/
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "EISPACK_PRB1\n" );
  printf ( "  Normal end of execution.\n" );

  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void test06 ( )

/******************************************************************************/
/*
  Purpose:

    TEST06 tests RS.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 November 2012

  Author:

    John Burkardt
*/
{
  double a[4*4] = {
    5.0, 4.0, 1.0, 1.0,
    4.0, 5.0, 1.0, 1.0,
    1.0, 1.0, 4.0, 2.0,
    1.0, 1.0, 2.0, 4.0 };
  double a2[4*4];
  int i;
  int ierr;
  int j;
  int k;
  int matz;
  int n = 4;
  double *r;
  double w[4];
  double x[4*4];
/*
  Save a copy of the matrix.
*/
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      a2[i+j*n] = a[i+j*n];
    }
  }
  printf ( "\n" );
  printf ( "TEST06\n" );
  printf ( "  RS computes the eigenvalues and eigenvectors\n" );
  printf ( "  of a real symmetric matrix.\n" );
  printf ( "\n" );
  printf ( "  Matrix order = %d\n", n );

  r8mat_print ( n, n, a, "  The matrix A:" );

  matz = 1;

  ierr = rs ( n, a, w, matz, x );

  if ( ierr != 0 )
  {
    printf ( "\n" );
    printf ( "TEST06 - Warning!\n" );
    printf ( "  The error return flag IERR = %d\n", ierr );
    return;
  }

  r8vec_print ( n, w, "  The eigenvalues Lambda:" );

  if ( matz != 0 )
  {
    r8mat_print ( n, n, x, "  The eigenvector matrix:" );

    r = r8mat_mm_new ( n, n, n, a2, x );

    for ( j = 0; j < n; j++ )
    {
      for ( i = 0; i < n; i++ )
      {
        r[i+j*n] = r[i+j*n] - w[j] * x[i+j*n];
      }
    }

    r8mat_print ( n, n, r, "  The residual (A-Lambda*I)*X:" );
  }

  free ( r );

  return;
}
/******************************************************************************/

void test065 ( )

/******************************************************************************/
/*
  Purpose:

    TEST065 tests RS.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 November 2012

  Author:

    John Burkardt
*/
{
  double *a;
  double a2[3*3];
  int i;
  int ierr;
  int j;
  int k;
  int matz;
  int n = 3;
  double *r;
  int seed;
  double t;
  double w[3];
  double x[3*3];

  printf ( "\n" );
  printf ( "TEST065\n" );
  printf ( "  RS computes the eigenvalues and eigenvectors\n" );
  printf ( "  of a real symmetric matrix.\n" );
  printf ( "\n" );
  printf ( "  Matrix order = %d\n", n );

  seed = 123456789;

  a = r8mat_uniform_01_new ( n, n, &seed );
  for ( i = 0; i < n - 1; i++ )
  {
    for ( j = i + 1; j < n; j++ )
    {
      t = 0.5 * ( a[i+j*n] + a[j+i*n] );
      a[i+j*n] = t;
      a[j+i*n] = t;
    }
  }
/*
  Save a copy of the matrix.
*/
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      a2[i+j*n] = a[i+j*n];
    }
  }

  r8mat_print ( n, n, a, "  The matrix A:" );

  matz = 1;

  ierr = rs ( n, a, w, matz, x );

  if ( ierr != 0 )
  {
    printf ( "\n" );
    printf ( "TEST065 - Warning!\n" );
    printf ( "  The error return flag IERR = %d\n", ierr );
    return;
  }

  r8vec_print ( n, w, "  The eigenvalues Lambda:" );

  if ( matz != 0 )
  {
    r8mat_print ( n, n, x, "  The eigenvector matrix:" );

    r = r8mat_mm_new ( n, n, n, a2, x );

    for ( j = 0; j < n; j++ )
    {
      for ( i = 0; i < n; i++ )
      {
        r[i+j*n] = r[i+j*n] - w[j] * x[i+j*n];
      }
    }

    r8mat_print ( n, n, r, "  The residual (A-Lambda*I)*X:" );

    free ( r );
  }

  free ( a );

  return;
}
/******************************************************************************/

void test07 ( )

/******************************************************************************/
/*
  Purpose:

    TEST07 tests RSB.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 November 2012

  Author:

    John Burkardt
*/
{
  double *a;
  double *a2;
  int i;
  int ierr;
  int j;
  int matz;
  int mb = 2;
  int n = 5;
  double *r;
  double *w;
  double *x;

  a = ( double * ) malloc ( n * mb * sizeof ( double ) );

  for ( j = 0; j < mb; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      a[i+j*n] = 0.0;
    }
  }
  j = mb - 1;
  for ( i = 0; i < n; i++ )
  {
    a[i+j*n] = 2.0;
  }
  j = 0;
  for ( i = 1; i < n; i++ )
  {
    a[i+j*n] = -1.0;
  }

  a2 = ( double * ) malloc ( n * n * sizeof ( double ) );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      if ( i == j )
      {
        a2[i+j*n] = 2.0;
      }
      else if ( abs ( i - j ) == 1 )
      {
        a2[i+j*n] = - 1.0;
      }
      else
      {
        a2[i+j*n] = 0.0;
      }
    }
  }

  printf ( "\n" );
  printf ( "TEST07\n" );
  printf ( "  RSB computes the eigenvalues and eigenvectors\n" );
  printf ( "  of a real symmetric band matrix.\n" );
  printf ( "\n" );
  printf ( "  Matrix order = %d\n", n );

  r8mat_print ( n, n, a2, "  The matrix A:" );

  w = ( double * ) malloc ( n * sizeof ( double ) );
  x = ( double * ) malloc ( n * n * sizeof ( double ) );
  matz = 1;

  ierr = rsb ( n, mb, a, w, matz, x );

  if ( ierr != 0 )
  {
    printf ( "\n" );
    printf ( "TEST07 - Warning!\n" );
    printf ( "  The error return flag IERR = %d\n", ierr );
    return;
  }

  r8vec_print ( n, w, "  The eigenvalues Lambda:" );

  if ( matz != 0 )
  {
    r8mat_print ( n, n, x, "  The eigenvector matrix X:" );

    r = r8mat_mm_new ( n, n, n, a2, x );

    for ( j = 0; j < n; j++ )
    {
      for ( i = 0; i < n; i++ )
      {
        r[i+j*n] = r[i+j*n] - x[i+j*n] * w[j];
      }
    }
    r8mat_print ( n, n, r, "  The residual (A-Lambda*I)*X:" );

    free ( r );
  }

  free ( a );
  free ( a2 );
  free ( w );
  free ( x );

  return;
}
