# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "asa006.h"

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for ASA006_PRB.

  Discussion:

    ASA006_PRB tests the ASA006 library.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    24 October 2010

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "ASA006_PRB:\n" );
  printf ( "  C version\n" );
  printf ( "  Test the ASA006 library.\n" );

  test01 ( );
  test02 ( );
  test03 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "ASA006_PRB:\n" );
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

    TEST01 demonstrates the use of CHOLESKY.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    01 February 2008

  Author:

    John Burkardt
*/
{
# define N_MAX 15

  double a[(N_MAX*(N_MAX+1))/2];
  double diff;
  int i;
  int ifault;
  int j;
  int k;
  int l;
  int n;
  int nn;
  int nullty;
  double r[N_MAX];
  double rmax;
  double u[(N_MAX*(N_MAX+1))/2];
  double ufull[N_MAX*N_MAX];
  double utu;

  printf ( "\n" );
  printf ( "TEST01:\n" );
  printf ( "  CHOLESKY computes the Cholesky factorization\n" );
  printf ( "  of a positive definite symmetric matrix.\n" );
  printf ( "  A compressed storage format is used\n" );
  printf ( "\n" );
  printf ( "  Here we look at the matrix A which is\n" );
  printf ( "  N+1 on the diagonal and\n" );
  printf ( "  N   on the off diagonals.\n" );

  for ( n = 1; n <= N_MAX; n++ )
  {
    nn = ( n * ( n + 1 ) ) / 2;
/*
  Set A to the lower triangle of the matrix which is N+1 on the diagonal
  and N on the off diagonals.
*/
    k = 0;
    for ( i = 1; i <= n; i++ )
    {
      for ( j = 1; j < i; j++ )
      {
        a[k] = ( double ) ( n );
        k = k + 1;
      }
      a[k] = ( double ) ( n + 1 );
      k = k + 1;
    }

    cholesky ( a, n, nn, u, &nullty, &ifault );

    printf ( "\n" );
    printf ( "  Matrix order N = %d\n", n );
    printf ( "  Maxtrix nullity NULLTY = %d\n", nullty );

    k = 0;
    for ( j = 1; j <= n; j++ )
    {
      for ( i = 1; i <= j; i++ )
      {
        ufull[i-1+(j-1)*n] = u[k];
        k = k + 1;
      }
      for ( i = j + 1; i <= n; i++ )
      {
        ufull[i-1+(j-1)*n] = 0.0;
      }
    }
/*
  Compute U' * U and compare to A.
*/
    k = 0;
    diff = 0.0;
    for ( i = 1; i <= n; i++ )
    {
      for ( j = 1; j <= i; j++ )
      {
        utu = 0.0;
        for ( l = 1; l <= n; l++ )
        {
          utu = utu + ufull[l-1+(i-1)*n] * ufull[l-1+(j-1)*n];
        }
        diff = diff + ( a[k] - utu ) * ( a[k] - utu );
        k = k + 1;
      }
    }

    diff = sqrt ( diff );

    printf ( "  RMS ( A - U'*U ) = %e\n", diff );
  }

  return;
# undef N_MAX
}
/******************************************************************************/

void test02 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST02 demonstrates the use of CHOLESKY.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    01 February 2008

  Author:

    John Burkardt
*/
{
# define N_MAX 15

  double a[(N_MAX*(N_MAX+1))/2];
  double diff;
  int i;
  int ifault;
  int j;
  int k;
  int l;
  int n;
  int nn;
  int nullty;
  double u[(N_MAX*(N_MAX+1))/2];
  double ufull[N_MAX*N_MAX];
  double utu;

  printf ( "\n" );
  printf ( "TEST02:\n" );
  printf ( "  CHOLESKY computes the Cholesky factorization\n" );
  printf ( "  of a positive definite symmetric matrix.\n" );
  printf ( "  A compressed storage format is used\n" );
  printf ( "\n" );
  printf ( "  Here we look at the Hilbert matrix\n" );
  printf ( "  A(I,J) = 1/(I+J-1)\n" );
  printf ( "\n" );
  printf ( "  For this matrix, we expect errors to grow quickly.\n" );

  for ( n = 1; n <= N_MAX; n++ )
  {
    nn = ( n * ( n + 1 ) ) / 2;
/*
  Set A to the Hilbert matrix.
*/
    k = 0;
    for ( i = 1; i <= n; i++ )
    {
      for ( j = 1; j <= i; j++ )
      {
        a[k] = 1.0 / ( double ) ( i + j - 1 );
        k = k + 1;
      }
    }

    cholesky ( a, n, nn, u, &nullty, &ifault );

    printf ( "\n" );
    printf ( "  Matrix order N = %d\n", n );
    printf ( "  Maxtrix nullity NULLTY = %d\n", nullty );

    k = 0;
    for ( j = 1; j <= n; j++ )
    {
      for ( i = 1; i <= j; i++ )
      {
        ufull[i-1+(j-1)*n] = u[k];
        k = k + 1;
      }
      for ( i = j + 1; i <= n; i++ )
      {
        ufull[i-1+(j-1)*n] = 0.0;
      }
    }
/*
  Compute U' * U and compare to A.
*/
    k = 0;
    diff = 0.0;
    for ( i = 1; i <= n; i++ )
    {
      for ( j = 1; j <= i; j++ )
      {
        utu = 0.0;
        for ( l = 1; l <= n; l++ )
        {
          utu = utu + ufull[l-1+(i-1)*n] * ufull[l-1+(j-1)*n];
        }
        diff = diff + ( a[k] - utu ) * ( a[k] - utu );
        k = k + 1;
      }
    }

    diff = sqrt ( diff );

    printf ( "  RMS ( A - U'*U ) = %e\n", diff );
  }

  return;
# undef N_MAX
}
/******************************************************************************/

void test03 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST03 demonstrates the use of SUBCHL.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    10 February 2008

  Author:

    John Burkardt
*/
{
# define N_MAX 15
# define NN_MAX ((N_MAX*(N_MAX+1))/2)

  double a[NN_MAX];
  int b[N_MAX];
  double det;
  double diff;
  int i;
  int ifault;
  int j;
  int k;
  int l;
  int n;
  int nullty;
  double u[NN_MAX];
  double ufull[N_MAX*N_MAX];
  double utu;

  printf ( "\n" );
  printf ( "TEST03:\n" );
  printf ( "  SUBCHL computes the Cholesky factor\n" );
  printf ( "  of a submatrix\n" );
  printf ( "  of a positive definite symmetric matrix.\n" );
  printf ( "  A compressed storage format is used.\n" );
  printf ( "\n" );
  printf ( "  Here we look at the Hilbert matrix\n" );
  printf ( "  A(I,J) = 1/(I+J-1).\n" );
  printf ( "\n" );
  printf ( "  For this particular matrix, we expect the\n" );
  printf ( "  errors to grow rapidly.\n" );
/*
  Set A to the N_MAX order Hilbert matrix.
*/
  k = 0;
  for ( i = 1; i <= N_MAX; i++ )
  {
    for ( j = 1; j <= i; j++ )
    {

      a[k] = 1.0 / ( double ) ( i + j - 1 );
      k = k + 1;
    }
  }

  for ( n = 1; n <= N_MAX; n++ )
  {
    for ( i = 1; i <= n; i++ )
    {
      b[i-1] = i;
    }

    subchl ( a, b, n, u, &nullty, &ifault, NN_MAX, &det );

    printf ( "\n" );
    printf ( "  Matrix order N = %d\n", n );
    printf ( "  Maxtrix nullity NULLTY = %d\n", nullty );
    printf ( "  Matrix determinant DET = %f\n", det );

    k = 0;
    for ( j = 1; j <= n; j++ )
    {
      for ( i = 1; i <= j; i++ )
      {
        k = k + 1;
        ufull[i-1+(j-1)*n] = u[k-1];
      }
      for ( i = j + 1; i <= n; i++ )
      {
        ufull[i-1+(j-1)*n] = 0.0;
      }
    }
/*
  Compute U' * U and compare to A.
*/
    k = 0;
    diff = 0.0;
    for ( i = 1; i <= n; i++ )
    {
      for ( j = 1; j <= i; j++ )
      {
        k = k + 1;
        utu = 0.0;
        for ( l = 1; l <= n; l++ )
        {
          utu = utu + ufull[l-1+(i-1)*n] * ufull[l-1+(j-1)*n];
        }
        diff = diff + pow ( a[k-1] - utu, 2 );
      }
    }
    diff = sqrt ( diff );
    printf ( "  RMS ( A - U'*U ) = %e\n", diff );
  }

  return;
# undef N_MAX
# undef NN_MAX
}
