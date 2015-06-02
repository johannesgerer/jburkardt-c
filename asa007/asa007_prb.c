# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "asa007.h"

int main ( );
void test01 ( );
void test02 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for ASA007_PRB.

  Discussion:

    ASA007_PRB tests the ASA007 library.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    25 October 2010

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "ASA007_PRB:\n" );
  printf ( "  C version\n" );
  printf ( "  Test the ASA007 library.\n" );

  test01 ( );
  test02 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "ASA007_PRB:\n" );
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

    TEST01 demonstrates the use of SYMINV.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    25 October 2010

  Author:

    John Burkardt
*/
{
# define N_MAX 15

  double a[(N_MAX*(N_MAX+1))/2];
  double afull[N_MAX*N_MAX];
  double c[(N_MAX*(N_MAX+1))/2];
  double cfull[N_MAX*N_MAX];
  double cta;
  double diff;
  int i;
  int ifault;
  int j;
  int k;
  int l;
  int n;
  int nullty;
  double w[N_MAX];

  printf ( "\n" );
  printf ( "TEST01:\n" );
  printf ( "  SYMINV computes the inverse of a positive\n" );
  printf ( "  definite symmetric matrix.\n" );
  printf ( "  A compressed storage format is used\n" );
  printf ( "\n" );
  printf ( "  Here we look at the matrix A which is\n" );
  printf ( "  N+1 on the diagonal and\n" );
  printf ( "  N   on the off diagonals.\n" );

  for ( n = 1; n <= N_MAX; n++ )
  {
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

    syminv ( a, n, c, w, &nullty, &ifault );

    printf ( "\n" );
    printf ( "  Matrix order N = %d\n", n );
    printf ( "  Maxtrix nullity NULLTY = %d\n", nullty );

    k = 0;
    for ( j = 1; j <= n; j++ )
    {
      for ( i = 1; i < j; i++ )
      {
        cfull[i-1+(j-1)*n] = c[k];
        cfull[j-1+(i-1)*n] = c[k];
        k = k + 1;
      }
      cfull[j-1+(j-1)*n] = c[k];
      k = k + 1;
    }

    k = 0;
    for ( j = 1; j <= n; j++ )
    {
      for ( i = 1; i < j; i++ )
      {
        afull[i-1+(j-1)*n] = a[k];
        afull[j-1+(i-1)*n] = a[k];
        k = k + 1;
      }
      afull[j-1+(j-1)*n] = a[k];
      k = k + 1;
    }
/*
  Compute C * A - I.
*/
    diff = 0.0;
    for ( i = 1; i <= n; i++ )
    {
      for ( j = 1; j <= i; j++ )
      {
        cta = 0.0;
        for ( k = 1; k <= n; k++ )
        {
          cta = cta + cfull[i-1+(k-1)*n] * afull[k-1+(j-1)*n];
        }
        if ( i == j )
        {
          diff = diff + pow ( 1.0 - cta, 2 );
        }
        else
        {
          diff = diff + cta * cta;
        }
      }
    }

    diff = sqrt ( diff );

    printf ( "  RMS ( C * A - I ) = %e\n", diff );
  }

  return;
# undef N_MAX
}
/******************************************************************************/

void test02 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST02 demonstrates the use of SYMINV.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    25 October 2010

  Author:

    John Burkardt
*/
{
# define N_MAX 15

  double a[(N_MAX*(N_MAX+1))/2];
  double afull[N_MAX*N_MAX];
  double c[(N_MAX*(N_MAX+1))/2];
  double cfull[N_MAX*N_MAX];
  double cta;
  double diff;
  int i;
  int ifault;
  int j;
  int k;
  int l;
  int n;
  int nullty;
  double w[N_MAX];

  printf ( "\n" );
  printf ( "TEST02:\n" );
  printf ( "  SYMINV computes the inverse of a positive\n" );
  printf ( "  definite symmetric matrix.\n" );
  printf ( "  A compressed storage format is used\n" );
  printf ( "\n" );
  printf ( "  Here we look at the Hilbert matrix\n" );
  printf ( "  A(I,J) = 1/(I+J-1)\n" );
  printf ( "\n" );
  printf ( "  For this matrix, we expect errors to grow quickly.\n" );

  for ( n = 1; n <= N_MAX; n++ )
  {
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

    syminv ( a, n, c, w, &nullty, &ifault );

    printf ( "\n" );
    printf ( "  Matrix order N = %d\n", n );
    printf ( "  Maxtrix nullity NULLTY = %d\n", nullty );

    k = 0;
    for ( j = 1; j <= n; j++ )
    {
      for ( i = 1; i < j; i++ )
      {
        cfull[i-1+(j-1)*n] = c[k];
        cfull[j-1+(i-1)*n] = c[k];
        k = k + 1;
      }
      cfull[j-1+(j-1)*n] = c[k];
      k = k + 1;
    }

    k = 0;
    for ( j = 1; j <= n; j++ )
    {
      for ( i = 1; i < j; i++ )
      {
        afull[i-1+(j-1)*n] = a[k];
        afull[j-1+(i-1)*n] = a[k];
        k = k + 1;
      }
      afull[j-1+(j-1)*n] = a[k];
      k = k + 1;
    }
/*
  Compute C * A - I.
*/
    diff = 0.0;
    for ( i = 1; i <= n; i++ )
    {
      for ( j = 1; j <= i; j++ )
      {
        cta = 0.0;
        for ( k = 1; k <= n; k++ )
        {
          cta = cta + cfull[i-1+(k-1)*n] * afull[k-1+(j-1)*n];
        }
        if ( i == j )
        {
          diff = diff + pow ( 1.0 - cta, 2 );
        }
        else
        {
          diff = diff + cta * cta;
        }
      }
    }

    diff = sqrt ( diff );

    printf ( "  RMS ( C * A - I ) = %e\n", diff );
  }

  return;
# undef N_MAX
}
