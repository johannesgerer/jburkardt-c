# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <string.h>

# include "test_eigen.h"

int main ( );
void test01 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for TEST_EIGEN_PRB.

  Discussion:

    TEST_EIGEN_PRB tests the TEST_EIGEN library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    23 February 2012

  Author:

    John Burkardt
*/
{
  timestamp ( );

  printf ( "\n" );
  printf ( "TEST_EIGEN_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the TEST_EIGEN library.\n" );

  test01 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "TEST_EIGEN_PRB\n" );
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

    TEST01 tests the use of R8SYMM_TEST to make symmetric test matrices.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    23 February 2012

  Author:

    John Burkardt
*/
{
  double *a;
  double *aq;
  int *bin;
  double *bin_limit;
  int bin_num = 10;
  int i;
  int j;
  double *lambda;
  double *lambda2;
  double lambda_dev = 1.0;
  double lambda_max;
  double lambda_mean = 1.0;
  double lambda_min;
  int n = 100;
  double *q;
  int seed = 123456789;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  R8SYMM_TEST makes an arbitrary size symmetric matrix\n" );
  printf ( "  with known eigenvalues and eigenvectors.\n" );

  a = ( double * ) malloc ( n * n * sizeof ( double ) );
  q = ( double * ) malloc ( n * n * sizeof ( double ) );
  lambda = ( double * ) malloc ( n * sizeof ( double ) );

  r8symm_test ( n, lambda_mean, lambda_dev, &seed, a, q, lambda );
/*
  Get the eigenvalue range.
*/
  lambda_min = r8vec_min ( n, lambda );
  lambda_max = r8vec_max ( n, lambda );

  printf ( "\n" );
  printf ( "  LAMBDA_MIN = %g\n", lambda_min );
  printf ( "  LAMBDA_MAX = %g\n", lambda_max );
/*
  Bin the eigenvalues.
*/
  bin = ( int * ) malloc ( ( bin_num + 2 ) * sizeof ( int ) );
  bin_limit = ( double * ) malloc ( ( bin_num + 1 ) * sizeof ( double ) );

  r8vec_bin ( n, lambda, bin_num, lambda_min, lambda_max, bin, bin_limit );

  r8bin_print ( bin_num, bin, bin_limit, "  Lambda bins:" );

  if ( 0 )
  {
    r8mat_print ( n, n, a, "  The matrix A:" );
  }

  if ( 0 )
  {
    r8mat_print ( n, n, q, "  The eigenvector matrix Q:" );
  }

  aq = r8mat_mm_new ( n, n, n, a, q );

  lambda2 = ( double * ) malloc ( n * sizeof ( double ) );

  for ( j = 0; j < n; j++ )
  {
    lambda2[j] = 0.0;
    for ( i = 0; i < n; i++ )
    {
      lambda2[j] = lambda2[j] + pow ( aq[i+j*n], 2 );
    }
    lambda2[j] = sqrt ( lambda2[j] );
  }

  r8vec2_print ( n, lambda, lambda2, "  LAMBDA versus the column norms of A*Q:" );

  free ( a );
  free ( aq );
  free ( bin );
  free ( bin_limit );
  free ( lambda );
  free ( lambda2 );
  free ( q );

  return;
}
