# include <stdlib.h>
# include <stdio.h>

# include "test_ls.h"
# include "r8lib.h"

int main ( );
void test01 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for TEST_LS_PRB.

  Discussion:

    TEST_LS_PRB tests the TEST_LS library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 September 2012

  Author:

    John Burkardt
*/
{
  timestamp ( );

  printf ( "\n" );
  printf ( "TEST_LS_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the TEST_LS library.\n" );
  printf ( "  This test also requires the R8LIB library.\n" );

  test01 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "TEST_LS_PRB\n" );
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

    TEST01 summarizes the test data.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 September 2012

  Author:

    John Burkardt
*/
{
  double *a;
  double *b;
  double b_norm;
  int i;
  int j;
  int m;
  int n;
  int prob;
  int prob_num;
  double *r;
  double r_norm;
  double *x;
  double x_norm;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  Get each least squares test and compute the maximum residual.\n" );
  printf ( "  The L2 norm of the residual MUST be no greater than\n" );
  printf ( "  the L2 norm of the right hand side, else 0 is a better solution.\n" );

  prob_num = p00_prob_num ( );

  printf ( "\n" );
  printf ( "  Number of problems = %d\n", prob_num );
  printf ( "\n" );
  printf ( "  Index     M     N     ||B||         ||X||         ||R||\n" );
  printf ( "\n" );

  for ( prob = 1; prob <= prob_num; prob++ )
  {
    m = p00_m ( prob );
    n = p00_n ( prob );

    a = p00_a ( prob, m, n );
    b = p00_b ( prob, m );
    x = p00_x ( prob, n );

    r = ( double * ) malloc ( m * sizeof ( double ) );

    for ( i = 0; i < m; i++ )
    {
      r[i] = - b[i];
      for ( j = 0; j < n; j++ )
      {
        r[i] = r[i] + a[i+j*m] * x[j];
      }
    }

    b_norm = r8vec_norm ( m, b );
    x_norm = r8vec_norm ( n, x );
    r_norm = r8vec_norm ( m, r );

    printf ( "  %5d  %4d  %4d  %12g  %12g  %12g\n", prob, m, n, b_norm, x_norm, r_norm );

    free ( a );
    free ( b );
    free ( r );
    free ( x );
  }
  return;
}
