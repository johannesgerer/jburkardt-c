# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "test_matrix_exponential.h"
# include "r8lib.h"

int main ( void );
void test_matrix_exponential_test01 ( void );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    TEST_MATRIX_EXPONENTIAL_TEST tests some matrix exponential algorithms.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    29 November 2011

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "TEST_MATRIX_EXPONENTIAL_TEST:\n" );
  printf ( "  C version\n" );
  printf ( "  Test the TEST_MATRIX_EXPONENTIAL library.\n" );
  printf ( "  The R8LIB library is needed.\n" );

  test_matrix_exponential_test01 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "TEST_MATRIX_EXPONENTIAL_TEST:\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void test_matrix_exponential_test01 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST_MATRIX_EXPONENTIAL_TEST01 retrieves the test data.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    29 November 2011

  Author:

    John Burkardt
*/
{
  double *a;
  double *expa;
  int n;
  int test;
  int test_num;

  printf ( "\n" );
  printf ( "TEST_MATRIX_EXPONENTIAL_TEST01:\n" );
  printf ( "  Retrieve the data for each matrix exponential test.\n" );

  test_num = mexp_test_num ( );

  for ( test = 1; test <= test_num; test++ )
  {
    printf ( "\n" );
    printf ( "  Test #%d\n", test );

    n = mexp_n ( test );

    mexp_story ( test );

    printf ( "\n" );
    printf ( "  Matrix order N = %d\n", n );

    a = mexp_a ( test, n );
    r8mat_print ( n, n, a, "  Matrix A:" );

    expa = mexp_expa ( test, n );
    r8mat_print ( n, n, expa, "  Exact Exponential exp(A):" );

    free ( a );
    free ( expa );
  }
  return;
}
