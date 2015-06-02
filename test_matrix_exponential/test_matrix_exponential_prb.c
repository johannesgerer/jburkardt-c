# include <stdlib.h>
# include <stdio.h>
# include <complex.h>
# include <math.h>

# include "test_matrix_exponential.h"
# include "c8lib.h"
# include "r8lib.h"

int main ( void );
void test_matrix_exponential_test01 ( void );
void test_matrix_exponential_test02 ( void );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for TEST_MATRIX_EXPONENTIAL_PRB.

  Discussion:

    TEST_MATRIX_EXPONENTIAL_PRB tests the TEST_MATRIX_EXPONENTIAL library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    03 March 2013

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "TEST_MATRIX_EXPONENTIAL_TEST:\n" );
  printf ( "  C version\n" );
  printf ( "  Test the TEST_MATRIX_EXPONENTIAL library.\n" );
  printf ( "  The C8LIB and R8LIB libraries are needed.\n" );

  test_matrix_exponential_test01 ( );
  test_matrix_exponential_test02 ( );
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

    TEST_MATRIX_EXPONENTIAL_TEST01 retrieves real test data.

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

  test_num = r8mat_exp_test_num ( );

  for ( test = 1; test <= test_num; test++ )
  {
    printf ( "\n" );
    printf ( "  Test #%d\n", test );

    n = r8mat_exp_n ( test );

    r8mat_exp_story ( test );

    printf ( "\n" );
    printf ( "  Matrix order N = %d\n", n );

    a = r8mat_exp_a ( test, n );
    r8mat_print ( n, n, a, "  Matrix A:" );

    expa = r8mat_exp_expa ( test, n );
    r8mat_print ( n, n, expa, "  Exact Exponential exp(A):" );

    free ( a );
    free ( expa );
  }
  return;
}
/******************************************************************************/

void test_matrix_exponential_test02 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST_MATRIX_EXPONENTIAL_TEST02 retrieves complex test data.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    03 March 2013

  Author:

    John Burkardt
*/
{
  double complex *a;
  double complex *expa;
  int n;
  int test;
  int test_num;

  printf ( "\n" );
  printf ( "TEST_MATRIX_EXPONENTIAL_TEST02:\n" );
  printf ( "  Retrieve the data for each matrix exponential test.\n" );

  test_num = c8mat_exp_test_num ( );

  for ( test = 1; test <= test_num; test++ )
  {
    printf ( "\n" );
    printf ( "  Test #%d\n", test );

    n = c8mat_exp_n ( test );

    c8mat_exp_story ( test );

    printf ( "\n" );
    printf ( "  Matrix order N = %d\n", n );

    a = c8mat_exp_a ( test, n );
    c8mat_print ( n, n, a, "  Matrix A:" );

    expa = c8mat_exp_expa ( test, n );
    c8mat_print ( n, n, expa, "  Exact Exponential exp(A):" );

    free ( a );
    free ( expa );
  }
  return;
}
