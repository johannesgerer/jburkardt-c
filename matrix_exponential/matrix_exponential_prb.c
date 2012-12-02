# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "matrix_exponential.h"
# include "test_matrix_exponential.h"
# include "r8lib.h"

int main ( void );
void matrix_exponential_test01 ( void );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    MATRIX_EXPONENTIAL_TEST tests some matrix exponential algorithms.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 December 2011

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "MATRIX_EXPONENTIAL_TEST:\n" );
  printf ( "  C version\n" );
  printf ( "  Test the MATRIX_EXPONENTIAL library.\n" );
  printf ( "  The R8LIB library is needed.\n" );
  printf ( "  This test needs the TEST_MATRIX_EXPONENTIAL library.\n" );

  matrix_exponential_test01 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "MATRIX_EXPONENTIAL_TEST:\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void matrix_exponential_test01 ( void )

/******************************************************************************/
/*
  Purpose:

    MATRIX_EXPONENTIAL_TEST01 compares matrix exponential algorithms.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 December 2011

  Author:

    John Burkardt
*/
{
  double *a;
  double *a_exp;
  int n;
  int test;
  int test_num;

  printf ( "\n" );
  printf ( "MATRIX_EXPONENTIAL_TEST01:\n" );
  printf ( "  EXPM is MATLAB's matrix exponential function\n" );
  printf ( "  EXPM11 is an equivalent to EXPM\n" );
  printf ( "  EXPM2 uses a Taylor series approach\n" );
  printf ( "  EXPM3 relies on an eigenvalue calculation.\n" );

  test_num = mexp_test_num ( );

  for ( test = 1; test <= test_num; test++ )
  {
    printf ( "\n" );
    printf ( "  Test #%d\n", test );

    mexp_story ( test );

    n = mexp_n ( test );

    printf ( "  Matrix order N = %d\n", n );

    a = mexp_a ( test, n );

    r8mat_print ( n, n, a, "  Matrix:" );

    a_exp = expm11 ( n, a );
    r8mat_print ( n, n, a_exp, "  EXPM1(A):" );
    free ( a_exp );

    a_exp = expm2 ( n, a );
    r8mat_print ( n, n, a_exp, "  EXPM2(A):" );
    free ( a_exp );
/*
    a_exp = expm3 ( n, a );
    r8mat_print ( n, n, a_exp, "  EXPM3(A):" );
    free ( a_exp );
*/
    a_exp = mexp_expa ( test, n );
    r8mat_print ( n, n, a_exp, "  Exact Exponential:" );
    free ( a_exp );

    free ( a );
  }

  return;
}
