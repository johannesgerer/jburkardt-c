# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <complex.h>

# include "matrix_exponential.h"
# include "test_matrix_exponential.h"
# include "c8lib.h"
# include "r8lib.h"

int main ( void );
void matrix_exponential_test01 ( void );
void matrix_exponential_test02 ( void );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for MATRIX_EXPONENTIAL_PRB.

  Discussion:

    MATRIX_EXPONENTIAL_PRB tests the MATRIX_EXPONENTIAL library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 March 2013

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "MATRIX_EXPONENTIAL_PRB:\n" );
  printf ( "  C version\n" );
  printf ( "  Test the MATRIX_EXPONENTIAL library.\n" );
  printf ( "  The C8LIB and R8LIB libraries are needed.\n" );
  printf ( "  This test needs the TEST_MATRIX_EXPONENTIAL library.\n" );

  matrix_exponential_test01 ( );
  matrix_exponential_test02 ( );

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

    MATRIX_EXPONENTIAL_TEST01 compares real matrix exponential algorithms.

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
  printf ( "  R8MAT_EXPM1 is an equivalent to EXPM\n" );
  printf ( "  R8MAT_EXPM2 uses a Taylor series approach\n" );
  printf ( "  R8MAT_EXPM3 relies on an eigenvalue calculation.\n" );

  test_num = r8mat_exp_test_num ( );

  for ( test = 1; test <= test_num; test++ )
  {
    printf ( "\n" );
    printf ( "  Test #%d\n", test );

    r8mat_exp_story ( test );

    n = r8mat_exp_n ( test );

    printf ( "  Matrix order N = %d\n", n );

    a = r8mat_exp_a ( test, n );

    r8mat_print ( n, n, a, "  Matrix:" );

    a_exp = r8mat_expm1 ( n, a );
    r8mat_print ( n, n, a_exp, "  R8MAT_EXPM1(A):" );
    free ( a_exp );

    a_exp = r8mat_expm2 ( n, a );
    r8mat_print ( n, n, a_exp, "  R8MAT_EXPM2(A):" );
    free ( a_exp );
/*
    a_exp = r8mat_expm3 ( n, a );
    r8mat_print ( n, n, a_exp, "  R8MAT_EXPM3(A):" );
    free ( a_exp );
*/
    a_exp = r8mat_exp_expa ( test, n );
    r8mat_print ( n, n, a_exp, "  Exact Exponential:" );
    free ( a_exp );

    free ( a );
  }

  return;
}
/******************************************************************************/

void matrix_exponential_test02 ( void )

/******************************************************************************/
/*
  Purpose:

    MATRIX_EXPONENTIAL_TEST02 compares complex matrix exponential algorithms.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 March 2013

  Author:

    John Burkardt
*/
{
  double complex *a;
  double complex *a_exp;
  int n;
  int test;
  int test_num;

  printf ( "\n" );
  printf ( "MATRIX_EXPONENTIAL_TEST02:\n" );
  printf ( "  EXPM is MATLAB's matrix exponential function\n" );
  printf ( "  C8MAT_EXPM1 is an equivalent to EXPM\n" );
  printf ( "  C8MAT_EXPM2 uses a Taylor series approach\n" );
  printf ( "  C8MAT_EXPM3 relies on an eigenvalue calculation.\n" );

  test_num = c8mat_exp_test_num ( );

  for ( test = 1; test <= test_num; test++ )
  {
    printf ( "\n" );
    printf ( "  Test #%d\n", test );

    c8mat_exp_story ( test );

    n = c8mat_exp_n ( test );

    printf ( "  Matrix order N = %d\n", n );

    a = c8mat_exp_a ( test, n );

    c8mat_print ( n, n, a, "  Matrix:" );
    a_exp = c8mat_expm1 ( n, a );
    c8mat_print ( n, n, a_exp, "  C8MAT_EXPM1(A):" );
    free ( a_exp );
/*
    a_exp = c8mat_expm2 ( n, a );
    c8mat_print ( n, n, a_exp, "  C8MAT_EXPM2(A):" );
    free ( a_exp );
*/
/*
    a_exp = c8mat_expm3 ( n, a );
    c8mat_print ( n, n, a_exp, "  C8MAT_EXPM3(A):" );
    free ( a_exp );
*/
    a_exp = c8mat_exp_expa ( test, n );
    c8mat_print ( n, n, a_exp, "  Exact Exponential:" );
    free ( a_exp );

    free ( a );
  }

  return;
}
