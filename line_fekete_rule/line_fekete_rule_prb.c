# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>

# include "line_fekete_rule.h"
# include "qr_solve.h"
# include "r8lib.h"

int main ( );
void test01 ( int m );
void test02 ( int m );
void test03 ( int m );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for LINE_FEKETE_RULE_PRB.

  Discussion:

    LINE_FEKETE_RULE_PRB tests the LINE_FEKETE_RULE library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 April 2014

  Author:

    John Burkardt
*/
{
# define TEST_NUM 3

  int m;
  int m_test[TEST_NUM] = { 5, 11, 21 };
  int test;
  int test_num = TEST_NUM;

  timestamp ( );
  printf ( "\n" );
  printf ( "LINE_FEKETE_RULE_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the LINE_FEKETE_RULE library.\n" );

  for ( test = 0; test < test_num; test++ )
  {
    m = m_test[test];
    test01 ( m );
  }

  for ( test = 0; test < test_num; test++ )
  {
    m = m_test[test];
    test02 ( m );
  }

  for ( test = 0; test < test_num; test++ )
  {
    m = m_test[test];
    test03 ( m );
  }
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "LINE_FEKETE_RULE_PRB\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
# undef TEST_NUM
}
/******************************************************************************/

void test01 ( int m )

/******************************************************************************/
/*
  Purpose:

    TEST01 seeks Fekete points in [-1,+1].

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 April 2014

  Author:

    John Burkardt

  Reference:

    Alvise Sommariva, Marco Vianello,
    Computing approximate Fekete points by QR factorizations of Vandermonde 
    matrices,
    Computers and Mathematics with Applications,
    Volume 57, 2009, pages 1324-1336.

  Parameters:

    Input, int M, the dimension of the polynomial space.
*/
{
# define N 5001

  double a;
  double b;
  int n = N;
  int nf;
  double *wf;
  double wf_sum;
  double *x;
  double *xf;

  a = -1.0;
  b = +1.0;
  x = r8vec_linspace_new ( n, a, b );

  printf ( "\n" );
  printf ( "TEST01:\n" );
  printf ( "  Seek Fekete points in [%f,%f]\n", a, b );
  printf ( "  using %d equally spaced sample points\n", n );
  printf ( "  for polynomials of degree M = %d\n", m );
  printf ( "  using the monomial basis and uniform weight.\n" );

  wf = ( double * ) malloc ( m * sizeof ( double ) );
  xf = ( double * ) malloc ( m * sizeof ( double ) );
  line_fekete_monomial ( m, a, b, n, x, &nf, xf, wf );

  printf ( "\n" );
  printf ( "  NF = %d\n", nf );
  r8vec_print ( nf, xf, "  Estimated Fekete points XF:" );

  wf_sum = r8vec_sum ( nf, wf );
  printf ( "\n" );
  printf ( "  Sum(WF) = %g\n", wf_sum );

  free ( wf );
  free ( x );
  free ( xf );

  return;
#  undef N
}
/******************************************************************************/

void test02 ( int m )

/******************************************************************************/
/*
  Purpose:

    TEST02 seeks Fekete points in [-1,+1].

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 April 2014

  Author:

    John Burkardt

  Reference:

    L Bos, N Levenberg,
    On the calculation of approximate Fekete points: the univariate case,
    Electronic Transactions on Numerical Analysis,
    Volume 30, pages 377-397, 2008.

  Parameters:

    Input, int M, the dimension of the polynomial space.
*/
{
# define N 5001

  double a;
  double b;
  int n = N;
  int nf;
  double *wf;
  double wf_sum;
  double *x;
  double *xf;

  a = -1.0;
  b = +1.0;
  x = r8vec_linspace_new ( n, a, b );

  printf ( "\n" );
  printf ( "TEST02:\n" );
  printf ( "  Seek Fekete points in [%f,%f]\n", a, b );
  printf ( "  using %d equally spaced sample points\n", n );
  printf ( "  for polynomials of degree M = %d\n", m );
  printf ( "  with the Chebyshev basis.\n" );

  wf = ( double * ) malloc ( m * sizeof ( double ) );
  xf = ( double * ) malloc ( m * sizeof ( double ) );
  line_fekete_chebyshev ( m, a, b, n, x, &nf, xf, wf );

  printf ( "\n" );
  printf ( "  NF = %d\n", nf );
  r8vec_print ( nf, xf, "  Estimated Fekete points XF:" );
  wf_sum = r8vec_sum ( nf, wf );
  printf ( "\n" );
  printf ( "  Sum(WF) = %g\n", wf_sum );

  free ( wf );
  free ( x );
  free ( xf );

  return;
#  undef N
}
/******************************************************************************/

void test03 ( int m )

/******************************************************************************/
/*
  Purpose:

    TEST03 seeks Fekete points in [-1,+1].

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 April 2014

  Author:

    John Burkardt

  Parameters:

    Input, int M, the dimension of the polynomial space.
*/
{
# define N 5001

  double a;
  double b;
  int n = N;
  int nf;
  double *wf;
  double wf_sum;
  double *x;
  double *xf;

  a = -1.0;
  b = +1.0;
  x = r8vec_linspace_new ( n, a, b );

  printf ( "\n" );
  printf ( "TEST03:\n" );
  printf ( "  Seek Fekete points in [%f,%f]\n", a, b );
  printf ( "  using %d equally spaced sample points\n", n );
  printf ( "  for polynomials of degree M = %d\n", m );
  printf ( "  with the Legendre basis and uniform weight.\n" );

  wf = ( double * ) malloc ( m * sizeof ( double ) );
  xf = ( double * ) malloc ( m * sizeof ( double ) );
  line_fekete_legendre ( m, a, b, n, x, &nf, xf, wf );

  printf ( "\n" );
  printf ( "  NF = %d\n", nf );
  r8vec_print ( nf, xf, "  Estimated Fekete points XF:" );
  wf_sum = r8vec_sum ( nf, wf );
  printf ( "\n" );
  printf ( "  Sum(WF) = %g\n", wf_sum );

  free ( wf );
  free ( x );
  free ( xf );

  return;
#  undef N
}
