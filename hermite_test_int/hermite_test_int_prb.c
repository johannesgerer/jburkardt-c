# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <string.h>

# include "hermite_test_int.h"

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );
void test04 ( );
void test05 ( );
void test06 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for HERMITE_TEST_INT_PRB.

  Discussion:

    HERMITE_TEST_INT_PRB tests the HERMITE_TEST_INT library.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 September 2012

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "HERMITE_TEST_INT_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the HERMITE_TEST_INT library.\n" );

  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
  test05 ( );
  test06 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "HERMITE_TEST_INT_PRB\n" );
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

    TEST01 tests P00_PROBLEM_NUM and P00_TITLE.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 September 2012

  Author:

    John Burkardt
*/
{
  int problem;
  int problem_num;
  char *title;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  P00_PROBLEM_NUM returns the number of problems.\n" );
  printf ( "  P00_TITLE returns the title of a problem.\n" );

  problem_num = p00_problem_num ( );

  printf ( "\n" );
  printf ( "  P00_PROBLEM_NUM: number of problems is %d\n", problem_num );
  printf ( "\n" );
  printf ( "   Problem       Title\n" );
  printf ( "\n" );

  for ( problem = 1; problem <= problem_num; problem++ )
  {
    title = p00_title ( problem );

    printf ( "  %8d  \"%s\".\n", problem, title );

    free ( title );
  }

  return;
}
/******************************************************************************/

void test02 ( )

/******************************************************************************/
/*
  Purpose:

    TEST02 tests P00_EXACT.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 September 2012

  Author:

    John Burkardt
*/
{
  double exact;
  int m;
  int problem;
  int problem_num;

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  P00_EXACT returns the \"exact\" integral.\n" );

  problem_num = p00_problem_num ( );

  m = 4;
  p06_param ( 'S', 'M', &m );

  printf ( "\n" );
  printf ( "   Problem       EXACT\n" );
  printf ( "\n" );

  for ( problem = 1; problem <= problem_num; problem++ )
  {
    exact = p00_exact ( problem );

    printf ( "  %8d  %24.16g\n", problem, exact );
  }

  return;
}
/******************************************************************************/

void test03 ( )

/******************************************************************************/
/*
  Purpose:

    TEST03 tests P00_GAUSS_HERMITE.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 September 2012

  Author:

    John Burkardt
*/
{
  double error;
  double estimate;
  double exact;
  int m;
  int order;
  int order_log;
  int problem;
  int problem_num;

  printf ( "\n" );
  printf ( "TEST03\n" );
  printf ( "  P00_GAUSS_HERMITE applies a Gauss-Hermite rule\n" );
  printf ( "  to estimate an integral on (-oo,+oo).\n" );

  problem_num = p00_problem_num ( );

  m = 4;
  p06_param ( 'S', 'M', &m );

  printf ( "\n" );
  printf ( "   Problem     Order          Estimate        Exact          Error\n" );

  for ( problem = 1; problem <= problem_num; problem++ )
  {
    exact = p00_exact ( problem );

    order = 1;

    printf ( "\n" );

    for ( order_log = 0; order_log <= 6; order_log++ )
    {
      estimate = p00_gauss_hermite ( problem, order );

      error = r8_abs ( exact - estimate );

      printf ( "  %8d  %8d  %14.6g  %14.6g  %14.6g\n", 
        problem, order, estimate, exact, error );

      order = order * 2;
    }
  }
  return;
}
/******************************************************************************/

void test04 ( )

/******************************************************************************/
/*
  Purpose:

    TEST04 tests P00_TURING.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 September 2012

  Author:

    John Burkardt
*/
{
  double error;
  double estimate;
  double exact;
  double h;
  int m;
  int n;
  int order_log;
  int problem;
  int problem_num;
  int test;
  double tol;

  printf ( "\n" );
  printf ( "TEST04\n" );
  printf ( "  P00_TURING applies a Turing procedure\n" );
  printf ( "  to estimate an integral on (-oo,+oo).\n" );

  problem_num = p00_problem_num ( );

  m = 4;
  p06_param ( 'S', 'M', &m );

  for ( test = 1; test <= 2; test++ )
  {
    if ( test == 1 )
    {
      tol = 1.0E-4;
    }
    else if ( test == 2 )
    {
      tol = 1.0E-07;
    }
    printf ( "\n" );
    printf ( "  Using a tolerance of TOL = %g\n", tol );
    printf ( "\n" );
    printf ( 
      "   Problem     Order          Estimate        Exact          Error\n" );

    for ( problem = 1; problem <= problem_num; problem++ )
    {
      exact = p00_exact ( problem );

      h = 1.0;

      printf ( "\n" );

      for ( order_log = 0; order_log <= 6; order_log++ )
      {
        estimate = p00_turing ( problem, h, tol, &n );

        error = r8_abs ( exact - estimate );

        printf ( "  %8d  %10g  %8d  %14.6g  %14.6g  %14.6g\n",
          problem, h, n, estimate, exact, error );

        h = h / 2.0;
      }
    }
  }
  return;
}
/******************************************************************************/

void test05 ( )

/******************************************************************************/
/*
  Purpose:

    TEST05 tests P00_GAUSS_HERMITE against the polynomials.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 September 2012

  Author:

    John Burkardt
*/
{
  double error;
  double estimate;
  double exact;
  int m;
  int order;
  int order_log;
  int problem;

  printf ( "\n" );
  printf ( "TEST05\n" );
  printf ( "  P00_GAUSS_HERMITE applies a Gauss-Hermite rule to\n" );
  printf ( "  estimate the integral x^m exp(-x*x) over (-oo,+oo).\n" );

  problem = 6;

  printf ( "\n" );
  printf ( "         M     Order      Estimate        Exact           Error\n" );

  for ( m = 0; m <= 6; m++ )
  {
    p06_param ( 'S', 'M', &m );

    exact = p00_exact ( problem );

    printf ( "\n" );

    for ( order = 1; order <= 3 + ( m / 2 ); order++ )
    {
      estimate = p00_gauss_hermite ( problem, order );

      error = r8_abs ( exact - estimate );

      printf ( "  %8d  %8d  %14g  %14g  %14g\n", 
        m, order, estimate, exact, error );
    }
  }

  return;
}
/******************************************************************************/

void test06 ( )

/******************************************************************************/
/*
  Purpose:

    TEST06 tests P00_MONTE_CARLO.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 September 2012

  Author:

    John Burkardt
*/
{
  double error;
  double estimate;
  double exact;
  int m;
  int order;
  int order_log;
  int problem;
  int problem_num;

  printf ( "\n" );
  printf ( "TEST06\n" );
  printf ( "  P00_MONTE_CARLO uses a weighted form of the Monte Carlo method\n" );
  printf ( "  to estimate a Hermite integral on (-oo,+oo).\n" );

  problem_num = p00_problem_num ( );

  m = 4;
  p06_param ( 'S', 'M', &m );

  printf ( "\n" );
  printf ( 
    "   Problem     Order          Estimate        Exact          Error\n" );

  for ( problem = 1; problem <= problem_num; problem++ )
  {
    exact = p00_exact ( problem );

    order = 128;

    printf ( "\n" );

    for ( order_log = 0; order_log <= 6; order_log++ )
    {
      estimate = p00_monte_carlo ( problem, order );

      error = r8_abs ( exact - estimate );

      printf ( "  %8d  %8d  %14.6g  %14.6g  %14.6g\n", 
        problem, order, estimate, exact, error );

      order = order * 4;
    }
  }
  return;
}
