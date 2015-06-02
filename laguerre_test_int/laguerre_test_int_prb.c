# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <string.h>

# include "laguerre_test_int.h"

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );
void test04 ( );
void test05 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for LAGUERRE_TEST_INT_PRB.

  Discussion:

    LAGUERRE_TEST_INT_PRB tests the LAGUERRE_TEST_INT library.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    14 September 2012

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "LAGUERRE_TEST_INT_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the LAGUERRE_TEST_INT library.\n" );

  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
  test05 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "LAGUERRE_TEST_INT_PRB\n" );
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

    14 September 2012

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

    TEST02 tests P00_ALPHA and P00_EXACT.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    14 September 2012

  Author:

    John Burkardt
*/
{
  double alpha;
  double exact;
  int problem;
  int problem_num;

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  P00_ALPHA returns the lower limit of integration.\n" );
  printf ( "  P00_EXACT returns the \"exact\" integral.\n" );

  problem_num = p00_problem_num ( );

  printf ( "\n" );
  printf ( "   Problem       ALPHA           EXACT\n" );
  printf ( "\n" );

  for ( problem = 1; problem <= problem_num; problem++ )
  {
    alpha = p00_alpha ( problem );

    exact = p00_exact ( problem );

    printf ( "  %8d  %14g  %24.16g\n", problem, alpha, exact );
  }

  return;
}
/******************************************************************************/

void test03 ( )

/******************************************************************************/
/*
  Purpose:

    TEST03 tests P00_GAUSS_LAGUERRE.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    14 September 2012

  Author:

    John Burkardt
*/
{
  double error;
  double estimate;
  double exact;
  int order;
  int order_log;
  int problem;
  int problem_num;

  printf ( "\n" );
  printf ( "TEST03\n" );
  printf ( "  P00_GAUSS_LAGUERRE applies a Gauss-Laguerre rule\n" );
  printf ( "  to estimate an integral on [ALPHA,+oo).\n" );

  problem_num = p00_problem_num ( );

  printf ( "\n" );
  printf ( "                              Exact\n" );
  printf ( "   Problem     Order          Estimate        Error\n" );

  for ( problem = 1; problem <= problem_num; problem++ )
  {
    exact = p00_exact ( problem );

    order = 1;

    printf ( "\n" );
    printf ( "  %8d            %14.6g\n", problem, exact );

    for ( order_log = 0; order_log <= 6; order_log++ )
    {
      estimate = p00_gauss_laguerre ( problem, order );

      error = r8_abs ( exact - estimate );

      printf ( "          %8d  %14.6g  %14.6g\n", order, estimate, error );

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

    TEST04 tests P00_EXP_TRANSFORM.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    14 September 2012

  Author:

    John Burkardt
*/
{
  double error;
  double estimate;
  double exact;
  int order;
  int order_log;
  int problem;
  int problem_num;

  printf ( "\n" );
  printf ( "TEST04\n" );
  printf ( "  P00_EXP_TRANSFORM applies an exponential transform\n" );
  printf ( "  to estimate an integral on [ALPHA,+oo)\n" );
  printf ( "  as a transformed integral on (0,exp(-ALPHA)],\n" );
  printf ( "  and applying a Gauss-Legendre rule.\n" );

  problem_num = p00_problem_num ( );

  printf ( "\n" );
  printf ( "                              Exact\n" );
  printf ( "   Problem     Order          Estimate        Error\n" );

  for ( problem = 1; problem <= problem_num; problem++ )
  {
    exact = p00_exact ( problem );

    order = 1;

    printf ( "\n" );
    printf ( "  %8d            %14.6g\n", problem, exact );

    for ( order_log = 0; order_log <= 6; order_log++ )
    {
      estimate = p00_exp_transform ( problem, order );

      error = r8_abs ( exact - estimate );

      printf ( "          %8d  %14.6g  %14.6g\n", order, estimate, error );

      order = order * 2;
    }
  }
  return;
}
/******************************************************************************/

void test05 ( )

/******************************************************************************/
/*
  Purpose:

    TEST05 tests P00_RAT_TRANSFORM.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    14 September 2012

  Author:

    John Burkardt
*/
{
  double error;
  double estimate;
  double exact;
  int order;
  int order_log;
  int problem;
  int problem_num;

  printf ( "\n" );
  printf ( "TEST05\n" );
  printf ( "  P00_RAT_TRANSFORM applies a rational transform\n" );
  printf ( "  to estimate an integral on [ALPHA,+oo)\n" );
  printf ( "  as a transformed integral on (0,1/(1+ALPHA)],\n" );
  printf ( "  and applying a Gauss-Legendre rule.\n" );

  problem_num = p00_problem_num ( );

  printf ( "\n" );
  printf ( "                              Exact\n" );
  printf ( "   Problem     Order          Estimate        Error\n" );

  for ( problem = 1; problem <= problem_num; problem++ )
  {
    exact = p00_exact ( problem );

    order = 1;

    printf ( "\n" );
    printf ( "  %8d            %14.6g\n", problem, exact );

    for ( order_log = 0; order_log <= 6; order_log++ )
    {
      estimate = p00_rat_transform ( problem, order );

      error = r8_abs ( exact - estimate );

      printf ( "          %8d  %14.6g  %14.6g\n", order, estimate, error );

      order = order * 2;
    }
  }
  return;
}
