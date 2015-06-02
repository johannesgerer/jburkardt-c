# include <stdlib.h>
# include <stdio.h>

# include "test_min.h"

int main ( void );
void test01 ( void );
void test02 ( void );
void test03 ( void );
void test04 ( void );
void test05 ( void );
void test06 ( void );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for TEST_MIN_PRB.

  Discussion:

    TEST_MIN_PRB tests the TEST_MIN library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    03 February 2012

  Author:

    John Burkardt
*/
{
  timestamp (  );
  printf ( "\n" );
  printf ( "TEST_MIN_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the TEST_MIN library.\n" );

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
  printf ( "TEST_MIN_PRB\n" );
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

    TEST01 prints the title of each problem.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    03 February 2012

  Author:

    John Burkardt
*/
{
  int problem_num;
  int problem;
  char title[50];

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  For each problem, print the title.\n" );
/*
  Get the number of problems.
*/
  problem_num = p00_problem_num ( );

  printf ( "\n" );
  printf ( "   Problem Title\n" );
  printf ( "\n" );

  for ( problem = 1; problem <= problem_num; problem++ )
  {
    p00_title ( problem, title );

    printf ( "  %2d  %s\n", problem, title );
  }

  return;
}
/******************************************************************************/

void test02 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST02 evaluates the objective function at each starting point.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    03 February 2012

  Author:

    John Burkardt
*/
{
  double f_sol;
  double f_start;
  int know;
  int problem_num;
  int problem;
  char title[50];
  double x;

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  For each problem, evaluate the function\n" );
  printf ( "  at the starting point and the solution.\n" );
/*
  Get the number of problems.
*/
  problem_num = p00_problem_num ( );

  for ( problem = 1; problem <= problem_num; problem++ )
  {
    p00_title ( problem, title );

    printf ( "\n" );
    printf ( "  Problem %d\n", problem );
    printf ( "  %s\n", title );
    printf ( "\n" );

    x = p00_start ( problem );

    f_start = p00_f ( problem, x );

    printf ( "    F(X_START) = %g\n", f_start );

    p00_sol ( problem, &know, &x );

    if ( 0 < know )
    {
      f_sol = p00_f ( problem, x );
      printf ( "    F(X_SOL) = %g\n", f_sol );
    }
  }

  return;
}
/******************************************************************************/

void test03 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST03 compares the exact and approximate first derivatives.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    03 February 2012

  Author:

    John Burkardt
*/
{
  double f1;
  double f1_dif;
  int problem_num;
  int problem;
  char title[50];
  double x;

  printf ( "\n" );
  printf ( "TEST03\n" );
  printf ( "  For each problem, compare the exact and\n" );
  printf ( "  approximate gradients at the starting point.\n" );
/*
  Get the number of problems.
*/
  problem_num = p00_problem_num ( );

  for ( problem = 1; problem <= problem_num; problem++ )
  {
    p00_title ( problem, title );

    printf ( "\n" );
    printf ( "  Problem %d\n", problem );
    printf ( "  %s\n", title );

    x = p00_start ( problem );

    f1 = p00_f1 ( problem, x );

    f1_dif = p00_f1_dif ( problem, x );

    printf ( "\n" );
    printf ( "  X\n" );
    printf ( "  %g\n", x );
    printf ( "  F'(X) (exact)\n" );
    printf ( "  %g\n", f1 );
    printf ( "  F'(X) (difference)\n" );
    printf ( "  %g\n", f1_dif );
  }

  return;
}
/******************************************************************************/

void test04 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST04 compares the exact and approximate second derivatives.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    03 February 2012

  Author:

    John Burkardt
*/
{
  double f2;
  double f2_dif;
  int problem_num;
  int problem;
  char title[50];
  double x;

  printf ( "\n" );
  printf ( "TEST04\n" );
  printf ( "  For each problem, compare the exact and\n" );
  printf ( "  approximate second derivatives at the starting point.\n" );
/*
  Get the number of problems.
*/
  problem_num = p00_problem_num ( );

  for ( problem = 1; problem <= problem_num; problem++ )
  {
    p00_title ( problem, title );

    printf ( "\n" );
    printf ( "  Problem %d\n", problem );
    printf ( "  %s\n", title );

    x = p00_start ( problem );

    printf ( "\n" );
    printf ( "  X:\n" );
    printf ( "  %g\n", x );

    f2 = p00_f2 ( problem, x );

    printf ( "  F\"(X) (exact):\n" );
    printf ( "  %g\n", f2 );

    f2_dif = p00_f2_dif ( problem, x );

    printf ( "  F\"(X) (difference):\n" );
    printf ( "  %g\n", f2_dif );
  }

  return;
}
/******************************************************************************/

void test05 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST05 carries out a simple bisection method.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    03 February 2012

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double c;
  double d;
  double e;
  double fa;
  double fb;
  double fc;
  double fd;
  double fe;
  int i;
  int max_step = 10;
  int problem_num;
  int problem;
  char title[50];

  printf ( "\n" );
  printf ( "TEST05\n" );
  printf ( "  For each problem, take a few steps of \n" );
  printf ( "  the bisection method.\n" );
/*
  Get the number of problems.
*/
  problem_num = p00_problem_num ( );

  for ( problem = 1; problem <= problem_num; problem++ )
  {
    p00_title ( problem, title );

    printf ( "\n" );
    printf ( "  Problem %d\n", problem );
    printf ( "  %s\n", title );

    p00_interval ( problem, &a, &c );
    b = 0.5 * ( a + c );
    fa = p00_f ( problem, a );
    fc = p00_f ( problem, c );
    fb = p00_f ( problem, b );

    i = 0;
    printf ( "\n" );
    printf ( "  %d\n", i );
    printf ( "  X:  %10g  %10g  %10g\n", a, b, c );
    printf ( "  F:  %10g  %10g  %10g\n", fa, fb, fc );

    for ( i = 1; i <= max_step; i++ )
    {
      d = 0.5 * ( a + b );
      fd = p00_f ( problem, d );

      e = 0.5 * ( b + c );
      fe = p00_f ( problem, e );

      if ( fd <= fb )
      {
        c = b;
        fc = fb;
        b = d;
        fb = fd;
      }
      else if ( fe <= fb )
      {
        a = b;
        fa = fb;
        b = e;
        fb = fe;
      }
      else
      {
        a = d;
        fa = fd;
        c = e;
        fc = fe;
      }
    printf ( "  %d\n", i );
    printf ( "  X:  %10g  %10g  %10g\n", a, b, c );
    printf ( "  F:  %10g  %10g  %10g\n", fa, fb, fc );
    }
  }

  return;
}
/******************************************************************************/

void test06 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST06 carries out a version of Brent's derivative-free minimizer.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    03 February 2012

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double fa;
  double fb;
  double fx;
  int problem_num;
  int problem;
  char title[50];
  double tol = 0.000001;
  double x;

  printf ( "\n" );
  printf ( "TEST06\n" );
  printf ( "  For each problem, use Brent's method.\n" );
/*
  Get the number of problems.
*/
  problem_num = p00_problem_num ( );

  for ( problem = 1; problem <= problem_num; problem++ )
  {
    p00_title ( problem, title );

    printf ( "\n" );
    printf ( "  Problem %d\n", problem );
    printf ( "  %s\n", title );

    p00_interval ( problem, &a, &b );

    fa = p00_f ( problem, a );
    fb = p00_f ( problem, b );

    printf ( "\n" );
    printf ( "  Initial interval [A,B]:\n" );
    printf ( "\n" );
    printf ( "   A,       B:  %16g                      %16g\n", a, b );
    printf ( "  FA,      FB:  %16g                      %16g\n", fa, fb );

    x = p00_fmin ( &a, &b, problem, tol );

    fa = p00_f ( problem, a );
    fb = p00_f ( problem, b );
    fx = p00_f ( problem, x );

    printf ( "\n" );
    printf ( "  Final interval [A,X*,B]:\n" );
    printf ( "\n" );
    printf ( "   A,  X*,  B:  %16g  %16g  %16g\n", a, x, b );
    printf ( "  FA, FX*, FB:  %16g  %16g  %16g\n", fa, fx, fb );
  }

  return;
}
