# include <stdlib.h>
# include <stdio.h>
# include <string.h>

# include "test_opt_con.h";

int main ( void );
void test01 ( void );
void test02 ( void );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for TEST_OPT_CON_PRB.

  Discussion:

    TEST_OPT_CON_PRB tests the TEST_OPT_CON library.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 February 2012

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "TEST_OPT_CON_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the TEST_OPT_CON library.\n" );

  test01 ( );
  test02 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "TEST_OPT_CON_PRB\n" );
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

    TEST01 simply prints the title of each problem.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 February 2012

  Author:

    John Burkardt
*/
{
  int problem_num;
  int problem;
  char title[100];

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  For each problem, print the title.\n" );
/*
  Get the number of problems.
*/
  problem_num = p00_problem_num ( );

  printf ( "\n" );
  printf ( "  Problem    Title\n" );
  printf ( "\n" );

  for ( problem = 1; problem <= problem_num; problem++ )
  {
    p00_title ( problem, title );

    printf ( "  %6d  %s\n", problem, title );
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

    17 February 2012

  Author:

    John Burkardt
*/
{
  double *a;
  double *b;
  double *f;
  double *fs;
  int i;
  int know;
  int m;
  int n = 100000;
  int problem;
  int problem_num;
  int seed;
  char title[100];
  double *x;
  double *xs;

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  For each problem, evaluate the function at many points.\n" );
  printf ( "  Number of sample points = %d\n", n );
/*
  Get the number of problems.
*/
  problem_num = p00_problem_num ( );

  for ( problem = 1; problem <= problem_num; problem++ )
  {
    printf ( "\n" );
    printf ( "  Problem %d\n", problem );

    p00_title ( problem, title );

    printf ( "  %s\n", title );

    m = p00_m ( problem );

    printf ( "  M =     %d\n", m );

    a = ( double * ) malloc ( m * sizeof ( double ) );
    b = ( double * ) malloc ( m * sizeof ( double ) );
 
    p00_ab ( problem, m, a, b );

    printf ( "\n" );
    printf ( "    I      A(i)      B(i)\n" );
    printf ( "\n" );

    for ( i = 0; i < m; i++ )
    {
      printf ( "  %4d  %10g  %10g\n", i, a[i], b[i] );
    }

    seed = 123456789;
    x = r8col_uniform_new ( m, n, a, b, &seed );
    f = p00_f ( problem, m, n, x );

    printf ( "\n" );
    printf ( "  Max(F) = %g\n", r8vec_max ( n, f ) );
    printf ( "  Min(F) = %g\n", r8vec_min ( n, f ) );

    know = 0;
    xs = p00_sol ( problem, m, &know );
    if ( know != 0 )
    {
      fs = p00_f ( problem, m, 1, xs );
      printf ( "  F(X*)  = %g\n", fs[0] );
      free ( fs );
      free ( xs );
    }
    else
    {
      printf ( "  X* is not given.\n" );
    }

    free ( a );
    free ( b );
    free ( f );
    free ( x );
  }
  return;
}
