# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "line_ncc_rule.h"

int main ( );
void test01 ( );
void test02 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for LINE_NCC_RULE_PRB.

  Discussion:

    LINE_NCC_RULE_PRB tests the LINE_NCC_RULE library.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    10 April 2014

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "LINE_NCC_RULE_PRB\n" );
  printf ( "  C version:\n" );
  printf ( "  Test the LINE_NCC_RULE library.\n" );

  test01 ( );
  test02 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "LINE_NCC_RULE_PRB\n" );
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

    TEST01 computes and prints NCC rules.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    10 April 2014

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  int i;
  int n;
  double *w;
  double w_sum;
  double *x;

  a = -1.0;
  b = +1.0;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  LINE_NCC_RULE computes the Newton-Cotes (closed) rule\n" );
  printf ( "  using N equally spaced points for an interval [A,B].\n" );

  for ( n = 1; n <= 12; n++ )
  {
    x = ( double * ) malloc ( n * sizeof ( double ) );
    w = ( double * ) malloc ( n * sizeof ( double ) );

    line_ncc_rule ( n, a, b, x, w );
    printf ( "\n" );
    printf ( "  Newton-Cotes (Closed) Rule #%d\n", n );
    printf ( "   I       X(I)            W(I)\n" );
    printf ( "\n" );
    w_sum = 0.0;
    for ( i = 0; i < n; i++ )
    {
      printf ( "  %2d  %14.6g  %14.6g\n", i, x[i], w[i] );
      w_sum = w_sum + fabs ( w[i] );
    }
    printf ( "        Sum(|W)|) =  %14.6g\n", w_sum );

    free ( x );
    free ( w );
  }
  return;
}
/******************************************************************************/

void test02 ( )

/******************************************************************************/
/*
  Purpose:

    TEST02 uses NCC rules to estimate the integral of exp(x) from 0 to 1.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    10 April 2014

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double error;
  double exact;
  int i;
  int n;
  double q;
  double *w;
  double *x;

  a =  0.0;
  b = +1.0;
  exact = exp ( b ) - exp ( a );

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  Use a sequence of NCC rules to compute an estimate Q\n" );
  printf ( "  of the integral:\n" );
  printf ( "    I = integral ( 0 <= x <= 1 ) exp(x) dx.\n" );
  printf ( "  The exact value is:\n" );
  printf ( "    I = %g\n", exact );

  printf ( "\n" );
  printf ( "   N       Q             |Q-I|\n" );
  printf ( "\n" );

  for ( n = 1; n <= 22; n++ )
  {
    x = ( double * ) malloc ( n * sizeof ( double ) );
    w = ( double * ) malloc ( n * sizeof ( double ) );

    line_ncc_rule ( n, a, b, x, w );

    q = 0.0;
    for ( i = 0; i < n; i++ )
    {
      q = q + w[i] * exp ( x[i] );
    }
    error = fabs ( exact - q );
    printf ( "  %2d  %14.6g  %14.6g\n", n, q, error );

    free ( x );
    free ( w );
  }
  return;
}