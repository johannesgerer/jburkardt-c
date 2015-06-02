# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "asa109.h"

int main ( );
void test01 ( );
void test02 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for ASA109_PRB.

  Discussion:

    ASA109_PRB tests the ASA109 library.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    25 September 2014

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "ASA109_PRB:\n" );
  printf ( "  C version\n" );
  printf ( "  Test the ASA109 library.\n" );

  test01 ( );
  test02 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "ASA109_PRB:\n" );
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

    TEST01 demonstrates the use of XINBTA.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    28 April 2013

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double beta_log;
  double fx;
  int ifault;
  int n_data;
  double x;
  double x2;

  printf ( "\n" );
  printf ( "TEST01:\n" );
  printf ( "  XINBTA inverts the incomplete Beta function.\n" );
  printf ( "  Given CDF, it computes an X.\n" );
  printf ( "\n" );
  printf ( "           A           B           CDF    " );
  printf ( "    X                         X\n" );
  printf ( "                                          " );
  printf ( "    (Tabulated)               (XINBTA)            DIFF\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    beta_inc_values ( &n_data, &a, &b, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    beta_log = lgamma ( a )
             + lgamma ( b )
             - lgamma ( a + b );

    x2 = xinbta ( a, b, beta_log, fx, &ifault );

    printf ( "  %10.4f  %10.4f  %10.4f  %24.16g  %24.16g  %10.4e\n",
      a, b, fx, x, x2, fabs ( x - x2 ) );
  }

  return;
}
/******************************************************************************/

void test02 ( )

/******************************************************************************/
/*
  Purpose:

    TEST02 tests BETA_INC_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    25 September 2014

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double beta_log;
  double e;
  double fx;
  double fx2;
  int ifault;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST02:\n" );
  printf ( "  BETA_INC_VALUES stores values of\n" );
  printf ( "  the incomplete Beta function.\n" );
  printf ( "\n" );
  printf ( "      A            B            X            BETA_INC(A,B)(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    beta_inc_values ( &n_data, &a, &b, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    beta_log = lgamma ( a )
             + lgamma ( b )
             - lgamma ( a + b );

    fx2 = betain ( x, a, b, beta_log, &ifault );

    e = fabs ( fx - fx2 );

    printf ( "  %12f  %12f  %12f  %24.16g  %24.16g  %8.4e\n", a, b, x, fx, fx2, e );
  }
  return;
}
