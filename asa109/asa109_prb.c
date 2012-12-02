# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "asa109.h"

int main ( void );
void test01 ( void );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for ASA109_PRB.

  Discussion:

    ASA109_PRB calls the ASA109 routines.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    05 November 2010

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

void test01 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST01 demonstrates the use of XINBTA.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    05 November 2010

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

    beta_log = alngam ( a, &ifault )
             + alngam ( b, &ifault )
             - alngam ( a + b, &ifault );

    x2 = xinbta ( a, b, beta_log, fx, &ifault );

    printf ( "  %10.4f  %10.4f  %10.4f  %24.16g  %24.16g  %10.4e\n",
      a, b, fx, x, x2, r8_abs ( x - x2 ) );
  }

  return;
}
