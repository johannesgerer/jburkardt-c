# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "asa032.h"

int main ( );
void test01 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for ASA032_PRB.

  Discussion:

    ASA032_PRB tests the ASA032 library.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    26 October 2010

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "ASA032_PRB:\n" );
  printf ( "  C version\n" );
  printf ( "  Test the ASA032 library.\n" );

  test01 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "ASA032_PRB:\n" );
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

    TEST01 demonstrates the use of ALNGAM.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 January 2008

  Author:

    John Burkardt
*/
{
  double a;
  double fx;
  double fx2;
  int ifault;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST01:\n" );
  printf ( "  GAMAIN computes the incomplete Gamma function.\n" );
  printf ( "  Compare the result to tabulated values.\n" );
  printf ( "\n" );
  printf ( "          A               X           " );
  printf ( "FX                        FX2\n" );
  printf ( "                                      " );
  printf ( "(Tabulated)               (GAMAIN)                DIFF\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    gamma_inc_values ( &n_data, &a, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = gamain ( x, a, &ifault );

    printf ( "  %12.8g  %12.8g  %24.16g  %24.16g  %10.4e\n", 
      a, x, fx, fx2, fabs ( fx - fx2 ) );
  }

  return;
}
