# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "asa121.h"

int main ( );
void test01 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for ASA121_PRB.

  Discussion:

    ASA121_PRB tests the ASA121 library.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    08 November 2010

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "ASA121_PRB:\n" );
  printf ( "  C version\n" );
  printf ( "  Test the ASA121 library.\n" );

  test01 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "ASA121_PRB:\n" );
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

    TEST01 demonstrates the use of TRIGAMMA.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    19 January 2008

  Author:

    John Burkardt
*/
{
  double fx;
  double fx2;
  int ifault;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST01:\n" );
  printf ( "  TRIGAMMA computes the trigamma function. \n" );
  printf ( "  We compare the result to tabulated values.\n" );
  printf ( "\n" );
  printf ( "          X                     " ); 
  printf ( "FX                        FX2\n" );
  printf ( "                                " );
  printf ( "(Tabulated)               (TRIGAMMA)                DIFF\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    trigamma_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = trigamma ( x, &ifault );

    printf ( "  %24.16f  %24.16f  %24.16f  %10.4g\n",
      x, fx, fx2, fabs ( fx - fx2 ) );
  }

  return;
}
