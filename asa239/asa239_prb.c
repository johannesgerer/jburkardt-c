# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "asa239.h"

int main ( void );
void test01 ( void );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for ASA239_PRB.

  Discussion:

    ASA239_PRB calls the ASA239 routines.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    20 January 2008

  Author:

    John Burkardt
*/
{
  timestamp ( );

  printf ( "\n" );
  printf ( "ASA239_PRB:\n" );
  printf ( "  C version\n" );
  printf ( "  Test the ASA239 library.\n" );

  test01 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "ASA239_PRB:\n" );
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

    TEST01 demonstrates the use of GAMMAD.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    13 November 2010

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
  printf ( "  GAMMAD computes the incomplete Gamma function.\n" );
  printf ( "  Compare to tabulated values.\n" );
  printf ( "\n" );
  printf ( "             A             X      " );
  printf ( "    FX                        FX2\n" );
  printf ( "                                  " );
  printf ( "    (Tabulated)               (GAMMAD)            DIFF\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    gamma_inc_values ( &n_data, &a, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = gammad ( x, a, &ifault );

    printf ( "  %12.4f  %12.4f  %24.16g  %24.16g  %10.4e\n",
      a, x, fx, fx2, fabs ( fx - fx2 ) );
  }

  return;
}
