# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "asa103.h"

int main ( void );
void test01 ( void );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for ASA103_PRB.

  Discussion:

    ASA103_PRB tests the ASA103 library.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    04 November 2010

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "ASA103_PRB:\n" );
  printf ( "  C version\n" );
  printf ( "  Test the ASA103 library.\n" );

  test01 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "ASA103_PRB:\n" );
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

    TEST01 demonstrates the use of DIGAMA.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    04 November 2010

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
  printf ( "  DIGAMA computes the Digama or Psi function. \n" );
  printf ( "  Compare the result to tabulated values.\n" );
  printf ( "\n" );
  printf ( "          X       " );
  printf ( "FX                        FX2\n" );
  printf ( "                  " );
  printf ( "(Tabulated)               (DIGAMA)                DIFF\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    psi_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = digama ( x, &ifault );

    printf ( "  %10.4f  %24.16g  %24.16g  %10.4e\n",
      x, fx, fx2, fabs ( fx - fx2 ) );
  }

  return;
}
