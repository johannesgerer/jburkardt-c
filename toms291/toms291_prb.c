# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "toms291.h"

int main ( void );
void test01 ( void );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for TOMS291_PRB.

  Discussion:

    TOMS291_PRB calls the TOMS291 routines.

  Modified:

    30 October 2010

  Author:

    John Burkardt
*/
{
  timestamp ( );

  printf ( "\n" );
  printf ( "TOMS291_PRB:\n" );
  printf ( "  C version\n" );
  printf ( "  Test the TOMS291 library.\n" );

  test01 ( );

  printf ( "\n" );
  printf ( "TOMS291_PRB:\n" );
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

    TEST01 demonstrates the use of ALOGAM.

  Modified:

    30 October 2010

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
  printf ( "  ALOGAM computes the logarithm of the \n" );
  printf ( "  Gamma function.  We compare the result\n" );
  printf ( "  to tabulated values.\n" );
  printf ( "\n" );
  printf ( "          X                     " );
  printf ( "FX                        FX2\n" );
  printf ( "                                " );
  printf ( "(Tabulated)               (ALOGAM)                DIFF\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    gamma_log_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = alogam ( x, &ifault );

    printf ( "  %24.16g  %24.16g  %24.16g  %10.4e\n",
      x, fx, fx2, r8_abs ( fx - fx2 ) );
  }

  return;
}
