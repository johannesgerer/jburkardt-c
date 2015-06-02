# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "toms179.h"

int main ( void );
void test01 ( void );
void test02 ( void );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for TOMS179_PRB.

  Discussion:

    TOMS179_PRB tests the TOMS179 library.

  Modified:

    30 January 2008

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "TOMS179_PRB:\n" );
  printf ( "  C version\n" );
  printf ( "  Test the TOMS179 library.\n" );

  test01 ( );
  test02 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "TOMS179_PRB:\n" );
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

    30 January 2008

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
  printf ( "          X                     FX                        FX2\n" );
  printf ( "                                (Tabulated)               (ALOGAM)                DIFF\n" );
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

    printf ( "  %24.16e  %24.16e  %24.16e  %10.4e\n",
      x, fx, fx2, r8_abs ( fx - fx2 ) );
  }

  return;
}
/******************************************************************************/

void test02 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST02 demonstrates the use of MDBETA.

  Modified:

    30 January 2008

  Author:

    John Burkardt
*/
{
  double fx;
  double fx2;
  int ier;
  int n_data;
  double p;
  double q;
  double x;

  printf ( "\n" );
  printf ( "TEST02:\n" );
  printf ( "  MDBETA estimates the value of th modified Beta function.\n" );
  printf ( "  Compare with tabulated values.\n" );
  printf ( "\n" );
  printf ( "         X         P         Q         Beta                       Beta                  DIFF\n" );
  printf ( "                                       (Tabulated)                (MDBETA)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    beta_cdf_values ( &n_data, &p, &q, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = mdbeta ( x, p, q, &ier );

    printf ( "  %8.4f  %8.4f  %8.4f  %24.16e  %24.16e %10.4e\n",
      x, p, q, fx, fx2, r8_abs ( fx - fx2 ) );
  }

  return;
}
