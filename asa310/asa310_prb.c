# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "asa310.h"

int main ( );
void test01 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for ASA310_PRB.

  Discussion:

    ASA310_PRB tests the ASA310 library.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    20 November 2010

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "ASA310_PRB:\n" );
  printf ( "  C version\n" );
  printf ( "  Test the ASA310 library.\n" );

  test01 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "ASA310_PRB:\n" );
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

    TEST01 demonstrates the use of NCBETA.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    20 November 2010

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double errmax = 1.0E-10;
  double fx;
  double fx2;
  int ifault;
  double lambda;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST01:\n" );
  printf ( "  NCBETA computes the noncentral incomplete Beta function.\n" );
  printf ( "  Compare to tabulated values.\n" );
  printf ( "\n" );
  printf ( "      A        B     LAMBDA        X      " );
  printf ( "    FX                        FX2\n" );
  printf ( "                                          " );
  printf ( "    (Tabulated)               (NCBETA)            DIFF\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    beta_noncentral_cdf_values ( &n_data, &a, &b, &lambda, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = ncbeta ( a, b, lambda, x, errmax, &ifault );

    printf ( "  %7.2f  %7.2f  %7.3f  %10.4f  %24.16f  %24.16f  %10.4e\n",
      a, b, lambda, x, fx, fx2, fabs ( fx - fx2 ) );
  }

  return;
}
