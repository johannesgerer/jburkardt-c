# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "beta_nc.h"

int main ( void );
void test01 ( void );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for BETA_NC_PRB.

  Discussion:

    BETA_NC_PRB tests the BETA_NC library.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    21 November 2010

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "BETA_NC_PRB:\n" );
  printf ( "  C version\n" );
  printf ( "  Test the BETA_NC library.\n" );

  test01 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "BETA_NC_PRB:\n" );
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

    TEST01 demonstrates the use of BETA_NONCENTRAL_CDF.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    21 November 2010

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double error_max;
  double fx;
  double fx2;
  double lambda;
  int n_data;
  double x;

  error_max = 1.0E-10;

  printf ( "\n" );
  printf ( "TEST01:\n" );
  printf ( "  BETA_NONCENTRAL_CDF computes the noncentral incomplete \n" );
  printf ( "  Beta function.\n" );
  printf ( "  Compare to tabulated values.\n" );
  printf ( "\n" );
  printf ( "      A      B     LAMBDA        X      " );
  printf ( "    FX                        FX2\n" );
  printf ( "                                        " );
  printf ( "    (Tabulated)               (computed)          DIFF\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    beta_noncentral_cdf_values ( &n_data, &a, &b, &lambda, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = beta_noncentral_cdf ( a, b, lambda, x, error_max );

    printf ( "  %5.2f  %5.2f  %7.3f  %10.4f  %24.16f  %24.16f  %10.4e\n",
      a, b, lambda, x, fx, fx2, r8_abs ( fx - fx2 ) );
  }

  return;
}
