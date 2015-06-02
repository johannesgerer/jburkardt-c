# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>

# include "asa243.h"

int main ( );
void test01 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for ASA243_PRB.

  Discussion:

    ASA243_PRB tests the ASA243 library.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    13 November 2010

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "ASA243_PRB:\n" );
  printf ( "  C version\n" );
  printf ( "  Test the ASA243 library.\n" );

  test01 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "ASA243_PRB:\n" );
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

    TEST01 demonstrates the use of TNC.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    13 November 2010

  Author:

    John Burkardt
*/
{
  double delta;
  int df;
  double df_real;
  double fx;
  double fx2;
  int ifault;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST01:\n" );
  printf ( "  TNC computes the noncentral Student T\n" );
  printf ( "  Cumulative Density Function.\n" );
  printf ( "  Compare with tabulated values.\n" );
  printf ( "\n" );
  printf ( "        X         LAMBDA        DF     " );
  printf ( " CDF             CDF           DIFF\n" );
  printf ( "                                       " );
  printf ( " Tabulated       PRNCST\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    student_noncentral_cdf_values ( &n_data, &df, &delta, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    df_real = ( double ) ( df );

    fx2 = tnc ( x, df_real, delta, &ifault );

    printf ( "  %10.4f  %10.4f  %8d  %24.16f  %24.16f  %10.4e\n",
      x, delta, df, fx, fx2, fabs ( fx - fx2 ) );
  }

  return;
}
