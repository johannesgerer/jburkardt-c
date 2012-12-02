# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>

# include "asa111.h"

int main ( void );
void test01 ( void );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    ASA111_PRB tests ASA111.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 November 2010

  Author:

    John Burkardt
*/
{
  printf ( "\n" );
  timestamp ( );

  printf ( "\n" );
  printf ( "ASA111_PRB:\n" );
  printf ( "  C version\n" );
  printf ( "  Test the ASA111 library.\n" );

  test01 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "ASA111_PRB:\n" );
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

    TEST01 compares PPND against tabulated values.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 November 2010

  Author:

    John Burkardt
*/
{
  double fx;
  int ifault;
  int n_data;
  double x;
  double x2;

  printf ( "\n" );
  printf ( "TEST01:\n" );
  printf ( "  PPND computes percentage points of the normal distribution.\n" );
  printf ( "  Compare against tabulated values.\n" );
  printf ( "\n" );
  printf ( "         CDF        X                           X  " );
  printf ( "                  DIFF\n" );
  printf ( "                 (tabulated)                   (PPND)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    normal_01_cdf_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    x2 = ppnd ( fx, &ifault );

    printf ( "  %10.4g  %24.16g  %24.16g  %10.4g\n",
      fx, x, x2, r8_abs ( x - x2 ) );
  }

  return;
}
