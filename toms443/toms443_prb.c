# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "toms443.h"

int main ( );
void test01 ( );
void test02 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for TOMS443_PRB.

  Discussion:

    TOMS443_PRB tests the TOMS443 library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 June 2014

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "TOMS443_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the TOMS443 library.\n" );

  test01 ( );
  test02 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "TOMS433_PRB\n" );
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

    TEST01 tests WEW_A

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 June 2014

  Author:

    John Burkardt
*/
{
  double en;
  int n_data;
  double w1;
  double w2;
  double x;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  Test WEW_A to evaluate\n" );
  printf ( "  Lambert's W function.\n" );
  printf ( "\n" );
  printf ( "          X             Exact             Computed      Error\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    lambert_w_values ( &n_data, &x, &w1 );

    if ( n_data <= 0 )
    {
      break;
    }

    if ( x == 0.0 )
    {
      w2 = 0.0;
    }
    else
    {
      w2 = wew_a ( x, &en );
    }

    printf ( "  %12.4f  %16.8g  %16.8g  %10.2e\n",
      x, w1, w2, fabs ( w1 - w2 ) );
  }

  return;
}
/******************************************************************************/

void test02 ( )

/******************************************************************************/
/*
  Purpose:

    TEST02 tests WEW_B

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 June 2014

  Author:

    John Burkardt
*/
{
  double en;
  int n_data;
  double w1;
  double w2;
  double x;

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  Test WEW_B to evaluate\n" );
  printf ( "  Lambert's W function.\n" );
  printf ( "\n" );
  printf ( "          X             Exact             Computed      Error\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    lambert_w_values ( &n_data, &x, &w1 );

    if ( n_data <= 0 )
    {
      break;
    }

    if ( x == 0.0 )
    {
      w2 = 0.0;
    }
    else
    {
      w2 = wew_b ( x, &en );
    }

    printf ( "  %12.4f  %16.8g  %16.8g  %10.2e\n",
      x, w1, w2, fabs ( w1 - w2 ) );
  }

  return;
}

