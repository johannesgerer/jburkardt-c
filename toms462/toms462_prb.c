# include <stdlib.h>
# include <stdio.h>

# include "toms462.h"

int main ( );
void test01 ( );
void test02 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    TOMS462_PRB tests BIVNOR.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    13 April 2012

  Author:

    John Burkardt
*/
{
  timestamp ( );

  printf ( "\n" );
  printf ( "TOMS462_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the TOMS462 library.\n" );

  test01 ( );
  test02 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "TOMS462_PRB\n" );
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

    TEST01 tests BIVNOR.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    13 April 2012

  Author:

    John Burkardt
*/
{
  double cdf;
  double expect;
  double r;
  double x;
  double y;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  Compare BIVNOR with some simple data\n" );
  printf ( "  with 3 digit accuracy.\n" );
  printf ( "\n" );
  printf ( "       X         Y          R          P               P\n" );
  printf ( "                                      (Tabulated)     (BIVNOR)\n" );
  printf ( "\n" );

  x =  0.8;
  y = -1.5;
  r =  -0.9;
  expect = 0.148;
  cdf = bivnor ( x, y, r );
  printf ( "  %8.3f  %8.3f  %8.3f  %14.6g  %14.6g\n", x, y, r, expect, cdf );

  x =  0.6;
  y = -1.4;
  r =  -0.7;
  expect = 0.208;
  cdf = bivnor ( x, y, r );
  printf ( "  %8.3f  %8.3f  %8.3f  %14.6g  %14.6g\n", x, y, r, expect, cdf );

  x =  0.2;
  y = -1.0;
  r =  -0.5;
  expect = 0.304;
  cdf = bivnor ( x, y, r );
  printf ( "  %8.3f  %8.3f  %8.3f  %14.6g  %14.6g\n", x, y, r, expect, cdf );

  x = -1.2;
  y =  0.1;
  r =   0.0;
  expect = 0.407;
  cdf = bivnor ( x, y, r );
  printf ( "  %8.3f  %8.3f  %8.3f  %14.6g  %14.6g\n", x, y, r, expect, cdf );

  x = -1.2;
  y = -0.1;
  r =   0.3;
  expect = 0.501;
  cdf = bivnor ( x, y, r );
  printf ( "  %8.3f  %8.3f  %8.3f  %14.6g  %14.6g\n", x, y, r, expect, cdf );

  x = -0.4;
  y = -0.9;
  r =   0.6;
  expect = 0.601;
  cdf = bivnor ( x, y, r );
  printf ( "  %8.3f  %8.3f  %8.3f  %14.6g  %14.6g\n", x, y, r, expect, cdf );

  return;
}
/******************************************************************************/

void test02 ( )

/******************************************************************************/
/*
  Purpose:

    TEST02 tests BIVNOR.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    13 April 2012

  Author:

    John Burkardt
*/
{
  double fxy1;
  double fxy2;
  int n_data;
  double r;
  double x;
  double y;

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  Compare BIVNOR with some tabulated data.\n" );
  printf ( "\n" );
  printf ( "      X          Y          " );
  printf ( "R           P                         P" );
  printf ( "                      DIFF\n" );
  printf ( "                                " );
  printf ( "       (Tabulated)               (BIVNOR)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    bivariate_normal_cdf_values ( &n_data, &x, &y, &r, &fxy1 );

    if ( n_data == 0 )
    {
      break;
    }
/*
  BIVNOR computes the "tail" of the probability, and we want the
  initial part.
*/
    fxy2 = bivnor ( - x, - y, r );

    printf ( "  %8.3f  %8.3f  %8.3f  %24.16g  %24.16g  %10.4g\n",
      x, y, r, fxy1, fxy2, r8_abs ( fxy1 - fxy2 ) );
  }
  return;
}
