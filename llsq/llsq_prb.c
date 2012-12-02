# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "llsq.h"

int main ( );
void test01 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    LLSQ_PRB tests LLSQ.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    07 March 2012

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "LLSQ_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the LLSQ library.\n" );

  test01 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "LLSQ_PRB\n" );
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

    TEST01 calls LLSQ to match 15 data values.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    07 March 2012

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double error;
  int i;
  int n = 15;
  double x[15] = { 
    1.47, 1.50, 1.52, 1.55, 1.57, 1.60, 1.63, 1.65, 1.68, 1.70, 
    1.73, 1.75, 1.78, 1.80, 1.83 };
  double y[15] = {
    52.21, 53.12, 54.48, 55.84, 57.20, 58.57, 59.93, 61.29, 63.11, 64.47,
    66.28, 68.10, 69.92, 72.19, 74.46 };

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  LLSQ can compute the formula for a line of the form\n" );
  printf ( "    y = A * x + B\n" );
  printf ( "  which minimizes the RMS error to a set of N data values.\n" );

  llsq ( n, x, y, &a, &b );

  printf ( "\n" );
  printf ( "  Estimated relationship is y = %g * x + %g\n", a, b );
  printf ( "  Expected value is         y = 61.272 * x - 39.062\n" );
  printf ( "\n" );
  printf ( "     I      X       Y      B+A*X    |error|\n" );
  printf ( "\n" );
  error = 0.0;
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %4d  %7g  %7g  %7g  %7g\n", i, x[i], y[i], b + a * x[i], b + a * x[i] - y[i] );
    error = error + pow ( b + a * x[i] - y[i], 2 );
  }
  error = sqrt ( error / ( double ) n );
  printf ( "\n" );
  printf ( "  RMS error =                      %g\n", error );

  return;
}
