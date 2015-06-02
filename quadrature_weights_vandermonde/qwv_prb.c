# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "qwv.h"

int main ( );
void test01 ( );
void test02 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for QWV_PRB.

  Discussion:

    QWV_PRB tests the QWV library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 February 2014

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "QWV_PRB:\n" );
  printf ( "  C version\n" );
  printf ( "  Test the QWV library.\n" );

  test01 ( );
  test02 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "QWV_PRB:\n" );
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

    TEST01 tests QWV for a Newton-Cotes rule.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 April 2011

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  int i;
  int n;
  double pi = 3.141592653589793;
  double theta;
  double *w;
  double *x;

  a =  0.0;
  b = +1.0;
  n = 5;

  printf ( "\n" );
  printf ( "TEST01:\n" );
  printf ( "  Use the Vandermonde procedure to compute the\n" );
  printf ( "  quadrature weights for a Newton-Cotes rule.\n" );
  printf ( "  Order N = %d\n", n );
  printf ( "  Interval = [%g,%g]\n", a, b );
/*
  Set the points.
*/
  x = r8vec_even_new ( n, a, b );
  r8vec_print ( n, x, "  Abscissas:" );
/*
  Compute the weights.
*/
  w = qwv ( n, a, b, x );

  r8vec_print ( n, w, "  Weights:" );
/*
  Free memory.
*/
  free ( w );
  free ( x );

  return;
}
/******************************************************************************/

void test02 ( )

/******************************************************************************/
/*
  Purpose:

    TEST02 tests QWV for a Clenshaw-Curtis rule.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 February 2014

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  int i;
  int n;
  double pi = 3.141592653589793;
  double theta;
  double *w;
  double *x;

  a = -1.0;
  b = +1.0;
  n = 5;

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  Use the Vandermonde procedure to compute the\n" );
  printf ( "  quadrature weights for a Clenshaw-Curtis rule.\n" );
  printf ( "  Order N = %d\n", n );
  printf ( "  Interval is [%g,%g]\n", a, b );
/*
  Set the points.
*/
  x = ( double * ) malloc ( n * sizeof ( double ) );

  for ( i = 0; i < n; i++ )
  {
    theta = ( double ) ( n - i - 1 ) * pi 
          / ( double ) ( n - 1 );

    x[i] = ( ( 1.0 - cos ( theta ) ) * a   
           + ( 1.0 + cos ( theta ) ) * b ) 
           /   2.0;
  }

  r8vec_print ( n, x, "  Abscissas:" );
/*
  Determine the corresponding weights.
*/
  w = qwv ( n, a, b, x );

  r8vec_print ( n, w, "  Weights:" );
/*
  Free memory.
*/
  free ( w );
  free ( x );

  return;
}
