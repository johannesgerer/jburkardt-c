# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "chebyshev.h"

int main ( void );
void test01 ( void );
double f1 ( double x );
double f2 ( double x );
double f3 ( double x );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for CHEBYSHEV_PRB.

  Discussion:

    CHEBYSHEV_PRB tests the CHEBYSHEV library.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 February 2012

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "CHEBYSHEV_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the CHEBYSHEV library.\n" );
 
  test01 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "CHEBYSHEV_PRB\n" );
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

    TEST01 tests CHEBYSHEV_COEFFICIENTS and CHEBYSHEV_INTERPOLANT.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 February 2012

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double *c;
  double *fc;
  int i;
  int m;
  int n;
  double *x;
  
  printf ( "\n" );
  printf ( "CHEBYSHEV_TEST01\n" );
  printf ( "  CHEBYSHEV_COEFFICIENTS computes the coefficients of the\n" );
  printf ( "  Chebyshev interpolant.\n" );
  printf ( "  CHEBYSHEV_INTERPOLANT evaluates the interpolant.\n" );

  n = 5;
  a = -1.0;
  b = +1.0;

  c = chebyshev_coefficients ( a, b, n, f1 );

  x = chebyshev_zeros ( n );
  for ( i = 0; i < n; i++ )
  {
    x[i] = 0.5 * ( a + b ) + x[i] * 0.5 * ( b - a );
  }

  m = n;
  fc = chebyshev_interpolant ( a, b, n, c, m, x );

  printf ( "\n" );
  printf ( "  F(X) is a trig function:\n" );
  printf ( "\n" );
  printf ( "          X               C(I)            F(X)           C(F)(X)\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %14g  %14g  %14g  %14g\n", x[i], c[i], f1 ( x[i] ), fc[i] );
  }

  free ( c );
  free ( fc );
  free ( x );
/*
  Try a variant interval.
*/
  n = 5;
  a = 0.0;
  b = +3.0;

  c = chebyshev_coefficients ( a, b, n, f1 );

  x = chebyshev_zeros ( n );
  for ( i = 0; i < n; i++ )
  {
    x[i] = 0.5 * ( a + b ) + x[i] * 0.5 * ( b - a );
  }

  m = n;
  fc = chebyshev_interpolant ( a, b, n, c, m, x );

  printf ( "\n" );
  printf ( "  Consider the same F(X), but now over [0,3]:\n" );
  printf ( "\n" );
  printf ( "          X               C(I)            F(X)           C(F)(X)\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %14g  %14g  %14g  %14g\n", x[i], c[i], f1 ( x[i] ), fc[i] );
  }

  free ( c );
  free ( fc );
  free ( x );
/*
  Try a higher order.
*/
  n = 10;
  a = -1.0;
  b = +1.0;

  c = chebyshev_coefficients ( a, b, n, f1 );

  x = chebyshev_zeros ( n );
  for ( i = 0; i < n; i++ )
  {
    x[i] = 0.5 * ( a + b ) + x[i] * 0.5 * ( b - a );
  }

  m = n;
  fc = chebyshev_interpolant ( a, b, n, c, m, x );

  printf ( "\n" );
  printf ( "  Consider the same F(X), but now with higher order:\n" );
  printf ( "\n" );
  printf ( "          X               C(I)            F(X)           C(F)(X)\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %14g  %14g  %14g  %14g\n", x[i], c[i], f1 ( x[i] ), fc[i] );
  }

  free ( c );
  free ( fc );
  free ( x );
/*
  Try a polynomial.
*/
  n = 10;
  a = -1.0;
  b = +1.0;

  c = chebyshev_coefficients ( a, b, n, f3 );

  x = chebyshev_zeros ( n );
  for ( i = 0; i < n; i++ )
  {
    x[i] = 0.5 * ( a + b ) + x[i] * 0.5 * ( b - a );
  }

  m = n;
  fc = chebyshev_interpolant ( a, b, n, c, m, x );

  printf ( "\n" );
  printf ( "  F(X) is a degree 4 polynomial:\n" );
  printf ( "\n" );
  printf ( "          X               C(I)            F(X)           C(F)(X)\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %14g  %14g  %14g  %14g\n", x[i], c[i], f3 ( x[i] ), fc[i] );
  }

  free ( c );
  free ( fc );
  free ( x );
/*
  Try a function with decaying behavior.
*/
  n = 10;
  a = -1.0;
  b = +1.0;

  c = chebyshev_coefficients ( a, b, n, f2 );

  x = chebyshev_zeros ( n );
  for ( i = 0; i < n; i++ )
  {
    x[i] = 0.5 * ( a + b ) + x[i] * 0.5 * ( b - a );
  }

  m = n;
  fc = chebyshev_interpolant ( a, b, n, c, m, x );

  printf ( "\n" );
  printf ( "  The polynomial approximation to F(X) decays:\n" );
  printf ( "\n" );
  printf ( "          X               C(I)            F(X)           C(F)(X)\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %14g  %14g  %14g  %14g\n", x[i], c[i], f2 ( x[i] ), fc[i] );
  }

  free ( c );
  free ( fc );
  free ( x );

  return;
}
/******************************************************************************/

double f1 ( double x )

/******************************************************************************/
/*
  Purpose:

    F1 evaluates a function that can be used for Chebyshev interpolation.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 February 2012

  Author:

    John Burkardt

  Parameters:

    Input, double X, a point where the function is to be evaluated.

    Output, double F1, the function value.
*/
{
  double pi = 3.141592653589793;
  double value;

  value = cos ( 2.0 * pi * x ) * sin ( 3.0 * pi * x );

  return value;
}
/******************************************************************************/

double f2 ( double x )

/******************************************************************************/
/*
  Purpose:

    F2 evaluates a function that can be used for Chebyshev interpolation.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 February 2012

  Author:

    John Burkardt

  Parameters:

    Input,  double X, a point where the function is to be evaluated.

    Output, double F2, the function value.
*/
{
  double value;

  value = exp ( x );

  return value;
}
/******************************************************************************/

double f3 ( double x )

/******************************************************************************/
/*
  Purpose:

    F3 evaluates a function that can be used for Chebyshev interpolation.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 February 2012

  Author:

    John Burkardt

  Parameters:

    Input, double X, a point where the function is to be evaluated.

    Output, double F3, the function values.
*/
{
  double value;

  value = ( x - 3.0 ) * ( x - 1.0 ) * ( x - 1.0 ) * ( x + 2.0 );

  return value;
}
