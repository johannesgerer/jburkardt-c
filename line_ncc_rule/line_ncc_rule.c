# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <string.h>

# include "line_ncc_rule.h"

/******************************************************************************/

void line_ncc_rule ( int n, double a, double b, double x[], double w[] )

/******************************************************************************/
/*
  Purpose:

    LINE_NCC_RULE computes a Newton-Cotes Closed (NCC) quadrature rule.

  Discussion:

    The integral:

      Integral ( A <= X <= B ) F(X) dx

    The quadrature rule:

      Sum ( 1 <= I <= N ) W(I) * F ( X(I) ).

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    10 April 2014

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order.

    Input, double A, B, the endpoints of the interval.

    Input, double X[N], the abscissas.

    Output, double W[N], the weights.
*/
{
  double *d;
  int i;
  int j;
  int k;
  double y_a;
  double y_b;
/*
  Define the points X.
*/
  r8vec_linspace ( n, a, b, x );

  d = ( double * ) malloc ( n * sizeof ( double ) );

  for ( i = 0; i < n; i++ )
  {
/*
  Compute the Lagrange basis polynomial which is 1 at XTAB(I),
  and zero at the other nodes.
*/
    for ( j = 0; j < n; j++ )
    {
      d[j] = 0.0;
    }
    d[i] = 1.0;

    for ( j = 2; j <= n; j++ )
    {
      for ( k = j; k <= n; k++ )
      {
        d[n+j-k-1] = ( d[n+j-k-2] - d[n+j-k-1] ) / ( x[n-k] - x[n+j-k-1] );
      }
    }
    for ( j = 1; j <= n - 1; j++ )
    {
      for ( k = 1; k <= n - j; k++ )
      {
        d[n-k-1] = d[n-k-1] - x[n-k-j] * d[n-k];
      }
    }
/*
  Evaluate the antiderivative of the polynomial at the endpoints.
*/
    y_a = d[n-1] / ( double ) ( n );
    for ( j = n - 2; 0 <= j; j-- )
    {
      y_a = y_a * a + d[j] / ( double ) ( j + 1 );
    }
    y_a = y_a * a;

    y_b = d[n-1] / ( double ) ( n );
    for ( j = n - 2; 0 <= j; j-- )
    {
      y_b = y_b * b + d[j] / ( double ) ( j + 1 );
    }
    y_b = y_b * b;

    w[i] = y_b - y_a;
  }

  free ( d );

  return;
}
/******************************************************************************/

void r8vec_linspace ( int n, double a, double b, double x[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_LINSPACE creates a vector of linearly spaced values.

  Discussion:

    An R8VEC is a vector of R8's.

    4 points evenly spaced between 0 and 12 will yield 0, 4, 8, 12.
 
    In other words, the interval is divided into N-1 even subintervals,
    and the endpoints of intervals are used as the points.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    29 March 2011

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the vector.

    Input, double A, B, the first and last entries.

    Output, double X[N], a vector of linearly spaced data.
*/
{
  int i;

  if ( n == 1 )
  {
    x[0] = ( a + b ) / 2.0;
  }
  else
  {
    for ( i = 0; i < n; i++ )
    {
      x[i] = ( ( double ) ( n - 1 - i ) * a 
             + ( double ) (         i ) * b ) 
             / ( double ) ( n - 1     );
    }
  }
  return;
}
/******************************************************************************/

void timestamp ( )

/******************************************************************************/
/*
  Purpose:

    TIMESTAMP prints the current YMDHMS date as a time stamp.

  Example:

    31 May 2001 09:45:54 AM

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    24 September 2003

  Author:

    John Burkardt

  Parameters:

    None
*/
{
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct tm *tm;
  size_t len;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  len = strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  fprintf ( stdout, "%s\n", time_buffer );

  return;
# undef TIME_SIZE
}
