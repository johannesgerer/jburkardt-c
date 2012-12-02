# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>

# include "llsq.h"

/******************************************************************************/

void llsq ( int n, double x[], double y[], double *a, double *b )

/******************************************************************************/
/*
  Purpose:

    LLSQ solves a linear least squares problem matching a line to data.

  Discussion:

    A formula for a line of the form Y = A * X + B is sought, which
    will minimize the root-mean-square error to N data points ( X[I], Y[I] );

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    07 March 2012

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of data values.

    Input, double X[N], Y[N], the coordinates of the data points.

    Output, double *A, *B, the slope and Y-intercept of the least-squares
    approximant to the data.
*/
{
  double bot;
  int i;
  double top;
  double xbar;
  double ybar;
/*
  Special case.
*/
  if ( n == 1 )
  {
    *a = 0.0;
    *b = y[0];
    return;
  }
/*
  Average X and Y.
*/
  xbar = 0.0;
  ybar = 0.0;
  for ( i = 0; i < n; i++ )
  {
    xbar = xbar + x[i];
    ybar = ybar + y[i];
  }
  xbar = xbar / ( double ) n;
  ybar = ybar / ( double ) n;
/*
  Compute Beta.
*/
  top = 0.0;
  bot = 0.0;
  for ( i = 0; i < n; i++ )
  {
    top = top + ( x[i] - xbar ) * ( y[i] - ybar );
    bot = bot + ( x[i] - xbar ) * ( x[i] - xbar );
  }
  *a = top / bot;

  *b = ybar - *a * xbar;

  return;
}
/******************************************************************************/

void timestamp ( void )

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
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  printf ( "%s\n", time_buffer );

  return;
# undef TIME_SIZE
}
