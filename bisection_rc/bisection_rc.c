# include <stdlib.h>
# include <stdio.h>
# include <time.h>

# include "bisection_rc.h"

/******************************************************************************/

double bisection_rc ( double *a, double *b, double fx, int *job )

/******************************************************************************/
/*
  Purpose:

    BISECTION_RC seeks a zero of f(x) in a change of sign interval.

  Discussion:

    The bisection method is used.

    This routine uses reverse communication, so that the function is always
    evaluated in the calling program.

    On the first call, the user sets JOB = 0, and the values of A and B.
    Thereafter, the user checks the returned value of JOB and follows 
    directions.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 January 2013

  Author:

    John Burkardt

  Parameters:

    Input/output, double *A, *B, the endpoints of the change of 
    sign interval.  These values will change from call to call as the
    interval size is reduced.

    Input, double FX, the function value at the point X returned on
    the previous call.  On the first call, FX is ignored.

    Input/output, int *JOB, a communication flag.
    The user sets JOB to 0 before the first call.  Thereafter, the program
    controls setting the value of JOB, whose output values mean:

    Output, double BISECTION_RC, a point X at which the function is to 
    be evaluated.
*/
{
  static double fa;
  static double fb;
  static int state;
  double x;
  static double x_local;

  if ( *job == 0 )
  {
    fa = 0.0;
    fb = 0.0;
    state = 1;
    x = *a;
    *job = 1;
  }
  else if ( state == 1 )
  {
    fa = fx;
    x = *b;
    state = 2;
  }
  else if ( state == 2 )
  {
    fb = fx;

    if ( r8_sign ( fa ) == r8_sign ( fb ) )
    {
      fprintf ( stderr, "\n" );
      fprintf ( stderr, "BISECTION_RC - Fatal error!\n" );
      fprintf ( stderr, "  F(A) and F(B) have the same sign.\n" );
      exit ( 1 );
    }

    x = ( *a + *b ) / 2.0;
    state = 3;
  }
  else
  {
    if ( r8_sign ( fx ) == r8_sign ( fa ) )
    {
      *a = x_local;
      fa = fx;
    }
    else
    {
      *b = x_local;
      fb = fx;
    }
    x = ( *a + *b ) / 2.0;
    state = 3;
  }

  x_local = x;

  return x;
}
/******************************************************************************/

double r8_sign ( double x )

/******************************************************************************/
/*
  Purpose:

    R8_SIGN returns the sign of an R8.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 May 2006

  Author:

    John Burkardt

  Parameters:

    Input, double X, the number whose sign is desired.

    Output, double R8_SIGN, the sign of X.
*/
{
  double value;

  if ( x < 0.0 )
  {
    value = - 1.0;
  }
  else
  {
    value = + 1.0;
  }
  return value;
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
  size_t len;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  len = strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  fprintf ( stdout, "%s\n", time_buffer );

  return;
# undef TIME_SIZE
}

