# include <stdlib.h>
# include <stdio.h>
# include <time.h>

# include "collatz_recursive.h"

/******************************************************************************/

void collatz_path ( int n )

/******************************************************************************/
/*
  Purpose:

    COLLATZ_PATH prints the members of a Collatz sequence.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    09 March 2012

  Author:

    John Burkardt

  Parameters:

    Input, int N, the current path member.
*/
{
  printf ( "  %d\n", n );

  if ( n == 1 )
  {
  }
  else if ( n % 2 == 0 )
  {
    collatz_path ( n / 2 );
  }
  else
  {
    collatz_path ( 3 * n + 1 );
  }
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
