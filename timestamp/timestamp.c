# include <stdlib.h>
# include <stdio.h>
# include <time.h>

# include "timestamp.h"

/******************************************************************************/

double cpu_time ( )

/******************************************************************************/
/*
  Purpose:

    CPU_TIME returns the current reading on the CPU clock.

  Discussion:

    The CPU time measurements available through this routine are often
    not very accurate.  In some cases, the accuracy is no better than
    a hundredth of a second.  

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 June 2005

  Author:

    John Burkardt

  Parameters:

    Output, double CPU_TIME, the current reading of the CPU clock, in seconds.
*/
{
  double value;

  value = ( double ) clock ( ) 
        / ( double ) CLOCKS_PER_SEC;

  return value;
}
/******************************************************************************/

int *time_numbers ( )

/******************************************************************************/
/*
  Purpose:

    TIME_NUMBERS returns the data as a string of integers.

  Example:

    2001  Year
    5     Month
    31    Day
    9     Hour (0-23)
    45    Minute
    12    Second

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 June 2011

  Author:

    John Burkardt

  Parameters:

    Output, int TIME_NUMBERS[6], the year, month, day, hour, minute and second.
*/
{
  const struct tm *tm_ptr;
  time_t now;
  int *value;

  now = time ( 0 );
  tm_ptr = localtime ( &now );

  value = ( int * ) malloc ( 6 * sizeof ( int ) );

  value[0] = 1900 + tm_ptr->tm_year;
  value[1] = 1 + tm_ptr->tm_mon;
  value[2] = tm_ptr->tm_mday;
  value[3] = tm_ptr->tm_hour;
  value[4] = tm_ptr->tm_min;
  value[5] = tm_ptr->tm_sec;

  return value;
}
/******************************************************************************/

void timestamp ( )

/******************************************************************************/
/*
  Purpose:

    TIMESTAMP prints the current YMDHMS date as a time stamp.

  Example:

    17 June 2014 09:45:54 AM

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 June 2014

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
/******************************************************************************/

char *timestring ( )

/******************************************************************************/
/*
  Purpose:

    TIMESTRING returns the current YMDHMS date as a string.

  Example:

    31 May 2001 09:45:54 AM

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    15 January 2014

  Author:

    John Burkardt

  Parameters:

    Output, char *TIMESTRING, a string containing the current YMDHMS date.
*/
{
# define TIME_SIZE 40

  const struct tm *tm;
  size_t len;
  time_t now;
  char *s;

  now = time ( NULL );
  tm = localtime ( &now );

  s = ( char * ) malloc ( TIME_SIZE * sizeof ( char ) );

  len = strftime ( s, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  return s;
# undef TIME_SIZE
}
