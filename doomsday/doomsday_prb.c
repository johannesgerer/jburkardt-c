# include <stdlib.h>
# include <stdbool.h>
# include <stdio.h>
# include <string.h>

# include "doomsday.h"

int main ( );
void test01 ( );
void test02 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for DOOMSDAY_PRB.
 
  Discussion:

    DOOMSDAY_PRB tests the DOOMSDAY library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    28 May 2012

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "DOOMSDAY_PRB:\n" );
  printf ( "  C version\n" );
  printf ( "  Test the DOOMSDAY library.\n" );

  test01 ( );
  test02 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "DOOMSDAY_PRB:\n" );
  printf ( "  Test the DOOMSDAY library.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void test01 ( )

/******************************************************************************/
/*
  Purpose:

    TEST01 tests DOOMSDAY against a couple of test dates.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    28 May 2012

  Author:

    John Burkardt
*/
{
  int d;
  int m;
  int n_data;
  int w;
  char *s1;
  char s2[10];
  int y;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  Try a couple selected dates.\n" );
  printf ( "\n" );
  printf ( "  YYYY  MM  DD  Weekday    Weekday\n" );
  printf ( "                Tabulated  Computed\n" );
  printf ( "\n" );

  y = 1989;
  m = 7;
  d = 13;
  w = doomsday_gregorian ( y, m, d );
  s1 = weekday_to_name_common ( w );
  strcpy ( s2, "Thursday" );
  printf ( "  %4d  %2d  %2d  %10s  %10s\n", y, m, d, s1, s2 );

  y = 2012;
  m = 5;
  d = 26;
  w = doomsday_gregorian ( y, m, d );
  s1 = weekday_to_name_common ( w );
  strcpy ( s2, "Saturday" );
  printf ( "  %4d  %2d  %2d  %10s  %10s\n", y, m, d, s1, s2 );

  return;
}
/******************************************************************************/

void test02 ( )

/******************************************************************************/
/*
  Purpose:

    TEST02 tests DOOMSDAY against a number of known values.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    28 May 2012

  Author:

    John Burkardt
*/
{
  int d;
  int m;
  int n_data;
  int w1;
  int w2;
  char *s1;
  char *s2;
  int y;

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  WEEKDAY_VALUES supplies a list of dates and weekdays.\n" );
  printf ( "\n" );
  printf ( "  YYYY  MM  DD  Weekday    Weekday\n" );
  printf ( "                Tabulated  Computed\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    weekday_values ( &n_data, &y, &m, &d, &w1 );

    if ( n_data <= 0 )
    {
      break;
    }
/*
  The transition from Julian to Gregorian calendars occurred in 1582
  (for some people).  The data in "WEEKDAY_VALUES" before the transition
  is stored in Julian format, which DOOMSDAY_GREGORIAN can't handle.
  So let's just refuse to handle 1582 or earlier
*/
    if ( y <= 1582 )
    {
      continue;
    }

    w2 = doomsday_gregorian ( y, m, d );

    s1 = weekday_to_name_common ( w1 );
    s2 = weekday_to_name_common ( w2 );

    printf ( "  %4d  %2d  %2d  %10s  %10s\n", y, m, d, s1, s2 );
  }

  return;
}
