# include <stdlib.h>
# include <stdio.h>
# include <string.h>

# include "weekday.h"

int main ( );
void test01 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for WEEKDAY_PRB.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 May 2012

  Author:

    John Burkardt
*/
{
  timestamp ( );

  printf ( "\n" );
  printf ( "WEEKDAY_PRB:\n" );
  printf ( "  C version\n" );
  printf ( "  Test the WEEKDAY library.\n" );

  test01 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "WEEKDAY_PRB:\n" );
  printf ( "  Noraml end of execution.\n" );

  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void test01 ( )

/******************************************************************************/
/*
  Purpose:

    TEST01 tests YMD_TO_WEEKDAY_COMMON.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 May 2012

  Author:

    John Burkardt
*/
{
  int d;
  int m;
  int n_data;
  char *s1;
  char *s2;
  char *s3;
  int w1;
  int w2;
  int y;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  For dates in the Common calendar:\n" );
  printf ( "  YMD_TO_WEEKDAY_COMMON returns the day of the week.\n" );
  printf ( "\n" );
  printf ( "  YMD                   Weekday    Weekday\n" );
  printf ( "                        Tabulated  Computed\n" );
  printf ( "\n" );

  for ( ; ; )
  {
    weekday_values ( &n_data, &y, &m, &d, &w1 );

    if ( n_data == 0 )
    {
      break;
    }

    s3 = ymd_to_s_common ( y, m, d );
    w2 = ymd_to_weekday_common ( y, m, d );
    s1 = weekday_to_name_common ( w1 );
    s2 = weekday_to_name_common ( w2 );

    printf ( "  %20s  %9s  %9s\n", s3, s1, s2 );

    free ( s1 );
    free ( s2 );
    free ( s3 );
  }
  return;
}
