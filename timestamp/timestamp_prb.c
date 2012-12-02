# include <stdlib.h>
# include <stdio.h>

# include "timestamp.h"

int main ( void );
void test01 ( void );
void test02 ( void );
void test03 ( void );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    TIMESTAMP_PRB demonstrates the use of TIMESTAMP.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 June 2011

  Author:

    John Burkardt
*/
{
  timestamp ( );

  printf ( "\n" );
  printf ( "TIMESTAMP_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the TIMESTAMP library.\n" );

  test01 ( );
  test02 ( );
  test03 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "TIMESTAMP_PRB\n" );
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

    TEST01 demonstrates the use of TIMESTAMP.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 June 2011

  Author:

    John Burkardt
*/
{
  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  TIMESTAMP prints out the current wallclock time,\n" );
  printf ( "  including the year, month, day, hours, minutes,\n" );
  printf ( "  seconds, thousandths of a second, and AM/PM.\n" );
  printf ( "\n" );
  printf ( "  This can be useful in keeping track of the date\n" );
  printf ( "  of execution of a particular program\n" );
  printf ( "  or to give a rough idea of the length of time\n" );
  printf ( "  required to run a program.\n" );

  printf ( "\n" );
  timestamp ( );

  return;
}
/******************************************************************************/

void test02 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST02 demonstrates the use of TIMESTRING.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 June 2011

  Author:

    John Burkardt
*/
{
  char *s;

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  TIMESTRING returns the current wallclock time,\n" );
  printf ( "  including the year, month, day, hours, minutes,\n" );
  printf ( "  seconds, thousandths of a second, and AM/PM\n" );
  printf ( "  in a string, which the user may print or manipulate.\n" );

  s = timestring ( );

  printf ( "\n" );
  printf ( "  TIMESTRING returned the value \"%s\"\n", s );

  return;
}
/******************************************************************************/

void test03 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST03 demonstrates the use of TIME_NUMBERS

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    05 June 2011

  Author:

    John Burkardt
*/
{
  int *time_vec;

  printf ( "\n" );
  printf ( "TEST03\n" );
  printf ( "  TIME_NUMBERS returns the date as a string of integers.\n" );

  time_vec = time_numbers ( );

  printf ( "\n" );
  printf ( "  Year =        %d\n", time_vec[0] );
  printf ( "  Month =       %d\n", time_vec[1] );
  printf ( "  Day =         %d\n", time_vec[2] );
  printf ( "  Hour =        %d\n", time_vec[3] );
  printf ( "  Minute =      %d\n", time_vec[4] );
  printf ( "  Second =      %d\n", time_vec[5] );

  free ( time_vec );

  return;
}
