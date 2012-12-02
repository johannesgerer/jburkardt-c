# include <stdlib.h>
# include <stdio.h>
# include <time.h>

int main ( int argc, char *argv[] );
void timestamp ( void );

/******************************************************************************/

int main ( int argc, char *argv[] )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for ARGS.

  Discussion:

    ARGS reports on the command line arguments.

    On input, ARGC is the number of command line arguments (including the
    program name itself), and ARGV is a pointer to an array of pointers, 
    the null-terminated strings that constitute the command line arguments.

    Thus, if we have invoked a program by:

      fred 1 alpha 3.7

    then on input

      ARGC is 4, and 

      ARGV --> ARGV(0) --> "fred"
               ARGV(1) --> "1"
               ARGV(2) --> "alpha"
               ARGV(3) --> "3.7"

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    26 April 2007

  Author:

    John Burkardt
*/
{
  int i;
  int VERBOSE = 1;

  if ( VERBOSE )
  {
    timestamp ( );

    printf ( "\n" );
    printf ( "ARGS\n" );
    printf ( "  C version\n" );
    printf ( "\n" );
    printf ( "  Compiled on %s at %s.\n", __DATE__, __TIME__ );
    printf ( "\n" );
    printf ( "  Print the command line arguments of a C program.\n" );
    printf ( "\n" );
    printf ( "  ARGC reports the number of arguments as %d\n", argc );
    printf ( "\n" );
  }

  for ( i = 0; i < argc; i++ ) 
  {
    printf ( "  %d  \"%s\"\n", i, *argv );
    argv++;
  }
/*
  Terminate.
*/
  if ( VERBOSE )
  {
    printf ( "\n" );
    printf ( "ARGS:\n" );
    printf ( "  Normal end of execution.\n" );
    printf ( "\n" );
    timestamp ( );
  }

  return 0;
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

  printf ( "%s\n", time_buffer );

  return;
# undef TIME_SIZE
}
