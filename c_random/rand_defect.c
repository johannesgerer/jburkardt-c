# include <stdlib.h>
# include <stdio.h>
# include <time.h>

int main ( void );
void timestamp ( void );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for RAND_DEFECT.

  Discussion:

    RAND_DEFECT illustrates a simple flaw in a scheme to compute random real
    values using TIME(0) to seed SRAND() before calling RAND().

    It is common to use the C random number generator RAND()
    to construct random real values, by dividing by RAND_MAX.
    It is also common to get a ``random'' seed by calling time(0),
    which returns the number of seconds since the Unix epoch.
    However this means that the value output by time(0) tends
    to be about the same number, and therefore the first value
    output by RAND() is about the same, and, more troubling,
    the first real number computed in this way will be about the same.

    Of course, two sequences that start close will quickly diverge.
    But the point is, I was trying to demonstrate a way to randomize
    a biased coin tossing game by checking whether a random real number
    was less than 0.6, and the game ended after one toss every time.

    Here is the first "random" real I got on successive runs:

      0.202452
      0.202507
      0.202593
      0.202655
      0.202710
      0.202781

    This problem could be avoided by having the user specify a seed,
    or by using a random number generator that works harder to shuffle
    the initial seed.  But the point is, the way we want RAND() to work
    and hence, the way we assume it works, is not the case, and this
    can easily cause serious problems.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    28 June 2011

  Author:

    John Burkardt
*/
{
  int i;
  double rd;
  int ri;
  unsigned int seed;

  timestamp ( );
  printf ( "\n" );
  printf ( "RAND_DEFECT:\n" );
  printf ( "  C version\n" );

  seed = time ( 0 );
  srand ( seed );

  printf ( "\n" );
  printf ( "  Initial seed = %d\n", seed );
  printf ( "\n" );
  printf ( "        RAND()   RAND()/RAND_MAX\n" );
  printf ( "\n" );

  for ( i = 0; i < 5; i++ )
  {
    ri = rand ( );
    rd = ( double ) ri / ( double ) RAND_MAX;
    printf ( "  %12d  %14g\n", ri, rd );
  }
/*
  Wait 2 seconds so we're guaranteed to get a different time 
  if we run the program again immediately.
*/
  sleep ( 2 );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "RAND_DEFECT:\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

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
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  printf ( "%s\n", time_buffer );

  return;
# undef TIME_SIZE
}
