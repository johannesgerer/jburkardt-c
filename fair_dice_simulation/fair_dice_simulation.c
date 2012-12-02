# include <stdlib.h>
# include <stdio.h>
# include <time.h>

int main ( int argc, char *argv[] );
int random_int ( int a, int b );
void timestamp ( );

/******************************************************************************/

int main ( int argc, char *argv[] )

/******************************************************************************/
/*
  Purpose:

    FAIR_DICE_SIMULATION simulates N throws of two fair dice.

  Usage:

    fair_dice n

    where 

    * n is the number of times the dice should be thrown.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    04 September 2012

  Author:

    John Burkardt

  Parameters:

    Command line, int N, the number of times the dice are thrown.
*/
{
  int die1;
  int die2;
  int i;
  int n;
  int seed;
  int score;
  int score_count[13];

  if ( 0 )
  {
    timestamp ( );
    printf ( "\n" );
    printf ( "FAIR_DICE_SIMULATION:\n" );
    printf ( "  C version\n" );
    printf ( "  Simulate N throws of a pair of fair dice.\n" );
  }

  if ( 1 < argc )
  {
    n = atoi ( argv[1] );
  }
  else
  {
    printf ( "\n" );
    printf ( "  Enter N, the number of times the dice are thrown: " );
    scanf ( "%d", &n );
  }
/*
  Initialize the random number generator.
*/
  seed = time ( 0 );
  srand ( seed );
/*
  For convenience, include slots for 0 and 1, even though they never occur.
*/
  for ( i = 0; i <= 12; i++ )
  {
    score_count[i] = 0;
  }
/*
  Roll N times.
*/
  for ( i = 1; i <= n; i++ )
  {
    die1 = random_int ( 1, 6 );
    die2 = random_int ( 1, 6 );
    score = die1 + die2;
    score_count[score] = score_count[score] + 1;
  }
/*
  Print a table, suitable for treatement by GNUPLOT.
*/
  for ( score = 2; score <= 12; score++ )
  {
    printf ( "  %d  %d\n", score, score_count[score] );
  }
/*
  Terminate.
*/
  if ( 0 )
  {
    printf ( "\n" );
    printf ( "FAIR_DICE_SIMULATION:\n" );
    printf ( "  Normal end of execution.\n" );
    printf ( "\n" );
    timestamp ( );
  }

  return 0;
}
/******************************************************************************/

int random_int ( int a, int b )

/******************************************************************************/
/*
  Purpose:

   RANDOM_INT returns a random integer between A and B.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    04 September 2012

  Author:

    John Burkardt

  Parameters:

    Input, int A, B, the range of integers.

    Output, int RANDOM_INT, the random integer.
*/
{
  int range;
  int value;
/*
  If we want integers between A and B, there are actually
  B - A + 1 values.
*/
  range = b - a + 1;

  value = a + rand ( ) % range;

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
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  printf ( "%s\n", time_buffer );

  return;
# undef TIME_SIZE
}
