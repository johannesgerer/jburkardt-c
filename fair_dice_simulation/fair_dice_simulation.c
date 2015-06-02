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

    25 May 2013

  Author:

    John Burkardt

  Parameters:

    Command line, int N, the number of times the dice are thrown.
*/
{
  char *command_filename = "fair_dice_commands.txt";
  FILE *command;
  char *data_filename = "fair_dice_data.txt";
  FILE *data;
  int die1;
  int die2;
  int i;
  int n;
  int seed;
  int score;
  int score_count[13];

  timestamp ( );
  printf ( "\n" );
  printf ( "FAIR_DICE_SIMULATION:\n" );
  printf ( "  C version\n" );
  printf ( "  Simulate N throws of a pair of fair dice.\n" );

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
  Create the graphics data file.
*/
  data = fopen ( data_filename, "wt" );
  for ( score = 2; score <= 12; score++ )
  {
    fprintf ( data, "  %d  %d\n", score, score_count[score] );
  }
  fclose ( data );
  printf ( "\n" );
  printf ( "  Created the graphics data file \"%s\".\n", data_filename );
/*
  Create the graphics command file.
*/
  command = fopen ( command_filename, "wt" );
  fprintf ( command, "# %s\n", command_filename );
  fprintf ( command, "#\n" );
  fprintf ( command, "# Usage:\n" );
  fprintf ( command, "#  gnuplot < %s\n", command_filename );
  fprintf ( command, "#\n" );
  fprintf ( command, "set term png\n" );
  fprintf ( command, "set output 'fair_dice.png'\n" );
  fprintf ( command, "set xlabel 'Score'\n" );
  fprintf ( command, "set ylabel 'Frequency'\n" );
  fprintf ( command, "set title 'Score frequency for a pair of fair dice'\n" );
  fprintf ( command, "set grid\n" );
  fprintf ( command, "set style fill solid\n" );
  fprintf ( command, "set yrange [0:*]\n" );
  fprintf ( command, "set timestamp\n" );
  fprintf ( command, "plot 'fair_dice_data.txt' using 1:2:(0.90):xtic(3) with boxes\n" );
  fprintf ( command, "quit\n" );

  fclose ( command );

  printf ( "  Created the graphics command file \"%s\".\n", command_filename );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "FAIR_DICE_SIMULATION:\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

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
