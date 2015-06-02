# include <stdlib.h>
# include <stdio.h>
# include <string.h>
# include <time.h>

# include "snakes.h"

/******************************************************************************/

double *snakes ( )

/******************************************************************************/
/*
  Purpose:

    SNAKES sets up the Snakes and Ladders matrix.

  Discussion:

    Snakes and Ladders, also known as Chutes and Ladders, is a game
    played on a 10x10 board of 100 squares.  A player can be said to
    start at square 0, that is, off the board.  The player repeatedly
    rolls a die, and advances between 1 and 6 steps accordingly.
    The game is won when the player reaches square 100.  In some versions,
    the player must reach 100 by exact die count, forfeiting the move
    if 100 is exceeded; in others, reaching or exceeding 100 counts as
    a win.

    Play is complicated by the existence of snakes and ladders.  On
    landing on a square that begins a snake or ladder, the player is
    immediately tranported to another square, which will be lower for
    a snake, or higher for a ladder.

    Typically, several players play, alternating turns.

    Given a vector V(0:100) which is initially all zero except for the
    first entry, the matrix-vector product A'*V represents the probabilities
    that a player starting on square 0 will be on any given square after one
    roll.  Correspondingly, (A')^2*V considers two moves, and so on.  Thus,
    repeatedly multiplying by A' reveals the probability distribution 
    associated with the likelihood of occupying any particular square at a 
    given turn in the game.  

    There is a single eigenvalue of value 1, whose corresponding eigenvector
    is all zero except for a final entry of 1, representing a player who
    has reached square 100.  All other eigenvalues have norm less than 1,
    corresponding to the fact that there are no other long term steady
    states or cycles in the game.

    Note that no two descriptions of the Snakes and Ladders board seem to
    agree.  This is the board described by Nick Berry.  The board described 
    by Higham and Higham is close to this one, but differs in the description 
    of two of the jumps.

    While most commentators elect to move immediately from a snake mouth or
    ladder foot, I have decide there are reasons to treat the game in such a
    way that when you land on a ladder foot or snake mouth, you stay there
    as though you had landed on an ordinary square; the difference arises on
    your next turn, when, instead of rolling a die, you move up the ladder
    or down the snake.  This allows the player to "register" a stop at the
    given square, may be suitable for certain applications, and makes for
    a transition matrix whose structure is more obvious to understand.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 September 2014

  Author:

    John Burkardt

  Reference:

    Steve Althoen, Larry King, Kenneth Schilling,
    How long is a game of Snakes and Ladders?,
    The Mathematical Gazette,
    Volume 77, Number 478, March 1993, pages 71-76.

    Nick Berry,
    Mathematical Analysis of Chutes and Ladders,
    http://www.datagenetics.com/blog/november12011/index.html

    Desmond Higham, Nicholas Higham,
    MATLAB Guide,
    SIAM, 2005,
    ISBN13: 9780898717891.

  Parameters:

    Output, double SNAKES[101*101], the matrix.
*/
{
  double *a;
  int d;
  int i;
  int j;
  int j1;
  int j2;
  int *jump;
  int k;

  a = ( double * ) malloc ( 101 * 101 * sizeof ( double ) );

  jump = ( int * ) malloc ( 101 * sizeof ( int ) );

  for ( i = 0; i < 101; i++ )
  {
    jump[i] = i;
  }

  jump[ 1] =  38;
  jump[ 4] =  14;
  jump[ 9] =  31;
  jump[16] =   6;
  jump[21] =  42;
  jump[28] =  84;
  jump[36] =  44;
  jump[48] =  26;
  jump[49] =  11;
  jump[51] =  67;
  jump[56] =  53;
  jump[62] =  19;
  jump[64] =  60;
  jump[71] =  91;
  jump[80] = 100;
  jump[87] =  24;
  jump[93] =  73;
  jump[95] =  75;
  jump[98] =  78;

  for ( j = 0; j < 101; j++ )
  {
    for ( i = 0; i < 101; i++ )
    {
      a[i+j*101] = 0.0;
    }
  }
/*
  A(I,J) represents the probablity that a dice roll will take you from
  square I to square J.

  Starting in square I...
*/
  for ( i = 0; i < 101; i++ )
  {
/*
  If I is a snake or ladder, go to the next spot.
*/
    if ( i != jump[i] )
    {
      j = jump[i];
      a[i+j*101] = 1.0;
    }
/*
  Otherwise, roll a die
*/
    else
    {
      for ( d = 1; d <= 6; d++ )
      {
/*
  so theoretically, our new location J will be I + D,
*/
        j = i + d;
/*
  but if J is greater than 100, move us back to J,
*/
        if ( 100 < j )
        {
          j = 100;
        }
  
        a[i+j*101] = a[i+j*101] + 1.0 / 6.0;
      }
    }
  }

  free ( jump );

  return a;
}
/******************************************************************************/

void spy_ge ( int m, int n, double a[], char *header )

/******************************************************************************/
/*
  Purpose:

    SPY_GE plots a sparsity pattern for a general storage (GE) matrix.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 September 2014

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns
    in the matrix.

    Input, double A[M*N], the matrix.

    Input, char *HEADER, the name to be used for the
    title of the plot, and as part of the names of the data, command
    and plot files.
*/
{
  char command_filename[255];
  FILE *command_unit;
  char data_filename[255];
  FILE *data_unit;
  int i;
  int j;
  int nz_num;
  char png_filename[255];
/*
  Create data file.
*/
  strcpy ( data_filename, header );
  strcat ( data_filename, "_data.txt" );
  data_unit = fopen ( data_filename, "wt" );
  nz_num = 0;
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      if ( a[i+j*m] != 0.0 )
      {
        fprintf ( data_unit, "%d  %d\n", j, i );
        nz_num = nz_num + 1;
      }
    }
  }
  fclose ( data_unit );
  printf ( "\n" );
  printf ( "  Created sparsity data file '%s'\n", data_filename );
/*
  Create command file.
*/
  strcpy ( command_filename, header );
  strcat ( command_filename, "_commands.txt" );
  command_unit = fopen ( command_filename, "wt" );

  fprintf ( command_unit, "# %s\n", command_filename );
  fprintf ( command_unit, "#\n" );
  fprintf ( command_unit, "# Usage:\n" );
  fprintf ( command_unit, "#  gnuplot < %s\n", command_filename );
  fprintf ( command_unit, "#\n" );
  fprintf ( command_unit, "unset key\n" );
  fprintf ( command_unit, "set term png\n" );

  strcpy ( png_filename, header );
  strcat ( png_filename, ".png" );
  fprintf ( command_unit, "set output '%s'\n", png_filename );
  fprintf ( command_unit, "set size ratio -1\n" );
  fprintf ( command_unit, "set xlabel '<--- J --->'\n" );
  fprintf ( command_unit, "set ylabel '<--- I --->'\n" );
  fprintf ( command_unit, "set title '%d nonzeros for \"%s\"'\n", 
    nz_num, header );
  fprintf ( command_unit, "set timestamp\n" );
  fprintf ( command_unit, 
    "plot [x=0:%d] [y=%d:0] '%s' with points pt 5\n",
    n-1, m-1, data_filename );

  fclose ( command_unit );
  printf ( "  Created graphics command file '%s'\n", command_filename );

  return;
}
/******************************************************************************/

void timestamp ( )

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
