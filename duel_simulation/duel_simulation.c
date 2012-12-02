# include <stdlib.h>
# include <stdio.h>
# include <time.h>

int main ( void );
int duel_result ( double a_accuracy, double b_accuracy );
double random_double ( void );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    DUEL_SIMULATION simulates a duel.

  Discussion:

    Player 1 fires at player 2, and hits with a probability of P(1).
    If Player 2 misses, then Player 2 fires at Player 1, hitting with
    a probability of P(2).

    The duel continues with alternating shots until only one player 
    survives.

    The simulation is intended to estimate the probabilities that a
    player will survive, and the number of turns required.

    The exact probability that player 1 will survive is

      P(1) / ( P(1) + P(2) - P(1) * P(2) )

    Player 2's chance is

     P(2) * ( 1 - P(1) ) / ( P(1) + P(2) - P(1) * P(2) )

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    04 September 2012

  Author:

    John Burkardt

  Reference:

    Paul Nahin,
    Duelling Idiots and Other Probability Puzzlers,
    Princeton University Press, 2000,
    ISBN13: 978-0691009797,
    LC: QA273.N29.

    Martin Shubik,
    "Does the Fittest Necessarily Survive?",
    in Readings in Game Theory and Political Behavior,
    edited by Martin Shubik,
    Doubleday, 1954,
    LC: H61.S53.

  Parameters:

    Input, double A_ACCURACY, B_ACCURACY, the probabilities that players A and
    B will hit their opponent in a single shot.

    Input, int DUEL_NUM, the number of duels to run.

    Output, double A_PROB, B_PROB, the estimated probablities that players
    A and B will survive.

    Output, double TURN_AVERAGE, the average number of turns required to
    complete the duel.
*/
{
  double a_accuracy;
  double a_prob;
  int a_wins;
  double b_accuracy;
  double b_prob;
  int b_wins;
  int duel;
  int duel_num;
  int seed;
  int winner;

  printf ( "\n" );
  printf ( "DUEL_SIMULATION:\n" );
  printf ( "  C version\n" );

  printf ( "Enter number of duels to run: " );
  scanf ( "%d", &duel_num );  

  printf ( "Enter player A's accuracy between 0.0 and 1.0: " );
  scanf ( "%lf", &a_accuracy );

  printf ( "Enter player B's accuracy between 0.0 and 1.0: " );
  scanf ( "%lf", &b_accuracy );

  seed = time ( 0 );
  srand ( seed );

  a_wins = 0;
  b_wins = 0;

  for ( duel = 1; duel <= duel_num; duel++ )
  {
    winner = duel_result ( a_accuracy, b_accuracy );

    if ( winner == 1 )
    {
      a_wins = a_wins + 1;
    }
    else
    {
      b_wins = b_wins + 1;
    }
  }

  a_prob = ( double ) a_wins / ( double ) duel_num;
  b_prob = ( double ) b_wins / ( double ) duel_num;

  printf ( "\n" );
  printf ( "  Player A wins with probability %g\n", a_prob );
  printf ( "  Player B wins with probability %g\n", b_prob );

  printf ( "\n" );
  printf ( "DUEL_SIMULATION:\n" );
  printf ( "  Normal end of execution\n" );

  return 0;
}
/******************************************************************************/

int duel_result ( double a_accuracy, double b_accuracy )

/******************************************************************************/
/*
  Purpose:

    DUEL_RESULT returns the outcome of a single duel.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 September 2012

  Author:

    John Burkardt

  Reference:

    Martin Shubik,
    “Does the Fittest Necessarily Survive?”,
    in Readings in Game Theory and Political Behavior,
    edited by Martin Shubik,
    Doubleday, 1954,
    LC: H61.S53.

  Parameters:

    Input, double A_ACCURACY, B_ACCURACY, the probabilities that player A
    and B will hit their opponent in a single shot.

    Output, int DUEL_RESULT, the survivor of the duel.
*/
{
  double r;
  int winner;

  while ( 1 )
  {
    r = random_double ( );

    if ( r <= a_accuracy )
    {
      winner = 1;
      break;
    }

    r = random_double ( );

    if ( r <= b_accuracy )
    {
      winner = 2;
      break;
    }
  }
  return winner;
}
/******************************************************************************/

double random_double ( )

/******************************************************************************/
/*
  Purpose:

    RANDOM_DOUBLE returns a random number between 0 and 1.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 November 2009

  Author:

    John Burkardt

  Parameters:

    Output, double RANDOM_DOUBLE, the random value.
*/
{
  double r;
  
  r = ( double ) rand ( ) / ( double ) RAND_MAX;

  return r;
}
