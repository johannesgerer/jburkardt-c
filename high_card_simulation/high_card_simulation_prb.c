# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>

# include "high_card_simulation.h"

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for HIGH_CARD_SIMULATION_PRB.

  Discussion:

    HIGH_CARD_SIMULATION_PRB tests the HIGH_CARD_SIMULATION library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    25 February 2014

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "HIGH_CARD_SIMULATION_PRB\n" );
  printf ( "  C version.\n" );
  printf ( "  Test the HIGH_CARD_SIMULATION library.\n" );

  test01 ( );
  test02 ( );
  test03 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "HIGH_CARD_SIMULATION_PRB\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void test01 ( )

/******************************************************************************/
/*
  Purpose:

    TEST01 varies the skip number.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    25 February 2014

  Author:

    John Burkardt
*/
{
  int deck_size;
  int i;
  double p;
  int seed;
  int skip_num;
  int trial_num;

  deck_size = 100;
  trial_num = 100;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  Estimate the chances of picking the high\n" );
  printf ( "  card by skipping a given number of initial cards,\n" );
  printf ( "  using a deck of %d cards\n", deck_size );
  printf ( "  and simulating %d trials.\n", trial_num );
  printf ( "\n" );
  printf ( "  Skip   Deck Size    Chance of Winning\n" );
  printf ( "\n" );

  seed = 123456789;

  for ( i = 0; i < 10; i++ )
  {
    skip_num = 1 + ( i * deck_size ) / 10;

    p = high_card_simulation ( deck_size, skip_num, trial_num, &seed );

    printf ( "  %3d  %3d  %14.6g\n", skip_num, deck_size, p );
  }
  return;
}
/******************************************************************************/

void test02 ( )

/******************************************************************************/
/*
  Purpose:

    TEST02 plots the results for a deck of 100 cards.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    25 February 2014

  Author:

    John Burkardt
*/
{
  char command_filename[] = "test02_commands.txt";
  FILE *command_unit;
  char data_filename[] = "test02_data.txt";
  FILE *data_unit;
  int deck_size = 100;
  int i;
  double *p;
  int seed;
  int skip_num;
  int trial_num;

  trial_num = 1000;
  seed = 123456789;

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  Compute the changes of picking the high card\n" );
  printf ( "  after skipping 0 through 99 cards,\n" );
  printf ( "  using a deck with %d cards\n", deck_size );
  printf ( "  and taking %d trials.\n", trial_num );

  p = ( double * ) malloc ( deck_size * sizeof ( double ) );

  for ( skip_num = 0; skip_num < deck_size; skip_num++ )
  {
    p[skip_num] = high_card_simulation ( deck_size, skip_num, 
      trial_num, &seed );
  }
/*
  Create the data file.
*/
  data_unit = fopen ( data_filename, "wt" );
  for ( i = 0; i < deck_size; i++ )
  {
    fprintf ( data_unit, "  %d  %g\n", i, p[i] );
  }
  fclose ( data_unit );
  printf ( "\n" );
  printf ( "  Created graphics data file '%s'\n", data_filename );
/*
  Write the command file.
*/
  command_unit = fopen ( command_filename, "wt" );

  fprintf ( command_unit, "# %s\n", command_filename );
  fprintf ( command_unit, "#\n" );
  fprintf ( command_unit, "# Usage:\n" );
  fprintf ( command_unit, "#  gnuplot < %d\n", command_filename );
  fprintf ( command_unit, "#\n" );
  fprintf ( command_unit, "set term png\n" );
  fprintf ( command_unit, "set output 'test02.png'\n" );
  fprintf ( command_unit, "set xlabel '<--- Skip --->'\n" );
  fprintf ( command_unit, "set ylabel '<--- P(Correct) --->'\n" );
  fprintf ( command_unit, 
    "set title 'Estimated Prob of Correct Guess after Skipping K Cards'\n" );
  fprintf ( command_unit, "set grid\n" );
  fprintf ( command_unit, "set style data lines\n" );
  fprintf ( command_unit, "plot '%s' using 1:2 lw 3 linecolor rgb 'red'\n",
    data_filename );

  fclose ( command_unit );
  printf ( "  Created graphics command file '%s'\n", command_filename );

  free ( p );

  return;
}
/******************************************************************************/

void test03 ( )

/******************************************************************************/
/*
  Purpose:

    TEST03 plots the results for a deck of 100 cards.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    25 February 2014

  Author:

    John Burkardt
*/
{
  char command_filename[] = "test03_commands.txt";
  FILE *command_unit;
  char data_filename[] = "test03_data.txt";
  FILE *data_unit;
  int deck_size = 100;
  int i;
  double *p;

  printf ( "\n" );
  printf ( "TEST03\n" );
  printf ( "  HIGH_CARD_PROBABILITY computes the exact probability of \n" );
  printf ( "  winning a high card game with a deck of %d cards\n", deck_size );
  printf ( "  assuming we skip the first K cards and select the next card\n" );
  printf ( "  that is bigger.\n" );

  p = high_card_probability ( deck_size );

  printf ( "\n" );
  printf ( "    K   Prob(K)\n" );
  printf ( "\n" );
  for ( i = 0; i < deck_size; i++ )
  {
    printf ( "  %3d  %8.4f\n", i, p[i] );
  }
/*
  Create data file.
*/
  data_unit = fopen ( data_filename, "wt" );
  for ( i = 0; i < deck_size; i++ )
  {
    fprintf ( data_unit, "  %d  %g\n", i, p[i] );
  }
  fclose ( data_unit );
  printf ( "\n" );
  printf ( "  Created graphics data file '%s'\n", data_filename );
/*
  Create the command file.
*/
  command_unit = fopen ( command_filename, "wt" );

  fprintf ( command_unit, "# %s\n", command_filename );
  fprintf ( command_unit, "#\n" );
  fprintf ( command_unit, "# Usage:\n" );
  fprintf ( command_unit, "#  gnuplot < %d\n", command_filename );
  fprintf ( command_unit, "#\n" );
  fprintf ( command_unit, "set term png\n" );
  fprintf ( command_unit, "set output 'test03.png'\n" );
  fprintf ( command_unit, "set xlabel '<--- Skip --->'\n" );
  fprintf ( command_unit, "set ylabel '<--- P(Correct) --->'\n" );
  fprintf ( command_unit, 
    "set title 'Probability of Correct Guess after Skipping K Cards'\n" );
  fprintf ( command_unit, "set grid\n" );
  fprintf ( command_unit, "set style data lines\n" );
  fprintf ( command_unit, "plot '%s' using 1:2 lw 3 linecolor rgb 'red'\n",
    data_filename );

  fclose ( command_unit );
  printf ( "  Created graphics command file '%s'\n", command_filename );

  free ( p );

  return;
}
