# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <string.h>

# include "poisson_simulation.h"

int main ( );
void test01 ( );
void test02 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    POISSON_SIMULATION_TEST tests POISSON_SIMULATION.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    28 September 2012

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "POISSON_SIMULATION_TEST\n" );
  printf ( "  C version.\n" );
  printf ( "  Test the POISSON_SIMULATION library.\n" );

  test01 ( );
  test02 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "POISSON_SIMULATION_TEST\n" );
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

    TEST01 simulates waiting for a given number of events.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    28 September 2012

  Author:

    John Burkardt
*/
{
  int bin_num = 30;
  char command_filename[80];
  FILE *command;
  char data_filename[80];
  FILE *data;
  int event_num = 1000;
  int *f_bin;
  int i;
  int j;
  double lambda;
  int seed;
  double *t;
  double *w;
  double w_ave;
  double *w_bin;
  double w_max;
  double w_min;
  double width;

  printf ( "\n" );
  printf ( "TEST01:\n" );
  printf ( "  POISSON_FIXED_EVENTS simulates a Poisson process\n" );
  printf ( "  until a given number of events have occurred.\n" );
  printf ( "\n" );
  printf ( "  Simulate a Poisson process, for which, on average,\n" );
  printf ( "  LAMBDA events occur per unit time.\n" );
  printf ( "  Run until you have observed EVENT_NUM events.\n" );
 
  lambda = 0.5;
  seed = 123456789;

  printf ( "\n" );
  printf ( "  LAMBDA = %g\n", lambda );
  printf ( "  EVENT_NUM = %d\n", event_num );

  t = ( double * ) malloc ( ( event_num + 1 ) * sizeof ( double ) );
  w = ( double * ) malloc ( ( event_num + 1 ) * sizeof ( double ) );
  poisson_fixed_events ( lambda, event_num, &seed, t, w );

  w_min = r8vec_min ( event_num + 1, w );
  w_max = r8vec_max ( event_num + 1, w );
  w_ave = r8vec_mean ( event_num + 1, w );

  printf ( "\n" );
  printf ( "  Minimum wait = %g\n", w_min );
  printf ( "  Average wait = %g\n", w_ave );
  printf ( "  Maximum wait = %g\n", w_max );

  printf ( "\n" );
  printf ( " Count            Time            Wait\n" );
  printf ( "\n" );
  for ( i = 0; i <= 5; i++ )
  {
    printf ( "  %d  %g  %g\n", i, t[i], w[i] );
  }
  printf ( "  ....  ..............  ..............\n" );
  for ( i = event_num - 5; i <= event_num; i++ )
  {
    printf ( "  %d  %g  %g\n", i, t[i], w[i] );
  }
/*
  Create the data file.
*/
  strcpy ( data_filename, "poisson_timeline_data.txt" );

  data = fopen ( data_filename, "wt" );

  for ( i = 0; i <= event_num; i++ )
  {
    fprintf ( data, "  %g  %d\n", t[i], i );
  }
  fclose ( data );

  printf ( " \n" );
  printf ( "  Data stored in \"%s\".\n", data_filename );
/*
  Create the command file.
*/
  strcpy ( command_filename, "poisson_timeline_commands.txt" );

  command = fopen ( command_filename, "wt" );

  fprintf ( command, "# poisson_timeline_commands.txt\n" );
  fprintf ( command, "#\n" );
  fprintf ( command, "# Usage:\n" );
  fprintf ( command, "#  gnuplot < poisson_timeline_commands.txt\n" );
  fprintf ( command, "#\n" );
  fprintf ( command, "set term png\n" );
  fprintf ( command, "set output 'poisson_timeline.png'\n" );
  fprintf ( command, "set style data lines\n" );
  fprintf ( command, "set xlabel 'Time'\n" );
  fprintf ( command, "set ylabel 'Number of Events'\n" );
  fprintf ( command, "set title 'Observation of Fixed Number of Poisson Events'\n" );
  fprintf ( command, "set grid\n" );
  fprintf ( command, "plot 'poisson_timeline_data.txt' using 1:2 lw 2\n" );
  fprintf ( command, "quit\n" );

  fclose ( command );

  printf ( "  Plot commands stored in \"%s\".\n", command_filename );
/*
  Determine bin information.
*/
  w_min = r8vec_min ( event_num + 1, w );
  w_max = r8vec_max ( event_num + 1, w );

  w_bin = r8vec_midspace_new ( bin_num, w_min, w_max );
  f_bin = ( int * ) malloc ( bin_num * sizeof ( int ) );

  for ( i = 0; i < bin_num; i++ )
  {
    f_bin[i] = 0;
  }

  for ( i = 0; i <= event_num; i++ )
  {
    j = 1 + ( int ) ( ( double ) ( bin_num ) * ( w[i] - w_min ) / ( w_max - w_min ) );
    j = i4_min ( j, bin_num );
    f_bin[j] = f_bin[j] + 1;
  }
/*
  Create the data file.
*/
  strcpy ( data_filename, "poisson_times_data.txt" );

  data = fopen ( data_filename, "wt" );

  for ( i = 0; i < bin_num; i++ )
  {
    fprintf ( data, "  %g  %d\n", w_bin[i], f_bin[i] );
  }
  fclose ( data );

  printf ( " \n" );
  printf ( "  Data stored in \"%s\".\n", data_filename );
/*
  Create the command file.
*/
  strcpy ( command_filename, "poisson_times_commands.txt" );

  command = fopen ( command_filename, "wt" );

  fprintf ( command, "# poisson_times_commands.txt\n" );
  fprintf ( command, "#\n" );
  fprintf ( command, "# Usage:\n" );
  fprintf ( command, "#  gnuplot < poisson_times_commands.txt\n" );
  fprintf ( command, "#\n" );
  fprintf ( command, "set term png\n" );
  fprintf ( command, "set output 'poisson_times.png'\n" );
  fprintf ( command, "set xlabel 'Waiting Time'\n" );
  fprintf ( command, "set ylabel 'Frequency'\n" );
  fprintf ( command, "set title 'Waiting Times Observed Over Fixed Time'\n" );
  fprintf ( command, "set grid\n" );
  fprintf ( command, "set style fill solid\n" );
  width = 0.85 * ( w_max - w_min ) / ( double ) ( bin_num );
  fprintf ( command, "plot 'poisson_times_data.txt' using 1:2:(%g) with boxes\n", width );
  fprintf ( command, "quit\n" );

  fclose ( command );

  printf ( "  Plot commands stored in \"%s\".\n", command_filename );

  free ( f_bin );
  free ( t );
  free ( w );
  free ( w_bin );

  return;
}
/******************************************************************************/

void test02 ( )

/******************************************************************************/
/*
  Purpose:

    TEST02 simulates waiting for a given length of time.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    28 September 2012

  Author:

    John Burkardt
*/
{
  int bin_num = 30;
  char command_filename[80];
  FILE *command;
  char data_filename[80];
  FILE *data;
  int *f_bin;
  int i;
  double lambda;
  int *n;
  double *n_bin;
  double n_max;
  double n_mean;
  double n_min;
  double n_var;
  int seed;
  double t;
  int test;
  int test_num = 20000;
  double w;

  lambda = 0.5;
  t = 1000.0;
  seed = 123456789;

  printf ( "\n" );
  printf ( "TEST02:\n" );
  printf ( "  POISSON_FIXED_EVENTS simulates a Poisson process\n" );
  printf ( "  counting the number of events that occur during\n" );
  printf ( "  a given time.\n" );
  printf ( "\n" );
  printf ( "  Simulate a Poisson process, for which, on average,\n" );
  printf ( "  LAMBDA events occur per unit time.\n" );
  printf ( "  Run for a total of %g time units.\n", t );
  printf ( "  LAMBDA = %g\n", lambda );

  n = ( int * ) malloc ( test_num * sizeof ( int ) );

  for ( test = 0; test < test_num; test++ )
  {
    n[test] = poisson_fixed_time ( lambda, t, &seed );
  }

  n_mean = i4vec_mean ( test_num, n );
  n_var = i4vec_variance ( test_num, n );
  printf ( "\n" );
  printf ( "  Mean number of events = %g\n", n_mean );
  printf ( "  Variance = %g\n", n_var );
  printf ( "  STD = %g\n", sqrt ( n_var ) );

  n_min = ( double ) ( i4vec_min ( test_num, n ) );
  n_max = ( double ) ( i4vec_max ( test_num, n ) );

  n_bin = r8vec_midspace_new ( bin_num, n_min, n_max );

  f_bin = ( int * ) malloc ( bin_num * sizeof ( int ) );
  for ( i = 0; i < bin_num; i++ )
  {
    f_bin[i] = 0;
  }
  for ( test = 0; test < test_num; test++ )
  {
    i = 1 + ( int ) ( ( double ) ( bin_num * ( n[test] - n_min ) ) 
      / ( double ) ( n_max - n_min ) );
    i = i4_min ( i, bin_num );
    f_bin[i] = f_bin[i] + 1;
  }
/*
  Create the data file.
*/
  strcpy ( data_filename, "poisson_events_data.txt" );

  data = fopen ( data_filename, "wt" );

  for ( i = 0; i < bin_num; i++ )
  {
    fprintf ( data, "  %g  %d\n", n_bin[i], f_bin[i] );
  }
  fclose ( data );

  printf ( " \n" );
  printf ( "  Data stored in \"%s\".\n", data_filename );
/*
  Create the command file.
*/
  strcpy ( command_filename, "poisson_events_commands.txt" );

  command = fopen ( command_filename, "wt" );

  fprintf ( command, "# poisson_events_commands.txt\n" );
  fprintf ( command, "#\n" );
  fprintf ( command, "# Usage:\n" );
  fprintf ( command, "#  gnuplot < poisson_events_commands.txt\n" );
  fprintf ( command, "#\n" );
  fprintf ( command, "set term png\n" );
  fprintf ( command, "set output 'poisson_events.png'\n" );
  fprintf ( command, "set xlabel 'Number of Events'\n" );
  fprintf ( command, "set ylabel 'Frequency'\n" );
  fprintf ( command, "set title 'Number of Poisson Events Over Fixed Time'\n" );
  fprintf ( command, "set grid\n" );
  fprintf ( command, "set style fill solid\n" );
  w = 0.85 * ( n_max - n_min ) / ( double ) ( bin_num );
  fprintf ( command, "plot 'poisson_events_data.txt' using 1:2:(%g) with boxes\n", w );
  fprintf ( command, "quit\n" );

  fclose ( command );

  printf ( "  Plot commands stored in \"%s\".\n", command_filename );

  free ( f_bin );
  free ( n );
  free ( n_bin );

  return;
}
