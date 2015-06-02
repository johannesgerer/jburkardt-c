# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>

int main ( );
void test01 ( );
double cpu_time ( );
void timestamp ( );

/*******************************************************************************/

int main ( )

/*******************************************************************************/
/*
  Purpose:

    MAIN is the main program for TIMER_CLOCK.

  Discussion:

    TIMER_CLOCK uses CLOCK as the timer.

    CLOCK is a timing utility accessible to C codes. It returns the number of
    "ticks" of processor time devoted to the user's job.  Dividing this by
    CLOCKS_PER_SEC (the number of clock ticks per second) gives a rough value for
    the elapsed CPU time.

    In this program, we call a routine CPU_TIME, which conveniently converts
    CLOCK's output to seconds for us.
  
  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    27 September 2005

  Author:

    John Burkardt
*/
{
  timestamp ( );

  printf ( "\n" );
  printf ( "TIMER_CLOCK\n" );
  printf ( "  C version\n" );
  printf ( "\n" );
  printf ( "  Demonstrate the use of the CLOCK timer.\n" );
  printf ( "\n" );
  printf ( "  CLOCK is a standard C library routine\n" );
  printf ( "  (defined in time.h)\n" );
  printf ( "\n" );
  printf ( "  It returns the processor time used by the program\n" );
  printf ( "  since the beginning of program execution.\n" );
  printf ( "  Divide this by CLOCKS_PER_SEC to convert to seconds.\n" );
  printf ( "\n" );
  printf ( "  CLOCK is a crude timer, and results less than\n" );
  printf ( "  a tenth of a second are probably not reliable.\n" );
  printf ( "\n" );
  printf ( "  The number of clock ticks per second is %d\n", CLOCKS_PER_SEC );

  test01 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "TIMER_CLOCK\n" );
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

    TEST01 times the C RAND routine.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    27 September 2005

  Author:

    John Burkardt
*/
{
# define N_MAX 1048756

  double cpu1;
  double cpu2;
  int i;
  int i_rep;
  int n;
  int n_log;
  int n_log_min = 10;
  int n_log_max = 20;
  int n_max = 1048756;
  int n_min = 1024;
  int n_rep = 5;
  int x[N_MAX];

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  Time the C RAND routine by computing N values.\n" );
  printf ( "  For a given N, repeat the computation 5 times.\n" );
  printf ( "\n" );
  printf ( "  Data vectors will be of minimum size %d\n", n_min );
  printf ( "  Data vectors will be of maximum size %d\n", n_max );
  printf ( "\n" );
  printf ( "  CPU times are in seconds.\n" );
  printf ( "\n" );
  printf ( "         N      Rep #1      Rep #2      Rep #2      Rep #4      Rep #5\n" );
  printf ( "\n" );

  for ( n_log = n_log_min; n_log <= n_log_max; n_log++ )
  {
    n = ( int ) pow ( 2, n_log );

    printf ( "  %8d", n );

    for ( i_rep = 0; i_rep < n_rep; i_rep++ )
    {
      cpu1 = cpu_time ( );

      for ( i = 0; i < n; i++ )
      {
        x[i] = rand ( );
      }

      cpu2 = cpu_time ( );

      printf ( "  %10f", cpu2 - cpu1 );
    }
    printf ( "\n" );
  }

  return;
}
/******************************************************************************/

double cpu_time ( )

/******************************************************************************/
/*
  Purpose:

    CPU_TIME returns the current reading on the CPU clock.

  Discussion:

    The CPU time measurements available through this routine are often
    not very accurate.  In some cases, the accuracy is no better than
    a hundredth of a second.  

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 June 2005

  Author:

    John Burkardt

  Parameters:

    Output, double CPU_TIME, the current reading of the CPU clock, in seconds.
*/
{
  double value;

  value = ( double ) clock ( ) 
        / ( double ) CLOCKS_PER_SEC;

  return value;
}
/*******************************************************************************/

void timestamp ( )

/*******************************************************************************/
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
