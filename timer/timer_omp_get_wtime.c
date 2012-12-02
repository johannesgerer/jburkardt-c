# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <omp.h>

int main ( void );
void test01 ( void );
void timestamp ( void );

/*******************************************************************************/

int main ( void )

/*******************************************************************************/
/*
  Purpose:

    MAIN is the main program for TIMER_OMP_GET_WTIME.

  Discussion:

    TIMER_OMP_GET_WTIME uses OMP_GET_WTIME as the timer.

    OMP_GET_WTIME is a timing utility accessible to C codes that
    support OpenMP.  It returns the elapsed wallclock time in seconds.

    Here, we run on as many threads as there are processors.  We could
    force the number of threads to be 1 to make a better comparison to
    timers that run on a single processor.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    10 July 2008

  Author:

    John Burkardt
*/
{
  int proc_num;
  int thread_num;

  timestamp ( );

  printf ( "\n" );
  printf ( "TIMER_OMP_GET_WTIME\n" );
  printf ( "  C version\n" );
  printf ( "\n" );
  printf ( "  Demonstrate the use of the OMP_GET_WTIME timer.\n" );
  printf ( "\n" );
  printf ( "  omp_get_wtime ( ) is an OpenMP library function.\n" );
  printf ( "\n" );
  printf ( "  It returns the elapsed wall clock time in seconds.\n" );
/*
  How many processors are available?
*/
  proc_num = omp_get_num_procs ( );

  printf ( "\n" );
  printf ( "  The number of processors available:\n" );
  printf ( "  OMP_GET_NUM_PROCS () = %d\n", proc_num );
  
  thread_num = proc_num;

  printf ( "\n" );
  printf ( "  OMP_SET_NUM_THREADS requests %d threads.\n", thread_num );

  omp_set_num_threads ( thread_num );

  test01 ( );

  printf ( "\n" );
  printf ( "TIMER_OMP_GET_WTIME\n" );
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

    TEST01 times the C RAND routine.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    10 July 2008

  Author:

    John Burkardt
*/
{
# define N_MAX 1048756

  int i;
  int i_rep;
  int n;
  int n_log;
  int n_log_min = 10;
  int n_log_max = 20;
  int n_max = 1048756;
  int n_min = 1024;
  int n_rep = 5;
  double wtime1;
  double wtime2;
  int x[N_MAX];

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  Time the C RAND routine by computing N values.\n" );
  printf ( "  For a given N, repeat the computation 5 times.\n" );
  printf ( "\n" );
  printf ( "  Data vectors will be of minimum size %d\n", n_min );
  printf ( "  Data vectors will be of maximum size %d\n", n_max );
  printf ( "\n" );
  printf ( "  Wall clock times are in seconds.\n" );
  printf ( "\n" );
  printf ( "         N      Rep #1      Rep #2      Rep #2      Rep #4      Rep #5\n" );
  printf ( "\n" );

  for ( n_log = n_log_min; n_log <= n_log_max; n_log++ )
  {
    n = ( int ) pow ( 2, n_log );

    printf ( "  %8d", n );

    for ( i_rep = 0; i_rep < n_rep; i_rep++ )
    {
      wtime1 = omp_get_wtime ( );

      for ( i = 0; i < n; i++ )
      {
        x[i] = rand ( );
      }

      wtime2 = omp_get_wtime ( );

      printf ( "  %10f", wtime2 - wtime1 );
    }
    printf ( "\n" );
  }

  return;
}
/*******************************************************************************/

void timestamp ( void )

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
