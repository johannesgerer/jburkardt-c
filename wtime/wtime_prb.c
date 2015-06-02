# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>

# include "wtime.h"

int main ( void );
void test01 ( void );
int i4_power ( int i, int j );
void timestamp ( void );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for WTIME_PRB.

  Discussion:

    WTIME_PRB tests the WTIME library.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    27 April 2009

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "WTIME_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the WTIME library.\n" );

  test01 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "WTIME_PRB\n" );
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

    TEST01 times the RAND routine.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    27 April 2009

  Author:

    John Burkardt
*/
{
  double cpu1;
  double cpu2;
  int i;
  int i_rep;
  int n;
  int n_log;
  int n_log_min = 10;
  int n_log_max = 20;
  int n_max;
  int n_min;
  int n_rep = 5;
  double seconds;
  int *x;

  n_min = i4_power ( 2, n_log_min );
  n_max = i4_power ( 2, n_log_max );

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  Time the RAND routine by computing N values.\n" );
  printf ( "  For a given N, repeat the computation 5 times.\n" );
  printf ( "\n" );
  printf ( "  Data vectors will be of minimum size %d\n", n_min );
  printf ( "  Data vectors will be of maximum size %d\n", n_max );
  printf ( "\n" );
  printf ( "  Times are measured in seconds.\n" );
  printf ( "\n" );
  printf ( "         N      Rep #1      Rep #2      Rep #2      Rep #4      Rep #5\n" );
  printf ( "\n" );

  for ( n_log = n_log_min; n_log <= n_log_max; n_log++ )
  {
    n = i4_power ( 2, n_log );

    x = ( int * ) malloc ( n * sizeof ( int ) );

    printf ( "  %8d", n );

    for ( i_rep = 0; i_rep < n_rep; i_rep++ )
    {
      seconds = wtime ( );

      for ( i = 0; i < n; i++ )
      {
        x[i] = rand ( );
      }

      seconds = wtime ( ) - seconds;

      printf ( "  %10f", seconds );
    }
    printf ( "\n" );
    free ( x );
  }

  return;
}
/******************************************************************************/

int i4_power ( int i, int j )

/******************************************************************************/
/*
  Purpose:

    I4_POWER returns the value of I^J.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    23 October 2007

  Author:

    John Burkardt

  Parameters:

    Input, int I, J, the base and the power.  J should be nonnegative.

    Output, int I4_POWER, the value of I^J.
*/
{
  int k;
  int value;

  if ( j < 0 )
  {
    if ( i == 1 )
    {
      value = 1;
    }
    else if ( i == 0 )
    {
      printf ( "\n" );
      printf ( "I4_POWER - Fatal error!\n" );
      printf ( "  I^J requested, with I = 0 and J negative.\n" );
      exit ( 1 );
    }
    else
    {
      value = 0;
    }
  }
  else if ( j == 0 )
  {
    if ( i == 0 )
    {
      printf ( "\n" );
      printf ( "I4_POWER - Fatal error!\n" );
      printf ( "  I^J requested, with I = 0 and J = 0.\n" );
      exit ( 1 );
    }
    else
    {
      value = 1;
    }
  }
  else if ( j == 1 )
  {
    value = i;
  }
  else
  {
    value = 1;
    for ( k = 1; k <= j; k++ )
    {
      value = value * i;
    }
  }
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
  size_t len;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  len = strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  printf ( "%s\n", time_buffer );

  return;
# undef TIME_SIZE
}
