# include <stdlib.h>
# include <stdio.h>
# include <time.h>
# include <math.h>

# include "task_division.h"

/******************************************************************************/

int i4_div_rounded ( int a, int b )

/******************************************************************************/
/*
  Purpose:

    I4_DIV_ROUNDED computes the rounded result of I4 division.

  Discussion:

    This routine computes C = A / B, where A, B and C are integers
    and C is the closest integer value to the exact real result.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    23 October 2011

  Author:

    John Burkardt

  Parameters:

    Input, int A, B, the number to be divided,
    and the divisor.

    Output, int I4_DIV_ROUNDED, the rounded result
    of the division.
*/
{
  int a_abs;
  int b_abs;
  static int i4_huge = 2147483647;
  int value;

  if ( a == 0 && b == 0 )
  {
    value = i4_huge;
  }
  else if ( a == 0 )
  {
    value = 0;
  }
  else if ( b == 0 )
  {
    if ( a < 0 )
    {
      value = - i4_huge;
    }
    else
    {
      value = + i4_huge;
    }
  }
  else
  {
    a_abs = abs ( a );
    b_abs = abs ( b );

    value = a_abs / b_abs;
/*
  Round the value.
*/
    if ( ( 2 * value + 1 ) * b_abs < 2 * a_abs )
    {
      value = value + 1;
    }
/*
  Set the sign.
*/
    if ( ( a < 0 && 0 < b ) || ( 0 < a && b < 0 ) )
    {
      value = - value;
    }
  }
  return value;
}
/******************************************************************************/

void task_division ( int task_number, int proc_first, int proc_last )

/******************************************************************************/
/*
  Purpose:

    TASK_DIVISION divides tasks among processors.

  Discussion:

    This routine assigns each of T tasks to P processors, assuming that 
    the assignment is to be beforehand.

    In that case, we just want to make sure that we assign each task
    to a processor, that we assign about the same number of tasks
    to each processor, and that we assign each processor a contiguous
    range of tasks, say tasks I_LO to I_HI.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    23 October 2011

  Author:

    John Burkardt

  Parameters:

    Input, int TASK_NUMBER, the number of tasks.

    Input, int PROC_FIRST, PROC_LAST, the first and last processors.
*/
{
  int i_hi;
  int i_lo;
  int proc;
  int proc_number;
  int proc_remain;
  int task_proc;
  int task_remain;

  proc_number = proc_last + 1 - proc_first;

  printf ( "\n" );
  printf ( "TASK_DIVISION\n" );
  printf ( "  Divide T tasks among P processors.\n" );
  printf ( "\n" );
  printf ( "  Number of tasks T = %d\n", task_number );
  printf ( "  Number of processors P = %d\n", proc_number );
  printf ( "\n"  );
  printf ( "  P_FIRST = %d\n", proc_first );
  printf ( "  P_LAST  = %d\n", proc_last );


  printf ( "\n" );
  printf ( "             Number of   First      Last\n" );
  printf ( " Processor     Tasks     Task       Task\n" );
  printf ( "\n" );

  i_hi = 0;

  task_remain = task_number;
  proc_remain = proc_number;

  for ( proc = proc_first; proc <= proc_last; proc++ )
  {
    task_proc = i4_div_rounded ( task_remain, proc_remain );

    proc_remain = proc_remain - 1;
    task_remain = task_remain - task_proc;

    i_lo = i_hi + 1;
    i_hi = i_hi + task_proc;

    printf ( "  %8d  %8d  %8d  %8d\n", proc, task_proc, i_lo, i_hi );
  }

  return;
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
