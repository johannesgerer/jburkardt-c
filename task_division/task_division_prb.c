# include <stdlib.h>
# include <stdio.h>
# include <time.h>
# include <math.h>

# include "task_division.h"

int main ( int argc, char *argv[] );

/******************************************************************************/

int main ( int argc, char *argv[] )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for TASK_DIVISION_PRB.

  Discussion:

    TASK_DIVISION_PRB tests the TASK_DIVISION library.

    This program simply demonstrates how one might automate the
    assignment of T tasks to P processors, assuming that the assignment
    is to be beforehand.

    In that case, we just want to make sure that we assign each task
    to a processor, that we assign about the same number of tasks
    to each processor, and that we assign each processor a contiguous
    range of tasks, say tasks I_LO to I_HI.

    The routine that is called simulates this process.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    23 October 2011

  Author:

    John Burkardt
*/
{
  int proc_first;
  int proc_last;
  int task_number;

  timestamp ( );

  printf ( "\n" );
  printf ( "TASK_DIVISION_PRB:\n" );
  printf ( "  C version\n" );
  printf ( "  Test the TASK_DIVISION library.\n" );
  printf ( "  Demonstrate how to automate the division of\n" );
  printf ( "  T tasks among a range of P processors\n" );
  printf ( "  indexed from PROC_FIRST to PROC_LAST.\n" );

  task_number = 23;
  proc_first = 0;
  proc_last = 3;
  task_division ( task_number, proc_first, proc_last );

  task_number = 17;
  proc_first = 1;
  proc_last = 6;
  task_division ( task_number, proc_first, proc_last );

  task_number = 17;
  proc_first = 4;
  proc_last = 6;
  task_division ( task_number, proc_first, proc_last );

  task_number = 5;
  proc_first = -2;
  proc_last = 6;
  task_division ( task_number, proc_first, proc_last );

  task_number = 5;
  proc_first = 0;
  proc_last = 4;
  task_division ( task_number, proc_first, proc_last );

  task_number = 5;
  proc_first = 0;
  proc_last = 0;
  task_division ( task_number, proc_first, proc_last );

  task_number = 1000;
  proc_first = 1;
  proc_last = 17;
  task_division ( task_number, proc_first, proc_last );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "TASK_DIVISION_PRB:\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
