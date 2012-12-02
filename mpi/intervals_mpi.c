# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>

# include "mpi.h"

int main ( int argc, char *argv[] );
double f ( double x );
void timestamp ( void );

/******************************************************************************/

int main ( int argc, char *argv[] )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for INTERVALS.

  Discussion:

    INTERVALS uses MPI routines to multiprocess a computational task.

    We have a function F(X), an interval [XMIN,XMAX], 
    and a value N.

    We define N equally spaced points in the interval,

      X(I) = ( ( N - I     ) * XMIN 
             + (     I - 1 ) * XMAX ) 
             / ( N     - 1 )

    We thus have N-1 subintervals.

    We assume we have N processors available.

    Processor 0 is designated the master processor, assigned
    to estimating the integral of F(X) over the entire
    interval [ X(1), X(N) ].

    For I = 1 to N-1, processor I is assigned the subinterval

      [ X(I), X(I+1) ]

    and then estimates the integral Q(I) of F(X) over that
    subinterval.

    COMMUNICATION:

    Processor 0 communicates to processor I the endpoints of 
    the interval it is assigned, and the number of sample points
    to use in that interval.

    Processor I communicates to processor 0 the computed value of Q(I).

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 October 2007

  Author:

    John Burkardt

  Reference:

    William Gropp, Ewing Lusk, Anthony Skjellum,
    Using MPI: Portable Parallel Programming with the
    Message-Passing Interface,
    Second Edition,
    MIT Press, 1999,
    ISBN: 0262571323.

    Snir, Otto, Huss-Lederman, Walker, Dongarra,
    MPI - The Complete Reference,
    Volume 1, The MPI Core,
    second edition,
    MIT Press, 1998.
*/
{
  double end_time;
  double h;
  int i;
  int ierr;
  int m;
  int master = 0;
  int n;
  double pi = 3.141592653589793238462643;
  int process;
  int process_id;
  int process_num;
  double q_global;
  double q_local;
  int received;
  int source;
  double start_time;
  MPI_Status status;
  int tag;
  int target;
  double x;
  double xb[2];
  double x_max = 1.0;
  double x_min = 0.0;
/*
  Establish the MPI environment.
*/
  ierr = MPI_Init ( &argc, &argv );
/*
  Determine this processes's rank.
*/
  ierr = MPI_Comm_rank ( MPI_COMM_WORLD, &process_id );
/*
  Get the number of processes.
*/
  ierr = MPI_Comm_size ( MPI_COMM_WORLD, &process_num );
/*
  Say hello (once), and shut down right away unless we
  have at least 2 processes available.
*/
  if ( process_id == master )
  {
    timestamp ( );
    printf ( "\n" );
    printf ( "INTERVALS - Master process:\n" );
    printf ( "  C version\n" );
    printf ( "\n" );
    printf ( "  An MPI example program,\n" );
    printf ( "  A quadrature over an interval is done by\n" );
    printf ( "  assigning subintervals to processes.\n" );
    printf ( "\n" );
    printf ( "  The number of processes is %d\n", process_num );

    start_time = MPI_Wtime ( );

    if ( process_num <= 1 )
    {
      printf ( "\n" );
      printf ( "INTERVALS - Master process:\n" );
      printf ( "  Need at least 2 processes!\n" );
      ierr = MPI_Finalize ( );
      printf ( "\n" );
      printf ( "INTERVALS - Master process:\n" );
      printf ( "  Abnormal end of execution.\n" );
      exit ( 1 );
    }
  }

  printf ( "\n" );
  printf ( "Process %d is active!\n", process_id );
/*
  Every process could figure out the endpoints of its interval
  on its own.  But we want to demonstrate communication.  So we
  assume that the assignment of processes to intervals is done
  only by the master process, which then tells each process
  what job it is to do.
*/
  if ( process_id == master )
  {
    for ( process = 1; process <= process_num-1; process++ )
    {
      xb[0] = ( ( double ) ( process_num - process     ) * x_min   
              + ( double ) (               process - 1 ) * x_max ) 
              / ( double ) ( process_num           - 1 );

      xb[1] = ( ( double ) ( process_num - process - 1 ) * x_min   
              + ( double ) (               process     ) * x_max ) 
              / ( double ) ( process_num           - 1 );
 
      target = process;
      tag = 1;

      ierr = MPI_Send ( xb, 2, MPI_DOUBLE, target, tag, MPI_COMM_WORLD );
    }
  }
  else
  {
    source = master;
    tag = 1;

    ierr = MPI_Recv ( xb, 2, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status );
  }
/*
  Wait here until everyone has gotten their assignment.
*/
  ierr = MPI_Barrier ( MPI_COMM_WORLD );

  if ( process_id == master )
  {
    printf ( "\n" );
    printf ( "INTERVALS - Master process:\n" );
    printf ( "  Subintervals have been assigned.\n" );
  }
/*
  Every process needs to be told the number of points to use.
  Since this is the same value for everybody, we use a broadcast.
  Again, we are doing it in this roundabout way to emphasize that
  the choice for M could really be made at runtime, by processor 0,
  and then sent out to the others.
*/
  m = 100;
  source = master;

  ierr = MPI_Bcast ( &m, 1, MPI_INT, source, MPI_COMM_WORLD );
/*
  Now, every process EXCEPT 0 computes its estimate of the 
  integral over its subinterval, and sends the result back
  to process 0.
*/
  if ( process_id != master )
  {  
    q_local = 0.0;

    for ( i = 1; i <= m; i++ )
    {
      x = ( ( double ) ( 2 * m - 2 * i + 1 ) * xb[0]   
          + ( double ) (         2 * i - 1 ) * xb[1] ) 
          / ( double ) ( 2 * m             );

      q_local = q_local + f ( x );
    }

    q_local = q_local * ( xb[1] - xb[0] ) / ( double ) ( m );

    target = master;
    tag = 2;

    ierr = MPI_Send ( &q_local, 1, MPI_DOUBLE, target, tag, MPI_COMM_WORLD );
  }
/*
  Process 0 expects to receive N-1 partial results.
*/
  else
  {
    received = 0;
    q_global = 0.0;

    while ( received < process_num - 1 )
    {
      source = MPI_ANY_SOURCE;
      tag = 2;

      ierr = MPI_Recv ( &q_local, 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, 
        &status );

      q_global = q_global + q_local;
      received = received + 1;
    }
  }
/*
  The master process prints the answer.
*/
  if ( process_id == master )
  {
    printf ( "\n" );
    printf ( "INTERVALS - Master process:\n" );
    printf ( "  Estimate for PI is %f\n", q_global );
    printf ( "  Error is           %f\n", q_global - pi );

    end_time = MPI_Wtime ( );

    printf ( "\n" );
    printf ( "  Elapsed wall clock seconds = %f\n", 
      end_time - start_time );
  }
/*
  Terminate MPI.
*/
  ierr = MPI_Finalize ( );
/*
  Termiante.
*/
  if ( process_id == master ) 
  {
    printf ( "\n" );
    printf ( "INTERVALS - Master process:\n" );         
    printf ( "  Normal end of execution.\n" );
    printf ( "\n" );
    timestamp ( );     
  }
  return 0;
}
/******************************************************************************/

double f ( double x )

/******************************************************************************/
/*
  Purpose:

    F evaluates the function F(X) which we are integrating.

  Discussion:

    Integral ( 0 <= X <= 1 ) 1/(1+X*X) dX = PI/4

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 October 2007

  Author:

    John Burkardt
*/
{
  double value;

  value = 4.0 / ( 1.0 + x * x );

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
