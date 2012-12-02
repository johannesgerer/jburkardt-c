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

    MAIN is the main program for QUADRATURE.

  Discussion:

    QUADRATURE estimates an integral using quadrature.

    The integral of F(X) = 4 / ( 1 + X * X ) from 0 to 1 is PI.

    We break up the interval [0,1] into N subintervals, evaluate
    F(X) at the midpoint of each subinterval, and multiply the
    sum of these values by N to get an estimate for the integral.

    If we have M processes available because we are using MPI, then
    we can ask processes 0, 1, 2, ... M-1 to handle the subintervals
    in the following order:

          0      1       2            M-1  <-- Process numbers begin at 0
     ------ ------  ------  -----  ------
          1      2       3    ...       M
        M+1    M+2     M+3    ...     2*M
      2*M+1    2*M+2 2*M+3    ...     3*M
                              
    and so on up to subinterval N.  The partial sums collected by 
    each process are then sent to the master process to be added 
    together to get the estimated integral.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    22 August 2008

  Author:

    John Burkardt

  Reference:

    William Gropp, Ewing Lusk, Anthony Skjellum,
    Using MPI: Portable Parallel Programming with the
    Message-Passing Interface,
    Second Edition,
    MIT Press, 1999,
    ISBN: 0262571323.
*/
{
  double h;
  int i;
  int id;
  int ierr;
  int master = 0;
  int n;
  int n_part;
  int p;
  double q;
  double q_diff;
  double q_exact = 3.141592653589793238462643;
  double q_part;
  double wtime_diff;
  double wtime_end;
  double wtime_start;
  double x;
/*
  Initialize MPI.
*/
  ierr = MPI_Init ( &argc, &argv );
/*
  Get the number of processes.
*/
  ierr = MPI_Comm_size ( MPI_COMM_WORLD, &p );
/*
  Determine this processes's rank.
*/
  ierr = MPI_Comm_rank ( MPI_COMM_WORLD, &id );
/*
  Say hello.
*/
  if ( id == master ) 
  {
    timestamp ( );
    printf ( "\n" );
    printf ( "QUADRATURE - Master process:\n" );
    printf ( "  C version\n" );
    printf ( "\n" );
    printf ( "  An MPI example program.\n" );
    printf ( "  Estimate the value of PI by approximating an integral.\n" );
    printf ( "  The integral is approximated by a sum,\n" );
    printf ( "  whose calculation is divided among a number of processes.\n" );
    printf ( "\n" );
    printf ( "  Compiled on %s at %s.\n", __DATE__, __TIME__ );
    printf ( "\n" );
    printf ( "  The number of processes available is %d.\n", p );
  }
/*
  Record the starting time.
*/
  if ( id == master ) 
  {
    wtime_start = MPI_Wtime ( );
  }
/*
  We could assume that the number of intervals N is read in by process 0.
  To keep this program simple, we'll just assign a value to
  N, but only on process 0.
*/
  if ( id == master )
  {
    n = 1000;
    printf ( "\n" );
    printf ( "QUADRATURE - Master process:\n" );
    printf ( "  Number of intervals used is %d.\n", n );
  }
/*
  The master process broadcasts, and the other processes receive,
  the number of intervals N.
*/
  ierr = MPI_Bcast ( &n, 1, MPI_INT, master, MPI_COMM_WORLD );
/*
  Every process integrates F(X) over a subinterval determined by its process ID.
*/
  h = 1.0 / ( double ) n;

  q_part = 0.0;
  n_part = 0;

  for ( i = id + 1; i <= n; i = i + p ) 
  {
    x = ( double ) ( 2 * i - 1 )
      / ( double ) ( 2 * n     );

    n_part = n_part + 1;
    q_part = q_part + f ( x );
  }
  q_part = q_part * h;

  printf ( "\n" );
  printf ( "QUADRATURE - Process %d:\n", id );
  printf ( "  Points used = %d\n", n_part );
  printf ( "  Estimate %f.\n", q_part );
/*
  Each process sends its local result Q_PART to the MASTER process, 
  to be added to the global result QI.
*/
  ierr = MPI_Reduce ( &q_part, &q, 1, MPI_DOUBLE, MPI_SUM, master, 
    MPI_COMM_WORLD );
/*
  The master process scales the sum and reports the results.
*/
  if ( id == master )
  {
    printf ( "\n");         
    printf ( "QUADRATURE - Master process:\n");         
    printf ( "  Integral estimate      %24.16f\n", q );
    printf ( "  The exact value is     %24.16f\n", q_exact );
    q_diff = fabs ( q - q_exact );
    printf ( "  The error is           %24.16f\n", q_diff );

    wtime_end = MPI_Wtime ( );
    wtime_diff = wtime_end - wtime_start;

    printf ( "\n");       
    printf ( "  Wall clock elapsed seconds = %f\n", wtime_diff );      
  }
/*
  Terminate MPI.
*/
  ierr = MPI_Finalize ( );
/*
  Terminate.
*/
  if ( id == master ) 
  {
    printf ( "\n");         
    printf ( "QUADRATURE - Master process:\n");         
    printf ( "  Normal end of execution.\n");
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

    Integral ( 0 <= X <= 1 ) 4/(1+X*X) dX = PI

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    15 October 2007

  Author:

    John Burkardt

  Parameters:

    Input, double X, the point at which to evaluate F.

    Output, double F, the value of F(X).
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
