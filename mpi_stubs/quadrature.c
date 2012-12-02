# include <stdio.h>
# include <math.h>
# include <time.h>

# include "mpi_stubs_c.h"

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

    QUADRATURE uses MPI to estimate an integral.

    This program demonstrates the use of the MPI library to carry out
    an algorithm in parallel.

    In this case, we intend to compute the value of PI by means of the
    formula:

      PI = Integral ( 0 <= X <= 1 ) 4 / ( 1 + X * X ).

    This routine carries out this integration in parallel by dividing up
    the integration interval into subintervals, and having each process
    compute a subinterval estimate.  The master process adds these
    estimates together to get the estimate for the total integral, 
    and hence for PI.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    10 October 2002

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
  double endwtime;
  double h;
  int i;
  int ierr;
  int master = 0;
  int my_id;
  double mypi;
  int n;
  int numprocs;
  double pi;
  double PI25DT = 3.141592653589793238462643;
  double startwtime;
  double sum;
  double x;
/*
  Initialize MPI.
*/
  ierr = MPI_Init ( &argc, &argv );
/*
  Get the number of processes.
*/
  ierr = MPI_Comm_size ( MPI_COMM_WORLD, &numprocs );
/*
  Determine this processes's rank.
*/
  ierr = MPI_Comm_rank ( MPI_COMM_WORLD, &my_id );
/*
  Say hello.
*/
  if ( my_id == master ) 
  {
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
    printf ( "  The number of processes available is %d.\n", numprocs );
  }
/*
  Record the starting time.
*/
  if ( my_id == master ) 
  {
    startwtime = MPI_Wtime();
  }
/*
  Normally, we will assume, the number of intervals N is read in
  by process 0.  We'll just simulate this by assigning a value to
  N, but only on process 0!
*/
  if ( my_id == master )
  {
    n = 100;
    printf ( "\n" );
    printf ( "QUADRATURE - Master process:\n" );
    printf ( "  Number of intervals used is %d.\n", n );
  }
/*
  Broadcast the number of intervals to the other processes.
*/
  ierr = MPI_Bcast ( &n, 1, MPI_INT, master, MPI_COMM_WORLD );
/*
  Integrate F(X) over a subinterval determined by your process ID.
*/
  h = 1.0 / (double) n;
  sum = 0.0;
  for ( i = my_id + 1; i <= n; i = i + numprocs ) 
  {
    x = h * ( (double) i - 0.5 );
    sum = sum + f(x);
  }
  mypi = h * sum;

  printf ( "\n" );
  printf ( "QUADRATURE - Process %d:\n", my_id );
  printf ( "  My contribution to integral is %f.\n", mypi );
/*
  Send your local result (MYPI) to the MASTER process, to be added to
  the global result (PI).
*/
  ierr = MPI_Reduce ( &mypi, &pi, 1, MPI_DOUBLE, MPI_SUM, master, 
    MPI_COMM_WORLD );

  if ( my_id == master )
  {
    printf ( "\n");         
    printf ( "QUADRATURE - Master process:\n");         
    printf ( "  The estimate for PI is %.16f\n", pi );
    printf ( "  The error of this estimate is %.16f\n", fabs ( pi - PI25DT) );

    endwtime = MPI_Wtime();
    printf ( "\n");       
    printf ( "  Wall clock elapsed seconds = %f\n", endwtime-startwtime );      
  }
/*
  Terminate MPI.
*/
  ierr = MPI_Finalize();
/*
  Terminate.
*/
  if ( my_id == master ) 
  {
    printf ( "\n");         
    printf ( "QUADRATURE - Master process:\n");         
    printf ( "  Normal end of execution.\n");         
  }

  return 0;
}
/******************************************************************************/

double f ( double x )

/******************************************************************************/
/*
  Purpose:

    F evaluates the function F(X) which we are integrating.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    21 March 2009

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
  size_t len;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  len = strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  printf ( "%s\n", time_buffer );

  return;
# undef TIME_SIZE
}
