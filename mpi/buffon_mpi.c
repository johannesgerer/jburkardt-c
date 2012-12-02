# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>

# include "mpi.h"

int main ( int argc, char *argv[] );
int buffon_laplace_simulate ( double a, double b, double l, int trial_num );
double r8_abs ( double x );
double r8_huge ( void );
void timestamp ( void );

/******************************************************************************/

int main ( int argc, char *argv[] )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for BUFFON.

  Discussion:

    This program uses MPI to do a Buffon-Laplace simulation in parallel.

    This is an example of an "embarassingly parallel" computation.  Each
    processor does the same randomized computation, counting "successes",
    and returning the number of successes to the master process.

    The particular point made by this example is that each process must
    use a different stream of random numbers, and that this can be done
    if each process uses a different seed to initialize the random number
    generator, and that this, in turn, can be done by simply adding the
    process's rank to a common seed value.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    28 February 2007

  Author:

    John Burkardt

  Reference:

    William Gropp, Ewing Lusk, Anthony Skjellum,
    Using MPI: Portable Parallel Programming with the
    Message-Passing Interface,
    Second Edition,
    MIT Press, 1999,
    ISBN: 0262571323.

    Sudarshan Raghunathan,
    Making a Supercomputer Do What You Want: High Level Tools for
    Parallel Programming,
    Computing in Science and Engineering,
    Volume 8, Number 5, September/October 2006, pages 70-80.
*/
{
  double a = 1.0;
  double b = 1.0;
  int hit_num;
  int hit_total;
  int ierr;
  double l = 1.0;
  int master = 0;
  double pdf_estimate;
  double pi = 3.141592653589793238462643;
  double pi_error;
  double pi_estimate;
  int process_num;
  int process_rank;
  double random_value;
  int seed;
  int trial_num = 100000;
  int trial_total;
/*
  Initialize MPI.
*/
  ierr = MPI_Init ( &argc, &argv );
/*
  Get the number of processes.
*/
  ierr = MPI_Comm_size ( MPI_COMM_WORLD, &process_num );
/*
  Get the rank of this process.
*/
  ierr = MPI_Comm_rank ( MPI_COMM_WORLD, &process_rank );
/*
  The master process prints a message.
*/
  if ( process_rank == master ) 
  {
    timestamp ( );
    printf ( "\n" );
    printf ( "BUFFON - Master process:\n" );
    printf ( "  C version\n" );
    printf ( "  An MPI example program to estimate PI\n" );
    printf ( "  using the Buffon-Laplace needle experiment.\n" );
    printf ( "  On a grid of cells of  width A and height B,\n" );
    printf ( "  a needle of length L is dropped at random.\n" );
    printf ( "  We count the number of times it crosses\n" );
    printf ( "  at least one grid line, and use this to estimate \n" );
    printf ( "  the value of PI.\n" );
    printf ( "\n" );
    printf ( "  The number of processes is %d\n", process_num );
    printf ( "\n" );
    printf ( "  Cell width A =    %f\n", a );
    printf ( "  Cell height B =   %f\n", b );
    printf ( "  Needle length L = %f\n", l );
  }
/*
  Each process sets a random number seed.
*/
  seed = 123456789 + process_rank * 100;
  srand ( seed );
/*
  Just to make sure that we're all doing different things, have each
  process print out its rank, seed value, and a first test random value.
*/
  random_value  = ( double ) rand ( ) / ( double ) RAND_MAX;

  printf ( "  %8d  %12d  %14f\n", process_rank, seed, random_value );
/*
  Each process now carries out TRIAL_NUM trials, and then
  sends the value back to the master process.
*/
  hit_num = buffon_laplace_simulate ( a, b, l, trial_num );

  ierr = MPI_Reduce ( &hit_num, &hit_total, 1, MPI_INT, MPI_SUM, master, 
    MPI_COMM_WORLD );
/*
  The master process can now estimate PI.
*/
  if ( process_rank == master )
  {
    trial_total = trial_num * process_num;

    pdf_estimate = ( double ) ( hit_total ) / ( double ) ( trial_total );

    if ( hit_total == 0 )
    {
      pi_estimate = r8_huge ( );
    }
    else
    {
      pi_estimate = l * ( 2.0 * ( a + b ) - l ) / ( a * b * pdf_estimate );
    }

    pi_error = r8_abs ( pi - pi_estimate );

    printf ( "\n" );
    printf ( 
      "    Trials      Hits    Estimated PDF       Estimated Pi        Error\n" );
    printf ( "\n" );
    printf ( "  %8d  %8d  %16.12f  %16.12f  %16.12f\n",
      trial_total, hit_total, pdf_estimate, pi_estimate, pi_error );
  }
/*
  Terminate MPI.
*/
  ierr = MPI_Finalize ( );
/*
  Terminate.
*/
  if ( process_rank == master )
  {
    printf ( "\n" );
    printf ( "BUFFON - Master process:\n" );
    printf ( "  Normal end of execution.\n" );
    printf ( "\n" );
    timestamp ( );
  }
  return 0;
}
/******************************************************************************/

int buffon_laplace_simulate ( double a, double b, double l, int trial_num )

/******************************************************************************/
/*
  Purpose:

    BUFFON_LAPLACE_SIMULATE simulates a Buffon-Laplace needle experiment.

  Discussion:

    In the Buffon-Laplace needle experiment, we suppose that the plane has
    been tiled into a grid of rectangles of width A and height B, and that a
    needle of length L is dropped "at random" onto this grid.  
 
    We may assume that one end, the "eye" of the needle falls at the point
    (X1,Y1), taken uniformly at random in the cell [0,A]x[0,B].

    ANGLE, the angle that the needle makes is taken to be uniformly random.
    The point of the needle, (X2,Y2), therefore lies at

      (X2,Y2) = ( X1+L*cos(ANGLE), Y1+L*sin(ANGLE) )

    The needle will have crossed at least one grid line if any of the 
    following are true:

      X2 <= 0, A <= X2, Y2 <= 0, B <= Y2.

    This routine simulates the tossing of the needle, and returns the number
    of times that the needle crossed at least one grid line.

    If L is larger than sqrt ( A*A + B*B ), then the needle will
    cross every time, and the computation is uninteresting.  However, if
    L is smaller than this limit, then the probability of a crossing on
    a single trial is

      P(L,A,B) = ( 2 * L * ( A + B ) - L * L ) / ( PI * A * B )

    and therefore, a record of the number of hits for a given number of
    trials can be used as a very roundabout way of estimating PI.  
    (Particularly roundabout, since we actually will use a good value of
    PI in order to pick the random angles//)

    Since this routine invokes the C random number generator,
    the user should initialize the random number generator, particularly
    if it is desired to control whether the sequence is to be varied
    or repeated.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    28 February 2007

  Author:

    John Burkardt

  Reference:

    Sudarshan Raghunathan,
    Making a Supercomputer Do What You Want: High Level Tools for 
    Parallel Programming,
    Computing in Science and Engineering,
    Volume 8, Number 5, September/October 2006, pages 70-80.

  Parameters:

    Input, double A, B, the horizontal and vertical dimensions
    of each cell of the grid.  0 <= A, 0 <= B.

    Input, double L, the length of the needle.
    0 <= L <= min ( A, B ).

    Input, int TRIAL_NUM, the number of times the needle is
    to be dropped onto the grid.

    Output, int BUFFON_LAPLACE_SIMULATE, the number of times the needle 
    crossed at least one line of the grid of cells.
*/
{
  double angle;
  int hits;
  double pi = 3.141592653589793238462643;
  int trial;
  double x1;
  double x2;
  double y1;
  double y2;

  hits = 0;

  for ( trial = 1; trial <= trial_num; trial++ )
  {
/*
  Randomly choose the location of the eye of the needle in [0,0]x[A,B],
  and the angle the needle makes.
*/
    x1 = a * ( double ) rand ( ) / ( double ) RAND_MAX;
    y1 = b * ( double ) rand ( ) / ( double ) RAND_MAX;
    angle = 2.0 * pi * ( double ) rand ( ) / ( double ) RAND_MAX;
/*
  Compute the location of the point of the needle.
*/
    x2 = x1 + l * cos ( angle );
    y2 = y1 + l * sin ( angle );
/*
  Count the end locations that lie outside the cell.
*/
    if ( x2 <= 0.0 || a <= x2 || y2 <= 0.0 || b <= y2 )
    {
      hits = hits + 1;
    }
  }
  return hits;
}
/******************************************************************************/

double r8_abs ( double x )

/******************************************************************************/
/*
  Purpose:

    R8_ABS returns the absolute value of an R8.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    07 May 2006

  Author:

    John Burkardt

  Parameters:

    Input, double X, the quantity whose absolute value is desired.

    Output, double R8_ABS, the absolute value of X.
*/
{
  if ( 0.0 <= x )
  {
    return x;
  } 
  else
  {
    return ( -x );
  }
}
/******************************************************************************/

double r8_huge ( void )

/******************************************************************************/
/*
  Purpose:

    R8_HUGE returns a "huge" R8.

  Discussion:

    The value returned by this function is NOT required to be the
    maximum representable R8.  This value varies from machine to machine,
    from compiler to compiler, and may cause problems when being printed.
    We simply want a "very large" but non-infinite number.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 October 2007

  Author:

    John Burkardt

  Parameters:

    Output, double R8_HUGE, a "huge" R8 value.
*/
{
  double value;

  value = 1.0E+30;

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
