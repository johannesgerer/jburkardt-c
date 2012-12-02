# include <stdio.h>
# include <math.h>
# include <time.h>

# include "mpi.h"


# define DEBUG             0
# define CHUNKSIZE      1000
# define RANDOM_SEED       0

/* 
  Message tags
*/
# define NEED_NUMBERS      1
# define RANDOM_NUMBERS    2

int main ( int argc, char *argv[] );
void timestamp ( void );

/******************************************************************************/

int main ( int argc, char *argv[] )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for MONTE_CARLO.

  Discussion:

    MONTE_CARLO uses Monte Carlo methods to estimate Pi.

    Generate N random points in the unit square.  Count M, the number
    of points that are in the quarter circle.  Then PI is approximately
    equal to the ratio 4 * M / N.

    It's important that each processor use DIFFERENT random numbers.
    One way to ensure this is to have a single master processor
    generate all the random numbers, and then divide them up.

    (A second way, not explored here, is simply to ensure that each
    processor uses a different seed, either chosen by a master processor,
    or generated from the processor ID.)

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    26 February 2007

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
  double calculatedPi;
  int done;
  double error;
  int i;
  int ierr;
  int in;
  int max;
  MPI_Status mesgStatus;
  int my_id;
  int numprocs;
  int out;
  int point_max = 1000000;
  int randServer;
  int randNums[CHUNKSIZE];
  int ranks[1];
  int request;
  int temp;
  double tolerance;
  int totalin;
  int totalout;
  MPI_Group worker_group;
  MPI_Comm workers;
  MPI_Group world_group;
  double wtime;
  double x;
  double y;
/*
  Initialize MPI.
*/
  ierr = MPI_Init ( &argc, &argv );
/*
  Get the number of processors.
*/
  ierr = MPI_Comm_size ( MPI_COMM_WORLD, &numprocs );
/*
  Get the rank of this processor.
*/
  ierr = MPI_Comm_rank ( MPI_COMM_WORLD, &my_id );

  if ( my_id == 0 ) 
  {
    timestamp ( );
    printf ( "\n" );
    printf ( "MONTE_CARLO - Master process:\n" );
    printf ( "  C version\n" );
    printf ( "  An MPI example program.\n" );
    printf ( "  Estimate pi by the Monte Carlo method, using MPI.\n" );
    printf ( "\n" );
    printf ( "  Compiled on %s at %s.\n", __DATE__, __TIME__ );
    printf ( "\n" );
    printf ( "  The number of processes is %d.\n", numprocs );
    printf ( "\n" );
    printf ( "  Points in the unit square will be tested\n" );
    printf ( "  to see if they lie in the unit quarter circle.\n" );
  }

  if ( my_id == 0 )
  {
    wtime = MPI_Wtime ( );
  }
/*
  Pretend that the tolerance TOLERANCE is supplied externally
  to the master process, which must then broadcast it to all
  other processes.
*/
  if ( my_id == 0 )
  {
    tolerance = 0.0001;

    printf ( "\n" );
    printf ( "  The method will continue to improve the estimate until:\n" );
    printf ( "  PI is computed to within a tolerance = %f,\n", tolerance );
    printf ( "  or the number of points examined reaches %d.\n", point_max );
  }

  ierr = MPI_Bcast ( &tolerance, 1, MPI_DOUBLE_PRECISION, 0,
    MPI_COMM_WORLD );

  printf ( "  Process %d is active.\n", my_id );
/*
  Start by getting the group corresponding to the world communicator.
*/
  ierr = MPI_Comm_group ( MPI_COMM_WORLD, &world_group );  
/*
  Put SERVER on the list of processes to exclude, and create the new
  worker group.
*/
  randServer = numprocs-1;
  ranks[0] = randServer;
  ierr = MPI_Group_excl ( world_group, 1, ranks, &worker_group );
/* 
  Use the worker group to create the new worker communicator.
*/
  ierr = MPI_Comm_create ( MPI_COMM_WORLD, worker_group, &workers );
/*
  Since we only needed the worker group to create the worker
  communicator, we can free the worker group now.
*/
  ierr = MPI_Group_free ( &worker_group );
/*
  Here is where the computation is carried out.
*/

/* 
  I am the rand server.
*/
  if ( my_id == randServer ) 
  {
# if RANDOM_SEED
    struct timeval time;
    gettimeofday( &time, 0 );
/* 
  Initialize the random number generator 
*/
    srandom ( (int)(time.tv_usec*1000000+time.tv_sec) );
# endif
    do
    {
      ierr = MPI_Recv ( &request, 1, MPI_INT, MPI_ANY_SOURCE, NEED_NUMBERS,
        MPI_COMM_WORLD, &mesgStatus );

      if ( request ) 
      {
        for ( i = 0; i < CHUNKSIZE; i++) 
        {
          randNums[i] = random();
        }
        ierr = MPI_Send ( randNums, CHUNKSIZE, MPI_INT, 
          mesgStatus.MPI_SOURCE, RANDOM_NUMBERS, MPI_COMM_WORLD );
      }
    } while ( 0 < request );

  }
/* 
  I am a worker process.
*/
  else  
  {
    request = 1;
    done = in = out = 0;
    max = 2147483647;

    ierr = MPI_Send ( &request, 1, MPI_INT, randServer, NEED_NUMBERS, 
      MPI_COMM_WORLD );
/* 
  Request a string of random numbers.
*/
    while (!done) 
    {
      request = 1;
      ierr = MPI_Recv ( randNums, CHUNKSIZE, MPI_INT, randServer,
        RANDOM_NUMBERS, MPI_COMM_WORLD, &mesgStatus );

      for ( i = 0; i < CHUNKSIZE; ) 
      {
        x = ( ( float ) randNums[i++] ) / max;
        y = ( ( float ) randNums[i++] ) / max;

        if ( x * x + y * y < 1.0 ) 
        {
          in++;
        } 
        else 
        {
          out++;
        }

      }

      temp = in;
      ierr = MPI_Reduce ( &temp, &totalin, 1, MPI_INT, MPI_SUM, 0, workers );
/* 
  Count total of ins.
*/ 
      temp = out;
      ierr = MPI_Reduce ( &temp, &totalout, 1, MPI_INT, MPI_SUM, 0, workers );
/* 
  Count total of outs.
*/
      if ( my_id == 0 ) 
      {
        calculatedPi = ( 4.0 * totalin ) / ( totalin + totalout );
        error = fabs ( calculatedPi - 3.141592653589793238462643 );
        done = ( error < tolerance ) || point_max <= ( totalin + totalout );
        printf( "pi = %23.20lf\n", calculatedPi );

        if ( done )
        {
          request = 0;
        }
        else
        {
          request = 1;
        }

        ierr = MPI_Send ( &request, 1, MPI_INT, randServer, NEED_NUMBERS,
          MPI_COMM_WORLD );

        ierr = MPI_Bcast ( &done, 1, MPI_INT, 0, workers );
      } 
      else
      {
        ierr = MPI_Bcast ( &done, 1, MPI_INT, 0, workers );

        if ( !done ) 
        {
          request = 1;
          ierr = MPI_Send ( &request, 1, MPI_INT, randServer, NEED_NUMBERS,
            MPI_COMM_WORLD );
        }
      }
    }
  }

  if ( my_id == 0 ) 
  {
    printf( "\npoints: %d\nin: %d, out: %d\n", totalin + totalout, totalin,
      totalout );

    wtime = MPI_Wtime ( ) - wtime;
    printf ( "\n" );
    printf ( "  Elapsed wallclock time = %g seconds.\n", wtime );
  }
/*
  Terminate MPI.
*/
  ierr = MPI_Finalize();
/*
  Terminate.
*/
  if ( my_id == 0 )
  {
    printf ( "\n" );
    printf ( "MONTE_CARLO - Master process:\n" );
    printf ( "  Normal end of execution.\n" );
    printf ( "\n" );
    timestamp ( );
  }
  return 0;
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

