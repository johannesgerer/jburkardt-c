# include <stdio.h>
# include <stdlib.h>
# include <time.h>

# include "mpi.h"

int main ( int argc, char *argv[] );
void timestamp ( void );

/******************************************************************************/

int main ( int argc, char *argv[] )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for DAY1.

  Discussion:

    DAY1 is exercise 3 for first day of the MPI workshop

    The instructions say:

    Process 1 computes the squares of the first 200 integers.
    It sends this data to process 3.

    Process 3 should divide the integers between 20 and 119 by 53,
    getting a real result, and passes this data back to process 1.

    * I presume the first 200 integers are the numbers 0 through 199.

    * The instructions literally mean that process 3 should look
      at integers whose VALUES are between 20 and 119.  I doubt that
      is what the instructor meant, but it's more interesting than
      simply picking the entries with index between 20 and 119,
      so that's what I'll do.

    * It is also not completely clear whether only the selected data
      should be sent back, or the entire array.  Again, it is more
      interesting to send back only part of the data.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    08 October 2002

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
  int i_dim = 200;
  int r_dim = 200;

  int count;
  int count2;
  int dest;
  int i;
  int i_buffer[i_dim];
  int ierr;
  int num_procs;
  float r_buffer[r_dim];
  int rank;
  int source;
  MPI_Status status;
  int tag;
/*
  Initialize MPI.
*/
  ierr = MPI_Init ( &argc, &argv );
/*
  Determine this process's rank.
*/
  ierr = MPI_Comm_rank ( MPI_COMM_WORLD, &rank );
/*
  Have Process 0 say hello.
*/
  if ( rank == 0 )
  {
    timestamp ( );
    printf ( "\n" );
    printf ( "DAY1:\n" );
    printf ( "  C version\n" );
    printf ( "  An MPI example program.\n" );
    printf ( "  MPI exercise #3 for day 1.\n" );
    printf ( "\n" );
    printf ( "  Compiled on %s at %s.\n", __DATE__, __TIME__ );
    printf ( "\n" );
    printf ( "  The number of processes available is %d.\n", num_procs );
  }
/*
  Get the number of processes.
*/
  ierr = MPI_Comm_size ( MPI_COMM_WORLD, &num_procs );
/*
  If we don't have at least 4 processes, then bail out now.
*/
  if ( num_procs < 4 )
  {
    printf ( "\n" );
    printf ( "DAY1 - Process %d.\n", rank );
    printf ( "  Not enough processes for this task!\n" );
    printf ( "  Bailing out now!\n" );
    ierr = MPI_Finalize ( );
    return 1;
  }
/*
  Process 1 knows that it will generate 200 integers, and may receive no more
  than 200 reals.
*/
  if ( rank == 1 )
  {
    count = 200;

    for ( i = 0; i < count; i++ ) 
    {
      i_buffer[i] = i * i;
    }

    dest = 3;
    tag = 1;

    ierr = MPI_Send ( i_buffer, count, MPI_INT, dest, tag, MPI_COMM_WORLD );

    printf ( "P:%d sent %d integers to process %d.\n", rank, count, dest );

    source = 3;
    tag = 2;

    ierr = MPI_Recv ( r_buffer, r_dim, MPI_FLOAT, source, tag, 
      MPI_COMM_WORLD, &status );

    printf ( "P:%d received real values from process 3.\n", rank );

    ierr = MPI_Get_count ( &status, MPI_FLOAT, &count );

    printf ( "P:%d Number of real values received is %d.\n", rank, count );

    printf ( "P:%d First 3 values = %f %f %f\n", rank, r_buffer[0], 
      r_buffer[1], r_buffer[2] );
  }
/*
  Process 3 receives the integer data from process 1, selects some of the data, does
  a real computation on it, and sends that part back to process 1.
*/
  else if ( rank == 3 ) 
  {
    source = 1;
    tag = 1;

    ierr = MPI_Recv ( i_buffer, i_dim, MPI_INT, source, tag, 
      MPI_COMM_WORLD, &status );

    printf ( "\n" );
    printf ( "P:%d received integer values from process 1.\n", rank );

    ierr = MPI_Get_count ( &status, MPI_INT, &count );

    printf ( "P:%d - Number of integers received is %d.\n", rank, count );
    printf ( "P:%d First 3 values = %d %d %d.\n", 
      rank, i_buffer[0], i_buffer[1], i_buffer[2] );

    count2 = 0;
     
    for ( i = 0; i < count; i++ ) 
    {
      if ( 20 <= i_buffer[i] && i_buffer[i] <= 119 ) 
      {

        r_buffer[count2] = ( float ) i_buffer[i] / 53.0E+00;
        count2 = count2 + 1;

        if ( count2 <= 3 ) 
        {
          printf ( "P:%d Input integer %d becomes %f.\n", 
            rank, i_buffer[i], r_buffer[count2-1] );
        }

      }
    }

    dest = 1;
    tag = 2;
  
    ierr = MPI_Send ( r_buffer, count2, MPI_FLOAT, dest, tag, MPI_COMM_WORLD );

    printf ( "P:%d sent %d reals to process %d.\n", rank, count2, dest );
  }
  else
  {
    printf ( "\n" );
    printf ( "P:%d - MPI has no work for me!\n", rank );
  }
/*
  Terminate MPI.
*/
  ierr = MPI_Finalize ( );
/*
  Terminate.
*/
  if ( rank == 0 )
  {
    printf ( "\n" );
    printf ( "DAY1:\n" );
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

