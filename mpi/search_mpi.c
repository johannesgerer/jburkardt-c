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

    MAIN is the main program for SEARCH.

  Discussion:

    SEARCH demonstrates the use of MPI routines to carry out a search

    An array of given size is to be searched for occurrences of a
    specific value.

    The search is done in parallel.  A master process generates the
    array and the target value, then distributes the information among
    a set of worker processes, and waits for them to communicate back
    the (global) index values at which occurrences of the target value
    were found.

    An interesting feature of this program is the use of allocatable
    arrays, which allows the master program to set aside just enough
    memory for the whole array, and for each worker program to set aside
    just enough memory for its own part of the array.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    09 October 2002

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
  int *a;
  int dest;
  float factor;
  int global;
  int i;
  int ierr;
  int master = 0;
  int my_id;
  int n;
  int npart;
  int num_procs;
  int source;
  int start;
  MPI_Status status;
  int tag;
  int tag_target = 1;
  int tag_size = 2;
  int tag_data = 3;
  int tag_found = 4;
  int tag_done = 5;
  int target;
  int workers_done;
  int x;
/*
  Initialize MPI.
*/
  ierr = MPI_Init ( &argc, &argv );
/*
  Get this processes's rank.
*/
  ierr = MPI_Comm_rank ( MPI_COMM_WORLD, &my_id );
/*
  Find out how many processes are available.
*/
  ierr = MPI_Comm_size ( MPI_COMM_WORLD, &num_procs );

  if ( my_id == master )
  {
    timestamp ( );
    printf ( "\n" );
    printf ( "SEARCH - Master process:\n" );
    printf ( "  C version\n" );
    printf ( "  An MPI example program to search an array.\n" );
    printf ( "\n" );
    printf ( "  Compiled on %s at %s.\n", __DATE__, __TIME__ );
    printf ( "\n" );
    printf ( "  The number of processes is %d,\n", num_procs );
  }

  printf ( "\n" );
  printf ( "Process %d is active.\n", my_id );
/*
  Have the master process generate the target and data.  In a more 
  realistic application, the data might be in a file which the master 
  process would read.  Here, the master process decides.
*/
  if ( my_id == master )
  {
/*
  Pick the number of data items per process, and set the total.
*/
    factor = ( float ) rand ( ) / ( float ) RAND_MAX;
    npart = 50 + ( int ) ( factor * 100.0E+00 );
    n = npart * num_procs;

    printf ( "\n" );
    printf ( "SEARCH - Master process:\n" );
    printf ( "  The number of data items per process is %d,\n", npart );
    printf ( "  The total number of data items is       %d.\n", n );
/*
  Now allocate the master copy of A, fill it with values, and pick 
  a value for the target.
*/
    a = malloc ( n * sizeof ( int ) );

    factor = ( float ) n / 10.0E+00 ;

    for ( i = 0; i < n; i++ ) 
    {
      a[i] = ( int ) ( factor * ( float ) rand ( ) / ( float ) RAND_MAX );
    }
    target = a[n/2];

    printf ( "  The target value is %d.\n", target );
/*
  The worker processes need to have the target value, the number of data items,
  and their individual chunk of the data vector.
*/
    for ( i = 1; i <= num_procs-1; i++ )
    {
      dest = i;
      tag = tag_target;

      ierr = MPI_Send ( &target, 1, MPI_INT, dest, tag, MPI_COMM_WORLD );

      tag = tag_size;

      ierr = MPI_Send ( &npart, 1, MPI_INT, dest, tag, MPI_COMM_WORLD );

      start = ( i - 1 ) * npart;
      tag = tag_data;

      ierr = MPI_Send ( a+start, npart, MPI_INT, dest, tag,
        MPI_COMM_WORLD );
    }
/*
  Now the master process simply waits for each worker process to report that 
  it is done.
*/
    workers_done = 0;

    while ( workers_done < num_procs-1 )
    {
      ierr = MPI_Recv ( &x, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG,
        MPI_COMM_WORLD, &status );

      source = status.MPI_SOURCE;
      tag = status.MPI_TAG;
    
      if ( tag == tag_done )
      {
        workers_done = workers_done + 1;
      }
      else if ( tag == tag_found )
      {
        printf ( "P%d %d %d\n", source, x, a[x] );
      }
      else
      {
        printf ( "  Master process received message with unknown tag = %d.\n", tag );
      }

    }
/*
  The master process can throw away A now.
*/
    free ( a );
  }
/*
  Each worker process expects to receive the target value, the number of data
  items, and the data vector.
*/
  else 
  {
    source = master;
    tag = tag_target;

    ierr = MPI_Recv ( &target, 1, MPI_INT, source, tag, MPI_COMM_WORLD,
      &status );
 
    source = master;
    tag = tag_size;

    ierr = MPI_Recv ( &npart, 1, MPI_INT, source, tag, MPI_COMM_WORLD, 
      &status );

    a = malloc ( npart * sizeof ( int ) );

    source = master;
    tag = tag_data;

    ierr = MPI_Recv ( a, npart, MPI_INT, source, tag, MPI_COMM_WORLD,
      &status );
/*
  The worker simply checks each entry to see if it is equal to the target
  value.
*/
    for ( i = 0; i < npart; i++ )
    {
      if ( a[i] == target )
      {
        global = ( my_id - 1 ) * npart + i;
        dest = master;
        tag = tag_found;

        ierr = MPI_Send ( &global, 1, MPI_INT, dest, tag, MPI_COMM_WORLD );
      }

    }
/*
  When the worker is finished with the loop, it sends a dummy data value with
  the tag "TAG_DONE" indicating that it is done.
*/
    dest = master;
    tag = tag_done;

    ierr = MPI_Send ( &target, 1, MPI_INT, dest, tag, MPI_COMM_WORLD );

    free ( a );
     
  }
/*
  Terminate MPI.
*/
  ierr = MPI_Finalize ( );
/*
  Terminate.
*/
  if ( my_id == master )
  {
    printf ( "\n" );
    printf ( "SEARCH - Master process:\n" );
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
