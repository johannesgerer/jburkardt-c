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

    MAIN is the main program for TYPE.

  Discussion:

    TYPE demonstrates the use of a user-defined MPI datatype.

    The datatype defined will be a structure that contains three
    integers.

    Process 0 will set up an example of this structure, and send it
    to Proces 1, which will alter it and send it back.

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
  int dest;
  int i;
  int ierr;
  int master = 0;
  int my_id;
  int num_procs;
  struct 
  {
    int x;
    int y;
    int z;
  } point;
  MPI_Datatype point_type;
  int source;
  MPI_Status status;
  int tag;
/*
  Initialize MPI.
*/
  ierr = MPI_Init ( &argc, &argv );
/*
  Get the number of processes.
*/
  ierr = MPI_Comm_size ( MPI_COMM_WORLD, &num_procs );
/*
  Get the individual process ID.
*/
  ierr = MPI_Comm_rank ( MPI_COMM_WORLD, &my_id );
/*
  Print a message.
*/
  if ( my_id == master ) 
  {
    timestamp ( );
    printf ( "\n" );
    printf ( "TYPE - Master process:\n" );
    printf ( "  C version\n" );
    printf ( "  An MPI example program that uses an MPI datatype.\n" );
    printf ( "\n" );
    printf ( "  Compiled on %s at %s.\n", __DATE__, __TIME__ );
    printf ( "\n" );
    printf ( "  The number of processes is %d\n", num_procs );
  }
  printf ( "\n" );
  printf ( "  Process %d says 'Hello, world!'\n", my_id );
/*
  Define and commit the new datatype.
*/
  ierr = MPI_Type_contiguous ( 3, MPI_INT, &point_type );
  ierr = MPI_Type_commit ( &point_type );

  if ( my_id == master ) 
  {
    point.x = 1;
    point.y = 2;
    point.z = 4;
    dest = 1;
    tag = 1;

    ierr = MPI_Send ( &point, 1, point_type, dest, tag, MPI_COMM_WORLD );

    printf ( "\n" );
    printf ( "  Process %d sent an item of type POINT_TYPE\n" );
    printf ( "  with value %d %d %d.\n", my_id, 
      point.x, point.y, point.z );

    source = 1;
    tag = 2;
    ierr = MPI_Recv ( &point, 1, point_type, source, tag, MPI_COMM_WORLD,
      &status );

    printf ( "  Process %d received a modified item of type POINT_TYPE, with value %d %d %d.\n", my_id, 
      point.x, point.y, point.z );

  }
  else if ( my_id == 1 )
  {
    source = 0;
    tag = 1;

    ierr = MPI_Recv ( &point, 1, point_type, source, tag, MPI_COMM_WORLD,
      &status );

    printf ( "\n" );
    printf ( "  Process %d received an item of type POINT_TYPE, with value %d %d %d.\n", my_id, 
      point.x, point.y, point.z );

    i = point.x;
    point.x = point.z * 100;
    point.y = point.y * 10;
    point.z = i;
    dest = 0;
    tag = 2;

    ierr = MPI_Send ( &point, 1, point_type, dest, tag, MPI_COMM_WORLD );

    printf ( "  Process %d sent a modified item of type POINT_TYPE, with value %d %d %d.\n", my_id, 
      point.x, point.y, point.z );
  }
  else
  {
    printf ( "\n" );
    printf ( "  Process %d: MPI has nothing for me to do!\n", my_id );
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
    printf ( "TYPE - Master process:\n" );
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
