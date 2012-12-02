# include <stdlib.h>
# include <stdio.h>
# include <time.h>

# include "mpi.h"

int main ( int argc, char *argv[] );
int search ( int a, int b, int c, int id, int p );
int f ( int i );
void timestamp ( void );

/******************************************************************************/

int main ( int argc, char *argv[] )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for SEARCH_MPI.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 November 2012

  Author:

    John Burkardt
*/
{
  int a;
  int b;
  int c;
  int fj;
  int i4_huge = 2147483647;
  int id;
  int j;
  int p;
  double wtime;
/*
  Initialize MPI.
*/
  MPI_Init ( &argc, &argv );
/*
  Get this processor's ID.
*/
  MPI_Comm_rank ( MPI_COMM_WORLD, &id );
/*
  Get the number of processes.
*/
  MPI_Comm_size ( MPI_COMM_WORLD, &p );

  a = 1;
  b = i4_huge;
  c = 45;

  if ( id == 0 )
  {
    timestamp ( );
    printf ( "\n" );
    printf ( "SEARCH_MPI:\n" );
    printf ( "  C/MPI version\n" );
    printf ( "  Search the integers from A to B\n" );
    printf ( "  for a value J such that F(J) = C.\n" );
    printf ( "\n" );
    printf ( "  A           = %d\n", a );
    printf ( "  B           = %d\n", b );
    printf ( "  C           = %d\n", c );
  }

  wtime = MPI_Wtime ( );

  j = search ( a, b, c, id, p );

  wtime = MPI_Wtime ( ) - wtime;

  if ( j != -1 )
  {
    printf ( "\n" );
    printf ( "  Process %d found     J = %d\n", id, j );
    printf ( "  Verify F(J) = %d\n", f ( j ) );
  }

  if ( id == 0 )
  {
    printf ( "  Elapsed wallclock time is %g\n", wtime );
  }
/*
  Terminate.
*/
  if ( id == 0 )
  {
    printf ( "\n" );
    printf ( "SEARCH_MPI:\n" );
    printf ( "  Normal end of execution.\n" );
    printf ( "\n" );
    timestamp ( );
  }

  MPI_Finalize ( );

  return 0;
}
/******************************************************************************/

int search ( int a, int b, int c, int id, int p )

/******************************************************************************/
/*
  Purpose:

    SEARCH searches integers in [A,B] for a J so that F(J) = C.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 November 2012

  Author:

    John Burkardt

  Parameters:

    Input, int A, B, the search range.

    Input, int C, the desired function value.

    Input, int ID, the process ID.

    Input, int P, the number of processes.

    Output, int SEARCH, the computed solution, or -1
    if no solution was found.
*/
{
  int fi;
  int i;
  int j;

  j = -1;
/*
  i = i + p can take us "over top" so that i becomes negative!
  So we have to be more careful here!
*/
  for ( i = a + id; 0 < i && i <= b; i = i + p )
  {
    fi = f ( i );

    if ( fi == c )
    {
      j = i;
      break;
    }
  }

  return j;
}
/******************************************************************************/

int f ( int i )

/******************************************************************************/
/*
  Purpose:

    F is the function we are analyzing.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    22 October 2012

  Author:

    John Burkardt

  Parameters:

    Input, int I, the argument.

    Input, int F, the value.
*/
{
  int i4_huge = 2147483647;
  int j;
  int k;
  int value;

  value = i;

  for ( j = 1; j <= 5; j++ )
  {
    k = value / 127773;

    value = 16807 * ( value - k * 127773 ) - k * 2836;

    if ( value <= 0 )
    {
      value = value + i4_huge;
    }
  }

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
/******************************************************************************/

double cpu_time ( void )

/******************************************************************************/
/*
  Purpose:

    CPU_TIME returns the current reading on the CPU clock.

  Discussion:

    The CPU time measurements available through this routine are often
    not very accurate.  In some cases, the accuracy is no better than
    a hundredth of a second.  

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 June 2005

  Author:

    John Burkardt

  Parameters:

    Output, double CPU_TIME, the current reading of the CPU clock, in seconds.
*/
{
  double value;

  value = ( double ) clock ( ) 
        / ( double ) CLOCKS_PER_SEC;

  return value;
}
