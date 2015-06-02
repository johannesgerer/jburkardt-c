# include <stdlib.h>
# include <stdio.h>
# include <time.h>

int main ( int argc, char *argv[] );
void timestamp ( void );

/******************************************************************************/

int main ( int argc, char *argv[] )

/******************************************************************************/
/*
 Purpose:

    NOT_ALLOCATED_ARRAYS shows why it's good to initialize array pointers.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    02 March 2006

  Author:

    John Burkardt
*/
{
  int *a;
  int *b = NULL;

  timestamp ( );

  printf ( "\n" );
  printf ( "NOT_ALLOCATED_ARRAYS\n" );
  printf ( "  C version\n" );
  printf ( "\n" );
  printf ( "  When an array starts out as a pointer, you have to use the\n" );
  printf ( "  MALLOC command to allocate memory.  You should always \n" );
  printf ( "  initialize such array pointers to NULL, so you can tell if\n" );
  printf ( "  they've been allocated or not!\n" );
  printf ( "\n" );
  printf ( "  Unfortunately, when you FREE an array, you also have to\n" );
  printf ( "  reset the pointer to NULL; that does not happen\n" );
  printf ( "  automatically either!\n" );

  printf ( "\n" );
  printf ( "  The pointer A is not preset to NULL.\n" );
  printf ( "  Before allocation, we check the value:\n" );
  printf ( "    a = %d\n", a );
  printf ( "  The test 'if ( !a )' is not guaranteed to return 1\n" );
  printf ( "  because we did not initialize A properly.\n" );
  printf ( "    !a = %d\n", !a );
  printf ( "  Now we allocate A.\n" );

  a = ( int * ) malloc ( 10 * sizeof ( int ) );

  printf ( "    a = %d\n", a );
  printf ( "    !a = %d\n", !a );
  printf ( "  Now we FREE A.\n" );
  free ( a );
  printf ( "    a = %d\n", a );
  printf ( "    !a = %d\n", !a );
  printf ( "  Now we RESET A to NULL!\n" );
  a = NULL;
  printf ( "    a = %d\n", a );
  printf ( "    !a = %d\n", !a );

  printf ( "\n" );
  printf ( "  The pointer B is preset to NULL.\n" );
  printf ( "  Before allocation, we check the value:\n" );
  printf ( "    b = %d\n", b );
  printf ( "  The test 'if ( !b )' is guaranteed to return 1\n" );
  printf ( "  because we initialized B properly.\n" );
  printf ( "    !b = %d\n", !b );
  printf ( "  Now we allocate B.\n" );

  b = ( int * ) malloc ( 10 * sizeof ( int ) );

  printf ( "    b = %d\n", b );
  printf ( "    !b = %d\n", !b );
  printf ( "  Now we FREE B.\n" );
  free ( b );
  printf ( "    b = %d\n", b );
  printf ( "    !b = %d\n", !b );
  printf ( "  Now we RESET B to NULL!\n" );
  b = NULL;
  printf ( "    b = %d\n", b );
  printf ( "    !b = %d\n", !b );

  printf ( "\n" );
  printf ( "NOT_ALLOCATED_ARRAYS:\n" );
  printf ( "  Normal end of execution.\n" );

  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void timestamp ( void )

/******************************************************************************/
/*
  Purpose:

    TIMESTAMP prints the current YMDHMS date as a time stamp.

  Example:

    May 31 2001 09:45:54 AM

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    03 October 2003

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
