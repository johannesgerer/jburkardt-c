# include <stdlib.h>
# include <stdio.h>
# include <time.h>

int main ( int argc, char *argv[] );
void timestamp ( void );

/**********************************************************************/

int main ( int argc, char *argv[] )

/**********************************************************************/
/*
  Purpose:

    Demonstrate the use of arrays.

  Discussion:

    Here is an example of how to set up a dynamic two dimensional array, 
    essentially by declaring it to be a one dimensional array, and doing 
    the double indexing yourself.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    07 October 2012

  Author:

    John Burkardt
*/
{
  int a[10];
  int *b;
  int c[5][5];
  int *d;
  int *e;
  int i;
  int j;
  int k;
  
  timestamp ( );
  printf ( "\n" );
  printf ( "ARRAYS\n" );
  printf ( "  C version\n" );
  printf ( "\n" );
  printf ( "  Examples of array use.\n" );
  printf ( "  You can create and manipulate multiple-dimension arrays\n" );
  printf ( "  as though they were vectors, as long as you're willing.\n" );
  printf ( "  to do the indexing yourself.\n" );

  for ( i = 0; i < 10; i++ )
  {
    a[i] = 500 + i;
  }
  
  printf ( "\n" );
  printf ( "  The 1D array A[]:\n" );
  printf ( "\n" );
  printf ( "       i    A[i]\n" );
  printf ( "\n" );

  for ( i = 0; i < 10; i++ )
  {
    printf ( "  %6d  %6d\n", i, a[i] );
  }

  b = malloc ( 10 * sizeof ( int ) );

  for ( i = 0; i < 10; i++ )
  {
    b[i] = 500 + i;
  }
  
  printf ( "\n" );
  printf ( "  The 1D vector B[], allocated dynamically:\n" );
  printf ( "\n" );
  printf ( "       i    B[i]\n" );
  printf ( "\n" );

  for ( i = 0; i < 10; i++ )
  {
    printf ( "  %6d  %6d\n", i, b[i] );
  }

  free ( b );

  for ( j = 0; j < 5; j++ )
  {
    for ( i = 0; i < 5; i++ )
    {
      c[i][j] = 10 * ( i + 1 ) + ( j + 1 );
    }
  }
  
  printf ( "\n" );
  printf ( "  The matrix C[][]:\n" );
  printf ( "\n" );
  printf ( "       k       i       j    C[i][j]\n" );
  printf ( "\n" );

  k = 0;
  for ( j = 0; j < 5; j++ )
  {
    for ( i = 0; i < 5; i++ )
    {
      printf ( "  %6d  %6d  %6d  %6d\n", k, i, j, c[i][j] );
      k = k + 1;
    }
  }

  d = malloc ( 5*7*sizeof ( int ) );

  for ( j = 0; j < 7; j++ )
  {
    for ( i = 0; i < 5; i++ )
    {
      d[i+j*5] = 10 * ( i + 1 ) + ( j + 1 );
    }
  }
  
  printf ( "\n" );
  printf ( "  The matrix D, stored as a column-major vector:\n" );
  printf ( "\n" );
  printf ( "       k       i       j    D[i+j*5]\n" );
  printf ( "\n" );

  k = 0;
  for ( j = 0; j < 7; j++ )
  {
    for ( i = 0; i < 5; i++ )
    {
      printf ( "  %6d  %6d  %6d  %6d\n", k, i, j, d[i+j*5] );
     k = k + 1;
    }
  }

  free ( d );
/*
  Same thing, but now row major.
*/
  e = malloc ( 5*7*sizeof ( int ) );

  for ( j = 0; j < 7; j++ )
  {
    for ( i = 0; i < 5; i++ )
    {
      e[i*7+j] = 10 * ( i + 1 ) + ( j + 1 );
    }
  }
  
  printf ( "\n" );
  printf ( "  The matrix E, stored as a row-major vector:\n" );
  printf ( "\n" );
  printf ( "       k       i       j    E[i*7+j]\n" );
  printf ( "\n" );

  k = 0;
  for ( i = 0; i < 5; i++ )
  {
    for ( j = 0; j < 7; j++ )
    {
      printf ( "  %6d  %6d  %6d  %6d\n", k, i, j, e[i*7+j] );
      k = k + 1;
    }
  }

  free ( e );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "ARRAYS:\n" );
  printf ( "  Normal end of execution.\n" );

  printf ( "\n" );
  timestamp ( );

  return 0;
}
/********************************************************************/

void timestamp ( void )

/********************************************************************/
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
