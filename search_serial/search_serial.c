# include <stdlib.h>
# include <stdio.h>
# include <time.h>

int main ( );
int search ( int a, int b, int c );
int f ( int i );
void timestamp ( void );
double cpu_time ( void );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for SEARCH_SERIAL.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    22 October 2012

  Author:

    John Burkardt
*/
{
  int a;
  int b;
  int c;
  int fj;
  int i4_huge = 2147483647;
  int j;
  double wtime;

  a = 1;
  b = i4_huge;
  c = 45;

  timestamp ( );
  printf ( "\n" );
  printf ( "SEARCH_SERIAL:\n" );
  printf ( "  C version\n" );
  printf ( "  Search the integers from A to B\n" );
  printf ( "  for a value J such that F(J) = C.\n" );
  printf ( "\n" );
  printf ( "  A           = %d\n", a );
  printf ( "  B           = %d\n", b );
  printf ( "  C           = %d\n", c );

  wtime = cpu_time ( );

  j = search ( a, b, c );

  wtime = cpu_time ( ) - wtime;

  if ( j == -1 )
  {
    printf ( "\n" );
    printf ( "  No solution was found.\n" );
  }
  else
  {
    printf ( "\n" );
    printf ( "  Found     J = %d\n", j );
    printf ( "  Verify F(J) = %d\n", f ( j ) );
  }

  printf ( "  Elapsed CPU time is %g\n", wtime );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "SEARCH_SERIAL:\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

int search ( int a, int b, int c )

/******************************************************************************/
/*
  Purpose:

    SEARCH searches integers in [A,B] for a J so that F(J) = C.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    22 October 2012

  Author:

    John Burkardt

  Parameters:

    Input, int A, B, the search range.

    Input, int C, the desired function value.

    Output, int SEARCH, the computed solution, or -1
    if no solution was found.
*/
{
  int fi;
  int i;
  int j;

  j = -1;

  for ( i = a; i <= b; i++ )
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
