# include <stdlib.h>
# include <stdio.h>
# include <float.h>
# include <limits.h>
# include <time.h>

int main ( );
void test01 ( );
void timestamp ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for BIG_INTS_REAL.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    28 September 2014

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "BIG_INTS_REAL:\n" );
  printf ( "  C version\n" );
  printf ( "  Examine the transfer of integer values into and out of real variables.\n" );

  test01 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "BIG_INTS_REAL:\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void test01 ( )

/******************************************************************************/
/*
  Purpose:

    TEST01 stores huge integers as reals.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    28 September 2014

  Author:

    John Burkardt
*/
{
  int i4;
  int i4r4i4;
  int i4r8i4;
  long long int i8;
  long long int i8r4i8;
  long long int i8r8i8;
  float r4;
  float r4i4;
  float r4i8;
  double r8;
  double r8i4;
  double r8i8;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  Compute the largest possible integers.\n" );
  printf ( "  Try to store them as real values.\n" );
  printf ( "  Then copy them back.\n" );

  printf ( "\n" );
  printf ( "  'Huge' integers and huge reals:\n" );
  printf ( "\n" );

  i4 = INT_MAX;
  i8 = LLONG_MAX;
  r4 = FLT_MAX;
  r8 = DBL_MAX;

  printf ( "  i4 = INT_MAX   = %d\n", i4 );
  printf ( "  i8 = LLONG_MAX = %d\n", i8 );
  printf ( "  r4 = FLT_MAX   = %g\n", r4 );
  printf ( "  r8 = DBL_MAX   = %g\n", r8 );

  printf ( "\n" );
  printf ( "  Convert huge integers to real values:\n" );
  printf ( "\n" );

  r4i4 = ( float ) ( i4 );
  r4i8 = ( float ) ( i8 );
  r8i4 = ( double ) ( i4 );
  r8i8 = ( double ) ( i8 );

  printf ( "  r4i4 = ( float ) ( i4 )  = %g\n", r4i4 );
  printf ( "  r4i8 = ( float ) ( i8 )  = %g\n", r4i8 );
  printf ( "  r8i4 = ( double ) ( i4 ) = %g\n", r8i4 );
  printf ( "  r8i8 = ( double ) ( i8 ) = %g\n", r8i8 );

  printf ( "\n" );
  printf ( "  Convert real values of integers back to integers:\n" );
  printf ( "\n" );

  i4r4i4 = ( int ) ( r4i4 );
  i4r8i4 = ( int ) ( r8i4 );
  i8r4i8 = ( long long int ) ( r4i8 );
  i8r8i8 = ( long long int ) ( r8i8 );

  printf ( "  i4r4i4 = ( int ) ( r4i4 )           = %d\n", i4r4i4 );
  printf ( "  i4r8i4 = ( int ) ( r8i4 )           = %d\n", i4r8i4 );
  printf ( "  i8r4i8 = ( long long int ) ( r4i8 ) = %ld\n", i8r4i8 );
  printf ( "  i8r8i8 = ( long long int ) ( r8i8 ) = %ld\n", i8r8i8 );

  return;
}
/******************************************************************************/

void timestamp ( )

/******************************************************************************/
/*
  Purpose:

    TIMESTAMP prints the current YMDHMS date as a time stamp.

  Example:

    17 June 2014 09:45:54 AM

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 June 2014

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
