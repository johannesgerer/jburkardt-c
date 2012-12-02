# include <stdlib.h>
# include <stdio.h>
# include <time.h>

int main ( int argc, char *argv[] );
void test01 ( );
void timestamp ( );

/******************************************************************************/

int main ( int argc, char *argv[] )

/******************************************************************************/
/*
  Purpose:

    DRAND48_TEST generates random numbers using DRAND48().

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 September 2012

  Author:

    John Burkardt
*/
{
  timestamp ( );

  printf ( "\n" );
  printf ( "DRAND48_TEST:\n" );
  printf ( "  C version\n" );
  printf ( "  Generate random numbers using\n" );
  printf ( "  SRAND48 to set the seed, and\n" );
  printf ( "  DRAND48 to return the random values.\n" );
  printf ( "\n" );
  printf ( "  Compiled on %s at %s\n", __DATE__, __TIME__ );

  test01 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "DRAND48_TEST:\n" );
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

    TEST01 simply calls DRAND48 a few times.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 September 2012

  Author:

    John Burkardt
*/
{
  int i;
  long long int seed = 123456789LL;
  double x;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  Call SRAND48 to initialize the seed.\n" );
  printf ( "  Call DRAND48 to generate some values.\n" );
  printf ( "\n" );
  printf ( "  Initial SEED = %d\n", seed );

  srand48 ( seed );
  printf ( "\n" );
  printf ( "      Step    DRAND48()\n" );
  printf ( "\n" );

  for ( i = 0; i < 10; i++ )
  {
    x = drand48 ( );
    printf ( "  %8d  %12g\n", i, x );
  }
  return;
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
  size_t len;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  len = strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  fprintf ( stdout, "%s\n", time_buffer );

  return;
# undef TIME_SIZE
}

