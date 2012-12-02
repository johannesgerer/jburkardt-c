# include <stdlib.h>
# include <stdio.h>
# include <time.h>

int main ( int argc, char *argv[] );
void test01 ( );
void test02 ( );
void timestamp ( void );

/********************************************************************/

int main ( int argc, char *argv[] )

/********************************************************************/
/*
  Purpose:

    RAND_TEST generates random numbers.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    24 September 2006

  Author:

    John Burkardt
*/
{
  timestamp ( );

  printf ( "\n" );
  printf ( "RAND_TEST:\n" );
  printf ( "  C version\n" );
  printf ( "  Generate random numbers using\n" );
  printf ( "  SRAND to set the seed, and\n" );
  printf ( "  RAND to return the random values.\n" );
  printf ( "\n" );
  printf ( "  Compiled on %s at %s\n", __DATE__, __TIME__ );

  test01 ( );
  test02 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "RAND_TEST:\n" );
  printf ( "  Normal end of execution.\n" );

  printf ( "\n" );
  timestamp ( );

  return 0;
}
/********************************************************************/

void test01 ( void )

/********************************************************************/
/*
  Purpose:

    TEST01 simply calls RAND a few times.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    24 September 2006

  Author:

    John Burkardt
*/
{
  int i;
  unsigned int seed = 123456789;
  long int x;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  Call SRAND to initialize the seed.\n" );
  printf ( "  Call RAND to generate some values.\n" );
  printf ( "\n" );
  printf ( "  Initial SEED = %d\n", seed );

  srand ( seed );

  printf ( "\n" );
  printf ( "      Step    RANDOM()\n" );
  printf ( "\n" );

  for ( i = 0; i < 10; i++ )
  {
    x = rand ( );
    printf ( "  %8d  %12d\n", i, x );
  }
}
/********************************************************************/

void test02 ( void )

/********************************************************************/
/*
  Purpose:

    TEST02 uses RAND to generate real values in [0,1];

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    19 September 2006

  Author:

    John Burkardt
*/
{
  int i;
  unsigned int seed = 123456789;
  long int x;
  double y;

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  Call SRAND to initialize the seed.\n" ); 
  printf ( "  Call RAND to generate some values.\n" );
  printf ( "  Set X = ( double ) RAND ( ) / ( double ) RAND_MAX\n" );
  printf ( "  so that X is a random real in [0,1].\n" );

  printf ( "\n" );
  printf ( "  RAND_MAX = %d\n", RAND_MAX );

  printf ( "\n" );
  printf ( "  Initial SEED = %d\n", seed );
  srand ( seed );

  printf ( "\n" );
  printf ( "      Step      RAND()        RAND()/RAND_MAX\n" );
  printf ( "\n" );

  for ( i = 0; i < 10; i++ )
  {
    x = rand ( );

    y = ( double ) x / ( double ) RAND_MAX;

    printf ( "  %8d  %12d  %12f\n", i, x, y );
  }
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
