# include <stdlib.h>
# include <stdbool.h>
# include <stdio.h>
# include <complex.h>
# include <time.h>

int main ( int argc, char *argv[] );
void test01 ( void );
void test02 ( void );
void test03 ( void );
void timestamp ( void );

/******************************************************************************/

int main ( int argc, char *argv[] )

/******************************************************************************/
/*
  Purpose:

    SIZES returns the size of various data types.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    26 June 2012

  Author:

    John Burkardt
*/
{
  timestamp ( );

  printf ( "\n" );
  printf ( "SIZES\n" );
  printf ( "  C version\n" );
  printf ( "  Return the size of various data types.\n" );
  printf ( "  Compiled on %s at %s\n", __DATE__, __TIME__ );

  test01 ( );
  test02 ( );
  test03 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "SIZES\n" );
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

    TEST01 shows that you can ask for the size of a datatype by name.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    26 June 2012

  Author:

    John Burkardt
*/
{
  printf ( "\n" );
  printf ( "TEST01:\n" );
  printf ( "  You can ask for the size of a datatype by name:\n" );
  printf ( "\n" );
  printf ( "  sizeof (            int ) = %i\n", sizeof ( int ) );
  printf ( "  sizeof (      short int ) = %i\n", sizeof ( short int ) );
  printf ( "  sizeof (       long int ) = %i\n", sizeof ( long int ) );
  printf ( "  sizeof (  long long int ) = %i\n", sizeof ( long long int ) );

  printf ( "  sizeof (          float ) = %i\n", sizeof ( float ) );
  printf ( "  sizeof (         double ) = %i\n", sizeof ( double ) );
  printf ( "  sizeof (    long double ) = %i\n", sizeof ( long double ) );

  printf ( "  sizeof (        complex ) = %i\n", sizeof ( complex ) );
  printf ( "  sizeof ( double complex ) = %i\n", sizeof ( double complex ) );
  printf ( "  sizeof (           char ) = %i\n", sizeof ( char ) );
  printf ( "  sizeof (           bool ) = %i\n", sizeof ( bool ) );

  return;
}
/******************************************************************************/

void test02 ( )

/******************************************************************************/
/*
  Purpose:

    TEST02 shows that you can ask for the size of a specific variable.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    26 June 2012

  Author:

    John Burkardt
*/
{
  int a;
  short int b;
  long int c;
  long long int d;
  
  float e;
  double f;
  long double g;

  complex h;
  double complex i;

  char j;
  bool k;

  printf ( "\n" );
  printf ( "TEST02:\n" );
  printf ( "  You can ask for the size of a specific variable:\n" );
  printf ( "\n" );
  printf ( "  sizeof ( a ) = %i\n", sizeof ( a ) );
  printf ( "  sizeof ( b ) = %i\n", sizeof ( b ) );
  printf ( "  sizeof ( c ) = %i\n", sizeof ( c ) );
  printf ( "  sizeof ( d ) = %i\n", sizeof ( d ) );

  printf ( "  sizeof ( e ) = %i\n", sizeof ( e ) );
  printf ( "  sizeof ( f ) = %i\n", sizeof ( f ) );
  printf ( "  sizeof ( g ) = %i\n", sizeof ( g ) );

  printf ( "  sizeof ( h ) = %i\n", sizeof ( h ) );
  printf ( "  sizeof ( i ) = %i\n", sizeof ( i ) );
  printf ( "  sizeof ( j ) = %i\n", sizeof ( j ) );
  printf ( "  sizeof ( k ) = %i\n", sizeof ( k ) );

  return;
}
/******************************************************************************/

void test03 ( )

/******************************************************************************/
/*
  Purpose:

    TEST03 shows that you can ask for the size of an array.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    26 June 2012

  Author:

    John Burkardt
*/
{
  int l[10];
  float m[10];
  double n[10];
  char o[10];

  int p[5][2];
  float q[5][2];
  double r[5][2];
  char s[5][2];

  printf ( "\n" );
  printf ( "TEST03:\n" );
  printf ( "  You can ask for the size of an array:\n" );
  printf ( "\n" );
  printf ( "  sizeof ( l ) = %i\n", sizeof ( l ) );
  printf ( "  sizeof ( m ) = %i\n", sizeof ( m ) );
  printf ( "  sizeof ( n ) = %i\n", sizeof ( n ) );
  printf ( "  sizeof ( o ) = %i\n", sizeof ( o ) );

  printf ( "  sizeof ( p ) = %i\n", sizeof ( p ) );
  printf ( "  sizeof ( q ) = %i\n", sizeof ( q ) );
  printf ( "  sizeof ( r ) = %i\n", sizeof ( r ) );
  printf ( "  sizeof ( s ) = %i\n", sizeof ( s ) );

  return;
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
