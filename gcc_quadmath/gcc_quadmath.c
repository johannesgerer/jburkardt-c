# include <stdlib.h>
# include <stdio.h>
# include <math.h>
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

    MAIN is the main program for GCC_QUADMATH.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    18 April 2011

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "GCC_QUADMATH:\n" );
  printf ( "  C version.\n" );
  printf ( "  Test the GCC quadmath library.\n" );

  test01 ( );
  test02 ( );
  test03 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "GCC_QUADMATH:\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void test01 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST01 uses flaot arithmetic.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    18 April 2011

  Author:

    John Burkardt
*/
{
  int divs;
  float x;
  float x_old;
  float y;
  float z;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  Using FLOAT arithmetic:\n" );
  printf ( "  Compute smallest 1/2^DIV that can be added to 1.\n" );

  x = 1.0;
  z = 1.0;
  divs = 0;

  for ( ; ; )
  {
    x_old = x;
    x = x / 2.0;
    y = 1.0 + x;
    if ( y <= z )
    {
      break;
    }
    divs = divs + 1;
  }

  printf ( "  Number of divisions DIV = %d\n", divs );
  printf ( "  1/2^DIV =         %e\n", x_old );

  return;
}
/******************************************************************************/

void test02 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST02 uses double arithmetic.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    18 April 2011

  Author:

    John Burkardt
*/
{
  int divs;
  double x;
  double x_old;
  double y;
  double z;

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  Using DOUBLE arithmetic:\n" );
  printf ( "  Compute smallest 1/2^DIV that can be added to 1.\n" );

  x = 1.0;
  z = 1.0;
  divs = 0;

  for ( ; ; )
  {
    x_old = x;
    x = x / 2.0;
    y = 1.0 + x;
    if ( y <= z )
    {
      break;
    }
    divs = divs + 1;
  }

  printf ( "  Number of divisions DIV = %d\n", divs );
  printf ( "  1/2^DIV =         %e\n", x_old );

  return;
}
/******************************************************************************/

void test03 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST03 uses __float80 arithmetic.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    18 April 2011

  Author:

    John Burkardt
*/
{
  int divs;
  __float80 x;
  __float80 x_old;
  __float80 y;
  __float80 z;

  printf ( "\n" );
  printf ( "TEST03\n" );
  printf ( "  Using __FLOAT80 arithmetic:\n" );
  printf ( "  Compute smallest 1/2^DIV that can be added to 1.\n" );

  x = 1.0;
  z = 1.0;
  divs = 0;

  for ( ; ; )
  {
    x_old = x;
    x = x / 2.0;
    y = 1.0 + x;
    if ( y <= z )
    {
      break;
    }
    divs = divs + 1;
  }

  printf ( "  Number of divisions DIV = %d\n", divs );
  printf ( "  1/2^DIV =         %Le\n", x_old );

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
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  printf ( "%s\n", time_buffer );

  return;
# undef TIME_SIZE
}

