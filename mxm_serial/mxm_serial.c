# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>

int main ( void );
void timestamp ( void );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for MXM_SERIAL.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    30 October 2011

  Author:

    John Burkardt
*/
{
  double a[500][500];
  double angle;
  double b[500][500];
  double c[500][500];
  int i;
  int j;
  int k;
  int n = 500;
  double pi = 3.141592653589793;
  double s;

  timestamp ( );

  printf ( "\n" );
  printf ( "MXM_SERIAL:\n" );
  printf ( "  C version\n" );
  printf ( "  Compute matrix product C = A * B.\n" );
  printf ( "\n" );
  printf ( "  The matrix order N                 = %d\n", n );
/*
  Loop 1: Evaluate A.
*/
  s = 1.0 / sqrt ( ( double ) ( n ) );

  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      angle = 2.0 * pi * i * j / ( double ) n;
      a[i][j] = s * ( sin ( angle ) + cos ( angle ) );
    }
  }
/*
  Loop 2: Copy A into B.
*/
  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      b[i][j] = a[i][j];
    }
  }
/*
  Loop 3: Compute C = A * B.
*/
  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      c[i][j] = 0.0;
      for ( k = 0; k < n; k++ )
      {
        c[i][j] = c[i][j] + a[i][k] * b[k][j];
      }
    }
  }

  printf ( "  C(100,100)  = %g\n", c[99][99] );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "MXM_SERIAL:\n" );
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
