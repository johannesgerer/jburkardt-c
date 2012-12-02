# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>

int main ( int argc, char *argv[] );
double f ( double x );
double cpu_time ( void );
void timestamp ( void );

/******************************************************************************/

int main ( int argc, char *argv[] )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for QUAD_SERIAL.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    25 October 2011

  Author:

    John Burkardt
*/
{
  double a = 0.0;
  double b = 10.0;
  double error;
  double exact = 0.49936338107645674464;
  int i;
  int n = 10000000;
  double total;
  double wtime;
  double wtime1;
  double wtime2;
  double x;

  timestamp ( );
  printf ( "\n" );
  printf ( "QUAD_SERIAL:\n" );
  printf ( "  C version\n" );
  printf ( "  Estimate the integral of f(x) from A to B.\n" );
  printf ( "  f(x) = 50 / ( pi * ( 2500 * x * x + 1 ) ).\n" );
  printf ( "\n" );
  printf ( "  A        = %f\n", a );
  printf ( "  B        = %f\n", b );
  printf ( "  N        = %d\n", n );
  printf ( "  Exact    = %24.16f\n", exact );

  wtime1 = cpu_time ( );

  total = 0.0;
  for ( i = 0; i < n; i++ )
  {
    x = ( ( n - i - 1 ) * a + ( i ) * b ) / ( n - 1 );
    total = total + f ( x );
  }

  wtime2 = cpu_time ( );

  total = ( b - a ) * total / ( double ) n;
  error = fabs ( total - exact );
  wtime = wtime2 - wtime1;

  printf ( "\n" );
  printf ( "  Estimate = %24.16f\n", total );
  printf ( "  Error    = %e\n", error );
  printf ( "  Time     = %f\n", wtime );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "QUAD_SERIAL:\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/*******************************************************************************/

double f ( double x )

/*******************************************************************************/
/*
  Purpose:
 
    F evaluates the function.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    18 July 2010

  Author:

    John Burkardt

  Parameters:

    Input, double X, the argument.

    Output, double F, the value of the function.
*/
{
  double pi = 3.141592653589793;
  double value;

  value = 50.0 / ( pi * ( 2500.0 * x * x + 1.0 ) );

  return value;
}
/*******************************************************************************/

double cpu_time ( void )

/*******************************************************************************/
/*
  Purpose:
 
    CPU_TIME reports the total CPU time for a program.

  Modified:

    27 September 2005

  Author:

    John Burkardt

  Parameters:

    Output, double CPU_TIME, the current total elapsed CPU time in second.
*/
{
  double value;

  value = ( double ) clock ( ) / ( double ) CLOCKS_PER_SEC;

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
