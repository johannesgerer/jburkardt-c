# include <stdlib.h>
# include <stdio.h>
# include <time.h>

int main ( void );
double cpu_time ( void );
double *set_up ( int n );
void sum_up ( int n, double x[], double *total, double *ctime );
void timestamp ( void );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for SUM_MILLION.

  Discussion:

    This code estimates the power of a computer by summing the integers
    from 1 to 1,000,000.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    24 July 2008

  Author:

    John Burkardt
*/
{
  double ctime;
  double error;
  double exact = 500000500000.0;
  int i;
  double mflops;
  int n = 1000000;
  double total;
  double *x;

  printf ( "\n" );
  timestamp ( );
  printf ( "\n" );
  printf ( "SUM_MILLION\n" );
  printf ( "  C version\n" );
  printf ( "  Sum the integers from 1 to 1,000,000.\n" );
  printf ( "  Correct answer is 500000500000.\n" );

  x = set_up ( n );

  printf ( "\n" );
  printf ( "         N      CPU time        MFLOPS          ERROR\n" );
  printf ( "                (seconds)\n" );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    sum_up ( n, x, &total, &ctime );

    mflops = ( double ) ( n ) / 1000000.0 / ctime;

    error = total - exact;

    printf ( "  %8d  %14f  %14f  %14e\n", n, ctime, mflops, error );
  }

  free ( x );

  printf ( "\n" );
  printf ( "SUM_MILLION:\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
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
/******************************************************************************/

double *set_up ( int n )

/******************************************************************************/
/*
  Purpose:

    SET_UP sets up the data for the SUM_MILLION program.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 July 2008

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of values to define.

    Output, double X[N], a vector which contains the values 1 through N.
*/
{
  int i;
  double *x;

  x = ( double * ) malloc ( n * sizeof ( double ) );

  for ( i = 0; i < n; i++ )
  {
    x[i] = ( double ) ( i + 1 );
  }
  return x;
}
/******************************************************************************/

void sum_up ( int n, double x[], double *total, double *ctime )

/******************************************************************************/
/*
  Purpose:

    SUM_UP carries out the sum for the SUM_MILLION program.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 July 2008

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of values to define.

    Input, double X[N], the data to be summed.

    Output, double TOTAL, the sum of the data.

    Output, double CTIME, the cpu time required to sum the data.
*/
{
  double ctime1;
  double ctime2;
  int i;

  ctime1 = cpu_time ( );

  *total = 0.0;
  for ( i = 0; i < n; i++ )
  {
    *total = *total + x[i];
  }

  ctime2 = cpu_time ( );

  *ctime = ctime2 - ctime1;

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

  printf ( "%s\n", time_buffer );

  return;
# undef TIME_SIZE
}
