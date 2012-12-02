# include <stdlib.h>
# include <stdio.h>
# include <time.h>

int main ( int argc, char *argv[] );
double cpu_time ( void );
double *matgen ( int m, int n );
double mxv_foriforj ( int m, int n, double a[], double x[], double y[] );
double mxv_forjfori ( int m, int n, double a[], double x[], double y[] );
void timestamp ( void );

/******************************************************************************/

int main ( int argc, char *argv[] )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for MXV.

  Discussion:

    MXV computes a matrix-vector product in a number of ways, and reports
    the elapsed CPU time.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    25 May 2008

  Author:

    John Burkardt

  Usage:

    mxv m n

  Parameters:

    Command line argument, int M, the number of rows in the matrix.

    Command line argument, int N, the number of columns in the matrix.
*/
{
  double *a;
  double cpu_seconds;
  int flop_count;
  int i;
  int m;
  double mflops;
  int n;
  double *x;
  double *y;

  timestamp ( );

  printf ( "\n" );
  printf ( "MXV:\n" );
  printf ( "  C version\n" );
  printf ( "  Compiled on %s at %s.\n", __DATE__, __TIME__ );
  printf ( "\n" );
  printf ( "  Compute matrix vector products y = A*x.\n" );
/*
  Get the number of rows, M.
*/
  if ( 2 <= argc )
  {
    m = atoi ( argv[1] );
  } 
  else
  {
    printf ( "\n" );
    printf ( "  Enter the number of rows, M\n" );
    scanf ( "%d", &m );
  }
/*
  Get the number of columns, N.
*/
  if ( 3 <= argc )
  {
    n = atoi ( argv[2] );
  }
  else
  {
    printf ( "\n" );
    printf ( "  Enter the number of columns, N\n" );
    scanf ( "%d", &n );
  }
/*
  Record the amount of work.
  Each of the M entries of Y requires N multiplies and N adds.
*/
  flop_count = 2 * m * n;

  printf ( "\n" );
  printf ( "  Number of matrix rows M =             %12d\n", m );
  printf ( "  Number of matrix columns N =          %12d\n", n );
  printf ( "  Number of floating point operations = %12d\n", flop_count );
/*
  Set A and X.
*/
  a = matgen ( m, n );

  x = ( double * ) malloc ( n * sizeof ( double ) );

  for ( i = 0; i < n; i++ )
  {
    x[i] = ( double ) ( i + 1 );
  }

  y = ( double * ) malloc ( m * sizeof ( double ) );

  printf ( "\n" );
  printf ( "  Method     Cpu Seconds       MegaFlopS\n" );
  printf ( "  ------  --------------  --------------\n" );
/*
  FORIFORJ
*/
  cpu_seconds = mxv_foriforj ( m, n, a, x, y );

  if ( 0.0 < cpu_seconds )
  {
    mflops = ( double ) ( flop_count ) / cpu_seconds / 1000000.0;
  }
  else
  {
    mflops = - 1.0;
  }

  printf ( "  FORIFORJ  %14f  %14f\n", cpu_seconds, mflops );
/*
  FORJFORI
*/
  cpu_seconds = mxv_forjfori ( m, n, a, x, y );

  if ( 0.0 < cpu_seconds )
  {
    mflops = ( double ) ( flop_count ) / cpu_seconds / 1000000.0;
  }
  else
  {
    mflops = - 1.0;
  }

  printf ( "  FORJFORI  %14f  %14f\n", cpu_seconds, mflops );
/*
  Deallocate arrays.
*/
  free ( a );
  free ( x );
  free ( y );

  printf ( "\n" );
  printf ( "MXV:\n" );
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

double *matgen ( int m, int n )

/******************************************************************************/
/*
  Purpose:

    MATGEN generates a random matrix A and vector X.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    23 May 2008

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns of the matrix.

    Output, double MATGEN[M*N], the matrix.
*/
{
  double *a;
  int i;
  int j;
  int k;
  int seed;

  a = ( double * ) malloc ( m * n * sizeof ( double ) );

  seed = 1325;
/*
 Set the matrix A.
*/
  k = 0;
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      seed = ( ( 3125 * seed ) % 65536 );
      a[k] = ( seed - 32768.0 ) / 16384.0;
      k = k + 1;
    }
  }
  return a;
}
/******************************************************************************/

double mxv_foriforj ( int m, int n, double a[], double x[], double y[] )

/******************************************************************************/
/*
  Purpose:

    MXV_FORIFORJ computes y = A * x, using FOR I, FOR J loops.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    23 May 2008

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns
    of the matrix.

    Input, double A[M*N], the matrix.

    Input, double X[N], the vector to be multiplied.

    Output, double Y[M], the product vector.

    Output, double MXV_FORIFORJ, the elapsed CPU time.
*/
{
  double cpu_seconds;
  int i;
  int j;
  double time1;
  double time2;

  time1 = cpu_time ( );

  for ( i = 0; i < m; i++ )
  {
    y[i] = 0.0;
    for ( j = 0; j < n; j++ )
    {
      y[i] = y[i] + a[i+j*m] * x[j];
    }
  }

  time2 = cpu_time ( );

  cpu_seconds = time2 - time1;

  return cpu_seconds;
}
/******************************************************************************/

double mxv_forjfori ( int m, int n, double a[], double x[], double y[] )

/******************************************************************************/
/*
  Purpose:

    MXV_FORJFORI computes y = A * x, using FOR J, FOR I loops.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 May 2008

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns
    of the matrix.

    Input, double A[M*N], the matrix.

    Input, double X[N], the vector to be multiplied.

    Output, double Y[M], the product vector.

    Output, double MXV_FORJFORI, the elapsed CPU time.
*/
{
  double cpu_seconds;
  int i;
  int j;;
  double time1;
  double time2;

  time1 = cpu_time ( );

  for ( i = 0; i < m; i++ )
  {
    y[i] = 0.0;
  }

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      y[i] = y[i] + a[i+j*m] * x[j];
    }
  }

  time2 = cpu_time ( );

  cpu_seconds = time2 - time1;

  return cpu_seconds;
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
