# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <omp.h>

int main ( void );
int *prime_table ( int prime_num );
double *sine_table ( int sine_num );
void timestamp ( void );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for MULTITASK_OPENMP.

  Discussion:

    This program demonstrates how OpenMP can be used for multitasking, that 
    is, a simple kind of parallel processing in which a certain number of 
    perhaps quite unrelated tasks must be done.

    The OpenMP SECTIONS directive identifies the portion of the program where
    the code for these tasks is given.

    The OpenMP SECTION directive is used repeatedly to divide this area of
    the program into independent tasks.

    The code will get the benefit of parallel processing up to the point where
    there are as many threads as there are tasks.

    The code will get a substantial speedup if the tasks take roughly the
    same amount of time.  However, if one task takes substantially more time
    than the others, this results in a limit to the parallel speedup that is
    possible.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 October 2011

  Author:

    John Burkardt

*/
{
  int prime_num;
  int *primes;
  int sine_num;
  double *sines;
  double wtime;
  double wtime1;
  double wtime2;

  timestamp ( );
  printf ( "\n" );
  printf ( "MULTITASK_OPENMP:\n" );
  printf ( "  C/OpenMP version\n" );
  printf ( "  Demonstrate how OpenMP can \"multitask\" by using the\n" );
  printf ( "  SECTIONS directive to carry out several tasks in parallel.\n" );

  prime_num = 20000;
  sine_num = 20000;

  wtime = omp_get_wtime ( );

# pragma omp parallel shared ( prime_num, primes, sine_num, sines )
{
  # pragma omp sections
  {
    # pragma omp section
    {
      wtime1 = omp_get_wtime ( );
      primes = prime_table ( prime_num );
      wtime1 = omp_get_wtime ( ) - wtime1;
    }
    # pragma omp section
    {
      wtime2 = omp_get_wtime ( );
      sines = sine_table ( sine_num );
      wtime2 = omp_get_wtime ( ) - wtime2;
    }
  }
}
  wtime = omp_get_wtime ( ) - wtime;

  printf ( "\n" );
  printf ( "  Number of primes computed was %d\n", prime_num );
  printf ( "  Last prime was %d\n", primes[prime_num-1] );
  printf ( "  Number of sines computed was %d\n", sine_num );
  printf ( "  Last sine computed was %g\n", sines[sine_num-1] );
  printf ( "\n" );
  printf ( "  Elapsed time = %g\n", wtime );
  printf ( "  Task 1 time = %g\n", wtime1 );
  printf ( "  Task 2 time = %g\n", wtime2 );

  free ( primes );
  free ( sines );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "MULTITASK_OPENMP:\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

int *prime_table ( int prime_num )

/******************************************************************************/
/*
  Purpose:

    PRIME_TABLE computes a table of the first PRIME_NUM prime numbers.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 October 2011

  Author:

    John Burkardt

  Parameters:

    Input, int PRIME_NUM, the number of primes to compute.

    Output, int PRIME_TABLE[PRIME_NUM], the computed primes.
*/
{
  int i;
  int j;
  int p;
  int prime;
  int *primes;

  primes = ( int * ) malloc ( prime_num * sizeof ( int ) );

  i = 2;
  p = 0;

  while ( p < prime_num )
  {
    prime = 1;

    for ( j = 2; j < i; j++ )
    {
      if ( ( i % j ) == 0 )
      {
        prime = 0;
        break;
      }
    }
      
    if ( prime )
    {
      primes[p] = i;
      p = p + 1;
    }
    i = i + 1;
  }

  return primes;
}
/******************************************************************************/

double *sine_table ( int sine_num )

/******************************************************************************/
/*
  Purpose:

    SINE_TABLE computes a table of sines.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 October 2011

  Author:

    John Burkardt

  Parameters:

    Input, int SINE_NUM, the number of sines to compute.

    Output, double SINE_TABLE[SINE_NUM], the sines.
*/
{
  double a;
  int i;
  int j;
  double pi = 3.141592653589793;
  double *sines;

  sines = ( double * ) malloc ( sine_num * sizeof ( double ) );

  for ( i = 0; i < sine_num; i++ )
  {
    sines[i] = 0.0;
    for ( j = 0; j <= i; j++ )
    {
      a = ( double ) ( j ) * pi / ( double ) ( sine_num - 1 );
      sines[i] = sines[i] + sin ( a );
    }
  }

  return sines;
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
