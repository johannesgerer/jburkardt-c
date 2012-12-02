# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <omp.h>
# include <time.h>

int main ( void );
void monte_carlo ( int n, int *seed );
double random_value ( int *seed );
void timestamp ( void );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for RANDOM_OPENMP.

  Discussion:

    This program simply explores one issue in the generation of random
    numbers in a parallel program.  If the random number generator uses
    an integer seed to determine the next entry, then it is not easy for
    a parallel program to reproduce the same exact sequence.

    But what is worse is that it might not be clear how the separate
    OpenMP threads should handle the SEED value - as a shared or private
    variable?  It seems clear that each thread should have a private
    seed that is initialized to a distinct value at the beginning of
    the computation.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    03 September 2012

  Author:

    John Burkardt
*/
{
  int n;
  int seed;

  timestamp ( );

  printf ( "\n" );
  printf ( "RANDOM_OPENMP\n" );
  printf ( "  C version\n" );
  printf ( "  An OpenMP program using random numbers.\n" );
  printf ( "  The random numbers depend on a seed.\n" );
  printf ( "  We need to insure that each OpenMP thread\n" );
  printf ( "  starts with a different seed.\n" );
  printf ( "\n" );
  printf ( "  Number of processors available = %d\n", omp_get_num_procs ( ) );
  printf ( "  Number of threads =              %d\n", omp_get_max_threads ( ) );

  n = 100;
  seed = 123456789;
  monte_carlo ( n, &seed );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "RANDOM_OPENMP\n" );
  printf ( "  Normal end of execution.\n" );

  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void monte_carlo ( int n, int *seed )

/******************************************************************************/
/*
  Purpose:

    MONTE_CARLO carries out a Monte Carlo calculation with random values.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    03 September 2012

  Author:

    John Burkardt

  Parameter:

    Input, int N, the number of values to generate.

    Input, int *SEED, a seed for the random number generator.
*/
{
  int i;
  int my_id;
  int my_seed;
  double *x;

  x = ( double * ) malloc ( n * sizeof ( double ) );

# pragma omp master
{
  printf ( "\n" );
  printf ( "  Thread   Seed  I   X(I)\n" );
  printf ( "\n" );
}

# pragma omp parallel private ( i, my_id, my_seed ) shared ( n, x )
{
  my_id = omp_get_thread_num ( );
  my_seed = *seed + my_id;
  printf ( "  %6d  %12d\n", my_id, my_seed );

# pragma omp for
  for ( i = 0; i < n; i++ )
  {
    x[i] = random_value ( &my_seed );
    printf ( "  %6d  %12d  %6d  %14.6g\n", my_id, my_seed, i, x[i] );
  }

}

  free ( x );

  return;
}
/******************************************************************************/

double random_value ( int *seed )

/******************************************************************************/
/*
  Purpose:

    RANDOM_VALUE generates a random value R.

  Discussion:

    This is not a good random number generator.  It is a SIMPLE one.
    It illustrates a model which works by accepting an integer seed value
    as input, performing some simple operation on the seed, and then
    producing a "random" real value using some simple transformation.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    03 September 2012

  Author:

    John Burkardt

  Parameters:

    Input/output, int *SEED, a seed for the random 
    number generator.

    Output, double RANDOM_VALUE, the random value.
*/
{
  double r;

  *seed = ( *seed % 65536 );
  *seed = ( ( 3125 * *seed ) % 65536 );
  r = ( double ) ( *seed ) / 65536.0;

  return r;
}
/******************************************************************************/

void timestamp ( void )

/******************************************************************************/
/*
  Purpose:

    TIMESTAMP prints the current YMDHMS date as a time stamp.

  Example:

    31 May 2001 09:45:54 AM

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
