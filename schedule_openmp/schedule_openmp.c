# include <stdlib.h>
# include <stdio.h>

# include <omp.h>

int main ( int argc, char *argv[] );
int prime_default ( int n );
int prime_static ( int n );
int prime_dynamic ( int n );

/******************************************************************************/

int main ( int argc, char *argv[] )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for SCHEDULE_OPENMP.

  Discussion:

    This program demonstrates the difference between default,
    static and dynamic scheduling for a loop parallelized in OpenMP.

    The purpose of scheduling is to deal with loops in which there is
    known or suspected imbalance in the work load.  In this example,
    if the work is divided in the default manner between two threads,
    the second thread has 3 times the work of the first.  

    Both static and dynamic scheduling, if used, even out the work
    so that both threads have about the same load.  This could be
    expected to decrease the run time of the loop by about 1/3.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    10 July 2010

  Author:

    John Burkardt
*/
{
  int n;
  int n_factor;
  int n_hi;
  int n_lo;
  int primes;
  double time1;
  double time2;
  double time3;

  printf ( "\n" );
  printf ( "SCHEDULE_OPENMP\n" );
  printf ( "  C/OpenMP version\n" );
  printf ( "  Count the primes from 1 to N.\n" );
  printf ( "  This is an unbalanced work load, particular for two threads.\n" );
  printf ( "  Demonstrate default, static and dynamic scheduling.\n" );
  printf ( "\n" );
  printf ( "  Number of processors available = %d\n", omp_get_num_procs ( )  );
  printf ( "  Number of threads =              %d\n", omp_get_max_threads ( )  );

  n_lo = 1;
  n_hi = 131072;
  n_factor = 2;

  printf ( "\n" );
  printf ( "                           Default        Static       Dynamic\n" );
  printf ( "         N     Pi(N)          Time          Time          Time\n" );
  printf ( "\n" );

  n = n_lo;

  while ( n <= n_hi )
  {
    time1 = omp_get_wtime ( );
    primes = prime_default ( n );
    time1 = omp_get_wtime ( ) - time1;

    time2 = omp_get_wtime ( );
    primes = prime_static ( n );
    time2 = omp_get_wtime ( ) - time2;

    time3 = omp_get_wtime ( );
    primes = prime_dynamic ( n );
    time3 = omp_get_wtime ( ) - time3;

    printf ( "  %8d  %8d  %12f  %12f  %12f\n", n, primes, time1, time2, time3 );

    n = n * n_factor;
  }
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "SCHEDULE_OPENMP\n" );
  printf ( "  Normal end of execution.\n" );

  return 0;
}
/******************************************************************************/

int prime_default ( int n )

/******************************************************************************/
/*
  Purpose:

    PRIME_DEFAULT counts primes, using default scheduling.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    10 July 2010

  Author:

    John Burkardt

  Parameters:

    Input, int N, the maximum number to check.

    Output, int PRIME_DEFAULT, the number of prime numbers up to N.
*/
{
  int i;
  int j;
  int prime;
  int total = 0;

# pragma omp parallel \
  shared ( n ) \
  private ( i, j, prime )

# pragma omp for reduction ( + : total )
  for ( i = 2; i <= n; i++ )
  {
    prime = 1;

    for ( j = 2; j < i; j++ )
    {
      if ( i % j == 0 )
      {
        prime = 0;
        break;
      }
    }
    total = total + prime;
  }

  return total;
}
/******************************************************************************/

int prime_static ( int n )

/******************************************************************************/
/*
  Purpose:

    PRIME_STATIC counts primes using static scheduling.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    10 July 2010

  Author:

    John Burkardt

  Parameters:

    Input, int N, the maximum number to check.

    Output, int PRIME_STATIC, the number of prime numbers up to N.
*/
{
  int i;
  int j;
  int prime;
  int total = 0;

# pragma omp parallel \
  shared ( n ) \
  private ( i, j, prime )

# pragma omp for reduction ( + : total ) schedule ( static, 100 )
  for ( i = 2; i <= n; i++ )
  {
    prime = 1;

    for ( j = 2; j < i; j++ )
    {
      if ( i % j == 0 )
      {
        prime = 0;
        break;
      }
    }
    total = total + prime;
  }

  return total;
}
/******************************************************************************/

int prime_dynamic ( int n )

/******************************************************************************/
/*
  Purpose:

    PRIME_DYNAMIC counts primes using dynamic scheduling.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    10 July 2010

  Author:

    John Burkardt

  Parameters:

    Input, int N, the maximum number to check.

    Output, int PRIME_DYNAMIC, the number of prime numbers up to N.
*/
{
  int i;
  int j;
  int prime;
  int total = 0;

# pragma omp parallel \
  shared ( n ) \
  private ( i, j, prime )

# pragma omp for reduction ( + : total ) schedule ( dynamic, 100 )
  for ( i = 2; i <= n; i++ )
  {
    prime = 1;

    for ( j = 2; j < i; j++ )
    {
      if ( i % j == 0 )
      {
        prime = 0;
        break;
      }
    }
    total = total + prime;
  }

  return total;
}
