#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>

# define NUMTHRDS 4
double sum;
double a[256], b[256];
int status;
int n = 256;
pthread_t thds[NUMTHRDS];
pthread_mutex_t mutex_sum;

int main ( int argc, char *argv[] );
void *dotprod ( void *arg );

/******************************************************************************/

int main ( int argc, char *argv[] )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for DOT_PRODUCT.

  Reference:

    Barbara Chapman, Gabriele Jost, Ruud vanderPas, David Kuck,
    Using OpenMP: Portable Shared Memory Parallel Processing,
    MIT Press, 2007,
    ISBN13: 978-0262533027,
    LC: QA76.642.C49.
*/
{
  int i;

  pthread_attr_t attr;

  printf ( "\n" );
  printf ( "DOT_PRODUCT:\n" );
  printf ( "  C++ version\n" );
  printf ( "\n" );
  printf ( "  Demonstrate the use of Posix Threads by computing the\n" );
  printf ( "  dot product of two vectors.\n" );
  printf ( "\n" );
  printf ( "  The length of the vectors N = %d\n", n );
  printf ( "  The number of threads used is %d\n", NUMTHRDS );

  for ( i = 0; i < n; i++ )
  {
    a[i] = i * 0.5;
    b[i] = i * 2.0;
  }

  pthread_mutex_init ( &mutex_sum, NULL );
  pthread_attr_init ( &attr );
  pthread_attr_setdetachstate ( &attr, PTHREAD_CREATE_JOINABLE );

  for ( i = 0; i < NUMTHRDS; i++ )
  {
    pthread_create ( &thds[i], &attr, dotprod, ( void * ) i );
  }

  pthread_attr_destroy ( &attr );

  for ( i = 0; i < NUMTHRDS; i++ )
  {
    pthread_join ( thds[i], ( void ** ) &status );
  }

  printf ( "\n" );
  printf ( "  Sum = %f\n", sum );

  pthread_mutex_destroy ( &mutex_sum );
  pthread_exit ( NULL );

  printf ( "\n" );
  printf ( "DOT_PRODUCT:\n" );
  printf ( "  Normal end of execution.\n" );

  return 0;
}
/******************************************************************************/

void *dotprod ( void *arg )

/******************************************************************************/
/*
  Purpose:

    DOT_PRODUCT is executed by each thread.

  Reference:

    Barbara Chapman, Gabriele Jost, Ruud vanderPas, David Kuck,
    Using OpenMP: Portable Shared Memory Parallel Processing,
    MIT Press, 2007,
    ISBN13: 978-0262533027,
    LC: QA76.642.C49.
*/
{
  int i;
  int my_first;
  int my_last;
  int myid;
  double sum_local;

  myid = ( int ) arg;
/*
  Determine the portion of the dot product to be computed by this thread.
*/
  my_first = myid * n / NUMTHRDS;
  my_last = ( myid + 1 ) * n / NUMTHRDS;
/*
  Compute a part of the dot product.
*/
  sum_local = 0;
  for ( i = my_first; i <= my_last; i++ )
  {
    sum_local = sum_local + a[i] * b[i];
  }
/*
  Lock the variable MUTEX_SUM, update it, and then unlock it.
*/
  pthread_mutex_lock ( &mutex_sum );
  sum = sum + sum_local;
  pthread_mutex_unlock ( &mutex_sum );

  pthread_exit ( ( void * ) 0 );
}
