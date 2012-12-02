# include <stdlib.h>
# include <stdio.h>
# include <omp.h>

int main ( int argc, char *argv[] );

/******************************************************************************/

int main ( int argc, char *argv[] )

/******************************************************************************/
/*
  Purpose:

    HELLO has each thread print out its ID.

  Discussion:

    HELLO is a "Hello, World" program for OpenMP.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    23 June 2010

  Author:

    John Burkardt
*/
{
  int id;
  double wtime;

  printf ( "\n" );
  printf ( "HELLO_OPENMP\n" );
  printf ( "  C/OpenMP version\n" );

  printf ( "\n" );
  printf ( "  Number of processors available = %d\n", omp_get_num_procs ( ) );
  printf ( "  Number of threads =              %d\n", omp_get_max_threads ( ) );

  wtime = omp_get_wtime ( );

  printf ( "\n" );
  printf ( "  OUTSIDE the parallel region.\n" );
  printf ( "\n" );

  id = omp_get_thread_num ( );
  printf ( "  HELLO from process %d\n", id ) ;

  printf ( "\n" );
  printf ( "  Going INSIDE the parallel region:\n" );
  printf ( "\n" );
/*
  INSIDE THE PARALLEL REGION, have each thread say hello.
*/
# pragma omp parallel \
  private ( id )
  {
    id = omp_get_thread_num ( );
    printf ("  Hello from process %d\n", id );
  }
/*
  Finish up by measuring the elapsed time.
*/
  wtime = omp_get_wtime ( ) - wtime;

  printf ( "\n" );
  printf ( "  Back OUTSIDE the parallel region.\n" );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "HELLO_OPENMP\n" );
  printf ( "  Normal end of execution.\n" );

  printf ( "\n" );
  printf ( "  Elapsed wall clock time = %f\n", wtime );

  return 0;
}
