# include <stdlib.h>
# include <stdio.h>
# include <time.h>
# include <omp.h>

int main ( int argc, char *argv[] );
int circuit_value ( int n, int bvec[] );
void i4_to_bvec ( int i4, int n, int bvec[] );
void timestamp ( void );

/******************************************************************************/

int main ( int argc, char *argv[] )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for SATISFY_OPENMP.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    24 March 2009

  Author:

    John Burkardt

  Reference:

    Michael Quinn,
    Parallel Programming in C with MPI and OpenMP,
    McGraw-Hill, 2004,
    ISBN13: 978-0071232654,
    LC: QA76.73.C15.Q55.
*/
{
# define N 23

  int bvec[N];
  int i;
  int id;
  int ihi;
  int ihi2;
  int ilo;
  int ilo2;
  int j;
  int n = N;
  int proc_num;
  int solution_num;
  int solution_num_local;
  int thread_num;
  int value;
  double wtime;

  printf ( "\n" );
  timestamp ( );
  printf ( "\n" );
  printf ( "SATISFY_OPENMP\n" );
  printf ( "  C + OpenMP version\n" );
  printf ( "  We have a logical function of N logical arguments.\n" );
  printf ( "  We do an exhaustive search of all 2^N possibilities,\n" );
  printf ( "  seeking those inputs that make the function TRUE.\n" );

  printf ( "\n" );
  printf ( "  Number of processors available = %d\n", omp_get_num_procs ( ) );
  printf ( "  Number of threads =              %d\n", omp_get_max_threads ( ) );
/*
  Compute the number of binary vectors to check.
*/
  ilo = 0;
  ihi = 1;
  for ( i = 1; i <= n; i++ )
  {
    ihi = ihi * 2;
  }
  printf ( "\n" );
  printf ( "  The number of logical variables is N = %d\n", n );
  printf ( "  The number of input vectors to check is %d\n", ihi );
  printf ( "\n" );
  printf ( "   # Processor       Index    ---------Input Values------------------------\n" );
  printf ( "\n" );
/*
  Processor ID takes the interval ILO2 <= I < IHI2.
  Using the formulas below yields a set of nonintersecting intervals
  which cover the original interval [ILO,IHI).
*/
  thread_num = omp_get_max_threads ( );

  solution_num = 0;

  wtime = omp_get_wtime ( );

# pragma omp parallel \
  shared ( ihi, ilo, n, thread_num ) \
  private ( bvec, i, id, ihi2, ilo2, j, solution_num_local, value ) \
  reduction ( + : solution_num )
  {
    id = omp_get_thread_num ( );

    ilo2 = ( ( thread_num - id     ) * ilo   
           + (              id     ) * ihi ) 
           / ( thread_num          );

    ihi2 = ( ( thread_num - id - 1 ) * ilo   
           + (              id + 1 ) * ihi ) 
           / ( thread_num          );

    printf ( "\n" );
    printf ( "  Processor %8d iterates from %8d <= I < %8d.\n", id, ilo2, ihi2 );
    printf ( "\n" );
/*
  Check every possible input vector.
*/
    solution_num_local = 0;

    for ( i = ilo2; i < ihi2; i++ )
    {
      i4_to_bvec ( i, n, bvec );

      value = circuit_value ( n, bvec );

      if ( value == 1 )
      {
        solution_num_local = solution_num_local + 1;
  
        printf ( "  %2d  %8d  %10d:  ", solution_num_local, id, i );
        for ( j = 0; j < n; j++ )
        {
          printf ( " %d", bvec[j] );
        }
        printf ( "\n" );
      }
    }
    solution_num = solution_num + solution_num_local;
  }
  wtime = omp_get_wtime ( ) - wtime;
  printf ( "\n" );
  printf ( "  Number of solutions found was %d\n", solution_num );
  printf ( "  Elapsed wall clock time (seconds) %f\n", wtime );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "SATISFY_OPENMP\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
# undef N
}
/******************************************************************************/

int circuit_value ( int n, int bvec[] )

/******************************************************************************/
/*
  Purpose:

    CIRCUIT_VALUE returns the value of a circuit for a given input set.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 March 2009

  Author:

    John Burkardt

  Reference:

    Michael Quinn,
    Parallel Programming in C with MPI and OpenMP,
    McGraw-Hill, 2004,
    ISBN13: 978-0071232654,
    LC: QA76.73.C15.Q55.

  Parameters:

    Input, int N, the length of the input vector.

    Input, int BVEC[N], the binary inputs.

    Output, int CIRCUIT_VALUE, the output of the circuit.
*/
{
  int value;

  value = 
       (  bvec[0]  ||  bvec[1]  )
    && ( !bvec[1]  || !bvec[3]  )
    && (  bvec[2]  ||  bvec[3]  )
    && ( !bvec[3]  || !bvec[4]  )
    && (  bvec[4]  || !bvec[5]  )
    && (  bvec[5]  || !bvec[6]  )
    && (  bvec[5]  ||  bvec[6]  )
    && (  bvec[6]  || !bvec[15] )
    && (  bvec[7]  || !bvec[8]  )
    && ( !bvec[7]  || !bvec[13] )
    && (  bvec[8]  ||  bvec[9]  )
    && (  bvec[8]  || !bvec[9]  )
    && ( !bvec[9]  || !bvec[10] )
    && (  bvec[9]  ||  bvec[11] )
    && (  bvec[10] ||  bvec[11] )
    && (  bvec[12] ||  bvec[13] )
    && (  bvec[13] || !bvec[14] )
    && (  bvec[14] ||  bvec[15] )
    && (  bvec[14] ||  bvec[16] )
    && (  bvec[17] ||  bvec[1]  )
    && (  bvec[18] || !bvec[0]  )
    && (  bvec[19] ||  bvec[1]  )
    && (  bvec[19] || !bvec[18] )
    && ( !bvec[19] || !bvec[9]  )
    && (  bvec[0]  ||  bvec[17] )
    && ( !bvec[1]  ||  bvec[20] )
    && ( !bvec[21] ||  bvec[20] )
    && ( !bvec[22] ||  bvec[20] )
    && ( !bvec[21] || !bvec[20] )
    && (  bvec[22] || !bvec[20] );

  return value;
}
/******************************************************************************/

void i4_to_bvec ( int i4, int n, int bvec[] )

/******************************************************************************/
/*
  Purpose:

    I4_TO_BVEC converts an integer into a binary vector.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 March 2009

  Author:

    John Burkardt

  Parameters:

    Input, int I4, the integer.

    Input, int N, the dimension of the vector.

    Output, int BVEC[N], the vector of binary remainders.
*/
{
  int i;

  for ( i = n - 1; 0 <= i; i-- )
  {
    bvec[i] = i4 % 2;
    i4 = i4 / 2;
  }

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
