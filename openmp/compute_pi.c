# include <stdlib.h>
# include <stdio.h>
# include <string.h>
# include <omp.h>

int main ( int argc, char *argv[] );
void r8_test ( int r8_logn_max );
double r8_abs ( double r8 );
double r8_pi_est_omp ( int n );
double r8_pi_est_seq ( int n );
double r8_pi ( void );

/******************************************************************************/

int main ( int argc, char *argv[] )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for COMPUTE_PI.

  Discussion:

    COMPUTE_PI estimates the value of PI.

    This program uses Open MP parallelization directives.  

    It should run properly whether parallelization is used or not.

    However, the parallel version computes the sum in a different
    order than the serial version; some of the quantities added are
    quite small, and so this will affect the accuracy of the results.

    The single precision code is noticeably less accurate than the
    double code.  Again, the amount of difference depends
    on whether the computation is done in parallel or not.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    18 April 2009

  Author:

    John Burkardt
*/
{
  int r8_logn_max = 10;

  printf ( "\n" );
  printf ( "COMPUTE_PI\n" );
  printf ( "  C/OpenMP version\n" );
  printf ( "\n" );
  printf ( "  Estimate the value of PI by summing a series.\n" );

  printf ( "\n" );
  printf ( "  Number of processors available = %d\n", omp_get_num_procs ( ) );
  printf ( "  Number of threads =              %d\n", omp_get_max_threads ( ) );

  r8_test ( r8_logn_max );

  printf ( "\n" );
  printf ( "COMPUTE_PI\n" );
  printf ( "  Normal end of execution.\n" );

  return 0;
}
/******************************************************************************/

void r8_test ( int logn_max )

/******************************************************************************/
/*
  Purpose:

    R8_TEST estimates the value of PI using double.

  Discussion:

    PI is estimated using N terms.  N is increased from 10^2 to 10^LOGN_MAX.
    The calculation is repeated using both sequential and Open MP enabled code.
    Wall clock time is measured by calling SYSTEM_CLOCK.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    13 November 2007

  Author:

    John Burkardt
*/
{
  double error;
  double estimate;
  int logn;
  char mode[4];
  int n;
  double r8_pi = 3.141592653589793;
  double wtime;

  printf ( "\n" );
  printf ( "R8_TEST:\n" );
  printf ( "  Estimate the value of PI,\n" );
  printf ( "  using double arithmetic.\n" );
  printf ( "\n" );
  printf ( "  N = number of terms computed and added;\n" );
  printf ( "\n" );
  printf ( "  MODE = SEQ for sequential code;\n" );
  printf ( "  MODE = OMP for Open MP enabled code;\n" );
  printf ( "  (performance depends on whether Open MP is used,\n" );
  printf ( "  and how many processes are available)\n" );
  printf ( "\n" );
  printf ( "  ESTIMATE = the computed estimate of PI;\n" );
  printf ( "\n" );
  printf ( "  ERROR = ( the computed estimate - PI );\n" );
  printf ( "\n" );
  printf ( "  TIME = elapsed wall clock time;\n" );
  printf ( "\n" );
  printf ( "  Note that you can''t increase N forever, because:\n" );
  printf ( "  A) ROUNDOFF starts to be a problem, and\n" );
  printf ( "  B) maximum integer size is a problem.\n" );
  printf ( "\n" );
  printf ( "             N Mode    Estimate        Error           Time\n" );
  printf ( "\n" );

  n = 1;

  for ( logn = 1; logn <= logn_max; logn++ )
  {
/*
  Note that when I set N = 10**LOGN directly, rather than using
  recursion, I got inaccurate values of N when LOGN was "large",
  that is, for LOGN = 10, despite the fact that N itself was
  a KIND = 8 integer!  
 
  Sequential calculation.
*/
    strcpy ( mode, "SEQ" );

    wtime = omp_get_wtime ( );

    estimate = r8_pi_est_seq ( n );

    wtime = omp_get_wtime ( ) - wtime;

    error = r8_abs ( estimate - r8_pi );

    printf ( "%14d  %s  %14f  %14g  %14f\n", n, mode, estimate, error, wtime );
/*
  Open MP enabled calculation.
*/
    strcpy ( mode, "OMP" );

    wtime = omp_get_wtime ( );

    estimate = r8_pi_est_omp ( n );

    wtime = omp_get_wtime ( ) - wtime;

    error = r8_abs ( estimate - r8_pi );

    printf ( "%14d  %s  %14f  %14g  %14f\n", n, mode, estimate, error, wtime );

    n = n * 10;
  }

  return;
}
/******************************************************************************/

double r8_abs ( double r8 )

/******************************************************************************/
/*
  Purpose:
 
    R8_ABS returns the absolute value of a double.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    14 November 2007

  Author:

    John Burkardt

  Parameters:

    Input, double R8, a number.

    Output, double R8_ABS, the absolute value of R4.
*/
{
  double value;

  if ( 0.0 <= r8 )
  {
    value = r8;
  }
  else
  {
    value = - r8;
  }
  return value;
}
/******************************************************************************/

double r8_pi_est_omp ( int n )

/******************************************************************************/
/*
  Purpose:

    R8_PI_EST_OMP estimates the value of PI, using Open MP.

  Discussion:

    The calculation is based on the formula for the indefinite integral:

      Integral 1 / ( 1 + X**2 ) dx = Arctan ( X ) 

    Hence, the definite integral

      Integral ( 0 <= X <= 1 ) 1 / ( 1 + X**2 ) dx 
      = Arctan ( 1 ) - Arctan ( 0 )
      = PI / 4.

    A standard way to approximate an integral uses the midpoint rule.
    If we create N equally spaced intervals of width 1/N, then the
    midpoint of the I-th interval is 

      X(I) = (2*I-1)/(2*N).  

    The approximation for the integral is then:

      Sum ( 1 <= I <= N ) (1/N) * 1 / ( 1 + X(I)**2 )

    In order to compute PI, we multiply this by 4; we also can pull out
    the factor of 1/N, so that the formula you see in the program looks like:

      ( 4 / N ) * Sum ( 1 <= I <= N ) 1 / ( 1 + X(I)**2 )

    Until roundoff becomes an issue, greater accuracy can be achieved by 
    increasing the value of N.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    13 November 2007

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of terms to add up.

    Output, double R8_PI_EST_OMP, the estimated value of pi.
*/
{
  double h;
  double estimate;
  int i;
  double sum2;
  double x;

  h = 1.0 / ( double ) ( 2 * n );

  sum2 = 0.0;

# pragma omp parallel \
  shared ( h, n ) \
  private ( i, x )

# pragma omp for reduction ( +: sum2 )

  for ( i = 1; i <= n; i++ )
  {
    x = h * ( double ) ( 2 * i - 1 );
    sum2 = sum2 + 1.0 / ( 1.0 + x * x );
  }

  estimate = 4.0 * sum2 / ( double ) ( n );

  return estimate;
}
/******************************************************************************/

double r8_pi_est_seq ( int n )

/******************************************************************************/
/*
  Purpose:

    R8_PI_EST_SEQ estimates the value of PI, using sequential execution.

  Discussion:

    The calculation is based on the formula for the indefinite integral:

      Integral 1 / ( 1 + X**2 ) dx = Arctan ( X ) 

    Hence, the definite integral

      Integral ( 0 <= X <= 1 ) 1 / ( 1 + X**2 ) dx 
      = Arctan ( 1 ) - Arctan ( 0 )
      = PI / 4.

    A standard way to approximate an integral uses the midpoint rule.
    If we create N equally spaced intervals of width 1/N, then the
    midpoint of the I-th interval is 

      X(I) = (2*I-1)/(2*N).  

    The approximation for the integral is then:

      Sum ( 1 <= I <= N ) (1/N) * 1 / ( 1 + X(I)**2 )

    In order to compute PI, we multiply this by 4; we also can pull out
    the factor of 1/N, so that the formula you see in the program looks like:

      ( 4 / N ) * Sum ( 1 <= I <= N ) 1 / ( 1 + X(I)**2 )

    Until roundoff becomes an issue, greater accuracy can be achieved by 
    increasing the value of N.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    13 November 2007

  Author:

    John Burkardt

  Parameters:

    Input, integer N, the number of terms to add up.

    Output, double R8_PI_EST_SEQ, the estimated value of pi.
*/
{
  double h;
  double estimate;
  int i;
  double sum2;
  double x;

  h = 1.0 / ( double ) ( 2 * n );

  sum2 = 0.0;

  for ( i = 1; i <= n; i++ )
  {
    x = h * ( double ) ( 2 * i - 1 );
    sum2 = sum2 + 1.0 / ( 1.0 + x * x );
  }

  estimate = 4.0 * sum2 / ( double ) ( n );

  return estimate;
}

