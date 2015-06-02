# include <stdio.h>
# include <stdlib.h>
# include <time.h>
# include <math.h>

# include "umfpack.h"

int main ( );
void timestamp ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for UMFPACK_SIMPLE.

  Discussion:

    This program uses UMFPACK to solve the 5x5 linear system A*X=B:

        2  3  0  0  0        1.0         8.0
        3  0  4  0  6        2.0        45.0
    A = 0 -1 -3  2  0    X = 3.0    B = -3.0
        0  0  1  0  0        4.0         3.0
        0  4  2  0  1        5.0        10.0

    The matrix contains 12 nonzero values.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 July 2014

  Author:

    John Burkardt

  Reference:

    Timothy Davis,
    UMFPACK User Guide,
    Version 5.6.2, 25 April 2013
    http://suitesparse.com
*/
{
# define N 5
# define NCC 12

  int Ai[NCC] = { 
    0, 1, 
    0, 2, 4, 
    1, 2, 3, 4, 
    2, 
    1, 4 };
  int Ap[N+1] = { 0, 2, 5, 9, 10, 12 };
  double Ax[NCC] = { 
    2.0,  3.0, 
    3.0, -1.0, 4.0, 
    4.0, -3.0, 1.0, 2.0, 
    2.0, 
    6.0, 1.0 };
  double b[N] = { 8.0, 45.0, -3.0, 3.0, 19.0 };
  int i;
  int n = N;
  double *null = ( double * ) NULL;
  void *Numeric;
  int status;
  void *Symbolic;
  double x[N];

  timestamp ( );
  printf ( "\n" );
  printf ( "UMFPACK_SIMPLE:\n" );
  printf ( "  C version\n" );
  printf ( "  Use UMFPACK to solve the sparse linear system A*x=b.\n" );
/*
  From the matrix data, create the symbolic factorization information.
*/
  status = umfpack_di_symbolic ( n, n, Ap, Ai, Ax, &Symbolic, null, null );
/*
  From the symbolic factorization information, carry out the numeric factorization.
*/
  status = umfpack_di_numeric ( Ap, Ai, Ax, Symbolic, &Numeric, null, null );
/*
  Free the symbolic factorization memory.
*/
  umfpack_di_free_symbolic ( &Symbolic );
/*
  Using the numeric factorization, solve the linear system.
*/
  status = umfpack_di_solve ( UMFPACK_A, Ap, Ai, Ax, x, b, Numeric, null, null );
/*
  Free the numeric factorization.
*/
  umfpack_di_free_numeric ( &Numeric );
/*
  Print the solution.
*/
  printf ( "\n" );
  printf ( "  Computed solution:\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ ) 
  {
    printf ( "  x[%d] = %g\n", i, x[i] );
  }
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "UMFPACK_SIMPLE:\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
# undef N
# undef NCC
}
/******************************************************************************/

void timestamp ( )

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

  fprintf ( stdout, "%s\n", time_buffer );

  return;
# undef TIME_SIZE
}
