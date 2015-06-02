# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>

# include "jacobi.h"

/******************************************************************************/

double *dif2 ( int m, int n )

/******************************************************************************/
/*
  Purpose:

    DIF2 returns the DIF2 matrix.

  Example:

    N = 5

    2 -1  .  .  .
   -1  2 -1  .  .
    . -1  2 -1  .
    .  . -1  2 -1
    .  .  . -1  2

  Properties:

    A is banded, with bandwidth 3.

    A is tridiagonal.

    Because A is tridiagonal, it has property A (bipartite).

    A is a special case of the TRIS or tridiagonal scalar matrix.

    A is integral, therefore det ( A ) is integral, and 
    det ( A ) * inverse ( A ) is integral.

    A is Toeplitz: constant along diagonals.

    A is symmetric: A' = A.

    Because A is symmetric, it is normal.

    Because A is normal, it is diagonalizable.

    A is persymmetric: A(I,J) = A(N+1-J,N+1-I).

    A is positive definite.

    A is an M matrix.

    A is weakly diagonally dominant, but not strictly diagonally dominant.

    A has an LU factorization A = L * U, without pivoting.

      The matrix L is lower bidiagonal with subdiagonal elements:

        L(I+1,I) = -I/(I+1)

      The matrix U is upper bidiagonal, with diagonal elements

        U(I,I) = (I+1)/I

      and superdiagonal elements which are all -1.

    A has a Cholesky factorization A = L * L', with L lower bidiagonal.

      L(I,I) =    sqrt ( (I+1) / I )
      L(I,I-1) = -sqrt ( (I-1) / I )

    The eigenvalues are

      LAMBDA(I) = 2 + 2 * COS(I*PI/(N+1))
                = 4 SIN^2(I*PI/(2*N+2))

    The corresponding eigenvector X(I) has entries

       X(I)(J) = sqrt(2/(N+1)) * sin ( I*J*PI/(N+1) ).

    Simple linear systems:

      x = (1,1,1,...,1,1),   A*x=(1,0,0,...,0,1)

      x = (1,2,3,...,n-1,n), A*x=(0,0,0,...,0,n+1)

    det ( A ) = N + 1.

    The value of the determinant can be seen by induction,
    and expanding the determinant across the first row:

      det ( A(N) ) = 2 * det ( A(N-1) ) - (-1) * (-1) * det ( A(N-2) )
                = 2 * N - (N-1)
                = N + 1

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 September 2010

  Author:

    John Burkardt

  Reference:

    Robert Gregory, David Karney,
    A Collection of Matrices for Testing Computational Algorithms,
    Wiley, 1969,
    ISBN: 0882756494,
    LC: QA263.68

    Morris Newman, John Todd,
    Example A8,
    The evaluation of matrix inversion programs,
    Journal of the Society for Industrial and Applied Mathematics,
    Volume 6, Number 4, pages 466-476, 1958.

    John Todd,
    Basic Numerical Mathematics,
    Volume 2: Numerical Algebra,
    Birkhauser, 1980,
    ISBN: 0817608117,
    LC: QA297.T58.

    Joan Westlake,
    A Handbook of Numerical Matrix Inversion and Solution of 
    Linear Equations,
    John Wiley, 1968,
    ISBN13: 978-0471936756,
    LC: QA263.W47.

  Parameters:

    Input, int M, N, the order of the matrix.

    Output, double DIF2[M*N], the matrix.
*/
{
  double *a;
  int i;
  int j;

  a = ( double * ) malloc ( m * n * sizeof ( double ) );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      if ( j == i - 1 )
      {
        a[i+j*m] = -1.0;
      }
      else if ( j == i )
      {
        a[i+j*m] = 2.0;
      }
      else if ( j == i + 1 )
      {
        a[i+j*m] = -1.0;
      }
      else
      {
        a[i+j*m] = 0.0;
      }
    }
  }

  return a;
}
/******************************************************************************/

double *jacobi1 ( int n, double a[], double b[], double x[] )

/******************************************************************************/
/*
  Purpose:

    JACOBI1 carries out one step of the Jacobi iteration.

  Discussion:

    The linear system A*x=b is to be solved.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    10 December 2013

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of the matrix.

    Input, double A[N,N], the matrix.

    Input, double B[N], the right hand side.

    Input, double X[N], the current solution estimate.

    Output, double JACOBI1[N], the solution estimate updated by
    one step of the Jacobi iteration.
*/
{
  int i;
  int j;
  double *x_new;

  x_new = ( double * ) malloc ( n * sizeof ( double ) );

  for ( i = 0; i < n; i++ )
  {
    x_new[i] = b[i];
    for ( j = 0; j < n; j++ )
    {
      if ( j != i )
      {
        x_new[i] = x_new[i] - a[i+j*n] * x[j];
      }
    }
    x_new[i] = x_new[i] / a[i+i*n];
  }

  return x_new;
}
/******************************************************************************/

double *r8mat_mv_new ( int m, int n, double a[], double x[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_MV_NEW multiplies a matrix times a vector.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

    For this routine, the result is returned as the function value.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 April 2007

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns of the matrix.

    Input, double A[M,N], the M by N matrix.

    Input, double X[N], the vector to be multiplied by A.

    Output, double R8MAT_MV_NEW[M], the product A*X.
*/
{
  int i;
  int j;
  double *y;

  y = ( double * ) malloc ( m * sizeof ( double ) );

  for ( i = 0; i < m; i++ )
  {
    y[i] = 0.0;
    for ( j = 0; j < n; j++ )
    {
      y[i] = y[i] + a[i+j*m] * x[j];
    }
  }

  return y;
}
/******************************************************************************/

double r8mat_residual_norm ( int m, int n, double a[], double x[], double b[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_RESIDUAL_NORM returns the norm of A*x-b.

  Discussion:

    A is an MxN R8MAT, a matrix of R8's.

    X is an N R8VEC, and B is an M R8VEC.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    24 June 2011

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns of the matrix.

    Input, double A[M,N], the M by N matrix.

    Input, double X[N], the vector to be multiplied by A.

    Input, double B[M], the right hand side vector.

    Output, double R8MAT_RESIDUAL_NORM, the norm of A*x-b.
*/
{
  int i;
  int j;
  double *r;
  double r_norm;

  r = ( double * ) malloc ( m * sizeof ( double ) );

  for ( i = 0; i < m; i++ )
  {
    r[i] = - b[i];
    for ( j = 0; j < n; j++ )
    {
      r[i] = r[i] + a[i+j*m] * x[j];
    }
  }

  r_norm = 0.0;
  for ( i = 0; i < m; i++ )
  {
    r_norm = r_norm + r[i] * r[i];
  }
  r_norm = sqrt ( r_norm );

  free ( r );

  return r_norm;
}
/******************************************************************************/

void r8vec_copy ( int n, double a1[], double a2[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_COPY copies an R8VEC.

  Discussion:

    An R8VEC is a vector of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    03 July 2005

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the vectors.

    Input, double A1[N], the vector to be copied.

    Input, double A2[N], the copy of A1.
*/
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    a2[i] = a1[i];
  }
  return;
}
/******************************************************************************/

double r8vec_diff_norm_squared ( int n, double a[], double b[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_DIFF_NORM_SQUARED: square of the L2 norm of the difference of R8VEC's.

  Discussion:

    An R8VEC is a vector of R8's.

    The square of the L2 norm of the difference of A and B is:

      R8VEC_DIFF_NORM_SQUARED = sum ( 1 <= I <= N ) ( A[I] - B[I] )^2.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    24 June 2011

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in A.

    Input, double A[N], B[N], the vectors.

    Output, double R8VEC_DIFF_NORM_SQUARED, the square of the L2 norm of A - B.
*/
{
  int i;
  double value;

  value = 0.0;

  for ( i = 0; i < n; i++ )
  {
    value = value + ( a[i] - b[i] ) * ( a[i] - b[i] );
  }

  return value;
}
/******************************************************************************/

void r8vec_print ( int n, double a[], char *title )

/******************************************************************************/
/*
  Purpose:

    R8VEC_PRINT prints an R8VEC.

  Discussion:

    An R8VEC is a vector of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 April 2009

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of components of the vector.

    Input, double A[N], the vector to be printed.

    Input, char *TITLE, a title.
*/
{
  int i;

  fprintf ( stdout, "\n" );
  fprintf ( stdout, "%s\n", title );
  fprintf ( stdout, "\n" );
  for ( i = 0; i < n; i++ )
  {
    fprintf ( stdout, "  %8d: %14g\n", i, a[i] );
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

  fprintf ( stdout, "%s\n", time_buffer );

  return;
# undef TIME_SIZE
}

