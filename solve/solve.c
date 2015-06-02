# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>

# include "solve.h"

/******************************************************************************/

double **r8rmat_copy_new ( int m, int n, double **a )

/******************************************************************************/
/*
  Purpose:

    R8RMAT_COPY_NEW makes a new copy of an R8RMAT .

  Discussion:

    An R8RMAT is a matrix stored in row major form, using M pointers
    to the beginnings of rows.

    A declaration of the form
      double **a;
    is necesary.  Then an assignment of the form:
      a = r8rmat_new ( m, n );
    allows the user to assign entries to the matrix using typical
    2D array notation:
      a[2][3] = 17.0;
      y = a[1][0];
    and so on.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    26 May 2014

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns.

    Input, double **A, the array to copy.

    Output, double **R8RMAT_COPY_NEW, the copied array.
*/
{
  double **b;
  int i;
  int j;

  b = r8rmat_new ( m, n );

  for ( i = 0; i < m; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      b[i][j] = a[i][j];
    }
  }  
  return b;
}
/******************************************************************************/

void r8rmat_delete ( int m, int n, double **a )

/******************************************************************************/
/*
  Purpose:

    R8RMAT_DELETE frees the memory set aside by R8RMAT_NEW.

  Discussion:

    This function releases the memory associated with a
    row-major array that was created by a command like:

      double **a;
      a = r8rmat_new ( m, n );

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    08 September 2013

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns.

    Input, double **A, the array.
*/
{
  int i;
/*
  Delete the pointers to rows.
*/
  for ( i = 0; i < m; i++ )
  {
    free ( a[i] );
  }
/*
  Delete the pointer to the pointers to rows.
*/
  free ( a );

  return;
}
/******************************************************************************/

double *r8rmat_fs_new ( int n, double **a, double b[] )

/******************************************************************************/
/*
  Purpose:

    R8RMAT_FS_NEW factors and solves an R8RMAT system with one right hand side.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    26 May 2014

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of the matrix.
    N must be positive.

    Input, double **A, the coefficient matrix of the linear system.

    Input, double B[N], the right hand side of the linear system.

    Output, double R8RMAT_FS_NEW[N], the solution of the linear system.
*/
{
  double **a2;
  int i;
  int j;
  int k;
  int p;
  double t;
  double *x;

  a2 = r8rmat_copy_new ( n, n, a );
  x = r8vec_copy_new ( n, b );

  for ( k = 0; k < n; k++ )
  {
/*
  Find the maximum element in column I.
*/
    p = k;

    for ( i = k + 1; i < n; i++ )
    {
      if ( fabs ( a2[p][k] ) < fabs ( a2[i][k] ) )
      {
        p = i;
      }
    }

    if ( a2[p][k] == 0.0 )
    {
      fprintf ( stderr, "\n" );
      fprintf ( stderr, "R8RMAT_FS_NEW - Fatal error!\n" );
      fprintf ( stderr, "  Zero pivot on step %d\n", k );
      exit ( 1 );
    }
/*
  Switch rows K and P.
*/
    if ( k != p )
    {
      for ( j = 0; j < n; j++ )
      {
        t        = a2[k][j];
        a2[k][j] = a2[p][j];
        a2[p][j] = t;
      }
      t    = x[k];
      x[k] = x[p];
      x[p] = t;
    }
/*
  Scale the pivot row.
*/
    t = a2[k][k];
    a2[k][k] = 1.0;
    for ( j = k + 1; j < n; j++ )
    {
      a2[k][j] = a2[k][j] / t;
    }
    x[k] = x[k] / t;
/*
  Use the pivot row to eliminate lower entries in that column.
*/
    for ( i = k + 1; i < n; i++ )
    {
      if ( a2[i][k] != 0.0 )
      {
        t = - a2[i][k];
        a2[i][k] = 0.0;
        for ( j = k + 1; j < n; j++ )
        {
          a2[i][j] = a2[i][j] + t * a2[k][j];
        }
        x[i] = x[i] + t * x[k];
      }
    }
  }
/*
  Back solve.
*/
  for ( j = n - 1; 1 <= j; j-- )
  {
    for ( i = 0; i < j; i++ )
    {
      x[i] = x[i] - a2[i][j] * x[j];
    }
  }

  r8rmat_delete ( n, n, a2 );

  return x;
}
/******************************************************************************/

double **r8rmat_new ( int m, int n )

/******************************************************************************/
/*
  Purpose:

    R8RMAT_NEW sets up an R8RMAT of the desired dimensions.

  Discussion:

    An R8RMAT is a matrix stored in row major form, using M pointers
    to the beginnings of rows.

    A declaration of the form
      double **a;
    is necesary.  Then an assignment of the form:
      a = r8rmat_new ( m, n );
    allows the user to assign entries to the matrix using typical
    2D array notation:
      a[2][3] = 17.0;
      y = a[1][0];
    and so on.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    08 September 2013

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns.

    Output, double **R8RMAT_NEW, the array.
*/
{
  double **a;
  int i;

  a = ( double ** ) malloc ( m * sizeof ( double * ) );

  for ( i = 0; i < m; i++ )
  {
    a[i] = ( double * ) malloc ( n * sizeof ( double ) );
  }
  return a;
}
/******************************************************************************/

double **r8rmat_zero ( int m, int n )

/******************************************************************************/
/*
  Purpose:

    R8RMAT_ZERO sets up and zeros an R8RMAT of the desired dimensions.

  Discussion:

    An R8RMAT is a matrix stored in row major form, using M pointers
    to the beginnings of rows.

    A declaration of the form
      double **a;
    is necesary.  Then an assignment of the form:
      a = r8rmat_new ( m, n );
    allows the user to assign entries to the matrix using typical
    2D array notation:
      a[2][3] = 17.0;
      y = a[1][0];
    and so on.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    26 May 2014

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns.

    Output, double **R8RMAT_ZERO, the array of zeros.
*/
{
  double **a;
  int i;
  int j;

  a = ( double ** ) malloc ( m * sizeof ( double * ) );

  for ( i = 0; i < m; i++ )
  {
    a[i] = ( double * ) malloc ( n * sizeof ( double ) );
  }
  for ( i = 0; i < m; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      a[i][j] = 0.0;
    }
  }
  return a;
}
/******************************************************************************/

double *r8vec_copy_new ( int n, double a1[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_COPY_NEW copies an R8VEC.

  Discussion:

    An R8VEC is a vector of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 August 2008

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the vectors.

    Input, double A1[N], the vector to be copied.

    Output, double R8VEC_COPY_NEW[N], the copy of A1.
*/
{
  double *a2;
  int i;

  a2 = ( double * ) malloc ( n * sizeof ( double ) );

  for ( i = 0; i < n; i++ )
  {
    a2[i] = a1[i];
  }
  return a2;
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
