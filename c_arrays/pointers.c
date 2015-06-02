# include <stdlib.h>
# include <stdio.h>
# include <time.h>

int main ( int argc, char *argv[] );
void test01 ( );
void test02 ( );
void r8cmat_delete ( double **a, int m, int n );
double **r8cmat_new ( int m, int n );
void r8rmat_delete ( double **a, int m, int n );
double **r8rmat_new ( int m, int n );
void timestamp ( void );

/******************************************************************************/

int main ( int argc, char *argv[] )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for POINTERS.

  Discussion:

    This program demonstrates the use of pointers to define objects
    that can be indexed like a two-dimensional array.  Either row-major
    or column-major storage can be handled.  However, the column major
    storage means that an expression like a[i][j] refers to the i-th
    column, j-th row, in other words, what we would naturally refer
    to as A(J,I) instead!

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    08 September 2013

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "POINTERS\n" );
  printf ( "  C version\n" );

  test01 ( );
  test02 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "POINTERS:\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void test01 ( )

/******************************************************************************/
/*
  Purpose:

    TEST01 demonstrates R8RMAT_NEW and R8RMAT_DELETE.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    08 September 2013

  Author:

    John Burkardt
*/
{
  double **a;
  double **b;
  int c;
  int ma = 4;
  int mb;
  int na = 3;
  int nb;
  int r;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  R8RMAT_NEW creates row-major doubly dimensioned arrays.\n" );
  printf ( "  R8RMAT_DELETE deletes them.\n" );
/*
  Create the matrix.
*/
  a = r8rmat_new ( ma, na );
/*
  Put some stuff in it.
*/
  for ( r = 0; r < ma; r++ )
  {
    for ( c = 0; c < na; c++ )
    {
      a[r][c] = 10.0 * r + c;
    }
  }
/*
  Print the matrix.
*/
  printf ( "\n" );
  printf ( "  Matrix A:\n" );
  printf ( "\n" );
  for ( r = 0; r < ma; r++ )
  {
    for ( c = 0; c < na; c++ )
    {
      printf ( "  %9f", a[r][c] );
    }
    printf ( "\n" );
  }
/*
  Create the transpose.
*/
  mb = na;
  nb = ma;

  b = r8rmat_new ( mb, nb );
  for ( r = 0; r < ma; r++ )
  {
    for ( c = 0; c < na; c++ )
    {
      b[c][r] = a[r][c];
    }
  }
/*
  Print the transpose;
*/
  printf ( "\n" );
  printf ( "  Matrix B = transpose ( A ):\n" );
  printf ( "\n" );
  for ( r = 0; r < mb; r++ )
  {
    for ( c = 0; c < nb; c++ )
    {
      printf ( "  %9f", b[r][c] );
    }
    printf ( "\n" );
  }
/*
  To print row R,
  we need to access row R, column C, so realize that 
  B[R] is a pointer to the first entry, so 
  B[R]+C is a pointer to the entry B(R,C) and 
  *(B[R]+C) is its value.
  These are CONSECUTIVE memory locations.
*/
  r = 1;
  printf ( "\n" );
  printf ( "  Row %d of matrix B:\n", r );
  printf ( "  (consecutive memory locations)\n" );
  printf ( "\n" );
  for ( c = 0; c < nb; c++ )
  {
    printf ( "  %9f", *(b[r]+c) );
  }
  printf ( "\n" );
/*
  To print column C, we end up using a similar formula,
  but because B[R] changes on each step, these are NOT consecutive
  memory locations.
*/
  c = 2;
  printf ( "\n" );
  printf ( "  Column %d of matrix B:\n", c );
  printf ( "\n" );
  for ( r = 0; r < mb; r++ )
  {
    printf ( "  %9f", *(b[r]+c) );
  }
  printf ( "\n" );
/*
  Free memory.
*/
  r8rmat_delete ( a, ma, na );
  r8rmat_delete ( b, mb, nb );

  return;
}
/******************************************************************************/

void test02 ( )

/******************************************************************************/
/*
  Purpose:

    TEST02 demonstrates R8CMAT_NEW and R8CMAT_DELETE.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    08 September 2013

  Author:

    John Burkardt
*/
{
  double **a;
  double **b;
  int c;
  int ma = 4;
  int mb;
  int na = 3;
  int nb;
  int r;

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  R8CMAT_NEW creates column-major doubly dimensioned arrays.\n" );
  printf ( "  R8CMAT_DELETE deletes them.\n" );
  printf ( "\n" );
  printf ( "  Unfortunately, A(i,j) is referenced as A[J][I]!\n" );
/*
  Create the matrix.
*/
  a = r8cmat_new ( ma, na );
/*
  Put some stuff in it.
*/
  for ( r = 0; r < ma; r++ )
  {
    for ( c = 0; c < na; c++ )
    {
      a[c][r] = 10.0 * r + c;
    }
  }
/*
  Print the matrix.
*/
  printf ( "\n" );
  printf ( "  Matrix A:\n" );
  printf ( "\n" );
  for ( r = 0; r < ma; r++ )
  {
    for ( c = 0; c < na; c++ )
    {
      printf ( "  %9f", a[c][r] );
    }
    printf ( "\n" );
  }
/*
  Create the transpose.
*/
  mb = na;
  nb = ma;

  b = r8cmat_new ( mb, nb );
  for ( r = 0; r < ma; r++ )
  {
    for ( c = 0; c < na; c++ )
    {
      b[r][c] = a[c][r];
    }
  }
/*
  Print the transpose;
*/
  printf ( "\n" );
  printf ( "  Matrix B = transpose ( A ):\n" );
  printf ( "\n" );
  for ( r = 0; r < mb; r++ )
  {
    for ( c = 0; c < nb; c++ )
    {
      printf ( "  %9f", b[c][r] );
    }
    printf ( "\n" );
  }
/*
  To print row R, B[C] changes on each step, these are NOT consecutive
  memory locations.
*/
  r = 1;
  printf ( "\n" );
  printf ( "  Row %d of matrix B:\n", r );
  printf ( "\n" );
  for ( c = 0; c < nb; c++ )
  {
    printf ( "  %9f", *(b[c]+r) );
  }
  printf ( "\n" );
/*
  To print column C,
  we need to access row R, column C, so realize that 
  B[C] is a pointer to the first entry of the column, so 
  B[C]+R is a pointer to the entry B(R,C) and 
  *(B[C]+R) is its value.
  These are CONSECUTIVE memory locations.

*/
  c = 2;
  printf ( "\n" );
  printf ( "  Column %d of matrix B:\n", c );
  printf ( "  (consecutive memory locations)\n" );
  printf ( "\n" );
  for ( r = 0; r < mb; r++ )
  {
    printf ( "  %9f", *(b[c]+r) );
  }
  printf ( "\n" );
/*
  Free memory.
*/
  r8cmat_delete ( a, ma, na );
  r8cmat_delete ( b, mb, nb );

  return;
}
/******************************************************************************/

void r8cmat_delete ( double **a, int m, int n )

/******************************************************************************/
/*
  Purpose:

    R8CMAT_DELETE frees the memory set aside by R8CMAT_NEW.

  Discussion:

    This function releases the memory associated with a
    column-major array that was created by a command like:

      double **a;
      a = r8cmat_new ( m, n );

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    08 September 2013

  Author:

    John Burkardt

  Parameters:

    Input, double **A, the array.

    Input, int M, N, the number of rows and columns.
*/
{
  int j;
/*
  Delete the pointers to columns.
*/
  for ( j = 0; j < n; j++ )
  {
    free ( a[j] );
  }
/*
  Delete the pointer to the pointers to rows.
*/
  free ( a );

  return;
}
/******************************************************************************/

double **r8cmat_new ( int m, int n )

/******************************************************************************/
/*
  Purpose:

    R8CMAT_NEW sets up an R8CMAT of the desired dimensions.

  Discussion:

    An R8CMAT is a matrix stored in column major form, using N pointers
    to the beginnings of columns.

    A declaration of the form
      double **a;
    is necesary.  Then an assignment of the form:
      a = r8cmat_new ( m, n );
    allows the user to assign entries to the matrix using typical
    2D array notation, except that the column index is listed FIRST:
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

    Output, double **R8CMAT_NEW, the array.
*/
{
  double **a;
  int j;

  a = ( double ** ) malloc ( n * sizeof ( double * ) );

  for ( j = 0; j < n; j++ )
  {
    a[j] = ( double * ) malloc ( m * sizeof ( double ) );
  }
  return a;
}
/******************************************************************************/

void r8rmat_delete ( double **a, int m, int n )

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

    Input, double **A, the array.

    Input, int M, N, the number of rows and columns.
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
/********************************************************************/

void timestamp ( void )

/********************************************************************/
/*
  Purpose:

    TIMESTAMP prints the current YMDHMS date as a time stamp.

  Example:

    May 31 2001 09:45:54 AM

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    03 October 2003

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
