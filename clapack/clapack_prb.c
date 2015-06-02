# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
/*
  I have merged the blaswrap.h, f2c.h and clapack.h include files into one.
*/
# include "clapack.h"

int main ( );

void dgesv_test ( );
void dgesvd_test ( );
void dgetrf_test ( );
void dgetri_test ( );
void dsyev_test ( );
void zgesv_test ( );

double *clement2 ( int n );
void r8mat_print ( int m, int n, double a[], char *title );
void r8mat_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi,
  int jhi, char *title );
void r8vec_print ( int n, double a[], char *title );
void timestamp ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for CLAPACK_PRB.

  Discussion:

    CLAPACK_PRB tests the CLAPACK library.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    30 July 2013

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "CLAPACK_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the CLAPACK library.\n" );

  dgesv_test ( );
  dgesvd_test ( );
  dgetrf_test ( );
  dgetri_test ( );
  dsyev_test ( );
  zgesv_test ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "CLAPACK_PRB:\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void dgesv_test ( )

/******************************************************************************/
/*
  Purpose:

    DGESV_TEST uses DGESV to solve a general linear system A*x=b.

  Discussion:

    Solve

      1.0 -1.0  2.0 -1.0   x(0)     -8.0
      2.0 -2.0  3.0 -3.0 * x(1)  = -20.0
      1.0  1.0  1.0  0.0   x(2)     -2.0
      1.0 -1.0  4.0  3.0   x(3)      4.0

    The solution is ( -7, 3, 2, 2 ).

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    07 January 2014

  Author:

    John Burkardt
*/
{
/*
  The matrix must be stored in column-major form, so what you
  see below looks like the transpose of the actual matrix.
*/
  double A[4*4] = {
     1.0,  2.0,  1.0,  1.0,
    -1.0, -2.0,  1.0, -1.0,
     2.0,  3.0,  1.0,  4.0,
    -1.0, -3.0,  0.0,  3.0 };
  double B[4] = {
    -8.0, -20.0, -2.0, 4.0 };
  int i;
  static long int INFO;
  int info2;
  static long int IPIV[4];
  int j;
  long int LDA;
  long int LDB;
  long int N = 4;
  long int NRHS;

  printf ( "\n" );
  printf ( "DGESV_TEST\n" );
  printf ( "  Demonstrate the use of DGESV to solve a linear system\n" );
  printf ( "  using double precision real arithmetic.\n" );
/*
  Print the coefficient matrix.
*/
  r8mat_print ( N, N, A, "  Coefficient matrix A:" );
/*
  Print the right hand side.
*/
  r8vec_print ( N, B, "  Right hand side B:\n" );
/*
  Call DGESV to compute the solution.
*/
  NRHS = 1;
  LDA = N;
  LDB = N;

  dgesv_ ( &N, &NRHS, A, &LDA, IPIV, B, &LDB, &INFO );

  printf ( "\n" );
  printf ( "  Return value of error flag INFO = %d\n", ( int ) INFO );
/*
  Print the solution.
*/
  r8vec_print ( N, B, "  Computed solution X:\n" );

  return;
}
/******************************************************************************/

void dgesvd_test ( )

/******************************************************************************/
/*
  Purpose:

    DGESVD_TEST demonstrates DGESVD.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    07 January 2014

  Author:

    John Burkardt
*/
{
# define MVAL 4
# define NVAL 4
/*
  The entries of A are listed by columns, not rows!
*/
  long int LWORK = 201; 

  double a[MVAL*NVAL] = { 
    16.0,  5.0, 9.0,  4.0, 
     2.0, 11.0, 7.0, 14.0, 
     3.0, 10.0, 6.0, 15.0,
    13.0, 8.0, 12.0,  1.0 };
  int i;
  long int INFO;
  int j;
  char JOBU = 'A';
  char JOBVT = 'A';
  long int LDA = MVAL;
  long int LDU = MVAL;
  long int LDVT = NVAL;
  long int M = MVAL;
  long int N = NVAL;
  long int mn = min ( MVAL, NVAL );  
  long int MN = max ( MVAL, NVAL );
  double s[MVAL];
  double uu[MVAL*MVAL];
  double vt[NVAL*NVAL];
  double wk[LWORK];

  printf ( "\n" );
  printf ( "DGESVD_TEST\n" );
  printf ( "  Demonstrate the use of DGESVD to compute the\n" );
  printf ( "  singular value decomposition A = U * S * V',\n" );
  printf ( "  using double precision real arithmetic.\n" );
/*
  Print the coefficient matrix.
*/
  r8mat_print ( M, N, a, "  Coefficient matrix A:" );
/*
  Call DGESVD for singular value decomposition A = U * S * V'.
*/
  dgesvd_ ( &JOBU, &JOBVT, &M, &N, a, &LDA, s, uu, &LDU, vt, &LDVT, wk, 
    &LWORK, &INFO );
     
  printf ( "\n" );     
  printf ( "  Error flag INFO = %d\n", ( int ) INFO );          
/*
  Print the singular values.
*/
  r8vec_print ( M, s, "  Singular values:\n" );

  return;

# undef MVAL
# undef NVAL
}
/******************************************************************************/

void dgetrf_test ( )

/******************************************************************************/
/*
  Purpose:

    DGETRF_TEST demonstrates DGETRF and DGETRS.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    07 January 2014

  Author:

    John Burkardt
*/
{
  double A[4*4] = {
     1.0,  2.0,  1.0,  1.0,
    -1.0, -2.0,  1.0, -1.0,
     2.0,  3.0,  1.0,  4.0,
    -1.0, -3.0,  0.0,  3.0 };
  double B[4] = {
    -8.0, -20.0, -2.0, 4.0 };
  int i;
  static long int INFO;
  int info2;
  static long int IPIV[4];
  int j;
  long int LDA;
  long int LDB;
  long int N = 4;
  long int NRHS;
  char TRANS;

  printf ( "\n" );
  printf ( "DGETRF_TEST\n" );
  printf ( "  Demonstrate the use of:\n" );
  printf ( "  DGETRF to factor a general matrix A,\n" );
  printf ( "  DGETRS to solve A*x=b after A has been factored,\n" );
  printf ( "  using double precision real arithmetic.\n" );
/*
  Print the coefficient matrix.
*/
  r8mat_print ( N, N, A, "  Coefficient matrix A:" );
/*
  Call DGETRF to factor the matrix.
*/
  LDA = N;
  dgetrf_ ( &N, &N, A, &LDA, IPIV, &INFO );

  printf ( "\n" );
  printf ( "  Return value of DGETRF error flag INFO = %d\n", ( int ) INFO );
/*
  Print the right hand side.
*/
  r8vec_print ( N, B, "  Right hand side B:\n" );
/*
  Call DGETRS to solve the linear system A*x=b.
*/
  TRANS = 'N';
  NRHS = 1;
  LDB = N;

  dgetrs_ ( &TRANS, &N, &NRHS, A, &LDA, IPIV, B, &LDB, &INFO );

  printf ( "\n" );
  printf ( "  Return value of DGETRS error flag INFO = %d\n", ( int ) INFO );
/*
  Solution X is returned in B.
*/
  r8vec_print ( N, B, "  Computed solution X:\n" );

  return;
}
/******************************************************************************/

void dgetri_test ( )

/******************************************************************************/
/*
  Purpose:

    DGETRI_TEST tests DGETRF and DGETRI.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 January 2014

  Author:

    John Burkardt
*/
{
  double A[3*3] = {
    1.0, 4.0, 7.0,
    2.0, 5.0, 8.0,
    3.0, 6.0, 0.0 };
  long int INFO;
  long int IPIV[3];
  long int LDA;
  long int LWORK;
  long int N = 3;
  double WORK[3];

  printf ( "\n" );
  printf ( "DGETRI_TEST\n" );
  printf ( "  For a double precision real matrix (D)\n" );
  printf ( "  in general storage mode (GE):\n" );
  printf ( "\n" );
  printf ( "  DGETRF factors a general matrix;\n" );
  printf ( "  DGETRI computes the inverse.\n" );

  r8mat_print ( N, N, A, "  The matrix A:" );
/*
  Factor the matrix.
*/
  LDA = N;
  dgetrf_ ( &N, &N, A, &LDA, IPIV, &INFO );

  if ( ( int ) INFO != 0 )
  {
    printf ( "\n" );
    printf ( "  DGETRF returned INFO = %d\n", INFO );
    printf ( "  The matrix is numerically singular.\n" );
    return;
  }
/*
  Compute the inverse matrix.
*/
  LWORK = N;
  dgetri_ ( &N, A, &LDA, IPIV, WORK, &LWORK, &INFO );

  if ( ( int ) INFO != 0 )
  {
    printf ( "\n" );
    printf ( "  The inversion procedure failed!\n" );
    printf ( " '  INFO = %d\n", INFO );
    return;
  }
/*
  Print the inverse matrix.
*/
  r8mat_print ( N, N, A, "  The inverse matrix:" );

  return;
}
/******************************************************************************/

void dsyev_test ( )

/******************************************************************************/
/*
  Purpose:

    DSYEV_TEST tests DSYEV.

  Discussion:

    For some reason, you can't use "int" variables as arguments to CLAPACK
    functions; you have to use "integer" variables, which are apparently
    a name for the standard "long int" datatype.  If you also want to
    use int variables here and there, you may need to declare two versions
    of the same quantity.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 July 2013

  Author:

    John Burkardt
*/
{
  double *a;
  int info;
  long int INFO;
  char jobz;
  double *lambda;
  int lwork;
  long int LWORK;
  int n;
  long int N = 7;
  char uplo;
  double *work;

  printf ( "\n" );
  printf ( "DSYEV_TEST\n" );
  printf ( "  For a double precision real matrix (D)\n" );
  printf ( "  in symmetric storage mode (SY):\n" );
  printf ( "\n" );
  printf ( "  For a symmetric matrix in general storage,\n" );
  printf ( "  DSYEV computes eigenvalues and eigenvectors;\n" );
/*
  Set A.
*/
  n = ( int ) N;
  a = clement2 ( n );

  r8mat_print ( n, n, a, "  The matrix A:" );
/*
  Compute the eigenvalues and eigenvectors.
*/
  jobz = 'V';
  uplo = 'U';
  lambda = ( double * ) malloc ( N * sizeof ( double ) );
  LWORK = 3 * N - 1;
  work = ( double * ) malloc ( LWORK * sizeof ( double ) );

  dsyev_ ( &jobz, &uplo, &N, a, &N, lambda, work, &LWORK, &INFO );

  info = ( int ) INFO;
  if ( info != 0 )
  {
    printf ( "\n" );
    printf ( "  DSYEV returned nonzero INFO = %d\n", info );
  }
  else
  {
    r8vec_print ( n, lambda, "  The eigenvalues:" );
 
    if ( jobz == 'V' )
    {
      r8mat_print ( n, n, a, "  The eigenvector matrix:" );
    }
  }

  free ( a );
  free ( lambda );
  free ( work );

  return;
}
/******************************************************************************/

void zgesv_test ( )

/******************************************************************************/
/*
  Purpose:

    ZGESV_TEST demonstrates ZGESV.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    14 May 2013

  Author:

    John Burkardt
*/
{
# define NDIM 2

  doublecomplex *A;
  doublecomplex *B;
  int i;
  static long int INFO;
  int info2;
  static long int IPIV[NDIM];
  int j;
  long int LDA;
  long int LDB;
  long int N;
  long int NRHS;
  const double pi = 3.141592653589793;

  printf ( "\n" );
  printf ( "ZGESV_TEST\n" );
  printf ( "  Demonstrate the use of ZGESV to solve a linear system\n" );
  printf ( "  using double precision complex arithmetic.\n" );

  A = ( doublecomplex * ) malloc ( NDIM * NDIM * sizeof ( doublecomplex ) );
  B = ( doublecomplex * ) malloc ( NDIM * sizeof ( doublecomplex ) );

  N = NDIM;
  NRHS = 1;
  LDA = NDIM;
  LDB = NDIM;
/*
  Print the coefficient matrix.
*/
  printf ( "\n" );
  printf ( "  Coefficient matrix A:\n" );
  printf ( "\n" );

  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      A[i+N*j].r = cos ( pi * ( double ) ( i + 1 ) * 3.0 / 4.0 );  
      A[i+N*j].i = sin ( pi * ( double ) ( j + 1 ) / 5.0 );  
      printf ( "  %f + %fi", A[i+N*j].r, A[i+N*j].i );
    }
    printf ( "\n" );
  }
/*
  Print the right hand side.
*/
  printf ( "\n" );
  printf ( "  Right hand side B:\n" );
  printf ( "\n" );

  B[0].r = 1.0; 
  B[0].i = 1.0; 
  B[1].r = 2.0;
  B[1].i = 3.0;
  for ( i = 0; i < N; i++ )
  {
    printf ( "  %f + %fi\n", B[i].r, B[i].i );
  }
/*
  Call ZGESV to compute the solution.
*/
  info2 = zgesv_ ( &N, &NRHS, A, &LDA, IPIV, B, &LDB, &INFO );

  printf ( "\n" );
  printf ( "  Return value of error flag INFO = %d\n", ( int ) INFO );

  printf ( "\n" );
  printf ( "  Computed solution X:\n" );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    printf ( "  %f + %fi\n", B[i].r, B[i].i );
  }

  free ( A );
  free ( B );

  return;
# undef NDIM
}
/******************************************************************************/

void r8mat_print ( int m, int n, double a[], char *title )

/******************************************************************************/
/*
  Purpose:

    R8MAT_PRINT prints an R8MAT.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

    Entry A(I,J) is stored as A[I+J*M]

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    28 May 2008

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows in A.

    Input, int N, the number of columns in A.

    Input, double A[M*N], the M by N matrix.

    Input, char *TITLE, a title.
*/
{
  r8mat_print_some ( m, n, a, 1, 1, m, n, title );

  return;
}
/******************************************************************************/

void r8mat_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi,
  int jhi, char *title )

/******************************************************************************/
/*
  Purpose:

    R8MAT_PRINT_SOME prints some of an R8MAT.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 June 2013

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows of the matrix.
    M must be positive.

    Input, int N, the number of columns of the matrix.
    N must be positive.

    Input, double A[M*N], the matrix.

    Input, int ILO, JLO, IHI, JHI, designate the first row and
    column, and the last row and column to be printed.

    Input, char *TITLE, a title.
*/
{
# define INCX 5

  int i;
  int i2hi;
  int i2lo;
  int j;
  int j2hi;
  int j2lo;

  fprintf ( stdout, "\n" );
  fprintf ( stdout, "%s\n", title );

  if ( m <= 0 || n <= 0 )
  {
    fprintf ( stdout, "\n" );
    fprintf ( stdout, "  (None)\n" );
    return;
  }
/*
  Print the columns of the matrix, in strips of 5.
*/
  for ( j2lo = jlo; j2lo <= jhi; j2lo = j2lo + INCX )
  {
    j2hi = j2lo + INCX - 1;
    if ( n < j2hi )
    {
      j2hi = n;
    }
    if ( jhi < j2hi )
    {
      j2hi = jhi;
    }

    fprintf ( stdout, "\n" );
/*
  For each column J in the current range...

  Write the header.
*/
    fprintf ( stdout, "  Col:  ");
    for ( j = j2lo; j <= j2hi; j++ )
    {
      fprintf ( stdout, "  %7d     ", j - 1 );
    }
    fprintf ( stdout, "\n" );
    fprintf ( stdout, "  Row\n" );
    fprintf ( stdout, "\n" );
/*
  Determine the range of the rows in this strip.
*/
    if ( 1 < ilo )
    {
      i2lo = ilo;
    }
    else
    {
      i2lo = 1;
    }
    if ( m < ihi )
    {
      i2hi = m;
    }
    else
    {
      i2hi = ihi;
    }

    for ( i = i2lo; i <= i2hi; i++ )
    {
/*
  Print out (up to) 5 entries in row I, that lie in the current strip.
*/
      fprintf ( stdout, "%5d:", i - 1 );
      for ( j = j2lo; j <= j2hi; j++ )
      {
        fprintf ( stdout, "  %14f", a[i-1+(j-1)*m] );
      }
      fprintf ( stdout, "\n" );
    }
  }

  return;
# undef INCX
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
    fprintf ( stdout, "  %8d: %14f\n", i, a[i] );
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
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  printf ( "%s\n", time_buffer );

  return;
# undef TIME_SIZE
}
/******************************************************************************/

double *clement2 ( int n )

/******************************************************************************/
/*
  Purpose:

    CLEMENT2 returns the CLEMENT2 matrix.

  Formula:

    if ( J = I + 1 )
      A(I,J) = sqrt(I*(N-I))
    else if ( I = J + 1 )
      A(I,J) = sqrt(J*(N-J))
    else
      A(I,J) = 0

  Example:

    N = 5

       .    sqrt(4)    .       .       .
    sqrt(4)    .    sqrt(6)    .       .
       .    sqrt(6)    .    sqrt(6)    .
       .       .    sqrt(6)    .    sqrt(4)
       .       .       .    sqrt(4)    .

  Properties:

    A is tridiagonal.

    A is banded, with bandwidth 3.

    Because A is tridiagonal, it has property A (bipartite).

    A is symmetric: A' = A.

    Because A is symmetric, it is normal.

    Because A is normal, it is diagonalizable.

    A is persymmetric: A(I,J) = A(N+1-J,N+1-I).

    The diagonal of A is zero.

    A is singular if N is odd.

    About 64 percent of the entries of the inverse of A are zero.

    The eigenvalues are plus and minus the numbers

      N-1, N-3, N-5, ..., (1 or 0).

    If N is even,

      det ( A ) = (-1)^(N/2) * (N-1) * (N+1)^(N/2)

    and if N is odd,

      det ( A ) = 0

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 June 2008

  Author:

    John Burkardt

  Reference:

    Paul Clement,
    A class of triple-diagonal matrices for test purposes,
    SIAM Review,
    Volume 1, 1959, pages 50-52.

  Parameters:

    Input, int N, the order of the matrix.

    Output, double CLEMENT2[N*N], the matrix.
*/
{
  double *a;
  int i;
  int j;

  a = ( double * ) malloc ( n * n * sizeof ( double ) );

  for ( i = 1; i <= n; i++ )
  {
    for ( j = 1; j <= n; j++ )
    {
      if ( j == i + 1 )
      {
        a[i-1+(j-1)*n] = sqrt ( ( double ) ( i * ( n - i ) ) );
      }
      else if ( i == j + 1 )
      {
        a[i-1+(j-1)*n] = sqrt ( ( double ) ( j * ( n - j ) ) );
      }
      else
      {
        a[i-1+(j-1)*n] = 0.0;
      }
    }
  }
  return a;
}
