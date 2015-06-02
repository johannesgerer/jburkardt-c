# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <time.h>
# include <math.h>

# include "umfpack.h"

int main ( );
void cc_data_read ( char *prefix, int ncc, int n, int icc[], int ccc[], 
  double acc[] );
void cc_header_read ( char *prefix, int *ncc, int *n );
double *cc_mv ( int m, int n, int ncc, int icc[], int ccc[], double acc[], 
  double x[] );
void cc_print ( int m, int n, int ncc, int icc[], int ccc[], double acc[], 
  char *title );
void cc_print_some ( int i_min, int i_max, int j_min, int j_max, int ncc, 
  int n, int icc[], int ccc[], double acc[], char *title );
int file_row_count ( char *input_filename );
void i4vec_data_read ( char *input_filename, int n, int a[] );
void r8vec_data_read ( char *input_filename, int n, double x[] );
double r8vec_diff_norm ( int n, double a[], double b[] );
void r8vec_print ( int n, double a[], char *title );
double *r8vec_uniform_01_new ( int n, int *seed );
int s_len_trim ( char *s );
void timestamp ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for UMFPACK_WEST.

  Discussion:

    This program uses UMFPACK to solve a linear system A*X=B for which the
    matrix is stored, in compressed column (CC) format, in three files.

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
  double *acc;
  double *b;
  int *ccc;
  int i;
  int *icc;
  int m;
  int n;
  int ncc;
  double *null = ( double * ) NULL;
  void *Numeric;
  char prefix[] = "west";
  double r;
  int seed;
  int status;
  void *Symbolic;
  double *x1;
  double *x2;

  timestamp ( );
  printf ( "\n" );
  printf ( "UMFPACK_WEST:\n" );
  printf ( "  C version\n" );
  printf ( "  Use UMFPACK to solve the sparse linear system A*x=b.\n" );
  printf ( "  The matrix A is stored, in CC format, in 3 files.\n" );
/*
  Get the matrix size.
*/
  cc_header_read ( prefix, &ncc, &n );
  printf ( "\n" );
  printf ( "  Number of rows and columns = %d\n", n );
  printf ( "  Number of nonzeros NCC = %d\n", ncc );
/*
  Allocate space.
*/
  acc = ( double * ) malloc ( ncc * sizeof ( double ) );
  ccc = ( int * ) malloc ( ( n + 1 ) * sizeof ( int ) );
  icc = ( int * ) malloc ( ncc * sizeof ( int ) );
/*
  Read the matrix data.
*/
  cc_data_read ( prefix, ncc, n, icc, ccc, acc );
/*
  Print the matrix.
*/
  m = n;
  cc_print ( m, n, ncc, icc, ccc, acc, "  The CC matrix:" ); 
/*
  Set up the solution.
*/
  seed = 123456789;
  x1 = r8vec_uniform_01_new ( n, &seed );
/*
  Set the right hand side.
*/
  b = cc_mv ( m, n, ncc, icc, ccc, acc, x1 );
/*
  From the matrix data, create the symbolic factorization information.
*/
  status = umfpack_di_symbolic ( n, n, ccc, icc, acc, &Symbolic, null, null );
/*
  From the symbolic factorization information, carry out the numeric factorization.
*/
  status = umfpack_di_numeric ( ccc, icc, acc, Symbolic, &Numeric, null, null );
/*
  Free the symbolic factorization memory.
*/
  umfpack_di_free_symbolic ( &Symbolic );
/*
  Using the numeric factorization, solve the linear system.
*/
  x2 = ( double * ) malloc ( n * sizeof ( double ) );
  status = umfpack_di_solve ( UMFPACK_A, ccc, icc, acc, x2, b, Numeric, null, null );
/*
  Free the numeric factorization.
*/
  umfpack_di_free_numeric ( &Numeric );
/*
  Compute the error:
*/
  r = r8vec_diff_norm ( n, x1, x2 );
  printf ( "\n" );
  printf ( "  Residual: ||A*x-b|| = %g\n", r );
/*
  Free memory.
*/
  free ( acc );
  free ( b );
  free ( ccc );
  free ( icc );
  free ( x1 );
  free ( x2 );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "UMFPACK_WEST:\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void cc_data_read ( char *prefix, int ncc, int n, int icc[], int ccc[], 
  double acc[] )

/******************************************************************************/
/*
  Purpose:

    CC_DATA_READ reads data about a sparse matrix in CC format.

  Discussion:

    Three files are presumed to exist:
    * prefix_icc.txt contains NCC ICC values;
    * prefix_ccc.txt contains N+1 CCC values;
    * prefix_acc.txt contains NCC ACC values.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 July 2014

  Author:

    John Burkardt

  Parameters:

    Input, char *PREFIX, a common prefix for the filenames.

    Input, int NCC, the number of CC elements.

    Input, int N, the number of columns in the matrix.

    Output, int ICC[NCC], the CC rows.

    Output, int CCC[N+1], the compressed CC columns.

    Output, double ACC[NCC], the CC values.
*/
{
  char filename_acc[255];
  char filename_ccc[255];
  char filename_icc[255];

  strcpy ( filename_icc, prefix );
  strcat ( filename_icc, "_icc.txt" );
  i4vec_data_read ( filename_icc, ncc, icc );

  strcpy ( filename_ccc, prefix );
  strcat ( filename_ccc, "_ccc.txt" );
  i4vec_data_read ( filename_ccc, n + 1, ccc );

  strcpy ( filename_acc, prefix );
  strcat ( filename_acc, "_acc.txt" );
  r8vec_data_read ( filename_acc, ncc, acc );

  return;
}
/******************************************************************************/

void cc_header_read ( char *prefix, int *ncc, int *n )

/******************************************************************************/
/*
  Purpose:

    CC_HEADER_READ reads header information about a sparse matrix in CC format.

  Discussion:

    Three files are presumed to exist:
    * prefix_icc.txt contains NCC ICC values;
    * prefix_ccc.txt contains N+1 CCC values;
    * prefix_acc.txt contains NCC ACC values.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 July 2014

  Author:

    John Burkardt

  Parameters:

    Input, char *PREFIX, a common prefix for the filenames.

    Output, int *NCC, the number of CC elements.

    Output, int *N, the number of columns in the matrix.
*/
{
  char filename_ccc[255];
  char filename_icc[255];

  strcpy ( filename_icc, prefix );
  strcat ( filename_icc, "_icc.txt" );
  *ncc = file_row_count ( filename_icc );

  strcpy ( filename_ccc, prefix );
  strcat ( filename_ccc, "_ccc.txt" );
  *n = file_row_count ( filename_ccc ) - 1;

  return;
}
/******************************************************************************/

double *cc_mv ( int m, int n, int ncc, int icc[], int ccc[], double acc[], 
  double x[] )

/******************************************************************************/
/*
  Purpose:

    CC_MV multiplies a CC matrix by a vector

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    15 July 2014

  Author:

    John Burkardt

  Reference:

    Iain Duff, Roger Grimes, John Lewis,
    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
    October 1992

  Parameters:

    Input, int M, the number of rows.

    Input, int N, the number of columns.

    Input, int NCC, the number of CC values.

    Input, int RCC[NCC], the CC rows.

    Input, int CCC[N+1], the compressed CC columns

    Input, double ACC[NCC], the CC values.

    Input, double X[N], the vector to be multiplied.

    Output, double CC_MV[M], the product A*X.
*/
{
  double *b;
  int i;
  int j;
  int k;

  b = ( double * ) malloc ( m * sizeof ( double ) );

  for ( i = 0; i < m; i++ )
  {
    b[i] = 0.0;
  }

  for ( j = 0; j < n; j++ )
  {
    for ( k = ccc[j]; k < ccc[j+1]; k++ )
    {
      i = icc[k];
      b[i] = b[i] + acc[k] * x[j];
    }
  }

  return b;
}
/******************************************************************************/

void cc_print ( int m, int n, int ncc, int icc[], int ccc[], double acc[], 
  char *title )

/******************************************************************************/
/*
  Purpose:

    CC_PRINT prints a sparse matrix in CC format.

  Discussion:

    The index data in ICC and CCC is assumed to be 0-based.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 July 2014

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows in the matrix.

    Input, int N, the number of columns in the matrix.

    Input, int NCC, the number of CC elements.

    Input, int ICC[NCC], the CC rows.

    Input, int CCC[N+1], the compressed CC columns.

    Input, double ACC[NCC], the CC values.

    Input, char *TITLE, a title.
*/
{
  cc_print_some ( 0, m - 1, 0, n - 1, ncc, n, icc, ccc, acc, title );

  return;
}
/******************************************************************************/

void cc_print_some ( int i_min, int i_max, int j_min, int j_max, int ncc, 
  int n, int icc[], int ccc[], double acc[], char *title )

/******************************************************************************/
/*
  Purpose:

    CC_PRINT_SOME prints some of a sparse matrix in CC format.

  Discussion:

    The index data in ICC and CCC is assumed to be 0-based.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 July 2014

  Author:

    John Burkardt

  Parameters:

    Input, int I_MIN, IMAX, the first and last rows to print.

    Input, int J_MIN, J_MAX, the first and last columns 
    to print.

    Input, int NCC, the number of CC elements.

    Input, int N, the number of columns.

    Input, int ICC[NCC], the CC rows.

    Input, int CCC[N+1], the compressed CC columns.

    Input, double ACC[NCC], the CC values.

    Input, char *TITLE, a title.
*/
{
  int i;
  int j;
  int jnext;
  int k;

  printf ( "\n" );
  printf ( "%s\n", title );
  printf ( "     #     I     J         A\n" );
  printf ( "  ----  ----  ----  ----------------\n" );
  printf ( "\n" );

  j = 0;
  jnext = ccc[1];

  for ( k = 0; k < ncc; k++ )
  {
    i = icc[k];
    while ( jnext <= k )
    {
      j = j + 1;
      jnext = ccc[j+1];
    }
 
    if ( i_min <= i && i <= i_max &&
         j_min <= j && j <= j_max )
    {
      printf ( "  %4d  %4d  %4d  %16.8g\n", k, i, j, acc[k] );
    }
  }
  return;
}
/******************************************************************************/

int file_row_count ( char *input_filename )

/******************************************************************************/
/*
  Purpose:

    FILE_ROW_COUNT counts the number of row records in a file.

  Discussion:

    It does not count lines that are blank, or that begin with a
    comment symbol '#'.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    13 June 2003

  Author:

    John Burkardt

  Parameters:

    Input, char *INPUT_FILENAME, the name of the input file.

    Output, int FILE_ROW_COUNT, the number of rows found.
*/
{
# define LINE_MX 255

  int bad_num;
  int comment_num;
  char *error;
  FILE *input;
  int i;
  char line[LINE_MX];
  int record_num;
  int row_num;

  row_num = 0;
  comment_num = 0;
  record_num = 0;
  bad_num = 0;

  input = fopen ( input_filename, "r" );

  if ( !input )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "FILE_ROW_COUNT - Fatal error!\n" );
    fprintf ( stderr, "  Could not open the input file: \"%s\"\n",
      input_filename );
    exit ( 1 );
  }

  for ( ; ; )
  {
    error = fgets ( line, LINE_MX, input );

    if ( !error )
    {
      break;
    }

    record_num = record_num + 1;

    if ( line[0] == '#' )
    {
      comment_num = comment_num + 1;
      continue;
    }

    if ( s_len_trim ( line ) == 0 )
    {
      comment_num = comment_num + 1;
      continue;
    }

    row_num = row_num + 1;
  }

  fclose ( input );

  return row_num;

# undef LINE_MX
}
/******************************************************************************/

void i4vec_data_read ( char *input_filename, int n, int a[] )

/******************************************************************************/
/*
  Purpose:

    I4VEC_DATA_READ reads the data from an I4VEC file.

  Discussion:

    An I4VEC is an array of I4's.

    The file is assumed to contain one record per line.

    Records beginning with the '#' character are comments, and are ignored.
    Blank lines are also ignored.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 July 2014

  Author:

    John Burkardt

  Parameters:

    Input, char *INPUT_FILENAME, the name of the input file.

    Input, int N, the number of points.

    Output, int A[N], the data.
*/
{
# define LINE_MX 255

  char *got_string;
  FILE *input;
  int j;
  int l;
  char line[255];;

  input = fopen ( input_filename, "r" );

  if ( !input )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "I4VEC_DATA_READ - Fatal error!\n" );
    fprintf ( stderr, 
      "  Could not open the input file: \"%s\"\n", input_filename );
    exit ( 1 );
  }

  j = 0;

  while ( j < n )
  {
    got_string = fgets ( line, LINE_MX, input );

    if ( !got_string )
    {
      break;
    }

    if ( line[0] == '#' || s_len_trim ( line ) == 0 )
    {
      continue;
    }

    a[j] = atoi ( line );
    j = j + 1;
  }

  fclose ( input );

  return;

# undef LINE_MX
}
/******************************************************************************/

void r8vec_data_read ( char *input_filename, int n, double x[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_DATA_READ reads the data from an R8VEC file.

  Discussion:

    An R8VEC is a vector of R8's.

    The file is assumed to contain one record per line.

    Records beginning with the '#' character are comments, and are ignored.
    Blank lines are also ignored.

    There are assumed to be exactly (or at least) N such records.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 July 2014

  Author:

    John Burkardt

  Parameters:

    Input, char *INPUT_FILENAME, the name of the input file.

    Input, int N, the number of points.

    Output, double X[N], the data.
*/
{
# define LINE_MX 255

  int error;
  char *got_string;
  FILE *input;
  int i;
  int j;
  int lchar;
  char line[255];

  input = fopen ( input_filename, "r" );

  if ( !input )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8VEC_DATA_READ - Fatal error!\n" );
    fprintf ( stderr, "  Could not open the input file: \"%s\"\n", input_filename );
    exit ( 1 );
  }

  j = 0;

  while ( j < n )
  {
    got_string = fgets ( line, LINE_MX, input );

    if ( !got_string )
    {
      break;
    }

    if ( line[0] == '#' || s_len_trim ( line ) == 0 )
    {
      continue;
    }

    x[j] = atof ( line );

    j = j + 1;
  }

  fclose ( input );

  return;

# undef LINE_MX
}
/******************************************************************************/

double r8vec_diff_norm ( int n, double a[], double b[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_DIFF_NORM returns the L2 norm of the difference of R8VEC's.

  Discussion:

    An R8VEC is a vector of R8's.

    The vector L2 norm is defined as:

      R8VEC_NORM_L2 = sqrt ( sum ( 1 <= I <= N ) A(I)^2 ).

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    24 June 2011

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in A.

    Input, double A[N], B[N], the vectors.

    Output, double R8VEC_DIFF_NORM, the L2 norm of A - B.
*/
{
  int i;
  double value;

  value = 0.0;

  for ( i = 0; i < n; i++ )
  {
    value = value + ( a[i] - b[i] ) * ( a[i] - b[i] );
  }
  value = sqrt ( value );

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

double *r8vec_uniform_01_new ( int n, int *seed )

/******************************************************************************/
/*
  Purpose:

    R8VEC_UNIFORM_01_NEW returns a unit pseudorandom R8VEC.

  Discussion:

    This routine implements the recursion

      seed = 16807 * seed mod ( 2^31 - 1 )
      unif = seed / ( 2^31 - 1 )

    The integer arithmetic never requires more than 32 bits,
    including a sign bit.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 August 2004

  Author:

    John Burkardt

  Reference:

    Paul Bratley, Bennett Fox, Linus Schrage,
    A Guide to Simulation,
    Second Edition,
    Springer, 1987,
    ISBN: 0387964673,
    LC: QA76.9.C65.B73.

    Bennett Fox,
    Algorithm 647:
    Implementation and Relative Efficiency of Quasirandom
    Sequence Generators,
    ACM Transactions on Mathematical Software,
    Volume 12, Number 4, December 1986, pages 362-376.

    Pierre L'Ecuyer,
    Random Number Generation,
    in Handbook of Simulation,
    edited by Jerry Banks,
    Wiley, 1998,
    ISBN: 0471134031,
    LC: T57.62.H37.

    Peter Lewis, Allen Goodman, James Miller,
    A Pseudo-Random Number Generator for the System/360,
    IBM Systems Journal,
    Volume 8, Number 2, 1969, pages 136-143.

  Parameters:

    Input, int N, the number of entries in the vector.

    Input/output, int *SEED, a seed for the random number generator.

    Output, double R8VEC_UNIFORM_01_NEW[N], the vector of pseudorandom values.
*/
{
  int i;
  const int i4_huge = 2147483647;
  int k;
  double *r;

  if ( *seed == 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8VEC_UNIFORM_01_NEW - Fatal error!\n" );
    fprintf ( stderr, "  Input value of SEED = 0.\n" );
    exit ( 1 );
  }

  r = ( double * ) malloc ( n * sizeof ( double ) );

  for ( i = 0; i < n; i++ )
  {
    k = *seed / 127773;

    *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

    if ( *seed < 0 )
    {
      *seed = *seed + i4_huge;
    }

    r[i] = ( double ) ( *seed ) * 4.656612875E-10;
  }

  return r;
}
/******************************************************************************/

int s_len_trim ( char *s )

/******************************************************************************/
/*
  Purpose:

    S_LEN_TRIM returns the length of a string to the last nonblank.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    26 April 2003

  Author:

    John Burkardt

  Parameters:

    Input, char *S, a pointer to a string.

    Output, int S_LEN_TRIM, the length of the string to the last nonblank.
    If S_LEN_TRIM is 0, then the string is entirely blank.
*/
{
  int n;
  char *t;

  n = strlen ( s );
  t = s + strlen ( s ) - 1;

  while ( 0 < n ) 
  {
    if ( *t != ' ' )
    {
      return n;
    }
    t--;
    n--;
  }

  return n;
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
