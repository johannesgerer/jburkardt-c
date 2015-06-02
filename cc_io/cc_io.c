# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <string.h>

# include "cc_io.h"

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

    Input, int ICC[NCC], the CC rows.

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

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    23 July 2014

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
  int i;
  int j;
  int k;

  printf ( "\n" );
  printf ( "%s\n", title );
  printf ( "     #     I     J           A\n" );
  printf ( "  ----  ----  ----  ----------------\n" );
  printf ( "\n" );

  if ( ccc[0] == 0 )
  {
    j = 0;
    for ( k = 0; k < ncc; k++ )
    {
      i = icc[k];
      while ( ccc[j+1] <= k )
      {
        j = j + 1;
      }
      printf ( "  %4d  %4d  %4d  %16.8g\n", k, i, j, acc[k] );
    }
  }
  else
  {
    j = 1;
    for ( k = 0; k < ncc; k++ )
    {
      i = icc[k];
      while ( ccc[j] <= k + 1 )
      {
        j = j + 1;
      }
      printf ( "  %4d  %4d  %4d  %16.8g\n", k + 1, i, j, acc[k] );
    }
  }

  return;
}
/******************************************************************************/

void cc_print_some ( int i_min, int i_max, int j_min, int j_max, int ncc, 
  int n, int icc[], int ccc[], double acc[], char *title )

/******************************************************************************/
/*
  Purpose:

    CC_PRINT_SOME prints some of a sparse matrix in CC format.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    23 July 2014

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
  printf ( "     #     I     J           A\n" );
  printf ( "  ----  ----  ----  ----------------\n" );
  printf ( "\n" );

  if ( ccc[0] == 0 )
  {
    j = 0;
    for ( k = 0; k < ncc; k++ )
    {
      i = icc[k];
      while ( ccc[j+1] <= k )
      {
        j = j + 1;
      }
      if ( i_min <= i && i <= i_max &&
           j_min <= j && j <= j_max )
      {
        printf ( "  %4d  %4d  %4d  %16.8g\n", k, i, j, acc[k] );
      }
    }
  }
  else
  {
    j = 1;
    for ( k = 0; k < ncc; k++ )
    {
      i = icc[k];
      while ( ccc[j] <= k + 1 )
      {
        j = j + 1;
      }
      if ( i_min <= i && i <= i_max &&
           j_min <= j && j <= j_max )
      {
        printf ( "  %4d  %4d  %4d  %16.8g\n", k + 1, i, j, acc[k] );
      }
    }
  }

  return;
}
/******************************************************************************/

void cc_write ( char *prefix, int ncc, int n, int icc[], int ccc[], 
  double acc[] )

/******************************************************************************/
/*
  Purpose:

    CC_WRITE writes a sparse matrix in CC format to 3 files.

  Discussion:

    Three files will be created:
    * prefix_icc.txt contains NCC ICC values;
    * prefix_ccc.txt contains N+1 CCC values;
    * prefix_acc.txt contains NCC ACC values.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 July 2014

  Author:

    John Burkardt

  Parameters:

    Input, char *PREFIX, a common prefix for the filenames.

    Input, int NCC, the number of CC elements.

    Input, int N, the number of columns in the matrix.

    Input, int ICC[NCC], the CC rows.

    Input, int CCC[N+1], the compressed CC columns.

    Input, double ACC[NCC], the CC values.
*/
{
  char filename_acc[255];
  char filename_ccc[255];
  char filename_icc[255];

  strcpy ( filename_icc, prefix );
  strcat ( filename_icc, "_icc.txt" );
  i4vec_write ( filename_icc, ncc, icc );

  strcpy ( filename_ccc, prefix );
  strcat ( filename_ccc, "_ccc.txt" );
  i4vec_write ( filename_ccc, n + 1, ccc );

  strcpy ( filename_acc, prefix );
  strcat ( filename_acc, "_acc.txt" );
  r8vec_write ( filename_acc, ncc, acc );

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
# define LINE_MAX 255

  int bad_num;
  int comment_num;
  char *error;
  FILE *input;
  int i;
  char line[LINE_MAX];
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
    error = fgets ( line, LINE_MAX, input );

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

# undef LINE_MAX
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
# define LINE_MAX 255

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
    got_string = fgets ( line, LINE_MAX, input );

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

# undef LINE_MAX
}
/******************************************************************************/

void i4vec_dec ( int n, int a[] )

/******************************************************************************/
/*
  Purpose:

    I4VEC_DEC decrements an I4VEC.

  Discussion:

    An I4VEC is a vector of I4's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 July 2014

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the vectors.

    Input/output, int A[N], the vector to be decremented.
*/
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    a[i] = a[i] - 1;
  }
  return;
}
/******************************************************************************/

void i4vec_inc ( int n, int a[] )

/******************************************************************************/
/*
  Purpose:

    I4VEC_INC increments an I4VEC.

  Discussion:

    An I4VEC is a vector of I4's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 July 2014

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the vectors.

    Input/output, int A[N], the vector to be incremented.
*/
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    a[i] = a[i] + 1;
  }
  return;
}
/******************************************************************************/

void i4vec_write ( char *output_filename, int n, int table[] )

/******************************************************************************/
/*
  Purpose:

    I4VEC_WRITE writes an I4VEC to a file.

  Discussion:

    An I4VEC is a vector of I4's.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 July 2014

  Author:

    John Burkardt

  Parameters:

    Input, char *OUTPUT_FILENAME, the output filename.

    Input, int N, the number of points.

    Input, int TABLE[N], the data.
*/
{
  int j;
  FILE *output;
/*
  Open the file.
*/
  output = fopen ( output_filename, "wt" );

  if ( !output )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "I4VEC_WRITE - Fatal error!\n" );
    fprintf ( stderr, "  Could not open the output file.\n" );
    exit ( 1 );
  }
/*
  Write the data.
*/
  for ( j = 0; j < n; j++ )
  {
    fprintf ( output, "%d\n", table[j] );
  }
/*
  Close the file.
*/
  fclose ( output );

  return;
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
# define LINE_MAX 255

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
    got_string = fgets ( line, LINE_MAX, input );

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

# undef LINE_MAX
}
/******************************************************************************/

void r8vec_write ( char *output_filename, int n, double x[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_WRITE writes an R8VEC file.

  Discussion:

    An R8VEC is a vector of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 July 2014

  Author:

    John Burkardt

  Parameters:

    Input, char *OUTPUT_FILENAME, the output filename.

    Input, int N, the number of points.

    Input, double X[N], the data.
*/
{
  int j;
  FILE *output;
/*
  Open the file.
*/
  output = fopen ( output_filename, "wt" );

  if ( !output )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8VEC_WRITE - Fatal error!\n" );
    fprintf ( stderr, "  Could not open the output file.\n" );
    exit ( 1 );
  }
/*
  Write the data.
*/
  for ( j = 0; j < n; j++ )
  {
    fprintf ( output, "%g\n", x[j] );
  }
/*
  Close the file.
*/
  fclose ( output );

  return;
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
