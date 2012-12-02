# include <stdlib.h>
# include <stdio.h>
# include <time.h>
# include <string.h>

# include "latinize.h"

/******************************************************************************/

char ch_cap ( char ch )

/******************************************************************************/
/*
  Purpose:

    CH_CAP capitalizes a single character.

  Discussion:

    This routine should be equivalent to the library "toupper" function.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    19 July 1998

  Author:

    John Burkardt

  Parameters:

    Input, char CH, the character to capitalize.

    Output, char CH_CAP, the capitalized character.
*/
{
  if ( 97 <= ch && ch <= 122 ) 
  {
    ch = ch - 32;
  }   

  return ch;
}
/******************************************************************************/

int ch_eqi ( char ch1, char ch2 )

/******************************************************************************/
/*
  Purpose:

    CH_EQI is TRUE (1) if two characters are equal, disregarding case.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    13 June 2003

  Author:

    John Burkardt

  Parameters:

    Input, char CH1, CH2, the characters to compare.

    Output, int CH_EQI, is TRUE (1) if the two characters are equal,
    disregarding case and FALSE (0) otherwise.
*/
{
  int value;

  if ( 97 <= ch1 && ch1 <= 122 ) 
  {
    ch1 = ch1 - 32;
  } 
  if ( 97 <= ch2 && ch2 <= 122 ) 
  {
    ch2 = ch2 - 32;
  }     
  if ( ch1 == ch2 )
  {
    value = 1;
  }
  else
  {
    value = 0;
  }
  return value;
}
/******************************************************************************/

int ch_to_digit ( char ch )

/******************************************************************************/
/*
  Purpose:

    CH_TO_DIGIT returns the integer value of a base 10 digit.

  Example:

     CH  DIGIT
    ---  -----
    '0'    0
    '1'    1
    ...  ...
    '9'    9
    ' '    0
    'X'   -1

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    13 June 2003

  Author:

    John Burkardt

  Parameters:

    Input, char CH, the decimal digit, '0' through '9' or blank are legal.

    Output, int CH_TO_DIGIT, the corresponding integer value.  If the 
    character was 'illegal', then DIGIT is -1.
*/
{
  int digit;

  if ( '0' <= ch && ch <= '9' )
  {
    digit = ch - '0';
  }
  else if ( ch == ' ' )
  {
    digit = 0;
  }
  else
  {
    digit = -1;
  }

  return digit;
}
/******************************************************************************/

int file_column_count ( char *input_filename )

/******************************************************************************/
/*
  Purpose:

    FILE_COLUMN_COUNT counts the number of columns in the first line of a file.

  Discussion:

    The file is assumed to be a simple text file.

    Most lines of the file is presumed to consist of COLUMN_NUM words, separated
    by spaces.  There may also be some blank lines, and some comment lines,
    which have a "#" in column 1.

    The routine tries to find the first non-comment non-blank line and
    counts the number of words in that line.

    If all lines are blanks or comments, it goes back and tries to analyze
    a comment line.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    13 June 2003

  Author:

    John Burkardt

  Parameters:

    Input, char *INPUT_FILENAME, the name of the file.

    Output, int FILE_COLUMN_COUNT, the number of columns assumed 
    to be in the file.
*/
{
# define LINE_MAX 256

  int column_num;
  char *error;
  FILE *input;
  int got_one;
  char line[LINE_MAX];
/*
  Open the file.
*/
  input = fopen ( input_filename, "r" );

  if ( !input )
  {
    column_num = -1;
    printf ( "\n" );
    printf ( "FILE_COLUMN_COUNT - Fatal error!\n" );
    printf ( "  Could not open the input file: \"%s\"\n", input_filename );
    return column_num;
  }
/*
  Read one line, but skip blank lines and comment lines.
*/
  got_one = 0;

  for ( ; ; )
  {
    error = fgets ( line, LINE_MAX, input );

    if ( !error )
    {
      break;
    }

    if ( s_len_trim ( line ) == 0 )
    {
      continue;
    }

    if ( line[0] == '#' )
    {
      continue;
    }

    got_one = 1;
    break;

  }

  if ( got_one == 0 )
  {
    fclose ( input );

    input = fopen ( input_filename, "r" );

    for ( ; ; )
    {
      error = fgets ( line, LINE_MAX, input );

      if ( !error )
      {
        break;
      }

      if ( s_len_trim ( line ) == 0 )
      {
        continue;
      }

      got_one = 1;
      break;
    }
  }

  fclose ( input );

  if ( got_one == 0 )
  {
    printf ( "\n" );
    printf ( "FILE_COLUMN_COUNT - Warning!\n" );
    printf ( "  The file does not seem to contain any data.\n" );
    return -1;
  }

  column_num = s_word_count ( line );

  return column_num;

# undef LINE_MAX
}
//****************************************************************************80

char *file_name_ext_swap ( char *filename, char *ext )

//****************************************************************************80
//
//  Purpose:
//
//    FILE_NAME_EXT_SWAP replaces the current "extension" of a file name.
//
//  Discussion:
//
//    The "extension" of a file name is the string of characters
//    that appears after the LAST period in the name.  A file
//    with no period, or with a period as the last character
//    in the name, has a "null" extension.
//
//  Example:
//
//          Input           Output
//    ================     ==================
//    FILENAME     EXT     FILE_NAME_EXT_SWAP
//
//    bob.for      obj     bob.obj
//    bob.bob.bob  txt     bob.bob.txt
//    bob          yak     bob.yak
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 February 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char *FILENAME, a file name.
//
//    Input, char *EXT, the extension to be added to the file name.
//
//    Output, char *FILE_NAME_EXT_SWAP, the file name with the new extension.
//
{
  char *filename2;
  int i;
  int l1;
  int l2;

  l1 = s_index_last_c ( filename, '.' );
  l2 = s_len_trim ( ext );

  filename2 = ( char * ) malloc ( ( l1 + l2 + 1 ) * sizeof ( char ) );

  if ( l1 == - 1 ) 
  {
    l1 = s_len_trim ( filename );
    filename2 = ( char * ) malloc ( ( l1 + 1 + l2 + 1 ) * sizeof ( char ) );
    for ( i = 0; i < l1; i++ )
    {
      filename2[i] = filename[i];
    }
    filename2[l1] = '.';
    for ( i = 0; i < l2; i++ )
    {
      filename2[l1+1+i] = ext[i];
    }
    filename2[l1+1+l2] = '\0';
  }
  else
  {
    filename2 = ( char * ) malloc ( ( l1 + l2 + 1 ) * sizeof ( char ) );
    for ( i = 0; i < l1; i++ )
    {
      filename2[i] = filename[i];
    }
    for ( i = 0; i < l2; i++ )
    {
      filename2[l1+i] = ext[i];
    }
    filename2[l1+l2] = '\0';
  }

  return filename2;
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
# define LINE_MAX 256

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
    printf ( "\n" );
    printf ( "FILE_ROW_COUNT - Fatal error!\n" );
    printf ( "  Could not open the input file: \"%s\"\n", input_filename );
    return (-1);
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

int i4_max ( int i1, int i2 )

/******************************************************************************/
/*
  Purpose:

    I4_MAX returns the maximum of two I4's.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    29 August 2006

  Author:

    John Burkardt

  Parameters:

    Input, int I1, I2, are two integers to be compared.

    Output, int I4_MAX, the larger of I1 and I2.
*/
{
  int value;

  if ( i2 < i1 )
  {
    value = i1;
  }
  else
  {
    value = i2;
  }
  return value;
}
/******************************************************************************/

int i4_min ( int i1, int i2 )

/******************************************************************************/
/*
  Purpose:

    I4_MIN returns the smaller of two I4's.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    29 August 2006

  Author:

    John Burkardt

  Parameters:

    Input, int I1, I2, two integers to be compared.

    Output, int I4_MIN, the smaller of I1 and I2.
*/
{
  int value;

  if ( i1 < i2 )
  {
    value = i1;
  }
  else
  {
    value = i2;
  }
  return value;
}
/******************************************************************************/

double *r8mat_data_read ( char *input_filename, int m, int n )

/******************************************************************************/
/*
  Purpose:

    R8MAT_DATA_READ reads the data from an R8MAT file.

  Discussion:

    An R8MAT is an array of R8's.

    The file is assumed to contain one record per line.

    Records beginning with the '#' character are comments, and are ignored.
    Blank lines are also ignored.

    Each line that is not ignored is assumed to contain exactly (or at least)
    M real numbers, representing the coordinates of a point.

    There are assumed to be exactly (or at least) N such records.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    27 January 2005

  Author:

    John Burkardt

  Parameters:

    Input, char *INPUT_FILENAME, the name of the input file.

    Input, int M, the number of spatial dimensions.

    Input, int N, the number of points.  The program
    will stop reading data once N values have been read.

    Output, double R8MAT_DATA_READ[M*N], the data.
*/
{
# define LINE_MAX 255

  int error;
  char *got_string;
  FILE *input;
  int i;
  int j;
  char line[255];
  double *table;
  double *x;

  input = fopen ( input_filename, "r" );

  if ( !input )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8MAT_DATA_READ - Fatal error!\n" );
    fprintf ( stderr, "  Could not open the input file: \"%s\"\n", input_filename );
    exit ( 1 );
  }

  table = ( double * ) malloc ( m * n * sizeof ( double ) );

  x = ( double * ) malloc ( m * sizeof ( double ) );

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

    error = s_to_r8vec ( line, m, x );

    if ( error == 1 )
    {
      continue;
    }

    for ( i = 0; i < m; i++ )
    {
      table[i+j*m] = x[i];
    }
    j = j + 1;

  }

  fclose ( input );

  free ( x );

  return table;

# undef LINE_MAX
}
/******************************************************************************/
 
void r8mat_header_read ( char *input_filename, int *m, int *n )
 
/******************************************************************************/
/*
  Purpose:

    R8MAT_HEADER_READ reads the header from an R8MAT file.

  Discussion:

    An R8MAT is an array of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    04 June 2004

  Author:

    John Burkardt

  Parameters:

    Input, char *INPUT_FILENAME, the name of the input file.

    Output, int *M, the number of spatial dimensions.

    Output, int *N, the number of points.
*/
{
  *m = file_column_count ( input_filename );

  if ( *m <= 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8MAT_HEADER_READ - Fatal error!\n" );
    fprintf ( stderr, "  FILE_COLUMN_COUNT failed.\n" );
    exit ( 1 );
  }

  *n = file_row_count ( input_filename );

  if ( *n <= 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8MAT_HEADER_READ - Fatal error!\n" );
    fprintf ( stderr, "  FILE_ROW_COUNT failed.\n" );
    exit ( 1 );
  }

  return;
}
/*******************************************************************************/

void r8mat_latinize ( int m, int n, double table[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_LATINIZE "latinizes" a real table.

  Discussion:

    On output, each row of the table will have the properties that:
    1) the minimum and maximum row values are the same as on input;
    2) the row contains N evenly spaced values between the
       minimum and maximum;
    3) in each row, the elements retain their ordering.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    04 February 2012

  Author:

    John Burkardt

  Parameters:

    Input, int M, the spatial dimension.

    Input, int N, the number of columns.

    Input/output, double TABLE[M*N].  On input, the dataset to
    be "latinized".  On output, the latinized dataset.
*/
{
  int i;
  int *indx;
  int j;
  double *v;
  double v_max;
  double v_min;

  v = ( double * ) malloc ( n * sizeof ( double ) );

  for ( i = 0; i < m; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      v[j] = table[i+j*m];
    }
    v_min = r8vec_min ( n, v );
    v_max = r8vec_max ( n, v );
    indx = r8vec_sort_heap_index_a_new ( n, v );

    for ( j = 0; j < n; j++ )
    {
      table[i+indx[j]*m] =  ( ( double ) ( n - j - 1 ) * v_min   
                            + ( double ) (     j     ) * v_max ) 
                            / ( double ) ( n     - 1 );

    }
    free ( indx );
  }
  free ( v );

  return;
}
/******************************************************************************/

void r8mat_print ( int m, int n, double a[], char *title )

/******************************************************************************/
/*
  Purpose:

    R8MAT_PRINT prints an R8MAT, with an optional title.

  Discussion:

    An R8MAT is an array of R8's.

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

    An R8MAT is an array of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    28 May 2008

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

  printf ( "\n" );
  printf ( "%s\n", title );
/*
  Print the columns of the matrix, in strips of 5.
*/
  for ( j2lo = jlo; j2lo <= jhi; j2lo = j2lo + INCX )
  {
    j2hi = j2lo + INCX - 1;
    j2hi = i4_min ( j2hi, n );
    j2hi = i4_min ( j2hi, jhi );

    printf ( "\n" );
/*
  For each column J in the current range...

  Write the header.
*/
    printf ( "  Col:  ");
    for ( j = j2lo; j <= j2hi; j++ )
    {
      printf ( "  %7d     ", j );
    }
    printf ( "\n" );
    printf ( "  Row\n" );
    printf ( "\n" );
/*
  Determine the range of the rows in this strip.
*/
    i2lo = i4_max ( ilo, 1 );
    i2hi = i4_min ( ihi, m );

    for ( i = i2lo; i <= i2hi; i++ )
    {
/*
  Print out (up to) 5 entries in row I, that lie in the current strip.
*/
      printf ( "%5d", i );
      for ( j = j2lo; j <= j2hi; j++ )
      {
        printf ( "  %12f", a[i-1+(j-1)*m] );
      }
      printf ( "\n" );
    }
  }

  return;
# undef INCX
}
/******************************************************************************/

void r8mat_transpose_print ( int m, int n, double a[], char *title )

/******************************************************************************/
/*
  Purpose:

    R8MAT_TRANSPOSE_PRINT prints an R8MAT, transposed.

  Discussion:

    An R8MAT is an array of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    28 May 2008

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns.

    Input, double A[M*N], an M by N matrix to be printed.

    Input, char *TITLE, a title.
*/
{
  r8mat_transpose_print_some ( m, n, a, 1, 1, m, n, title );

  return;
}
/******************************************************************************/

void r8mat_transpose_print_some ( int m, int n, double a[], int ilo, int jlo, 
  int ihi, int jhi, char *title )

/******************************************************************************/
/*
  Purpose:

    R8MAT_TRANSPOSE_PRINT_SOME prints some of an R8MAT, transposed.

  Discussion:

    An R8MAT is an array of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    28 May 2008

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns.

    Input, double A[M*N], an M by N matrix to be printed.

    Input, int ILO, JLO, the first row and column to print.

    Input, int IHI, JHI, the last row and column to print.

    Input, char *TITLE, a title.
*/
{
# define INCX 5

  int i;
  int i2;
  int i2hi;
  int i2lo;
  int inc;
  int j;
  int j2hi;
  int j2lo;

  printf ( "\n" );
  printf ( "%s\n", title );
 
  for ( i2lo = i4_max ( ilo, 1 ); i2lo <= i4_min ( ihi, m ); i2lo = i2lo + INCX )
  {
    i2hi = i2lo + INCX - 1;
    i2hi = i4_min ( i2hi, m );
    i2hi = i4_min ( i2hi, ihi );

    inc = i2hi + 1 - i2lo;

    printf ( "\n" );
    printf ( "  Row:" );
    for ( i = i2lo; i <= i2hi; i++ )
    {
      printf ( "  %7d     ", i );
    }
    printf ( "\n" );
    printf ( "  Col\n" );
    printf ( "\n" );

    j2lo = i4_max ( jlo, 1 );
    j2hi = i4_min ( jhi, n );

    for ( j = j2lo; j <= j2hi; j++ )
    {
      printf ( "%5d", j );
      for ( i2 = 1; i2 <= inc; i2++ )
      {
        i = i2lo - 1 + i2;
        printf ( "  %14f", a[(i-1)+(j-1)*m] );
      }
      printf ( "\n" );
    }
  }
  printf ( "\n" );

  return;
# undef INCX
}
/******************************************************************************/

void r8mat_write ( char *output_filename, int m, int n, double table[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_WRITE writes an R8MAT file.

  Discussion:

    An R8MAT is an array of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    01 June 2009

  Author:

    John Burkardt

  Parameters:

    Input, char *OUTPUT_FILENAME, the output filename.

    Input, int M, the spatial dimension.

    Input, int N, the number of points.

    Input, double TABLE[M*N], the data.
*/
{
  int i;
  int j;
  FILE *output;
/*
  Open the file.
*/
  output = fopen ( output_filename, "wt" );

  if ( !output )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8MAT_WRITE - Fatal error!\n" );
    fprintf ( stderr, "  Could not open the output file.\n" );
    exit ( 1 );
  }
/*
  Write the data.
*/
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      fprintf ( output, "  %24.16g", table[i+j*m] );
    }
    fprintf ( output, "\n" );
  }
/*
  Close the file.
*/
  fclose ( output );

  return;
}
/******************************************************************************/

double r8vec_max ( int n, double r8vec[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_MAX returns the value of the maximum element in a R8VEC.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 May 2006

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the array.

    Input, double R8VEC[N], a pointer to the first entry of the array.

    Output, double R8VEC_MAX, the value of the maximum element.  This
    is set to 0.0 if N <= 0.
*/
{
  int i;
  double *r8vec_pointer;
  double value;

  if ( n <= 0 )
  {
    value = 0.0;
    return value;
  }

  value = r8vec[0];

  for ( i = 1; i < n; i++ )
  {
    if ( value < r8vec[i] )
    {
      value = r8vec[i];
    }
  }
  return value;
}
/******************************************************************************/

double r8vec_min ( int n, double r8vec[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_MIN returns the value of the minimum element in a R8VEC.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 May 2006

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the array.

    Input, double R8VEC[N], the array to be checked.

    Output, double R8VEC_MIN, the value of the minimum element.
*/
{
  int i;
  double value;

  if ( n <= 0 )
  {
    value = 0.0;
    return;
  }

  value = r8vec[0];

  for ( i = 1; i < n; i++ )
  {
    if ( r8vec[i] < value )
    {
      value = r8vec[i];
    }
  }
  return value;
}
/******************************************************************************/

int *r8vec_sort_heap_index_a_new ( int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_SORT_HEAP_INDEX_A_NEW does an indexed heap ascending sort of an R8VEC

  Discussion:

    An R8VEC is a vector of R8's.

    The sorting is not actually carried out.  Rather an index array is
    created which defines the sorting.  This array may be used to sort
    or index the array, or to sort or index related arrays keyed on the
    original array.

    Once the index array is computed, the sorting can be carried out
    "implicitly:

      a(indx(*))

    or explicitly, by the call

      r8vec_permute ( n, indx, 0, a )

    after which a(*) is sorted.

    Note that the index vector is 0-based.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 October 2010

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the array.

    Input, double A[N], an array to be index-sorted.

    Output, int R8VEC_SORT_HEAP_INDEX_A_NEW[N], contains the sort index.  The
    I-th element of the sorted array is A(INDX(I)).
*/
{
  double aval;
  int i;
  int *indx;
  int indxt;
  int ir;
  int j;
  int l;

  if ( n < 1 )
  {
    return NULL;
  }

  indx = ( int * ) malloc ( n * sizeof ( int ) );

  for ( i = 0; i < n; i++ )
  {
    indx[i] = i;
  }

  if ( n == 1 )
  {
    return indx;
  }

  l = n / 2 + 1;
  ir = n;

  for ( ; ; )
  {
    if ( 1 < l )
    {
      l = l - 1;
      indxt = indx[l-1];
      aval = a[indxt];
    }
    else
    {
      indxt = indx[ir-1];
      aval = a[indxt];
      indx[ir-1] = indx[0];
      ir = ir - 1;

      if ( ir == 1 )
      {
        indx[0] = indxt;
        break;
      }
    }

    i = l;
    j = l + l;

    while ( j <= ir )
    {
      if ( j < ir )
      {
        if ( a[indx[j-1]] < a[indx[j]] )
        {
          j = j + 1;
        }
      }

      if ( aval < a[indx[j-1]] )
      {
        indx[i-1] = indx[j-1];
        i = j;
        j = j + j;
      }
      else
      {
        j = ir + 1;
      }
    }
    indx[i-1] = indxt;
  }

  return indx;
}
/******************************************************************************/

int s_index_last_c ( char *s, char c )

/******************************************************************************/
/*
  Purpose:

    S_INDEX_LAST_C points to the last occurrence of a given character.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 January 2009

  Author:

    John Burkardt

  Parameters:

    Input, char *S, a pointer to a string.

    Input, char C, the character to search for.

    Output, int S_INDEX_LAST_C, the index in S of the last occurrence
    of the character, or -1 if it does not occur.
*/
{
  int n;
  char *t;

  n = strlen ( s ) - 1;
  t = s + strlen ( s ) - 1;

  while ( 0 <= n )
  {
    if ( *t == c )
    {
      return n;
    }
    t--;
    n--;
  }

  return (-1);
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

double s_to_r8 ( char *s, int *lchar, int *error )

/******************************************************************************/
/*
  Purpose:

    S_TO_R8 reads an R8 value from a string.

  Discussion:

    We have had some trouble with input of the form 1.0E-312.
    For now, let's assume anything less than 1.0E-20 is zero.

    This routine will read as many characters as possible until it reaches
    the end of the string, or encounters a character which cannot be
    part of the real number.

    Legal input is:

       1 blanks,
       2 '+' or '-' sign,
       2.5 spaces
       3 integer part,
       4 decimal point,
       5 fraction part,
       6 'E' or 'e' or 'D' or 'd', exponent marker,
       7 exponent sign,
       8 exponent integer part,
       9 exponent decimal point,
      10 exponent fraction part,
      11 blanks,
      12 final comma or semicolon.

    with most quantities optional.

  Example:

    S                 R

    '1'               1.0
    '     1   '       1.0
    '1A'              1.0
    '12,34,56'        12.0
    '  34 7'          34.0
    '-1E2ABCD'        -100.0
    '-1X2ABCD'        -1.0
    ' 2E-1'           0.2
    '23.45'           23.45
    '-4.2E+2'         -420.0
    '17d2'            1700.0
    '-14e-2'         -0.14
    'e2'              100.0
    '-12.73e-9.23'   -12.73 * 10.0**(-9.23)

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    24 June 2005

  Author:

    John Burkardt

  Parameters:

    Input, char *S, the string containing the
    data to be read.  Reading will begin at position 1 and
    terminate at the end of the string, or when no more
    characters can be read to form a legal real.  Blanks,
    commas, or other nonnumeric data will, in particular,
    cause the conversion to halt.

    Output, int *LCHAR, the number of characters read from
    the string to form the number, including any terminating
    characters such as a trailing comma or blanks.

    Output, int *ERROR, is TRUE (1) if an error occurred and FALSE (0)
    otherwise.

    Output, double S_TO_R8, the value that was read from the string.
*/
{
  char c;
  int ihave;
  int isgn;
  int iterm;
  int jbot;
  int jsgn;
  int jtop;
  int nchar;
  int ndig;
  double r;
  double rbot;
  double rexp;
  double rtop;
  char TAB = 9;

  nchar = s_len_trim ( s );
  *error = 0;
  r = 0.0;
  *lchar = -1;
  isgn = 1;
  rtop = 0.0;
  rbot = 1.0;
  jsgn = 1;
  jtop = 0;
  jbot = 1;
  ihave = 1;
  iterm = 0;

  for ( ; ; )
  {
    c = s[*lchar+1];
    *lchar = *lchar + 1;
/*
  Blank or TAB character.
*/
    if ( c == ' ' || c == TAB )
    {
      if ( ihave == 2 )
      {
      }
      else if ( ihave == 6 || ihave == 7 )
      {
        iterm = 1;
      }
      else if ( 1 < ihave )
      {
        ihave = 11;
      }
    }
/*
  Comma.
*/
    else if ( c == ',' || c == ';' )
    {
      if ( ihave != 1 )
      {
        iterm = 1;
        ihave = 12;
        *lchar = *lchar + 1;
      }
    }
/*
  Minus sign.
*/
    else if ( c == '-' )
    {
      if ( ihave == 1 )
      {
        ihave = 2;
        isgn = -1;
      }
      else if ( ihave == 6 )
      {
        ihave = 7;
        jsgn = -1;
      }
      else
      {
        iterm = 1;
      }
    }
/*
  Plus sign.
*/
    else if ( c == '+' )
    {
      if ( ihave == 1 )
      {
        ihave = 2;
      }
      else if ( ihave == 6 )
      {
        ihave = 7;
      }
      else
      {
        iterm = 1;
      }
    }
/*
  Decimal point.
*/
    else if ( c == '.' )
    {
      if ( ihave < 4 )
      {
        ihave = 4;
      }
      else if ( 6 <= ihave && ihave <= 8 )
      {
        ihave = 9;
      }
      else
      {
        iterm = 1;
      }
    }
/*
  Exponent marker.
*/
    else if ( ch_eqi ( c, 'E' ) || ch_eqi ( c, 'D' ) )
    {
      if ( ihave < 6 )
      {
        ihave = 6;
      }
      else
      {
        iterm = 1;
      }
    }
/*
  Digit.
*/
    else if ( ihave < 11 && '0' <= c && c <= '9' )
    {
      if ( ihave <= 2 )
      {
        ihave = 3;
      }
      else if ( ihave == 4 )
      {
        ihave = 5;
      }
      else if ( ihave == 6 || ihave == 7 )
      {
        ihave = 8;
      }
      else if ( ihave == 9 )
      {
        ihave = 10;
      }

      ndig = ch_to_digit ( c );

      if ( ihave == 3 )
      {
        rtop = 10.0 * rtop + ( double ) ndig;
      }
      else if ( ihave == 5 )
      {
        rtop = 10.0 * rtop + ( double ) ndig;
        rbot = 10.0 * rbot;
      }
      else if ( ihave == 8 )
      {
        jtop = 10 * jtop + ndig;
      }
      else if ( ihave == 10 )
      {
        jtop = 10 * jtop + ndig;
        jbot = 10 * jbot;
      }
    }
/*
  Anything else is regarded as a terminator.
*/
    else
    {
      iterm = 1;
    }
/*
  If we haven't seen a terminator, and we haven't examined the
  entire string, go get the next character.
*/
    if ( iterm == 1 || nchar <= *lchar + 1 )
    {
      break;
    }

  }
/*
  If we haven't seen a terminator, and we have examined the
  entire string, then we're done, and LCHAR is equal to NCHAR.
*/
  if ( iterm != 1 && (*lchar) + 1 == nchar )
  {
    *lchar = nchar;
  }
/*
  Number seems to have terminated.  Have we got a legal number?
  Not if we terminated in states 1, 2, 6 or 7!
*/
  if ( ihave == 1 || ihave == 2 || ihave == 6 || ihave == 7 )
  {
    *error = 1;
    return r;
  }
/*
  Number seems OK.  Form it.

  We have had some trouble with input of the form 1.0E-312.
  For now, let's assume anything less than 1.0E-20 is zero.
*/
  if ( jtop == 0 )
  {
    rexp = 1.0;
  }
  else
  {
    if ( jbot == 1 )
    {
      if ( jsgn * jtop < -20 )
      {
        rexp = 0.0;
      }
      else
      {
        rexp = pow ( ( double ) 10.0, ( double ) ( jsgn * jtop ) );
      }
    }
    else
    {
      if ( jsgn * jtop < -20 * jbot )
      {
        rexp = 0.0;
      }
      else
      {
        rexp = jsgn * jtop;
        rexp = rexp / jbot;
        rexp = pow ( ( double ) 10.0, ( double ) rexp );
      }
    }

  }

  r = isgn * rexp * rtop / rbot;

  return r;
}
/******************************************************************************/

int s_to_r8vec ( char *s, int n, double rvec[] )

/******************************************************************************/
/*
  Purpose:

    S_TO_R8VEC reads an R8VEC from a string.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    19 February 2001

  Author:

    John Burkardt

  Parameters:

    Input, char *S, the string to be read.

    Input, int N, the number of values expected.

    Output, double RVEC[N], the values read from the string.

    Output, int S_TO_R8VEC, is TRUE (1) if an error occurred and FALSE (0)
    otherwise.
*/
{
  int error;
  int i;
  int lchar;

  error = 0;

  for ( i = 0; i < n; i++ )
  {
    rvec[i] = s_to_r8 ( s, &lchar, &error );

    if ( error )
    {
      return error;
    }

    s = s + lchar;

  }

  return error;
}
/******************************************************************************/

int s_word_count ( char *s )

/******************************************************************************/
/*
  Purpose:

    S_WORD_COUNT counts the number of "words" in a string.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    30 January 2006

  Author:

    John Burkardt

  Parameters:

    Input, char *S, the string to be examined.

    Output, int S_WORD_COUNT, the number of "words" in the string.
    Words are presumed to be separated by one or more blanks.
*/
{
  int blank;
  int i;
  int word_num;

  word_num = 0;
  blank = 1;

  while ( *s ) 
  {
    if ( *s == ' ' || *s == '\n' )
    {
      blank = 1;
    }
    else if ( blank )
    {
      word_num = word_num + 1;
      blank = 0;
    }
    *s++;
  }

  return word_num;
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
