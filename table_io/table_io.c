# include <stdlib.h>
# include <stdio.h>
# include <string.h>
# include <math.h>
# include <time.h>

# include "table_io.h"

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
# define LINE_MAX 255

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
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "FILE_COLUMN_COUNT - Fatal error!\n" );
    fprintf ( stderr, "  Could not open the input file: \"%s\"\n", input_filename );
    exit ( 1 );
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
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "FILE_COLUMN_COUNT - Warning!\n" );
    fprintf ( stderr, "  The file does not seem to contain any data.\n" );
    exit ( 1 );
  }

  column_num = s_word_count ( line );

  return column_num;

# undef LINE_MAX
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
    fprintf ( stderr, "  Could not open the input file: \"%s\"\n", input_filename );
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

int i4_log_10 ( int i )

/******************************************************************************/
/*
  Purpose:

    I4_LOG_10 returns the integer part of the logarithm base 10 of an I4.

  Example:

        I  I4_LOG_10
    -----  --------
        0    0
        1    0
        2    0
        9    0
       10    1
       11    1
       99    1
      100    2
      101    2
      999    2
     1000    3
     1001    3
     9999    3
    10000    4

  Discussion:

    I4_LOG_10 ( I ) + 1 is the number of decimal digits in I.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    23 October 2007

  Author:

    John Burkardt

  Parameters:

    Input, int I, the number whose logarithm base 10 is desired.

    Output, int I4_LOG_10, the integer part of the logarithm base 10 of
    the absolute value of X.
*/
{
  int i_abs;
  int ten_pow;
  int value;

  if ( i == 0 )
  {
    value = 0;
  }
  else
  {
    value = 0;
    ten_pow = 10;

    i_abs = abs ( i );

    while ( ten_pow <= i_abs )
    {
      value = value + 1;
      ten_pow = ten_pow * 10;
    }

  }

  return value;
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

int i4_power ( int i, int j )

/******************************************************************************/
/*
  Purpose:

    I4_POWER returns the value of I^J.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    23 October 2007

  Author:

    John Burkardt

  Parameters:

    Input, int I, J, the base and the power.  J should be nonnegative.

    Output, int I4_POWER, the value of I^J.
*/
{
  int k;
  int value;

  if ( j < 0 )
  {
    if ( i == 1 )
    {
      value = 1;
    }
    else if ( i == 0 )
    {
      fprintf ( stderr, "\n" );
      fprintf ( stderr, "I4_POWER - Fatal error!\n" );
      fprintf ( stderr, "  I^J requested, with I = 0 and J negative.\n" );
      exit ( 1 );
    }
    else
    {
      value = 0;
    }
  }
  else if ( j == 0 )
  {
    if ( i == 0 )
    {
      fprintf ( stderr, "\n" );
      fprintf ( stderr, "I4_POWER - Fatal error!\n" );
      fprintf ( stderr, "  I^J requested, with I = 0 and J = 0.\n" );
      exit ( 1 );
    }
    else
    {
      value = 1;
    }
  }
  else if ( j == 1 )
  {
    value = i;
  }
  else
  {
    value = 1;
    for ( k = 1; k <= j; k++ )
    {
      value = value * i;
    }
  }
  return value;
}
/******************************************************************************/

int *i4mat_border_add ( int m, int n, int table[] )

/******************************************************************************/
/*
  Purpose:

    I4MAT_BORDER_ADD adds a "border" to an I4MAT.

  Discussion:

    An I4MAT is an array of I4's.

    We suppose the input data gives values of a quantity on nodes
    in the interior of a 2D grid, and we wish to create a new table
    with additional positions for the nodes that would be on the
    border of the 2D grid.

                  0 0 0 0 0 0
      * * * *     0 * * * * 0
      * * * * --> 0 * * * * 0
      * * * *     0 * * * * 0
                  0 0 0 0 0 0

    The illustration suggests the situation in which a 3 by 4 array
    is input, and a 5 by 6 array is to be output.

    The old data is shifted to its correct positions in the new array.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    28 May 2008

  Author:

    John Burkardt

  Parameters:

    Input, int M, the spatial dimension.

    Input, int N, the number of points.

    Input, int TABLE[M*N], the data.

    Output, int TABLE2[(M+2)*(N+2)], the augmented data.
*/
{
  int i;
  int j;
  int *table2;

  table2 = ( int * ) malloc ( ( m + 2 ) * ( n + 2 ) * sizeof ( int ) );

  for ( j = 0; j < n+2; j++ )
  {
    for ( i = 0; i < m+2; i++ )
    {
      if ( i == 0 || i == m+1 || j == 0 || j == n+1 )
      {
        table2[i+j*(m+2)] = 0;
      }
      else
      {
        table2[i+j*(m+2)] = table[(i-1)+(j-1)*m];
      }
    }
  }

  return table2;
}
/******************************************************************************/

int *i4mat_border_cut ( int m, int n, int table[] )

/******************************************************************************/
/*
  Purpose:

    I4MAT_BORDER_CUT cuts the "border" of an I4MAT.

  Discussion:

    An I4MAT is an array of I4's.

    We suppose the input data gives values of a quantity on nodes
    on a 2D grid, and we wish to create a new table corresponding only
    to those nodes in the interior of the 2D grid.

      0 0 0 0 0 0
      0 * * * * 0    * * * *
      0 * * * * 0 -> * * * *
      0 * * * * 0    * * * *
      0 0 0 0 0 0

    The illustration suggests the situation in which a 5 by 6 array
    is input, and a 3 by 4 array is to be output.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    28 May 2008

  Author:

    John Burkardt

  Parameters:

    Input, int M, the spatial dimension.

    Input, int N, the number of points.

    Input, int TABLE[M*N], the data.

    Output, int TABLE2[(M-2)*(N-2)], the "interior" data.
*/
{
  int i;
  int j;
  int *table2;

  if ( m <= 2 || n <= 2 )
  {
    return NULL;
  }

  table2 = ( int * ) malloc ( ( m - 2 ) * ( n - 2 ) * sizeof ( int ) );

  for ( j = 0; j < n-2; j++ )
  {
    for ( i = 0; i < m-2; i++ )
    {
      table2[i+j*(m-2)] = table[(i+1)+(j+1)*m];
    }
  }

  return table2;
}
/******************************************************************************/

int *i4mat_data_read ( char *input_filename, int m, int n )

/******************************************************************************/
/*
  Purpose:

    I4MAT_DATA_READ reads the data from an I4MAT file.

  Discussion:

    An I4MAT is an array of I4's.

    The file is assumed to contain one record per line.

    Records beginning with the '#' character are comments, and are ignored.
    Blank lines are also ignored.

    Each line that is not ignored is assumed to contain exactly (or at least)
    M real numbers, representing the coordinates of a point.

    There are assumed to be exactly (or at least) N such records.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    28 May 2008

  Author:

    John Burkardt

  Parameters:

    Input, char *INPUT_FILENAME, the name of the input file.

    Input, int M, the number of spatial dimensions.

    Input, int N, the number of points.  The program
    will stop reading data once N values have been read.

    Output, int I4MAT_DATA_READ[M*N], the data.
*/
{
# define LINE_MAX 255

  int error;
  char *got_string;
  FILE *input;
  int i;
  int j;
  char line[255];
  int *table;
  int *x;

  input = fopen ( input_filename, "r" );

  if ( !input )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "I4MAT_DATA_READ - Fatal error!\n" );
    fprintf ( stderr, "  Could not open the input file: \"%s\"\n", input_filename );
    exit ( 1 );
  }

  table = ( int * ) malloc ( m * n * sizeof ( int ) );

  x = ( int * ) malloc ( m * sizeof ( int ) );

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

    error = s_to_i4vec ( line, m, x );

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
 
void i4mat_header_read ( char *input_filename, int *m, int *n )
 
/******************************************************************************/
/*
  Purpose:

    I4MAT_HEADER_READ reads the header from an I4MAT file.

  Discussion:

    An I4MAT is an array of I4's.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    28 May 2008

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
    fprintf ( stderr, "I4MAT_HEADER_READ - Fatal error!\n" );
    fprintf ( stderr, "  FILE_COLUMN_COUNT failed.\n" );
    exit ( 1 );
  }

  *n = file_row_count ( input_filename );

  if ( *n <= 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "I4MAT_HEADER_READ - Fatal error!\n" );
    fprintf ( stderr, "  FILE_ROW_COUNT failed.\n" );
    exit ( 1 );
  }

  return;
}
/******************************************************************************/

int *i4mat_indicator_new ( int m, int n )

/******************************************************************************/
/*
  Purpose:

    I4MAT_INDICATOR_NEW sets up an "indicator" I4MAT.

  Discussion:

    An I4MAT is an array of I4's.

    The value of each entry suggests its location, as in:

      11  12  13  14
      21  22  23  24
      31  32  33  34

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

    Output, int I4MAT_INDICATOR_NEW[M*N], the table.
*/
{
  int *a;
  int fac;
  int i;
  int j;

  a = ( int * ) malloc ( m * n * sizeof ( int ) );

  fac = i4_power ( 10, i4_log_10 ( n ) + 1 );

  for ( i = 1; i <= m; i++ )
  {
    for ( j = 1; j <= n; j++ )
    {
      a[i-1+(j-1)*m] = fac * i + j;
    }
  }
  return a;
}
/******************************************************************************/

void i4mat_print ( int m, int n, int a[], char *title )

/******************************************************************************/
/*
  Purpose:

    I4MAT_PRINT prints an I4MAT, with an optional title.

  Discussion:

    An I4MAT is an array of I4's.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    28 May 2008

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows in A.

    Input, int N, the number of columns in A.

    Input, int A[M*N], the M by N matrix.

    Input, char *TITLE, a title.
*/
{
  int i;
  int j;
  int jhi;
  int jlo;

  i4mat_print_some ( m, n, a, 1, 1, m, n, title );

  return;
}
/******************************************************************************/

void i4mat_print_some ( int m, int n, int a[], int ilo, int jlo, int ihi, 
  int jhi, char *title )

/******************************************************************************/
/*
  Purpose:

    I4MAT_PRINT_SOME prints some of an I4MAT.

  Discussion:

    An I4MAT is an array of I4's.


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

    Input, int A[M*N], the matrix.

    Input, int ILO, JLO, IHI, JHI, designate the first row and
    column, and the last row and column to be printed.

    Input, char *TITLE, a title.
*/
{
# define INCX 10

  int i;
  int i2hi;
  int i2lo;
  int j;
  int j2hi;
  int j2lo;

  printf ( "\n" );
  printf ( "%s\n", title );
/*
  Print the columns of the matrix, in strips of INCX.
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
    printf ( "  Col:" );
    for ( j = j2lo; j <= j2hi; j++ )
    {
      printf ( "  %6d", j );
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
  Print out (up to INCX) entries in row I, that lie in the current strip.
*/
      printf ( "%5d", i );
      for ( j = j2lo; j <= j2hi; j++ )
      {
        printf ( "  %6d", a[i-1+(j-1)*m] );
      }
      printf ( "\n" );
    }
  }

  return;
# undef INCX
}
/******************************************************************************/

int *i4mat_read ( char *input_filename, int *m, int *n )

/******************************************************************************/
/*
  Purpose:

    I4MAT_READ reads information from an I4MAT file.

  Discussion:

    An I4MAT is an array of I4's.

    The file is assumed to contain one record per line.

    Records beginning with the '#' character are comments, and are ignored.
    Blank lines are also ignored.

    Each line that is not ignored is assumed to contain exactly (or at least)
    M real numbers, representing the coordinates of a point.

    There are assumed to be exactly (or at least) N such records.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    28 May 2008

  Author:

    John Burkardt

  Parameters:

    Input, char *INPUT_FILENAME, the name of the input file.

    Output, int *M, the number of spatial dimensions.

    Output, int *N, the number of points.

    Output, int I4MAT_READ[M*N], the data.
*/
{
  int *table;

  i4mat_header_read ( input_filename, m, n );

  table = i4mat_data_read ( input_filename, *m, *n );

  return table;
}
/******************************************************************************/

void i4mat_write ( char *output_filename, int m, int n, int table[] )

/******************************************************************************/
/*
  Purpose:

    I4MAT_WRITE writes an I4MAT file.

  Discussion:

    An I4MAT is an array of I4's.

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

    Input, int TABLE[M*N], the data.
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
    fprintf ( stderr, "I4MAT_WRITE - Fatal error!\n" );
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
      fprintf ( output, "  %d", table[i+j*m] );
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

int *lvec_data_read ( char *input_filename, int n )

/******************************************************************************/
/*
  Purpose:

    LVEC_DATA_READ reads the data from an LVEC file.

  Discussion:

    An LVEC is an array of L's.

    The file is assumed to contain one record per line.

    Records beginning with the '#' character are comments, and are ignored.
    Blank lines are also ignored.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    03 December 2010

  Author:

    John Burkardt

  Parameters:

    Input, char *INPUT_FILENAME, the name of the input file.

    Input, int N, the number of points.

    Output, int LVEC_DATA_READ[N], the data.
*/
{
# define LINE_MAX 255

  char *got_string;
  FILE *input;
  int j;
  int l;
  char line[255];
  int *table;

  input = fopen ( input_filename, "r" );

  if ( !input )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "LVEC_DATA_READ - Fatal error!\n" );
    fprintf ( stderr, "  Could not open the input file: \"%s\"\n", input_filename );
    exit ( 1 );
  }

  table = ( int * ) malloc ( n * sizeof ( int ) );

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

    table[j] = s_to_l ( line );
    j = j + 1;
  }

  fclose ( input );

  return table;

# undef LINE_MAX
}
/******************************************************************************/
 
void lvec_header_read ( char *input_filename, int *n )
 
/******************************************************************************/
/*
  Purpose:

    LVEC_HEADER_READ reads the header from an LVEC file.

  Discussion:

    An LVEC is a vector of L's.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    03 December 2010

  Author:

    John Burkardt

  Parameters:

    Input, char *INPUT_FILENAME, the name of the input file.

    Output, int *N, the number of points.
*/
{
  *n = file_row_count ( input_filename );

  if ( *n <= 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "LVEC_HEADER_READ - Fatal error!\n" );
    fprintf ( stderr, "  FILE_ROW_COUNT failed.\n" );
    exit ( 1 );
  }

  return;
}
/******************************************************************************/

void lvec_write ( char *output_filename, int n, int table[] )

/******************************************************************************/
/*
  Purpose:

    LVEC_WRITE writes an LVEC to a file.

  Discussion:

    An LVEC is a vector of L's.

    An L is an integer value such that 0 represents FALSE and 
    any nonzero value represents TRUE.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    02 December 2010

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
    fprintf ( stderr, "LVEC_WRITE - Fatal error!\n" );
    fprintf ( stderr, "  Could not open the output file.\n" );
    exit ( 1 );
  }
/*
  Write the data.
*/
  for ( j = 0; j < n; j++ )
  {
    if ( !table[j] )
    {
      fprintf ( output, "0\n" );
    }
    else
    {
      fprintf ( output, "1\n" );
    }
  }
/*
  Close the file.
*/
  fclose ( output );

  return;
}
/******************************************************************************/

float *r4mat_data_read ( char *input_filename, int m, int n )

/******************************************************************************/
/*
  Purpose:

    R4MAT_DATA_READ reads the data from an R4MAT file.

  Discussion:

    An R4MAT is an array of R4's.

    The file is assumed to contain one record per line.

    Records beginning with the '#' character are comments, and are ignored.
    Blank lines are also ignored.

    Each line that is not ignored is assumed to contain exactly (or at least)
    M real numbers, representing the coordinates of a point.

    There are assumed to be exactly (or at least) N such records.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    27 April 2011

  Author:

    John Burkardt

  Parameters:

    Input, char *INPUT_FILENAME, the name of the input file.

    Input, int M, the number of spatial dimensions.

    Input, int N, the number of points.  The program
    will stop reading data once N values have been read.

    Output, float R4MAT_DATA_READ[M*N], the data.
*/
{
# define LINE_MAX 255

  int error;
  char *got_string;
  FILE *input;
  int i;
  int j;
  char line[255];
  float *table;
  float *x;

  input = fopen ( input_filename, "r" );

  if ( !input )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R4MAT_DATA_READ - Fatal error!\n" );
    fprintf ( stderr, "  Could not open the input file: \"%s\"\n", input_filename );
    exit ( 1 );
  }

  table = ( float * ) malloc ( m * n * sizeof ( float ) );

  x = ( float * ) malloc ( m * sizeof ( float ) );

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

    error = s_to_r4vec ( line, m, x );

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
 
void r4mat_header_read ( char *input_filename, int *m, int *n )
 
/******************************************************************************/
/*
  Purpose:

    R4MAT_HEADER_READ reads the header from an R4MAT file.

  Discussion:

    An R4MAT is an array of R4's.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    27 April 2011

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
    fprintf ( stderr, "R4MAT_HEADER_READ - Fatal error!\n" );
    fprintf ( stderr, "  FILE_COLUMN_COUNT failed.\n" );
    exit ( 1 );
  }

  *n = file_row_count ( input_filename );

  if ( *n <= 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R4MAT_HEADER_READ - Fatal error!\n" );
    fprintf ( stderr, "  FILE_ROW_COUNT failed.\n" );
    exit ( 1 );
  }

  return;
}
/******************************************************************************/

float *r4mat_indicator_new ( int m, int n )

/******************************************************************************/
/*
  Purpose:

    R4MAT_INDICATOR_NEW sets up an "indicator" R4MAT.

  Discussion:

    An R8MAT is an array of R8's.

    The value of each entry suggests its location, as in:

      11  12  13  14
      21  22  23  24
      31  32  33  34

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    27 April 2011

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows of the matrix.
    M must be positive.

    Input, int N, the number of columns of the matrix.
    N must be positive.

    Output, float R4MAT_INDICATOR_NEW[M*N], the table.
*/
{
  double *a;
  int fac;
  int i;
  int j;

  a = ( float * ) malloc ( m * n * sizeof ( float ) );

  fac = i4_power ( 10, i4_log_10 ( n ) + 1 );

  for ( i = 1; i <= m; i++ )
  {
    for ( j = 1; j <= n; j++ )
    {
      a[i-1+(j-1)*m] = ( float ) ( fac * i + j );
    }
  }
  return a;
}
/******************************************************************************/

void r4mat_print ( int m, int n, float a[], char *title )

/******************************************************************************/
/*
  Purpose:

    R4MAT_PRINT prints an R4MAT, with an optional title.

  Discussion:

    An R4MAT is an array of R4's.

    Entry A(I,J) is stored as A[I+J*M]

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    27 April 2011

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows in A.

    Input, int N, the number of columns in A.

    Input, float A[M*N], the M by N matrix.

    Input, char *TITLE, a title.
*/
{
  r4mat_print_some ( m, n, a, 1, 1, m, n, title );

  return;
}
/******************************************************************************/

void r4mat_print_some ( int m, int n, float a[], int ilo, int jlo, int ihi, 
  int jhi, char *title )

/******************************************************************************/
/*
  Purpose:

    R4MAT_PRINT_SOME prints some of an R4MAT.

  Discussion:

    An R4MAT is an array of R4's.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    27 April 2011

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows of the matrix.
    M must be positive.

    Input, int N, the number of columns of the matrix.
    N must be positive.

    Input, float A[M*N], the matrix.

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

float *r4mat_read ( char *input_filename, int *m, int *n )

/******************************************************************************/
/*
  Purpose:

    R4MAT_READ reads information from an R4MAT file.

  Discussion:

    An R4MAT is an array of R4's.

    The file is assumed to contain one record per line.

    Records beginning with the '#' character are comments, and are ignored.
    Blank lines are also ignored.

    Each line that is not ignored is assumed to contain exactly (or at least)
    M real numbers, representing the coordinates of a point.

    There are assumed to be exactly (or at least) N such records.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    27 April 2011

  Author:

    John Burkardt

  Parameters:

    Input, char *INPUT_FILENAME, the name of the input file.

    Output, int *M, the number of spatial dimensions.

    Output, int *N, the number of points.  The program
    will stop reading data once N values have been read.

    Output, float R4MAT_READ[M*N], the data.
*/
{
  float *table;

  r4mat_header_read ( input_filename, m, n );

  table = r4mat_data_read ( input_filename, *m, *n );

  return table;
}
/******************************************************************************/

void r4mat_transpose_print ( int m, int n, float a[], char *title )

/******************************************************************************/
/*
  Purpose:

    R4MAT_TRANSPOSE_PRINT prints an R4MAT, transposed.

  Discussion:

    An R4MAT is an array of R4's.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    27 April 2011

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns.

    Input, float A[M*N], an M by N matrix to be printed.

    Input, char *TITLE, a title.
*/
{
  r4mat_transpose_print_some ( m, n, a, 1, 1, m, n, title );

  return;
}
/******************************************************************************/

void r4mat_transpose_print_some ( int m, int n, float a[], int ilo, int jlo, 
  int ihi, int jhi, char *title )

/******************************************************************************/
/*
  Purpose:

    R4MAT_TRANSPOSE_PRINT_SOME prints some of an R4MAT, transposed.

  Discussion:

    An R4MAT is an array of R4's.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    27 April 2011

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns.

    Input, float A[M*N], an M by N matrix to be printed.

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

float *r4mat_uniform_01 ( int m, int n, int *seed )

/******************************************************************************/
/*
  Purpose:

    R4MAT_UNIFORM_01 returns a unit pseudorandom R4MAT.

  Discussion:

    An R4MAT is an array of R4's.

    This routine implements the recursion

      seed = 16807 * seed mod ( 2**31 - 1 )
      unif = seed / ( 2**31 - 1 )

    The integer arithmetic never requires more than 32 bits,
    including a sign bit.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    27 April 2011

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

    Input, int M, N, the number of rows and columns.

    Input/output, int *SEED, the "seed" value.  Normally, this
    value should not be 0.  On output, SEED has 
    been updated.

    Output, float R4MAT_UNIFORM_01[M*N], a matrix of pseudorandom values.
*/
{
  int i;
  int i4_huge = 2147483647;
  int j;
  int k;
  float *r;

  if ( *seed == 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R4MAT_UNIFORM_01 - Fatal error!\n" );
    fprintf ( stderr, "  Input value of SEED = 0.\n" );
    exit ( 1 );
  }

  r = malloc ( m * n * sizeof ( float ) );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      k = *seed / 127773;

      *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

      if ( *seed < 0 )
      {
        *seed = *seed + i4_huge;
      }
/*
  Although SEED can be represented exactly as a 32 bit integer,
  it generally cannot be represented exactly as a 32 bit real number!
*/
      r[i+j*m] = ( float ) ( *seed ) * 4.656612875E-10;
    }
  }

  return r;
}
/******************************************************************************/

void r4mat_write ( char *output_filename, int m, int n, float table[] )

/******************************************************************************/
/*
  Purpose:

    R4MAT_WRITE writes an R4MAT file.

  Discussion:

    An R4MAT is an array of R4's.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    27 April 2011

  Author:

    John Burkardt

  Parameters:

    Input, char *OUTPUT_FILENAME, the output filename.

    Input, int M, the spatial dimension.

    Input, int N, the number of points.

    Input, float TABLE[M*N], the data.
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
    fprintf ( stderr, "R4MAT_WRITE - Fatal error!\n" );
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

double *r8mat_border_add ( int m, int n, double table[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_BORDER_ADD adds a "border" to an R8MAT.

  Discussion:

    An R8MAT is an array of R8's.

    We suppose the input data gives values of a quantity on nodes
    in the interior of a 2D grid, and we wish to create a new table
    with additional positions for the nodes that would be on the
    border of the 2D grid.

                  0 0 0 0 0 0
      * * * *     0 * * * * 0
      * * * * --> 0 * * * * 0
      * * * *     0 * * * * 0
                  0 0 0 0 0 0

    The illustration suggests the situation in which a 3 by 4 array
    is input, and a 5 by 6 array is to be output.

    The old data is shifted to its correct positions in the new array.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    28 May 2008

  Author:

    John Burkardt

  Parameters:

    Input, int M, the spatial dimension.

    Input, int N, the number of points.

    Input, double TABLE[M*N], the data.

    Output, double TABLE2[(M+2)*(N+2)], the augmented data.
*/
{
  int i;
  int j;
  double *table2;

  table2 = ( double * ) malloc ( ( m + 2 ) * ( n + 2 ) * sizeof ( double ) );

  for ( j = 0; j < n+2; j++ )
  {
    for ( i = 0; i < m+2; i++ )
    {
      if ( i == 0 || i == m+1 || j == 0 || j == n+1 )
      {
        table2[i+j*(m+2)] = 0.0;
      }
      else
      {
        table2[i+j*(m+2)] = table[(i-1)+(j-1)*m];
      }
    }
  }

  return table2;
}
/******************************************************************************/

double *r8mat_border_cut ( int m, int n, double table[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_BORDER_CUT cuts the "border" of an R8MAT.

  Discussion:

    An R8MAT is an array of R8's.

    We suppose the input data gives values of a quantity on nodes
    on a 2D grid, and we wish to create a new table corresponding only
    to those nodes in the interior of the 2D grid.

      0 0 0 0 0 0
      0 * * * * 0    * * * *
      0 * * * * 0 -> * * * *
      0 * * * * 0    * * * *
      0 0 0 0 0 0

    The illustration suggests the situation in which a 5 by 6 array
    is input, and a 3 by 4 array is to be output.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    28 May 2008

  Author:

    John Burkardt

  Parameters:

    Input, int M, the spatial dimension.

    Input, int N, the number of points.

    Input, double TABLE[M*N], the data.

    Output, double TABLE2[(M-2)*(N-2)], the "interior" data.
*/
{
  int i;
  int j;
  double *table2;

  if ( m <= 2 || n <= 2 )
  {
    return NULL;
  }

  table2 = ( double * ) malloc ( ( m - 2 ) * ( n - 2 ) * sizeof ( double ) );

  for ( j = 0; j < n-2; j++ )
  {
    for ( i = 0; i < m-2; i++ )
    {
      table2[i+j*(m-2)] = table[(i+1)+(j+1)*m];
    }
  }

  return table2;
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
/******************************************************************************/

double *r8mat_indicator_new ( int m, int n )

/******************************************************************************/
/*
  Purpose:

    R8MAT_INDICATOR_NEW sets up an "indicator" R8MAT.

  Discussion:

    An R8MAT is an array of R8's.

    The value of each entry suggests its location, as in:

      11  12  13  14
      21  22  23  24
      31  32  33  34

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    25 January 2005

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows of the matrix.
    M must be positive.

    Input, int N, the number of columns of the matrix.
    N must be positive.

    Output, double R8MAT_INDICATOR_NEW[M*N], the table.
*/
{
  double *a;
  int fac;
  int i;
  int j;

  a = ( double * ) malloc ( m * n * sizeof ( double ) );

  fac = i4_power ( 10, i4_log_10 ( n ) + 1 );

  for ( i = 1; i <= m; i++ )
  {
    for ( j = 1; j <= n; j++ )
    {
      a[i-1+(j-1)*m] = ( double ) ( fac * i + j );
    }
  }
  return a;
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

double *r8mat_read ( char *input_filename, int *m, int *n )

/******************************************************************************/
/*
  Purpose:

    R8MAT_READ reads information from an R8MAT file.

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

    04 June 2004

  Author:

    John Burkardt

  Parameters:

    Input, char *INPUT_FILENAME, the name of the input file.

    Output, int *M, the number of spatial dimensions.

    Output, int *N, the number of points.  The program
    will stop reading data once N values have been read.

    Output, double R8MAT_READ[M*N], the data.
*/
{
  double *table;

  r8mat_header_read ( input_filename, m, n );

  table = r8mat_data_read ( input_filename, *m, *n );

  return table;
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

double *r8mat_uniform_01 ( int m, int n, int *seed )

/******************************************************************************/
/*
  Purpose:

    R8MAT_UNIFORM_01 returns a unit pseudorandom R8MAT.

  Discussion:

    An R8MAT is an array of R8's.

    This routine implements the recursion

      seed = 16807 * seed mod ( 2^31 - 1 )
      unif = seed / ( 2^31 - 1 )

    The integer arithmetic never requires more than 32 bits,
    including a sign bit.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    03 October 2005

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

    Input, int M, N, the number of rows and columns.

    Input/output, int *SEED, the "seed" value.  Normally, this
    value should not be 0.  On output, SEED has 
    been updated.

    Output, double R8MAT_UNIFORM_01[M*N], a matrix of pseudorandom values.
*/
{
  int i;
  int i4_huge = 2147483647;
  int j;
  int k;
  double *r;

  if ( *seed == 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8MAT_UNIFORM_01 - Fatal error!\n" );
    fprintf ( stderr, "  Input value of SEED = 0.\n" );
    exit ( 1 );
  }

  r = malloc ( m * n * sizeof ( double ) );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      k = *seed / 127773;

      *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

      if ( *seed < 0 )
      {
        *seed = *seed + i4_huge;
      }
/*
  Although SEED can be represented exactly as a 32 bit integer,
  it generally cannot be represented exactly as a 32 bit real number!
*/
      r[i+j*m] = ( double ) ( *seed ) * 4.656612875E-10;
    }
  }

  return r;
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

double *r8vec_data_read ( char *input_filename, int n )

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

    11 June 2012

  Author:

    John Burkardt

  Parameters:

    Input, char *INPUT_FILENAME, the name of the input file.

    Input, int N, the number of points.

    Output, double R8VEC_DATA_READ[N], the data.
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
  double *table;
  double x;

  input = fopen ( input_filename, "r" );

  if ( !input )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8VEC_DATA_READ - Fatal error!\n" );
    fprintf ( stderr, "  Could not open the input file: \"%s\"\n", input_filename );
    exit ( 1 );
  }

  table = ( double * ) malloc ( n * sizeof ( double ) );

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

    x = s_to_r8 ( line, &lchar, &error );

    if ( error == 1 )
    {
      continue;
    }

    table[j] = x;
    j = j + 1;
  }

  fclose ( input );

  return table;

# undef LINE_MAX
}
/******************************************************************************/
 
void r8vec_header_read ( char *input_filename, int *n )
 
/******************************************************************************/
/*
  Purpose:

    R8VEC_HEADER_READ reads the header from an R8VEC file.

  Discussion:

    An R8VEC is a vector of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 June 2012

  Author:

    John Burkardt

  Parameters:

    Input, char *INPUT_FILENAME, the name of the input file.

    Output, int *N, the number of points.
*/
{
  *n = file_row_count ( input_filename );

  if ( *n <= 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8VEC_HEADER_READ - Fatal error!\n" );
    fprintf ( stderr, "  FILE_ROW_COUNT failed.\n" );
    exit ( 1 );
  }

  return;
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

    10 July 2011

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
    fprintf ( output, "  %24.16g\n", x[j] );
  }
/*
  Close the file.
*/
  fclose ( output );

  return;
}
/******************************************************************************/

void r8vec2_write ( char *output_filename, int n, double x[], double y[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC2_WRITE writes an R8VEC2 file.

  Discussion:

    An R8VEC2 is a pair of vectors of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    10 July 2011

  Author:

    John Burkardt

  Parameters:

    Input, char *OUTPUT_FILENAME, the output filename.

    Input, int N, the number of points.

    Input, double X[N], Y[N], the data.
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
    fprintf ( stderr, "R8VEC2_WRITE - Fatal error!\n" );
    fprintf ( stderr, "  Could not open the output file.\n" );
    exit ( 1 );
  }
/*
  Write the data.
*/
  for ( j = 0; j < n; j++ )
  {
    fprintf ( output, "  %24.16g  %24.16g\n", x[j], y[j] );
  }
/*
  Close the file.
*/
  fclose ( output );

  return;
}
/******************************************************************************/

void r8vla2_write ( char *output_filename, int m, int n, double a[m][n] )

/******************************************************************************/
/*
  Purpose:

    R8VLA2_WRITE writes an R8VLA2 file.

  Discussion:

    An R8VLA2 is a 2D variable length array (VLA).

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    28 June 2012

  Author:

    John Burkardt

  Parameters:

    Input, char *OUTPUT_FILENAME, the output filename.

    Input, int M, the spatial dimension.

    Input, int N, the number of points.

    Input, double TABLE[M][N], the table data.
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
    printf ( "\n" );
    printf ( "R8MAT_WRITE - Fatal error!\n" );
    printf ( "  Could not open the output file.\n" );
    return;
  }
/*
  Write the data.
*/
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      fprintf ( output, "  %24.16e", a[i][j] );
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

int s_to_i4 ( char *s, int *last, int *error )

/******************************************************************************/
/*
  Purpose:

    S_TO_I4 reads an I4 from a string.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    13 June 2003

  Author:

    John Burkardt

  Parameters:

    Input, char *S, a string to be examined.

    Output, int *LAST, the last character of S used to make IVAL.

    Output, int *ERROR is TRUE (1) if an error occurred and FALSE (0) otherwise.

    Output, int *S_TO_I4, the integer value read from the string.
    If the string is blank, then IVAL will be returned 0.
*/
{
  char c;
  int i;
  int isgn;
  int istate;
  int ival;

  *error = 0;
  istate = 0;
  isgn = 1;
  i = 0;
  ival = 0;

  while ( *s ) 
  {
    c = s[i];
    i = i + 1;
/*
  Haven't read anything.
*/
    if ( istate == 0 )
    {
      if ( c == ' ' )
      {
      }
      else if ( c == '-' )
      {
        istate = 1;
        isgn = -1;
      }
      else if ( c == '+' )
      {
        istate = 1;
        isgn = + 1;
      }
      else if ( '0' <= c && c <= '9' )
      {
        istate = 2;
        ival = c - '0';
      }
      else
      {
        *error = 1;
        return ival;
      }
    }
/*
  Have read the sign, expecting digits.
*/
    else if ( istate == 1 )
    {
      if ( c == ' ' )
      {
      }
      else if ( '0' <= c && c <= '9' )
      {
        istate = 2;
        ival = c - '0';
      }
      else
      {
        *error = 1;
        return ival;
      }
    }
/*
  Have read at least one digit, expecting more.
*/
    else if ( istate == 2 )
    {
      if ( '0' <= c && c <= '9' )
      {
        ival = 10 * (ival) + c - '0';
      }
      else
      {
        ival = isgn * ival;
        *last = i - 1;
        return ival;
      }

    }
  }
/*
  If we read all the characters in the string, see if we're OK.
*/
  if ( istate == 2 )
  {
    ival = isgn * ival;
    *last = s_len_trim ( s );
  }
  else
  {
    *error = 1;
    *last = 0;
  }

  return ival;
}
/******************************************************************************/

int s_to_i4vec ( char *s, int n, int ivec[] )

/******************************************************************************/
/*
  Purpose:

    S_TO_I4VEC reads an I4VEC from a string.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    19 February 2001

  Author:

    John Burkardt

  Parameters:

    Input, char *S, the string to be read.

    Input, int N, the number of values expected.

    Output, int IVEC[N], the values read from the string.

    Output, int S_TO_I4VEC, is TRUE (1) if an error occurred and
    FALSE (0) otherwise.
*/
{
  int error;
  int i;
  int lchar;

  error = 0;

  for ( i = 0; i < n; i++ )
  {
    ivec[i] = s_to_i4 ( s, &lchar, &error );

    if ( error )
    {
      return error;
    }

    s = s + lchar;

  }

  return error;
}
/******************************************************************************/

int s_to_l ( char *s )

/******************************************************************************/
/*
  Purpose:

    S_TO_L reads an L from a string.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    03 December 2010

  Author:

    John Burkardt

  Parameters:

    Input, char *S, the string to be read.

    Output, int S_TO_L, the logical value.
*/
{
  int i;
  int l;
  int length;

  length = strlen ( s );

  if ( length < 1 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "S_TO_L - Fatal error!\n" );
    fprintf ( stderr, "  Input string is empty.\n" );
    exit ( 1 );
  }

  for ( i = 0; i < length; i++ )
  {
    if ( s[i] == '0' ||
         s[i] == 'f' ||
         s[i] == 'F' )
    {
      l = 0;
      return l;
    }
    else if ( s[i] == '1' ||
              s[i] == 't' ||
              s[i] == 'T' )
    {
      l = 1;
      return 1;
    }
  }

  fprintf ( stderr, "\n" );
  fprintf ( stderr, "S_TO_L - Fatal error!\n" );
  fprintf ( stderr, "  Input did not contain boolean data.\n" );
  exit ( 1 );
}
/******************************************************************************/

float s_to_r4 ( char *s, int *lchar, int *error )

/******************************************************************************/
/*
  Purpose:

    S_TO_R4 reads an R4 value from a string.

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

    27 April 2011

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

    Output, float S_TO_R4, the value that was read from the string.
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
  float r;
  float rbot;
  float rexp;
  float rtop;
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
        rtop = 10.0 * rtop + ( float ) ndig;
      }
      else if ( ihave == 5 )
      {
        rtop = 10.0 * rtop + ( float ) ndig;
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
        rexp = pow ( ( float ) 10.0, ( float ) ( jsgn * jtop ) );
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
        rexp = pow ( ( float ) 10.0, ( float ) rexp );
      }
    }
  }

  r = isgn * rexp * rtop / rbot;

  return r;
}
/******************************************************************************/

int s_to_r4vec ( char *s, int n, float rvec[] )

/******************************************************************************/
/*
  Purpose:

    S_TO_R4VEC reads an R4VEC from a string.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    27 April 2011

  Author:

    John Burkardt

  Parameters:

    Input, char *S, the string to be read.

    Input, int N, the number of values expected.

    Output, float RVEC[N], the values read from the string.

    Output, int S_TO_R4VEC, is TRUE (1) if an error occurred and FALSE (0)
    otherwise.
*/
{
  int error;
  int i;
  int lchar;

  error = 0;

  for ( i = 0; i < n; i++ )
  {
    rvec[i] = s_to_r4 ( s, &lchar, &error );

    if ( error )
    {
      return error;
    }

    s = s + lchar;
  }

  return error;
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
