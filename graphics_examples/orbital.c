# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <string.h>

# include "dislin.h"

int main ( int argc, char *argv[] );
char ch_cap ( char ch );
int ch_eqi ( char ch1, char ch2 );
int ch_to_digit ( char ch );
int file_column_count ( char *input_filename );
int file_row_count ( char *input_filename );
float r4_huge ( void );
float r4_max ( float x, float y );
float r4_min ( float x, float y );
float *r4mat_data_read ( char *input_filename, int m, int n );
void r4mat_header_read ( char *input_filename, int *m, int *n );
float r4mat_max ( int m, int n, float a[] );
float r4mat_min ( int m, int n, float a[] );
float r4vec_max ( int n, float r4vec[] );
float r4vec_min ( int n, float r4vec[] );
int s_len_trim ( char *s );
float s_to_r4 ( char *s, int *lchar, int *error );
int s_to_r4vec ( char *s, int n, float rvec[] );
int s_word_count ( char *s );

/******************************************************************************/

int main ( int argc, char *argv[] )

/******************************************************************************/
/*
  Purpose:

    ORBITAL uses DISLIN to display a contour plot of Z(X,Y) data.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    20 May 2011

  Author:

    John Burkardt

  Reference:

    Helmut Michels,
    The Data Plotting Software DISLIN - version 10.4,
    Shaker Media GmbH, January 2010,
    ISBN13: 978-3-86858-517-9.
*/
{
  int i;
  int j;
  int k;
  int level;
  int level_num;
  float level_value;
  int m;
  int n;
  int nn;
  float *x;
  float xmax;
  float xmin;
  float *xyz;
  float *y;
  float ymax;
  float ymin;
  float *z;
  float zmax;
  float zmin;

  printf ( "\n" );
  printf ( "ORBITAL:\n" );
  printf ( "  C version:\n" );
  printf ( "  Use DISLIN to make a contour plot of Z(X,Y) data.\n" );
/*
  Read the data.
*/
  r4mat_header_read ( "orbital.txt", &m, &nn );

  xyz = r4mat_data_read ( "orbital.txt", m, nn );
/*
  Split the data.
  The contouring routine expects that data is along fixed X and Y coordinates,
  and so the X and Y data is to be given as vectors, not arrays.
*/
  n = 101;
  x = ( float * ) malloc ( n * sizeof ( float ) );
  y = ( float * ) malloc ( n * sizeof ( float ) );
  z = ( float * ) malloc ( n * n * sizeof ( float ) );

  k = 0;
  for ( i = 0; i < n; i++ )
  {
    x[i] = xyz[0+k*3];
    k = k + 1;
  }
  xmin = r4vec_min ( n, x );
  xmax = r4vec_max ( n, x );

  k = 0;
  for ( i = 0; i < n; i++ )
  {
    y[i] = xyz[1+k*3];
    k = k + n;
  }
  ymin = r4vec_min ( n, y );
  ymax = r4vec_max ( n, y );
/*
  Z is a table.
  The first dimension should contain values for constant Y.
*/
  k = 0;
  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      z[i+j*n] = xyz[2+k*3];
      k = k + 1;
    }
  }
  zmax = r4mat_max ( n, n, z );
  zmin = r4mat_min ( n, n, z );
/*
  Specify the format of the output file.
*/
  metafl ( "png" );
/*
  Indicate that new data overwrites old data.
*/
  filmod ( "delete" );
/*
  Specify the name of the output graphics file.
*/
  setfil ( "orbital.png" );
/*
  Choose the page size and orientation.
  'USA' is 2160 plot units wide and 2790 plot units high.
  'P' requests PORTRAIT orientation.
*/
  setpag ( "usap" );
/*
  For PNG output, reverse the default black background to white.
*/
  scrmod ( "reverse" );
/*
  Open DISLIN.
*/
  disini ( );
/*
  Plot a border around the page.
*/
  pagera ( );
/*
  Use the SIMPLX font.
*/
  simplx ( );
/*
  Set the axis origin in plot units to the right, and plot units DOWN.
*/
  axspos ( 230, 2500 );
/*
  Define the X and Y sizes of the axis system in plot units.
*/
  axslen ( 1700, 1700 );
/*
  Label the X and Y axes.
*/
  name ( "X axis", "X" );
  name ( "Y axis", "Y" );
/*
  Relate the physical coordinates to the axes, and specify tick marks.
*/
  graf ( xmin, xmax, xmin, 1.0, ymin, ymax, ymin, 1.0 );
/*
  BEGIN LEVEL 2 COMMANDS.
*/

/*
  Define the title.
*/
  titlin ( "Orbital contour plot", 1 );
  title ( );
/*
  Set color to "blue".
*/
  color ( "blue" );
/*
  Draw the contour plot.
*/
  level_num = 10;

  for ( level = 1; level <= level_num; level++ )
  {
    level_value = ( ( float ) ( level_num + 1 - level ) * zmin   
                  + ( float ) (                 level ) * zmax ) 
                  / ( float ) ( level_num + 1         );

    contur ( x, n, y, n, z, level_value );
  }
/*
  End this graph.
*/
  endgrf ( );
/*
  RETURN FROM LEVEL 2 TO LEVEL 1.
*/

/*
  Close DISLIN.
*/
  disfin ( );
/*
  Free memory.
*/
  free ( x );
  free ( xyz );
  free ( y );
  free ( z );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "ORBITAL:\n" );
  printf ( "  Normal end of execution.\n" );

  return 0;
}
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

float r4_huge ( void )

/******************************************************************************/
/*
  Purpose:

    R4_HUGE returns a "huge" R4.

  Discussion:

    The value returned by this function is NOT required to be the
    maximum representable R4.  This value varies from machine to machine,
    from compiler to compiler, and may cause problems when being printed.
    We simply want a "very large" but non-infinite number.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 August 2010

  Author:

    John Burkardt

  Parameters:

    Output, float R4_HUGE, a "huge" R4 value.
*/
{
  float value;

  value = 1.0E+30;

  return value;
}
/******************************************************************************/

float r4_max ( float x, float y )

/******************************************************************************/
/*
  Purpose:

    R4_MAX returns the maximum of two R4's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 May 2006

  Author:

    John Burkardt

  Parameters:

    Input, float X, Y, the quantities to compare.

    Output, float R4_MAX, the maximum of X and Y.
*/
{
  float value;

  if ( y < x )
  {
    value = x;
  }
  else
  {
    value = y;
  }
  return value;
}
/******************************************************************************/

float r4_min ( float x, float y )

/******************************************************************************/
/*
  Purpose:

    R4_MIN returns the minimum of two R4's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 May 2006

  Author:

    John Burkardt

  Parameters:

    Input, float X, Y, the quantities to compare.

    Output, float R4_MIN, the minimum of X and Y.
*/
{
  float value;

  if ( y < x )
  {
    value = y;
  }
  else
  {
    value = x;
  }
  return value;
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
    printf ( "\n" );
    printf ( "R4MAT_DATA_READ - Fatal error!\n" );
    printf ( "  Could not open the input file: \"%s\"\n", input_filename );
    return NULL;
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
    printf ( "\n" );
    printf ( "R4MAT_HEADER_READ - Fatal error!\n" );
    printf ( "  FILE_COLUMN_COUNT failed.\n" );
    *n = -1;
    return;
  }

  *n = file_row_count ( input_filename );

  if ( *n <= 0 )
  {
    printf ( "\n" );
    printf ( "R4MAT_HEADER_READ - Fatal error!\n" );
    printf ( "  FILE_ROW_COUNT failed.\n" );
    return;
  }

  return;
}
/******************************************************************************/

float r4mat_max ( int m, int n, float a[] )

/******************************************************************************/
/*
  Purpose:

    R4MAT_MAX returns the maximum entry of an R4MAT.

  Discussion:

    An R4MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 May 2010

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows in A.

    Input, int N, the number of columns in A.

    Input, float A[M*N], the M by N matrix.

    Output, float R4MAT_MAX, the maximum entry of A.
*/
{
  int i;
  int j;
  float value;

  value = - r4_huge ( );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      value = r4_max ( value, a[i+j*m] );
    }
  }

  return value;
}
/******************************************************************************/

float r4mat_min ( int m, int n, float a[] )

/******************************************************************************/
/*
  Purpose:

    R4MAT_MIN returns the minimum entry of an R4MAT.

  Discussion:

    An R4MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 May 2010

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows in A.

    Input, int N, the number of columns in A.

    Input, float A[M*N], the M by N matrix.

    Output, float R4MAT_MIN, the minimum entry of A.
*/
{
  int i;
  int j;
  float value;

  value = r4_huge ( );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      value = r4_min ( value, a[i+j*m] );
    }
  }

  return value;
}
/******************************************************************************/

float r4vec_max ( int n, float r4vec[] )

/******************************************************************************/
/*
  Purpose:

    R4VEC_MAX returns the value of the maximum element in a R4VEC.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 May 2006

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the array.

    Input, float R4VEC[N], a pointer to the first entry of the array.

    Output, float R4VEC_MAX, the value of the maximum element.  This
    is set to 0.0 if N <= 0.
*/
{
  int i;
  float value;

  value = - r4_huge ( );

  if ( n <= 0 )
  {
    return value;
  }

  for ( i = 0; i < n; i++ )
  {
    if ( value < r4vec[i] )
    {
      value = r4vec[i];
    }
  }
  return value;
}
/******************************************************************************/

float r4vec_min ( int n, float r4vec[] )

/******************************************************************************/
/*
  Purpose:

    R4VEC_MIN returns the value of the minimum element in a R4VEC.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 May 2006

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the array.

    Input, float R4VEC[N], the array to be checked.

    Output, float R4VEC_MIN, the value of the minimum element.
*/
{
  int i;
  float value;

  value = r4_huge ( );

  if ( n <= 0 )
  {
    return value;
  }

  for ( i = 0; i < n; i++ )
  {
    if ( r4vec[i] < value )
    {
      value = r4vec[i];
    }
  }
  return value;
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
