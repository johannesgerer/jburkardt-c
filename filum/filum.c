# include <stdlib.h>
# include <stdio.h>
# include <string.h>
# include <math.h>
# include <time.h>

# include "filum.h"

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

int ch_index_last ( char *s, char c )

/******************************************************************************/
/*
  Purpose:

    CH_INDEX_LAST finds the last occurrence of a character in a string.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    01 July 2008

  Author:

    John Burkardt

  Parameters:

    Input, char *S, a pointer to a string to be searched.

    Input, char C, the character to be searched for in s.

    Output, int CH_INDEX_LAST, the index in s of the last occurrence 
    of C, or -1 if c does not occur in s.
*/
{
  int i;
  int j;
  int nchar;

  j = -1;

  nchar = strlen ( s );

  for ( i = 0; i < nchar; i++ ) 
  {
    if ( s[i] == c ) 
    {
      j = i;
    }
  }

  return j;
}
/******************************************************************************/

int ch_is_digit ( char c )

/******************************************************************************/
/*
  Purpose:

    CH_IS_DIGIT returns TRUE if a character is a decimal digit.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    05 December 2003

  Author:

    John Burkardt

  Parameters:

    Input, char C, the character to be analyzed.

    Output, bool CH_IS_DIGIT, is TRUE if C is a digit.
*/
{
  if ( '0' <= c && c <= '9' )
  {
    return 1;
  }
  else
  {
    return 0;
  }
}
/******************************************************************************/

char ch_low ( char ch )

/******************************************************************************/
/*
  Purpose:

    CH_LOW lowercases a single character.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    30 November 2009

  Author:

    John Burkardt

  Parameters:

    Input, char CH, the character.

    Output, char CH_LOW, the lowercase character.
*/
{
  if ( 65 <= ch && ch <= 90 ) 
  {
    ch = ch + 32;
  }   

  return ch;
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

char ch_to_rot13 ( char ch )

/******************************************************************************/
/*
  Purpose:

    CH_TO_ROT13 converts a character to its ROT13 equivalent.

  Discussion:

    Two applications of CH_TO_ROT13 to a character will return the original.!

    As a further scrambling, digits are similarly rotated using
    a "ROT5" scheme.

  Example:

    Input:  Output:

    a       n
    C       P
    J       W
    1       6
    5       0

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    03 June 2008

  Author:

    John Burkardt

  Parameters:

    Input, character CH, the character to be converted.

    Output, character CH_TO_ROT13, the ROT13 equivalent of the character.
*/
{
  char rot13;
/*
  [0:4] -> [5:9]
*/
  if ( '0' <= ch && ch <= '4' )
  {
    rot13 = ch + 5;
  }
/*
  [5:9] -> [0:4]
*/
  else if ( '5' <= ch && ch <= '9' )
  {
    rot13 = ch - 5;
  }
/*
  [A:M] -> [N:Z]
*/
  else if ( 'A' <= ch && ch <= 'M' )
  {
    rot13 = ch + 13;
  }
/*
  [N:Z] -> [A:M]
*/
  else if ( 'N' <= ch && ch <= 'Z' )
  {
    rot13 = ch - 13;
  }
/*
  [a:m] -> [n:z]
*/
  else if ( 'a' <= ch && ch <= 'm' )
  {
    rot13 = ch + 13;
  }
/*
  [n:z] -> [a:m]
*/
  else if ( 'n' <= ch && ch <= 'z' )
  {
    rot13 = ch - 13;
  }
  else
  {
    rot13 = ch;
  }

  return rot13;
}
/******************************************************************************/

char digit_inc ( char c )

/******************************************************************************/
/*
  Purpose:

    DIGIT_INC increments a decimal digit.

  Example:

    Input  Output
    -----  ------
    '0'    '1'
    '1'    '2'
    ...
    '8'    '9'
    '9'    '0'
    'A'    'A'

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    01 July 2008

  Author:

    John Burkardt

  Parameters:

    Input, char C, a digit to be incremented.

    Output, char DIGIT_INC, the incremented digit.
*/
{
  if ( '0' <= c && c <= '8' )
  {
    return ( c + 1 );
  }
  else if ( c == '9' ) 
  {
    return '0';
  }
  else
  {
    return c;
  }
}
/******************************************************************************/

char digit_to_ch ( int digit )

/******************************************************************************/
/*
  Purpose:

    DIGIT_TO_CH returns the character representation of a decimal digit.

  Example:

    DIGIT   C
    -----  ---
      0    '0'
      1    '1'
    ...    ...
      9    '9'
     17    '*'

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    28 May 2003

  Author:

    John Burkardt

  Parameters:

    Input, int DIGIT, the digit value between 0 and 9.

    Output, char DIGIT_TO_CH, the corresponding character, or '*' if DIGIT
    was illegal.
*/
{
  if ( 0 <= digit && digit <= 9 )
  {
    return ( digit + 48 );
  }
  else
  {
    return '*';
  }
}
/******************************************************************************/

int file_byte_count ( char *input_filename )

/******************************************************************************/
/*
  Purpose:

    FILE_BYTE_COUNT counts the number of bytes in a file.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 September 2009

  Author:

    John Burkardt

  Reference:

    Bob Stout,
    C Snippet #1,
    Dr Dobb's,
    August 5 2009.

  Parameters:

    Input, char *INPUT_FILENAME, the name of the input file.

    Output, int FILE_BYTE_COUNT, the number of bytes in the file.
*/
{
  FILE *input_filepointer;
  long offset;
  int value;

  value = 0;

  input_filepointer = fopen ( input_filename, "rb" );

  if ( input_filepointer == NULL )
  {
    return value;
  }

  offset = 0L;

  fseek ( input_filepointer, offset, SEEK_END );
  value = ftell ( input_filepointer );
  fclose ( input_filepointer );

  return value;
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
/******************************************************************************/

int file_exist ( char *filename )

/******************************************************************************/
/*
  Purpose:

    FILE_EXIST reports whether a file exists.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    28 November 2009

  Author:

    John Burkardt

  Parameters:

    Input, char *FILENAME, the name of the file.

    Output, int FILE_EXIST, is TRUE (1) if the file exists.
*/
{
  FILE *filepointer;
  int value;

  filepointer = fopen ( filename, "r" );

  if ( filepointer == NULL )
  {
    value = 0;
  }
  else
  {
    value = 1;
  }
  return value;
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

void file_name_inc ( char *filename )

/******************************************************************************/
/*
  Purpose:

    FILE_NAME_INC increments a partially numeric file name.

  Discussion:

    It is assumed that the digits in the name, whether scattered or
    connected, represent a number that is to be increased by 1 on
    each call.  If this number is all 9's on input, the output number
    is all 0's.  Non-numeric letters of the name are unaffected.

    If the name is empty, then the routine stops.

    If the name contains no digits, the empty string is returned.

  Example:

      Input            Output
      -----            ------
      "a7to11.txt"     "a7to12.txt"  (typical case.  Last digit incremented)
      "a7to99.txt"     "a8to00.txt"  (last digit incremented, with carry.)
      "a9to99.txt"     "a0to00.txt"  (wrap around)
      "cat.txt"        " "           (no digits to increment)
      " "              STOP!         (error)

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    22 November 2011

  Author:

    John Burkardt

  Parameters:

    Input/output, char *FILENAME, the filename to be incremented.
*/
{
  char c;
  int change;
  int i;
  int n;
  char *t;

  n = s_len_trim ( filename );

  if ( n <= 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "FILE_NAME_INC - Fatal error!\n" );
    fprintf ( stderr, "  The input string is empty.\n" );
    exit ( 1 );
  }

  change = 0;

  t = filename + n - 1;
  
  while ( 0 < n )
  {
    if ( '0' <= *t && *t <= '9' )
    {
      change = change + 1;

      if ( *t == '9' )
      {
        *t = '0';
      }
      else
      {
        *t = *t + 1;
        return;
      }
    }
    t--;
    n--;
  }
/*
  No digits were found.  Return blank.
*/
  if ( change == 0 )
  {
    n = s_len_trim ( filename );
    t = filename + n - 1;
    while ( 0 < n )
    {
      *t = ' ';
      t--;
      n--;
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

int i4_uniform ( int a, int b, int *seed )

/******************************************************************************/
/*
  Purpose:

    I4_UNIFORM returns a scaled pseudorandom I4.

  Discussion:

    The pseudorandom number should be uniformly distributed
    between A and B.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 November 2006

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

    Input, int A, B, the limits of the interval.

    Input/output, int *SEED, the "seed" value, which should NOT be 0.
    On output, SEED has been updated.

    Output, int I4_UNIFORM, a number between A and B.
*/
{
  int i4_huge = 2147483647;
  int k;
  float r;
  int value;

  if ( *seed == 0 )
  {
    printf ( "\n" );
    printf ( "I4_UNIFORM - Fatal error!\n" );
    printf ( "  Input value of SEED = 0.\n" );
    exit ( 1 );
  }

  k = *seed / 127773;

  *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

  if ( *seed < 0 )
  {
    *seed = *seed + i4_huge;
  }

  r = ( float ) ( *seed ) * 4.656612875E-10;
/*
  Scale R to lie between A-0.5 and B+0.5.
*/
  r = ( 1.0 - r ) * ( ( float ) ( i4_min ( a, b ) ) - 0.5 ) 
    +         r   * ( ( float ) ( i4_max ( a, b ) ) + 0.5 );
/*
  Use rounding to convert R to an integer between A and B.
*/
  value = r4_nint ( r );

  value = i4_max ( value, i4_min ( a, b ) );
  value = i4_min ( value, i4_max ( a, b ) );

  return value;
}
/******************************************************************************/

int r4_nint ( float x )

/******************************************************************************/
/*
  Purpose:

    R4_NINT returns the nearest integer to an R4.

  Example:

        X         R4_NINT

      1.3         1
      1.4         1
      1.5         1 or 2
      1.6         2
      0.0         0
     -0.7        -1
     -1.1        -1
     -1.6        -2

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    05 May 2006

  Author:

    John Burkardt

  Parameters:

    Input, float X, the value.

    Output, int R4_NINT, the nearest integer to X.
*/
{
  int s;
  int value;

  if ( x < 0.0 )
  {
    s = -1;
  }
  else
  {
    s = 1;
  }
  value = s * ( int ) ( fabs ( x ) + 0.5 );

  return value;
}
/******************************************************************************/

double r8_uniform_01 ( int *seed )

/******************************************************************************/
/*
  Purpose:

    R8_UNIFORM_01 returns a unit pseudorandom R8.

  Discussion:

    This routine implements the recursion

      seed = 16807 * seed mod ( 2**31 - 1 )
      r8_uniform_01 = seed / ( 2**31 - 1 )

    The integer arithmetic never requires more than 32 bits,
    including a sign bit.

    If the initial seed is 12345, then the first three computations are

      Input     Output      R8_UNIFORM_01
      SEED      SEED

         12345   207482415  0.096616
     207482415  1790989824  0.833995
    1790989824  2035175616  0.947702

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 August 2004

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

    Input/output, int *SEED, the "seed" value.  Normally, this
    value should not be 0.  On output, SEED has been updated.

    Output, double R8_UNIFORM_01, a new pseudorandom variate, strictly between
    0 and 1.
*/
{
  int i4_huge = 2147483647;
  int k;
  double r;

  if ( *seed == 0 )
  {
    printf ( "\n" );
    printf ( "R8_UNIFORM_01 - Fatal error!\n" );
    printf ( "  Input value of SEED = 0.\n" );
    exit ( 1 );
  }

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
  r = ( double ) ( *seed ) * 4.656612875E-10;

  return r;
}
/******************************************************************************/

int s_eqi ( char *s1, char *s2 )

/******************************************************************************/
/*
  Purpose:

    S_EQI reports whether two strings are equal, ignoring case.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    03 June 2008

  Author:

    John Burkardt

  Parameters:

    Input, char *S1, char *S2, pointers to two strings.

    Output, int S_EQI, is true if the strings are equal. 
*/
{
  int i;
  int nchar;
  int nchar1;
  int nchar2;

  nchar1 = strlen ( s1 );
  nchar2 = strlen ( s2 );
  if ( nchar1 < nchar2 )
  {
    nchar = nchar1;
  }
  else
  {
    nchar = nchar2;
  }
/*
  The strings are not equal if they differ over their common length.
*/
  for ( i = 0; i < nchar; i++ ) 
  {

    if ( ch_cap ( s1[i] ) != ch_cap ( s2[i] ) ) 
    {
      return 0;
    }
  }
/*
  The strings are not equal if the longer one includes nonblanks
  in the tail.
*/
  if ( nchar < nchar1 ) 
  {
    for ( i = nchar; i < nchar1; i++ ) 
    {
      if ( s1[i] != ' ' ) 
      {
        return 0;
      }
    } 
  }
  else if ( nchar < nchar2 ) 
  {
    for ( i = nchar; i < nchar2; i++ )
    {
      if ( s2[i] != ' ' ) 
      {
        return 0;
      }
    } 
  }

  return 1;
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

void s_low ( char *s )

/******************************************************************************/
/*
  Purpose:

    S_LOW replaces all uppercase characters by lowercase ones.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    30 November 2009

  Author:

    John Burkardt

  Parameters:

    Input/output, char *S, a pointer to a string.  On output, all the
    characters in the string are uppercase.
*/
{
  char ch;

  while ( *s != '\0' ) 
  {
    ch = *s;
    *s = ch_low ( ch );
    s++;
  }

  return;
}
/******************************************************************************/

char *s_reverse ( char *s )

/******************************************************************************/
/*
  Purpose:

    S_REVERSE reverses the characters in a string.

  Example:

    Input        Output

    ' Cat'       'taC '
    'Goo gol  '  'log ooG  '

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    25 January 2010

  Author:

    John Burkardt

  Parameters:

    Input, char *S, the string to reverse.

    Output, char *S_REVERSE, the reversed string.
*/
{
  int i;
  int j;
  int n;
  char *s2;

  n = strlen ( s );

  s2 = ( char * ) malloc ( ( n + 1 ) * sizeof ( char ) );

  for ( j = 0; j < n; j++ )
  {
    *(s2+j) = *(s+n-j-1);
  }
  
  return s2;
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

    01 February 2012

  Author:

    John Burkardt

  Parameters:

    Input, char *S, the string to be read.

    Input, int N, the number of values expected.

    Output, int IVEC[N], the values read from the string.

    Output, int S_TO_I4VEC, is TRUE (1) if an error occurred and FALSE (0)
    otherwise.
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

void s_to_rot13 ( char *s )

/******************************************************************************/
/*
  Purpose:

    S_TO_ROT13 "rotates" the alphabetical characters in a string by 13 positions.

  Discussion:

    Two applications of the routine will return the original string.

  Example:

    Input:                      Output:

    abcdefghijklmnopqrstuvwxyz  nopqrstuvwxyzabcdefghijklm
    Cher                        Pure
    James Thurston Howell       Wnzrf Guhefgba Ubjryy
    0123456789                  5678901234

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

   02 January 2009

  Author:

    John Burkardt

  Parameters:

    Input/output, char *S, a string to be "rotated".
*/
{
  while ( *s != 0 )
  {
    *s = ch_to_rot13 ( *s );
    s++;
  }
  return;
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
