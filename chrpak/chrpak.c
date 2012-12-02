# include <stdlib.h>
# include <stdio.h>
# include <string.h>
# include <math.h>
# include <time.h>

# include "chrpak.h"

/******************************************************************************/

int a_to_i4 ( char ch )

/******************************************************************************/
/*
  Purpose:

    A_TO_I4 returns the index of an alphabetic character.

  Example:

    CH  A_TO_I4

    'A'   1
    'B'   2
    ...
    'Z'  26
    'a'  27
    'b'  28
    ...
    'z'  52
    '$'   0

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 July 2008

  Author:

    John Burkardt

  Parameters:

    Input, char CH, a character.

    Output, int A_TO_I4, is the alphabetic index of the character,
    between 1 and 26 if the character is a capital letter,
    between 27 and 52 if it is lower case, and 0 otherwise.
*/
{
  if ( 'A' <= ch && ch <= 'Z' )
  {
    return ( ( int ) ( ch - 'A' + 1 ) );
  }
  else if ( 'a' <= ch && ch <= 'z' )
  {
    return ( ( int ) ( ch - 'a' + 26 + 1 ) );
  }
  else
  {
    return 0;
  }
}
/******************************************************************************/

int base_to_i4 ( char *s, int base )

/******************************************************************************/
/*
  Purpose:

    BASE_TO_I4 returns the value of an integer represented in some base.

  Discussion:

    BASE = 1 is allowed, in which case we allow the digits '1' and '0',
    and we simply count the '1' digits for the result.

    Negative bases between -16 and -2 are allowed.

    The base -1 is allowed, and essentially does a parity check on
    a string of 1's.

  Example:

        Input      Output
    -------------  ------
         S   BASE       I
    ------  -----  ------
      '101'     2       5
    '-1000'     3     -27
      '100'     4      16
   '111111'     2      63
   '111111'    -2      21
   '111111'     1       6
   '111111'    -1       0

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    27 July 2000

  Author:

    John Burkardt

  Parameters:

    Input, char *S, the string.  The elements of S are
    blanks, a plus or minus sign, and digits.  Normally, the digits
    are representations of integers between 0 and |BASE-1|.  In the
    special case of base 1 or base -1, we allow both 0 and 1 as digits.

    Input, int BASE, the base in which the representation is given.
    Normally, 2 <= BASE <= 16.  However, there are two exceptions.

    Output, int BASE_TO_I4, the integer.
*/
{
  char c;
  int i;
  int ichr;
  int idig;
  int isgn;
  int istate;
  int nchar;

  nchar = s_len_trim ( s );

  if ( base == 0 )
  {
    printf ( "\n" );
    printf ( "BASE_TO_I4 - Serious error!\n" );
    printf ( "  The input base is zero.\n" );
    i = -1;
    return i;
  }

  if ( 16 < abs ( base ) )
  {
    printf ( "\n" );
    printf ( "BASE_TO_I4 - Serious error!\n" );
    printf ( "  The input base is greater than 16!\n" );
    i = -1;
    return i;
  }

  i = 0;
  istate = 0;
  isgn = 1 ;
  ichr = 1;

  while ( ichr <= nchar )
  {
    c = s[ichr-1];
//
//  Blank.
//
    if ( c == ' ' )
    {
      if ( istate == 2 )
      {
        break;
      }
    }
//
//  Sign, + or -.
//
    else if ( c == '-' )
    {
      if ( istate != 0 )
      {
        break;
      }
      istate = 1;
      isgn = -1;
    }
    else if ( c == '+' )
    {
      if ( istate != 0 )
      {
        break;
      }
      istate = 1;
    }
    else
/*
  Digit?
*/
    {
      idig = hex_digit_to_i4 ( c );

      if ( abs ( base ) == 1 && ( idig == 0 || idig == 1 ) )
      {
        i = base * i + idig;
        istate = 2;
      }
      else if ( 0 <= idig && idig < abs ( base ) )
      {
        i = base * i + idig;
        istate = 2;
      }
      else
      {
        printf ( "\n" );
        printf ( "BASE_TO_I4 - Serious error!\n" );
        printf ( "  Illegal digit = \"%c\"\n", c );
        printf ( "  Conversion halted prematurely!\n" );
        i = -1;
        return i;
      }
    }
    ichr = ichr + 1;
  }
/*
  Once we're done reading information, we expect to be in state 2.
*/
  if ( istate != 2 )
  {
    printf ( "\n" );
    printf ( "BASE_TO_I4 - Serious error!\n" );
    printf ( "  Unable to decipher input!\n" );
    i = -1;
    return i;
  }
/*
  Account for the sign.
*/
  i = isgn * i;

  return i;
}
/******************************************************************************/

void byte_to_int ( unsigned char *bvec, unsigned int *ival )

/******************************************************************************/
/*
  Purpose:

    BYTE_TO_INT converts 4 bytes into an unsigned integer.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    22 November 2010

  Author:

    John Burkardt

  Parameters:

    Input, unsigned char *BVEC, is a pointer to a character string.
    The contents of BVEC through BVEC+3 are the bytes of IVAL,
    from high order to low.

    Output, unsigned int IVAL, the integer represented by the bytes.
*/
{
  int i;

  *ival = 0;

  for ( i = 0; i < 4; i++ )
  {
    *ival = *ival << 8;
    *ival = *ival + *bvec;
    bvec = bvec + 1;
  }
  return;
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

void ch_count_cvec_add ( int n, unsigned char cvec[], int count[256] )

/******************************************************************************/
/*
  Purpose:

    CH_COUNT_CVEC_ADD adds a character vector to a character count.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 July 2008

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the vector.

    Input, unsigned char CVEC[n], a vector of characters.

    Input/output, int COUNT[256], the character counts.
*/
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    count[cvec[i]] = count[cvec[i]] + 1;
  }

  return;
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

int ch_index_first ( char *s, char c )

/******************************************************************************/
/*
  Purpose:

    CH_INDEX_FIRST finds the first occurrence of a character in a string.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    03 June 2008

  Author:

    John Burkardt

  Parameters:

    Input, char *S, a pointer to a string to be searched.

    Input, char C, the character to be searched for in the string.

    Output, int CH_INDEX_FIRST, the index in S of the first occurrence
    of c, or -1 if C does not occur in S.
*/
{
  int i;
  int nchar;

  nchar = strlen ( s );

  for ( i = 0; i < nchar; i++ )
  {
    if ( s[i] == c )
    {
      return i;
    }
  }

  return -1;
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

int ch_is_alpha ( char ch )

/******************************************************************************/
/*
  Purpose:

    CH_IS_ALPHA is TRUE if a charaacter is alphabetic.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 July 2008

  Author:

    John Burkardt

  Parameters:

    Input, char CH, a character to check.

    Output, int CH_IS_ALPHA is TRUE if the character is alphabetic.
*/
{
  int value;

  if ( ( 'a' <= ch && ch <= 'z' ) ||
       ( 'A' <= ch && ch <= 'Z' ) )
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

int ch_is_alphanumeric ( char ch )

/******************************************************************************/
/*
  Purpose:

    CH_IS_ALPHANUMERIC is TRUE if a character is alphanumeric.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 July 2008

  Author:

    John Burkardt

  Parameters:

    Input, char CH, a character to check.

    Output, int CH_IS_ALPHANUMERIC is TRUE if the character is alphanumeric.
*/
{
  int value;

  if ( ( 'a' <= ch && ch <= 'z' ) ||
       ( 'A' <= ch && ch <= 'Z' ) ||
       ( '0' <= ch && ch <= '9' ) )
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

int ch_is_control ( char ch )

/******************************************************************************/
/*
  Purpose:

    CH_IS_CONTROL is TRUE if a character is a control character.

  Discussion:

    A "control character" has ASCII code <= 31 or 127 <= ASCII code.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    27 September 2008

  Author:

    John Burkardt

  Parameters:

    Input, char CH, the character to be tested.

    Output, int CH_IS_CONTROL, TRUE if the character is a control
    character, and FALSE otherwise.
*/
{
  int value;

  if ( ch <= 31 || 127 <= ch )
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

    Output, int CH_IS_DIGIT, is TRUE if C is a digit.
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

int ch_is_lower ( char c )

/******************************************************************************/
/*
  Purpose:

    CH_IS_LOWER is TRUE if C is a lowercase alphabetic character.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    22 November 2010

  Author:

    John Burkardt

  Parameters:

    Input, char C, a character to check.

    Output, int CH_IS_LOWER is TRUE if C is a lowercase alphabetic character.
*/
{
  int value;

  if ( ( 'a' <= c && c <= 'z' ) )
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

int ch_is_space ( char c )

/******************************************************************************/
/*
  Purpose:

    CH_IS_SPACE is TRUE if a character represents "white space".

  Discussion:

    A white space character is a space, a form feed, a newline, a carriage
    return, a horizontal tab, or a vertical tab.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    03 June 2008

  Author:

    John Burkardt

  Parameters:

    Input, char C, the character to be analyzed.

    Output, int CH_IS_SPACE, is TRUE if C is a whitespace character.
*/
{
  int value;

  if ( c == ' ' )
  {
    value = 1;
  }
  else if ( c == '\f' )
  {
    value = 1;
  }
  else if ( c == '\n' )
  {
    value = 1;
  }
  else if ( c == '\r' )
  {
    value = 1;
  }
  else if ( c == '\t' )
  {
    value = 1;
  }
  else if ( c == '\v' )
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

char ch_read ( FILE *filein )

/******************************************************************************/
/*
  Purpose:

    CH_READ reads one character from a binary file.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    24 May 1999

  Author:

    John Burkardt

  Parameters:

    Input, FILE *FILEIN, a pointer to the file.

    Output, char CH_READ, the character that was read.
*/
{
  char c;

  c = ( char ) fgetc ( filein );

  return c;
}
/******************************************************************************/

int ch_roman_to_i4 ( char ch )

/******************************************************************************/
/*
  Purpose:

    CH_ROMAN_TO_I4 returns the integer value of a single Roman digit.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    03 June 2008

  Author:

    John Burkardt

  Parameters:

    Input, char CH, a Roman digit.

    Output, int CH_ROMAN_TO_I4, the value of the Roman
    numeral.  If the Roman numeral was not recognized, 0 is returned.
*/
{
  int value;

  if ( ch == 'M' || ch == 'm' )
  {
    value = 1000;
  }
  else if ( ch == 'D' || ch == 'd' )
  {
    value = 500;
  }
  else if ( ch == 'C' || ch == 'c' )
  {
    value = 100;
  }
  else if ( ch == 'L' || ch == 'l' )
  {
    value = 50;
  }
  else if ( ch == 'X' || ch == 'x' )
  {
    value = 10;
  }
  else if ( ch == 'V' || ch == 'v' )
  {
    value = 5;
  }
  else if ( ch == 'I' || ch == 'i' || ch == 'J' || ch == 'j' )
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

char ch_scrabble ( int tile )

/******************************************************************************/
/*
  Purpose:

    CH_SCRABBLE returns the character on a given Scrabble tile.

  Discussion:

    The tiles are numbered 1 to 100, and are labeled 'A' through 'Z',
    plus two blanks.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 January 2009

  Author:

    John Burkardt

  Parameters:

    Input, int TILE, the index of the desired Scrabble tile, between 1
    and 100.

    Output, char CH_SCRABBLE, the character on the given tile.
*/
{
  char scrabble[100] = {
    'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'B',
    'B', 'C', 'C', 'D', 'D', 'D', 'D', 'E', 'E', 'E',
    'E', 'E', 'E', 'E', 'E', 'E', 'E', 'E', 'E', 'F',
    'F', 'G', 'G', 'G', 'H', 'H', 'I', 'I', 'I', 'I',
    'I', 'I', 'I', 'I', 'I', 'J', 'K', 'L', 'L', 'L',
    'L', 'M', 'M', 'N', 'N', 'N', 'N', 'N', 'N', 'O',
    'O', 'O', 'O', 'O', 'O', 'O', 'O', 'P', 'P', 'Q',
    'R', 'R', 'R', 'R', 'R', 'R', 'S', 'S', 'S', 'S',
    'T', 'T', 'T', 'T', 'T', 'T', 'U', 'U', 'U', 'U',
    'V', 'V', 'W', 'W', 'X', 'X', 'Y', 'Z', ' ', ' ' };
  int value;

  if ( 1 <= tile && tile <= 100 )
  {
    value = scrabble[tile-1];
  }
  else
  {
    value = '?';
  }

  return value;
}
/******************************************************************************/

void ch_swap ( char *ch1, char *ch2 )

/******************************************************************************/
/*
  Purpose:

    CH_SWAP swaps two characters.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    06 December 2003

  Author:

    John Burkardt

  Parameters:

    Input/output, char *CH1, *CH2.  On output, the values have been
    interchanged.
*/
{
  char ch3;

   ch3 = *ch1;
  *ch1 = *ch2;
  *ch2 =  ch3;

  return;
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

int ch_to_digit_bin ( char c )

/******************************************************************************/
/*
  Purpose:

    CH_TO_DIGIT_BIN returns the integer value of a binary digit.

  Discussion:

    This routine handles other traditional binary pairs of "digits"
    besides '0' and '1'.

  Example:

     C   DIGIT
    ---  -----
    '0'    0
    '1'    1
    'T'    1
    'F'    0
    'Y'    1
    'N'    0
    '+'    1
    '-'    0

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    22 November 2010

  Author:

    John Burkardt

  Parameters:

    Input, char C, the binary digit.

    Output, int CH_TO_DIGIT_BIN, the corresponding integer value.  If C was
    'illegal', then DIGIT is -1.
*/
{
  int digit;

  if ( c == '0' ||
      c == 'F' ||
      c == 'f' ||
      c == '-' ||
      c == 'N' ||
      c == 'n' )
  {
    digit = 0;
  }
  else if ( c == '1' ||
            c == 'T' ||
            c == 't' ||
            c == '+' ||
            c == 'Y' ||
            c == 'y' )
  {
    digit = 1;
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

char ch_uniform ( char clo, char chi, int *seed )

/******************************************************************************/
/*
  Purpose:

    CH_UNIFORM returns a random character in a given range.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    04 May 2005

  Author:

    John Burkardt

  Parameters:

    Input, char CLO, CHI, the minimum and maximum acceptable characters.

    Input/output, int *SEED, a seed for the random number generator.

    Output, char CH_UNIFORM, the randomly chosen character.
*/
{
  char c;
  double d;

  d = r8_uniform_01 ( seed );

  c = ( char ) ( ( 1.0 - d ) * ( double ) clo + d * ( double ) chi );

  return c;
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

int hex_digit_to_i4 ( char c )

/******************************************************************************/
/*
  Purpose:

    HEX_DIGIT_TO_I4 converts a hexadecimal digit to an I4.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    22 November 2010

  Author:

    John Burkardt

  Parameters:

    Input, char C, the hexadecimal digit, '0'
    through '9', or 'A' through 'F', or also 'a' through 'f'
    are allowed.

    Output, int HEX_DIGIT_TO_I4, the corresponding integer,
    or -1 if C was illegal.
*/
{
  int i;

  if ( '0' <= c && c <= '9' )
  {
    i = ( int ) ( c - '0' );
  }
  else if ( 'A' <= c && c <= 'F' )
  {
    i = 10 + ( int ) ( c - 'A' );
  }
  else if ( 'a' <= c && c <= 'f' )
  {
    i = 10 + ( int ) ( c - 'a' );
  }
  else if ( c == ' ' )
  {
    i = 0;
  }
  else
  {
    i = -1;
  }

  return i;
}
/******************************************************************************/

int i4_huge ( void )

/******************************************************************************/
/*
  Purpose:

    I4_HUGE returns a "huge" I4.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    29 August 2006

  Author:

    John Burkardt

  Parameters:

    Output, int I4_HUGE, a "huge" integer.
*/
{
  return 2147483647;
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

void i4_swap ( int *i, int *j )

/******************************************************************************/
/*
  Purpose:

    I4_SWAP switches two I4's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 January 2002

  Author:

    John Burkardt

  Parameters:

    Input/output, int *I, *J.  On output, the values of I and
    J have been interchanged.
*/
{
  int k;

  k = *i;
  *i = *j;
  *j = k;

  return;
}
/******************************************************************************/

char i4_to_a ( int i )

/******************************************************************************/
/*
  Purpose:

    I4_TO_A returns the I-th alphabetic character.

  Example:

    I  I4_TO_A

   -8  ' '
    0  ' '
    1  'A'
    2  'B'
   ..
   26  'Z'
   27  'a'
   52  'z'
   53  ' '
   99  ' '

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 July 2008

  Author:

    John Burkardt

  Parameters:

    Input, int I, the index of the letter to be returned.
    0 is a space;
    1 through 26 requests 'A' through 'Z', (ASCII 65:90);
    27 through 52 requests 'a' through 'z', (ASCII 97:122);

    Output, char I4_TO_A, the requested alphabetic letter.
*/
{
  if ( i <= 0 )
  {
    return ( ' ' );
  }
  else if ( 1 <= i && i <= 26 )
  {
    return ( 'A' + i - 1 );
  }
  else if ( 27 <= i && i <= 52 )
  {
    return ( 'a' + i - 27 );
  }
  else
  {
    return ( ' ' );
  }
}
/******************************************************************************/

char i4_to_amino_code ( int i )

/******************************************************************************/
/*
  Purpose:

    I4_TO_AMINO_CODE converts an integer to an amino code.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 January 2009

  Author:

    John Burkardt

  Reference:

    Carl Branden, John Tooze,
    Introduction to Protein Structure,
    Garland Publishing, 1991.

  Parameters:

    Input, int I, the index of an amino acid, between 1 and 23.

    Output, char I4_TO_AMINO_CODE, the one letter code for an amino acid.
*/
{
# define N 23

  char c;
  static char ch_table[N] = {
    'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K',
    'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W',
    'X', 'Y', 'Z' };

  if ( 1 <= i & i <= N )
  {
    c = ch_table[i-1];
  }
  else
  {
    c = '?';
  }

  return c;
# undef N
}
/******************************************************************************/

char i4_to_hex_digit ( int i )

/******************************************************************************/
/*
  Purpose:

    I4_TO_HEX_DIGIT converts a (small) I4 to a hexadecimal digit.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 August 2009

  Author:

    John Burkardt

  Parameters:

    Input, int I, the integer, between 0 and 15.

    Output, char DI4_TO_HEX_DIGIT, the hexadecimal digit corresponding
    to the integer.
*/
{
  char c;

  if ( 0 <= i && i <= 9 )
  {
    c = '0' + i;
  }
  else if ( 10 <= i && i <= 15 )
  {
    c = 'a' + ( i - 10 );
  }
  else
  {
    c = '*';
  }

  return c;
}
/******************************************************************************/

char i4_to_isbn ( int i )

/******************************************************************************/
/*
  Purpose:

    I4_TO_ISBN converts an I4 to an ISBN digit.

  Discussion:

    Only the integers 0 through 10 can be input.  The representation
    of 10 is 'X'.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 August 2008

  Author:

    John Burkardt

  Reference:

    Book Industry Study Group,
    The Evolution in Product Identification:
    Sunrise 2005 and the ISBN-13,
    http://www.bisg.org/docs/The_Evolution_in_Product_ID.pdf

  Parameters:

    Input, int I, an integer between 0 and 10.

    Output, char I4_TO_ISBN, the ISBN character code of the integer.
    If I is illegal, then I4_TO_ISBN is set to '?'.
*/
{
       if ( i == 0 )
  {
    return '0';
  }
  else if ( i == 1 )
  {
    return '1';
  }
  else if ( i == 2 )
  {
    return '2';
  }
  else if ( i == 3 )
  {
    return '3';
  }
  else if ( i == 4 )
  {
    return '4';
  }
  else if ( i == 5 )
  {
    return '5';
  }
  else if ( i == 6 )
  {
    return '6';
  }
  else if ( i == 7 )
  {
    return '7';
  }
  else if ( i == 8 )
  {
    return '8';
  }
  else if ( i == 9 )
  {
    return '9';
  }
  else if ( i == 10 )
  {
    return 'X';
  }
  else
  {
    return '?';
  }
}
/******************************************************************************/

char *i4_to_month_abb ( int i )

/******************************************************************************/
/*
  Purpose:

    I4_TO_MONTH_ABB returns an abbreviated month name.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 January 2009

  Author:

    John Burkardt

  Parameters:

    Input, int I, the number of the desired month.

    Output, character *I4_TO_MONTH_ABB, a 3 character abbreviation for
    the month, such as 'Jan', 'Feb', 'Mar', and so on.
*/
{
  static char *month_list[13] =
  {
    "???",
    "Jan", "Feb", "Mar", "Apr", "May", "Jun",
    "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"
  };

  if ( 1 <= i && i <= 12 )
  {
    return month_list[i];
  }
  else
  {
    return month_list[0];
  }
}
/******************************************************************************/

char *i4_to_s ( int i )

/******************************************************************************/
/*
  Purpose:

    I4_TO_S converts an I4 to a string.

  Example:

    INTVAL  S

         1  1
        -1  -1
         0  0
      1952  1952
    123456  123456
   1234567  1234567

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 May 2011

  Author:

    John Burkardt

  Parameters:

    Input, int I, an integer to be converted.

    Output, char *I4_TO_S, the representation of the integer.
*/
{
  int digit;
  int j;
  int length;
  int ten_power;
  char *s;
  static double ten = 10.0;

  length = i4_log_10 ( i );

  ten_power = ( int ) ( pow ( ten, length ) );

  if ( i < 0 )
  {
    length = length + 1;
  }
/*
  Add one position for the trailing null.
*/
  length = length + 1;

  s = malloc ( length * sizeof ( char ) );

  if ( i == 0 )
  {
    s[0] = '0';
    s[1] = '\0';
    return s;
  }
/*
  Now take care of the sign.
*/
  j = 0;
  if ( i < 0 )
  {
    s[j] = '-';
    j = j + 1;
    i = abs ( i );
  }
/*
  Find the leading digit of I, strip it off, and stick it into the string.
*/
  while ( 0 < ten_power )
  {
    digit = i / ten_power;
    s[j] = digit_to_ch ( digit );
    j = j + 1;
    i = i - digit * ten_power;
    ten_power = ten_power / 10;
  }
/*
  Tack on the trailing NULL.
*/
  s[j] = '\0';
  j = j + 1;

  return s;
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

int isbn_to_i4 ( char c )

/******************************************************************************/
/*
  Purpose:

    ISBN_TO_I4 converts an ISBN character into an I4.

  Discussion:

    The characters '0' through '9' stand for themselves, but
    the character 'X' or 'x' stands for 10.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    24 July 2010

  Author:

    John Burkardt

  Reference:

    Book Industry Study Group,
    The Evolution in Product Identification:
    Sunrise 2005 and the ISBN-13,
    http://www.bisg.org/docs/The_Evolution_in_Product_ID.pdf

  Parameters:

    Input, char C, the ISBN character code to be converted.

    Output, int ISBN_TO_I4, the numeric value of the character
    code, between 0 and 10.  This value is returned as -1 if C is
    not a valid character code.
*/
{
  int value;

  if ( '0' <= c && c <= '9' )
  {
    value = c - '0';
  }
  else if ( c == 'X' || c == 'x' )
  {
    value = 10;
  }
  else
  {
    value = -1;
  }

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
    s = - 1;
  }
  else
  {
    s = + 1;
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

      seed = 16807 * seed mod ( 2^31 - 1 )
      r8_uniform_01 = seed / ( 2^31 - 1 )

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
    Springer Verlag, pages 201-202, 1983.

    Pierre L'Ecuyer,
    Random Number Generation,
    in Handbook of Simulation
    edited by Jerry Banks,
    Wiley Interscience, page 95, 1998.

    Bennett Fox,
    Algorithm 647:
    Implementation and Relative Efficiency of Quasirandom
    Sequence Generators,
    ACM Transactions on Mathematical Software,
    Volume 12, Number 4, pages 362-376, 1986.

    P A Lewis, A S Goodman, J M Miller,
    A Pseudo-Random Number Generator for the System/360,
    IBM Systems Journal,
    Volume 8, pages 136-143, 1969.

  Parameters:

    Input/output, int *SEED, the "seed" value.  Normally, this
    value should not be 0.  On output, SEED has been updated.

    Output, double R8_UNIFORM_01, a new pseudorandom variate, strictly between
    0 and 1.
*/
{
  int k;
  double r;

  k = *seed / 127773;

  *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

  if ( *seed < 0 )
  {
    *seed = *seed + 2147483647;
  }
/*
  Although SEED can be represented exactly as a 32 bit integer,
  it generally cannot be represented exactly as a 32 bit real number!
*/
  r = ( ( double ) ( *seed ) ) * 4.656612875E-10;

  return r;
}
/******************************************************************************/

void s_adjustl ( char *s )

/******************************************************************************/
/*
  Purpose:

    S_ADJUSTL flushes a string left.

  Discussion:

    Both blanks and tabs are treated as "white space".

    This routine is similar to the FORTRAN90 ADJUSTL routine.

  Example:

    Input             Output

    '     Hello'      'Hello'
    ' Hi there!  '    'Hi there!'
    'Fred  '          'Fred'

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 June 2003

  Author:

    John Burkardt

  Parameters:

    Input/output, char *S.
    On input, S is a string of characters.
    On output, any initial blank or tab characters have been cut.
*/
{
  int i;
  int length;
  int nonb;
  char TAB = 9;
/*
  Check the length of the string to the last nonblank.
  If nonpositive, return.
*/
  length = s_len_trim ( s );

  if ( length <= 0 )
  {
    return;
  }
/*
  Find NONB, the location of the first nonblank, nontab.
*/
  nonb = 0;

  for ( i = 0; i < length; i++ )
  {
    if ( s[i] != ' ' && s[i] != TAB )
    {
      nonb = i;
      break;
    }
  }

  if ( 0 < nonb )
  {
    for ( i = nonb; i < length; i++ )
    {
      s[i-nonb] = s[i];
    }

    s[length-nonb] = '\0';
  }
  return;
}
/******************************************************************************/

int s_begin ( char *s1, char *s2 )

/******************************************************************************/
/*
  Purpose:

    S_BEGIN reports whether string 1 begins with string 2.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    22 November 2010

  Author:

    John Burkardt

  Parameters:

    Input, char *S1, char *S2, two strings.

    Output, int S_BEGIN, is true if S1 is the same as S2 up to
    the end of S2, and false otherwise.
*/
{
  int i;
  int n;
  int n1;
  int n2;

  n1 = strlen ( s1 );
  n2 = strlen ( s2 );

  if ( n1 < n2 )
  {
    return 0;
  }

  for ( i = 0; i < n2; i++ )
  {
    if ( ch_cap ( s1[i] ) != ch_cap ( s2[i] ) )
    {
      return 0;
    }
  }
  return 1;
}
/******************************************************************************/

void s_blank_delete ( char *s )

/******************************************************************************/
/*
  Purpose:

    S_BLANK_DELETE removes blanks and left justifies the remainder.

  Discussion:

    All TAB characters are also removed.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 July 2008

  Author:

    John Burkardt

  Parameters:

    Input/output, char *S, the string to be transformed.
*/
{
  char *get;
  char *put;
  char TAB = 9;

  put = s;
  get = s;

  while ( *get != '\0' )
  {
    if ( *get != ' ' && *get != TAB )
    {
      *put = *get;
      put = put + 1;
    }
    get = get + 1;
  }

  *put = *get;

  return;
}
/******************************************************************************/

void s_blanks_delete ( char *s )

/******************************************************************************/
/*
  Purpose:

    S_BLANKS_DELETE replaces consecutive blanks by one blank.

  Discussion:

    The remaining characters are left justified and right padded with blanks.
    TAB characters are converted to spaces.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    03 June 2008

  Author:

    John Burkardt

  Parameters:

    Input/output, char *S, the string to be transformed.
*/
{
  int blank;
  char *get;
  char *put;

  put = s;
  get = s;

  blank = 1;

  while ( *get != '\0' )
  {
    if ( *get != ' ' )
    {
      *put = *get;
      put = put + 1;
      blank = 0;
    }
    else if ( !blank )
    {
      *put = *get;
      put = put + 1;
      blank = 1;
    }
    get = get + 1;
  }

  *put = '\0';

  return;
}
/******************************************************************************/

void s_cap ( char *s )

/******************************************************************************/
/*
  Purpose:

    S_CAP capitalizes all the characters in a string.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    03 June 2008

  Author:

    John Burkardt

  Parameters:

    Input/output, char *S, is a pointer to a string to be capitalized.  On output,
    all alphabetic characters in the string have been capitalized.
*/
{
  int i;
  int nchar;

  nchar = strlen ( s );

  for ( i = 0; i < nchar; i++ )
  {
    s[i] = ch_cap ( s[i] );
  }

  return;
}
/******************************************************************************/

char *s_cat ( char *s1, char *s2 )

/******************************************************************************/
/*
  Purpose:

    S_CAT concatenates two strings to make a third string.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 January 2009

  Author:

    John Burkardt

  Parameters:

    Input, char *S1, the "prefix" string.

    Input, char *S2, the "postfix" string.

    Output, char *S_CAT, the string made by
    concatenating S1 and S2, ignoring any trailing blanks.
*/
{
  int i;
  int l1;
  int l2;
  char *s3;

  l1 = s_len_trim ( s1 );
  l2 = s_len_trim ( s2 );

  if ( l1 == 0 && l2 == 0 )
  {
    s3 = NULL;
    return s3;
  }

  s3 = ( char * ) malloc ( ( l1 + l2 + 1 ) * sizeof ( char ) );

  for ( i = 0; i < l1; i++ )
  {
    s3[i] = s1[i];
  }

  for ( i = 0; i < l2; i++ )
  {
    s3[l1+i] = s2[i];
  }

  s3[l1+l2] = '\0';

  return s3;
}
/******************************************************************************/

int s_ch_count ( char *s, char ch )

/******************************************************************************/
/*
  Purpose:

    S_CH_COUNT counts occurrences of a particular character in a string.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 November 2011

  Author:

    John Burkardt

  Parameters:

    Input, char *S, the string.

    Input, char CH, the character to be counted.

    Output, int S_CH_COUNT, the number of occurrences.
*/
{
  int ch_count;
  int i;
  int s_length;

  ch_count = 0;

  s_length = strlen ( s );

  for ( i = 0; i < s_length; i++ )
  {
    if ( s[i] == ch )
    {
      ch_count = ch_count + 1;
    }
  }

  return ch_count;
}
/******************************************************************************/

char *s_ch_delete ( char *s, char ch )

/******************************************************************************/
/*
  Purpose:

    S_CH_DELETE removes all occurrences of a character from a string.

  Discussion:

    Each time the given character is found in the string, the characters
    to the right of the string are shifted over one position.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 November 2011

  Author:

    John Burkardt

  Parameters:

    Input, char *S, the string to be transformed.

    Input, character CH, the character to be removed.

    Output, char *S_CH_DELETE, a copy of the string with the character removed.
*/
{
  int ch_num;
  int get;
  int put;
  int s_length;
  char *s2;

  s_length = strlen ( s );

  ch_num = s_ch_count ( s, ch );

  s2 = ( char * ) malloc ( ( s_length - ch_num + 1 ) * sizeof ( char ) );

  put = 0;

  for ( get = 0; get < s_length; get++ )
  {
    if ( s2[get] != ch )
    {
      s2[put] = s[get];
      put = put + 1;
    }
  }
  return s2;
}
/******************************************************************************/

void s_control_blank ( char *s )

/******************************************************************************/
/*
  Purpose:

    S_CONTROL_BLANK replaces control characters with blanks.

  Discussion:

    A "control character" has ASCII code <= 31 or 127 <= ASCII code.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    27 September 2008

  Author:

    John Burkardt

  Parameters:

    Input/output, char *S, the string to be transformed.
*/
{
  char *get;

  get = s;

  while ( *get != '\0' )
  {
    if ( ch_is_control ( *get ) )
    {
      *get = ' ';
    }
    get = get + 1;
  }
  return;
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

char *s_escape_tex ( char *s1 )

/******************************************************************************/
/*
  Purpose:

    S_ESCAPE_TEX de-escapes TeX escape sequences.

  Discussion:

    In particular, every occurrence of the characters '\', '_',
    '^', '{' and '}' will be replaced by '\\', '\_', '\^',
    '\{' and '\}'.  A TeX interpreter, on seeing these character
    strings, is then likely to return the original characters.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 January 2007

  Author:

    John Burkardt

  Parameters:

    Input, char *S1, the string to be de-escaped.

    Output, char *S_ESCAPE_TEX, a copy of the string,
    modified to avoid TeX escapes.
*/
{
  char ch;
  int s1_length;
  int s1_pos;
  char *s2;
  int s2_length;
  int s2_pos;
  int slash_count;

  s1_length = strlen ( s1 );
/*
 We need to know how many slashes occur in S1, so we
  can allocate S2.

  Note that, alas, the backslash is also the escape in C++,
  so we have to say '\\' when we mean '\'!
*/
  slash_count = 0;
  for ( s1_pos = 0; s1_pos < s1_length; s1_pos++ )
  {
    ch = s1[s1_pos];

    if ( ch == '\\' ||
         ch == '_' ||
         ch == '^' ||
         ch == '{' ||
         ch == '}' )
    {
      slash_count = slash_count + 1;
    }
  }

  s2_length = s1_length + slash_count;
  s2 = ( char * ) malloc ( ( s2_length + 1 ) * sizeof ( char ) );
/*
  Now copy S1 into S2.
*/
  s1_pos = 0;
  s2_pos = 0;

  for ( s1_pos = 0; s1_pos < s1_length; s1_pos++ )
  {
    ch = s1[s1_pos];

    if ( ch == '\\' ||
         ch == '_' ||
         ch == '^' ||
         ch == '{' ||
         ch == '}' )
    {
      s2[s2_pos] = '\\';
      s2_pos = s2_pos + 1;
    }

    s2[s2_pos] = ch;
    s2_pos = s2_pos + 1;
  }

  s2[s2_pos] = '\0';
  s2_pos = s2_pos + 1;

  return s2;
}
/******************************************************************************/

char *s_first_ch ( char *s, char ch )

/******************************************************************************/
/*
  Purpose:

    S_FIRST_CH points to the first occurrence of a character in a string.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 January 2009

  Author:

    John Burkardt

  Parameters:

    Input, char *S, a pointer to a string.

    Input, char CH, a character.

    Output, char *S_FIRST_CH, a pointer to the first occurrence of the
    character in the string, or NULL if the character does not occur.
*/
{
  while ( *s != ch )
  {
    if ( *s == '\0' )
    {
      return NULL;
    }
    s++;
  }

  return s;
}
/******************************************************************************/

char *s_first_nonblank ( char *s )

/******************************************************************************/
/*
  Purpose:

    S_FIRST_NONBLANK points to the first nonblank character in a string.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 November 2002

  Author:

    John Burkardt

  Parameters:

    Input, char *S, a pointer to a string.

    Output, char *S_FIRST_NONBLANK, a pointer to the first nonblank character
    in the string, or NULL if the entire string is blank.
*/
{
  char *t = NULL;

  while ( *s != '\0' )
  {
    if ( *s != ' ' )
    {
      t = s;
      break;
    }
    s++;
  }

  return t;
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

void s_newline_to_null ( char *s )

/******************************************************************************/
/*
  Purpose:

    S_NEWLINE_TO_NULL replaces carriage returns or newlines by nulls.

  Discussion:

    The function FGETS will read a string containing a line of text read from
    input.  However, the string will include the linefeed character '/n', or,
    for a PC-formatted file, the carriage return and linefeed pair '/r' + '/n'.

    It may be desirable that the string not contain these characters.  The
    easiest way to deal with this is simply to replace the first instance of
    '/r' or '/n' by a null character, which terminates the string.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    29 November 2010

  Author:

    John Burkardt

  Parameters:

    Input/output, char *S, a pointer to a string.  On output, the first
    carriage return or line feed encountered has been replaced by a null.
*/
{
  int i;
  int n;

  n = strlen ( s );

  for ( i = 0; i < n; i++ )
  {
/*
  Handle carriage return.
*/
    if ( s[i] == '\r' )
    {
      s[i] = '\0';
      return;
    }
/*
  Handle linefeed.
*/
    if ( s[i] == '\n' )
    {
      s[i] = '\0';
      return;
    }
  }

  return;
}
/******************************************************************************/

void s_replace_ch ( char *s, char c1, char c2 )

/******************************************************************************/
/*
  Purpose:

    S_REPLACE_CH replaces all occurrences of one character by another.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    22 November 2010

  Author:

    John Burkardt

  Parameters:

    Input/output, char *S, the string.

    Input, char C1, C2, the character to be replaced, and the
    replacement character.
*/
{
  int i;
  int s_length;

  s_length = strlen ( s );

  for ( i = 0; i < s_length; i++ )
  {
    if ( s[i] == c1 )
    {
      s[i] = c2;
    }
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

void s_sort_a ( char *s )

/******************************************************************************/
/*
  Purpose:

    S_SORT_A sorts a string into ascending order.

  Discussion:

    The string is assumed to be short, and so a simple bubble sort is used.

    ALL the characters are sorted, including blanks and punctuation.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 July 2008

  Author:

    John Burkardt

  Parameters:

    Input/output, char *S, the string to be sorted.
*/
{
  char c;
  int i;
  int j;
  int k;
  int s_length;

  s_length = strlen ( s );

  for ( i = 0; i < s_length - 1; i++ )
  {
    c = s[i];
    j = i;

    for ( k = i + 1; k < s_length; k++ )
    {
      if ( s[k] < s[j] )
      {
        j = k;
      }
    }

    if ( i != j )
    {
      s[i] = s[j];
      s[j] = c;
    }
  }

  return;
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
  static float ten = 10.0;

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
        rexp = pow ( ten, jsgn * jtop );
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
        rexp = pow ( ten, rexp );
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

    02 May 2011

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
  static double ten = 10.0;

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
        rexp = pow ( ten, jsgn * jtop );
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
        rexp = pow ( ten, rexp );
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

void s_word_cap ( char *s )

/******************************************************************************/
/*
  Purpose:

    S_WORD_CAP capitalizes the first character of each word in a string.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    22 November 2010

  Author:

    John Burkardt

  Parameters:

    Input/output, char *S, the string to be capitalized.
*/
{
  int blank;
  int i;
  int s_length;

  s_length = strlen ( s );

  blank = 1;

  for ( i = 0; i < s_length; i++ )
  {
    if ( blank )
    {
      s[i] = ch_cap ( s[i] );
    }
    else
    {
      s[i] = ch_low ( s[i] );
    }
    blank = ( s[i] == ' ' );
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

    28 May 20080

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

char *s_word_extract_first ( char *s )

/******************************************************************************/
/*
  Purpose:

    S_WORD_EXTRACT_FIRST extracts the first word from a string.

  Discussion:

    A "word" is a string of characters terminated by a blank or
    the end of the string.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 January 2009

  Author:

    John Burkardt

  Parameters:

    Input/output, char *S, the string.  On output, the first
    word has been removed, and the remaining string has been shifted left.

    Output, char *S_WORD_EXTRACT_FIRST, the leading word of the string.
    NULL is returned if there are no more words to read.
*/
{
  int get1;
  int get2;
  int i;
  int s_len;
  char *w;

  s_len = s_len_trim ( s );

  if ( s_len < 1 )
  {
    return NULL;
  }
/*
  Find the first nonblank.
*/
  get1 = 0;

  for ( ; ; )
  {
    if ( s_len <= get1 )
    {
      return NULL;
    }

    if ( *(s+get1) != ' ' )
    {
      break;
    }
    get1 = get1 + 1;
  }
/*
  Look for the last contiguous nonblank.
*/
  get2 = get1;

  for ( ; ; )
  {
    if ( s_len <= get2 + 1 )
    {
      break;
    }

    if ( *(s+get2+1) == ' ' )
    {
      break;
    }
    get2 = get2 + 1;
  }
/*
  Copy the word.
*/
  w = ( char * ) malloc ( ( get2 + 2 - get1 ) * sizeof ( char ) );
  for ( i = 0; i < get2+1-get1; i++ )
  {
    *(w+i) = *(s+get1+i);
    *(s+get1+i) = ' ';
  }
  *(w+get2+2-get1) = '\0';
/*
  Shift the string.
*/
  s_adjustl ( s );

  return w;
}
/******************************************************************************/

void sort_heap_external ( int n, int *indx, int *i, int *j, int isgn )

/******************************************************************************/
/*
  Purpose:

    SORT_HEAP_EXTERNAL externally sorts a list of items into ascending order.

  Discussion:

    The actual list is not passed to the routine.  Hence it may
    consist of integers, reals, numbers, names, etc.  The user,
    after each return from the routine, will be asked to compare or
    interchange two items.

    The current version of this code mimics the FORTRAN version,
    so the values of I and J, in particular, are FORTRAN indices.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 February 2004

  Author:

    FORTRAN77 original version by Albert Nijenhuis, Herbert Wilf.
    C version by John Burkardt

  Reference:

    Albert Nijenhuis, Herbert Wilf,
    Combinatorial Algorithms,
    Academic Press, 1978, second edition,
    ISBN 0-12-519260-6.

  Parameters:

    Input, int N, the length of the input list.

    Input/output, int *INDX.
    The user must set INDX to 0 before the first call.
    On return,
      if INDX is greater than 0, the user must interchange
      items I and J and recall the routine.
      If INDX is less than 0, the user is to compare items I
      and J and return in ISGN a negative value if I is to
      precede J, and a positive value otherwise.
      If INDX is 0, the sorting is done.

    Output, int *I, *J.  On return with INDX positive,
    elements I and J of the user's list should be
    interchanged.  On return with INDX negative, elements I
    and J are to be compared by the user.

    Input, int ISGN. On return with INDX negative, the
    user should compare elements I and J of the list.  If
    item I is to precede item J, set ISGN negative,
    otherwise set ISGN positive.
*/
{
  static int i_save = 0;
  static int j_save = 0;
  static int k = 0;
  static int k1 = 0;
  static int n1 = 0;
/*
  INDX = 0: This is the first call.
*/
  if ( *indx == 0 )
  {

    i_save = 0;
    j_save = 0;
    k = n / 2;
    k1 = k;
    n1 = n;
  }
/*
  INDX < 0: The user is returning the results of a comparison.
*/
  else if ( *indx < 0 )
  {
    if ( *indx == -2 )
    {
      if ( isgn < 0 )
      {
        i_save = i_save + 1;
      }
      j_save = k1;
      k1 = i_save;
      *indx = -1;
      *i = i_save;
      *j = j_save;
      return;
    }

    if ( 0 < isgn )
    {
      *indx = 2;
      *i = i_save;
      *j = j_save;
      return;
    }

    if ( k <= 1 )
    {
      if ( n1 == 1 )
      {
        i_save = 0;
        j_save = 0;
        *indx = 0;
      }
      else
      {
        i_save = n1;
        j_save = 1;
        n1 = n1 - 1;
        *indx = 1;
      }
      *i = i_save;
      *j = j_save;
      return;
    }
    k = k - 1;
    k1 = k;
  }
/*
  0 < INDX: the user was asked to make an interchange.
*/
  else if ( *indx == 1 )
  {
    k1 = k;
  }

  for ( ; ; )
  {

    i_save = 2 * k1;

    if ( i_save == n1 )
    {
      j_save = k1;
      k1 = i_save;
      *indx = -1;
      *i = i_save;
      *j = j_save;
      return;
    }
    else if ( i_save <= n1 )
    {
      j_save = i_save + 1;
      *indx = -2;
      *i = i_save;
      *j = j_save;
      return;
    }

    if ( k <= 1 )
    {
      break;
    }

    k = k - 1;
    k1 = k;
  }

  if ( n1 == 1 )
  {
    i_save = 0;
    j_save = 0;
    *indx = 0;
    *i = i_save;
    *j = j_save;
  }
  else
  {
    i_save = n1;
    j_save = 1;
    n1 = n1 - 1;
    *indx = 1;
    *i = i_save;
    *j = j_save;
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

  printf ( "%s\n", time_buffer );

  return;
# undef TIME_SIZE
}
