# include <stdlib.h>
# include <stdio.h>
# include <stdbool.h>
# include <string.h>
# include <ctype.h>
# include <math.h>
# include <time.h>

# include "pbmb_io.h"

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

bool pbmb_check_data ( int xsize, int ysize, int *barray )

/******************************************************************************/
/*
  Purpose:

    PBMB_CHECK_DATA checks the data for a binary portable bit map file.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 April 2003

  Author:

    John Burkardt

  Parameters:

    Input, int XSIZE, YSIZE, the number of rows and columns of data.

    Input, int *BARRAY, the array of XSIZE by YSIZE bits.

    Output, bool PBMB_CHECK_DATA, is true if an error occurred.
*/
{
  int i;
  int *indexb;
  int j;

  if ( xsize <= 0 )
  {
    printf ( "\n" );
    printf ( "PBMB_CHECK_DATA - Fatal error!\n" );
    printf ( "  XSIZE <= 0\n" );
    printf ( "  XSIZE = %d\n", xsize );
    return true;
  }

  if ( ysize <= 0 )
  {
    printf ( "\n" );
    printf ( "PBMB_CHECK_DATA - Fatal error!\n" );
    printf ( "  YSIZE <= 0\n" );
    printf ( "  YSIZE = %d\n", ysize );
    return true;
  }

  if ( barray == NULL )
  {
    printf ( "\n" );
    printf ( "PBMB_CHECK_DATA - Fatal error!\n" );
    printf ( "  Null pointer to data.\n" );
    return true;
  }

  indexb = barray;

  for ( j = 0; j < ysize; j++ )
  {
    for ( i = 0; i < xsize; i++ )
    {
      if ( *indexb != 0 && *indexb != 1 )
      {
        printf ( "\n" );
        printf ( "PBMB_CHECK_DATA - Fatal error!\n" );
        printf ( "  b[%d][%d] = %d.\n", i, j, *indexb );
        return true;
      }

      indexb = indexb + 1;
    }
  }

  return false;
}
/******************************************************************************/

bool pbmb_example ( int xsize, int ysize, int *barray )

/******************************************************************************/
/*
  Purpose:

    PBMB_EXAMPLE sets up some sample PBMB data.

  Discussion:

    The data represents an ellipse.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    15 June 2012

  Author:

    John Burkardt

  Parameters:

    Input, int XSIZE, YSIZE, the number of rows and columns of data.
    Values of 200 would be reasonable.

    Output, int *BARRAY, the array of XSIZE by YSIZE data values.

    Output, bool PBMB_EXAMPLE, is true if an error occurred.
*/
{
  int i;
  int *indexb;
  int j;
  float r;
  float test;
  float x;
  float xc;
  float y;
  float yc;
 
  indexb = barray;
  if ( xsize < ysize )
  {
    r = ( float ) xsize / 3.0;
  }
  else
  {
    r = ( float ) ysize / 3.0;
  }
  xc = ( xsize ) / 2.0;
  yc = ( ysize ) / 2.0;

  for ( i = 0; i < ysize; i++ )
  {
    y = ( float ) i;
    for ( j = 0; j < xsize; j++ )
    {
      x = ( float ) j;
      test = r - sqrt ( ( x - xc ) * ( x - xc ) 
               + 0.75 * ( y - yc ) * ( y - yc ) );
      if ( fabs ( test ) <= 3.0 )
      {
        *indexb = 1;
      }
      else
      {
        *indexb = 0;
      }
      indexb = indexb + 1;
    }
  }

  return false;
}
/******************************************************************************/

bool pbmb_read ( char *file_name, int *xsize, int *ysize, int **barray )

/******************************************************************************/
/*
  Purpose:

    PBMB_READ reads the header and data from a binary portable bit map file.
 
  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:
 
    04 October 1998
 
  Author:
 
    John Burkardt

  Parameters:

    Input, char *FILE_NAME, the name of the file containing the binary
    portable bit map data.

    Output, int *XSIZE, *YSIZE, the number of rows and columns of data.

    Output, int **BARRAY, the array of XSIZE by YSIZE data values.

    Output, bool PBMB_READ, is
    true, if an error was detected, or
    false, if no error occurred.
*/
{
  FILE *file_pointer;
  int numbytes;
  int result;

  file_pointer = fopen ( file_name, "rb" );

  if ( file_pointer == NULL )
  {
    printf ( "\n" );
    printf ( "PBMB_READ: Fatal error!\n" );
    printf ( "  Cannot open the input file %s.\n", file_name );
    return true;
  }
/*
  Read the header.
*/
  result = pbmb_read_header ( file_pointer, xsize, ysize );

  if ( result != 0 )
  {
    printf ( "\n" );
    printf ( "PBMB_READ: Fatal error!\n" );
    printf ( "  PBMB_READ_HEADER failed.\n" );
    return true;
  }
/*
  Allocate storage for the data.
*/
  numbytes = (*xsize) * (*ysize) * sizeof ( int );

  *barray = ( int * ) malloc ( numbytes );

  if ( *barray == NULL )
  {
    printf ( "\n" );
    printf ( "PBMB_READ: Fatal error!\n" );
    printf ( "  Unable to allocate memory for data.\n" );
    printf ( "  Seeking %d bytes.\n", numbytes );
    return true;
  }
/*
  Read the data.
*/
  result = pbmb_read_data ( file_pointer, *xsize, *ysize, *barray );

  if ( result != 0 )
  {
    printf ( "\n" );
    printf ( "PBMB_READ: Fatal error!\n" );
    printf ( "  PBMB_READ_DATA failed.\n" );
    return true;
  }
/*
  Close the file.
*/
  fclose ( file_pointer );

  return false;
}
/******************************************************************************/

bool pbmb_read_data ( FILE *file_pointer, int xsize, int ysize, int *barray )

/******************************************************************************/
/*
  Purpose:

    PBMB_READ_DATA reads the data in a binary portable bit map file.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    29 April 1999

  Author:

    John Burkardt

  Parameters:

    Input, FILE *FILE_POINTER, a pointer to the file containing the binary
    portable bit map data.

    Input, int XSIZE, YSIZE, the number of rows and columns of data.

    Output, int *BARRAY, the array of XSIZE by YSIZE data values.

    Output, bool PBMB_READ_DATA, is
    true, if an error was detected, or
    false, if the data was read.
*/
{
  int bit;
  int c;
  unsigned char c2;
  int i;
  int *indexb;
  int j;
  int k;
  int numbyte;

  indexb = barray;
  numbyte = 0;

  for ( j = 0; j < ysize; j++ )
  {
    for ( i = 0; i < xsize; i++ )
    {
      if ( i%8 == 0 )
      {
        c = fgetc ( file_pointer );

        if ( c == EOF )
        {
          printf ( "\n" );
          printf ( "PBMB_READ_DATA: Failed reading data byte %d.\n", numbyte );
          return true;
        }
        c2 = ( unsigned char ) c;
        numbyte = numbyte + 1;
      }

      k = 7 - i%8;
      bit = ( c2 >> k )%2;

      *indexb = bit;
      indexb = indexb + 1;
    }
  }
  return false;
}
/******************************************************************************/

bool pbmb_read_header ( FILE *file_pointer, int *xsize, int *ysize )

/******************************************************************************/
/*
  Purpose:

    PBMB_READ_HEADER reads the header of a binary portable bit map file.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    29 April 1999

  Author:

    John Burkardt

  Parameters:

    Input, FILE *FILE_POINTER, a pointer to the file containing the binary
    portable bit map data.

    Output, int *XSIZE, *YSIZE, the number of rows and columns of data.

    Output, bool PBMB_READ_HEADER, is
    true, if an error was detected, or
    false, if the header was read.
*/
{
  int c_val;
  int count;
  int flag;
  int nchar;
  int state;
  char string[80];

  state = 0;
  nchar = 0;

  for ( ; ; )
  {

    c_val = fgetc ( file_pointer );

    if ( c_val == EOF )
    {
      return true;
    }
/*
  If not whitespace, add the character to the current string.
*/
    flag = isspace ( c_val );

    if ( !flag )
    {
      string[nchar] = c_val;
      nchar = nchar + 1;
    }
/*
  See if we have finished an old item, or begun a new one.
*/
    if ( state == 0 )
    {
      if ( !flag )
      {
        state = 1;
      }
      else 
      {
        return true;
      }
    }
    else if ( state == 1 )
    {
      if ( flag )
      {
        string[nchar] = 0;
        nchar = nchar + 1;
        if ( strcmp ( string, "P4" ) != 0 && strcmp ( string, "p4" ) != 0 )
        {
          printf ( "\n" );
          printf ( "PBMB_READ_HEADER: Fatal error.\n" );
          printf ( "  Bad magic number = %s.\n", string );
          return true;
        }
        nchar = 0;
        state = 2;
      }
    }
    else if ( state == 2 )
    {
      if ( !flag )
      {
        state = 3;
      }
    }
    else if ( state == 3 )
    {
      if ( flag )
      {
        string[nchar] = 0;
        nchar = nchar + 1;
        count = sscanf ( string, "%d", xsize );
        if ( count == EOF )
        {
          return true;
        }
        nchar = 0;
        state = 4;
      }
    }
    else if ( state == 4 )
    {
      if ( !flag )
      {
        state = 5;
      }
    }
    else if ( state == 5 )
    {
      if ( flag )
      {
        string[nchar] = 0;
        nchar = nchar + 1;
        count = sscanf ( string, "%d", ysize );
        if ( count == EOF )
        {
          return true;
        }
        nchar = 0;
        return false;
      }
    }
  }
}
/******************************************************************************/

bool pbmb_read_test ( char *file_name )

/******************************************************************************/
/*
  Purpose:

    PBMB_READ_TEST tests the binary portable bit map read routines.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    04 October 1998

  Author:

    John Burkardt

  Parameters:

    Input, char *FILE_NAME, the name of the file containing the binary
    portable bit map data.

    Output, int PBMB_READ_TEST, is
    true, if an error was detected, or
    false, if the test was carried out.
*/
{
  int *barray;
  bool error;
  int xsize;
  int ysize;

  barray = NULL;
/*
  Read the data.
*/
  error = pbmb_read ( file_name, &xsize, &ysize, &barray );

  if ( error )
  {
    printf ( "\n" );
    printf ( "PBMB_READ_TEST: Fatal error!\n" );
    printf ( "  PBMB_READ failed.\n" );
    if ( barray != NULL )
    {
      free ( barray );
    }
    return true;
  }
/*
  Check the data.
*/
  error = pbmb_check_data ( xsize, ysize, barray );

  if ( barray != NULL )
  {
    free ( barray );
  }

  if ( error )
  {
    printf ( "\n" );
    printf ( "  PBM_CHECK_DATA reports bad data from the file.\n" );
    return true;
  }

  printf ( "\n" );
  printf ( "  PBM_CHECK_DATA passes the data from the file.\n" );

  return false;
}
/******************************************************************************/

bool pbmb_write ( char *file_name, int xsize, int ysize, int *barray )

/******************************************************************************/
/*
  Purpose:

    PBMB_WRITE writes the header and data for a binary portable bit map file.
 
  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:
 
    04 October 1998
 
  Author:
 
    John Burkardt

  Parameters:

    Input, char *FILE_NAME, the name of the file to contain the binary
    portable bit map data.

    Input, int XSIZE, YSIZE, the number of rows and columns of data.

    Input, int *BARRAY, the array of XSIZE by YSIZE data values.

    Output, bool PBMB_WRITE, is
    true, if an error was detected, or
    false, if the file was written.
*/
{
  bool error;
  FILE *file_pointer;

  file_pointer = fopen ( file_name, "wb" );

  if ( file_pointer == NULL )
  {
    printf ( "\n" );
    printf ( "PBMB_WRITE: Fatal error!\n" );
    printf ( "  Cannot open the output file %s.\n", file_name );
    return 1;
  }
/*
  Write the header.
*/
  error = pbmb_write_header ( file_pointer, xsize, ysize );

  if ( error )
  {
    printf ( "\n" );
    printf ( "PBMB_WRITE: Fatal error!\n" );
    printf ( "  PBMB_WRITE_HEADER failed.\n" );
    return true;
  }
/*
  Write the data.
*/
  error = pbmb_write_data ( file_pointer, xsize, ysize, barray );

  if ( error )
  {
    printf ( "\n" );
    printf ( "PBMB_WRITE: Fatal error!\n" );
    printf ( "  PBMB_WRITE_DATA failed.\n" );
    return true;
  }
/*
  Close the file.
*/
  fclose ( file_pointer );

  return false;
}
/******************************************************************************/

bool pbmb_write_data ( FILE *file_pointer, int xsize, int ysize, int *barray )

/******************************************************************************/
/*
  Purpose:

    PBMB_WRITE_DATA writes the data for a binary portable bit map file.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    04 October 1998

  Author:

    John Burkardt

  Parameters:

    Input, FILE *FILE_POINTER, a pointer to the file to contain the binary
    portable bit map data.

    Input, int XSIZE, YSIZE, the number of rows and columns of data.

    Input, int *BARRAY, the array of XSIZE by YSIZE data values.

    Output, bool PBMB_WRITE_DATA, is
    true, if an error was detected, or
    false, if the data was written.
*/
{
  int bit;
  unsigned char c;
  int i;
  int *indexb;
  int j;
  int k;

  indexb = barray;
  c = 0;

  for ( j = 0; j < ysize; j++ )
  {
    for ( i = 0; i < xsize; i++ )
    {
      k = 7 - i%8;
      bit = (*indexb)%2;
      c = c | ( bit << k );

      indexb = indexb + 1;

      if ( (i+1)%8 == 0 || i == ( xsize - 1 ) )
      {
        fputc ( c, file_pointer );
        c = 0;
      }

    }
  }
  return false;
}
/******************************************************************************/

bool pbmb_write_header ( FILE *file_pointer, int xsize, int ysize )

/******************************************************************************/
/*
  Purpose:

    PBMB_WRITE_HEADER writes the header of a binary portable bit map file.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    04 October 1998

  Author:

    John Burkardt

  Parameters:

    Input, FILE *FILE_POINTER, a pointer to the file to contain the binary
    portable bit map data.

    Input, int XSIZE, YSIZE, the number of rows and columns of data.

    Output, bool PBMB_WRITE_HEADER, is
    true, if an error was detected, or
    false, if the header was written.
*/
{
  fprintf ( file_pointer, "P4 %d %d ", xsize, ysize );

  return false;
}
/******************************************************************************/

bool pbmb_write_test ( char *file_name )

/******************************************************************************/
/*
  Purpose:

    PBMB_WRITE_TEST tests the binary portable bit map write routines.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 December 2002

  Author:

    John Burkardt

  Parameters:

    Input, char *FILE_NAME, the name of the file to contain the binary
    portable bit map data.

    Output, bool PBMB_WRITE_TEST, is
    true, if an error was detected, or
    false, if the test was carried out.
*/
{
  int *barray;
  bool error;
  int xsize;
  int ysize;
/*
  Set the data.
*/
  xsize = 250;
  ysize = 150;
 
  barray = ( int * ) malloc ( xsize * ysize * sizeof ( int ) );

  if ( barray == NULL )
  {
    printf ( "\n" );
    printf ( "PBMB_WRITE_TEST: Fatal error!\n" );
    printf ( "  Unable to allocate memory for data.\n" );
    return true;
  }

  error = pbmb_example ( xsize, ysize, barray );

  if ( error )
  {
    printf ( "\n" );
    printf ( "PBMB_WRITE_TEST: Fatal error!\n" );
    printf ( "  PBM_EXAMPLE failed.\n" );
    return true;
  }

  error = pbmb_write ( file_name, xsize, ysize, barray );

  if ( barray != NULL )
  {
    free ( barray );
  }

  if ( error )
  {
    printf ( "\n" );
    printf ( "PBMB_WRITE_TEST: Fatal error!\n" );
    printf ( "  PBMB_WRITE failed.\n" );
    return true;
  }

  return false;
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