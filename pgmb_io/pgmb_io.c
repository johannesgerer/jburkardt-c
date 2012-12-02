# include <stdlib.h>
# include <stdio.h>
# include <stdbool.h>
# include <string.h>
# include <ctype.h>
# include <math.h>
# include <time.h>

# include "pgmb_io.h"

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

bool pgmb_check_data ( int xsize, int ysize, unsigned char maxgray, 
  unsigned char *garray )

/******************************************************************************/
/*
  Purpose:

    PGMB_CHECK_DATA checks the data for a binary portable gray map file.

  Example:

    P2
    # feep.pgm
    24 7
    15
    0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
    0  3  3  3  3  0  0  7  7  7  7  0  0 11 11 11 11  0  0 15 15 15 15  0
    0  3  0  0  0  0  0  7  0  0  0  0  0 11  0  0  0  0  0 15  0  0 15  0
    0  3  3  3  0  0  0  7  7  7  0  0  0 11 11 11  0  0  0 15 15 15 15  0
    0  3  0  0  0  0  0  7  0  0  0  0  0 11  0  0  0  0  0 15  0  0  0  0
    0  3  0  0  0  0  0  7  7  7  7  0  0 11 11 11 11  0  0 15  0  0  0  0
    0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 June 2012

  Author:

    John Burkardt

  Parameters:

    Input, int XSIZE, YSIZE, the number of rows and columns of data.

    Input, unsigned chart MAXGRAY, the maximum gray value.

    Input, unsigned char *GARRAY, the array of XSIZE by YSIZE data values.

    Output, bool PGMB_CHECK_DATA, is
    true, if an error was detected, or
    false, if the data was legal.
*/
{
  int i;
  unsigned char *indexg;
  int j;

  if ( xsize <= 0 )
  {
    printf ( "\n" );
    printf ( "PGMB_CHECK_DATA: 0 >= XSIZE = %d.\n", xsize );
    return true;
  }
  if ( ysize <= 0 )
  {
    printf ( "\n" );
    printf ( "PGMB_CHECK_DATA: 0 >= YSIZE = %d.\n", ysize );
    return true;
  }

  if ( garray == NULL )
  {
    printf ( "\n" );
    printf ( "PGMB_CHECK_DATA: Null pointer to data.\n" );
    return true;
  }

  indexg = garray;

  for ( j = 0; j < ysize; j++ )
  {
    for ( i = 0; i < xsize; i++ )
    {
      if ( *indexg < 0  )
      {
        printf ( "\n" );
        printf ( "PGMB_CHECK_DATA: G(%d,%d) = %d < 0.\n", 
          i, j, *indexg );
        return true;
      }
      else if ( maxgray < *indexg )
      {
        printf ( "\n" );
        printf ( "PGMB_CHECK_DATA: G(%d,%d) = %d > %d.\n", 
          i, j, *indexg, maxgray );
        return true;
      }
      indexg = indexg + 1;
    }
  }

  return false;
}
/******************************************************************************/

bool pgmb_example ( int xsize, int ysize, unsigned char *garray )

/******************************************************************************/
/*
  Purpose:

    PGMB_EXAMPLE sets up some PGMB data.

  Discussion:

    The data is based on three periods of a sine curve.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 June 2012

  Author:

    John Burkardt

  Parameters:

    Input, int XSIZE, YSIZE, the number of rows and columns of data.
    Values of 200 would be reasonable.

    Output, unsigned char *GARRAY, the array of XSIZE by YSIZE data values.

    Output, bool PGMB_EXAMPLE, is
    true, if an error occurred,
    false, if no error occurred.
*/
{
  int i;
  unsigned char *indexg;
  int j;
  int periods = 3;
  float pi = 3.14159265;
  float x;
  float y;

  indexg = garray;

  for ( i = 0; i < ysize; i++ )
  {
    y = 2.0 * ( float ) ( i ) / ( float ) ( ysize - 1 ) - 1.0;
    for ( j = 0; j < xsize; j++ )
    {
      x = 2.0 * pi * ( float ) ( periods * ( j ) ) / ( float ) ( xsize - 1 );
      *indexg = ( unsigned char ) ( 20.0 * ( sin ( x ) - y + 2.0 ) );
      indexg = indexg + 1;
    }
  }

  return false;
}
/******************************************************************************/

bool pgmb_read ( char *file_name, int *xsize, int *ysize, unsigned char *maxgray,
  unsigned char **garray )

/******************************************************************************/
/*
  Purpose:

    PGMB_READ reads the header and data from a binary portable gray map file.
 
  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:
 
    16 June 2012
 
  Author:
 
    John Burkardt

  Parameters:

    Input, char *FILE_NAME, the name of the file containing the binary
    portable gray map data.

    Output, int *XSIZE, *YSIZE, the number of rows and columns of data.

    Output, unsigned char *MAXGRAY, the maximum gray value.

    Output, unsigned char **GARRAY, the array of XSIZE by YSIZE data values.

    Output, bool PGMB_READ, is
    true, if an error was detected, or
    false, if the file was read.
*/
{
  FILE *file_pointer;
  int numbytes;
  int result;

  file_pointer = fopen ( file_name, "rb" );

  if ( file_pointer == NULL )
  {
    printf ( "\n" );
    printf ( "PGMB_READ: Fatal error!\n" );
    printf ( "  Cannot open the input file %s.\n", file_name );
    return true;
  }
/*
  Read the header.
*/
  result = pgmb_read_header ( file_pointer, xsize, ysize, maxgray );

  if ( result != 0 )
  {
    printf ( "\n" );
    printf ( "PGMB_READ: Fatal error!\n" );
    printf ( "  PGMB_READ_HEADER failed.\n" );
    return true;
  }
/*
  Allocate storage for the data.
*/
  numbytes = ( *xsize ) * ( *ysize ) * sizeof ( unsigned char );

  *garray = ( unsigned char * ) malloc ( numbytes );

  if ( *garray == NULL )
  {
    printf ( "\n" );
    printf ( "PGMB_READ: Fatal error!\n" );
    printf ( "  Unable to allocate memory for data.\n" );
    printf ( "  Seeking %d bytes.\n", numbytes );
    return true;
  }
/*
  Read the data.
*/
  result = pgmb_read_data ( file_pointer, *xsize, *ysize, *garray );

  if ( result != 0 )
  {
    printf ( "\n" );
    printf ( "PGMB_READ: Fatal error!\n" );
    printf ( "  PGMB_READ_DATA failed.\n" );
    return true;
  }
/*
  Close the file.
*/
  fclose ( file_pointer );

  return false;
}
/******************************************************************************/

bool pgmb_read_data ( FILE *file_pointer, int xsize, int ysize, 
  unsigned char *garray )

/******************************************************************************/
/*
  Purpose:

    PGMB_READ_DATA reads the data in a binary portable gray map file.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 June 2012

  Author:

    John Burkardt

  Parameters:

    Input, FILE *FILE_POINTER, a pointer to the file containing the binary
    portable gray map data.

    Input, int XSIZE, YSIZE, the number of rows and columns of data.

    Output, unsigned char *GARRAY, the array of XSIZE by YSIZE data values.

    Output, bool PGMB_READ_DATA, is
    true, if an error was detected, or
    false, if the data was read.
*/
{
  char c;
  int i;
  unsigned char *indexg;
  int j;

  indexg = garray;

  for ( j = 0; j < ysize; j++ )
  {
    for ( i = 0; i < xsize; i++ )
    {
      c = fgetc ( file_pointer );

      if ( c == EOF )
      {
        printf ( "\n" );
        printf ( "PGMB_READ_DATA: End of file reading pixed[%d][%d].\n", i,j );
        return true;
      }

      *indexg = ( unsigned char ) c;
      indexg = indexg + 1;
    }
  }
  return false;
}
/******************************************************************************/

bool pgmb_read_header ( FILE *file_pointer, int *xsize, int *ysize, 
  unsigned char *maxgray )

/******************************************************************************/
/*
  Purpose:

    PGMB_READ_HEADER reads the header of a binary portable gray map file.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 June 2012

  Author:

    John Burkardt

  Parameters:

    Input, FILE *FILE_POINTER, a pointer to the file containing the binary
    portable gray map data.

    Output, int *XSIZE, *YSIZE, the number of rows and columns of data.

    Output, unsigned char *MAXGRAY, the maximum gray value.

    Output, bool PGMB_READ_HEADER, is
    true, if an error was detected, or
    false, if the header was read.
*/
{
  char c;
  int count;
  int flag;
  int nchar;
  int state;
  char string[255];

  state = 0;
  nchar = 0;

  for ( ; ; )
  {

    c = fgetc ( file_pointer );

    if ( c == EOF )
    {
      return true;
    }
/*
  If not whitespace, add the character to the current string.
*/
    flag = isspace ( c );

    if ( !flag )
    {
      string[nchar] = c;
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
        if ( strcmp ( string, "P5" ) != 0 && strcmp ( string, "p5" ) != 0 )
        {
          printf ( "\n" );
          printf ( "PGMB_READ_HEADER: Fatal error.\n" );
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
        state = 6;
      }
    }
    else if ( state == 6 )
    {
      if ( !flag )
      {
        state = 7;
      }
    }
    else if ( state == 7 )
    {
      if ( flag )
      {
        string[nchar] = 0;
        nchar = nchar + 1;
        count = sscanf ( string, "%d", maxgray );
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

bool pgmb_read_test ( char *file_name )

/******************************************************************************/
/*
  Purpose:

    PGMB_READ_TEST tests the binary portable gray map read routines.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 June 2012

  Author:

    John Burkardt

  Parameters:

    Input, char *FILE_NAME, the name of the file containing the binary
    portable gray map data.

    Output, bool PGMB_TEST, is
    true, if an error was detected, or
    false, if the test was carried out.
*/
{
  unsigned char *garray;
  unsigned char maxgray;
  int result;
  int xsize;
  int ysize;

  garray = NULL;
/*
  Read the data.
*/
  result = pgmb_read ( file_name, &xsize, &ysize, &maxgray, &garray );

  if ( result != 0 )
  {
    printf ( "\n" );
    printf ( "PGMB_READ_TEST: Fatal error!\n" );
    printf ( "  PGMB_READ failed.\n" );
    if ( garray != NULL )
    {
      free ( garray );
    }
    return true;
  }
/*
  Check the data.
*/
  result = pgmb_check_data ( xsize, ysize, maxgray, garray );

  if ( garray != NULL )
  {
    free ( garray );
  }

  if ( result != 0 )
  {
    printf ( "\n" );
    printf ( "  PGM_CHECK_DATA reports bad data from the file.\n" );
    return true;
  }

  printf ( "\n" );
  printf ( "  PGM_CHECK_DATA passes the data from the file.\n" );

  return false;
}
/******************************************************************************/

bool pgmb_write ( char *file_name, int xsize, int ysize, unsigned char *garray )

/******************************************************************************/
/*
  Purpose:

    PGMB_WRITE writes the header and data for a binary portable gray map file.
 
  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:
 
    16 June 2012
 
  Author:
 
    John Burkardt

  Parameters:

    Input, char *FILE_NAME, the name of the file to contain the binary
    portable gray map data.

    Input, int XSIZE, YSIZE, the number of rows and columns of data.

    Input, unsigned char *GARRAY, the array of XSIZE by YSIZE data values.

    Output, bool PGMB_WRITE, is
    true, if an error was detected, or
    false, if the file was written.
*/
{
  FILE *file_pointer;
  int i;
  unsigned char *indexg;
  int j;
  unsigned char maxgray;
  int result;

  maxgray = 0;
  indexg = garray;

  for ( i = 0; i < xsize; i++ )
  {
    for ( j = 0; j < ysize; j++ )
    {
      if ( maxgray < *indexg )
      {
        maxgray = *indexg;
      }
      indexg = indexg + 1;
    }
  }

  file_pointer = fopen ( file_name, "wb" );

  if ( file_pointer == NULL )
  {
    printf ( "\n" );
    printf ( "PGMB_WRITE: Fatal error!\n" );
    printf ( "  Cannot open the output file %s.\n", file_name );
    return true;
  }
/*
  Write the header.
*/
  result = pgmb_write_header ( file_pointer, xsize, ysize, maxgray );

  if ( result != 0 )
  {
    printf ( "\n" );
    printf ( "PGMB_WRITE: Fatal error!\n" );
    printf ( "  PGMB_WRITE_HEADER failed.\n" );
    return true;
  }
/*
  Write the data.
*/
  result = pgmb_write_data ( file_pointer, xsize, ysize, garray );

  if ( result != 0 )
  {
    printf ( "\n" );
    printf ( "PGMB_WRITE: Fatal error!\n" );
    printf ( "  PGMB_WRITE_DATA failed.\n" );
    return true;
  }
/*
  Close the file.
*/
  fclose ( file_pointer );

  return false;
}
/******************************************************************************/

bool pgmb_write_data ( FILE *file_pointer, int xsize, int ysize, 
  unsigned char *garray )

/******************************************************************************/
/*
  Purpose:

    PGMB_WRITE_DATA writes the data for a binary portable gray map file.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 June 2012

  Author:

    John Burkardt

  Parameters:

    Input, FILE *FILE_POINTER, a pointer to the file to contain the binary
    portable gray map data.

    Input, int XSIZE, YSIZE, the number of rows and columns of data.

    Input, unsigned char *GARRAY, the array of XSIZE by YSIZE data values.

    Output, bool PGMB_WRITE_DATA, is
    true, if an error was detected, or
    false, if the data was written.
*/
{
  int i;
  unsigned char *indexg;
  int j;

  indexg = garray;

  for ( j = 0; j < ysize; j++ )
  {
    for ( i = 0; i < xsize; i++ )
    {
      fputc ( *indexg, file_pointer );
      indexg = indexg + 1;
    }
  }

  return false;
}
/******************************************************************************/

bool pgmb_write_header ( FILE *file_pointer, int xsize, int ysize, 
  unsigned char maxgray )

/******************************************************************************/
/*
  Purpose:

    PGMB_WRITE_HEADER writes the header of a binary portable gray map file.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 June 2012

  Author:

    John Burkardt

  Parameters:

    Input, FILE *FILE_POINTER, a pointer to the file to contain the binary
    portable gray map data.

    Input, int XSIZE, YSIZE, the number of rows and columns of data.

    Input, unsigned char MAXGRAY, the maximum gray value.

    Output, bool PGMB_WRITE_HEADER, is
    true, if an error was detected, or
    false, if the header was written.
*/
{
  fprintf ( file_pointer, "P5 %d %d %d ", xsize, ysize, maxgray );

  return false;
}
/******************************************************************************/

bool pgmb_write_test ( char *file_name )

/******************************************************************************/
/*
  Purpose:

    PGMB_WRITE_TEST tests the binary portable gray map write routines.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 June 2012

  Author:

    John Burkardt

  Parameters:

    Input, char *FILE_NAME, the name of the file to contain the binary
    portable gray map data.

    Output, bool PGMB_WRITE_TEST, is
    true, if an error was detected, or
    false, if the test was carried out.
*/
{
  unsigned char *garray;
  int result;
  int xsize;
  int ysize;
/*
  Set the data.
*/
  xsize = 300;
  ysize = 200;

  garray = ( unsigned char * ) malloc ( xsize * ysize * sizeof ( unsigned char ) );

  if ( garray == NULL )
  {
    printf ( "\n" );
    printf ( "PGMB_WRITE_TEST: Fatal error!\n" );
    printf ( "  Unable to allocate memory for data.\n" );
    return true;
  }
  result = pgmb_example ( xsize, ysize, garray );

  if ( result != 0 )
  {
    printf ( "\n" );
    printf ( "PGMB_WRITE_TEST: Fatal error!\n" );
    printf ( "  PGMB_EXAMPLE failed.\n" );
    return true;
  }

  result = pgmb_write ( file_name, xsize, ysize, garray );

  if ( garray != NULL )
  {
    free ( garray );
  }

  if ( result != 0 )
  {
    printf ( "\n" );
    printf ( "PGMB_WRITE_TEST: Fatal error!\n" );
    printf ( "  PGMB_WRITE failed.\n" );
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
