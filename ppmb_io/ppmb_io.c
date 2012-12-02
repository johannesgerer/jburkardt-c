# include <stdlib.h>
# include <stdio.h>
# include <stdbool.h>
# include <string.h>
# include <ctype.h>
# include <math.h>
# include <time.h>

# include "ppmb_io.h"

/******************************************************************************/

bool ppmb_check_data ( int xsize, int ysize, int maxrgb, unsigned char *rarray,
  unsigned char *garray, unsigned char *barray )

/******************************************************************************/
/*
  Purpose:

    PPMB_CHECK_DATA checks the data for an ASCII portable pixel map file.

  Example:

    P3
    # feep.ppm
    4 4
    15
     0  0  0    0  0  0    0  0  0   15  0 15
     0  0  0    0 15  7    0  0  0    0  0  0
     0  0  0    0  0  0    0 15  7    0  0  0
    15  0 15    0  0  0    0  0  0    0  0  0

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 June 2012

  Author:

    John Burkardt

  Parameters:

    Input, int XSIZE, YSIZE, the number of rows and columns of data.

    Input, int MAXRGB, the maximum RGB value.

    Input, unsigned char *RARRAY, *GARRAY, *BARRAY, the arrays of XSIZE by YSIZE 
    data values.

    Output, bool PPMB_CHECK_DATA, is
    true, if an error was detected, or
    false, if the data was legal.
*/
{
  int i;
  unsigned char *index;
  int j;
  int k;

  if ( xsize <= 0 )
  {
    printf ( "\n" );
    printf ( "PPMb_CHECK_DATA: 0 >= XSIZE = %d.\n", xsize );
    return true;
  }

  if ( ysize <= 0 )
  {
    printf ( "\n" );
    printf ( "PPMB_CHECK_DATA: 0 >= YSIZE = %d.\n", ysize );
    return true;
  }

  if ( rarray == NULL || garray == NULL || barray == NULL )
  {
    printf ( "\n" );
    printf ( "PPMB_CHECK_DATA: Null pointer to data.\n" );
    return true;
  }

  for ( k = 0; k < 3; k++ )
  {
    if ( k == 0 )
    {
      index = rarray;
    }
    else if ( k == 1 )
    {
      index = garray;
    }
    else if ( k == 2 )
    {
      index = barray;
    }

    for ( j = 0; j < ysize; j++ )
    {
      for ( i = 0; i < xsize; i++ )
      {
        if ( *index < 0 )
        {
          if ( k == 0 )
          {
            printf ( "\n" );
            printf ( "PPMB_CHECK_DATA: R(%d,%d) = %d < 0.\n", i, j, *index );
          }
          else if ( k == 1 )
          {
            printf ( "\n" );
            printf ( "PPMB_CHECK_DATA: G(%d,%d) = %d < 0.\n", i, j, *index );
          }
          else if ( k == 2 )
          {
            printf ( "\n" );
            printf ( "PPMB_CHECK_DATA: B(%d,%d) = %d < 0.\n", i, j, *index );
          }
          return true;
        }
        else if ( *index > maxrgb )
        {
          if ( k == 0 )
          {
            printf ( "\n" );
            printf ( "PPMB_CHECK_DATA: R(%d,%d) = %d > %d.\n", i, j, *index, 
              maxrgb );
          }
          else if ( k == 1 )
          {
            printf ( "\n" );
            printf ( "PPMB_CHECK_DATA: G(%d,%d) = %d > %d.\n", i, j, *index, 
              maxrgb );
          }
          else if ( k == 2 )
          {
            printf ( "\n" );
            printf ( "PPMB_CHECK_DATA: B(%d,%d) = %d > %d.\n", i, j, *index, 
              maxrgb );
          }
          return true;
        }

        index = index + 1;
      }
    } 
  }

  return false;
}
/******************************************************************************/

bool ppmb_example ( int xsize, int ysize, unsigned char *rarray, 
  unsigned char *garray, unsigned char *barray )

/******************************************************************************/
/*
  Purpose:

    PPMB_EXAMPLE sets up some PPM data.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 June 2012

  Author:

    John Burkardt

  Parameters:

    Input, int XSIZE, YSIZE, the number of rows and columns of data.
    Values of 200 would be reasonable.

    Output, unsigned char *RARRAY, *GARRAY, *BARRAY, the arrays of XSIZE by YSIZE 
    RGB values.

    Output, bool PPMB_EXAMPLE, is
    true, if an error occurred,
    false, if no error occurred.
*/
{
  float f1;
  float f2;
  float f3;
  int i;
  unsigned char *indexr;
  unsigned char *indexg;
  unsigned char *indexb;
  int j;
  float x;
  float y;

  indexr = rarray;
  indexg = garray;
  indexb = barray;

  for ( i = 0; i < ysize; i++ )
  {
    y = ( float ) ( ysize + 1 - i ) / ( float ) ( ysize - 1 );
    for ( j = 0; j < xsize; j++ )
    {
      x = ( float ) ( j ) / ( float ) ( xsize - 1 );

      f1 = 4.0 * ( x - 0.5 ) * ( x - 0.5 );
      f2 = sin ( 3.14159265 * x );
      f3 = x;

      if ( y <= f1 )
      {
        *indexr = ( int ) ( 255.0 * f1 );
      }
      else
      {
        *indexr = 50;
      }

      if ( y <= f2 )
      {
        *indexg = ( int ) ( 255.0 * f2 );
      }
      else
      {
        *indexg = 150;
      }

      if ( y <= f3 )
      {
        *indexb = ( int ) ( 255.0 * f3 );
      }
      else
      {
        *indexb = 250;
      }

      indexr = indexr + 1;
      indexg = indexg + 1;
      indexb = indexb + 1;
    }
  }

  return false;
}
/******************************************************************************/

bool ppmb_read ( char *file_name, int *xsize, int *ysize, int *maxrgb,
  unsigned char **rarray, unsigned char **garray, unsigned char **barray )

/******************************************************************************/
/*
  Purpose:

    PPMB_READ reads the header and data from a binary portable pixel map file.
 
  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:
 
    16 June 2012
 
  Author:
 
    John Burkardt

  Parameters:

    Input, char *FILE_NAME, the name of the file containing the binary
    portable pixel map data.

    Output, int *XSIZE, *YSIZE, the number of rows and columns of data.

    Output, int *MAXRGB, the maximum RGB value.

    Output, unsigned char **RARRAY, **GARRAY, **BARRAY, the arrays of XSIZE by YSIZE 
    data values.

    Output, bool PPMB_READ, equals
    true, if the file could not be read,
    false, if the file was read.
*/
{
  FILE *file_pointer;
  int numbytes;
  bool result;

  file_pointer = fopen ( file_name, "rb" );

  if ( file_pointer == NULL )
  {
    printf ( "\n" );
    printf ( "PPMB_READ: Fatal error!\n" );
    printf ( "  Cannot open the input file %s.\n", file_name );
    return true;
  }
/*
  Read the header.
*/
  result = ppmb_read_header ( file_pointer, xsize, ysize, maxrgb );

  if ( result )
  {
    printf ( "\n" );
    printf ( "PPMB_READ: Fatal error!\n" );
    printf ( "  PPMB_READ_HEADER failed.\n" );
    return true;
  }
/*
  Allocate storage for the data.
*/
  numbytes = ( *xsize ) * ( *ysize ) * sizeof ( unsigned char );

  *rarray = ( unsigned char * ) malloc ( numbytes );

  if ( *rarray == NULL )
  {
    printf ( "\n" );
    printf ( "PPMB_READ: Fatal error!\n" );
    printf ( "  Unable to allocate memory for data.\n" );
    printf ( "  Seeking %d bytes.\n", numbytes );
    return true;
  }

  *garray = ( unsigned char * ) malloc ( numbytes );

  if ( *garray == NULL )
  {
    printf ( "\n" );
    printf ( "PPMB_READ: Fatal error!\n" );
    printf ( "  Unable to allocate memory for data.\n" );
    printf ( "  Seeking %d bytes.\n", numbytes );
    return true;
  }

  *barray = ( unsigned char * ) malloc ( numbytes );

  if ( *barray == NULL )
  {
    printf ( "\n" );
    printf ( "PPMB_READ: Fatal error!\n" );
    printf ( "  Unable to allocate memory for data.\n" );
    printf ( "  Seeking %d bytes.\n", numbytes );
    return true;
  }
/*
  Read the data.
*/
  result = ppmb_read_data ( file_pointer, *xsize, *ysize, *rarray, 
    *garray, *barray );

  if ( result )
  {
    printf ( "\n" );
    printf ( "PPMB_READ: Fatal error!\n" );
    printf ( "  PPMB_READ_DATA failed.\n" );
    return true;
  }
/*
  Close the file.
*/
  fclose ( file_pointer );

  return false;
}
/******************************************************************************/

bool ppmb_read_data ( FILE *file_pointer, int xsize, int ysize, 
  unsigned char *rarray, unsigned char *garray, unsigned char *barray )

/******************************************************************************/
/*
  Purpose:

    PPMB_READ_DATA reads the data in a binary portable pixel map file.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 June 2012

  Author:

    John Burkardt

  Parameters:

    Input, FILE *FILE_POINTER, a pointer to the file containing the binary
    portable pixel map data.

    Input, int XSIZE, YSIZE, the number of rows and columns of data.

    Input, unsigned char *RARRAY, *GARRAY, *BARRAY, the arrays of XSIZE by YSIZE 
    data values.

    Output, bool PPMB_READ_DATA, equals
    true, if the data could not be read,
    false, if the data was read.
*/
{
  int i;
  char c;
  unsigned char *indexb;
  unsigned char *indexg;
  unsigned char *indexr;
  int j;
  int k;
  int numval;

  indexr = rarray;
  indexg = garray;
  indexb = barray;
  numval = 0;

  for ( j = 0; j < ysize; j++ )
  {
    for ( i = 0; i < xsize; i++ )
    {
      for ( k = 0; k < 3; k++ )
      {
        c = fgetc ( file_pointer );

        if ( c == EOF )
        {
          printf ( "\n" );
          printf ( "PPMB_READ_DATA: Failed reading data byte %d.\n", numval );
          return true;
        }
        else
        {
          if ( k == 0 )
          {
            *indexr = ( unsigned char ) c;
            indexr = indexr + 1;
          }
          else if ( k == 1 )
          {
            *indexg = ( unsigned char ) c;
            indexg = indexg + 1;
          }
          else if ( k == 2 )
          {
            *indexb = ( unsigned char ) c;
            indexb = indexb + 1;
          }
        }
        numval = numval + 1;
      }
    }
  }
  return false;
}
/******************************************************************************/

bool ppmb_read_header ( FILE *file_pointer, int *xsize, int *ysize, int *maxrgb )

/******************************************************************************/
/*
  Purpose:

    PPMB_READ_HEADER reads the header of a binary portable pixel map file.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 June 2012

  Author:

    John Burkardt

  Parameters:

    Input, FILE *FILE_POINTER, a pointer to the file containing the binary
    portable pixel map data.

    Output, int *XSIZE, *YSIZE, the number of rows and columns of data.

    Output, int *MAXRGB, the maximum RGB value.

    Output, bool PPMB_READ_HEADER, equals
    true, if the header could not be read,
    false, if the header was read.
*/
{
  int c_val;
  int count;
  int flag;
  int nchar;
  int state;
  char string[255];

  state = 0;
  nchar = 0;

  for ( ;; )
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
        if ( strcmp ( string, "P6" ) != 0 && strcmp ( string, "p6" ) != 0 )
        {
          printf ( "\n" );
          printf ( "PPMB_READ_HEADER: Fatal error.\n" );
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
        count = sscanf ( string, "%d", maxrgb );
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

bool ppmb_read_test ( char *file_name )

/******************************************************************************/
/*
  Purpose:

    PPMB_READ_TEST tests the binary portable pixel map read routines.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 June 2012

  Author:

    John Burkardt

  Parameters:

    Input, char *FILE_NAME, the name of the file containing the binary
    portable pixel map data.

    Output, bool PPMB_READ_TEST, equals
    true, if the test could not be carried out,
    false, if the test was carried out.
*/
{
  unsigned char *barray;
  unsigned char *garray;
  int maxrgb;
  unsigned char *rarray;
  bool result;
  int xsize;
  int ysize;

  rarray = NULL;
  garray = NULL;
  barray = NULL;
/*
  Read the data.
*/
  result = ppmb_read ( file_name, &xsize, &ysize, &maxrgb, &rarray,
    &garray, &barray );

  if ( result )
  {
    printf ( "\n" );
    printf ( "PPMB_READ_TEST: Fatal error!\n" );
    printf ( "  PPMB_READ failed.\n" );
    if ( rarray != NULL )
    {
      free ( rarray );
    }
    if ( garray != NULL )
    {
      free ( garray );
    }
    if ( barray != NULL )
    {
      free ( barray );
    }
    return true;
  }
/*
  Check the data.
*/
  result = ppmb_check_data ( xsize, ysize, maxrgb, rarray, garray, barray );

  if ( rarray != NULL )
  {
    free ( rarray );
  }
  if ( garray != NULL )
  {
    free ( garray );
  }
  if ( barray != NULL )
  {
    free ( barray );
  }

  if ( result )
  {
    printf ( "\n" );
    printf ( "  PPMB_CHECK_DATA reports bad data from the file.\n" );
    return true;
  }

  printf ( "\n" );
  printf ( "  PPMB_CHECK_DATA passes the data from the file.\n" );

  return false;
}
/******************************************************************************/

bool ppmb_write ( char *file_name, int xsize, int ysize, unsigned char *rarray, 
  unsigned char *garray, unsigned char *barray )

/******************************************************************************/
/*
  Purpose:

    PPMB_WRITE writes the header and data for a binary portable pixel map file.
 
   Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:
 
    16 June 2012
 
  Author:
 
    John Burkardt

  Parameters:

    Input, char *FILE_NAME, the name of the file to contain the binary
    portable pixel map data.

    Input, int XSIZE, YSIZE, the number of rows and columns of data.

    Input, unsigned char *RARRAY, *GARRAY, *BARRAY, the arrays of XSIZE by YSIZE 
    data values.

    Output, bool PPMB_WRITE, equals
    true, if the file could not be written,
    false, if the file was written.
*/
{
  FILE *file_pointer;
  int i;
  unsigned char *indexb;
  unsigned char *indexg;
  unsigned char *indexr;
  int j;
  int maxrgb;
  bool result;
/*
  Open the output file.
*/
  file_pointer = fopen ( file_name, "wb" );

  if ( file_pointer == NULL )
  {
    printf ( "\n" );
    printf ( "PPMB_WRITE: Fatal error!\n" );
    printf ( "  Cannot open the output file %s.\n", file_name );
    return true;
  }
/*
  Compute the maximum.
*/
  maxrgb = 0;
  indexr = rarray;
  indexg = garray;
  indexb = barray;

  for ( j = 0; j < ysize; j++ )
  {
    for ( i = 0; i < xsize; i++ )
    {
      if ( maxrgb < *indexr )
      {
        maxrgb = *indexr;
      }
      if ( maxrgb < *indexg )
      {
        maxrgb = *indexg;
      }
      if ( maxrgb < *indexb )
      {
        maxrgb = *indexb;
      }
      indexr = indexr + 1;
      indexg = indexg + 1;
      indexb = indexb + 1;
    }
  }
/*
  Write the header.
*/
  result = ppmb_write_header ( file_pointer, xsize, ysize, maxrgb );

  if ( result )
  {
    printf ( "\n" );
    printf ( "PPMB_WRITE: Fatal error!\n" );
    printf ( "  PPMB_WRITE_HEADER failed.\n" );
    return true;
  }
/*
  Write the data.
*/
  result = ppmb_write_data ( file_pointer, xsize, ysize, rarray, garray, 
    barray );

  if ( result )
  {
    printf ( "\n" );
    printf ( "PPMB_WRITE: Fatal error!\n" );
    printf ( "  PPMB_WRITE_DATA failed.\n" );
    return true;
  }
/*
  Close the file.
*/
  fclose ( file_pointer );

  return false;
}
/******************************************************************************/

bool ppmb_write_data ( FILE *file_pointer, int xsize, int ysize, 
  unsigned char *rarray, unsigned char *garray, unsigned char *barray )

/******************************************************************************/
/*
  Purpose:

    PPMB_WRITE_DATA writes the data for a binary portable pixel map file.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 June 2012

  Author:

    John Burkardt

  Parameters:

    Input, FILE *FILE_POINTER, a pointer to the file to contain the binary
    portable pixel map data.

    Input, int XSIZE, YSIZE, the number of rows and columns of data.

    Input, unsigned char *RARRAY, *GARRAY, *BARRAY, the arrays of XSIZE by YSIZE 
    data values.

    Output, bool PPMB_WRITE_DATA, equals
    true, if the data could not be written,
    false, if the data was written.
*/
{
  int i;
  unsigned char *indexb;
  unsigned char *indexg;
  unsigned char *indexr;
  int j;

  indexr = rarray;
  indexg = garray;
  indexb = barray;

  for ( j = 0; j < ysize; j++ )
  {
    for ( i = 0; i < xsize; i++ )
    {
      fputc ( *indexr, file_pointer );
      fputc ( *indexg, file_pointer );
      fputc ( *indexb, file_pointer );
      indexr = indexr + 1;
      indexg = indexg + 1;
      indexb = indexb + 1;
    }
  }
  return false;
}
/******************************************************************************/

bool ppmb_write_header ( FILE *file_pointer, int xsize, int ysize, int maxrgb )

/******************************************************************************/
/*
  Purpose:

    PPMB_WRITE_HEADER writes the header of a binary portable pixel map file.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 June 2012

  Author:

    John Burkardt

  Parameters:

    Input, FILE *FILE_POINTER, a pointer to the file to contain the binary
    portable pixel map data.

    Input, int XSIZE, YSIZE, the number of rows and columns of data.

    Input, int MAXRGB, the maximum RGB value.

    Output, bool PPMB_WRITE_HEADER, equals
    true, if the header could not be written,
    false, if the header was written.
*/
{
  fprintf ( file_pointer, "P6 %d %d %d ", xsize, ysize, maxrgb );

  return false;
}
/******************************************************************************/

bool ppmb_write_test ( char *file_name )

/******************************************************************************/
/*
  Purpose:

    PPMB_WRITE_TEST tests the binary portable pixel map write routines.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 June 2012

  Author:

    John Burkardt

  Parameters:

    Input, char *FILE_NAME, the name of the file to contain the binary
    portable pixel map data.

    Output, bool PPMB_WRITE_TEST equals
    true, if the test could not be carried out,
    false, if the test was carried out.
*/
{
  unsigned char *barray;
  unsigned char *garray;
  int  maxrgb;
  unsigned char *rarray;
  bool result;
  int xsize;
  int ysize;

  xsize = 200;
  ysize = 200;
/*
  Allocate memory.
*/ 
  rarray = ( unsigned char * ) malloc ( xsize * ysize * sizeof ( unsigned char ) );

  if ( rarray == NULL )
  {
    printf ( "\n" );
    printf ( "PPMB_WRITE_TEST: Fatal error!\n" );
    printf ( "  Unable to allocate memory for data.\n" );
    return true;
  }

  garray = ( unsigned char * ) malloc ( xsize * ysize * sizeof ( unsigned char ) );

  if ( garray == NULL )
  {
    printf ( "\n" );
    printf ( "PPMB_WRITE_TEST: Fatal error!\n" );
    printf ( "  Unable to allocate memory for data.\n" );
    return true;
  }

  barray = ( unsigned char * ) malloc ( xsize * ysize * sizeof ( unsigned char ) );

  if ( barray == NULL )
  {
    printf ( "\n" );
    printf ( "PPMB_WRITE_TEST: Fatal error!\n" );
    printf ( "  Unable to allocate memory for data.\n" );
    return true;
  }
/*
  Set the data.
*/
  result = ppmb_example ( xsize, ysize, rarray, garray, barray );

  if ( result )
  {
    printf ( "\n" );
    printf ( "PPMB_WRITE_TEST: Fatal error!\n" );
    printf ( "  PPM_EXAMPLE failed.\n" );
    return true;
  }
/*
  Write the data to the file.
*/
  result = ppmb_write ( file_name, xsize, ysize, rarray, garray, barray );

  if ( rarray != NULL )
  {
    free ( rarray );
  }

  if ( garray != NULL )
  {
    free ( garray );
  }

  if ( barray != NULL )
  {
    free ( barray );
  }

  if ( result )
  {
    printf ( "\n" );
    printf ( "PPMB_WRITE_TEST: Fatal error!\n" );
    printf ( "  PPMB_WRITE failed.\n" );
    return true;
  }

  return false;
}
