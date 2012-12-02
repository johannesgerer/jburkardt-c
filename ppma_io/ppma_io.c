# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>

# include "ppma_io.h"

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

int ppma_check_data ( int xsize, int ysize, int rgb_max, int *r, int *g, 
  int *b )

/******************************************************************************/
/*
  Purpose:

    PPMA_CHECK_DATA checks the data for an ASCII portable pixel map file.

  Discussion:

    XSIZE and YSIZE must be positive, the pointers must not be null,
    and the data must be nonnegative and no greater than RGB_MAX.

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

    27 May 2011

  Author:

    John Burkardt

  Parameters:

    Input, int XSIZE, YSIZE, the number of rows and columns of data.

    Input, int RGB_MAX, the maximum RGB value.

    Input, int *R, *G, *B, the arrays of XSIZE by YSIZE data values.

    Output, int PPMA_CHECK_DATA, is
    1, if an error was detected, or
    0, if the data was legal.
*/
{
  char c;
  int i;
  int *index;
  int j;
  int k;

  if ( xsize <= 0 )
  {
    printf ( "\n" );
    printf ( "PPMA_CHECK_DATA: Fatal error!\n" );
    printf ( "  xsize <= 0.\n" );
    printf ( "  xsize = %d\n", xsize );
    return 1;
  }

  if ( ysize <= 0 )
  {
    printf ( "\n" );
    printf ( "PPMA_CHECK_DATA: Fatal error!\n" );
    printf ( "  ysize <= 0.\n" );
    printf ( "  ysize = %d\n", ysize );
    return 1;
  }

  if ( r == NULL )
  {
    printf ( "\n" );
    printf ( "PPMA_CHECK_DATA: Fatal error!\n" );
    printf ( "  Null pointer to R.\n" );
    return 1;
  }

  if ( g == NULL )
  {
    printf ( "\n" );
    printf ( "PPMA_CHECK_DATA: Fatal error!\n" );
    printf ( "  Null pointer to G.\n" );
    return 1;
  }

  if ( b == NULL )
  {
    printf ( "\n" );
    printf ( "PPMA_CHECK_DATA: Fatal error!\n" );
    printf ( "  Null pointer to B.\n" );
    return 1;
  }

  for ( k = 0; k < 3; k++ )
  {

    if ( k == 0 )
    {
      index = r;
      c = 'R';
    }
    else if ( k == 1 )
    {
      index = g;
      c = 'G';
    }
    else if ( k == 2 )
    {
      index = b;
      c = 'B';
    }

    for ( j = 0; j < ysize; j++ )
    {
      for ( i = 0; i < xsize; i++ )
      {
        if ( *index < 0 )
        {
          printf ( "\n" );
          printf ( "PPMA_CHECK_DATA - Fatal error!\n" );
          printf ( "  Negative data.\n" );
          printf ( "  %c(%d,%d) = %d\n", c, i, j, *index );
          return 1;
        }
        else if ( rgb_max < *index )
        {
          printf ( "\n" );
          printf ( "PPMA_CHECK_DATA - Fatal error!\n" );
          printf ( "  Data exceeds RGB_MAX = %d\n", rgb_max );
          printf ( "  %c(%d,%d) = %d\n", c, i, j, *index );
          return 1;
        }

        index = index + 1;
      }
    } 
  }

  return 0;
}
/******************************************************************************/

void ppma_example ( int xsize, int ysize, int *r, int *g, int *b )

/******************************************************************************/
/*
  Purpose:

    PPMA_EXAMPLE sets up some RGB data.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    27 May 2011

  Author:

    John Burkardt

  Parameters:

    Input, int XSIZE, YSIZE, the number of rows and columns of data.

    Output, int *R, *G, *B, the arrays of XSIZE by YSIZE RGB values.
*/
{
  int *b_index;
  float f1;
  float f2;
  float f3;
  int *g_index;
  int i;
  int j;
  int *r_index;
  float x;
  float y;

  r_index = r;
  g_index = g;
  b_index = b;

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
        *r_index = ( int ) ( 255.0 * f1 );
      }
      else
      {
        *r_index = 50;
      }

      if ( y <= f2 )
      {
        *g_index = ( int ) ( 255.0 * f2 );
      }
      else
      {
        *g_index = 150;
      }

      if ( y <= f3 )
      {
        *b_index = ( int ) ( 255.0 * f3 );
      }
      else
      {
        *b_index = 250;
      }

      r_index = r_index + 1;
      g_index = g_index + 1;
      b_index = b_index + 1;
    }
  }

  return;
}
/******************************************************************************/

void ppma_read ( char *input_name, int *xsize, int *ysize, int *rgb_max,
  int **r, int **g, int **b )

/******************************************************************************/
/*
  Purpose:

    PPMA_READ reads the header and data from an ASCII portable pixel map file.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:
 
    27 May 2011
 
  Author:
 
    John Burkardt

  Parameters:

    Input, char *INPUT_NAME, the name of the file containing the ASCII
    portable pixel map data.

    Output, int *XSIZE, *YSIZE, the number of rows and columns of data.

    Output, int *RGB_MAX, the maximum RGB value.

    Output, int **R, **G, **B, the arrays of XSIZE by YSIZE data values.
*/
{
  int error;
  FILE *input;
  int numbytes;

  input = fopen ( input_name, "rt" );

  if ( !input )
  {
    printf ( "\n" );
    printf ( "PPMA_READ - Fatal error!\n" );
    printf ( "  Cannot open the input file \"%s\".\n", input_name );
    exit ( 1 );
  }
/*
  Read the header.
*/
  ppma_read_header ( input, xsize, ysize, rgb_max );
/*
  Allocate storage for the data.
*/
  numbytes = (*xsize) * (*ysize) * sizeof ( int );

  *r = ( int * ) malloc ( numbytes * sizeof ( int ) );
  *g = ( int * ) malloc ( numbytes * sizeof ( int ) );
  *b = ( int * ) malloc ( numbytes * sizeof ( int ) );
/*
  Read the data.
*/
  ppma_read_data ( input, *xsize, *ysize, *r, *g, *b );
/*
  Close the file.
*/
  fclose ( input );

  return;
}
/******************************************************************************/

void ppma_read_header ( FILE *input, int *xsize, int *ysize, int *rgb_max )

/******************************************************************************/
/*
  Purpose:

    PPMA_READ_HEADER reads the header of an ASCII portable pixel map file.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    27 May 2011

  Author:

    John Burkardt

  Parameters:

    Input, FILE *INPUT, a pointer to the file containing the ASCII
    portable pixel map data.

    Output, int *XSIZE, *YSIZE, the number of rows and columns of data.

    Output, int *RGB_MAX, the maximum RGB value.
*/
{
# define LINE_MAX 255

  int count;
  char *error;
  char line[LINE_MAX];
  char *next;
  int step;
  int width;
  char word[LINE_MAX];

  step = 0;

  while ( 1 )
  {
    error = fgets ( line, LINE_MAX, input );

    if ( !error )
    {
      fprintf ( stderr, "\n" );
      fprintf ( stderr, "PGMA_READ_HEADER - Fatal error!\n" );
      fprintf ( stderr, "  End of file.\n" );
      exit ( 1 );
    }

    next = line;

    if ( line[0] == '#' )
    {
      continue;
    }

    if ( step == 0 )
    {
      count = sscanf ( next, "%s%n", word, &width );
      if ( count == EOF )
      {
        continue;
      }
      next = next + width;
      if ( strcmp ( word, "P3" ) != 0 && strcmp ( word, "p3" ) != 0 )
      {
        fprintf ( stderr, "\n" );
        fprintf ( stderr, "PPMA_READ_HEADER - Fatal error.\n" );
        fprintf ( stderr, "  Bad magic number = \"%s\".\n", word );
        exit ( 1 );
      }
      step = 1;
    }

    if ( step == 1 )
    {

      count = sscanf ( next, "%d%n", xsize, &width );
      next = next + width;
      if ( count == EOF )
      {
        continue;
      }
      step = 2;
    }

    if ( step == 2 )
    {
      count = sscanf ( next, "%d%n", ysize, &width );
      next = next + width;
      if ( count == EOF )
      {
        continue;
      }
      step = 3;
    }

    if ( step == 3 )
    {
      count = sscanf ( next, "%d%n", rgb_max, &width );
      next = next + width;
      if ( count == EOF )
      {
        continue;
      }
      break;
    }

  }

  return;
# undef LINE_MAX
}

/******************************************************************************/

void ppma_read_data ( FILE *input, int xsize, int ysize, int *r,
  int *g, int *b )

/******************************************************************************/
/*
  Purpose:

    PPMA_READ_DATA reads the data in an ASCII portable pixel map file.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    27 May 2011

  Author:

    John Burkardt

  Parameters:

    Input, FILE *INPUT, a pointer to the file containing the ASCII
    portable pixel map data.

    Input, int XSIZE, YSIZE, the number of rows and columns of data.

    Output, int *R, *G, *B, the arrays of XSIZE by YSIZE data values.
*/
{
  int i;
  int j;
  int n;

  for ( j = 0; j < ysize; j++ )
  {
    for ( i = 0; i < xsize; i++ )
    {
      n = fscanf ( input, "%d %d %d", r, g, b );

      if ( n != 3 )
      {
        fprintf ( stderr, "\n" );
        fprintf ( stderr, "PPMA_READ_DATA - Fatal error!\n" );
        fprintf ( stderr, "  Unable to read data.\n" );
        exit ( 1 );
      }
      r = r + 1;
      g = g + 1;
      b = b + 1;
    }
  }
  return;
}
/******************************************************************************/

void ppma_read_test ( char *input_name )

/******************************************************************************/
/*
  Purpose:

    PPMA_READ_TEST tests the ASCII portable pixel map read routines.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    27 May 2011

  Author:

    John Burkardt

  Parameters:

    Input, char *INPUT_NAME, the name of the file containing the ASCII
    portable pixel map data.
*/
{
  int *b;
  int error;
  int *g;
  int *r;
  int rgb_max;
  int xsize;
  int ysize;

  r = NULL;
  g = NULL;
  b = NULL;
/*
  Read the data.
*/
  ppma_read ( input_name, &xsize, &ysize, &rgb_max, &r, &g, &b );

  printf ( "\n" );
  printf ( "PPMA_READ_TEST:\n" );
  printf ( "  PPMA_READ was able to read \"%s\".\n", input_name );
/*
  Check the data.
*/
  error = ppma_check_data ( xsize, ysize, rgb_max, r, g, b );

  free ( r );
  free ( g );
  free ( b );

  if ( error == 1 )
  {
    printf ( "\n" );
    printf ( "PPMA_READ_TEST - Fatal error!\n" );
    printf ( "  PPMA_CHECK_DATA reports bad data in the file.\n" );
  }
  else
  {
    printf ( "\n" );
    printf ( "PPMA_READ_TEST:\n" );
    printf ( "  PPMA_CHECK_DATA has approved the data from the file.\n" );
  }
  return;
}
/******************************************************************************/

int ppma_write ( char *file_out_name, int xsize, int ysize, int *r,
  int *g, int *b )

/******************************************************************************/
/*
  Purpose:

    PPMA_WRITE writes the header and data for an ASCII portable pixel map file.

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

    28 February 2003

  Author:

    John Burkardt

  Parameters:

    Input, char *FILE_OUT_NAME, the name of the file to contain the ASCII
    portable pixel map data.

    Input, int XSIZE, YSIZE, the number of rows and columns of data.

    Input, int *R, *G, *B, the arrays of XSIZE by YSIZE data values.

    Output, int PPMA_WRITE, is
    true, if an error was detected, or
    false, if the file was written.
*/
{
  int *b_index;
  int error;
  FILE *file_out;
  int *g_index;
  int i;
  int j;
  int *r_index;
  int rgb_max;
/*
  Open the output file.
*/
  file_out = fopen ( file_out_name, "wt" );

  if ( !file_out )
  {
    printf ( "\n" );
    printf ( "PPMA_WRITE - Fatal error!\n" );
    printf ( "  Cannot open the output file \"%s\".\n", file_out_name );
    return 1;
  }
/*
  Compute the maximum.
*/
  rgb_max = 0;
  r_index = r;
  g_index = g;
  b_index = b;

  for ( j = 0; j < ysize; j++ )
  {
    for ( i = 0; i < xsize; i++ )
    {
      if ( rgb_max < *r_index )
      {
        rgb_max = *r_index;
      }
      r_index = r_index + 1;

      if ( rgb_max < *g_index )
      {
        rgb_max = *g_index;
      }
      g_index = g_index + 1;

      if ( rgb_max < *b_index )
      {
        rgb_max = *b_index;
      }
      b_index = b_index + 1;
    }
  }
/*
  Write the header.
*/
  error = ppma_write_header ( file_out, file_out_name, xsize, ysize, rgb_max );

  if ( error )
  {
    printf ( "\n" );
    printf ( "PPMA_WRITE - Fatal error!\n" );
    printf ( "  PPMA_WRITE_HEADER failed.\n" );
    return 1;
  }
/*
  Write the data.
*/
  error = ppma_write_data ( file_out, xsize, ysize, r, g, b );

  if ( error )
  {
    printf ( "\n" );
    printf ( "PPMA_WRITE - Fatal error!\n" );
    printf ( "  PPMA_WRITE_DATA failed.\n" );
    return 1;
  }
/*
  Close the file.
*/
  fclose ( file_out );

  return 0;
}
/******************************************************************************/

int ppma_write_data ( FILE *file_out, int xsize, int ysize, int *r,
  int *g, int *b )

/******************************************************************************/
/*
  Purpose:

    PPMA_WRITE_DATA writes the data for an ASCII portable pixel map file.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    28 February 2003

  Author:

    John Burkardt

  Parameters:

    Input, ofstream &FILE_OUT, a pointer to the file to contain the ASCII
    portable pixel map data.

    Input, int XSIZE, YSIZE, the number of rows and columns of data.

    Input, int *R, *G, *B, the arrays of XSIZE by YSIZE data values.

    Output, int PPMA_WRITE_DATA, is
    true, if an error was detected, or
    false, if the data was written.
*/
{
  int *b_index;
  int *g_index;
  int i;
  int j;
  int *r_index;
  int rgb_num;

  r_index = r;
  g_index = g;
  b_index = b;
  rgb_num = 0;

  for ( j = 0; j < ysize; j++ )
  {
    for ( i = 0; i < xsize; i++ )
    {
      fprintf ( file_out, "%d  %d  %d", *r_index, *g_index, *b_index );
      rgb_num = rgb_num + 3;
      r_index = r_index + 1;
      g_index = g_index + 1;
      b_index = b_index + 1;

      if ( rgb_num % 12 == 0 || i == xsize - 1 || rgb_num == 3 * xsize * ysize )
      {
        fprintf ( file_out, "\n" );
      }
      else
      {
        fprintf ( file_out, " " );
      }
    }
  }
  return 0;
}
/******************************************************************************/

int ppma_write_header ( FILE *file_out, char *file_out_name, int xsize,
  int ysize, int rgb_max )

/******************************************************************************/
/*
  Purpose:

    PPMA_WRITE_HEADER writes the header of an ASCII portable pixel map file.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    28 February 2003

  Author:

    John Burkardt

  Parameters:

    Input, FILE *FILE_OUT, a pointer to the file to contain the ASCII
    portable pixel map data.

    Input, char *FILE_OUT_NAME, the name of the file.

    Input, int XSIZE, YSIZE, the number of rows and columns of data.

    Input, int RGB_MAX, the maximum RGB value.

    Output, int PPMA_WRITE_HEADER, is
    true, if an error was detected, or
    false, if the header was written.
*/
{
  fprintf ( file_out, "P3\n" );
  fprintf ( file_out, "# %s created by ppma_write.c.\n", file_out_name );
  fprintf ( file_out, "%d  %d\n", xsize, ysize );
  fprintf ( file_out, "%d\n", rgb_max );

  return 0;
}
/******************************************************************************/

int ppma_write_test ( char *output_name )

/******************************************************************************/
/*
  Purpose:

    PPMA_WRITE_TEST tests the ASCII portable pixel map write routines.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    27 May 2011

  Author:

    John Burkardt

  Parameters:

    Input, char *OUTPUT_NAME, the name of the file to contain the ASCII
    portable pixel map data.

    Output, int PPMA_WRITE_TEST, equals
    true, if the test could not be carried out,
    false, if the test was carried out.
*/
{
  int *b;
  int error;
  int *g;
  int *r;
  int xsize;
  int ysize;

  xsize = 300;
  ysize = 300;
/*
  Allocate memory.
*/
  r = ( int * ) malloc ( xsize * ysize * sizeof ( int ) );
  g = ( int * ) malloc ( xsize * ysize * sizeof ( int ) );
  b = ( int * ) malloc ( xsize * ysize * sizeof ( int ) );
/*
  Set the data.
*/
  ppma_example ( xsize, ysize, r, g, b );
/*
  Write the data to the file.
*/
  error = ppma_write ( output_name, xsize, ysize, r, g, b );

  free ( r );
  free ( g );
  free ( b );

  if ( error )
  {
    printf ( "\n" );
    printf ( "PPMA_WRITE_TEST - Fatal error!\n" );
    printf ( "  PPMA_WRITE failed.\n" );
    return 1;
  }

  return 0;
}

