# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "pbma_io.h"

/******************************************************************************/

void pbma_check_data ( int xsize, int ysize, int *b )

/******************************************************************************/
/*
  Purpose:

    PBMA_CHECK_DATA checks the data for an ASCII PBM file.

  Discussion:

    XSIZE and YSIZE must be positive, the pointers must not be null,
    and the data must be 0 or 1.

  Example:

    P1
    # feep.pbm
    24 7
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 1 1 1 1 0 0 1 1 1 1 0 0 1 1 1 1 0 0 1 1 1 1 0
    0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 1 0
    0 1 1 1 0 0 0 1 1 1 0 0 0 1 1 1 0 0 0 1 1 1 1 0
    0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0
    0 1 0 0 0 0 0 1 1 1 1 0 0 1 1 1 1 0 0 1 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    05 June 2010

  Author:

    John Burkardt

  Parameters:

    Input, int XSIZE, YSIZE, the number of rows and columns of data.

    Input, int *B, the array of XSIZE by YSIZE data values.
*/
{
  int i;
  int *index;
  int j;
  int k;

  if ( xsize <= 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "PBMA_CHECK_DATA: Error!\n" );
    fprintf ( stderr, "  XSIZE <= 0.\n" );
    fprintf ( stderr, "  XSIZE = %d\n", xsize );
    exit ( 1 );
  }

  if ( ysize <= 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "PBMA_CHECK_DATA: Error!\n" );
    fprintf ( stderr, "  YSIZE <= 0.\n" );
    fprintf ( stderr, "  YSIZE = %d\n", ysize );
    exit ( 1 );
  }

  if ( b == NULL )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "PBMA_CHECK_DATA: Error!\n" );
    fprintf ( stderr, "  Null pointer to B.\n" );
    exit ( 1 );
  }

  index = b;

  for ( j = 0; j < ysize; j++ )
  {
    for ( i = 0; i < xsize; i++ )
    {
      if ( *index < 0 )
      {
        fprintf ( stderr, "\n" );
        fprintf ( stderr, "PBMA_CHECK_DATA - Fatal error!\n" );
        fprintf ( stderr, "  Negative data.\n" );
        fprintf ( stderr, "  B(%d,%d)=%d\n", i, j, *index );
        exit ( 1 );
      }
      else if ( 1 < *index )
      {
        fprintf ( stderr, "\n" );
        fprintf ( stderr, "PBMA_CHECK_DATA - Fatal error!\n" );
        fprintf ( stderr, "  Data exceeds 1\n" );
        fprintf ( stderr, "  B(%d,%d)=%d\n", i, j, *index );
        exit ( 1 );
      }

      index = index + 1;
    }
  } 
  return;
}
/******************************************************************************/

void pbma_example ( int xsize, int ysize, int *b )

/******************************************************************************/
/*
  Purpose:

    PBMA_EXAMPLE sets up some ASCII PBM data.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    05 June 2010

  Author:

    John Burkardt

  Parameters:

    Input, int XSIZE, YSIZE, the number of rows and columns of data.

    Output, int *B, the array of XSIZE by YSIZE gray values.
*/
{
  int i;
  int *indexb;
  int j;
  double r;
  double test;
  double x;
  double xc;
  double y;
  double yc;
 
  indexb = b;

  if ( xsize < ysize )
  {
    r = ( double ) xsize / 3.0;
  }
  else
  {
    r = ( double ) ysize / 3.0;
  }

  xc = ( ( double ) xsize ) / 2.0;
  yc = ( ( double ) ysize ) / 2.0;

  for ( i = 0; i < ysize; i++ )
  {
    y = ( ( double ) i );
    for ( j = 0; j < xsize; j++ )
    {
      x = ( ( double ) j );
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

  return;
}
/******************************************************************************/

void pbma_read ( char *file_in_name, int *xsize, int *ysize, int **b )

/******************************************************************************/
/*
  Purpose:

    PBMA_READ reads the header and data from an ASCII PBM file.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:
 
    05 June 2010
 
  Author:
 
    John Burkardt

  Parameters:

    Input, char *FILE_IN_NAME, the name of the file.

    Output, int *XSIZE, *YSIZE, the number of rows and columns of data.

    Output, int **B, the array of XSIZE by YSIZE data values.
*/
{
  FILE *file_in;
  int numbytes;

  file_in = fopen ( file_in_name, "rt" );

  if ( !file_in )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "PBMA_READ - Fatal error!\n" );
    fprintf ( stderr, "  Cannot open the input file \"%s\".\n", file_in_name );
    exit ( 1 );
  }
/*
  Read the header.
*/
  pbma_read_header ( file_in, xsize, ysize );
/*
  Allocate storage for the data.
*/
  *b = ( int * ) malloc ( ( *xsize ) * ( *ysize ) * sizeof ( int ) );
/*
  Read the data.
*/
  pbma_read_data ( file_in, *xsize, *ysize, *b );
/*
  Close the file.
*/
  fclose ( file_in );

  return;
}
/******************************************************************************/

void pbma_read_data ( FILE *file_in, int xsize, int ysize, int *b )

/******************************************************************************/
/*
  Purpose:

    PBMA_READ_DATA reads the data in an ASCII PBM file.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    05 June 2010

  Author:

    John Burkardt

  Parameters:

    Input, FILE *FILE_IN, a pointer to the file.

    Input, int XSIZE, YSIZE, the number of rows and columns of data.

    Output, int *B, the array of XSIZE by YSIZE data values.
*/
{
  int i;
  int j;
  int n;

  for ( j = 0; j < ysize; j++ )
  {
    for ( i = 0; i < xsize; i++ )
    {
      n = fscanf ( file_in, "%d", b );

      if ( n != 1 )
      {
        fprintf ( stdout, "\n" );
        fprintf ( stdout, "PBMA_READ_DATA - Fatal error!\n" );
        fprintf ( stdout, "  Unable to read data.\n" );
        exit ( 1 );
      }

      b = b + 1;
    }
  }
  return;
}
/******************************************************************************/

void pbma_read_header ( FILE *file_in, int *xsize, int *ysize )

/******************************************************************************/
/*
  Purpose:

    PBMA_READ_HEADER reads the header of an ASCII PBM file.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    05 June 2010

  Author:

    John Burkardt

  Parameters:

    Input, FILE *FILE_IN, a pointer to the file.

    Output, int *XSIZE, *YSIZE, the number of rows and columns of data.
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
    error = fgets ( line, LINE_MAX, file_in );

    if ( !error )
    {
      fprintf ( stdout, "\n" );
      fprintf ( stdout, "PBMA_READ_HEADER - Fatal error!\n" );
      fprintf ( stdout, "  End of file.\n" );
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
      if ( strcmp ( word, "P1" ) != 0 && strcmp ( word, "p1" ) != 0 )
      {
        fprintf ( stderr, "\n" );
        fprintf ( stderr, "PBMA_READ_HEADER - Fatal error.\n" );
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
      break;
    }
  }
  return;
}
/******************************************************************************/

void pbma_read_test ( char *file_in_name )

/******************************************************************************/
/*
  Purpose:

    PBMA_READ_TEST tests the ASCII PBM read routines.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    05 June 2010

  Author:

    John Burkardt

  Parameters:

    Input, char *FILE_IN_NAME, the name of the file.
*/
{
  int *b;
  int xsize;
  int ysize;

  b = NULL;
/*
  Read the data.
*/
  pbma_read ( file_in_name, &xsize, &ysize, &b );

  fprintf ( stdout, "\n" );
  fprintf ( stdout, "  PBMA_READ was able to read \"%s\".\n", file_in_name );
/*
  Check the data.
*/
  pbma_check_data ( xsize, ysize, b );

  free ( b );

  fprintf ( stdout, "\n" );
  fprintf ( stdout, "  PBMA_CHECK_DATA approved the data from the file.\n" );

  return;
}
/******************************************************************************/

void pbma_write ( char *file_out_name, int xsize, int ysize, int *b )

/******************************************************************************/
/*
  Purpose:

    PBMA_WRITE writes the header and data for an ASCII PBM file.
 
  Example:

    P1
    # feep.pbm
    24 7
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 1 1 1 1 0 0 1 1 1 1 0 0 1 1 1 1 0 0 1 1 1 1 0
    0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 1 0
    0 1 1 1 0 0 0 1 1 1 0 0 0 1 1 1 0 0 0 1 1 1 1 0
    0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0
    0 1 0 0 0 0 0 1 1 1 1 0 0 1 1 1 1 0 0 1 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:
 
    05 June 2010
 
  Author:
 
    John Burkardt

  Parameters:

    Input, char *FILE_OUT_NAME, the name of the file.

    Input, int XSIZE, YSIZE, the number of rows and columns of data.

    Input, int *B, the array of XSIZE by YSIZE data values.
*/
{
  FILE *file_out;
  int i;
  int *indexg;
  int j;
/*
  Open the output file.
*/
  file_out = fopen ( file_out_name, "wt" );

  if ( !file_out )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "PBMA_WRITE - Fatal error!\n" );
    fprintf ( stderr, "  Cannot open the output file \"%d\".\n", file_out_name );
    exit ( 1 );
  }
/*
  Write the header.
*/
  pbma_write_header ( file_out, file_out_name, xsize, ysize );
/*
  Write the data.
*/
  pbma_write_data ( file_out, xsize, ysize, b );
/*
  Close the file.
*/
  fclose ( file_out );

  return;
}
/******************************************************************************/

void pbma_write_data ( FILE *file_out, int xsize, int ysize, int *b )

/******************************************************************************/
/*
  Purpose:

    PBMA_WRITE_DATA writes the data for an ASCII PBM file.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    05 June 2010

  Author:

    John Burkardt

  Parameters:

    Input, ofstream &FILE_OUT, a pointer to the file.

    Input, int XSIZE, YSIZE, the number of rows and columns of data.

    Input, int *B, the arrays of XSIZE by YSIZE data values.
*/
{
  int i;
  int *indexb;
  int j;
  int numval;

  indexb = b;
  numval = 0;

  for ( j = 0; j < ysize; j++ )
  {
    for ( i = 0; i < xsize; i++ )
    {
      fprintf ( file_out, " %d", *indexb );
      numval = numval + 1;
      indexb = indexb + 1;

      if ( numval % 12 == 0 || i == xsize - 1 || numval == xsize * ysize )
      {
        fprintf ( file_out, "\n" );
      }
    }
  }
  return;
}
/******************************************************************************/

void pbma_write_header ( FILE *file_out, char *file_out_name, int xsize, 
  int ysize )

/******************************************************************************/
/*
  Purpose:

    PBMA_WRITE_HEADER writes the header of an ASCII PBM file.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    05 June 2010

  Author:

    John Burkardt

  Parameters:

    Input, ofstream &FILE_OUT, a pointer to the file.

    Input, char *FILE_OUT_NAME, the name of the file.

    Input, int XSIZE, YSIZE, the number of rows and columns of data.
*/
{
  fprintf ( file_out, "%s\n", "P1" );
  fprintf ( file_out, "# %s created by PBMA_IO::PBMA_WRITE.\n", file_out_name );
  fprintf ( file_out, "%d  %d\n", xsize, ysize );

  return;
}
/******************************************************************************/

void pbma_write_test ( char *file_out_name )

/******************************************************************************/
/*
  Purpose:

    PBMA_WRITE_TEST tests the ASCII PBM routines.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    05 June 2010

  Author:

    John Burkardt

  Parameters:

    Input, char *FILE_OUT_NAME, the name of the file.
*/
{
  int *b;
  int xsize;
  int ysize;

  xsize = 200;
  ysize = 200;
/*
  Allocate memory.
*/
  b = ( int * ) malloc ( xsize * ysize * sizeof ( int ) );
/*
  Set the data.
*/
  pbma_example ( xsize, ysize, b );
/*
  Write the data to the file.
*/
  pbma_write ( file_out_name, xsize, ysize, b );

  free ( b );
 
  return;
}
