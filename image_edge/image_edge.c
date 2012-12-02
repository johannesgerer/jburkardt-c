# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>

# include "image_edge.h"

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
  static int value = 2147483647;

  return value;
}
/******************************************************************************/

int *i4mat_histogram ( int m, int n, int a[], int histo_num )

/******************************************************************************/
/*
  Purpose:

    I4MAT_HISTOGRAM computes a histogram of the elements of an I4MAT.

  Discussion:

    An I4MAT is an array of I4's.

    It is assumed that the entries in the vector A are nonnegative.
    Only values between 0 and HISTO_NUM will be histogrammed.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    04 June 2010

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the order of A.

    Input, int A[M*N], the array to examine.

    Input, int HISTO_NUM, the maximum value for which a
    histogram entry will be computed.

    Output, int I4MAT_HISTOGRAM[HISTO_NUM+1], contains the number of
    entries of A with the values of 0 through HISTO_NUM.
*/
{
  int *histo_gram;
  int i;
  int j;

  histo_gram = ( int * ) malloc ( ( histo_num + 1 ) * sizeof ( int ) );

  for ( i = 0; i <= histo_num; i++ )
  {
    histo_gram[i] = 0;
  }

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      if ( 0 <= a[i+j*m] && a[i+j*m] <= histo_num )
      {
        histo_gram[a[i+j*m]] = histo_gram[a[i+j*m]] + 1;
      }
    }
  }

  return histo_gram;
}
/******************************************************************************/

int i4mat_max ( int m, int n, int a[] )

/******************************************************************************/
/*
  Purpose:

    I4MAT_MAX returns the maximum of an I4MAT.

  Discussion:

    An I4MAT is an MxN array of I4's, stored by (I,J) -> [I+J*M].

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 June 2010

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows in A.

    Input, int N, the number of columns in A.

    Input, int A[M*N], the M by N matrix.

    Output, int I4MAT_MAX, the maximum entry of A.
*/
{
  int i;
  int j;
  int value;

  value = - i4_huge ( );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      if ( value < a[i+j*m] )
      {
        value = a[i+j*m];
      }
    }
  }
  return value;
}
/******************************************************************************/

int *news ( int m, int n, int a[] )

/******************************************************************************/
/*
  Purpose:

    NEWS demonstrates the NEWS stencil for edge detection.

  Discussion:

    Given a black and white image A, which we regard as an M by N array
    of pixels, we want to produce an array E of the same shape, which
    contains information describing the location of edges.

    A simple algorithm for trying to detect edges in an array that
    represents an image is the NEWS scheme.  For each pixel A(C),
    we consider its North, East, West, and South pixel neighbors.  The
    indexing of arrays and images do not correspond, so we will use
    these directions instead:

             A(N)
              |
              |
      A(W)---A(C)---A(E)
              |
              |
             A(S)

    Entry E(C) of the edge array will be computed by

      E(C) = abs ( A(N) - A(S) ) + abs ( A(E) - A(W) )

    Pixels of A that represent edges will tend to have high values
    of E, while pixels that are interior to a region of roughly the
    same shade will tend to have low values.

    Thus, an edge detection scheme would use the NEWS stencil to
    compute the E array, determine E_MAX, the maximum entry in E,
    choose some threshold value E_THRESH, and declare pixel A(I,J)
    to be associated with an edge whenever E(I,J) is greater than E_THRESH.

    In this program, we demonstrate the NEWS stencil using a PGM
    grayscale image of coins.  At the end, we use the edge information
    to produce a color image in which the edges of the coins have been
    outlined in red.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 June 2010

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns in the image.

    Input, int A[M*N], the gray scale image data, presumably
    integers between 0 and 255.

    Output, int NEWS[M*N], is 1 for each pixel that is part of an
    edge, and 0 otherwise.
*/
{
  int *b;
  int *e;
  int e_max;
  int i;
  int j;
  int thresh;
/*
  For neatness, we add a border of zeros to the image,
  then fill in the border by copying the nearby original values.
  This will be our M+2 by N+2 data array B.
*/
  b = ( int * ) malloc ( ( m + 2 ) * ( n + 2 ) * sizeof ( int ) );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      b[i+1+(j+1)*(m+2)] = a[i+j*m];
    }
  }

  for ( j = 1; j < n + 1; j++ )
  {
    b[  0+j*(m+2)] = b[1+j*(m+2)];
    b[m+1+j*(m+2)] = b[m+j*(m+2)];
  }

  for ( i = 1; i < m + 1; i++ )
  {
    b[i+   0 *(m+2)] = b[i+1*(m+2)];
    b[i+(n+1)*(m+2)] = b[i+n*(m+2)];
  }

  b[  0+   0 *(m+2)] = ( b[  0+   1 *(m+2)] + b[  1+   0 *(m+2)] ) / 2;
  b[m+1+   0 *(m+2)] = ( b[m+1+   1 *(m+2)] + b[m+0+   0 *(m+2)] ) / 2;
  b[  0+(n+1)*(m+2)] = ( b[  0+(n+0)*(m+2)] + b[  1+(n+1)*(m+2)] ) / 2;
  b[m+1+(n+1)*(m+2)] = ( b[m+1+(n+0)*(m+2)] + b[m+0+(n+1)*(m+2)] ) / 2;
/*
  Apply the NEWS Operator.  We do not process the boundary pixels.

  The picture is:

   |  0 +1  0 |     |  0  0   0 |
   |  0  0  0 |  +  | -1  0  +1 |
   |  0 -1  0 |     |  0  0   0 |
*/
  e = ( int * ) malloc ( m * n * sizeof ( int ) );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      e[i+j*m] = abs ( - b[i  +(j+1)*(m+2)] + b[i+2+(j+1)*(m+2)] )
               + abs ( - b[i+1+ j   *(m+2)] + b[i+1+(j+2)*(m+2)] );
    }
  }
  free ( b );
/*
  Remap E so the largest value is 255.
*/
  e_max = i4mat_max ( m, n, e );
/*
  Threshold the data.  Set the threshold to give enough detail
  to guess the coin denominations.
*/
  thresh = e_max / 5;

  fprintf ( stdout, "\n" );
  fprintf ( stdout, "NEWS:\n" );
  fprintf ( stdout, "  E_MAX = %d\n", e_max );
  fprintf ( stdout, "  Using threshold value THRESH = %d\n", thresh );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      if ( e[i+j*m] < thresh )
      {
        e[i+j*m] = 0;
      }
      else
      {
        e[i+j*m] = 1;
      }
    }
  }
  return e;
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

void pgma_read_data ( FILE *file_in, int xsize, int ysize, int *g )

/******************************************************************************/
/*
  Purpose:

    PGMA_READ_DATA reads the data in an ASCII PGM file.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    04 June 2010

  Author:

    John Burkardt

  Parameters:

    Input, FILE *FILE_IN, a pointer to the file containing the data.

    Input, int XSIZE, YSIZE, the number of rows and columns of data.

    Output, int *G, the array of XSIZE by YSIZE data values.
*/
{
  int i;
  int j;
  int n;

  for ( j = 0; j < ysize; j++ )
  {
    for ( i = 0; i < xsize; i++ )
    {
      n = fscanf ( file_in, "%d", g );

      if ( n != 1 )
      {
        fprintf ( stderr, "\n" );
        fprintf ( stderr, "PGMA_READ_DATA - Fatal error!\n" );
        fprintf ( stderr, "  Unable to read data.\n" );
        exit ( 1 );
      }
      g = g + 1;
    }
  }

  return;
}
/******************************************************************************/

void pgma_read_header ( FILE *file_in, int *xsize, int *ysize, int *maxg )

/******************************************************************************/
/*
  Purpose:

    PGMA_READ_HEADER reads the header of an ASCII PGM file.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 June 2010

  Author:

    John Burkardt

  Parameters:

    Input, FILE *FILE_IN, a pointer to the file.

    Output, int *XSIZE, *YSIZE, the number of rows and columns of data.

    Output, int *MAXG, the maximum gray value.
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
      if ( strcmp ( word, "P2" ) != 0 && strcmp ( word, "p2" ) != 0 )
      {
        fprintf ( stderr, "\n" );
        fprintf ( stderr, "PGMA_READ_HEADER - Fatal error.\n" );
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
      count = sscanf ( next, "%d%n", maxg, &width );
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

  fprintf ( stdout, "%s\n", time_buffer );

  return;
# undef TIME_SIZE
}
