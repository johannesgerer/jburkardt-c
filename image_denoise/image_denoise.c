# include <stdlib.h>
# include <stdio.h>
# include <time.h>
# include <string.h>

# include "image_denoise.h"

/******************************************************************************/

int *gray_median_news ( int m, int n, int gray[] )

/******************************************************************************/
/*
  Purpose:

    GRAY_MEDIAN_NEWS uses a median NEWS filter on a gray scale image to remove noise.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 February 2012

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns of pixels.

    Input, int GRAY[M*N], the noisy grayscale data.

    Output, int GRAY_MEDIAN_NEWS[M*N], the grayscale data for the filtered image.
*/
{
  int *gray2;
  int i;
  int j;
  int p[5];

  gray2 = ( int * ) malloc ( m * n * sizeof ( int ) );
/*
  Process the main part of the image:
*/
  for ( i = 1; i < m - 1; i++ )
  {
    for ( j = 1; j < n - 1; j++ )
    {
      p[0] = gray[i-1+ j   *m];
      p[1] = gray[i+1+ j   *m];
      p[2] = gray[i  +(j+1)*m];
      p[3] = gray[i  +(j-1)*m];
      p[4] = gray[i  + j   *m];

      gray2[i+j*m] = i4vec_median ( 5, p );
    }
  }
/*
  Process the four borders.
  Get an odd number of data points, 
*/
  for ( i = 1; i < m - 1; i++ )
  {
    j = 0;
    p[0] = gray[i-1+ j   *m];
    p[1] = gray[i+1+ j   *m];
    p[2] = gray[i  + j   *m];
    p[3] = gray[i  +(j+1)*m];
    p[4] = gray[i  +(j+2)*m];
    gray2[i+j*m] = i4vec_median ( 5, p );

    j = n - 1;
    p[0] = gray[i-1+ j   *m];
    p[1] = gray[i+1+ j   *m];
    p[2] = gray[i  +(j-2)*m];
    p[3] = gray[i  +(j-1)*m];
    p[4] = gray[i  + j   *m];
    gray2[i+j*m] = i4vec_median ( 5, p );
  }

  for ( j = 1; j < n - 1; j++ )
  {
    i = 0;
    p[0] = gray[i  + j   *m];
    p[1] = gray[i+1+ j   *m];
    p[2] = gray[i+2+ j   *m];
    p[3] = gray[i  +(j-1)*m];
    p[4] = gray[i  +(j+1)*m];
    gray2[i+j*m] = i4vec_median ( 5, p );

    i = m - 1;
    p[0] = gray[i-2+ j   *m];
    p[1] = gray[i-1+ j   *m];
    p[2] = gray[i  + j   *m];
    p[3] = gray[i  +(j-1)*m];
    p[4] = gray[i  +(j+1)*m];
    gray2[i+j*m] = i4vec_median ( 5, p );
  }
/*
  Process the four corners.
*/
  i = 0;
  j = 0;
  p[0] = gray[i+1+ j   *m];
  p[1] = gray[i  + j   *m];
  p[2] = gray[i  +(j+1)*m];
  gray2[i+j*m] = i4vec_median ( 3, p );

  i = 0;
  j = n - 1;
  p[0] = gray[i+1+ j   *m];
  p[1] = gray[i  + j   *m];
  p[2] = gray[i  +(j-1)*m];
  gray2[i+j*m] = i4vec_median ( 3, p );

  i = m - 1;
  j = 0;
  p[0] = gray[i-1+ j   *m];
  p[1] = gray[i  + j   *m];
  p[2] = gray[i  +(j+1)*m];
  gray2[i+j*m] = i4vec_median ( 3, p );

  i = m - 1;
  j = n - 1;
  p[0] = gray[i-1+ j   *m];
  p[1] = gray[i  + j   *m];
  p[2] = gray[i  +(j-1)*m];
  gray2[i+j*m] = i4vec_median ( 3, p );

  return gray2;
}
/******************************************************************************/

int i4vec_frac ( int n, int a[], int k )

/******************************************************************************/
/*
  Purpose:

    I4VEC_FRAC searches for the K-th smallest entry in an I4VEC.

  Discussion:

    An I4VEC is a vector of I4's.

    Hoare's algorithm is used.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    18 September 2005

  Parameters:

    Input, int N, the number of elements of A.

    Input/output, int A[N].
    On input, A is the array to search.
    On output, the elements of A have been somewhat rearranged.

    Input, int K, the fractile to be sought.  If K = 1, the minimum
    entry is sought.  If K = N, the maximum is sought.  Other values
    of K search for the entry which is K-th in size.  K must be at
    least 1, and no greater than N.

    Output, double I4VEC_FRAC, the value of the K-th fractile of A.
*/
{
  int frac;
  int i;
  int iryt;
  int j;
  int left;
  int temp;
  int w;
  int x;

  if ( n <= 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "I4VEC_FRAC - Fatal error!\n" );
    fprintf ( stderr, "  Illegal nonpositive value of N = %d\n", n );
    exit ( 1 );
  }

  if ( k <= 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "I4VEC_FRAC - Fatal error!\n" );
    fprintf ( stderr, "  Illegal nonpositive value of K = %d\n", k );
    exit ( 1 );
  }

  if ( n < k )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "I4VEC_FRAC - Fatal error!\n" );
    fprintf ( stderr, "  Illegal N < K, K = %d\n", k );
    exit ( 1 );
  }

  left = 1;
  iryt = n;

  for ( ; ; )
  {
    if ( iryt <= left )
    {
      frac = a[k-1];
      break;
    }

    x = a[k-1];
    i = left;
    j = iryt;

    for ( ; ; )
    {
      if ( j < i )
      {
        if ( j < k )
        {
          left = i;
        }
        if ( k < i )
        {
          iryt = j;
        }
        break;
      }
/*
  Find I so that X <= A(I).
*/
      while ( a[i-1] < x )
      {
        i = i + 1;
      }
/*
  Find J so that A(J) <= X.
*/
      while ( x < a[j-1] )
      {
        j = j - 1;
      }

      if ( i <= j )
      {
        temp   = a[i-1];
        a[i-1] = a[j-1];
        a[j-1] = temp;
        i = i + 1;
        j = j - 1;
      }
    }
  }

  return frac;
}
/******************************************************************************/

int i4vec_median ( int n, int a[] )

/******************************************************************************/
/*
  Purpose:

    I4VEC_MEDIAN returns the median of an unsorted I4VEC.

  Discussion:

    An I4VEC is a vector of I4's.

    Hoare's algorithm is used.  The values of the vector are
    rearranged by this routine.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    18 September 2005

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of elements of A.

    Input/output, int A[N], the array to search.  On output,
    the order of the elements of A has been somewhat changed.

    Output, int I4VEC_MEDIAN, the value of the median of A.
*/
{
  int k;
  int median;

  k = ( n + 1 ) / 2;

  median = i4vec_frac ( n, a, k );

  return median;
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

void pgma_write ( char *file_out_name, int xsize, int ysize, int *g )

/******************************************************************************/
/*
  Purpose:

    PGMA_WRITE writes the header and data for an ASCII PGM file.
 
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
 
    05 June 2010
 
  Author:
 
    John Burkardt

  Parameters:

    Input, char *FILE_OUT_NAME, the name of the file.

    Input, int XSIZE, YSIZE, the number of rows and columns of data.

    Input, int *G, the array of XSIZE by YSIZE data values.
*/
{
  FILE *file_out;
  int i;
  int *indexg;
  int j;
  int maxg;
/*
  Open the output file.
*/
  file_out = fopen ( file_out_name, "wt" );

  if ( !file_out )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "PGMA_WRITE - Fatal error!\n" );
    fprintf ( stderr, "  Cannot open the output file \"%s\".\n", file_out_name );
    exit ( 1 );
  }
/*
  Compute the maximum.
*/
  maxg = 0;
  indexg = g;

  for ( j = 0; j < ysize; j++ )
  {
    for ( i = 0; i < xsize; i++ )
    {
      if ( maxg < *indexg )
      {
        maxg = *indexg;
      }
      indexg = indexg + 1;

    }
  }
/*
  Write the header.
*/
  pgma_write_header ( file_out, file_out_name, xsize, ysize, maxg );
/*
  Write the data.
*/
  pgma_write_data ( file_out, xsize, ysize, g );
/*
  Close the file.
*/
  fclose ( file_out );

  return;
}
/******************************************************************************/

void pgma_write_data ( FILE *file_out, int xsize, int ysize, int *g )

/******************************************************************************/
/*
  Purpose:

    PGMA_WRITE_DATA writes the data for an ASCII PGM file.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    05 June 2010

  Author:

    John Burkardt

  Parameters:

    Input, FILE *FILE_OUT, a pointer to the file.

    Input, int XSIZE, YSIZE, the number of rows and columns of data.

    Input, int *G, the array of XSIZE by YSIZE data.
*/
{
  int i;
  int *indexg;
  int j;
  int numval;

  indexg = g;
  numval = 0;

  for ( j = 0; j < ysize; j++ )
  {
    for ( i = 0; i < xsize; i++ )
    {
      fprintf ( file_out, "%d", *indexg );
      numval = numval + 1;
      indexg = indexg + 1;

      if ( numval % 12 == 0 || i == xsize - 1 || numval == xsize * ysize )
      {
        fprintf ( file_out, "\n" );
      }
      else
      {
        fprintf ( file_out, " " );
      }

    }
  }
  return;
}
/******************************************************************************/

void pgma_write_header ( FILE *file_out, char *file_out_name, int xsize, 
  int ysize, int maxg )

/******************************************************************************/
/*
  Purpose:

    PGMA_WRITE_HEADER writes the header of an ASCII PGM file.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    05 June 2010

  Author:

    John Burkardt

  Parameters:

    Input, FILE *FILE_OUT, a pointer to the file.

    Input, char *FILE_OUT_NAME, the name of the file.

    Input, int XSIZE, YSIZE, the number of rows and columns of data.

    Input, int MAXG, the maximum gray value.
*/
{
  fprintf ( file_out, "P2\n" );
  fprintf ( file_out, "# %s created by PGMA_IO::PGMA_WRITE.\n",
    file_out_name );
  fprintf ( file_out, "%d %d\n", xsize, ysize );
  fprintf ( file_out, "%d\n", maxg );

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

  fprintf ( stdout, "%s\n", time_buffer );

  return;
# undef TIME_SIZE
}