# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <string.h>

# include "ellipse_grid.h"

/******************************************************************************/

double *ellipse_grid ( int n, double r[2], double c[2], int ng )

/******************************************************************************/
/*
  Purpose:

    ELLIPSE_GRID generates grid points inside an ellipse.

  Discussion:

    The ellipse is specified as

      ( ( X - C1 ) / R1 )^2 + ( ( Y - C2 ) / R2 )^2 = 1

    The user supplies a number N.  There will be N+1 grid points along
    the shorter axis.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 November 2011

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of subintervals.

    Input, double R[2], the half axis lengths.

    Input, double C[2], the center of the ellipse.

    Input, int NG, the number of grid points inside the ellipse.

    Output, double ELLIPSE_GRID[2*NG], the grid points.
*/
{
  double h;
  int i;
  int j;
  int ni;
  int nj;
  int p;
  double x;
  double *xy;
  double y;

  xy = ( double * ) malloc ( 2 * ng * sizeof ( double ) );

  if ( r[0] < r[1] )
  {
    h = 2.0 * r[0] / ( double ) ( 2 * n + 1 );
    ni = n;
    nj = i4_ceiling ( r[1] / r[0] ) * ( double ) ( n );
  }
  else
  {
    h = 2.0 * r[1] / ( double ) ( 2 * n + 1 );
    nj = n;
    ni = i4_ceiling ( r[0] / r[1] ) * ( double ) ( n );
  }

  p = 0;

  for ( j = 0; j <= nj; j++ )
  {
    i = 0;
    x = c[0];
    y = c[1] + ( double ) ( j ) * h;

    xy[0+p*2] = x;
    xy[1+p*2] = y;
    p = p + 1;

    if ( 0 < j )
    {
      xy[0+p*2] = x;
      xy[1+p*2] = 2.0 * c[1] - y;
      p = p + 1;
    }

    for ( ; ; )
    {
      i = i + 1;
      x = c[0] + ( double ) ( i ) * h;

      if ( 1.0 < pow ( ( x - c[0] ) / r[0], 2 ) 
               + pow ( ( y - c[1] ) / r[1], 2 ) )
      {
        break;
      }

      xy[0+p*2] = x;
      xy[1+p*2] = y;
      p = p + 1;
      xy[0+p*2] = 2.0 * c[0] - x;
      xy[1+p*2] = y;
      p = p + 1;

      if ( 0 < j )
      {
        xy[0+p*2] = x;
        xy[1+p*2] = 2.0 * c[1] - y;
        p = p + 1;
        xy[0+p*2] = 2.0 * c[0] - x;
        xy[1+p*2] = 2.0 * c[1] - y;
        p = p + 1;
      }
    }
  }
  return xy;
}
/******************************************************************************/

int ellipse_grid_count ( int n, double r[2], double c[2] )

/******************************************************************************/
/*
  Purpose:

    ELLIPSE_GRID_COUNT counts the grid points inside an ellipse.

  Discussion:

    The ellipse is specified as

      ( ( X - C1 ) / R1 )^2 + ( ( Y - C2 ) / R2 )^2 = 1

    The user supplies a number N.  There will be N+1 grid points along
    the shorter axis.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 November 2011

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of subintervals.

    Input, double R[2], the half axis lengths.

    Input, double C[2], the center of the ellipse.

    Output, int ELLIPSE_GRID)_COUNT, the number of grid points inside 
    the ellipse.
*/
{
  double h;
  int i;
  int j;
  int ni;
  int nj;
  int p;
  double x;
  double y;

  if ( r[0] < r[1] )
  {
    h = 2.0 * r[0] / ( double ) ( 2 * n + 1 );
    ni = n;
    nj = i4_ceiling ( r[1] / r[0] ) * ( double ) ( n );
  }
  else
  {
    h = 2.0 * r[1] / ( double ) ( 2 * n + 1 );
    nj = n;
    ni = i4_ceiling ( r[0] / r[1] ) * ( double ) ( n );
  }

  p = 0;

  for ( j = 0; j <= nj; j++ )
  {
    i = 0;
    x = c[0];
    y = c[1] + ( double ) ( j ) * h;

    p = p + 1;

    if ( 0 < j )
    {
      p = p + 1;
    }

    for ( ; ; )
    {
      i = i + 1;
      x = c[0] + ( double ) ( i ) * h;

      if ( 1.0 < pow ( ( x - c[0] ) / r[0], 2 ) 
               + pow ( ( y - c[1] ) / r[1], 2 ) )
      {
        break;
      }

      p = p + 1;
      p = p + 1;

      if ( 0 < j )
      {
        p = p + 1;
        p = p + 1;
      }
    }
  }

  return p;
}
/******************************************************************************/

int i4_ceiling ( double x )

/******************************************************************************/
/*
  Purpose:

    I4_CEILING rounds an R8 up to the nearest I4.

  Discussion:

    The "ceiling" of X is the value of X rounded towards plus infinity.

  Example:

    X        I4_CEILING(X)

   -1.1      -1
   -1.0      -1
   -0.9       0
   -0.1       0
    0.0       0
    0.1       1
    0.9       1
    1.0       1
    1.1       2
    2.9       3
    3.0       3
    3.14159   4

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    10 November 2011

  Author:

    John Burkardt

  Parameters:

    Input, double X, the number whose ceiling is desired.

    Output, int I4_CEILING, the ceiling of X.
*/
{
  int value;

  value = ( int ) x;

  if ( value < x )
  {
    value = value + 1;
  }

  return value;
}
/******************************************************************************/

void r82vec_print_part ( int n, double a[], int max_print, char *title )

/******************************************************************************/
/*
  Purpose:

    R82VEC_PRINT_PART prints "part" of an R82VEC.

  Discussion:

    The user specifies MAX_PRINT, the maximum number of lines to print.

    If N, the size of the vector, is no more than MAX_PRINT, then
    the entire vector is printed, one entry per line.

    Otherwise, if possible, the first MAX_PRINT-2 entries are printed,
    followed by a line of periods suggesting an omission,
    and the last entry.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 November 2011

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries of the vector.

    Input, double A[2*N], the vector to be printed.

    Input, int MAX_PRINT, the maximum number of lines
    to print.

    Input, char *TITLE, a title.
*/
{
  int i;

  if ( max_print <= 0 )
  {
    return;
  }

  if ( n <= 0 )
  {
    return;
  }

  fprintf ( stdout, "\n" );
  fprintf ( stdout, "%s\n", title );
  fprintf ( stdout, "\n" );

  if ( n <= max_print )
  {
    for ( i = 0; i < n; i++ )
    {
      fprintf ( stdout, "  %8d: %14g  %14g\n", i, a[0+i*2], a[1+i*2] );
    }
  }
  else if ( 3 <= max_print )
  {
    for ( i = 0; i < max_print - 2; i++ )
    {
      fprintf ( stdout, "  %8d: %14g  %14g\n", i, a[0+i*2], a[1+i*2] );
    }
    fprintf ( stdout, "  ........  ..............  ..............\n" );
    i = n - 1;
    fprintf ( stdout, "  %8d: %14g  %14g\n", i, a[0+i*2], a[1+i*2] );
  }
  else
  {
    for ( i = 0; i < max_print - 1; i++ )
    {
      fprintf ( stdout, "  %8d: %14g  %14g\n", i, a[0+i*2], a[1+i*2] );
    }
    i = max_print - 1;
    fprintf ( stdout, "  %8d: %14g  %14g  ...more entries...\n", 
      i, a[0+i*2], a[1+i*2] );
  }

  return;
}
/******************************************************************************/

void r8mat_write ( char *output_filename, int m, int n, double table[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_WRITE writes an R8MAT file.

  Discussion:

    An R8MAT is an array of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    01 June 2009

  Author:

    John Burkardt

  Parameters:

    Input, char *OUTPUT_FILENAME, the output filename.

    Input, int M, the spatial dimension.

    Input, int N, the number of points.

    Input, double TABLE[M*N], the data.
*/
{
  int i;
  int j;
  FILE *output;
/*
  Open the file.
*/
  output = fopen ( output_filename, "wt" );

  if ( !output )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8MAT_WRITE - Fatal error!\n" );
    fprintf ( stderr, "  Could not open the output file.\n" );
    exit ( 1 );
  }
/*
  Write the data.
*/
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      fprintf ( output, "  %24.16g", table[i+j*m] );
    }
    fprintf ( output, "\n" );
  }
/*
  Close the file.
*/
  fclose ( output );

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
