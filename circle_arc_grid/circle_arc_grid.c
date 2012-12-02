# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <string.h>

# include "circle_arc_grid.h"

/******************************************************************************/

double *circle_arc_grid ( double r, double c[2], double a[2], int n )

/******************************************************************************/
/*
  Purpose:

    CIRCLE_ARC_GRID computes grid points along a circular arc.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    13 November 2011

  Author:

    John Burkardt

  Parameters:

    Input, double R, the radius of the circle.

    Input, double C[2], the coordinates of the center.

    Input, double A[2], the angle of the first and last
    points, in DEGREES.

    Input, int N, the number of points.

    Output, double CIRCLE_ARC_GRID[2*N], the grid points.
*/
{
  double aj;
  int j;
  double pi = 3.141592653589793;
  double *xy;

  xy = ( double * ) malloc ( 2 * n * sizeof ( double ) );

  for ( j = 0; j < n; j++ )
  {
    aj = ( ( double ) ( n - j - 1 ) * a[0]
         + ( double ) (     j     ) * a[1] ) 
         / ( double ) ( n     - 1 );

    xy[0+j*2] = c[0] + r * cos ( aj * pi / 180.0 );
    xy[1+j*2] = c[1] + r * sin ( aj * pi / 180.0 );
  }

  return xy;
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
