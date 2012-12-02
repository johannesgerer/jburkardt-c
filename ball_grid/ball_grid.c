# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>

# include "ball_grid.h"

/******************************************************************************/

double *ball_grid ( int n, double r, double c[3], int ng )

/******************************************************************************/
/*
  Purpose:

    BALL_GRID computes grid points inside a ball.

  Discussion:

    The grid is defined by specifying the radius and center of the ball,
    and the number of subintervals N into which the horizontal radius
    should be divided.  Thus, a value of N = 2 will result in 5 points
    along that horizontal line.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 November 2011

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of subintervals.

    Input, double R, the radius of the ball.

    Input, double C[3], the coordinates of the center of the ball.

    Input, int NG, the number of grid points, as determined by
    BALL_GRID_COUNT.

    Output, double BALL_GRID[3*NG], the grid points inside the ball.
*/
{
  double *bg;
  int i;
  int j;
  int k;
  int p;
  double x;
  double y;
  double z;

  bg = ( double * ) malloc ( 3 * ng * sizeof ( double ) );

  p = 0;

  for ( i = 0; i <= n; i++ )
  {
    x = c[0] + r * ( double ) ( 2 * i ) / ( double ) ( 2 * n + 1 );
    for ( j = 0; j <= n; j++ )
    {    
      y = c[1] + r * ( double ) ( 2 * j ) / ( double ) ( 2 * n + 1 );
      for ( k = 0; k <= n; k++ )
      {
        z = c[2] + r * ( double ) ( 2 * k ) / ( double ) ( 2 * n + 1 );

        if ( r * r < pow ( x - c[0], 2 ) 
                   + pow ( y - c[1], 2 )
                   + pow ( z - c[2], 2 ) )
        {
          break;
        }


        bg[0+p*3] = x;
        bg[1+p*3] = y;
        bg[2+p*3] = z;
        p = p + 1;

        if ( 0 < i )
        {
          bg[0+p*3] = 2.0 * c[0] - x;
          bg[1+p*3] = y;
          bg[2+p*3] = z;
          p = p + 1;
        }

        if ( 0 < j )
        {
          bg[0+p*3] = x;
          bg[1+p*3] = 2.0 * c[1] - y;
          bg[2+p*3] = z;
          p = p + 1;
        }

        if ( 0 < k )
        {
          bg[0+p*3] = x;
          bg[1+p*3] = y;
          bg[2+p*3] = 2.0 * c[2] - z;
          p = p + 1;
        }

        if ( 0 < i && 0 < j )
        {
          bg[0+p*3] = 2.0 * c[0] - x;
          bg[1+p*3] = 2.0 * c[1] - y;
          bg[2+p*3] = z;
          p = p + 1;
        }

        if ( 0 < i && 0 < k )
        {
          bg[0+p*3] = 2.0 * c[0] - x;
          bg[1+p*3] = y;
          bg[2+p*3] = 2.0 * c[2] - z;
          p = p + 1;
        }

        if ( 0 < j && 0 < k )
        {
          bg[0+p*3] = x;
          bg[1+p*3] = 2.0 * c[1] - y;
          bg[2+p*3] = 2.0 * c[2] - z;
          p = p + 1;
        }

        if ( 0 < i && 0 < j && 0 < k )
        {
          bg[0+p*3] = 2.0 * c[0] - x;
          bg[1+p*3] = 2.0 * c[1] - y;
          bg[2+p*3] = 2.0 * c[2] - z;
          p = p + 1;
        }
      }
    }
  }

  return bg;
}
/******************************************************************************/

int ball_grid_count ( int n, double r, double c[3] )

/******************************************************************************/
/*
  Purpose:

    BALL_GRID computes grid points inside a ball.

  Discussion:

    The grid is defined by specifying the radius and center of the ball,
    and the number of subintervals N into which the horizontal radius
    should be divided.  Thus, a value of N = 2 will result in 5 points
    along that horizontal line.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 November 2011

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of subintervals.

    Input, double R, the radius of the ball.

    Input, double C[3], the coordinates of the center of the ball.

    Output, int BALL_GRID_COUNT, the number of grid points inside the ball.
*/
{
  int i;
  int j;
  int k;
  int ng;
  double x;
  double y;
  double z;

  ng = 0;

  for ( i = 0; i <= n; i++ )
  {
    x = c[0] + r * ( double ) ( 2 * i ) / ( double ) ( 2 * n + 1 );
    for ( j = 0; j <= n; j++ )
    {    
      y = c[1] + r * ( double ) ( 2 * j ) / ( double ) ( 2 * n + 1 );
      for ( k = 0; k <= n; k++ )
      {
        z = c[2] + r * ( double ) ( 2 * k ) / ( double ) ( 2 * n + 1 );

        if ( r * r < pow ( x - c[0], 2 ) 
                   + pow ( y - c[1], 2 )
                   + pow ( z - c[2], 2 ) )
        {
          break;
        }

        ng = ng + 1;

        if ( 0 < i )
        {
          ng = ng + 1;
        }

        if ( 0 < j )
        {
          ng = ng + 1;
        }

        if ( 0 < k )
        {
          ng = ng + 1;
        }

        if ( 0 < i && 0 < j )
        {
          ng = ng + 1;
        }

        if ( 0 < i && 0 < k )
        {
          ng = ng + 1;
        }

        if ( 0 < j && 0 < k )
        {
          ng = ng + 1;
        }

        if ( 0 < i && 0 < j && 0 < k )
        {
          ng = ng + 1;
        }

      }
    }
  }

  return ng;
}
/******************************************************************************/

void r83vec_print_part ( int n, double a[], int max_print, char *title )

/******************************************************************************/
/*
  Purpose:

    R83VEC_PRINT_PART prints "part" of an R83VEC.

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

    11 November 2011

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries of the vector.

    Input, double A[3*N], the vector to be printed.

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
      fprintf ( stdout, "  %8d: %14g  %14g  %14g\n", i, a[0+i*3], a[1+i*3], a[2+i*3] );
    }
  }
  else if ( 3 <= max_print )
  {
    for ( i = 0; i < max_print - 2; i++ )
    {
      fprintf ( stdout, "  %8d: %14g  %14g  %14g\n", i, a[0+i*3], a[1+i*3], a[2+i*3] );
    }
    fprintf ( stdout, "  ........  ..............  ..............  ..............\n" );
    i = n - 1;
    fprintf ( stdout, "  %8d: %14g  %14g  %14g\n", i, a[0+i*3], a[1+i*3], a[2+i*3] );
  }
  else
  {
    for ( i = 0; i < max_print - 1; i++ )
    {
      fprintf ( stdout, "  %8d: %14g  %14g  %14g\n", i, a[0+i*3], a[1+i*3], a[2+i*3] );
    }
    i = max_print - 1;
    fprintf ( stdout, "  %8d: %14g  %14g  %14g  ...more entries...\n", 
      i, a[0+i*3], a[1+i*3], a[2+i*3] );
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
    fprintf ( stderr, "R8MAT_WRITE - Fatal error\n" );
    fprintf ( stderr, "  Could not open the output file \"%s\".\n", output_filename );
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
