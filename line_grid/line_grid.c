# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>

# include "line_grid.h"

/******************************************************************************/

double *line_grid ( int n, double a, double b, int c )

/******************************************************************************/
/*
  Purpose:

    LINE_GRID: grid points over the interior of a line segment in 1D.

  Discussion:

    In 1D, a grid is to be created using N points.

    Over the interval [A(i),B(i)], we have 5 choices for grid centering:
      1: 0,   1/3, 2/3, 1
      2: 1/5, 2/5, 3/5, 4/5
      3: 0,   1/4, 2/4, 3/4
      4: 1/4, 2/4, 3/4, 1
      5: 1/8, 3/8, 5/8, 7/8

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    02 September 2014

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of points.
 
    Input, double A, B, the endpoints.

    Input, int C, the grid centering.
    1 <= C <= 5.

    Output, double LINE_GRID[N], the points.
*/
{
  int j;
  double *x;

  x = ( double * ) malloc ( n * sizeof ( double ) );
/*
  Create the 1D grids in each dimension.
*/
  for ( j = 0; j < n; j++ )
  {
    if ( c == 1 )
    {
      if ( n == 1 )
      {
        x[j] = 0.5 * ( a + b );
      }
      else
      {
        x[j] = (   ( double ) ( n - j - 1 ) * a   
                 + ( double ) (     j     ) * b ) 
                 / ( double ) ( n    - 1 );
      }
    }
    else if ( c == 2 )
    {
      x[j] = (   ( double ) ( n - j     ) * a   
               + ( double ) (     j + 1 ) * b ) 
               / ( double ) ( n     + 1 );
    }
    else if ( c == 3 )
    {
      x[j] = (   ( double ) ( n - j     ) * a   
               + ( double ) (     j - 2 ) * b ) 
               / ( double ) ( n         );
    }
    else if ( c == 4 )
    {
      x[j] = (   ( double ) ( n - j - 1 ) * a   
               + ( double ) (     j + 1 ) * b ) 
               / ( double ) ( n         );
    }
    else if ( c == 5 )
    {
      x[j] = (   ( double ) ( 2 * n - 2 * j - 1 ) * a   
               + ( double ) (         2 * j + 1 ) * b ) 
               / ( double ) ( 2 * n             );
    }
  }

  return x;
}
/******************************************************************************/

void r8vec_print ( int n, double a[], char *title )

/******************************************************************************/
/*
  Purpose:

    R8VEC_PRINT prints an R8VEC.

  Discussion:

    An R8VEC is a vector of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 April 2009

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of components of the vector.

    Input, double A[N], the vector to be printed.

    Input, char *TITLE, a title.
*/
{
  int i;

  fprintf ( stdout, "\n" );
  fprintf ( stdout, "%s\n", title );
  fprintf ( stdout, "\n" );
  for ( i = 0; i < n; i++ )
  {
    fprintf ( stdout, "  %8d: %14g\n", i, a[i] );
  }

  return;
}
/******************************************************************************/

void timestamp ( )

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

