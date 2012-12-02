# include <stdlib.h>
# include <stdio.h>
# include <time.h>

# include "triangle_grid.h"

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

  printf ( "%s\n", time_buffer );

  return;
# undef TIME_SIZE
}
/******************************************************************************/

double *triangle_grid ( int n, double t[] )

/******************************************************************************/
/*
  Purpose:

    TRIANGLE_GRID computes points on a triangular grid.

  Discussion:

    The grid is defined by specifying the coordinates of an enclosing
    triangle T, and the number of subintervals each side of the triangle
    should be divided into.

    Choosing N = 10, for instance, breaks each side into 10 subintervals,
    and produces a grid of ((10+1)*(10+2))/2 = 66 points.

              X
             9 X
            8 9 X
           7 8 9 X
          6 7 8 9 X
         5 6 7 8 9 X
        4 5 6 7 8 9 X
       3 4 5 6 7 8 9 X
      2 3 4 5 6 7 8 9 X
     1 2 3 4 5 6 7 8 9 X
    0 1 2 3 4 5 6 7 8 9 X

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 September 2010

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of subintervals.

    Input, double T[2*3], the coordinates of the points
    defining the triangle.

    Output, double *TRIANGLE_GRID[2*((N+1)*(N+2))/2], the coordinates
    of the points in the triangle.
*/
{
  int i;
  double ir;
  int j;
  double jr;
  int k;
  double kr;
  int l;
  int ng;
  double nr;
  int p;
  double *tg;

  ng = ( ( n + 1 ) * ( n + 2 ) ) / 2;
  tg = ( double * ) malloc ( 2 * ng * sizeof ( double ) );

  p = 0;
  nr = ( double ) ( n );

  for ( i = 0; i <= n; i++ )
  {
    ir = ( double ) ( i );
    for ( j = 0; j <= n - i; j++ )
    {
      jr = ( double ) ( j );
      k = n - i - j;
      kr = ( double ) ( k );
      for ( l = 0; l < 2; l++ )
      {
        tg[l+p*2] = ( ir * t[l+0*2] + jr * t[l+1*2] + kr * t[l+2*2] ) / nr;
      }
      p = p + 1;
    }
  }

  return tg;
}
