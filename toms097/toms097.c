# include <stdlib.h>
# include <stdio.h>
# include <time.h>

# include "toms097.h"

/******************************************************************************/

int i4_huge ( )

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
  const int value = 2147483647;

  return value;
}
/******************************************************************************/

void i4mat_shortest_path ( int n, int m[] )

/******************************************************************************/
/*
  Purpose:

    I4MAT_SHORTEST_PATH computes the shortest distance between all pairs of points.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 March 2014

  Author:

    John Burkardt

  Reference:

    Robert Floyd,
    Algorithm 97, Shortest Path,
    Communications of the ACM,
    Volume 5, Number 6, June 1962, page 345.

  Parameters:

    Input, int N, the number of points.

    Input/output, int M[N*N].
    On input, M(I,J) contains the length of the direct link between 
    nodes I and J, or HUGE if there is no direct link.
    On output, M(I,J) contains the distance between nodes I and J,
    that is, the length of the shortest path between them.  If there
    is no such path, then M(I,J) will remain HUGE.
*/
{
  int i;
  const int i4_inf = 2147483647;
  int j;
  int k;
  int s;

  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      if ( m[j+i*n] < i4_inf )
      {
        for ( k = 0; k < n; k++ )
        {
          if ( m[i+k*n] < i4_inf )
          {
            s = m[j+i*n] + m[i+k*n];
            if ( s < m[j+k*n] )
            {
              m[j+k*n] = s;
            }
          }
        }
      }
    }
  }
  return;
}
/******************************************************************************/

double r8_huge ( )

/******************************************************************************/
/*
  Purpose:

    R8_HUGE returns a "huge" R8.

  Discussion:

    The value returned by this function is NOT required to be the
    maximum representable R8.  This value varies from machine to machine,
    from compiler to compiler, and may cause problems when being printed.
    We simply want a "very large" but non-infinite number.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    06 October 2007

  Author:

    John Burkardt

  Parameters:

    Output, double R8_HUGE, a "huge" R8 value.
*/
{
  double value;

  value = 1.0E+30;

  return value;
}
/******************************************************************************/

void r8mat_shortest_path ( int n, double m[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_SHORTEST_PATH computes the shortest distance between all pairs of points.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 March 2014

  Author:

    John Burkardt

  Reference:

    Robert Floyd,
    Algorithm 97, Shortest Path,
    Communications of the ACM,
    Volume 5, Number 6, June 1962, page 345.

  Parameters:

    Input, int N, the number of points.

    Input/output, double M[N*N].
    On input, M(I,J) contains the length of the direct link between 
    nodes I and J, or HUGE if there is no direct link.
    On output, M(I,J) contains the distance between nodes I and J,
    that is, the length of the shortest path between them.  If there
    is no such path, then M(I,J) will remain HUGE.
*/
{
  int i;
  int j;
  int k;
  const double r8_inf = 1.0E+30;
  double s;

  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      if ( m[j+i*n] < r8_inf )
      {
        for ( k = 0; k < n; k++ )
        {
          if ( m[i+k*n] < r8_inf )
          {
            s = m[j+i*n] + m[i+k*n];
            if ( s < m[j+k*n] )
            {
              m[j+k*n] = s;
            }
          }
        }
      }
    }
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
