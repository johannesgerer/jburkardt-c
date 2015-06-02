# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <string.h>

# include "naca.h"

/******************************************************************************/

void naca4_cambered ( double m, double p, double t, double c, int n, 
  double xc[], double xu[], double yu[], double xl[], double yl[] )

/******************************************************************************/
/*
  Purpose:

    NACA4_CAMBERED: (xu,yu), (xl,yl) for a NACA cambered 4-digit airfoil.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    21 May 2014

  Author:

    John Burkardt

  Reference:

    Eastman Jacobs, Kenneth Ward, Robert Pinkerton,
    "The characteristics of 78 related airfoil sections from tests in
    the variable-density wind tunnel",
    NACA Report 460, 1933.

  Parameters:

    Input, double M, the maximum camber.
    0.0 < M.

    Input, double P, the location of maximum camber.
    0.0 < P < 1.0

    Input, double T, the maximum relative thickness.
    0.0 < T <= 1.0

    Input, double C, the chord length.
    0.0 < C.

    Input, int N, the number of sample points.

    Input, double XC[N], points along the chord length.  
    0.0 <= XC(*) <= C.

    Output, double XU[N], YU[N], XL[N], YL[N], for each value of 
    XC, measured along the camber line, the corresponding values (XU,YU) 
    on the upper airfoil surface and (XL,YL) on the lower airfoil surface.
*/
{
  double divisor;
  double dycdx;
  int i;
  double theta;
  double yc;
  double yt;

  for ( i = 0; i < n; i++ )
  {
    if ( 0.0 <= xc[i] / c && xc[i] / c <= p )
    {
      divisor = p * p;
    }
    else if ( p <= xc[i] / c && xc[i] / c <= 1.0 )
    {
      divisor = pow ( 1.0 - p, 2 );
    }
    else
    {
      divisor = 1.0;
    }

    dycdx = 2.0 * m * ( p - xc[i] / c ) / divisor;

    theta = atan ( dycdx );
   
    yt = 5.0 * t * c * ( 
       0.2969 * sqrt ( xc[i] / c ) 
       + (((( 
         - 0.1015 ) * ( xc[i] / c ) 
         + 0.2843 ) * ( xc[i] / c ) 
         - 0.3516 ) * ( xc[i] / c ) 
         - 0.1260 ) * ( xc[i] / c ) );

    if ( 0.0 <= xc[i] / c && xc[i] / c <= p )
    {
      yc = m * xc[i] * ( 2.0 * p - xc[i] / c ) / p / p;
    }
    else if ( p <= xc[i] / c && xc[i] / c <= 1.0 )
    {
      yc = m * ( xc[i] - c ) * ( 2.0 * p - xc[i] / c - 1.0 )
        / ( 1.0 - p ) / ( 1.0 - p );
    }
    else
    {
      yc = 0.0;
    }

    xu[i] = xc[i] - yt * sin ( theta );
    yu[i] = yc + yt * cos ( theta );
    xl[i] = xc[i] + yt * sin ( theta );
    yl[i] = yc - yt * cos ( theta );
  }
  return;
}
/******************************************************************************/

double *naca4_symmetric ( double t, double c, int n, double x[] )

/******************************************************************************/
/*
  Purpose:

    NACA4_SYMMETRIC evaluates y(x) for a NACA symmetric 4-digit airfoil.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    22 May 2014

  Author:

    John Burkardt

  Reference:

    Eastman Jacobs, Kenneth Ward, Robert Pinkerton,
    "The characteristics of 78 related airfoil sections from tests in
    the variable-density wind tunnel",
    NACA Report 460, 1933.

  Parameters:

    Input, double T, the maximum relative thickness.

    Input, double C, the chord length.

    Input, int N, the number of sample points.

    Input, double X[N], points along the chord length.  
    0.0 <= X(*) <= C.

    Output, double NACA4_SYMMETRIC[N], for each value of X, the corresponding
    value of Y so that (X,Y) is on the upper wing surface, and (X,-Y) is on the
    lower wing surface.
*/
{
  int i;
  double *y;

  y = ( double * ) malloc ( n * sizeof ( double ) );

  for ( i = 0; i < n; i++ )
  {
    y[i] = 5.0 * t * c * ( 
      0.2969 * sqrt ( x[i] / c ) 
      + (((( 
      - 0.1015 ) * ( x[i] / c ) 
      + 0.2843 ) * ( x[i] / c ) 
      - 0.3516 ) * ( x[i] / c ) 
      - 0.1260 ) * ( x[i] / c ) );
  }

  return y;
}
/******************************************************************************/

double r8_max ( double x, double y )

/******************************************************************************/
/*
  Purpose:

    R8_MAX returns the maximum of two R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 May 2006

  Author:

    John Burkardt

  Parameters:

    Input, double X, Y, the quantities to compare.

    Output, double R8_MAX, the maximum of X and Y.
*/
{
  double value;

  if ( y < x )
  {
    value = x;
  }
  else
  {
    value = y;
  }
  return value;
}
/******************************************************************************/

double r8_min ( double x, double y )

/******************************************************************************/
/*
  Purpose:

    R8_MIN returns the minimum of two R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 May 2006

  Author:

    John Burkardt

  Parameters:

    Input, double X, Y, the quantities to compare.

    Output, double R8_MIN, the minimum of X and Y.
*/
{
  double value;

  if ( y < x )
  {
    value = y;
  }
  else
  {
    value = x;
  }
  return value;
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

double *r8vec_linspace_new ( int n, double a, double b )

/******************************************************************************/
/*
  Purpose:

    R8VEC_LINSPACE_NEW creates a vector of linearly spaced values.

  Discussion:

    An R8VEC is a vector of R8's.

    4 points evenly spaced between 0 and 12 will yield 0, 4, 8, 12.
 
    In other words, the interval is divided into N-1 even subintervals,
    and the endpoints of intervals are used as the points.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    29 March 2011

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the vector.

    Input, double A, B, the first and last entries.

    Output, double R8VEC_LINSPACE_NEW[N], a vector of linearly spaced data.
*/
{
  int i;
  double *x;

  x = ( double * ) malloc ( n * sizeof ( double ) );

  if ( n == 1 )
  {
    x[0] = ( a + b ) / 2.0;
  }
  else
  {
    for ( i = 0; i < n; i++ )
    {
      x[i] = ( ( double ) ( n - 1 - i ) * a 
             + ( double ) (         i ) * b ) 
             / ( double ) ( n - 1     );
    }
  }
  return x;
}
/******************************************************************************/

double r8vec_max ( int n, double r8vec[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_MAX returns the value of the maximum element in a R8VEC.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 May 2006

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the array.

    Input, double R8VEC[N], a pointer to the first entry of the array.

    Output, double R8VEC_MAX, the value of the maximum element.  This
    is set to 0.0 if N <= 0.
*/
{
  int i;
  double value;

  if ( n <= 0 )
  {
    value = 0.0;
    return value;
  }

  value = r8vec[0];

  for ( i = 1; i < n; i++ )
  {
    if ( value < r8vec[i] )
    {
      value = r8vec[i];
    }
  }
  return value;
}
/******************************************************************************/

double r8vec_min ( int n, double r8vec[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_MIN returns the value of the minimum element in a R8VEC.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 May 2006

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the array.

    Input, double R8VEC[N], the array to be checked.

    Output, double R8VEC_MIN, the value of the minimum element.
*/
{
  int i;
  double value;

  value = r8vec[0];

  for ( i = 1; i < n; i++ )
  {
    if ( r8vec[i] < value )
    {
      value = r8vec[i];
    }
  }
  return value;
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

