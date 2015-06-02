# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <string.h>

# include "triangle_svg.h"

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

int r8_to_i4 ( double xmin, double xmax, double x, int ixmin, int ixmax )

/******************************************************************************/
/*
  Purpose:

    R8_TO_I4 maps real X in [XMIN, XMAX] to integer IX in [IXMIN, IXMAX].

  Discussion:

    IX := IXMIN + ( IXMAX - IXMIN ) * ( X - XMIN ) / ( XMAX - XMIN )
    IX := min ( IX, max ( IXMIN, IXMAX ) )
    IX := max ( IX, min ( IXMIN, IXMAX ) )

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 April 2014

  Author:

    John Burkardt

  Parameters:

    Input, double XMIN, XMAX, the real range.  XMAX and XMIN must not be
    equal.  It is not necessary that XMIN be less than XMAX.

    Input, double X, the real number to be converted.

    Input, int IXMIN, IXMAX, the allowed range of the output
    variable.  IXMAX corresponds to XMAX, and IXMIN to XMIN.
    It is not necessary that IXMIN be less than IXMAX.

    Output, int R8_TO_I4, the value in the range [IXMIN,IXMAX] that
    corresponds to X.
*/
{
  int ix;
  double temp;

  if ( xmax == xmin )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8_TO_I4 - Fatal error!\n" );
    fprintf ( stderr, "  XMAX = XMIN, making a zero divisor.\n" );
    fprintf ( stderr, "  XMAX = %g\n", xmax );
    fprintf ( stderr, "  XMIN = %g\n", xmin );
    exit ( 1 );
  }

  temp =
      ( ( xmax - x        ) * ( double ) ixmin
      + (        x - xmin ) * ( double ) ixmax )
      / ( xmax     - xmin );

  if ( 0.0 <= temp )
  {
    temp = temp + 0.5;
  }
  else
  {
    temp = temp - 0.5;
  }

  ix = ( int ) temp;

  return ix;
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
/******************************************************************************/

void triangle_svg ( char *plot_filename, double t[], int p_num, double p[] )

/******************************************************************************/
/*
  Purpose:

    TRIANGLE_SVG plots a triangle and points in SVG format.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    21 April 2014

  Author:

    John Burkardt

  Parameters:

    Input, char *PLOT_FILENAME, the name of the output file.

    Input, double T[2*3], points forming a triangle.

    Input, int P_NUM, the number of points.

    Input, double P[2*P_NUM], the points.
*/
{
  int i;
  int i4;
  int i4_max;
  int i4_min;
  int ii;
  int j;
  int j4;
  int j4_max;
  int j4_min;
  int node;
  FILE *output;
  int r;
  char string;
  double x;
  double x_max;
  double x_min;
  double x_scale;
  double y;
  double y_max;
  double y_min;
  double y_scale;
/*
  Determine SCALE, the maximum data range.
*/
  x_max = p[0+0*2];
  x_min = p[0+0*2];
  for ( j = 0; j < p_num; j++ )
  {
    x_max = r8_max ( x_max, p[0+j*2] );
    x_min = r8_min ( x_min, p[0+j*2] );
  }
  for ( j = 0; j < 3; j++ )
  {
    x_max = r8_max ( x_max, t[0+j*2] );
    x_min = r8_min ( x_min, t[0+j*2] );
  }
  x_scale = x_max - x_min;
  x_max = x_max + 0.05 * x_scale;
  x_min = x_min - 0.05 * x_scale;
  x_scale = x_max - x_min;

  y_max = p[1+0*2];
  y_min = p[1+0*2];
  for ( j = 0; j < p_num; j++ )
  {
    y_max = r8_max ( y_max, p[1+j*2] );
    y_min = r8_min ( y_min, p[1+j*2] );
  }
  for ( j = 0; j < 3; j++ )
  {
    y_max = r8_max ( y_max, t[1+j*2] );
    y_min = r8_min ( y_min, t[1+j*2] );
  }
  y_scale = y_max - y_min;
  y_max = y_max + 0.05 * y_scale;
  y_min = y_min - 0.05 * y_scale;
  y_scale = y_max - y_min;

  i4_min = 1;
  j4_min = 1;
  if ( x_scale < y_scale )
  {
    i4_max = ( int ) ( 0.5 + 500.0 * x_scale / y_scale );
    j4_max = 500;
  }
  else
  {
    i4_max = 500;
    j4_max = ( int ) ( 0.5 + 500.0 * y_scale / x_scale );
  }
/*
  Open the file.
*/
  output = fopen ( plot_filename, "wt");
/*
  Write that junk.
*/
  fprintf ( output, "<?xml version = \"1.0\" standalone=\"no\"?>\n" );
  fprintf ( output, "\n" );
  fprintf ( output, "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\"\n" );
  fprintf ( output, "  \"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n" );
  fprintf ( output, "\n" );
  fprintf ( output, "<svg\n" );
  fprintf ( output, "  width=\"%d\"\n", i4_max );
  fprintf ( output, "  height=\"%d\"\n", j4_max );
  fprintf ( output, "  viewbox=\"%d %d %d %d\"\n", 
    i4_min, j4_min, i4_max, j4_max);
  fprintf ( output, "  xmlns=\"http://www.w3.org/2000/svg\"\n" );
  fprintf ( output, "  version=\"1.1\">\n" );
  fprintf ( output, "  <desc>\n" );
  fprintf ( output, "    Triangulation created by triangle_svg.c\n" );
  fprintf ( output, "  </desc>\n" );
/*
  Draw the triangle.
*/
  fprintf ( output, "  <polygon\n" );
  fprintf ( output, "    fill=\"pink\"\n" );
  fprintf ( output, "    stroke=\"black\"\n" );
  fprintf ( output, "    stroke-width=\"2\"\n" );
  fprintf ( output, "    points=\"\n" );

  for ( j = 0; j < 3; j++ )
  {
    i4 = r8_to_i4 ( x_min, x_max, t[0+j*2], i4_min, i4_max );
    j4 = r8_to_i4 ( y_max, y_min, t[1+j*2], j4_min, j4_max );
    fprintf ( output, "      %d,%d\n", i4, j4 );
  }
  fprintf ( output, "  \" />\n" );
/*
  Draw points.
*/
  for ( j = 0; j < p_num; j++ )
  {
    i4 = r8_to_i4 ( x_min, x_max, p[0+j*2], i4_min, i4_max );
    j4 = r8_to_i4 ( y_max, y_min, p[1+j*2], j4_min, j4_max );
    r = 5;

    fprintf ( output, "  <circle\n" );
    fprintf ( output, "    cx=\"%d\"\n", i4 );
    fprintf ( output, "    cy=\"%d\"\n", j4 );
    fprintf ( output, "    r=\"%d\"\n", r );
    fprintf ( output, "    fill=\"blue\"\n" );
    fprintf ( output, "    stroke=\"black\"\n" );
    fprintf ( output, "    stroke-width=\"2\"\n" );
    fprintf ( output, "  />\n" );
  }
/*
  End of plot.
*/
  fprintf ( output, "</svg>\n" );

  fclose ( output );

  printf ( "\n" );
  printf ( "  Graphics data written to file \"%s\"\n", plot_filename );

  return;
}
