# include <math.h>
# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <time.h>

# include "sphere_llt_grid.h"

/******************************************************************************/

void i4mat_transpose_print ( int m, int n, int a[], char *title )

/******************************************************************************/
/*
  Purpose:

    I4MAT_TRANSPOSE_PRINT prints an I4MAT, transposed.

  Discussion:

    An I4MAT is an MxN array of I4's, stored by (I,J) -> [I+J*M].

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 January 2005

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows in A.

    Input, int N, the number of columns in A.

    Input, int A[M*N], the M by N matrix.

    Input, char *TITLE, a title.
*/
{
  i4mat_transpose_print_some ( m, n, a, 1, 1, m, n, title );

  return;
}
/******************************************************************************/

void i4mat_transpose_print_some ( int m, int n, int a[], int ilo, int jlo,
  int ihi, int jhi, char *title )

/******************************************************************************/
/*
  Purpose:

    I4MAT_TRANSPOSE_PRINT_SOME prints some of an I4MAT, transposed.

  Discussion:

    An I4MAT is an MxN array of I4's, stored by (I,J) -> [I+J*M].

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 August 2010

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows of the matrix.
    M must be positive.

    Input, int N, the number of columns of the matrix.
    N must be positive.

    Input, int A[M*N], the matrix.

    Input, int ILO, JLO, IHI, JHI, designate the first row and
    column, and the last row and column to be printed.

    Input, char *TITLE, a title.
*/
{
# define INCX 10

  int i;
  int i2hi;
  int i2lo;
  int j;
  int j2hi;
  int j2lo;

  fprintf ( stdout, "\n" );
  fprintf ( stdout, "%s\n", title );

  if ( m <= 0 || n <= 0 )
  {
    fprintf ( stdout, "\n" );
    fprintf ( stdout, "  (None)\n" );
    return;
  }
/*
  Print the columns of the matrix, in strips of INCX.
*/
  for ( i2lo = ilo; i2lo <= ihi; i2lo = i2lo + INCX )
  {
    i2hi = i2lo + INCX - 1;
    if ( m < i2hi )
    {
      i2hi = m;
    }
    if ( ihi < i2hi )
    {
      i2hi = ihi;
    }

    fprintf ( stdout, "\n" );
/*
  For each row I in the current range...

  Write the header.
*/
    fprintf ( stdout, "  Row: " );
    for ( i = i2lo; i <= i2hi; i++ )
    {
      fprintf ( stdout, "%6d  ", i - 1 );
    }
    fprintf ( stdout, "\n" );
    fprintf ( stdout, "  Col\n" );
    fprintf ( stdout, "\n" );
/*
  Determine the range of the rows in this strip.
*/
    j2lo = jlo;
    if ( j2lo < 1 )
    {
      j2lo = 1;
    }
    j2hi = jhi;
    if ( n < jhi )
    {
      j2hi = n;
    }

    for ( j = j2lo; j <= j2hi; j++ )
    {
/*
  Print out (up to INCX) entries in column J, that lie in the current strip.
*/
      fprintf ( stdout, "%5d: ", j - 1 );
      for ( i = i2lo; i <= i2hi; i++ )
      {
        fprintf ( stdout, "%6d  ", a[i-1+(j-1)*m] );
      }
      fprintf ( stdout, "\n" );
    }
  }

  return;
# undef INCX
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

void sphere_llt_grid_display ( int ng, double xg[], int line_num, 
  int line_data[], char *prefix )

/******************************************************************************/
/*
  Purpose:

    SPHERE_LLT_GRID_DISPLAY displays a latitude/longitude triangle grid.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 May 2015

  Author:

    John Burkardt

  Parameters:

    Input, int NG, the number of points.

    Input, double XG[3*NG], the points.

    Input, int LINE_NUM, the number of grid lines.

    Input, inte LINE_DATA[2*LINE_NUM], contains pairs of 
    point indices for line segments that make up the grid.

    Input, char *PREFIX, a prefix for the filenames.
*/
{
  char command_filename[255];
  FILE *command_unit;
  int i;
  int j;
  int j1;
  int j2;
  int l;
  char line_filename[255];
  FILE *line_unit;
  char node_filename[255];
  FILE *node_unit;
  char plot_filename[255];
/*
  Create graphics data files.
*/
  strcpy ( node_filename, prefix );
  strcat ( node_filename, "_nodes.txt" );

  node_unit = fopen ( node_filename, "wt" );
  for ( j = 0; j < ng; j++ )
  {
    for ( i = 0; i < 3; i++ )
    {
      fprintf ( node_unit, "  %g", xg[i+j*3] );
    }
    fprintf ( node_unit, "\n" );
  }
  fclose ( node_unit );
  printf ( "\n" );
  printf ( "  Created node file '%s'\n", node_filename );

  strcpy ( line_filename, prefix );
  strcat ( line_filename, "_lines.txt" );

  line_unit = fopen ( line_filename, "wt" );
  for ( l = 0; l < line_num; l++ )
  {
    if ( 0 < l )
    {
      fprintf ( line_unit, "\n" );
      fprintf ( line_unit, "\n" );
    }
    j1 = line_data[0+l*2];
    j2 = line_data[1+l*2];
    fprintf ( line_unit, "  %g  %g  %g\n", xg[0+j1*3], xg[1+j1*3], xg[2+j1*3] );
    fprintf ( line_unit, "  %g  %g  %g\n", xg[0+j2*3], xg[1+j2*3], xg[2+j2*3] );
  }
  fclose ( line_unit );
  printf ( "\n" );
  printf ( "  Created line file '%s'\n", line_filename );
/*
  Create graphics command file.
*/
  strcpy ( command_filename, prefix );
  strcat ( command_filename, "_commands.txt" );

  command_unit = fopen ( command_filename, "wt" );
  fprintf ( command_unit, "# %s\n", command_filename );
  fprintf ( command_unit, "#\n" );
  fprintf ( command_unit, "# Usage:\n" );
  fprintf ( command_unit, "#  gnuplot < %s\n", command_filename );
  fprintf ( command_unit, "#\n" );
  fprintf ( command_unit, "set term png\n" );

  strcpy ( plot_filename, prefix );
  strcat ( plot_filename, ".png" );

  fprintf ( command_unit, "set output '%s'\n", plot_filename );
  fprintf ( command_unit, "set xlabel '<--- X --->'\n" );
  fprintf ( command_unit, "set ylabel '<--- Y --->'\n" );
  fprintf ( command_unit, "set zlabel '<--- Z --->'\n" );
  fprintf ( command_unit, "set title '%s'\n", prefix );
  fprintf ( command_unit, "set grid\n" );
  fprintf ( command_unit, "set key off\n" );
  fprintf ( command_unit, "set style data points\n" );
  fprintf ( command_unit, "set timestamp\n" );
  fprintf ( command_unit, "set view equal xyz\n" );
  fprintf ( command_unit, "splot '%s' with lines lw 3, \\\n", line_filename );
  fprintf ( command_unit, "      '%s' with points pt 7 lt 0\n", node_filename );
  fprintf ( command_unit, "quit\n" );
  fclose ( command_unit );

  printf ( "  Created command file '%s'\n", command_filename );

  return;
}
/******************************************************************************/

int sphere_llt_grid_line_count ( int lat_num, int long_num )

/******************************************************************************/
/*
  Purpose:

    SPHERE_LLT_GRID_LINE_COUNT counts latitude/longitude triangle grid lines.

  Discussion:

    A SPHERE LLT grid imposes a grid of triangles on a sphere,
    using latitude and longitude lines.

    The number returned is the number of pairs of points to be connected.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    20 May 2015

  Author:

    John Burkardt

  Parameters:

    Input, int LAT_NUM, LONG_NUM, the number of latitude and
    longitude lines to draw.  The latitudes do not include the North and South
    poles, which will be included automatically, so LAT_NUM = 5, for instance,
    will result in points along 7 lines of latitude.

    Output, int SPHERE_LLT_GRID_LINE_COUNT, the number of grid lines.
*/
{
  int line_num;

  line_num = long_num * ( lat_num + 1 ) 
           + long_num *   lat_num 
           + long_num * ( lat_num - 1 );

  return line_num;
}
/******************************************************************************/

int *sphere_llt_grid_lines ( int nlat, int nlong, int line_num )

/******************************************************************************/
/*
  Purpose:

    SPHERE_LLT_GRID_LINES: latitude/longitude triangle grid lines.

  Discussion:

    A SPHERE LLT grid imposes a grid of triangles on a sphere,
    using latitude and longitude lines.

    The point numbering system is the same used in SPHERE_LLT_POINTS,
    and that routine may be used to compute the coordinates of the points.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    20 May 2015

  Author:

    John Burkardt

  Parameters:

    Input, int NLAT, NLONG, the number of latitude and longitude
    lines to draw.  The latitudes do not include the North and South
    poles, which will be included automatically, so NLAT = 5, for instance,
    will result in points along 7 lines of latitude.

    Input, int LINE_NUM, the number of grid lines.

    Output, int SPHERE_LLT_GRID_LINES[2*LINE_NUM], contains pairs of point 
    indices for line segments that make up the grid.
*/
{
  int i;
  int j;
  int l;
  int *line;
  int next;
  int newcol;
  int old;

  line = ( int * ) malloc ( 2 * line_num * sizeof ( int ) );

  l = 0;
/*
  "Vertical" lines.
*/
  for ( j = 0; j <= nlong - 1; j++ )
  {
    old = 0;
    next = j + 1;
    line[0+l*2] = old;
    line[1+l*2] = next;
    l = l + 1;

    for ( i = 1; i <= nlat-1; i++ )
    {
      old = next;
      next = old + nlong;
      line[0+l*2] = old;
      line[1+l*2] = next;
      l = l + 1;
    }

    old = next;
    line[0+l*2] = old;
    line[1+l*2] = 1 + nlat * nlong;
    l = l + 1;
  }
/*
  "Horizontal" lines.
*/
  for ( i = 1; i <= nlat; i++ )
  {
    next = ( i - 1 ) * nlong + 1;

    for ( j = 0; j <= nlong - 2; j++ )
    {
      old = next;
      next = old + 1;
      line[0+l*2] = old;
      line[1+l*2] = next;
      l = l + 1;
    }

    old = next;
    next = ( i - 1 ) * nlong + 1;
    line[0+l*2] = old;
    line[1+l*2] = next;
    l = l + 1;
  }
/*
  "Diagonal" lines.
*/
  for ( j = 0; j < nlong; j++ )
  {
    old = 0;
    next = j + 1;
    newcol = j;

    for ( i = 1; i < nlat; i++ )
    {
      old = next;
      next = old + nlong + 1;
      newcol = newcol + 1;
      if ( nlong - 1 < newcol )
      {
        newcol = 0;
        next = next - nlong;
      }

      line[0+l*2] = old;
      line[1+l*2] = next;
      l = l + 1;
    }
  }

  return line;
}
/******************************************************************************/

int sphere_llt_grid_point_count ( int lat_num, int long_num )

/******************************************************************************/
/*
  Purpose:

    SPHERE_LLT_GRID_POINT_COUNT counts points for a latitude/longitude grid.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    20 May 2015

  Author:

    John Burkardt

  Parameters:

    Input, int LAT_NUM, LONG_NUM, the number of latitude 
    and longitude lines to draw.  The latitudes do not include the North and 
    South poles, which will be included automatically, so LAT_NUM = 5, for 
    instance, will result in points along 7 lines of latitude.

    Output, int SPHERE_LLT_GRID_POINT_COUNT, the number of grid points.
*/
{
  int point_num;

  point_num = 2 + lat_num * long_num;

  return point_num;
}
/******************************************************************************/

double *sphere_llt_grid_points ( double r, double pc[3], int lat_num, 
  int lon_num, int point_num )

/******************************************************************************/
/*
  Purpose:

    SPHERE_LLT_GRID_POINTS produces points on a latitude/longitude grid.

  Discussion:

    A SPHERE LLT grid imposes a grid of triangles on a sphere,
    using latitude and longitude lines.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    20 May 2015

  Author:

    John Burkardt

  Parameters:

    Input, double R, the radius of the sphere.

    Input, double PC[3], the coordinates of the center of the sphere.

    Input, int LAT_NUM, LON_NUM, the number of latitude and longitude
    lines to draw.  The latitudes do not include the North and South
    poles, which will be included automatically, so LAT_NUM = 5, for instance,
    will result in points along 7 lines of latitude.

    Input, int POINT_NUM, the number of points.

    Output, double SPHERE_LLT_GRID_POINTS[3*POINT_NUM], the coordinates 
    of the grid points.
*/
{
  int lat;
  int lon;
  int n;
  double *p;
  double phi;
  double r8_pi = 3.141592653589793;
  double theta;

  p = ( double * ) malloc ( 3 * point_num * sizeof ( double ) );
  n = 0;
/*
  The north pole.
*/
  theta = 0.0;
  phi = 0.0;

  p[0+n*3] = pc[0] + r * sin ( phi ) * cos ( theta );
  p[1+n*3] = pc[1] + r * sin ( phi ) * sin ( theta );
  p[2+n*3] = pc[2] + r * cos ( phi );
  n = n + 1;
/*
  Do each intermediate ring of latitude.
*/
  for ( lat = 1; lat <= lat_num; lat++ )
  {
    phi = ( double ) ( lat ) * r8_pi / ( double ) ( lat_num + 1 );
/*
  Along that ring of latitude, compute points at various longitudes.
*/
    for ( lon = 0; lon < lon_num; lon++ )
    {
      theta = ( double ) ( lon ) * 2.0 * r8_pi / ( double ) ( lon_num );

      p[0+n*3] = pc[0] + r * sin ( phi ) * cos ( theta );
      p[1+n*3] = pc[1] + r * sin ( phi ) * sin ( theta );
      p[2+n*3] = pc[2] + r * cos ( phi );
      n = n + 1;
    }
  }
/*
  The south pole.
*/
  theta = 0.0;
  phi = r8_pi;
  p[0+n*3] = pc[0] + r * sin ( phi ) * cos ( theta );
  p[1+n*3] = pc[1] + r * sin ( phi ) * sin ( theta );
  p[2+n*3] = pc[2] + r * cos ( phi );
  n = n + 1;

  return p;
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
