# include <stdlib.h>
# include <stdio.h>
# include <string.h>
# include <time.h>

# include "pyramid_grid.h"

/******************************************************************************/

int pyramid_grid_count ( int n )

/******************************************************************************/
/*
  Purpose:

    PYRAMID_GRID_COUNT counts the points in a pyramid grid.

  Discussion:

    0:  x

    1:  x  x
        x  x

    2:  x  x  x
        x  x  x
        x  x  x

    3:  x  x  x  x
        x  x  x  x
        x  x  x  x
        x  x  x  x

    N  Size

    0     1
    1     5 = 1 + 4
    2    14 = 1 + 4 + 9
    3    30 = 1 + 4 + 9 + 16
    4    55 = 1 + 4 + 9 + 16 + 25
    5    91 = 1 + 4 + 9 + 16 + 25 + 36

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 August 2014

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of subintervals.

    Output, int PYRAMID_GRID_COUNT, the number of
    nodes in the grid of size N.
*/
{
  int np1;
  int value;

  np1 = n + 1;

  value = ( np1 * ( np1 + 1 ) * ( 2 * np1 + 1 ) ) / 6;

  return value;
}
/******************************************************************************/

double *pyramid_unit_grid ( int n, int ng )

/******************************************************************************/
/*
  Purpose:

    PYRAMID_UNIT_GRID computes grid points in the unit pyramid.

  Discussion:

    The unit pyramid has base (-1,-1,0), (+1,1,0), (+1,+1,0), (-1,+1,0)
    and vertex (0,0,1).

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 August 2014

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of subintervals.

    Input, int NG, the number of nodes to generate,
    as determined by pyramid_grid_count().

    Output, double PYRAMID_UNIT_GRID[3*NG], the grid point coordinates.
*/
{
  int g;
  int hi;
  int i;
  int j;
  int k;
  int lo;
  double *pg;

  pg = ( double * ) malloc ( 3 * ng * sizeof ( double ) );

  g = 0;

  for ( k = n; 0 <= k; k-- )
  {
    hi = n - k;
    lo = - hi;
    for ( j = lo; j <= hi; j = j + 2 )
    {
      for ( i = lo; i <= hi; i = i + 2 )
      {
        pg[0+g*3] = ( double ) ( i ) / ( double ) ( n );
        pg[1+g*3] = ( double ) ( j ) / ( double ) ( n );
        pg[2+g*3] = ( double ) ( k ) / ( double ) ( n );
        g = g + 1;
      }
    }
  }

  return pg;
}
/******************************************************************************/

void pyramid_unit_grid_plot ( int n, int ng, double pg[], char *header )

/******************************************************************************/
/*
  Purpose:

    PYRAMID_UNIT_GRID_PLOT sets up a GNUPLOT plot of a unit pyramid grid.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 August 2014

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of subintervals.

    Input, int NG, the number of nodes to generate,
    as determined by pyramid_grid_count().

    Input, double PG[3*NG], the grid point coordinates.

    Input, char *HEADER, the header for the files.
*/
{
  char command_filename[255];
  FILE *command_unit;
  int j;
  char node_filename[255];
  FILE *node_unit;
  char plot_filename[255];
  double v1[3];
  double v2[3];
  double v3[3];
  double v4[3];
  double v5[3];
  char vertex_filename[255];
  FILE *vertex_unit;
/*
  Create the vertex file.
*/
  pyramid_unit_vertices ( v1, v2, v3, v4, v5 );

  strcpy ( vertex_filename, header );
  strcat ( vertex_filename, "_vertices.txt" );

  vertex_unit = fopen ( vertex_filename, "wt" );

  fprintf ( vertex_unit, "%g  %g  %g\n", v2[0], v2[1], v2[2] );
  fprintf ( vertex_unit, "%g  %g  %g\n", v3[0], v3[1], v3[2] );
  fprintf ( vertex_unit, "%g  %g  %g\n", v4[0], v4[1], v4[2] );
  fprintf ( vertex_unit, "%g  %g  %g\n", v5[0], v5[1], v5[2] );
  fprintf ( vertex_unit, "%g  %g  %g\n", v2[0], v2[1], v2[2] );
  fprintf ( vertex_unit, "\n" );

  fprintf ( vertex_unit, "%g  %g  %g\n", v1[0], v1[1], v1[2] );
  fprintf ( vertex_unit, "%g  %g  %g\n", v2[0], v2[1], v2[2] );
  fprintf ( vertex_unit, "\n" );

  fprintf ( vertex_unit, "%g  %g  %g\n", v1[0], v1[1], v1[2] );
  fprintf ( vertex_unit, "%g  %g  %g\n", v3[0], v3[1], v3[2] );
  fprintf ( vertex_unit, "\n" );

  fprintf ( vertex_unit, "%g  %g  %g\n", v1[0], v1[1], v1[2] );
  fprintf ( vertex_unit, "%g  %g  %g\n", v4[0], v4[1], v4[2] );
  fprintf ( vertex_unit, "\n" );

  fprintf ( vertex_unit, "%g  %g  %g\n", v1[0], v1[1], v1[2] );
  fprintf ( vertex_unit, "%g  %g  %g\n", v5[0], v5[1], v5[2] );

  fclose ( vertex_unit );

  printf ( "\n" );
  printf ( "  Created vertex file '%s'\n", vertex_filename );
/*
  Create the node file.
*/
  strcpy ( node_filename, header );
  strcat ( node_filename, "_nodes.txt" );

  node_unit = fopen ( node_filename, "wt" );

  for ( j = 0; j < ng; j++ )
  {
    fprintf ( node_unit, "%g  %g  %g\n", pg[0+j*3], pg[1+j*3], pg[2+j*3] );
  }
  fclose ( node_unit );

  printf ( " Created node file '%s'\n", node_filename );
/*
  Create the command file.
*/
  strcpy ( command_filename, header );
  strcat ( command_filename, "_commands.txt" );

  command_unit = fopen ( command_filename, "wt" );

  fprintf ( command_unit, "# %s\n", command_filename );
  fprintf ( command_unit, "#\n" );
  fprintf ( command_unit, "# Usage:\n" );
  fprintf ( command_unit, "#  gnuplot < %s\n", command_filename );
  fprintf ( command_unit, "#\n" );
  fprintf ( command_unit, "set term png\n" );

  strcpy ( plot_filename, header );
  strcat ( plot_filename, ".png" );

  fprintf ( command_unit, "set output '%s'\n", plot_filename );
  fprintf ( command_unit, "set xlabel '<--- X --->'\n" );
  fprintf ( command_unit, "set ylabel '<--- Y --->'\n" );
  fprintf ( command_unit, "set zlabel '<--- Z --->'\n" );
  fprintf ( command_unit, "set title '%s'\n", header );
  fprintf ( command_unit, "set grid\n" );
  fprintf ( command_unit, "set key off\n" );
  fprintf ( command_unit, "set view equal xyz\n" );
  fprintf ( command_unit, "set view 80, 40\n" );
  fprintf ( command_unit, "set style data lines\n" );
  fprintf ( command_unit, "set timestamp\n" );
  fprintf ( command_unit, "splot '%s' with lines lw 3, \\\n", vertex_filename );
  fprintf ( command_unit, "      '%s' with points pt 7 lt 0\n", node_filename );

  fclose ( command_unit );

  printf ( "  Created command file '%s'\n", command_filename );

  return;
}
/******************************************************************************/

void pyramid_unit_vertices ( double v1[], double v2[], double v3[], 
  double v4[], double v5[] )

/******************************************************************************/
/*
  Purpose:

    PYRAMID_UNIT_VERTICES returns the vertices of the unit pyramid.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 August 2014

  Author:

    John Burkardt

  Parameters:

    Output, double V1[3], V2[3], V3[3], V4[3], V5[3], the vertices.
*/
{
  v1[0] =  0.0;
  v1[1] =  0.0;
  v1[2] = +1.0;

  v2[0] = -1.0;
  v2[1] = -1.0;
  v2[2] =  0.0;

  v3[0] = +1.0;
  v3[1] = -1.0;
  v3[2] =  0.0;

  v4[0] = +1.0;
  v4[1] = +1.0;
  v4[2] =  0.0;

  v5[0] = -1.0;
  v5[1] = +1.0;
  v5[2] =  0.0;

  return;
}
/******************************************************************************/

void r8_print ( double r, char *title )

/******************************************************************************/
/*
  Purpose:

    R8_PRINT prints an R8.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 August 2014

  Author:

    John Burkardt

  Parameters:

    Input, double R, the value to print.

    Input, char *TITLE, a title.
*/
{
  printf ( "%s  %g\n", title, r );

  return;
}
/******************************************************************************/

void r8mat_transpose_print ( int m, int n, double a[], char *title )

/******************************************************************************/
/*
  Purpose:

    R8MAT_TRANSPOSE_PRINT prints an R8MAT, transposed.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    28 May 2008

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns.

    Input, double A[M*N], an M by N matrix to be printed.

    Input, char *TITLE, a title.
*/
{
  r8mat_transpose_print_some ( m, n, a, 1, 1, m, n, title );

  return;
}
/******************************************************************************/

void r8mat_transpose_print_some ( int m, int n, double a[], int ilo, int jlo,
  int ihi, int jhi, char *title )

/******************************************************************************/
/*
  Purpose:

    R8MAT_TRANSPOSE_PRINT_SOME prints some of an R8MAT, transposed.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    10 September 2013

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns.

    Input, double A[M*N], an M by N matrix to be printed.

    Input, int ILO, JLO, the first row and column to print.

    Input, int IHI, JHI, the last row and column to print.

    Input, char *TITLE, a title.
*/
{
# define INCX 5

  int i;
  int i2;
  int i2hi;
  int i2lo;
  int i2lo_hi;
  int i2lo_lo;
  int inc;
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

  if ( ilo < 1 )
  {
    i2lo_lo = 1;
  }
  else
  {
    i2lo_lo = ilo;
  }

  if ( ihi < m )
  {
    i2lo_hi = m;
  }
  else
  {
    i2lo_hi = ihi;
  }

  for ( i2lo = i2lo_lo; i2lo <= i2lo_hi; i2lo = i2lo + INCX )
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

    inc = i2hi + 1 - i2lo;

    fprintf ( stdout, "\n" );
    fprintf ( stdout, "  Row:" );
    for ( i = i2lo; i <= i2hi; i++ )
    {
      fprintf ( stdout, "  %7d     ", i - 1 );
    }
    fprintf ( stdout, "\n" );
    fprintf ( stdout, "  Col\n" );
    fprintf ( stdout, "\n" );

    if ( jlo < 1 )
    {
      j2lo = 1;
    }
    else
    {
      j2lo = jlo;
    }
    if ( n < jhi )
    {
      j2hi = n;
    }
    else
    {
      j2hi = jhi;
    }
    for ( j = j2lo; j <= j2hi; j++ )
    {
      fprintf ( stdout, "%5d:", j - 1 );
      for ( i2 = 1; i2 <= inc; i2++ )
      {
        i = i2lo - 1 + i2;
        fprintf ( stdout, "  %14g", a[(i-1)+(j-1)*m] );
      }
      fprintf ( stdout, "\n" );
    }
  }

  return;
# undef INCX
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

