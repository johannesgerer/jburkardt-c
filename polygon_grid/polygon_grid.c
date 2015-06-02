# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <string.h>

# include "polygon_grid.h"

/******************************************************************************/

int polygon_grid_count ( int n, int nv )

/******************************************************************************/
/*
  Purpose:

    POLYGON_GRID_COUNT counts the grid points inside a polygon.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 May 2015

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of subintervals on a side.

    Input, int NV, the number of vertices.
    3 <= NV.

    Output, int POLYGON_GRID_COUNT, the number of grid points.
*/
{
  int ng;

  ng = 1 + nv * ( n * ( n + 1 ) ) / 2;

  return ng;
}

/******************************************************************************/

void polygon_grid_display ( int n, int nv, double v[], int ng, double xg[], 
  char *prefix )

/******************************************************************************/
/*
  Purpose:

    POLYGON_GRID_DISPLAY displays grid points inside a polygon.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    10 May 2015

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of subintervals.

    Input, int NV, the number of vertices in the polygon.

    Input, double V[2*NV], the coordinates of the vertices.

    Input, int NG, the number of grid points.

    Input, double XG[2*NG], the grid points.

    Input, char *PREFIX, a string used to name the files.
*/
{
  FILE *command_unit;
  char command_filename[255];
  FILE *grid_unit;
  char grid_filename[255];
  int j;
  char plot_filename[255];
  double vc[2];
  FILE *vertex_unit;
  char vertex_filename[255];
/*
  Determine the centroid.
*/
  vc[0] = 0.0;
  vc[1] = 0.0;
  for ( j = 0; j < nv; j++ )
  {
    vc[0] = vc[0] + v[0+j*2];
    vc[1] = vc[1] + v[1+j*2];
  }
  vc[0] = vc[0] / ( double ) ( nv );
  vc[1] = vc[1] / ( double ) ( nv );
/*
  Write the vertex file.
*/
  strcpy ( vertex_filename, prefix );
  strcat ( vertex_filename, "_vertex.txt" );
  vertex_unit = fopen ( vertex_filename, "wt" );

  for ( j = 0; j < nv; j++ )
  {
    fprintf ( vertex_unit, "  %14.6g  %14.6g\n", v[0+2*j], v[1+2*j] );
  }
  fprintf ( vertex_unit, "  %14.6g  %14.6g\n", v[0+0*j], v[1+0*j] );
  for ( j = 0; j < nv; j++ )
  {
    fprintf ( vertex_unit, "\n" );
    fprintf ( vertex_unit, "  %14.6g  %14.6g\n", v[0+j*2], v[1+j*2] );
    fprintf ( vertex_unit, "  %14.6g  %14.6g\n", vc[0], vc[1] );
  }
  fclose ( vertex_unit );
  printf ( "\n" );
  printf ( "  Created vertex file '%s'\n", vertex_filename );
/*
  Write the gridpoint file.
*/
  strcpy ( grid_filename, prefix );
  strcat ( grid_filename, "_grid.txt" );
  grid_unit = fopen ( grid_filename, "wt" );
  for ( j = 0; j < ng; j++ )
  {
    fprintf ( grid_unit, "  %14.6g  %14.6g\n", xg[0+j*2], xg[1+j*2] );
  }
  fclose ( grid_unit );
  printf ( "\n" );
  printf ( "  Created grid file '%s'\n", grid_filename );
/*
  Write the command file.
*/
  strcpy ( plot_filename, prefix );
  strcat ( plot_filename, ".png" );

  strcpy ( command_filename, prefix );
  strcat ( command_filename, "_commands.txt" );

  command_unit = fopen ( command_filename, "wt" );
  fprintf ( command_unit, "# %s\n", command_filename );
  fprintf ( command_unit, "#\n" );
  fprintf ( command_unit, "# Usage:\n" );
  fprintf ( command_unit, "#  gnuplot < %s\n", command_filename );
  fprintf ( command_unit, "#\n" );
  fprintf ( command_unit, "set term png\n" );
  fprintf ( command_unit, "set output '%s'\n", plot_filename );
  fprintf ( command_unit, "set xlabel '<--- X --->'\n" );
  fprintf ( command_unit, "set ylabel '<--- Y --->'\n" );
  fprintf ( command_unit, "set title '%s'\n", prefix );
  fprintf ( command_unit, "set grid\n" );
  fprintf ( command_unit, "set key off\n" );
  fprintf ( command_unit, "set size ratio -1\n" );
  fprintf ( command_unit, "set style data lines\n" );

  fprintf ( command_unit, 
    "plot '%s' using 1:2 with points lt 3 pt 3,\\\n", grid_filename );
  fprintf ( command_unit, 
    "     '%s' using 1:2 lw 3 linecolor rgb 'black'\n", vertex_filename );
  fprintf ( command_unit, "quit\n" );
  fclose ( command_unit );

  printf ( "  Created command file '%s'\n", command_filename );

  return;
}
/******************************************************************************/

double *polygon_grid_points ( int n, int nv, double v[], int ng )

/******************************************************************************/
/*
  Purpose:

    POLYGON_GRID_POINTS computes points on a polygonal grid.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 May 2015

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of subintervals.

    Input, int NV, the number of vertices in the polygon.

    Input, double V[2*NV], the coordinates of the vertices.

    Input, int NG, the number of grid points.

    Output, double POLYGON_GRID_POINTS[2*NG], the coordinates of the 
    grid points.
*/
{
  int i;
  int j;
  int k;
  int l;
  int lp1;
  int p;
  double vc[2];
  double *xg;

  xg = ( double * ) malloc ( 2 * ng * sizeof ( double ) );
  p = 0;
/*
  Determine the centroid.
*/
  vc[0] = 0.0;
  vc[1] = 0.0;
  for ( j = 0; j < nv; j++ )
  {
    vc[0] = vc[0] + v[0+j*2];
    vc[1] = vc[1] + v[1+j*2];
  }
  vc[0] = vc[0] / ( double ) ( nv );
  vc[1] = vc[1] / ( double ) ( nv );
/*
  The centroid is the first point.
*/
  xg[0+p*2] = vc[0];
  xg[1+p*2] = vc[1];
  p = p + 1;
/*
  Consider each triangle formed by two consecutive vertices and the centroid,
  but skip the first line of points.
*/
  for ( l = 0; l < nv; l++ )
  {
    lp1 = ( ( l + 1 ) % nv );
    for ( i = 1; i <= n; i++ )
    {
      for ( j = 0; j <= n - i; j++ )
      {
        k = n - i - j;
        xg[0+p*2] = ( ( double ) ( i ) * v[0+l*2]   
                    + ( double ) ( j ) * v[0+lp1*2] 
                    + ( double ) ( k ) * vc[0] )  
                    / ( double ) ( n );
        xg[1+p*2] = ( ( double ) ( i ) * v[1+l*2]   
                    + ( double ) ( j ) * v[1+lp1*2] 
                    + ( double ) ( k ) * vc[1] )  
                    / ( double ) ( n );
        p = p + 1;
      }
    }
  }

  return xg;
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
    fprintf ( stderr, "  Could not open the file '%s'.\n", output_filename );
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
