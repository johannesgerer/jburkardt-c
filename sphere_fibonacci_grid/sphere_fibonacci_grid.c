# include <math.h>
# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <time.h>

# include "sphere_fibonacci_grid.h"

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

void sphere_fibonacci_grid_display ( int ng, double xg[], char *prefix )

/******************************************************************************/
/*
  Purpose:

    SPHERE_FIBONACCI_GRID_DISPLAY displays sphere points on a Fibonacci spiral.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 May 2015

  Author:

    John Burkardt

  Parameters:

    Input, int NG, the number of points.

    Input, double XG[3*NG], the Fibonacci spiral points.

    Input, char *PREFIX, a prefix for the filenames.
*/
{
  char command_filename[255];
  FILE *command_unit;
  char data_filename[255];
  FILE *data_unit;
  int i;
  int j;
  char plot_filename[255];
/*
  Create graphics data file.
*/
  strcpy ( data_filename, prefix );
  strcat ( data_filename, "_data.txt" );

  data_unit = fopen ( data_filename, "wt" );
  for ( j = 0; j < ng; j++ )
  {
    for ( i = 0; i < 3; i++ )
    {
      fprintf ( data_unit, "  %g", xg[i+j*3] );
    }
    fprintf ( data_unit, "\n" );
  }
  fclose ( data_unit );
  printf ( "\n" );
  printf ( "  Created data file '%s'\n", data_filename );
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
  fprintf ( command_unit, "splot '%s'\n", data_filename );
  fprintf ( command_unit, "quit\n" );
  fclose ( command_unit );

  printf ( "  Created command file '%s'\n", command_filename );

  return;
}
/******************************************************************************/

double *sphere_fibonacci_grid_points ( int ng )

/******************************************************************************/
/*
  Purpose:

    SPHERE_FIBONACCI_GRID_POINTS computes sphere points on a Fibonacci spiral.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 May 2015

  Author:

    John Burkardt

  Reference:

    Richard Swinbank, James Purser,
    Fibonacci grids: A novel approach to global modelling,
    Quarterly Journal of the Royal Meteorological Society,
    Volume 132, Number 619, July 2006 Part B, pages 1769-1793.

  Parameters:

    Input, int NG, the number of points.

    Output, double SPHERE_FIBONACCI_GRID_POINTS[3*NG], the Fibonacci 
  spiral points.
*/
{
  double cphi;
  int i;
  double i_r8;
  int j;
  double ng_r8;
  double r8_phi;
  const double r8_pi = 3.141592653589793;
  double sphi;
  double theta;
  double *xyz;

  xyz = ( double * ) malloc ( 3 * ng * sizeof ( double ) );

  r8_phi = ( 1.0 + sqrt ( 5.0 ) ) / 2.0;
  ng_r8 = ( double ) ( ng );

  for ( j = 0; j < ng; j++ )
  {
    i_r8 = ( double ) ( - ng + 1 + 2 * j );
    theta = 2.0 * r8_pi * i_r8 / r8_phi;
    sphi = i_r8 / ng_r8;
    cphi = sqrt ( ( ng_r8 + i_r8 ) * ( ng_r8 - i_r8 ) ) / ng_r8;
    xyz[0+j*3] = cphi * sin ( theta );
    xyz[1+j*3] = cphi * cos ( theta );
    xyz[2+j*3] = sphi;
  }

  return xyz;
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

