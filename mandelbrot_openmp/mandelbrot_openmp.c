# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>

# include <omp.h>

int main ( void );
int i4_min ( int i1, int i2 );
void timestamp ( void );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose

    MAIN is the main program for MANDELBROT_OPENMP.

  Discussion:

    MANDELBROT_OPENMP computes an image of the Mandelbrot set.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    03 September 2012

  Author:

    John Burkardt

  Local Parameters:

    Local, int COUNT_MAX, the maximum number of iterations taken
    for a particular pixel.
*/
{
  int m = 500;
  int n = 500;

  int b[m][n];
  int c;
  int c_max;
  int count[m][n];
  int count_max = 2000;
  int g[m][n];
  int i;
  int ierror;
  int j;
  int jhi;
  int jlo;
  int k;
  char *output_filename = "mandelbrot.ppm";
  FILE *output_unit;
  int r[m][n];
  double wtime;
  double wtime_total;
  double x_max =   1.25;
  double x_min = - 2.25;
  double x;
  double x1;
  double x2;
  double y_max =   1.75;
  double y_min = - 1.75;
  double y;
  double y1;
  double y2;

  timestamp ( );
  printf ( "\n" );
  printf ( "MANDELBROT_OPENMP\n" );
  printf ( "  C/OpenMP version\n" );
  printf ( "\n" );
  printf ( "  Create an ASCII PPM image of the Mandelbrot set.\n" );
  printf ( "\n" );
  printf ( "  For each point C = X + i*Y\n" );
  printf ( "  with X range [%g,%g]\n", x_min, x_max );
  printf ( "  and  Y range [%g,%g]\n", y_min, y_max );
  printf ( "  carry out %d iterations of the map\n", count_max );
  printf ( "  Z(n+1) = Z(n)^2 + C.\n" );
  printf ( "  If the iterates stay bounded (norm less than 2)\n" );
  printf ( "  then C is taken to be a member of the set.\n" );
  printf ( "\n" );
  printf ( "  An ASCII PPM image of the set is created using\n" );
  printf ( "    M = %d pixels in the X direction and\n", m );
  printf ( "    N = %d pixels in the Y direction.\n", n );

  wtime = omp_get_wtime ( );
/*
  Carry out the iteration for each pixel, determining COUNT.
*/
# pragma omp parallel \
  shared ( b, count, count_max, g, r, x_max, x_min, y_max, y_min ) \
  private ( i, j, k, x, x1, x2, y, y1, y2 )
{
# pragma omp for

  for ( i = 0; i < m; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      x = ( ( double ) (     j - 1 ) * x_max   
          + ( double ) ( m - j     ) * x_min ) 
          / ( double ) ( m     - 1 );

      y = ( ( double ) (     i - 1 ) * y_max   
          + ( double ) ( n - i     ) * y_min ) 
          / ( double ) ( n     - 1 );

      count[i][j] = 0;

      x1 = x;
      y1 = y;

      for ( k = 1; k <= count_max; k++ )
      {
        x2 = x1 * x1 - y1 * y1 + x;
        y2 = 2 * x1 * y1 + y;

        if ( x2 < -2.0 || 2.0 < x2 || y2 < -2.0 || 2.0 < y2 )
        {
          count[i][j] = k;
          break;
        }
        x1 = x2;
        y1 = y2;
      }

      if ( ( count[i][j] % 2 ) == 1 )
      {
        r[i][j] = 255;
        g[i][j] = 255;
        b[i][j] = 255;
      }
      else
      {
        c = ( int ) ( 255.0 * sqrt ( sqrt ( sqrt ( 
          ( ( double ) ( count[i][j] ) / ( double ) ( count_max ) ) ) ) ) );
        r[i][j] = 3 * c / 5;
        g[i][j] = 3 * c / 5;
        b[i][j] = c;
      }
    }
  }
}

  wtime = omp_get_wtime ( ) - wtime;
  printf ( "\n" );
  printf ( "  Time = %g seconds.\n", wtime );
/*
  Write data to an ASCII PPM file.
*/
  output_unit = fopen ( output_filename, "wt" );

  fprintf ( output_unit, "P3\n" );
  fprintf ( output_unit, "%d  %d\n", n, m );
  fprintf ( output_unit, "%d\n", 255 );
  for ( i = 0; i < m; i++ )
  {
    for ( jlo = 0; jlo < n; jlo = jlo + 4 )
    {
      jhi = i4_min ( jlo + 4, n );
      for ( j = jlo; j < jhi; j++ )
      {
        fprintf ( output_unit, "  %d  %d  %d", r[i][j], g[i][j], b[i][j] );
      }
      fprintf ( output_unit, "\n" );
    }
  }

  fclose ( output_unit );
  printf ( "\n" );
  printf ( "  Graphics data written to \"%s\".\n", output_filename );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "MANDELBROT_OPENMP\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

int i4_min ( int i1, int i2 )

/******************************************************************************/
/*
  Purpose:

    I4_MIN returns the smaller of two I4's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    29 August 2006

  Author:

    John Burkardt

  Parameters:

    Input, int I1, I2, two integers to be compared.

    Output, int I4_MIN, the smaller of I1 and I2.
*/
{
  int value;

  if ( i1 < i2 )
  {
    value = i1;
  }
  else
  {
    value = i2;
  }
  return value;
}
/******************************************************************************/

void timestamp ( void )

/******************************************************************************/
/*
  Purpose:

    TIMESTAMP prints the current YMDHMS date as a time stamp.

  Example:

    31 May 2001 09:45:54 AM

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
