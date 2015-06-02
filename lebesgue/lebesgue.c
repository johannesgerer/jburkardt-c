# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>

# include "lebesgue.h"

/******************************************************************************/

double *chebyshev1 ( int n )

/******************************************************************************/
/*
  Purpose:

    CHEBYSHEV1 returns the Type 1 Chebyshev points.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    03 March 2018

  Author:

    John Burkardt.

  Parameters:

    Input, int N, the number of points.

    Input, double CHEBYSHEV1[N], the points.
*/
{
  double angle;
  int i;
  const double r8_pi = 3.141592653589793;
  double *x;

  x = ( double * ) malloc ( n * sizeof ( double ) );

  for ( i = 0; i < n; i++ )
  {
    angle = r8_pi * ( double ) ( 2 * i + 1 ) / ( double ) ( 2 * n );
    x[i] = cos ( angle );
  }
  return x;
}
/******************************************************************************/

double *chebyshev2 ( int n )

/******************************************************************************/
/*
  Purpose:

    CHEBYSHEV2 returns the Type 2 Chebyshev points.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    03 March 2018

  Author:

    John Burkardt.

  Parameters:

    Input, int N, the number of points.

    Input, double CHEBYSHEV2[N], the points.
*/
{
  double angle;
  int i;
  const double r8_pi = 3.141592653589793;
  double *x;

  x = ( double * ) malloc ( n * sizeof ( double ) );

  if ( n == 1 )
  {
    x[0] = 0.0;
  }
  else
  {
    for ( i = 0; i < n; i++ )
    {
      angle = r8_pi * ( double ) ( n - i - 1 ) / ( double ) ( n - 1 );
      x[i] = cos ( angle );
    }
  }

  return x;
}
/******************************************************************************/

double *chebyshev3 ( int n )

/******************************************************************************/
/*
  Purpose:

    CHEBYSHEV3 returns the Type 3 Chebyshev points.

  Discussion:

    Note that this point set is NOT symmetric in [-1,+1].
    It is sometimes augmented by the value -1.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    03 March 2018

  Author:

    John Burkardt.

  Parameters:

    Input, int N, the number of points.

    Input, double CHEBYSHEV3[N], the points.
*/
{
  double angle;
  int i;
  const double r8_pi = 3.141592653589793;
  double *x;

  x = ( double * ) malloc ( n * sizeof ( double ) );

  for ( i = 0; i < n; i++ )
  {
    angle = r8_pi * ( double ) ( 2 * n - 2 * i - 1 ) 
                  / ( double ) ( 2 * n         + 1 );
    x[i] = cos ( angle );
  }

  return x;
}
/******************************************************************************/

double *chebyshev4 ( int n )

/******************************************************************************/
/*
  Purpose:

    CHEBYSHEV4 returns the Type 4 Chebyshev points.

  Discussion:

    Note that this point set is NOT symmetric in [-1,+1].
    It is sometimes augmented by the value +1.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    03 March 2018

  Author:

    John Burkardt.

  Parameters:

    Input, int N, the number of points.

    Input, double CHEBYSHEV4[N], the points.
*/
{
  double angle;
  int i;
  const double r8_pi = 3.141592653589793;
  double *x;

  x = ( double * ) malloc ( n * sizeof ( double ) );

  for ( i = 0; i < n; i++ )
  {
    angle = r8_pi * ( double ) ( 2 * n - 2 * i )
                  / ( double ) ( 2 * n + 1 );
    x[i] = cos ( angle );
  }

  return x;
}
/******************************************************************************/

double *equidistant1 ( int n )

/******************************************************************************/
/*
  Purpose:

    EQUIDISTANT1 returns the Type 1 Equidistant points.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    03 March 2018

  Author:

    John Burkardt.

  Parameters:

    Input, int N, the number of points.

    Input, double EQUIDISTANT1[N], the points.
*/
{
  int i;
  double *x;

  x = ( double * ) malloc ( n * sizeof ( double ) );

  for ( i = 0; i < n; i++ )
  {
    x[i] = ( double ) ( - n + 1 + 2 * i ) / ( double ) ( n + 1 );
  }

  return x;
}
/******************************************************************************/

double *equidistant2 ( int n )

/******************************************************************************/
/*
  Purpose:

    EQUIDISTANT2 returns the Type 2 Equidistant points.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    03 March 2018

  Author:

    John Burkardt.

  Parameters:

    Input, int N, the number of points.

    Input, double EQUIDISTANT2[N], the points.
*/
{
  int i;
  double *x;

  x = ( double * ) malloc ( n * sizeof ( double ) );

  if ( n == 1 )
  {
    x[0] = 0.0;
  }
  else
  {
    for ( i = 0; i < n; i++ )
    {
      x[i] = ( double ) ( - n + 1 + 2 * i ) / ( double ) ( n - 1 );
    }
  }

  return x;
}
/******************************************************************************/

double *equidistant3 ( int n )

/******************************************************************************/
/*
  Purpose:

    EQUIDISTANT3 returns the Type 3 Equidistant points.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    03 March 2018

  Author:

    John Burkardt.

  Parameters:

    Input, int N, the number of points.

    Input, double EQUIDISTANT3[N], the points.
*/
{
  int i;
  double *x;

  x = ( double * ) malloc ( n * sizeof ( double ) );

  for ( i = 0; i < n; i++ )
  {
    x[i] = ( double ) ( - n + 1 + 2 * i ) / ( double ) ( n );
  }

  return x;
}
/******************************************************************************/

double *fejer1 ( int n )

/******************************************************************************/
/*
  Purpose:

    FEJER1 returns the Type 1 Fejer points.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    03 March 2018

  Author:

    John Burkardt.

  Parameters:

    Input, int N, the number of points.

    Input, double FEJER1[N], the points.
*/
{
  int i;
  const double r8_pi = 3.141592653589793;
  double theta;
  double *x;

  x = ( double * ) malloc ( n * sizeof ( double ) );

  for ( i = 0; i < n; i++ )
  {
    theta = r8_pi * ( double ) ( 2 * n - 1 - 2 * i ) 
                  / ( double ) ( 2 * n );
    x[i] = cos ( theta );
  }
  return x;
}
/******************************************************************************/

double *fejer2 ( int n )

/******************************************************************************/
/*
  Purpose:

    FEJER2 returns the Type 2 Fejer points.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    03 March 2018

  Author:

    John Burkardt.

  Parameters:

    Input, int N, the number of points.

    Input, double FEJER2[N], the points.
*/
{
  int i;
  const double r8_pi = 3.141592653589793;
  double theta;
  double *x;

  x = ( double * ) malloc ( n * sizeof ( double ) );

  for ( i = 0; i < n; i++ )
  {
    theta = r8_pi * ( double ) ( n - i ) 
                  / ( double ) ( n + 1 );
    x[i] = cos ( theta );
  }

  return x;
}
/******************************************************************************/

double *lagrange_value ( int data_num, double t_data[], int interp_num, 
  double t_interp[] )

/******************************************************************************/
/*
  Purpose:

    LAGRANGE_VALUE evaluates the Lagrange polynomials.

  Discussion:

    Given DATA_NUM distinct abscissas, T_DATA(1:DATA_NUM),
    the I-th Lagrange polynomial L(I)(T) is defined as the polynomial of
    degree DATA_NUM - 1 which is 1 at T_DATA(I) and 0 at the DATA_NUM - 1
    other abscissas.

    A formal representation is:

      L(I)(T) = Product ( 1 <= J <= DATA_NUM, I /= J )
       ( T - T(J) ) / ( T(I) - T(J) )

    This routine accepts a set of INTERP_NUM values at which all the Lagrange
    polynomials should be evaluated.

    Given data values P_DATA at each of the abscissas, the value of the
    Lagrange interpolating polynomial at each of the interpolation points
    is then simple to compute by matrix multiplication:

      P_INTERP(1:INTERP_NUM) =
        P_DATA(1:DATA_NUM) * L_INTERP(1:DATA_NUM,1:INTERP_NUM)

    or, in the case where P is multidimensional:

      P_INTERP(1:M,1:INTERP_NUM) =
        P_DATA(1:M,1:DATA_NUM) * L_INTERP(1:DATA_NUM,1:INTERP_NUM)

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    03 December 2007

  Author:

    John Burkardt

  Parameters:

    Input, int DATA_NUM, the number of data points.
    DATA_NUM must be at least 1.

    Input, double T_DATA[DATA_NUM], the data points.

    Input, int INTERP_NUM, the number of
    interpolation points.

    Input, double T_INTERP[INTERP_NUM], the
    interpolation points.

    Output, double LAGRANGE_VALUE[DATA_NUM*INTERP_NUM], the values
    of the Lagrange polynomials at the interpolation points.
*/
{
  int i;
  int i1;
  int i2;
  int j;
  double *l_interp;

  l_interp = ( double * ) malloc ( data_num * interp_num * sizeof ( double ) );
/*
  Evaluate the polynomial.
*/
  for ( j = 0; j < interp_num; j++ )
  {
    for ( i = 0; i < data_num; i++ )
    {
      l_interp[i+j*data_num] = 1.0;
    }
  }

  for ( i1 = 0; i1 < data_num; i1++ )
  {
    for ( i2 = 0; i2 < data_num; i2++ )
    {
      if ( i1 != i2 )
      {
        for ( j = 0; j < interp_num; j++ )
        {
          l_interp[i1+j*data_num] = l_interp[i1+j*data_num] 
            * ( t_interp[j] - t_data[i2] ) / ( t_data[i1] - t_data[i2] );
        }
      }
    }
  }

  return l_interp;
}
/******************************************************************************/

double lebesgue_constant ( int n, double x[], int nfun, double xfun[] )

/******************************************************************************/
/*
  Purpose:

    LEBESGUE_CONSTANT estimates the Lebesgue constant for a set of points.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    03 March 2014

  Author:

    John Burkardt.

  Parameters:

    Jean-Paul Berrut, Lloyd Trefethen,
    Barycentric Lagrange Interpolation,
    SIAM Review,
    Volume 46, Number 3, September 2004, pages 501-517.

  Parameters:

    Input, int N, the number of interpolation points.

    Input, double X[N], the interpolation points.

    Input, int NFUN, the number of evaluation points.

    Input, double XFUN[NFUN], the evaluation points.

    Output, double LEBESGUE_CONSTANT, an estimate of the Lebesgue constant 
    for the points.
*/
{
  double *lfun;
  double lmax;

  lfun = lebesgue_function ( n, x, nfun, xfun );

  lmax = r8vec_max ( nfun, lfun );

  free ( lfun );

  return lmax;
}
/******************************************************************************/

double *lebesgue_function ( int n, double x[], int nfun, double xfun[] )

/******************************************************************************/
/*
  Purpose:

    LEBESGUE_FUNCTION evaluates the Lebesgue function for a set of points.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    03 March 2014

  Author:

    Original MATLAB version by Greg von Winckel.
    C version by John Burkardt.

  Parameters:

    Jean-Paul Berrut, Lloyd Trefethen,
    Barycentric Lagrange Interpolation,
    SIAM Review,
    Volume 46, Number 3, September 2004, pages 501-517.

  Parameters:

    Input, int N, the number of interpolation points.

    Input, double X[N], the interpolation points.

    Input, int NFUN, the number of evaluation points.

    Input, double XFUN[NFUN], the evaluation points.

    Output, double LEBESGUE_FUNCTION[NFUN], the Lebesgue function values.
*/
{
  int i;
  int j;
  double *lfun;
  double *llfun;
  double t;

  lfun = ( double * ) malloc ( nfun * sizeof ( double ) );
/*
  Handle special case.
*/
  if ( n == 1 )
  {
    for ( j = 0; j < nfun; j++ )
    {
      lfun[j] = 1.0;
    }
    return lfun;
  }

  llfun = lagrange_value ( n, x, nfun, xfun );

  for ( j = 0; j < nfun; j++ )
  {
    t = 0.0;
    for ( i = 0; i < n; i++ )
    {
      t = t + fabs ( llfun[i+j*n] );
    }
    lfun[j] = t;
  }

  free ( llfun );

  return lfun;
}
/******************************************************************************/

void lebesgue_plot ( int n, double x[], int nfun, double xfun[], 
  char *label, char *filename )

/******************************************************************************/
/*
  Purpose:

    LEBESGUE_PLOT plots the Lebesgue function for a set of points.

  Discussion:

    The interpolation interval is assumed to be [min(XFUN), max(XFUN)].

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    03 March 2014

  Author:

    John Burkardt.

  Parameters:

    Jean-Paul Berrut, Lloyd Trefethen,
    Barycentric Lagrange Interpolation,
    SIAM Review,
    Volume 46, Number 3, September 2004, pages 501-517.

  Parameters:

    Input, int N, the number of interpolation points.

    Input, double X[N], the interpolation points.

    Input, int NFUN, the number of evaluation points.

    Input, double XFUN[NFUN], the evaluation points.  

    Input, char *LABEL, a title for the plot.

    Input, char *FILENAME, a partial filename.
    The program will create "filename_commands.txt', 'filename_data.txt',
    and 'filename.png'.
*/
{
  char command_filename[255];
  FILE *command_unit;
  char data_filename[255];
  FILE *data_unit;
  int i;
  double *lfun;
  char png_filename[255];

  lfun = lebesgue_function ( n, x, nfun, xfun );
/*
  Create data file.
*/
  sprintf ( data_filename, "%s_data.txt", filename );
  data_unit = fopen ( data_filename, "wt" );
  for ( i = 0; i < nfun; i++ )
  {
    fprintf ( data_unit, "%g  %g\n", xfun[i], lfun[i] );
  }
  fclose ( data_unit );
  printf ( "\n" );
  printf ( "  Created graphics data file '%s'\n", data_filename );
/*
  Create command file.
*/
  sprintf ( command_filename, "%s_commands.txt", filename );
  command_unit = fopen ( command_filename, "wt" );

  fprintf ( command_unit, "# %s\n", command_filename );
  fprintf ( command_unit, "#\n" );
  fprintf ( command_unit, "# Usage:\n" );
  fprintf ( command_unit, "#  gnuplot < %s\n", command_filename );
  fprintf ( command_unit, "#\n" );
  fprintf ( command_unit, "set term png\n" );

  sprintf ( png_filename, "%s.png", filename );
  fprintf ( command_unit, "set output '%s'\n", png_filename );
  fprintf ( command_unit, "set xlabel '<--- X --->'\n" );
  fprintf ( command_unit, "set ylabel '<--- Lebesgue(X) --->'\n" );
  fprintf ( command_unit, "set title '%s'\n", label );
  fprintf ( command_unit, "set grid\n" );
  fprintf ( command_unit, "set style data lines\n" );
  fprintf ( command_unit, "set timestamp\n" );
  fprintf ( command_unit, "plot '%s' using 1:2 lw 3 linecolor rgb 'red'\n", 
    data_filename );

  fclose ( command_unit );
  printf ( "  Created graphics command file '%s'\n", command_filename );

  free ( lfun );

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
  double *r8vec_pointer;
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
