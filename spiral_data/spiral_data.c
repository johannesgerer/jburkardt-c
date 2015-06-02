# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <string.h>
# include <time.h>

# include "spiral_data.h"

/******************************************************************************/

void grid_2d ( int x_num, double x_lo, double x_hi, int y_num, double y_lo, 
  double y_hi, double x[], double y[] )

/******************************************************************************/
/*
  Purpose:

    GRID_2D returns a regular 2D grid.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 January 2015

  Author:

    John Burkardt

  Parameters:

    Input, int X_NUM, the number of X values to use.

    Input, double X_LO, X_HI, the range of X values.

    Input, int Y_NUM, the number of Y values to use.

    Input, double Y_LO, Y_HI, the range of Y values.

    Output, double X[X_NUM*Y_NUM], Y[X_NUM*Y_NUM], 
    the coordinates of the grid.
*/
{
  int i;
  int j;
  double xi;
  double yj;

  if ( x_num == 1 )
  {
    for ( j = 0; j < y_num; j++ )
    {
      for ( i = 0; i < x_num; i++ )
      {
        x[i+j*x_num] = ( x_lo + x_hi ) / 2.0;
      }
    }
  }
  else
  {
    for ( i = 0; i < x_num; i++ )
    {
      xi = ( ( double ) ( x_num - i - 1 ) * x_lo   
           + ( double ) (         i     ) * x_hi ) 
           / ( double ) ( x_num     - 1 );
      for ( j = 0; j < y_num; j++ )
      {
        x[i+j*x_num] = xi;
      }
    }
  }

  if ( y_num == 1 )
  {
    for ( j = 0; j < y_num; j++ )
    {
      for ( i = 0; i < x_num; i++ )
      {
        y[i+j*x_num] = ( y_lo + y_hi ) / 2.0;
      }
    }
  }
  else
  {
    for ( j = 0; j < y_num; j++ )
    {
      yj = ( ( double ) ( y_num - j - 1 ) * y_lo   
           + ( double ) (         j     ) * y_hi ) 
           / ( double ) ( y_num     - 1 );
      for ( i = 0; i < x_num; i++ )
      {
        y[i+j*x_num] = yj;
      }
    }
  }

  return;
}
/******************************************************************************/

double r8vec_amax ( int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_AMAX returns the maximum absolute value in an R8VEC.

  Discussion:

    An R8VEC is a vector of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 May 2009

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the array.

    Input, double A[N], the array.

    Output, double AMAX, the value of the entry
    of largest magnitude.
*/
{
  double amax;
  int i;

  amax = 0.0;
  for ( i = 0; i < n; i++ )
  {
    if ( amax < fabs ( a[i] ) )
    {
      amax = fabs ( a[i] );
    }
  }
  return amax;
}
/******************************************************************************/

double r8vec_amin ( int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_AMIN returns the minimum absolute value in an R8VEC.

  Discussion:

    An R8VEC is a vector of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 May 2009

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the array.

    Input, double A[N], the array.

    Output, double R8VEC_AMIN, the value of the entry
    of smallest magnitude.
*/
{
  double amin;
  int i;
  const double r8_huge = 1.79769313486231571E+308;

  amin = r8_huge;
  for ( i = 0; i < n; i++ )
  {
    if ( fabs ( a[i] ) < amin )
    {
      amin = fabs ( a[i] );
    }
  }

  return amin;
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

double *r8vec_uniform_ab_new ( int n, double a, double b, int *seed )

/******************************************************************************/
/*
  Purpose:

    R8VEC_UNIFORM_AB_NEW returns a scaled pseudorandom R8VEC.

  Discussion:

    Each dimension ranges from A to B.

    This routine implements the recursion

      seed = 16807 * seed mod ( 2^31 - 1 )
      unif = seed / ( 2^31 - 1 )

    The integer arithmetic never requires more than 32 bits,
    including a sign bit.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    30 January 2005

  Author:

    John Burkardt

  Reference:

    Paul Bratley, Bennett Fox, Linus Schrage,
    A Guide to Simulation,
    Second Edition,
    Springer, 1987,
    ISBN: 0387964673,
    LC: QA76.9.C65.B73.

    Bennett Fox,
    Algorithm 647:
    Implementation and Relative Efficiency of Quasirandom
    Sequence Generators,
    ACM Transactions on Mathematical Software,
    Volume 12, Number 4, December 1986, pages 362-376.

    Pierre L'Ecuyer,
    Random Number Generation,
    in Handbook of Simulation,
    edited by Jerry Banks,
    Wiley, 1998,
    ISBN: 0471134031,
    LC: T57.62.H37.

    Peter Lewis, Allen Goodman, James Miller,
    A Pseudo-Random Number Generator for the System/360,
    IBM Systems Journal,
    Volume 8, Number 2, 1969, pages 136-143.

  Parameters:

    Input, int N, the number of entries in the vector.

    Input, double A, B, the lower and upper limits of the pseudorandom values.

    Input/output, int *SEED, a seed for the random number generator.

    Output, double R8VEC_UNIFORM_AB_NEW[N], the vector of pseudorandom values.
*/
{
  int i;
  const int i4_huge = 2147483647;
  int k;
  double *r;

  if ( *seed == 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8VEC_UNIFORM_AB_NEW - Fatal error!\n" );
    fprintf ( stderr, "  Input value of SEED = 0.\n" );
    exit ( 1 );
  }

  r = ( double * ) malloc ( n * sizeof ( double ) );

  for ( i = 0; i < n; i++ )
  {
    k = *seed / 127773;

    *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

    if ( *seed < 0 )
    {
      *seed = *seed + i4_huge;
    }

    r[i] = a + ( b - a ) * ( double ) ( *seed ) * 4.656612875E-10;
  }

  return r;
}
/******************************************************************************/

void resid_spiral ( int n, double x[], double y[], double c, double pr[] )

/******************************************************************************/
/*
  Purpose:

    RESID_SPIRAL computes the residual for a spiral velocity vector field.

  Discussion:

    Note that the continuous velocity field (U,V)(X,Y) that is discretely
    sampled here satisfies the homogeneous continuity equation, that is,
    it has zero divergence.  In other words:

      dU/dX + dV/dY = 0.

    This is by construction, since we have

      U(X,Y) =  10 * d/dY ( PHI(X) * PHI(Y) )
      V(X,Y) = -10 * d/dX ( PHI(X) * PHI(Y) )

    which guarantees zero divergence.

    The underlying function PHI is defined by

      PHI(Z) = ( 1 - cos ( C * pi * Z ) ) * ( 1 - Z )^2

    where C is a parameter.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 January 2015

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of evaluation points.

    Input, double X[N], Y[N], the coordinates of the 
    evaluation points.

    Input, double C, a parameter, typically between 0 and 2 * PI.

    Output, double PR[N], the residual in the continuity equation.
*/
{
  int i;
  const double r8_pi = 3.141592653589793;
  double u;
  double ux;
  double v;
  double vy;

  for ( i = 0; i < n; i++ )
  {
    u =   10.0 * ( 1.0 - cos ( c * r8_pi * x[i] ) ) 
           * pow ( 1.0 - x[i], 2 )
           * ( 
               c * r8_pi * sin ( c * r8_pi * y[i] ) * pow ( 1.0 - y[i], 2 )
             - ( 1.0 - cos ( c * r8_pi * y[i] ) ) 
             * 2.0 * ( 1.0 - y[i] ) 
             );

    ux =   10.0 * 
      ( 
        c * r8_pi * sin ( c * r8_pi * x[i] ) * pow ( 1.0 - x[i], 2 )
        - ( 1.0 - cos ( c * r8_pi * x[i] ) ) 
        * 2.0 * ( 1.0 - x[i] ) 
      ) 
      * 
      ( 
        c * r8_pi * sin ( c * r8_pi * y[i] ) * pow ( 1.0 - y[i], 2 )
        - ( 1.0 - cos ( c * r8_pi * y[i] ) ) 
        * 2.0 * ( 1.0 - y[i] ) 
      );

    v = - 10.0 * ( 1.0 - cos ( c * r8_pi * y[i] ) ) 
      * pow ( 1.0 - y[i], 2 )
      * ( 
          c * r8_pi * sin ( c * r8_pi * x[i] ) * pow ( 1.0 - x[i], 2 )
        - ( 1.0 - cos ( c * r8_pi * x[i] ) ) 
        * 2.0 * ( 1.0 - x[i] ) 
        );

    vy =  - 10.0 * 
      ( 
        c * r8_pi * sin ( c * r8_pi * x[i] ) * pow ( 1.0 - x[i], 2 )
        - ( 1.0 - cos ( c * r8_pi * x[i] ) ) 
        * 2.0 * ( 1.0 - x[i] ) 
      ) 
      * 
      ( 
        c * r8_pi * sin ( c * r8_pi * y[i] ) * pow ( 1.0 - y[i], 2 )
        - ( 1.0 - cos ( c * r8_pi * y[i] ) ) 
        * 2.0 * ( 1.0 - y[i] ) 
      );

    pr[i] = ux + vy;
  }

  return;
}
/******************************************************************************/

void spiral_gnuplot ( char *header, int n, double x[], double y[], double u[], 
  double v[], double s )

/******************************************************************************/
/*
  Purpose:

    SPIRAL_GNUPLOT writes the spiral vector field to files for GNUPLOT.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 January 2015

  Author:

    John Burkardt

  Parameters:

    Input, char *HEADER, a header to be used to name the files.

    Input, int N, the number of evaluation points.

    Input, double X[N], Y[N], the coordinates of the evaluation points.

    Input, double U[N], V[N], the velocity components.

    Input, double S, a scale factor for the velocity vectors.
*/
{
  char command_filename[255];
  FILE *command_unit;
  char data_filename[255];
  FILE *data_unit;
  int i;
  char plot_filename[255];
/*
  Write the data file.
*/
  strcpy ( data_filename, header );
  strcat ( data_filename, "_data.txt" );

  data_unit = fopen ( data_filename, "wt" );

  for ( i = 0; i < n; i++ )
  {
    fprintf ( data_unit, "  %g  %g  %g  %g\n", x[i], y[i], s * u[i], s * v[i] );
  }

  fclose ( data_unit );

  printf ( "\n" );
  printf ( "  Data written to '%s'\n", data_filename );
/*
  Write the command file.
*/
  strcpy ( command_filename, header );
  strcat ( command_filename, "_commands.txt" );

  strcpy ( plot_filename, header );
  strcat ( plot_filename, ".png" );

  command_unit = fopen ( command_filename, "wt" );

  fprintf ( command_unit, "#  %s\n", command_filename );
  fprintf ( command_unit, "#\n" );
  fprintf ( command_unit, "set term png\n" );
  fprintf ( command_unit, "set output '%s'\n", plot_filename );
  fprintf ( command_unit, "#\n" );
  fprintf ( command_unit, "#  Add titles and labels.\n" );
  fprintf ( command_unit, "#\n" );
  fprintf ( command_unit, "set xlabel '<--- X --->'\n" );
  fprintf ( command_unit, "set ylabel '<--- Y --->'\n" );
  fprintf ( command_unit, "set title 'Spiral velocity flow'\n" );
  fprintf ( command_unit, "unset key\n" );
  fprintf ( command_unit, "#\n" );
  fprintf ( command_unit, "#  Add grid lines.\n" );
  fprintf ( command_unit, "#\n" );
  fprintf ( command_unit, "set grid\n" );
  fprintf ( command_unit, "set size ratio -1\n" );
  fprintf ( command_unit, "#\n" );
  fprintf ( command_unit, "#  Timestamp the plot.\n" );
  fprintf ( command_unit, "#\n" );
  fprintf ( command_unit, "set timestamp\n" );
  fprintf ( command_unit, 
    "plot '%s' using 1:2:3:4 with vectors \\\n", data_filename );
  fprintf ( command_unit, "  head filled lt 2 linecolor rgb 'blue'\n" );
  fprintf ( command_unit, "quit\n" );

  fclose ( command_unit );

  printf ( "  Commands written to '%s'\n", command_filename );


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
/******************************************************************************/

void uv_spiral ( int n, double x[], double y[], double c, double u[], 
  double v[] )

/******************************************************************************/
/*
  Purpose:

    UV_SPIRAL computes a spiral velocity vector field.

  Discussion:

    Note that the continuous velocity field (U,V)(X,Y) that is discretely
    sampled here satisfies the homogeneous continuity equation, that is,
    it has zero divergence.  In other words:

      dU/dX + dV/dY = 0.

    This is by construction, since we have

      U(X,Y) =  10 * d/dY ( PHI(X) * PHI(Y) )
      V(X,Y) = -10 * d/dX ( PHI(X) * PHI(Y) )

    which guarantees zero divergence.

    The underlying function PHI is defined by

      PHI(Z) = ( 1 - cos ( C * pi * Z ) ) * ( 1 - Z )^2

    where C is a parameter.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    18 January 2015

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of evaluation points.

    Input, double X[N], Y[N], the coordinates of the 
    evaluation points.

    Input, double C, a parameter, typically between 0 and 2 * PI.

    Output, double U[N], V[N], the velocity components.
*/
{
  int i;
  const double r8_pi = 3.141592653589793;

  for ( i = 0; i < n; i++ )
  {
    u[i] =   10.0 * ( 1.0 - cos ( c * r8_pi * x[i] ) ) 
      * pow ( 1.0 - x[i], 2 )
      * ( 
          c * r8_pi * sin ( c * r8_pi * y[i] ) * pow ( 1.0 - y[i], 2 )
        - ( 1.0 - cos ( c * r8_pi * y[i] ) ) 
          * 2.0 * ( 1.0 - y[i] ) 
        );

    v[i] = - 10.0 * ( 1.0 - cos ( c * r8_pi * y[i] ) ) 
      * pow ( 1.0 - y[i], 2 )
      * ( 
          c * r8_pi * sin ( c * r8_pi * x[i] ) * pow ( 1.0 - x[i], 2 )
        - ( 1.0 - cos ( c * r8_pi * x[i] ) ) 
        * 2.0 * ( 1.0 - x[i] ) 
        );

  }

  return;
}
