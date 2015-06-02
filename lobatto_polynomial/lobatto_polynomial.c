# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <string.h>

# include "lobatto_polynomial.h"

/******************************************************************************/

double *lobatto_polynomial_derivative ( int m, int n, double x[] )

/******************************************************************************/
/*
  Purpose:

    LOBATTO_POLYNOMIAL_DERIVATIVE: derivative of completed Lobatto polynomial.

  Discussion:

    L(N,X)  =  N * ( P(N-1,X) - X * P(N,X) ) 
    L'(N,X) =  N * ( P'(N-1,X) - P(N,X) - X * P'(N,X) )

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    20 November 2014

  Author:

    John Burkardt

  Reference:

    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.

    Larry Andrews,
    Special Functions of Mathematics for Engineers,
    Second Edition,
    Oxford University Press, 1998,
    ISBN: 0819426164,
    LC: QA351.A75.

    Daniel Zwillinger, editor,
    CRC Standard Mathematical Tables and Formulae,
    30th Edition,
    CRC Press, 1996.

  Parameters:

    Input, int M, the number of evaluation points.

    Input, int N, the highest order polynomial to evaluate.
    Note that polynomials 0 through N will be evaluated.

    Input, double X[M], the evaluation points.

    Output, double LOBATTO_POLYNOMIAL_DERIVATIVE[M*N], the derivative of 
    the completed Lobatto polynomials of order 1 through N at the point X.
*/
{
  int i;
  int j;
  double *lp;
  double *p;
  double *pp;

  lp = ( double * ) malloc ( m * n * sizeof ( double ) );
  p = ( double * ) malloc ( m * ( n + 1 ) * sizeof ( double ) );
  pp = ( double * ) malloc ( m * ( n + 1 ) * sizeof ( double ) );

  for ( i = 0; i < m; i++ )
  {
    if ( 1 <= n )
    {
      lp[i+0*m] = - 2.0 * x[i];

      if ( 2 <= n )
      {
        p[i+0*m] = 1.0;
        p[i+1*m] = x[i];
        for ( j = 1; j < n; j++ )
        {
          p[i+(j+1)*m] = 
            ( ( double ) ( 2 * j + 1 ) * x[i] * p[i+j*m]     
            - ( double ) (     j     ) *        p[i+(j-1)*m] ) 
            / ( double ) (     j + 1 );
        }

        pp[i+0*m] = 0.0;
        pp[i+1*m] = 1.0;
        for ( j = 1; j < n; j++ )
        {
          pp[i+(j+1)*m] = 
            ( ( double ) ( 2 * j + 1 ) * ( p[i+j*m] + x[i] * pp[i+j*m] )   
            - ( double ) (     j     ) *                     pp[i+(j-1)*m] ) 
            / ( double ) (     j + 1 );
        }

        for ( j = 1; j < n; j++ )
        {
          lp[i+j*m] = 
            ( double ) ( j + 1 ) 
            * ( pp[i+j*m] - p[i+(j+1)*m] - x[i] * pp[i+(j+1)*m] );
        }
      }
    }
  }

  free ( p );
  free ( pp );

  return lp;
}
/******************************************************************************/

void lobatto_polynomial_derivatives ( int *n_data, int *n, double *x, 
  double *fx )

/******************************************************************************/
/*
  Purpose:

    LOBATTO_POLYNOMIAL_DERIVATIVES: derivatives of completed Lobatto polynomials.

  Discussion:

    In Mathematica, the function can be evaluated by:

      n * LegendreP [ n - 1, x ] - n * x * LegendreP [ n, x ]

     In Mathematica, the completed Lobatto polynomial can be evaluated by:
 
       n * LegendreP [ n - 1, x ] - n * x * LegendreP [ n, x ]
 
     The derivative is:
 
         n * D[LegendreP [ n - 1, x ], {x} ] 
       - n * LegendreP [ n, x ] 
       - n * x * D[LegendreP [ n, x ], {x}]

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    18 November 2014

  Author:

    John Burkardt

  Parameters:

    Input/output, int *N_DATA.  The user sets N_DATA to 0
    before the first call.  On each call, the routine increments N_DATA by 1,
    and returns the corresponding data; when there is no more data, the
    output value of N_DATA will be 0 again.

    Output, int *N, the order of the function.

    Output, double *X, the point where the function is evaluated.

    Output, double *FX, the value of the function.
*/
{
# define N_MAX 31

  static double fx_vec[N_MAX] = {
     -0.5, 
      2.437500000000000, 
      4.031250000000000, 
     -3.154296875000000, 
    -10.19165039062500, 
     -1.019622802734375, 
     15.67544555664063, 
     10.97668933868408, 
    -15.91419786214828, 
    -24.33202382177114, 
     12.00000000000000, 
      5.670000000000000, 
      0.9600000000000000, 
     -2.310000000000000, 
     -4.320000000000000, 
     -5.250000000000000, 
     -5.280000000000000, 
     -4.590000000000000, 
     -3.360000000000000, 
     -1.770000000000000, 
      0.0, 
      1.770000000000000, 
      3.360000000000000, 
      4.590000000000000, 
      5.280000000000000, 
      5.250000000000000, 
      4.320000000000000, 
      2.310000000000000, 
     -0.9600000000000000, 
     -5.670000000000000, 
    -12.00000000000000 };

  static int n_vec[N_MAX] = {
     1,  2, 
     3,  4,  5, 
     6,  7,  8, 
     9, 10,  3, 
     3,  3,  3, 
     3,  3,  3, 
     3,  3,  3, 
     3,  3,  3, 
     3,  3,  3, 
     3,  3,  3, 
     3,  3 };

  static double x_vec[N_MAX] = {
    0.25, 
    0.25, 
    0.25, 
    0.25, 
    0.25, 
    0.25, 
    0.25, 
    0.25, 
    0.25, 
    0.25, 
   -1.00, 
   -0.90, 
   -0.80, 
   -0.70, 
   -0.60, 
   -0.50, 
   -0.40, 
   -0.30, 
   -0.20, 
   -0.10, 
    0.00, 
    0.10, 
    0.20, 
    0.30, 
    0.40, 
    0.50, 
    0.60, 
    0.70, 
    0.80, 
    0.90, 
    1.00 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  *n_data = *n_data + 1;

  if ( N_MAX < *n_data )
  {
    *n_data = 0;
    *n = 0;
    *x = 0.0;
    *fx = 0.0;
  }
  else
  {
    *n  = n_vec[*n_data-1];
    *x  = x_vec[*n_data-1];
    *fx = fx_vec[*n_data-1];
  }

  return;
# undef N_MAX
}
/******************************************************************************/

void lobatto_polynomial_plot ( int ndx_num, int ndx[], char *prefix )

/******************************************************************************/
/*
  Purpose:

    LOBATTO_POLYNOMIAL_PLOT plots one or more completed Lobatto polynomials.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 November 2014

  Author:

    John Burkardt

  Parameters:

    Input, int NDX_NUM, the number of polynomials to plot.

    Input, int NDX[NDX_NUM], the orders of 1 or more 
    Legendre polynomials to be plotted together.

    Input, char *PREFIX. the filename prefix.
*/
{
  char command_filename[255];
  FILE *command_unit;
  char data_filename[255];
  FILE *data_unit;
  int i;
  int j;
  double *l;
  double *lp;
  int n;
  char plot_filename[255];
  double *x;
  double x_hi;
  double x_lo;
  static int x_num = 501;
  double *y;
  double *yp;

  x_lo = -1.0;
  x_hi = +1.0;
  x = r8vec_linspace_new ( x_num, x_lo, x_hi );
/*
  Collect the data.
*/
  y = ( double * ) malloc ( x_num * ndx_num * sizeof ( double ) );
  yp = ( double * ) malloc ( x_num * ndx_num * sizeof ( double ) );

  for ( j = 0; j < ndx_num; j++ )
  {
    n = ndx[j];

    l = lobatto_polynomial_value ( x_num, n, x );
    for ( i = 0; i < x_num; i++ )
    {
      y[i+j*x_num] = l[i+(n-1)*x_num];
    }
    free ( l );

    lp = lobatto_polynomial_derivative ( x_num, n, x );
    for ( i = 0; i < x_num; i++ )
    {
      yp[i+j*x_num] = lp[i+(n-1)*x_num];
    }
    free ( lp );
  }

  printf ( "\n" );
/*
  Make data file for values.
*/
  strcpy ( data_filename, prefix );
  strcat ( data_filename, "_value_data.txt" );

  data_unit = fopen ( data_filename, "wt" );

  for ( i = 0; i < x_num; i++ )
  {
    fprintf ( data_unit, "%g", x[i] );
    for ( j = 0; j < ndx_num; j++ )
    {
      fprintf ( data_unit, "  %g", y[i+j*x_num] );
    }
    fprintf ( data_unit, "\n" );
  }
  fclose ( data_unit );

  printf ( "  Lobatto value data in '%s'\n", data_filename );
/*
  Make command file for values.
*/
  strcpy ( command_filename, prefix );
  strcat ( command_filename, "_value_commands.txt" );

  command_unit = fopen ( command_filename, "wt" );

  fprintf ( command_unit, "# %s\n", command_filename );
  fprintf ( command_unit, "#\n" );
  fprintf ( command_unit, "# Usage:\n" );
  fprintf ( command_unit, "#  gnuplot < %s\n", command_filename );
  fprintf ( command_unit, "#\n" );
  fprintf ( command_unit, "unset key\n" );
  fprintf ( command_unit, "set term png\n" );
  fprintf ( command_unit, "set timestamp\n" );

  strcpy ( plot_filename, prefix );
  strcat ( plot_filename, "_value.png" );

  fprintf ( command_unit, "set output '%s'\n", plot_filename );
  fprintf ( command_unit, "set xlabel 'x'\n" );
  fprintf ( command_unit, "set ylabel 'L(n,x)'\n" );
  fprintf ( command_unit, "set title 'Lobatto values'\n" );
  fprintf ( command_unit, "set grid\n" );
  fprintf ( command_unit, "set style data lines\n" );
  for ( j = 0; j < ndx_num; j++ )
  {
    if ( j == 0 )
    {
      fprintf ( command_unit, "plot '" );
    }
    else
    {
      fprintf ( command_unit, "     '" );
    }
    fprintf ( command_unit, "%s' using 1:%d", data_filename, j + 2 );
    if ( j < ndx_num - 1 )
    {
      fprintf ( command_unit, ", \\" );
    }
    fprintf ( command_unit, "\n" );
  }
  fprintf ( command_unit, "quit\n" );

  fclose ( command_unit );

  printf ( "  Lobatto value commands in '%s'\n", command_filename );
/*
  Make data file for derivatives.
*/
  strcpy ( data_filename, prefix );
  strcat ( data_filename, "_derivative_data.txt" );

  data_unit = fopen ( data_filename, "wt" );

  for ( i = 0; i < x_num; i++ )
  {
    fprintf ( data_unit, "%g", x[i] );
    for ( j = 0; j < ndx_num; j++ )
    {
      fprintf ( data_unit, "  %g", yp[i+j*x_num] );
    }
    fprintf ( data_unit, "\n" );
  }
  fclose ( data_unit );

  printf ( "  Lobatto derivative data stored in '%s'\n", data_filename );
/*
  Make command file for derivatives.
*/
  strcpy ( command_filename, prefix );
  strcat ( command_filename, "_derivative_commands.txt" );

  command_unit = fopen ( command_filename, "wt" );

  fprintf ( command_unit, "# %s\n", command_filename );
  fprintf ( command_unit, "#\n" );
  fprintf ( command_unit, "# Usage:\n" );
  fprintf ( command_unit, "#  gnuplot < %s\n", command_filename );
  fprintf ( command_unit, "#\n" );
  fprintf ( command_unit, "unset key\n" );
  fprintf ( command_unit, "set term png\n" );
  fprintf ( command_unit, "set timestamp\n" );

  strcpy ( plot_filename, prefix );
  strcat ( plot_filename, "_derivative.png" );

  fprintf ( command_unit, "set output '%s'\n", plot_filename );
  fprintf ( command_unit, "set xlabel 'x'\n" );
  fprintf ( command_unit, "set ylabel 'L(n,x)'\n" );
  fprintf ( command_unit, "set title 'Lobatto derivatives'\n" );
  fprintf ( command_unit, "set grid\n" );
  fprintf ( command_unit, "set style data lines\n" );
  for ( j = 0; j < ndx_num; j++ )
  {
    if ( j == 0 )
    {
      fprintf ( command_unit, "plot '" );
    }
    else
    {
      fprintf ( command_unit, "     '" );
    }
    fprintf ( command_unit, "%s' using 1:%d", data_filename, j + 2 );
    if ( j < ndx_num - 1 )
    {
      fprintf ( command_unit, ", \\" );
    }
    fprintf ( command_unit, "\n" );
  }
  fprintf ( command_unit, "quit\n" );

  fclose ( command_unit );

  printf ( "  Lobatto derivative commands in '%s'\n", command_filename );

  free ( x );
  free ( y );
  free ( yp );

  return;
}
/******************************************************************************/

double *lobatto_polynomial_value ( int m, int n, double x[] )

/******************************************************************************/
/*
  Purpose:

    LOBATTO_POLYNOMIAL_VALUE evaluates the completed Lobatto polynomials Lo(n,x).

  Discussion:

    L(N,X) = ( 1 - X^2 ) * P'(N,X)
           = N * ( P(N-1,X) - X * P(N,X) ) 

    The Lobatto polynomials are 0 at -1 and +1.

      (1-x^2) * 1
      (1-x^2) * 3X
      (1-x^2) * ( -3 + 15x^2 ) / 2
      (1-x^2) * ( -60x + 140x^3 ) / 8
      (1-x^2) * ( -15 - 210x^2 + 315x^4 ) / 8
      (1-x^2) * ( 210x - 1260x^3 + 1386x^5 ) / 16
      (1-x^2) * ( -35 + 945x^2 - 3465x^4 + 3003x^6 ) / 16
      (1-x^2) * ( -2520x + 27720x^3 - 72072x^5 + 51480x^7 ) / 128
      (1-x^2) * ( 315 - 13860x^2 + 90090x^4 - 180180x^6 + 109395x^8 ) / 128
      (1-x^2) * ( 6930x - 120120x^3 + 540540x^5 - 875160x^7 + 461890x^9 ) / 256

    Mathematica: (replacing "n" by desired index):

      Expand [ ( 1-x^2) * D [ LegendreP[n,x], {x} ] ]

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    20 November 2014

  Author:

    John Burkardt

  Reference:

    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.

    Larry Andrews,
    Special Functions of Mathematics for Engineers,
    Second Edition,
    Oxford University Press, 1998,
    ISBN: 0819426164,
    LC: QA351.A75.

    Daniel Zwillinger, editor,
    CRC Standard Mathematical Tables and Formulae,
    30th Edition,
    CRC Press, 1996.

  Parameters:

    Input, int M, the number of evaluation points.

    Input, int N, the highest order polynomial to evaluate.
    Note that polynomials 0 through N will be evaluated.

    Input, double X[M], the evaluation points.

    Output, double LOBATTO_POLYNOMIAL_VALUE[M*N], the values of the 
    completed Lobatto polynomials of order 1 through N at the point X.
*/
{
  int i;
  int j;
  double *l;
  double *p;

  l = ( double * ) malloc ( m * n * sizeof ( double ) );
  p = ( double * ) malloc ( m * ( n + 1 ) * sizeof ( double ) );

  for ( i = 0; i < m; i++ )
  {
    if ( 1 <= n )
    {
      l[i+0*m] = 1.0 - x[i] * x[i];

      if ( 2 <= n )
      {
        p[i+0*m] = 1.0;
        p[i+1*m] = x[i];

        for ( j = 1; j < n; j++ )
        {
          p[i+(j+1)*m] = 
            ( ( double ) ( 2 * j + 1 ) * x[i] * p[i+j*m]     
            - ( double ) (     j ) *            p[i+(j-1)*m] ) 
            / ( double ) (     j + 1 );
        }

        for ( j = 1; j < n; j++ )
        {
          l[i+j*m] = ( double ) ( j + 1 ) * ( p[i+j*m] - x[i] * p[i+(j+1)*m] );
        }
      }
    }
  }

  free ( p );

  return l;
}
/******************************************************************************/

void lobatto_polynomial_values ( int *n_data, int *n, double *x, double *fx )

/******************************************************************************/
/*
  Purpose:

    LOBATTO_POLYNOMIAL_VALUES: values of the completed Lobatto polynomials.

  Discussion:

    In Mathematica, the function can be evaluated by:

      n * LegendreP [ n - 1, x ] - n * x * LegendreP [ n, x ]

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 May 2013

  Author:

    John Burkardt

  Parameters:

    Input/output, int *N_DATA.  The user sets N_DATA to 0
    before the first call.  On each call, the routine increments N_DATA by 1,
    and returns the corresponding data; when there is no more data, the
    output value of N_DATA will be 0 again.

    Output, int *N, the order of the function.

    Output, double *X, the point where the function is evaluated.

    Output, double *FX, the value of the function.
*/
{
# define N_MAX 31

  static double fx_vec[N_MAX] = {
    0.9375000000000000, 
    0.7031250000000000, 
   -0.9667968750000000, 
   -1.501464843750000, 
    0.3639221191406250, 
    2.001914978027344, 
    0.6597948074340820, 
   -1.934441328048706, 
   -1.769941113889217, 
    1.215243665501475, 
    0.000000000000000, 
    0.8692500000000000, 
    1.188000000000000, 
    1.109250000000000, 
    0.7680000000000000, 
    0.2812500000000000, 
   -0.2520000000000000, 
   -0.7507500000000000, 
   -1.152000000000000, 
   -1.410750000000000, 
   -1.500000000000000, 
   -1.410750000000000, 
   -1.152000000000000, 
   -0.7507500000000000, 
   -0.2520000000000000, 
    0.2812500000000000, 
    0.7680000000000000, 
    1.109250000000000, 
    1.188000000000000, 
    0.8692500000000000, 
    0.000000000000000 };

  static int n_vec[N_MAX] = {
     1,  2, 
     3,  4,  5, 
     6,  7,  8, 
     9, 10,  3, 
     3,  3,  3, 
     3,  3,  3, 
     3,  3,  3, 
     3,  3,  3, 
     3,  3,  3, 
     3,  3,  3, 
     3,  3 };

  static double x_vec[N_MAX] = {
    0.25, 
    0.25, 
    0.25, 
    0.25, 
    0.25, 
    0.25, 
    0.25, 
    0.25, 
    0.25, 
    0.25, 
   -1.00, 
   -0.90, 
   -0.80, 
   -0.70, 
   -0.60, 
   -0.50, 
   -0.40, 
   -0.30, 
   -0.20, 
   -0.10, 
    0.00, 
    0.10, 
    0.20, 
    0.30, 
    0.40, 
    0.50, 
    0.60, 
    0.70, 
    0.80, 
    0.90, 
    1.00 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  *n_data = *n_data + 1;

  if ( N_MAX < *n_data )
  {
    *n_data = 0;
    *n = 0;
    *x = 0.0;
    *fx = 0.0;
  }
  else
  {
    *n  = n_vec[*n_data-1];
    *x  = x_vec[*n_data-1];
    *fx = fx_vec[*n_data-1];
  }

  return;
# undef N_MAX
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

