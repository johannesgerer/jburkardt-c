# include <stdlib.h>
# include <stdio.h>
# include <time.h>

# include "compass_search.h"

/******************************************************************************/

double *compass_search ( double function_handle ( int m, double x[] ), int m, 
  double x0[], double delta_tol, double delta_init, int k_max, double *fx, 
  int *k )

/******************************************************************************/
/*
  Purpose:

    COMPASS_SEARCH carries out a direct search minimization algorithm.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 January 2012

  Author:

    John Burkardt

  Reference:

    Tamara Kolda, Robert Michael Lewis, Virginia Torczon,
    Optimization by Direct Search: New Perspectives on Some Classical 
    and Modern Methods,
    SIAM Review,
    Volume 45, Number 3, 2003, pages 385-482. 

  Parameters:

    Input, double FUNCTION_HANDLE ( int m, double x[] ), the name of
    a function which evaluates the function to be minimized.

    Input, int M, the number of variables.

    Input, double X0[M], a starting estimate for the minimizer.

    Input, double DELTA_TOL, the smallest step size that is allowed.

    Input, double DELTA_INIT, the starting stepsize.  

    Input, int K_MAX, the maximum number of steps allowed.

    Output, double COMPASS_SEARCH[M], the estimated minimizer.

    Output, double *FX, the function value at X.

    Output, int *K, the number of steps taken.
*/
{
  int decrease;
  double delta;
  double fxd;
  int i;
  int ii;
  double s;
  double *x;
  double *xd;

  *k = 0;
  x = ( double * ) malloc ( m * sizeof ( double ) );
  xd = ( double * ) malloc ( m * sizeof ( double ) );
  r8vec_copy ( m, x0, x );
  *fx = function_handle ( m, x );

  if ( delta_tol <= 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "COMPASS_SEARCH - Fatal error!\n" );
    fprintf ( stderr, "  DELTA_TOL <= 0.0.\n" );
    fprintf ( stderr, "  DELTA_TOL = %g\n", delta_tol );
    exit ( 1 );
  }

  if ( delta_init <= delta_tol )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "COMPASS_SEARCH - Fatal error!\n" );
    fprintf ( stderr, "  DELTA_INIT < DELTA_TOL.\n" );
    fprintf ( stderr, "  DELTA_INIT = %g\n", delta_init );
    fprintf ( stderr, "  DELTA_TOL = %g\n", delta_tol );
    exit ( 1 );
  }

  delta = delta_init;

  while ( *k < k_max )
  {
    *k = *k + 1;
/*
  For each coordinate direction I, seek a lower function value
  by increasing or decreasing X(I) by DELTA.
*/
    decrease = 0;
    s = + 1.0;
    i = 0;

    for ( ii = 1; ii <= 2 * m; ii++ )
    {
      r8vec_copy ( m, x, xd );
      xd[i] = xd[i] + s * delta;
      fxd = function_handle ( m, xd );
/*
  As soon as a decrease is noticed, accept the new point.
*/
      if ( fxd < *fx )
      {
        r8vec_copy ( m, xd, x );
        *fx = fxd;
        decrease = 1;
        break;
      }

      s = - s;
      if ( s == + 1.0 )
      {
        i = i + 1;
      }
    }
/*
  If no decrease occurred, reduce DELTA.
*/
    if ( !decrease )
    {
      delta = delta / 2.0;
      if ( delta < delta_tol )
      {
        break;
      }
    }
  }

  free ( xd );

  return x;
}
/******************************************************************************/

double r8_abs ( double x )

/******************************************************************************/
/*
  Purpose:

    R8_ABS returns the absolute value of an R8.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 May 2006

  Author:

    John Burkardt

  Parameters:

    Input, double X, the quantity whose absolute value is desired.

    Output, double R8_ABS, the absolute value of X.
*/
{
  double value;

  if ( 0.0 <= x )
  {
    value = + x;
  }
  else
  {
    value = - x;
  }
  return value;
}
/******************************************************************************/

void r8vec_copy ( int n, double a1[], double a2[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_COPY copies an R8VEC.

  Discussion:

    An R8VEC is a vector of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    03 July 2005

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the vectors.

    Input, double A1[N], the vector to be copied.

    Input, double A2[N], the copy of A1.
*/
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    a2[i] = a1[i];
  }
  return;
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
    fprintf ( stdout, "  %8d: %14f\n", i, a[i] );
  }

  return;
}
/******************************************************************************/

void timestamp ( void )

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
