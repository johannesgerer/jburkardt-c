# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>

# include "square_exactness.h"

/******************************************************************************/

void legendre_2d_exactness ( double a[], double b[], int n, double x[], 
  double y[], double w[], int t )

/******************************************************************************/
/*
  Purpose:

    LEGENDRE_2D_EXACTNESS: monomial exactness for the 2D Legendre integral.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    31 May 2014

  Author:

    John Burkardt

  Parameters:

    Input, double A[2], the lower limits of integration.

    Input, double B[2], the upper limits of integration.

    Input, int N, the number of points in the rule.

    Input, double X[N], Y[N], the quadrature points.

    Input, double W[N], the quadrature weights.

    Input, int T, the maximum total degree.
    0 <= T.
*/
{
  double e;
  int i;
  int j;
  int p[2];
  double q;
  double s;
  int tt;
  double *v;

  printf ( "\n" );
  printf ( "  Quadrature rule for the 2D Legendre integral.\n" );
  printf ( "  Number of points in rule is %d\n", n );
  printf ( "\n" );
  printf ( "   D   I       J          Relative Error\n" );

  v = ( double * ) malloc ( n * sizeof ( double ) );

  for ( tt = 0; tt <= t; tt++ )
  {
    printf ( "  %2d\n", tt );
    for ( j = 0; j <= tt; j++ )
    {
      i = tt - j;

      p[0] = i;
      p[1] = j;

      s = legendre_2d_monomial_integral ( a, b, p );

      for ( i = 0; i < n; i++ )
      {
        v[i] = pow ( x[i], p[0] ) * pow ( y[i], p[1] );
      }
      q = r8vec_dot_product ( n, w, v );

      if ( s == 0.0 )
      {
        e = fabs ( q );
      }
      else
      {
        e = fabs ( q - s ) / fabs ( s );
      }
      printf ( "  %6d  %6d  %24.16f\n", p[0], p[1], e );
    }
  }

  free ( v );

  return;
}
/******************************************************************************/

double legendre_2d_monomial_integral ( double a[], double b[], int p[] )

/******************************************************************************/
/*
  Purpose:

    LEGENDRE_2D_MONOMIAL_INTEGRAL the Legendre integral of a monomial.

  Discussion:

    The Legendre integral to be evaluated has the form

      I(f) = integral ( y1 <= y <= y2 ) 
             integral ( x1 <= x <= x2 ) x^i y^j dx dy

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    31 May 2014

  Author:

    John Burkardt

  Parameters:

    Input, double A[2], the lower limits of integration.

    Input, double B[2], the upper limits of integration.

    Input, int P[2], the exponents of X and Y.

    Output, double LEGENDRE_2D_MONOMIAL_INTEGRAL, the value of the 
    exact integral.
*/
{
  double exact;

  exact = ( pow ( b[0], p[0] + 1 ) - pow ( a[0], p[0] + 1 ) ) 
        / ( double ) ( p[0] + 1 ) 
        * ( pow ( b[1], p[1] + 1 ) - pow ( a[1], p[1] + 1 ) ) 
        / ( double ) ( p[1] + 1 );

  return exact;
}
/******************************************************************************/

double r8vec_dot_product ( int n, double a1[], double a2[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_DOT_PRODUCT computes the dot product of a pair of R8VEC's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 July 2007

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the vectors.

    Input, double A1[N], A2[N], the two vectors to be considered.

    Output, double R8VEC_DOT_PRODUCT, the dot product of the vectors.
*/
{
  int i;
  double value;

  value = 0.0;
  for ( i = 0; i < n; i++ )
  {
    value = value + a1[i] * a2[i];
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

void r8vec2_print ( int n, double a1[], double a2[], char *title )

/******************************************************************************/
/*
  Purpose:

    R8VEC2_PRINT prints an R8VEC2.

  Discussion:

    An R8VEC2 is a dataset consisting of N pairs of real values, stored
    as two separate vectors A1 and A2.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 March 2009

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of components of the vector.

    Input, double A1[N], double A2[N], the vectors to be printed.

    Input, char *TITLE, a title.
*/
{
  int i;

  fprintf ( stdout, "\n" );
  fprintf ( stdout, "%s\n", title );
  fprintf ( stdout, "\n" );
  for ( i = 0; i < n; i++ )
  {
    fprintf ( stdout, "  %4d: %14g  %14g\n", i, a1[i], a2[i] );
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

