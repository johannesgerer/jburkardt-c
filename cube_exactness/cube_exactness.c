# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>

# include "cube_exactness.h"

/******************************************************************************/

void legendre_3d_exactness ( double a[], double b[], int n, double x[], 
  double y[], double z[], double w[], int t )

/******************************************************************************/
/*
  Purpose:

    LEGENDRE_3D_EXACTNESS: monomial exactness for the 3D Legendre integral.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 August 2014

  Author:

    John Burkardt

  Parameters:

    Input, double A[3], the lower limits of integration.

    Input, double B[3], the upper limits of integration.

    Input, int N, the number of points in the rule.

    Input, double X[N], Y[N], Z[N], the quadrature points.

    Input, double W[N], the quadrature weights.

    Input, int T, the maximum total degree.
    0 <= T.
*/
{
  double e;
  int i;
  int j;
  int k;
  int l;
  int p[3];
  double q;
  double s;
  int tt;
  double *v;

  v = ( double * ) malloc ( n * sizeof ( double ) );

  printf ( "\n" );
  printf ( "  Quadrature rule for the 3D Legendre integral.\n" );
  printf ( "  Number of points in rule is %d\n", n );
  printf ( "\n" );
  printf ( "   D   I       J       K          Relative Error\n" );

  for ( tt = 0; tt <= t; tt++ )
  {
    printf ( "  %2d\n", tt );

    for ( k = 0; k <= tt; k++ )
    {
      for ( j = 0; j <= tt - k; j++ )
      {
        i = tt - j - k;

        p[0] = i;
        p[1] = j;
        p[2] = k;

        s = legendre_3d_monomial_integral ( a, b, p );

        for ( l = 0; l < n; l++ )
        {
          v[l] = pow ( x[l], p[0] ) * pow ( y[l], p[1] ) * pow ( z[l], p[2] );
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
        printf ( "  %6d  %6d  %6d  %24.16f\n", p[0], p[1], p[2], e );
      }
    }
  }

  free ( v );

  return;
}
/******************************************************************************/

double legendre_3d_monomial_integral ( double a[], double b[], int p[] )

/******************************************************************************/
/*
  Purpose:

    LEGENDRE_3D_MONOMIAL_INTEGRAL the Legendre integral of a monomial.

  Discussion:

    The Legendre integral to be evaluated has the form

      I(f) = integral ( z1 <= z <= z2 )
             integral ( y1 <= y <= y2 ) 
             integral ( x1 <= x <= x2 ) x^i y^j z^k dx dy dz

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 August 2014

  Author:

    John Burkardt

  Parameters:

    Input, double A[3], the lower limits of integration.

    Input, double B[3], the upper limits of integration.

    Input, int P[3], the exponents of X and Y.

    Output, double LEGENDRE_3D_MONOMIAL_INTEGRAL, the value of the exact integral.
*/
{
  double value;

  value = ( pow ( b[0], p[0] + 1 ) - pow ( a[0], p[0] + 1 ) ) / ( double ) ( p[0] + 1 ) 
        * ( pow ( b[1], p[1] + 1 ) - pow ( a[1], p[1] + 1 ) ) / ( double ) ( p[1] + 1 ) 
        * ( pow ( b[2], p[2] + 1 ) - pow ( a[2], p[2] + 1 ) ) / ( double ) ( p[2] + 1 );

  return value;
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

