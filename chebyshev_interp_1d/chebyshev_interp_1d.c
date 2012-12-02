# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>

# include "chebyshev_interp_1d.h"
# include "qr_solve.h"
# include "r8lib.h"

/******************************************************************************/

double *chebyshev_coef_1d ( int nd, double xd[], double yd[], double *xmin, 
  double *xmax )

/******************************************************************************/
/*
  Purpose:

    CHEBYSHEV_COEF_1D determines the Chebyshev interpolant coefficients.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    27 September 2012

  Author:

    John Burkardt

  Parameters:

    Input, int ND, the number of data points.
    ND must be at least 1.

    Input, double XD[ND], the data locations.

    Input, double YD[ND], the data values.

    Output, double *XMIN, *XMAX, the interpolation interval.

    Output, double CHEBYSHEV_COEF_1D[ND], the Chebyshev coefficients.
*/
{
  double *a;
  double *c;
  int i;
  int j;
  double *x;

  if ( nd == 1 )
  {
    *xmin = xd[0];
    *xmax = xd[0];
    c = ( double * ) malloc ( nd * sizeof ( double ) );
    c[0] = 1.0;
    return c;
  }

  *xmin = r8vec_min ( nd, xd );
  *xmax = r8vec_max ( nd, xd );
/*
  Map XD to [-1,+1].
*/
  x = ( double * ) malloc ( nd * sizeof ( double ) );
  for ( i = 0; i < nd; i++ )
  {
    x[i] = ( 2.0 * xd[i] - *xmin - *xmax ) / ( *xmax - *xmin );
  }
/*
  Form the Chebyshev Vandermonde matrix.
*/
  a = ( double * ) malloc ( nd * nd * sizeof ( double ) );
  for ( j = 0; j < nd; j++ )
  {
    for ( i = 0; i < nd; i++ )
    {
      a[i+j*nd] = cos ( acos ( x[i] ) * ( double ) ( j ) );
    }
  }
/*
  Solve for the expansion coefficients.
*/
  c = qr_solve ( nd, nd, a, yd );

  free ( a );
  free ( x );

  return c;
}
/******************************************************************************/

double *chebyshev_interp_1d ( int nd, double xd[], double yd[], int ni, 
  double xi[] )

/******************************************************************************/
/*
  Purpose:

    CHEBYSHEV_INTERP_1D determines and evaluates the Chebyshev interpolant.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    27 September 2012

  Author:

    John Burkardt

  Parameters:

    Input, int ND, the number of data points.
    ND must be at least 1.

    Input, double XD[ND], the data locations.

    Input, double YD[ND], the data values.

    Input, int NI, the number of interpolation points.

    Input, double XI[NI], the interpolation points, which
    must be each be in the interval [ min(XD), max(XD)].

    Output, double YI[NI], the interpolated values.
*/
{
  double *c;
  double xmax;
  double xmin;
  double *yi;

  c = chebyshev_coef_1d ( nd, xd, yd, &xmin, &xmax );

  yi = chebyshev_value_1d ( nd, c, xmin, xmax, ni, xi );

  free ( c );

  return yi;
}
/******************************************************************************/

double *chebyshev_value_1d ( int nd, double c[], double xmin, double xmax, 
  int ni, double xi[] )

/******************************************************************************/
/*
  Purpose:

    CHEBYSHEV_VALUE_1D evaluates a Chebyshev interpolant, given its coefficients.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    27 September 2012

  Author:

    John Burkardt

  Parameters:

    Input, int ND, the number of data points.
    ND must be at least 1.

    Input, double C[ND], the Chebyshev coefficients.

    Input, double XMIN, XMAX, the interpolation interval.

    Input, int NI, the number of interpolation points.

    Input, double XI[NI], the interpolation points, which
    must be each be in the interval [XMIN,XMAX].

    Output, double YI[NI], the interpolated values.
*/
{
  double *a;
  int i;
  int j;
  double *x;
  double *yi;

  if ( nd == 1 )
  {
    yi = ( double * ) malloc ( nd * sizeof ( double ) );
    yi[0] = c[0];
    return yi;
  }
/*
  Map XI to [-1,+1].
*/
  x = ( double * ) malloc ( ni * sizeof ( double ) );
  for ( i = 0; i < ni; i++ )
  {
    x[i] = ( 2.0 * xi[i] - xmin - xmax ) / ( xmax - xmin );
  }

  a = ( double * ) malloc ( ni * nd * sizeof ( double ) );
  for ( j = 0; j < nd; j++ )
  {
    for ( i = 0; i < ni; i++ )
    {
      a[i+j*ni] = cos ( acos ( x[i] ) * ( double ) ( j ) );
    }
  }

  yi = r8mat_mv_new ( ni, nd, a, c );

  free ( a );
  free ( x );

  return yi;
}
