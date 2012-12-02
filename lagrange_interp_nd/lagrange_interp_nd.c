# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <string.h>

# include "lagrange_interp_nd.h"
# include "r8lib.h"

/******************************************************************************/

double *cc_compute_points ( int n )

/******************************************************************************/
/*
  Purpose:

    CC_COMPUTE_POINTS computes Clenshaw Curtis quadrature points.

  Discussion:

    Our convention is that the abscissas are numbered from left to right.

    This rule is defined on [-1,1].

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    29 September 2012

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of the rule.

    Output, double CC_COMPUTE_POINTS[N], the abscissas.
*/
{
  int i;
  double pi = 3.141592653589793;
  double *x;

  if ( n < 1 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "CC_COMPUTE_POINTS - Fatal error!\n" );
    fprintf ( stderr, "  N < 1.\n" );
    exit ( 1 );
  }

  x = ( double * ) malloc ( n * sizeof ( double ) );

  if ( n == 1 )
  {
    x[0] = 0.0;
  }
  else
  {
    for ( i = 1; i <= n; i++ )
    {
      x[i-1] =  cos ( ( double ) ( n - i ) * pi 
                    / ( double ) ( n - 1     ) );
    }
    x[0] = -1.0;
    if ( ( n % 2 ) == 1 )
    {
      x[(n-1)/2] = 0.0;
    }
    x[n-1] = +1.0;
  }
  return x;
}
/******************************************************************************/

int i4vec_product ( int n, int a[] )

/******************************************************************************/
/*
  Purpose:

    I4VEC_PRODUCT multiplies the entries of an I4VEC.

  Discussion:

    An I4VEC is a vector of I4's.

  Example:

    Input:

      A = ( 1, 2, 3, 4 )

    Output:

      I4VEC_PRODUCT = 24

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 August 2008

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the vector.

    Input, int A[N], the vector

    Output, int I4VEC_PRODUCT, the product of the entries of A.
*/
{
  int i;
  int product;

  product = 1;
  for ( i = 0; i < n; i++ )
  {
    product = product * a[i];
  }

  return product;
}
/******************************************************************************/

double *lagrange_basis_1d ( int nd, double xd[], int ni, double xi[] )

/******************************************************************************/
/*
  Purpose:

    LAGRANGE_BASIS_1D evaluates the Lagrange basis polynomials.

  Discussion:

    Given ND distinct abscissas, XD(1:ND),
    the I-th Lagrange basis polynomial LB(I)(T) is defined as the polynomial of
    degree ND - 1 which is 1 at XD(I) and 0 at the ND - 1
    other abscissas.

    A formal representation is:

      LB(I)(T) = Product ( 1 <= J <= ND, I /= J )
       ( T - T(J) ) / ( T(I) - T(J) )

    This routine accepts a set of NI values at which all the Lagrange
    basis polynomials should be evaluated.

    Given data values YD at each of the abscissas, the value of the
    Lagrange interpolating polynomial at each of the interpolation points
    is then simple to compute by matrix multiplication:

      YI(1:NI) = LB(1:NI,1:ND) * YD(1:ND)

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 September 2012

  Author:

    John Burkardt

  Parameters:

    Input, int ND, the number of data points.
    ND must be at least 1.

    Input, double XD[ND], the data points.

    Input, int NI, the number of interpolation points.

    Input, double XI[NI], the interpolation points.

    Output, double LAGRANGE_BASIS_1D[NI*ND], the values
    of the Lagrange basis polynomials at the interpolation points.
*/
{
  int i;
  int j;
  int k;
  double *lb;
/*
  Evaluate the polynomial.
*/
  lb = ( double * ) malloc ( ni * nd * sizeof ( double ) );

  for ( j = 0; j < nd; j++ )
  {
    for ( i = 0; i < ni; i++ )
    {
      lb[i+j*ni] = 1.0;
    }
  }

  for ( i = 0; i < nd; i++ )
  {
    for ( j = 0; j < nd; j++ )
    {
      if ( j != i )
      {
        for ( k = 0; k < ni; k++ )
        {
          lb[k+i*ni] = lb[k+i*ni] * ( xi[k] - xd[j] ) / ( xd[i] - xd[j] );
        }
      }
    }
  }

  return lb;
}
/******************************************************************************/

double *lagrange_interp_nd_grid ( int m, int n_1d[], double a[], double b[], 
  int nd )

/******************************************************************************/
/*
  Purpose:

    LAGRANGE_INTERP_ND_GRID sets an M-dimensional Lagrange interpolant grid.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    29 September 2012

  Author:

    John Burkardt

  Parameters:

    Input, int M, the spatial dimension.

    Input, int N_1D[M], the order of the 1D rule to be used
    in each dimension.

    Input, double A[M], B[M], the upper and lower limits.

    Input, int ND, the number of points in the product grid.

    Output, double LAGRANGE_INTERP_ND_GRID[M*ND], the points at which data 
    is to be sampled.
*/
{
  int i;
  int j;
  int n;
  double *x_1d;
  double *xd;
/*
  Compute the data points.
*/
  xd = ( double * ) malloc ( m * nd * sizeof ( double ) );

  for ( j = 0; j < nd; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      xd[i+j*m] = 0.0;
    }
  }

  for ( i = 0; i < m; i++ )
  {
    n = n_1d[i];
    x_1d = cc_compute_points ( n );
    for ( j = 0; j < n; j++ )
    {
      x_1d[j] = 0.5 * ( ( 1.0 - x_1d[j] ) * a[i] 
                      + ( 1.0 + x_1d[j] ) * b[i] );
    }
    r8vec_direct_product ( i, n, x_1d, m, nd, xd );
    free ( x_1d );
  }

  return xd;
}
/******************************************************************************/

double *lagrange_interp_nd_grid2 ( int m, int ind[], double a[], double b[], 
  int nd )

/******************************************************************************/
/*
  Purpose:

    LAGRANGE_INTERP_ND_GRID2 sets an M-dimensional Lagrange interpolant grid.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    29 September 2012

  Author:

    John Burkardt

  Parameters:

    Input, int M, the spatial dimension.

    Input, int IND[M], the index or level of the 1D rule 
    to be used in each dimension.

    Input, double A[M], B[M], the upper and lower limits.

    Input, int ND, the number of points in the product grid.

    Output, double LAGRANGE_INTERP_ND_GRID2[M*ND], the points at which data 
    was sampled.
*/
{
  int i;
  int j;
  int n;
  double *x_1d;
  double *xd;
/*
  Compute the data points.
*/
  xd = ( double * ) malloc ( m * nd * sizeof ( double ) );

  for ( j = 0; j < nd; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      xd[i+j*m] = 0.0;
    }
  }

  for ( i = 0; i < m; i++ )
  {
    n = order_from_level_135 ( ind[i] );
    x_1d = cc_compute_points ( n );
    for ( j = 0; j < n; j++ )
    {
      x_1d[j] = 0.5 * ( ( 1.0 - x_1d[j] ) * a[i] 
                      + ( 1.0 + x_1d[j] ) * b[i] );
    }
    r8vec_direct_product ( i, n, x_1d, m, nd, xd );
    free ( x_1d );
  }

  return xd;
}
/******************************************************************************/

int lagrange_interp_nd_size ( int m, int n_1d[] )

/******************************************************************************/
/*
  Purpose:

    LAGRANGE_INTERP_ND_SIZE sizes an M-dimensional Lagrange interpolant.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    29 September 2012

  Author:

    John Burkardt

  Parameters:

    Input, int M, the spatial dimension.

    Input, int N_1D[M], the order of the 1D rule to be used
    in each dimension.

    Output, int LAGRANGE_INTERP_ND_SIZE, the number of points in the product grid.
*/
{
  int nd;
/*
  Determine the number of data points.
*/
  nd = i4vec_product ( m, n_1d );

  return nd;
}
/******************************************************************************/

int lagrange_interp_nd_size2 ( int m, int ind[] )

/******************************************************************************/
/*
  Purpose:

    LAGRANGE_INTERP_ND_SIZE2 sizes an M-dimensional Lagrange interpolant.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    29 September 2012

  Author:

    John Burkardt

  Parameters:

    Input, int M, the spatial dimension.

    Input, int IND[M], the index or level of the 1D rule 
    to be used in each dimension.

    Output, int LAGRANGE_INTERP_ND_SIZE2, the number of points in the product grid.
*/
{
  int i;
  int n;
  int nd;
/*
  Determine the number of data points.
*/
  nd = 1;
  for ( i = 0; i < m; i++ )
  {
    n = order_from_level_135 ( ind[i] );
    nd = nd * n;
  }

  return nd;
}
/******************************************************************************/

double *lagrange_interp_nd_value ( int m, int n_1d[], double a[], double b[], 
  int nd, double zd[], int ni, double xi[] )

/******************************************************************************/
/*
  Purpose:

    LAGRANGE_INTERP_ND_VALUE evaluates an ND Lagrange interpolant.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 September 2012

  Author:

    John Burkardt

  Parameters:

    Input, int M, the spatial dimension.

    Input, int N_1D[M], the order of the 1D rule to be used
    in each dimension.

    Input, double A[M], B[M], the upper and lower limits.

    Input, int ND, the number of points in the product grid.

    Input, double ZD[ND], the function evaluated at the points XD.

    Input, int NI, the number of points at which the 
    interpolant is to be evaluated.

    Input, double XI[M*NI], the points at which the interpolant is 
    to be evaluated.

    Output, double LAGRANGE_INTERP_ND_VALUE[NI], the interpolant evaluated 
    at the points XI.
*/
{
  int i;
  int j;
  int k;
  int n;
  double *value;
  double *w;
  double *x_1d;
  double *zi;

  w = ( double * ) malloc ( nd * sizeof ( double ) );
  zi = ( double * ) malloc ( ni * sizeof ( double ) );

  for ( j = 0; j < ni; j++ )
  {
    for ( i = 0; i < nd; i++ )
    {
      w[i] = 1.0;
    }
    for ( i = 0; i < m; i++ )
    {
      n = n_1d[i];
      x_1d = cc_compute_points ( n );
      for ( k = 0; k < n; k++ )
      {
        x_1d[k] = 0.5 * ( ( 1.0 - x_1d[k] ) * a[i] 
                        + ( 1.0 + x_1d[k] ) * b[i] );
      }
      value = lagrange_basis_1d ( n, x_1d, 1, xi+i+j*m );
      r8vec_direct_product2 ( i, n, value, m, nd, w );
      free ( value );
      free ( x_1d );
    }
    zi[j] = r8vec_dot_product ( nd, w, zd );
  }

  free ( w );

  return zi;
}
/******************************************************************************/

double *lagrange_interp_nd_value2 ( int m, int ind[], double a[], double b[], 
  int nd, double zd[], int ni, double xi[] )

/******************************************************************************/
/*
  Purpose:

    LAGRANGE_INTERP_ND_VALUE2 evaluates an ND Lagrange interpolant.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 September 2012

  Author:

    John Burkardt

  Parameters:

    Input, int M, the spatial dimension.

    Input, int IND[M], the index or level of the 1D rule 
    to be used in each dimension.

    Input, double A[M], B[M], the upper and lower limits.

    Input, int ND, the number of points in the product grid.

    Input, double ZD[ND], the function evaluated at the points XD.

    Input, int NI, the number of points at which the 
    interpolant is to be evaluated.

    Input, double XI[M*NI], the points at which the interpolant 
    is to be evaluated.

    Output, double ZI[NI], the interpolant evaluated at the 
    points XI.
*/
{
  int i;
  int j;
  int k;
  int n;
  double *value;
  double *w;
  double *x_1d;
  double *zi;

  w = ( double * ) malloc ( nd * sizeof ( double ) );
  zi = ( double * ) malloc ( ni * sizeof ( double ) );

  for ( j = 0; j < ni; j++ )
  {
    for ( i = 0; i < nd; i++ )
    {
      w[i] = 1.0;
    }

    for ( i = 0; i < m; i++ )
    {
      n = order_from_level_135 ( ind[i] );
      x_1d = cc_compute_points ( n );
      for ( k = 0; k < n; k++ )
      {
        x_1d[k] = 0.5 * ( ( 1.0 - x_1d[k] ) * a[i] 
                        + ( 1.0 + x_1d[k] ) * b[i] );
      }
      value = lagrange_basis_1d ( n, x_1d, 1, xi+i+j*m );
      r8vec_direct_product2 ( i, n, value, m, nd, w );
      free ( value );
      free ( x_1d );
    }
    zi[j] = r8vec_dot_product ( nd, w, zd );
  }

  free ( w );

  return zi;
}
/******************************************************************************/

int order_from_level_135 ( int l )

/******************************************************************************/
/*
  Purpose:

    ORDER_FROM_LEVEL_135 evaluates the 135 level-to-order relationship.

  Discussion:

    Clenshaw Curtis rules, and some others, often use the following
    scheme:

    L: 0  1  2  3   4   5
    N: 1  3  5  9  17  33 ... 2^L+1

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    29 September 2012

  Author:

    John Burkardt

  Parameters:

    Input, int L, the level, which should be 0 or greater.

    Output, int ORDER_FROM_LEVEL_135, the order.
*/
{
  int n;

  if ( l < 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "ORDER_FROM_LEVEL_135 - Fatal error!\n" );
    fprintf ( stderr, "  Illegal input value of L!\n" );
    exit ( 1 );
  }
  else if ( l == 0 )
  {
    n = 1;
  }
  else
  {
    n = i4_power ( 2, l ) + 1;
  }
  return n;
}
