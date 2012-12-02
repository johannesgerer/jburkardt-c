# include <stdlib.h>
# include <stdio.h>
# include <string.h>
# include <math.h>
# include <time.h>

# include "sparse_interp_nd.h"
# include "r8lib.h"

int main ( );
void test01 ( int m, int sparse_max );
double *f_sinr ( int m, int n, double x[] );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    SPARSE_INTERP_ND_PRB tests SPARSE_INTERP_ND.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    06 October 2012

  Author:

    John Burkardt
*/
{
  int m;
  int sparse_max;

  timestamp ( );
  printf ( " \n" );
  printf ( "SPARSE_INTERP_ND_PRB\n" );
  printf ( "  C version.\n" );
  printf ( "  Test the SPARSE_INTERP_ND library.\n" );
  printf ( "  The R8LIB library is also required.\n" );

  m = 1;
  sparse_max = 9;
  test01 ( m, sparse_max );

  m = 2;
  sparse_max = 9;
  test01 ( m, sparse_max );

  m = 3;
  sparse_max = 9;
  test01 ( m, sparse_max );

  m = 4;
  sparse_max = 7;
  test01 ( m, sparse_max );
/*
  Terminate.
*/
  printf ( " \n" );
  printf ( "SPARSE_INTERP_ND_PRB\n" );
  printf ( "  Normal end of execution.\n" );

  printf ( " \n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void test01 ( int m, int sparse_max )

/******************************************************************************/
/*
  Purpose:

    TEST01: sequence of sparse interpolants to an M-dimensional function.

  Discussion:

    We have functions that can generate a Lagrange interpolant to data
    in M dimensions, with specified order or level in each dimension.

    We use the Lagrange function as the inner evaluator for a sparse
    grid procedure. 

    The procedure computes sparse interpolants of levels 0 to SPARSE_MAX
    to a given function.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    06 October 2012

  Author:

    John Burkardt

  Parameters:

    Local, int M, the spatial dimension.

    Input, int SPARSE_MAX, the maximum sparse grid level to try.

  Local Parameters:

    Local, double A[M], B[M], the upper and lower variable limits 
    in each dimension.

    Local, double APP_ERROR, the averaged Euclidean norm of the 
    difference between the sparse interpolant and the exact function at 
    the interpolation points.

    Local, int C[L_MAX+1], the sparse grid coefficient vector.
    Results at level L are weighted by C(L).

    Local, int IND[M], the 1D indices defining a Lagrange grid.
    Each entry is a 1d "level" that specifies the order of a 
    Clenshaw Curtis 1D grid.

    Local, int L, the current Lagrange grid level.

    Local, int L_MAX, the current sparse grid level.

    Local, int MORE, used to control the enumeration of all the
    Lagrange grids at a current grid level.

    Local, int ND, the number of points used in a Lagrange grid.

    Local, int ND_TOTAL, the total number of points used in all the
    Lagrange interpolants for a given sparse interpolant points that occur
    in more than one Lagrange grid are counted multiple times.

    Local, int NI, the number of interpolant evaluation points.

    Local, int SPARSE_MIN, the minimum sparse grid level to try.

    Local, double XD[M*ND], the data points for a Lagrange grid.

    Local, double XI[M*NI], the interpolant evaluation points.

    Local, double ZD[ND], the data values for a Lagrange grid.

    Local, double ZE[NI], the exact function values at XI.

    Local, double ZI[NI], the sparse interpolant values at XI.

    Local, double ZPI[NI], one set of Lagrange interpolant values at XI.
*/
{
  double *a;
  double app_error;
  double *b;
  int *c;
  int h;
  int i;
  int *ind;
  int l;
  int l_max;
  int l_min;
  int more;
  int nd;
  int nd_total;
  int ni;
  int seed;
  int sparse_min;
  int t;
  int *w;
  double *xd;
  double *xi;
  double *zd;
  double *ze;
  double *zi;
  double *zpi;

  printf ( "\n" );
  printf ( "TEST01:\n" );
  printf ( "  Sparse interpolation for a function f(x) of M-dimensional argument.\n" );
  printf ( "  Use a sequence of sparse grids of levels 0 through SPARSE_MAX.\n" );
  printf ( "  Invoke a general Lagrange interpolant function to do this.\n" );
  printf ( "\n" );
  printf ( "  Compare the exact function and the interpolants at a grid of points.\n" );
  printf ( "\n" );
  printf ( "  The \"order\" is the sum of the orders of all the product grids\n" );
  printf ( "  used to make a particular sparse grid.\n" );
/*
  User input.
*/
  printf ( "\n" );
  printf ( "  Spatial dimension M = %d\n", m );
  printf ( "  Maximum sparse grid level = %d\n", sparse_max );
/*
  Define the region.
*/
  a = ( double * ) malloc ( m * sizeof ( double ) );
  b = ( double * ) malloc ( m * sizeof ( double ) );

  for ( i = 0; i < m; i++ )
  {
    a[i] = 0.0;
    b[i] = 1.0;
  }
/*
  Define the interpolation evaluation information.
*/
  ni = 100;
  seed = 123456789;
  xi = r8mat_uniform_abvec_new ( m, ni, a, b, &seed );

  printf ( "  Number of interpolation points is NI = %d\n", ni );

  ze = f_sinr ( m, ni, xi );
/*
  Compute a sequence of sparse grid interpolants of increasing level.
*/
  printf ( "\n" );
  printf ( "   L     Order    ApproxError\n" );
  printf ( "\n" );

  ind = ( int * ) malloc ( m * sizeof ( int ) );
  zi = ( double * ) malloc ( ni * sizeof ( double ) );

  sparse_min = 0;

  for ( l_max = sparse_min; l_max <= sparse_max; l_max++ )
  {
    c = ( int * ) malloc ( ( l_max + 1 ) * sizeof ( int ) );
    w = ( int * ) malloc ( ( l_max + 1 ) * sizeof ( int ) );
    smolyak_coefficients ( l_max, m, c, w );
    
    for ( i = 0; i < ni; i++ )
    {
      zi[i] = 0.0;
    }
    nd_total = 0;

    l_min = i4_max ( l_max + 1 - m, 0 );

    for ( l = l_min; l <= l_max; l++ )
    {
      more = 0;
      while ( 1 )
      {
/*
  Define the next product grid.
*/
        comp_next ( l, m, ind, &more, &h, &t );
/*
  Count the grid, find its points, evaluate the data there.
*/
        nd = lagrange_interp_nd_size2 ( m, ind );
        xd = lagrange_interp_nd_grid2 ( m, ind, a, b, nd );
        zd = f_sinr ( m, nd, xd );
/*
  Use the grid to evaluate the interpolant.
*/
        zpi = lagrange_interp_nd_value2 ( m, ind, a, b, nd, zd, ni, xi );
/*
  Weighted the interpolant values and add to the sparse grid interpolant.
*/
        nd_total = nd_total + nd;
        for ( i = 0; i < ni; i++ )
        {
          zi[i] = zi[i] + c[l] * zpi[i];
        }

        free ( xd );
        free ( zd );
        free ( zpi );

        if ( !more )
        {
          break;
        }
      }
    }
/*
  Compare sparse interpolant and exact function at interpolation points.
*/
    app_error = r8vec_norm_affine ( ni, zi, ze ) / ( double ) ( ni );

    printf ( "  %2d  %8d  %8.2e\n", l, nd_total, app_error );

    free ( c );
    free ( w );

  }

  free ( a );
  free ( b );
  free ( ind );
  free ( xi );
  free ( ze );
  free ( zi );

  return;
}
/******************************************************************************/

double *f_sinr ( int m, int n, double x[] )

/******************************************************************************/
/*
  Purpose:

    F_SINR is a scalar function of an M-dimensional argument, to be interpolated.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    06 October 2012

  Author:

    John Burkardt

  Parameters:

    Input, int M, the spatial dimension.

    Input, int N, the number of points.

    Input, double X(M,N), the points.

    Output, double F_SINR[N], the value of the function at each point.
*/
{
  int i;
  int j;
  double r;
  double *z;

  z = ( double * ) malloc ( n * sizeof ( double ) );

  for ( j = 0; j < n; j++ )
  {
    r = 0.0;
    for ( i = 0; i < m; i++ )
    {
      r = r + x[i+j*m] * x[i+j*m];
    }
    r = sqrt ( r );
    z[j] = sin ( r );
  }

  return z;
}
