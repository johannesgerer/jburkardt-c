# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <string.h>

# include "sparse_grid_cc.h"

int main ( );

void test01 ( int dim_min, int dim_max, int level_max_min, int level_max_max );
void test015 ( int dim_min, int dim_max, int level_max_min, int level_max_max );
void test02 ( int dim_num, int level_max );
void test03 ( int dim_num, int level_max );
void test04 ( int dim_num, int level_max );
void test05 ( int dim_num, int level_max, int degree_max );
void test06 ( int dim_num, int level_max );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for SPARSE_GRID_CC_PRB.

  Discussion:

    SPARSE_GRID_CC_PRB tests the SPARSE_GRID_CC library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 March 2013

  Author:

    John Burkardt

  Reference:

    Fabio Nobile, Raul Tempone, Clayton Webster,
    A Sparse Grid Stochastic Collocation Method for Partial Differential
    Equations with Random Input Data,
    SIAM Journal on Numerical Analysis,
    Volume 46, Number 5, 2008, pages 2309-2345.
*/
{
  int dim_max;
  int dim_min;
  int dim_num;
  int level_max;
  int level_max_max;
  int level_max_min;

  timestamp ( );
  printf ( "\n" );
  printf ( "SPARSE_GRID_CC_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the SPARSE_GRID_CC library.\n" );
/*
  Count number of points in sparse rule from DIM_MIN to DIM_MAX, LEVEL_MAX_MAX.
*/
  dim_min = 1;
  dim_max = 5;
  level_max_min = 0;
  level_max_max = 10;
  test01 ( dim_min, dim_max, level_max_min, level_max_max );

  printf ( "\n" );
  timestamp ( );

  dim_min = 6;
  dim_max = 10;
  level_max_min = 0;
  level_max_max = 10;
  test01 ( dim_min, dim_max, level_max_min, level_max_max );

  printf ( "\n" );
  timestamp ( );

  dim_min = 100;
  dim_max = 100;
  level_max_min = 0;
  level_max_max = 4;
  test01 ( dim_min, dim_max, level_max_min, level_max_max );

  printf ( "\n" );
  timestamp ( );
/*
  Count number of points in sparse rule from DIM_MIN to DIM_MAX, LEVEL_MAX_MAX.
*/
  dim_min = 1;
  dim_max = 5;
  level_max_min = 0;
  level_max_max = 10;
  test015 ( dim_min, dim_max, level_max_min, level_max_max );

  printf ( "\n" );
  timestamp ( );

  dim_min = 6;
  dim_max = 10;
  level_max_min = 0;
  level_max_max = 10;
  test015 ( dim_min, dim_max, level_max_min, level_max_max );

  printf ( "\n" );
  timestamp ( );

  dim_min = 100;
  dim_max = 100;
  level_max_min = 0;
  level_max_max = 4;
  test015 ( dim_min, dim_max, level_max_min, level_max_max );

  printf ( "\n" );
  timestamp ( );
/*
  Compute abstract grid indices of sparse grid points as selected from product grid
  for DIMENSION, LEVEL_MAX.
*/
  test02 ( 2, 3 );
  test02 ( 2, 4 );
  test02 ( 3, 0 );
  test02 ( 3, 2 );
  test02 ( 6, 2 );
/*
  Compute sparse Clenshaw-Curtis rule for DIMENSION, LEVEL_MAX.
*/
  test03 ( 2, 3 );
  test03 ( 3, 0 );
  test03 ( 3, 1 );
/*
  Test sum of weights for DIMENSION, LEVEL_MAX.
*/
  test04 ( 2, 4 );
  test04 ( 3, 0 );
  test04 ( 3, 1 );
  test04 ( 3, 6 );
  test04 ( 10, 3 );
/*
  Test monomial exactness for DIMENSION, LEVEL_MAX, DEGREE_MAX.
*/ 
  test05 ( 2, 0, 3 );
  test05 ( 2, 1, 5 );
  test05 ( 2, 2, 7 );
  test05 ( 2, 3, 9 );
  test05 ( 2, 4, 11 );
  test05 ( 2, 5, 13 );
  
  test05 ( 3, 0, 2 );
  test05 ( 3, 1, 4 );
  test05 ( 3, 2, 6 );
  test05 ( 3, 3, 8 );
/*
  Show how to write a rule to a file.
*/
  dim_num = 2;
  level_max = 3;

  test06 ( dim_num, level_max );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "SPARSE_GRID_CC_PRB\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );
  
  return 0;
}
/******************************************************************************/

void test01 ( int dim_min, int dim_max, int level_max_min, int level_max_max )

/******************************************************************************/
/*
  Purpose:

    TEST01 tests SPARSE_GRID_CFN_SIZE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 March 2013

  Author:

    John Burkardt

  Parameters:

    Input, int DIM_MIN, the minimum spatial dimension to consider.

    Input, int DIM_MAX, the maximum spatial dimension to consider.

    Input, int LEVEL_MAX_MIN, the minimum value of LEVEL_MAX to consider.

    Input, int LEVEL_MAX_MAX, the maximum value of LEVEL_MAX to consider.
*/
{
  int dim_num;
  int level_max;
  int point_num;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  SPARSE_GRID_CFN_SIZE returns the number of distinct\n" );
  printf ( "  points in a sparse grid of Closed Fully Nested rules.\n" );
  printf ( "\n" );
  printf ( "  Each sparse grid is of spatial dimension DIM,\n" );
  printf ( "  and is made up of all product grids of levels up to LEVEL_MAX.\n" );
  printf ( "\n" );
  printf ( "   DIM: " );
  for ( dim_num = dim_min; dim_num <= dim_max; dim_num++)
  {
    printf ( "  %8d", dim_num );
  }
  printf ( "\n" );
  printf ( "\n" );
  printf ( "   LEVEL_MAX\n" );
  printf ( "\n" );

  for ( level_max = level_max_min; level_max <= level_max_max; level_max++ )
  {
    printf ( "    %4d", level_max );
    for ( dim_num = dim_min; dim_num <= dim_max; dim_num++ )
    {
      point_num = sparse_grid_cfn_size ( dim_num, level_max );
      printf ( "  %8d", point_num );
    }
    printf ( "\n" );
  }

  return;
}
/******************************************************************************/

void test015 ( int dim_min, int dim_max, int level_max_min, int level_max_max )

/******************************************************************************/
/*
  Purpose:

    TEST015 tests SPARSE_GRID_CCS_SIZE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 March 2013

  Author:

    John Burkardt

  Parameters:

    Input, int DIM_MIN, the minimum spatial dimension to consider.

    Input, int DIM_MAX, the maximum spatial dimension to consider.

    Input, int LEVEL_MAX_MIN, the minimum value of LEVEL_MAX to consider.

    Input, int LEVEL_MAX_MAX, the maximum value of LEVEL_MAX to consider.
*/
{
  int dim_num;
  int level_max;
  int point_num;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  SPARSE_GRID_CCS_SIZE returns the number of distinct\n" );
  printf ( "  points in a Clenshaw Curtis Slow-Growth sparse grid.\n" );
  printf ( "\n" );
  printf ( "  Each sparse grid is of spatial dimension DIM,\n" );
  printf ( "  and is made up of all product grids of levels up to LEVEL_MAX.\n" );
  printf ( "\n" );
  printf ( "   DIM: " );
  for ( dim_num = dim_min; dim_num <= dim_max; dim_num++)
  {
    printf ( "  %8d", dim_num );
  }
  printf ( "\n" );
  printf ( "\n" );
  printf ( "   LEVEL_MAX\n" );
  printf ( "\n" );

  for ( level_max = level_max_min; level_max <= level_max_max; level_max++ )
  {
    printf ( "    %4d", level_max );
    for ( dim_num = dim_min; dim_num <= dim_max; dim_num++ )
    {
      point_num = sparse_grid_ccs_size ( dim_num, level_max );
      printf ( "  %8d", point_num );
    }
    printf ( "\n" );
  }

  return;
}
/******************************************************************************/

void test02 ( int dim_num, int level_max )

/******************************************************************************/
/*
  Purpose:

    TEST02 tests SPARSE_GRID_CC_INDEX.

  Discussion:

    The routine computes the indices of the unique points used in a sparse 
    multidimensional grid whose size is controlled by a parameter LEVEL_MAX.

    Once these indices are returned, they can be converted into
    Clenshaw Curtis points.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 March 2013

  Author:

    John Burkardt

  Parameters:

    Input, int DIM_NUM, the spatial dimension.

    Input, int LEVEL_MAX, the level.
*/
{
  int dim;
  int *grid_index;
  int point;
  int point_num;

  printf ( "\n" );
  printf ( "TEST02:\n" );
  printf ( "  SPARSE_GRID_CC_INDEX returns all grid indexes\n" );
  printf ( "  whose level value satisfies\n" );
  printf ( "    0 <= LEVEL <= LEVEL_MAX.\n" );
  printf ( "  Here, LEVEL is the sum of the levels of the 1D rules,\n" );
  printf ( "  and the order of the rule is 2**LEVEL + 1.\n" );

  printf ( "\n" );
  printf ( "  LEVEL_MAX = %d\n", level_max );
  printf ( "  Spatial dimension DIM_NUM = %d\n", dim_num );

  point_num = sparse_grid_cfn_size ( dim_num, level_max );

  printf ( "\n" );
  printf ( "  Number of unique points in the grid = %d\n", point_num );
/*
  Compute the orders and points.
*/
  grid_index = sparse_grid_cc_index ( dim_num, level_max, point_num );
/*
  Now we're done.  Print the merged grid data.
*/
  printf ( "\n" );
  printf ( "  Grid index:\n" );
  printf ( "\n" );
  for ( point = 0; point < point_num; point++ )
  {
    printf ( "  %4d  ", point );
    for ( dim = 0; dim < dim_num; dim++ )
    {
      printf ( "%6d", grid_index[dim+point*dim_num] );
    }
    printf ( "\n" );
  }

  free ( grid_index );

  return;
}
/******************************************************************************/

void test03 ( int dim_num, int level_max )

/******************************************************************************/
/*
  Purpose:

    TEST03 call SPARSE_GRID_CC to create a Clenshaw Curtis grid.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 March 2013

  Author:

    John Burkardt

  Parameters:

    Input, int DIM_NUM, the spatial dimension.

    Input, int LEVEL_MAX, the level.
*/
{
  int dim;
  double *grid_point;
  double *grid_weight;
  int point;
  int point_num;

  printf ( "\n" );
  printf ( "TEST03:\n" );
  printf ( "  SPARSE_GRID_CC makes a sparse Clenshaw Curtis grid.\n" );
  printf ( "\n" );
  printf ( "  LEVEL_MAX = %d\n", level_max );
  printf ( "  Spatial dimension DIM_NUM = %d\n", dim_num );
/*
  Determine the number of points.
*/
  point_num = sparse_grid_cfn_size ( dim_num, level_max );

  printf ( "\n" );
  printf ( "  Number of unique points in the grid = %d\n", point_num );
/*
  Allocate space for the weights and points.
*/
  grid_weight = ( double * ) malloc ( point_num * sizeof ( double ) );
  grid_point = ( double * ) malloc ( dim_num * point_num * sizeof ( double ) );
/*
  Compute the weights and points.
*/
  sparse_grid_cc ( dim_num, level_max, point_num, grid_weight, grid_point );
/*
  Print them out.
*/
  printf ( "\n" );
  printf ( "  Grid weights:\n" );
  printf ( "\n" );
  for ( point = 0; point < point_num; point++ )
  {
    printf ( "  %4d  %10f\n", point, grid_weight[point] );
  }
  
  printf ( "\n" );
  printf ( "  Grid points:\n" );
  printf ( "\n" );
  for ( point = 0; point < point_num; point++ )
  {
    printf ( "  %4d", point );
    for ( dim = 0; dim < dim_num; dim++ )
    {
      printf ( "%10f", grid_point[dim+point*dim_num] );
    }
    printf ( "\n" );
  }

  free ( grid_point );
  free ( grid_weight );

  return;
}
/******************************************************************************/

void test04 ( int dim_num, int level_max )

/******************************************************************************/
/*
  Purpose:

    TEST04 sums the weights and compares them to 2^DIM_NUM.

  Discussion:

    This routine gets the sparse grid indices and determines the 
    corresponding sparse grid abscissas.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 March 2013

  Author:

    John Burkardt

  Parameters:

    Input, int DIM_NUM, the spatial dimension.

    Input, int LEVEL_MAX, the level.
*/
{
  double *grid_point;
  double *grid_weight;
  int point;
  int point_num;
  double weight_sum;
  double weight_sum_error;
  double weight_sum_exact;

  printf ( "\n" );
  printf ( "TEST04:\n" );
  printf ( "  Compute the weights of a Clenshaw Curtis sparse grid .\n" );
  printf ( "\n" );
  printf ( "  As a simple test, sum these weights.\n" );
  printf ( "  They should sum to exactly 2^DIM_NUM.\n" );
  printf ( "\n" );
  printf ( "  LEVEL_MAX = %d\n", level_max );
  printf ( "  Spatial dimension DIM_NUM = %d\n", dim_num );
/*
  Determine the number of points.
*/
  point_num = sparse_grid_cfn_size ( dim_num, level_max );

  printf ( "\n" );
  printf ( "  Number of unique points in the grid = %d\n", point_num );
/*
  Allocate space for the weights and points.
*/
  grid_weight = ( double * ) malloc ( point_num * sizeof ( double ) );
  grid_point = ( double * ) malloc ( dim_num * point_num * sizeof ( double ) );
/*
  Compute the weights and points.
*/
  sparse_grid_cc ( dim_num, level_max, point_num, grid_weight, grid_point );
/*
  Sum the weights.
*/
  weight_sum = 0.0;
  for ( point = 0; point < point_num; point++ )
  {
    weight_sum = weight_sum + grid_weight[point];  }
  
  weight_sum_exact = pow ( 2.0, dim_num );
  
  weight_sum_error = fabs ( weight_sum - weight_sum_exact );
  
  printf ( "\n" );
  printf ( "    Weight sum     Exact sum    Difference\n" );
  printf ( "\n" );
  printf ( "  %12g  %12g  %12g\n", weight_sum, weight_sum_exact, weight_sum_error );

  free ( grid_point );
  free ( grid_weight );

  return;
}
/******************************************************************************/

void test05 ( int dim_num, int level_max, int degree_max )

/******************************************************************************/
/*
  Purpose:

    TEST05 tests a Clenshaw Curtis sparse grid rule for monomial exactness.

  Discussion:

    This test is going to check EVERY monomial of total degree DEGREE_MAX
    or less.  Even for a moderately high dimension of DIM_NUM = 10, you
    do NOT want to use a large value of DEGREE_MAX, since there are

      1         monomials of total degree 0,
      DIM_NUM   monomials of total degree 1,
      DIM_NUM^2 monomials of total degree 2,
      DIM_NUM^3 monomials of total degree 3, and so on.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 March 2013

  Author:

    John Burkardt

  Parameters:

    Input, int DIM_NUM, the spatial dimension.

    Input, int LEVEL_MAX, the level.

    Input, int DEGREE_MAX, the maximum monomial total degree to check.
*/
{
  int degree;
  int dim;
  int error;
  int *expon;
  double *grid_point;
  double *grid_weight;
  int h;
  int last;
  int more;
  int point;
  int point_num;
  double quad_error;
  int t;
  double volume;

  printf ( "\n" );
  printf ( "TEST05\n" );
  printf ( "  Check the exactness of a Clenshaw Curtis sparse grid quadrature rule,\n" );
  printf ( "  applied to all monomials of orders 0 to DEGREE_MAX.\n" );
  printf ( "\n" );
  printf ( "  LEVEL_MAX = %d\n", level_max );
  printf ( "  Spatial dimension DIM_NUM = %d\n", dim_num );
  printf ( "\n" );
  printf ( "  The maximum total degree to be checked is DEGREE_MAX = %d\n", degree_max );
  printf ( "\n" );
  printf ( "  We expect this rule to be accurate up to and including total degree %d\n", 
	   2 * level_max + 1 );
/*
  Determine the number of points in the rule.
*/
  point_num = sparse_grid_cfn_size ( dim_num, level_max );

  printf ( "\n" );
  printf ( "  Number of unique points in the grid = %d\n", point_num );
/*
  Allocate space for the weights and points.
*/
  grid_weight = ( double * ) malloc ( point_num * sizeof ( double ) );
  grid_point = ( double * ) malloc ( dim_num * point_num * sizeof ( double ) );
/*
  Compute the weights and points.
*/
  sparse_grid_cc ( dim_num, level_max, point_num, grid_weight, 
  grid_point );
/*
  Rescale the weights, and translate the abscissas.
*/
  volume = pow ( 2.0, dim_num );

  for ( point = 0; point < point_num; point++ )
  {
    grid_weight[point] = grid_weight[point] / volume;
  }

  for ( dim = 0; dim < dim_num; dim++ )
  {
    for ( point = 0; point < point_num; point++ )
    {
      grid_point[dim+point*dim_num] = ( grid_point[dim+point*dim_num] + 1.0 ) 
      / 2.0; 
    }
  }
/*
  Explore the monomials.
*/
  expon = ( int * ) malloc ( dim_num * sizeof ( int ) );

  printf ( "\n" );
  printf ( "      Error      Total   Monomial\n" );
  printf ( "                 Degree  Exponents\n" );
  printf ( "\n" );

  for ( degree = 0; degree <= degree_max; degree++ )
  {
    more = 0;

    for ( ; ; )
    {
      comp_next ( degree, dim_num, expon, &more, &h, &t );

      quad_error = monomial_quadrature ( dim_num, expon, point_num, 
        grid_weight, grid_point );

      printf ( "  %12g     %2d      ", quad_error, degree );

      for ( dim = 0; dim < dim_num; dim++ )
      {
        printf ( "%2d", expon[dim] );
      }
      printf ( "\n" );

      if ( !more )
      {
        break;
      }
    }
    printf ( "\n" );
  }

  free ( expon );
  free ( grid_point );
  free ( grid_weight );

  return;
}
/*****************************************************************************/

void test06 ( int dim_num, int level_max )

/*****************************************************************************/
/*
  Purpose:

    TEST06 creates a sparse Clenshaw-Curtis grid and writes it to a file.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 March 2013

  Author:

    John Burkardt

  Parameters:

    Input, int DIM_NUM, the spatial dimension.

    Input, int LEVEL_MAX, the level.
*/
{
  int dim;
  int point;
  int point_num;
  double *r;
  char r_filename[80];
  double *w;
  char w_filename[80];
  double *x;
  char x_filename[80];

  printf ( "\n" );
  printf ( "TEST06:\n" );
  printf ( "  Call SPARSE_GRID_CC to make a sparse Clenshaw-Curtis grid.\n" );
  printf ( "  Write the data to a set of quadrature files.\n" );
    
  printf ( "\n" );
  printf ( "  LEVEL_MAX = %d\n", level_max );
  printf ( "  Spatial dimension DIM_NUM = %d\n", dim_num );
/*
  Determine the number of points.
*/
  point_num = sparse_grid_cfn_size ( dim_num, level_max );
/*
  Allocate space for the weights and points.
*/
  r = ( double * ) malloc ( dim_num * 2 * sizeof ( double ) );
  w = ( double * ) malloc ( point_num * sizeof ( double ) );
  x = ( double * ) malloc ( dim_num * point_num * sizeof ( double ) );
/*
  Compute the weights and points.
*/
  for ( dim = 0; dim < dim_num; dim++ )
  {
    r[dim+0*dim_num] = -1.0;
    r[dim+1*dim_num] = +1.0;
  }

  sparse_grid_cc ( dim_num, level_max, point_num, w, x );
/*
  Write the data out.
*/
  sprintf ( r_filename, "cc_d%d_level%d_r.txt", dim_num, level_max );
  sprintf ( w_filename, "cc_d%d_level%d_w.txt", dim_num, level_max );
  sprintf ( x_filename, "cc_d%d_level%d_x.txt", dim_num, level_max );

  r8mat_write ( r_filename, dim_num, 2,         r );
  r8mat_write ( w_filename, 1,       point_num, w );
  r8mat_write ( x_filename, dim_num, point_num, x );

  printf ( "\n" );
  printf ( "  R data written to \"%s\".\n", r_filename );
  printf ( "  W data written to \"%s\".\n", w_filename);
  printf ( "  X data written to \"%s\".\n", x_filename );

  free ( r );
  free ( w );
  free ( x );

  return;
}

