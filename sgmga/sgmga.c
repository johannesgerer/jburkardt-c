# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "sandia_rules.h"
# include "sgmga.h"

/******************************************************************************/

void sgmga_aniso_normalize ( int option, int dim_num, double level_weight[] )

/******************************************************************************/
/*
  Purpose:

    SGMGA_ANISO_NORMALIZE normalizes the SGMGA anisotropic weight vector.

  Discussion:

    It is convenient for the user to initialize the anisotropic weight
    vector with any set of positive values.  These values are to be used
    as coefficients of the 1D levels, to evaluate an expression which 
    determines which 1D levels will be included in a given rule.

    This means that a relatively LARGE coefficient forces the corresponding 
    level to be relatively SMALL.  This is perhaps the opposite of what
    a user might expect.  If a user wishes to use an importance vector,
    so that a relatively large importance should correspond to more levels,
    and hence more points, in that dimension, then the function
    SGMGA_IMPORTANCE_TO_ANISO should be called first!

    Since the weights only represent the relative importance of the
    components, they may be multiplied by any (positive) scale factor.
    Nonetheless, it may be convenient ot choose a particular normalization
    for the weights.  

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    27 November 2009

  Author:

    John Burkardt

  Reference:

    Fabio Nobile, Raul Tempone, Clayton Webster,
    A Sparse Grid Stochastic Collocation Method for Partial Differential
    Equations with Random Input Data,
    SIAM Journal on Numerical Analysis,
    Volume 46, Number 5, 2008, pages 2309-2345.

    Fabio Nobile, Raul Tempone, Clayton Webster,
    An Anisotropic Sparse Grid Stochastic Collocation Method for Partial 
    Differential Equations with Random Input Data,
    SIAM Journal on Numerical Analysis,
    Volume 46, Number 5, 2008, pages 2411-2442.

  Parameters:

    Input, int OPTION, the normalization option.
    0, no scaling is applied.
    1, the weights are scaled so that the minimum nonzero entry is 1.
    2, the weights are scaled so that they sum to DIM_NUM.

    Input, int DIM_NUM, the spatial dimension.

    Input/output, double LEVEL_WEIGHT[DIM_NUM], the anisotropic
    weights.  The input values must be strictly positive.  
    On output, these have been normalized.
*/
{
  int dim;
  int found;
  double level_weight_min;
  double level_weight_sum;
/*
  Option 0, no normalization.
*/
  if ( option == 0 )
  {
  }
/*
  Option 1, the minimum nonzero entry is 1.
*/
  else if ( option == 1 )
  {
    level_weight_min = r8_huge ( );
    found = 0;
    for ( dim = 0; dim < dim_num; dim++ )
    {
      if ( 0.0 < level_weight[dim] )
      {
        if ( level_weight[dim] < level_weight_min )
        {
          level_weight_min = level_weight[dim];
          found = found + 1;
        }
      }
    }

    if ( found == 0 )
    {
      printf ( "\n" );
      printf ( "SGMGA_ANISO_NORMALIZE - Fatal error!\n" );
      printf ( "  Could not find a positive entry in LEVEL_WEIGHT.\n" );
      exit ( 1 );
    }

    for ( dim = 0; dim < dim_num; dim++ )
    {
      level_weight[dim] = level_weight[dim] / level_weight_min;
    }
  }
/*
  Option 2, rescale so sum of weights is DIM_NUM.
*/
  else if ( option == 2 )
  {
    level_weight_sum = r8vec_sum ( dim_num, level_weight );

    if ( level_weight_sum <= 0.0 )
    {
      printf ( "\n" );
      printf ( "SGMGA_ANISO_NORMALIZE - Fatal error!\n" );
      printf ( "  Sum of level weights is not positive.\n" );
      exit ( 1 );
    }
    for ( dim = 0; dim < dim_num; dim++ )
    {
      level_weight[dim] = ( ( double ) ( dim_num ) * level_weight[dim] )
        / level_weight_sum;
    }
  }

  return;
}
/******************************************************************************/

void sgmga_importance_to_aniso ( int dim_num, double importance[], 
  double level_weight[] )

/******************************************************************************/
/*
  Purpose:

    SGMGA_IMPORTANCE_TO_ANISO: importance vector to anisotropic weight vector.

  Discussion:

    To specify the anisotropy of a multidimensional problem, the user is
    allowed to specify an "importance vector".  This vector can contain
    any set of positive values.  These values represent the relative
    importance of each dimension.  These values, with a suitable normalization,
    will be used to evaluate a constraint of the following form:

      QMIN < Level(1) / Importance(1) + Level(2) / Importance(2) + ...
             Level(N) / Importance(N) <= QMAX

    and a set of levels that satisfies this constraint will then be included
    in a given anistotropic sparse grid rule.  Thus, increasing the
    importance value of a particular dimension allows larger level values
    in that dimension to satisfy the constraint.

    The program actually works with coefficients LEVEL_WEIGHT that are
    the inverse of the importance vector entries, with a suitable
    normalization.  This function is supplied to convert between the
    more natural "importance vector" and the internally useful 
    "level_weight" vector.

    This function converts the importance vector to an unnormalized 
    anisotropy weight vector.

    Note that some (but not all) of the IMPORTANCE vector entries may be zero.
    This indicates that the corresponding dimension is of "zero" or
    rather "minimal" importance.  In such a case, only a one-point quadrature
    rule will be applied for that dimension, no matter what sparse grid
    level is requested for the overall problem.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 November 2009

  Author:

    John Burkardt

  Reference:

    Fabio Nobile, Raul Tempone, Clayton Webster,
    A Sparse Grid Stochastic Collocation Method for Partial Differential
    Equations with Random Input Data,
    SIAM Journal on Numerical Analysis,
    Volume 46, Number 5, 2008, pages 2309-2345.

    Fabio Nobile, Raul Tempone, Clayton Webster,
    An Anisotropic Sparse Grid Stochastic Collocation Method for Partial 
    Differential Equations with Random Input Data,
    SIAM Journal on Numerical Analysis,
    Volume 46, Number 5, 2008, pages 2411-2442.

  Parameters:

    Input, int DIM_NUM, the spatial dimension.

    Input, double IMPORTANCE[DIM_NUM], the importance vector.
    All entries must be nonnegative, and at least one must be positive.

    Output, double LEVEL_WEIGHT[DIM_NUM], the anisotropic
    weights.
*/
{
  int dim;
  int found;
  double level_weight_norm;

  for ( dim = 0; dim < dim_num; dim++ )
  {
    if ( importance[dim] < 0.0 )
    {
      printf ( "\n" );
      printf ( "SGMGA_IMPORTANCE_TO_ANISO - Fatal error!\n" );
      printf ( "  Some IMPORTANCE entries are not positive.\n" );
      exit ( 1 );
    }
  }

  found = 0;

  for ( dim = 0; dim < dim_num; dim++ )
  {
    if ( 0.0 < importance[dim] )
    {
      level_weight[dim] = 1.0 / importance[dim];
      found = found + 1;
    }
    else
    {
      level_weight[dim] = 0.0;
    }
  }

  if ( found == 0 )
  {
    printf ( "\n" );
    printf ( "SGMGA_IMPORTANCE_TO_ANISO - Fatal error!\n" );
    printf ( "  No importance entry is positive.\n" );
    exit ( 1 );
  }

  return;
}
/******************************************************************************/

void sgmga_index ( int dim_num, double level_weight[], int level_max, 
  int rule[], int point_num, int point_total_num, int sparse_unique_index[],
  void level_to_order ( int dim_num, int level[], int rule[], int order[] ),
  int sparse_order[], int sparse_index[] )

/******************************************************************************/
/*
  Purpose:

    SGMGA_INDEX indexes an SGMGA grid.

  Discussion:

    For each "unique" point in the sparse grid, we return its INDEX and ORDER.

    That is, for the I-th unique point P, we determine the product grid which
    first generated this point, and  and we return in SPARSE_ORDER the orders 
    of the 1D rules in that grid, and  and in SPARSE_INDEX the component 
    indexes in those rules that generated this specific point.

    For instance, say P was first generated by a rule which was a 3D product
    of a 9th order CC rule and  and a 15th order GL rule, and  and that to 
    generate P, we used the 7-th point of the CC rule and  and the 3rh point 
    of the GL rule.  Then the SPARSE_ORDER information would be (9,15) and
    the SPARSE_INDEX information would be (7,3).  This, combined with the 
    information in RULE, is enough to regenerate the value of P.

    The user must preallocate space for the output arrays SPARSE_ORDER and
    SPARSE_INDEX.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 November 2009

  Author:

    John Burkardt

  Reference:

    Fabio Nobile, Raul Tempone, Clayton Webster,
    A Sparse Grid Stochastic Collocation Method for Partial Differential
    Equations with Random Input Data,
    SIAM Journal on Numerical Analysis,
    Volume 46, Number 5, 2008, pages 2309-2345.

    Fabio Nobile, Raul Tempone, Clayton Webster,
    An Anisotropic Sparse Grid Stochastic Collocation Method for Partial 
    Differential Equations with Random Input Data,
    SIAM Journal on Numerical Analysis,
    Volume 46, Number 5, 2008, pages 2411-2442.

  Parameters:

    Input, int DIM_NUM, the spatial dimension.

    Input, double LEVEL_WEIGHT[DIM_NUM], the anisotropic weights.

    Input, int LEVEL_MAX, the maximum value of LEVEL.

    Input, int RULE[DIM_NUM], the rule in each dimension.
     1, "CC",  Clenshaw Curtis, Closed Fully Nested rule.
     2, "F2",  Fejer Type 2, Open Fully Nested rule.
     3, "GP",  Gauss Patterson, Open Fully Nested rule.
     4, "GL",  Gauss Legendre, Open Weakly Nested rule.
     5, "GH",  Gauss Hermite, Open Weakly Nested rule.
     6, "GGH", Generalized Gauss Hermite, Open Weakly Nested rule.
     7, "LG",  Gauss Laguerre, Open Non Nested rule.
     8, "GLG", Generalized Gauss Laguerre, Open Non Nested rule.
     9, "GJ",  Gauss Jacobi, Open Non Nested rule.
    10, "GW",  Golub Welsch, (presumed) Open Non Nested rule.
    11, "CCS", Clenshaw Curtis "Slow", Closed Fully Nested rule.
    12, "F2S", Fejer Type 2 Slow, Closed Fully Nested rule.
    13, "GPS", Gauss Patterson Slow, Closed Fully Nested rule.

    Input, int POINT_NUM, the number of unique points 
    in the grid. 

    Input, int POINT_TOTAL_NUM, the total number of points in the grid.

    Input, int SPARSE_UNIQUE_INDEX[POINT_TOTAL_NUM], associates each
    point in the grid with its unique representative.

    Input, void LEVEL_TO_ORDER ( int dim_num, int level[], int rule[], 
    int order[] ), the function converting levels to orders.
    The choices are "level_to_order_default", "level_to_order_exponential" or 
    "level_to_order_linear".

    Output, int SPARSE_ORDER[DIM_NUM*POINT_NUM], lists, 
    for each point, the order of the 1D rules used in the grid that 
    generated it.

    Output, int SPARSE_INDEX[DIM_NUM*POINT_NUM)] lists, for 
    each point, its index in each of the 1D rules in the grid that generated 
    it.  The indices are 1-based.
*/
{
  double coef;
  int dim;
  int *level_1d;
  int *level_1d_max;
  double level_weight_min_pos;
  int more_grids;
  int more_points;
  int *order_1d;
  int point;
  int point_count;
  int *point_index;
  int point_unique;
  double q_max;
  double q_min;
/*
  Special cases.
*/
  if ( level_max < 0 )
  {
    return;
  }

  if ( level_max == 0 )
  {
    point = 0;
    for ( dim = 0; dim < dim_num; dim++ )
    {
      sparse_order[dim+point*dim_num] = 1;
      sparse_index[dim+point*dim_num] = 1;
    }
    return;
  }
/*
  Initialize the INDEX and ORDER arrays to -1 to help catch errors.
*/
  for ( point = 0; point < point_num; point++ )
  {
    for ( dim = 0; dim < dim_num; dim++ )
    {
      sparse_order[dim+point*dim_num] = -1;
      sparse_index[dim+point*dim_num] = -1;
    }
  }

  point_count = 0;

  level_1d = ( int * ) malloc ( dim_num * sizeof ( int ) );
  level_1d_max = ( int * ) malloc ( dim_num * sizeof ( int ) );
  order_1d = ( int * ) malloc ( dim_num * sizeof ( int ) );
  point_index = ( int * ) malloc ( dim_num * sizeof ( int ) );
/*
  Initialization for SGMGA_VCN_ORDERED.
*/
  level_weight_min_pos = r8vec_min_pos ( dim_num, level_weight );
  q_min = ( double ) ( level_max ) * level_weight_min_pos 
    - r8vec_sum ( dim_num, level_weight );
  q_max = ( double ) ( level_max ) * level_weight_min_pos;
  for ( dim = 0; dim < dim_num; dim++ )
  {
    if ( 0.0 < level_weight[dim] )
    {
      level_1d_max[dim] = r8_floor ( q_max / level_weight[dim] ) + 1;
    }
    else
    {
      level_1d_max[dim] = 0;
    }
  }
  more_grids = 0;
/*
  Repeatedly call SGMGA_VCN_ORDERED, seeking all vectors LEVEL_1D 
  which satisfy the constraint:

    LEVEL_MAX * LEVEL_WEIGHT_MIN_POS - sum ( LEVEL_WEIGHT ) 
      < sum ( 0 <= I < DIM_NUM ) LEVEL_WEIGHT[I] * LEVEL_1D[I]
      <= LEVEL_MAX * LEVEL_WEIGHT_MIN_POS.
*/
  for ( ; ; )
  {
    sgmga_vcn_ordered ( dim_num, level_weight, level_1d_max, level_1d, 
      q_min, q_max, &more_grids );

    if ( !more_grids )
    {
      break;
    }
/*
  Compute the combinatorial coefficient.
*/
    coef = sgmga_vcn_coef ( dim_num, level_weight, level_1d_max, level_1d, 
      q_min, q_max );

    if ( coef == 0.0 )
    {
      continue;
    }
/*
  Transform each 1D level to a corresponding 1D order.
*/
    level_to_order ( dim_num, level_1d, rule, order_1d );
/*
  The inner loop generates a POINT of the GRID of the LEVEL.
*/
    more_points = 0;

    for ( ; ; )
    {
      vec_colex_next3 ( dim_num, order_1d, point_index, &more_points );

      if ( !more_points )
      {
        break;
      }
      point_unique = sparse_unique_index[point_count];
      for ( dim = 0; dim < dim_num; dim++ )
      {
        sparse_order[dim+point_unique*dim_num] = order_1d[dim];
      }
      for ( dim = 0; dim < dim_num; dim++ )
      {
        sparse_index[dim+point_unique*dim_num] = point_index[dim];
      }
      point_count = point_count + 1;
    }
  }

  free ( level_1d );
  free ( level_1d_max );
  free ( order_1d );
  free ( point_index );

  return;
}
/******************************************************************************/

void sgmga_point ( int dim_num, double level_weight[], int level_max, 
  int rule[], int np[], double p[], 
  void ( *gw_compute_points[] ) ( int order, int np, double p[], double x[] ),
  int point_num, int sparse_order[], int sparse_index[], 
  void level_to_order ( int dim_num, int level[], int rule[], int order[] ),
  double sparse_point[] )

/******************************************************************************/
/*
  Purpose:

    SGMGA_POINT computes the points of an SGMGA rule.

  Discussion:

    The sparse grid is the logical sum of low degree product rules.

    Each product rule is the product of 1D factor rules.

    The user specifies:
    * the spatial dimension of the quadrature region,
    * the level that defines the Smolyak grid.
    * the quadrature rules.
    * the number of points.

    The user must preallocate space for the output array SPARSE_POINT.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 December 2009

  Author:

    John Burkardt

  Reference:

    Fabio Nobile, Raul Tempone, Clayton Webster,
    A Sparse Grid Stochastic Collocation Method for Partial Differential
    Equations with Random Input Data,
    SIAM Journal on Numerical Analysis,
    Volume 46, Number 5, 2008, pages 2309-2345.

    Fabio Nobile, Raul Tempone, Clayton Webster,
    An Anisotropic Sparse Grid Stochastic Collocation Method for Partial 
    Differential Equations with Random Input Data,
    SIAM Journal on Numerical Analysis,
    Volume 46, Number 5, 2008, pages 2411-2442.

  Parameters:

    Input, int DIM_NUM, the spatial dimension.

    Input, double LEVEL_WEIGHT[DIM_NUM], the anisotropic weights.

    Input, int LEVEL_MAX, controls the size of the final
    sparse grid.

    Input, int RULE[DIM_NUM], the rule in each dimension.
     1, "CC",  Clenshaw Curtis, Closed Fully Nested rule.
     2, "F2",  Fejer Type 2, Open Fully Nested rule.
     3, "GP",  Gauss Patterson, Open Fully Nested rule.
     4, "GL",  Gauss Legendre, Open Weakly Nested rule.
     5, "GH",  Gauss Hermite, Open Weakly Nested rule.
     6, "GGH", Generalized Gauss Hermite, Open Weakly Nested rule.
     7, "LG",  Gauss Laguerre, Open Non Nested rule.
     8, "GLG", Generalized Gauss Laguerre, Open Non Nested rule.
     9, "GJ",  Gauss Jacobi, Open Non Nested rule.
    10, "GW",  Golub Welsch, (presumed) Open Non Nested rule.
    11, "CCS", Clenshaw Curtis "Slow", Closed Fully Nested rule.
    12, "F2S", Fejer Type 2 Slow, Closed Fully Nested rule.
    13, "GPS", Gauss Patterson Slow, Closed Fully Nested rule.

    Input, int NP[DIM_NUM], the number of parameters used by each rule.

    Input, double P[sum(NP[*])], the parameters needed by each rule.

    Input, void ( *GW_COMPUTE_POINTS[] ) ( int order, int np, double p[], double x[] ),
    an array of pointers to functions which return the 1D quadrature points 
    associated with each spatial dimension for which a Golub Welsch rule 
    is used.

    Input, int POINT_NUM, the number of points in the grid,
    as determined by SGMGA_SIZE.

    Input, int SPARSE_ORDER[DIM_NUM*POINT_NUM], lists, for each point,
    the order of the 1D rules used in the grid that generated it.

    Input, int SPARSE_INDEX[DIM_NUM*POINT_NUM], lists, for each point,
    its index in each of the 1D rules in the grid that generated it.
    The indices are 1-based.

    Input, void LEVEL_TO_ORDER ( int dim_num, int level[], int rule[], 
    int order[] ), the function converting levels to orders.
    The choices are "level_to_order_default", "level_to_order_exponential" or 
    "level_to_order_linear".

    Output, double SPARSE_POINT[DIM_NUM*POINT_NUM], the points.
*/
{
  int dim;
  int level;
  int *level_1d_max;
  double level_weight_min_pos;
  int order;
  int p_index;
  int point;
  double *points;
  double q_max;

  for ( point = 0; point < point_num; point++ )
  {
    for ( dim = 0; dim < dim_num; dim++ )
    {
      sparse_point[dim+point*dim_num] = - r8_huge ( );
    }
  }
/*
  Compute the point coordinates.
*/
  level_1d_max = ( int * ) malloc ( dim_num * sizeof ( int ) );
  level_weight_min_pos = r8vec_min_pos ( dim_num, level_weight );
  q_max = ( double ) ( level_max ) * level_weight_min_pos;

  p_index = 0;

  for ( dim = 0; dim < dim_num; dim++ )
  {
    if ( 0.0 < level_weight[dim] )
    {
      level_1d_max[dim] = r8_floor ( q_max / level_weight[dim] ) + 1;
    }
    else
    {
      level_1d_max[dim] = 0;
    }

    for ( level = 0; level <= level_1d_max[dim]; level++ )
    {
      level_to_order ( 1, &level, rule+dim, &order );

      points = ( double * ) malloc ( order * sizeof ( double ) );

      if ( rule[dim] == 1 )
      {
        clenshaw_curtis_compute_points_np ( 
          order, np[dim], p+p_index, points );
      }
      else if ( rule[dim] == 2 )
      {
        fejer2_compute_points_np ( 
          order, np[dim], p+p_index, points );
      }
      else if ( rule[dim] == 3 )
      {
        patterson_lookup_points_np ( 
          order, np[dim], p+p_index, points );
      }
      else if ( rule[dim] == 4 )
      {
        legendre_compute_points_np ( 
          order, np[dim], p+p_index, points );
      }
      else if ( rule[dim] == 5 )
      {
        hermite_compute_points_np ( 
          order, np[dim], p+p_index, points );
      }
      else if ( rule[dim] == 6 )
      {
        gen_hermite_compute_points_np ( 
          order, np[dim], p+p_index, points );
      }
      else if ( rule[dim] == 7 )
      {
        laguerre_compute_points_np ( 
          order, np[dim], p+p_index, points );
      }
      else if ( rule[dim] == 8 )
      {
        gen_laguerre_compute_points_np ( 
          order, np[dim], p+p_index, points );
      }
      else if ( rule[dim] == 9 )
      {
        jacobi_compute_points_np ( 
          order, np[dim], p+p_index, points );
      }
      else if ( rule[dim] == 10 )
      {
        gw_compute_points[dim] (
          order, np[dim], p+p_index, points );
      }
      else if ( rule[dim] == 11 )
      {
        clenshaw_curtis_compute_points_np ( 
          order, np[dim], p+p_index, points );
      }
      else if ( rule[dim] == 12 )
      {
        fejer2_compute_points_np ( 
          order, np[dim], p+p_index, points );
      }
      else if ( rule[dim] == 13 )
      {
        patterson_lookup_points_np ( 
          order, np[dim], p+p_index, points );
      }
      else
      {
        printf ( "\n" );
        printf ( "SGMGA_POINT - Fatal error!\n" );
        printf ( "  Unexpected value of RULE[%d] = %d.\n", dim, rule[dim] );
        exit ( 1 );
      }

      for ( point = 0; point < point_num; point++ )
      {
        if ( sparse_order[dim+point*dim_num] == order )
        {
          sparse_point[dim+point*dim_num] = 
            points[sparse_index[dim+point*dim_num]-1];
        }
      }
      free ( points );
    }
    p_index = p_index + np[dim];
  }
/*
  Check to see if we missed any points.
*/
  for ( point = 0; point < point_num; point++ )
  {
    for ( dim = 0; dim < dim_num; dim++ )
    {
      if ( sparse_point[dim+point*dim_num] == -  r8_huge ( ) )
      {
        printf ( "\n" );
        printf ( "SGMGA_POINT - Fatal error!\n" );
        printf ( "  At least one point component was not assigned.\n" );
        printf ( "  POINT = %d\n", point );
        printf ( "  DIM = %d\n", dim );
        printf ( "  SPARSE_ORDER(DIM,POINT) = %d\n", 
          sparse_order[dim+point*dim_num] );
        printf ( "  LEVEL_WEIGHT(DIM) = %d\n", level_weight[dim] );
        exit ( 1 );
      }
    }
  }

  free ( level_1d_max );

  return;
}
/******************************************************************************/

void sgmga_product_weight ( int dim_num, int order_1d[], int order_nd, 
  int rule[], int np[], double p[], 
  void ( *gw_compute_weights[] ) ( int order, int np, double p[], double w[] ),
  double weight_nd[] )

/******************************************************************************/
/*
  Purpose:

    SGMGA_PRODUCT_WEIGHT computes the weights of a mixed product rule.

  Discussion:

    This routine computes the weights for a quadrature rule which is
    a product of 1D rules of varying order and kind.

    The user must preallocate space for the output array WEIGHT_ND.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 December 2009

  Author:

    John Burkardt

  Reference:

    Fabio Nobile, Raul Tempone, Clayton Webster,
    A Sparse Grid Stochastic Collocation Method for Partial Differential
    Equations with Random Input Data,
    SIAM Journal on Numerical Analysis,
    Volume 46, Number 5, 2008, pages 2309-2345.

    Fabio Nobile, Raul Tempone, Clayton Webster,
    An Anisotropic Sparse Grid Stochastic Collocation Method for Partial 
    Differential Equations with Random Input Data,
    SIAM Journal on Numerical Analysis,
    Volume 46, Number 5, 2008, pages 2411-2442.

  Parameters:

    Input, int DIM_NUM, the spatial dimension.

    Input, int ORDER_1D[DIM_NUM], the order of the 1D rules.

    Input, int ORDER_ND, the order of the product rule.

    Input, int RULE[DIM_NUM], the rule in each dimension.
     1, "CC",  Clenshaw Curtis, Closed Fully Nested rule.
     2, "F2",  Fejer Type 2, Open Fully Nested rule.
     3, "GP",  Gauss Patterson, Open Fully Nested rule.
     4, "GL",  Gauss Legendre, Open Weakly Nested rule.
     5, "GH",  Gauss Hermite, Open Weakly Nested rule.
     6, "GGH", Generalized Gauss Hermite, Open Weakly Nested rule.
     7, "LG",  Gauss Laguerre, Open Non Nested rule.
     8, "GLG", Generalized Gauss Laguerre, Open Non Nested rule.
     9, "GJ",  Gauss Jacobi, Open Non Nested rule.
    10, "GW",  Golub Welsch, (presumed) Open Non Nested rule.
    11, "CCS", Clenshaw Curtis "Slow", Closed Fully Nested rule.
    12, "F2S", Fejer Type 2 Slow, Closed Fully Nested rule.
    13, "GPS", Gauss Patterson Slow, Closed Fully Nested rule.

    Input, int NP[DIM_NUM], the number of parameters used by each rule.

    Input, double P[sum(NP[*])], the parameters needed by each rule.

    Input, void ( *GW_COMPUTE_WEIGHTS[] ) ( int order, int np, double p[], double w[] ),
    an array of pointers to functions which return the 1D quadrature weights 
    associated with each spatial dimension for which a Golub Welsch rule 
    is used.

    Output, double WEIGHT_ND[ORDER_ND], the product rule weights.
*/
{
  int dim;
  int i;
  int p_index;
  double *weight_1d;

  for ( i = 0; i < order_nd; i++ )
  {
    weight_nd[i] = 1.0;
  }

  p_index = 0;

  for ( dim = 0; dim < dim_num; dim++ )
  {
    weight_1d = ( double * ) malloc ( order_1d[dim] * sizeof ( double ) );

    if ( rule[dim] == 1 )
    {
      clenshaw_curtis_compute_weights_np ( 
        order_1d[dim], np[dim], p+p_index, weight_1d );
    }
    else if ( rule[dim] == 2 )
    {
      fejer2_compute_weights_np ( 
        order_1d[dim], np[dim], p+p_index, weight_1d );
    }
    else if ( rule[dim] == 3 )
    {
      patterson_lookup_weights_np ( 
        order_1d[dim], np[dim], p+p_index, weight_1d );
    }
    else if ( rule[dim] == 4 )
    {
      legendre_compute_weights_np ( 
        order_1d[dim], np[dim], p+p_index, weight_1d );
    }
    else if ( rule[dim] == 5 )
    {
      hermite_compute_weights_np ( 
        order_1d[dim], np[dim], p+p_index, weight_1d );
    }
    else if ( rule[dim] == 6 )
    {
      gen_hermite_compute_weights_np ( 
        order_1d[dim], np[dim], p+p_index, weight_1d );
    }
    else if ( rule[dim] == 7 )
    {
      laguerre_compute_weights_np ( 
        order_1d[dim], np[dim], p+p_index, weight_1d );
    }
    else if ( rule[dim] == 8 )
    {
      gen_laguerre_compute_weights_np ( 
        order_1d[dim], np[dim], p+p_index, weight_1d );
    }
    else if ( rule[dim] == 9 )
    {
      jacobi_compute_weights_np ( 
        order_1d[dim], np[dim], p+p_index, weight_1d );
    }
    else if ( rule[dim] == 10 )
    {
      gw_compute_weights[dim] ( 
        order_1d[dim], np[dim], p+p_index, weight_1d );
    }
    else if ( rule[dim] == 11 )
    {
      clenshaw_curtis_compute_weights_np ( 
        order_1d[dim], np[dim], p+p_index, weight_1d );
    }
    else if ( rule[dim] == 12 )
    {
      fejer2_compute_weights_np ( 
        order_1d[dim], np[dim], p+p_index, weight_1d );
    }
    else if ( rule[dim] == 13 )
    {
      patterson_lookup_weights_np ( 
        order_1d[dim], np[dim], p+p_index, weight_1d );
    }
    else
    {
      printf ( "\n" );
      printf ( "SGMGA_PRODUCT_WEIGHT - Fatal error!\n" );
      printf ( "  Unexpected value of RULE[%d] = %d.\n", dim, rule[dim] );
      exit ( 1 );
    }

    p_index = p_index + np[dim];

    r8vec_direct_product2 ( dim, order_1d[dim], weight_1d, 
      dim_num, order_nd, weight_nd );

    free ( weight_1d );
  }
  return;
}
/******************************************************************************/

int sgmga_size ( int dim_num, double level_weight[], int level_max, int rule[], 
  int np[], double p[], 
  void ( *gw_compute_points[] ) ( int order, int np, double p[], double x[] ),
  double tol,
  void level_to_order ( int dim_num, int level[], int rule[], int order[] ) )

/******************************************************************************/
/*
  Purpose:

    SGMGA_SIZE sizes an SGMGA grid, discounting duplicates.

  Discussion:

    The sparse grid is the logical sum of product grids that satisfy
    a particular constraint.

    Depending on the 1D rules involved, there may be many duplicate points
    in the sparse grid.

    This function counts the unique points in the sparse grid.  It does this
    in a straightforward way, by actually generating all the points, and
    comparing them, with a tolerance for equality.

    This function has been modified to automatically omit points for which
    the "combinatorial coefficient" is zero, since such points would have
    a weight of zero in the grid.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 July 2010

  Author:

    John Burkardt

  Reference:

    Fabio Nobile, Raul Tempone, Clayton Webster,
    A Sparse Grid Stochastic Collocation Method for Partial Differential
    Equations with Random Input Data,
    SIAM Journal on Numerical Analysis,
    Volume 46, Number 5, 2008, pages 2309-2345.

    Fabio Nobile, Raul Tempone, Clayton Webster,
    An Anisotropic Sparse Grid Stochastic Collocation Method for Partial 
    Differential Equations with Random Input Data,
    SIAM Journal on Numerical Analysis,
    Volume 46, Number 5, 2008, pages 2411-2442.

  Parameters:

    Input, int DIM_NUM, the spatial dimension.

    Input, double LEVEL_WEIGHT[DIM_NUM], the anisotropic weights.

    Input, int LEVEL_MAX, the maximum value of LEVEL.

    Input, int RULE[DIM_NUM], the rule in each dimension.
     1, "CC",  Clenshaw Curtis, Closed Fully Nested rule.
     2, "F2",  Fejer Type 2, Open Fully Nested rule.
     3, "GP",  Gauss Patterson, Open Fully Nested rule.
     4, "GL",  Gauss Legendre, Open Weakly Nested rule.
     5, "GH",  Gauss Hermite, Open Weakly Nested rule.
     6, "GGH", Generalized Gauss Hermite, Open Weakly Nested rule.
     7, "LG",  Gauss Laguerre, Open Non Nested rule.
     8, "GLG", Generalized Gauss Laguerre, Open Non Nested rule.
     9, "GJ",  Gauss Jacobi, Open Non Nested rule.
    10, "GW",  Golub Welsch, (presumed) Open Non Nested rule.
    11, "CCS", Clenshaw Curtis "Slow", Closed Fully Nested rule.
    12, "F2S", Fejer Type 2 Slow, Closed Fully Nested rule.
    13, "GPS", Gauss Patterson Slow, Closed Fully Nested rule.

    Input, int NP[DIM_NUM], the number of parameters used by each rule.

    Input, double P[sum(NP[*])], the parameters needed by each rule.

    Input, void ( *GW_COMPUTE_POINTS[] ) ( int order, int np, double p[], double x[] ),
    an array of pointers to functions which return the 1D quadrature points 
    associated with each spatial dimension for which a Golub Welsch rule 
    is used.

    Input, double TOL, a tolerance for point equality.

    Input, void LEVEL_TO_ORDER ( int dim_num, int level[], int rule[], 
    int order[] ), the function converting levels to orders.
    The choices are "level_to_order_default", "level_to_order_exponential" or 
    "level_to_order_linear".

    Output, int SGMGA_SIZE, the number of unique points.
*/
{
  double coef;
  int dim;
  int level;
  int *level_1d;
  int *level_1d_max;
  double level_weight_min_pos;
  int more_grids;
  int more_points;
  int order;
  int *order_1d;
  int p_index;
  int point;
  int *point_index;
  int point_num;
  int point_total_num;
  int point_total_num2;
  double *points;
  double q_max;
  double q_min;
  int seed;
  int *sparse_total_index;
  int *sparse_total_order;
  double *sparse_total_point;
/*
  Special cases.
*/
  if ( level_max < 0 )
  {
    point_num = -1;
    return point_num;
  }

  if ( level_max == 0 )
  {
    point_num = 1;
    return point_num;
  }
/*
  Get total number of points, including duplicates.
*/
  point_total_num = sgmga_size_total ( dim_num, level_weight, level_max,
    rule, level_to_order );
/*
  Generate SPARSE_TOTAL_ORDER and SPARSE_TOTAL_INDEX arrays 
  for the TOTAL set of points.
*/
  sparse_total_order = ( int * ) malloc ( dim_num * point_total_num * sizeof ( int ) );
  sparse_total_index = ( int * ) malloc ( dim_num * point_total_num * sizeof ( int ) );

  point_total_num2 = 0;

  level_1d = ( int * ) malloc ( dim_num * sizeof ( int ) );
  level_1d_max = ( int * ) malloc ( dim_num * sizeof ( int ) );
  order_1d = ( int * ) malloc ( dim_num * sizeof ( int ) );
  point_index = ( int * ) malloc ( dim_num * sizeof ( int ) );
/*
  Initialization for SGMGA_VCN_ORDERED.
*/
  level_weight_min_pos = r8vec_min_pos ( dim_num, level_weight );
  q_min = ( double ) ( level_max ) * level_weight_min_pos 
    - r8vec_sum ( dim_num, level_weight );
  q_max = ( double ) ( level_max ) * level_weight_min_pos;
  for ( dim = 0; dim < dim_num; dim++ )
  {
    if ( 0.0 < level_weight[dim] )
    {
      level_1d_max[dim] = r8_floor ( q_max / level_weight[dim] ) + 1;
    }
    else
    {
      level_1d_max[dim] = 0;
    }
  }
  more_grids = 0;
/*
  Repeatedly call SGMGA_VCN_ORDERED, seeking all vectors LEVEL_1D 
  which satisfy the constraint:

    LEVEL_MAX * LEVEL_WEIGHT_MIN_POS - sum ( LEVEL_WEIGHT ) 
      < sum ( 0 <= I < DIM_NUM ) LEVEL_WEIGHT[I] * LEVEL_1D[I]
      <= LEVEL_MAX * LEVEL_WEIGHT_MIN_POS.
*/
  for ( ; ; )
  {
    sgmga_vcn_ordered ( dim_num, level_weight, level_1d_max, level_1d, 
      q_min, q_max, &more_grids );

    if ( !more_grids )
    {
      break;
    }
/*
  Compute the combinatorial coefficient.
*/
    coef = sgmga_vcn_coef ( dim_num, level_weight, level_1d_max, level_1d, 
      q_min, q_max );

    if ( coef == 0.0 )
    {
      continue;
    }
/*
  Transform each 1D level to a corresponding 1D order.
*/
    level_to_order ( dim_num, level_1d, rule, order_1d );
/*
  The inner loop generates a POINT of the GRID of the LEVEL.
*/
    more_points = 0;

    for ( ; ; )
    {
      vec_colex_next3 ( dim_num, order_1d, point_index, &more_points );

      if ( !more_points )
      {
        break;
      }
      for ( dim = 0; dim < dim_num; dim++ )
      {
        sparse_total_order[dim+point_total_num2*dim_num] = order_1d[dim];
      }
      for ( dim = 0; dim < dim_num; dim++ )
      {
        sparse_total_index[dim+point_total_num2*dim_num] = point_index[dim];
      }
      point_total_num2 = point_total_num2 + 1;
    }
  }
  free ( level_1d );
  free ( order_1d );
  free ( point_index );
/*
  Now compute the coordinates of the TOTAL set of points.
*/
  sparse_total_point = ( double * ) malloc ( dim_num * point_total_num * sizeof ( double ) );

  for ( point = 0; point < point_total_num; point++ )
  {
    for ( dim = 0; dim < dim_num; dim++ )
    {
      sparse_total_point[dim+point*dim_num] = r8_huge ( );
    }
  }
/*
  Compute the point coordinates.
*/
  level_1d_max = ( int * ) malloc ( dim_num * sizeof ( int ) );
  level_weight_min_pos = r8vec_min_pos ( dim_num, level_weight );
  q_max = ( double ) ( level_max ) * level_weight_min_pos;

  p_index = 0;
  for ( dim = 0; dim < dim_num; dim++ )
  {
    if ( 0.0 < level_weight[dim] )
    {
      level_1d_max[dim] = r8_floor ( q_max / level_weight[dim] ) + 1;
    }
    else
    {
      level_1d_max[dim] = 0;
    }

    for ( level = 0; level <= level_1d_max[dim]; level++ )
    {
      level_to_order ( 1, &level, rule+dim, &order );

      points = ( double * ) malloc ( order * sizeof ( double ) );

      if ( rule[dim] == 1 )
      {
        clenshaw_curtis_compute_points_np ( 
          order, np[dim], p+p_index, points );
      }
      else if ( rule[dim] == 2 )
      {
        fejer2_compute_points_np ( 
          order, np[dim], p+p_index, points );
      }
      else if ( rule[dim] == 3 )
      {
        patterson_lookup_points_np ( 
          order, np[dim], p+p_index, points );
      }
      else if ( rule[dim] == 4 )
      {
        legendre_compute_points_np ( 
          order, np[dim], p+p_index, points );
      }
      else if ( rule[dim] == 5 )
      {
        hermite_compute_points_np ( 
          order, np[dim], p+p_index, points );
      }
      else if ( rule[dim] == 6 )
      {
        gen_hermite_compute_points_np ( 
          order, np[dim], p+p_index, points );
      }
      else if ( rule[dim] == 7 )
      {
        laguerre_compute_points_np ( 
          order, np[dim], p+p_index, points );
      }
      else if ( rule[dim] == 8 )
      {
        gen_laguerre_compute_points_np ( 
          order, np[dim], p+p_index, points );
      }
      else if ( rule[dim] == 9 )
      {
        jacobi_compute_points_np ( 
          order, np[dim], p+p_index, points );
      }
      else if ( rule[dim] == 10 )
      {
        gw_compute_points[dim] (
          order, np[dim], p+p_index, points );
      }
      else if ( rule[dim] == 11 )
      {
        clenshaw_curtis_compute_points_np ( 
          order, np[dim], p+p_index, points );
      }
      else if ( rule[dim] == 12 )
      {
        fejer2_compute_points_np ( 
          order, np[dim], p+p_index, points );
      }
      else if ( rule[dim] == 13 )
      {
        patterson_lookup_points_np ( 
          order, np[dim], p+p_index, points );
      }
      else
      {
        printf ( "\n" );
        printf ( "SGMGA_SIZE - Fatal error!\n" );
        printf ( "  Unexpected value of RULE[%d] = %d.\n", dim, rule[dim] );
        exit ( 1 );
      }

      for ( point = 0; point < point_total_num; point++ )
      {
        if ( sparse_total_order[dim+point*dim_num] == order )
        {
          sparse_total_point[dim+point*dim_num] = 
            points[sparse_total_index[dim+point*dim_num]-1];
        }
      }
      free ( points );
    }
    p_index = p_index + np[dim];
  }
/*
  Count the tolerably unique points. 
*/
  seed = 123456789;

  point_num = point_radial_tol_unique_count ( dim_num, point_total_num,
    sparse_total_point, tol, &seed );

  free ( level_1d_max );
  free ( sparse_total_index );
  free ( sparse_total_order );
  free ( sparse_total_point );

  return point_num;
}

/******************************************************************************/

int sgmga_size_total ( int dim_num, double level_weight[], int level_max, 
  int rule[], 
  void level_to_order ( int dim_num, int level[], int rule[], int order[] ) )

/******************************************************************************/
/*
  Purpose:

    SGMGA_SIZE_TOTAL sizes an SGMGA grid, counting duplicates.

  Discussion:

    This routine returns the total point count for an SGMGA
    ( Sparse Grid of Mixed type with Growth rule and Anisotropic weights).

    The sparse grid is the logical sum of product grids.

    The sparse grid has an associated integer index LEVEL_MAX, whose lowest 
    value is 0.  LEVEL_MAX = 0 indicates the sparse grid made up of one product 
    grid, which in turn is the product of 1D factor grids of the lowest level.
    This usually means the sparse grid with LEVEL_MAX equal to 0 is a
    one point grid.

    We can assign a level to each factor grid, and hence a LEVEL vector
    to the corresponding product grid, and a weighted index
    LEVEL_GRID (which will in general be a real number):

      LEVEL_GRID = sum ( 1 <= I <= DIM_NUM ) LEVEL_WEIGHT(I) * LEVEL(I)

    The product grid will participate in the formation of the sparse grid
    if it satisfies the following weighted constraint:

      LEVEL_MAX - DIM_NUM < LEVEL_GRID <= LEVEL_MAX

    This routine determines the total number of abscissas in all the 
    product rules used to form the SGMGA associated with the index LEVEL_MAX.
    The count disregards duplication.  If the same multidimensional abcsissa
    occurs in two different product rules that are part of the SGMGA, then
    that single abcissa is counted twice. 

    This computation is useful in cases where the entire set of abscissas
    is going to be generated, preparatory to compression to finding, indexing
    and merging the duplicate abcissass.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 November 2009

  Author:

    John Burkardt

  Reference:

    Fabio Nobile, Raul Tempone, Clayton Webster,
    A Sparse Grid Stochastic Collocation Method for Partial Differential
    Equations with Random Input Data,
    SIAM Journal on Numerical Analysis,
    Volume 46, Number 5, 2008, pages 2309-2345.

    Fabio Nobile, Raul Tempone, Clayton Webster,
    An Anisotropic Sparse Grid Stochastic Collocation Method for Partial 
    Differential Equations with Random Input Data,
    SIAM Journal on Numerical Analysis,
    Volume 46, Number 5, 2008, pages 2411-2442.

  Parameters:

    Input, int DIM_NUM, the spatial dimension.

    Input, double LEVEL_WEIGHT[DIM_NUM], the anisotropic weights.

    Input, int LEVEL_MAX, the maximum value of LEVEL.

    Input, int RULE[DIM_NUM], the rule in each dimension.
     1, "CC",  Clenshaw Curtis, Closed Fully Nested rule.
     2, "F2",  Fejer Type 2, Open Fully Nested rule.
     3, "GP",  Gauss Patterson, Open Fully Nested rule.
     4, "GL",  Gauss Legendre, Open Weakly Nested rule.
     5, "GH",  Gauss Hermite, Open Weakly Nested rule.
     6, "GGH", Generalized Gauss Hermite, Open Weakly Nested rule.
     7, "LG",  Gauss Laguerre, Open Non Nested rule.
     8, "GLG", Generalized Gauss Laguerre, Open Non Nested rule.
     9, "GJ",  Gauss Jacobi, Open Non Nested rule.
    10, "GW",  Golub Welsch, (presumed) Open Non Nested rule.
    11, "CCS", Clenshaw Curtis "Slow", Closed Fully Nested rule.
    12, "F2S", Fejer Type 2 Slow, Closed Fully Nested rule.
    13, "GPS", Gauss Patterson Slow, Closed Fully Nested rule.

    Input, void LEVEL_TO_ORDER ( int dim_num, int level[], int rule[], 
    int order[] ), the function converting levels to orders.
    The choices are "level_to_order_default", "level_to_order_exponential" or 
    "level_to_order_linear".

    Output, int SGMGA_SIZE_TOTAL, the number of points
    including repetitions.
*/
{
  double coef;
  int dim;
  int *level_1d;
  int *level_1d_max;
  double level_weight_min_pos;
  int more_grids;
  int *order_1d;
  int point_total_num;
  double q_max;
  double q_min;
/*
  Special case.
*/
  if ( level_max == 0 )
  {
    point_total_num = 1;
    return point_total_num;
  }

  point_total_num = 0;

  level_1d = ( int * ) malloc ( dim_num * sizeof ( int ) );
  level_1d_max = ( int * ) malloc ( dim_num * sizeof ( int ) );
  order_1d = ( int * ) malloc ( dim_num * sizeof ( int ) );
/*
  Initialization for SGMGA_VCN_ORDERED.
*/
  level_weight_min_pos = r8vec_min_pos ( dim_num, level_weight );
  q_min = ( double ) ( level_max ) * level_weight_min_pos 
    - r8vec_sum ( dim_num, level_weight );
  q_max = ( double ) ( level_max ) * level_weight_min_pos;
  for ( dim = 0; dim < dim_num; dim++ )
  {
    if ( 0.0 < level_weight[dim] )
    {
      level_1d_max[dim] = r8_floor ( q_max / level_weight[dim] ) + 1;
    }
    else
    {
      level_1d_max[dim] = 0;
    }
  }
  more_grids = 0;
/*
  Repeatedly call SGMGA_VCN_ORDERED, seeking all vectors LEVEL_1D 
  which satisfy the constraint:

    LEVEL_MAX * LEVEL_WEIGHT_MIN_POS - sum ( LEVEL_WEIGHT ) 
      < sum ( 0 <= I < DIM_NUM ) LEVEL_WEIGHT[I] * LEVEL_1D[I]
      <= LEVEL_MAX * LEVEL_WEIGHT_MIN_POS.
*/
  for ( ; ; )
  {
    sgmga_vcn_ordered ( dim_num, level_weight, level_1d_max, level_1d, 
      q_min, q_max, &more_grids );

    if ( !more_grids )
    {
      break;
    }
/*
  Compute the combinatorlal coefficient.
*/
    coef = sgmga_vcn_coef ( dim_num, level_weight, level_1d_max, level_1d, 
      q_min, q_max );

    if ( coef == 0.0 )
    {
      continue;
    }
/*
  Transform each 1D level to a corresponding 1D order.
*/
    level_to_order ( dim_num, level_1d, rule, order_1d );

    point_total_num = point_total_num + i4vec_product ( dim_num, 
      order_1d );
  }
  free ( level_1d );
  free ( level_1d_max );
  free ( order_1d );

  return point_total_num;
}
/******************************************************************************/

void sgmga_unique_index ( int dim_num, double level_weight[], int level_max, 
  int rule[], int np[], double p[], 
  void ( *gw_compute_points[] ) ( int order, int np, double p[], double x[] ),
  double tol, int point_num, int point_total_num, 
  void level_to_order ( int dim_num, int level[], int rule[], int order[] ), 
  int sparse_unique_index[] )

/******************************************************************************/
/*
  Purpose:

    SGMGA_UNIQUE_INDEX maps nonunique to unique points.

  Discussion:

    The sparse grid usually contains many points that occur in more
    than one product grid.

    When generating the point locations, it is easy to realize that a point
    has already been generated.

    But when it's time to compute the weights of the sparse grids, it is
    necessary to handle situations in which weights corresponding to 
    the same point generated in multiple grids must be collected together.

    This routine generates ALL the points, including their multiplicities,
    and figures out a mapping from them to the collapsed set of unique points.

    This mapping can then be used during the weight calculation so that
    a contribution to the weight gets to the right place.

    The user must preallocate space for the output array SPARSE_UNIQUE_INDEX.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    29 July 2010

  Author:

    John Burkardt

  Reference:

    Fabio Nobile, Raul Tempone, Clayton Webster,
    A Sparse Grid Stochastic Collocation Method for Partial Differential
    Equations with Random Input Data,
    SIAM Journal on Numerical Analysis,
    Volume 46, Number 5, 2008, pages 2309-2345.

    Fabio Nobile, Raul Tempone, Clayton Webster,
    An Anisotropic Sparse Grid Stochastic Collocation Method for Partial 
    Differential Equations with Random Input Data,
    SIAM Journal on Numerical Analysis,
    Volume 46, Number 5, 2008, pages 2411-2442.

  Parameters:

    Input, int DIM_NUM, the spatial dimension.

    Input, int LEVEL_MAX, the maximum value of LEVEL.

    Input, double LEVEL_WEIGHT[DIM_NUM], the anisotropic weights.

    Input, int RULE[DIM_NUM], the rule in each dimension.
     1, "CC",  Clenshaw Curtis, Closed Fully Nested rule.
     2, "F2",  Fejer Type 2, Open Fully Nested rule.
     3, "GP",  Gauss Patterson, Open Fully Nested rule.
     4, "GL",  Gauss Legendre, Open Weakly Nested rule.
     5, "GH",  Gauss Hermite, Open Weakly Nested rule.
     6, "GGH", Generalized Gauss Hermite, Open Weakly Nested rule.
     7, "LG",  Gauss Laguerre, Open Non Nested rule.
     8, "GLG", Generalized Gauss Laguerre, Open Non Nested rule.
     9, "GJ",  Gauss Jacobi, Open Non Nested rule.
    10, "GW",  Golub Welsch, (presumed) Open Non Nested rule.
    11, "CCS", Clenshaw Curtis "Slow", Closed Fully Nested rule.
    12, "F2S", Fejer Type 2 Slow, Closed Fully Nested rule.
    13, "GPS", Gauss Patterson Slow, Closed Fully Nested rule.

    Input, int NP[DIM_NUM], the number of parameters used by each rule.

    Input, double P[sum(NP[*])], the parameters needed by each rule.

    Input, void ( *GW_COMPUTE_POINTS[] ) ( int order, int np, double p[], double x[] ),
    an array of pointers to functions which return the 1D quadrature points 
    associated with each spatial dimension for which a Golub Welsch rule 
    is used.

    Input, double TOL, a tolerance for point equality.

    Input, int POINT_NUM, the number of unique points 
    in the grid. 

    Input, int POINT_TOTAL_NUM, the total number of points 
    in the grid. 

    Input, void LEVEL_TO_ORDER ( int dim_num, int level[], int rule[], 
    int order[] ), the function converting levels to orders.
    The choices are "level_to_order_default", "level_to_order_exponential" or 
    "level_to_order_linear".

    Output, int SPARSE UNIQUE_INDEX[POINT_TOTAL_NUM], lists, 
    for each (nonunique) point, the corresponding index of the same point in 
    the unique listing.
*/
{
  double coef;
  int dim;
  int level;
  int *level_1d;
  int *level_1d_max;
  double level_weight_min_pos;
  int more_grids;
  int more_points;
  int order;
  int *order_1d;
  int p_index;
  int point;
  int *point_index;
  int point_num2;
  int point_total_num2;
  double *points;
  double q_max;
  double q_min;
  int rep;
  int seed;
  int *sparse_total_index;
  int *sparse_total_order;
  double *sparse_total_point;
  int *undx;
/*
  Special cases.
*/
  if ( level_max < 0 )
  {
    return;
  }

  if ( level_max == 0 )
  {
    sparse_unique_index[0] = 0;
    return;
  }
/*
  Generate SPARSE_TOTAL_ORDER and SPARSE_TOTAL_INDEX arrays 
  for the TOTAL set of points.
*/
  sparse_total_order = ( int * ) malloc ( dim_num * point_total_num * sizeof ( int ) );
  sparse_total_index = ( int * ) malloc ( dim_num * point_total_num * sizeof ( int ) );

  level_1d = ( int * ) malloc ( dim_num * sizeof ( int ) );
  level_1d_max = ( int * ) malloc ( dim_num * sizeof ( int ) );
  order_1d = ( int * ) malloc ( dim_num * sizeof ( int ) );
  point_index = ( int * ) malloc ( dim_num * sizeof ( int ) );

  point_total_num2 = 0;
/*
  Initialization for SGMGA_VCN_ORDERED.
*/
  level_weight_min_pos = r8vec_min_pos ( dim_num, level_weight );
  q_min = ( double ) ( level_max ) * level_weight_min_pos 
    - r8vec_sum ( dim_num, level_weight );
  q_max = ( double ) ( level_max ) * level_weight_min_pos;
  for ( dim = 0; dim < dim_num; dim++ )
  {
    if ( 0.0 < level_weight[dim] )
    {
      level_1d_max[dim] = r8_floor ( q_max / level_weight[dim] ) + 1;
    }
    else
    {
      level_1d_max[dim] = 0;
    }
  }
  more_grids = 0;
/*
  Repeatedly call SGMGA_VCN_ORDERED, seeking all vectors LEVEL_1D 
  which satisfy the constraint:

    LEVEL_MAX * LEVEL_WEIGHT_MIN_POS - sum ( LEVEL_WEIGHT ) 
      < sum ( 0 <= I < DIM_NUM ) LEVEL_WEIGHT[I] * LEVEL_1D[I]
      <= LEVEL_MAX * LEVEL_WEIGHT_MIN_POS.
*/
  for ( ; ; )
  {
    sgmga_vcn_ordered ( dim_num, level_weight, level_1d_max, level_1d, 
      q_min, q_max, &more_grids );

    if ( !more_grids )
    {
      break;
    }
/*
  Compute the combinatorial coefficient.
*/
    coef = sgmga_vcn_coef ( dim_num, level_weight, level_1d_max, level_1d, 
      q_min, q_max );

    if ( coef == 0.0 )
    {
      continue;
    }
/*
  Transform each 1D level to a corresponding 1D order.
*/
    level_to_order ( dim_num, level_1d, rule, order_1d );
/*
  The inner loop generates a POINT of the GRID of the LEVEL.
*/
    more_points = 0;

    for ( ; ; )
    {
      vec_colex_next3 ( dim_num, order_1d, point_index, &more_points );

      if ( !more_points )
      {
        break;
      }
      for ( dim = 0; dim < dim_num; dim++ )
      {
        sparse_total_order[dim+point_total_num2*dim_num] = order_1d[dim];
      }
      for ( dim = 0; dim < dim_num; dim++ )
      {
        sparse_total_index[dim+point_total_num2*dim_num] = point_index[dim];
      }
      point_total_num2 = point_total_num2 + 1;
    }
  }
  free ( level_1d );
  free ( level_1d_max );
  free ( order_1d );
  free ( point_index );
/*
  Now compute the coordinates of the TOTAL set of points.
*/
  sparse_total_point = ( double * ) 
    malloc ( dim_num * point_total_num * sizeof ( double ) );

  for ( point = 0; point < point_total_num; point++ )
  {
    for ( dim = 0; dim < dim_num; dim++ )
    {
      sparse_total_point[dim+point*dim_num] = r8_huge ( );
    }
  }
/*
  Compute the point coordinates.
*/
  level_1d_max = ( int * ) malloc ( dim_num * sizeof ( int ) );
  level_weight_min_pos = r8vec_min_pos ( dim_num, level_weight );
  q_max = ( double ) ( level_max ) * level_weight_min_pos;

  p_index = 0;
  for ( dim = 0; dim < dim_num; dim++ )
  {
    if ( 0.0 < level_weight[dim] )
    {
      level_1d_max[dim] = r8_floor ( q_max / level_weight[dim] ) + 1;
    }
    else
    {
      level_1d_max[dim] = 0;
    }

    for ( level = 0; level <= level_1d_max[dim]; level++ )
    {
      level_to_order ( 1, &level, rule+dim, &order );

      points = ( double * ) malloc ( order * sizeof ( double ) );

      if ( rule[dim] == 1 )
      {
        clenshaw_curtis_compute_points_np ( 
          order, np[dim], p+p_index, points );
      }
      else if ( rule[dim] == 2 )
      {
        fejer2_compute_points_np ( 
          order, np[dim], p+p_index, points );
      }
      else if ( rule[dim] == 3 )
      {
        patterson_lookup_points_np ( 
          order, np[dim], p+p_index, points );
      }
      else if ( rule[dim] == 4 )
      {
        legendre_compute_points_np ( 
          order, np[dim], p+p_index, points );
      }
      else if ( rule[dim] == 5 )
      {
        hermite_compute_points_np ( 
          order, np[dim], p+p_index, points );
      }
      else if ( rule[dim] == 6 )
      {
        gen_hermite_compute_points_np ( 
          order, np[dim], p+p_index, points );
      }
      else if ( rule[dim] == 7 )
      {
        laguerre_compute_points_np ( 
          order, np[dim], p+p_index, points );
      }
      else if ( rule[dim] == 8 )
      {
        gen_laguerre_compute_points_np ( 
          order, np[dim], p+p_index, points );
      }
      else if ( rule[dim] == 9 )
      {
        jacobi_compute_points_np ( 
          order, np[dim], p+p_index, points );
      }
      else if ( rule[dim] == 10 )
      {
        gw_compute_points[dim] (
          order, np[dim], p+p_index, points );
      }
      else if ( rule[dim] == 11 )
      {
        clenshaw_curtis_compute_points_np ( 
          order, np[dim], p+p_index, points );
      }
      else if ( rule[dim] == 12 )
      {
        fejer2_compute_points_np ( 
          order, np[dim], p+p_index, points );
      }
      else if ( rule[dim] == 13 )
      {
        patterson_lookup_points_np ( 
          order, np[dim], p+p_index, points );
      }
      else
      {
        printf ( "\n" );
        printf ( "SGMGA_UNIQUE_INDEX - Fatal error!\n" );
        printf ( "  Unexpected value of RULE[%d] = %d.\n", dim, rule[dim] );
        exit ( 1 );
      }

      for ( point = 0; point < point_total_num; point++ )
      {
        if ( sparse_total_order[dim+point*dim_num] == order )
        {
          sparse_total_point[dim+point*dim_num] = 
            points[sparse_total_index[dim+point*dim_num]-1];
        }
      }
      free ( points );
    }
    p_index = p_index + np[dim];
  }
/*
  Merge points that are too close.
*/
  seed = 123456789;

  undx = ( int * ) malloc ( point_num * sizeof ( int ) ) ;

  point_num2 = point_radial_tol_unique_index ( dim_num, point_total_num, 
    sparse_total_point, tol, &seed, undx, sparse_unique_index );

  for ( point = 0; point < point_total_num; point++ )
  {
    rep = undx[sparse_unique_index[point]];
    if ( point != rep )
    {
      for ( dim = 0; dim < dim_num; dim++ )
      {
        sparse_total_point[dim+point*dim_num] = sparse_total_point[dim+rep*dim_num];
      }
    }
  }
/*
  Construct an index that indicates the "rank" of the unique points.
*/
  point_unique_index ( dim_num, point_total_num, sparse_total_point,
    point_num, undx, sparse_unique_index );

  free ( undx );

  free ( sparse_total_index );
  free ( sparse_total_order );
  free ( sparse_total_point );

  return;
}
/******************************************************************************/

void sgmga_vcn ( int dim_num, double level_weight[], int x_max[], int x[], 
  double q_min, double q_max, int *more )

/******************************************************************************/
/*
  Purpose:

    SGMGA_VCN returns the next constrained vector.

  Discussion:

    We consider vectors X of dimension DIM_NUM satisfying:

      0 <= X(1:DIM_NUM) <= X_MAX(1:DIM_NUM).

    and define

      Q = sum ( 1 <= I <= DIM_NUM ) LEVEL_WEIGHT(I) * X(I)

    and seek X satisfying the constraint:

      Q_MIN < Q <= Q_MAX

    For sparse grid applications, we compute

      LEVEL_WEIGHT_MIN_POS = minimum positive entry in LEVEL_WEIGHT

    and assume there is an underlying LEVEL used to index the sets of 
    constrained vectors, and that 

      Q_MAX = LEVEL * LEVEL_WEIGHT_MIN_POS
      Q_MIN = LEVEL - LEVEL_WEIGHT_MIN_POS * sum ( LEVEL_WEIGHT(:) )
      X_MAX(I) = LEVEL * LEVEL_WEIGHT_MIN_POS / LEVEL_WEIGHT(I)

    This routine returns, one at a time exactly those X which satisfy
    the constraint.  No attempt is made to return the X values in 
    any particular order as far as Q goes.  

  Example:

    LEVEL_WEIGHT:          1.000000        1.000000

    Q_MIN:        0.000000
    Q_MAX:        2.000000
    X_MAX:                         2         2

         1        1.000000         1         0
         2        2.000000         2         0
         3        1.000000         0         1
         4        2.000000         1         1
         5        2.000000         0         2

    LEVEL_WEIGHT:          1.000000        2.000000

    Q_MIN:       -1.000000
    Q_MAX:        2.000000
    X_MAX:                         2         1

         1        0.000000         0         0
         2        1.000000         1         0
         3        2.000000         2         0
         4        2.000000         0         1

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 November 2009

  Author:

    John Burkardt

  Reference:

    Fabio Nobile, Raul Tempone, Clayton Webster,
    A Sparse Grid Stochastic Collocation Method for Partial Differential
    Equations with Random Input Data,
    SIAM Journal on Numerical Analysis,
    Volume 46, Number 5, 2008, pages 2309-2345.

    Fabio Nobile, Raul Tempone, Clayton Webster,
    An Anisotropic Sparse Grid Stochastic Collocation Method for Partial 
    Differential Equations with Random Input Data,
    SIAM Journal on Numerical Analysis,
    Volume 46, Number 5, 2008, pages 2411-2442.

  Parameters:

    Input, int DIM_NUM, the number of components in the vector.

    Input, double LEVEL_WEIGHT[DIM_NUM], the anisotropic weights.

    Input, int X_MAX[DIM_NUM], the maximum values allowed in each component.

    Input/output, int X[DIM_NUM].  On first call (with MORE = FALSE),
    the input value of X is not important.  On subsequent calls, the
    input value of X should be the output value from the previous call.
    On output, (with MORE = TRUE), the value of X will be the "next"
    vector in the reverse lexicographical list of vectors that satisfy
    the condition.  However, on output with MORE = FALSE, the vector
    X is meaningless, because there are no more vectors in the list.

    Input, double Q_MIN, Q_MAX, the lower and upper
    limits on the sum.

    Input/output, int *MORE.  On input, if the user has set MORE
    FALSE, the user is requesting the initiation of a new sequence
    of values.  If MORE is TRUE, then the user is requesting "more"
    values in the current sequence.  On output, if MORE is TRUE,
    then another value was found and returned in X, but if MORE is
    FALSE, then there are no more values in the sequence, and X is
    NOT the next value.
*/
{
  int i;
  int j;
  double q;

  if ( ! ( *more ) )
  {
    *more = 1;
    for ( i = 0; i < dim_num; i++ )
    {
      x[i] = 0;
    }

    q = 0.0;
    for ( i = 0; i < dim_num; i++ )
    {
      q = q + level_weight[i] * ( double ) ( x[i] );
    }

    if ( q_min < q && q <= q_max )
    {
      return;
    }
  }

  for ( ; ; )
  {
    j = 0;

    for ( ; ; )
    {
      if ( x[j] < x_max[j] )
      {
        break;
      }

      if ( dim_num - 1 <= j )
      {
        *more = 0;
        return;
      }
      j = j + 1;
    }

    x[j] = x[j] + 1;
    for ( i = 0; i < j; i++ )
    {
      x[i] = 0;
    }

    q = 0.0;
    for ( i = 0; i < dim_num; i++ )
    {
      q = q + level_weight[i] * ( double ) ( x[i] );
    }

    if ( q_min < q && q <= q_max )
    {
      break;
    }
  }

  return;
}
/******************************************************************************/

double sgmga_vcn_coef ( int dim_num, double level_weight[], int x_max[], 
  int x[], double q_min, double q_max )

/******************************************************************************/
/*
  Purpose:

    SGMGA_VCN_COEF returns the "next" constrained vector's coefficient.

  Discussion:

    We are given a vector X of dimension DIM_NUM which satisfies:

      0 <= X(1:DIM_NUM) <= X_MAX(1:DIM_NUM).

    and the following constraint:

      Q_MIN < sum ( 1 <= I <= DIM_NUM ) LEVEL_WEIGHT(I) * X(I) <= Q_MAX

    This routine computes the appropriate coefficient for X in the
    anisotropic sparse grid scheme.

    The coefficient is calculated as follows:

      Let B be a binary vector of length DIM_NUM, and let ||B|| represent
      the sum of the entries of B.

      Coef = sum ( all B such that X+B satisfies constraints ) (-1)^||B||

    Since X+0 satisfies the constraint, there is always at least one 
    summand.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 November 2009

  Author:

    John Burkardt

  Reference:

    Fabio Nobile, Raul Tempone, Clayton Webster,
    A Sparse Grid Stochastic Collocation Method for Partial Differential
    Equations with Random Input Data,
    SIAM Journal on Numerical Analysis,
    Volume 46, Number 5, 2008, pages 2309-2345.

    Fabio Nobile, Raul Tempone, Clayton Webster,
    An Anisotropic Sparse Grid Stochastic Collocation Method for Partial 
    Differential Equations with Random Input Data,
    SIAM Journal on Numerical Analysis,
    Volume 46, Number 5, 2008, pages 2411-2442.

  Parameters:

    Input, int DIM_NUM, the number of components in the vector.

    Input, double LEVEL_WEIGHT[DIM_NUM], the anisotropic weights.

    Input, int X_MAX[DIM_NUM], the maximum
    values allowed in each component.

    Input, int X[DIM_NUM], a point which satisifies the constraints.

    Input, double Q_MIN, Q_MAX, the lower and upper
    limits on the sum.

    Output, double SGMGA_VCN_COEF, the combinatorial coefficient.
*/
{
  int *b;
  int b_sum;
  double coef;
  int i;
  int legal;
  double q;
  int *x2;

  b = ( int * ) malloc ( dim_num * sizeof ( int ) );
  x2 = ( int * ) malloc ( dim_num * sizeof ( int ) );

  for ( i = 0; i < dim_num; i++ )
  {
    b[i] = 0;
  }
  coef = 1.0;

  for ( ; ; )
  {
/*
  Generate the next binary perturbation.
*/
    binary_vector_next ( dim_num, b );
    b_sum = i4vec_sum ( dim_num, b );
/*
  We're done if we've got back to 0.
*/
    if ( b_sum == 0 )
    {
      break;
    }
/*
  Perturb the vector.
*/
    for ( i = 0; i < dim_num; i++ )
    {
      x2[i] = x[i] + b[i];
    }
/*
  Does it satisfy the XMAX constraint?
*/
    legal = 1;
    for ( i = 0; i < dim_num; i++ )
    {
      if ( x_max[i] < x2[i] )
      {
        legal = 0;
        break;
      }
    }
    if ( !legal )
    {
      continue;
    }
/*
  Does it satisfy the Q_MIN, Q_MAX constraint?
*/
    q = 0.0;
    for ( i = 0; i < dim_num; i++ )
    {
      q = q + level_weight[i] * ( double ) ( x2[i] );
    }

    if ( q_min < q && q <= q_max )
    {
      coef = coef + r8_mop ( b_sum );
    }
  }

  free ( b );
  free ( x2 );

  return coef;
}
/******************************************************************************/

void sgmga_vcn_ordered ( int dim_num, double level_weight[], int x_max[], 
  int x[], double q_min, double q_max, int *more )

/******************************************************************************/
/*
  Purpose:

    SGMGA_VCN_ORDERED returns the "next" constrained vector, with ordering.

  Discussion:

    We consider vectors X of dimension DIM_NUM satisfying:

      0 <= X(1:DIM_NUM) <= X_MAX(1:DIM_NUM).

    and define

      Q = sum ( 1 <= I <= DIM_NUM ) LEVEL_WEIGHT(I) * X(I)

    and seek X's satisfying the constraint:

      Q_MIN < Q <= Q_MAX

    For sparse grid applications, we compute

      LEVEL_WEIGHT_MIN_POS = minimum positive entry in LEVEL_WEIGHT

    and assume there is an underlying LEVEL used to index the sets of 
    constrained vectors, and that 

      Q_MAX = LEVEL * LEVEL_WEIGHT_MIN_POS
      Q_MIN = LEVEL - LEVEL_WEIGHT_MIN_POS * sum ( LEVEL_WEIGHT(:) )
      X_MAX(I) = LEVEL * LEVEL_WEIGHT_MIN_POS / LEVEL_WEIGHT(I)

    This function returns, one at a time exactly those X which satisfy
    the constraint.

    A weak ordering is imposed on the solution vectors.  This function 
    subdivides the range Q_MIN through Q_MAX into subintervals of width 1, so 
    that the X vectors returned are roughly sorted (or at least binned) 
    by Q value.

  Example:

    If the weights are also integral, then the X vectors are in fact SORTED 
    by Q value:

    LEVEL_WEIGHT:          1.000000        1.000000
    Q_MIN:        0.000000
    Q_MAX:        2.000000
    X_MAX:                         2         2

         1        1.000000         1         0
         2        1.000000         0         1
         3        2.000000         2         0
         4        2.000000         1         1
         5        2.000000         0         2

    When the weights are not integral, then the X values are only BINNED
    by Q value, that is, we first get all X's with Q values between Q_MIN
    and Q_MIN+1, then Q_MIN+1 to Q_MIN+2 and so on, as demonstrated here:

    LEVEL_WEIGHT:             1.5               1
    Q_MIN:  0.5
    Q_MAX:  3
    X_MAX:                           2         3

           1             1.5         1         0
           2               1         0         1
           3             2.5         1         1
           4               2         0         2
           5               3         2         0
           6               3         0         3

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 November 2009

  Author:

    John Burkardt

  Reference:

    Fabio Nobile, Raul Tempone, Clayton Webster,
    A Sparse Grid Stochastic Collocation Method for Partial Differential
    Equations with Random Input Data,
    SIAM Journal on Numerical Analysis,
    Volume 46, Number 5, 2008, pages 2309-2345.

    Fabio Nobile, Raul Tempone, Clayton Webster,
    An Anisotropic Sparse Grid Stochastic Collocation Method for Partial 
    Differential Equations with Random Input Data,
    SIAM Journal on Numerical Analysis,
    Volume 46, Number 5, 2008, pages 2411-2442.

  Parameters:

    Input, int DIM_NUM, the number of components in the vector.

    Input, double LEVEL_WEIGHT[DIM_NUM], the anisotropic weights.

    Input, int X_MAX[DIM_NUM], the maximum values allowed in each component.

    Input/output, int X[DIM_NUM].  On first call (with MORE = FALSE),
    the input value of X is not important.  On subsequent calls, the
    input value of X should be the output value from the previous call.
    On output, (with MORE = TRUE), the value of X will be the "next"
    vector in the reverse lexicographical list of vectors that satisfy
    the condition.  However, on output with MORE = FALSE, the vector
    X is meaningless, because there are no more vectors in the list.

    Input, double Q_MIN, Q_MAX, the lower and upper
    limits on the sum.

    Input/output, int *MORE.  On input, if the user has set MORE
    FALSE, the user is requesting the initiation of a new sequence
    of values.  If MORE is TRUE, then the user is requesting "more"
    values in the current sequence.  On output, if MORE is TRUE,
    then another value was found and returned in X, but if MORE is
    FALSE, then there are no more values in the sequence, and X is
    NOT the next value.
*/
{
  double q;
  static double q_max2;
  static double q_min2;
/*
  On first call, initialize the subrange.
*/
  if ( !(*more) )
  {
    q_min2 = q_min;
    q_max2 = r8_min ( q_min + 1.0, q_max );
  }
/*
  Call a lower level function to search the subrange.
*/
  for ( ; ; )
  {
    sgmga_vcn ( dim_num, level_weight, x_max, x, q_min2, q_max2, more );
/*
  If another solution was found, return it.
*/
    if ( *more )
    {
      return;
    }
/*
  If the current subrange is exhausted, try to move to the next one.
*/
    if ( q_max2 < q_max )
    {
      q_min2 = q_max2;
      q_max2 = r8_min ( q_max2 + 1.0, q_max );
    }
/*
  If there are no more subranges, we're done.
*/
    else
    {
      break;
    }
  }

  return;
}
/******************************************************************************/

void sgmga_weight ( int dim_num, double level_weight[], int level_max, 
  int rule[], int np[], double p[], 
  void ( *gw_compute_weights[] ) ( int order, int np, double p[], double w[] ),
  int point_num, int point_total_num, int sparse_unique_index[], 
  void level_to_order ( int dim_num, int level[], int rule[], int order[] ),
  double sparse_weight[] )

/******************************************************************************/
/*
  Purpose:

    SGMGA_WEIGHT computes weights for an SGMGA grid.

  Discussion:

    The user must preallocate space for the output array SPARSE_WEIGHT.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 November 2009

  Author:

    John Burkardt

  Reference:

    Fabio Nobile, Raul Tempone, Clayton Webster,
    A Sparse Grid Stochastic Collocation Method for Partial Differential
    Equations with Random Input Data,
    SIAM Journal on Numerical Analysis,
    Volume 46, Number 5, 2008, pages 2309-2345.

    Fabio Nobile, Raul Tempone, Clayton Webster,
    An Anisotropic Sparse Grid Stochastic Collocation Method for Partial 
    Differential Equations with Random Input Data,
    SIAM Journal on Numerical Analysis,
    Volume 46, Number 5, 2008, pages 2411-2442.

  Parameters:

    Input, int DIM_NUM, the spatial dimension.

    Input, double LEVEL_WEIGHT[DIM_NUM], the anisotropic weights.

    Input, int LEVEL_MAX, the maximum value of LEVEL.

    Input, int RULE[DIM_NUM], the rule in each dimension.
     1, "CC",  Clenshaw Curtis, Closed Fully Nested rule.
     2, "F2",  Fejer Type 2, Open Fully Nested rule.
     3, "GP",  Gauss Patterson, Open Fully Nested rule.
     4, "GL",  Gauss Legendre, Open Weakly Nested rule.
     5, "GH",  Gauss Hermite, Open Weakly Nested rule.
     6, "GGH", Generalized Gauss Hermite, Open Weakly Nested rule.
     7, "LG",  Gauss Laguerre, Open Non Nested rule.
     8, "GLG", Generalized Gauss Laguerre, Open Non Nested rule.
     9, "GJ",  Gauss Jacobi, Open Non Nested rule.
    10, "GW",  Golub Welsch, (presumed) Open Non Nested rule.
    11, "CCS", Clenshaw Curtis "Slow", Closed Fully Nested rule.
    12, "F2S", Fejer Type 2 Slow, Closed Fully Nested rule.
    13, "GPS", Gauss Patterson Slow, Closed Fully Nested rule.

    Input, int NP[DIM_NUM], the number of parameters used by each rule.

    Input, double P[sum(NP[*])], the parameters needed by each rule.

    Input, void ( *GW_COMPUTE_WEIGHTS[] ) ( int order, int np, double p[], double w[] ),
    an array of pointers to functions which return the 1D quadrature weights 
    associated with each spatial dimension for which a Golub Welsch rule 
    is used.

    Input, int POINT_NUM, the number of unique points 
    in the grid. 

    Input, int POINT_TOTAL_NUM, the total number of points 
    in the grid. 

    Input, int SPARSE UNIQUE_INDEX[POINT_TOTAL_NUM], lists, 
    for each (nonunique) point, the corresponding index of the same point in 
    the unique listing.

    Input, void LEVEL_TO_ORDER ( int dim_num, int level[], int rule[], 
    int order[] ), the function converting levels to orders.
    The choices are "level_to_order_default", "level_to_order_exponential" or 
    "level_to_order_linear".

    Output, double SPARSE_WEIGHT[POINT_NUM], the weights
    associated with the sparse grid points.
*/
{
  double coef;
  int dim;
  double *grid_weight;
  int level;
  int *level_1d;
  int *level_1d_max;
  double level_weight_min_pos;
  int more_grids;
  int order;
  int *order_1d;
  int order_nd;
  int point;
  int point_total;
  int point_unique;
  double q_max;
  double q_min;

  for ( point = 0; point < point_num; point++ )
  {
    sparse_weight[point] = 0.0;
  }

  point_total = 0;

  level_1d = ( int * ) malloc ( dim_num * sizeof ( int ) );
  order_1d = ( int * ) malloc ( dim_num * sizeof ( int ) );
  level_1d_max = ( int * ) malloc ( dim_num * sizeof ( int ) );
/*
  Initialization for SGMGA_VCN_ORDERED.
*/
  level_weight_min_pos = r8vec_min_pos ( dim_num, level_weight );
  q_min = ( double ) ( level_max ) * level_weight_min_pos 
    - r8vec_sum ( dim_num, level_weight );
  q_max = ( double ) ( level_max ) * level_weight_min_pos;
  for ( dim = 0; dim < dim_num; dim++ )
  {
    if ( 0.0 < level_weight[dim] )
    {
      level_1d_max[dim] = r8_floor ( q_max / level_weight[dim] ) + 1;
    }
    else
    {
      level_1d_max[dim] = 0;
    }
  }
  more_grids = 0;
/*
  Repeatedly call SGMGA_VCN_ORDERED, seeking all vectors LEVEL_1D 
  which satisfy the constraint:

    LEVEL_MAX * LEVEL_WEIGHT_MIN_POS - sum ( LEVEL_WEIGHT ) 
      < sum ( 0 <= I < DIM_NUM ) LEVEL_WEIGHT[I] * LEVEL_1D[I]
      <= LEVEL_MAX * LEVEL_WEIGHT_MIN_POS.
*/
  for ( ; ; )
  {
    sgmga_vcn_ordered ( dim_num, level_weight, level_1d_max, level_1d, 
      q_min, q_max, &more_grids );

    if ( !more_grids )
    {
      break;
    }
/*
  Compute the combinatorial coefficient.
*/
    coef = sgmga_vcn_coef ( dim_num, level_weight, level_1d_max, level_1d, 
      q_min, q_max );

    if ( coef == 0.0 )
    {
      continue;
    }
/*
  Transform each 1D level to a corresponding 1D order.
*/
    level_to_order ( dim_num, level_1d, rule, order_1d );
/*
  The product of the 1D orders gives us the number of points in this grid.
*/
    order_nd = i4vec_product ( dim_num, order_1d );
/*
  Compute the weights for this grid.

  The correct transfer of data from the product grid to the sparse grid
  depends on the fact that the product rule weights are stored under colex
  order of the points, and this is the same ordering implicitly used in
  generating the SPARSE_UNIQUE_INDEX array.
*/
    grid_weight = ( double * ) malloc ( order_nd * sizeof ( double ) );

    sgmga_product_weight ( dim_num, order_1d, order_nd, rule, 
      np, p, gw_compute_weights, grid_weight );
/*
  Add these weights to the rule.
*/
    for ( order = 0; order < order_nd; order++ )
    {
      point_unique = sparse_unique_index[point_total];

      point_total = point_total + 1;

      sparse_weight[point_unique] = sparse_weight[point_unique] 
        + coef * grid_weight[order];
    }

    free ( grid_weight );
  }

  free ( level_1d );
  free ( level_1d_max );
  free ( order_1d );

  return;
}
/******************************************************************************/

void sgmga_write ( int dim_num, double level_weight[], int rule[], int np[],
  double p[], int point_num, double sparse_weight[], double sparse_point[],
  char *file_name )

/******************************************************************************/
/*
  Purpose:

    SGMGA_WRITE writes an SGMGA rule to six files.

  Discussion:

    The files are:
    * the "A" file stores the anisotropic weights, as a DIM_NUM x 1 list.
    * the "N" file stores the NP values, as a DIM_NUM x 1 list.
    * the "P" file stores the P values, as a sum(NP[*]) x 1 list.
    * the "R" file stores the region, as a DIM_NUM x 2 list;
    * the "W" file stores the weights as a POINT_NUM list;
    * the "X" file stores the abscissas as a DIM_NUM x POINT_NUM list.

    The entries in the "R" file are the two corners of the DIM_NUM dimensional
    rectangle that constitutes the integration region.  Coordinates that
    should be infinite are set to 1.0E+30.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 December 2009

  Author:

    John Burkardt

  Reference:

    Fabio Nobile, Raul Tempone, Clayton Webster,
    A Sparse Grid Stochastic Collocation Method for Partial Differential
    Equations with Random Input Data,
    SIAM Journal on Numerical Analysis,
    Volume 46, Number 5, 2008, pages 2309-2345.

    Fabio Nobile, Raul Tempone, Clayton Webster,
    An Anisotropic Sparse Grid Stochastic Collocation Method for Partial 
    Differential Equations with Random Input Data,
    SIAM Journal on Numerical Analysis,
    Volume 46, Number 5, 2008, pages 2411-2442.

  Parameters:

    Input, int DIM_NUM, the spatial dimension.

    Input, double LEVEL_WEIGHT[DIM_NUM], the anisotropic weights.

    Input, int RULE[DIM_NUM], the rule in each dimension.
     1, "CC",  Clenshaw Curtis, Closed Fully Nested rule.
     2, "F2",  Fejer Type 2, Open Fully Nested rule.
     3, "GP",  Gauss Patterson, Open Fully Nested rule.
     4, "GL",  Gauss Legendre, Open Weakly Nested rule.
     5, "GH",  Gauss Hermite, Open Weakly Nested rule.
     6, "GGH", Generalized Gauss Hermite, Open Weakly Nested rule.
     7, "LG",  Gauss Laguerre, Open Non Nested rule.
     8, "GLG", Generalized Gauss Laguerre, Open Non Nested rule.
     9, "GJ",  Gauss Jacobi, Open Non Nested rule.
    10, "GW",  Golub Welsch, (presumed) Open Non Nested rule.
    11, "CCS", Clenshaw Curtis "Slow", Closed Fully Nested rule.
    12, "F2S", Fejer Type 2 Slow, Closed Fully Nested rule.
    13, "GPS", Gauss Patterson Slow, Closed Fully Nested rule.

    Input, int NP[DIM_NUM], the number of parameters used by each rule.

    Input, double P[sum(NP[*])], the parameters needed by each rule.

    Input, int POINT_NUM, the number of unique points 
    in the grid. 

    Input, double SPARSE_WEIGHT[POINT_NUM], the weights.

    Input, double SPARSE_POINT[DIM_NUM*POINT_NUM], the points.

    Input, string FILE_NAME, the main part of the file name.
*/
{
  int dim;
  char file_name_a[255];
  char file_name_n[255];
  char file_name_p[255];
  char file_name_r[255];
  char file_name_w[255];
  char file_name_x[255];
  int np_sum;
  int point;
  double *sparse_region;
  double t1;
  double t2;

  sparse_region = ( double * ) malloc ( dim_num * 2 * sizeof ( double ) );

  for ( dim = 0; dim < dim_num; dim++ )
  {
    if ( rule[dim] == 1 )
    {
      sparse_region[dim+0*dim_num] = -1.0;
      sparse_region[dim+1*dim_num] = +1.0;
    }
    else if ( rule[dim] == 2 )
    {
      sparse_region[dim+0*dim_num] = -1.0;
      sparse_region[dim+1*dim_num] = +1.0;
    }
    else if ( rule[dim] == 3 )
    {
      sparse_region[dim+0*dim_num] = -1.0;
      sparse_region[dim+1*dim_num] = +1.0;
    }
    else if ( rule[dim] == 4 )
    {
      sparse_region[dim+0*dim_num] = -1.0;
      sparse_region[dim+1*dim_num] = +1.0;
    }
    else if ( rule[dim] == 5 )
    {
      sparse_region[dim+0*dim_num] = - r8_huge ( );
      sparse_region[dim+1*dim_num] = + r8_huge ( );
    }
    else if ( rule[dim] == 6 )
    {
      sparse_region[dim+0*dim_num] = - r8_huge ( );
      sparse_region[dim+1*dim_num] = + r8_huge ( );
    }
    else if ( rule[dim] == 7 )
    {
      sparse_region[dim+0*dim_num] = 0.0;
      sparse_region[dim+1*dim_num] = r8_huge ( );
    }
    else if ( rule[dim] == 8 )
    {
      sparse_region[dim+0*dim_num] = 0.0;
      sparse_region[dim+1*dim_num] = r8_huge ( );
    }
    else if ( rule[dim] == 9 )
    {
      sparse_region[dim+0*dim_num] = -1.0;
      sparse_region[dim+1*dim_num] = +1.0;
    }
/*
  Best guess as to region extent for rules of type 10.
*/
    else if ( rule[dim] == 10 )
    {
      t1 =   r8_huge ( );
      t2 = - r8_huge ( );
      for ( point = 0; point < point_num; point++ )
      {
        t1 = r8_min ( t1, sparse_point[dim+point*dim_num] );
        t2 = r8_max ( t2, sparse_point[dim+point*dim_num] );
      }
      sparse_region[dim+0*dim_num] = t1;
      sparse_region[dim+1*dim_num] = t2;
    }
    else if ( rule[dim] == 11 )
    {
      sparse_region[dim+0*dim_num] = -1.0;
      sparse_region[dim+1*dim_num] = +1.0;
    }
    else if ( rule[dim] == 12 )
    {
      sparse_region[dim+0*dim_num] = -1.0;
      sparse_region[dim+1*dim_num] = +1.0;
    }
    else if ( rule[dim] == 13 )
    {
      sparse_region[dim+0*dim_num] = -1.0;
      sparse_region[dim+1*dim_num] = +1.0;
    }
    else 
    {
      printf ( "\n" );
      printf ( "SGMGA_WRITE - Fatal error!\n" );
      printf ( "  Unexpected value of RULE[%d] = %d.\n", dim, rule[dim] );
      exit ( 1 );
    }
  }
  printf ( "\n" );
  printf ( "SGMGA_WRITE:\n" );

  sprintf ( file_name_a, "%s_a.txt", file_name );
  r8mat_write ( file_name_a, 1, dim_num, level_weight );
  printf ( "  Wrote the A file = \"%s\".\n", file_name_a );

  sprintf ( file_name_n, "%s_n.txt", file_name );
  i4mat_write ( file_name_n, 1, dim_num, np );
  printf ( "  Wrote the N file = \"%s\".\n", file_name_n );

  np_sum = i4vec_sum ( dim_num, np );
  sprintf ( file_name_p, "%s_p.txt", file_name );
  r8mat_write ( file_name_p, 1, np_sum, p );
  printf ( "  Wrote the P file = \"%s\".\n", file_name_p );

  sprintf ( file_name_r, "%s_r.txt", file_name );
  r8mat_write ( file_name_r, dim_num, 2, sparse_region );
  printf ( "  Wrote the R file = \"%s\".\n", file_name_r );

  sprintf ( file_name_w, "%s_w.txt", file_name );
  r8mat_write ( file_name_w, 1, point_num, sparse_weight );
  printf ( "  Wrote the W file = \"%s\".\n", file_name_w );

  sprintf ( file_name_x, "%s_x.txt", file_name );
  r8mat_write ( file_name_x, dim_num, point_num, sparse_point );
  printf ( "  Wrote the X file = \"%s\".\n", file_name_x );

  free ( sparse_region );

  return;
}
