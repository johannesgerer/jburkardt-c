# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "sandia_rules.h"
# include "sgmga.h"

int main ( void );

void sgmga_size_tests ( void );
void sgmga_size_test ( int dim_num, double importance[], double level_weight[],
  int level_max_min, int level_max_max, int rule[], int np[], double p[], 
  void ( *gw_compute_points[] ) ( int order, int np, double p[], double w[] ),
  double tol );

typedef void ( *GWPointer ) ( int order, int np, double p[], double w[] );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for SGMGA_SIZE_PRB.

  Discussion:

    SGMGA_SIZE_PRB tests the SGMGA_SIZE and SGMGA_SIZE_TOTAL functions.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    27 November 2009

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "SGMGA_SIZE_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the SGMGA_SIZE and SGMGA_SIZE_TOTAL functions.\n" );

  sgmga_size_tests ( );

  printf ( "\n" );
  printf ( "SGMGA__SIZE_PRB\n" );
  printf ( "  Normal end of execution.\n" );

  printf ( "\n" );
  timestamp ( );
  
  return 0;
}
/******************************************************************************/

void sgmga_size_tests ( void )

/******************************************************************************/
/*
  Purpose:

    SGMGA_SIZE_TESTS calls SGMGA_SIZE_TEST.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    27 November 2009

  Author:

    John Burkardt

  Local Parameters:

    Local, double TOL, a tolerance for point equality.
    A value of sqrt ( eps ) is reasonable, and will allow the code to
    consolidate points which are equal, or very nearly so.  A value of
    -1.0, on the other hand, will force the code to use every point, 
    regardless of duplication.
*/
{
  int dim;
  int dim_num;
  GWPointer *gw_compute_points;
  double *importance;
  int level_max_max;
  int level_max_min;
  double *level_weight;
  int *np;
  int np_sum;
  int *order_1d;
  int order_nd;
  double *p;
  int *rule;
  double tol;

  printf ( "\n" );
  printf ( "SGMGA_SIZE_TESTS\n" );
  printf ( "  Call SGMGA_SIZE_TEST with various arguments.\n" );
/*
  Set the point equality tolerance.
*/
  tol = sqrt ( r8_epsilon ( ) );
  printf ( "\n" );
  printf ( "  Point equality tolerance = %e\n", tol );

  dim_num = 2;
  importance = ( double * ) malloc ( dim_num * sizeof ( double ) );
  for ( dim = 0; dim < dim_num; dim++ )
  {
    importance[dim] = 1.0;
  }
  level_weight = ( double * ) malloc ( dim_num * sizeof ( double ) );
  sgmga_importance_to_aniso ( dim_num, importance, level_weight );
  level_max_min = 0;
  level_max_max = 5;
  rule = ( int * ) malloc ( dim_num * sizeof ( int ) );
  rule[0] = 1;
  rule[1] = 1;
  np = ( int * ) malloc ( dim_num * sizeof ( int ) );
  np[0] = 0;
  np[1] = 0;
  np_sum = i4vec_sum ( dim_num, np );
  p = ( double * ) malloc ( np_sum * sizeof ( double ) );
  gw_compute_points = ( GWPointer * ) malloc ( dim_num * sizeof ( GWPointer ) );
  gw_compute_points[0] = clenshaw_curtis_compute_points_np;
  gw_compute_points[1] = clenshaw_curtis_compute_points_np;
  sgmga_size_test ( dim_num, importance, level_weight, level_max_min, 
    level_max_max, rule, np, p, gw_compute_points, tol );
  free ( gw_compute_points );
  free ( importance );
  free ( level_weight );
  free ( np );
  free ( p );
  free ( rule );

  dim_num = 2;
  importance = ( double * ) malloc ( dim_num * sizeof ( double ) );
  for ( dim = 0; dim < dim_num; dim++ )
  {
    importance[dim] = ( double ) ( dim + 1 );
  }
  level_weight = ( double * ) malloc ( dim_num * sizeof ( double ) );
  sgmga_importance_to_aniso ( dim_num, importance, level_weight );
  level_max_min = 0;
  level_max_max = 5;
  rule = ( int * ) malloc ( dim_num * sizeof ( int ) );
  rule[0] = 1;
  rule[1] = 1;
  np = ( int * ) malloc ( dim_num * sizeof ( int ) );
  np[0] = 0;
  np[1] = 0;
  np_sum = i4vec_sum ( dim_num, np );
  p = ( double * ) malloc ( np_sum * sizeof ( double ) );
  gw_compute_points = ( GWPointer * ) malloc ( dim_num * sizeof ( GWPointer ) );
  gw_compute_points[0] = clenshaw_curtis_compute_points_np;
  gw_compute_points[1] = clenshaw_curtis_compute_points_np;
  sgmga_size_test ( dim_num, importance, level_weight, level_max_min, 
    level_max_max, rule, np, p, gw_compute_points, tol );
  free ( gw_compute_points );
  free ( importance );
  free ( level_weight );
  free ( np );
  free ( p );
  free ( rule );

  dim_num = 3;
  importance = ( double * ) malloc ( dim_num * sizeof ( double ) );
  for ( dim = 0; dim < dim_num; dim++ )
  {
    importance[dim] = 1.0;
  }
  level_weight = ( double * ) malloc ( dim_num * sizeof ( double ) );
  sgmga_importance_to_aniso ( dim_num, importance, level_weight );
  level_max_min = 0;
  level_max_max = 5;
  rule = ( int * ) malloc ( dim_num * sizeof ( int ) );
  rule[0] = 1;
  rule[1] = 1;
  rule[2] = 1;
  np = ( int * ) malloc ( dim_num * sizeof ( int ) );
  np[0] = 0;
  np[1] = 0;
  np[2] = 0;
  np_sum = i4vec_sum ( dim_num, np );
  p = ( double * ) malloc ( np_sum * sizeof ( double ) );
  gw_compute_points = ( GWPointer * ) malloc ( dim_num * sizeof ( GWPointer ) );
  gw_compute_points[0] = clenshaw_curtis_compute_points_np;
  gw_compute_points[1] = clenshaw_curtis_compute_points_np;
  sgmga_size_test ( dim_num, importance, level_weight, level_max_min, 
    level_max_max, rule, np, p, gw_compute_points, tol );
  free ( gw_compute_points );
  free ( importance );
  free ( level_weight );
  free ( np );
  free ( p );
  free ( rule );

  dim_num = 3;
  importance = ( double * ) malloc ( dim_num * sizeof ( double ) );
  for ( dim = 0; dim < dim_num; dim++ )
  {
    importance[dim] = ( double ) ( dim + 1 );
  }
  level_weight = ( double * ) malloc ( dim_num * sizeof ( double ) );
  sgmga_importance_to_aniso ( dim_num, importance, level_weight );
  level_max_min = 0;
  level_max_max = 5;
  rule = ( int * ) malloc ( dim_num * sizeof ( int ) );
  rule[0] = 1;
  rule[1] = 1;
  rule[2] = 1;
  np = ( int * ) malloc ( dim_num * sizeof ( int ) );
  np[0] = 0;
  np[1] = 0;
  np[2] = 0;
  np_sum = i4vec_sum ( dim_num, np );
  p = ( double * ) malloc ( np_sum * sizeof ( double ) );
  gw_compute_points = ( GWPointer * ) malloc ( dim_num * sizeof ( GWPointer ) );
  gw_compute_points[0] = clenshaw_curtis_compute_points_np;
  gw_compute_points[1] = clenshaw_curtis_compute_points_np;
  sgmga_size_test ( dim_num, importance, level_weight, level_max_min, 
    level_max_max, rule, np, p, gw_compute_points, tol );
  free ( gw_compute_points );
  free ( importance );
  free ( level_weight );
  free ( np );
  free ( p );
  free ( rule );

  dim_num = 2;
  importance = ( double * ) malloc ( dim_num * sizeof ( double ) );
  for ( dim = 0; dim < dim_num; dim++ )
  {
    importance[dim] = ( double ) ( dim + 1 );
  }
  level_weight = ( double * ) malloc ( dim_num * sizeof ( double ) );
  sgmga_importance_to_aniso ( dim_num, importance, level_weight );
  level_max_min = 0;
  level_max_max = 2;
  rule = ( int * ) malloc ( dim_num * sizeof ( int ) );
  rule[0] = 1;
  rule[1] = 3;
  np = ( int * ) malloc ( dim_num * sizeof ( int ) );
  np[0] = 0;
  np[1] = 0;
  np_sum = i4vec_sum ( dim_num, np );
  p = ( double * ) malloc ( np_sum * sizeof ( double ) );
  gw_compute_points = ( GWPointer * ) malloc ( dim_num * sizeof ( GWPointer ) );
  gw_compute_points[0] = clenshaw_curtis_compute_points_np;
  gw_compute_points[1] = patterson_lookup_points_np;
  sgmga_size_test ( dim_num, importance, level_weight, level_max_min, 
    level_max_max, rule, np, p, gw_compute_points, tol );
  free ( gw_compute_points );
  free ( importance );
  free ( level_weight );
  free ( np );
  free ( p );
  free ( rule );

  dim_num = 2;
  importance = ( double * ) malloc ( dim_num * sizeof ( double ) );
  for ( dim = 0; dim < dim_num; dim++ )
  {
    importance[dim] = ( double ) ( dim + 1 );
  }
  level_weight = ( double * ) malloc ( dim_num * sizeof ( double ) );
  sgmga_importance_to_aniso ( dim_num, importance, level_weight );
  level_max_min = 0;
  level_max_max = 2;
  rule = ( int * ) malloc ( dim_num * sizeof ( int ) );
  rule[0] = 1;
  rule[1] = 4;
  np = ( int * ) malloc ( dim_num * sizeof ( int ) );
  np[0] = 0;
  np[1] = 0;
  np_sum = i4vec_sum ( dim_num, np );
  p = ( double * ) malloc ( np_sum * sizeof ( double ) );
  gw_compute_points = ( GWPointer * ) malloc ( dim_num * sizeof ( GWPointer ) );
  gw_compute_points[0] = clenshaw_curtis_compute_points_np;
  gw_compute_points[1] = legendre_compute_points_np;
  sgmga_size_test ( dim_num, importance, level_weight, level_max_min, 
    level_max_max, rule, np, p, gw_compute_points, tol );
  free ( gw_compute_points );
  free ( importance );
  free ( level_weight );
  free ( np );
  free ( p );
  free ( rule );

  dim_num = 2;
  importance = ( double * ) malloc ( dim_num * sizeof ( double ) );
  for ( dim = 0; dim < dim_num; dim++ )
  {
    importance[dim] = ( double ) ( dim + 1 );
  }
  level_weight = ( double * ) malloc ( dim_num * sizeof ( double ) );
  sgmga_importance_to_aniso ( dim_num, importance, level_weight );
  level_max_min = 0;
  level_max_max = 2;
  rule = ( int * ) malloc ( dim_num * sizeof ( int ) );
  rule[0] = 1;
  rule[1] = 7;
  np = ( int * ) malloc ( dim_num * sizeof ( int ) );
  np[0] = 0;
  np[1] = 0;
  np_sum = i4vec_sum ( dim_num, np );
  p = ( double * ) malloc ( np_sum * sizeof ( double ) );
  gw_compute_points = ( GWPointer * ) malloc ( dim_num * sizeof ( GWPointer ) );
  gw_compute_points[0] = clenshaw_curtis_compute_points_np;
  gw_compute_points[1] = laguerre_compute_points_np;
  sgmga_size_test ( dim_num, importance, level_weight, level_max_min, 
    level_max_max, rule, np, p, gw_compute_points, tol );
  free ( gw_compute_points );
  free ( importance );
  free ( level_weight );
  free ( np );
  free ( p );
  free ( rule );

  dim_num = 2;
  importance = ( double * ) malloc ( dim_num * sizeof ( double ) );
  for ( dim = 0; dim < dim_num; dim++ )
  {
    importance[dim] = ( double ) ( dim + 1 );
  }
  level_weight = ( double * ) malloc ( dim_num * sizeof ( double ) );
  sgmga_importance_to_aniso ( dim_num, importance, level_weight );
  level_max_min = 0;
  level_max_max = 2;
  rule = ( int * ) malloc ( dim_num * sizeof ( int ) );
  rule[0] = 1;
  rule[1] = 8;
  np = ( int * ) malloc ( dim_num * sizeof ( int ) );
  np[0] = 0;
  np[1] = 1;
  np_sum = i4vec_sum ( dim_num, np );
  p = ( double * ) malloc ( np_sum * sizeof ( double ) );
  p[0] = 1.5;
  gw_compute_points = ( GWPointer * ) malloc ( dim_num * sizeof ( GWPointer ) );
  gw_compute_points[0] = clenshaw_curtis_compute_points_np;
  gw_compute_points[1] = gen_laguerre_compute_points_np;
  sgmga_size_test ( dim_num, importance, level_weight, level_max_min, 
    level_max_max, rule, np, p, gw_compute_points, tol );
  free ( gw_compute_points );
  free ( importance );
  free ( level_weight );
  free ( np );
  free ( p );
  free ( rule );

  dim_num = 2;
  importance = ( double * ) malloc ( dim_num * sizeof ( double ) );
  for ( dim = 0; dim < dim_num; dim++ )
  {
    importance[dim] = ( double ) ( dim + 1 );
  }
  level_weight = ( double * ) malloc ( dim_num * sizeof ( double ) );
  sgmga_importance_to_aniso ( dim_num, importance, level_weight );
  level_max_min = 0;
  level_max_max = 2;
  rule = ( int * ) malloc ( dim_num * sizeof ( int ) );
  rule[0] = 2;
  rule[1] = 9;
  np = ( int * ) malloc ( dim_num * sizeof ( int ) );
  np[0] = 0;
  np[1] = 2;
  np_sum = i4vec_sum ( dim_num, np );
  p = ( double * ) malloc ( np_sum * sizeof ( double ) );
  p[0] = 0.5;
  p[1] = 1.5;
  gw_compute_points = ( GWPointer * ) malloc ( dim_num * sizeof ( GWPointer ) );
  gw_compute_points[0] = fejer2_compute_points_np;
  gw_compute_points[1] = jacobi_compute_points_np;
  sgmga_size_test ( dim_num, importance, level_weight, level_max_min, 
    level_max_max, rule, np, p, gw_compute_points, tol );
  free ( gw_compute_points );
  free ( importance );
  free ( level_weight );
  free ( np );
  free ( p );
  free ( rule );

  dim_num = 2;
  importance = ( double * ) malloc ( dim_num * sizeof ( double ) );
  for ( dim = 0; dim < dim_num; dim++ )
  {
    importance[dim] = ( double ) ( dim + 1 );
  }
  level_weight = ( double * ) malloc ( dim_num * sizeof ( double ) );
  sgmga_importance_to_aniso ( dim_num, importance, level_weight );
  level_max_min = 0;
  level_max_max = 2;
  rule = ( int * ) malloc ( dim_num * sizeof ( int ) );
  rule[0] = 6;
  rule[1] = 4;
  np = ( int * ) malloc ( dim_num * sizeof ( int ) );
  np[0] = 1;
  np[1] = 0;
  np_sum = i4vec_sum ( dim_num, np );
  p = ( double * ) malloc ( np_sum * sizeof ( double ) );
  p[0] = 2.0;
  gw_compute_points = ( GWPointer * ) malloc ( dim_num * sizeof ( GWPointer ) );
  gw_compute_points[0] = gen_hermite_compute_points_np;
  gw_compute_points[1] = legendre_compute_points_np;
  sgmga_size_test ( dim_num, importance, level_weight, level_max_min, 
    level_max_max, rule, np, p, gw_compute_points, tol );
  free ( gw_compute_points );
  free ( importance );
  free ( level_weight );
  free ( np );
  free ( p );
  free ( rule );

  dim_num = 3;
  importance = ( double * ) malloc ( dim_num * sizeof ( double ) );
  for ( dim = 0; dim < dim_num; dim++ )
  {
    importance[dim] = ( double ) ( dim + 1 );
  }
  level_weight = ( double * ) malloc ( dim_num * sizeof ( double ) );
  sgmga_importance_to_aniso ( dim_num, importance, level_weight );
  level_max_min = 0;
  level_max_max = 2;
  rule = ( int * ) malloc ( dim_num * sizeof ( int ) );
  rule[0] = 1;
  rule[1] = 4;
  rule[2] = 5;
  np = ( int * ) malloc ( dim_num * sizeof ( int ) );
  np[0] = 0;
  np[1] = 0;
  np[2] = 0;
  np_sum = i4vec_sum ( dim_num, np );
  p = ( double * ) malloc ( np_sum * sizeof ( double ) );
  gw_compute_points = ( GWPointer * ) malloc ( dim_num * sizeof ( GWPointer ) );
  gw_compute_points[0] = clenshaw_curtis_compute_points_np;
  gw_compute_points[1] = legendre_compute_points_np;
  gw_compute_points[2] = hermite_compute_points_np;
  sgmga_size_test ( dim_num, importance, level_weight, level_max_min, 
    level_max_max, rule, np, p, gw_compute_points, tol );
  free ( gw_compute_points );
  free ( importance );
  free ( level_weight );
  free ( np );
  free ( p );
  free ( rule );
/*
  Repeat, treating rules #2 and #3 as Golub Welsch rules.
*/
  dim_num = 3;
  importance = ( double * ) malloc ( dim_num * sizeof ( double ) );
  for ( dim = 0; dim < dim_num; dim++ )
  {
    importance[dim] = ( double ) ( dim + 1 );
  }
  level_weight = ( double * ) malloc ( dim_num * sizeof ( double ) );
  sgmga_importance_to_aniso ( dim_num, importance, level_weight );
  level_max_min = 0;
  level_max_max = 2;
  rule = ( int * ) malloc ( dim_num * sizeof ( int ) );
  rule[0] = 1;
  rule[1] = 10;
  rule[2] = 10;
  np = ( int * ) malloc ( dim_num * sizeof ( int ) );
  np[0] = 0;
  np[1] = 0;
  np[2] = 0;
  np_sum = i4vec_sum ( dim_num, np );
  p = ( double * ) malloc ( np_sum * sizeof ( double ) );
  gw_compute_points = ( GWPointer * ) malloc ( dim_num * sizeof ( GWPointer ) );
  gw_compute_points[0] = clenshaw_curtis_compute_points_np;
  gw_compute_points[1] = legendre_compute_points_np;
  gw_compute_points[2] = hermite_compute_points_np;
  sgmga_size_test ( dim_num, importance, level_weight, level_max_min, 
    level_max_max, rule, np, p, gw_compute_points, tol );
  free ( gw_compute_points );
  free ( importance );
  free ( level_weight );
  free ( np );
  free ( p );
  free ( rule );
/*
  Try a case involving a dimension of "0" importance.
*/
  dim_num = 3;
  importance = ( double * ) malloc ( dim_num * sizeof ( double ) );
  importance[0] = 1.0;
  importance[1] = 0.0;
  importance[2] = 1.0;
  level_weight = ( double * ) malloc ( dim_num * sizeof ( double ) );
  sgmga_importance_to_aniso ( dim_num, importance, level_weight );
  level_max_min = 0;
  level_max_max = 3;
  rule = ( int * ) malloc ( dim_num * sizeof ( int ) );
  rule[0] = 1;
  rule[1] = 1;
  rule[2] = 1;
  np = ( int * ) malloc ( dim_num * sizeof ( int ) );
  np[0] = 0;
  np[1] = 0;
  np[2] = 0;
  np_sum = i4vec_sum ( dim_num, np );
  p = ( double * ) malloc ( np_sum * sizeof ( double ) );
  gw_compute_points = ( GWPointer * ) malloc ( dim_num * sizeof ( GWPointer ) );
  gw_compute_points[0] = clenshaw_curtis_compute_points_np;
  gw_compute_points[1] = clenshaw_curtis_compute_points_np;
  gw_compute_points[2] = clenshaw_curtis_compute_points_np;
  sgmga_size_test ( dim_num, importance, level_weight, level_max_min, 
    level_max_max, rule, np, p, gw_compute_points, tol );
  free ( gw_compute_points );
  free ( importance );
  free ( level_weight );
  free ( np );
  free ( p );
  free ( rule );

  return;
}
/******************************************************************************/

void sgmga_size_test ( int dim_num, double importance[], double level_weight[], 
  int level_max_min, int level_max_max, int rule[], int np[], double p[], 
  void ( *gw_compute_points[] ) ( int order, int np, double p[], double w[] ),
  double tol )

/******************************************************************************/
/*
  Purpose:

    SGMGA_SIZE_TEST tests SGMG_SIZE, SGMG_SIZE_TOTAL.

  Discussion:

    The "level_to_order" argument for the SGMG_SIZE and SGMG_SIZE_TOTAL
    functions is set to "LEVEL_TO_ORDER_DEFAULT", which uses exponential
    growth for fully nested rules, and linear otherwise.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 November 2009

  Author:

    John Burkardt

  Parameters:

    Input, int DIM_NUM, the spatial dimension.

    Input, double IMPORTANCE[DIM_NUM], the importance for each dimension.

    Input, double LEVEL_WEIGHT[DIM_NUM], the weights for each dimension.

    Input, int LEVEL_MAX_MIN, LEVEL_MAX_MAX, the minimum and
    maximum values of LEVEL_MAX.

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

    Input, int NP[RULE_NUM], the number of parameters used by each rule.

    Input, double P[sum(NP[*])], the parameters needed by each rule.

    Input, void ( *GW_COMPUTE_WEIGHTS[] ) ( int order, int np, double p[], double w[] ),
    an array of pointers to functions which return the 1D quadrature weights 
    associated with each spatial dimension for which a Golub Welsch rule 
    is used.

    Input, double TOL, a tolerance for point equality.
*/
{
  double alpha;
  double beta;
  int dim;
  int i;
  int level_max;
  int p_index;
  int point_num;
  int point_total_num;

  printf ( "\n" );
  printf ( "SGMGA_SIZE_TEST\n" );
  printf ( "  SGMGA_SIZE_TOTAL counts the total number of points,\n" );
  printf ( "  including duplications, in an SGMGA sparse grid.\n" );
  printf ( "  SGMGA_SIZE counts the total number of points,\n" );
  printf ( "  excluding duplications, in an SGMGA sparse grid.\n" );
  printf ( "\n" );
  printf ( "  IMPORTANCE:  " );
  for ( dim = 0; dim < dim_num; dim++ )
  {
    printf ( "  %14f", importance[dim] );
  }
  printf ( "\n" );
  printf ( "  LEVEL_WEIGHT:" );
  for ( dim = 0; dim < dim_num; dim++ )
  {
    printf ( "  %14f", level_weight[dim] );
  }
  printf ( "\n" );
  printf ( "\n" );
  printf ( " Dimension      Rule  Growth rate       Parameters\n" );
  printf ( "\n" );

  p_index = 0;

  for ( dim = 0; dim < dim_num; dim++ )
  {
    if ( rule[dim] == 1 )
    {
      printf ( "  %8d  %8d   Exponential\n", dim, rule[dim] );
    }
    else if ( rule[dim] == 2 )
    {
      printf ( "  %8d  %8d   Exponential\n", dim, rule[dim] );
    }
    else if ( rule[dim] == 3 )
    {
      printf ( "  %8d  %8d   Exponential\n", dim, rule[dim] );
    }
    else if ( rule[dim] == 4 )
    {
      printf ( "  %8d  %8d   Linear\n", dim, rule[dim] );
    }
    else if ( rule[dim] == 5 )
    {
      printf ( "  %8d  %8d   Linear\n", dim, rule[dim] );
    }
    else if ( rule[dim] == 6 )
    {
      alpha = p[p_index];
      p_index = p_index + 1;
      printf ( "  %8d  %8d   Linear     %14f\n", dim, rule[dim], alpha );
    }
    else if ( rule[dim] == 7 )
    {
      printf ( "  %8d  %8d   Linear\n", dim, rule[dim] );
    }
    else if ( rule[dim] == 8 )
    {
      alpha = p[p_index];
      p_index = p_index + 1;
      printf ( "  %8d  %8d   Linear     %14f\n", dim, rule[dim], alpha );
    }
    else if ( rule[dim] == 9 )
    {
      alpha = p[p_index];
      p_index = p_index + 1;
      beta = p[p_index];
      p_index = p_index + 1;
      printf ( "  %8d  %8d   Linear     %14f  %14f\n", dim, rule[dim], alpha, beta );
    }
    else if ( rule[dim] == 10 )
    {
      printf ( "  %8d  %8d   Linear     ", dim, rule[dim] );
      for ( i = 0; i < np[dim]; i++ )
      {
        alpha = p[p_index];
        p_index = p_index + 1;
        printf ( "  %14f", alpha );
      }
      printf ( "\n" );
    }
    else if ( rule[dim] == 11 )
    {
      printf ( "  %8d  %8d   Slow Expo  \n", dim, rule[dim] );
    }
    else
    {
      printf ( "\n" );
      printf ( "SGMGA_SIZE_TEST - Fatal error!\n" );
      printf ( "  Unexpected value of RULE = %d\n", rule[dim] );
      exit ( 1 );
    }
  }

  printf ( "\n" );
  printf ( "   DIM_NUM LEVEL_MAX POINT_NUM POINT_NUM\n" );
  printf ( "                        Unique     Total\n" );
  printf ( "\n" );

  for ( level_max = level_max_min; level_max <= level_max_max; level_max++ )
  {
    point_total_num = sgmga_size_total ( dim_num, level_weight,
      level_max, rule, level_to_order_default );

    point_num = sgmga_size ( dim_num, level_weight, level_max, 
      rule, np, p, gw_compute_points, tol, 
      level_to_order_default );

    printf ( "  %8d  %8d  %8d  %8d\n", 
      dim_num, level_max, point_num, point_total_num );
  }

  return;
}
