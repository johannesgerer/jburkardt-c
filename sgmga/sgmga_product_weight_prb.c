# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "sandia_rules.h"
# include "sgmga.h"

int main ( void );

void sgmga_product_weight_tests ( );
void sgmga_product_weight_test ( int dim_num, int order_1d[], 
  int order_nd, int rule[], int np[], double p[],
  void ( *gw_compute_weights[] ) ( int order, int np, double p[], double w[] ) );

typedef void ( *GWPointer ) ( int order, int np, double p[], double w[] );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for SGMGA_PRODUCT_WEIGHT_PRB.

  Discussion:

    SGMGA_PRODUCT_WEIGHT_PRB tests the SGMGA_PRODUCT_WEIGHT function.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 November 2009

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "SGMGA_PRODUCT_WEIGHT_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the SGMGA_PRODUCT_WEIGHT function.\n" );
/*
  Make sure the individual product rule weights are computed correctly.
*/
  sgmga_product_weight_tests ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "SGMGA_PRODUCT_WEIGHT_PRB\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );
  
  return 0;
}
/******************************************************************************/

void sgmga_product_weight_tests ( )

/******************************************************************************/
/*
  Purpose:

    SGMGA_PRODUCT_WEIGHT_TESTS calls SGMGA_PRODUCT_WEIGHT_TEST.

  Discussion:

    To test Golub Welsch rules for a spatial dimension DIM, we can
    set RULE[DIM] = 10, and set the corresponding entry of 
    GW_COMPUTE_WEIGHTS to the name of a function that we know is already
    available, such as "clenshaw_curtis_compute_weights".

    Note that, for ALL the tests, we set every entry of the GW_COMPUTE_WEIGHTS
    array.  However, a particular entry is only inspected if the corresponding
    entry of RULE is 10.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 November 2009

  Author:

    John Burkardt
*/
{
  int dim_num;
  GWPointer *gw_compute_weights;
  int *np;
  int np_sum;
  int *order_1d;
  int order_nd;
  double *p;
  int *rule;

  printf ( "\n" );
  printf ( "SGMGA_PRODUCT_WEIGHT_TESTS\n" );
  printf ( "  Call SGMGA_PRODUCT_WEIGHT_TEST with various arguments.\n" );

  dim_num = 2;
  order_1d = ( int * ) malloc ( dim_num * sizeof ( int ) );
  order_1d[0] = 3;
  order_1d[1] = 5;
  order_nd = i4vec_product ( dim_num, order_1d );
  rule = ( int * ) malloc ( dim_num * sizeof ( int ) );
  rule[0] = 1;
  rule[1] = 1;
  np = ( int * ) malloc ( dim_num * sizeof ( int ) );
  np[0] = 0;
  np[1] = 0;
  np_sum = i4vec_sum ( dim_num, np );
  p = ( double * ) malloc ( np_sum * sizeof ( double ) );
  gw_compute_weights = ( GWPointer * ) malloc ( dim_num * sizeof ( GWPointer ) );
  gw_compute_weights[0] = clenshaw_curtis_compute_weights_np;
  gw_compute_weights[1] = clenshaw_curtis_compute_weights_np;
  sgmga_product_weight_test ( dim_num, order_1d, order_nd, rule, 
    np, p, gw_compute_weights );
  free ( gw_compute_weights );
  free ( np );
  free ( order_1d );
  free ( p );
  free ( rule );

  dim_num = 2;
  order_1d = ( int * ) malloc ( dim_num * sizeof ( int ) );
  order_1d[0] = 3;
  order_1d[1] = 7;
  order_nd = i4vec_product ( dim_num, order_1d );
  rule = ( int * ) malloc ( dim_num * sizeof ( int ) );
  rule[0] = 1;
  rule[1] = 5;
  np = ( int * ) malloc ( dim_num * sizeof ( int ) );
  np[0] = 0;
  np[1] = 0;
  np_sum = i4vec_sum ( dim_num, np );
  p = ( double * ) malloc ( np_sum * sizeof ( double ) );
  gw_compute_weights = ( GWPointer * ) malloc ( dim_num * sizeof ( GWPointer ) );
  gw_compute_weights[0] = clenshaw_curtis_compute_weights_np;
  gw_compute_weights[1] = hermite_compute_weights_np;
  sgmga_product_weight_test ( dim_num, order_1d, order_nd, rule, 
    np, p, gw_compute_weights );
  free ( gw_compute_weights );
  free ( np );
  free ( order_1d );
  free ( p );
  free ( rule );

  dim_num = 2;
  order_1d = ( int * ) malloc ( dim_num * sizeof ( int ) );
  order_1d[0] = 3;
  order_1d[1] = 3;
  order_nd = i4vec_product ( dim_num, order_1d );
  rule = ( int * ) malloc ( dim_num * sizeof ( int ) );
  rule[0] = 3;
  rule[1] = 7;
  np = ( int * ) malloc ( dim_num * sizeof ( int ) );
  np[0] = 0;
  np[1] = 0;
  np_sum = i4vec_sum ( dim_num, np );
  p = ( double * ) malloc ( np_sum * sizeof ( double ) );
  gw_compute_weights = ( GWPointer * ) malloc ( dim_num * sizeof ( GWPointer ) );
  gw_compute_weights[0] = patterson_lookup_weights_np;
  gw_compute_weights[1] = laguerre_compute_weights_np;
  sgmga_product_weight_test ( dim_num, order_1d, order_nd, rule, 
    np, p, gw_compute_weights );
  free ( gw_compute_weights );
  free ( np );
  free ( order_1d );
  free ( p );
  free ( rule );

  dim_num = 2;
  order_1d = ( int * ) malloc ( dim_num * sizeof ( int ) );
  order_1d[0] = 5;
  order_1d[1] = 5;
  order_nd = i4vec_product ( dim_num, order_1d );
  rule = ( int * ) malloc ( dim_num * sizeof ( int ) );
  rule[0] = 1;
  rule[1] = 8;
  np = ( int * ) malloc ( dim_num * sizeof ( int ) );
  np[0] = 0;
  np[1] = 1;
  np_sum = i4vec_sum ( dim_num, np );
  p = ( double * ) malloc ( np_sum * sizeof ( double ) );
  p[0] = 1.5;
  gw_compute_weights = ( GWPointer * ) malloc ( dim_num * sizeof ( GWPointer ) );
  gw_compute_weights[0] = clenshaw_curtis_compute_weights_np;
  gw_compute_weights[1] = gen_laguerre_compute_weights_np;
  sgmga_product_weight_test ( dim_num, order_1d, order_nd, rule, 
    np, p, gw_compute_weights );
  free ( gw_compute_weights );
  free ( np );
  free ( order_1d );
  free ( p );
  free ( rule );

  dim_num = 2;
  order_1d = ( int * ) malloc ( dim_num * sizeof ( int ) );
  order_1d[0] = 5;
  order_1d[1] = 5;
  order_nd = i4vec_product ( dim_num, order_1d );
  rule = ( int * ) malloc ( dim_num * sizeof ( int ) );
  rule[0] = 2;
  rule[1] = 9;
  np = ( int * ) malloc ( dim_num * sizeof ( int ) );
  np[0] = 0;
  np[1] = 0;
  np_sum = i4vec_sum ( dim_num, np );
  p = ( double * ) malloc ( np_sum * sizeof ( double ) );
  p[0] = 0.5;
  p[1] = 1.5;
  gw_compute_weights = ( GWPointer * ) malloc ( dim_num * sizeof ( GWPointer ) );
  gw_compute_weights[0] = fejer2_compute_weights_np;
  gw_compute_weights[1] = jacobi_compute_weights_np;
  sgmga_product_weight_test ( dim_num, order_1d, order_nd, rule, 
    np, p, gw_compute_weights );
  free ( gw_compute_weights );
  free ( np );
  free ( order_1d );
  free ( p );
  free ( rule );

  dim_num = 2;
  order_1d = ( int * ) malloc ( dim_num * sizeof ( int ) );
  order_1d[0] = 7;
  order_1d[1] = 7;
  order_nd = i4vec_product ( dim_num, order_1d );
  rule = ( int * ) malloc ( dim_num * sizeof ( int ) );
  rule[0] = 6;
  rule[1] = 4;
  np = ( int * ) malloc ( dim_num * sizeof ( int ) );
  np[0] = 0;
  np[1] = 0;
  np_sum = i4vec_sum ( dim_num, np );
  p = ( double * ) malloc ( np_sum * sizeof ( double ) );
  p[0] = 2.0;
  gw_compute_weights = ( GWPointer * ) malloc ( dim_num * sizeof ( GWPointer ) );
  gw_compute_weights[0] = gen_hermite_compute_weights_np;
  gw_compute_weights[1] = legendre_compute_weights_np;
  sgmga_product_weight_test ( dim_num, order_1d, order_nd, rule, 
    np, p, gw_compute_weights );
  free ( gw_compute_weights );
  free ( np );
  free ( order_1d );
  free ( p );
  free ( rule );

  dim_num = 3;
  order_1d = ( int * ) malloc ( dim_num * sizeof ( int ) );
  order_1d[0] = 2;
  order_1d[1] = 3;
  order_1d[2] = 3;
  order_nd = i4vec_product ( dim_num, order_1d );
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
  gw_compute_weights = ( GWPointer * ) malloc ( dim_num * sizeof ( GWPointer ) );
  gw_compute_weights[0] = clenshaw_curtis_compute_weights_np;
  gw_compute_weights[1] = legendre_compute_weights_np;
  gw_compute_weights[2] = hermite_compute_weights_np;
  sgmga_product_weight_test ( dim_num, order_1d, order_nd, rule, 
    np, p, gw_compute_weights );
  free ( gw_compute_weights );
  free ( np );
  free ( order_1d );
  free ( p );
  free ( rule );
/*
  Repeat, treating  rules #2 and #3 as Golub Welsch rules.
*/
  dim_num = 3;
  order_1d = ( int * ) malloc ( dim_num * sizeof ( int ) );
  order_1d[0] = 2;
  order_1d[1] = 3;
  order_1d[2] = 3;
  order_nd = i4vec_product ( dim_num, order_1d );
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
  gw_compute_weights = ( GWPointer * ) malloc ( dim_num * sizeof ( GWPointer ) );
  gw_compute_weights[0] = clenshaw_curtis_compute_weights_np;
  gw_compute_weights[1] = legendre_compute_weights_np;
  gw_compute_weights[2] = hermite_compute_weights_np;
  sgmga_product_weight_test ( dim_num, order_1d, order_nd, rule, 
    np, p, gw_compute_weights );
  free ( gw_compute_weights );
  free ( np );
  free ( order_1d );
  free ( p );
  free ( rule );

  return;
}
/******************************************************************************/

void sgmga_product_weight_test ( int dim_num, int order_1d[], 
  int order_nd, int rule[], int np[], double p[],
  void ( *gw_compute_weights[] ) ( int order, int np, double p[], double w[] ) )

/******************************************************************************/
/*
  Purpose:

    SGMGA_PRODUCT_WEIGHT_TEST: weights of a mixed factor product rule.

  Discussion:

    This routine computes a sparse grid and compares the sum of the weights
    to the expected exact value.

    The routine cannot produce a result for rules that include one or more
    component rules of type 10, that is, Golub-Welsch rules.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 November 2009

  Author:

    John Burkardt

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

    Input, int NP[RULE_NUM], the number of parameters used by each rule.

    Input, double P[sum(NP[*])], the parameters needed by each rule.

    Input, void ( *GW_COMPUTE_WEIGHTS[] ) ( int order, int np, double p[], double w[] ),
    an array of pointers to functions which return the 1D quadrature weights 
    associated with each spatial dimension for which a Golub Welsch rule 
    is used.
*/
{
  double alpha;
  double beta;
  double arg1;
  double arg2;
  double arg3;
  double arg4;
  int dim;
  int i;
  int p_index;
  double pi = 3.141592653589793;
  double value1;
  double value2;
  double *weight;
  double weight_sum;
  double weight_sum_error;
  double weight_sum_exact;
/*
  Determine the integral of 1 over the multidimensional weighted region.
*/
  p_index = 0;

  weight_sum_exact = 1.0;

  for ( dim = 0; dim < dim_num; dim++ )
  {
    if ( rule[dim] == 1 )
    {
      weight_sum_exact = weight_sum_exact * 2.0;
    }
    else if ( rule[dim] == 2 )
    {
      weight_sum_exact = weight_sum_exact * 2.0;
    }
    else if ( rule[dim] == 3 )
    {
      weight_sum_exact = weight_sum_exact * 2.0;
    }
    else if ( rule[dim] == 4 )
    {
      weight_sum_exact = weight_sum_exact * 2.0;
    }
    else if ( rule[dim] == 5 )
    {
      weight_sum_exact = weight_sum_exact * sqrt ( pi );
    }
    else if ( rule[dim] == 6 )
    {
      alpha = p[p_index];
      p_index = p_index + 1;

      weight_sum_exact = weight_sum_exact 
        * r8_gamma ( 0.5 * ( alpha + 1.0 ) );
    }
    else if ( rule[dim] == 7 )
    {
      weight_sum_exact = weight_sum_exact * 1.0;
    }
    else if ( rule[dim] == 8 )
    {
      alpha = p[p_index];
      p_index = p_index + 1;

      weight_sum_exact = weight_sum_exact * r8_gamma ( alpha + 1.0 );
    }
    else if ( rule[dim] == 9 )
    {
      alpha = p[p_index];
      p_index = p_index + 1;
      beta = p[p_index];
      p_index = p_index + 1;
      arg1 = - alpha;
      arg2 = 1.0;
      arg3 = beta + 2.0;
      arg4 = - 1.0;
      value1 = r8_hyper_2f1 ( arg1, arg2, arg3, arg4 );
      arg1 = - beta;
      arg2 = 1.0;
      arg3 = alpha + 2.0;
      arg4 = - 1.0;
      value2 = r8_hyper_2f1 ( arg1, arg2, arg3, arg4 );
      weight_sum_exact = weight_sum_exact * ( 
        value1 / ( beta + 1.0 ) + value2 / ( alpha + 1.0 ) );
    }
    else if ( rule[dim] == 10 )
    {
      for ( i = 0; i < np[dim]; i++ )
      {
        alpha = p[p_index];
        p_index = p_index + 1;
      } 
      weight_sum_exact = 0.0;
    }
    else if ( rule[dim] == 11 )
    {
      weight_sum_exact = weight_sum_exact * 2.0;
    }
    else
    {
      printf ( "\n" );
      printf ( "SGMGA_PRODUCT_WEIGHT_TEST - Fatal error!\n" );
      printf ( "  Unexpected value of RULE[%d] = %d.\n", dim, rule[dim] );
      exit ( 1 );
    }
  }

  printf ( "\n" );
  printf ( "SGMGA_PRODUCT_WEIGHT_TEST:\n" );
  printf ( "  Compute the weights of a mixed factor product grid.\n" );
  if ( weight_sum_exact != 0.0 )
  {
    printf ( "\n" );
    printf ( "  As a simple test, sum these weights.\n" );
    printf ( "  They should sum to exactly %f\n", weight_sum_exact );
  }
  printf ( "\n" );
  printf ( "  Spatial dimension DIM_NUM = %d\n", dim_num );

  printf ( "\n" );
  printf ( " Dimension      Rule     Order        Parameters\n" );
  printf ( "\n" );

  p_index = 0;

  for ( dim = 0; dim < dim_num; dim++ )
  {
    if ( rule[dim] == 1 )
    {
      printf ( "  %8d  %8d  %8d\n", dim, rule[dim], order_1d[dim] );
    }
    else if ( rule[dim] == 2 )
    {
      printf ( "  %8d  %8d  %8d\n", dim, rule[dim], order_1d[dim] );
    }
    else if ( rule[dim] == 3 )
    {
      printf ( "  %8d  %8d  %8d\n", dim, rule[dim], order_1d[dim] );
    }
    else if ( rule[dim] == 4 )
    {
      printf ( "  %8d  %8d  %8d\n", dim, rule[dim], order_1d[dim] );
    }
    else if ( rule[dim] == 5 )
    {
      printf ( "  %8d  %8d  %8d\n", dim, rule[dim], order_1d[dim] );
    }
    else if ( rule[dim] == 6 )
    {
      alpha = p[p_index];
      p_index = p_index + 1;
      printf ( "  %8d  %8d  %8d  %14f\n", dim, rule[dim], order_1d[dim], alpha );
    }
    else if ( rule[dim] == 7 )
    {
      printf ( "  %8d  %8d  %8d\n", dim, rule[dim], order_1d[dim] );
    }
    else if ( rule[dim] == 8 )
    {
      alpha = p[p_index];
      p_index = p_index + 1;
      printf ( "  %8d  %8d  %8d  %14f\n", dim, rule[dim], order_1d[dim], alpha );
    }
    else if ( rule[dim] == 9 )
    {
      alpha = p[p_index];
      p_index = p_index + 1;
      beta = p[p_index];
      p_index = p_index + 1;
      printf ( "  %8d  %8d  %8d  %14f  %14f\n", dim, rule[dim], order_1d[dim], alpha, beta );
    }
    else if ( rule[dim] == 10 )
    {
      printf ( "  %8d  %8d  %8d  ", dim, rule[dim], order_1d[dim] );
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
      printf ( "  %8d  %8d  %8d\n", dim, rule[dim], order_1d[dim] );
    }
    else
    {
      printf ( "\n" );
      printf ( "SGMGA_PRODUCT_WEIGHT_TEST - Fatal error!\n" );
      printf ( "  Unexpected value of RULE = %d\n", rule[dim] );
      exit ( 1 );
    }
  }
/*
  Compute the weights.
*/
  weight = ( double * ) malloc ( order_nd * sizeof ( double ) );

  sgmga_product_weight ( dim_num, order_1d, order_nd, rule, 
    np, p, gw_compute_weights, weight );
/*
  Sum the weights to get the approximation to the integral of 1.
*/
  weight_sum = r8vec_sum ( order_nd, weight );
/*
  Compare the exact and estimated integrals.
*/
  weight_sum_error = r8_abs ( weight_sum - weight_sum_exact );

  if ( weight_sum_exact != 0.0 )
  {
    printf ( "\n" );
    printf ( "    Weight sum  Expected sum    Difference\n" );
    printf ( "\n" );
    printf ( "  %14f  %14f  %14e\n", weight_sum, weight_sum_exact, weight_sum_error );
  }
  else
  {
    printf ( "\n" );
    printf ( "    Weight sum\n" );
    printf ( "\n" );
    printf ( "  %14f\n", weight_sum );
  }
  free ( weight );

  return;
}

