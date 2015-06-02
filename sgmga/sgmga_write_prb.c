# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <string.h>

# include "sandia_rules.h"
# include "sgmga.h"

int main ( void );
void sgmga_write_tests ( void );
void sgmga_write_test ( int dim_num, double level_weight[], int level_max, 
  int rule[], int np[], double p[], 
  void ( *gw_compute_points[] ) ( int order, int np, double p[], double w[] ),
  void ( *gw_compute_weights[] ) ( int order, int np, double p[], double w[] ),
  double tol, char *file_name );

typedef void ( *GWPointer ) ( int order, int np, double p[], double w[] );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for SGMGA_WRITE_PRB.

  Discussion:

    SGMGA_WRITE_PRB tests the SGMGA_WRITE function.

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
  printf ( "SGMGA_WRITE_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the SGMGA_WRITE function.\n" );
/*
  Generate sparse grid rules and write them to files.
*/
  sgmga_write_tests ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "SGMGA_WRITE_PRB\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );
  
  return 0;
}
/******************************************************************************/

void sgmga_write_tests ( void )

/******************************************************************************/
/*
  Purpose:

    SGMGA_WRITE_TESTS calls SGMGA_WRITE_TEST.

  Discussion:
  
    We can't test Golub-Welsch rules in this routine, because the program
    that writes out the files needs to know the integration region for each
    component, and we have not specified how that would be done with 
    Golub Welsch rules.

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
  char file_name[255];
  GWPointer *gw_compute_points;
  GWPointer *gw_compute_weights;
  double *importance;
  int level_max;
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
  printf ( "SGMGA_WRITE_TESTS\n" );
  printf ( "  Call SGMGA_WRITE_TEST with various arguments.\n" );
/*
  Set the point equality tolerance.
*/
  tol = sqrt ( r8_epsilon ( ) );
  printf ( "\n" );
  printf ( "  All tests will use a point equality tolerance of %e\n\n", tol );

  dim_num = 2;
  importance = ( double * ) malloc ( dim_num * sizeof ( double ) );
  for ( dim = 0; dim < dim_num; dim++ )
  {
    importance[dim] = 1.0;
  }
  level_weight = ( double * ) malloc ( dim_num * sizeof ( double ) );
  sgmga_importance_to_aniso ( dim_num, importance, level_weight );
  level_max = 2;
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
  gw_compute_weights = ( GWPointer * ) malloc ( dim_num * sizeof ( GWPointer ) );
  gw_compute_weights[0] = clenshaw_curtis_compute_weights_np;
  gw_compute_weights[1] = clenshaw_curtis_compute_weights_np;
  strcpy ( file_name, "sgmga_d2_l2_ccxcc_iso" );
  sgmga_write_test ( dim_num, level_weight, level_max, rule, np, 
    p, gw_compute_points, gw_compute_weights, tol, file_name );
  free ( gw_compute_points );
  free ( gw_compute_weights );
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
  level_max = 2;
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
  gw_compute_weights = ( GWPointer * ) malloc ( dim_num * sizeof ( GWPointer ) );
  gw_compute_weights[0] = clenshaw_curtis_compute_weights_np;
  gw_compute_weights[1] = clenshaw_curtis_compute_weights_np;
  strcpy ( file_name, "sgmga_d2_l2_ccxcc_aniso" );
  sgmga_write_test ( dim_num, level_weight, level_max, rule, np, 
    p, gw_compute_points, gw_compute_weights, tol, file_name );
  free ( gw_compute_points );
  free ( gw_compute_weights );
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
  level_max = 2;
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
  gw_compute_weights = ( GWPointer * ) malloc ( dim_num * sizeof ( GWPointer ) );
  gw_compute_weights[0] = clenshaw_curtis_compute_weights_np;
  gw_compute_weights[1] = clenshaw_curtis_compute_weights_np;
  gw_compute_weights[2] = clenshaw_curtis_compute_weights_np;
  strcpy ( file_name, "sgmga_d3_l2_ccxccxcc_iso" );
  sgmga_write_test ( dim_num, level_weight, level_max, rule, np, 
    p, gw_compute_points, gw_compute_weights, tol, file_name );
  free ( gw_compute_points );
  free ( gw_compute_weights );
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
  level_max = 2;
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
  gw_compute_weights = ( GWPointer * ) malloc ( dim_num * sizeof ( GWPointer ) );
  gw_compute_weights[0] = clenshaw_curtis_compute_weights_np;
  gw_compute_weights[1] = clenshaw_curtis_compute_weights_np;
  gw_compute_weights[2] = clenshaw_curtis_compute_weights_np;
  strcpy ( file_name, "sgmga_d3_l2_ccxccxcc_aniso" );
  sgmga_write_test ( dim_num, level_weight, level_max, rule, np, 
    p, gw_compute_points, gw_compute_weights, tol, file_name );
  free ( gw_compute_points );
  free ( gw_compute_weights );
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
  level_max = 3;
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
  gw_compute_weights = ( GWPointer * ) malloc ( dim_num * sizeof ( GWPointer ) );
  gw_compute_weights[0] = clenshaw_curtis_compute_weights_np;
  gw_compute_weights[1] = patterson_lookup_weights_np;
  strcpy ( file_name, "sgmga_d2_l3_ccxgp_aniso" );
  sgmga_write_test ( dim_num, level_weight, level_max, rule, np,
    p, gw_compute_points, gw_compute_weights, tol, file_name );
  free ( gw_compute_points );
  free ( gw_compute_weights );
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
  level_max = 2;
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
  gw_compute_weights = ( GWPointer * ) malloc ( dim_num * sizeof ( GWPointer ) );
  gw_compute_weights[0] = clenshaw_curtis_compute_weights_np;
  gw_compute_weights[1] = legendre_compute_weights_np;
  strcpy ( file_name, "sgmga_d2_l2_ccxgl_aniso" );
  sgmga_write_test ( dim_num, level_weight, level_max, rule, np, 
    p, gw_compute_points, gw_compute_weights, tol, file_name );
  free ( gw_compute_points );
  free ( gw_compute_weights );
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
  level_max = 2;
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
  gw_compute_weights = ( GWPointer * ) malloc ( dim_num * sizeof ( GWPointer ) );
  gw_compute_weights[0] = clenshaw_curtis_compute_weights_np;
  gw_compute_weights[1] = laguerre_compute_weights_np;
  strcpy ( file_name, "sgmga_d2_l2_ccxlg_aniso" );
  sgmga_write_test ( dim_num, level_weight, level_max, rule, np, 
    p, gw_compute_points, gw_compute_weights, tol, file_name );
  free ( gw_compute_points );
  free ( gw_compute_weights );
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
  level_max = 2;
  rule = ( int * ) malloc ( dim_num * sizeof ( int ) );
  rule[0] = 1;
  rule[1] = 8;
  np = ( int * ) malloc ( dim_num * sizeof ( int ) );
  np[0] = 1;
  np[1] = 0;
  np_sum = i4vec_sum ( dim_num, np );
  p = ( double * ) malloc ( np_sum * sizeof ( double ) );
  p[0] = 1.5;
  gw_compute_points = ( GWPointer * ) malloc ( dim_num * sizeof ( GWPointer ) );
  gw_compute_points[0] = clenshaw_curtis_compute_points_np;
  gw_compute_points[1] = gen_laguerre_compute_points_np;
  gw_compute_weights = ( GWPointer * ) malloc ( dim_num * sizeof ( GWPointer ) );
  gw_compute_weights[0] = clenshaw_curtis_compute_weights_np;
  gw_compute_weights[1] = gen_laguerre_compute_weights_np;
  strcpy ( file_name, "sgmga_d2_l2_ccxglg_aniso" );
  sgmga_write_test ( dim_num, level_weight, level_max, rule, np, 
    p, gw_compute_points, gw_compute_weights, tol, file_name );
  free ( gw_compute_points );
  free ( gw_compute_weights );
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
  level_max = 2;
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
  gw_compute_weights = ( GWPointer * ) malloc ( dim_num * sizeof ( GWPointer ) );
  gw_compute_weights[0] = fejer2_compute_weights_np;
  gw_compute_weights[1] = jacobi_compute_weights_np;
  strcpy ( file_name, "sgmga_d2_l2_f2xgj_aniso" );
  sgmga_write_test ( dim_num, level_weight, level_max, rule, np, 
    p, gw_compute_points, gw_compute_weights, tol, file_name );
  free ( gw_compute_points );
  free ( gw_compute_weights );
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
  level_max = 2;
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
  gw_compute_weights = ( GWPointer * ) malloc ( dim_num * sizeof ( GWPointer ) );
  gw_compute_weights[0] = gen_hermite_compute_weights_np;
  gw_compute_weights[1] = legendre_compute_weights_np;
  strcpy ( file_name, "sgmga_d2_l2_gghxgl_aniso" );
  sgmga_write_test ( dim_num, level_weight, level_max, rule, np, 
    p, gw_compute_points, gw_compute_weights, tol, file_name );
  free ( gw_compute_points );
  free ( gw_compute_weights );
  free ( importance );
  free ( level_weight );
  free ( np );
  free ( p );
  free ( rule );
/*
  LEVEL_MAX = 1
*/
  dim_num = 2;
  importance = ( double * ) malloc ( dim_num * sizeof ( double ) );
  for ( dim = 0; dim < dim_num; dim++ )
  {
    importance[dim] = ( double ) ( dim + 1 );
  }
  level_weight = ( double * ) malloc ( dim_num * sizeof ( double ) );
  sgmga_importance_to_aniso ( dim_num, importance, level_weight );
  level_max = 1;
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
  gw_compute_weights = ( GWPointer * ) malloc ( dim_num * sizeof ( GWPointer ) );
  gw_compute_weights[0] = clenshaw_curtis_compute_weights_np;
  gw_compute_weights[1] = clenshaw_curtis_compute_weights_np;
  strcpy ( file_name, "sgmga_d2_l1_ccxcc_aniso" );
  sgmga_write_test ( dim_num, level_weight, level_max, rule, np, 
    p, gw_compute_points, gw_compute_weights, tol, file_name );
  free ( gw_compute_points );
  free ( gw_compute_weights );
  free ( importance );
  free ( level_weight );
  free ( np );
  free ( p );
  free ( rule );
/*
  LEVEL_MAX = 2 (already done)

  LEVEL_MAX = 3
*/
  dim_num = 2;
  importance = ( double * ) malloc ( dim_num * sizeof ( double ) );
  for ( dim = 0; dim < dim_num; dim++ )
  {
    importance[dim] = ( double ) ( dim + 1 );
  }
  level_weight = ( double * ) malloc ( dim_num * sizeof ( double ) );
  sgmga_importance_to_aniso ( dim_num, importance, level_weight );
  level_max = 3;
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
  gw_compute_weights = ( GWPointer * ) malloc ( dim_num * sizeof ( GWPointer ) );
  gw_compute_weights[0] = clenshaw_curtis_compute_weights_np;
  gw_compute_weights[1] = clenshaw_curtis_compute_weights_np;
  strcpy ( file_name, "sgmga_d2_l3_ccxcc_aniso" );
  sgmga_write_test ( dim_num, level_weight, level_max, rule, np, 
    p, gw_compute_points, gw_compute_weights, tol, file_name );
  free ( gw_compute_points );
  free ( gw_compute_weights );
  free ( importance );
  free ( level_weight );
  free ( np );
  free ( p );
  free ( rule );
/*
  LEVEL_MAX = 4
*/
  dim_num = 2;
  importance = ( double * ) malloc ( dim_num * sizeof ( double ) );
  for ( dim = 0; dim < dim_num; dim++ )
  {
    importance[dim] = ( double ) ( dim + 1 );
  }
  level_weight = ( double * ) malloc ( dim_num * sizeof ( double ) );
  sgmga_importance_to_aniso ( dim_num, importance, level_weight );
  level_max = 4;
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
  gw_compute_weights = ( GWPointer * ) malloc ( dim_num * sizeof ( GWPointer ) );
  gw_compute_weights[0] = clenshaw_curtis_compute_weights_np;
  gw_compute_weights[1] = clenshaw_curtis_compute_weights_np;
  strcpy ( file_name, "sgmga_d2_l4_ccxcc_aniso" );
  sgmga_write_test ( dim_num, level_weight, level_max, rule, np, 
    p, gw_compute_points, gw_compute_weights, tol, file_name );
  free ( gw_compute_points );
  free ( gw_compute_weights );
  free ( importance );
  free ( level_weight );
  free ( np );
  free ( p );
  free ( rule );
/*
  LEVEL_MAX = 5
*/
  dim_num = 2;
  importance = ( double * ) malloc ( dim_num * sizeof ( double ) );
  for ( dim = 0; dim < dim_num; dim++ )
  {
    importance[dim] = ( double ) ( dim + 1 );
  }
  level_weight = ( double * ) malloc ( dim_num * sizeof ( double ) );
  sgmga_importance_to_aniso ( dim_num, importance, level_weight );
  level_max = 5;
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
  gw_compute_weights = ( GWPointer * ) malloc ( dim_num * sizeof ( GWPointer ) );
  gw_compute_weights[0] = clenshaw_curtis_compute_weights_np;
  gw_compute_weights[1] = clenshaw_curtis_compute_weights_np;
  strcpy ( file_name, "sgmga_d2_l5_ccxcc_aniso" );
  sgmga_write_test ( dim_num, level_weight, level_max, rule, np, 
    p, gw_compute_points, gw_compute_weights, tol, file_name );
  free ( gw_compute_points );
  free ( gw_compute_weights );
  free ( importance );
  free ( level_weight );
  free ( np );
  free ( p );
  free ( rule );
/*
  Dimension 3
*/
  dim_num = 3;
  importance = ( double * ) malloc ( dim_num * sizeof ( double ) );
  for ( dim = 0; dim < dim_num; dim++ )
  {
    importance[dim] = ( double ) ( dim + 1 );
  }
  level_weight = ( double * ) malloc ( dim_num * sizeof ( double ) );
  sgmga_importance_to_aniso ( dim_num, importance, level_weight );
  level_max = 2;
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
  gw_compute_weights = ( GWPointer * ) malloc ( dim_num * sizeof ( GWPointer ) );
  gw_compute_weights[0] = clenshaw_curtis_compute_weights_np;
  gw_compute_weights[1] = legendre_compute_weights_np;
  gw_compute_weights[2] = hermite_compute_weights_np;
  strcpy ( file_name, "sgmga_d3_l2_ccxglxgh_aniso" );
  sgmga_write_test ( dim_num, level_weight, level_max, rule, np, 
    p, gw_compute_points, gw_compute_weights, tol, file_name );
  free ( gw_compute_points );
  free ( gw_compute_weights );
  free ( importance );
  free ( level_weight );
  free ( np );
  free ( p );
  free ( rule );

  return;
}
/******************************************************************************/

void sgmga_write_test ( int dim_num, double level_weight[], int level_max, 
  int rule[], int np[], double p[], 
  void ( *gw_compute_points[] ) ( int order, int np, double p[], double w[] ),
  void ( *gw_compute_weights[] ) ( int order, int np, double p[], double w[] ),
  double tol, char *file_name )

/******************************************************************************/
/*
  Purpose:

    SGMGA_WRITE_TEST tests SGMGA_WRITE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    18 November 2009

  Author:

    John Burkardt

  Parameters:

    Input, integer DIM_NUM, the spatial dimension.

    Input, double LEVEL_WEIGHT[DIM_NUM], the weights for each dimension.

    Input, integer LEVEL_MAX, the level that defines the grid.

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

    Input, void ( *GW_COMPUTE_POINTS[] ) ( int order, int np, double p[], double x[] ),
    an array of pointers to functions which return the 1D quadrature points 
    associated with each spatial dimension for which a Golub Welsch rule 
    is used.

    Input, void ( *GW_COMPUTE_WEIGHTS[] ) ( int order, int np, double p[], double w[] ),
    an array of pointers to functions which return the 1D quadrature weights 
    associated with each spatial dimension for which a Golub Welsch rule 
    is used.

    Input, double TOL, a tolerance for point equality.

    Input, char *FILE_NAME, the main name of the output files.
*/
{
  int point_num;
  int point_total_num;
  int *sparse_index;
  int *sparse_order;
  double *sparse_point;
  int *sparse_unique_index;
  double *sparse_weight;

  printf ( "\n" );
  printf ( "SGMGA_WRITE_TEST\n" );
  printf ( "  SGMGA_WRITE writes a sparse grid rule to files.\n" );
/*
  Compute necessary data.
*/
  point_total_num = sgmga_size_total ( dim_num, level_weight,
    level_max, rule, level_to_order_default );

  point_num = sgmga_size ( dim_num, level_weight, level_max, 
    rule, np, p, gw_compute_points, tol, 
    level_to_order_default );

  sparse_unique_index = ( int * ) malloc ( point_total_num * sizeof ( int ) );

  sgmga_unique_index ( dim_num, level_weight, level_max, rule, 
    np, p, gw_compute_points, tol, point_num, point_total_num, 
    level_to_order_default, sparse_unique_index );

  sparse_order = ( int * ) malloc ( dim_num * point_num * sizeof ( int ) );
  sparse_index = ( int * ) malloc ( dim_num * point_num * sizeof ( int ) );

  sgmga_index ( dim_num, level_weight, level_max, rule, point_num, 
    point_total_num, sparse_unique_index, level_to_order_default,
    sparse_order, sparse_index );
/*
  Compute points and weights.
*/
  sparse_point = ( double * ) malloc ( dim_num * point_num * sizeof ( double ) );

  sgmga_point ( dim_num, level_weight, level_max, rule, np, 
    p, gw_compute_points, point_num, sparse_order, sparse_index, 
    level_to_order_default, sparse_point );

  sparse_weight = ( double * ) malloc ( point_num * sizeof ( double ) );

  sgmga_weight ( dim_num, level_weight, level_max, rule, np, 
    p, gw_compute_weights, point_num, point_total_num, sparse_unique_index, 
    level_to_order_default, sparse_weight );
/*
  Write points and weights to files.
*/
  sgmga_write ( dim_num, level_weight, rule, np, p, 
    point_num, sparse_weight, sparse_point, file_name );

  free ( sparse_index );
  free ( sparse_order );
  free ( sparse_point );
  free ( sparse_unique_index );
  free ( sparse_weight );

  return;
}
