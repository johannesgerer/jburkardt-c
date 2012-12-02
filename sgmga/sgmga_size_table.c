# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "sandia_rules.h"
# include "sgmga.h"


int main ( void );

void sgmga_size_table ( int rule_1d, int pd_1d, double p_1d[], int dim_min, 
  int dim_max, int level_max_min, int level_max_max, 
  void gw_compute_points_1d ( int order, int np, double p[], double x[] ) );

typedef void ( *GWPointer ) ( int order, int np, double p[], double w[] );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for SGMGA_SIZE_TABLE.

  Discussion:

    SGMGA_INDEX_PRB tests the SGMGA_INDEX function.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 December 2009

  Author:

    John Burkardt
*/
{
  double ctime;
  int dim_max;
  int dim_min;
  int level_max_max;
  int level_max_min;
  int np_1d;
  double *p_1d;
  int rule_1d;

  timestamp ( );
  printf ( "\n" );
  printf ( "SGMGA_SIZE_TABLE_TESTS\n" );
  printf ( "  C version\n" );
  printf ( "  Make tables of point counts.\n" );
  printf ( "  Print the CPU time required for each table.\n" );
/*
  For the Clenshaw-Curtis Grid (rule 1), compute some point counts
  for dimensions DIM_MIN to DIM_MAX, levels LEVE_MAX_MIN to LEVEL_MAX_MAX.
*/
  rule_1d = 1;
  np_1d = 0;
  p_1d = ( double * ) malloc ( np_1d * sizeof ( double ) );
  dim_min = 1;
  dim_max = 5;
  level_max_min = 0;
  level_max_max = 7;
  ctime = cpu_time ( );
  sgmga_size_table ( rule_1d, np_1d, p_1d, dim_min, dim_max, level_max_min, 
    level_max_max, clenshaw_curtis_compute_points_np );
  ctime = cpu_time ( ) - ctime;
  printf ( "\n" );
  printf ( "  CPU Time = %f\n", ctime );
  free ( p_1d );
/*
  Code is not efficient for this case.
*/
  if ( 0 )
  {
    rule_1d = 1;
    np_1d = 0;
    p_1d = ( double * ) malloc ( np_1d * sizeof ( double ) );
    dim_min = 6;
    dim_max = 10;
    level_max_min = 0;
    level_max_max = 7;
    ctime = cpu_time ( );
    sgmga_size_table ( rule_1d, np_1d, p_1d, dim_min, dim_max, level_max_min, 
      level_max_max, clenshaw_curtis_compute_points_np );
    ctime = cpu_time ( ) - ctime;
    printf ( "\n" );
    printf ( "  CPU Time = %f\n", ctime );
    free ( p_1d );
  }
/*
  Code is not efficient for this case.
*/
  if ( 0 ) 
  {
    rule_1d = 1;
    np_1d = 0;
    p_1d = ( double * ) malloc ( np_1d * sizeof ( double ) );
    dim_min = 100;
    dim_max = 100;
    level_max_min = 0;
    level_max_max = 2;
    ctime = cpu_time ( );
    sgmga_size_table ( rule_1d, np_1d, p_1d, dim_min, dim_max, level_max_min, 
      level_max_max, clenshaw_curtis_compute_points_np );
    ctime = cpu_time ( ) - ctime;
    printf ( "\n" );
    printf ( "  CPU Time = %f\n", ctime );
    free ( p_1d );
  }
/*
  For the Clenshaw-Curtis Grid with "slow exponential growth", (rule 11), 
  compute some point counts
  for dimensions DIM_MIN to DIM_MAX, levels LEVE_MAX_MIN to LEVEL_MAX_MAX.
*/
  rule_1d = 11;
  np_1d = 0;
  p_1d = ( double * ) malloc ( np_1d * sizeof ( double ) );
  dim_min = 1;
  dim_max = 5;
  level_max_min = 0;
  level_max_max = 7;
  ctime = cpu_time ( );
  sgmga_size_table ( rule_1d, np_1d, p_1d, dim_min, dim_max, level_max_min, 
    level_max_max, clenshaw_curtis_compute_points_np );
  ctime = cpu_time ( ) - ctime;
  printf ( "\n" );
  printf ( "  CPU Time = %f\n", ctime );
  free ( p_1d );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "SGMGA_SIZE_TABLE_TESTS\n" );
  printf ( "  Normal end of execution.\n" );

  printf ( "\n" );
  timestamp ( );
  
  return 0;
}
/******************************************************************************/

void sgmga_size_table ( int rule_1d, int np_1d, double p_1d[], int dim_min, 
  int dim_max, int level_max_min, int level_max_max,
  void gw_compute_points_1d ( int order, int np, double p[], double x[] ) )

/******************************************************************************/
/*
  Purpose:

    SGMGA_SIZE_TABLE tests SGMGA_SIZE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 December 2009

  Author:

    John Burkardt

  Parameters:

    Input, int RULE_1D, the 1D rule.

    Input, int NP_1D, the number of parameters in the 1D rule.

    Input, double P_1D[NP_1D], the parameters.

    Input, int DIM_MIN, the minimum spatial dimension to consider.

    Input, int DIM_MAX, the maximum spatial dimension to consider.

    Input, int LEVEL_MAX_MIN, the minimum value of LEVEL_MAX to consider.

    Input, int LEVEL_MAX_MAX, the maximum value of LEVEL_MAX to consider.

    Input, GW_COMPUTE_POINTS_1D ( int order, int np, double p[], double x[] ),
    a function which return the 1D quadrature points.
*/
{
  int dim;
  int dim_num;
  GWPointer *gw_compute_points;
  int i;
  int level_max;
  double *level_weight;
  int *np;
  int np_sum;
  double *p;
  int point_num;
  int *rule;
  double tol;

  printf ( "\n" );
  printf ( "SGMGA_SIZE_TABLE\n" );
  printf ( "  SGMGA_SIZE returns the number of distinct\n" );
  printf ( "  points in a sparse grid.\n" );
  printf ( "\n" );
  printf ( "  We use the same rule in all dimensions, and count the points,\n" );
  printf ( "  for a range of dimensions and levels.\n" );
  printf ( "\n" );
  printf ( "  1D rule index = %d\n", rule_1d );
  printf ( "\n" );

  tol = sqrt ( r8_epsilon ( ) );

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
      level_weight = ( double * ) malloc ( dim_num * sizeof ( double ) );
      rule = ( int * ) malloc ( dim_num * sizeof ( int ) );
      np = ( int * ) malloc ( dim_num * sizeof ( int ) );
      np_sum = dim_num * np_1d;
      p = ( double * ) malloc ( np_sum * sizeof ( double ) );
      gw_compute_points = ( GWPointer * ) malloc ( dim_num * sizeof ( GWPointer ) );

      for ( dim = 0; dim < dim_num; dim++ )
      {
        level_weight[dim] = 1.0;
        rule[dim] = rule_1d;
        np[dim] = np_1d;
        for ( i = 0; i < np_1d; i++ )
        {
          p[i+dim*np_1d] = p_1d[i];
        }
        gw_compute_points[i] = gw_compute_points_1d;
      }

      point_num = sgmga_size ( dim_num, level_weight, level_max, rule,
        np, p, gw_compute_points, tol, level_to_order_default );

      printf ( "  %8d", point_num );

      free ( gw_compute_points );
      free ( level_weight );
      free ( np );
      free ( p );
      free ( rule );
    }
    printf ( "\n" );
  }

  return;
}
