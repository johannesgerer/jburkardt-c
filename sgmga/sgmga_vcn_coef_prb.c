# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "sandia_rules.h"
# include "sgmga.h"

int main ( void );

void sgmga_vcn_coef_tests ( void );
void sgmga_vcn_coef_test ( int dim_num, double importance[], 
  double level_weight[], int level_max_min, int level_max_max );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for SGMGA_VCN_COEF_PRB.

  Discussion:

    SGMGA_VCN_COEF_PRB tests the SGMGA_VCN_COEF function.

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
  printf ( "SGMGA_VCN_COEF_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the SGMGA_VCN_COEF function.\n" );
/*
  Compute examples of the combinatorial coefficent.
*/
  sgmga_vcn_coef_tests ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "SGMGA_VCN_COEF_PRB\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );
  
  return 0;
}
/******************************************************************************/

void sgmga_vcn_coef_tests ( void )

/******************************************************************************/
/*
  Purpose:

    SGMGA_VCN_COEF_TESTS calls SGMGA_VCN_COEF_TEST.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    27 November 2009

  Author:

    John Burkardt
*/
{
  int dim;
  int dim_num;
  double *importance;
  int level_max;
  int level_max_max;
  int level_max_min;
  double *level_weight;

  printf ( "\n" );
  printf ( "SGMGA_VCN_COEF_TESTS\n" );
  printf ( "  calls SGMGA_VCN_COEF_TEST.\n" );

  dim_num = 2;
  importance = ( double * ) malloc ( dim_num * sizeof ( double ) );
  for ( dim = 0; dim < dim_num; dim++ )
  {
    importance[dim] = 1.0;
  }
  level_weight = ( double * ) malloc ( dim_num * sizeof ( double ) );
  sgmga_importance_to_aniso ( dim_num, importance, level_weight );
  level_max_min = 0;
  level_max_max = 4;
  sgmga_vcn_coef_test ( dim_num, importance, level_weight, level_max_min, 
    level_max_max );
  free ( importance );
  free ( level_weight );

  dim_num = 2;
  importance = ( double * ) malloc ( dim_num * sizeof ( double ) );
  for ( dim = 0; dim < dim_num; dim++ )
  {
    importance[dim] = ( double ) ( dim + 1 );
  }
  level_weight = ( double * ) malloc ( dim_num * sizeof ( double ) );
  sgmga_importance_to_aniso ( dim_num, importance, level_weight );
  level_max_min = 0;
  level_max_max = 4;
  sgmga_vcn_coef_test ( dim_num, importance, level_weight, level_max_min, 
    level_max_max );
  free ( importance );
  free ( level_weight );

  dim_num = 3;
  importance = ( double * ) malloc ( dim_num * sizeof ( double ) );
  for ( dim = 0; dim < dim_num; dim++ )
  {
    importance[dim] = 1.0;
  }
  level_weight = ( double * ) malloc ( dim_num * sizeof ( double ) );
  sgmga_importance_to_aniso ( dim_num, importance, level_weight );
  level_max_min = 0;
  level_max_max = 4;
  sgmga_vcn_coef_test ( dim_num, importance, level_weight, level_max_min, 
    level_max_max );
  free ( importance );
  free ( level_weight );

  dim_num = 3;
  importance = ( double * ) malloc ( dim_num * sizeof ( double ) );
  for ( dim = 0; dim < dim_num; dim++ )
  {
    importance[dim] = ( double ) ( dim + 1 );
  }
  level_weight = ( double * ) malloc ( dim_num * sizeof ( double ) );
  sgmga_importance_to_aniso ( dim_num, importance, level_weight );
  level_max_min = 0;
  level_max_max = 4;
  sgmga_vcn_coef_test ( dim_num, importance, level_weight, level_max_min, 
    level_max_max );
  free ( importance );
  free ( level_weight );

  dim_num = 4;
  importance = ( double * ) malloc ( dim_num * sizeof ( double ) );
  for ( dim = 0; dim < dim_num; dim++ )
  {
    importance[dim] = ( double ) ( dim + 1 );
  }
  level_weight = ( double * ) malloc ( dim_num * sizeof ( double ) );
  sgmga_importance_to_aniso ( dim_num, importance, level_weight );
  level_max_min = 0;
  level_max_max = 3;
  sgmga_vcn_coef_test ( dim_num, importance, level_weight, level_max_min, 
    level_max_max );
  free ( importance );
  free ( level_weight );
/*
  Try a case with a dimension of "0 importance".
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
  sgmga_vcn_coef_test ( dim_num, importance, level_weight, level_max_min, 
    level_max_max );
  free ( importance );
  free ( level_weight );

  return;
}
/******************************************************************************/

void sgmga_vcn_coef_test ( int dim_num, double importance[], 
  double level_weight[], int level_max_min, int level_max_max )

/******************************************************************************/
/*
  Purpose:

    SGMGA_VCN_COEF_TEST tests SGMGA_VCN_COEF.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 November 2009

  Author:

    John Burkardt
*/
{
  double coef;
  double coef_sum;
  int dim;
  int i;
  int *level_1d;
  int *level_1d_max;
  int level_max;
  double level_weight_min_pos;
  double level_weight_norm;
  int more_grids;
  double q;
  double q_max;
  double q_min;

  level_1d = ( int * ) malloc ( dim_num * sizeof ( int ) );
  level_1d_max = ( int * ) malloc ( dim_num * sizeof ( int ) );

  printf ( "\n" );
  printf ( "SGMGA_VCN_COEF_TEST\n" );
  printf ( "  SGMGA_VCN_COEF computes a \"combinatorial coefficennt\"\n" );
  printf ( "  for vectors that are produced by SGMGA_VCN.\n" );
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

  for ( level_max = level_max_min; level_max <= level_max_max; level_max++ )
  {
    i = 0;
    coef_sum = 0.0;
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

    printf ( "\n" );
    printf ( "  Q_MIN:    %14f\n", q_min );
    printf ( "  Q_MAX:    %14f\n", q_max );
    printf ( "  LEVEL_1D_MAX:                           " );
    for ( dim = 0; dim < dim_num; dim++ )
    {
      printf ( "  %8d", level_1d_max[dim] );
    }
    printf ( "\n" );
    printf ( "\n" );
    printf ( "         I               Q            Coef         X\n" );
    printf ( "\n" );
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
      coef = sgmga_vcn_coef ( dim_num, level_weight, level_1d_max, 
        level_1d, q_min, q_max );

      i = i + 1;

      q = 0.0;
      for ( dim = 0; dim < dim_num; dim++ )
      {
        q = q + level_weight[dim] * ( double ) level_1d[dim];
      }

      coef_sum = coef_sum + coef;

      printf ( "  %8d  %14f  %14f", i, q, coef );
      for ( dim = 0; dim < dim_num; dim++ )
      {
        printf ( "  %8d", level_1d[dim] );
      }
      printf ( "\n" );
    }
    printf ( "\n" );
    printf ( "  Sum of COEFs =            %14f\n", coef_sum );
  }
  free ( level_1d );
  free ( level_1d_max );

  return;
}

