# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "sandia_rules.h"
# include "sgmga.h"

int main ( void );
void sgmga_vcn_tests ( void );
void sgmga_vcn_test ( int dim_num, double importance[], double level_weight[], 
  double q_min, double q_max );
void sgmga_vcn_ordered_tests ( void );
void sgmga_vcn_ordered_test ( int dim_num, double importance[], 
  double level_weight[], double q_min, double q_max );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for SGMGA_VCN_PRB.

  Discussion:

    SGMGA_VCN_PRB tests SGMGA_VCN.

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
  printf ( "SGMGA_VCN_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the SGMGA_VCN and SGMGA_VCN_ORDERED functions.\n" );

  sgmga_vcn_tests ( );

  sgmga_vcn_ordered_tests ( );

  printf ( "\n" );
  printf ( "SGMGA_VCN_PRB\n" );
  printf ( "  Normal end of execution.\n" );

  printf ( "\n" );
  timestamp ( );
  
  return 0;
}
/******************************************************************************/

void sgmga_vcn_tests ( void )

/******************************************************************************/
/*
  Purpose:

    SGMGA_VCN_TESTS calls SGMGA_VCN_TEST.

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
  int dim_num_array[12] = {
    2, 2, 2, 2, 2, 
    3, 3, 3, 3, 3, 
    4, 4 };
  double *importance;
  int level_max;
  int level_max_array[12] = {
    0, 1, 2, 3, 4, 
    0, 1, 2, 3, 4, 
    2, 3 };
  double *level_weight;
  double q_max;
  double q_min;
  int test;
  int test_num = 12;

  printf ( "\n" );
  printf ( "SGMGA_VCN_TESTS\n" );
  printf ( "  calls SGMGA_VCN_TEST.\n" );
/*
  Isotropic examples.
*/
  for ( test = 0; test < test_num; test++ )
  {
    dim_num = dim_num_array[test];
    importance = ( double * ) malloc ( dim_num * sizeof ( double ) );
    for ( dim = 0; dim < dim_num; dim++ )
    {
      importance[dim] = 1.0;
    }
    level_weight = ( double * ) malloc ( dim_num * sizeof ( double ) );
    sgmga_importance_to_aniso ( dim_num, importance, level_weight );
    level_max = level_max_array[test];
    q_min = ( double ) ( level_max ) - r8vec_sum ( dim_num, level_weight );
    q_max = ( double ) ( level_max );

    sgmga_vcn_test ( dim_num, importance, level_weight, q_min, q_max );

    free ( importance );
    free ( level_weight );
  }
/*
  Anisotropic examples.
*/
  for ( test = 0; test < test_num; test++ )
  {
    dim_num = dim_num_array[test];
    importance = ( double * ) malloc ( dim_num * sizeof ( double ) );
    for ( dim = 0; dim < dim_num; dim++ )
    {
      importance[dim] = ( double ) ( dim + 1 );
    }
    level_weight = ( double * ) malloc ( dim_num * sizeof ( double ) );
    sgmga_importance_to_aniso ( dim_num, importance, level_weight );
    level_max = level_max_array[test];
    q_min = ( double ) ( level_max ) - r8vec_sum ( dim_num, level_weight );
    q_max = ( double ) ( level_max );

    sgmga_vcn_test ( dim_num, importance, level_weight, q_min, q_max );

    free ( importance );
    free ( level_weight );
  }

  return;
}
/******************************************************************************/

void sgmga_vcn_test ( int dim_num, double importance[], double level_weight[], 
  double q_min, double q_max )

/******************************************************************************/
/*
  Purpose:

    SGMGA_VCN_TEST tests SGMGA_VCN.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 November 2009

  Author:

    John Burkardt
*/
{
  int dim;
  int i;
  int *level_1d;
  int *level_1d_max;
  double level_weight_norm;
  int more_grids;
  double q;
  int test;

  level_1d = ( int * ) malloc ( dim_num * sizeof ( int ) );
  level_1d_max = ( int * ) malloc ( dim_num * sizeof ( int ) );

  printf ( "\n" );
  printf ( "SGMGA_VCN_TEST\n" );
  printf ( "  SGMGA_VCN considers vectors:\n" );
  printf ( "    0 <= LEVEL_1D(1:N) <= LEVEL_1D_MAX(1:N),\n" );
  printf ( "  Set\n" );
  printf ( "    Q = sum ( LEVEL_WEIGHT(1:N) * LEVEL_1D(1:N) )\n" );
  printf ( "  Accept only vectors for which:\n" );
  printf ( "    Q_MIN < Q <= Q_MAX\n" );
  printf ( "\n" );
  printf ( "  No attempt is made to order the LEVEL_1D values.\n" );

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
  printf ( "  Q_MIN:    %14f\n", q_min );
  printf ( "  Q_MAX:    %14f\n", q_max );
  printf ( "  LEVEL_1D_MAX:           " );
  for ( dim = 0; dim < dim_num; dim++ )
  {
    printf ( "  %8d", level_1d_max[dim] );
  }
  printf ( "\n" );
  printf ( "\n" );

  i = 0;

  for ( ; ; )
  {
    sgmga_vcn ( dim_num, level_weight, level_1d_max, level_1d, q_min, q_max, 
      &more_grids );

    if ( !more_grids )
    {
      printf ( "\n" );
      printf ( "  End of solutions.\n" );
      break;
    }

    q = 0.0;
    for ( dim = 0; dim < dim_num; dim++ )
    {
      q = q + level_weight[dim] * ( double ) level_1d[dim];
    }
    i = i + 1;
    printf ( "  %8d  %14f", i, q );
    for ( dim = 0; dim < dim_num; dim++ )
    {
      printf ( "  %8d", level_1d[dim] );
    }
    printf ( "\n" );
  }

  free ( level_1d );
  free ( level_1d_max );

  return;
}
/******************************************************************************/

void sgmga_vcn_ordered_tests ( )

/******************************************************************************/
/*
  Purpose:

    SGMGA_VCN_ORDERED_TESTS calls SGMGA_VCN_ORDERED_TEST.

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
  int dim_num_array[12] = {
    2, 2, 2, 2, 2, 
    3, 3, 3, 3, 3, 
    4, 4 };
  double *importance;
  int level_max;
  int level_max_array[12] = {
    0, 1, 2, 3, 4, 
    0, 1, 2, 3, 4, 
    2, 3 };
  double *level_weight;
  double q_max;
  double q_min;
  int test;
  int test_num = 12;

  printf ( "\n" );
  printf ( "SGMGA_VCN_ORDERED_TESTS\n" );
  printf ( "  calls SGMGA_VCN_ORDERED_TEST.\n" );
/*
  Isotropic examples.
*/
  for ( test = 0; test < test_num; test++ )
  {
    dim_num = dim_num_array[test];
    importance = ( double * ) malloc ( dim_num * sizeof ( double ) );
    for ( dim = 0; dim < dim_num; dim++ )
    {
      importance[dim] = 1.0;
    }
    level_weight = ( double * ) malloc ( dim_num * sizeof ( double ) );
    sgmga_importance_to_aniso ( dim_num, importance, level_weight );
    level_max = level_max_array[test];
    q_min = ( double ) ( level_max ) - r8vec_sum ( dim_num, level_weight );
    q_max = ( double ) ( level_max );

    sgmga_vcn_ordered_test ( dim_num, importance, level_weight, q_min, q_max );

    free ( importance );
    free ( level_weight );
  }
/*
  Anisotropic examples.
*/
  for ( test = 0; test < test_num; test++ )
  {
    dim_num = dim_num_array[test];
    importance = ( double * ) malloc ( dim_num * sizeof ( double ) );
    for ( dim = 0; dim < dim_num; dim++ )
    {
      importance[dim] = ( double ) ( dim + 1 );
    }
    level_weight = ( double * ) malloc ( dim_num * sizeof ( double ) );
    sgmga_importance_to_aniso ( dim_num, importance, level_weight );
    level_max = level_max_array[test];
    q_min = ( double ) ( level_max ) - r8vec_sum ( dim_num, level_weight );
    q_max = ( double ) ( level_max );

    sgmga_vcn_ordered_test ( dim_num, importance, level_weight, q_min, q_max );

    free ( importance );
    free ( level_weight );
  }

  return;
}
/******************************************************************************/

void sgmga_vcn_ordered_test ( int dim_num, double importance[], 
  double level_weight[], double q_min, double q_max )

/******************************************************************************/
/*
  Purpose:

    SGMGA_VCN_ORDERED_TEST tests SGMGA_VCN_ORDERED.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 November 2009

  Author:

    John Burkardt
*/
{
  int dim;
  int i;
  int *level_1d;
  int *level_1d_max;
  double level_weight_norm;
  int more_grids;
  double q;
  int test;

  level_1d = ( int * ) malloc ( dim_num * sizeof ( int ) );
  level_1d_max = ( int * ) malloc ( dim_num * sizeof ( int ) );

  printf ( "\n" );
  printf ( "SGMGA_VCN_ORDERED_TEST\n" );
  printf ( "  SGMGA_VCN_ORDERED considers vectors:\n" );
  printf ( "    0 <= LEVEL_1D(1:N) <= LEVEL_1D_MAX(1:N),\n" );
  printf ( "  Set\n" );
  printf ( "    Q = sum ( LEVEL_WEIGHT(1:N) * LEVEL_1D(1:N) )\n" );
  printf ( "  Accept only vectors for which:\n" );
  printf ( "    Q_MIN < Q <= Q_MAX\n" );
  printf ( "\n" );
  printf ( "  The solutions are weakly ordered by the value of Q.\n" );

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
  printf ( "  Q_MIN:    %14f\n", q_min );
  printf ( "  Q_MAX:    %14f\n", q_max );
  printf ( "  LEVEL_1D_MAX:           " );
  for ( dim = 0; dim < dim_num; dim++ )
  {
    printf ( "  %8d", level_1d_max[dim] );
  }
  printf ( "\n" );
  printf ( "\n" );

  i = 0;

  for ( ; ; )
  {
    sgmga_vcn_ordered ( dim_num, level_weight, level_1d_max, level_1d, 
      q_min, q_max, &more_grids );

    if ( !more_grids )
    {
      printf ( "\n" );
      printf ( "  End of solutions.\n" );
      break;
    }

    q = 0.0;
    for ( dim = 0; dim < dim_num; dim++ )
    {
      q = q + level_weight[dim] * ( double ) level_1d[dim];
    }
    i = i + 1;
    printf ( "  %8d  %14f", i, q );
    for ( dim = 0; dim < dim_num; dim++ )
    {
      printf ( "  %8d", level_1d[dim] );
    }
    printf ( "\n" );
  }

  free ( level_1d );
  free ( level_1d_max );

  return;
}

