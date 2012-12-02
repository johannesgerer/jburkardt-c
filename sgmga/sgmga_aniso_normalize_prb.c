# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "sandia_rules.h"
# include "sgmga.h"

int main ( void );
void sgmga_aniso_normalize_tests ( void );
void sgmga_aniso_normalize_test ( int dim_num, double level_weight[] );
void sgmga_importance_to_aniso_tests ( void );
void sgmga_importance_to_aniso_test ( int dim_num, double importance[], 
  double level_weight[] );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for SGMGA_ANISO_NORMALIZE_PRB.

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
  printf ( "SGMGA_ANISO_NORMALIZE_PRB:\n" );
  printf ( "  C version\n" );
  printf ( "  Test the SGMGA_ANISO_NORMALIZE and\n" );
  printf ( "  SGMGA_IMPORTANCE_TO_ANISO functions.\n" );

  sgmga_aniso_normalize_tests ( );

  sgmga_importance_to_aniso_tests ( );

  printf ( "\n" );
  printf ( "SGMGA_ANISO_NORMALIZE_PRB:\n" );
  printf ( "  Normal end of execution.\n" );

  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void sgmga_aniso_normalize_tests ( )

/******************************************************************************/
/*
  Purpose:

    SGMGA_ANISO_NORMALIZE_TESTS call SGMGA_ANISO_NORMALIZE_TEST.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 November 2009

  Author:

    John Burkardt
*/
{
  int dim_num;
  double *level_weight;

  printf ( "\n" );
  printf ( "SGMGA_ANISO_NORMALIZE_TESTS\n" );
  printf ( "  Call SGMGA_ANISO_NORMALIZE_TEST with various arguments.\n" );

  dim_num = 2;
  level_weight = ( double * ) malloc ( dim_num * sizeof ( double ) );
  level_weight[0] = 1.0;
  level_weight[1] = 1.0;
  sgmga_aniso_normalize_test ( dim_num, level_weight );
  free ( level_weight );

  dim_num = 2;
  level_weight = ( double * ) malloc ( dim_num * sizeof ( double ) );
  level_weight[0] = 10.0;
  level_weight[1] = 10.0;
  sgmga_aniso_normalize_test ( dim_num, level_weight );
  free ( level_weight );

  dim_num = 2;
  level_weight = ( double * ) malloc ( dim_num * sizeof ( double ) );
  level_weight[0] = 10.0;
  level_weight[1] = 2.0;
  sgmga_aniso_normalize_test ( dim_num, level_weight );
  free ( level_weight );

  dim_num = 2;
  level_weight = ( double * ) malloc ( dim_num * sizeof ( double ) );
  level_weight[0] = 1.0;
  level_weight[1] = 2.0;
  sgmga_aniso_normalize_test ( dim_num, level_weight );
  free ( level_weight );

  dim_num = 3;
  level_weight = ( double * ) malloc ( dim_num * sizeof ( double ) );
  level_weight[0] = 1.0;
  level_weight[1] = 2.0;
  level_weight[2] = 3.0;
  sgmga_aniso_normalize_test ( dim_num, level_weight );
  free ( level_weight );
/*
  Try a case in which one variable has 0 weight.
*/
  dim_num = 3;
  level_weight = ( double * ) malloc ( dim_num * sizeof ( double ) );
  level_weight[0] = 2.0;
  level_weight[1] = 0.0;
  level_weight[2] = 1.5;
  sgmga_aniso_normalize_test ( dim_num, level_weight );
  free ( level_weight );

  dim_num = 4;
  level_weight = ( double * ) malloc ( dim_num * sizeof ( double ) );
  level_weight[0] = 1.0;
  level_weight[1] = 2.0;
  level_weight[2] = 3.0;
  level_weight[3] = 4.0;
  sgmga_aniso_normalize_test ( dim_num, level_weight );
  free ( level_weight );

  return;
}
/******************************************************************************/

void sgmga_aniso_normalize_test ( int dim_num, double level_weight[] )

/******************************************************************************/
/*
  Purpose:

    SGMGA_ANISO_NORMALIZE_TEST calls SGMGA_ANISO_NORMALIZE.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    27 November 2009

  Author:

    John Burkardt
*/
{
  int dim;
  int option;

  printf ( "\n" );
  printf ( "SGMGA_ANISO_NORMALIZE_TEST\n" );
  printf ( "  Input weight sum: %f\n", r8vec_sum ( dim_num, level_weight ) );
  for ( dim = 0; dim < dim_num; dim++ )
  {
    printf ( "  %12f", level_weight[dim] );
  }
  printf ( "\n" );

  for ( option = 0; option <= 2; option++ )
  {
    sgmga_aniso_normalize ( option, dim_num, level_weight );

    printf ( "  For OPTION = %d,  Normalized weight sum: %f\n",  
      option, r8vec_sum ( dim_num, level_weight ) );
    for ( dim = 0; dim < dim_num; dim++ )
    {
      printf ( "  %12f", level_weight[dim] );
    }
    printf ( "\n" );
  }

  return;
}
/******************************************************************************/

void sgmga_importance_to_aniso_tests ( void )

/******************************************************************************/
/*
  Purpose:

    SGMGA_IMPORTANCE_TO_ANISO_TESTS call SGMGA_IMPORTANCE_TO_ANISO_TEST.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 November 2009

  Author:

    John Burkardt
*/
{
  int dim_num;
  double *importance;
  double *level_weight;

  printf ( "\n" );
  printf ( "SGMGA_IMPORTANCE_TO_ANISO_TESTS\n" );
  printf ( "  Call SGMGA_IMPORTANCE_TO_ANISO_TEST with various arguments.\n" );

  dim_num = 2;
  importance = ( double * ) malloc ( dim_num * sizeof ( double ) );
  level_weight = ( double * ) malloc ( dim_num * sizeof ( double ) );
  importance[0] = 1.0;
  importance[1] = 1.0;
  sgmga_importance_to_aniso_test ( dim_num, importance, level_weight );
  free ( importance );
  free ( level_weight );

  dim_num = 2;
  importance = ( double * ) malloc ( dim_num * sizeof ( double ) );
  level_weight = ( double * ) malloc ( dim_num * sizeof ( double ) );
  importance[0] = 10.0;
  importance[1] = 10.0;
  sgmga_importance_to_aniso_test ( dim_num, importance, level_weight );
  free ( importance );
  free ( level_weight );

  dim_num = 2;
  importance = ( double * ) malloc ( dim_num * sizeof ( double ) );
  level_weight = ( double * ) malloc ( dim_num * sizeof ( double ) );
  importance[0] = 10.0;
  importance[1] = 2.0;
  sgmga_importance_to_aniso_test ( dim_num, importance, level_weight );
  free ( importance );
  free ( level_weight );

  dim_num = 2;
  importance = ( double * ) malloc ( dim_num * sizeof ( double ) );
  level_weight = ( double * ) malloc ( dim_num * sizeof ( double ) );
  importance[0] = 1.0;
  importance[1] = 2.0;
  sgmga_importance_to_aniso_test ( dim_num, importance, level_weight );
  free ( importance );
  free ( level_weight );

  dim_num = 3;
  importance = ( double * ) malloc ( dim_num * sizeof ( double ) );
  level_weight = ( double * ) malloc ( dim_num * sizeof ( double ) );
  importance[0] = 1.0;
  importance[1] = 2.0;
  importance[2] = 3.0;
  sgmga_importance_to_aniso_test ( dim_num, importance, level_weight );
  free ( importance );
  free ( level_weight );
/*
  Try a case in which one variable has 0 importance.
*/
  dim_num = 3;
  importance = ( double * ) malloc ( dim_num * sizeof ( double ) );
  level_weight = ( double * ) malloc ( dim_num * sizeof ( double ) );
  importance[0] = 2.0;
  importance[1] = 0.0;
  importance[2] = 1.5;
  sgmga_importance_to_aniso_test ( dim_num, importance, level_weight );
  free ( importance );
  free ( level_weight );

  dim_num = 4;
  importance = ( double * ) malloc ( dim_num * sizeof ( double ) );
  level_weight = ( double * ) malloc ( dim_num * sizeof ( double ) );
  importance[0] = 1.0;
  importance[1] = 2.0;
  importance[2] = 3.0;
  importance[3] = 4.0;
  sgmga_importance_to_aniso_test ( dim_num, importance, level_weight );
  free ( importance );
  free ( level_weight );

  return;
}
/******************************************************************************/

void sgmga_importance_to_aniso_test ( int dim_num, double importance[], 
  double level_weight[] )

/******************************************************************************/
/*
  Purpose:

    SGMGA_IMPORTANCE_TO_ANISO_TEST calls SGMGA_IMPORTANCE_TO_ANISO.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    27 November 2009

  Author:

    John Burkardt
*/
{
  int dim;

  printf ( "\n" );
  printf ( "SGMGA_IMPORTANCE_TO_ANISO_TEST\n" );
  printf ( "  Importances:\n" );
  for ( dim = 0; dim < dim_num; dim++ )
  {
    printf ( "  %12f", importance[dim] );
  }
  printf ( "\n" );

  sgmga_importance_to_aniso ( dim_num, importance, level_weight );

  printf ( "  Anisotropic coefficients:\n" );
  for ( dim = 0; dim < dim_num; dim++ )
  {
    printf ( "  %12f", level_weight[dim] );
  }
  printf ( "\n" );

  return;
}
