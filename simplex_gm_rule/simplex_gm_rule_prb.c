# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "simplex_gm_rule.h"

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );
void test04 ( );
void test05 ( );
void test06 ( );
void test07 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for SIMPLEX_GM_RULE_PRB.

  Discussion:

    SIMPLEX_GM_RULE_PRB tests the SIMPLEX_GM_RULE library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 August 2012

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "SIMPLEX_GM_RULE_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the SIMPLEX_GM_RULE library.\n" );

  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
  test05 ( );
  test06 ( );
  test07 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "SIMPLEX_GM_RULE_PRB\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void test01 ( )

/******************************************************************************/
/*
  Purpose:

    TEST01 tests SIMPLEX_UNIT_TO_GENERAL.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 August 2012

  Author:

    John Burkardt
*/
{
# define DIM_NUM 2

  int dim;
  int dim_num = DIM_NUM;
  int j;
  double *phy;
  double *phy_unit;
  int point_num = 10;
  double *ref;
  int seed = 123456789;
  double t[DIM_NUM*(DIM_NUM + 1 )] = {
    1.0, 1.0,
    3.0, 1.0,
    2.0, 5.0 };
  double t_unit[DIM_NUM*(DIM_NUM + 1 )] = {
    0.0, 0.0,
    1.0, 0.0,
    0.0, 1.0 };
  int vertex_num = DIM_NUM + 1;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  SIMPLEX_UNIT_TO_GENERAL\n" );
  printf ( "  maps points in the unit simplex to a general simplex.\n" );
  printf ( "\n" );
  printf ( "  Here we consider a simplex in 2D, a triangle.\n" );
  printf ( "\n" );
  printf ( "  The vertices of the general triangle are:\n" );
  printf ( "\n" );
  for ( j = 0; j < vertex_num; j++ )
  {
    for ( dim = 0; dim < dim_num; dim++ )
    {
      printf ( "  %8g", t[dim+j*dim_num] );
    }
    printf ( "\n" );
  }

  printf ( "\n" );
  printf ( "   (  XSI     ETA )   ( X       Y  )\n" );
  printf ( "\n" );

  phy_unit = simplex_unit_to_general ( dim_num, dim_num+1, t, t_unit );

  for ( j = 0; j < dim_num + 1; j++ )
  {
    for ( dim = 0; dim < dim_num; dim++ )
    {
      printf ( "  %9g", t_unit[dim+j*dim_num] );
    }
    for ( dim = 0; dim < dim_num; dim++ )
    {
      printf ( "  %9g", phy_unit[dim+j*dim_num] );
    }
    printf ( "\n" );
  }
  ref = simplex_unit_sample ( dim_num, point_num, &seed );

  phy = simplex_unit_to_general ( dim_num, point_num, t, ref );

  for ( j = 0; j < point_num; j++ )
  {
    for ( dim = 0; dim < dim_num; dim++ )
    {
      printf ( "  %9g", ref[dim+j*dim_num] );
    }
    for ( dim = 0; dim < dim_num; dim++ )
    {
      printf ( "  %9g", phy[dim+j*dim_num] );
    }
    printf ( "\n" );
  }

  free ( phy );
  free ( phy_unit );
  free ( ref );

  return;
# undef DIM_NUM
}
/******************************************************************************/

void test02 ( )

/******************************************************************************/
/*
  Purpose:

    TEST02 tests SIMPLEX_UNIT_TO_GENERAL.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 August 2012

  Author:

    John Burkardt
*/
{
# define DIM_NUM 3

  int dim;
  int dim_num = DIM_NUM;
  int vertex_num = DIM_NUM + 1;
  int j;
  double *phy;
  double *phy_unit;
  int point_num = 10;
  double *ref;
  int seed = 123456789;
  double t[DIM_NUM*(DIM_NUM + 1 )] = {
    1.0, 1.0, 1.0,
    3.0, 1.0, 1.0,
    1.0, 4.0, 1.0,
    1.0, 1.0, 5.0 };
  double t_unit[DIM_NUM*(DIM_NUM + 1 )] = {
    0.0, 0.0, 0.0,
    1.0, 0.0, 0.0,
    0.0, 1.0, 0.0,
    0.0, 0.0, 1.0 };

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  SIMPLEX_UNIT_TO_GENERAL\n" );
  printf ( "  maps points in the unit simplex to a general simplex.\n" );
  printf ( "\n" );
  printf ( "  Here we consider a simplex in 3D, a tetrahedron.\n" );
  printf ( "\n" );
  printf ( "  The vertices of the general tetrahedron are:\n" );
  printf ( "\n" );
  for ( j = 0; j < vertex_num; j++ )
  {
    for ( dim = 0; dim < dim_num; dim++ )
    {
      printf ( "  %8g", t[dim+j*dim_num] );
    }
    printf ( "\n" );
  }

  printf ( "\n" );
  printf ( "   (  XSI     ETA     MU )    ( X       Y       Z )\n" );
  printf ( "\n" );

  phy_unit = simplex_unit_to_general ( dim_num, dim_num+1, t, t_unit );

  for ( j = 0; j < dim_num + 1; j++ )
  {
    for ( dim = 0; dim < dim_num; dim++ )
    {
      printf ( "  %9g", t_unit[dim+j*dim_num] );
    }
    for ( dim = 0; dim < dim_num; dim++ )
    {
      printf ( "  %9g", phy_unit[dim+j*dim_num] );
    }
    printf ( "\n" );
  }

  ref = simplex_unit_sample ( dim_num, point_num, &seed );

  phy = simplex_unit_to_general ( dim_num, point_num, t, ref );

  for ( j = 0; j < point_num; j++ )
  {
    for ( dim = 0; dim < dim_num; dim++ )
    {
      printf ( "  %9g", ref[dim+j*dim_num] );
    }
    for ( dim = 0; dim < dim_num; dim++ )
    {
      printf ( "  %9g", phy[dim+j*dim_num] );
    }
    printf ( "\n" );
  }

  free ( phy );
  free ( phy_unit );
  free ( ref );

  return;
# undef DIM_NUM
}
/******************************************************************************/

void test03 ( )

/******************************************************************************/
/*
  Purpose:

    TEST03 tests GM_RULE_SIZE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 August 2012

  Author:

    John Burkardt
*/
{
# define TEST_NUM 4

  int dim_num;
  int dim_num_test[TEST_NUM] = { 2, 3, 5, 10 };
  int degree;
  int point_num;
  int rule;
  int test;
  int test_num = TEST_NUM;

  printf ( "\n" );
  printf ( "TEST03\n" );
  printf ( "  GM_RULE_SIZE returns POINT_NUM, the number of points\n" );
  printf ( "  associated with a Grundmann-Moeller quadrature rule\n" );
  printf ( "  for the unit simplex of dimension DIM_NUM\n" );
  printf ( "  with rule index RULE\n" );
  printf ( "  and degree of exactness DEGREE = 2*RULE+1.\n" );

  printf ( "\n" );
  printf ( "   DIM_NUM      RULE    DEGREE POINT_NUM\n" );

  for ( test = 0; test < test_num; test++ )
  {
    dim_num = dim_num_test[test];

    printf ( "\n" );

    for ( rule = 0; rule <= 5; rule++ )
    {
      point_num = gm_rule_size ( rule, dim_num );
      degree = 2 * rule + 1;

      printf ( "  %8d  %8d  %8d  %8d\n", dim_num, rule, degree, point_num );
    }
  }

  return;
# undef TEST_NUM
}
/******************************************************************************/

void test04 ( )

/******************************************************************************/
/*
  Purpose:

    TEST04 tests GM_RULE_SET.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 August 2012

  Author:

    John Burkardt
*/
{
  int dim;
  int dim_num;
  int point;
  int point_num;
  int rule;
  double *w;
  double *x;

  printf ( "\n" );
  printf ( "TEST04\n" );
  printf ( "  GM_RULE_SET determines the weights and abscissas\n" );
  printf ( "  of a Grundmann-Moeller quadrature rule for\n" );
  printf ( "  the DIM_NUM dimensional simplex,\n" );
  printf ( "  using a rule of in index RULE,\n" );
  printf ( "  which will have degree of exactness 2*RULE+1.\n" );

  dim_num = 3;
  rule = 2;

  printf ( "\n" );
  printf ( "  Here we use DIM_NUM = %d\n",dim_num );
  printf ( "  RULE = %d\n", rule );
  printf ( "  DEGREE = %d\n", 2 * rule + 1 );

  point_num = gm_rule_size ( rule, dim_num );

  w = ( double * ) malloc ( point_num * sizeof ( double ) );
  x = ( double * ) malloc ( dim_num * point_num * sizeof ( double ) );

  gm_rule_set ( rule, dim_num, point_num, w, x );

  printf ( "\n" );
  printf ( "     POINT        W             X             Y             Z\n" );
  printf ( "\n" );

  for ( point = 0; point < point_num; point++ )
  {
    printf ( "  %8d  %12g", point + 1, w[point] );
    for ( dim = 0; dim < dim_num; dim++ )
    {
      printf ( "  %12g", x[dim+point*dim_num] );
    }
    printf ( "\n" );
  }

  free ( w );
  free ( x );

  return;
}
/******************************************************************************/

void test05 ( )

/******************************************************************************/
/*
  Purpose:

    TEST05 tests GM_RULE_SET.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 August 2012

  Author:

    John Burkardt
*/
{
# define TEST_NUM 4

  int dim_num;
  int dim_num_test[TEST_NUM] = { 2, 3, 5, 10 };
  int point;
  int point_num;
  int rule;
  int test;
  int test_num = TEST_NUM;
  double *w;
  double w_sum;
  double *x;

  printf ( "\n" );
  printf ( "TEST05\n" );
  printf ( "  GM_RULE_SET determines the weights and abscissas\n" );
  printf ( "  of a Grundmann-Moeller quadrature rule for\n" );
  printf ( "  the DIM_NUM dimensional simplex,\n" );
  printf ( "  using a rule of in index RULE,\n" );
  printf ( "  which will have degree of exactness 2*RULE+1.\n" );
  printf ( "\n" );
  printf ( "  In this test, we compute various rules, and simply\n" );
  printf ( "  report the number of points, and the sum of weights.\n" );

  printf ( "\n" );
  printf ( "   DIM_NUM      RULE    POINT_NUM  WEIGHT SUM\n" );

  for ( test = 0; test < test_num; test++ )
  {
    dim_num = dim_num_test[test];

    printf ( "\n" );

    for ( rule = 0; rule <= 5; rule++ )
    {
      point_num = gm_rule_size ( rule, dim_num );

      w = ( double * ) malloc ( point_num * sizeof ( double ) );
      x = ( double * ) malloc ( dim_num * point_num * sizeof ( double ) );

      gm_rule_set ( rule, dim_num, point_num, w, x );

      w_sum = 0.0;
      for ( point = 0; point < point_num; point++ )
      {
        w_sum = w_sum + w[point];
      }

      printf ( "  %8d  %8d  %8d  %24.16g\n", dim_num, rule, point_num, w_sum );

      free ( w );
      free ( x );
    }
  }

  return;
# undef TEST_NUM
}
/******************************************************************************/

void test06 ( )

/******************************************************************************/
/*
  Purpose:

    TEST06 tests GM_RULE_SET.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 August 2012

  Author:

    John Burkardt
*/
{
  int degree;
  int dim;
  int dim_num;
  int point;
  int point_num;
  int rule;
  double *w;
  char w_file[127];
  FILE *w_unit;
  double *x;
  char x_file[127];
  FILE *x_unit;

  printf ( "\n" );
  printf ( "TEST06\n" );
  printf ( "  GM_RULE_SET determines the weights and abscissas\n" );
  printf ( "  of a Grundmann-Moeller quadrature rule for\n" );
  printf ( "  the DIM_NUM dimensional simplex,\n" );
  printf ( "  using a rule of in index RULE,\n" );
  printf ( "  which will have degree of exactness 2*RULE+1.\n" );
  printf ( "\n" );
  printf ( "  In this test, we write a rule to a file.\n" );

  dim_num = 3;
  rule = 2;

  printf ( "\n" );
  printf ( "  Here we use DIM_NUM = %d\n", dim_num  );
  printf ( "  RULE = %d\n", rule );
  printf ( "  DEGREE = %d\n", 2 * rule + 1 );

  point_num = gm_rule_size ( rule, dim_num );

  w = ( double * ) malloc ( point_num * sizeof ( double ) );
  x = ( double * ) malloc ( dim_num * point_num * sizeof ( double ) );

  gm_rule_set ( rule, dim_num, point_num, w, x );

  sprintf ( w_file, "gm%d_%dd_w.txt", rule, dim_num );

  w_unit = fopen ( w_file, "wt" );

  for ( point = 0; point < point_num; point++ )
  {
    fprintf ( w_unit, "%20.16g\n", w[point] );
  }

  fclose ( w_unit );

  sprintf ( x_file, "gm%d_%dd_x.txt", rule, dim_num );

  x_unit = fopen ( x_file, "wt" );

  for ( point = 0; point < point_num; point++ )
  {
    for ( dim = 0; dim < dim_num; dim++ )
    {
      fprintf ( x_unit, "%20.16", x[dim+point*dim_num] );
    }
    fprintf ( x_unit, "\n" );
  }

  fclose ( x_unit );

  printf ( "\n" );
  printf ( "  Wrote rule %d to \"%s\" and \"%s\".\n", rule, w_file, x_file );

  free ( w );
  free ( x );

  return;
}
/******************************************************************************/

void test07 ( )

/******************************************************************************/
/*
  Purpose:

    TEST07 tests GM_RULE_SET.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 August 2012

  Author:

    John Burkardt
*/
{
  int degree;
  int degree_max = 4;
  int dim_num = 5;
  int *expon;
  int h;
  int t;
  double *mono;
  int more;
  int point;
  int point_num;
  double quad;
  double quad_error;
  int rule;
  int rule_max = 3;
  double *w;
  double *x;

  printf ( "\n" );
  printf ( "TEST07\n" );
  printf ( "  GM_RULE_SET determines the weights and abscissas\n" );
  printf ( "  of a Grundmann-Moeller quadrature rule for\n" );
  printf ( "  the DIM_NUM dimensional simplex,\n" );
  printf ( "  using a rule of in index RULE,\n" );
  printf ( "  which will have degree of exactness 2*RULE+1.\n" );
  printf ( "\n" );
  printf ( "  In this test, look at all the monomials up to\n" );
  printf ( "  some maximum degree, choose a few low order rules\n" );
  printf ( "  and determine the quadrature error for each.\n" );
  printf ( "\n" );
  printf ( "  Here we use DIM_NUM = %d\n", dim_num );

  printf ( "\n" );
  printf ( "      Rule     Order     Quad_Error\n" );
  printf ( "\n" );

  expon = ( int * ) malloc ( dim_num * sizeof ( int ) );

  for ( degree = 0; degree <= degree_max; degree++ )
  {
    more = 0;

    for ( ; ; )
    {
      comp_next ( degree, dim_num, expon, &more, &h, &t );

      printf ( "\n" );
      printf ( "  F(X) = X1^%d * X2^%d * X3^%d * X4^%d * X5^%d\n", 
        expon[0], expon[1], expon[2], expon[3], expon[4] );
      printf ( "\n" );

      for ( rule = 0; rule <= rule_max; rule++ )
      {
        point_num = gm_rule_size ( rule, dim_num );

        mono = ( double * ) malloc ( point_num * sizeof ( double ) );
        w = ( double * ) malloc ( point_num * sizeof ( double ) );
        x = ( double * ) malloc ( dim_num * point_num * sizeof ( double ) );

        gm_rule_set ( rule, dim_num, point_num, w, x );

        quad_error = simplex_unit_monomial_quadrature ( dim_num, expon,
          point_num, x, w );

        printf ( "  %8d  %8d  %14g\n", rule, point_num, quad_error );

        free ( mono );
        free ( w );
        free ( x );
      }

      if ( !more )
      {
        break;
      }
    }
  }

  free ( expon );

  return;
}
