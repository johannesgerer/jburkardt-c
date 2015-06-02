# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "triangle_ncc_rule.h"

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );
void test04 ( );
void test05 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for TRIANGLE_NCC_RULE_PRB.

  Discussion:

    TRIANGLE_NCC_RULE_PRB tests the TRIANGLE_NCC_RULE library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 June 2014

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "TRIANGLE_NCC_RULE_PRB:\n" );
  printf ( "  C version\n" );
  printf ( "  Test the TRIANGLE_NCC_RULE library.\n" );

  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
  test05 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "TRIANGLE_NCC_RULE_PRB:\n" );
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

    TEST01 tests TRIANGLE_NCC_RULE_NUM, TRIANGLE_NCC_DEGREE, TRIANGLE_NCC_ORDER_NUM.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 January 2007

  Author:

    John Burkardt
*/
{
  int degree;
  int order_num;
  int rule;
  int rule_num;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  TRIANGLE_NCC_RULE_NUM returns the number of rules;\n" );
  printf ( "  TRIANGLE_NCC_DEGREE returns the degree of a rule;\n" );
  printf ( "  TRIANGLE_NCC_ORDER_NUM returns the order of a rule.\n" );

  rule_num = triangle_ncc_rule_num ( );

  printf ( "\n" );
  printf ( "  Number of available rules = %d\n", rule_num );
  printf ( "\n" );
  printf ( "      Rule    Degree     Order\n" );
  printf ( "\n" );

  for ( rule = 1; rule <= rule_num; rule++ )
  {
    order_num = triangle_ncc_order_num ( rule );
    degree = triangle_ncc_degree ( rule );
    printf ( "  %8d  %8d  %8d\n", rule, degree, order_num );
  }

  return;
}
/******************************************************************************/

void test02 ( )

/******************************************************************************/
/*
  Purpose:

    TEST02 tests TRIANGLE_NCC_RULE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 January 2007

  Author:

    John Burkardt
*/
{
  int order;
  int order_num;
  int rule;
  int rule_num;
  double *wtab;
  double wtab_sum;
  double *xytab;

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  TRIANGLE_NCC_RULE returns the points and weights\n" );
  printf ( "  of an NCC rule for the triangle.\n" );
  printf ( "\n" );
  printf ( "  In this test, we simply check that the weights\n" );
  printf ( "  sum to 1.\n" );

  rule_num = triangle_ncc_rule_num ( );

  printf ( "\n" );
  printf ( "  Number of available rules = %d\n", rule_num );

  printf ( "\n" );
  printf ( "      Rule      Order     Sum of weights\n" );
  printf ( "\n" );

  for ( rule = 1; rule <= rule_num; rule++ )
  {
    order_num = triangle_ncc_order_num ( rule );

    xytab = ( double * ) malloc ( 2 * order_num * sizeof ( double ) );
    wtab = ( double * ) malloc ( order_num * sizeof ( double ) );

    triangle_ncc_rule ( rule, order_num, xytab, wtab );

    wtab_sum = 0.0;
    for ( order = 0; order < order_num; order++ )
    {
      wtab_sum = wtab_sum + wtab[order];
    }

    printf ( "  %8d  %8d  %14g\n", rule, order_num, wtab_sum );

    free ( wtab );
    free ( xytab );
  }
  return;
}
/******************************************************************************/

void test03 ( )

/******************************************************************************/
/*
  Purpose:

    TEST03 tests TRIANGLE_NCC_RULE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 June 2014

  Author:

    John Burkardt
*/
{
  int rule;
  int rule_num;
  int suborder;
  int suborder_num;
  double *suborder_w;
  double *suborder_xyz;
  double xyz_sum;

  printf ( "\n" );
  printf ( "TEST03\n" );
  printf ( "  TRIANGLE_NCC_RULE returns the points and weights\n" );
  printf ( "  of an NCC rule for the triangle.\n" );
  printf ( "\n" );
  printf ( "  In this test, we simply check that, for each\n" );
  printf ( "  quadrature point, the barycentric coordinates\n" );
  printf ( "  add up to 1.\n" );

  rule_num = triangle_ncc_rule_num ( );

  printf ( "\n" );
  printf ( "      Rule    Suborder    Sum of coordinates\n" );
  printf ( "\n" );

  for ( rule = 1; rule <= rule_num; rule++ )
  {
    suborder_num = triangle_ncc_suborder_num ( rule );

    suborder_xyz = ( double * ) malloc ( 3 * suborder_num * sizeof ( double ) );
    suborder_w = ( double * ) malloc ( suborder_num * sizeof ( double ) );

    triangle_ncc_subrule ( rule, suborder_num, suborder_xyz, suborder_w );

    printf ( "\n" );
    printf ( "  %8d  %8d\n", rule, suborder_num );

    for ( suborder = 0; suborder < suborder_num; suborder++ )
    {
      xyz_sum = suborder_xyz[0+suborder*3]
              + suborder_xyz[1+suborder*3]
              + suborder_xyz[2+suborder*3];
     printf ( "                      %25.16g\n", xyz_sum );
    }

    free ( suborder_w );
    free ( suborder_xyz );
  }
  return;
}
/******************************************************************************/

void test04 ( )

/******************************************************************************/
/*
  Purpose:

    TEST04 tests TRIANGLE_NCC_RULE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 January 2007

  Author:

    John Burkardt
*/
{
  int a;
  double area = 0.5;
  int b;
  double coef;
  double err;
  double exact;
  int i;
  int order;
  int order_num;
  double quad;
  int rule;
  int rule_num;
  double value;
  double *wtab;
  double x;
  double *xytab;
  double y;

  printf ( "\n" );
  printf ( "TEST04\n" );
  printf ( "  TRIANGLE_NCC_RULE returns the points and weights of\n" );
  printf ( "  an NCC rule for the unit triangle.\n" );
  printf ( "\n" );
  printf ( "  This routine uses those rules to estimate the\n" );
  printf ( "  integral of monomomials in the unit triangle.\n" );

  rule_num = triangle_ncc_rule_num ( );

  for ( a = 0; a <= 10; a++ )
  {
    for ( b = 0; b <= 10 - a; b++ )
    {
/*
  Multiplying X^A * Y^B by COEF will give us an integrand
  whose integral is exactly 1.  This makes the error calculations easy.
*/
      coef = ( double ) ( a + b + 2 ) * ( double ) ( a + b + 1 );
      for ( i = 1; i <= b; i++ )
      {
        coef = coef * ( double ) ( a + i ) / ( double ) ( i );
      }

      printf ( "\n" );
      printf ( "  Integrate %g * X^%d * Y^%d\n", coef, a, b );
      printf ( "\n" );
      printf ( "      Rule       QUAD           ERROR\n" );
      printf ( "\n" );

      for ( rule = 1; rule <= rule_num; rule++ )
      {
        order_num = triangle_ncc_order_num ( rule );

        xytab = ( double * ) malloc ( 2 * order_num * sizeof ( double ) );
        wtab = ( double * ) malloc ( order_num * sizeof ( double ) );

        triangle_ncc_rule ( rule, order_num, xytab, wtab );

        quad = 0.0;

        for ( order = 0; order < order_num; order++ )
        {
          x = xytab[0+order*2];
          y = xytab[1+order*2];

          if ( a == 0 && b == 0 )
          {
            value = coef;
          }
          else if ( a == 0 && b != 0 )
          {
            value = coef * pow ( y, b );
          }
          else if ( a != 0 && b == 0 )
          {
            value = coef * pow ( x, a );
          }
          else if ( a != 0 && b != 0 )
          {
            value = coef * pow ( x, a ) * pow ( y, b );
          }
          quad = quad + wtab[order] * value;
        }
        quad = area * quad;

        exact = 1.0;
        err = fabs ( exact - quad );

        printf ( "  %8d  %14g  %14g\n", rule, quad, err );

        free ( wtab );
        free ( xytab ); 
      }
    }
  }
  return;
}
/******************************************************************************/

void test05 ( )

/******************************************************************************/
/*
  Purpose:

    TEST05 demonstrates REFERENCE_TO_PHYSICAL_T3.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 January 2007

  Author:

    John Burkardt
*/
{
# define DIM_NUM 2
# define NODE_NUM 3

  double area;
  double area2;
  int i;
  int node;
  double node_xy[2*NODE_NUM] = {
    0.0, 0.0,
    1.0, 0.0,
    0.0, 1.0 };
  double node_xy2[2*NODE_NUM] = {
    1.0, 2.0,
    1.0, 1.0,
    3.0, 2.0 };
  int order;
  int order_num;
  int rule;
  double *w;
  double *xy;
  double *xy2;

  printf ( "\n" );
  printf ( "TEST05\n" );
  printf ( "  REFERENCE_TO_PHYSICAL_T3 transforms a rule\n" );
  printf ( "  on the unit (reference) triangle to a rule on \n" );
  printf ( "  an arbitrary (physical) triangle.\n" );

  rule = 3;

  order_num = triangle_ncc_order_num ( rule );

  xy = ( double * ) malloc ( 2 * order_num * sizeof ( double ) );
  xy2 = ( double * ) malloc ( 2 * order_num * sizeof ( double ) );
  w = ( double * ) malloc ( order_num * sizeof ( double ) );

  triangle_ncc_rule ( rule, order_num, xy, w );
/*
  Here is the reference triangle, and its rule.
*/
  printf ( "\n" );
  printf ( "  The reference triangle:\n" );
  printf ( "\n" );

  for ( node = 0; node < NODE_NUM; node++ )
  {
    printf ( "  %8d  %14g  %14g\n", node+1, node_xy[0+node*2], node_xy[1+node*2] );
  }

  area = triangle_area ( node_xy );

  printf ( "\n" );
  printf ( "  Rule %d for reference triangle\n", rule );
  printf ( "  with area = %g\n", area);
  printf ( "\n" );
  printf ( "                X               Y               W\n" );
  printf ( "\n" );

  for ( order = 0; order < order_num; order++ )
  {
    printf ( "  %8d  %14g  %14g  %14g\n", order, xy[0+order*2], xy[1+order*2], w[order] );
  }
/*
  Transform the rule.
*/
  reference_to_physical_t3 ( node_xy2, order_num, xy, xy2 );
/*
  Here is the physical triangle, and its transformed rule.
*/
  printf ( "\n" );
  printf ( "  The physical triangle:\n" );
  printf ( "\n" );

  for ( node = 0; node < NODE_NUM; node++ )
  {
    printf ( "  %8d  %14g  %14g\n", node+1, node_xy2[0+node*2], node_xy2[1+node*2] );
  }

  area2 = triangle_area ( node_xy2 );

  printf ( "\n" );
  printf ( "  Rule %d for physical triangle\n", rule );
  printf ( "  with area = %g\n", area2 );
  printf ( "\n" );
  printf ( "                X               Y               W\n" );
  printf ( "\n" );

  for ( order = 0; order < order_num; order++ )
  {
    printf ( "  %8d  %14g  %14g  %14g\n",
      order, xy2[0+order*2], xy2[1+order*2], w[order] );
  }

  free ( w );
  free ( xy );
  free ( xy2 );

  return;
# undef DIM_NUM
# undef NODE_NUM
}
