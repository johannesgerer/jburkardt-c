# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "tetrahedron_ncc_rule.h"

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );
void test04 ( );
void test05 ( );
void test06 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for TETRAHEDRON_NCC_RULE_PRB.

  Discussion:

    TETRAHEDRON_NCC_RULE_PRB tests the TETRAHEDRON_NCC_RULE library.

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
  printf ( "TETRAHEDRON_NCC_RULE_PRB:\n" );
  printf ( "  C version\n" );
  printf ( "  Test the TETRAHEDRON_NCC_RULE library.\n" );

  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
  test05 ( );
  test06 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "TETRAHEDRON_NCC_RULE_PRB:\n" );
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

    TEST01 tests TETRAHEDRON_NCC_RULE_NUM, TETRAHEDRON_NCC_DEGREE, TETRAHEDRON_NCC_ORDER_NUM.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 January 2007

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
  printf ( "  TETRAHEDRON_NCC_RULE_NUM returns the number of rules;\n" );
  printf ( "  TETRAHEDRON_NCC_DEGREE returns the degree of a rule;\n" );
  printf ( "  TETRAHEDRON_NCC_ORDER_NUM returns the order of a rule.\n" );

  rule_num = tetrahedron_ncc_rule_num ( );

  printf ( "\n" );
  printf ( "  Number of available rules = %d\n", rule_num );
  printf ( "\n" );
  printf ( "      Rule    Degree     Order\n" );
  printf ( "\n" );

  for ( rule = 1; rule <= rule_num; rule++ )
  {
    order_num = tetrahedron_ncc_order_num ( rule );
    degree = tetrahedron_ncc_degree ( rule );
    printf ( "  %8d  %8d  %8d\n", rule, degree, order_num );
  }

  return;
}
/******************************************************************************/

void test02 ( )

/******************************************************************************/
/*
  Purpose:

    TEST02 tests TETRAHEDRON_NCC_RULE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 January 2007

  Author:

    John Burkardt
*/
{
  int dim_num = 3;
  int order;
  int order_num;
  int rule;
  int rule_num;
  double *wtab;
  double wtab_sum;
  double *xyztab;

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  TETRAHEDRON_NCC_RULE returns the points and weights\n" );
  printf ( "  of an NCC rule for the tetrahedron.\n" );
  printf ( "\n" );
  printf ( "  In this test, we simply check that the weights\n" );
  printf ( "  sum to 1.\n" );

  rule_num = tetrahedron_ncc_rule_num ( );

  printf ( "\n" );
  printf ( "  Number of available rules = %d\n", rule_num );

  printf ( "\n" );
  printf ( "      Rule      Sum of weights\n" );
  printf ( "\n" );

  for ( rule = 1; rule <= rule_num; rule++ )
  {
    order_num = tetrahedron_ncc_order_num ( rule );

    xyztab = ( double * ) malloc ( dim_num * order_num * sizeof ( double ) );
    wtab = ( double * ) malloc ( order_num * sizeof ( double ) );

    tetrahedron_ncc_rule ( rule, order_num, xyztab, wtab );

    wtab_sum = 0.0;
    for ( order = 0; order < order_num; order++ )
    {
      wtab_sum = wtab_sum + wtab[order];
    }

    printf ( "  %8d  %14g\n", rule, wtab_sum );

    free ( wtab );
    free ( xyztab );
  }
  return;
}
/******************************************************************************/

void test03 ( )

/******************************************************************************/
/*
  Purpose:

    TEST03 tests TETRAHEDRON_NCC_RULE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 January 2007

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
  printf ( "  TETRAHEDRON_NCC_RULE returns the points and weights\n" );
  printf ( "  of an NCC rule for the tetrahedron.\n" );
  printf ( "\n" );
  printf ( "  In this test, we simply check that, for each\n" );
  printf ( "  quadrature point, the barycentric coordinates\n" );
  printf ( "  add up to 1.\n" );

  rule_num = tetrahedron_ncc_rule_num ( );

  printf ( "\n" );
  printf ( "      Rule    Suborder    Sum of coordinates\n" );
  printf ( "\n" );

  for ( rule = 1; rule <= rule_num; rule++ )
  {
    suborder_num = tetrahedron_ncc_suborder_num ( rule );

    suborder_xyz = ( double * ) malloc ( 4 * suborder_num * sizeof ( double ) );
    suborder_w = ( double * ) malloc ( suborder_num * sizeof ( double ) );

    tetrahedron_ncc_subrule ( rule, suborder_num, suborder_xyz, suborder_w );

    printf ( "\n" );
    printf ( "  %8d  %8d\n", rule, suborder_num );

    for ( suborder = 0; suborder < suborder_num; suborder++ )
    {
      xyz_sum = suborder_xyz[0+suborder*4]
              + suborder_xyz[1+suborder*4]
              + suborder_xyz[2+suborder*4]
              + suborder_xyz[3+suborder*4];
     printf ( "                      %25.16g\n", xyz_sum  );
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

    TEST04 tests TETRAHEDRON_NCC_RULE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 January 2007

  Author:

    John Burkardt
*/
{
  int a;
  int b;
  int c;
  double coef;
  int dim_num = 3;
  double err;
  double exact;
  int i;
  int order;
  int order_num;
  double quad;
  int rule;
  int rule_num;
  double value;
  double volume;
  double *wtab;
  double x;
  double *xyztab;
  double y;
  double z;

  printf ( "\n" );
  printf ( "TEST04\n" );
  printf ( "  TETRAHEDRON_NCC_RULE returns the points and weights of\n" );
  printf ( "  an NCC rule for the unit tetrahedron.\n" );
  printf ( "\n" );
  printf ( "  This routine uses those rules to estimate the\n" );
  printf ( "  integral of monomomials in the unit tetrahedron.\n" );

  rule_num = tetrahedron_ncc_rule_num ( );

  volume = 1.0 / 6.0;

  for ( a = 0; a <= 6; a++ )
  {
    for ( b = 0; b <= 6 - a; b++ )
    {
      for ( c = 0; c <= 6 - a - b; c++ )
      {
/*
  Multiplying X^A * Y^B * Z^C by COEF will give us an integrand
  whose integral is exactly 1.  This makes the error calculations easy.
*/
        coef = 1.0;

        for ( i = a + 1; i <= a + b; i++ )
        {
          coef = coef * ( double ) ( i ) / ( double ) ( i - a );
        }
        for ( i = a + b + 1; i <= a + b + c; i++ )
        {
          coef = coef * ( double ) ( i ) / ( double ) ( i - a - b );
        }
        for ( i = a + b + c + 1; i <= a + b + c + 3; i++ )
        {
          coef = coef * ( double ) ( i );
        }

        printf ( "\n" );
        printf ( "  Integrate %g * X^%d * Y^%d * Z^%d\n", coef, a, b, c );
        printf ( "\n" );
        printf ( "      Rule       QUAD           ERROR\n" );
        printf ( "\n" );

        for ( rule = 1; rule <= rule_num; rule++ )
        {
          order_num = tetrahedron_ncc_order_num ( rule );

          xyztab = ( double * ) malloc ( dim_num * order_num * sizeof ( double ) );
          wtab = ( double * ) malloc ( order_num * sizeof ( double ) );

          tetrahedron_ncc_rule ( rule, order_num, xyztab, wtab );

          quad = 0.0;

          for ( order = 0; order < order_num; order++ )
          {
            x = xyztab[0+order*dim_num];
            y = xyztab[1+order*dim_num];
            z = xyztab[2+order*dim_num];
/*
  Some tedious calculations to avoid 0^0 complaints.
*/
            value = coef;

            if ( a != 0 )
            {
              value = value * pow ( x, a );
            }

            if ( b != 0 )
            {
              value = value * pow ( y, b );
            }

            if ( c != 0 )
            {
              value = value * pow ( z, c );
            }

            quad = quad + wtab[order] * value;
          }
          quad = volume * quad;

          exact = 1.0;
          err = fabs ( exact - quad );

          printf ( "  %8d  %14g  %14g\n", rule, quad, err  );

          free ( wtab );
          free ( xyztab );
        }
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

    TEST05 demonstrates REFERENCE_TO_PHYSICAL_T4.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 January 2007

  Author:

    John Burkardt
*/
{
# define DIM_NUM 3
# define NODE_NUM 4

  int i;
  int node;
  double node_xyz[DIM_NUM*NODE_NUM] = {
    0.0, 0.0, 0.0,
    1.0, 0.0, 0.0,
    0.0, 1.0, 0.0,
    0.0, 0.0, 1.0 };
  double node_xyz2[DIM_NUM*NODE_NUM] = {
    4.0, 5.0, 1.0,
    6.0, 5.0, 1.0,
    4.0, 8.0, 1.0,
    4.0, 5.0, 5.0 };
  int order;
  int order_num;
  int rule;
  double volume;
  double volume2;
  double *w;
  double *xyz;
  double *xyz2;

  printf ( "\n" );
  printf ( "TEST06\n" );
  printf ( "  REFERENCE_TO_PHYSICAL_T4 transforms a rule\n" );
  printf ( "  on the unit (reference) tetrahedron to a rule on \n" );
  printf ( "  an arbitrary (physical) tetrahedron.\n" );

  rule = 3;

  order_num = tetrahedron_ncc_order_num ( rule );

  xyz = ( double * ) malloc ( 3 * order_num * sizeof ( double ) );
  xyz2 = ( double * ) malloc ( 3 * order_num * sizeof ( double ) );
  w = ( double * ) malloc ( order_num * sizeof ( double ) );

  tetrahedron_ncc_rule ( rule, order_num, xyz, w );
/*
  Here is the reference tetrahedron, and its rule.
*/
  printf ( "\n" );
  printf ( "  The reference tetrahedron:\n" );
  printf ( "\n" );

  for ( node = 0; node < NODE_NUM; node++ )
  {
    printf ( "  %8d  %14g  %14g  %14g\n",
      node+1, node_xyz[0+node*DIM_NUM], node_xyz[1+node*DIM_NUM], node_xyz[2+node*DIM_NUM] );
  }

  volume = tetrahedron_volume ( node_xyz );

  printf ( "\n" );
  printf ( "  Rule %d for reference tetrahedron\n", rule );
  printf ( "  with volume = %g\n", volume );
  printf ( "\n" );
  printf ( "                X               Y               Z               W\n" );
  printf ( "\n" );

  for ( order = 0; order < order_num; order++ )
  {
    printf ( "  %8d  %14g  %14g  %14g  %14g\n",
      order, xyz[0+order*DIM_NUM], xyz[1+order*DIM_NUM], xyz[2+order*DIM_NUM], w[order] );
  }
/*
  Transform the rule.
*/
  reference_to_physical_t4 ( node_xyz2, order_num, xyz, xyz2 );
/*
  Here is the physical tetrahedron, and its transformed rule.
*/
  printf ( "\n" );
  printf ( "  The physical tetrahedron:\n" );
  printf ( "\n" );

  for ( node = 0; node < NODE_NUM; node++ )
  {
    printf ( "  %8d  %14g  %14g  %14g\n",
      node+1, node_xyz2[0+node*DIM_NUM], node_xyz2[1+node*DIM_NUM], 
      node_xyz2[2+node*DIM_NUM] );
  }

  volume2 = tetrahedron_volume ( node_xyz2 );

  printf ( "\n" );
  printf ( "  Rule %d for physical tetrahedron\n", rule );
  printf ( "  with volume = %g\n", volume2 );
  printf ( "\n" );
  printf ( "                X               Y               Z               W\n" );
  printf ( "\n" );

  for ( order = 0; order < order_num; order++ )
  {
    printf ( "  %8d  %14g  %14g  %14g  %14g\n",
      order, xyz2[0+order*DIM_NUM], xyz2[1+order*DIM_NUM], 
      xyz2[2+order*DIM_NUM], w[order] );
  }

  free ( w );
  free ( xyz );
  free ( xyz2 );

  return;
# undef DIM_NUM
# undef NODE_NUM
}
/******************************************************************************/

void test06 ( )

/******************************************************************************/
/*
  Purpose:

    TEST06 tests TETRAHEDRON_NCC_RULE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 January 2007

  Author:

    John Burkardt
*/
{
  int dim_num = 3;
  int o;
  int order_num;
  int rule;
  int s;
  int *suborder;
  int suborder_num;
  double *suborder_w;
  double *suborder_xyz;
  double *w;
  double *xyz;

  printf ( "\n" );
  printf ( "TEST06\n" );
  printf ( "  TETRAHEDRON_NCC_RULE returns the points and weights\n" );
  printf ( "  of an NCC rule for the tetrahedron.\n" );

  rule = 4;

  printf ( "\n" );
  printf ( "  In this test, we simply print rule %d\n", rule );

  suborder_num = tetrahedron_ncc_suborder_num ( rule );

  suborder = tetrahedron_ncc_suborder ( rule, suborder_num );

  suborder_w = ( double * ) malloc ( suborder_num * sizeof ( double ) );
  suborder_xyz = ( double * ) malloc ( 4 * suborder_num * sizeof ( double ) );

  tetrahedron_ncc_subrule ( rule, suborder_num, suborder_xyz, suborder_w );

  printf ( "\n" );
  printf ( "  The compressed rule:\n" );
  printf ( "\n" );
  printf ( "  Number of suborders = %d\n", suborder_num );
  printf ( "\n" );
  printf ( "     S   Sub     Weight     Xsi1      Xsi2      Xsi3      Xsi4\n" );
  printf ( "\n" );

  for ( s = 0; s < suborder_num; s++ )
  {
    printf ( "  %4d  %8d  %8d  %8g  %8g  %8g\n",
      s+1, suborder[s], suborder_w[s], suborder_xyz[0+s*4],
      suborder_xyz[1+s*4], suborder_xyz[2+s*4], suborder_xyz[3+s*4] );
  }

  order_num = tetrahedron_ncc_order_num ( rule );

  xyz = ( double * ) malloc ( dim_num * order_num * sizeof ( double ) );
  w = ( double * ) malloc ( order_num * sizeof ( double ) );

  tetrahedron_ncc_rule ( rule, order_num, xyz, w );

  printf ( "\n" );
  printf ( "  The full rule:\n" );
  printf ( "\n" );
  printf ( "  Order = %d\n", order_num );
  printf ( "\n" );
  printf ( "     O    Weight        X           Y           Z\n" );
  printf ( "\n" );

  for ( o = o; o < order_num; o++ )
  {
    printf ( "  %4d  %8g  %8g  %8g  %8g\n",
      o+1, w[o], xyz[0+o*3], xyz[1+o*3], xyz[2+o*3] );
  }

  free ( suborder );
  free ( suborder_w );
  free ( suborder_xyz );
  free ( w );
  free ( xyz );

  return;
}

