# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "triangle_felippa_rule.h"

int main ( );
void triangle_unit_monomial_test ( int degree_max );
void triangle_unit_quad_test ( int degree_max );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for TRIANGLE_FELIPPA_RULE_PRB.

  Discussion:

    TRIANGLE_FELIPPA_RULE_PRB tests the TRIANGLE_FELIPPA_RULE library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 August 2014

  Author:

    John Burkardt
*/
{
  int degree_max;

  timestamp ( );
  printf ( "\n" );
  printf ( "TRIANGLE_FELIPPA_RULE_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the TRIANGLE_FELIPPA_RULE library.\n" );

  degree_max = 4;
  triangle_unit_monomial_test ( degree_max );

  degree_max = 7;
  triangle_unit_quad_test ( degree_max );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "TRIANGLE_FELIPPA_RULE_PRB\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void triangle_unit_monomial_test ( int degree_max )

/******************************************************************************/
/*
  Purpose:

    TRIANGLE_UNIT_MONOMIAL_TEST tests TRIANGLE_UNIT_MONOMIAL.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 April 2009

  Author:

    John Burkardt

  Parameters:

    Input, int DEGREE_MAX, the maximum total degree of the
    monomials to check.
*/
{
  int alpha;
  int beta;
  int expon[2];
  double value;

  printf ( "\n" );
  printf ( "TRIANGLE_UNIT_MONOMIAL_TEST\n" );
  printf ( "  For the unit triangle,\n" );
  printf ( "  TRIANGLE_UNIT_MONOMIAL returns the exact value of the\n" );
  printf ( "  integral of X^ALPHA Y^BETA\n" );
  printf ( "\n" );
  printf ( "  Volume = %g\n", triangle_unit_volume ( ) );
  printf ( "\n" );
  printf ( "     ALPHA      BETA      INTEGRAL\n" );
  printf ( "\n" );

  for ( alpha = 0; alpha <= degree_max; alpha++ )
  {
    expon[0] = alpha;
    for ( beta = 0; beta <= degree_max - alpha; beta++ )
    {
      expon[1] = beta;

      value = triangle_unit_monomial ( expon );

      printf ( "  %8d  %8d  %14.6g\n", expon[0], expon[1], value );
    }
  }

  return;
}
/******************************************************************************/

void triangle_unit_quad_test ( int degree_max )

/******************************************************************************/
/*
  Purpose:

    TRIANGLE_UNIT_QUAD_TEST tests the rules for the unit triangle.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    18 April 2008

  Author:

    John Burkardt

  Parameters:

    Input, int DEGREE_MAX, the maximum total degree of the
    monomials to check.
*/
{
# define DIM_NUM 2

  int dim;
  int dim_num = DIM_NUM;
  int expon[DIM_NUM];
  int h;
  int more;
  int order;
  double quad;
  int t;
  double *v;
  double *w;
  double *xy;

  printf ( "\n" );
  printf ( "TRIANGLE_UNIT_QUAD_TEST\n" );
  printf ( "  For the unit triangle,\n" );
  printf ( "  we approximate monomial integrals with:\n" );
  printf ( "  TRIANGLE_UNIT_O01,\n" );
  printf ( "  TRIANGLE_UNIT_O03,\n" );
  printf ( "  TRIANGLE_UNIT_O03b,\n" );
  printf ( "  TRIANGLE_UNIT_O06,\n" );
  printf ( "  TRIANGLE_UNIT_O06b,\n" );
  printf ( "  TRIANGLE_UNIT_O07,\n" );
  printf ( "  TRIANGLE_UNIT_O12,\n" );

  more = 0;

  for ( ; ; )
  {
    subcomp_next ( degree_max, dim_num, expon, &more, &h, &t );

    printf ( "\n" );
    printf ( "  Monomial exponents: " );
    for ( dim = 0; dim < dim_num; dim++ )
    {
      printf ( "  %2d", expon[dim] );
    }
    printf ( "\n" );
    printf ( "\n" );

    order = 1;
    w = ( double * ) malloc ( order * sizeof ( double ) );
    xy = ( double * ) malloc ( dim_num * order * sizeof ( double ) );
    triangle_unit_o01 ( w, xy );
    v = monomial_value ( dim_num, order, expon, xy );
    quad = triangle_unit_volume ( ) * r8vec_dot_product ( order, w, v );
    printf ( "  %6d  %14.6g\n", order, quad );
    free ( v );
    free ( w );
    free ( xy );

    order = 3;
    w = ( double * ) malloc ( order * sizeof ( double ) );
    xy = ( double * ) malloc ( dim_num * order * sizeof ( double ) );
    triangle_unit_o03 ( w, xy );
    v = monomial_value ( dim_num, order, expon, xy );
    quad = triangle_unit_volume ( ) * r8vec_dot_product ( order, w, v );
    printf ( "  %6d  %14.6g\n", order, quad );
    free ( v );
    free ( w );
    free ( xy );

    order = 3;
    w = ( double * ) malloc ( order * sizeof ( double ) );
    xy = ( double * ) malloc ( dim_num * order * sizeof ( double ) );
    triangle_unit_o03b ( w, xy );
    v = monomial_value ( dim_num, order, expon, xy );
    quad = triangle_unit_volume ( ) * r8vec_dot_product ( order, w, v );
    printf ( "  %6d  %14.6g\n", order, quad );
    free ( v );
    free ( w );
    free ( xy );

    order = 6;
    w = ( double * ) malloc ( order * sizeof ( double ) );
    xy = ( double * ) malloc ( dim_num * order * sizeof ( double ) );
    triangle_unit_o06 ( w, xy );
    v = monomial_value ( dim_num, order, expon, xy );
    quad = triangle_unit_volume ( ) * r8vec_dot_product ( order, w, v );
    printf ( "  %6d  %14.6g\n", order, quad );
    free ( v );
    free ( w );
    free ( xy );

    order = 6;
    w = ( double * ) malloc ( order * sizeof ( double ) );
    xy = ( double * ) malloc ( dim_num * order * sizeof ( double ) );
    triangle_unit_o06b ( w, xy );
    v = monomial_value ( dim_num, order, expon, xy );
    quad = triangle_unit_volume ( ) * r8vec_dot_product ( order, w, v );
    printf ( "  %6d  %14.6g\n", order, quad );
    free ( v );
    free ( w );
    free ( xy );

    order = 7;
    w = ( double * ) malloc ( order * sizeof ( double ) );
    xy = ( double * ) malloc ( dim_num * order * sizeof ( double ) );
    triangle_unit_o07 ( w, xy );
    v = monomial_value ( dim_num, order, expon, xy );
    quad = triangle_unit_volume ( ) * r8vec_dot_product ( order, w, v );
    printf ( "  %6d  %14.6g\n", order, quad );
    free ( v );
    free ( w );
    free ( xy );

    order = 12;
    w = ( double * ) malloc ( order * sizeof ( double ) );
    xy = ( double * ) malloc ( dim_num * order * sizeof ( double ) );
    triangle_unit_o12 ( w, xy );
    v = monomial_value ( dim_num, order, expon, xy );
    quad = triangle_unit_volume ( ) * r8vec_dot_product ( order, w, v );
    printf ( "  %6d  %14.6g\n", order, quad );
    free ( v );
    free ( w );
    free ( xy );

    printf ( "\n" );
    quad = triangle_unit_monomial ( expon );
    printf ( "   Exact  %14.6g\n", quad );

    if ( !more )
    {
      break;
    }
  }

  return;
# undef DIM_NUM
}