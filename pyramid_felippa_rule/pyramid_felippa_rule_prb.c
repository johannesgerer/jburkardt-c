# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "pyramid_felippa_rule.h"

int main ( );
void pyramid_unit_monomial_test ( int degree_max );
void pyramid_unit_quad_test ( int degree_max );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for PYRAMID_FELIPPA_RULE_PRB.

  Discussion:

    PYRAMID_FELIPPA_RULE_PRB tests the PYRAMID_FELIPPA_RULE library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    28 August 2014

  Author:

    John Burkardt
*/
{
  int degree_max;

  timestamp ( );
  printf ( "\n" );
  printf ( "PYRAMID_FELIPPA_RULE_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the PYRAMID_FELIPPA_RULE library.\n" );

  degree_max = 4;
  pyramid_unit_monomial_test ( degree_max );

  degree_max = 5;
  pyramid_unit_quad_test ( degree_max );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "PYRAMID_FELIPPA_RULE_PRB\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void pyramid_unit_monomial_test ( int degree_max )

/******************************************************************************/
/*
  Purpose:

    PYRAMID_UNIT_MONOMIAL_TEST tests PYRAMID_UNIT_MONOMIAL.

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
  int expon[3];
  int gamma;
  double value;

  printf ( "\n" );
  printf ( "PYRAMID_UNIT_MONOMIAL_TEST\n" );
  printf ( "  For the unit pyramid,\n" );
  printf ( "  PYRAMID_UNIT_MONOMIAL returns the exact value of the\n" );
  printf ( "  integral of X^ALPHA Y^BETA Z^GAMMA\n" );
  printf ( "\n" );
  printf ( "  Volume = %g\n", pyramid_unit_volume ( ) );
  printf ( "\n" );
  printf ( "     ALPHA      BETA     GAMMA      INTEGRAL\n" );
  printf ( "\n" );

  for ( alpha = 0; alpha <= degree_max; alpha++ )
  {
    expon[0] = alpha;
    for ( beta = 0; beta <= degree_max - alpha; beta++ )
    {
      expon[1] = beta;
      for ( gamma = 0; gamma <= degree_max - alpha - beta; gamma++ )
      {
        expon[2] = gamma;

        value = pyramid_unit_monomial ( expon );

        printf ( "  %8d  %8d  %8d  %14.6g\n",
          expon[0], expon[1], expon[2], value );
      }
    }
  }

  return;
}
/******************************************************************************/

void pyramid_unit_quad_test ( int degree_max )

/******************************************************************************/
/*
  Purpose:

    PYRAMID_UNIT_QUAD_TEST tests the rules for the unit pyramid.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 April 2008

  Author:

    John Burkardt

  Parameters:

    Input, int DEGREE_MAX, the maximum total degree of the
    monomials to check.
*/
{
# define DIM_NUM 3

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
  double *xyz;

  printf ( "\n" );
  printf ( "PYRAMID_UNIT_QUAD_TEST\n" );
  printf ( "  For the unit pyramid,\n" );
  printf ( "  we approximate monomial integrals with:\n" );
  printf ( "  PYRAMID_UNIT_O01,\n" );
  printf ( "  PYRAMID_UNIT_O05,\n" );
  printf ( "  PYRAMID_UNIT_O06,\n" );
  printf ( "  PYRAMID_UNIT_O08,\n" );
  printf ( "  PYRAMID_UNIT_O08b,\n" );
  printf ( "  PYRAMID_UNIT_O09,\n" );
  printf ( "  PYRAMID_UNIT_O13,\n" );
  printf ( "  PYRAMID_UNIT_O18,\n" );
  printf ( "  PYRAMID_UNIT_O27,\n" );
  printf ( "  PYRAMID_UNIT_O48,\n" );

  more = 0;

  for ( ; ; )
  {
    subcomp_next ( degree_max, dim_num, expon, &more, &h, &t );

    if ( ( expon[0] % 2 ) == 1 || ( expon[1] % 2 ) == 1 )
    {
      continue;
    }

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
    xyz = ( double * ) malloc ( dim_num * order * sizeof ( double ) );
    pyramid_unit_o01 ( w, xyz );
    v = monomial_value ( dim_num, order, expon, xyz );
    quad = pyramid_unit_volume ( ) * r8vec_dot_product ( order, w, v );
    printf ( "  %6d  %14.6g\n", order, quad );
    free ( v );
    free ( w );
    free ( xyz );

    order = 5;
    w = ( double * ) malloc ( order * sizeof ( double ) );
    xyz = ( double * ) malloc ( dim_num * order * sizeof ( double ) );
    pyramid_unit_o05 ( w, xyz );
    v = monomial_value ( dim_num, order, expon, xyz );
    quad = pyramid_unit_volume ( ) * r8vec_dot_product ( order, w, v );
    printf ( "  %6d  %14.6g\n", order, quad );
    free ( v );
    free ( w );
    free ( xyz );

    order = 6;
    w = ( double * ) malloc ( order * sizeof ( double ) );
    xyz = ( double * ) malloc ( dim_num * order * sizeof ( double ) );
    pyramid_unit_o06 ( w, xyz );
    v = monomial_value ( dim_num, order, expon, xyz );
    quad = pyramid_unit_volume ( ) * r8vec_dot_product ( order, w, v );
    printf ( "  %6d  %14.6g\n", order, quad );
    free ( v );
    free ( w );
    free ( xyz );

    order = 8;
    w = ( double * ) malloc ( order * sizeof ( double ) );
    xyz = ( double * ) malloc ( dim_num * order * sizeof ( double ) );
    pyramid_unit_o08 ( w, xyz );
    v = monomial_value ( dim_num, order, expon, xyz );
    quad = pyramid_unit_volume ( ) * r8vec_dot_product ( order, w, v );
    printf ( "  %6d  %14.6g\n", order, quad );
    free ( v );
    free ( w );
    free ( xyz );

    order = 8;
    w = ( double * ) malloc ( order * sizeof ( double ) );
    xyz = ( double * ) malloc ( dim_num * order * sizeof ( double ) );
    pyramid_unit_o08b ( w, xyz );
    v = monomial_value ( dim_num, order, expon, xyz );
    quad = pyramid_unit_volume ( ) * r8vec_dot_product ( order, w, v );
    printf ( "  %6d  %14.6g\n", order, quad );
    free ( v );
    free ( w );
    free ( xyz );

    order = 9;
    w = ( double * ) malloc ( order * sizeof ( double ) );
    xyz = ( double * ) malloc ( dim_num * order * sizeof ( double ) );
    pyramid_unit_o09 ( w, xyz );
    v = monomial_value ( dim_num, order, expon, xyz );
    quad = pyramid_unit_volume ( ) * r8vec_dot_product ( order, w, v );
    printf ( "  %6d  %14.6g\n", order, quad );
    free ( v );
    free ( w );
    free ( xyz );

    order = 13;
    w = ( double * ) malloc ( order * sizeof ( double ) );
    xyz = ( double * ) malloc ( dim_num * order * sizeof ( double ) );
    pyramid_unit_o13 ( w, xyz );
    v = monomial_value ( dim_num, order, expon, xyz );
    quad = pyramid_unit_volume ( ) * r8vec_dot_product ( order, w, v );
    printf ( "  %6d  %14.6g\n", order, quad );
    free ( v );
    free ( w );
    free ( xyz );

    order = 18;
    w = ( double * ) malloc ( order * sizeof ( double ) );
    xyz = ( double * ) malloc ( dim_num * order * sizeof ( double ) );
    pyramid_unit_o18 ( w, xyz );
    v = monomial_value ( dim_num, order, expon, xyz );
    quad = pyramid_unit_volume ( ) * r8vec_dot_product ( order, w, v );
    printf ( "  %6d  %14.6g\n", order, quad );
    free ( v );
    free ( w );
    free ( xyz );

    order = 27;
    w = ( double * ) malloc ( order * sizeof ( double ) );
    xyz = ( double * ) malloc ( dim_num * order * sizeof ( double ) );
    pyramid_unit_o27 ( w, xyz );
    v = monomial_value ( dim_num, order, expon, xyz );
    quad = pyramid_unit_volume ( ) * r8vec_dot_product ( order, w, v );
    printf ( "  %6d  %14.6g\n", order, quad );
    free ( v );
    free ( w );
    free ( xyz );

    order = 48;
    w = ( double * ) malloc ( order * sizeof ( double ) );
    xyz = ( double * ) malloc ( dim_num * order * sizeof ( double ) );
    pyramid_unit_o48 ( w, xyz );
    v = monomial_value ( dim_num, order, expon, xyz );
    quad = pyramid_unit_volume ( ) * r8vec_dot_product ( order, w, v );
    printf ( "  %6d  %14.6g\n", order, quad );
    free ( v );
    free ( w );
    free ( xyz );

    printf ( "\n" );
    quad = pyramid_unit_monomial ( expon );
    printf ( "   Exact  %14.6g\n", quad );

    if ( !more )
    {
      break;
    }
  }

  return;
# undef DIM_NUM
}