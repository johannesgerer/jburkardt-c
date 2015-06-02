# include <stdlib.h>
# include <stdio.h>
# include <time.h>
# include <string.h>
# include <math.h>

# include "wedge_felippa_rule.h"

int main ( );
void test01 ( int degree_max );
void test02 ( int degree_max );
void test03 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for WEDGE_FELIPPA_RULE_PRB.

  Discussion:

    FELIPPA_PRB tests the FELIPPA library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 August 2014

  Author:

    John Burkardt
*/
{
  int degree_max;

  timestamp ( );
  printf ( "\n" );
  printf ( "WEDGE_FELIPPA_RULE_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the WEDGE_FELIPPA_RULE library.\n" );

  degree_max = 4;

  test01 ( degree_max );
  test02 ( degree_max );

  test03 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "WEDGE_FELIPPA_RULE_PRB\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void test01 ( int degree_max )

/******************************************************************************/
/*
  Purpose:

    TEST01 tests WEDGE_INTEGRAL.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 August 2014

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
  printf ( "TEST01\n" );
  printf ( "  For the unit wedge,\n" );
  printf ( "  WEDGE_INTEGRAL returns the exact value of the\n" );
  printf ( "  integral of X^ALPHA Y^BETA Z^GAMMA\n" );
  printf ( "\n" );
  printf ( "  Volume = %g\n", wedge_volume ( ) );
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
        value = wedge_integral ( expon );
        printf ( "  %8d  %8d  %8d  %14.6g\n",
          expon[0], expon[1], expon[2], value );
      }
    }
  }

  return;
}
/******************************************************************************/

void test02 ( int degree_max )

/******************************************************************************/
/*
  Purpose:

    TEST02 tests the rules for the unit wedge.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 August 2014

  Author:

    John Burkardt

  Parameters:

    Input, int DEGREE_MAX, the maximum total degree of the
    monomials to check.
*/
{
  int dim_num = 3;
  int expon[3];
  int h;
  int line_order;
  int line_order_array[7] = {
    1, 2, 2, 3, 2, 3, 4 };
  int more;
  int order;
  double quad;
  int t;
  int test;
  int test_num = 7;
  int triangle_order;
  int triangle_order_index;
  int triangle_order_array[7] = {
    1, 3, -3, 6, -6, 7, 12 };
  double *v;
  double *w;
  double *xyz;

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  For the unit wedge,\n" );
  printf ( "  we approximate monomial integrals with WEDG_UNIT_RULE.\n" );

  more = 0;

  for ( ; ; )
  {
    subcomp_next ( degree_max, dim_num, expon, &more, &h, &t );

    if ( ( expon[2] % 2 ) == 1 )
    {
      if ( !more )
      {
        break;
      }
      else
      {
        continue;
      }
    }

    printf ( "\n" );
    printf ( "  Monomial exponents:   %2d  %2d  %2d\n", expon[0], expon[1], expon[2] );
    printf ( "\n" );

    for ( test = 0; test < test_num; test++ )
    {
      line_order = line_order_array[test];
      triangle_order = triangle_order_array[test];

      order = line_order * fabs ( triangle_order );

      w = ( double * ) malloc ( order * sizeof ( double ) );
      xyz = ( double * ) malloc ( dim_num * order * sizeof ( double ) );
      wedge_rule ( line_order, triangle_order, w, xyz );
      v = monomial_value ( dim_num, order, expon, xyz );
      quad = wedge_volume ( ) * r8vec_dot_product ( order, w, v );
      printf ( "  %6d  %6d  %6d  %14.6g\n",
        triangle_order, line_order, order, quad );
      free ( v );
      free ( w );
      free ( xyz );
    }

    printf ( "\n" );
    quad = wedge_integral ( expon );
    printf ( "   Exact                  %14.6g\n", quad );

    if ( !more )
    {
      break;
    }

  }

  return;
}
/******************************************************************************/

void test03 ( )

/******************************************************************************/
/*
  Purpose:

    TEST03 writes out some rules for the unit wedge.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 August 2014

  Author:

    John Burkardt

  Parameters:
*/
{
  int dim_num = 3;
  int line_order;
  int line_order_array[7] = {
    1, 2, 2, 3, 2, 3, 4 };
  int order;
  int rule;
  int rule_num = 7;
  int triangle_order;
  int triangle_order_array[7] = {
    1, 3, -3, 6, -6, 7, 12 };
  double *w;
  char w_filename[255];
  double *x;
  char x_filename[255];

  printf ( "\n" );
  printf ( "TEST03\n" );
  printf ( "  For the unit wedge,\n" );
  printf ( "  write some rules to a file\n" );
  printf ( "\n" );
  printf ( "   Rule  Trig    Line   Total  W_File X_File\n" );
  printf ( "         Order   Order  Order\n" );
  printf ( "\n" );

  for ( rule = 0; rule < rule_num; rule++ )
  {
    if ( rule == 0 )
    {
      strcpy ( w_filename, "wedge_felippa_1x1_w.txt" );
      strcpy ( x_filename, "wedge_felippa_1x1_x.txt" );
    }
    else if ( rule == 1 )
    {
      strcpy ( w_filename, "wedge_felippa_3x2_w.txt" );
      strcpy ( x_filename, "wedge_felippa_3x2_x.txt" );
    }
    else if ( rule == 2 )
    {
      strcpy ( w_filename, "wedge_felippa_3bx2_w.txt" );
      strcpy ( x_filename, "wedge_felippa_3bx2_x.txt" );
    }
    else if ( rule == 3 )
    {
      strcpy ( w_filename, "wedge_felippa_6x3_w.txt" );
      strcpy ( x_filename, "wedge_felippa_6x3_x.txt" );
    }
    else if ( rule == 4 )
    {
      strcpy ( w_filename, "wedge_felippa_6bx2_w.txt" );
      strcpy ( x_filename, "wedge_felippa_6bx2_x.txt" );
    }
    else if ( rule == 5 )
    {
      strcpy ( w_filename, "wedge_felippa_7x3_w.txt" );
      strcpy ( x_filename, "wedge_felippa_7x3_x.txt" );
    }
    else if ( rule == 6 )
    {
      strcpy ( w_filename, "wedge_felippa_12x4_w.txt" );
      strcpy ( x_filename, "wedge_felippa_12x4_x.txt" );
    }

    line_order = line_order_array[rule];
    triangle_order = triangle_order_array[rule];

    order = line_order * fabs ( triangle_order );

    w = ( double * ) malloc ( order * sizeof ( double ) );
    x = ( double * ) malloc ( dim_num * order * sizeof ( double ) );
    wedge_rule ( line_order, triangle_order, w, x );
    r8mat_write ( w_filename, 1, order, w );
    r8mat_write ( x_filename, dim_num, order, x );
    printf ( "  %6d  %6d  %6d  %6d  %s  %s\n",
      rule, triangle_order, line_order, order, w_filename, x_filename );

    free ( w );
    free ( x );
  }

  return;
}
