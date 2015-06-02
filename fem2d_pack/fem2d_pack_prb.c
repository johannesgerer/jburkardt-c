# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "fem2d_pack.h"

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );
void test04 ( );
void test05 ( );
void test07 ( );
void test08 ( );
void test09 ( );
void test10 ( );
void test105 ( );
void test11 ( );
void test12 ( );
void test13 ( );
void test135 ( );
void test14 ( );
void test15 ( );
void test16 ( );
void test18 ( );
void test19 ( );
void test20 ( );
void test21 ( );
void test22 ( );
void test23 ( );
void test24 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for FEM2D_PACK_PRB.

  Discussion:

    FEM2D_PACK_PRB tests the FEM2D_PACK library.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    22 January 2013

  Author:

    John Burkardt
*/
{
  int i;

  timestamp ( );
  printf ( "\n" );
  printf ( "FEM2D_PACK_PRB:\n" );
  printf ( "  C version\n" );
  printf ( "  Test the FEM2D_PACK library.\n" );
/*
  test01 ( );
  test02 ( );
*/
  test03 ( );
  test04 ( );
/*
  test05 ( );
  test07 ( );
  test08 ( );
  test09 ( );

  test10 ( );
*/
  test105 ( );
/*
  test11 ( );
*/
  test12 ( );
  test13 ( );
  test135 ( );
  test14 ( );
  test15 ( );
  test16 ( );
  test18 ( );
/*
  test19 ( );

  test20 ( );
  test21 ( );
  test22 ( );
  test23 ( );
*/
  test24 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "FEM2D_PACK_PRB:\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void test03 ( )

/******************************************************************************/
/*
  Purpose:

    TEST03 tests BASIS_11_**_TEST.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 January 2013

  Author:

    John Burkardt
*/
{
  printf ( "\n" );
  printf ( "TEST03\n" );
  printf ( "  BASIS_11_T3_TEST - Test the T3 basis functions.\n" );
  printf ( "  BASIS_11_T4_TEST - Test the T4 basis functions.\n" );
  printf ( "  BASIS_11_T6_TEST - Test the T6 basis functions.\n" );

  basis_11_t3_test ( );

  basis_11_t4_test ( );

  basis_11_t6_test ( );

  return;
}
/******************************************************************************/

void test04 ( )

/******************************************************************************/
/*
  Purpose:

    TEST04 tests BASIS_MN_**_TEST.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 January 2013

  Author:

    John Burkardt
*/
{
  printf ( "\n" );
  printf ( "TEST04\n" );
  printf ( "  Test the computation of basis functions by evaluating them\n" );
  printf ( "  at the nodes that define the basis functions.\n" );
  printf ( "\n" );
  printf ( "  BASIS_MN_Q4_TEST - for the Q4 element.\n" );
  printf ( "  BASIS_MN_T3_TEST - for the T3 element.\n" );
  printf ( "  BASIS_MN_T4_TEST - for the T4 element.\n" );
  printf ( "  BASIS_MN_T6_TEST - for the T6 element.\n" );

  basis_mn_q4_test ( );

  basis_mn_t3_test ( );

  basis_mn_t4_test ( );

  basis_mn_t6_test ( );

  return;
}
/******************************************************************************/

void test105 ( )

/******************************************************************************/
/*
  Purpose:

    TEST105 tests GRID_NODES_01.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    24 January 2013

  Author:

    John Burkardt
*/
{
  int node;
  int node_num;
  double *node_xy;
  int num_x = 5;
  int num_y = 3;

  printf ( "\n" );
  printf ( "TEST105\n" );
  printf ( "  GRID_NODES_01 computes a regular grid in the unit square.\n" );
  printf ( "\n" );
  printf ( "  NUM_X =    %d\n", num_x );
  printf ( "  NUM_Y =    %d\n", num_y );
  node_num = num_x * num_y;
  printf ( "  NODE_NUM = %d\n", node_num );
  printf ( "\n" );

  node_xy = grid_nodes_01 ( num_x, num_y );

  for ( node = 0; node < node_num; node++ )
  {
    printf ( "  %8d  %14g  %14g\n", node, node_xy[0+node*2], node_xy[1+node*2] );
  }

  free ( node_xy );

  return;
}
/******************************************************************************/

void test12 ( )

/******************************************************************************/
/*
  Purpose:

    TEST12 tests INTERP_TEST.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    22 January 2013

  Author:

    John Burkardt
*/
{
  printf ( "\n" );
  printf ( "TEST12\n" );
  printf ( "  INTERP_TEST tests the interpolating power\n" );
  printf ( "  of the element.\n" );

  interp_test ( "Q4" );

  interp_test ( "Q8" );

  interp_test ( "Q9" );

  interp_test ( "Q12" );

  interp_test ( "Q16" );

  interp_test ( "QL" );

  interp_test ( "T3" );

  interp_test ( "T4" );

  interp_test ( "T6" );

  interp_test ( "T10" );

  return;
}
/******************************************************************************/

void test13 ( )

/******************************************************************************/
/*
  Purpose:

    TEST13 tests MAP_TEST.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    22 January 2013

  Author:

    John Burkardt
*/
{
  printf ( "\n" );
  printf ( "TEST13\n" );
  printf ( "  MAP_TEST tests the map routines.\n" );

  map_test ( "Q4" );

  map_test ( "Q8" );

  map_test ( "Q9" );

  map_test ( "Q12" );

  map_test ( "Q16" );

  map_test ( "QL" );

  map_test ( "T3" );

  map_test ( "T6" );

  map_test ( "T10" );

  return;
}
/******************************************************************************/

void test135 ( )

/******************************************************************************/
/*
  Purpose:

    TEST135 tests MASS_MATRIX_T3.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 January 2013

  Author:

    John Burkardt
*/
{
# define ELEMENT_NUM 8
# define NODE_NUM 9

  double *a;
  int element_node[3*ELEMENT_NUM] = {
    1, 4, 2, 
    5, 2, 4, 
    4, 7, 5, 
    8, 5, 7, 
    2, 5, 3, 
    6, 3, 5, 
    5, 8, 6, 
    9, 6, 8 };
  double node_xy[2*NODE_NUM] = {
    0.0, 0.0,
    0.0, 0.5,
    0.0, 1.0,
    0.5, 0.0,
    0.5, 0.5,
    0.5, 1.0,
    1.0, 0.0,
    1.0, 0.5,
    1.0, 1.0 };

  printf ( "\n" );
  printf ( "TEST135\n" );
  printf ( "  MASS_MATRIX_T3 computes the mass matrix for\n" );
  printf ( "  a finite element system using T3 elements\n" );
  printf ( "  (linear triangles).\n" );

  a = mass_matrix_t3 ( NODE_NUM, ELEMENT_NUM, element_node, node_xy );

  r8mat_print ( NODE_NUM, NODE_NUM, a, "  The T3 mass matrix:" );

  free ( a );

  return;
# undef ELEMENT_NUM
# undef NODE_NUM
}
/******************************************************************************/

void test14 ( )

/******************************************************************************/
/*
  Purpose:

    TEST14 tests MASS_MATRIX_T6.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 January 2013

  Author:

    John Burkardt
*/
{
# define ELEMENT_NUM 2
# define NODE_NUM 9

  double *a;
  int element_node[6*ELEMENT_NUM] = {
    1, 3, 7, 2, 5, 4,
    9, 7, 3, 8, 5, 6 };
  double node_xy[2*NODE_NUM] = {
    0.0, 0.0,
    0.0, 0.5,
    0.0, 1.0,
    0.5, 0.0,
    0.5, 0.5,
    0.5, 1.0,
    1.0, 0.0,
    1.0, 0.5,
    1.0, 1.0 };

  printf ( "\n" );
  printf ( "TEST14\n" );
  printf ( "  MASS_MATRIX_T6 computes the mass matrix for\n" );
  printf ( "  a finite element system using T6 elements\n" );
  printf ( "  (quadratic triangles).\n" );

  a = mass_matrix_t6 ( NODE_NUM, ELEMENT_NUM, element_node, node_xy );

  r8mat_print ( NODE_NUM, NODE_NUM, a, "  The T6 mass matrix:" );

  free ( a );

  return;
# undef ELEMENT_NUM
# undef NODE_NUM
}
/******************************************************************************/

void test15 ( )

/******************************************************************************/
/*
  Purpose:

    TEST15 tests PHYSICAL_TO_REFERENCE_T3 and REFERENCE_TO_PHYSICAL_T3.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    15 January 2013

  Author:

    John Burkardt
*/
{
# define N 10

  int i;
  int j;
  double phy[2*N];
  double ref[2*N];
  double ref2[2*N];
  int seed;
  double t[2*3] = {
    1.0, 1.0, 
    3.0, 1.0, 
    2.0, 5.0 };

  seed = 123456789;

  printf ( "\n" );
  printf ( "TEST15\n" );
  printf ( "  For an order 3 triangle,\n" );
  printf ( "  PHYSICAL_TO_REFERENCE_T3 maps a physical point to\n" );
  printf ( "  a reference point.\n" );
  printf ( "  REFERENCE_TO_PHYSICAL_T3 maps a reference point to\n" );
  printf ( "  a physical point.\n" );
  printf ( "\n" );
  printf ( "      XSI     ETA  ==>  X       Y    ==>  XSI2    ETA2\n" );
  printf ( "\n" );

  for ( j = 0; j < N; j++ )
  {
    for ( i = 0; i < 2; i++ )
    {
      ref[i+j*2] = r8_uniform_01 ( &seed );
    }

    if ( 1.0 < ref[0+j*2] + ref[1+j*2] )
    {
      ref[0+j*2] = 1.0 - ref[0+j*2];
      ref[1+j*2] = 1.0 - ref[1+j*2];
    }
  }

  reference_to_physical_t3 ( t, N, ref, phy );

  physical_to_reference_t3 ( t, N, phy, ref2 );

  for ( j = 0; j < N; j++ )
  {
    printf ( "  %10f  %10f  %10f  %10f  %10f  %10f\n",
     ref[0+j*2], ref[1+j*2], phy[0+j*2], phy[1+j*2], ref2[0+j*2], ref2[1+j*2] );
  }

  return;
# undef N
}
/******************************************************************************/

void test16 ( )

/******************************************************************************/
/*
  Purpose:

    TEST16 tests REFERENCE_TO_PHYSICAL_T6.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    08 February 2006

  Author:

    John Burkardt
*/
{
# define N 16

  int i;
  int j;
  double phy[2*N];
  double ref[2*N] = {
    0.00, 0.00, 
    1.00, 0.00, 
    0.00, 1.00, 
    0.50, 0.00, 
    0.50, 0.50, 
    0.00, 0.50, 
    0.25, 0.75, 
    0.75, 0.25, 
    0.40, 0.10, 
    0.30, 0.20, 
    0.20, 0.30, 
    0.10, 0.40, 
    0.10, 0.10, 
    0.20, 0.20, 
    0.30, 0.30, 
    0.40, 0.40 };
  double t[2*6] = {
    0.0, 0.0, 
    2.0, 0.0, 
    0.0, 4.0, 
    1.0, 0.0, 
    1.0, 1.0, 
    0.0, 2.0 };

  printf ( "\n" );
  printf ( "TEST16\n" );
  printf ( "  For an order 6 triangle,\n" );
  printf ( "  REFERENCE_TO_PHYSICAL_T6 maps a reference point to\n" );
  printf ( "  a physical point.\n" );
  printf ( "\n" );
  printf ( "      XSI     ETA  ==>  X       Y\n" );
  printf ( "\n" );

  reference_to_physical_t6 ( t, N, ref, phy );

  for ( j = 0; j < N; j++ )
  {
    printf ( "  %8f  %8f  %8f  %8f\n", ref[0+j*2], ref[1+j*2], phy[0+j*2], phy[1+j*2] );
  }

  return;
# undef N
}
/******************************************************************************/

void test18 ( )

/******************************************************************************/
/*
  Purpose:

    TEST18 tests the shape routines.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 January 2013

  Author:

    John Burkardt
*/
{
  printf ( "\n" );
  printf ( "TEST18\n" );
  printf ( "  SHAPE_TEST tests the shape routines.\n" );

  shape_test ( "Q4" );

  shape_test ( "Q8" );

  shape_test ( "Q9" );

  shape_test ( "Q12" );

  shape_test ( "Q16" );

  shape_test ( "QL" );

  shape_test ( "T3" );

  shape_test ( "T6" );

  shape_test ( "T10" );

  return;
}
/******************************************************************************/

void test24 ( )

/******************************************************************************/
/*
  Purpose:

    TEST24 tests TRIANGLE_UNIT_SET.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    18 February 2013

  Author:

    John Burkardt
*/
{
# define ORDER_MAX 64

  int a;
  int b;
  double coef;
  double err;
  double exact;
  int i;
  int order;
  double quad;
  int rule;
  int rule_max = 20;
  double value;
  double weight[ORDER_MAX];
  double x;
  double xtab[ORDER_MAX];
  double y;
  double ytab[ORDER_MAX];

  printf ( "\n" );
  printf ( "TEST24\n" );
  printf ( "  TRIANGLE_UNIT_SET sets up a quadrature\n" );
  printf ( "  in the unit triangle.\n" );
  printf ( "\n" );

  for ( a = 0; a <= 10; a++ )
  {
    for ( b = 0; b <= 10 - a; b++ )
    {
      coef = ( double ) ( a + b + 2 ) * ( double ) ( a + b + 1 );
      for ( i = 1; i <= b; i++ )
      {
        coef = coef * ( double ) ( a + i ) / ( double ) ( i );
      }

      printf ( "\n" );
      printf ( "  A = %d  B = %d\n", a, b );
      printf ( "\n" );
      printf ( "  Rule       QUAD           ERROR\n" );
      printf ( "\n" );

      for ( rule = 1; rule <= rule_max; rule++ )
      {
        order = triangle_unit_size ( rule );

        triangle_unit_set ( rule, order, xtab, ytab, weight );

        quad = 0.0;

        for ( i = 0; i < order; i++ )
        {
          x = xtab[i];
          y = ytab[i];

          if ( a == 0 && b == 0 )
          {
            value = coef;
          }
          else if ( a == 0 && b != 0 )
          {
            value = coef * pow ( ytab[i], b );
          }
          else if ( a != 0 && b == 0 )
          {
            value = coef * pow ( xtab[i], a );
          }
          else if ( a != 0 && b != 0 )
          {
            value = coef * pow ( xtab[i], a ) * pow ( ytab[i], b );
          }

          quad = quad + 0.5 * weight[i] * value;

        }

        exact = 1.0;
        err = fabs ( exact - quad );

        printf ( "  %4d  %14g  %14g\n", rule, quad, err );
      }
    }
  }

  return;
# undef ORDER_MAX
}
