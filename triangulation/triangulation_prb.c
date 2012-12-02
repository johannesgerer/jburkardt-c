# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <string.h>

# include "triangulation.h"

int main ( );

void test01 ( );
void test02 ( );
void test025 ( );
void test026 ( );
void test03 ( );
void test04 ( );
void test05 ( );
void test06 ( );
void test07 ( );
void test08 ( );
void test09 ( );

void test10 ( );
void test11 ( );
void test12 ( );
void test125 ( );
void test13 ( );
void test14 ( );
void test15 ( );
void test16 ( );
void test17 ( );
void test18 ( );
void test19 ( );

void test20 ( );
void test21 ( );
void test213 ( );
void quad_fun ( int n, double xy_vec[], double f_vec[] );
void test215 ( );
void test217 ( );
void test219 ( );
void test22 ( );
void test23 ( );
void test24 ( );
void test25 ( );
void test26 ( );
void test265 ( );
void test27 ( );

void test31 ( );
void test32 ( );
void test33 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for TRIANGULATION_PRB.

  Discussion:

    TRIANGULATION_PRB tests routines from the TRIANGULATION library.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    19 October 2012

  Author:

    John Burkardt
*/
{
  timestamp ( );

  printf ( "\n" );
  printf ( "TRIANGULATION_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the TRIANGULATION library.\n" );

  test01 ( );
  test02 ( );
  test025 ( );
  test026 ( );
  test03 ( );
  test04 ( );
  test05 ( );
  test06 ( );
  test07 ( );
  test08 ( );
  test09 ( );

  test10 ( );
  test11 ( );
  test12 ( );
  test125 ( );
  test13 ( );
  test14 ( );
  test15 ( );
  test16 ( );
  test17 ( );
  test18 ( );
  test19 ( );

  test20 ( );
  test21 ( );
  test213 ( );
  test215 ( );
  test217 ( );
  test219 ( );
  test22 ( );
  test23 ( );
  test24 ( );
  test25 ( );
  test26 ( );
  test265 ( );
  test27 ( );

  test31 ( );
  test32 ( );
  test33 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "TRIANGULATION_PRB\n" );
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

    TEST01 tests ALPHA_MEASURE.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    18 October 2012

  Author:

    John Burkardt
*/
{
  double alpha_ave;
  double alpha_area;
  double alpha_min;
  int hole_num;
  int node_num;
  double *node_xy;
  double quality;
  int triangle_num;
  int *triangle_node;
  int *triangle_neighbor;
  int triangle_order = 3;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  ALPHA_MEASURE returns the ALPHA measure of\n" );
  printf ( "  quality of a triangulation.\n" );
/*
  Get the sizes.
*/
  triangulation_order3_example1_size ( &node_num, &triangle_num, &hole_num );
/*
  Allocate space.
*/
  node_xy = ( double * ) malloc ( 2 * node_num * sizeof ( double ) );
  triangle_node = ( int * ) malloc ( triangle_order * triangle_num * sizeof ( int ) );
  triangle_neighbor = ( int * ) malloc ( 3 * triangle_num * sizeof ( int ) );
/*
  Get the triangulation data.
*/
  triangulation_order3_example1 ( node_num, triangle_num, node_xy, 
    triangle_node, triangle_neighbor );
/*
  Compute the triangulation quality.
*/
  alpha_measure ( node_num, node_xy, triangle_order, triangle_num,
    triangle_node, &alpha_min, &alpha_ave, &alpha_area );

  printf ( "\n" );
  printf ( "  ALPHA_MIN  = %g\n", alpha_min );
  printf ( "  ALPHA_AVE  = %g\n", alpha_ave );
  printf ( "  ALPHA_AREA = %g\n", alpha_area );
/*
  Free the memory.
*/
  free ( node_xy );
  free ( triangle_node );
  free ( triangle_neighbor );

  return;
}
/******************************************************************************/

void test02 ( )

/******************************************************************************/
/*
  Purpose:

    TEST02 tests AREA_MEASURE.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    21 June 2009

  Author:

    John Burkardt
*/
{
  double area_ave;
  double area_max;
  double area_min;
  double area_ratio;
  double area_std;
  int hole_num;
  int node_num;
  double *node_xy;
  double quality;
  int triangle_num;
  int *triangle_node;
  int *triangle_neighbor;
  int triangle_order = 3;

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  AREA_MEASURE returns the AREA measure of\n" );
  printf ( "  quality of a triangulation.\n" );
/*
  Get the sizes.
*/
  triangulation_order3_example1_size ( &node_num, &triangle_num, &hole_num );
/*
  Allocate space.
*/
  node_xy = ( double * ) malloc ( 2 * node_num * sizeof ( double ) );
  triangle_node = ( int * ) malloc ( triangle_order * triangle_num * sizeof ( int ) );
  triangle_neighbor = ( int * ) malloc ( 3 * triangle_num * sizeof ( int ) );
/*
  Get the triangulation data.
*/
  triangulation_order3_example1 ( node_num, triangle_num, node_xy, 
    triangle_node, triangle_neighbor );
/*
  Compute the triangulation quality.
*/
  area_measure ( node_num, node_xy, triangle_order, triangle_num,
    triangle_node, &area_min, &area_max, &area_ratio, &area_ave, &area_std );

  printf ( "\n" );
  printf ( "  AREA_MIN   = %g\n", area_min );
  printf ( "  AREA_MAX   = %g\n", area_max );
  printf ( "  AREA_RATIO = %g\n", area_ratio );
  printf ( "  AREA_AVE   = %g\n", area_ave );
  printf ( "  AREA_STD   = %g\n", area_std );
/*
  Free the memory.
*/
  free ( node_xy );
  free ( triangle_node );
  free ( triangle_neighbor );

  return;
}
/******************************************************************************/

void test025 ( )

/******************************************************************************/
/*
  Purpose:

    TEST025 tests DELAUNAY_SWAP_TEST.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    18 October 2012

  Author:

    John Burkardt
*/
{
  int node_num = 4;
  int triangle_num = 2;
  int triangle_order = 3;

  double alpha_area;
  double alpha_ave;
  double alpha_min_swapped;
  double alpha_min_unswapped;
  double node_xy[2*4];
  int seed = 123456789;
  int swap;
  int test;
  int test_num = 10;
  int triangle_node[3*2];
  int value;

  printf ( "\n" );
  printf ( "TEST025\n" );
  printf ( "  DELAUNAY_SWAP_TEST determines whether two triangles\n" );
  printf ( "  with a common edge need to \"swap\" diagonals.\n" );
  printf ( "  If swapping is indicated, then ALPHA_MIN should increase.\n" );
  printf ( "\n" );
  printf ( "  Swap   ALPHA_MIN   ALPHA_MIN\n" );
  printf ( "         Unswapped   Swapped\n" );
  printf ( "\n" );

  for ( test = 1; test <= test_num; test++ )
  {
/*
  Generate a random quadrilateral (1,2,3,4).
*/
    quad_convex_random ( &seed, node_xy );
/*
  Does it need swapping?
*/
    swap = delaunay_swap_test ( node_xy );
/*
  Compute ALPHA_MIN unswapped.
*/
    triangle_node[0+0*3] = 1;
    triangle_node[1+0*3] = 2;
    triangle_node[2+0*3] = 3;
    triangle_node[0+1*3] = 1;
    triangle_node[1+1*3] = 3;
    triangle_node[2+1*3] = 4;

    alpha_measure ( node_num, node_xy, triangle_order, triangle_num, 
      triangle_node, &alpha_min_unswapped, &alpha_ave, &alpha_area );
/*
  Compute ALPHA_MIN swapped.
*/
    triangle_node[0+0*3] = 1;
    triangle_node[1+0*3] = 2;
    triangle_node[2+0*3] = 4;
    triangle_node[0+1*3] = 2;
    triangle_node[1+1*3] = 3;
    triangle_node[2+1*3] = 4;

    alpha_measure ( node_num, node_xy, triangle_order, triangle_num, 
      triangle_node, &alpha_min_swapped, &alpha_ave, &alpha_area );

    if ( 0 )
    {
      r8mat_transpose_print ( 2, node_num, node_xy, "  Quadrilateral" );
    }

    printf ( "     %d  %g10  %g10\n", swap, alpha_min_unswapped, alpha_min_swapped );
  }

  return;
}
/******************************************************************************/

void test026 ( )

/******************************************************************************/
/*
  Purpose:

    TEST026 tests DIAEDG.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    25 June 2009

  Author:

    John Burkardt
*/
{
  int node_num = 4;
  int triangle_num = 2;
  int triangle_order = 3;

  double alpha_area;
  double alpha_ave;
  double alpha_min_swapped;
  double alpha_min_unswapped;
  double node_xy[2*4];
  int seed = 123456789;
  int swap;
  int test;
  int test_num = 10;
  int triangle_node[3*2];
  int value;

  printf ( "\n" );
  printf ( "TEST026\n" );
  printf ( "  DIAEDG determines whether two triangles\n" );
  printf ( "  with a common edge need to \"swap\" diagonals.\n" );
  printf ( "  If swapping is indicated, then ALPHA_MIN should increase.\n" );
  printf ( "\n" );
  printf ( "  Swap   ALPHA_MIN   ALPHA_MIN\n" );
  printf ( "         Unswapped   Swapped\n" );
  printf ( "\n" );

  for ( test = 1; test <= test_num; test++ )
  {
/*
  Generate a random quadrilateral (1,2,3,4).
*/
    quad_convex_random ( &seed, node_xy );
/*
  Does it need swapping?
*/
    value = diaedg ( 
      node_xy[0+0*2], node_xy[1+0*2], 
      node_xy[0+1*2], node_xy[1+1*2], 
      node_xy[0+2*2], node_xy[1+2*2], 
      node_xy[0+3*2], node_xy[1+3*2] );

    if ( value == 1 )
    {
      swap = 0;
    }
    else
    {
      swap = 1;
    }
/*
  Compute ALPHA_MIN unswapped.
*/
    triangle_node[0+0*3] = 1;
    triangle_node[1+0*3] = 2;
    triangle_node[2+0*3] = 3;
    triangle_node[0+1*3] = 1;
    triangle_node[1+1*3] = 3;
    triangle_node[2+1*3] = 4;

    alpha_measure ( node_num, node_xy, triangle_order, triangle_num, 
      triangle_node, &alpha_min_unswapped, &alpha_ave, &alpha_area );
/*
  Compute ALPHA_MIN swapped.
*/
    triangle_node[0+0*3] = 1;
    triangle_node[1+0*3] = 2;
    triangle_node[2+0*3] = 4;
    triangle_node[0+1*3] = 2;
    triangle_node[1+1*3] = 3;
    triangle_node[2+1*3] = 4;

    alpha_measure ( node_num, node_xy, triangle_order, triangle_num, 
      triangle_node, &alpha_min_swapped, &alpha_ave, &alpha_area );

    if ( 0 )
    {
      r8mat_transpose_print ( 2, node_num, node_xy, "  Quadrilateral" );
    }

    printf ( "    %d  %10g  %10g\n", swap, alpha_min_unswapped, alpha_min_swapped );
  }

  return;
}
/******************************************************************************/

void test03 ( )

/******************************************************************************/
/*
  Purpose:

    TEST03 tests NODE_MERGE.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    24 August 2006

  Author:

    John Burkardt
*/
{
# define DIM_NUM 2
# define NODE_NUM 15
# define TEST_NUM 4

  int node;
  int node_rep[NODE_NUM];
  double node_xy[DIM_NUM*NODE_NUM] = {
       0.0, 0.0, 
       1.0, 0.0,
       3.0, 0.0,
       4.0, 0.0,
       1.0, 1.0,
       4.0, 1.0,
       2.0, 2.0,
       3.0, 3.0,
       2.0, 3.5,
       0.5, 4.0,
       1.0, 4.0,
       1.5, 4.0,
       4.0, 4.0,
       1.0, 4.5,
       1.0, 4.5 };
  int rep;
  int rep_num;
  int test;
  double tolerance;
  double tolerance_test[TEST_NUM] = { 0.01, 0.75, 1.2, 1.5 };

  printf ( "\n" );
  printf ( "TEST03\n" );
  printf ( "  NODE_MERGE identifies groups of nodes\n" );
  printf ( "  that can be merged, with a given tolerance.\n" );

  r8mat_transpose_print ( DIM_NUM, NODE_NUM, node_xy, "  Node coordinates:" );

  for ( test = 0; test < TEST_NUM; test++ )
  {
    tolerance = tolerance_test[test];

    node_merge ( DIM_NUM, NODE_NUM, node_xy, tolerance, node_rep );

    printf ( "\n" );
    printf ( "  TOLERANCE = %g\n", tolerance );
    printf ( "\n" );
    printf ( "      Node  Representatives:\n" );
    printf ( "\n" );

    for ( node = 0; node < NODE_NUM; node++ )
    {
      printf ( "  %8d  %8d\n", node, node_rep[node] );
    }
/*
  Make a list of the node representatives.
*/
    printf ( "\n" );
    printf ( "      Rep   Coordinates:\n" );
    printf ( "\n" );

    i4vec_sort_heap_a ( NODE_NUM, node_rep );

    rep_num = 0;

    for ( node = 0; node < NODE_NUM; node++ )
    {
      if ( 1 <= node )
      {
        if ( node_rep[node-1] == node_rep[node] )
        {
          continue;
        }
      }

      rep = node_rep[node];

      printf ( "  %8d  %12g  %12g\n", 
        rep_num, node_xy[0+rep*2], node_xy[1+rep*2] );

      rep_num = rep_num + 1;
    }
  }

  return;
# undef DIM_NUM
# undef NODE_NUM
# undef TEST_NUM
}
/******************************************************************************/

void test04 ( )

/******************************************************************************/
/*
  Purpose:

    TEST04 tests NS_ADJ_COL_SET, NS_ADJ_COUNT and NS_ADJ_ROW_SET.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    26 September 2006

  Author:

    John Burkardt
*/
{
# define NODE_NUM 15
# define TRIANGLE_NUM 4
# define TRIANGLE_ORDER 6
# define VARIABLE_NUM 36

  int adj_col[VARIABLE_NUM+1];
  int adj_num;
  int *adj_row;
  char file_name[80] = "ns_triangulation.eps";
  int node;
  int node_show;
  int node_p_variable[NODE_NUM] = {
    3, -1,  8, -1, 13, 
   -1, -1, -1, -1, 
   24, -1, 29, 
   -1, -1, 
   36 };
  int node_u_variable[NODE_NUM] = {
    1,  4,  6,  9, 11, 
   14, 16, 18, 20, 
   22, 25, 27, 
   30, 32, 
   34 };
  int node_v_variable[NODE_NUM] = {
    2,  5,  7, 10, 12, 
   15, 17, 19, 21, 
   23, 26, 28, 
   31, 33, 
   35 };
  double node_xy[2*NODE_NUM] = {
   0.0, 0.0, 
   0.0, 1.0, 
   0.0, 2.0, 
   0.0, 3.0, 
   0.0, 4.0, 
   1.0, 0.0, 
   1.0, 1.0, 
   1.0, 2.0, 
   1.0, 3.0, 
   2.0, 0.0, 
   2.0, 1.0, 
   2.0, 2.0, 
   3.0, 0.0, 
   3.0, 1.0, 
   4.0, 0.0 };
  int num;
  int r;
  int rhi;
  int rlo;
  int triangle_neighbor[3*TRIANGLE_NUM] = {
    -1,  2, -1, 
     3,  1,  4, 
     2, -1, -1, 
    -1, -1,  2 };
  int triangle_node[TRIANGLE_ORDER*TRIANGLE_NUM] = {
     1, 10,  3,  6,  7,  2, 
    12,  3, 10,  8,  7, 11, 
     3, 12,  5,  8,  9,  4, 
    10, 15, 12, 13, 14, 11 };
  int triangle_show;
  int variable;

  printf ( "\n" );
  printf ( "TEST04\n" );
  printf ( "  For an order 3/order 6 Taylor Hood triangulation\n" );
  printf ( "  for Navier Stokes velocity and pressure,\n" );
  printf ( "  NS_ADJ_COUNT counts variable adjacencies\n" );
  printf ( "    and sets up the sparse compressed column\n" );
  printf ( "    column pointer array.\n" );
  printf ( "  NS_ADJ_COL_SET sets up the sparse compressed column\n" );
  printf ( "    COL vector.\n" );
  printf ( "  NS_ADJ_ROW_SET sets up the sparse compressed column\n" );
  printf ( "    ROW vector.\n" );
/*
  Plot the example.
*/
  node_show = 2;
  triangle_show = 2;

  triangulation_order6_plot ( file_name, NODE_NUM, node_xy, 
    TRIANGLE_NUM, triangle_node, node_show, triangle_show );
/*
  Get the count of the variable adjacencies.
  We don't really need to make this call, since the next
  call does the calculation as part of getting ADJ_COL.
*/
  printf ( "\n" );
  printf ( "  Number of variables is %d\n", VARIABLE_NUM );

  adj_num = ns_adj_count ( NODE_NUM, TRIANGLE_NUM, VARIABLE_NUM, triangle_node,
    triangle_neighbor, node_u_variable, node_v_variable, node_p_variable );

  printf ( "\n" );
  printf ( "  As computed by NS_ADJ_COUNT,\n" );
  printf ( "  Number of variable adjacency entries is %d\n", adj_num );
/*
  Get the count of the variable adjacencies and the COL vector.
*/
  adj_num = ns_adj_col_set ( NODE_NUM, TRIANGLE_NUM, VARIABLE_NUM, triangle_node,
    triangle_neighbor, node_u_variable, node_v_variable, node_p_variable, 
    adj_col );

  printf ( "\n" );
  printf ( "  As computed by NS_ADJ_COL_SET,\n" );
  printf ( "  Number of variable adjacency entries is %d\n", adj_num );

  printf ( "\n" );
  printf ( "  Variable adjacency pointers:\n" );
  printf ( "\n" );
  printf ( "  Variable     First      Last    Number\n" );
  printf ( "\n" );

  for ( variable = 0; variable < VARIABLE_NUM; variable++ )
  {
    num = adj_col[variable+1] - adj_col[variable];

    printf ( "  %8d  %8d  %8d  %8d\n", 
      variable+1, adj_col[variable], adj_col[variable+1]-1, num );
  }
/*
  Get the variable adjacencies.
*/
  adj_row = ( int * ) malloc ( adj_num * sizeof ( int ) );

  ns_adj_row_set ( NODE_NUM, TRIANGLE_NUM, VARIABLE_NUM, triangle_node, 
    triangle_neighbor, node_u_variable, node_v_variable, node_p_variable, 
    adj_num, adj_col, adj_row );
/*
  This is a huge array.  We only print out the beginning and end.
*/
  printf ( "\n" );
  printf ( "  Variable adjacency row entries:\n" );
  printf ( "  (Partial printout only)\n" );
  printf ( "\n" );
  printf ( "     Entry     Row       Col\n" );
  printf ( "\n" );

  for ( variable = 0; variable < VARIABLE_NUM; variable++ )
  {
    rlo = adj_col[variable]-1;
    rhi = adj_col[variable+1]-2;

    if ( variable <= 2 || VARIABLE_NUM - 4 <= variable )
    {
      printf ( "\n" );

      for ( r = rlo; r <= rhi; r++ )
      {
        printf ( "  %8d  %8d  %8d\n", r+1, adj_row[r], variable+1 );
      }
    }

    if ( variable == 2 )
    {
      printf ( "\n" );
      printf ( "  (SKIPPING MANY MANY ENTRIES...)\n" );
      printf ( "\n" );
    }

  }

  free ( adj_row );

  return;
# undef NODE_NUM
# undef TRIANGLE_NUM
# undef TRIANGLE_ORDER
# undef VARIABLE_NUM
}
/******************************************************************************/

void test05 ( )

/******************************************************************************/
/*
  Purpose:  

    TEST04 tests POINTS_DELAUNAY_NAIVE_2D.

  Diagram:

    !....3&11....
    !............
    !............
    X..9.........
    !.....5......
    !...........6
    !.4.2...10...
    !.....8...12.
    V............
    !..7.........
    !......1.....
    !............
    !............
    !----V----X--

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    24 August 2006

  Author:

    John Burkardt
*/
{
# define NODE_NUM 12
# define DIM_NUM 2

  int i;
  int triangle_num;
  int *triangle_node;
  double node_xy[DIM_NUM*NODE_NUM] = {
     7.0, 3.0,
     4.0,  7.0,
     5.0, 13.0,
     2.0,  7.0,
     6.0,  9.0,
    12.0, 10.0,
     3.0,  4.0,
     6.0,  6.0,
     3.0, 10.0,
     8.0,  7.0,
     5.0, 13.0,
    10.0,  6.0 };

  printf ( "\n" );
  printf ( "TEST05\n" );
  printf ( "  POINTS_DELAUNAY_NAIVE_2D computes the Delaunay\n" );
  printf ( "  triangulation of a set of nodes.\n" );

  r8mat_transpose_print ( DIM_NUM, NODE_NUM, node_xy, "  The nodes:" );

  triangle_node = points_delaunay_naive_2d ( NODE_NUM, node_xy, &triangle_num );

  printf ( "\n" );
  printf ( "  Number of triangles is TRIANGLE_NUM = %d\n", triangle_num );

  i4mat_transpose_print ( 3, triangle_num, triangle_node, 
    "  The Delaunay triangles:" );

  free ( triangle_node );

  return;
# undef NODE_NUM
# undef DIM_NUM
}
/******************************************************************************/

void test06 ( )

/******************************************************************************/
/*
  Purpose:

    TEST06 tests POINTS_HULL_2D.

  Diagram:

    !....3.......
    !............
    !..9.........
    !.....5......
    !...........6
    !.4.2...10...
    !.....8......
    !.........12.
    !..7.........
    !......1.....
    !............
    !............
    !-----------

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    24 August 2006

  Author:

    John Burkardt
*/
{
# define NODE_NUM 12
# define DIM_NUM 2

  int i;
  int j;
  int ival[NODE_NUM];
  int nval;
  double node_xy[DIM_NUM*NODE_NUM] = {
       7.0,  3.0, 
       4.0,  7.0, 
       5.0, 13.0, 
       2.0,  7.0, 
       6.0,  9.0, 
      12.0,  8.0, 
       3.0,  4.0, 
       6.0,  6.0, 
       3.0, 10.0, 
       8.0,  7.0, 
       5.0, 13.0, 
      10.0,  6.0 };

  printf ( "\n" );
  printf ( "TEST06\n" );
  printf ( "  POINTS_HULL_2D computes the convex hull\n" );
  printf ( "    of a set of nodes.\n" );

  r8mat_transpose_print ( DIM_NUM, NODE_NUM, node_xy, "  The nodes:" );

  points_hull_2d ( NODE_NUM, node_xy, &nval, ival );

  printf ( "\n" );
  printf ( "  The convex hull is formed by connecting:\n" );
  printf ( "\n" );
  for ( j = 0; j < nval; j++ )
  {
    printf ( "  %3d  %3d", j, ival[j] );
    for ( i = 0; i < DIM_NUM; i++ )
    {
      printf ( "  %14g", node_xy[i+(ival[j]-1)*DIM_NUM] );
    }
    printf ( "\n" );    
  }

  printf ( "\n" );
  printf ( "  The correct sequence of nodes is:\n" );
  printf ( "  4, 9, 3, 6, 12, 1, 7, (4).\n" );

  return;
# undef NODE_NUM
# undef DIM_NUM
}
/******************************************************************************/

void test07 ( )

/******************************************************************************/
/*
  Purpose:

    TEST07 tests Q_MEASURE.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    21 June 2009

  Author:

    John Burkardt
*/
{
  int hole_num;
  int node_num;
  double *node_xy;
  double q_area;
  double q_ave;
  double q_max;
  double q_min;
  int triangle_num;
  int *triangle_node;
  int *triangle_neighbor;
  int triangle_order = 3;

  printf ( "\n" );
  printf ( "TEST07\n" );
  printf ( "  Q_MEASURE returns the Q measure of\n" );
  printf ( "  quality of a triangulation.\n" );
/*
  Get the sizes.
*/
  triangulation_order3_example1_size ( &node_num, &triangle_num, &hole_num );
/*
  Allocate space.
*/
  node_xy = ( double * ) malloc ( 2 * node_num * sizeof ( double ) );
  triangle_node = ( int * ) malloc ( triangle_order * triangle_num * sizeof ( int ) );
  triangle_neighbor = ( int * ) malloc ( 3*triangle_num * sizeof ( int ) );
/*
  Get the triangulation data.
*/
  triangulation_order3_example1 ( node_num, triangle_num, node_xy, 
    triangle_node, triangle_neighbor );
/*
  Compute the triangulation quality.
*/
  q_measure ( node_num, node_xy, triangle_order, triangle_num,
    triangle_node, &q_min, &q_max, &q_ave, &q_area );

  printf ( "\n" );
  printf ( "  Q_MIN  = %g\n", q_min );
  printf ( "  Q_MAX  = %g\n", q_max );
  printf ( "  Q_AVE  = %g\n", q_ave );
  printf ( "  Q_AREA = %g\n", q_area );
/*
  Free the memory.
*/
  free ( node_xy );
  free ( triangle_node );
  free ( triangle_neighbor );

  return;
}
/******************************************************************************/

void test08 ( )

/******************************************************************************/
/*
  Purpose:

    TEST08 tests R8TRIS2.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    24 August 2006

  Author:

    John Burkardt
*/
{
# define NODE_NUM 9

  int error;
  double node_xy[NODE_NUM*2] = {
       0.0, 0.0,
       0.0, 1.0,
       0.2, 0.5,
       0.3, 0.6,
       0.4, 0.5,
       0.6, 0.4,
       0.6, 0.5,
       1.0, 0.0,
       1.0, 1.0 };
  int triangle_node[2*NODE_NUM*3];
  int triangle_neighbor[2*NODE_NUM*3];
  int triangle_num;

  printf ( "\n" );
  printf ( "TEST08\n" );
  printf ( "  R8RIS2 computes the Delaunay triangulation of a\n" );
  printf ( "  pointset in 2D.\n" );
/*
  Set up the Delaunay triangulation.
*/
  error = r8tris2 ( NODE_NUM, node_xy, &triangle_num, triangle_node,
    triangle_neighbor );

  if ( error == 0 ) 
  {
    printf ( "\n" );
    printf ( "  R8TRIS2 computed the Delaunay triangulation with no\n" );
    printf ( "  errors detected.\n" );
  }
  else
  {
    printf ( "\n" );
    printf ( "  R8TRIS2 detected an error condition of index %d\n", error );
    return;
  }

  triangulation_order3_print ( NODE_NUM, triangle_num, node_xy,
    triangle_node, triangle_neighbor );

  return;
# undef NODE_NUM
}

/******************************************************************************/

void test09 ( )

/******************************************************************************/
/*
  Purpose:

    TEST09 tests TRIANGLE_ORDER3_PHYSICAL_TO_REFERENCE.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    07 December 2006

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
  printf ( "TEST09\n" );
  printf ( "  For an order 3 triangle,\n" );
  printf ( "  TRIANGLE_ORDER3_PHYSICAL_TO_REFERENCE\n" );
  printf ( "    maps a physical point to a reference point.\n" );
  printf ( "  TRIANGLE_ORDER3_REFERENCE_TO_PHYSICAL \n" );
  printf ( "    maps a reference point to a physical point.\n" );
  printf ( "\n" );
  printf ( "   (  XSI    ETA ) ==> ( X      Y  ) ==> ( XSI2    ETA2 )\n" );
  printf ( "\n" );

  triangle_reference_sample ( N, &seed, ref );

  triangle_order3_reference_to_physical ( t, N, ref, phy );

  triangle_order3_physical_to_reference ( t, N, phy, ref2 );

  for ( j = 0; j < N; j++ )
  {
    printf ( "  %10g  %10g    %10g  %10g    %10g  %10g\n",
      ref[0+j*2], ref[1+j*2], phy[0+j*2], phy[1+j*2], ref2[0+j*2], ref2[1+j*2] );
  }

  return;
# undef N
}


/******************************************************************************/

void test10 ( )

/******************************************************************************/
/*
  Purpose:

    TEST10 tests TRIANGLE_ORDER6_PHYSICAL_TO_REFERENCE.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    07 December 2006

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
  double t[2*6] = {
    7.0, 2.0, 
    9.0, 2.0, 
    7.0, 3.0, 
    8.0, 2.0, 
    8.0, 2.5, 
    7.0, 2.5 };

  seed = 123456789;

  printf ( "\n" );
  printf ( "TEST10\n" );
  printf ( "  For an order 6 triangle,\n" );
  printf ( "  TRIANGLE_ORDER6_PHYSICAL_TO_REFERENCE\n" );
  printf ( "    maps a physical point to a reference point\n" );
  printf ( "  TRIANGLE_ORDER6_REFERENCE_TO_PHYSICAL\n" );
  printf ( "    maps a reference point to a physical point.\n" );
  printf ( "\n" );
  printf ( "   (  XSI    ETA ) ==> ( X      Y  ) ==> ( XSI2    ETA2 )\n" );
  printf ( "\n" );

  triangle_reference_sample ( N, &seed, ref );

  triangle_order6_reference_to_physical ( t, N, ref, phy );

  triangle_order6_physical_to_reference ( t, N, phy, ref2 );

  for ( j = 0; j < N; j++ )
  {
    printf ( "  %10g  %10g    %10g  %10g    %10g  %10g\n",
      ref[0+j*2], ref[1+j*2], phy[0+j*2], phy[1+j*2], ref2[0+j*2], ref2[1+j*2] );
  }

  return;
# undef N
}
/******************************************************************************/

void test11 ( )

/******************************************************************************/
/*
  Purpose:

    TEST11 tests TRIANGULATION_NODE_ORDER.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    29 August 2005

  Author:

    John Burkardt
*/
{
# define NODE_NUM 36
# define TRIANGLE_NUM 41
# define TRIANGLE_ORDER 3

  int *node_order;
  int triangle_node[3*TRIANGLE_NUM] = {
     1,  8,  7, 
     1,  2,  8, 
     2,  9,  8, 
     2,  3,  9, 
     3, 10,  9, 
     3,  4, 10, 
     4, 11, 10, 
     4,  5, 11, 
     5, 12, 11, 
     5,  6, 12, 
     7, 14, 13, 
     7,  8, 14, 
     8, 15, 14, 
     8,  9, 15, 
    11, 18, 17, 
    11, 12, 18, 
    13, 20, 19, 
    13, 14, 20, 
    14, 21, 20, 
    14, 15, 21, 
    15, 22, 21, 
    15, 16, 22, 
    16, 23, 22, 
    16, 17, 23, 
    17, 24, 23, 
    17, 18, 24, 
    19, 26, 25, 
    19, 20, 26, 
    21, 28, 27, 
    21, 22, 28, 
    25, 30, 29, 
    25, 26, 30, 
    26, 31, 30, 
    27, 32, 31, 
    27, 28, 32, 
    29, 34, 33, 
    29, 30, 34, 
    30, 35, 34, 
    30, 31, 35, 
    31, 36, 35, 
    31, 32, 36 };

  printf ( "\n" );
  printf ( "TEST11\n" );
  printf ( "  TRIANGULATION_NODE_ORDER computes the order\n" );
  printf ( "  of the nodes in a triangulation.\n" );

  node_order = triangulation_node_order ( TRIANGLE_ORDER, TRIANGLE_NUM,
    triangle_node, NODE_NUM );

  i4vec_print ( NODE_NUM, node_order, "  NODE ORDER:" );

  free ( node_order );

  return;
# undef NODE_NUM
# undef TRIANGLE_NUM
# undef TRIANGLE_ORDER
}
/******************************************************************************/

void test12 ( )

/******************************************************************************/
/*
  Purpose:

    TEST12 tests TRIANGULATION_ORDER3_ADJ_SET.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    24 August 2006

  Author:

    John Burkardt
*/
{
  int *adj;
  int adj_num;
  int *adj_col;
  int hole_num;
  int k;
  int node;
  int node_num;
  double *node_xy;
  int *triangle_neighbor;
  int *triangle_node;
  int triangle_order = 3;
  int triangle_num;

  printf ( "\n" );
  printf ( "TEST12\n" );
  printf ( "  For an order3 triangulation:\n" );
  printf ( "  TRIANGULATION_ORDER3_ADJ_COUNT counts adjacencies\n" );
  printf ( "  TRIANGULATION_ORDER3_ADJ_SET sets adjacencies.\n" );
/*
  Get the sizes.
*/
  triangulation_order3_example1_size ( &node_num, &triangle_num, &hole_num );

  adj_col = ( int * ) malloc ( ( node_num + 1 ) * sizeof ( int ) );
  node_xy = ( double * ) malloc ( 2 * node_num * sizeof ( double ) );
  triangle_neighbor = ( int * ) malloc ( 3*triangle_num * sizeof ( int ) );
  triangle_node = ( int * ) malloc ( triangle_order*triangle_num * sizeof ( int ) );
/*
  Get the example data.
*/
  triangulation_order3_example1 ( node_num, triangle_num, node_xy, 
    triangle_node, triangle_neighbor );
/*
  Get the count of the adjacencies.
*/
  adj_num = triangulation_order3_adj_count ( node_num, triangle_num, 
    triangle_node, triangle_neighbor, adj_col );

  printf ( "\n" );
  printf ( "  Number of adjacency entries is %d\n", adj_num );

  printf ( "\n" );
  printf ( "  Adjacency pointers:\n" );
  printf ( "\n" );
  for ( node = 1; node <= node_num; node++ )
  {
    printf ( "  %8d  %8d  %8d\n", node, adj_col[node-1], adj_col[node]-1 );
  }
/*
  Get the adjacencies.
*/
  adj = triangulation_order3_adj_set ( node_num, triangle_num, triangle_node,
    triangle_neighbor, adj_num, adj_col );
/*
  Print the adjacencies.
*/
  for ( node = 1; node <= node_num; node++ )
  {
    printf ( "\n" );
    printf ( "  Nodes adjacent to node %d\n", node );
    printf ( "\n" );

    for ( k = adj_col[node-1]; k <= adj_col[node]-1; k++ )
    {
      printf ( "  %8d\n", adj[k-1] );
    }
  }

  free ( adj );
  free ( adj_col );
  free ( node_xy );
  free ( triangle_neighbor );
  free ( triangle_node );

  return;
}
/******************************************************************************/

void test125 ( )

/******************************************************************************/
/*
  Purpose:

    TEST125 tests TRIANGULATION_ORDER3_ADJ_SET2.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    24 August 2006

  Author:

    John Burkardt
*/
{
  int adj;
  int adj_num;
  int *adj_col;
  int hole_num;
  int *ia;
  int *ja;
  int k;
  int node;
  int node_num;
  double *node_xy;
  int *triangle_neighbor;
  int *triangle_node;
  int triangle_order = 3;
  int triangle_num;

  printf ( "\n" );
  printf ( "TEST125\n" );
  printf ( "  For an order3 triangulation:\n" );
  printf ( "  TRIANGULATION_ORDER3_ADJ_COUNT counts adjacencies\n" );
  printf ( "  TRIANGULATION_ORDER3_ADJ_SET2 sets adjacencies\n" );
  printf ( "  as a pair of vectors IA(*), JA(*).\n" );
/*
  Get the sizes.
*/
  triangulation_order3_example1_size ( &node_num, &triangle_num, &hole_num );

  adj_col = ( int * ) malloc ( ( node_num + 1 ) * sizeof ( int ) );
  node_xy = ( double * ) malloc ( 2 * node_num * sizeof ( double ) );
  triangle_neighbor = ( int * ) malloc ( 3 * triangle_num * sizeof ( int ) );
  triangle_node = ( int * ) malloc ( triangle_order * triangle_num * sizeof ( int ) );
/*
  Get the example data.
*/
  triangulation_order3_example1 ( node_num, triangle_num, node_xy, 
    triangle_node, triangle_neighbor );
/*
  Get the count of the adjacencies.
*/
  adj_num = triangulation_order3_adj_count ( node_num, triangle_num, 
    triangle_node, triangle_neighbor, adj_col );

  printf ( "\n" );
  printf ( "  Number of adjacency entries is %d\n", adj_num );

  printf ( "\n" );
  printf ( "  Adjacency pointers:\n" );
  printf ( "\n" );
  for ( node = 1; node <= node_num; node++ )
  {
    printf ( "  %8d  %8d  %8d\n", node, adj_col[node-1], adj_col[node]-1 );
  }
/*
  Get the adjacencies.
*/
  ia = ( int * ) malloc ( adj_num * sizeof ( int ) );
  ja = ( int * ) malloc ( adj_num * sizeof ( int ) );

  triangulation_order3_adj_set2 ( node_num, triangle_num, triangle_node,
    triangle_neighbor, adj_num, adj_col, ia, ja );
/*
  Print the adjacencies.
*/
  printf ( "\n" );
  printf ( "  Adjacency list:\n" );
  printf ( "\n" );
  for ( adj = 0; adj < adj_num; adj++ )
  {
    printf ( "  %8d  (%2d,%2d)\n", adj+1, ia[adj], ja[adj] );
  }

  free ( adj_col );
  free ( ia );
  free ( ja );
  free ( node_xy );
  free ( triangle_neighbor );
  free ( triangle_node );

  return;
}
/******************************************************************************/

void test13 ( )

/******************************************************************************/
/*
  Purpose:

    TEST13 tests TRIANGULATION_ORDER3_BOUNDARY_EDGE_COUNT.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    24 August 2006

  Author:

    John Burkardt
*/
{
# define DIM_NUM 2
# define NODE_NUM 36
# define TRIANGLE_NUM 41
# define TRIANGLE_ORDER 3

  int boundary_edge_num;
  char file_name[80] = "triangulation_order3_plot2.eps";
  int node_show = 2;
  double node_xy[DIM_NUM*NODE_NUM] = {
    0.0, 0.0, 
    1.0, 0.0, 
    2.0, 0.0, 
    3.0, 0.0, 
    4.0, 0.0, 
    5.0, 0.0, 
    0.0, 1.0, 
    1.0, 1.0, 
    2.0, 1.0, 
    3.0, 1.0, 
    4.0, 1.0, 
    5.0, 1.0, 
    0.0, 2.0, 
    1.0, 2.0, 
    2.0, 2.0, 
    3.0, 2.0, 
    4.0, 2.0, 
    5.0, 2.0, 
    0.0, 3.0, 
    1.0, 3.0, 
    2.0, 3.0, 
    3.0, 3.0, 
    4.0, 3.0, 
    5.0, 3.0, 
    0.0, 4.0, 
    1.0, 4.0, 
    2.0, 4.0, 
    3.0, 4.0, 
    0.0, 5.0, 
    1.0, 5.0, 
    2.0, 5.0, 
    3.0, 5.0, 
    0.0, 6.0, 
    1.0, 6.0, 
    2.0, 6.0, 
    3.0, 6.0 };
  int triangle_node[TRIANGLE_ORDER*TRIANGLE_NUM] = {
     1,  8,  7, 
     1,  2,  8, 
     2,  9,  8, 
     2,  3,  9, 
     3, 10,  9, 
     3,  4, 10, 
     4, 11, 10, 
     4,  5, 11, 
     5, 12, 11, 
     5,  6, 12, 
     7, 14, 13, 
     7,  8, 14, 
     8, 15, 14, 
     8,  9, 15, 
    11, 18, 17, 
    11, 12, 18, 
    13, 20, 19, 
    13, 14, 20, 
    14, 21, 20, 
    14, 15, 21, 
    15, 22, 21, 
    15, 16, 22, 
    16, 23, 22, 
    16, 17, 23, 
    17, 24, 23, 
    17, 18, 24, 
    19, 26, 25, 
    19, 20, 26, 
    21, 28, 27, 
    21, 22, 28, 
    25, 30, 29, 
    25, 26, 30, 
    26, 31, 30, 
    27, 32, 31, 
    27, 28, 32, 
    29, 34, 33, 
    29, 30, 34, 
    30, 35, 34, 
    30, 31, 35, 
    31, 36, 35, 
    31, 32, 36 };
  int triangle_show = 2;

  printf ( "\n" );
  printf ( "TEST13\n" );
  printf ( "  For an order3 triangulation:\n" );
  printf ( "  TRIANGULATION_ORDER3_BOUNDARY_EDGE_COUNT counts the\n" );
  printf ( "    boundary edges;\n" );
  printf ( "  TRIANGULATION_ORDER3_PLOT plots the triangulation.\n" );

  triangulation_order3_plot ( file_name, NODE_NUM, node_xy, 
    TRIANGLE_NUM, triangle_node, node_show, triangle_show );

  printf ( "\n" );
  printf ( "  An Encapsulated PostScript image of this\n" );
  printf ( "  triangulation is in \"%s\".\n", file_name );

  boundary_edge_num = triangulation_order3_boundary_edge_count ( 
    TRIANGLE_NUM, triangle_node );

  printf ( "\n" );
  printf ( "  Number of boundary edges = %d\n", boundary_edge_num );
  printf ( "  Correct number =           33\n" );

  return;
# undef DIM_NUM
# undef NODE_NUM
# undef TRIANGLE_NUM
# undef TRIANGLE_ORDER
}
/******************************************************************************/

void test14 ( )

/******************************************************************************/
/*
  Purpose:

    TEST14 tests TRIANGULATION_ORDER3_BOUNDARY_EDGE_COUNT_EULER.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    24 August 2006

  Author:

    John Burkardt
*/
{
  int boundary_num;
  int hole_num = 2;
  int node_num = 36;
  int triangle_num = 41;

  printf ( "\n" );
  printf ( "TEST14\n" );
  printf ( "  For an order3 triangulation:\n" );
  printf ( "  TRIANGULATION_ORDER3_BOUNDARY_EDGE_COUNT_EULER\n" );
  printf ( "  determines the number of edges that lie on the\n" );
  printf ( "  boundary of a region that has been triangulated.\n" );
  printf ( "\n" );
  printf ( "  Number of points =         %d\n", node_num );
  printf ( "  Number of triangles =      %d\n", triangle_num );
  printf ( "  Number of holes =          %d\n", hole_num );

  boundary_num = triangulation_order3_boundary_edge_count_euler ( node_num, 
    triangle_num, hole_num );

  printf ( "  Number of boundary edges = %d\n", boundary_num );

  return;
}
/******************************************************************************/

void test15 ( )

/******************************************************************************/
/*
  Purpose:

    TEST15 tests TRIANGULATION_ORDER3_BOUNDARY_NODE.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    24 August 2006

  Author:

    John Burkardt
*/
{
# define DIM_NUM 2
# define NODE_NUM 36
# define TRIANGLE_NUM 41
# define TRIANGLE_ORDER 3

  int i;
  int *node_boundary;
  int triangle_node[TRIANGLE_ORDER*TRIANGLE_NUM] = {
     1,  8,  7, 
     1,  2,  8, 
     2,  9,  8, 
     2,  3,  9, 
     3, 10,  9, 
     3,  4, 10, 
     4, 11, 10, 
     4,  5, 11, 
     5, 12, 11, 
     5,  6, 12, 
     7, 14, 13, 
     7,  8, 14, 
     8, 15, 14, 
     8,  9, 15, 
    11, 18, 17, 
    11, 12, 18, 
    13, 20, 19, 
    13, 14, 20, 
    14, 21, 20, 
    14, 15, 21, 
    15, 22, 21, 
    15, 16, 22, 
    16, 23, 22, 
    16, 17, 23, 
    17, 24, 23,
    17, 18, 24, 
    19, 26, 25, 
    19, 20, 26, 
    21, 28, 27, 
    21, 22, 28, 
    25, 30, 29, 
    25, 26, 30, 
    26, 31, 30, 
    27, 32, 31, 
    27, 28, 32, 
    29, 34, 33, 
    29, 30, 34, 
    30, 35, 34, 
    30, 31, 35, 
    31, 36, 35, 
    31, 32, 36 };

  printf ( "\n" );
  printf ( "TEST15\n" );
  printf ( "  For an order3 triangulation:\n" );
  printf ( "  TRIANGULATION_ORDER3_BOUNDARY_NODE determines which\n" );
  printf ( "  nodes lie on the boundary of a triangulation.\n" );

  node_boundary = triangulation_order3_boundary_node ( NODE_NUM, TRIANGLE_NUM, 
    triangle_node );

  lvec_print ( NODE_NUM, node_boundary, "    Node  BN?" );

  free ( node_boundary );

  return;
# undef DIM_NUM
# undef NODE_NUM
# undef TRIANGLE_NUM
# undef TRIANGLE_ORDER
}
/******************************************************************************/

void test16 ( )

/******************************************************************************/
/*
  Purpose:

    TEST16 tests TRIANGULATION_ORDER3_CHECK.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    24 August 2006

  Author:

    John Burkardt
*/
{
# define NODE_NUM 13
# define TRIANGLE_NUM 16
# define TRIANGLE_ORDER 3

  int i;
  int ierror;
  int isave;
  int node_num2;
  int triangle_num2;
  int triangle_node[TRIANGLE_ORDER*TRIANGLE_NUM] = {
     3,   4,   1, 
     3,   1,   2, 
     3,   2,   8, 
     2,   1,   5, 
     8,   2,  13, 
     8,  13,   9, 
     3,   8,   9, 
    13,   2,   5, 
     9,  13,   7, 
     7,  13,   5, 
     6,   7,   5, 
     9,   7,   6, 
    10,   9,   6, 
     6,   5,  12, 
    11,   6,  12, 
    10,   6,  11 };

  printf ( "\n" );
  printf ( "TEST16\n" );
  printf ( "  For an order3 triangulation:\n" );
  printf ( "  TRIANGULATION_ORDER3_CHECK checks the triangulation.\n" );

  i4mat_transpose_print ( TRIANGLE_ORDER, TRIANGLE_NUM, triangle_node, 
    "  Triangles:" );
/*
  Pass all tests.
*/
  ierror = triangulation_order3_check ( NODE_NUM, TRIANGLE_NUM, triangle_node );

  printf ( "  Error code = %d\n", ierror );
/*
  Fail test 1.
*/
  node_num2 = 2;

  ierror = triangulation_order3_check ( node_num2, TRIANGLE_NUM, 
    triangle_node );

  printf ( "  Error code = %d\n", ierror );
/*
  Fail test 2.
*/
  triangle_num2 = 0;

  ierror = triangulation_order3_check ( NODE_NUM, triangle_num2, 
    triangle_node );

  printf ( "  Error code = %d\n", ierror );
/*
  Fail test 3.
*/
  isave = triangle_node[1+4*3];
  triangle_node[1+4*3] = 0;

  ierror = triangulation_order3_check ( NODE_NUM, TRIANGLE_NUM, triangle_node );

  printf ( "  Error code = %d\n", ierror );
  triangle_node[1+4*3] = isave;
/*
  Fail test 4.
*/
  isave = triangle_node[2+9*3];
  triangle_node[2+9*3] = 2 * NODE_NUM + 1;

  ierror = triangulation_order3_check ( NODE_NUM, TRIANGLE_NUM, triangle_node );

  printf ( "  Error code = %d\n", ierror );
  triangle_node[2+9*3] = isave;
/*
  Fail test 5.
*/
  triangle_node[2+3*3] = 3;
  triangle_node[2+7*3] = 3;
  triangle_node[2+9*3] = 3;
  triangle_node[2+10*3] = 3;
  triangle_node[1+13*3] = 3;

  ierror = triangulation_order3_check ( NODE_NUM, TRIANGLE_NUM, triangle_node );

  printf ( "  Error code = %d\n", ierror );

  triangle_node[2+3*3] = 5;
  triangle_node[2+7*3] = 5;
  triangle_node[2+9*3] = 5;
  triangle_node[2+10*3] = 5;
  triangle_node[1+13*3] = 5;
/*
  Fail test 6.
*/
  triangle_node[0+8*3] = 7;

  ierror = triangulation_order3_check ( NODE_NUM, TRIANGLE_NUM, triangle_node );

  printf ( "  Error code = %d\n", ierror );
  triangle_node[0+8*3] = 9;
/*
  Fail test 7.
*/
  triangle_node[2+6*3] = 2;

  ierror = triangulation_order3_check ( NODE_NUM, TRIANGLE_NUM, triangle_node );

  printf ( "  Error code = %d\n", ierror );

  triangle_node[2+6*3] = 9;

  return;
# undef NODE_NUM
# undef TRIANGLE_NUM
# undef TRIANGLE_ORDER
}
/******************************************************************************/

void test17 ( )

/******************************************************************************/
/*
  Purpose:

    TEST17 tests TRIANGULATION_ORDER3_EXAMPLE1.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    24 August 2006

  Author:

    John Burkardt
*/
{
  int hole_num;
  int node_num;
  double *node_xy;
  int *triangle_node;
  int *triangle_neighbor;
  int triangle_order = 3;
  int triangle_num;

  printf ( "\n" );
  printf ( "TEST17\n" );
  printf ( "  For an order3 triangulation:\n" );
  printf ( "  TRIANGULATION_ORDER3_EXAMPLE1_SIZE gives the sizes\n" );
  printf ( "    for an example triangulation;\n" );
  printf ( "  TRIANGULATION_ORDER3_EXAMPLE1 returns the information\n" );
  printf ( "    for an example triangulation;\n" );
  printf ( "  TRIANGULATION_ORDER3_PRINT prints a triangulation.\n" );
/*
  Get the sizes.
*/
  triangulation_order3_example1_size ( &node_num, &triangle_num, &hole_num );

  node_xy = ( double * ) malloc ( 2 * node_num * sizeof ( double ) );
  triangle_node = ( int * ) malloc ( triangle_order * triangle_num * sizeof ( int ) );
  triangle_neighbor = ( int * ) malloc ( 3 * triangle_num * sizeof ( int ) );
/*
  Get the data.
*/
  triangulation_order3_example1 ( node_num, triangle_num, node_xy, 
    triangle_node, triangle_neighbor );

  triangulation_order3_print ( node_num, triangle_num, node_xy, 
    triangle_node, triangle_neighbor );

  free ( node_xy );
  free ( triangle_node );
  free ( triangle_neighbor );

  return;
}
/******************************************************************************/

void test18 ( )

/******************************************************************************/
/*
  Purpose:

    TEST18 tests TRIANGULATION_ORDER3_NEIGHBOR.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    24 August 2006

  Author:

    John Burkardt
*/
{
# define NODE_NUM 13
# define TRIANGLE_NUM 16
# define TRIANGLE_ORDER 3

  int s1;
  int s2;
  int t1;
  int t2;
  int triangle_node[TRIANGLE_ORDER*TRIANGLE_NUM] = {
     3,   4,   1, 
     3,   1,   2, 
     3,   2,   8, 
     2,   1,   5, 
     8,   2,  13, 
     8,  13,   9, 
     3,   8,   9, 
    13,   2,   5, 
     9,  13,   7, 
     7,  13,   5, 
     6,   7,   5, 
     9,   7,   6, 
    10,   9,   6, 
     6,   5,  12, 
    11,   6,  12, 
    10,   6,  11 };

  printf ( "\n" );
  printf ( "TEST18\n" );
  printf ( "  For an order3 triangulation:\n" );
  printf ( "  TRIANGULATION_ORDER3_NEIGHBOR determines the\n" );
  printf ( "  triangle neighbors.\n" );
  printf ( "\n" );
  printf ( "    T1    S1    T2    S2\n" );
  printf ( "\n" );

  for ( t1 = 1; t1 <= TRIANGLE_NUM; t1++ )
  {
    for ( s1 = 1; s1 <= 3; s1++ )
    {
      triangulation_order3_neighbor ( TRIANGLE_NUM, triangle_node, 
        t1, s1, &t2, &s2 );

      printf ( "  %4d  %4d  %4d  %4d\n", t1, s1, t2, s2 );
    }
  }

  return;
# undef NODE_NUM
# undef TRIANGLE_NUM
# undef TRIANGLE_ORDER
}
/******************************************************************************/

void test19 ( )

/******************************************************************************/
/*
  Purpose:

    TEST19 tests TRIANGULATION_NEIGHBOR_ELEMENTS.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    07 September 2009

  Author:

    John Burkardt
*/
{
# define TRIANGLE_NUM 16
# define TRIANGLE_ORDER 3

  int i;
  int triangle_node[TRIANGLE_ORDER*TRIANGLE_NUM] = {
     3,   4,   1,
     3,   1,   2,
     3,   2,   8,
     2,   1,   5,
     8,   2,  13,
     8,  13,   9,
     3,   8,   9,
    13,   2,   5,
     9,  13,   7,
     7,  13,   5,
     6,   7,   5,
     9,   7,   6,
    10,   9,   6,
     6,   5,  12,
    11,   6,  12,
    10,   6,  11 };
  int *triangle_neighbor;

  printf ( "\n" );
  printf ( "TEST19\n" );
  printf ( "  TRIANGULATION_NEIGHBOR_ELEMENTS determines the\n" );
  printf ( "  adjacency relationships between elements.\n" );

  i4mat_transpose_print ( TRIANGLE_ORDER, TRIANGLE_NUM, triangle_node, 
    "  Elements:" );

  triangle_neighbor = triangulation_neighbor_elements ( TRIANGLE_ORDER, 
    TRIANGLE_NUM, triangle_node );

  i4mat_transpose_print ( 3, TRIANGLE_NUM, triangle_neighbor, 
    "  Element neighbors:" );

  free ( triangle_neighbor );

  return;
# undef TRIANGLE_NUM
# undef TRIANGLE_ORDER
}
/******************************************************************************/

void test20 ( )

/******************************************************************************/
/*
  Purpose:

    TEST20 tests TRIANGULATION_ORDER3_PLOT.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    24 August 2006

  Author:

    John Burkardt
*/
{
  char file_name[80] = "triangulation_order3_plot.eps";
  int hole_num;
  int node_show = 0;
  int node_num;
  double *node_xy;
  int triangle_show = 2;
  int *triangle_neighbor;
  int *triangle_node;
  int triangle_num;
  int triangle_order = 3;

  printf ( "\n" );
  printf ( "TEST20\n" );
  printf ( "  For an order3 triangulation:\n" );
  printf ( "  TRIANGULATION_ORDER3_PLOT can plot a triangulation.\n" );
/*
  Get the sizes.
*/
  triangulation_order3_example1_size ( &node_num, &triangle_num, &hole_num );

  printf ( "\n" );
  printf ( "  Example data has %d points,\n", node_num );
  printf ( "  organized into %d triangles.\n", triangle_num );
/*
  Allocate space.
*/
  node_xy = ( double * ) malloc ( 2 * node_num * sizeof ( double ) );
  triangle_node = ( int * ) malloc ( triangle_order * triangle_num * sizeof ( int ) );
  triangle_neighbor = ( int * ) malloc ( 3 * triangle_num * sizeof ( int ) );
/*
  Get the example data.
*/
  triangulation_order3_example1 ( node_num, triangle_num, node_xy, 
    triangle_node, triangle_neighbor );
/*
  Make the plot.
*/
  triangulation_order3_plot ( file_name, node_num, node_xy, triangle_num, 
    triangle_node, node_show, triangle_show );

  printf ( "\n" );
  printf ( "  TRIANGULATION_ORDER3_PLOT has created an\n" );
  printf ( "  Encapsulated PostScript file (EPS) containing\n" );
  printf ( "  an image of the triangulation.\n" );
  printf ( "\n" );
  printf ( "  This file is called \"%s\".\n", file_name );

  free ( node_xy );
  free ( triangle_node );
  free ( triangle_neighbor );

  return;
}
/******************************************************************************/

void test21 ( )

/******************************************************************************/
/*
  Purpose:

    TEST21 tests TRIANGULATION_ORDER3_PRINT.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    24 August 2006

  Author:

    John Burkardt
*/
{
# define NODE_NUM 9
# define TRIANGLE_NUM 12
# define TRIANGLE_ORDER 3

  double node_xy[2*NODE_NUM] = {
       0.0, 0.0, 
       0.0, 1.0, 
       0.2, 0.5,
       0.3, 0.6, 
       0.4, 0.5, 
       0.6, 0.4, 
       0.6, 0.5, 
       1.0, 0.0, 
       1.0, 1.0 };
  int triangle_node[TRIANGLE_ORDER*TRIANGLE_NUM] = {
       2, 1, 3, 
       3, 1, 6, 
       2, 3, 4, 
       4, 3, 5, 
       7, 4, 5, 
       5, 3, 6, 
       7, 5, 6, 
       9, 4, 7, 
       6, 1, 8, 
       7, 6, 8, 
       7, 8, 9, 
       2, 4, 9 };
  int triangle_neighbor[3*TRIANGLE_NUM] = {
       -28,   2,  3, 
         1,   9,  6, 
         1,   4, 12, 
         3,   6,  5, 
         8,   4,  7, 
         4,   2,  7, 
         5,   6, 10, 
        12,   5, 11, 
         2, -34, 10, 
         7,   9, 11, 
        10, -38,  8, 
         3,   8, -3 };

  printf ( "\n" );
  printf ( "TEST21\n" );
  printf ( "  For an order3 triangulation:\n" );
  printf ( "  TRIANGULATION_ORDER3_PRINT prints out a triangulation.\n" );

  triangulation_order3_print ( NODE_NUM, TRIANGLE_NUM, node_xy,
    triangle_node, triangle_neighbor );

  return;
# undef NODE_NUM
# undef TRIANGLE_NUM
# undef TRIANGLE_ORDER
}
/******************************************************************************/

void test213 ( )

/******************************************************************************/
/*
  Purpose:

    TEST213 tests TRIANGULATION_ORDER3_QUAD.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    22 January 2007

  Author:

    John Burkardt
*/
{
# define QUAD_NUM 6

  int i;
  int j;
  int k;
  int n;
  int n11;
  int n12;
  int n21;
  int n22;
  int node_num;
  double *node_xy;
  double quad_value;
  double quad_w[QUAD_NUM] = {
    0.1666666666666666, 
    0.1666666666666666, 
    0.1666666666666666, 
    0.1666666666666666, 
    0.1666666666666666, 
    0.16666666666666660 };
  double quad_xy[2*QUAD_NUM] = {
    0.659027622374092,  0.231933368553031, 
    0.659027622374092,  0.109039009072877, 
    0.231933368553031,  0.659027622374092, 
    0.231933368553031,  0.109039009072877, 
    0.109039009072877,  0.659027622374092, 
    0.109039009072877,  0.231933368553031 };
  double region_area;
  int test;
  int *triangle_node;
  int test_num = 4;
  int triangle_order = 3;
  int triangle_num;
  double x;
  double y;

  printf ( "\n" );
  printf ( "TEST213\n" );
  printf ( "  TRIANGULATION_ORDER3_QUAD can apply a quadrature rule\n" );
  printf ( "  to every triangle in a triangulated region,\n" );
  printf ( "  and estimate the integral of a function over\n" );
  printf ( "  that region.\n" );
  printf ( "\n" );
  printf ( "  NODE_NUM   TRI_NUM  Integral estim  Area of Region\n" );
  printf ( "\n" );

  for ( test = 1; test <= test_num; test++ )
  {
/*
  Set up the grid.
*/
    n = i4_power ( 2, test - 1 );
    node_num = ( n + 1 ) * ( n + 1 );

    node_xy = ( double * ) malloc ( 2 * node_num * sizeof ( double ) );

    k = 0;
    for ( j = 1; j <= n + 1; j++ )
    {
      y = ( double ) ( j - 1 ) / ( double ) ( n + 1 - 1 );
      for ( i = 1; i <= n + 1; i++ )
      {
        x = ( double ) ( i - 1 ) / ( double ) ( n + 1 - 1 );
        node_xy[0+k*2] = x;
        node_xy[1+k*2] = y;
        k = k + 1;
      }
    }
/*
  Set up the triangulation.
*/
    triangle_num = 2 * n * n;

    triangle_node = ( int * ) malloc ( triangle_order * triangle_num * sizeof ( int ) );

    k = 0;
    for ( j = 1; j <= n; j++ )
    {
      for ( i = 1; i <= n; i++ )
      {
        n11 = i     + ( j     - 1 ) * ( n + 1 );
        n12 = i     + ( j + 1 - 1 ) * ( n + 1 );
        n21 = i + 1 + ( j     - 1 ) * ( n + 1 );
        n22 = i + 1 + ( j + 1 - 1 ) * ( n + 1 );

        triangle_node[0+k*3] = n11;
        triangle_node[1+k*3] = n21;
        triangle_node[2+k*3] = n12;
        k = k + 1;

        triangle_node[0+k*3] = n22;
        triangle_node[1+k*3] = n12;
        triangle_node[2+k*3] = n21;
        k = k + 1;
      }
    }
/*
  Estimate the integral.
*/
    triangulation_order3_quad ( node_num, node_xy, triangle_order, 
      triangle_num, triangle_node, &quad_fun, QUAD_NUM, quad_xy, quad_w, 
      &quad_value, &region_area );

    printf ( "  %8d  %8d  %14g  %14g\n", node_num, triangle_num, quad_value, region_area );
/*
  Delete allocatables.
*/
    free ( node_xy );
    free ( triangle_node );
  }

  return;;
# undef QUAD_NUM
}
/******************************************************************************/

void quad_fun ( int n, double xy_vec[], double f_vec[] )

/******************************************************************************/
/*
  Purpose:

    QUAD_FUN is a sample integrand function for TRIANGULATION_QUAD.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    22 January 2007

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of evaluation points.

    Input, double XY_VEC[2*N], the evaluation points.

    Output, double F_VEC[N], the value of the integrand
    function at the evaluation points.
*/
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    f_vec[i] = exp ( pow ( xy_vec[0+i*2], 2 ) 
                   + pow ( xy_vec[1+i*2], 2 ) );
  }
  return;
}
/******************************************************************************/

void test215 ( )

/******************************************************************************/
/*
  Purpose:

    TEST215 tests TRIANGULATION_ORDER3_REFINE_COMPUTE.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    28 January 2007

  Author:

    John Burkardt
*/
{
# define DIM_NUM 2
# define NODE_NUM1 5
# define TRIANGLE_NUM1 3
# define TRIANGLE_ORDER 3

  int *edge_data;
  int node_num2;
  double node_xy1[DIM_NUM*NODE_NUM1] = {
       0.0, 0.0, 
       1.0, 0.0, 
       0.0, 1.0, 
       1.0, 1.0, 
       0.5, 1.5 };
  double *node_xy2;
  int triangle_node1[TRIANGLE_ORDER*TRIANGLE_NUM1] = {
       1, 2, 3, 
       4, 3, 2, 
       3, 4, 5 };
  int *triangle_node2;
  int triangle_num2;

  printf ( "\n" );
  printf ( "TEST215\n" );
  printf ( "  For an order3 triangulation:\n" );
  printf ( "  TRIANGULATION_ORDER3_REFINE_SIZE determines the\n" );
  printf ( "  size of a refined triangulation.\n" );
  printf ( "  TRIANGULATION_ORDER3_REFINE_COMPUTES computes the\n" );
  printf ( "  refined triangulation.\n" );

  printf ( "\n" );
  printf ( "  The number of nodes is %d\n", NODE_NUM1 );
  printf ( "  The number of triangles is %d\n", TRIANGLE_NUM1 );

  r8mat_transpose_print ( DIM_NUM, NODE_NUM1, node_xy1, 
    "  The nodes" );

  i4mat_transpose_print ( TRIANGLE_ORDER, TRIANGLE_NUM1, triangle_node1, 
    "  The triangles:" );

  edge_data = ( int * ) malloc ( 5 * 3 * TRIANGLE_NUM1 * sizeof ( int ) );

  printf ( "\n" );
  printf ( "  Sizing the refined mesh:\n" );

  triangulation_order3_refine_size ( NODE_NUM1, TRIANGLE_NUM1, 
    triangle_node1, &node_num2, &triangle_num2, edge_data );

  printf ( "\n" );
  printf ( "  Information about the refined mesh:\n" );
  printf ( "\n" );
  printf ( "  The number of nodes is %d\n", node_num2 );
  printf ( "  The number of triangles is %d\n", triangle_num2 );

  printf ( "\n" );
  printf ( "  Computing the refined mesh:\n" );

  node_xy2 = ( double * ) malloc ( DIM_NUM * node_num2 * sizeof ( double ) );
  triangle_node2 = ( int * ) malloc ( TRIANGLE_ORDER * triangle_num2 * sizeof ( int ) );

  triangulation_order3_refine_compute ( NODE_NUM1, TRIANGLE_NUM1, 
    node_xy1, triangle_node1, node_num2, triangle_num2, edge_data, node_xy2, 
    triangle_node2 );

  r8mat_transpose_print ( DIM_NUM, node_num2, node_xy2, 
    "  The refined nodes" );

  i4mat_transpose_print ( TRIANGLE_ORDER, triangle_num2, triangle_node2, 
    "  The refined triangles:" );

  free ( edge_data );
  free ( node_xy2 );
  free ( triangle_node2 );

  return;
# undef DIM_NUM
# undef NODE_NUM1
# undef TRIANGLE_NUM1
# undef TRIANGLE_ORDER
}
/******************************************************************************/

void test217 ( )

/******************************************************************************/
/*
  Purpose:

    TEST217 tests TRIANGULATION_SEARCH_DELAUNAY.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    18 August 2009

  Author:

    John Burkardt
*/
{
# define DIM_NUM 2
# define NODE_NUM 13
# define TEST_NUM 10
# define TRIANGLE_ORDER 3

  double alpha;
  double beta;
  double d1;
  double d2;
  double d3;
  double dist;
  double dnear;
  int edge;
  int error;
  double gamma;
  int i;
  int i1;
  int i2;
  int i3;
  int nnear;
  double node_xy[DIM_NUM*NODE_NUM] = {
       0.0, 0.0,
       2.0, 2.0,
      -1.0, 3.0,
      -2.0, 2.0,
       8.0, 2.0,
       9.0, 5.0,
       7.0, 4.0,
       5.0, 6.0,
       6.0, 7.0,
       8.0, 8.0,
      11.0, 7.0,
      10.0, 4.0,
       6.0, 4.0 };
  double p[DIM_NUM];
  int seed;
  int step_num;
  int td[TEST_NUM];
  int test;
  int triangle_index;
  int triangle_neighbor[3*2*NODE_NUM];
  int triangle_node[TRIANGLE_ORDER*2*NODE_NUM];
  int triangle_num;
  double xd[DIM_NUM*TEST_NUM];

  printf ( "\n" );
  printf ( "TEST217\n" );
  printf ( "  Given a set of nodes NODE_XY, and a single point XD,\n" );
  printf ( "  find the nearest node in NODE_XY to XD.\n" );
  printf ( "\n" );
  printf ( "  POINTS_POINT_NEAR_NAIVE_ND uses a naive method.\n" );
  printf ( "  TRIANGULATION_SEARCH_DELAUNAY finds a triangle\n" );
  printf ( "    containing the point.  Often, one of these vertices\n" );
  printf ( "    is the closest point.\n" );
/*
  Set up the Delaunay triangulation.
*/
  error = r8tris2 ( NODE_NUM, node_xy, &triangle_num, triangle_node,
    triangle_neighbor );

  if ( error == 0 )
  {
    printf ( "\n" );
    printf ( "  R8TRIS2 computed the Delaunay triangulation.\n" );
  }
  else
  {
    printf ( "\n" );
    printf ( "  R8TRIS2 returned an error condition.\n" );
    exit ( 1 );
  }
/*
  Get the test points.
*/
  seed = 123456789;

  triangulation_order3_sample ( NODE_NUM, node_xy, triangle_num,
    triangle_node, TEST_NUM, &seed, xd, td );

  printf ( "\n" );
  printf ( "              X         Y     Distance    Index     Steps\n" );
  printf ( "\n" );

  for ( test = 0; test < TEST_NUM; test++ )
  {
    for ( i = 0; i < DIM_NUM; i++ )
    {
      p[i] = xd[i+test*DIM_NUM];
    }

    nnear = points_point_near_naive_nd ( DIM_NUM, NODE_NUM, node_xy, 
      p, &dnear );

    printf ( "\n" );
    printf ( "  XD       %8d  %8d\n", p[0], p[1] );
    printf ( "  Naive    %8f  %8f  %8f  %6d\n", 
      node_xy[0+(nnear-1)*DIM_NUM], node_xy[1+(nnear-1)*DIM_NUM], dnear, nnear );

    triangulation_search_delaunay ( NODE_NUM, node_xy, TRIANGLE_ORDER, triangle_num,
      triangle_node, triangle_neighbor, p, &triangle_index, &alpha, &beta,
      &gamma, &edge, &step_num );
   
    if ( triangle_index < 1 )
    {
      printf ( "\n" );
      printf ( "  Error: the search failed.\n" );
      continue;
    }

    i1 = triangle_node[0+(triangle_index-1)*TRIANGLE_ORDER];
    d1 = sqrt ( pow ( p[0] - node_xy[0+i1*2], 2 ) 
              + pow ( p[1] - node_xy[1+i1*2], 2 ) );

    dist = d1;
    nnear = i1;

    i2 = triangle_node[1+(triangle_index-1)*TRIANGLE_ORDER];
    d2 = sqrt ( pow ( p[0] - node_xy[0+i2*2], 2 ) 
              + pow ( p[1] - node_xy[1+i2*2], 2 ) );

    if ( d2 < dist )
    {
      dnear = d2;
      nnear = i2;
    }

    i3 = triangle_node[2+(triangle_index-1)*TRIANGLE_ORDER];
    d3 = sqrt ( pow ( p[0] - node_xy[0+i3*2], 2 ) 
              + pow ( p[1] - node_xy[1+i3*2], 2 ) );

    if ( d3 < dist )
    {
      dnear = d3;
      nnear = i3;
    }

    printf ( "  Delaunay %8g  %8g  %8g  %8d  %8d\n",
      node_xy[0+nnear*2], node_xy[1+nnear*2], dnear, nnear + 1, step_num );
  }

  return;
# undef DIM_NUM
# undef NODE_NUM
# undef TEST_NUM
# undef TRIANGLE_ORDER
}
/******************************************************************************/

void test219 ( )

/******************************************************************************/
/*
  Purpose:

    TEST219 tests TRIANGULATION_SEARCH_DELAUNAY, TRIANGULATION_SEARCH_NAIVE.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    01 August 2009

  Author:

    John Burkardt
*/
{
# define DIM_NUM 2
# define NODE_NUM 13
# define TEST_NUM 10
# define TRIANGLE_ORDER 3

  double alpha;
  double beta;
  int edge;
  int error;
  double gamma;
  int i;
  int nnear;
  double node_xy[DIM_NUM*NODE_NUM] = {
       0.0, 0.0,
       2.0, 2.0,
      -1.0, 3.0,
      -2.0, 2.0,
       8.0, 2.0,
       9.0, 5.0,
       7.0, 4.0,
       5.0, 6.0,
       6.0, 7.0,
       8.0, 8.0,
      11.0, 7.0,
      10.0, 4.0,
       6.0, 4.0 };
  double p_test[DIM_NUM*TEST_NUM];
  int seed = 123456789;
  int step_num;
  int t_test[TEST_NUM];
  int test;
  int triangle_index1;
  int triangle_index2;
  int triangle_neighbor[3*2*NODE_NUM];
  int triangle_node[TRIANGLE_ORDER*2*NODE_NUM];
  int triangle_num;

  printf ( "\n" );
  printf ( "TEST219\n" );
  printf ( "  Given a triangulation, and a point P,\n" );
  printf ( "  find the triangle T containing to P.\n" );
  printf ( "\n" );
  printf ( "  TRIANGULATION_SEARCH_NAIVE uses a naive method.\n" );
  printf ( "  TRIANGULATION_SEARCH_DELAUNAY uses a method that will work\n" );
  printf ( "    fast if the triangulation is Delaunay.\n" );
/*
  Set up the Delaunay triangulation.
*/
  error = r8tris2 ( NODE_NUM, node_xy, &triangle_num, triangle_node,
    triangle_neighbor );

  if ( error == 0 )
  {
    printf ( "\n" );
    printf ( "  R8TRIS2 computed the Delaunay triangulation.\n" );
  }
  else
  {
    printf ( "\n" );
    printf ( "  R8TRIS2 returned an error condition.\n" );
    exit ( 1 );
  }
/*
  Get the test points.
*/
  triangulation_order3_sample ( NODE_NUM, node_xy, triangle_num,
    triangle_node, TEST_NUM, &seed, p_test, t_test );

  printf ( "\n" );
  printf ( "         X           Y     Naive   Delaunay  Steps\n" );
  printf ( "\n" );

  for ( test = 0; test < TEST_NUM; test++ )
  {
    triangle_index1 = triangulation_search_naive ( NODE_NUM, node_xy, 
      TRIANGLE_ORDER, triangle_num, triangle_node, p_test+DIM_NUM*test );

    triangulation_search_delaunay ( NODE_NUM, node_xy, TRIANGLE_ORDER, 
      triangle_num, triangle_node, triangle_neighbor, p_test+DIM_NUM*test, 
      &triangle_index2, &alpha, &beta, &gamma, &edge, &step_num );

    printf ( "  %10g  %10g  %8d  %8d  %8d\n",
      p_test[0+test*DIM_NUM], p_test[1+test*DIM_NUM], triangle_index1, 
      triangle_index2, step_num );

  }

  return;
# undef DIM_NUM
# undef NODE_NUM
# undef TEST_NUM
# undef TRIANGLE_ORDER
}
/******************************************************************************/

void test22 ( )

/******************************************************************************/
/*
  Purpose:

    TEST22 tests TRIANGULATION_ORDER6_ADJ_SET.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    24 August 2006

  Author:

    John Burkardt
*/
{
  int *adj;
  int adj_num;
  int *adj_col;
  int hole_num;
  int k;
  int node;
  int node_num;
  double *node_xy;
  int *triangle_neighbor;
  int *triangle_node;
  int triangle_order = 6;
  int triangle_num;

  printf ( "\n" );
  printf ( "TEST22\n" );
  printf ( "  For an order6 triangulation:\n" );
  printf ( "  TRIANGULATION_ORDER6_ADJ_COUNT counts adjacencies\n" );
  printf ( "  TRIANGULATION_ORDER6_ADJ_SET sets adjacencies.\n" );
/*
  Get the sizes.
*/
  triangulation_order6_example1_size ( &node_num, &triangle_num, &hole_num );

  adj_col = ( int * ) malloc ( ( node_num + 1 ) * sizeof ( int ) );
  node_xy = ( double * ) malloc ( 2 * node_num * sizeof ( double ) );
  triangle_neighbor = ( int * ) malloc ( 3*triangle_num * sizeof ( int ) );
  triangle_node = ( int * ) malloc ( triangle_order * triangle_num * sizeof ( int ) );
/*
  Get the example data.
*/
  triangulation_order6_example1 ( node_num, triangle_num, node_xy, 
    triangle_node, triangle_neighbor );
/*
  Get the count of the adjacencies.
*/
  adj_num = triangulation_order6_adj_count ( node_num, triangle_num, 
    triangle_node, triangle_neighbor, adj_col );

  printf ( "\n" );
  printf ( "  Number of adjacency entries is %d\n", adj_num );

  printf ( "\n" );
  printf ( "  Adjacency pointers:\n" );
  printf ( "\n" );
  for ( node = 1; node <= node_num; node++ )
  {
    printf ( "  %8d  %8d  %8d\n", node, adj_col[node-1], adj_col[node]-1 );
  }
/*
  Get the adjacencies.
*/
adj = triangulation_order6_adj_set ( node_num, triangle_num, triangle_node,
  triangle_neighbor, adj_num, adj_col );
/*
  Print the adjacencies.
*/
  for ( node = 1; node <= node_num; node++ )
  {
    printf ( "\n" );
    printf ( "  Nodes adjacent to node %d\n", node );
    printf ( "\n" );

    for ( k = adj_col[node-1]; k <= adj_col[node]-1; k++ )
    {
      printf ( "  %8d\n", adj[k-1] );
    }
  }

  free ( adj );
  free ( adj_col );
  free ( node_xy );
  free ( triangle_neighbor );
  free ( triangle_node );

  return;
}
/******************************************************************************/

void test23 ( )

/******************************************************************************/
/*
  Purpose:

    TEST23 tests TRIANGULATION_ORDER6_BOUNDARY_EDGE_COUNT.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    24 August 2006

  Author:

    John Burkardt
*/
{
  int boundary_edge_num;
  int dim_num = 2;
  int hole_num;
  int node_num;
  double *node_xy;
  int *triangle_node;
  int *triangle_neighbor;
  int triangle_num;
  int triangle_order = 6;

  printf ( "\n" );
  printf ( "TEST23\n" );
  printf ( "  For an order6 triangulation:\n" );
  printf ( "  TRIANGULATION_ORDER6_BOUNDARY_EDGE_COUNT counts the\n" );
  printf ( "    boundary edges.\n" );

  triangulation_order6_example1_size ( &node_num, &triangle_num, &hole_num );

  node_xy = ( double * ) malloc ( dim_num * node_num * sizeof ( double ) );
  triangle_node = ( int * ) malloc ( triangle_order * triangle_num * sizeof ( int ) );
  triangle_neighbor = ( int * ) malloc ( 3 * triangle_num * sizeof ( int ) );

  triangulation_order6_example1 ( node_num, triangle_num, node_xy, 
    triangle_node, triangle_neighbor );

  boundary_edge_num = triangulation_order6_boundary_edge_count ( triangle_num, 
    triangle_node );

  printf ( "\n" );
  printf ( "  Number of boundary edges = %d\n", boundary_edge_num );
  printf ( "  Correct number =           16\n" );

  free ( node_xy );
  free ( triangle_node );
  free ( triangle_neighbor );

  return;
}
/******************************************************************************/

void test24 ( )

/******************************************************************************/
/*
  Purpose:

    TEST24 tests TRIANGULATION_ORDER6_BOUNDARY_EDGE_COUNT_EULER.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    24 August 2006

  Author:

    John Burkardt
*/
{
  int boundary_num;
  int hole_num;
  int node_num;
  int triangle_num;

  printf ( "\n" );
  printf ( "TEST24\n" );
  printf ( "  For an order6 triangulation:\n" );
  printf ( "  TRIANGULATION_ORDER6_BOUNDARY_EDGE_COUNT_EULER\n" );
  printf ( "  determines the number of edges that lie on the\n" );
  printf ( "  boundary of a region that has been triangulated.\n" );

  triangulation_order6_example1_size ( &node_num, &triangle_num, &hole_num );

  printf ( "\n" );
  printf ( "  Number of nodes =          %d\n", node_num );
  printf ( "  Number of triangles =      %d\n", triangle_num );
  printf ( "  Number of holes =          %d\n", hole_num );

  boundary_num = triangulation_order6_boundary_edge_count_euler ( node_num, 
    triangle_num, hole_num );

  printf ( "  Number of boundary edges = %d\n", boundary_num );
  printf ( "  Correct number =           16\n" );

  return;
}
/******************************************************************************/

void test25 ( )

/******************************************************************************/
/*
  Purpose:

    TEST25 tests TRIANGULATION_ORDER6_BOUNDARY_NODE.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    24 August 2006

  Author:

    John Burkardt
*/
{
  char file_name[80] = "triangulation_order6_plot.eps";
  int i;
  int dim_num = 2;
  int hole_num;
  int *node_boundary;
  int node_num;
  int node_show = 2;
  double *node_xy;
  int triangle_num;
  int *triangle_node;
  int *triangle_neighbor;
  int triangle_order = 6;
  int triangle_show = 2;

  printf ( "\n" );
  printf ( "TEST25\n" );
  printf ( "  For an order6 triangulation:\n" );
  printf ( "  TRIANGULATION_ORDER6_BOUNDARY_COUNT counts the boundary\n" );
  printf ( "    edges.\n" );
  printf ( "  TRIANGULATION_ORDER6_PLOT plots the triangulation.\n" );

  triangulation_order6_example1_size ( &node_num, &triangle_num, &hole_num );

  node_xy = ( double * ) malloc ( dim_num * node_num * sizeof ( double ) );
  triangle_node = ( int * ) malloc ( triangle_order * triangle_num * sizeof ( int ) );
  triangle_neighbor = ( int * ) malloc ( 3 * triangle_num * sizeof ( int ) );

  triangulation_order6_example1 ( node_num, triangle_num, node_xy, 
    triangle_node, triangle_neighbor );
/*
  Make the plot.
*/
  triangulation_order6_plot ( file_name, node_num, node_xy, triangle_num, 
    triangle_node, node_show, triangle_show );

  printf ( "\n" );
  printf ( "  An Encapsulated PostScript image of this\n" );
  printf ( "  triangulation is in \"%s\"\n", file_name );

  node_boundary = triangulation_order6_boundary_node ( node_num, triangle_num, 
    triangle_node );

  printf ( "\n" );
  printf ( "    Node  BN?\n" );
  printf ( "\n" );

  for ( i = 1; i <= node_num; i++ )
  {
    printf ( "  %6d  %d\n", i, node_boundary[i-1] );
  }

  free ( node_boundary );
  free ( node_xy );
  free ( triangle_node );
  free ( triangle_neighbor );

  return;
}
/******************************************************************************/

void test26 ( )

/******************************************************************************/
/*
  Purpose:

    TEST26 tests TRIANGULATION_ORDER6_PRINT.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    24 August 2006

  Author:

    John Burkardt
*/
{
  int hole_num;
  int node_num;
  double *node_xy;
  int *triangle_neighbor;
  int *triangle_node;
  int triangle_num;
  int triangle_order = 6;

  printf ( "\n" );
  printf ( "TEST26\n" );
  printf ( "  For an order6 triangulation:\n" );
  printf ( "  TRIANGULATION_ORDER6_PRINT prints the data.\n" );

  triangulation_order6_example1_size ( &node_num, &triangle_num, &hole_num );

  node_xy = ( double * ) malloc ( 2 * node_num * sizeof ( double ) );
  triangle_neighbor = ( int * ) malloc ( 3 * triangle_num * sizeof ( int ) );
  triangle_node = ( int * ) malloc ( triangle_order * triangle_num * sizeof ( int ) );

  triangulation_order6_example1 ( node_num, triangle_num, node_xy, 
    triangle_node, triangle_neighbor );

  triangulation_order6_print ( node_num, triangle_num, node_xy, 
    triangle_node, triangle_neighbor );

  free ( node_xy );
  free ( triangle_neighbor );
  free ( triangle_node );

  return;
}
/******************************************************************************/

void test265 ( )

/******************************************************************************/
/*
  Purpose:

    TEST265 tests TRIANGULATION_ORDER6_REFINE_COMPUTE.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    10 February 2007

  Author:

    John Burkardt
*/
{
# define DIM_NUM 2
# define NODE_NUM1 12
# define TRIANGLE_NUM1 3
# define TRIANGLE_ORDER 6

  int *edge_data;
  int node_num2;
  double node_xy1[DIM_NUM*NODE_NUM1] = {
       0.0, 0.0, 
       2.0, 0.0, 
       0.0, 2.0, 
       2.0, 2.0, 
       1.0, 3.0, 
       1.0, 0.0, 
       0.0, 1.0, 
       1.0, 1.0, 
       2.0, 1.0, 
       1.0, 2.0, 
       0.5, 2.5, 
       1.5, 2.5 };
  double *node_xy2;
  int triangle_node1[TRIANGLE_ORDER*TRIANGLE_NUM1] = {
       1,  2,  3,  6,  8,  7, 
       4,  3,  2,  9, 10,  8, 
       3,  4,  5, 10, 12, 11 };
  int *triangle_node2;
  int triangle_num2;

  printf ( "\n" );
  printf ( "TEST265\n" );
  printf ( "  For an order6 triangulation:\n" );
  printf ( "  TRIANGULATION_ORDER6_REFINE_SIZE determines the\n" );
  printf ( "  size of a refined triangulation.\n" );
  printf ( "  TRIANGULATION_ORDER6_REFINE_COMPUTES computes the\n" );
  printf ( "  refined triangulation.\n" );

  printf ( "\n" );
  printf ( "  The number of nodes is %d\n", NODE_NUM1 );
  printf ( "  The number of triangles is %d\n", TRIANGLE_NUM1 );

  r8mat_transpose_print ( DIM_NUM, NODE_NUM1, node_xy1, 
    "  The nodes" );

  i4mat_transpose_print ( TRIANGLE_ORDER, TRIANGLE_NUM1, triangle_node1, 
    "  The triangles:" );

  edge_data = ( int * ) malloc ( 5 *  3 * TRIANGLE_NUM1 * sizeof ( int ) );

  printf ( "\n" );
  printf ( "  Sizing the refined mesh:\n" );

  triangulation_order6_refine_size ( NODE_NUM1, TRIANGLE_NUM1, 
    triangle_node1, &node_num2, &triangle_num2, edge_data );

  printf ( "\n" );
  printf ( "  Information about the refined mesh:\n" );
  printf ( "\n" );
  printf ( "  The number of nodes is %d\n", node_num2 );
  printf ( "  The number of triangles is %d\n", triangle_num2 );

  printf ( "\n" );
  printf ( "  Computing the refined mesh:\n" );

  node_xy2 = ( double * ) malloc ( DIM_NUM * node_num2 * sizeof ( double ) );
  triangle_node2 = ( int * ) malloc ( TRIANGLE_ORDER * triangle_num2 * sizeof ( int ) );

  triangulation_order6_refine_compute ( NODE_NUM1, TRIANGLE_NUM1, 
    node_xy1, triangle_node1, node_num2, triangle_num2, edge_data, node_xy2, 
    triangle_node2 );

  r8mat_transpose_print ( DIM_NUM, node_num2, node_xy2, 
    "  The refined nodes" );

  i4mat_transpose_print ( TRIANGLE_ORDER, triangle_num2, triangle_node2, 
    "  The refined triangles:" );

  free ( edge_data );
  free ( node_xy2 );
  free ( triangle_node2 );

  return;
# undef DIM_NUM
# undef NODE_NUM1
# undef TRIANGLE_NUM1
# undef TRIANGLE_ORDER
}
/******************************************************************************/

void test27 ( )

/******************************************************************************/
/*
  Purpose:

    TEST27 tests TRIANGULATION_ORDER6_VERTEX_COUNT.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    24 August 2006

  Author:

    John Burkardt
*/
{
  int hole_num;
  int midside_num;
  int node_num;
  double *node_xy;
  int *triangle_neighbor;
  int *triangle_node;
  int triangle_num;
  int triangle_order = 6;
  int vertex_num;

  printf ( "\n" );
  printf ( "TEST27\n" );
  printf ( "  For an order6 triangulation:\n" );
  printf ( "  TRIANGULATION_ORDER6_VERTEX_COUNT counts the \n" );
  printf ( "  vertex nodes and midside nodes.\n" );

  triangulation_order6_example1_size ( &node_num, &triangle_num, &hole_num );

  node_xy = ( double * ) malloc ( 2 * node_num * sizeof ( double ) );
  triangle_neighbor = ( int * ) malloc ( 3 * triangle_num * sizeof ( int ) );
  triangle_node = ( int * ) malloc ( triangle_order * triangle_num * sizeof ( int ) );

  triangulation_order6_example1 ( node_num, triangle_num, node_xy, 
    triangle_node, triangle_neighbor );

  vertex_num = triangulation_order6_vertex_count ( triangle_num, 
    triangle_node );

  midside_num = node_num - vertex_num;

  printf ( "\n" );
  printf ( "  Number of nodes =         %d\n", node_num );
  printf ( "  Number of vertex nodes =  %d\n", vertex_num );
  printf ( "  Number of midside nodes = %d\n", midside_num );

  free ( node_xy );
  free ( triangle_neighbor );
  free ( triangle_node );

  return;
}
/******************************************************************************/

void test31 ( )

/******************************************************************************/
/*
  Purpose:

    TEST31 tests VORONOI_POLYGON_AREA.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    24 August 2006

  Author:

    John Burkardt
*/
{
# define DIM_NUM 2
# define NEIGHBOR_NUM 4
# define NODE_NUM 5

  double area;
  double area_correct = 0.5;
  int neighbor_index[NEIGHBOR_NUM] = { 0, 1, 2, 3 };
  int node = 4;
  double node_xy[DIM_NUM*NODE_NUM] = { 
    0.0, 0.0,
    1.0, 0.0,
    1.0, 1.0,
    0.0, 1.0,
    0.5, 0.5 };

  printf ( "\n" );
  printf ( "TEST31\n" );
  printf ( "  VORONOI_POLYGON_AREA computes the area of\n" );
  printf ( "  a finite Voronoi polygon.\n" );

  area = voronoi_polygon_area ( node, NEIGHBOR_NUM, neighbor_index,
    NODE_NUM, node_xy );

  printf ( "\n" );
  printf ( "  The computed area is %g\n", area );
  printf ( "  The correct area is  %g\n", area_correct );

  return;
# undef DIM_NUM
# undef NEIGHBOR_NUM
# undef NODE_NUM
}
/******************************************************************************/

void test32 ( )

/******************************************************************************/
/*
  Purpose:

    TEST32 tests VORONOI_POLYGON_CENTROID.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    24 August 2006

  Author:

    John Burkardt
*/
{
# define DIM_NUM 2
# define NEIGHBOR_NUM 4
# define NODE_NUM 5

  double *centroid;
  double centroid_exact[2] = { 0.5, 0.5 };
  int neighbor_index[NEIGHBOR_NUM] = { 0, 1, 2, 3 };
  int node = 4;
  double node_xy[DIM_NUM*NODE_NUM] = { 
    0.0, 0.0,
    1.0, 0.0,
    1.0, 1.0,
    0.0, 1.0,
    0.5, 0.5 };

  printf ( "\n" );
  printf ( "TEST32\n" );
  printf ( "  VORONOI_POLYGON_CENTROID computes the centroid of\n" );
  printf ( "  a finite Voronoi polygon.\n" );

  centroid = voronoi_polygon_centroid ( node, NEIGHBOR_NUM, 
    neighbor_index, NODE_NUM, node_xy );

  printf ( "\n" );
  printf ( "  The computed centroid is %10g  %10g\n", centroid[0], centroid[1] );
  printf ( "  The correct centroid is  %10g  %10g\n", centroid_exact[0], centroid_exact[1] );

  return;
# undef DIM_NUM
# undef NEIGHBOR_NUM
# undef NODE_NUM
}
/******************************************************************************/

void test33 ( )

/******************************************************************************/
/*
  Purpose:

    TEST33 tests VORONOI_POLYGON_VERTICES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    24 August 2006

  Author:

    John Burkardt
*/
{
# define DIM_NUM 2
# define NEIGHBOR_NUM 4
# define NODE_NUM 5

  int neighbor_index[NEIGHBOR_NUM] = { 0, 1, 2, 3 };
  int node = 4;
  double v[DIM_NUM*NEIGHBOR_NUM];
  double v_y[NEIGHBOR_NUM];
  double node_xy[DIM_NUM*NODE_NUM] = { 
    0.0, 0.0,
    1.0, 0.0,
    1.0, 1.0,
    0.0, 1.0,
    0.5, 0.5 };
 
  printf ( "\n" );
  printf ( "TEST33\n" );
  printf ( "  VORONOI_POLYGON_VERTICES computes the vertices of\n" );
  printf ( "  a finite Voronoi polygon.\n" );

  voronoi_polygon_vertices ( node, NEIGHBOR_NUM, neighbor_index, 
    NODE_NUM, node_xy, v );

  r8mat_transpose_print ( DIM_NUM, NEIGHBOR_NUM, v, "  Vertices:" );

  return;
# undef DIM_NUM
# undef NEIGHBOR_NUM
# undef NODE_NUM
}
