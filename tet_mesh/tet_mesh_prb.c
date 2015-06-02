/*  Output from c_comment */

# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "tet_mesh.h"

int main ( );
void test001 ( );
void test002 ( );
void test003 ( );
void test004 ( );
void test005 ( );
void test006 ( );
void test007 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for TET_MESH_PRB

  Discussion:

    TET_MESH_PRB tests the TET_MESH library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 October 2014

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "TET_MESH_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the TET_MESH library.\n" );

  test001 ( );
  test002 ( );
  test003 ( );
  test004 ( );
  test005 ( );
  test006 ( );
  test007 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "TET_MESH_PRB\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void test001 ( )

/******************************************************************************/
/*
  Purpose:

    TEST001 tests R8MAT_SOLVE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    06 December 2006

  Author:

    John Burkardt
*/
{
# define N 3
# define RHS_NUM 2

  double a[N*(N+RHS_NUM)] = {
     1.0,  4.0,  7.0,
     2.0,  5.0,  8.0,
     3.0,  6.0,  0.0,
    14.0, 32.0, 23.0,
     7.0, 16.0,  7.0 };
  int i;
  int info;
  int j;

  printf ( "\n" );
  printf ( "TEST001\n" );
  printf ( "  R8MAT_SOLVE solves linear systems.\n" );
/*
  Print out the matrix to be inverted.
*/
  r8mat_print ( N, N+RHS_NUM, a, "  The linear system:" );
/*
  Solve the systems.
*/
  info = r8mat_solve ( N, RHS_NUM, a );

  if ( info != 0 )
  {
    printf ( "\n" );
    printf ( "  The input matrix was singular.\n" );
    printf ( "  The solutions could not be computed.\n" );
    return;
  }

  printf ( "\n" );
  printf ( "  The computed solutions:\n" );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    for ( j = N; j < N+RHS_NUM; j++ )
    {
      printf ( "%10g.4  ", a[i+j*N] );
    }
    printf ( "\n" );
  }

  return;
# undef N
# undef RHS_NUM
}
/******************************************************************************/

void test002 ( )

/******************************************************************************/
/*
  Purpose:

    TEST002 tests TETRAHEDRON_ORDER4_PHYSICAL_TO_REFERENCE, and
    TETRAHEDRON_ORDER4_REFERENCE_TO_PHYSICAL.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 December 2006

  Author:

    John Burkardt
*/
{
# define N 10

  int i;
  int j;
  double phy[3*N];
  double ref[3*N];
  double ref2[3*N];
  int seed;
  double t[3*4] = {
    5.0, 0.0, 0.0,
    8.0, 0.0, 0.0,
    5.0, 2.0, 0.0,
    6.0, 1.0, 2.0 };

  seed = 123456789;

  printf ( "\n" );
  printf ( "TEST002\n" );
  printf ( "  For an order 4 tetrahedron,\n" );
  printf ( "  TETRAHEDRON_ORDER4_PHYSICAL_TO_REFERENCE\n" );
  printf ( "  maps a physical point to a reference point.\n" );
  printf ( "  TETRAHEDRON_ORDER4_REFERENCE_TO_PHYSICAL \n" );
  printf ( "  maps a reference point to a physical point.\n" );
  printf ( "\n" );
  printf ( "     ( R, S, T )          ==>  ( X, Y, Z )           ==> ( R2, S2, T2 )\n" );
  printf ( "\n" );

  tetrahedron_reference_sample ( N, &seed, ref );

  tetrahedron_order4_reference_to_physical ( t, N, ref, phy );
  tetrahedron_order4_physical_to_reference ( t, N, phy, ref2 );

  for ( j = 0; j < N; j++ )
  {
    printf ( "  %8g  %8g  %8g    %8g  %8g  %8g    %8g  %8g  %8g\n",
      ref[0+j*3], ref[1+j*3], ref[2+j*3],
      phy[0+j*3], phy[1+j*3], phy[2+j*3],
      ref2[0+j*3], ref2[1+j*3], ref2[2+j*3] );
  }

  return;
# undef N
}
/******************************************************************************/

void test003 ( )

/******************************************************************************/
/*
  Purpose:

    TEST003 tests TETRAHEDRON_ORDER10_TO_ORDER4.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    22 January 2007

  Author:

    John Burkardt
*/
{
  int node_num1;
  int node_num2;
  double *node_xyz;
  int *tet_node1;
  int *tet_node2;
  int tet_num1;
  int tet_num2;
  int tet_order1 = 10;
  int tet_order2 = 4;

  printf ( "\n" );
  printf ( "TEST003\n" );
  printf ( "  For an order 10 tet mesh,\n" );
  printf ( "  TETRAHEDRON_ORDER10_TO_ORDER4\n" );
  printf ( "  makes a linear (order 4) tet mesh by using\n" );
  printf ( "  the existing nodes, and replacing each\n" );
  printf ( "  quadratic tetrahedron by 8 linear tetrahedrons.\n" );

  tet_mesh_order10_example_size ( &node_num1, &tet_num1 );

  node_xyz = ( double * ) malloc ( 3 * node_num1 * sizeof ( double ) );
  tet_node1 = ( int * ) malloc ( tet_order1 * tet_num1 * sizeof ( int ) );

  tet_mesh_order10_example_set ( node_num1, tet_num1,
    node_xyz, tet_node1 );

  i4mat_transpose_print_some ( tet_order1, tet_num1, tet_node1,
    1, 1, tet_order1, 5, "  First 5 quadratic tetrahedrons:" );

  tet_mesh_order10_to_order4_size ( node_num1, tet_num1,
    &node_num2, &tet_num2 );

  printf ( "\n" );
  printf ( "  Quadratic mesh size is       %d\n", tet_num1 );
  printf ( "  Linearized mesh size will be %d\n", tet_num2 );

  tet_node2 = ( int * ) malloc ( tet_order2 * tet_num2 * sizeof ( int ) );

  tet_mesh_order10_to_order4_compute ( tet_num1, tet_node1,
    tet_num2, tet_node2 );

  i4mat_transpose_print_some ( tet_order2, tet_num2, tet_node2,
    1, 1, tet_order2, 5, "  First 5 linear tetrahedrons:" );

  free ( node_xyz );
  free ( tet_node1 );
  free ( tet_node2 );

  return;
}
/******************************************************************************/

void test004 ( )

/******************************************************************************/
/*
  Purpose:

    TEST004 tests TETRAHEDRON_ORDER10_TO_ORDER4.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 July 2009

  Author:

    John Burkardt
*/
{
  int node_num;
  int *node_order;
  double *node_xyz;
  int *tet_node;
  int tet_num;
  int tet_order = 10;

  printf ( "\n" );
  printf ( "TEST004\n" );
  printf ( "  TET_MESH_NODE_ORDER determines the order of \n" );
  printf ( "  each node in a tet mesh.\n" );
  printf ( "\n" );
  printf ( "  The order of a node is the number of tetrahedrons\n" );
  printf ( "  that use the node as part of their definition.\n" );

  tet_mesh_order10_example_size ( &node_num, &tet_num );

  printf ( "\n" );
  printf ( "  This mesh has tetrahedron order %d\n", tet_order );
  printf ( "  The number of tetrahedrons is   %d\n", tet_num );

  node_xyz = ( double * ) malloc ( 3 * node_num * sizeof ( double ) );
  tet_node = ( int * ) malloc ( tet_order * tet_num * sizeof ( int ) );

  tet_mesh_order10_example_set ( node_num, tet_num,
    node_xyz, tet_node );

  i4mat_transpose_print ( tet_order, tet_num, tet_node,
    "  The tet mesh:" );

  node_order = tet_mesh_node_order ( tet_order, tet_num, tet_node, node_num );

  i4vec_print ( node_num, node_order, "  Node orders:" );

  printf ( "\n" );
  printf ( "  Check that the following are equal:\n" );
  printf ( "\n" );
  printf ( "  Number of tetrahedrons * order = %d\n", tet_num * tet_order );
  printf ( "  Sum of node orders             = %d\n", i4vec_sum ( node_num, node_order ) );

  free ( node_order );
  free ( node_xyz );
  free ( tet_node );

  return;
}
/******************************************************************************/

void test005 ( )

/******************************************************************************/
/*
  Purpose:

    TEST005 tests TETRAHEDRON_BARYCENTRIC.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 August 2009

  Author:

    John Burkardt
*/
{
  double *c1;
  double c1_sum;
  double *c2;
  int i;
  double *p;
  int seed;
  int test1;
  int test1_num = 3;
  int test2;
  int test2_num = 5;
  double *tet_xyz;

  seed = 123456789;

  printf ( "\n" );
  printf ( "TEST005\n" );
  printf ( "  TETRAHEDRON_BARYCENTRIC computes the barycentric\n" );
  printf ( "  coordinates of a point.\n" );
/*
  Choose a random tetrahedron.
*/
  for ( test1 = 1; test1 <= test1_num; test1++ )
  {
    tet_xyz = r8mat_uniform_01_new ( 3, 4, &seed );

    r8mat_transpose_print ( 3, 4, tet_xyz, "  Random tetrahedron:" );
/*
  Choose barycentric coordinates C1 at random.

  Define a point P.

  Have TETRAHEDRON_BARYCENTRIC compute C2, the barycentric coordinates of P.
*/
    for ( test2 = 1; test2 <= test2_num; test2++ )
    {
      c1 = r8vec_uniform_01_new ( 4, &seed );
      c1_sum = r8vec_sum ( 4, c1 );
      for ( i = 0; i < 4; i++ )
      {
        c1[i] = c1[i] / c1_sum;
      }

      p = r8mat_mv_new ( 3, 4, tet_xyz, c1 );

      c2 = tetrahedron_barycentric ( tet_xyz, p );

      printf ( "\n" );
      printf ( "  C1 = " );
      for ( i = 0; i < 4; i++ )
      {
        printf ( "  %14.6g", c1[i] );
      }
      printf ( "\n" );
      printf ( "  C2 = " );
      for ( i = 0; i < 4; i++ )
      {
        printf ( "  %14.6g", c2[i] );
      }
      printf ( "\n" );

      free ( c1 );
      free ( c2 );
      free ( p );
    }
    free ( tet_xyz );
  }

  return;
}
/******************************************************************************/

void test006 ( )

/******************************************************************************/
/*
  Purpose:

    TEST006 tests TET_MESH_TET_NEIGHBORS.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    18 August 2009

  Author:

    John Burkardt
*/
{
  int node_num;
  double *node_xyz;
  int *tet_neighbor;
  int *tet_node;
  int tet_num;
  int tet_order = 4;

  printf ( "\n" );
  printf ( "TEST006\n" );
  printf ( "  TET_MESH_TET_NEIGHBORS computes the 4 neighboring\n" );
  printf ( "  tetrahedrons of each tetrahedron in a tet mesh.\n" );
  printf ( "  containing a point.\n" );
/*
  Set up the example tetrahedron mesh.
*/
  tet_mesh_order4_example_size ( &node_num, &tet_num );

  printf ( "\n" );
  printf ( "  This mesh has tetrahedron order %d\n", tet_order );
  printf ( "  The number of tetrahedrons is   %d\n", tet_num );

  node_xyz = ( double * ) malloc ( 3 * node_num * sizeof ( double ) );
  tet_node = ( int * ) malloc ( tet_order * tet_num * sizeof ( int ) );

  tet_mesh_order4_example_set ( node_num, tet_num, node_xyz, tet_node );
/*
  Print the tets.
*/
  i4mat_transpose_print_some ( tet_order, tet_num, tet_node,
    1, 1, tet_order, 10, "  First 10 Tets:" );
/*
  The TET_NEIGHBOR array is needed by TET_MESH_DELAUNAY_SEARCH.
*/
  tet_neighbor = tet_mesh_neighbor_tets ( tet_order, tet_num, tet_node );

  i4mat_transpose_print_some ( 4, tet_num, tet_neighbor,
    1, 1, 4, 10, "  First 10 Tet Neighbors:" );

  free ( node_xyz );
  free ( tet_neighbor );
  free ( tet_node );

  return;
}
/******************************************************************************/

void test007 ( )

/******************************************************************************/
/*
  Purpose:

    TEST007 tests TET_MESH_SEARCH_NAIVE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 August 2009

  Author:

    John Burkardt
*/
{
  int face;
  int i;
  int j;
  int k;
  int node_num;
  double *node_xyz;
  double p[3];
  int seed;
  int step_num;
  int test;
  int test_num = 5;
  int *tet_neighbor;
  int *tet_node;
  int tet_num;
  int tet_order = 4;
  double tet_xyz[3*4];
  int tet1;
  int tet2;
  int tet3;

  seed = 123456789;

  printf ( "\n" );
  printf ( "TEST007\n" );
  printf ( "  TET_MESH_SEARCH_NAIVE uses a naive algorithm\n" );
  printf ( "  to search a tetrahedral mesh for the tetrahedron\n" );
  printf ( "  containing a point.\n" );
/*
  Set up the example tetrahedron mesh.
*/
  tet_mesh_order4_example_size ( &node_num, &tet_num );

  printf ( "\n" );
  printf ( "  This mesh has tetrahedron order %d\n", tet_order );
  printf ( "  The number of tetrahedrons is   %d\n", tet_num );

  node_xyz = ( double * ) malloc ( 3 * node_num * sizeof ( double ) );
  tet_node = ( int * ) malloc ( tet_order*tet_num * sizeof ( int ) );

  tet_mesh_order4_example_set ( node_num, tet_num, node_xyz, tet_node );
/*
  The TET_NEIGHBOR array is needed for the Delaunay search.
*/
  tet_neighbor = tet_mesh_neighbor_tets ( tet_order, tet_num, tet_node );

  for ( test = 1; test <= test_num; test++ )
  {
/*
  Choose a tetrahedron at random.
*/
    tet1 = i4_uniform_ab ( 0, tet_num - 1, &seed );

    printf ( "\n" );
    printf ( "  Point was chosen from tetrahedron    %8d\n", tet1 );

    for ( j = 0; j < 4; j++ )
    {
      k = tet_node[j+tet1*4];
      for ( i = 0; i < 3; i++ )
      {
        tet_xyz[i+j*3] = node_xyz[i+k*3];
      }
    }
/*
  Choose a point in the tetrahedron at random.
*/
    tetrahedron_sample ( tet_xyz, 1, &seed, p );
/*
  Naive search.
*/
    tet2 = tet_mesh_search_naive ( node_num, node_xyz, tet_order, tet_num,
      tet_node, p, &step_num );

    printf ( "  Naive search ended in tetrahedron    %8d, number of steps = %d\n", 
      tet2, step_num );
/*
  Delaunay search.
*/
    tet3 = tet_mesh_search_delaunay ( node_num, node_xyz, tet_order,
      tet_num, tet_node, tet_neighbor, p, &face, &step_num );

    printf ( "  Delaunay search ended in tetrahedron %d, number of steps = %d\n",
      tet3, step_num );
  }

  free ( node_xyz );
  free ( tet_neighbor );
  free ( tet_node );

  return;
}
