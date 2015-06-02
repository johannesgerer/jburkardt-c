# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "fem_io.h"

int main ( );
void test01 ( );
void test02 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for FEM_IO_PRB.

  Discussion:

    FEM_IO_PRB tests the FEM_IO library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    22 March 2006

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "FEM_IO_PRB:\n" );
  printf ( "  C version\n" );
  printf ( "  Test the FEM_IO library.\n" );

  test01 ( );
  test02 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "FEM_IO_PRB:\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/*****************************************************************************80*/

void test01 ( )

/*****************************************************************************80

  Purpose:

    TEST01 tests FEM_READ.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    22 March 2006

  Author:

    John Burkardt
*/
{
  int dim_num;
  char element_file_name[] = "ell_elements.txt";
  int *element_node;
  int element_num;
  int element_order;
  double *node_coord;
  char node_coord_file_name[] = "ell_nodes.txt";
  double *node_data;
  char node_data_file_name[] = "ell_values.txt";
  int node_data_num;
  int node_num;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  FEM_READ reads finite element data from files.\n" );

  printf ( "\n" );
  printf ( "  The node coordinate file name is \"%s\"\n", node_coord_file_name );
  printf ( "  The element file name is \"%s\"\n", element_file_name );
  printf ( "  The node data file name is\"%s\"\n", node_data_file_name );

  fem_header_read ( node_coord_file_name, element_file_name,
    node_data_file_name, &dim_num, &node_num, &element_num,
    &element_order, &node_data_num );

  fem_header_print ( dim_num, node_num, element_order, element_num,
    node_data_num );

  fem_data_read ( node_coord_file_name, element_file_name,
    node_data_file_name, dim_num, node_num, element_num,
    element_order, node_data_num, &node_coord, &element_node, &node_data );

  r8mat_transpose_print_some ( dim_num, node_num, node_coord, 1, 1,
    dim_num,  10, "  First 10 node coordindates:" );

  i4mat_transpose_print_some ( element_order, element_num,
    element_node, 1, 1, element_order, 10, "  First 10 elements" );

  r8mat_transpose_print_some ( node_data_num, node_num, node_data,
    1, 1, node_data_num, 10, "  First 10 node data sets:" );

  free ( element_node );
  free ( node_coord );
  free ( node_data );

  return;
}
/*****************************************************************************80*/

void test02 ( )

/*****************************************************************************80

  Purpose:

    FEM_IO_TEST02 tests FEM_WRITE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    22 March 2006

  Author:

    John Burkardt
*/
{
# define DIM_NUM 2
# define NODE_NUM 5
# define ELEMENT_NUM 3
# define ELEMENT_ORDER 3
# define NODE_DATA_NUM 2

  char element_file_name[] = "tiny_elements.txt";
  int element_node[ELEMENT_ORDER*ELEMENT_NUM] = {
    1, 2, 4,
    5, 4, 2,
    2, 3, 5 };
  double node_coord[DIM_NUM*NODE_NUM] = {
    0.0, 0.0,
    1.0, 0.0,
    2.0, 0.0,
    0.0, 1.0,
    1.0, 1.0 };
  char node_coord_file_name[] = "tiny_nodes.txt";
  double node_data[NODE_DATA_NUM*NODE_NUM] = {
    1.0, 0.0,
    0.8, 0.2,
    0.6, 0.4,
    0.9, 0.1,
    0.5, 0.5 };
  char node_data_file_name[] = "tiny_values.txt";

  printf ( "\n" );
  printf ( "FEM_TEST02\n" );
  printf ( "  Demonstrate the use of FEM_WRITE to write finite\n" );
  printf ( "  element data to files.\n" );

  printf ( "\n" );
  printf ( "  The node coordinate file name is \"%s\"\n", node_coord_file_name );
  printf ( "  The element file name is \"%s\"\n", element_file_name );
  printf ( "  The node data file name is\"%s\"\n", node_data_file_name );

  fem_header_print ( DIM_NUM, NODE_NUM, ELEMENT_ORDER, ELEMENT_NUM,
    NODE_DATA_NUM );

  r8mat_transpose_print ( DIM_NUM, NODE_NUM, node_coord,
    "  Node coordindates:" );

  i4mat_transpose_print ( ELEMENT_ORDER, ELEMENT_NUM,
    element_node, "  Elements" );

  r8mat_transpose_print ( NODE_DATA_NUM, NODE_NUM, node_data,
    "  Node data sets:" );

  fem_write ( node_coord_file_name, element_file_name,
    node_data_file_name, DIM_NUM, NODE_NUM, ELEMENT_NUM,
    ELEMENT_ORDER, NODE_DATA_NUM, node_coord, element_node, node_data );

  return;
# undef DIM_NUM
# undef NODE_NUM
# undef ELEMENT_NUM
# undef ELEMENT_ORDER
# undef NODE_DATA_NUM
}
