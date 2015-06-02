# include <stdlib.h>
# include <stdio.h>
# include <time.h>

# include "triangle_io.h"

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );
void test04 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for TRIANGLE_IO_PRB.

  Discussion:

    TRIANGLE_IO_PRB tests the TRIANGLE_IO library.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    09 November 2010

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "TRIANGLE_IO_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the TRIANGLE_IO library.\n" );

  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "TRIANGLE_IO_PRB\n" );
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

    TEST01 gets the example node data and writes it to a file.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    09 November 2010

  Author:

    John Burkardt
*/
{
  double *node_att;
  int node_att_num;
  double *node_coord;
  int node_dim;
  char *node_file = "example.node";
  int *node_marker;
  int node_marker_num;
  int node_num;

  printf ( "\n" );
  printf ( "TEST01:\n" );
  printf ( "  Get example node data, write to a triangle node file.\n" );
/*
  Get node example size.
*/
  triangle_node_size_example ( &node_num, &node_dim, &node_att_num, 
    &node_marker_num );
/*
  Print the sizes.
*/
  printf ( "\n" );
  printf ( "  Number of nodes = %d\n", node_num );
  printf ( "  Spatial dimension = %d\n", node_dim );
  printf ( "  Number of node attributes = %d\n", node_att_num );
  printf ( "  Number of node markers = %d\n", node_marker_num );
/*
  Allocate memory for node data.
*/
  node_coord = ( double * ) malloc ( node_dim * node_num * sizeof ( double ) );
  node_att = ( double * ) malloc ( node_att_num * node_num * sizeof ( double ) );
  node_marker = ( int * ) malloc ( node_marker_num * node_num * sizeof ( int ) );
/*
  Get the node data.
*/
  triangle_node_data_example ( node_num, node_dim, node_att_num, 
    node_marker_num, node_coord, node_att, node_marker );
/*
  Print some of the data.
*/
  r8mat_transpose_print_some ( node_dim, node_num, node_coord, 
    1, 1, node_dim, 10, "  Coordinates for first 10 nodes:" );

  r8mat_transpose_print_some ( node_att_num, node_num, node_att,
    1, 1, node_att_num, 10, "  Attributes for first 10 nodes:" );

  i4mat_transpose_print_some ( node_marker_num, node_num, node_marker,
    1, 1, node_marker_num, 10, "  Markers for first 10 nodes:" ); 
/*
  Write the node information to node file.
*/
  triangle_node_write ( node_file, node_num, node_dim, node_att_num, 
    node_marker_num, node_coord, node_att, node_marker );

  printf ( "\n" );
  printf ( "  Node data written to file \"%s\"\n", node_file );
/*
  Clean up.
*/
  free ( node_att );
  free ( node_coord );
  free ( node_marker );

  return;
}
/******************************************************************************/

void test02 ( )

/******************************************************************************/
/*
  Purpose:

    TEST02 gets the example element data and writes it to a file.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    09 November 2010

  Author:

    John Burkardt
*/
{
  double *element_att;
  int element_att_num;
  char *element_file = "example.ele";
  int *element_node;
  int element_num;
  int element_order;

  printf ( "\n" );
  printf ( "TEST02:\n" );
  printf ( "  Get example element data, write to a triangle element file.\n" );
/*
  Get element example size.
*/
  triangle_element_size_example ( &element_num, &element_order, 
    &element_att_num );
/*
  Print the sizes.
*/
  printf ( "\n" );
  printf ( "  Number of elements = %d\n", element_num );
  printf ( "  Order of elements = %d\n", element_order );
  printf ( "  Number of element attributes = %d\n", element_att_num );
/*
  Allocate memory.
*/
  element_node = ( int * ) malloc ( element_order * element_num * sizeof ( int ) );
  element_att = ( double * ) 
    malloc ( element_att_num * element_num * sizeof ( double ) );
/*
  Get the data.
*/
  triangle_element_data_example ( element_num, element_order, element_att_num, 
    element_node, element_att );
/*
  Print some of the data.
*/
  i4mat_transpose_print_some ( element_order, element_num, element_node,
    1, 1, element_order, 10, "  Node connectivity of first 10 elements:" );

  r8mat_transpose_print_some ( element_att_num, element_num, element_att,
    1, 1, element_att_num, 10, "  Attributes for first 10 elements:" ); 
/*
  Write the node information to node file.
*/
  triangle_element_write ( element_file, element_num, element_order, 
    element_att_num, element_node, element_att );

  printf ( "\n" );
  printf ( "  Element data written to file \"%s\"\n", element_file );
/*
  Clean up.
*/
  free ( element_att );
  free ( element_node );

  return;
}
/******************************************************************************/

void test03 ( )

/******************************************************************************/
/*
  Purpose:

    TEST03 reads the example node data from a file.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    09 November 2010

  Author:

    John Burkardt
*/
{
  double *node_att;
  int node_att_num;
  double *node_coord;
  int node_dim;
  char *node_file = "example.node";
  int *node_marker;
  int node_marker_num;
  int node_num;

  printf ( "\n" );
  printf ( "TEST03:\n" );
  printf ( "  Read node data from a node file.\n" );
/*
  Get the data size.
*/
  triangle_node_size_read ( node_file, &node_num, &node_dim, &node_att_num, 
    &node_marker_num );
/*
  Print the sizes.
*/
  printf ( "\n" );
  printf ( "  Node data read from file \"%s\"\n", node_file );
  printf ( "\n" );
  printf ( "  Number of nodes = %d\n", node_num );
  printf ( "  Spatial dimension = %d\n", node_dim );
  printf ( "  Number of node attributes = %d\n", node_att_num );
  printf ( "  Number of node markers = %d\n", node_marker_num );
/*
  Allocate memory.
*/
  node_coord = ( double * ) malloc ( node_dim * node_num * sizeof ( double ) );
  node_att = ( double * ) malloc ( node_att_num * node_num * sizeof ( double ) );
  node_marker = ( int * ) malloc ( node_marker_num * node_num * sizeof ( int ) );
/*
  Get the data.
*/
  triangle_node_data_read ( node_file, node_num, node_dim, node_att_num, 
    node_marker_num, node_coord, node_att, node_marker );
/*
  Print some of the data.
*/
  r8mat_transpose_print_some ( node_dim, node_num, node_coord, 
    1, 1, node_dim, 10, "  Coordinates for first 10 nodes:" );

  r8mat_transpose_print_some ( node_att_num, node_num, node_att,
    1, 1, node_att_num, 10, "  Attributes for first 10 nodes:" );

  i4mat_transpose_print_some ( node_marker_num, node_num, node_marker,
    1, 1, node_marker_num, 10, "  Markers for first 10 nodes:" ); 
/*
  Clean up.
*/
  free ( node_att );
  free ( node_coord );
  free ( node_marker );

  return;
}
/******************************************************************************/

void test04 ( )

/******************************************************************************/
/*
  Purpose:

    TEST04 reads the example element data from a file.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    09 November 2010

  Author:

    John Burkardt
*/
{
  double *element_att;
  int element_att_num;
  char *element_file = "example.ele";
  int *element_node;
  int element_num;
  int element_order;

  printf ( "\n" );
  printf ( "TEST04:\n" );
  printf ( "  Read element data from an element file.\n" );
/*
  Get data size.
*/
  triangle_element_size_read ( element_file, &element_num, &element_order, 
    &element_att_num );
/*
  Print the sizes.
*/
  printf ( "\n" );
  printf ( "  Element data read from file \"%s\"\n", element_file );
  printf ( "\n" );
  printf ( "  Number of elements = %d\n", element_num );
  printf ( "  Element order = %d\n", element_order );
  printf ( "  Number of element attributes = %d\n", element_att_num );
/*
  Allocate memory.
*/
  element_node = ( int * ) malloc ( element_order * element_num * sizeof ( int ) );
  element_att = ( double * ) 
    malloc ( element_att_num * element_num * sizeof ( double ) );  
/*
  Get the data.
*/
  triangle_element_data_read ( element_file, element_num, element_order, 
    element_att_num, element_node, element_att );
/*
  Print some of the data.
*/
  i4mat_transpose_print_some ( element_order, element_num, element_node, 
    1, 1, element_order, 10, "  Connectivity for first 10 elements:" );

  r8mat_transpose_print_some ( element_att_num, element_num, element_att,
    1, 1, element_att_num, 10, "  Attributes for first 10 elements:" );
/*
  Clean up.
*/
  free ( element_att );
  free ( element_node );

  return;
}
