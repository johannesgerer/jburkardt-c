# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <string.h>

# include "gmsh_io.h"

int main ( );
void test01 ( );
void test02 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for GMSH_IO_PRB.

  Discussion:

    GMSH_IO_PRB tests the GMSH_IO library.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    18 October 2014

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "GMSH_IO_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the GMSH_IO library.\n" );

  test01 ( );
  test02 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "GMSH_IO_PRB\n" );
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

    TEST01 gets the example 2D data and writes it to a file.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    18 October 2014

  Author:

    John Burkardt
*/
{
  int *element_node;
  int element_num;
  int element_order;
  char gmsh_filename[] = "example_2d.msh";
  int m;
  int node_num;
  double *node_x;

  printf ( "\n" );
  printf ( "TEST01:\n" );
  printf ( "  Get example 2D data, write to a file.\n" );
/*
  Get sizes.
*/
  gmsh_mesh2d_node_size_example ( &node_num, &m );

  gmsh_mesh2d_element_size_example ( &element_num, &element_order );
/*
  Print the sizes.
*/
  printf ( "\n" );
  printf ( "  Number of nodes = %d\n", node_num );
  printf ( "  Spatial dimension = %d\n", m );
  printf ( "  Number of elements = %d\n", element_num );
  printf ( "  Order of elements = %d\n", element_order );
/*
  Get the data.
*/
  node_x = gmsh_mesh2d_node_data_example ( node_num, m );

  element_node = gmsh_mesh2d_element_data_example ( element_num, element_order );
/*
  Print some of the data.
*/
  r8mat_transpose_print_some ( m, node_num, node_x, 
    1, 1, m, 10, "  Coordinates for first 10 nodes:" );

  i4mat_transpose_print_some ( element_order, element_num, element_node, 
    1, 1, element_order, 10, "  Node connectivity of first 10 elements:" );
/*
  Write the GMSH file.
*/
  gmsh_mesh2d_write ( gmsh_filename, m, node_num, node_x, 
    element_order, element_num, element_node );

  printf ( "\n" );
  printf ( "  Wrote example data to file \"%s\"\n", gmsh_filename );
/*
  Clean up.
*/
  free ( element_node );
  free ( node_x );

  return;
}
/******************************************************************************/

void test02 ( )

/******************************************************************************/
/*
  Purpose:

    TEST02 reads the example data from a file.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    20 October 2014

  Author:

   John Burkardt
*/
{
  int *element_node;
  int element_num;
  int element_order;
  char gmsh_filename[] = "example_2d.msh";
  int m;
  int node_num;
  double *node_x;

  printf ( "\n" );
  printf ( "TEST02:\n" );
  printf ( "  Read data from a file.\n" );
/*
  Get the data size.
*/
  gmsh_size_read ( gmsh_filename, &node_num, &m, &element_num, 
    &element_order );
/*
  Print the sizes.
*/
  printf ( "\n" );
  printf ( "  Node data read from file \"%s\"\n", gmsh_filename );
  printf ( "\n" );
  printf ( "  Number of nodes = %d\n", node_num );
  printf ( "  Spatial dimension = %d\n", m );
  printf ( "  Number of elements = %d\n", element_num );
  printf ( "  Element order = %d\n", element_order );
/*
  Allocate memory.
*/
  node_x = ( double * ) malloc ( m * node_num * sizeof ( double ) );
  element_node = ( int * ) 
    malloc ( element_order * element_num * sizeof ( int ) );
/*
  Get the data.
*/
  gmsh_data_read ( gmsh_filename, m, node_num, node_x, 
    element_order, element_num, element_node );
/*
  Print some of the data.
*/
  r8mat_transpose_print_some ( m, node_num, node_x, 
    1, 1, m, 10, "  Coordinates for first 10 nodes:" );

  i4mat_transpose_print_some ( element_order, element_num, element_node, 
    1, 1, element_order, 10, "  Connectivity for first 10 elements:" );
/*
  Clean up.
*/
  free ( element_node );
  free ( node_x );

  return;
}
