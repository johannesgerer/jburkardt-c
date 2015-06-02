# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <string.h>

# include "medit_io.h"

int main ( );
void test01 ( );
void test02 ( );
void test03 ( char *filename );
void test04 ( char *filename );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for MEDIT_IO_PRB.

  Discussion:

    MEDIT_IO_PRB tests the MEDIT_IO library.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    22 November 2010

  Author:

    John Burkardt
*/
{
  char filename[255];

  timestamp ( );
  printf ( "\n" );
  printf ( "MEDIT_IO_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the MEDIT_IO library.\n" );
/*
  Create the file hexahexa_2x2x2.mesh
*/
  test01 ( );
/*
  Read and print the sizes of file hexahexa_2x2x2.mesh.
*/
  strcpy ( filename, "hexahexa_2x2x2.mesh" );
  test03 ( filename );
/*
  Create the file cyl248.mesh
*/
  test02 ( );
/*
  Read and print the sizes of file cyl248.mesh.
*/
  strcpy ( filename, "cyl248.mesh" );
  test03 ( filename );
/*
  Read and print the data in file cyl248.mesh.
*/
  test04 ( filename );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "MEDIT_IO_PRB\n" );
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

    TEST01 creates a MESH dataset and writes it to a file.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    22 November 2010

  Author:

    John Burkardt
*/
{
  int dim;
  int *edge_label;
  int *edge_vertex;
  int edges;
  char filename[255];
  int *hexahedron_label;
  int *hexahedron_vertex;
  int hexahedrons;
  int *quadrilateral_label;
  int *quadrilateral_vertex;
  int quadrilaterals;
  int *tetrahedron_label;
  int *tetrahedron_vertex;
  int tetrahedrons;
  int *triangle_label;
  int *triangle_vertex;
  int triangles;
  double *vertex_coordinate;
  int *vertex_label;
  int vertices;

  printf ( "\n" );
  printf ( "TEST01:\n" );
  printf ( "  Create a hexahedral mesh and write it to a file.\n" );
/*
  Get sizes.
*/
  hexahexa_2x2x2_size ( &dim, &vertices, &edges, &triangles, &quadrilaterals, 
    &tetrahedrons, &hexahedrons );
/*
  Allocate memory.
*/
  edge_label = ( int * ) malloc ( edges * sizeof ( int ) );
  edge_vertex = ( int * ) malloc ( 2 * edges * sizeof ( int ) );
  hexahedron_label = ( int * ) malloc ( hexahedrons * sizeof ( int ) );
  hexahedron_vertex = ( int * ) malloc ( 8 * hexahedrons * sizeof ( int ) );
  quadrilateral_label = ( int * ) malloc ( quadrilaterals * sizeof ( int ) );
  quadrilateral_vertex = ( int * ) malloc ( 4 * quadrilaterals * sizeof ( int ) );
  tetrahedron_label = ( int * ) malloc ( tetrahedrons * sizeof ( int ) );
  tetrahedron_vertex = ( int * ) malloc ( 4 * tetrahedrons * sizeof ( int ) );
  triangle_label = ( int * ) malloc ( triangles * sizeof ( int ) );
  triangle_vertex = ( int * ) malloc ( 3 * triangles * sizeof ( int ) );
  vertex_coordinate = ( double * ) malloc ( dim * vertices * sizeof ( double ) );
  vertex_label = ( int * ) malloc ( vertices * sizeof ( int ) );
/*
  Get the data.
*/
  hexahexa_2x2x2_data ( dim, vertices, edges, triangles, 
    quadrilaterals, tetrahedrons, hexahedrons, vertex_coordinate, vertex_label, 
    edge_vertex, edge_label, triangle_vertex, triangle_label, 
    quadrilateral_vertex, quadrilateral_label, tetrahedron_vertex, 
    tetrahedron_label, hexahedron_vertex, hexahedron_label );
/*
  Write the data.
*/
  strcpy ( filename, "hexahexa_2x2x2.mesh" );

  mesh_write ( filename, dim, vertices, edges, triangles, 
    quadrilaterals, tetrahedrons, hexahedrons, vertex_coordinate, 
    vertex_label, edge_vertex, edge_label, triangle_vertex, triangle_label, 
    quadrilateral_vertex, quadrilateral_label, tetrahedron_vertex, 
    tetrahedron_label, hexahedron_vertex, hexahedron_label );

  printf ( "\n" );
  printf ( "  Created the file \"%s\".\n", filename );
/*
  Deallocate memory.
*/
  free ( edge_label );
  free ( edge_vertex );
  free ( hexahedron_label );
  free ( hexahedron_vertex );
  free ( quadrilateral_label );
  free ( quadrilateral_vertex );
  free ( tetrahedron_label );
  free ( tetrahedron_vertex );
  free ( triangle_label );
  free ( triangle_vertex );
  free ( vertex_coordinate );
  free ( vertex_label );

  return;
}
/******************************************************************************/

void test02 ( )

/******************************************************************************/
/*
  Purpose:

    TEST02 creates a MESH dataset and writes it to a file.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    22 November 2010

  Author:

    John Burkardt
*/
{
  int dim;
  int *edge_label;
  int *edge_vertex;
  int edges;
  char filename[255];
  int *hexahedron_label;
  int *hexahedron_vertex;
  int hexahedrons;
  int *quadrilateral_label;
  int *quadrilateral_vertex;
  int quadrilaterals;
  int *tetrahedron_label;
  int *tetrahedron_vertex;
  int tetrahedrons;
  int *triangle_label;
  int *triangle_vertex;
  int triangles;
  double *vertex_coordinate;
  int *vertex_label;
  int vertices;

  printf ( "\n" );
  printf ( "TEST02:\n" );
  printf ( "  Create a tetrahedral mesh and write it to a file.\n" );
/*
  Get sizes.
*/
  cyl248_size ( &dim, &vertices, &edges, &triangles, &quadrilaterals, 
    &tetrahedrons, &hexahedrons );
/*
  Allocate memory.
*/
  edge_label = ( int * ) malloc ( edges * sizeof ( int ) );
  edge_vertex = ( int * ) malloc ( 2 * edges * sizeof ( int ) );
  hexahedron_label = ( int * ) malloc ( hexahedrons * sizeof ( int ) );
  hexahedron_vertex = ( int * ) malloc ( 8 * hexahedrons * sizeof ( int ) );
  quadrilateral_label = ( int * ) malloc ( quadrilaterals * sizeof ( int ) );
  quadrilateral_vertex = ( int * ) malloc ( 4 * quadrilaterals * sizeof ( int ) );
  tetrahedron_label = ( int * ) malloc ( tetrahedrons * sizeof ( int ) );
  tetrahedron_vertex = ( int * ) malloc ( 4 * tetrahedrons * sizeof ( int ) );
  triangle_label = ( int * ) malloc ( triangles * sizeof ( int ) );
  triangle_vertex = ( int * ) malloc ( 3 * triangles * sizeof ( int ) );
  vertex_coordinate = ( double * ) malloc ( dim * vertices * sizeof ( double ) );
  vertex_label = ( int * ) malloc ( vertices * sizeof ( int ) );
/*
  Get the data.
*/
  cyl248_data ( dim, vertices, edges, triangles, 
    quadrilaterals, tetrahedrons, hexahedrons, vertex_coordinate, vertex_label, 
    edge_vertex, edge_label, triangle_vertex, triangle_label, 
    quadrilateral_vertex, quadrilateral_label, tetrahedron_vertex, 
    tetrahedron_label, hexahedron_vertex, hexahedron_label );
/*
  Write the data.
*/
  strcpy ( filename, "cyl248.mesh" );

  mesh_write ( filename, dim, vertices, edges, triangles, 
    quadrilaterals, tetrahedrons, hexahedrons, vertex_coordinate, 
    vertex_label, edge_vertex, edge_label, triangle_vertex, triangle_label, 
    quadrilateral_vertex, quadrilateral_label, tetrahedron_vertex, 
    tetrahedron_label, hexahedron_vertex, hexahedron_label );

  printf ( "\n" );
  printf ( "  Created the file \"%s\".\n", filename );
/*
  Deallocate memory.
*/
  free ( edge_label );
  free ( edge_vertex );
  free ( hexahedron_label );
  free ( hexahedron_vertex );
  free ( quadrilateral_label );
  free ( quadrilateral_vertex );
  free ( tetrahedron_label );
  free ( tetrahedron_vertex );
  free ( triangle_label );
  free ( triangle_vertex );
  free ( vertex_coordinate );
  free ( vertex_label );

  return;
}
/******************************************************************************/

void test03 ( char *filename )

/******************************************************************************/
/*
  Purpose:

    MESH_IO_TEST03 reads and prints the sizes in a MESH dataset.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    22 November 2010

  Author:

    John Burkardt
*/
{
  int dim;
  int edges;
  int hexahedrons;
  int quadrilaterals;
  int tetrahedrons;
  int triangles;
  int vertices;

  printf ( "\n" );
  printf ( "MESH_IO_TEST03\n" );
  printf ( "  Read a mesh file and print its sizes.\n" );
/*
  Read sizes.
*/
  mesh_size_read ( filename, &dim, &vertices, &edges, &triangles, 
    &quadrilaterals, &tetrahedrons, &hexahedrons );
/*
  Print sizes.
*/
  printf ( "\n" );
  printf ( "  Header information for \"%s\"\n", filename );

  mesh_size_print ( dim, vertices, edges, triangles, quadrilaterals, 
    tetrahedrons, hexahedrons );

  return;
}
/******************************************************************************/

void test04 ( char *filename )

/******************************************************************************/
/*
  Purpose:

    MESH_IO_TEST04 reads a MESH dataset and prints its data.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    22 November 2010

  Author:

    John Burkardt
*/
{
  int dim;
  int *edge_label;
  int *edge_vertex;
  int edges;
  int *hexahedron_label;
  int *hexahedron_vertex;
  int hexahedrons;
  int *quadrilateral_label;
  int *quadrilateral_vertex;
  int quadrilaterals;
  int *tetrahedron_label;
  int *tetrahedron_vertex;
  int tetrahedrons;
  int *triangle_label;
  int *triangle_vertex;
  int triangles;
  double *vertex_coordinate;
  int *vertex_label;
  int vertices;

  printf ( "\n" );
  printf ( "MESH_IO_TEST04\n" );
  printf ( "  Read a mesh file and print its data.\n" );
/*
  Read sizes.
*/
  mesh_size_read ( filename, &dim, &vertices, &edges, &triangles, 
    &quadrilaterals, &tetrahedrons, &hexahedrons );
/*
  Allocate memory.
*/
  edge_label = ( int * ) malloc ( edges * sizeof ( int ) );
  edge_vertex = ( int * ) malloc ( 2 * edges * sizeof ( int ) );
  hexahedron_label = ( int * ) malloc ( hexahedrons * sizeof ( int ) );
  hexahedron_vertex = ( int * ) malloc ( 8 * hexahedrons * sizeof ( int ) );
  quadrilateral_label = ( int * ) malloc ( quadrilaterals * sizeof ( int ) );
  quadrilateral_vertex = ( int * ) malloc ( 4 * quadrilaterals * sizeof ( int ) );
  tetrahedron_label = ( int * ) malloc ( tetrahedrons * sizeof ( int ) );
  tetrahedron_vertex = ( int * ) malloc ( 4 * tetrahedrons * sizeof ( int ) );
  triangle_label = ( int * ) malloc ( triangles * sizeof ( int ) );
  triangle_vertex = ( int * ) malloc ( 3 * triangles * sizeof ( int ) );
  vertex_coordinate = ( double * ) malloc ( dim * vertices * sizeof ( double ) );
  vertex_label = ( int * ) malloc ( vertices * sizeof ( int ) );
/*
  Read the data.
*/
  mesh_data_read ( filename, dim, vertices, edges, triangles, 
    quadrilaterals, tetrahedrons, hexahedrons, vertex_coordinate, 
    vertex_label, edge_vertex, edge_label, triangle_vertex, triangle_label, 
    quadrilateral_vertex, quadrilateral_label, tetrahedron_vertex, 
    tetrahedron_label, hexahedron_vertex, hexahedron_label );
/*
  Print the data.
*/
  printf ( "\n" );
  printf ( "  Data for file \"%s\".\n", filename );

  mesh_data_print ( dim, vertices, edges, triangles, quadrilaterals, 
    tetrahedrons, hexahedrons, vertex_coordinate, vertex_label,  edge_vertex, 
    edge_label, triangle_vertex, triangle_label, quadrilateral_vertex, 
    quadrilateral_label, tetrahedron_vertex, tetrahedron_label, 
    hexahedron_vertex, hexahedron_label );
/*
  Deallocate memory.
*/
  free ( edge_label );
  free ( edge_vertex );
  free ( hexahedron_label );
  free ( hexahedron_vertex );
  free ( quadrilateral_label );
  free ( quadrilateral_vertex );
  free ( tetrahedron_label );
  free ( tetrahedron_vertex );
  free ( triangle_label );
  free ( triangle_vertex );
  free ( vertex_coordinate );
  free ( vertex_label );

  return;
}
