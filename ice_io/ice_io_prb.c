# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "ice_io.h"

int main ( void );
void test01 ( 
  void ice_size ( int *dim, int *vertices, int *edges, int *triangles, 
    int *quadrilaterals, int *tetrahedrons, int *hexahedrons ),
  void ice_data ( int dim, int vertices, int edges, int triangles, 
    int quadrilaterals, int tetrahedrons, int hexahedrons, 
    double vertex_coordinate[], int vertex_label[], int edge_vertex[], 
    int edge_label[], int triangle_vertex[], int triangle_label[], 
    int quadrilateral_vertex[], int quadrilateral_label[], int tetrahedron_vertex[], 
    int tetrahedron_label[], int hexahedron_vertex[], int hexahedron_label[] ),
  char* filename );
void test02 ( char* filename );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    ICE_IO_PRB tests the ICE_IO library.

  Discussion:

    We begin by creating a file.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    22 October 2010

  Author:

    John Burkardt

  Reference:

    Russ Rew, Glenn Davis, Steve Emmerson, Harvey Davies, Ed Hartne,
    The NETCDF User"s Guide,
    Unidata Program Center, March 2009.
*/
{
  char filename_hexa[] = "hexahexa_2x2x2.nc";
  char filename_cyl[] = "cyl248.nc";

  timestamp ( );
  printf ( "\n" );
  printf ( "ICE_IO_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the ICE_IO library.\n" );
/*
  Create "hexahexa_2x2x2.nc"
*/
  test01 ( hexahexa_2x2x2_size, hexahexa_2x2x2_data, filename_hexa );
/*
  Read "hexahexa_2x2x2.nc"
*/
  test02 ( filename_hexa );
/*
  Create "cyl248.nc"
*/
  test01 ( cyl248_size, cyl248_data, filename_cyl );
/*
  Read "cyl248.nc"
*/
  test02 ( filename_cyl );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "ICE_IO_PRB\n" );
  printf ( "  Normal end of execution.\n"  );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void test01 ( 
  void ice_size ( int *dim, int *vertices, int *edges, 
    int *triangles, int *quadrilaterals, int *tetrahedrons, int *hexahedrons ),
  void ice_data ( int dim, int vertices, int edges, int triangles, 
    int quadrilaterals, int tetrahedrons, int hexahedrons, 
    double vertex_coordinate[], int vertex_label[], int edge_vertex[], 
    int edge_label[], int triangle_vertex[], int triangle_label[], 
    int quadrilateral_vertex[], int quadrilateral_label[], 
    int tetrahedron_vertex[], int tetrahedron_label[], int hexahedron_vertex[], 
    int hexahedron_label[] ),
  char* filename )

/******************************************************************************/
/*
  Purpose:

    TEST01 writes an ICE grid dataset to a NETCDF file.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    22 October 2010

  Author:

    John Burkardt

  Reference:

    Russ Rew, Glenn Davis, Steve Emmerson, Harvey Davies, Ed Hartne,
    The NETCDF User"s Guide,
    Unidata Program Center, March 2009.

  Parameters:

    Input, void ICE_SIZE(), a function which determines sizes.

    Input, void ICE_DATA(), a function which determines data.

    Input, char *FILENAME, the name of the file to be created.
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
  printf ( "TEST01:\n" );
  printf ( "  Create an ICE grid dataset, print it,\n" );
  printf ( "  and write it to an NETCDF file.\n" );
/*
  Get sizes.
*/
  ice_size ( &dim, &vertices, &edges, &triangles, 
    &quadrilaterals, &tetrahedrons, &hexahedrons );
/*
  Print sizes;
*/
  size_print ( dim, vertices, edges, triangles, quadrilaterals, 
    tetrahedrons, hexahedrons );
/*
  Allocate memory.
*/
  vertex_coordinate = ( double * ) malloc ( 3 * vertices * sizeof ( double ) );
  vertex_label = ( int * ) malloc ( vertices * sizeof ( int ) );
  edge_vertex = ( int * ) malloc ( 2 * edges * sizeof ( int ) );
  edge_label = ( int * ) malloc ( edges * sizeof ( int ) );
  triangle_vertex = ( int * ) malloc ( 3 * triangles * sizeof ( int ) );
  triangle_label = ( int * ) malloc ( triangles * sizeof ( int ) );
  quadrilateral_vertex = ( int * ) malloc ( 4 * quadrilaterals * sizeof ( int ) );
  quadrilateral_label = ( int * ) malloc ( quadrilaterals * sizeof ( int ) );
  tetrahedron_vertex = ( int * ) malloc ( 4 * tetrahedrons * sizeof ( int ) );
  tetrahedron_label = ( int * ) malloc ( tetrahedrons * sizeof ( int ) );
  hexahedron_vertex = ( int * ) malloc ( 8 * hexahedrons * sizeof ( int ) );
  hexahedron_label = ( int * ) malloc ( hexahedrons * sizeof ( int ) );
/*
  Get data.
*/
  ice_data ( dim, vertices, edges, triangles, quadrilaterals, tetrahedrons, 
    hexahedrons, vertex_coordinate, vertex_label, edge_vertex, edge_label, 
    triangle_vertex, triangle_label, quadrilateral_vertex, quadrilateral_label, 
    tetrahedron_vertex, tetrahedron_label, hexahedron_vertex, hexahedron_label );
/*
  Print the data.
*/
  printf ( "\n" );
  printf ( "  Data to be written to \"%s\":\n", filename );

  data_print ( dim, vertices, edges, triangles, quadrilaterals, tetrahedrons,
    hexahedrons, vertex_coordinate, vertex_label, edge_vertex, edge_label, 
    triangle_vertex, triangle_label, quadrilateral_vertex, quadrilateral_label,
    tetrahedron_vertex, tetrahedron_label, hexahedron_vertex, hexahedron_label );
/*
  Create the file.
*/
  ice_write ( filename, dim, vertices, edges, triangles, quadrilaterals, 
    tetrahedrons, hexahedrons, vertex_coordinate, vertex_label, edge_vertex, 
    edge_label, triangle_vertex, triangle_label, quadrilateral_vertex, 
    quadrilateral_label, tetrahedron_vertex, tetrahedron_label, 
    hexahedron_vertex, hexahedron_label );

  printf ( "\n" );
  printf ( "  Created the file \"%s\"\n", filename );
/*
  Free memory.
*/
  free ( vertex_coordinate );
  free ( vertex_label );
  free ( edge_vertex );
  free ( edge_label );
  free ( triangle_vertex );
  free ( triangle_label );
  free ( quadrilateral_vertex );
  free ( quadrilateral_label );
  free ( tetrahedron_vertex );
  free ( tetrahedron_label );
  free ( hexahedron_vertex );
  free ( hexahedron_label );

  return;
}
/******************************************************************************/

void test02 ( char *filename )

/******************************************************************************/
/*
  Purpose:

    TEST02 reads an ICE grid dataset from a NETCDF file.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    22 October 2010

  Author:

    John Burkardt

  Reference:

    Russ Rew, Glenn Davis, Steve Emmerson, Harvey Davies, Ed Hartne,
    The NETCDF User"s Guide,
    Unidata Program Center, March 2009.

  Parameters:

    Input, char *FILENAME, the name of the file to be read.
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
  printf ( "TEST02:\n" );
  printf ( "  Read an ICE grid dataset from a NETCDF file,\n");
  printf ( "  and print the data.\n" );
/*
  Read sizes;
*/
  size_read ( filename, &dim, &vertices, &edges, &triangles, &quadrilaterals, 
    &tetrahedrons, &hexahedrons );
/*
  Print sizes;
*/
  size_print ( dim, vertices, edges, triangles, quadrilaterals, 
    tetrahedrons, hexahedrons );
/*
  Allocate memory.
*/
  vertex_coordinate = ( double * ) malloc ( 3 * vertices * sizeof ( double ) );
  vertex_label = ( int * ) malloc ( vertices * sizeof ( int ) );
  edge_vertex = ( int * ) malloc ( 2 * edges * sizeof ( int ) );
  edge_label = ( int * ) malloc ( edges * sizeof ( int ) );
  triangle_vertex = ( int * ) malloc ( 3 * triangles * sizeof ( int ) );
  triangle_label = ( int * ) malloc ( triangles * sizeof ( int ) );
  quadrilateral_vertex = ( int * ) malloc ( 4 * quadrilaterals * sizeof ( int ) );
  quadrilateral_label = ( int * ) malloc ( quadrilaterals * sizeof ( int ) );
  tetrahedron_vertex = ( int * ) malloc ( 4 * tetrahedrons * sizeof ( int ) );
  tetrahedron_label = ( int * ) malloc ( tetrahedrons * sizeof ( int ) );
  hexahedron_vertex = ( int * ) malloc ( 8 * hexahedrons * sizeof ( int ) );
  hexahedron_label = ( int * ) malloc ( hexahedrons * sizeof ( int ) );
/*
  Read the file
*/
  data_read ( filename, dim, vertices, edges, triangles, quadrilaterals,
    tetrahedrons, hexahedrons, vertex_coordinate, vertex_label, edge_vertex, 
    edge_label, triangle_vertex, triangle_label, quadrilateral_vertex, 
    quadrilateral_label, tetrahedron_vertex, tetrahedron_label, 
    hexahedron_vertex, hexahedron_label );
/*
  Print the data.
*/
  printf ( "\n" );
  printf ( "  Data from file \"%s\"\n", filename );

  data_print ( dim, vertices, edges, triangles, quadrilaterals, tetrahedrons,
    hexahedrons, vertex_coordinate, vertex_label, edge_vertex, edge_label, 
    triangle_vertex, triangle_label, quadrilateral_vertex, quadrilateral_label,
    tetrahedron_vertex, tetrahedron_label, hexahedron_vertex, hexahedron_label );
/*
  Free memory.
*/
  free ( vertex_coordinate );
  free ( vertex_label );
  free ( edge_vertex );
  free ( edge_label );
  free ( triangle_vertex );
  free ( triangle_label );
  free ( quadrilateral_vertex );
  free ( quadrilateral_label );
  free ( tetrahedron_vertex );
  free ( tetrahedron_label );
  free ( hexahedron_vertex );
  free ( hexahedron_label );

  return;
}
