# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <string.h>

# include "netcdf.h"

int main ( int argc, char **argv );
void data_print ( int dim, int vertices, int edges, int triangles,
  int quadrilaterals, int tetrahedrons, int hexahedrons,
  double vertex_coordinate[], int vertex_label[], int edge_vertex[],
  int edge_label[], int triangle_vertex[], int triangle_label[],
  int quadrilateral_vertex[], int quadrilateral_label[],
  int tetrahedron_vertex[], int tetrahedron_label[],
  int hexahedron_vertex[], int hexahedron_label[] );
void data_read ( char *filename, int dim, int vertices, int edges,
  int triangles, int quadrilaterals, int tetrahedrons, int hexahedrons,
  double vertex_coordinate[], int vertex_label[], int edge_vertex[],
  int edge_label[], int triangle_vertex[], int triangle_label[],
  int quadrilateral_vertex[], int quadrilateral_label[],
  int tetrahedron_vertex[], int tetrahedron_label[], int hexahedron_vertex[],
  int hexahedron_label[] );
void mesh_write ( char *filename, int dim, int vertices, int edges,
  int triangles, int quadrilaterals, int tetrahedrons, int hexahedrons,
  double vertex_coordinate[], int vertex_label[], int edge_vertex[],
  int edge_label[], int triangle_vertex[], int triangle_label[],
  int quadrilateral_vertex[], int quadrilateral_label[],
  int tetrahedron_vertex[], int tetrahedron_label[],
  int hexahedron_vertex[], int hexahedron_label[] );
void size_print ( int dim, int vertices, int edges, int triangles,
  int quadrilaterals, int tetrahedrons, int hexahedrons );
void size_read ( char *filename, int *dim, int *vertices, int *edges,
  int *triangles, int *quadrilaterals, int *tetrahedrons, int *hexahedrons );
void timestamp ( void );

/******************************************************************************/

int main ( int argc, char **argv )

/******************************************************************************/
/*
  Purpose:

    ICE_TO_MESH reads ICE data from a NETCDF file and writes to a MESH file.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    23 November 2010

  Author:

    John Burkardt

  Parameters:
*/
{
  int dim;
  int *edge_label;
  int *edge_vertex;
  int edges;
  char filename_mesh[255];
  char filename_nc[255];
  int *hexahedron_label;
  int *hexahedron_vertex;
  int hexahedrons;
  char prefix[255];
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

  timestamp ( );
  printf ( "\n" );
  printf ( "ICE_TO_MESH:\n" );
  printf ( "  C version\n" );
  printf ( "  Read ICE data from NETCDF file, write to MESH file.\n" );
/*
  Check the input argument.
*/
  if ( 2 <= argc )
  {
    strcpy ( prefix, argv[1] );
  }
  else
  {
    printf ( "\n" );
    printf ( "  Enter the filename prefix:\n" );
    scanf ( "%s", prefix );
  }
/*
  Create the file names.
*/
  strcpy ( filename_nc, prefix );
  strcat ( filename_nc, ".nc" );
  strcpy ( filename_mesh, prefix );
  strcat ( filename_mesh, ".mesh" );
/*
  Read sizes;
*/
  size_read ( filename_nc, &dim, &vertices, &edges, &triangles,
    &quadrilaterals, &tetrahedrons, &hexahedrons );
/*
  Print sizes.
*/
  size_print ( dim, vertices, edges, triangles, quadrilaterals,
    tetrahedrons, hexahedrons );
/*
  Allocate memory.
*/
  vertex_coordinate = ( double * ) malloc ( dim * vertices * sizeof ( double ) );
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
  Read the data.
*/
  printf ( "\n" );
  printf ( "  Reading \"%s\".\n", filename_nc );

  data_read ( filename_nc, dim, vertices, edges, triangles,
    quadrilaterals, tetrahedrons, hexahedrons, vertex_coordinate,
    vertex_label, edge_vertex, edge_label, triangle_vertex,
    triangle_label, quadrilateral_vertex, quadrilateral_label,
    tetrahedron_vertex, tetrahedron_label, hexahedron_vertex,
    hexahedron_label );
/*
  Print the data.
*/
  if ( vertices < 250 )
  {
    data_print ( dim, vertices, edges, triangles, quadrilaterals,
      tetrahedrons, hexahedrons, vertex_coordinate, vertex_label,
      edge_vertex, edge_label, triangle_vertex, triangle_label,
      quadrilateral_vertex, quadrilateral_label, tetrahedron_vertex,
      tetrahedron_label, hexahedron_vertex, hexahedron_label );
  }
/*
  Write the data.
*/
  printf ( "\n" );
  printf ( "  Writing \"%s\".\n", filename_mesh );

  mesh_write ( filename_mesh, dim, vertices, edges, triangles,
    quadrilaterals, tetrahedrons, hexahedrons, vertex_coordinate,
    vertex_label, edge_vertex, edge_label,  triangle_vertex, triangle_label,
    quadrilateral_vertex, quadrilateral_label, tetrahedron_vertex,
    tetrahedron_label, hexahedron_vertex, hexahedron_label );

  printf ( "  Conversion completed.\n" );
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
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "ICE_TO_MESH:\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void data_print ( int dim, int vertices, int edges, int triangles,
  int quadrilaterals, int tetrahedrons, int hexahedrons,
  double vertex_coordinate[], int vertex_label[], int edge_vertex[],
  int edge_label[], int triangle_vertex[], int triangle_label[],
  int quadrilateral_vertex[], int quadrilateral_label[],
  int tetrahedron_vertex[], int tetrahedron_label[],
  int hexahedron_vertex[], int hexahedron_label[] )

/******************************************************************************/
/*
  Purpose:

    DATA_PRINT prints the data of an ICE grid dataset.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    21 October 2010

  Author:

    John Burkardt

  Reference:

    Pascal Frey,
    MEDIT: An interactive mesh visualization software,
    Technical Report RT-0253,
    Institut National de Recherche en Informatique et en Automatique,
    03 December 2001.

  Parameters:

    Input, int DIM, the spatial dimension, which should be 2 or 3.

    Input, int VERTICES, the number of vertices.

    Input, int EDGES, the number of edges (may be 0).

    Input, int TRIANGLES, the number of triangles (may be 0).

    Input, int QUADRILATERALS, the number of quadrilaterals (may be 0).

    Input, int TETRAHEDRONS, the number of tetrahedrons (may be 0).

    Input, int HEXAHEDRONS, the number of hexahedrons (may be 0).

    Input, double VERTEX_COORDINATE[DIM*VERTICES], the coordinates
    of each vertex.

    Input, int VERTEX_LABEL[VERTICES], a label for each vertex.

    Input, int EDGE_VERTEX[2*EDGES], the vertices that form each edge.

    Input, int EDGE_LABEL[EDGES], a label for each edge.

    Input, int TRIANGLE_VERTEX[3*TRIANGLES], the vertices that form
    each triangle.

    Input, int TRIANGLE_LABEL[TRIANGLES], a label for each triangle.

    Input, int QUADRILATERAL_VERTEX[4*QUADRILATERALS], the vertices that
    form each quadrilateral.

    Input, int QUADRILATERAL_LABEL[QUADRILATERALS], a label for
    each quadrilateral.

    Input, int TETRAHEDRON_VERTEX[4*TETRAHEDRONS], the vertices that
    form each tetrahedron.

    Input, int TETRAHEDRON_LABEL[TETRAHEDRONS], a label for
    each tetrahedron.

    Input, int HEXAHEDRON_VERTEX[8*HEXAHEDRONS], the vertices that form
    each hexahedron.

    Input, int HEXAHEDRON_LABEL[HEXAHEDRONS], a label for each hexahedron.
*/
{
  int i;
  int j;

  printf ( "\n" );
  printf ( "  Vertices:\n" );
  printf ( "\n" );
  for ( j = 0; j < vertices; j++ )
  {
    for ( i = 0; i < dim; i++ )
    {
      printf ( "  %f", vertex_coordinate[i+j*dim] );
    }
    printf ( "  (%d)\n", vertex_label[j] );
  }

  if ( 0 < edges )
  {
    printf ( "\n" );
    printf ( "  Edges:\n" );
    printf ( "\n" );
    for ( j = 0; j < edges; j++ )
    {
      for ( i = 0; i < 2; i++ )
      {
        printf ( "  %d", edge_vertex[i+j*2] );
    }
    printf ( "  (%d)\n", edge_label[j] );
    }
  }

  if ( 0 < triangles )
  {
    printf ( "\n" );
    printf ( "  Triangles:\n" );
    printf ( "\n" );
    for ( j = 0; j < triangles; j++ )
    {
      for ( i = 0; i < 3; i++ )
      {
        printf ( "  %d", triangle_vertex[i+j*3] );
      }
      printf ( "  (%d)\n", triangle_label[j] );
    }
  }

  if ( 0 < quadrilaterals )
  {
    printf ( "\n" );
    printf ( "  Quadrilaterals:\n" );
    printf ( "\n" );
    for ( j = 0; j < quadrilaterals; j++ )
    {
      for ( i = 0; i < 4; i++ )
      {
        printf ( "  %d", quadrilateral_vertex[i+j*4] );
      }
      printf ( "  (%d)\n", quadrilateral_label[j] );
    }
  }

  if ( 0 < tetrahedrons )
  {
    printf ( "\n" );
    printf ( "  Tetrahedrons:\n" );
    printf ( "\n" );
    for ( j = 0; j < tetrahedrons; j++ )
    {
      for ( i = 0; i < 4; i++ )
      {
        printf ( "  %d", tetrahedron_vertex[i+j*4] );
      }
      printf ( "  (%d)\n", tetrahedron_label[j] );
    }
  }

  if ( 0 < hexahedrons )
  {
    printf ( "\n" );
    printf ( "  Hexahedrons:\n" );
    printf ( "\n" );
    for ( j = 0; j < hexahedrons; j++ )
    {
      for ( i = 0; i < 8; i++ )
      {
        printf ( "  %d", hexahedron_vertex[i+j*8] );
      }
      printf ( "  (%d)\n", hexahedron_label[j] );
    }
  }
  return;
}
/******************************************************************************/

void data_read ( char *filename, int dim, int vertices, int edges,
  int triangles, int quadrilaterals, int tetrahedrons, int hexahedrons,
  double vertex_coordinate[], int vertex_label[], int edge_vertex[],
  int edge_label[], int triangle_vertex[], int triangle_label[],
  int quadrilateral_vertex[], int quadrilateral_label[],
  int tetrahedron_vertex[], int tetrahedron_label[], int hexahedron_vertex[],
  int hexahedron_label[] )

/******************************************************************************/
/*
  Purpose:

    DATA_READ reads ICE data from a NETCDF file.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    23 October 2010

  Author:

    John Burkardt

  Reference:

    Pascal Frey,
    MEDIT: An interactive mesh visualization software,
    Technical Report RT-0253,
    Institut National de Recherche en Informatique et en Automatique,
    03 December 2001.

  Parameters:

    Input, char *FILENAME, the name of the file to be read.
    Ordinarily, the name should include the extension ".nc".

    Input, int DIM, the spatial dimension, which should be 2 or 3.

    Input, int VERTICES, the number of vertices.

    Input, int EDGES, the number of edges (may be 0).

    Input, int TRIANGLES, the number of triangles (may be 0).

    Input, int QUADRILATERALS, the number of quadrilaterals (may be 0).

    Input, int TETRAHEDRONS, the number of tetrahedrons (may be 0).

    Input, int HEXAHEDRONS, the number of hexahedrons (may be 0).

    Output, double VERTEX_COORDINATE[DIM*VERTICES], the coordinates
    of each vertex.

    Output, int VERTEX_LABEL[VERTICES], a label for each vertex.

    Output, int EDGE_VERTEX[2*EDGES], the vertices that form each edge.

    Output, int EDGE_LABEL[EDGES], a label for each edge.

    Output, int TRIANGLE_VERTEX[3*TRIANGLES], the vertices that form
    each triangle.

    Output, int TRIANGLE_LABEL[TRIANGLES], a label for each triangle.

    Output, int QUADRILATERAL_VERTEX[4*QUADRILATERALS], the vertices that
    form each quadrilateral.

    Output, int QUADRILATERAL_LABEL[QUADRILATERALS], a label for
    each quadrilateral.

    Output, int TETRAHEDRON_VERTEX[4*TETRAHEDRONS], the vertices that
    form each tetrahedron.

    Output, int TETRAHEDRON_LABEL[TETRAHEDRONS], a label for
    each tetrahedron.

    Output, int HEXAHEDRON_VERTEX[8*HEXAHEDRONS], the vertices that form
    each hexahedron.

    Output, int HEXAHEDRON_LABEL[HEXAHEDRONS], a label for each hexahedron.
*/
{
  int i;
  int ncid;
  int dimid;
  int mode;
  int varid;
/*
  Open the file.
*/
  mode = NC_NOCLOBBER;
  nc_open ( filename, mode, &ncid );
/*
  Vertices.
*/
  nc_inq_varid ( ncid, "Vertex_Coordinate", &varid );
  nc_get_var_double ( ncid, varid, vertex_coordinate );

  nc_inq_varid ( ncid, "Vertex_Label", &varid );
  nc_get_var_int ( ncid, varid, vertex_label );
/*
  Edges.
*/
  if ( 0 < edges )
  {
    nc_inq_varid ( ncid, "Edge_Vertex", &varid );
    nc_get_var_int ( ncid, varid, edge_vertex );

    nc_inq_varid ( ncid, "Edge_Label", &varid );
    nc_get_var_int ( ncid, varid, edge_label );
  }
/*
  Triangles.
*/
  if ( 0 < triangles )
  {
    nc_inq_varid ( ncid, "Triangle_Vertex", &varid );
    nc_get_var_int ( ncid, varid, triangle_vertex );

    nc_inq_varid ( ncid, "Triangle_Label", &varid );
    nc_get_var_int ( ncid, varid, triangle_label );
  }
/*
  Quadrilaterals.
*/
  if ( 0 < quadrilaterals )
  {
    nc_inq_varid ( ncid, "Quadrilateral_Vertex", &varid );
    nc_get_var_int ( ncid, varid, quadrilateral_vertex );

    nc_inq_varid ( ncid, "Quadrilateral_Label", &varid );
    nc_get_var_int ( ncid, varid, quadrilateral_label );
  }
/*
  Tetrahedrons.
*/
  if ( 0 < tetrahedrons )
  {
    nc_inq_varid ( ncid, "Tetrahedron_Vertex", &varid );
    nc_get_var_int ( ncid, varid, tetrahedron_vertex );

    nc_inq_varid ( ncid, "Tetrahedron_Label", &varid );
    nc_get_var_int ( ncid, varid, tetrahedron_label );
  }
/*
  Hexahedrons.
*/
  if ( 0 < hexahedrons )
  {
    nc_inq_varid ( ncid, "Hexahedron_Vertex", &varid );
    nc_get_var_int ( ncid, varid, hexahedron_vertex );

    nc_inq_varid ( ncid, "Hexahedron_Label", &varid );
    nc_get_var_int ( ncid, varid, hexahedron_label );
  }
/*
  Close the file.
*/
  nc_close ( ncid );

  return;
}
/******************************************************************************/

void mesh_write ( char *filename, int dim, int vertices, int edges,
  int triangles, int quadrilaterals, int tetrahedrons, int hexahedrons,
  double vertex_coordinate[], int vertex_label[], int edge_vertex[],
  int edge_label[], int triangle_vertex[], int triangle_label[],
  int quadrilateral_vertex[], int quadrilateral_label[],
  int tetrahedron_vertex[], int tetrahedron_label[],
  int hexahedron_vertex[], int hexahedron_label[] )

/******************************************************************************/
/*
  Purpose:

    MESH_WRITE writes mesh data to a MESH file.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    03 December 2010

  Author:

    John Burkardt

  Reference:

    Pascal Frey,
    MEDIT: An interactive mesh visualization software,
    Technical Report RT-0253,
    Institut National de Recherche en Informatique et en Automatique,
    03 December 2001.

  Parameters:

    Input, char *FILENAME, the name of the file to be created.
    Ordinarily, the name should include the extension ".mesh".

    Input, int DIM, the spatial dimension, which should be 2 or 3.

    Input, int VERTICES, the number of vertices.

    Input, int EDGES, the number of edges (may be 0).

    Input, int TRIANGLES, the number of triangles (may be 0).

    Input, int QUADRILATERALS, the number of quadrilaterals (may be 0).

    Input, int TETRAHEDRONS, the number of tetrahedrons (may be 0).

    Input, int HEXAHEDRONS, the number of hexahedrons (may be 0).

    Input, double VERTEX_COORDINATE[DIM*VERTICES], the coordinates
    of each vertex.

    Input, int VERTEX_LABEL[VERTICES], a label for each vertex.

    Input, int EDGE_VERTEX[2*EDGES], the vertices that form each edge.

    Input, int EDGE_LABEL[EDGES], a label for each edge.

    Input, int TRIANGLE_VERTEX[3*TRIANGLES], the vertices that form
    each triangle.

    Input, int TRIANGLE_LABEL[TRIANGLES], a label for each triangle.

    Input, int QUADRILATERAL_VERTEX[4*QUADRILATERALS], the vertices that
    form each quadrilateral.

    Input, int QUADRILATERAL_LABEL[QUADRILATERALS], a label for
    each quadrilateral.

    Input, int TETRAHEDRON_VERTEX[4*TETRAHEDRONS], the vertices that
    form each tetrahedron.

    Input, int TETRAHEDRON_LABEL[TETRAHEDRONS], a label for
    each tetrahedron.

    Input, int HEXAHEDRON_VERTEX[8*HEXAHEDRONS], the vertices that form
    each hexahedron.

    Input, int HEXAHEDRON_LABEL[HEXAHEDRONS], a label for each hexahedron.
*/
{
  int i;
  int j;
  FILE *output;

  output = fopen ( filename, "wt" );

  if ( !output )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "MESH_WRITE - Fatal error!\n" );
    fprintf ( stderr, "  Unable to open output file.\n" );
    exit ( 1 );
  }

  fprintf ( output, "MeshVersionFormatted 1\n" );
  fprintf ( output, "#  Created by mesh_write.c\n" );
/*
  Vertices.
*/
  fprintf ( output, "\n" );
  fprintf ( output, "Vertices\n" );
  fprintf ( output, "%d\n", vertices );
  for ( j = 0; j < vertices; j++ )
  {
    for ( i = 0; i < dim; i++ )
    {
      fprintf ( output, "  %f", vertex_coordinate[i+j*dim] );
    }
    fprintf ( output, "  %d\n", vertex_label[j] );
  }
/*
  Edges.
*/
  if ( 0 < edges )
  {
    fprintf ( output, "\n" );
    fprintf ( output, "Edges\n" );
    fprintf ( output, "%d\n", edges );
    for ( j = 0; j < edges; j++ )
    {
      for ( i = 0; i < 2; i++ )
      {
        fprintf ( output, "  %d", edge_vertex[i+j*2] );
    }
    fprintf ( output, "  %d\n", edge_label[j] );
    }
  }
/*
  Triangles.
*/
  if ( 0 < triangles )
  {
    fprintf ( output, "\n" );
    fprintf ( output, "Triangles\n" );
    fprintf ( output, "%d\n", triangles );
    for ( j = 0; j < triangles; j++ )
    {
      for ( i = 0; i < 3; i++ )
      {
        fprintf ( output, "  %d", triangle_vertex[i+j*3] );
      }
      fprintf ( output, "  %d\n", triangle_label[j] );
    }
  }
/*
  Quadrilaterals.
*/
  if ( 0 < quadrilaterals )
  {
    fprintf ( output, "\n" );
    fprintf ( output, "Quadrilaterals\n" );
    fprintf ( output, "%d\n", quadrilaterals );
    for ( j = 0; j < quadrilaterals; j++ )
    {
      for ( i = 0; i < 4; i++ )
      {
        fprintf ( output, "  %d", quadrilateral_vertex[i+j*4] );
      }
      fprintf ( output, "  %d\n", quadrilateral_label[j] );
    }
  }
/*
  Tetrahedra.
*/
  if ( 0 < tetrahedrons )
  {
    fprintf ( output, "\n" );
    fprintf ( output, "Tetrahedra\n" );
    fprintf ( output, "%d\n", tetrahedrons );
    for ( j = 0; j < tetrahedrons; j++ )
    {
      for ( i = 0; i < 4; i++ )
      {
        fprintf ( output, "  %d", tetrahedron_vertex[i+j*4] );
      }
      fprintf ( output, "  %d\n", tetrahedron_label[j] );
    }
  }
/*
  Hexahedra.
*/
  if ( 0 < hexahedrons )
  {
    fprintf ( output, "\n" );
    fprintf ( output, "Hexahedra\n" );
    fprintf ( output, "%d\n", hexahedrons );
    for ( j = 0; j < hexahedrons; j++ )
    {
      for ( i = 0; i < 8; i++ )
      {
        fprintf ( output, "  %d", hexahedron_vertex[i+j*8] );
      }
      fprintf ( output, "  %d\n", hexahedron_label[j] );
    }
  }
/*
  End
*/
  fprintf ( output, "\n" );
  fprintf ( output, "End\n" );

  fclose ( output );

  return;
}
/******************************************************************************/

void size_print ( int dim, int vertices, int edges, int triangles,
  int quadrilaterals, int tetrahedrons, int hexahedrons )

/******************************************************************************/
/*
  Purpose:

    SIZE_PRINT prints the sizes of an ICE grid dataset.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    22 October 2010

  Author:

    John Burkardt

  Reference:

    Pascal Frey,
    MEDIT: An interactive mesh visualization software,
    Technical Report RT-0253,
    Institut National de Recherche en Informatique et en Automatique,
    03 December 2001.

  Parameters:

    Input, int DIM, the spatial dimension, which should be 2 or 3.

    Input, int VERTICES, the number of vertices.

    Input, int EDGES, the number of edges (may be 0).

    Input, int TRIANGLES, the number of triangles (may be 0).

    Input, int QUADRILATERALS, the number of quadrilaterals (may be 0).

    Input, int TETRAHEDRONS, the number of tetrahedrons (may be 0).

    Input, int HEXAHEDRONS, the number of hexahedrons (may be 0).
*/
{
  printf ( "\n" );
  printf ( "  Number of dimensions = %d\n", dim );
  printf ( "  Number of vertices = %d\n", vertices );
  printf ( "  Number of edges = %d\n", edges );
  printf ( "  Number of triangles = %d\n", triangles );
  printf ( "  Number of quadrilaterals = %d\n", quadrilaterals );
  printf ( "  Number of tetrahedrons = %d\n", tetrahedrons );
  printf ( "  Number of hexahedrons = %d\n", hexahedrons );

  return;
}
/******************************************************************************/

void size_read ( char *filename, int *dim, int *vertices, int *edges,
  int *triangles, int *quadrilaterals, int *tetrahedrons, int *hexahedrons )

/******************************************************************************/
/*
  Purpose:

    SIZE_READ reads ICE sizes from a NETCDF file.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    22 October 2010

  Author:

    John Burkardt

  Reference:

    Pascal Frey,
    MEDIT: An interactive mesh visualization software,
    Technical Report RT-0253,
    Institut National de Recherche en Informatique et en Automatique,
    03 December 2001.

  Parameters:

    Input, char *FILENAME, the name of the file to be read.
    Ordinarily, the name should include the extension ".nc".

    Output, int *DIM, the spatial dimension, which should be 2 or 3.

    Output, int *VERTICES, the number of vertices.

    Output, int *EDGES, the number of edges (may be 0).

    Output, int *TRIANGLES, the number of triangles (may be 0).

    Output, int *QUADRILATERALS, the number of quadrilaterals (may be 0).

    Output, int *TETRAHEDRONS, the number of tetrahedrons (may be 0).

    Output, int *HEXAHEDRONS, the number of hexahedrons (may be 0).
*/
{
  int ncid;
  int dimid;
  int mode;
  size_t value;
/*
  Initialize everything to nothing.
*/
  *dim = 0;
  *vertices = 0;
  *edges = 0;
  *triangles = 0;
  *quadrilaterals = 0;
  *tetrahedrons = 0;
  *hexahedrons = 0;
/*
  Open the file.
*/
  mode = NC_NOCLOBBER;
  nc_open ( filename, mode, &ncid );
/*
  Get the dimension information.
*/
  if ( nc_inq_dimid ( ncid, "Dimension", &dimid ) == NC_NOERR )
  {
    nc_inq_dimlen ( ncid, dimid, &value );
    *dim = value;
  }
  if ( nc_inq_dimid ( ncid, "Vertices", &dimid ) == NC_NOERR )
  {
    nc_inq_dimlen ( ncid, dimid, &value );
    *vertices = value;
  }
  if ( nc_inq_dimid ( ncid, "Edges", &dimid ) == NC_NOERR )
  {
    nc_inq_dimlen ( ncid, dimid, &value );
    *edges = value;
  }
  if ( nc_inq_dimid ( ncid, "Triangles", &dimid ) == NC_NOERR )
  {
    nc_inq_dimlen ( ncid, dimid, &value );
    *triangles = value;
  }
  if ( nc_inq_dimid ( ncid, "Quadrilaterals", &dimid ) == NC_NOERR )
  {
    nc_inq_dimlen ( ncid, dimid, &value );
    *quadrilaterals = value;
  }
  if ( nc_inq_dimid ( ncid, "Tetrahedra", &dimid ) == NC_NOERR )
  {
    nc_inq_dimlen ( ncid, dimid, &value );
    *tetrahedrons = value;
  }
  if ( nc_inq_dimid ( ncid, "Tetrahedrons", &dimid ) == NC_NOERR )
  {
    nc_inq_dimlen ( ncid, dimid, &value );
    *tetrahedrons = value;
  }
  if ( nc_inq_dimid ( ncid, "Hexahedra", &dimid ) == NC_NOERR )
  {
    nc_inq_dimlen ( ncid, dimid, &value );
    *hexahedrons = value;
  }
  if ( nc_inq_dimid ( ncid, "Hexahedrons", &dimid ) == NC_NOERR )
  {
    nc_inq_dimlen ( ncid, dimid, &value );
    *hexahedrons = value;
  }
/*
  Close the file.
*/
  nc_close ( ncid );

  return;
}
/******************************************************************************/

void timestamp ( void )

/******************************************************************************/
/*
  Purpose:

    TIMESTAMP prints the current YMDHMS date as a time stamp.

  Example:

    31 May 2001 09:45:54 AM

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    24 September 2003

  Author:

    John Burkardt

  Parameters:

    None
*/
{
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct tm *tm;
  size_t len;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  len = strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  printf ( "%s\n", time_buffer );

  return;
# undef TIME_SIZE
}

