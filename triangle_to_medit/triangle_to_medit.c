# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <string.h>

int main ( int argc, char *argv[] );
void element_data_read ( char *element_file, int element_num, int element_order,
  int element_att_num, int element_node[], double element_att[] );
void element_size_read ( char *element_file, int *element_num,
  int *element_order, int *element_att_num );
void mesh_write ( char *filename, int dim, int vertices, int edges,
  int triangles, int quadrilaterals, int tetrahedrons, int hexahedrons,
  double vertex_coordinate[], int vertex_label[], int edge_vertex[],
  int edge_label[], int triangle_vertex[], int triangle_label[],
  int quadrilateral_vertex[], int quadrilateral_label[],
  int tetrahedron_vertex[], int tetrahedron_label[],
  int hexahedron_vertex[], int hexahedron_label[] );
void node_data_read ( char *node_file, int node_num, int node_dim,
  int node_att_num, int node_marker_num, double node_coord[], double node_att[],
  int node_marker[] );
void node_size_read ( char *node_file, int *node_num, int *node_dim,
  int *node_att_num, int *node_marker_num );
void timestamp ( );

/******************************************************************************/

int main ( int argc, char *argv[] )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for TRIANGLE_TO_MEDIT.

  Discussion:

    The TRIANGLE program creates "node" and "element" files that define
    a triangular mesh.  A typical pair of such files might have the names
    "suv.node" and "suv.ele".

    This program reads this pair of files and creates a MEDIT mesh file, whose
    name might be "suv.mesh".

  Usage:

    triangle_to_medit prefix

    where 'prefix' is the common filename prefix so that:

    * prefix.node contains the coordinates of the nodes;
    * prefix.ele contains the indices of nodes forming each element.
    * prefix.mesh will be the MESH file created by the program.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 December 2010

  Author:

    John Burkardt
*/
{
  int dim;
  int *edge_label;
  int *edge_vertex;
  int edges;
  double *element_att;
  int element_att_num;
  char element_filename[255];
  int *element_node;
  int element_num;
  int element_order;
  int *hexahedron_label;
  int *hexahedron_vertex;
  int hexahedrons;
  int i;
  int j;
  char mesh_filename[255];
  double *node_att;
  int node_att_num;
  char node_filename[255];
  int *node_marker;
  double *node_coord;
  int node_dim;
  int node_marker_num;
  int node_num;
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
  printf ( "TRIANGLE_TO_MEDIT:\n" );
  printf ( "  C version\n" );
  printf ( "  Read a pair of NODE and ELE files created by TRIANGLE.\n" );
  printf ( "  Write a corresponding MEDIT mesh file.\n" );
/*
  Get the filename prefix.
*/
  if ( 1 <= argc )
  {
    strcpy ( prefix, argv[1] );
  }
  else
  {
    printf ( "\n" );
    printf ( "  Please enter the filename prefix:\n" );

    scanf ( "%s", prefix );
  }
/*
  Create the file names.
*/
  strcpy ( node_filename, prefix );
  strcat ( node_filename, ".node" );
  strcpy ( element_filename, prefix );
  strcat ( element_filename, ".ele" );
  strcpy ( mesh_filename, prefix );
  strcat ( mesh_filename, ".mesh" );

  printf ( "\n" );
  printf ( "  Read:\n" );
  printf ( "  * TRIANGLE node file \"%s\"\n", node_filename );
  printf ( "  * TRIANGLE element file \"%s\"\n", element_filename );
  printf ( "  Create:\n" );
  printf ( "  * MEDIT mesh file \"%s\"\n", mesh_filename );
/*
  Read the TRIANGLE NODE data.
*/
  node_size_read ( node_filename, &node_num, &node_dim, &node_att_num,
    &node_marker_num );

  node_coord = ( double * ) malloc ( 2 * node_num * sizeof ( double ) );
  node_att = ( double * ) malloc ( node_att_num * node_num * sizeof ( double ) );
  node_marker = ( int * ) malloc ( node_num * sizeof ( int ) );

  node_data_read ( node_filename, node_num, node_dim, node_att_num,
    node_marker_num, node_coord, node_att, node_marker );
/*
  Read the TRIANGLE ELE data.
*/
  element_size_read ( element_filename, &element_num, &element_order,
    &element_att_num );

  element_node = ( int * ) malloc ( element_order * element_num * sizeof ( int ) );
  element_att = ( double * ) malloc ( element_att_num * element_num * sizeof ( double ) );

  element_data_read ( element_filename, element_num, element_order,
    element_att_num, element_node, element_att );
/*
  Write the MEDIT data.
*/
  dim = 2;
  vertices = node_num;
  edges = 0;
  triangles = element_num;
  quadrilaterals = 0;
  tetrahedrons = 0;
  hexahedrons = 0;
  vertex_coordinate = ( double * ) malloc ( 2 * vertices * sizeof ( double ) );
  for ( j = 0; j < vertices; j++ )
  {
    for ( i = 0; i < 2; i++ )
    {
      vertex_coordinate[i+j*2] = node_coord[i+j*2];
    }
  }
  vertex_label = ( int * ) malloc ( vertices * sizeof ( int ) );
  for ( j = 0; j < vertices; j++ )
  {
    vertex_label[j] = node_marker[j];
  }
  edge_vertex = NULL;
  edge_label = NULL;
  triangle_vertex = ( int * ) malloc ( 3 * triangles * sizeof ( int ) );
  for ( j = 0; j < triangles; j++ )
  {
    for ( i = 0; i < 3; i++ )
    {
      triangle_vertex[i+j*3] = element_node[i+j*3];
    }
  }
  triangle_label = ( int * ) malloc ( triangles * sizeof ( int ) );
  for ( j = 0; j < triangles; j++ )
  {
    triangle_label[j] = 0;
  }
  quadrilateral_vertex = NULL;
  quadrilateral_label = NULL;
  tetrahedron_vertex = NULL;
  tetrahedron_label = NULL;
  hexahedron_vertex = NULL;
  hexahedron_label = NULL;

  mesh_write ( mesh_filename, dim, vertices, edges, triangles, quadrilaterals,
    tetrahedrons, hexahedrons, vertex_coordinate, vertex_label, edge_vertex,
    edge_label, triangle_vertex, triangle_label, quadrilateral_vertex,
    quadrilateral_label, tetrahedron_vertex, tetrahedron_label,
    hexahedron_vertex, hexahedron_label );
/*
  Free memory.
*/
  free ( element_att );
  free ( element_node );
  free ( node_att );
  free ( node_coord );
  free ( node_marker );
  free ( triangle_label );
  free ( triangle_vertex );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "TRIANGLE_TO_MEDIT:\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void element_data_read ( char *element_file, int element_num, int element_order,
  int element_att_num, int element_node[], double element_att[] )

/******************************************************************************/
/*
  Purpose:

    ELEMENT_DATA_READ reads the data from an element file.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 November 2010

  Author:

    John Burkardt

  Parameters:

    Input, char *ELEMENT_FILE, the name of the file to be read.

    Input, int ELEMENT_NUM, the number of elements.

    Input, int ELEMENT_ORDER, the order of the elements.

    Input, int ELEMENT_ATT_NUM, number of element attributes listed on each
    node record.

    Output, int ELEMENT_NODE[ELEMENT_ORDER*ELEMENT_NUM], the indices of the
    nodes that make up each element.

    Output, double ELEMENT_ATT[ELEMENT_ATT_NUM*ELEMENT_NUM], the attributes
    of each element.
*/
{
  int element;
  char *error;
  int i;
  int i1;
  int i2;
  int i3;
  FILE *input;
  int ival;
  double value;

  element = -1;

  input = fopen ( element_file, "rt" );

  if ( !input )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "ELEMENT_DATA_READ - Fatal error!\n" );
    fprintf ( stderr, "  Could not open file.\n" );
    exit ( 1 );
  }

  for ( ; ; )
  {
/*
  Read, but ignore, dimension line.
*/
    if ( element == -1 )
    {
      fscanf ( input, "%d  %d  %d", &i1, &i2, &i3 );
    }
    else
    {
      fscanf ( input, "%d", &ival );

      for ( i = 0; i < element_order; i++ )
      {
        fscanf ( input, "%d", &ival );
        element_node[i+element*element_order] = ival;
      }
      for ( i = 0; i < element_att_num; i++ )
      {
        fscanf ( input, "%lf", &value );
        element_att[i+element*element_att_num] = value;
      }
    }

    element = element + 1;

    if ( element_num <= element )
    {
      break;
    }
  }

  fclose ( input );

  return;
}
/******************************************************************************/

void element_size_read ( char *element_file, int *element_num,
  int *element_order, int *element_att_num )

/******************************************************************************/
/*
  Purpose:

    ELEMENT_SIZE_READ reads the header information from an element file.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 November 2010

  Author:

    John Burkardt

  Parameters:

    Input, char *ELEMENT_FILE, the name of the file to be read.

    Output, int *ELEMENT_NUM, the number of elements.

    Output, int *ELEMENT_ORDER, the order of the elements.

    Output, int *ELEMENT_ATT_NUM, the number of element attributes.
*/
{
  char *error;
  FILE *input;

  input = fopen ( element_file, "rt" );

  if ( !input )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "ELEMENT_SIZE_READ - Fatal error!\n" );
    fprintf ( stderr, "  Could not open file.\n" );
    exit ( 1 );
  }

  fscanf ( input, "%d  %d  %d",
    element_num, element_order, element_att_num );

  fclose ( input );

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

    24 October 2010

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

  output = fopen ( filename, "wr" );

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
  Dimension.
*/
  fprintf ( output, "\n" );
  fprintf ( output, "Dimension\n" );
  fprintf ( output, "%d\n", dim );
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
      fprintf ( output, "  %g", vertex_coordinate[i+j*dim] );
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

void node_data_read ( char *node_file, int node_num, int node_dim,
  int node_att_num, int node_marker_num, double node_coord[], double node_att[],
  int node_marker[] )

/******************************************************************************/
/*
  Purpose:

    NODE_DATA_READ reads the data from a node file.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    24 September 2003

  Author:

    John Burkardt

  Parameters:

    Input, char *NODE_FILE, the name of the node file to be read.

    Input, int NODE_NUM, the number of nodes.

    Input, int NODE_DIM, the spatial dimension.

    Input, int NODE_ATT_NUM, number of node attributes listed on each
    node record.

    Input, int NODE_MARKER_NUM, 1 if every node record includes a final
    boundary marker value.

    Output, double NODE_COORD[NODE_DIM*NODE_NUM], the nodal coordinates.

    Output, double NODE_ATT[NODE_ATT_NUM*NODE_NUM], the nodal attributes.

    Output, int NODE_MARKER[NODE_MARKER_NUM*NODE_NUM], the node markers.
*/
{
  char *error;
  int i;
  int i1;
  int i2;
  int i3;
  int i4;
  FILE *input;
  int ival;
  int node;
  double value;

  node = -1;

  input = fopen ( node_file, "rt" );

  if ( !input )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "NODE_DATA_READ - Fatal error!\n" );
    fprintf ( stderr, "  Could not open file.\n" );
    exit ( 1 );
  }

  for ( ; ; )
  {
/*
  Read, but ignore, dimension line.
*/
    if ( node == -1 )
    {
      fscanf ( input, "%d  %d  %d  %d", &i1, &i2, &i3, &i4 );
    }
    else
    {
      fscanf ( input, "%d", &ival );

      for ( i = 0; i < node_dim; i++ )
      {
        fscanf ( input, "%lf", &value );
        node_coord[i+node*node_dim] = value;
      }
      for ( i = 0; i < node_att_num; i++ )
      {
        fscanf ( input, "%lf", &value );
        node_att[i+node*node_att_num] = value;
      }
      for ( i = 0; i < node_marker_num; i++ )
      {
        fscanf ( input, "%d", &ival );
        node_marker[i+node*node_marker_num] = ival;
      }
    }

    node = node + 1;

    if ( node_num <= node )
    {
      break;
    }
  }

  fclose ( input );

  return;
}
/******************************************************************************/

void node_size_read ( char *node_file, int *node_num, int *node_dim,
  int *node_att_num, int *node_marker_num )

/******************************************************************************/
/*
  Purpose:

    NODE_SIZE_READ reads the header information from a node file.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    24 September 2003

  Author:

    John Burkardt

  Parameters:

    Input, char *NODE_FILE, the name of the node file to be read.

    Output, int *NODE_NUM, the number of nodes.

    Output, int *NODE_DIM, the spatial dimension.

    Output, int *NODE_ATT_NUM, number of node attributes listed on each
    node record.

    Output, int *NODE_MARKER_NUM, 1 if every node record includes a final
    boundary marker value.
*/
{
  char *error;
  FILE *input;

  input = fopen ( node_file, "rt" );

  fscanf ( input, "%d  %d  %d  %d",
    node_num, node_dim, node_att_num, node_marker_num );

  fclose ( input );

  return;
}
/******************************************************************************/

void timestamp ( )

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
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  printf ( "%s\n", time_buffer );

  return;
# undef TIME_SIZE
}
