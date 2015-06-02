# include <stdlib.h>
# include <stdio.h>
# include <string.h>
# include <time.h>

int main ( int argc, char *argv[] );
void mesh_base_zero ( int node_num, int element_order, int element_num, 
  int element_node[] );
void timestamp ( );
void triangle_element_data_read ( char *element_file, int element_num, int element_order,
  int element_att_num, int element_node[], double element_att[] );
void triangle_element_size_read ( char *element_file, int *element_num,
  int *element_order, int *element_att_num );
void triangle_node_data_read ( char *node_file, int node_num, int node_dim,
  int node_att_num, int node_marker_num, double node_coord[], double node_att[],
  int node_marker[] );
void triangle_node_size_read ( char *node_file, int *node_num, int *node_dim,
  int *node_att_num, int *node_marker_num );
void xml_mesh2d_write ( char *xml_filename, int m, int node_num, double node_x[],
  int element_order, int element_num, int element_node[] );

/******************************************************************************/

int main ( int argc, char *argv[] )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for TRIANGLE_TO_XML.

  Discussion:

    TRIANGLE_TO_XML converts mesh data from TRIANGLE format to XML format.

  Usage:

    triangle_to_xml prefix

    where 'prefix' is the common filename prefix:

    * 'prefix'.node contains the triangle node coordinates,
    * 'prefix'.ele contains the triangle element node connectivity.
    * 'prefix'.xml will contain the mesh as an XML file.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    15 October 2014

  Author:

    John Burkardt
*/
{
  double *element_att;
  int element_att_num;
  int *element_node;
  int element_num;
  int element_order;

  int m;
  double *node_att;
  int node_att_num;
  int *node_marker;
  int node_marker_num;
  int node_num;
  double *node_x;
  char prefix[255];
  char triangle_element_filename[255];
  char triangle_node_filename[255];
  char xml_filename[255];

  timestamp ( );
  printf ( "\n" );
  printf ( "TRIANGLE_TO_XML\n" );
  printf ( "  C version:\n" );
  printf ( "  Read a mesh description created by the TRIANGLE program:\n" );
  printf ( "  * \"prefix\".node, node coordinates.\n" );
  printf ( "  * \"prefix\".ele, element connectivity.\n" );
  printf ( "  Write out an XML file for use by FENICS or DOLFIN.\n" );
  printf ( "  * \"prefix\".xml, the XML file.\n" );
/*
  Get the filename prefix.
*/
  if ( argc <= 1 ) 
  {
    printf ( "\n" );
    printf ( "TRIANGLE_TO_XML:\n" );
    printf ( "  Please enter the filename prefix.\n" );

    scanf ( "%s", prefix );
  }
  else 
  {
    strcpy ( prefix, argv[1] );
  }
/*
  Create the filenames.
*/
  strcpy ( triangle_node_filename, prefix );
  strcat ( triangle_node_filename, ".node" );
  strcpy ( triangle_element_filename, prefix );
  strcat ( triangle_element_filename, ".ele" );
  strcpy ( xml_filename, prefix );
  strcat ( xml_filename, ".xml" );
/*
  Read the triangle node size information.
*/
  triangle_node_size_read ( triangle_node_filename, &node_num, &m, &node_att_num, 
    &node_marker_num );
/*
  Read the triangle element size information.
*/
  triangle_element_size_read ( triangle_element_filename, &element_num, &element_order, 
    &element_att_num );
/*
  Report sizes.
*/
  printf ( "\n" );
  printf ( "  Size information from TRIANGLE files:\n" );
  printf ( "  Spatial dimension M = %d\n", m );
  printf ( "  Number of nodes NODE_NUM = %d\n", node_num );
  printf ( "  NODE_ATT_NUM = %d\n", node_att_num );
  printf ( "  NODE_MARKER_NUM = %d\n", node_marker_num );
  printf ( "  Number of elements ELEMENT_NUM = %d\n", element_num );
  printf ( "  Element order ELEMENT_ORDER = %d\n", element_order );
  printf ( "  ELEMENT_ATT_NUM = %d\n", element_att_num );
/*
  Allocate memory.
*/
  node_att = ( double * ) malloc ( node_att_num * node_num * sizeof ( double ) );
  node_marker = ( int * ) malloc ( node_marker_num * node_num * sizeof ( int ) );
  node_x = ( double * ) malloc ( m * node_num * sizeof ( double ) );
  element_node = ( int * ) malloc ( element_order * element_num * sizeof ( int ) );
  element_att = ( double * ) malloc ( element_att_num * element_num * sizeof ( double ) );
/*
  Read TRIANGLE data.
*/
  triangle_node_data_read ( triangle_node_filename, node_num, m, node_att_num, 
    node_marker_num, node_x, node_att, node_marker );

  triangle_element_data_read ( triangle_element_filename, element_num, element_order, 
    element_att_num, element_node, element_att );
/*
  Write XML data.
*/
  xml_mesh2d_write ( xml_filename, m, node_num, node_x, element_order,
    element_num, element_node );

  printf ( "\n" );
  printf ( "  Created XML file \"%s\".\n", xml_filename );
/*
  Free memory.
*/
  free ( element_att );
  free ( element_node );
  free ( node_att );
  free ( node_marker );
  free ( node_x );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "TRIANGLE_TO_XML:\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void mesh_base_zero ( int node_num, int element_order, int element_num, 
  int element_node[] )

/******************************************************************************/
/*
  Purpose:

    MESH_BASE_ZERO ensures that the element definition is zero-based.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    15 October 2014

  Author:

    John Burkardt

  Parameters:

    Input, int NODE_NUM, the number of nodes.

    Input, int ELEMENT_ORDER, the order of the elements.

    Input, int ELEMENT_NUM, the number of elements.

    Input/output, int ELEMENT_NODE[ELEMENT_ORDER*ELEMENT_NUM], the element
    definitions.
*/
{
  int i;
  const int i4_huge = 2147483647;
  int j;
  int node;
  int node_max;
  int node_min;
  int t;

  node_min = + i4_huge;
  node_max = - i4_huge;

  for ( j = 0; j < element_num; j++ )
  {
    for ( i = 0; i < element_order; i++ )
    {
      t = element_node[i+j*element_order];
      if ( t < node_min )
      {
        node_min = t;
      }
      if ( node_max < t )
      {
        node_max = t;
      }
    }
  }

  if ( node_min == 0 && node_max == node_num - 1 )
  {
    printf ( "\n" );
    printf ( "MESH_BASE_ZERO:\n" );
    printf ( "  The element indexing appears to be 0-based!\n" );
    printf ( "  No conversion is necessary.\n" );

  }
  else if ( node_min == 1 && node_max == node_num )
  {
    printf ( "\n" );
    printf ( "MESH_BASE_ZERO:\n" );
    printf ( "  The element indexing appears to be 1-based!\n" );
    printf ( "  This will be converted to 0-based.\n" );
    for ( j = 0; j < element_num; j++ )
    {
      for ( i = 0; i < element_order; i++ )
      {
        element_node[i+j*element_order] = element_node[i+j*element_order] - 1;
      }
    }
  }
  else
  {
    printf ( "\n" );
    printf ( "MESH_BASE_ZERO - Warning!\n" );
    printf ( "  The element indexing is not of a recognized type.\n" );
    printf ( "  NODE_MIN = %d\n", node_min );
    printf ( "  NODE_MAX = %d\n", node_max );
    printf ( "  NODE_NUM = %d\n", node_num );
  }
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
/******************************************************************************/

void triangle_element_data_read ( char *element_file, int element_num, 
  int element_order, int element_att_num, int element_node[], double element_att[] )

/******************************************************************************/
/*
  Purpose:

    TRIANGLE_ELEMENT_DATA_READ reads the data from an element file.

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
    fprintf ( stderr, "TRIANGLE_ELEMENT_DATA_READ - Fatal error!\n" );
    fprintf ( stderr, "  Could not open file \"%s\".\n", element_file );
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

void triangle_element_size_read ( char *element_file, int *element_num,
  int *element_order, int *element_att_num )

/******************************************************************************/
/*
  Purpose:

    TRIANGLE_ELEMENT_SIZE_READ reads the header information from an element file.

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
    fprintf ( stderr, "TRIANGLE_ELEMENT_SIZE_READ - Fatal error!\n" );
    fprintf ( stderr, "  Could not open file \"%s\".\n", element_file );
    exit ( 1 );
  }

  fscanf ( input, "%d  %d  %d",
    element_num, element_order, element_att_num );

  fclose ( input );

  return;
}
/******************************************************************************/

void triangle_node_data_read ( char *node_file, int node_num, int node_dim,
  int node_att_num, int node_marker_num, double node_coord[], double node_att[],
  int node_marker[] )

/******************************************************************************/
/*
  Purpose:

    TRIANGLE_NODE_DATA_READ reads the data from a node file.

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
    fprintf ( stderr, "TRIANGLE_NODE_DATA_READ - Fatal error!\n" );
    fprintf ( stderr, "  Could not open file \"%s\".\n", node_file );
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

void triangle_node_size_read ( char *node_file, int *node_num, int *node_dim,
  int *node_att_num, int *node_marker_num )

/******************************************************************************/
/*
  Purpose:

    TRIANGLE_NODE_SIZE_READ reads the header information from a node file.

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

  if ( !input )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "TRIANGLE_NODE_SIZE_READ - Fatal error!\n" );
    fprintf ( stderr, "  Could not open file \"%s\".\n", node_file );
    exit ( 1 );
  }

  fscanf ( input, "%d  %d  %d  %d",
    node_num, node_dim, node_att_num, node_marker_num );

  fclose ( input );

  return;
}
/******************************************************************************/

void xml_mesh2d_write ( char *xml_filename, int m, int node_num, 
  double node_x[], int element_order, int element_num, int element_node[] )

/******************************************************************************/
/*
  Purpose:

    XML_MESH2D_WRITE writes a 2D mesh as a DOLFIN XML file.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 May 2014

  Author:

    John Burkardt

  Reference:

    Anders Logg, Kent-Andre Mardal, Garth Wells,
    Automated Solution of Differential Equations by the Finite Element
    Method: The FEniCS Book,
    Lecture Notes in Computational Science and Engineering,
    Springer, 2011,
    ISBN13: 978-364223098

  Parameters:

    Input, char *XML_FILENAME, the name of the XML file 
    to create.

    Input, int M, the spatial dimension.

    Input, int NODE_NUM, the number of nodes.

    Input, double NODE_X[M*NODE_NUM], the node coordinates.

    Input, int ELEMENT_ORDER, the order of the elements.

    Input, int ELEMENT_NUM, the number of elements.

    Input, int ELEMENT_NODE[ELEMENT_ORDER*ELEMENT_NUM], 
    the nodes that make up each element.
*/
{
  int element;
  int element_type;
  FILE *xml;
  int i;
  int node;
/*
  Force 0-based indexing.
*/
  mesh_base_zero ( node_num, element_order, element_num, element_node );
/*
  Open the file.
*/
  xml = fopen ( xml_filename, "wt" );
/*
  Write the data.
*/
  fprintf ( xml, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n" );
  fprintf ( xml, "\n" );
  fprintf ( xml, "<dolfin xmlns:dolfin=\"http://www.fenics.org/dolfin/\">\n" );
  fprintf ( xml, "  <mesh celltype=\"triangle\" dim=\"%d\">\n", m );

  fprintf ( xml, "    <vertices size=\"%d\">\n", node_num );
  for ( node = 0; node < node_num; node++ )
  {
    fprintf ( xml, 
      "      <vertex index =\"%d\" x =\"%g\" y =\"%g\"/>\n", 
      node, node_x[0+node*m], node_x[1+node*m] );
  }
  fprintf ( xml, "    </vertices>\n" );

  fprintf ( xml, "    <cells size=\"%d\">\n", element_num );
  for ( element = 0; element < element_num; element++ )
  {
    fprintf ( xml, 
      "      <triangle index =\"%d\" v0 =\"%d\" v1 =\"%d\" v2 =\"%d\"/>\n",
      element, element_node[0+element*element_order], 
      element_node[1+element*element_order], 
      element_node[2+element*element_order] );
  }
  fprintf ( xml, "    </cells>\n" );
  fprintf ( xml, "  </mesh>\n" );
  fprintf ( xml, "</dolfin>\n" );

  fclose ( xml );

  return;
}
