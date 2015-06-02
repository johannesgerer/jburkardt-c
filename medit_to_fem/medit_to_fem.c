# include <math.h>
# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <time.h>

int main ( int argc, char *argv[] );
char ch_cap ( char ch );
int ch_eqi ( char ch1, char ch2 );
int ch_to_digit ( char ch );
void i4mat_write ( char *output_filename, int m, int n, int table[] );
void i4vec_zero ( int n, int a[] );
void mesh_data_read ( char *filename, int dim, int vertices, int edges,
  int triangles, int quadrilaterals, int tetrahedrons, int hexahedrons,
  double vertex_coordinate[], int vertex_label[], int edge_vertex[],
  int edge_label[], int triangle_vertex[], int triangle_label[],
  int quadrilateral_vertex[], int quadrilateral_label[],
  int tetrahedron_vertex[], int tetrahedron_label[],
  int hexahedron_vertex[], int hexahedron_label[] );
void mesh_size_read ( char *filename, int *dim, int *vertices, int *edges,
  int *triangles, int *quadrilaterals, int *tetrahedrons, int *hexahedrons );
void r8mat_write ( char *output_filename, int m, int n, double table[] );
void r8vec_zero ( int n, double a[] );
int s_begin ( char *s1, char *s2 );
int s_eqi ( char *s1, char *s2 );
int s_len_trim ( char *s );
void s_newline_to_null ( char *s );
int s_to_i4 ( char *s, int *last, int *error );
int s_to_i4vec ( char *s, int n, int ivec[] );
double s_to_r8 ( char *s, int *lchar, int *error );
int s_to_r8vec ( char *s, int n, double rvec[] );
void timestamp ( );

/******************************************************************************/

int main ( int argc, char *argv[] )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for MEDIT_TO_FEM.

  Discussion:

    MEDIT_TO_FEM converts mesh data from MEDIT to FEM format.

  Usage:

    medit_to_fem prefix

    where 'prefix' is the common filename prefix:

    * 'prefix'.mesh is the MEDIT mesh file.
    * 'prefix'_nodes.txt will contain the node coordinates.
    * 'prefix'_elements.txt will contain the element node connectivity.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    08 November 2014

  Author:

    John Burkardt
*/
{
  int dim;
  int *edge_label;
  int *edge_vertex;
  int edges;
  int element_num;
  int element_order;
  char fem_element_filename[255];
  char fem_node_filename[255];
  int *hexahedron_label;
  int *hexahedron_vertex;
  int hexahedrons;
  int m;
  char medit_filename[255];
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
  printf ( "MEDIT_TO_FEM\n" );
  printf ( "  C version:\n" );
  printf ( "  Read a mesh description created by the MEDIT program:\n" );
  printf ( "  * \"prefix\".mesh, the MEDIT mesh file.\n" );
  printf ( "  Write out two simple FEM files listing nodes and elements.\n" );
  printf ( "  * \"prefix\"_nodes.txt, node coordinates.\n" );
  printf ( "  * \"prefix\"_elements.txt, element connectivity.\n" );
/*
  Get the filename prefix.
*/
  if ( argc <= 1 ) 
  {
    printf ( "\n" );
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
  strcpy ( medit_filename, prefix );
  strcat ( medit_filename, ".mesh" );
  strcpy ( fem_node_filename, prefix );
  strcat ( fem_node_filename, "_nodes.txt" );
  strcpy ( fem_element_filename, prefix );
  strcat ( fem_element_filename, "_elements.txt" );
/*
  Read MEDIT sizes.
*/
  mesh_size_read ( medit_filename, &dim, &vertices, &edges, &triangles, 
    &quadrilaterals, &tetrahedrons, &hexahedrons );
/*
  Report sizes.
*/
  printf ( "\n" );
  printf ( "  Number of dimensions = %d\n", dim );
  printf ( "  Number of vertices = %d\n", vertices );
  printf ( "  Number of edges = %d\n", edges );
  printf ( "  Number of triangles = %d\n", triangles );
  printf ( "  Number of quadrilaterals = %d\n", quadrilaterals );
  printf ( "  Number of tetrahedrons = %d\n", tetrahedrons );
  printf ( "  Number of hexahedrons = %d\n", hexahedrons );
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
  Read MEDIT data.
*/
  mesh_data_read ( medit_filename, dim, vertices, edges, triangles, 
    quadrilaterals, tetrahedrons, hexahedrons, vertex_coordinate, 
    vertex_label, edge_vertex, edge_label, triangle_vertex, triangle_label, 
    quadrilateral_vertex, quadrilateral_label, tetrahedron_vertex, 
    tetrahedron_label, hexahedron_vertex, hexahedron_label );
/*
  Choose the FEM data.

  We need to assume that there is only one element type.
  If there are elements of multiple dimension, take the highest.
*/
  m = dim;
  node_num = vertices;
  r8mat_write ( fem_node_filename, dim, vertices, vertex_coordinate );

  printf ( "\n" );
  printf ( "  Created node coordinate file '%s'\n", fem_node_filename );

  if ( 0 < hexahedrons && dim == 3 )
  {
    element_order = 8;
    element_num = hexahedrons;
    i4mat_write ( fem_element_filename, element_order, element_num, 
      hexahedron_vertex );
  }
  else if ( 0 < tetrahedrons && dim == 3 )
  {
    element_order = 4;
    element_num = tetrahedrons;
    i4mat_write ( fem_element_filename, element_order, element_num, 
      tetrahedron_vertex );
  }
  else if ( 0 < quadrilaterals && dim == 2 )
  {
    element_order = 4;
    element_num = quadrilaterals;
    i4mat_write ( fem_element_filename, element_order, element_num, 
      quadrilateral_vertex );
  }
  else if ( 0 < triangles && dim == 2 )
  {
    element_order = 3;
    element_num = triangles;
    i4mat_write ( fem_element_filename, element_order, element_num, 
      triangle_vertex );
  }
  else
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "MEDIT_TO_FEM - Fatal error!\n" );
    fprintf ( stderr, "  Unexpected values for spatial dimension\n" );
    fprintf ( stderr, "  and number of nonzero objects.\n" );
    exit ( 1 );
  }

  printf ( "  Created element connectivity file '%s'\n", fem_element_filename );
/*
  Free memory.
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
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "MEDIT_TO_FEM:\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

char ch_cap ( char ch )

/******************************************************************************/
/*
  Purpose:

    CH_CAP capitalizes a single character.

  Discussion:

    This routine should be equivalent to the library "toupper" function.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 July 1998

  Author:

    John Burkardt

  Parameters:

    Input, char CH, the character to capitalize.

    Output, char CH_CAP, the capitalized character.
*/
{
  if ( 97 <= ch && ch <= 122 )
  {
    ch = ch - 32;
  }

  return ch;
}
/******************************************************************************/

int ch_eqi ( char ch1, char ch2 )

/******************************************************************************/
/*
  Purpose:

    CH_EQI is TRUE (1) if two characters are equal, disregarding case.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 June 2003

  Author:

    John Burkardt

  Parameters:

    Input, char CH1, CH2, the characters to compare.

    Output, int CH_EQI, is TRUE (1) if the two characters are equal,
    disregarding case and FALSE (0) otherwise.
*/
{
  int value;

  if ( 97 <= ch1 && ch1 <= 122 )
  {
    ch1 = ch1 - 32;
  }
  if ( 97 <= ch2 && ch2 <= 122 )
  {
    ch2 = ch2 - 32;
  }
  if ( ch1 == ch2 )
  {
    value = 1;
  }
  else
  {
    value = 0;
  }
  return value;
}
/******************************************************************************/

int ch_to_digit ( char ch )

/******************************************************************************/
/*
  Purpose:

    CH_TO_DIGIT returns the value of a base 10 digit.

  Example:

     CH  DIGIT
    ---  -----
    '0'    0
    '1'    1
    ...  ...
    '9'    9
    ' '    0
    'X'   -1

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 June 2003

  Author:

    John Burkardt

  Parameters:

    Input, char CH, the decimal digit, '0' through '9' or blank are legal.

    Output, int CH_TO_DIGIT, the corresponding value.  If the
    character was 'illegal', then DIGIT is -1.
*/
{
  int digit;

  if ( '0' <= ch && ch <= '9' )
  {
    digit = ch - '0';
  }
  else if ( ch == ' ' )
  {
    digit = 0;
  }
  else
  {
    digit = -1;
  }

  return digit;
}
/******************************************************************************/

void i4mat_write ( char *output_filename, int m, int n, int table[] )

/******************************************************************************/
/*
  Purpose:

    I4MAT_WRITE writes an I4MAT file.

  Discussion:

    An I4MAT is an array of I4's.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    01 June 2009

  Author:

    John Burkardt

  Parameters:

    Input, char *OUTPUT_FILENAME, the output filename.

    Input, int M, the spatial dimension.

    Input, int N, the number of points.

    Input, int TABLE[M*N], the data.
*/
{
  int i;
  int j;
  FILE *output;
/*
  Open the file.
*/
  output = fopen ( output_filename, "wt" );

  if ( !output )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "I4MAT_WRITE - Fatal error!\n" );
    fprintf ( stderr, "  Could not open the file \"%s\".\n", output_filename );
    exit ( 1 );
  }
/*
  Write the data.
*/
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      fprintf ( output, "  %d", table[i+j*m] );
    }
    fprintf ( output, "\n" );
  }
/*
  Close the file.
*/
  fclose ( output );

  return;
}
/******************************************************************************/

void i4vec_zero ( int n, int a[] )

/******************************************************************************/
/*
  Purpose:

    I4VEC_ZERO zeroes an I4VEC.

  Discussion:

    An I4VEC is a vector of I4's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    03 July 2005

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the vector.

    Output, int A[N], a vector of zeroes.
*/
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    a[i] = 0;
  }
  return;
}
/******************************************************************************/

void mesh_data_read ( char *filename, int dim, int vertices, int edges,
  int triangles, int quadrilaterals, int tetrahedrons, int hexahedrons,
  double vertex_coordinate[], int vertex_label[], int edge_vertex[],
  int edge_label[], int triangle_vertex[], int triangle_label[],
  int quadrilateral_vertex[], int quadrilateral_label[],
  int tetrahedron_vertex[], int tetrahedron_label[],
  int hexahedron_vertex[], int hexahedron_label[] )

/******************************************************************************/
/*
  Purpose:

    MESH_DATA_READ reads data from a MESH file.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 November 2010

  Author:

    John Burkardt

  Reference:

    Pascal Frey,
    MEDIT: An interactive mesh visualization software,
    Technical Report RT-0253,
    Institut National de Recherche en Informatique et en Automatique,
    03 December 2001.

  Parameters:

    Input, char *FILENAME, the name of the MESH file.

    Input, int DIM, the spatial dimension, which should be 2 or 3.

    Input, int VERTICES, the number of vertices.

    Input, int EDGES, the number of edges (may be 0).

    Input, int TRIANGLES, the number of triangles (may be 0).

    Input, int QUADRILATERALS, the number of quadrilaterals
    (may be 0).

    Input, int TETRAHEDRONS, the number of tetrahedrons
    (may be 0).

    Input, int HEXAHEDRONS, the number of hexahedrons
    (may be 0).

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
  int dim2;
  int edge;
  int edges2;
  char *error;
  int hexahedron;
  int hexahedrons2;
  int i;
  int i4vec[9];
  int ierror;
  FILE *input;
  char keyword[255];
  int length;
  int line_num;
  int quadrilateral;
  int quadrilaterals2;
  double r8vec[9];
  int tetrahedron;
  int tetrahedrons2;
  char text[255];
  int triangle;
  int triangles2;
  int vertex;
  int vertices2;
/*
  Initialize everything to nothing.
*/
  i4vec_zero ( edges, edge_label );
  i4vec_zero ( 2 * edges, edge_vertex );
  i4vec_zero ( hexahedrons, hexahedron_label );
  i4vec_zero ( 8 * hexahedrons, hexahedron_vertex );
  i4vec_zero ( quadrilaterals, quadrilateral_label );;
  i4vec_zero ( 4 * quadrilaterals, quadrilateral_vertex );
  i4vec_zero ( tetrahedrons, tetrahedron_label );
  i4vec_zero ( 4 * tetrahedrons, tetrahedron_vertex );
  i4vec_zero ( triangles, triangle_label );
  i4vec_zero ( 3 * triangles, triangle_vertex );
  r8vec_zero ( dim * vertices, vertex_coordinate );
  i4vec_zero ( vertices, vertex_label );
/*
  Open the file.
*/
  input = fopen ( filename, "rt" );

  if ( !input )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "MESH_DATA_READ - Fatal error!\n" );
    fprintf ( stderr, "  Could not open file.\n" );
    exit ( 1 );
  }
/*
  Read lines til you get alphanumerics and determine a "mode"
*/
  line_num = 0;
  strcpy ( keyword, "NONE" );

  for ( ; ; )
  {
    error = fgets ( text, 255, input );

    if ( !error )
    {
      break;
    }

    s_newline_to_null ( text );

    line_num = line_num + 1;

    if ( s_len_trim ( text ) == 0 )
    {
      strcpy ( keyword, "NONE" );
      continue;
    }

    if ( text[0] == '#' )
    {
      continue;
    }
/*
  Remove initial blanks.
*/

/*
  Expecting a keyword.
*/
        if ( s_eqi ( text, "CORNERS" ) )
    {
      strcpy ( keyword, "SKIP" );
    }
    else if ( s_eqi ( text, "DIMENSION" ) )
    {
      strcpy ( keyword, "DIMENSION" );
    }
    else if ( s_eqi ( text, "EDGES" ) )
    {
      strcpy ( keyword, "EDGES" );
    }
    else if ( s_eqi ( text, "END" ) )
    {
      printf ( "\n" );
      printf ( "  END statement encountered.\n" );
      break;
    }
    else if ( s_eqi ( text, "HEXAHEDRA" ) ||
              s_eqi ( text, "HEXAHEDRONS" ) )
    {
      strcpy ( keyword, "HEXAHEDRONS" );
    }
    else if ( s_begin ( text, "MESHVERSIONFORMATTED" ) )
    {
    }
    else if ( s_eqi ( text, "NORMALATQUADRILATERALVERTICES" ) )
    {
      strcpy ( keyword, "SKIP" );
    }
    else if ( s_eqi ( text, "NORMALATTRIANGLEVERTICES" ) )
    {
      strcpy ( keyword, "SKIP" );
    }
    else if ( s_eqi ( text, "NORMALATVERTICES" ) )
    {
      strcpy ( keyword, "SKIP" );
    }
    else if ( s_eqi ( text, "NORMALS" ) )
    {
      strcpy ( keyword, "SKIP" );
    }
    else if ( s_eqi ( text, "QUADRILATERALS" ) )
    {
      strcpy ( keyword, "QUADRILATERALS" );
    }
    else if ( s_eqi ( text, "REQUIREDEDGES" ) )
    {
      strcpy ( keyword, "SKIP" );
    }
    else if ( s_eqi ( text, "REQUIREDVERTICES" ) )
    {
      strcpy ( keyword, "SKIP" );
    }
    else if ( s_eqi ( text, "RIDGES" ) )
    {
      strcpy ( keyword, "SKIP" );
    }
    else if ( s_eqi ( text, "TANGENTATEDGES" ) )
    {
      strcpy ( keyword, "SKIP" );
    }
    else if ( s_eqi ( text, "TANGENTS" ) )
    {
      strcpy ( keyword, "SKIP" );
    }
    else if ( s_eqi ( text, "TETRAHEDRA" ) ||
              s_eqi ( text, "TETRAHEDRONS" ) )
    {
      strcpy ( keyword, "TETRAHEDRONS" );
    }
    else if ( s_eqi ( text, "TRIANGLES" ) )
    {
      strcpy ( keyword, "TRIANGLES" );
    }
    else if ( s_eqi ( text, "VERTICES" ) )
    {
      strcpy ( keyword, "VERTICES" );
    }
/*
  Presumably, numeric data to be processed by keyword.
*/
    else if ( s_eqi ( keyword, "DIMENSION" ) )
    {
      dim2 = atoi ( text );
      strcpy ( keyword, "NONE" );
    }
    else if ( s_eqi ( keyword, "EDGES" ) )
    {
      edges2 = atoi ( text );
      strcpy ( keyword, "EDGE_VERTEX" );
      edge = 0;
    }
    else if ( s_eqi ( keyword, "EDGE_VERTEX" ) )
    {
      s_to_i4vec ( text, 3, i4vec );
      for ( i = 0; i < 2; i++ )
      {
        edge_vertex[i+edge*2] = i4vec[i];
      }
      edge_label[edge] = i4vec[2];
      edge = edge + 1;
    }
    else if ( s_eqi ( keyword, "HEXAHEDRONS" ) )
    {
      hexahedrons2 = atoi ( text );
      strcpy ( keyword, "HEXAHEDRON_VERTEX" );
      hexahedron = 0;
    }
    else if ( s_eqi ( keyword, "HEXAHEDRON_VERTEX" ) )
    {
      s_to_i4vec ( text, 9, i4vec );
      for ( i = 0; i < 8; i++ )
      {
        hexahedron_vertex[i+hexahedron*8] = i4vec[i];
      }
      hexahedron_label[hexahedron] = i4vec[8];
      hexahedron = hexahedron + 1;
    }
    else if ( s_eqi ( keyword, "QUADRILATERALS" ) )
    {
      quadrilaterals2 = atoi ( text );
      strcpy ( keyword, "QUADRILATERAL_VERTEX" );
      quadrilateral = 0;
    }
    else if ( s_eqi ( keyword, "QUADRILATERAL_VERTEX" ) )
    {
      s_to_i4vec ( text, 5, i4vec );
      for ( i = 0; i < 4; i++ )
      {
        quadrilateral_vertex[i+quadrilateral*4] = i4vec[i];
      }
      quadrilateral_label[quadrilateral] = i4vec[4];
      quadrilateral = quadrilateral + 1;
    }
    else if ( s_eqi ( keyword, "TETRAHEDRONS" ) )
    {
      tetrahedrons2 = atoi ( text );
      strcpy ( keyword, "TETRAHEDRON_VERTEX" );
      tetrahedron = 0;
    }
    else if ( s_eqi ( keyword, "TETRAHEDRON_VERTEX" ) )
    {
      s_to_i4vec ( text, 5, i4vec );
      for ( i = 0; i < 4; i++ )
      {
        tetrahedron_vertex[i+tetrahedron*4] = i4vec[i];
      }
      tetrahedron_label[tetrahedron] = i4vec[4];
      tetrahedron = tetrahedron + 1;
    }
    else if ( s_eqi ( keyword, "TRIANGLES" ) )
    {
      triangles2 = atoi ( text );
      strcpy ( keyword, "TRIANGLE_VERTEX" );
      triangle = 0;
    }
    else if ( s_eqi ( keyword, "TRIANGLE_VERTEX" ) )
    {
      s_to_i4vec ( text, 4, i4vec );
      for ( i = 0; i < 3; i++ )
      {
        triangle_vertex[i+triangle*3] = i4vec[i];
      }
      triangle_label[triangle] = i4vec[3];
      triangle = triangle + 1;
    }
    else if ( s_eqi ( keyword, "VERTICES" ) )
    {
      vertices2 = atoi ( text );
      strcpy ( keyword, "VERTEX_COORDINATE" );
      vertex = 0;
    }
    else if ( s_eqi ( keyword, "VERTEX_COORDINATE" ) )
    {
      s_to_r8vec ( text, dim + 1, r8vec );
      for ( i = 0; i < dim; i++ )
      {
        vertex_coordinate[i+vertex*dim] = r8vec[i];
      }
      vertex_label[vertex] = ( int ) r8vec[dim];
      vertex = vertex + 1;
    }
    else if ( s_eqi ( keyword, "SKIP" ) )
    {
    }
    else
    {
      fprintf ( stderr, "\n" );
      fprintf ( stderr, "MESH_DATA_READ - Fatal error!\n" );
      fprintf ( stderr, "  Could not find keyword while reading line %d\n",
        line_num );
      fprintf ( stderr, "\"%s\".\n", text );
      exit ( 1 );
    }
  }
/*
  Close the file.
*/
  fclose ( input );

  printf ( "\n" );
  printf ( "  Read %d lines from \"%s\".\n", line_num, filename );

  return;
}
/******************************************************************************/

void mesh_size_read ( char *filename, int *dim, int *vertices, int *edges,
  int *triangles, int *quadrilaterals, int *tetrahedrons, int *hexahedrons )

/******************************************************************************/
/*
  Purpose:

    MESH_SIZE_READ reads sizes from a MESH file.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    22 November 2010

  Author:

    John Burkardt

  Reference:

    Pascal Frey,
    MEDIT: An interactive mesh visualization software,
    Technical Report RT-0253,
    Institut National de Recherche en Informatique et en Automatique,
    03 December 2001.

  Parameters:

    Input, char *FILENAME, the name of the MESH file.

    Output, int *DIM, the spatial dimension, which should be 2 or 3.

    Output, int *VERTICES, the number of vertices.

    Output, int *EDGES, the number of edges (may be 0).

    Output, int *TRIANGLES, the number of triangles (may be 0).

    Output, int *QUADRILATERALS, the number of quadrilaterals
    (may be 0).

    Output, int *TETRAHEDRONS, the number of tetrahedrons
    (may be 0).

    Output, int *HEXAHEDRONS, the number of hexahedrons
    (may be 0).
*/
{
  char *error;
  int ierror;
  FILE *input;
  char keyword[255];
  int length;
  int line_num;

  char text[255];
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
  input = fopen ( filename, "rt" );

  if ( !input )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "MESH_SIZE_READ - Fatal error!\n" );
    fprintf ( stderr, "  Could not open file.\n" );
    exit ( 1 );
  }
/*
  Read lines til you get alphanumerics and determine a "mode"
*/
  line_num = 0;
  strcpy ( keyword, "NONE" );

  for ( ; ; )
  {
    error = fgets ( text, 255, input );

    if ( !error )
    {
      break;
    }

    s_newline_to_null ( text );

    line_num = line_num + 1;

    if ( s_len_trim ( text ) == 0 )
    {
      strcpy ( keyword, "NONE" );
      continue;
    }

    if ( text[0] == '#' )
    {
      continue;
    }
/*
  Remove initial blanks.
*/

/*
  Expecting a keyword.
*/
        if ( s_eqi ( text, "CORNERS" ) )
    {
      strcpy ( keyword, "SKIP" );
    }
    else if ( s_eqi ( text, "DIMENSION" ) )
    {
      strcpy ( keyword, "DIMENSION" );
    }
    else if ( s_eqi ( text, "EDGES" ) )
    {
      strcpy ( keyword, "EDGES" );
    }
    else if ( s_eqi ( text, "END" ) )
    {
      printf ( "\n" );
      printf ( "  END statement encountered.\n" );
      break;
    }
    else if ( s_eqi ( text, "HEXAHEDRA" ) ||
              s_eqi ( text, "HEXAHEDRONS" ) )
    {
      strcpy ( keyword, "HEXAHEDRONS" );
    }
    else if ( s_begin ( text, "MESHVERSIONFORMATTED" ) )
    {
    }
    else if ( s_eqi ( text, "NORMALATQUADRILATERALVERTICES" ) )
    {
      strcpy ( keyword, "SKIP" );
    }
    else if ( s_eqi ( text, "NORMALATTRIANGLEVERTICES" ) )
    {
      strcpy ( keyword, "SKIP" );
    }
    else if ( s_eqi ( text, "NORMALATVERTICES" ) )
    {
      strcpy ( keyword, "SKIP" );
    }
    else if ( s_eqi ( text, "NORMALS" ) )
    {
      strcpy ( keyword, "SKIP" );
    }
    else if ( s_eqi ( text, "QUADRILATERALS" ) )
    {
      strcpy ( keyword, "QUADRILATERALS" );
    }
    else if ( s_eqi ( text, "REQUIREDEDGES" ) )
    {
      strcpy ( keyword, "SKIP" );
    }
    else if ( s_eqi ( text, "REQUIREDVERTICES" ) )
    {
      strcpy ( keyword, "SKIP" );
    }
    else if ( s_eqi ( text, "RIDGES" ) )
    {
      strcpy ( keyword, "SKIP" );
    }
    else if ( s_eqi ( text, "TANGENTATEDGES" ) )
    {
      strcpy ( keyword, "SKIP" );
    }
    else if ( s_eqi ( text, "TANGENTS" ) )
    {
      strcpy ( keyword, "SKIP" );
    }
    else if ( s_eqi ( text, "TETRAHEDRA" ) ||
              s_eqi ( text, "TETRAHEDRONS" ) )
    {
      strcpy ( keyword, "TETRAHEDRONS" );
    }
    else if ( s_eqi ( text, "TRIANGLES" ) )
    {
      strcpy ( keyword, "TRIANGLES" );
    }
    else if ( s_eqi ( text, "VERTICES" ) )
    {
      strcpy ( keyword, "VERTICES" );
    }
/*
  Presumably, numeric data to be processed by keyword.
*/
    else if ( s_eqi ( keyword, "DIMENSION" ) )
    {
      *dim = atoi ( text );
      strcpy ( keyword, "NONE" );
    }
    else if ( s_eqi ( keyword, "EDGES" ) )
    {
      *edges = atoi ( text );
      strcpy ( keyword, "EDGE_VERTEX" );
    }
    else if ( s_eqi ( keyword, "EDGE_VERTEX" ) )
    {
    }
    else if ( s_eqi ( keyword, "HEXAHEDRONS" ) )
    {
      *hexahedrons = atoi ( text );
      strcpy ( keyword, "HEXAHEDRON_VERTEX" );
    }
    else if ( s_eqi ( keyword, "HEXAHEDRON_VERTEX" ) )
    {
    }
    else if ( s_eqi ( keyword, "QUADRILATERALS" ) )
    {
      *quadrilaterals = atoi ( text );
      strcpy ( keyword, "QUADRILATERAL_VERTEX" );
    }
    else if ( s_eqi ( keyword, "QUADRILATERAL_VERTEX" ) )
    {
    }
    else if ( s_eqi ( keyword, "TETRAHEDRONS" ) )
    {
      *tetrahedrons = atoi ( text );
      strcpy ( keyword, "TETRAHEDRON_VERTEX" );
    }
    else if ( s_eqi ( keyword, "TETRAHEDRON_VERTEX" ) )
    {
    }
    else if ( s_eqi ( keyword, "TRIANGLES" ) )
    {
      *triangles = atoi ( text );
      strcpy ( keyword, "TRIANGLE_VERTEX" );
    }
    else if ( s_eqi ( keyword, "TRIANGLE_VERTEX" ) )
    {
    }
    else if ( s_eqi ( keyword, "VERTICES" ) )
    {
      *vertices = atoi ( text );
      strcpy ( keyword, "VERTEX_COORDINATE" );
    }
    else if ( s_eqi ( keyword, "VERTEX_COORDINATE" ) )
    {
    }
    else if ( s_eqi ( keyword, "SKIP" ) )
    {
    }
    else
    {
      fprintf ( stderr, "\n" );
      fprintf ( stderr, "MESH_SIZE_READ - Fatal error!\n" );
      fprintf ( stderr, "  Could not find keyword while reading line %d\n",
        line_num );
      fprintf ( stderr, "\"%s\".\n", text );
      exit ( 1 );
    }
  }
/*
  Close the file.
*/
  fclose ( input );

  printf ( "\n" );
  printf ( "  Read %d lines from \"%s\".\n", line_num, filename );

  return;
}
/******************************************************************************/

void r8mat_write ( char *output_filename, int m, int n, double table[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_WRITE writes an R8MAT file.

  Discussion:

    An R8MAT is an array of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    01 June 2009

  Author:

    John Burkardt

  Parameters:

    Input, char *OUTPUT_FILENAME, the output filename.

    Input, int M, the spatial dimension.

    Input, int N, the number of points.

    Input, double TABLE[M*N], the data.
*/
{
  int i;
  int j;
  FILE *output;
/*
  Open the file.
*/
  output = fopen ( output_filename, "wt" );

  if ( !output )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8MAT_WRITE - Fatal error!\n" );
    fprintf ( stderr, 
      "  Could not open the output file \"%s\".\n", output_filename );
    exit ( 1 );
  }
/*
  Write the data.
*/
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      fprintf ( output, "  %g", table[i+j*m] );
    }
    fprintf ( output, "\n" );
  }
/*
  Close the file.
*/
  fclose ( output );

  return;
}
/******************************************************************************/

void r8vec_zero ( int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_ZERO zeroes an R8VEC.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 August 2008

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the vector.

    Output, double A[N], a vector of zeroes.
*/
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    a[i] = 0.0;
  }
  return;
}
/******************************************************************************/

int s_begin ( char *s1, char *s2 )

/******************************************************************************/
/*
  Purpose:

    S_BEGIN reports whether string 1 begins with string 2.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    22 November 2010

  Author:

    John Burkardt

  Parameters:

    Input, char *S1, char *S2, two strings.

    Output, int S_BEGIN, is true if S1 is the same as S2 up to
    the end of S2, and false otherwise.
*/
{
  int i;
  int n;
  int n1;
  int n2;

  n1 = strlen ( s1 );
  n2 = strlen ( s2 );

  if ( n1 < n2 )
  {
    return 0;
  }

  for ( i = 0; i < n2; i++ )
  {
    if ( ch_cap ( s1[i] ) != ch_cap ( s2[i] ) )
    {
      return 0;
    }
  }
  return 1;
}
/******************************************************************************/

int s_eqi ( char *s1, char *s2 )

/******************************************************************************/
/*
  Purpose:

    S_EQI reports whether two strings are equal, ignoring case.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    03 June 2008

  Author:

    John Burkardt

  Parameters:

    Input, char *S1, char *S2, pointers to two strings.

    Output, int S_EQI, is true if the strings are equal.
*/
{
  int i;
  int nchar;
  int nchar1;
  int nchar2;

  nchar1 = strlen ( s1 );
  nchar2 = strlen ( s2 );
  if ( nchar1 < nchar2 )
  {
    nchar = nchar1;
  }
  else
  {
    nchar = nchar2;
  }
/*
  The strings are not equal if they differ over their common length.
*/
  for ( i = 0; i < nchar; i++ )
  {

    if ( ch_cap ( s1[i] ) != ch_cap ( s2[i] ) )
    {
      return 0;
    }
  }
/*
  The strings are not equal if the longer one includes nonblanks
  in the tail.
*/
  if ( nchar < nchar1 )
  {
    for ( i = nchar; i < nchar1; i++ )
    {
      if ( s1[i] != ' ' )
      {
        return 0;
      }
    }
  }
  else if ( nchar < nchar2 )
  {
    for ( i = nchar; i < nchar2; i++ )
    {
      if ( s2[i] != ' ' )
      {
        return 0;
      }
    }
  }

  return 1;
}
/******************************************************************************/

int s_len_trim ( char *s )

/******************************************************************************/
/*
  Purpose:

    S_LEN_TRIM returns the length of a string to the last nonblank.

  Discussion:

    It turns out that I also want to ignore the '\n' character!

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 October 2014

  Author:

    John Burkardt

  Parameters:

    Input, char *S, a pointer to a string.

    Output, int S_LEN_TRIM, the length of the string to the last nonblank.
    If S_LEN_TRIM is 0, then the string is entirely blank.
*/
{
  int n;
  char *t;

  n = strlen ( s );
  t = s + strlen ( s ) - 1;

  while ( 0 < n )
  {
    if ( *t != ' ' && *t != '\n' )
    {
      return n;
    }
    t--;
    n--;
  }

  return n;
}
/******************************************************************************/

void s_newline_to_null ( char *s )

/******************************************************************************/
/*
  Purpose:

    S_NEWLINE_TO_NULL replaces carriage returns or newlines by nulls.

  Discussion:

    The function FGETS will read a string containing a line of text read from
    input.  However, the string will include the linefeed character '/n', or,
    for a PC-formatted file, the carriage return and linefeed pair '/r' + '/n'.

    It may be desirable that the string not contain these characters.  The
    easiest way to deal with this is simply to replace the first instance of
    '/r' or '/n' by a null character, which terminates the string.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    29 November 2010

  Author:

    John Burkardt

  Parameters:

    Input/output, char *S, a pointer to a string.  On output, the first
    carriage return or line feed encountered has been replaced by a null.
*/
{
  int i;
  int n;

  n = strlen ( s );

  for ( i = 0; i < n; i++ )
  {
/*
  Handle carriage return.
*/
    if ( s[i] == '\r' )
    {
      s[i] = '\0';
      return;
    }
/*
  Handle linefeed.
*/
    if ( s[i] == '\n' )
    {
      s[i] = '\0';
      return;
    }
  }

  return;
}
/******************************************************************************/

int s_to_i4 ( char *s, int *last, int *error )

/******************************************************************************/
/*
  Purpose:

    S_TO_I4 reads an I4 from a string.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 June 2003

  Author:

    John Burkardt

  Parameters:

    Input, char *S, a string to be examined.

    Output, int *LAST, the last character of S used to make IVAL.

    Output, int *ERROR is TRUE (1) if an error occurred and FALSE (0) otherwise.

    Output, int *S_TO_I4, the integer value read from the string.
    If the string is blank, then IVAL will be returned 0.
*/
{
  char c;
  int i;
  int isgn;
  int istate;
  int ival;

  *error = 0;
  istate = 0;
  isgn = 1;
  i = 0;
  ival = 0;

  while ( *s )
  {
    c = s[i];
    i = i + 1;
/*
  Haven't read anything.
*/
    if ( istate == 0 )
    {
      if ( c == ' ' )
      {
      }
      else if ( c == '-' )
      {
        istate = 1;
        isgn = -1;
      }
      else if ( c == '+' )
      {
        istate = 1;
        isgn = + 1;
      }
      else if ( '0' <= c && c <= '9' )
      {
        istate = 2;
        ival = c - '0';
      }
      else
      {
        *error = 1;
        return ival;
      }
    }
/*
  Have read the sign, expecting digits.
*/
    else if ( istate == 1 )
    {
      if ( c == ' ' )
      {
      }
      else if ( '0' <= c && c <= '9' )
      {
        istate = 2;
        ival = c - '0';
      }
      else
      {
        *error = 1;
        return ival;
      }
    }
/*
  Have read at least one digit, expecting more.
*/
    else if ( istate == 2 )
    {
      if ( '0' <= c && c <= '9' )
      {
        ival = 10 * (ival) + c - '0';
      }
      else
      {
        ival = isgn * ival;
        *last = i - 1;
        return ival;
      }
    }
  }
/*
  If we read all the characters in the string, see if we're OK.
*/
  if ( istate == 2 )
  {
    ival = isgn * ival;
    *last = s_len_trim ( s );
  }
  else
  {
    *error = 1;
    *last = 0;
  }

  return ival;
}
/******************************************************************************/

int s_to_i4vec ( char *s, int n, int ivec[] )

/******************************************************************************/
/*
  Purpose:

    S_TO_I4VEC reads an I4VEC from a string.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 February 2001

  Author:

    John Burkardt

  Parameters:

    Input, char *S, the string to be read.

    Input, int N, the number of values expected.

    Output, int IVEC[N], the values read from the string.

    Output, int S_TO_I4VEC, is TRUE (1) if an error occurred and
    FALSE (0) otherwise.
*/
{
  int error;
  int i;
  int lchar;

  error = 0;

  for ( i = 0; i < n; i++ )
  {
    ivec[i] = s_to_i4 ( s, &lchar, &error );

    if ( error )
    {
      return error;
    }
    s = s + lchar;
  }

  return error;
}
/******************************************************************************/

double s_to_r8 ( char *s, int *lchar, int *error )

/******************************************************************************/
/*
  Purpose:

    S_TO_R8 reads an R8 value from a string.

  Discussion:

    We have had some trouble with input of the form 1.0E-312.
    For now, let's assume anything less than 1.0E-20 is zero.

    This routine will read as many characters as possible until it reaches
    the end of the string, or encounters a character which cannot be
    part of the real number.

    Legal input is:

       1 blanks,
       2 '+' or '-' sign,
       2.5 spaces
       3 integer part,
       4 decimal point,
       5 fraction part,
       6 'E' or 'e' or 'D' or 'd', exponent marker,
       7 exponent sign,
       8 exponent integer part,
       9 exponent decimal point,
      10 exponent fraction part,
      11 blanks,
      12 final comma or semicolon.

    with most quantities optional.

  Example:

    S                 R

    '1'               1.0
    '     1   '       1.0
    '1A'              1.0
    '12,34,56'        12.0
    '  34 7'          34.0
    '-1E2ABCD'        -100.0
    '-1X2ABCD'        -1.0
    ' 2E-1'           0.2
    '23.45'           23.45
    '-4.2E+2'         -420.0
    '17d2'            1700.0
    '-14e-2'         -0.14
    'e2'              100.0
    '-12.73e-9.23'   -12.73 * 10.0^(-9.23)

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    24 June 2005

  Author:

    John Burkardt

  Parameters:

    Input, char *S, the string containing the
    data to be read.  Reading will begin at position 1 and
    terminate at the end of the string, or when no more
    characters can be read to form a legal real.  Blanks,
    commas, or other nonnumeric data will, in particular,
    cause the conversion to halt.

    Output, int *LCHAR, the number of characters read from
    the string to form the number, including any terminating
    characters such as a trailing comma or blanks.

    Output, int *ERROR, is TRUE (1) if an error occurred and FALSE (0)
    otherwise.

    Output, double S_TO_R8, the value that was read from the string.
*/
{
  char c;
  int ihave;
  int isgn;
  int iterm;
  int jbot;
  int jsgn;
  int jtop;
  int nchar;
  int ndig;
  double r;
  double rbot;
  double rexp;
  double rtop;
  char TAB = 9;

  nchar = s_len_trim ( s );
  *error = 0;
  r = 0.0;
  *lchar = -1;
  isgn = 1;
  rtop = 0.0;
  rbot = 1.0;
  jsgn = 1;
  jtop = 0;
  jbot = 1;
  ihave = 1;
  iterm = 0;

  for ( ; ; )
  {
    c = s[*lchar+1];
    *lchar = *lchar + 1;
/*
  Blank or TAB character.
*/
    if ( c == ' ' || c == TAB )
    {
      if ( ihave == 2 )
      {
      }
      else if ( ihave == 6 || ihave == 7 )
      {
        iterm = 1;
      }
      else if ( 1 < ihave )
      {
        ihave = 11;
      }
    }
/*
  Comma.
*/
    else if ( c == ',' || c == ';' )
    {
      if ( ihave != 1 )
      {
        iterm = 1;
        ihave = 12;
        *lchar = *lchar + 1;
      }
    }
/*
  Minus sign.
*/
    else if ( c == '-' )
    {
      if ( ihave == 1 )
      {
        ihave = 2;
        isgn = -1;
      }
      else if ( ihave == 6 )
      {
        ihave = 7;
        jsgn = -1;
      }
      else
      {
        iterm = 1;
      }
    }
/*
  Plus sign.
*/
    else if ( c == '+' )
    {
      if ( ihave == 1 )
      {
        ihave = 2;
      }
      else if ( ihave == 6 )
      {
        ihave = 7;
      }
      else
      {
        iterm = 1;
      }
    }
/*
  Decimal point.
*/
    else if ( c == '.' )
    {
      if ( ihave < 4 )
      {
        ihave = 4;
      }
      else if ( 6 <= ihave && ihave <= 8 )
      {
        ihave = 9;
      }
      else
      {
        iterm = 1;
      }
    }
/*
  Exponent marker.
*/
    else if ( ch_eqi ( c, 'E' ) || ch_eqi ( c, 'D' ) )
    {
      if ( ihave < 6 )
      {
        ihave = 6;
      }
      else
      {
        iterm = 1;
      }
    }
/*
  Digit.
*/
    else if ( ihave < 11 && '0' <= c && c <= '9' )
    {
      if ( ihave <= 2 )
      {
        ihave = 3;
      }
      else if ( ihave == 4 )
      {
        ihave = 5;
      }
      else if ( ihave == 6 || ihave == 7 )
      {
        ihave = 8;
      }
      else if ( ihave == 9 )
      {
        ihave = 10;
      }

      ndig = ch_to_digit ( c );

      if ( ihave == 3 )
      {
        rtop = 10.0 * rtop + ( double ) ndig;
      }
      else if ( ihave == 5 )
      {
        rtop = 10.0 * rtop + ( double ) ndig;
        rbot = 10.0 * rbot;
      }
      else if ( ihave == 8 )
      {
        jtop = 10 * jtop + ndig;
      }
      else if ( ihave == 10 )
      {
        jtop = 10 * jtop + ndig;
        jbot = 10 * jbot;
      }
    }
/*
  Anything else is regarded as a terminator.
*/
    else
    {
      iterm = 1;
    }
/*
  If we haven't seen a terminator, and we haven't examined the
  entire string, go get the next character.
*/
    if ( iterm == 1 || nchar <= *lchar + 1 )
    {
      break;
    }

  }
/*
  If we haven't seen a terminator, and we have examined the
  entire string, then we're done, and LCHAR is equal to NCHAR.
*/
  if ( iterm != 1 && (*lchar) + 1 == nchar )
  {
    *lchar = nchar;
  }
/*
  Number seems to have terminated.  Have we got a legal number?
  Not if we terminated in states 1, 2, 6 or 7!
*/
  if ( ihave == 1 || ihave == 2 || ihave == 6 || ihave == 7 )
  {
    *error = 1;
    return r;
  }
/*
  Number seems OK.  Form it.

  We have had some trouble with input of the form 1.0E-312.
  For now, let's assume anything less than 1.0E-20 is zero.
*/
  if ( jtop == 0 )
  {
    rexp = 1.0;
  }
  else
  {
    if ( jbot == 1 )
    {
      if ( jsgn * jtop < -20 )
      {
        rexp = 0.0;
      }
      else
      {
        rexp = pow ( ( double ) 10.0, ( double ) ( jsgn * jtop ) );
      }
    }
    else
    {
      if ( jsgn * jtop < -20 * jbot )
      {
        rexp = 0.0;
      }
      else
      {
        rexp = jsgn * jtop;
        rexp = rexp / jbot;
        rexp = pow ( ( double ) 10.0, ( double ) rexp );
      }
    }
  }

  r = isgn * rexp * rtop / rbot;

  return r;
}
/******************************************************************************/

int s_to_r8vec ( char *s, int n, double rvec[] )

/******************************************************************************/
/*
  Purpose:

    S_TO_R8VEC reads an R8VEC from a string.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 February 2001

  Author:

    John Burkardt

  Parameters:

    Input, char *S, the string to be read.

    Input, int N, the number of values expected.

    Output, double RVEC[N], the values read from the string.

    Output, int S_TO_R8VEC, is TRUE (1) if an error occurred and FALSE (0)
    otherwise.
*/
{
  int error;
  int i;
  int lchar;

  error = 0;

  for ( i = 0; i < n; i++ )
  {
    rvec[i] = s_to_r8 ( s, &lchar, &error );

    if ( error )
    {
      return error;
    }

    s = s + lchar;
  }

  return error;
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

