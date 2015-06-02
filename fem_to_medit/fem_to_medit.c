/*  Output from c_comment */

# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <string.h>

int main ( int argc, char *argv[] );
char ch_cap ( char c );
int ch_eqi ( char c1, char c2 );
int ch_to_digit ( char c );
int file_column_count ( char *input_filename );
int file_row_count ( char *input_filename );
int *i4mat_data_read ( char *input_filename, int m, int n );
void i4mat_header_read ( char *input_filename, int *m, int *n );
void medit_write ( char *filename, int dim, int vertices, int edges,
  int triangles, int quadrilaterals, int tetrahedrons, int hexahedrons,
  double vertex_coordinate[], int vertex_label[], int edge_vertex[],
  int edge_label[], int triangle_vertex[], int triangle_label[],
  int quadrilateral_vertex[], int quadrilateral_label[],
  int tetrahedron_vertex[], int tetrahedron_label[],
  int hexahedron_vertex[], int hexahedron_label[] );
double *r8mat_data_read ( char *input_filename, int m, int n );
void r8mat_header_read ( char *input_filename, int *m, int *n );
int s_len_trim ( char *s );
int s_to_i4 ( char *s, int *last, int *error );
int s_to_i4vec ( char *s, int n, int ivec[] );
double s_to_r8 ( char *s, int *lchar, int *error );
int s_to_r8vec ( char *s, int n, double rvec[] );
int s_word_count ( char *s );
void timestamp ( );

/******************************************************************************/

int main ( int argc, char *argv[] )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for FEM_TO_MEDIT.

  Discussion:

    FEM_TO_MEDIT converts mesh data from FEM format to MEDIT format.

    The FEM format defines "node", "element", and "boundary_node_mask",
    files for a mesh.  A typical set of such files might have the names
    "suv_nodes.txt", "suv_elements.txt" and "suv_boundary_node_mask.txt".

    This program reads these files and creates a MEDIT mesh file, whose
    name might be "suv.mesh".

  Usage:

    fem_to_medit prefix

    where 'prefix' is the common filename prefix so that:

    * prefix_nodes.txt contains the coordinates of the nodes;
    * prefix_elements.txt contains the indices of nodes forming each element;
    * prefix_boundary_node_mask.txt is 0 for interior nodes, 1 for boundary nodes;
    * prefix.mesh will be the MEDIT mesh file created by the program.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    10 October 2014

  Author:

    John Burkardt
*/
{
  int *boundary_node;
  char boundary_node_mask_filename[255];
  int dim;
  int *edge_label;
  int *edge_vertex;
  int edges;
  char element_filename[255];
  int *element_node;
  int element_num;
  int element_order;
  int *hexahedron_label;
  int *hexahedron_vertex;
  int hexahedrons;
  int i;
  int j;
  char medit_filename[255];
  double *node_coord;
  int node_dim;
  char node_filename[255];
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
  printf ( "FEM_TO_MEDIT:\n" );
  printf ( "  C version\n" );
  printf ( "  Read a set of FEM files.\n" );
  printf ( "  Write a corresponding MEDIT mesh file.\n" );
/*
  Get the filename prefix.
*/
  if ( argc <= 1 ) 
  {
    printf ( "\n" );
    printf ( "MESH_TO_XML:\n" );
    printf ( "  Please enter the filename prefix.\n" );

    scanf ( "%s", prefix );
  }
  else 
  {
    strcpy ( prefix, argv[1] );
  }
/*
  Create the file names.
*/
  strcpy ( node_filename, prefix );
  strcat ( node_filename, "_nodes.txt" );
  strcpy ( boundary_node_mask_filename, prefix );
  strcat ( boundary_node_mask_filename, "_boundary_node_mask.txt" );
  strcpy ( element_filename, prefix );
  strcat ( element_filename, "_elements.txt" );
  strcpy ( medit_filename, prefix );
  strcat ( medit_filename, ".mesh" );

  printf ( "\n" );
  printf ( "  Read FEM node file \"%s\"\n", node_filename );
  printf ( "    and FEM element file \"%s\"\n", element_filename );
  printf ( "    and FEM boundary node mask file \"%s\"\n", boundary_node_mask_filename );
  printf ( "  Create MEDIT file \"%s\"\n", medit_filename );
/*
  Read the FEM node data.
*/
  r8mat_header_read ( node_filename, &node_dim, &node_num );

  node_coord = r8mat_data_read ( node_filename, node_dim, node_num );

  printf ( "\n" );
  printf ( "  The node dimension is        %d\n", node_dim );
  printf ( "  The node number is           %d\n", node_num );
/*
  Read the FEM boundary node data.
*/
  boundary_node = i4mat_data_read ( boundary_node_mask_filename, 1, node_num );
/*
  Read the FEM element data.
*/
  i4mat_header_read ( element_filename, &element_order, &element_num );

  element_node = i4mat_data_read ( element_filename, element_order,
    element_num );

  printf ( "  The FEM element order is         %d\n", element_order );
  printf ( "  The FEM element number is        %d\n", element_num );
/*
  Write the MEDIT mesh data.
*/
  dim = node_dim;
  vertices = node_num;
  edges = 0;
  triangles = element_num;
  quadrilaterals = 0;
  tetrahedrons = 0;
  hexahedrons = 0;
  vertex_coordinate = ( double * ) malloc ( dim * vertices * sizeof ( double ) );
  for ( j = 0; j < vertices; j++ )
  {
    for ( i = 0; i < dim; i++ )
    {
      vertex_coordinate[i+j*dim] = node_coord[i+j*dim];
    }
  }
  vertex_label = ( int * ) malloc ( vertices * sizeof ( int ) );
  for ( j = 0; j < vertices; j++ )
  {
    vertex_label[j] = boundary_node[j];
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

  medit_write ( medit_filename, dim, vertices, edges, triangles, quadrilaterals,
    tetrahedrons, hexahedrons, vertex_coordinate, vertex_label, edge_vertex,
    edge_label, triangle_vertex, triangle_label, quadrilateral_vertex,
    quadrilateral_label, tetrahedron_vertex, tetrahedron_label,
    hexahedron_vertex, hexahedron_label );
/*
  Free memory.
*/
  free ( element_node );
  free ( node_coord );
  free ( boundary_node );
  free ( triangle_label );
  free ( triangle_vertex );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "FEM_TO_MEDIT:\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );
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

    CH_TO_DIGIT returns the integer value of a base 10 digit.

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

    Output, int CH_TO_DIGIT, the corresponding integer value.  If the 
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

int file_column_count ( char *input_filename )

/******************************************************************************/
/*
  Purpose:

    FILE_COLUMN_COUNT counts the number of columns in the first line of a file.

  Discussion:

    The file is assumed to be a simple text file.

    Most lines of the file is presumed to consist of COLUMN_NUM words, separated
    by spaces.  There may also be some blank lines, and some comment lines,
    which have a "#" in column 1.

    The routine tries to find the first non-comment non-blank line and
    counts the number of words in that line.

    If all lines are blanks or comments, it goes back and tries to analyze
    a comment line.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    13 June 2003

  Author:

    John Burkardt

  Parameters:

    Input, char *INPUT_FILENAME, the name of the file.

    Output, int FILE_COLUMN_COUNT, the number of columns assumed 
    to be in the file.
*/
{
# define MY_LINE_MAX 255

  int column_num;
  char *error;
  FILE *input;
  int got_one;
  char line[MY_LINE_MAX];
/*
  Open the file.
*/
  input = fopen ( input_filename, "r" );

  if ( !input )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "FILE_COLUMN_COUNT - Fatal error!\n" );
    fprintf ( stderr, "  Could not open the input file: \"%s\"\n",
      input_filename );
    exit ( 1 );
  }
/*
  Read one line, but skip blank lines and comment lines.
*/
  got_one = 0;

  for ( ; ; )
  {
    error = fgets ( line, MY_LINE_MAX, input );

    if ( !error )
    {
      break;
    }

    if ( s_len_trim ( line ) == 0 )
    {
      continue;
    }

    if ( line[0] == '#' )
    {
      continue;
    }

    got_one = 1;
    break;

  }

  if ( got_one == 0 )
  {
    fclose ( input );

    input = fopen ( input_filename, "r" );

    for ( ; ; )
    {
      error = fgets ( line, MY_LINE_MAX, input );

      if ( !error )
      {
        break;
      }

      if ( s_len_trim ( line ) == 0 )
      {
        continue;
      }

      got_one = 1;
      break;
    }
  }

  fclose ( input );

  if ( got_one == 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "FILE_COLUMN_COUNT - Warning!\n" );
    fprintf ( stderr, "  The file does not seem to contain any data.\n" );
    exit ( 1 );
  }

  column_num = s_word_count ( line );

  return column_num;

# undef MY_LINE_MAX
}
/******************************************************************************/

int file_row_count ( char *input_filename )

/******************************************************************************/
/*
  Purpose:

    FILE_ROW_COUNT counts the number of row records in a file.

  Discussion:

    It does not count lines that are blank, or that begin with a
    comment symbol '#'.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    13 June 2003

  Author:

    John Burkardt

  Parameters:

    Input, char *INPUT_FILENAME, the name of the input file.

    Output, int FILE_ROW_COUNT, the number of rows found.
*/
{
# define MY_LINE_MAX 255

  int bad_num;
  int comment_num;
  char *error;
  FILE *input;
  int i;
  char line[MY_LINE_MAX];
  int record_num;
  int row_num;

  row_num = 0;
  comment_num = 0;
  record_num = 0;
  bad_num = 0;

  input = fopen ( input_filename, "r" );

  if ( !input )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "FILE_ROW_COUNT - Fatal error!\n" );
    fprintf ( stderr, "  Could not open the input file: \"%s\"\n",
      input_filename );
    exit ( 1 );
  }

  for ( ; ; )
  {
    error = fgets ( line, MY_LINE_MAX, input );

    if ( !error )
    {
      break;
    }

    record_num = record_num + 1;

    if ( line[0] == '#' )
    {
      comment_num = comment_num + 1;
      continue;
    }

    if ( s_len_trim ( line ) == 0 )
    {
      comment_num = comment_num + 1;
      continue;
    }

    row_num = row_num + 1;
  }

  fclose ( input );

  return row_num;

# undef MY_LINE_MAX
}
/******************************************************************************/

int *i4mat_data_read ( char *input_filename, int m, int n )

/******************************************************************************/
/*
  Purpose:

    I4MAT_DATA_READ reads the data from an I4MAT file.

  Discussion:

    An I4MAT is an array of I4's.

    The file is assumed to contain one record per line.

    Records beginning with the '#' character are comments, and are ignored.
    Blank lines are also ignored.

    Each line that is not ignored is assumed to contain exactly (or at least)
    M real numbers, representing the coordinates of a point.

    There are assumed to be exactly (or at least) N such records.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    28 May 2008

  Author:

    John Burkardt

  Parameters:

    Input, char *INPUT_FILENAME, the name of the input file.

    Input, int M, the number of spatial dimensions.

    Input, int N, the number of points.  The program
    will stop reading data once N values have been read.

    Output, int I4MAT_DATA_READ[M*N], the data.
*/
{
# define MY_LINE_MAX 255

  int error;
  char *got_string;
  FILE *input;
  int i;
  int j;
  char line[255];
  int *table;
  int *x;

  input = fopen ( input_filename, "r" );

  if ( !input )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "I4MAT_DATA_READ - Fatal error!\n" );
    fprintf ( stderr, 
      "  Could not open the input file: \"%s\"\n", input_filename );
    exit ( 1 );
  }

  table = ( int * ) malloc ( m * n * sizeof ( int ) );

  x = ( int * ) malloc ( m * sizeof ( int ) );

  j = 0;

  while ( j < n )
  {
    got_string = fgets ( line, MY_LINE_MAX, input );

    if ( !got_string )
    {
      break;
    }

    if ( line[0] == '#' || s_len_trim ( line ) == 0 )
    {
      continue;
    }

    error = s_to_i4vec ( line, m, x );

    if ( error == 1 )
    {
      continue;
    }

    for ( i = 0; i < m; i++ )
    {
      table[i+j*m] = x[i];
    }
    j = j + 1;

  }

  fclose ( input );

  free ( x );

  return table;

# undef MY_LINE_MAX
}
/******************************************************************************/
 
void i4mat_header_read ( char *input_filename, int *m, int *n )
 
/******************************************************************************/
/*
  Purpose:

    I4MAT_HEADER_READ reads the header from an I4MAT file.

  Discussion:

    An I4MAT is an array of I4's.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    28 May 2008

  Author:

    John Burkardt

  Parameters:

    Input, char *INPUT_FILENAME, the name of the input file.

    Output, int *M, the number of spatial dimensions.

    Output, int *N, the number of points.
*/
{
  *m = file_column_count ( input_filename );

  if ( *m <= 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "I4MAT_HEADER_READ - Fatal error!\n" );
    fprintf ( stderr, "  FILE_COLUMN_COUNT failed.\n" );
    exit ( 1 );
  }

  *n = file_row_count ( input_filename );

  if ( *n <= 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "I4MAT_HEADER_READ - Fatal error!\n" );
    fprintf ( stderr, "  FILE_ROW_COUNT failed.\n" );
    exit ( 1 );
  }

  return;
}
/******************************************************************************/

void medit_write ( char *filename, int dim, int vertices, int edges,
  int triangles, int quadrilaterals, int tetrahedrons, int hexahedrons,
  double vertex_coordinate[], int vertex_label[], int edge_vertex[],
  int edge_label[], int triangle_vertex[], int triangle_label[],
  int quadrilateral_vertex[], int quadrilateral_label[],
  int tetrahedron_vertex[], int tetrahedron_label[],
  int hexahedron_vertex[], int hexahedron_label[] )

/******************************************************************************/
/*
  Purpose:

    MEDIT_WRITE writes mesh data to a MEDIT mesh file.

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
    fprintf ( stderr, "MEDIT_WRITE - Fatal error!\n" );
    fprintf ( stderr, "  Unable to open output file.\n" );
    exit ( 1 );
  }

  fprintf ( output, "MeshVersionFormatted 1\n" );
  fprintf ( output, "#  Created by medit_write.c\n" );
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

double *r8mat_data_read ( char *input_filename, int m, int n )

/******************************************************************************/
/*
  Purpose:

    R8MAT_DATA_READ reads the data from an R8MAT file.

  Discussion:

    An R8MAT is an array of R8's.

    The file is assumed to contain one record per line.

    Records beginning with the '#' character are comments, and are ignored.
    Blank lines are also ignored.

    Each line that is not ignored is assumed to contain exactly (or at least)
    M real numbers, representing the coordinates of a point.

    There are assumed to be exactly (or at least) N such records.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    27 January 2005

  Author:

    John Burkardt

  Parameters:

    Input, char *INPUT_FILENAME, the name of the input file.

    Input, int M, the number of spatial dimensions.

    Input, int N, the number of points.  The program
    will stop reading data once N values have been read.

    Output, double R8MAT_DATA_READ[M*N], the data.
*/
{
# define MY_LINE_MAX 255

  int error;
  char *got_string;
  FILE *input;
  int i;
  int j;
  char line[255];
  double *table;
  double *x;

  input = fopen ( input_filename, "r" );

  if ( !input )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8MAT_DATA_READ - Fatal error!\n" );
    fprintf ( stderr, "  Could not open the input file: \"%s\"\n", input_filename );
    exit ( 1 );
  }

  table = ( double * ) malloc ( m * n * sizeof ( double ) );

  x = ( double * ) malloc ( m * sizeof ( double ) );

  j = 0;

  while ( j < n )
  {
    got_string = fgets ( line, MY_LINE_MAX, input );

    if ( !got_string )
    {
      break;
    }

    if ( line[0] == '#' || s_len_trim ( line ) == 0 )
    {
      continue;
    }

    error = s_to_r8vec ( line, m, x );

    if ( error == 1 )
    {
      continue;
    }

    for ( i = 0; i < m; i++ )
    {
      table[i+j*m] = x[i];
    }
    j = j + 1;

  }

  fclose ( input );

  free ( x );

  return table;

# undef MY_LINE_MAX
}
/******************************************************************************/
 
void r8mat_header_read ( char *input_filename, int *m, int *n )
 
/******************************************************************************/
/*
  Purpose:

    R8MAT_HEADER_READ reads the header from an R8MAT file.

  Discussion:

    An R8MAT is an array of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    04 June 2004

  Author:

    John Burkardt

  Parameters:

    Input, char *INPUT_FILENAME, the name of the input file.

    Output, int *M, the number of spatial dimensions.

    Output, int *N, the number of points.
*/
{
  *m = file_column_count ( input_filename );

  if ( *m <= 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8MAT_HEADER_READ - Fatal error!\n" );
    fprintf ( stderr, "  FILE_COLUMN_COUNT failed.\n" );
    exit ( 1 );
  }

  *n = file_row_count ( input_filename );

  if ( *n <= 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8MAT_HEADER_READ - Fatal error!\n" );
    fprintf ( stderr, "  FILE_ROW_COUNT failed.\n" );
    exit ( 1 );
  }

  return;
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

int s_word_count ( char *s )

/******************************************************************************/
/*
  Purpose:

    S_WORD_COUNT counts the number of "words" in a string.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    30 January 2006

  Author:

    John Burkardt

  Parameters:

    Input, char *S, the string to be examined.

    Output, int S_WORD_COUNT, the number of "words" in the string.
    Words are presumed to be separated by one or more blanks.
*/
{
  int blank;
  int i;
  int word_num;

  word_num = 0;
  blank = 1;

  while ( *s ) 
  {
    if ( *s == ' ' || *s == '\n' )
    {
      blank = 1;
    }
    else if ( blank )
    {
      word_num = word_num + 1;
      blank = 0;
    }
    *s++;
  }
  return word_num;
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
  size_t len;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  len = strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  printf ( "%s\n", time_buffer );

  return;
# undef TIME_SIZE
}
