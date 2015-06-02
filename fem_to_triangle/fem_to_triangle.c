# include <stdlib.h>
# include <stdio.h>
# include <time.h>
# include <math.h>
# include <string.h>

int main ( int argc, char *argv[] );
char ch_cap ( char ch );
int ch_eqi ( char ch1, char ch2 );
int ch_to_digit ( char ch );
int file_column_count ( char *filename );
int file_row_count ( char *filename );
int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
int *i4mat_data_read ( char *input_filename, int m, int n );
void i4mat_header_read ( char *input_filename, int *m, int *n );
int i4mat_max ( int m, int n, int a[] );
int i4mat_min ( int m, int n, int a[] );
void i4mat_transpose_print ( int m, int n, int a[], char *title );
void i4mat_transpose_print_some ( int m, int n, int a[], int ilo, int jlo, 
  int ihi, int jhi, char *title );
void mesh_base_one ( int node_num, int element_order, int element_num, 
  int element_node[] );
double *r8mat_data_read ( char *input_filename, int m, int n );
void r8mat_header_read ( char *input_filename, int *m, int *n );
void r8mat_transpose_print ( int m, int n, double a[], char *title );
void r8mat_transpose_print_some ( int m, int n, double a[], int ilo, int jlo, 
  int ihi, int jhi, char *title );
int s_len_trim ( char *s );
int s_to_i4 ( char *s, int *last, int *error );
int s_to_i4vec ( char *s, int n, int ivec[] );
double s_to_r8 ( char *s, int *lchar, int *error );
int s_to_r8vec ( char *s, int n, double rvec[] );
int s_word_count ( char *s );
void timestamp ( );
void triangle_element_write ( char *element_file, int element_num, 
  int element_order,int element_att_num, int element_node[], 
  double element_att[] );
void triangle_node_write ( char *node_file, int node_num, int node_dim,
  int node_att_num, int node_marker_num, double node_coord[],
  double node_att[], int node_marker[] );

/******************************************************************************/

int main ( int argc, char *argv[] )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for FEM_TO_TRIANGLE.

  Discussion:

    FEM_TO_TRIANGLE converts a mesh of triangles from FEM to TRIANGLE format.

  Usage:

    fem_to_triangle prefix

    where 'prefix' is the common filename prefix:

    * 'prefix'_nodes.txt contains the FEM node coordinates,
    * 'prefix'_elements.txt contains the FEM element connectivities.
    * 'prefix'.node will contains the TRIANGLE node coordinates,
    * 'prefix'.ele will contains the TRIANGLE element connectivities.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    14 October 2014

  Author:

    John Burkardt
*/
{
  int element;
  double *element_att;
  int element_att_num;
  int *element_node;
  int element_num;
  int element_order;
  char fem_element_filename[255];
  char fem_node_filename[255];
  int i;
  int m;
  int node;
  double *node_att;
  int node_att_num;
  int *node_marker;
  int node_marker_num;
  int node_num;
  double *node_x;
  char prefix[255];
  char triangle_element_filename[255];
  char triangle_node_filename[255];

  timestamp ( );
  printf ( "\n" );
  printf ( "FEM_TO_TRIANGLE\n" );
  printf ( "  C version:\n" );
  printf ( "  Convert a 2D mesh from FEM to TRIANGLE format.\n" );
  printf ( "\n" );
  printf ( "  Read:\n" );
  printf ( "  * \"prefix\"_nodes.txt, FEM node coordinates.\n" );
  printf ( "  * \"prefix\"_elements.txt, FEM element connectivities.\n" );
  printf ( "\n" );
  printf ( "  Create:\n" );
  printf ( "  * \"prefix\".node, TRIANGLE node coordinates.\n" );
  printf ( "  * \"prefix\".ele, TRIANGLE element connectivities.\n" );
/*
  Get the filename prefix.
*/
  if ( argc <= 1 ) 
  {
    printf ( "\n" );
    printf ( "FEM_TO_TRIANGLE:\n" );
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
  strcpy ( fem_node_filename, prefix );
  strcat ( fem_node_filename, "_nodes.txt" );
  strcpy ( fem_element_filename, prefix );
  strcat ( fem_element_filename, "_elements.txt" );
  strcpy ( triangle_node_filename, prefix );
  strcat ( triangle_node_filename, ".node" );
  strcpy ( triangle_element_filename, prefix );
  strcat ( triangle_element_filename, ".ele" );
/*
  Read the node data.
*/
  r8mat_header_read ( fem_node_filename, &m, &node_num );

  printf ( "\n" );
  printf ( "  Read the header of \"%s\".\n", fem_node_filename );
  printf ( "\n" );
  printf ( "  Spatial dimension = %d\n", m );
  printf ( "  Number of nodes  = %d\n", node_num );

  if ( m != 2 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "FEM_TO_TRIANGLE - Fatal error!\n" );
    fprintf ( stderr, "  Spatial dimension must be 2.\n" );
    exit ( 1 );
  }
  node_x = r8mat_data_read ( fem_node_filename, m, node_num );

  printf ( "\n" );
  printf ( "  Read the data in \"%s\".\n", fem_node_filename );

  r8mat_transpose_print_some ( m, node_num, node_x, 1, 1, m, 5, 
    "  Portion of node coordinate data:" );
/*
  Read the element data.
*/
  i4mat_header_read ( fem_element_filename, &element_order, &element_num );


  printf ( "\n" );
  printf ( "  Read the header of \"%s\".\n", fem_element_filename );
  printf ( "\n" );
  printf ( "  Element order = %d\n", element_order );
  printf ( "  Number of elements  = %d\n", element_num );

  if ( element_order != 3 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "FEM_TO_TRIANGLE - Fatal error!\n" );
    fprintf ( stderr, "  Element order must be 3.\n" );
    exit ( 1 );
  }

  element_node = i4mat_data_read ( fem_element_filename, element_order,
    element_num );

  printf ( "\n" );
  printf ( "  Read the data in \"%s\".\n", fem_element_filename );

  i4mat_transpose_print_some ( element_order, element_num, element_node, 
    1, 1, element_order, 10, "  Initial portion of element data:" );
/*
  Force 1-based indexing.
*/
  mesh_base_one ( node_num, element_order, element_num, element_node );
/*
  Write out the TRIANGLE version of the data.
*/
  element_att_num = 0;
  element_att = ( double * ) NULL;

  triangle_element_write ( triangle_element_filename, element_num, 
    element_order, element_att_num, element_node, element_att );

  printf ( "\n" );
  printf ( "  Created the TRIANGLE element file \"%s\".\n", triangle_element_filename );

  node_att_num = 0;
  node_att = ( double * ) NULL;
  node_marker_num = 0;
  node_marker = ( int * ) NULL;

  triangle_node_write ( triangle_node_filename, node_num, m,
    node_att_num, node_marker_num, node_x, node_att, node_marker );

  printf ( "  Created the TRIANGLE node file \"%s\".\n", triangle_node_filename );
/*
  Free memory.
*/
  free ( element_node );
  free ( node_x );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "FEM_TO_TRIANGLE:\n" );
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
# define LINE_MAX 255

  int column_num;
  char *error;
  FILE *input;
  int got_one;
  char line[LINE_MAX];
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
    error = fgets ( line, LINE_MAX, input );

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
      error = fgets ( line, LINE_MAX, input );

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

# undef LINE_MAX
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
# define LINE_MAX 255

  int bad_num;
  int comment_num;
  char *error;
  FILE *input;
  int i;
  char line[LINE_MAX];
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
    error = fgets ( line, LINE_MAX, input );

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

# undef LINE_MAX
}
/******************************************************************************/

void gmsh_mesh1d_write ( char *gmsh_filename, int m, int node_num, 
  double node_x[], int element_order, int element_num, int element_node[] )

/******************************************************************************/
/*
  Purpose:

    GMSH_MESH1D_WRITE writes 1D mesh data as a Gmsh mesh file.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    06 October 2014

  Author:

    John Burkardt

  Reference:

    Christophe Geuzaine, Jean-Francois Remacle,
    Gmsh: a three-dimensional finite element mesh generator with
    built-in pre- and post-processing facilities,
    International Journal for Numerical Methods in Engineering,
    Volume 79, Number 11, pages 1309-1331, 2009.

  Parameters:

    Input, char *GMSH_FILENAME, the name of the Gmsh file.

    Input, int M, the spatial dimension.

    Input, inte NODE_NUM, the number of nodes.

    Input, double NODE_X[M*NODE_NUM], the node coordinates.

    Input, int ELEMENT_ORDER, the order of the elements.

    Input, int ELEMENT_NUM, the number of elements.

    Input, int ELEMENT_NODE[ELEMENT_ORDER*ELEMENT_NUM], 
    the nodes that make up each element.
*/
{
  int element;
  int element_type;
  FILE *gmsh;
  int i;
  int node;
  int tag_num;
  int tag1;
/*
  Open the file.
*/
  gmsh = fopen ( gmsh_filename, "wt" );
/*
  Write the data.
*/
  fprintf ( gmsh, "$MeshFormat\n" );
  fprintf ( gmsh, "2.2 0 8\n" );
  fprintf ( gmsh, "$EndMeshFormat\n" );

  fprintf ( gmsh, "$Nodes\n" );
  fprintf ( gmsh, "%d\n", node_num );
  for ( node = 0; node < node_num; node++ )
  {
    fprintf ( gmsh, "  %d  %g  0.0  0.0\n", 
      node + 1, node_x[0+node*m] );
  }
  fprintf ( gmsh, "$EndNodes\n" );

  element_type = 1;

  tag_num = 2;
  tag1 = 0;
  fprintf ( gmsh, "$Elements\n" );
  fprintf ( gmsh, "%d\n", element_num );
  for ( element = 0; element < element_num; element++ )
  {
    fprintf ( gmsh, "  %d  %d  %d  %d  %d",
      element + 1, element_type, tag_num, tag1, element + 1 );
    for ( i = 0; i < element_order; i++ )
    {
      fprintf ( gmsh, "  %d", element_node[i+element*element_order] );
    }
    fprintf ( gmsh, "\n" );
  }
  fprintf ( gmsh, "$EndElements\n" );

  fclose ( gmsh );

  return;
}
/******************************************************************************/

void gmsh_mesh2d_write ( char *gmsh_filename, int m, int node_num, 
  double node_x[], int element_order, int element_num, int element_node[] )

/******************************************************************************/
/*
  Purpose:

    GMSH_MESH2D_WRITE writes 2D mesh data as a Gmsh mesh file.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    06 October 2014

  Author:

    John Burkardt

  Reference:

    Christophe Geuzaine, Jean-Francois Remacle,
    Gmsh: a three-dimensional finite element mesh generator with
    built-in pre- and post-processing facilities,
    International Journal for Numerical Methods in Engineering,
    Volume 79, Number 11, pages 1309-1331, 2009.

  Parameters:

    Input, char *GMSH_FILENAME, the name of the Gmsh file.

    Input, int M, the spatial dimension.

    Input, inte NODE_NUM, the number of nodes.

    Input, double NODE_X[M*NODE_NUM], the node coordinates.

    Input, int ELEMENT_ORDER, the order of the elements.

    Input, int ELEMENT_NUM, the number of elements.

    Input, int ELEMENT_NODE[ELEMENT_ORDER*ELEMENT_NUM], 
    the nodes that make up each element.
*/
{
  int element;
  int element_type;
  FILE *gmsh;
  int i;
  int node;
  int tag_num;
  int tag1;
/*
  Open the file.
*/
  gmsh = fopen ( gmsh_filename, "wt" );
/*
  Write the data.
*/
  fprintf ( gmsh, "$MeshFormat\n" );
  fprintf ( gmsh, "2.2 0 8\n" );
  fprintf ( gmsh, "$EndMeshFormat\n" );

  fprintf ( gmsh, "$Nodes\n" );
  fprintf ( gmsh, "%d\n", node_num );
  for ( node = 0; node < node_num; node++ )
  {
    fprintf ( gmsh, "  %d  %g  %g  0.0\n", 
      node + 1, node_x[0+node*2], node_x[1+node*2] );
  }
  fprintf ( gmsh, "$EndNodes\n" );

  if ( element_order == 3 )
  {
    element_type = 2;
  }
  else if ( element_order == 6 )
  {
    element_type = 9;
  }

  tag_num = 2;
  tag1 = 0;
  fprintf ( gmsh, "$Elements\n" );
  fprintf ( gmsh, "%d\n", element_num );
  for ( element = 0; element < element_num; element++ )
  {
    fprintf ( gmsh, "  %d  %d  %d  %d  %d",
      element + 1, element_type, tag_num, tag1, element + 1 );
    for ( i = 0; i < element_order; i++ )
    {
      fprintf ( gmsh, "  %d", element_node[i+element*element_order] );
    }
    fprintf ( gmsh, "\n" );
  }
  fprintf ( gmsh, "$EndElements\n" );

  fclose ( gmsh );

  return;
}
/******************************************************************************/

void gmsh_mesh3d_write ( char *gmsh_filename, int m, int node_num, 
  double node_x[], int element_order, int element_num, int element_node[] )

/******************************************************************************/
/*
  Purpose:

    GMSH_MESH3D_WRITE writes 3D mesh data as a Gmsh mesh file.

  Discussion:

    The node ordering for the 20 node element is not standard.

    Assuming the vertices are A, B, C and D, Gmsh uses the following ordering:

    1:    a
    2:        b
    3:            c
    4:                d
    5: (2*a  +b        )/3
    6: (  a+2*b        )/3
    7: (    2*b+  c    )/3
    8: (      b+2*c    )/3
    9: (  a    +2*c    )/3
   10: (2*a    +  c    )/3
   11: (2*a        +  d)/3
   12: (  a        +2*d)/3
   13: (     b     +2*d)/3
   14: (   2*b     +  d)/3
   15: (       +  c+2*d)/3
   16: (       +2*c+  d)/3
   17: (  a+  b+  c    )/3
   18: (  a+  b    +  d)/3
   19: (      b+  c+  d)/3
   20: (  a+      c+  d)/3

    Leo Rebholz used the following ordering:

    1:    a
    2:        b
    3:            c
    4:                d
    5: (2*a  +b        )/3
    6: (2*a    +  c    )/3
    7: (  a+2*b        )/3
    8: (  a    +2*c    )/3
    9: (  a+  b+  c    )/3
   10: (    2*b+  c    )/3
   11: (      b+2*c    )/3
   12: (2*a        +  d)/3
   13: (   2*b     +  d)/3
   14: (       +2*c+  d)/3
   15: (  a+  b    +  d)/3
   16: (      b+  c+  d)/3
   17: (  a+      c+  d)/3
   18: (  a        +2*d)/3
   19: (     b     +2*d)/3
   20: (       +  c+2*d)/3

    Since the only 20 node data we have is from Leo, we will assume that
    all 20 node input data is in Leo's format, and needs to be converted
    to the Gmsh convention.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    06 October 2014

  Author:

    John Burkardt

  Reference:

    Christophe Geuzaine, Jean-Francois Remacle,
    Gmsh: a three-dimensional finite element mesh generator with
    built-in pre- and post-processing facilities,
    International Journal for Numerical Methods in Engineering,
    Volume 79, Number 11, pages 1309-1331, 2009.

  Parameters:

    Input, char *GMSH_FILENAME, the name of the Gmsh file.

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
  FILE *gmsh;
  int i;
  int i2;
  int leo_to_gmsh[20] = {
     0,  1,  2,  3,  4, 
     6,  9, 10,  7,  5, 
    11, 17, 18, 12, 19, 
    13,  8, 14, 15, 16 };
  int node;
  int tag_num;
  int tag1;
/*
  Open the file.
*/
  gmsh = fopen ( gmsh_filename, "wt" );
/*
  Write the data.
*/
  fprintf ( gmsh, "$MeshFormat\n" );
  fprintf ( gmsh, "2.2 0 8\n" );
  fprintf ( gmsh, "$EndMeshFormat\n" );

  fprintf ( gmsh, "$Nodes\n" );
  fprintf ( gmsh, "%d\n", node_num );
  for ( node = 0; node < node_num; node++ )
  {
    fprintf ( gmsh, "%d  %g  %g  %g\n", 
      node + 1, node_x[0+node*3], node_x[1+node*3], node_x[2+node*3] );
  }
  fprintf ( gmsh, "$EndNodes\n" );

  if ( element_order == 4 )
  {
    element_type = 4;
  }
  else if ( element_order == 10 )
  {
    element_type = 11;
  }
  else if ( element_order == 20 )
  {
    element_type = 29;
  }

  tag_num = 2;
  tag1 = 0;
  fprintf ( gmsh, "$Elements\n" );
  fprintf ( gmsh, "%d\n", element_num );
  for ( element = 0; element < element_num; element++ )
  {
    fprintf ( gmsh, "%d  %d  %d  %d  %d",
      element + 1, element_type, tag_num, tag1, element + 1 );
    if ( element_order == 20 )
    {
      for ( i = 0; i < element_order; i++ )
      {
        i2 = leo_to_gmsh[i];
        fprintf ( gmsh, "  %d", element_node[i2+element*element_order] );
      }
      fprintf ( gmsh, "\n" );
    }
    else
    {
      for ( i = 0; i < element_order; i++ )
      {
        fprintf ( gmsh, "  %d", element_node[i+element*element_order] );
      }
      fprintf ( gmsh, "\n" );
    }
  }
  fprintf ( gmsh, "$EndElements\n" );

  fclose ( gmsh );

  return;
}
/******************************************************************************/

int i4_max ( int i1, int i2 )

/******************************************************************************/
/*
  Purpose:

    I4_MAX returns the maximum of two I4's.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    29 August 2006

  Author:

    John Burkardt

  Parameters:

    Input, int I1, I2, are two integers to be compared.

    Output, int I4_MAX, the larger of I1 and I2.
*/
{
  int value;

  if ( i2 < i1 )
  {
    value = i1;
  }
  else
  {
    value = i2;
  }
  return value;
}
/******************************************************************************/

int i4_min ( int i1, int i2 )

/******************************************************************************/
/*
  Purpose:

    I4_MIN returns the smaller of two I4's.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    29 August 2006

  Author:

    John Burkardt

  Parameters:

    Input, int I1, I2, two integers to be compared.

    Output, int I4_MIN, the smaller of I1 and I2.
*/
{
  int value;

  if ( i1 < i2 )
  {
    value = i1;
  }
  else
  {
    value = i2;
  }
  return value;
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
# define LINE_MAX 255

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
    got_string = fgets ( line, LINE_MAX, input );

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

# undef LINE_MAX
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

int i4mat_max ( int m, int n, int a[] )

/******************************************************************************/
/*
  Purpose:

    I4MAT_MAX returns the maximum of an I4MAT.

  Discussion:

    An I4MAT is an MxN array of I4's, stored by (I,J) -> [I+J*M].

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 June 2010

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows in A.

    Input, int N, the number of columns in A.

    Input, int A[M*N], the M by N matrix.

    Output, int I4MAT_MAX, the maximum entry of A.
*/
{
  int i;
  const int i4_huge = 2147483647;
  int j;
  int value;

  value = - i4_huge;

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      if ( value < a[i+j*m] )
      {
        value = a[i+j*m];
      }
    }
  }
  return value;
}
/******************************************************************************/

int i4mat_min ( int m, int n, int a[] )

/******************************************************************************/
/*
  Purpose:

    I4MAT_MIN returns the minimum of an I4MAT.

  Discussion:

    An I4MAT is an MxN array of I4's, stored by (I,J) -> [I+J*M].

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 June 2010

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows in A.

    Input, int N, the number of columns in A.

    Input, int A[M*N], the M by N matrix.

    Output, int I4MAT_MIN, the minimum entry of A.
*/
{
  int i;
  const int i4_huge = 2147483647;
  int j;
  int value;

  value = i4_huge;

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      if ( a[i+j*m] < value )
      {
        value = a[i+j*m];
      }
    }
  }
  return value;
}
/******************************************************************************/

void i4mat_transpose_print ( int m, int n, int a[], char *title )

/******************************************************************************/
/*
  Purpose:

    I4MAT_TRANSPOSE_PRINT prints an I4MAT, transposed.

  Discussion:

    An I4MAT is an MxN array of I4's, stored by (I,J) -> [I+J*M].

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 January 2005

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows in A.

    Input, int N, the number of columns in A.

    Input, int A[M*N], the M by N matrix.

    Input, char *TITLE, a title.
*/
{
  int i;
  int j;
  int jhi;
  int jlo;

  i4mat_transpose_print_some ( m, n, a, 1, 1, m, n, title );

  return;
}
/******************************************************************************/

void i4mat_transpose_print_some ( int m, int n, int a[], int ilo, int jlo,
  int ihi, int jhi, char *title )

/******************************************************************************/
/*
  Purpose:

    I4MAT_TRANSPOSE_PRINT_SOME prints some of an I4MAT, transposed.

  Discussion:

    An I4MAT is an MxN array of I4's, stored by (I,J) -> [I+J*M].

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 August 2010

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows of the matrix.
    M must be positive.

    Input, int N, the number of columns of the matrix.
    N must be positive.

    Input, int A[M*N], the matrix.

    Input, int ILO, JLO, IHI, JHI, designate the first row and
    column, and the last row and column to be printed.

    Input, char *TITLE, a title.
*/
{
# define INCX 10

  int i;
  int i2hi;
  int i2lo;
  int j;
  int j2hi;
  int j2lo;

  fprintf ( stdout, "\n" );
  fprintf ( stdout, "%s\n", title );

  if ( m <= 0 || n <= 0 )
  {
    fprintf ( stdout, "\n" );
    fprintf ( stdout, "  (None)\n" );
    return;
  }
/*
  Print the columns of the matrix, in strips of INCX.
*/
  for ( i2lo = ilo; i2lo <= ihi; i2lo = i2lo + INCX )
  {
    i2hi = i2lo + INCX - 1;
    i2hi = i4_min ( i2hi, m );
    i2hi = i4_min ( i2hi, ihi );

    fprintf ( stdout, "\n" );
/*
  For each row I in the current range...

  Write the header.
*/
    fprintf ( stdout, "  Row: " );
    for ( i = i2lo; i <= i2hi; i++ )
    {
      fprintf ( stdout, "%6d  ", i - 1 );
    }
    fprintf ( stdout, "\n" );
    fprintf ( stdout, "  Col\n" );
    fprintf ( stdout, "\n" );
/*
  Determine the range of the rows in this strip.
*/
    j2lo = i4_max ( jlo, 1 );
    j2hi = i4_min ( jhi, n );

    for ( j = j2lo; j <= j2hi; j++ )
    {
/*
  Print out (up to INCX) entries in column J, that lie in the current strip.
*/
      fprintf ( stdout, "%5d: ", j - 1 );
      for ( i = i2lo; i <= i2hi; i++ )
      {
        fprintf ( stdout, "%6d  ", a[i-1+(j-1)*m] );
      }
      fprintf ( stdout, "\n" );
    }
  }

  return;
# undef INCX
}
/******************************************************************************/

void mesh_base_one ( int node_num, int element_order, int element_num, 
  int element_node[] )

/******************************************************************************/
/*
  Purpose:

    MESH_BASE_ONE ensures that the element definition is one-based.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    02 October 2009

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
  int element;
  int node;
  int node_max;
  int node_min;
  int order;

  node_min = i4mat_min ( element_order, element_num, element_node );
  node_max = i4mat_max ( element_order, element_num, element_node );

  if ( node_min == 1 && node_max == node_num )
  {
    printf ( "\n" );
    printf ( "MESH_BASE_ONE:\n" );
    printf ( "  The element indexing appears to be 1-based!\n" );
    printf ( "  No conversion is necessary.\n" );

  }
  else if ( node_min == 0 && node_max == node_num - 1 )
  {
    printf ( "\n" );
    printf ( "MESH_BASE_ONE:\n" );
    printf ( "  The element indexing appears to be 0-based!\n" );
    printf ( "  This will be converted to 1-based.\n" );
    for ( element = 0; element < element_num; element++ )
    {
      for ( order = 0; order < element_order; order++ )
      {
        element_node[order+element*element_order] =
          element_node[order+element*element_order] + 1;
      }
    }
  }
  else
  {
    printf ( "\n" );
    printf ( "MESH_BASE_ONE - Warning!\n" );
    printf ( "  The element indexing is not of a recognized type.\n" );
    printf ( "  NODE_MIN = %d\n", node_min );
    printf ( "  NODE_MAX = %d\n", node_max );
    printf ( "  NODE_NUM = %d\n", node_num );
  }
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
# define LINE_MAX 255

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
    got_string = fgets ( line, LINE_MAX, input );

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

# undef LINE_MAX
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

void r8mat_transpose_print ( int m, int n, double a[], char *title )

/******************************************************************************/
/*
  Purpose:

    R8MAT_TRANSPOSE_PRINT prints an R8MAT, transposed.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    28 May 2008

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns.

    Input, double A[M*N], an M by N matrix to be printed.

    Input, char *TITLE, a title.
*/
{
  r8mat_transpose_print_some ( m, n, a, 1, 1, m, n, title );

  return;
}
/******************************************************************************/

void r8mat_transpose_print_some ( int m, int n, double a[], int ilo, int jlo,
  int ihi, int jhi, char *title )

/******************************************************************************/
/*
  Purpose:

    R8MAT_TRANSPOSE_PRINT_SOME prints some of an R8MAT, transposed.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    10 September 2013

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns.

    Input, double A[M*N], an M by N matrix to be printed.

    Input, int ILO, JLO, the first row and column to print.

    Input, int IHI, JHI, the last row and column to print.

    Input, char *TITLE, a title.
*/
{
# define INCX 5

  int i;
  int i2;
  int i2hi;
  int i2lo;
  int i2lo_hi;
  int i2lo_lo;
  int inc;
  int j;
  int j2hi;
  int j2lo;

  fprintf ( stdout, "\n" );
  fprintf ( stdout, "%s\n", title );

  if ( m <= 0 || n <= 0 )
  {
    fprintf ( stdout, "\n" );
    fprintf ( stdout, "  (None)\n" );
    return;
  }

  if ( ilo < 1 )
  {
    i2lo_lo = 1;
  }
  else
  {
    i2lo_lo = ilo;
  }

  if ( ihi < m )
  {
    i2lo_hi = m;
  }
  else
  {
    i2lo_hi = ihi;
  }

  for ( i2lo = i2lo_lo; i2lo <= i2lo_hi; i2lo = i2lo + INCX )
  {
    i2hi = i2lo + INCX - 1;

    if ( m < i2hi )
    {
      i2hi = m;
    }
    if ( ihi < i2hi )
    {
      i2hi = ihi;
    }

    inc = i2hi + 1 - i2lo;

    fprintf ( stdout, "\n" );
    fprintf ( stdout, "  Row:" );
    for ( i = i2lo; i <= i2hi; i++ )
    {
      fprintf ( stdout, "  %7d     ", i - 1 );
    }
    fprintf ( stdout, "\n" );
    fprintf ( stdout, "  Col\n" );
    fprintf ( stdout, "\n" );

    if ( jlo < 1 )
    {
      j2lo = 1;
    }
    else
    {
      j2lo = jlo;
    }
    if ( n < jhi )
    {
      j2hi = n;
    }
    else
    {
      j2hi = jhi;
    }
    for ( j = j2lo; j <= j2hi; j++ )
    {
      fprintf ( stdout, "%5d:", j - 1 );
      for ( i2 = 1; i2 <= inc; i2++ )
      {
        i = i2lo - 1 + i2;
        fprintf ( stdout, "  %14g", a[(i-1)+(j-1)*m] );
      }
      fprintf ( stdout, "\n" );
    }
  }

  return;
# undef INCX
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
/******************************************************************************/

void triangle_element_write ( char *element_file, int element_num, 
  int element_order, int element_att_num, int element_node[], 
  double element_att[] )

/******************************************************************************/
/*
  Purpose:

    TRIANGLE_ELEMENT_WRITE writes a TRIANGLE ".ele" file.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 November 2010

  Author:

    John Burkardt

  Parameters:

    Input, char *ELEMENT_FILE, the name of the file to be written.

    Input, int ELEMENT_NUM, the number of elements.

    Input, int ELEMENT_ORDER, the order of the elements.

    Input, int ELEMENT_ATT_NUM, the number of element attributes.

    Input, int ELEMENT_NODE[ELEMENT_ORDER*ELEMENT_NUM], the indices of the
    nodes that make up each element.

    Input, double ELEMENT_ATT[ELEMENT_ATT_NUM*ELEMENT_NUM], the attributes
    of each element.
*/
{
  int att;
  int dim;
  int element;
  int order;
  FILE *output;

  output = fopen ( element_file, "w" );

  fprintf ( output, "%d  %d  %d\n",
    element_num, element_order, element_att_num );

  for ( element = 0; element < element_num; element++ )
  {
    fprintf ( output, "%d", element + 1 );
    for ( order = 0; order < element_order; order++ )
    {
      fprintf ( output, "  %d", element_node[order+element*element_order] );
    }
    for ( att = 0; att < element_att_num; att++ )
    {
      fprintf ( output, "  %f", element_att[att+element*element_att_num] );
    }
    fprintf ( output, "\n" );
  }

  fclose ( output );

  return;
}
/******************************************************************************/

void triangle_node_write ( char *node_file, int node_num, int node_dim,
  int node_att_num, int node_marker_num, double node_coord[],
  double node_att[], int node_marker[] )

/******************************************************************************/
/*
  Purpose:

    TRIANGLE_NODE_WRITE writes a TRIANGLE ".node" file.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 November 2010

  Author:

    John Burkardt

  Parameters:

    Input, char *NODE_FILE, the name of the node file to be written.

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
  int att;
  int dim;
  int node;
  FILE *output;

  output = fopen ( node_file, "w" );

  fprintf ( output, "%d  %d  %d  %d\n",
    node_num, node_dim, node_att_num, node_marker_num );

  for ( node = 0; node < node_num; node++ )
  {
    fprintf ( output, "%d", node + 1 );
    for ( dim = 0; dim < node_dim; dim++ )
    {
      fprintf ( output, "  %f", node_coord[dim+node*node_dim] );
    }
    for ( att = 0; att < node_att_num; att++ )
    {
      fprintf ( output, "  %f", node_att[att+node*node_att_num] );
    }
    if ( node_marker_num == 1 )
    {
      fprintf ( output, "  %d", node_marker[node] );
    }
    fprintf ( output, "\n" );
  }

  fclose ( output );

  return;
}
