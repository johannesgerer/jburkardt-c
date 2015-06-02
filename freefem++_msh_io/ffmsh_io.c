# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <string.h>

# include "ffmsh_io.h"

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

void ffmsh_2d_data_example ( int v_num, int e_num, int t_num, double v_xy[], 
  int v_l[], int e_v[], int e_l[], int t_v[], int t_l[] )

/******************************************************************************/
/*
  Purpose:

    FFMSH_2D_DATA_EXAMPLE returns example FFMSH data.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    21 December 2014

  Author:

    John Burkardt

  Parameters:

    Input, int V_NUM, the number of vertices.

    Input, int E_NUM, the number of boundary edges.

    Input, int T_NUM, the number of triangles.

    Output, double V_XY[2*V_NUM], vertex coordinates.

    Output, int V_L[V_NUM], vertex labels.

    Output, int E_V[2*E_NUM], edge vertices.

    Output, int E_L[E_NUM], vertex labels.

    Output, int T_V[3*T_NUM], triangle vertices.

    Output, int T_L[T_NUM], triangle labels.
*/
{
  int e_l_save[10] = {
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
  int e_v_save[2*10] = {
  11,  6, 
   6,  4, 
   4,  1, 
   1,  2, 
   2,  5, 
   5,  9, 
   9, 13, 
  13, 15, 
  15, 14, 
  14, 11 };
  int t_l_save[18] = {
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
  int t_v_save[3*18] = {
     1,  3,  4, 
     7,  2,  5, 
     9,  7,  5, 
     8,  6,  4, 
    12,  8,  7, 
    12, 11,  8, 
     3,  1,  2, 
     7,  3,  2, 
     7,  8,  3, 
     4,  3,  8, 
     6,  8, 11, 
    12,  7, 10, 
    11, 12, 14, 
    10,  9, 13, 
    12, 10, 13, 
     7,  9, 10, 
    12, 13, 15, 
    14, 12, 15 };
  int v_l_save[15] = {
    1, 1, 0, 1, 1, 1, 0, 0, 1, 0, 1, 0, 1, 1, 1 };
  double v_xy_save[2*15] = {
    -0.309016994375,  0.951056516295, 
    -0.809016994375,  0.587785252292, 
    -0.321175165867,  0.475528256720, 
     0.309016994375,  0.951056516295, 
    -1.000000000000,  0.000000000000, 
     0.809016994375,  0.587785252292, 
    -0.333333334358,  0.000000000000, 
     0.237841829972,  0.293892623813, 
    -0.809016994375, -0.587785252292, 
    -0.321175165867, -0.475528259963, 
     1.000000000000,  0.000000000000, 
     0.206011327827, -0.391856835534, 
    -0.309016994375, -0.951056516295, 
     0.809016994375, -0.587785252292, 
     0.309016994375, -0.951056516295 };

  i4vec_copy (    v_num, v_l_save,  v_l );
  r8mat_copy ( 2, v_num, v_xy_save, v_xy );
  i4vec_copy (    e_num, e_l_save,  e_l );
  i4mat_copy ( 2, e_num, e_v_save,  e_v );
  i4vec_copy (    t_num, t_l_save,  t_l );
  i4mat_copy ( 3, t_num, t_v_save,  t_v );

  return;
}
/******************************************************************************/

void ffmsh_2d_data_print ( char *title, int v_num, int e_num, int t_num, 
  double v_xy[], int v_l[], int e_v[], int e_l[], int t_v[], int t_l[] )

/******************************************************************************/
/*
  Purpose:

    FFMSH_2D_DATA_PRINT prints FFMSH data.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    21 December 2014

  Author:

    John Burkardt

  Parameters:

    Input, char *TITLE, a title.

    Input, int V_NUM, the number of vertices.

    Input, int E_NUM, the number of boundary edges.

    Input, int T_NUM, the number of triangles.

    Input, double V_XY[2*V_NUM], vertex coordinates.

    Input, int V_L[V_NUM], vertex labels.

    Input, int E_V[2*E_NUM], edge vertices.

    Input, int E_L[E_NUM], vertex labels.

    Input, int T_V[3*T_NUM], triangle vertices.

    Input, int T_L[T_NUM], triangle labels.
*/
{
  printf ( "\n" );
  printf ( "%s\n", title );
  
  i4vec_print (              v_num, v_l,  "  Vertex labels:" );
  r8mat_transpose_print ( 2, v_num, v_xy, "  Vertex coordinates:" );
  i4vec_print (              e_num, e_l,  "  Edge labels:" );
  i4mat_transpose_print ( 2, e_num, e_v,  "  Edge vertices:" );
  i4vec_print (              t_num, t_l,  "  Triangle labels:" );
  i4mat_transpose_print ( 3, t_num, t_v,  "  Triangle vertices:" );

  return;
}
/******************************************************************************/

void ffmsh_2d_data_read ( char *ffmsh_filename, int v_num, int e_num, int t_num, 
  double v_xy[], int v_l[], int e_v[], int e_l[], int t_v[], int t_l[] )

/******************************************************************************/
/*
  Purpose:

    FFMSH_2D_DATA_READ reads data from an FFMSH file.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    22 December 2014

  Author:

    John Burkardt

  Parameters:

    Input, char *FFMSH_FILENAME, the FFMSH filename.

    Input, int V_NUM, the number of vertices.

    Input, int E_NUM, the number of boundary edges.

    Input, int T_NUM, the number of triangles.

    Output, double V_XY[2*V_NUM], vertex coordinates.

    Output, int V_L[V_NUM], vertex labels.

    Output, int E_V[2*E_NUM], edge vertices.

    Output, int E_L[E_NUM], vertex labels.!

    Output, int T_V[3*T_NUM], triangle vertices.

    Output, int T_L[T_NUM], triangle labels.
*/
{
  int e_num2;
  char *error;
  FILE *ffmsh_unit;
  int i1;
  int i2;
  int i3;
  int i4;
  int ierror;
  int j;
  int length;
  double r1;
  double r2;
  int t_num2;
  char text[255];
  char *text_pointer;
  int v_num2;

  ffmsh_unit = fopen ( ffmsh_filename, "rt" );

  if ( ! ffmsh_unit )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "FFMSH_2D_DATA_READ - Fatal error!\n" );
    fprintf ( stderr, "  Could not open input file \"%s\"\n", ffmsh_filename );
    exit ( 1 );
  }
/*
  Read the sizes (again).
*/
  text_pointer = text;
  error = fgets ( text_pointer, 255, ffmsh_unit );

  v_num2 = s_to_i4 ( text_pointer, &length, &ierror );
  text_pointer = text_pointer + length;
  t_num2 = s_to_i4 ( text_pointer, &length, &ierror );
  text_pointer = text_pointer + length;
  e_num2 = s_to_i4 ( text_pointer, &length, &ierror );
  text_pointer = text_pointer + length;
/*
  Read Vertex X, Y, Label
*/
  for ( j = 0; j < v_num; j++ )
  {
    text_pointer = text;
    error = fgets ( text_pointer, 255, ffmsh_unit );

    r1 = s_to_r8 ( text_pointer, &length, &ierror );
    text_pointer = text_pointer + length;
    r2 = s_to_r8 ( text_pointer, &length, &ierror );
    text_pointer = text_pointer + length;
    i1 = s_to_i4 ( text_pointer, &length, &ierror );
    text_pointer = text_pointer + length;
    v_xy[0+j*2] = r1;
    v_xy[1+j*2] = r2;
    v_l[j] = i1;
  }
/*
  Read Triangle V1, V2, V3, Label
*/
  for ( j = 0; j < t_num; j++ )
  {
    text_pointer = text;
    error = fgets ( text_pointer, 255, ffmsh_unit );

    i1 = s_to_i4 ( text_pointer, &length, &ierror );
    text_pointer = text_pointer + length;
    i2 = s_to_i4 ( text_pointer, &length, &ierror );
    text_pointer = text_pointer + length;
    i3 = s_to_i4 ( text_pointer, &length, &ierror );
    text_pointer = text_pointer + length;
    i4 = s_to_i4 ( text_pointer, &length, &ierror );
    text_pointer = text_pointer + length;
    t_v[0+j*3] = i1;
    t_v[1+j*3] = i2;
    t_v[2+j*3] = i3;
    t_l[j] = i4;
  }
/*
  Read Edge V1, V2, Label
*/
  for ( j = 0; j < e_num; j++ )
  {
    text_pointer = text;
    error = fgets ( text_pointer, 255, ffmsh_unit );

    i1 = s_to_i4 ( text_pointer, &length, &ierror );
    text_pointer = text_pointer + length;
    i2 = s_to_i4 ( text_pointer, &length, &ierror );
    text_pointer = text_pointer + length;
    i3 = s_to_i4 ( text_pointer, &length, &ierror );
    text_pointer = text_pointer + length;
    e_v[0+j*2] = i1;
    e_v[1+j*2] = i2;
    e_l[j] = i3;
  }

  fclose ( ffmsh_unit );

  return;
}
/******************************************************************************/

void ffmsh_2d_size_example ( int *v_num, int *e_num, int *t_num )

/******************************************************************************/
/*
  Purpose:

    FFMSH_2D_SIZE_EXAMPLE returns sizes for the 2D example.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    21 December 2014

  Author:

    John Burkardt

  Parameters:

    Output, int *V_NUM, the number of vertices.

    Output, int *E_NUM, the number of boundary edges.

    Output, int *T_NUM, the number of triangles.
*/
{
  *e_num = 10;
  *t_num = 18;
  *v_num = 15;

  return;
}
/******************************************************************************/

void ffmsh_2d_size_print ( char *title, int v_num, int e_num, int t_num )

/******************************************************************************/
/*
  Purpose:

    FFMSH_2D_SIZE_PRINT prints the sizes of an FFMSH.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    21 December 2014

  Author:

    John Burkardt

  Parameters:

    Input, char *TITLE, a title.

    Input, int V_NUM, the number of vertices.

    Input, int E_NUM, the number of boundary edges.

    Input, int T_NUM, the number of triangles.
*/
{
  printf ( "\n" );
  printf ( "%s\n", title );
  printf ( "\n" );
  printf ( "  Number of vertices = %d\n", v_num );
  printf ( "  Number of boundary edges = %d\n", e_num );
  printf ( "  Number of triangles = %d\n", t_num );

  return;
}
/******************************************************************************/

void ffmsh_2d_size_read ( char *ffmsh_filename, int *v_num, int *e_num, 
  int *t_num )

/******************************************************************************/
/*
  Purpose:

    FFMSH_2D_SIZE_READ reads sizes from a FFMSH file.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    22 December 2014

  Author:

    John Burkardt

  Parameters:

    Input, char *FFMSH_FILENAME, the FFMSH filename.

    Output, int *V_NUM, the number of vertices.

    Output, int *E_NUM, the number of boundary edges.

    Output, int *T_NUM, the number of triangles.
*/
{
  char *error;
  FILE *ffmsh_unit;
  int ierror;
  int length;
  char text[255];
  char *text_pointer;

  ffmsh_unit = fopen ( ffmsh_filename, "rt" );

  if ( ! ffmsh_unit )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "FFMSH_SIZE_READ - Fatal error!\n" );
    fprintf ( stderr, "  Could not open the input file \"%s\"\n", ffmsh_filename );
    exit ( 1 );
  }

  text_pointer = text;
  error = fgets ( text_pointer, 255, ffmsh_unit );

  *v_num = s_to_i4 ( text_pointer, &length, &ierror );
  text_pointer = text_pointer + length;
  *t_num = s_to_i4 ( text_pointer, &length, &ierror );
  text_pointer = text_pointer + length;
  *e_num = s_to_i4 ( text_pointer, &length, &ierror );
  text_pointer = text_pointer + length;

  fclose ( ffmsh_unit );

  return;
}
/******************************************************************************/

void ffmsh_2d_write ( char *ffmsh_filename, int v_num, int e_num, int t_num, 
  double v_xy[], int v_l[], int e_v[], int e_l[], int t_v[], int t_l[] )

/******************************************************************************/
/*
  Purpose:

    FFMSH_2D_WRITE writes FFMSH data to a file.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    22 December 2014

  Author:

    John Burkardt

  Parameters:

    Input, char *FFMSH_FILENAME, the name of the file.

    Input, int V_NUM, the number of vertices.

    Input, int E_NUM, the number of boundary edges.

    Input, int T_NUM, the number of triangles.

    Input, double V_XY[2*V_NUM], vertex coordinates.

    Input, int V_L[V_NUM], vertex labels.

    Input, int E_V[2*E_NUM], edge vertices.

    Input, int E_L[E_NUM], vertex labels.

    Input, int T_V[3*T_NUM], triangle vertices.

    Input, int T_L[T_NUM], triangle labels.
*/
{
  FILE *ffmsh_unit;
  int j;
/*
  Open the file.
*/
  ffmsh_unit = fopen ( ffmsh_filename, "wt" );
/*
  Write the data.
*/
  fprintf ( ffmsh_unit, "%d  %d  %d\n", v_num, t_num, e_num );

  for ( j = 0; j < v_num; j++ )
  {
    fprintf ( ffmsh_unit, "%g  %g  %d\n", v_xy[0+j*2], v_xy[1+j*2], v_l[j] );
  }

  for ( j = 0; j < t_num; j++ )
  {
    fprintf ( ffmsh_unit, "%d  %d  %d  %d\n", t_v[0+j*3], t_v[1+j*3], t_v[2+j*3], t_l[j] );
  }

  for ( j = 0; j < e_num; j++ )
  {
    fprintf ( ffmsh_unit, "%d  %d  %d\n", e_v[0+j*2], e_v[1+j*2], e_l[j] );
  }

  fclose ( ffmsh_unit );

  return;
}
/******************************************************************************/

void i4mat_copy ( int m, int n, int a1[], int a2[] )

/******************************************************************************/
/*
  Purpose:

    I4MAT_COPY copies one I4MAT to another.

  Discussion:

    An I4MAT is an MxN array of I4's, stored by (I,J) -> [I+J*M].

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    27 August 2008

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns.

    Input, int A1[M*N], the matrix to be copied.

    Output, int A2[M*N], the copy of A1.
*/
{
  int i;
  int j;

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      a2[i+j*m] = a1[i+j*m];
    }
  }
  return;
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
    if ( m < i2hi )
    {
      i2hi = m;
    }
    if ( ihi < i2hi )
    {
      i2hi = ihi;
    }

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
    j2lo = jlo;
    if ( j2lo < 1 )
    {
      j2lo = 1;
    }
    j2hi = jhi;
    if ( n < jhi )
    {
      j2hi = n;
    }

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

void i4vec_copy ( int n, int a1[], int a2[] )

/******************************************************************************/
/*
  Purpose:

    I4VEC_COPY copies an I4VEC.

  Discussion:

    An I4VEC is a vector of I4's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    25 April 2007

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the vectors.

    Input, int A1[N], the vector to be copied.

    Input, int A2[N], the copy of A1.
*/
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    a2[i] = a1[i];
  }
  return;
}
/******************************************************************************/

void i4vec_print ( int n, int a[], char *title )

/******************************************************************************/
/*
  Purpose:

    I4VEC_PRINT prints an I4VEC.

  Discussion:

    An I4VEC is a vector of I4's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 November 2003

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of components of the vector.

    Input, int A[N], the vector to be printed.

    Input, char *TITLE, a title.
*/
{
  int i;

  fprintf ( stdout, "\n" );
  fprintf ( stdout, "%s\n", title );
  fprintf ( stdout, "\n" );

  for ( i = 0; i < n; i++ )
  {
    fprintf ( stdout, "  %6d: %8d\n", i, a[i] );
  }
  return;
}
/******************************************************************************/

void mesh_base_one ( int node_num, int element_order, int element_num,
  int element_node[] )

/******************************************************************************/
/*
  Purpose:

    MESH_BASE_ONE ensures that the element definition is one-based.

  Discussion:

    The ELEMENT_NODE array contains nodes indices that form elements.
    The convention for node indexing might start at 0 or at 1.

    If this function detects 0-based indexing, it converts to 1-based.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    18 October 2014

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
  const int i4_huge = 2147483647;
  int node;
  int node_max;
  int node_min;
  int order;

  node_min = + i4_huge;
  node_max = - i4_huge;
  for ( element = 0; element < element_num; element++ )
  {
    for ( order = 0; order < element_order; order++ )
    {
      node = element_node[order+element*element_order];
      if ( node < node_min )
      {
        node_min = node;
      }
      if ( node_max < node )
      {
        node_max = node;
      }
    }
  }

  if ( node_min == 0 && node_max == node_num - 1 )
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
  else if ( node_min == 1 && node_max == node_num )
  {
    printf ( "\n" );
    printf ( "MESH_BASE_ONE:\n" );
    printf ( "  The element indexing appears to be 1-based!\n" );
    printf ( "  No conversion is necessary.\n" );
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

void r8mat_copy ( int m, int n, double a1[], double a2[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_COPY copies one R8MAT to another.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 July 2008

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns.

    Input, double A1[M*N], the matrix to be copied.

    Output, double A2[M*N], the copy of A1.
*/
{
  int i;
  int j;

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      a2[i+j*m] = a1[i+j*m];
    }
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

  fprintf ( stdout, "%s\n", time_buffer );

  return;
# undef TIME_SIZE
}
