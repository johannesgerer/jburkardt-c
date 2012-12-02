# include <stdlib.h>
# include <stdio.h>
# include <time.h>
# include <math.h>
# include <string.h>

int main ( int argc, char *argv[] );
void bandwidth_mesh ( int element_order, int element_num, int element_node[],
  int *ml, int *mu, int *m );
int file_column_count ( char *input_filename );
int file_row_count ( char *input_filename );
int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
int *i4mat_data_read ( char *input_filename, int m, int n );
void i4mat_header_read ( char *input_filename, int *m, int *n );
void i4mat_transpose_print_some ( int m, int n, int a[], int ilo, int jlo, 
  int ihi, int jhi, char *title );
int s_len_trim ( char *s );
int s_to_i4 ( char *s, int *last, int *error );
int s_to_i4vec ( char *s, int n, int ivec[] );
int s_word_count ( char *s );
void timestamp ( );

/******************************************************************************/

int main ( int argc, char *argv[] )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for MESH_BANDWIDTH.

  Discussion:

    MESH_BANDWIDTH determines the geometric bandwidth of a mesh of elements.

    The user supplies an element file, listing the indices of the nodes that
    make up each element.

    The program computes the geometric bandwidth associated with the mesh.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 September 2012

  Author:

    John Burkardt

  Usage:

    mesh_bandwidth element_file
*/
{
  char element_file_name[81];
  int *element_node;
  int element_num;
  int element_order;
  int i;
  int j;
  int m;
  int ml;
  int mu;

  timestamp ( );

  printf ( "\n" );
  printf ( "MESH_BANDWIDTH\n" );
  printf ( "  C version\n" );
  printf ( "  Read a mesh file which defines\n" );
  printf ( "  a \"triangulation\" of a region in the plane,\n" );
  printf ( "  or a \"tetrahedronization\" of a region in space,\n" );
  printf ( "  or any division of a regino in ND space into elements,\n" );
  printf ( "  using a mesh of elements of uniform order.\n" );
  printf ( "\n" );
  printf ( "  Determine the geometric mesh bandwidth.\n" );
  printf ( "\n" );
  printf ( "    M = ML + 1 + MU.\n" );
  printf ( "\n" );
  printf ( "  which is the bandwidth of the vertex connectivity\n" );
  printf ( "  matrix.\n" );
  printf ( "\n" );
  printf ( "  Note that a matrix associated with variables defined\n" );
  printf ( "  at the  nodes could have a greater bandwidth than M,\n" );
  printf ( "  since you might have multiple variables at a vertex,\n" );
  printf ( "  or the variable might be a vector quantity,\n" );
  printf ( "  or physical effects might link two variables that are\n" );
  printf ( "  not associated with vertices that are connected.\n" );
/*
  If at least one command line argument, it's the element file.
*/
  if ( argc <= 1 ) 
  {
    printf ( "\n" );
    printf ( "MESH_BANDWIDTH:\n" );
    printf ( "  Please enter the name of the element file.\n" );

    scanf ( "%s", element_file_name );
  }
  else 
  {
    strcpy ( element_file_name, argv[1] );
  } 
/*
  Read the triangulation data.
*/
  i4mat_header_read ( element_file_name, &element_order, &element_num );

  printf ( "\n" );
  printf ( "  Read the header of \"%s\".\n", element_file_name );
  printf ( "\n" );
  printf ( "  Element order ELEMENT_ORDER =    %d\n", element_order );
  printf ( "  Number of element ELEMENT_NUM  = %d\n", element_num );

  element_node = i4mat_data_read ( element_file_name, element_order,
    element_num );

  printf ( "\n" );
  printf ( "  Read the data in \"%s\".\n", element_file_name );

  i4mat_transpose_print_some ( element_order, element_num, element_node, 1, 1, 
    element_order, 10, "  Portion of data read from file:" );
/*
  Compute the bandwidth.
*/
  bandwidth_mesh ( element_order, element_num, element_node, &ml, &mu, &m );

  printf ( "\n" );
  printf ( "  Lower bandwidth ML = %d\n", ml );
  printf ( "  Upper bandwidth MU = %d\n", mu );
  printf ( "  Total bandwidth M  = %d\n", m );
/*
  Free memory.
*/
  free ( element_node );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "MESH_BANDWIDTH\n" );
  printf ( "  Normal end of execution.\n" );

  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void bandwidth_mesh ( int element_order, int element_num, int element_node[],
  int *ml, int *mu, int *m )

/******************************************************************************/
/*
  Purpose:

    BANDWIDTH_MESH determines the bandwidth of the coefficient matrix.

  Discussion:

    The quantity computed here is the "geometric" bandwidth determined
    by the finite element mesh alone.

    If a single finite element variable is associated with each node
    of the mesh, and if the nodes and variables are numbered in the
    same way, then the geometric bandwidth is the same as the bandwidth
    of a typical finite element matrix.

    The bandwidth M is defined in terms of the lower and upper bandwidths:

      M = ML + 1 + MU

    where 

      ML = maximum distance from any diagonal entry to a nonzero
      entry in the same row, but earlier column,

      MU = maximum distance from any diagonal entry to a nonzero
      entry in the same row, but later column.

    Because the finite element node adjacency relationship is symmetric,
    we are guaranteed that ML = MU.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 September 2012

  Author:

    John Burkardt

  Parameters:

    Input, int ELEMENT_ORDER, the order of the elements.

    Input, int ELEMENT_NUM, the number of elements.

    Input,  ELEMENT_NODE[ELEMENT_ORDER*ELEMENT_NUM];
    ELEMENT_NODE(I,J) is the global index of local node I in element J.

    Output, int *ML, *MU, the lower and upper bandwidths of the matrix.

    Output, int *M, the bandwidth of the matrix.
*/
{
  int element;
  int global_i;
  int global_j;
  int local_i;
  int local_j;

  *ml = 0;
  *mu = 0;

  for ( element = 0; element < element_num; element++ )
  {
    for ( local_i = 0; local_i < element_order; local_i++ )
    {
      global_i = element_node[local_i+element*element_order];

      for ( local_j = 0; local_j < element_order; local_j++ )
      {
        global_j = element_node[local_j+element*element_order];

        *mu = i4_max ( *mu, global_j - global_i );
        *ml = i4_max ( *ml, global_i - global_j );
      }
    }
  }

  *m = *ml + 1 + *mu;

  return;
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
    fprintf ( stderr, "  Could not open the input file: \"%s\"\n", input_filename );
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
    fprintf ( stderr, "  Could not open the input file: \"%s\"\n", input_filename );
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
    fprintf ( stderr, "  Could not open the input file: \"%s\"\n", input_filename );
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

int s_len_trim ( char *s )

/******************************************************************************/
/*
  Purpose:

    S_LEN_TRIM returns the length of a string to the last nonblank.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    26 April 2003

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
    if ( *t != ' ' )
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
