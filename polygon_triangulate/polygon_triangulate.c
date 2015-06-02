# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <string.h>

# include "polygon_triangulate.h"

/******************************************************************************/

int between ( double xa, double ya, double xb, double yb, double xc, 
  double yc )

/******************************************************************************/
/*
  Purpose:

    BETWEEN is TRUE if vertex C is between vertices A and B.

  Discussion:

    The points must be (numerically) collinear.

    Given that condition, we take the greater of XA - XB and YA - YB
    as a "scale" and check where C's value lies.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 May 2014

  Author:

    Original C version by Joseph ORourke.
    This C version by John Burkardt.

  Reference:

    Joseph ORourke,
    Computational Geometry in C,
    Cambridge, 1998,
    ISBN: 0521649765,
    LC: QA448.D38.

  Parameters:

    Input, double XA, YA, XB, YB, XC, YC, the coordinates of 
    the vertices.

    Output, int BETWEEN, is TRUE if C is between A and B.
*/
{
  int value;
  double xmax;
  double xmin;
  double ymax;
  double ymin;

  if ( ! collinear ( xa, ya, xb, yb, xc, yc ) )
  {
    value = 0;
  }
  else if ( fabs ( ya - yb ) < fabs ( xa - xb ) )
  {
    xmax = r8_max ( xa, xb );
    xmin = r8_min ( xa, xb );
    value = ( xmin <= xc && xc <= xmax );
  }
  else
  {
    ymax = r8_max ( ya, yb );
    ymin = r8_min ( ya, yb );
    value = ( ymin <= yc && yc <= ymax );
  }

  return value;
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

int collinear ( double xa, double ya, double xb, double yb, double xc, 
  double yc )

/******************************************************************************/
/*
  Purpose:

    COLLINEAR returns a measure of collinearity for three points.

  Discussion:

    In order to deal with collinear points whose coordinates are not
    numerically exact, we compare the area of the largest square
    that can be created by the line segment between two of the points
    to (twice) the area of the triangle formed by the points.

    If the points are collinear, their triangle has zero area.
    If the points are close to collinear, then the area of this triangle
    will be small relative to the square of the longest segment.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 May 2014

  Author:

    Original C version by Joseph ORourke.
    This C version by John Burkardt.

  Reference:

    Joseph ORourke,
    Computational Geometry in C,
    Cambridge, 1998,
    ISBN: 0521649765,
    LC: QA448.D38.

  Parameters:

    Input, double XA, YA, XB, YB, XC, YC, the coordinates of 
    the vertices.

    Output, int COLLINEAR, is TRUE if the points are judged 
    to be collinear.
*/
{
  double area;
  const double r8_eps = 2.220446049250313E-016;
  double side_ab_sq;
  double side_bc_sq;
  double side_ca_sq;
  double side_max_sq;
  int value;

  area = 0.5 * ( 
      ( xb - xa ) * ( yc - ya ) 
    - ( xc - xa ) * ( yb - ya ) );

  side_ab_sq = pow ( xa - xb, 2 ) + pow ( ya - yb, 2 );
  side_bc_sq = pow ( xb - xc, 2 ) + pow ( yb - yc, 2 );
  side_ca_sq = pow ( xc - xa, 2 ) + pow ( yc - ya, 2 );

  side_max_sq = r8_max ( side_ab_sq, r8_max ( side_bc_sq, side_ca_sq ) );

  if ( side_max_sq <= r8_eps )
  {
    value = 1;
  }
  else if ( 2.0 * fabs ( area ) <= r8_eps * side_max_sq )
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

int diagonal ( int im1, int ip1, int n, int prev[], int next[], double x[], 
  double y[] )

/******************************************************************************/
/*
  Purpose:

    DIAGONAL: VERTEX(IM1) to VERTEX(IP1) is a proper internal diagonal.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 May 2014

  Author:

    Original C version by Joseph ORourke.
    This C version by John Burkardt.

  Reference:

    Joseph ORourke,
    Computational Geometry in C,
    Cambridge, 1998,
    ISBN: 0521649765,
    LC: QA448.D38.

  Parameters:

    Input, int IM1, IP1, the indices of two vertices.

    Input, int N, the number of vertices.

    Input, int PREV[N], the previous neighbor of each vertex.

    Input, int NEXT[N], the next neighbor of each vertex.

    Input, double X[N], Y[N], the coordinates of each vertex.

    Output, int DIAGONAL, the value of the test.
*/
{
  int value;
  int value1;
  int value2;
  int value3;

  value1 = in_cone ( im1, ip1, n, prev, next, x, y );
  value2 = in_cone ( ip1, im1, n, prev, next, x, y );
  value3 = diagonalie ( im1, ip1, n, next, x, y );

  value = ( value1 && value2 && value3 );

  return value;
}
/******************************************************************************/

int diagonalie ( int im1, int ip1, int n, int next[], double x[], double y[] )

/******************************************************************************/
/*
  Purpose:

    DIAGONALIE is true if VERTEX(IM1):VERTEX(IP1) is a proper diagonal.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 May 2014

  Author:

    Original C version by Joseph ORourke.
    This C version by John Burkardt.

  Reference:

    Joseph ORourke,
    Computational Geometry in C,
    Cambridge, 1998,
    ISBN: 0521649765,
    LC: QA448.D38.

  Parameters:

    Input, int IM1, IP1, the indices of two vertices.

    Input, int N, the number of vertices.

    Input, int NEXT[N], the next neighbor of each vertex.

    Input, double X[N], Y[N], the coordinates of each vertex.

    Output, int DIAGONALIE, the value of the test.
*/
{
  int first;
  int j;
  int jp1;
  int value;
  int value2;

  first = im1;
  j = first;
  jp1 = next[first];

  value = 1;
/*
  For each edge VERTEX(J):VERTEX(JP1) of the polygon:
*/
  while ( 1 )
  {
/*
  Skip any edge that includes vertex IM1 or IP1.
*/
    if ( j == im1 || j == ip1 || jp1 == im1 || jp1 == ip1 )
    {
    }
    else
    {
      value2 = intersect ( x[im1], y[im1], x[ip1], y[ip1], x[j], y[j], 
        x[jp1], y[jp1] );

      if ( value2 )
      {
        value = 0;
        break;
      }
    }
    j = jp1;
    jp1 = next[j];

    if ( j == first )
    {
      break;
    }
  }

  return value;
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

void i4mat_print ( int m, int n, int a[], char *title )

/******************************************************************************/
/*
  Purpose:

    I4MAT_PRINT prints an I4MAT.

  Discussion:

    An I4MAT is an MxN array of I4's, stored by (I,J) -> [I+J*M].

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    28 May 2008

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

  i4mat_print_some ( m, n, a, 1, 1, m, n, title );

  return;
}
/******************************************************************************/

void i4mat_print_some ( int m, int n, int a[], int ilo, int jlo, int ihi,
  int jhi, char *title )

/******************************************************************************/
/*
  Purpose:

    I4MAT_PRINT_SOME prints some of an I4MAT.

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
  for ( j2lo = jlo; j2lo <= jhi; j2lo = j2lo + INCX )
  {
    j2hi = j2lo + INCX - 1;
    if ( n < j2hi )
    {
      j2hi = n;
    }
    if ( jhi < j2hi )
    {
      j2hi = jhi;
    }

    fprintf ( stdout, "\n" );
/*
  For each column J in the current range...

  Write the header.
*/
    fprintf ( stdout, "  Col:" );
    for ( j = j2lo; j <= j2hi; j++ )
    {
      fprintf ( stdout, "  %6d", j - 1 );
    }
    fprintf ( stdout, "\n" );
    fprintf ( stdout, "  Row\n" );
    fprintf ( stdout, "\n" );
/*
  Determine the range of the rows in this strip.
*/
    if ( 1 < ilo )
    {
      i2lo = ilo;
    }
    else
    {
      i2lo = 1;
    }
    if ( m < ihi )
    {
      i2hi = m;
    }
    else
    {
      i2hi = ihi;
    }

    for ( i = i2lo; i <= i2hi; i++ )
    {
/*
  Print out (up to INCX) entries in row I, that lie in the current strip.
*/
      fprintf ( stdout, "%5d:", i - 1 );
      for ( j = j2lo; j <= j2hi; j++ )
      {
        fprintf ( stdout, "  %6d", a[i-1+(j-1)*m] );
      }
      fprintf ( stdout, "\n" );
    }
  }

  return;
# undef INCX
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
    fprintf ( stderr, "  Could not open the output file.\n" );
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

int in_cone ( int im1, int ip1, int n, int prev[], int next[], double x[], 
  double y[] )

/******************************************************************************/
/*
  Purpose:

    IN_CONE is TRUE if the diagonal VERTEX(IM1):VERTEX(IP1) is strictly internal.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 May 2014

  Author:

    Original C version by Joseph ORourke.
    This C version by John Burkardt.

  Reference:

    Joseph ORourke,
    Computational Geometry in C,
    Cambridge, 1998,
    ISBN: 0521649765,
    LC: QA448.D38.

  Parameters:

    Input, int IM1, IP1, the indices of two vertices.

    Input, int N, the number of vertices.

    Input, int PREV[N], the previous neighbor of each vertex.

    Input, int NEXT[N], the next neighbor of each vertex.

    Input, double X[N], Y[N], the coordinates of each vertex.

    Output, int IN_CONE, the value of the test.
*/
{
  int i;
  int im2;
  double t1;
  double t2;
  double t3;
  double t4;
  double t5;
  int value;

  im2 = prev[im1];
  i = next[im1];

  t1 = triangle_area ( x[im1], y[im1], x[i], y[i], x[im2], y[im2] );

  if ( 0.0 <= t1 )
  {
    t2 = triangle_area ( x[im1], y[im1], x[ip1], y[ip1], x[im2], y[im2] );
    t3 = triangle_area ( x[ip1], y[ip1], x[im1], y[im1], x[i], y[i] );
    value = ( ( 0.0 < t2 ) && ( 0.0 < t3 ) );
  }
  else
  {
    t4 = triangle_area ( x[im1], y[im1], x[ip1], y[ip1], x[i], y[i] );
    t5 = triangle_area ( x[ip1], y[ip1], x[im1], y[im1], x[im2], y[im2] );
    value = ! ( ( 0.0 <= t4 ) && ( 0.0 <= t5 ) );
  }
  return value;
}
/******************************************************************************/

int intersect ( double xa, double ya, double xb, double yb, double xc, 
  double yc, double xd, double yd )

/******************************************************************************/
/*
  Purpose:

    INTERSECT is true if lines VA:VB and VC:VD intersect.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 May 2014

  Author:

    Original C version by Joseph ORourke.
    This C version by John Burkardt.

  Reference:

    Joseph ORourke,
    Computational Geometry in C,
    Cambridge, 1998,
    ISBN: 0521649765,
    LC: QA448.D38.

  Parameters:

    Input, double XA, YA, XB, YB, XC, YC, XD, YD, the X and Y 
    coordinates of the four vertices.

    Output, int INTERSECT, the value of the test.
*/
{
  int value;

  if ( intersect_prop ( xa, ya, xb, yb, xc, yc, xc, yd ) )
  {
    value = 1;
  }
  else if ( between ( xa, ya, xb, yb, xc, yc ) )
  {
    value = 1;
  }
  else if ( between ( xa, ya, xb, yb, xd, yd ) )
  {
    value = 1;
  }
  else if ( between ( xc, yc, xd, yd, xa, ya ) )
  {
    value = 1;
  }
  else if ( between ( xc, yc, xd, yd, xb, yb ) )
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

int intersect_prop ( double xa, double ya, double xb, double yb, double xc, 
  double yc, double xd, double yd )

/******************************************************************************/
/*
  Purpose:

    INTERSECT_PROP is TRUE if lines VA:VB and VC:VD have a proper intersection.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    04 May 2014

  Author:

    Original C version by Joseph ORourke.
    This C version by John Burkardt.

  Reference:

    Joseph ORourke,
    Computational Geometry in C,
    Cambridge, 1998,
    ISBN: 0521649765,
    LC: QA448.D38.

  Parameters:

    Input, double XA, YA, XB, YB, XC, YC, XD, YD, the X and Y 
    coordinates of the four vertices.

    Output, int INTERSECT_PROP, the result of the test.
*/
{
  double t1;
  double t2;
  double t3;
  double t4;
  int value;
  int value1;
  int value2;
  int value3;
  int value4;

  if ( collinear ( xa, ya, xb, yb, xc, yc ) )
  {
    value = 0;
  }
  else if ( collinear ( xa, ya, xb, yb, xd, yd ) )
  {
    value = 0;
  }
  else if ( collinear ( xc, yc, xd, yd, xa, ya ) )
  {
    value = 0;
  }
  else if ( collinear ( xc, yc, xd, yd, xb, yb ) )
  {
    value = 0;
  }
  else
  {
    t1 = triangle_area ( xa, ya, xb, yb, xc, yc );
    t2 = triangle_area ( xa, ya, xb, yb, xd, yd );
    t3 = triangle_area ( xc, yc, xd, yd, xa, ya );
    t4 = triangle_area ( xc, yc, xd, yd, xb, yb );

    value1 = ( 0.0 < t1 );
    value2 = ( 0.0 < t2 );
    value3 = ( 0.0 < t3 );
    value4 = ( 0.0 < t4 );

    value = ( l4_xor ( value1, value2 ) ) && ( l4_xor ( value3, value4 ) );
  }
  return value;
}
/******************************************************************************/

int l4_xor ( int l1, int l2 )

/******************************************************************************/
/*
  Purpose:

    L4_XOR returns the exclusive OR of two L4's.

  Discussion:

    An L4 is a logical value.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    04 May 2014

  Author:

   John Burkardt

  Parameters:

    Input, int L1, L2, two values whose exclusive OR is needed.

    Output, int L4_XOR, the exclusive OR of L1 and L2.
*/
{
  int value;
  int value1;
  int value2;

  value1 = (     l1   && ( ! l2 ) );
  value2 = ( ( ! l1 ) &&     l2   );

  value = ( value1 || value2 );

  return value;
}
/******************************************************************************/

void l4vec_print ( int n, int a[], char *title )

/******************************************************************************/
/*
  Purpose:

    L4VEC_PRINT prints an L4VEC.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 May 2014

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of components of the vector.

    Input, int A[N], the (logical) vector to be printed.

    Input, char *TITLE, a title.
*/
{
  int i;

  printf ( "\n" );
  printf ( "%s\n", title );
  printf ( "\n" );
  for ( i = 0; i < n; i++ ) 
  {
    if ( a[i] == 0 )
    {
      printf ( "  %8d: F\n", i );
    }
    else
    {
      printf ( "  %8d: T\n", i );
    }
  }

  return;
}
/******************************************************************************/

int *polygon_triangulate ( int n, double x[], double y[] )

/******************************************************************************/
/*
  Purpose:

    POLYGON_TRIANGULATE determines a triangulation of a polygon.

  Discussion:

    There are N-3 triangles in the triangulation.

    For the first N-2 triangles, the first edge listed is always an
    internal diagonal.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 May 2014

  Author:

    Original C version by Joseph ORourke.
    This C version by John Burkardt.

  Reference:

    Joseph ORourke,
    Computational Geometry in C,
    Cambridge, 1998,
    ISBN: 0521649765,
    LC: QA448.D38.

  Parameters:

    Input, int N, the number of vertices.

    Input, double X[N], Y[N], the coordinates of each vertex.

    Output, int TRIANGLES[3*(N-2)], the triangles of the 
    triangulation.
*/
{
  double area;
  int *ear;
  int first;
  int i;
  int i0;
  int i1;
  int i2;
  int i3;
  int i4;
  int *next;
  int node;
  int node_m1;
  int *prev;
  int triangle_num;
  int *triangles;
/*
  We must have at least 3 vertices.
*/
  if ( n < 3 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "POLYGON_TRIANGULATE - Fatal error!\n" );
    fprintf ( stderr, "  N < 3.\n" );
    exit ( 1 );
  }
/*
  Consecutive vertices cannot be equal.
*/
  node_m1 = n - 1;
  for ( node = 0; node < n; node++ )
  {
    if ( x[node_m1] == x[node] && y[node_m1] == y[node] )
    {
      fprintf ( stderr, "\n" );
      fprintf ( stderr, "POLYGON_TRIANGULATE - Fatal error!\n" );
      fprintf ( stderr, "  Two consecutive nodes are identical.\n" );
      exit ( 1 );
    }
    node_m1 = node;
  }
/*
  Area must be positive.
*/
  area = 0.0;
  for ( node = 0; node < n - 2; node++ )
  {
    area = area + 0.5 * 
    ( 
        ( x[node+1] - x[node] ) * ( y[node+2] - y[node] ) 
      - ( x[node+2] - x[node] ) * ( y[node+1] - y[node] )
    );
  }

  if ( area <= 0.0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "POLYGON_TRIANGULATE - Fatal error!\n" );
    fprintf ( stderr, "  Polygon has zero or negative area.\n" );
    exit ( 1 );
  }

  triangles = ( int * ) malloc ( 3 * ( n - 2 ) * sizeof ( int ) );
/*
  PREV and NEXT point to the previous and next nodes.
*/
  prev = ( int * ) malloc ( n * sizeof ( int ) );
  next = ( int * ) malloc ( n * sizeof ( int ) );

  i = 0;
  prev[i] = n - 1;
  next[i] = i + 1;

  for ( i = 1; i < n - 1; i++ )
  {
    prev[i] = i - 1;
    next[i] = i + 1;
  }

  i = n - 1;
  prev[i] = i - 1;
  next[i] = 0;
/*
  EAR indicates whether the node and its immediate neighbors form an ear
  that can be sliced off immediately.
*/
  ear = ( int * ) malloc ( n * sizeof ( int ) );
  for ( i = 0; i < n; i++ )
  {
    ear[i] = diagonal ( prev[i], next[i], n, prev, next, x, y );
  }

  triangle_num = 0;

  i2 = 0;

  while ( triangle_num < n - 3 )
  {
/*
  If I2 is an ear, gather information necessary to carry out
  the slicing operation and subsequent "healing".
*/
    if ( ear[i2] )
    {
      i3 = next[i2];
      i4 = next[i3];
      i1 = prev[i2];
      i0 = prev[i1];
/*
  Make vertex I2 disappear.
*/
      next[i1] = i3;
      prev[i3] = i1;
/*
  Update the earity of I1 and I3, because I2 disappeared.
*/
      ear[i1] = diagonal ( i0, i3, n, prev, next, x, y );
      ear[i3] = diagonal ( i1, i4, n, prev, next, x, y );
/*
  Add the diagonal [I3, I1, I2] to the list.
*/
      triangles[0+triangle_num*3] = i3;
      triangles[1+triangle_num*3] = i1;
      triangles[2+triangle_num*3] = i2;
      triangle_num = triangle_num + 1;
    }
/*
  Try the next vertex.
*/
    i2 = next[i2];
  }
/*
  The last triangle is formed from the three remaining vertices.
*/
  i3 = next[i2];
  i1 = prev[i2];

  triangles[0+triangle_num*3] = i3;
  triangles[1+triangle_num*3] = i1;
  triangles[2+triangle_num*3] = i2;
  triangle_num = triangle_num + 1;

  free ( ear );
  free ( next );
  free ( prev );

  return triangles;
}
/******************************************************************************/

double r8_max ( double x, double y )

/******************************************************************************/
/*
  Purpose:

    R8_MAX returns the maximum of two R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 May 2006

  Author:

    John Burkardt

  Parameters:

    Input, double X, Y, the quantities to compare.

    Output, double R8_MAX, the maximum of X and Y.
*/
{
  double value;

  if ( y < x )
  {
    value = x;
  }
  else
  {
    value = y;
  }
  return value;
}
/******************************************************************************/

double r8_min ( double x, double y )

/******************************************************************************/
/*
  Purpose:

    R8_MIN returns the minimum of two R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 May 2006

  Author:

    John Burkardt

  Parameters:

    Input, double X, Y, the quantities to compare.

    Output, double R8_MIN, the minimum of X and Y.
*/
{
  double value;

  if ( y < x )
  {
    value = y;
  }
  else
  {
    value = x;
  }
  return value;
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

double triangle_area ( double xa, double ya, double xb, double yb, double xc, 
  double yc )

/******************************************************************************/
/*
  Purpose:

    TRIANGLE_AREA computes the signed area of a triangle.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 May 2014

  Author:

    John Burkardt

  Parameters:

    Input, double XA, YA, XB, YB, XC, YC, the coordinates of
    the vertices of the triangle, given in counterclockwise order.

    Output, double TRIANGLE_AREA, the signed area of the triangle.
*/
{
  double value;

  value = 0.5 * ( 
      ( xb - xa ) * ( yc - ya ) 
    - ( xc - xa ) * ( yb - ya ) );

  return value;
}
