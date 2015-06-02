# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <string.h>
# include <time.h>

int main ( int argc, char *argv[] );
char ch_cap ( char ch );
int ch_eqi ( char ch1, char ch2 );
int ch_to_digit ( char ch );
void comp_next ( int n, int k, int a[], int *more, int *h, int *t );
int file_column_count ( char *filename );
int file_row_count ( char *filename );
double *monomial_value ( int dim_num, int point_num, int expon[], double x[] );
double r8_gamma ( double x );
double *r8mat_data_read ( char *input_filename, int m, int n );
void r8mat_header_read ( char *input_filename, int *m, int *n );
double r8vec_dot_product ( int n, double a1[], double a2[] );
double r8vec_sum ( int n, double a[] );
int s_eqi ( char *s1, char *s2 );
int s_len_trim ( char *s );
int s_to_i4 ( char *s, int *last, int *error );
double s_to_r8 ( char *s, int *lchar, int *error );
int s_to_r8vec ( char *s, int n, double rvec[] );
int s_word_count ( char *s );
double sphere01_monomial_integral ( int e[] );
double sphere01_monomial_quadrature ( int expon[], int point_num, double xyz[], 
  double w[] );
void timestamp ( );

/******************************************************************************/

int main ( int argc, char *argv[] )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for SPHERE_EXACTNESS.

  Discussion:

    This program investigates the polynomial exactness of a quadrature
    rule for the unit sphere

  Usage:

    sphere_exactness files prefix degree_max

    where

    * files explains how the quadrature rule is stored:
      'XYZW'  for file 'prefix.xyzw' containing (X,Y,Z,Weight);
      'RTPW'  for file 'prefix.rtpw' containing  (Theta, Phi, Weight) (radians);
      'DTPW'  for file 'prefix.dtpw' containing  (Theta, Phi, Weight) (degrees);
      'XYZ+W' for file 'prefix.xyz' containing (X,Y,Z)
              and file 'prefix.w' containing Weight;
      'RTP+W' for file 'prefix.rtp' containing (Theta, Phi ) in radians,
              and file 'prefix.w' containing Weight;
      'DTP+W' for file 'prefix.dtp' containing (Theta, Phi ) in degrees,
              and file 'prefix.w' containing Weight;
      'XYZ1'  for file 'prefix.xyz' containing (X,Y,Z), 
              and equal weights, which do not need to be read in.
      'RTP1'  for file 'prefix.rtp' containing (Theta, Phi ) in radians,
              and equal weights, which do not need to be read in.
      'DTP1'  for file 'prefix.dtp' containing (Theta, Phi ) in degrees,'
              and equal weights, which do not need to be read in.
    * prefix is the common file prefix;
    * degree_max is the maximum monomial degree to check.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    10 September 2010

  Author:

    John Burkardt
*/
{
  int degree;
  int degree_max;
  int dim;
  int dim_num;
  double *dtp;
  double *dtpw;
  int error;
  int *expon;
  char filename[255];
  char files[255];
  int h;
  int i;
  int j;
  int last;
  int more;
  const double pi = 3.141592653589793;
  int point_num;
  char prefix[255];
  double quad_error;
  double *rtp;
  double *rtpw;
  int t;
  double *w;
  double w_sum;
  double *xyz;
  double *xyzw;

  timestamp ( );
  printf ( "\n" );
  printf ( "SPHERE_EXACTNESS\n" );
  printf ( "  C version\n" );
  printf ( "\n" );
  printf ( "  Investigate the polynomial exactness of a quadrature\n" );
  printf ( "  rule for the unit sphere by integrating all monomials\n" );
  printf ( "  of a given degree.\n" );
//
//  Get the file structure;
//
  if ( 1 < argc ) 
  {
    strcpy ( files, argv[1] );
  }
  else
  {
    printf ( "\n" );
    printf ( "SPHERE_EXACTNESS:\n" );
    printf ( "  Describe the files to be read:\n" );
    printf ( "\n" );
    printf ( "  For coordinates and weights in one file:\n" );
    printf ( "    XYZW     (X,Y,Z,Weight)\n" );
    printf ( "    RTPW     (Theta, Phi, Weight) (radians)\n" );
    printf ( "    DTPW     (Theta, Phi, Weight) (degrees)\n" );
    printf ( "  For coordinates in one file and weights in another:\n" );
    printf ( "    XYZ+W    (X,Y,Z)       + Weight\n" );
    printf ( "    RTP+W    (Theta, Phi ) + Weight\n" );
    printf ( "    DTP+W    (Theta, Phi ) + Weight\n" );
    printf ( "  For coordinates in one file, and equal weights:\\n" );
    printf ( "    XYZ1     (X,Y,Z)\n" );
    printf ( "    RTP1     (Theta, Phi ) (radians)\n" );
    printf ( "    DTP1     (Theta, Phi ) (degrees)\n" );

    scanf ( "%s", files );
  }
//
//  Get the file prefix.
//
  if ( 2 < argc ) 
  {
    strcpy ( prefix, argv[2] );
  }
  else
  {
    printf ( "\n" );
    printf ( "SPHERE_EXACTNESS:\n" );
    printf ( "  Enter the filename prefix.\n" );

    scanf ( "%s", prefix );
  }
//
//  Get the maximum degree.
//
  if ( 3 < argc )
  {
    degree_max = s_to_i4 ( argv[3], &last, &error );
  }
  else
  {
    printf ( "\n" );
    printf ( "SPHERE_EXACTNESS:\n" );
    printf ( "  Please enter the maximum total degree to check.\n" );

    scanf ( "%d", degree_max );
  }
//
//  Summarize the input.
//
  printf ( "\n" );
  printf ( "SPHERE_EXACTNESS: User input:\n" );
  printf ( "  File structure = \"%s\".\n", files );
  printf ( "  Filename prefix = \"%s\".\n", prefix );
  printf ( "  Maximum total degree to check = %d\n", degree_max );
//
//  Read data needed to create XYZ and W arrays.
//
  if ( s_eqi ( files, "xyzw" ) )
  {
    strcpy ( filename, prefix );
    strcat ( filename, ".xyzw" );

    r8mat_header_read ( filename, &dim_num, &point_num );

    xyzw = r8mat_data_read ( filename, 4, point_num );

    xyz = ( double * ) malloc ( 3 * point_num * sizeof ( double ) );
    w = ( double * ) malloc ( point_num * sizeof ( double ) );

    for ( j = 0; j < point_num; j++ )
    {
      for ( i = 0; i < 3; i++ )
      {
         xyz[i+j*3] = xyzw[i+j*4];
      }
      w[j] = xyzw[3+j*4];
    }
    free ( xyzw );
  }
  else if ( s_eqi ( files, "xyz+w" ) )
  { 
    strcpy ( filename, prefix );
    strcat ( filename, ".xyz" );

    r8mat_header_read ( filename, &dim_num, &point_num );

    xyz = r8mat_data_read ( filename, 3, point_num );

    strcpy ( filename, prefix );
    strcat ( filename, ".w" );

    w = r8mat_data_read ( filename, 1, point_num );
  }
  else if ( s_eqi ( files, "xyz1" ) )
  {
    strcpy ( filename, prefix );
    strcat ( filename, ".xyz" );

    r8mat_header_read ( filename, &dim_num, &point_num );

    xyz = r8mat_data_read ( filename, 3, point_num );

    w = ( double * ) malloc ( point_num * sizeof ( double ) );
    for ( i = 0; i < point_num; i++ )
    {
      w[i] = 4.0 * pi / point_num;
    }
  }
  else if ( s_eqi ( files, "rtpw" ) )
  {
    strcpy ( filename, prefix );
    strcat ( filename, ".rtpw" );

    r8mat_header_read ( filename, &dim_num, &point_num );

    rtpw = r8mat_data_read ( filename, 3, point_num );
 
    xyz = ( double * ) malloc ( 3 * point_num * sizeof ( double ) );
    w = ( double * ) malloc ( point_num * sizeof ( double ) );

    for ( j = 0; j < point_num; j++ )
    {
      xyz[0+j*3] = cos ( rtpw[0+j*3] ) * sin ( rtpw[1+j*3] );
      xyz[1+j*3] = sin ( rtpw[0+j*3] ) * sin ( rtpw[1+j*3] );
      xyz[2+j*3] =                       cos ( rtpw[1+j*3] );
      w[j] = rtpw[2+j*3];
    }
    free ( rtpw );
  }
  else if ( s_eqi ( files, "rtp+w" ) )
  {
    strcpy ( filename, prefix );
    strcat ( filename, ".rtp" );

    r8mat_header_read ( filename, &dim_num, &point_num );

    rtp = r8mat_data_read ( filename, 2, point_num );

    strcpy ( filename, prefix );
    strcat ( filename, ".w" );

    w = r8mat_data_read ( filename, 1, point_num );

    xyz = ( double * ) malloc ( 3 * point_num * sizeof ( double ) );

    for ( j = 0; j < point_num; j++ )
    {
      xyz[0+j*3] = cos ( rtp[0+j*2] ) * sin ( rtp[1+j*2] );
      xyz[1+j*3] = sin ( rtp[0+j*2] ) * sin ( rtp[1+j*2] );
      xyz[2+j*3] =                      cos ( rtp[1+j*2] );
    }
    free ( rtp );
  }
  else if ( s_eqi ( files, "rtp1" ) )
  {
    strcpy ( filename, prefix );
    strcat ( filename, ".rtp" );

    r8mat_header_read ( filename, &dim_num, &point_num );

    rtp = r8mat_data_read ( filename, 2, point_num );

    xyz = ( double * ) malloc ( 3 * point_num * sizeof ( double ) );

    for ( j = 0; j < point_num; j++ )
    {
      xyz[0+j*3] = cos ( rtp[0+j*2] ) * sin ( rtp[1+j*2] );
      xyz[1+j*3] = sin ( rtp[0+j*2] ) * sin ( rtp[1+j*2] );
      xyz[2+j*3] =                      cos ( rtp[1+j*2] );
    }

    free ( rtp );

    w = ( double * ) malloc ( point_num * sizeof ( double ) );
    for ( i = 0; i < point_num; i++ )
    {
      w[i] = 4.0 * pi / point_num;
    }
  }
  else if ( s_eqi ( files, "dtpw" ) )
  { 
    strcpy ( filename, prefix );
    strcat ( filename, ".dtpw" );

    r8mat_header_read ( filename, &dim_num, &point_num );

    dtpw = r8mat_data_read ( filename, 3, point_num );

    xyz = ( double * ) malloc ( 3 * point_num * sizeof ( double ) );
    w = ( double * ) malloc ( point_num * sizeof ( double ) );

    for ( j = 0; j < point_num; j++ )
    {
      dtpw[0+j*3] = dtpw[0+j*3] * pi / 180.0;
      dtpw[1+j*3] = dtpw[1+j*3] * pi / 180.0;
    }
    for ( j = 0; j < point_num; j++ )
    {
      xyz[0+j*3] = cos ( dtpw[0+j*3] ) * sin ( dtpw[1+j*3] );
      xyz[1+j*3] = sin ( dtpw[0+j*3] ) * sin ( dtpw[1+j*3] );
      xyz[2+j*3] =                       cos ( dtpw[1+j*3] );
      w[j] = dtpw[2+j*3];
    }

    free ( dtpw );
  }
  else if ( s_eqi ( files, "dtp+w" ) )
  {
    strcpy ( filename, prefix );
    strcat ( filename, ".dtp" );

    r8mat_header_read ( filename, &dim_num, &point_num );

    dtp = r8mat_data_read ( filename, 2, point_num );

    strcpy ( filename, prefix );
    strcat ( filename, ".w" );

    w = r8mat_data_read ( filename, 1, point_num );

    xyz = ( double * ) malloc ( 3 * point_num * sizeof ( double ) );

    for ( j = 0; j < point_num; j++ )
    {
      dtp[0+j*2] = dtp[0+j*2] * pi / 180.0;
      dtp[1+j*2] = dtp[1+j*2] * pi / 180.0;
    }

    for ( j = 0; j < point_num; j++ )
    {
      xyz[0+j*3] = cos ( dtp[0+j*2] ) * sin ( dtp[1+j*2] );
      xyz[1+j*3] = sin ( dtp[0+j*2] ) * sin ( dtp[1+j*2] );
      xyz[2+j*3] =                      cos ( dtp[1+j*2] );
    }
    free ( dtp );
  }
  else if ( s_eqi ( files, "dtp1" ) )
  { 
    strcpy ( filename, prefix );
    strcat ( filename, ".dtp" );

    r8mat_header_read ( filename, &dim_num, &point_num );

    dtp = r8mat_data_read ( filename, 2, point_num );

    xyz = ( double * ) malloc ( 3 * point_num * sizeof ( double ) );

    for ( j = 0; j < point_num; j++ )
    {
      dtp[0+j*2] = dtp[0+j*2] * pi / 180.0;
      dtp[1+j*2] = dtp[1+j*2] * pi / 180.0;
    }

    for ( j = 0; j < point_num; j++ )
    {
      xyz[0+j*3] = cos ( dtp[0+j*2] ) * sin ( dtp[1+j*2] );
      xyz[1+j*3] = sin ( dtp[0+j*2] ) * sin ( dtp[1+j*2] );
      xyz[2+j*3] =                      cos ( dtp[1+j*2] );
    }
    free ( dtp );

    w = ( double * ) malloc ( point_num * sizeof ( double ) );
    for ( i = 0; i < point_num; i++ )
    {
      w[i] = 4.0 * pi / point_num;
    }
  }
  else
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "SPHERE_EXACTNESS - Fatal error!\n" );
    fprintf ( stderr, "  Unrecognized file structure choice\n" );
    exit ( 1 );
  }

  printf ( "\n" );
  printf ( "  Number of points  = %d\n", point_num );
//
//  The W's should sum to 4 * PI.
//
  w_sum = r8vec_sum ( point_num, w );

  for ( i = 0; i < point_num; i++ )
  {
    w[i] = 4.0 * pi * w[i] / w_sum;
  }
//
//  Explore the monomials.
//
  expon = ( int * ) malloc ( dim_num * sizeof ( int ) );

  printf ( "\n" );
  printf ( "      Error    Degree  Exponents\n" );

  for ( degree = 0; degree <= degree_max; degree++ )
  {
    printf ( "\n" );
    more = 0;
    h = 0;
    t = 0;

    for ( ; ; )
    {
      comp_next ( degree, dim_num, expon, &more, &h, &t );

      quad_error = sphere01_monomial_quadrature ( expon, point_num, xyz, w ); 

      printf ( "  %12g    %2d  ", quad_error, degree );

      for ( dim = 0; dim < dim_num; dim++ )
      {
        printf ( "%3d", expon[dim] );
      }
      printf ( "\n" );

      if ( !more )
      {
        break;
      }
    }

  }

  free ( expon );
  free ( w );
  free ( xyz );
//
//  Terminate.
//
  printf ( "\n" );
  printf ( "SPHERE_EXACTNESS:\n" );
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

void comp_next ( int n, int k, int a[], int *more, int *h,  int *t )

/******************************************************************************/
/*
  Purpose:

    COMP_NEXT computes the compositions of the integer N into K parts.

  Discussion:

    A composition of the integer N into K parts is an ordered sequence
    of K nonnegative integers which sum to N.  The compositions (1,2,1)
    and (1,1,2) are considered to be distinct.

    The routine computes one composition on each call until there are no more.
    For instance, one composition of 6 into 3 parts is
    3+2+1, another would be 6+0+0.

    On the first call to this routine, set MORE = FALSE.  The routine
    will compute the first element in the sequence of compositions, and
    return it, as well as setting MORE = TRUE.  If more compositions
    are desired, call again, and again.  Each time, the routine will
    return with a new composition.

    However, when the LAST composition in the sequence is computed 
    and returned, the routine will reset MORE to FALSE, signaling that
    the end of the sequence has been reached.

    This routine originally used a STATICE statement to maintain the
    variables H and T.  I have decided (based on an wasting an
    entire morning trying to track down a problem) that it is safer
    to pass these variables as arguments, even though the user should
    never alter them.  This allows this routine to safely shuffle
    between several ongoing calculations.


    There are 28 compositions of 6 into three parts.  This routine will
    produce those compositions in the following order:

     I         A
     -     ---------
     1     6   0   0
     2     5   1   0
     3     4   2   0
     4     3   3   0
     5     2   4   0
     6     1   5   0
     7     0   6   0
     8     5   0   1
     9     4   1   1
    10     3   2   1
    11     2   3   1
    12     1   4   1
    13     0   5   1
    14     4   0   2
    15     3   1   2
    16     2   2   2
    17     1   3   2
    18     0   4   2
    19     3   0   3
    20     2   1   3
    21     1   2   3
    22     0   3   3
    23     2   0   4
    24     1   1   4
    25     0   2   4
    26     1   0   5
    27     0   1   5
    28     0   0   6

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    02 July 2008

  Author:

    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
    C version by John Burkardt.

  Reference:

    Albert Nijenhuis, Herbert Wilf,
    Combinatorial Algorithms for Computers and Calculators,
    Second Edition,
    Academic Press, 1978,
    ISBN: 0-12-519260-6,
    LC: QA164.N54.

  Parameters:

    Input, int N, the integer whose compositions are desired.

    Input, int K, the number of parts in the composition.

    Input/output, int A[K], the parts of the composition.

    Input/output, int *MORE.
    Set MORE = FALSE on first call.  It will be reset to TRUE on return
    with a new composition.  Each new call returns another composition until
    MORE is set to FALSE when the last composition has been computed
    and returned.

    Input/output, int *H, *T, two internal parameters needed for the
    computation.  The user should allocate space for these in the calling
    program, include them in the calling sequence, but never alter them!
*/
{
  int i;

  if ( !( *more ) )
  {
    *t = n;
    *h = 0;
    a[0] = n;
    for ( i = 1; i < k; i++ )
    {
       a[i] = 0;
    }
  }
  else
  {
    if ( 1 < *t )
    {
      *h = 0;
    }
    *h = *h + 1;
    *t = a[*h-1];
    a[*h-1] = 0;
    a[0] = *t - 1;
    a[*h] = a[*h] + 1;
  }

  *more = ( a[k-1] != n );

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
# define MY_LINE_MAX 256

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
    column_num = -1;
    printf ( "\n" );
    printf ( "FILE_COLUMN_COUNT - Fatal error!\n" );
    printf ( "  Could not open the input file: \"%s\"\n", input_filename );
    return column_num;
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
    printf ( "\n" );
    printf ( "FILE_COLUMN_COUNT - Warning!\n" );
    printf ( "  The file does not seem to contain any data.\n" );
    return -1;
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
# define MY_LINE_MAX 256

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
    printf ( "\n" );
    printf ( "FILE_ROW_COUNT - Fatal error!\n" );
    printf ( "  Could not open the input file: \"%s\"\n", input_filename );
    return (-1);
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

double *monomial_value ( int m, int n, int expon[], double x[] )

/******************************************************************************/
/*
  Purpose:

    MONOMIAL_VALUE evaluates a monomial.

  Discussion:

    F(X) = product ( 1 <= I <= M ) X(I)^EXPON(I)

    with the convention that 0^0 = 1.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 April 2009

  Author:

    John Burkardt

  Parameters:

    Input, int M, the spatial dimension.

    Input, int N, the number of points.

    Input, int EXPON[M], the exponents.

    Input, double X[M*N], the evaluation points.

    Output, double MONOMIAL_VALUE[N], the monomial values.
*/
{
  int i;
  int j;
  double *v;

  v = ( double * ) malloc ( n * sizeof ( double ) );

  for ( j = 0; j < n; j++ )
  {
    v[j] = 1.0;
  }

  for ( i = 0; i < m; i++ )
  {
    if ( expon[i] != 0.0 )
    {
      for ( j = 0; j < n; j++ )
      {
        v[j] = v[j] * pow ( x[i+j*m], expon[i] );
      }
    }
  }

  return v;
}
/******************************************************************************/

double r8_gamma ( double x )

/******************************************************************************/
/*
  Purpose:

    R8_GAMMA evaluates Gamma(X) for a real argument.

  Discussion:

    The C math library includes the GAMMA ( X ) function which should generally
    be used instead of this function.

    This routine calculates the gamma function for a real argument X.

    Computation is based on an algorithm outlined in reference 1.
    The program uses rational functions that approximate the gamma
    function to at least 20 significant decimal digits.  Coefficients
    for the approximation over the interval (1,2) are unpublished.
    Those for the approximation for 12 <= X are from reference 2.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 January 2010

  Author:

    Original FORTRAN77 version by William Cody, Laura Stoltz.
    C version by John Burkardt.

  Reference:

    William Cody,
    An Overview of Software Development for Special Functions,
    in Numerical Analysis Dundee, 1975,
    edited by GA Watson,
    Lecture Notes in Mathematics 506,
    Springer, 1976.

    John Hart, Ward Cheney, Charles Lawson, Hans Maehly,
    Charles Mesztenyi, John Rice, Henry Thatcher,
    Christoph Witzgall,
    Computer Approximations,
    Wiley, 1968,
    LC: QA297.C64.

  Parameters:

    Input, double X, the argument of the function.

    Output, double R8_GAMMA, the value of the function.
*/
{
  double c[7] = {
   -1.910444077728E-03, 
    8.4171387781295E-04, 
   -5.952379913043012E-04, 
    7.93650793500350248E-04, 
   -2.777777777777681622553E-03, 
    8.333333333333333331554247E-02, 
    5.7083835261E-03 };
  double eps = 2.22E-16;
  double fact;
  int i;
  int n;
  double p[8] = {
  -1.71618513886549492533811E+00,
   2.47656508055759199108314E+01, 
  -3.79804256470945635097577E+02,
   6.29331155312818442661052E+02, 
   8.66966202790413211295064E+02,
  -3.14512729688483675254357E+04, 
  -3.61444134186911729807069E+04,
   6.64561438202405440627855E+04 };
  int parity;
  double q[8] = {
  -3.08402300119738975254353E+01,
   3.15350626979604161529144E+02, 
  -1.01515636749021914166146E+03,
  -3.10777167157231109440444E+03, 
   2.25381184209801510330112E+04,
   4.75584627752788110767815E+03, 
  -1.34659959864969306392456E+05,
  -1.15132259675553483497211E+05 };
  const double r8_pi = 3.1415926535897932384626434;
  double res;
  const double sqrtpi = 0.9189385332046727417803297;
  double sum;
  double value;
  double xbig = 171.624;
  double xden;
  double xinf = 1.79E+308;
  double xminin = 2.23E-308;
  double xnum;
  double y;
  double y1;
  double ysq;
  double z;

  parity = 0;
  fact = 1.0;
  n = 0;
  y = x;
/*
  Argument is negative.
*/
  if ( y <= 0.0 )
  {
    y = - x;
    y1 = ( double ) ( int ) ( y );
    res = y - y1;

    if ( res != 0.0 )
    {
      if ( y1 != ( double ) ( int ) ( y1 * 0.5 ) * 2.0 )
      {
        parity = 1;
      }

      fact = - r8_pi / sin ( r8_pi * res );
      y = y + 1.0;
    }
    else
    {
      res = xinf;
      value = res;
      return value;
    }
  }
/*
  Argument is positive.
*/
  if ( y < eps )
  {
/*
  Argument < EPS.
*/
    if ( xminin <= y )
    {
      res = 1.0 / y;
    }
    else
    {
      res = xinf;
      value = res;
      return value;
    }
  }
  else if ( y < 12.0 )
  {
    y1 = y;
/*
  0.0 < argument < 1.0.
*/
    if ( y < 1.0 )
    {
      z = y;
      y = y + 1.0;
    }
/*
  1.0 < argument < 12.0.
  Reduce argument if necessary.
*/
    else
    {
      n = ( int ) ( y ) - 1;
      y = y - ( double ) ( n );
      z = y - 1.0;
    }
/*
  Evaluate approximation for 1.0 < argument < 2.0.
*/
    xnum = 0.0;
    xden = 1.0;
    for ( i = 0; i < 8; i++ )
    {
      xnum = ( xnum + p[i] ) * z;
      xden = xden * z + q[i];
    }
    res = xnum / xden + 1.0;
/*
  Adjust result for case  0.0 < argument < 1.0.
*/
    if ( y1 < y )
    {
      res = res / y1;
    }
/*
  Adjust result for case 2.0 < argument < 12.0.
*/
    else if ( y < y1 )
    {
      for ( i = 1; i <= n; i++ )
      {
        res = res * y;
        y = y + 1.0;
      }
    }
  }
  else
  {
/*
  Evaluate for 12.0 <= argument.
*/
    if ( y <= xbig )
    {
      ysq = y * y;
      sum = c[6];
      for ( i = 0; i < 6; i++ )
      {
        sum = sum / ysq + c[i];
      }
      sum = sum / y - y + sqrtpi;
      sum = sum + ( y - 0.5 ) * log ( y );
      res = exp ( sum );
    }
    else
    {
      res = xinf;
      value = res;
      return value;
    }
  }
/*
  Final adjustments and return.
*/
  if ( parity )
  {
    res = - res;
  }

  if ( fact != 1.0 )
  {
    res = fact / res;
  }

  value = res;

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

double r8vec_dot_product ( int n, double a1[], double a2[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_DOT_PRODUCT computes the dot product of a pair of R8VEC's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 July 2007

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the vectors.

    Input, double A1[N], A2[N], the two vectors to be considered.

    Output, double R8VEC_DOT_PRODUCT, the dot product of the vectors.
*/
{
  int i;
  double value;

  value = 0.0;
  for ( i = 0; i < n; i++ )
  {
    value = value + a1[i] * a2[i];
  }
  return value;
}
/******************************************************************************/

double r8vec_sum ( int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_SUM returns the sum of an R8VEC.

  Discussion:

    An R8VEC is a vector of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 August 2008

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the vector.

    Input, double A[N], the vector.

    Output, double R8VEC_SUM, the sum of the vector.
*/
{
  int i;
  double value;

  value = 0.0;
  for ( i = 0; i < n; i++ )
  {
    value = value + a[i];
  }

  return value;
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

double sphere01_monomial_integral ( int e[] )

/******************************************************************************/
/*
  Purpose:

    SPHERE01_MONOMIAL_INT returns monomial integrals on the unit sphere.

  Discussion:

    The integration region is 

      X^2 + Y^2 + Z^2 = 1.

    The monomial is F(X,Y,Z) = X^E(1) * Y^E(2) * Z^E(3).

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    10 September 2010

  Author:

    John Burkardt

  Reference:

    Philip Davis, Philip Rabinowitz,
    Methods of Numerical Integration,
    Second Edition,
    Academic Press, 1984, page 263.

  Parameters:

    Input, int E[3], the exponents of X, Y and Z in the 
    monomial.  Each exponent must be nonnegative.

    Output, double SPHERE01_MONOMIAL_INTEGRAL, the integral.
*/
{
  int i;
  double integral;
  const double r8_pi = 3.141592653589793;

  if ( e[0] < 0 || e[1] < 0 || e[2] < 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "SPHERE01_MONOMIAL_INTEGRAL - Fatal error!\n" );
    fprintf ( stderr, "  All exponents must be nonnegative.\n" );
    exit ( 1 );
  }

  if ( e[0] == 0 && e[1] == 0 && e[2] == 0 )
  {
    integral = 2.0 * sqrt ( r8_pi * r8_pi * r8_pi ) / r8_gamma ( 1.5 );
  }
  else if ( ( e[0] % 2 == 1 ) || ( e[1] % 2 == 1 ) || ( e[2] % 2 == 1 ) )
  {
    integral = 0.0;
  }
  else
  {
    integral = 2.0;
    for ( i = 0; i < 3; i++ )
    {
      integral = integral * r8_gamma ( 0.5 * ( double ) ( e[i] + 1 ) );
    }
    integral = integral 
      / r8_gamma ( 0.5 * ( double ) ( e[0] + e[1] + e[2] + 3 ) );
  }
  return integral;
}
/******************************************************************************/

double sphere01_monomial_quadrature ( int expon[], int point_num, double xyz[], 
  double w[] )

/******************************************************************************/
/*
  Purpose:

    SPHERE01_MONOMIAL_QUADRATURE applies quadrature to a monomial in a sphere.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    05 July 2007

  Author:

    John Burkardt

  Parameters:

    Input, int DIM_NUM, the spatial dimension.

    Input, int EXPON[DIM_NUM], the exponents.

    Input, int POINT_NUM, the number of points in the rule.

    Input, double XYZ[DIM_NUM*POINT_NUM], the quadrature points.

    Input, double W[POINT_NUM], the quadrature weights.

    Output, double SPHERE01_MONOMIAL_QUADRATURE, the quadrature error.
*/
{
  double exact;
  double quad;
  double quad_error;
  double *value;
  double volume;
/*
  Get the exact value of the integral.
*/
  exact = sphere01_monomial_integral ( expon );
/*
  Evaluate the monomial at the quadrature points.
*/
  value = monomial_value ( 3, point_num, expon, xyz );
/*
  Compute the weighted sum.
*/
  quad = r8vec_dot_product ( point_num, w, value );
/*
  Error:
*/
  quad_error = fabs ( quad - exact );

  free ( value );

  return quad_error;
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
