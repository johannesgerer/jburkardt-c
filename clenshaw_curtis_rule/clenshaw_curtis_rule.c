# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <string.h>

int main ( int argc, char *argv[] );
void clenshaw_curtis_compute ( int order, double xtab[], double weight[] );
void r8mat_write ( char *output_filename, int m, int n, double table[] );
void rescale ( double a, double b, int n, double x[], double w[] );
void rule_write ( int order, char *filename, double x[], double w[],
  double r[] );
void timestamp ( );

/******************************************************************************/

int main ( int argc, char *argv[] )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for CLENSHAW_CURTIS_RULE.

  Discussion:

    This program computes a standard Clenshaw Curtis quadrature rule
    and writes it to a file.

    The user specifies:
    * the ORDER (number of points) in the rule;
    * A, the left endpoint;
    * B, the right endpoint;
    * FILENAME, which defines the output filenames.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 October 2012

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  char filename[80];
  int order;
  double *r;
  double *w;
  double *x;

  timestamp ( );
  printf ( "\n" );
  printf ( "CLENSHAW_CURTIS_RULE\n" );
  printf ( "  C version\n" );
  printf ( "\n" );
  printf ( "  Compiled on %s at %s\n", __DATE__, __TIME__ );
  printf ( "\n" );
  printf ( "  Compute a Clenshaw Curtis rule for approximating\n" );
  printf ( "\n" );
  printf ( "    Integral ( -1 <= x <= +1 ) f(x) dx\n" );
  printf ( "\n" );
  printf ( "  of order ORDER.\n" );
  printf ( "\n" );
  printf ( "  The user specifies ORDER, A, B and FILENAME.\n" );
  printf ( "\n" );
  printf ( "  ORDER is the number of points.\n" );
  printf ( "\n" );
  printf ( "  A is the left endpoint.\n" );
  printf ( "\n" );
  printf ( "  B is the right endpoint.\n" );
  printf ( "\n" );
  printf ( "  FILENAME is used to generate 3 files:\n" );
  printf ( "\n" );
  printf ( "    filename_w.txt - the weight file\n" );
  printf ( "    filename_x.txt - the abscissa file.\n" );
  printf ( "    filename_r.txt - the region file.\n" );
/*
  Get ORDER.
*/
  if ( 1 < argc )
  {
    order = atoi ( argv[1] );
  }
  else
  {
    printf ( "\n" );
    printf ( "  Enter the value of ORDER (1 or greater)\n" );
    scanf ( "%d", &order );
  }
/*
  Get A.
*/
  if ( 2 < argc )
  {
    a = atof ( argv[2] );
  }
  else
  {
    printf ( "\n" );
    printf ( "  Enter the left endpoint A:\n" );
    scanf ( "%lf", &a );
  }
/*
  Get B.
*/
  if ( 3 < argc )
  {
    b = atof ( argv[3] );
  }
  else
  {
    printf ( "\n" );
    printf ( "  Enter the right endpoint B:\n" );
    scanf ( "%lf", &b );
  }
/*
  Get FILENAME:
*/
  if ( 4 < argc )
  {
    strcpy ( filename, argv[4] );
  }
  else
  {
    printf ( "\n" );
    printf ( "  Enter FILENAME, the \"root name\" of the quadrature files).\n" );
    scanf ( "%s", filename );
  }
/*
  Input summary.
*/
  printf ( "\n" );
  printf ( "  ORDER = %d\n", order );
  printf ( "  A = %g\n", a );
  printf ( "  B = %g\n", b );
  printf ( "  FILENAME = \"%s\".\n", filename );
/*
  Construct the rule.
*/
  r = ( double * ) malloc ( 2 * sizeof ( double ) );
  w = ( double * ) malloc ( order * sizeof ( double ) );
  x = ( double * ) malloc ( order * sizeof ( double ) );

  r[0] = a;
  r[1] = b;

  clenshaw_curtis_compute ( order, x, w );
/*
  Rescale the rule.
*/
  rescale ( a, b, order, x, w );
/*
  Output the rule.
*/
  rule_write ( order, filename, x, w, r );
/*
  Free memory.
*/
  free ( r );
  free ( w );
  free ( x );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "CLENSHAW_CURTIS_RULE:\n" );
  printf ( "  Normal end of execution.\n" );

  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void clenshaw_curtis_compute ( int order, double x[], double w[] )

/******************************************************************************/
/*
  Purpose:

    CLENSHAW_CURTIS_COMPUTE computes a Clenshaw Curtis quadrature rule.

  Discussion:

    The integration interval is [ -1, 1 ].

    The weight function is w(x) = 1.0.

    The integral to approximate:

      Integral ( -1 <= X <= 1 ) F(X) dX

    The quadrature rule:

      Sum ( 1 <= I <= ORDER ) W(I) * F ( X(I) )

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 October 2012

  Author:

    John Burkardt

  Parameters:

    Input, int ORDER, the order of the rule.
    1 <= ORDER.

    Output, double X[ORDER], the abscissas.

    Output, double W[ORDER], the weights.
*/
{
  double b;
  int i;
  int j;
  double pi = 3.141592653589793;
  double theta;

  if ( order < 1 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "CLENSHAW_CURTIS_COMPUTE - Fatal error!\n" );
    fprintf ( stderr, "  Illegal value of ORDER = %d\n", order );
    exit ( 1 );
  }
  else if ( order == 1 )
  {
    x[0] = 0.0;
    w[0] = 2.0;
  }
  else
  {
    for ( i = 0; i < order; i++ )
    {
      x[i] =  cos ( ( double ) ( order - 1 - i ) * pi
                  / ( double ) ( order - 1     ) );
    }
    x[0] = -1.0;
    if ( ( order % 2 ) == 1 )
    {
      x[(order-1)/2] = 0.0;
    }
    x[order-1] = +1.0;

    for ( i = 0; i < order; i++ )
    {
      theta = ( double ) ( i ) * pi / ( double ) ( order - 1 );

      w[i] = 1.0;

      for ( j = 1; j <= ( order - 1 ) / 2; j++ )
      {
        if ( 2 * j == ( order - 1 ) )
        {
          b = 1.0;
        }
        else
        {
          b = 2.0;
        }

        w[i] = w[i] - b * cos ( 2.0 * ( double ) ( j ) * theta )
          / ( double ) ( 4 * j * j - 1 );
      }
    }

    w[0] = w[0] / ( double ) ( order - 1 );
    for ( i = 1; i < order - 1; i++ )
    {
      w[i] = 2.0 * w[i] / ( double ) ( order - 1 );
    }
    w[order-1] = w[order-1] / ( double ) ( order - 1 );
  }

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
      fprintf ( output, "  %24.16g", table[i+j*m] );
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

void rescale ( double a, double b, int n, double x[], double w[] )

/******************************************************************************/
/*
  Purpose:

    RESCALE rescales a quadrature rule from [-1,+1] to [A,B].

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 October 2012

  Author:

    John Burkardt.

  Parameters:

    Input, double A, B, the endpoints of the new interval.

    Input, int N, the order.

    Input/output, double X[N], on input, the abscissas for [-1,+1].
    On output, the abscissas for [A,B].

    Input/output, double W[N], on input, the weights for [-1,+1].
    On output, the weights for [A,B].
*/
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    x[i] = ( ( a + b ) + ( b - a ) * x[i] ) / 2.0;
  }
  for ( i = 0; i < n; i++ )
  {
    w[i] = ( b - a ) * w[i] / 2.0;
  }
  return;
}
/******************************************************************************/

void rule_write ( int order, char *filename, double x[], double w[],
  double r[] )

/******************************************************************************/
/*
  Purpose:

    RULE_WRITE writes a quadrature rule to three files.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 October 2012

  Author:

    John Burkardt

  Parameters:

    Input, int ORDER, the order of the rule.

    Input, double A, the left endpoint.

    Input, double B, the right endpoint.

    Input, char *FILENAME, specifies the output filenames.
    "filename_w.txt", "filename_x.txt", "filename_r.txt"
    defining weights, abscissas, and region.
*/
{
  char filename_r[80];
  char filename_w[80];
  char filename_x[80];

  strcpy ( filename_r, filename );
  strcat ( filename_r, "_r.txt" );
  strcpy ( filename_w, filename );
  strcat ( filename_w, "_w.txt" );
  strcpy ( filename_x, filename );
  strcat ( filename_x, "_x.txt" );

  printf ( "\n" );
  printf ( "  Creating quadrature files.\n" );
  printf ( "\n" );
  printf ( "  Root file name is     \"%s\".\n", filename );
  printf ( "\n" );
  printf ( "  Weight file will be   \"%s\".\n", filename_w );
  printf ( "  Abscissa file will be \"%s\".\n", filename_x );
  printf ( "  Region file will be   \"%s\".\n", filename_r );

  r8mat_write ( filename_w, 1, order, w );
  r8mat_write ( filename_x, 1, order, x );
  r8mat_write ( filename_r, 1, 2,     r );

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

  fprintf ( stdout, "%s\n", time_buffer );

  return;
# undef TIME_SIZE
}

