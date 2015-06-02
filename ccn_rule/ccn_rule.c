# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <string.h>
# include <time.h>

int main ( int argc, char *argv[] );
double *ccn_compute_points_new ( int n );
int i4_min ( int i1, int i2 );
double *nc_compute_new ( int n, double x_min, double x_max, double x[] );
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

    MAIN is the main program for CCN_RULE.

  Discussion:

    This program computes a nested Clenshaw Curtis quadrature rule
    and writes it to a file.

    The user specifies:
    * N, the number of points in the rule;
    * A, the left endpoint;
    * B, the right endpoint;
    * FILENAME, which defines the output filenames.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 April 2014

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  char filename[255];
  int n;
  double *r;
  double *w;
  double *x;
  double x_max;
  double x_min;

  timestamp ( );
  printf ( "\n" );
  printf ( "CCN_RULE\n" );
  printf ( "  C version\n" );
  printf ( "\n" );
  printf ( "  Compute one of a family of nested Clenshaw Curtis rules\n" );
  printf ( "  for approximating\n" );
  printf ( "    Integral ( -1 <= x <= +1 ) f(x) dx\n" );
  printf ( "  of order N.\n" );
  printf ( "\n" );
  printf ( "  The user specifies N, A, B and FILENAME.\n" );
  printf ( "\n" );
  printf ( "  N is the number of points.\n" );
  printf ( "  A is the left endpoint.\n" );
  printf ( "  B is the right endpoint.\n" );
  printf ( "  FILENAME is used to generate 3 files:\n" );
  printf ( "    filename_w.txt - the weight file\n" );
  printf ( "    filename_x.txt - the abscissa file.\n" );
  printf ( "    filename_r.txt - the region file.\n" );
/*
  Get N.
*/
  if ( 1 < argc )
  {
    n = atoi ( argv[1] );
  }
  else
  {
    printf ( "\n" );
    printf ( "  Enter the value of N (1 or greater)\n" );
    scanf ( "%d", &n );
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
    printf ( "  Enter FILENAME, the \"root name\" of the quadrature files.\n" );
    scanf ( "%s", filename );
  }
/*
  Input summary.
*/
  printf ( "\n" );
  printf ( "  N = %d\n", n );
  printf ( "  A = %g\n", a );
  printf ( "  B = %g\n", b );
  printf ( "  FILENAME = \"%s\".\n", filename );
/*
  Construct the rule.
*/
  r = ( double * ) malloc ( 2 * sizeof ( double ) );

  r[0] = a;
  r[1] = b;

  x = ccn_compute_points_new ( n );

  x_min = -1.0;
  x_max = +1.0;
  w = nc_compute_new ( n, x_min, x_max, x );
/*
  Rescale the rule.
*/
  rescale ( a, b, n, x, w );
/*
  Output the rule.
*/
  rule_write ( n, filename, x, w, r );
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
  printf ( "CCN_RULE:\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

double *ccn_compute_points_new ( int n )

/******************************************************************************/
/*
  Purpose:

    CCN_COMPUTE_POINTS: compute Clenshaw Curtis Nested points.

  Discussion:

    We want to compute the following sequence:

    1/2,
    0, 1
    1/4, 3/4
    1/8, 3/8, 5/8, 7/8,
    1/16, 3/16, 5/16, 7/16, 9/16, 11/16, 13/16, 15/16, and so on.

    But we would prefer that the numbers in each row be regrouped in pairs
    that are symmetric about 1/2, with the number above 1/2 coming first.
    Thus, the last row might become:
    (9/16, 7/16), (11/16, 5/16), ..., (15/16, 1/16).

    Once we have our sequence, we apply the Chebyshev transformation
    which maps [0,1] to [-1,+1].

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    06 March 2011

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of elements to compute.

    Output, double CCN_COMPUTE_POINTS_NEW[N], the elements of the sequence.
*/
{
  int d;
  int i;
  int k;
  int m;
  double r8_pi = 3.141592653589793;
  int td;
  int tu;
  double *x;

  x = ( double * ) malloc ( n * sizeof ( double ) );
/*
  Handle first three entries specially.
*/
  if ( 1 <= n )
  {
    x[0] = 0.5;
  }

  if ( 2 <= n )
  {
    x[1] = 1.0;
  }

  if ( 3 <= n )
  {
    x[2] = 0.0;
  }

  m = 3;
  d = 2;

  while ( m < n )
  {
    tu = d + 1;
    td = d - 1;

    k = i4_min ( d, n - m );

    for ( i = 1; i <= k; i++ )
    {
      if ( ( i % 2 ) == 1 )
      {
        x[m+i-1] = tu / 2.0 / ( double ) ( k );
        tu = tu + 2;
      }
      else
      {
        x[m+i-1] = td / 2.0 / ( double ) ( k );
        td = td - 2;
      }
    }
    m = m + k;
    d = d * 2;
  }
/*
  Apply the Chebyshev transformation.
*/
  for ( i = 0; i < n; i++ )
  {
    x[i] = cos ( x[i] * r8_pi );
  }
  x[0] = 0.0;

  if ( 2 <= n )
  {
    x[1] = -1.0;
  }

  if ( 3 <= n )
  {
    x[2] = +1.0;
  }

  return x;
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

double *nc_compute_new ( int n, double x_min, double x_max, double x[] )

/******************************************************************************/
/*
  Purpose:

    NC_COMPUTE_NEW computes a Newton-Cotes quadrature rule.

  Discussion:

    For the interval [X_MIN,X_MAX], the Newton-Cotes quadrature rule
    estimates

      Integral ( X_MIN <= X <= X_MAX ) F(X) dX

    using N abscissas X and weights W:

      Sum ( 1 <= I <= N ) W(I) * F ( X(I) ).

    For the CLOSED rule, the abscissas include the end points.
    For the OPEN rule, the abscissas do not include the end points.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 November 2009

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order.

    Input, double X_MIN, X_MAX, the endpoints of the interval.

    Input, double X[N], the abscissas.

    Output, double NC_COMPUTE_NEW[N], the weights.
*/
{
  double *d;
  int i;
  int j;
  int k;
  double *w;
  double yvala;
  double yvalb;

  d = ( double * ) malloc ( n * sizeof ( double ) );
  w = ( double * ) malloc ( n * sizeof ( double ) );

  for ( i = 0; i < n; i++ )
  {
/*
  Compute the Lagrange basis polynomial which is 1 at XTAB(I),
  and zero at the other nodes.
*/
    for ( j = 0; j < n; j++ )
    {
      d[j] = 0.0;
    }
    d[i] = 1.0;

    for ( j = 2; j <= n; j++ )
    {
      for ( k = j; k <= n; k++ )
      {
        d[n+j-k-1] = ( d[n+j-k-1-1] - d[n+j-k-1] ) / ( x[n+1-k-1] - x[n+j-k-1] );
      }
    }

    for ( j = 1; j <= n - 1; j++ )
    {
      for ( k = 1; k <= n - j; k++ )
      {
        d[n-k-1] = d[n-k-1] - x[n-k-j] * d[n-k];
      }
    }
/*
  Evaluate the antiderivative of the polynomial at the left and
  right endpoints.
*/
    yvala = d[n-1] / ( double ) ( n );
    for ( j = n - 2; 0 <= j; j-- )
    {
      yvala = yvala * x_min + d[j] / ( double ) ( j + 1 );
    }
    yvala = yvala * x_min;

    yvalb = d[n-1] / ( double ) ( n );
    for ( j = n - 2; 0 <= j; j-- )
    {
      yvalb = yvalb * x_max + d[j] / ( double ) ( j + 1 );
    }
    yvalb = yvalb * x_max;

    w[i] = yvalb - yvala;
  }

  free ( d );

  return w;
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

