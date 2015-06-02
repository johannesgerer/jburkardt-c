# include <stdlib.h>
# include <math.h>
# include <stdio.h>
# include <time.h>
# include <string.h>

int main ( int argc, char *argv[] );
double *r8mat_uniform_01_new ( int m, int n, int *seed );
void r8mat_write ( char *output_filename, int m, int n, double table[] );
void timestamp ( );

/******************************************************************************/

int main ( int argc, char *argv[] )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for MONTE_CARLO_RULE.

  Discussion:

    MONTE_CARLO_RULE generates N points in the M-dimensional unit hypercube,
    and writes out files so that the data can be regarded as a quadrature rule.

  Usage:

    monte_carlo_rule m n seed

    where

    * M, the spatial dimension,
    * N, the number of points to generate,
    * SEED, the seed, a positive integer.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 March 2013

  Author:

    John Burkardt
*/
{
  char filename_r[255];
  char filename_w[255];
  char filename_x[255];
  int i;
  int m;
  int n;
  double *r;
  int s;
  int seed;
  double *w;
  double *x;

  timestamp ( );
  printf ( "\n" );
  printf ( "MONTE_CARLO_RULE\n" );
  printf ( "  C++ version\n" );
  printf ( "  Compiled on %s at %s.\n", __DATE__, __TIME__ );
  printf ( "\n" );
  printf ( "  Compute the abscissas and weights of a quadrature rule\n" );
  printf ( "  that is simply a Monte Carlo sampling.\n" );
  printf ( "\n" );
  printf ( "  The program requests input values from the user:\n" );
  printf ( "\n" );
  printf ( "  * M, the spatial dimension,\n" );
  printf ( "  * N, the number of points to generate,\n" );
  printf ( "  * SEED, a positive integer.\n" );
  printf ( "\n" );
  printf ( "  Output from the program includes\n" );
  printf ( "  a set of 3 files that define the quadrature rule.\n" );
  printf ( "\n" );
  printf ( "    (1) \"mc_m?_n?_s?_r.txt\", the ranges;\n" );
  printf ( "    (2) \"mc_m?_n?_s?_w.txt\", the weights;\n" );
  printf ( "    (3) \"mc_m?_n?_s?_x.txt\", the abscissas.\n" );
/*
  Get the spatial dimension M.
*/
  if ( 1 < argc )
  {
    m = atoi ( argv[1] );
  }
  else
  {
    printf ( "\n" );
    printf ( "  Enter the spatial dimension M (1 or greater)\n" );
    scanf ( "%d", &m );
  }
/*
  Get the number of points N.
*/
  if ( 2 < argc )
  {
    n = atoi ( argv[2] );
  }
  else
  {
    printf ( "\n" );
    printf ( "  Enter the number of points N (1 or greater):\n" );
    scanf ( "%d", &n );
  }
/*
  Get the seed S.
*/
  if ( 3 < argc )
  {
    s = atoi ( argv[3] );
  }
  else
  {
    printf ( "\n" );
    printf ( "  Enter the seed S (1 or greater):\n" );
    scanf ( "%d", &s );
  }
/*
  Input summary.
*/
  printf ( "\n" );
  printf ( "  M = %d\n", m );
  printf ( "  N = %d\n", n );
  printf ( "  S = %d\n", s );
/*
  Construct the rule.
*/
  r = ( double * ) malloc ( m * 2 * sizeof ( double ) );

  for ( i = 0; i < m; i++ )
  {
    r[i+0*m] = 0.0;
    r[i+1*m] = 1.0;
  }

  w = ( double * ) malloc ( n * sizeof ( double ) );

  for ( i = 0; i < n; i++ )
  {
    w[i] = 1.0 / ( double ) n;
  }

  seed = s;
  x = r8mat_uniform_01_new ( m, n, &seed );
/*
  Output the rule.
*/
  sprintf ( filename_r, "mc_d%d_n%d_s%d_r.txt", m, n, s );
  sprintf ( filename_w, "mc_d%d_n%d_s%d_w.txt", m, n, s );
  sprintf ( filename_x, "mc_d%d_n%d_s%d_x.txt", m, n, s );

  printf ( "\n" );
  printf ( "  Region file will be   \"%s\".\n", filename_r );
  printf ( "  Weight file will be   \"%s\".\n", filename_w );
  printf ( "  Abscissa file will be \"%s\".\n", filename_x );

  r8mat_write ( filename_r, m, 2, r );
  r8mat_write ( filename_w, 1, n, w );
  r8mat_write ( filename_x, m, n, x );
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
  printf ( "MONTE_CARLO_RULE:\n" );
  printf ( "  Normal end of execution.\n" );

  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

double *r8mat_uniform_01_new ( int m, int n, int *seed )

/******************************************************************************/
/*
  Purpose:

    R8MAT_UNIFORM_01_NEW fills an R8MAT with pseudorandom values scaled to [0,1].

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

    This routine implements the recursion

      seed = 16807 * seed mod ( 2^31 - 1 )
      unif = seed / ( 2^31 - 1 )

    The integer arithmetic never requires more than 32 bits,
    including a sign bit.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 June 2009

  Author:

    John Burkardt

  Reference:

    Paul Bratley, Bennett Fox, Linus Schrage,
    A Guide to Simulation,
    Springer Verlag, pages 201-202, 1983.

    Bennett Fox,
    Algorithm 647:
    Implementation and Relative Efficiency of Quasirandom
    Sequence Generators,
    ACM Transactions on Mathematical Software,
    Volume 12, Number 4, pages 362-376, 1986.

    Philip Lewis, Allen Goodman, James Miller,
    A Pseudo-Random Number Generator for the System/360,
    IBM Systems Journal,
    Volume 8, pages 136-143, 1969.

  Parameters:

    Input, int M, N, the number of rows and columns.

    Input/output, int *SEED, the "seed" value.  Normally, this
    value should not be 0, otherwise the output value of SEED
    will still be 0, and R8_UNIFORM will be 0.  On output, SEED has
    been updated.

    Output, double R8MAT_UNIFORM_01_NEW[M*N], a matrix of pseudorandom values.
*/
{
  int i;
  int j;
  int k;
  double *r;

  r = ( double * ) malloc ( m * n * sizeof ( double ) );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      k = *seed / 127773;

      *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

      if ( *seed < 0 )
      {
        *seed = *seed + 2147483647;
      }
      r[i+j*m] = ( double ) ( *seed ) * 4.656612875E-10;
    }
  }

  return r;
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
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  printf ( "%s\n", time_buffer );

  return;
# undef TIME_SIZE
}
