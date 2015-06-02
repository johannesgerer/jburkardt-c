# include <math.h>
# include <stdio.h>
# include <stdlib.h>
# include <time.h>

int main ( int argc, char *argv[] );
int i4_uniform_ab ( int ilo, int ihi, int *seed );
double *latin_random_new ( int dim_num, int point_num, int *seed );
int *perm_uniform_new ( int n, int *seed );
double *r8mat_uniform_01_new ( int m, int n, int *seed );
void r8mat_write ( char *output_filename, int m, int n, double table[] );
void timestamp ( );

/******************************************************************************/

int main ( int argc, char *argv[] )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for LATIN_RANDOM_DATASET.

  Discussion:

    LATIN_RANDOM_DATASET generates a Latin Random Square dataset 
    and writes it to a file.

  Usage:

    latin_random_dataset m n seed

    where

    * M, the spatial dimension,
    * N, the number of points to generate,
    * SEED, the seed, a positive integer.

    creates an M by N dataset and writes it to the
    file "latin_random_M_N.txt".

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    13 November 2014

  Author:

    John Burkardt
*/
{
  int i;
  int m;
  int n;
  char output_filename[255];
  double *r;
  int seed;

  timestamp ( );

  printf ( "\n" );
  printf ( "LATIN_RANDOM_DATASET\n" );
  printf ( "  C version\n" );
  printf ( "  Generate a Latin Random Square dataset.\n" );
/*
  Get the spatial dimension.
*/
  if ( 1 < argc )
  {
    m = atoi ( argv[1] );
  }
  else
  {
    printf ( "\n" );
    printf ( "  Enter the value of M\n" );
    scanf ( "%d", &m );
  }

  printf ( "\n" );
  printf ( "  Spatial dimension M = %d\n", m );
/*
  Get the number of points.
*/
  if ( 2 < argc )
  {
    n = atoi ( argv[2] );
  }
  else
  {
    printf ( "\n" );
    printf ( "  Enter the number of points N\n" );
    scanf ( "%d", &n );
  }

  printf ( "  Number of points N = %d\n", n );
/*
  Get the seed.
*/
  if ( 3 < argc )
  {
    seed = atoi ( argv[3] );
  }
  else
  {
    printf ( "\n" );
    printf ( "  Enter the value of SEED\n" );
    scanf ( "%d", &seed );
  }

  printf ( "  The seed is = %d\n", seed );
/*
  Compute the data.
*/
  r = latin_random_new ( m, n, &seed );
/*
  Write it to a file.
*/
  sprintf ( output_filename, "latin_random_%d_%d.txt", m, n );

  r8mat_write ( output_filename, m, n, r );

  printf ( "\n" );
  printf ( "  The data was written to the file \"%s\"\n", output_filename );
/*
  Free memory.
*/
  free ( r );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "LATIN_RANDOM_DATASET:\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

int i4_uniform_ab ( int a, int b, int *seed )

/******************************************************************************/
/*
  Purpose:

    I4_UNIFORM_AB returns a scaled pseudorandom I4 between A and B.

  Discussion:

    The pseudorandom number should be uniformly distributed
    between A and B.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    24 May 2012

  Author:

    John Burkardt

  Reference:

    Paul Bratley, Bennett Fox, Linus Schrage,
    A Guide to Simulation,
    Second Edition,
    Springer, 1987,
    ISBN: 0387964673,
    LC: QA76.9.C65.B73.

    Bennett Fox,
    Algorithm 647:
    Implementation and Relative Efficiency of Quasirandom
    Sequence Generators,
    ACM Transactions on Mathematical Software,
    Volume 12, Number 4, December 1986, pages 362-376.

    Pierre L'Ecuyer,
    Random Number Generation,
    in Handbook of Simulation,
    edited by Jerry Banks,
    Wiley, 1998,
    ISBN: 0471134031,
    LC: T57.62.H37.

    Peter Lewis, Allen Goodman, James Miller,
    A Pseudo-Random Number Generator for the System/360,
    IBM Systems Journal,
    Volume 8, Number 2, 1969, pages 136-143.

  Parameters:

    Input, int A, B, the limits of the interval.

    Input/output, int *SEED, the "seed" value, which should NOT be 0.
    On output, SEED has been updated.

    Output, int I4_UNIFORM_AB, a number between A and B.
*/
{
  int c;
  const int i4_huge = 2147483647;
  int k;
  float r;
  int value;

  if ( *seed == 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "I4_UNIFORM_AB - Fatal error!\n" );
    fprintf ( stderr, "  Input value of SEED = 0.\n" );
    exit ( 1 );
  }
/*
  Guaranteee A <= B.
*/
  if ( b < a )
  {
    c = a;
    a = b;
    b = c;
  }

  k = *seed / 127773;

  *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

  if ( *seed < 0 )
  {
    *seed = *seed + i4_huge;
  }

  r = ( float ) ( *seed ) * 4.656612875E-10;
/*
  Scale R to lie between A-0.5 and B+0.5.
*/
  r = ( 1.0 - r ) * ( ( float ) ( a ) - 0.5 ) 
    +         r   * ( ( float ) ( b ) + 0.5 );
/*
  Round R to the nearest integer.
*/
  value = round ( r );
/*
  Guarantee that A <= VALUE <= B.
*/
  if ( value < a )
  {
    value = a;
  }
  if ( b < value )
  {
    value = b;
  }

  return value;
}
/******************************************************************************/

double *latin_random_new ( int dim_num, int point_num, int *seed )

/******************************************************************************/
/*
  Purpose:

    LATIN_RANDOM_NEW returns points in a Latin Random square.

  Discussion:

    In each spatial dimension, there will be exactly one
    point whose coordinate value lies between consecutive
    values in the list:

      ( 0, 1, 2, ..., point_num ) / point_num

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    08 April 2003

  Author:

    John Burkardt

  Parameters:

    Input, int DIM_NUM, the spatial dimension.

    Input, int POINT_NUM, the number of points.

    Input/output, int *SEED, a seed for UNIFORM.

    Output, double LATIN_RANDOM_NEW[DIM_NUM,POINT_NUM], the points.
*/
{
  int i;
  int j;
  int *perm;
  double r;
  double *x;

  x = r8mat_uniform_01_new ( dim_num, point_num, seed );
/*
  For spatial dimension I, 
    pick a random permutation of 1 to POINT_NUM,
    force the corresponding I-th components of X to lie in the
    interval ( PERM[J]-1, PERM[J] ) / POINT_NUM.
*/
  for ( i = 0; i < dim_num; i++ )
  {
    perm = perm_uniform_new ( point_num, seed );

    for ( j = 0; j < point_num; j++ )
    {
      x[i+j*dim_num] = ( ( ( double ) perm[j] ) + x[i+j*dim_num] ) 
                       / ( ( double ) point_num );
    }
    free ( perm );
  }
  return x;
}
/******************************************************************************/

int *perm_uniform_new ( int n, int *seed )

/******************************************************************************/
/*
  Purpose:

    PERM_UNIFORM_NEW selects a random permutation of N objects.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    25 February 2014

  Author:

    John Burkardt

  Reference:

    Albert Nijenhuis, Herbert Wilf,
    Combinatorial Algorithms,
    Academic Press, 1978, second edition,
    ISBN 0-12-519260-6.

  Parameters:

    Input, int N, the number of objects to be permuted.

    Input/output, int *SEED, a seed for the random number generator.

    Output, int PERM_UNIFORM_NEW[N], a permutation of
    (BASE, BASE+1, ..., BASE+N-1).
*/
{
  int i;
  int j;
  int k;
  int *p;

  p = ( int * ) malloc ( n * sizeof ( int ) );

  for ( i = 0; i < n; i++ )
  {
    p[i] = i;
  }

  for ( i = 0; i < n - 1; i++ )
  {
    j = i4_uniform_ab ( i, n - 1, seed );
    k    = p[i];
    p[i] = p[j];
    p[j] = k;
  }

  return p;
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
  const int i4_huge = 2147483647;
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
        *seed = *seed + i4_huge;
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
    fprintf ( stderr, "  Could not open the file '%s'.\n", output_filename );
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
