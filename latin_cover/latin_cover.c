# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <string.h>

# include "latin_cover.h"

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

int i4_modp ( int i, int j )

/******************************************************************************/
/*
  Purpose:

    I4_MODP returns the nonnegative remainder of I4 division.

  Discussion:

    If
      NREM = I4_MODP ( I, J )
      NMULT = ( I - NREM ) / J
    then
      I = J * NMULT + NREM
    where NREM is always nonnegative.

    The MOD function computes a result with the same sign as the
    quantity being divided.  Thus, suppose you had an angle A,
    and you wanted to ensure that it was between 0 and 360.
    Then mod(A,360) would do, if A was positive, but if A
    was negative, your result would be between -360 and 0.

    On the other hand, I4_MODP(A,360) is between 0 and 360, always.

  Example:

        I         J     MOD  I4_MODP   I4_MODP Factorization

      107        50       7       7    107 =  2 *  50 + 7
      107       -50       7       7    107 = -2 * -50 + 7
     -107        50      -7      43   -107 = -3 *  50 + 43
     -107       -50      -7      43   -107 =  3 * -50 + 43

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 January 2007

  Author:

    John Burkardt

  Parameters:

    Input, int I, the number to be divided.

    Input, int J, the number that divides I.

    Output, int I4_MODP, the nonnegative remainder when I is
    divided by J.
*/
{
  int value;

  if ( j == 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "I4_MODP - Fatal error!\n" );
    fprintf ( stderr, "  I4_MODP ( I, J ) called with J = %d\n", j );
    exit ( 1 );
  }

  value = i % j;

  if ( value < 0 )
  {
    value = value + abs ( j );
  }

  return value;
}
/******************************************************************************/

int i4_uniform ( int a, int b, int *seed )

/******************************************************************************/
/*
  Purpose:

    I4_UNIFORM returns a scaled pseudorandom I4.

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

    Output, int I4_UNIFORM, a number between A and B.
*/
{
  int c;
  int i4_huge = 2147483647;
  int k;
  float r;
  int value;

  if ( *seed == 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "I4_UNIFORM - Fatal error!\n" );
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
  r = ( 1.0 - r ) * ( ( float ) a - 0.5 ) 
    +         r   * ( ( float ) b + 0.5 );
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

int i4_wrap ( int ival, int ilo, int ihi )

/******************************************************************************/
/*
  Purpose:

    I4_WRAP forces an I4 to lie between given limits by wrapping.

  Example:

    ILO = 4, IHI = 8

    I   Value

    -2     8
    -1     4
     0     5
     1     6
     2     7
     3     8
     4     4
     5     5
     6     6
     7     7
     8     8
     9     4
    10     5
    11     6
    12     7
    13     8
    14     4

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 July 2008

  Author:

    John Burkardt

  Parameters:

    Input, int IVAL, an integer value.

    Input, int ILO, IHI, the desired bounds for the integer value.

    Output, int I4_WRAP, a "wrapped" version of IVAL.
*/
{
  int jhi;
  int jlo;
  int value;
  int wide;

  jlo = i4_min ( ilo, ihi );
  jhi = i4_max ( ilo, ihi );

  wide = jhi + 1 - jlo;

  if ( wide == 1 )
  {
    value = jlo;
  }
  else
  {
    value = jlo + i4_modp ( ival - jlo, wide );
  }

  return value;
}
/******************************************************************************/

void i4block_print ( int l, int m, int n, int a[], char *title )

/******************************************************************************/
/*
  Purpose:

    I4BLOCK_PRINT prints an I4BLOCK.

  Discussion:

    An I4BLOCK is a 3D array of I4 values.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    25 June 2012

  Author:

    John Burkardt

  Parameters:

    Input, int L, M, N, the dimensions of the block.

    Input, int A[L*M*N], the matrix to be printed.

    Input, char *TITLE, a title.
*/
{
  int i;
  int j;
  int jhi;
  int jlo;
  int k;

  printf ( "\n" );
  printf ( "%s\n", title );

  for ( k = 0; k < n; k++ )
  {
    printf ( "\n" );
    printf ( "  K = %d\n", k );
    for ( jlo = 0; jlo < m; jlo = jlo + 10 )
    {
      jhi = i4_min ( jlo + 10, m );
      printf ( "\n" );
      printf ( "        J:" );
      for ( j = jlo; j < jhi; j++ )
      {
        printf ( "  %6d", j );
      }
      printf ( "\n" );
      printf ( "       I:\n" );
      for ( i = 0; i < l; i++ )
      {
        printf ( "  %6d: ", i );
        for ( j = jlo; j < jhi; j++ )
        {
          printf ( "  %6d", a[i+j*l+k*l*m] );
        }
        printf ( "\n" );
      }
    }
  }

  return;
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
    j2hi = i4_min ( j2hi, n );
    j2hi = i4_min ( j2hi, jhi );

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
    i2lo = i4_max ( ilo, 1 );
    i2hi = i4_min ( ihi, m );

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

int *latin_cover ( int n, int p[] )

/******************************************************************************/
/*
  Purpose:

    LATIN_COVER returns a 2D Latin Square Covering.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    05 August 2012

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of points.

    Input, int P[N], a permutation which describes the
    first Latin square.

    Output, int LATIN_COVER[N*N], the Latin cover.  A(I,J) = K
    means that (I,J) is one element of the K-th Latin square.
*/
{
  int *a;
  int i;
  int ik;
  int j;
  int k;

  a = ( int * ) malloc ( n * n * sizeof ( int ) );

  perm_check ( n, p );

  for ( i = 0; i < n; i++ )
  {
    for ( k = 0; k < n; k++ )
    {
      ik = i4_wrap ( i + k, 0, n - 1 );
      j = p[ik] - 1;
      a[i+j*n] = k + 1;
    }
  }
  return a;
}
/******************************************************************************/

int *latin_cover_2d ( int n, int p1[], int p2[] )

/******************************************************************************/
/*
  Purpose:

    LATIN_COVER_2D returns a 2D Latin Square Covering.

  Discussion:

    This procedure has a chance of being extended to M dimensions.

    A basic solution is computed, and the user is permitted to permute
    both the I and J coordinates.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    05 August 2012

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of points.

    Input, int P1[N], P2[N], permutations to be applied
    to the spatial dimensions.

    Output, int LATIN_COVER_2D[N*N], the Latin cover.  A(I,J) = K
    means that (I,J) is one element of the K-th Latin square.
*/
{
  int *a;
  int *b;
  int i;
  int i1;
  int j;
  int j1;

  perm_check ( n, p1 );
  perm_check ( n, p2 );

  a = ( int * ) malloc ( n * n * sizeof ( int ) );
  b = ( int * ) malloc ( n * n * sizeof ( int ) );
/*
  Set up the basic solution.
*/
  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      a[i+j*n] = i4_wrap ( i - j + 1, 1, n );
    }
  }
/*
  Apply permutation to dimension I.
*/
  for ( i = 0; i < n; i++ )
  {
    i1 = p1[i] - 1;
    for ( j = 0; j < n; j++ )
    {
      b[i1+j*n] = a[i+j*n];
    }
  }
/*
  Apply permutation to dimension J.
*/
  for ( j = 0; j < n; j++ )
  {
    j1 = p2[j] - 1;
    for ( i = 0; i < n; i++ )
    {
      a[i+j1*n] = b[i+j*n];
    }
  }

  free ( b );

  return a;
}
/******************************************************************************/

int *latin_cover_3d ( int n, int p1[], int p2[], int p3[] )

/******************************************************************************/
/*
  Purpose:

    LATIN_COVER_3D returns a 3D Latin Square Covering.

  Discussion:

    A basic solution is computed, and the user is permitted to permute
    I, J and K coordinates.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    24 June 2012

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of points.

    Input, int P1[N], P2[N], P3[N], permutations to be applied
    to the spatial dimensions.

    Output, int LATIN_COVER_3D[N*N*N], the Latin cover.  A(I,J,K) = L
    means that (I,J,K) is one element of the L-th Latin square.
*/
{
  int *a;
  int *b;
  int i;
  int i1;
  int ik;
  int j;
  int j1;
  int jk;
  int k;
  int k1;

  perm_check ( n, p1 );
  perm_check ( n, p2 );
  perm_check ( n, p3 );

  a = ( int * ) malloc ( n * n * n * sizeof ( int ) );
  b = ( int * ) malloc ( n * n * n * sizeof ( int ) );
/*
  Set up the basic solution.
*/
  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      for ( k = 0; k < n; k++ )
      {
        ik = i4_wrap ( i + 1 - k, 1, n );
        jk = i4_wrap ( j + 1 - k, 1, n );
        b[i+j*n+k*n*n] = ik + ( jk - 1 ) * n;
      }
    }
  }
/*
  Apply permutation to dimension I.
*/
  for ( i = 0; i < n; i++ )
  {
    i1 = p1[i] - 1;
    for ( j = 0; j < n; j++ )
    {
      for ( k = 0; k < n; k++ )
      {
        a[i1+j*n+k*n*n] = b[i+j*n+k*n*n];
      }
    }
  }
/*
  Apply permutation to dimension J.
*/
  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      j1 = p2[j] - 1;
      for ( k = 0; k < n; k++ )
      {
        b[i+j1*n+k*n*n] = a[i+j*n+k*n*n];
      }
    }
  }
/*
  Apply permutation to dimension K.
*/
  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      for ( k = 0; k < n; k++ )
      {
        k1 = p3[k] - 1;
        a[i+j*n+k1*n*n] = b[i+j*n+k*n*n];
      }
    }
  }

  free ( b );

  return a;
}
/******************************************************************************/

int perm_check ( int n, int p[] )

/******************************************************************************/
/*
  Purpose:

    PERM_CHECK checks that a vector represents a permutation.

  Discussion:

    The routine verifies that each of the integers from 1
    to N occurs among the N entries of the permutation.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    28 March 2009

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries.

    Input, int P[N], the permutation, in standard index form.

    Output, int PERM_CHECK, is 1 if the array is NOT a permutation.
*/
{
  int error;
  int find;
  int seek;

  for ( seek = 1; seek <= n; seek++ )
  {
    error = 1;

    for ( find = 1; find <= n; find++ )
    {
      if ( p[find-1] == seek )
      {
        error = 0;
        break;
      }
    }

    if ( error )
    {
      fprintf ( stderr, "\n" );
      fprintf ( stderr, "PERM_CHECK - Fatal error!\n" );
      fprintf ( stderr, "  The permutation was missing the entry %d\n", seek );
      i4vec_print ( n, p, "  The permutation:" );
      exit ( 1 );
    }
  }

  return 0;
}
/******************************************************************************/

void perm_print ( int n, int p[], char *title )

/******************************************************************************/
/*
  Purpose:

    PERM_PRINT prints a permutation.

  Discussion:

    The permutation is assumed to be 0-based.

  Example:

    Input:

      P = 6 1 3 0 4 2 5

    Printed output:

      "This is the permutation:"

      0 1 2 3 4 5 6
      6 1 3 0 4 2 5

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    14 May 2011

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of objects permuted.

    Input, int P[N], the permutation, in standard index form.

    Input, char *TITLE, an optional title.
    If no title is supplied, then only the permutation is printed.
*/
{
  int i;
  int ihi;
  int ilo;
  int inc = 20;

  if ( s_len_trim ( title ) != 0 )
  {
    printf ( "\n" );
    printf ( "%s\n", title );

    for ( ilo = 0; ilo < n; ilo = ilo + inc )
    {
      ihi = ilo + inc;
      if ( n < ihi ) 
      {
        ihi = n;
      }
      printf ( "\n" );
      printf ( "  " );
      for ( i = ilo; i < ihi; i++ )
      {
        printf ( "%4d", i );
      }
      printf ( "\n" );
      printf ( "  " );
      for ( i = ilo; i < ihi; i++ )
      {
        printf ( "%4d", p[i] );
      }
      printf ( "\n" );
    }
  }
  else
  {
    for ( ilo = 0; ilo < n; ilo = ilo + inc )
    {
      ihi = ilo + inc;
      if ( n < ihi ) 
      {
        ihi = n;
      }
      printf ( "  " );
      for ( i = ilo; i < ihi; i++ )
      {
        printf ( "%4d", p[i] );
      }
      printf ( "\n" );
    }
  }

  return;
}
/******************************************************************************/

int *perm_uniform_new ( int n, int *seed )

/******************************************************************************/
/*
  Purpose:

    PERM_UNIFORM selects a random permutation of N objects.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    05 August 2012

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

    Input, int N, the number of objects to be permuted.

    Input/output, int *SEED, a seed for the random number generator.

    Output, int PERM_UNIFORM_NEW[N], a permutation of (1, 2, ..., N).
*/
{
  int i;
  int j;
  int k;
  int *p;

  p = ( int * ) malloc ( n * sizeof ( int ) );

  for ( i = 0; i < n; i++ )
  {
    p[i] = i + 1;
  }

  for ( i = 0; i < n; i++ )
  {
    j = i4_uniform ( i, n - 1, seed );
    k    = p[i];
    p[i] = p[j];
    p[j] = k;
  }
 
  return p;
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
