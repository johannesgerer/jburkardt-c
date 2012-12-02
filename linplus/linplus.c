# include <stdlib.h>
# include <stdio.h>
# include <string.h>
# include <time.h>
# include <math.h>

# include "linplus.h"

/******************************************************************************/

int get_seed ( void )

/******************************************************************************/
/*
  Purpose:

    GET_SEED returns a random seed for the random number generator.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 November 2004

  Author:

    John Burkardt

  Parameters:

    Output, int GET_SEED, a random seed value.
*/
{
  time_t clock;
  int i;
  int i4_huge = 2147483647;
  int ihour;
  int imin;
  int isec;
  int seed;
  struct tm *lt;
  time_t tloc;
/*
  If the internal seed is 0, generate a value based on the time.
*/
  clock = time ( &tloc );
  lt = localtime ( &clock );
/*
  Hours is 1, 2, ..., 12.
*/
  ihour = lt->tm_hour;

  if ( 12 < ihour )
  {
    ihour = ihour - 12;
  }
/*
  Move Hours to 0, 1, ..., 11
*/
  ihour = ihour - 1;

  imin = lt->tm_min;

  isec = lt->tm_sec;

  seed = isec + 60 * ( imin + 60 * ihour );
/*
  We want values in [1,43200], not [0,43199].
*/
  seed = seed + 1;
/*
  Remap SEED from [1,43200] to [1,HUGE].
*/
  seed = ( int ) 
    ( ( ( double ) seed )
    * ( ( double ) i4_huge ) / ( 60.0 * 60.0 * 12.0 ) );
/*
  Never use a seed of 0.
*/
  if ( seed == 0 )
  {
    seed = 1;
  }

  return seed;
}
/******************************************************************************/

double *hilbert_inverse ( int n )

/******************************************************************************/
/*
  Purpose:

    HILBERT_INVERSE returns the inverse of the Hilbert matrix.

  Formula:

    A(I,J) =  (-1)**(I+J) * (N+I-1)! * (N+J-1)! /
           [ (I+J-1) * ((I-1)!*(J-1)!)**2 * (N-I)! * (N-J)! ]

  Example:

    N = 5

       25    -300     1050    -1400     630
     -300    4800   -18900    26880  -12600
     1050  -18900    79380  -117600   56700
    -1400   26880  -117600   179200  -88200
      630  -12600    56700   -88200   44100

  Properties:

    A is symmetric.

    Because A is symmetric, it is normal, so diagonalizable.

    A is almost impossible to compute accurately by general routines
    that compute the inverse.

    A is integral.

    The sum of the entries of A is N**2.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    08 May 2010

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of A.

    Output, double HILBERT_INVERSE[N*N], the inverse Hilbert matrix.
*/
{
  double *a;
  int i;
  int j;

  a = ( double * ) malloc ( n * n * sizeof ( double ) );
/*
  Set the (1,1) entry.
*/
  a[0+0*n] = ( double ) ( n * n );
/*
  Define Row 1, Column J by recursion on Row 1 Column J-1
*/
  i = 1;

  for ( j = 2; j <= n; j++ )
  {
    a[i-1+(j-1)*n] = -a[i-1+(j-2)*n] 
      * ( double ) ( ( n + j - 1 ) * ( i + j - 2 ) * ( n + 1 - j ) ) 
      / ( double ) ( ( i + j - 1 ) * ( j - 1 ) * ( j - 1 ) );
  }
/*
  Define Row I by recursion on row I-1
*/
  for ( i = 2; i <= n; i++ )
  {
    for (  j = 1; j <= n; j++ )
    {
      a[i-1+(j-1)*n] = -a[i-2+(j-1)*n] 
        * ( double ) ( ( n + i - 1 ) * ( i + j - 2 ) * ( n + 1 - i ) ) 
        / ( double ) ( ( i + j - 1 ) * ( i - 1 ) * ( i - 1 ) );
    }
  }

  return a;
}
/******************************************************************************/

int i4_huge ( void )

/******************************************************************************/
/*
  Purpose:

    I4_HUGE returns a "huge" I4.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    29 August 2006

  Author:

    John Burkardt

  Parameters:

    Output, int I4_HUGE, a "huge" integer.
*/
{
  static int value = 2147483647;

  return value;
}
/******************************************************************************/

int i4_log_10 ( int i )

/******************************************************************************/
/*
  Purpose:

    I4_LOG_10 returns the integer part of the logarithm base 10 of an I4.

  Example:

        I  I4_LOG_10
    -----  --------
        0    0
        1    0
        2    0
        9    0
       10    1
       11    1
       99    1
      100    2
      101    2
      999    2
     1000    3
     1001    3
     9999    3
    10000    4

  Discussion:

    I4_LOG_10 ( I ) + 1 is the number of decimal digits in I.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    23 October 2007

  Author:

    John Burkardt

  Parameters:

    Input, int I, the number whose logarithm base 10 is desired.

    Output, int I4_LOG_10, the integer part of the logarithm base 10 of
    the absolute value of X.
*/
{
  int i_abs;
  int ten_pow;
  int value;

  if ( i == 0 )
  {
    value = 0;
  }
  else
  {
    value = 0;
    ten_pow = 10;

    i_abs = abs ( i );

    while ( ten_pow <= i_abs )
    {
      value = value + 1;
      ten_pow = ten_pow * 10;
    }
  }
  return value;
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

int i4_power ( int i, int j )

/******************************************************************************/
/*
  Purpose:

    I4_POWER returns the value of I^J.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    23 October 2007

  Author:

    John Burkardt

  Parameters:

    Input, int I, J, the base and the power.  J should be nonnegative.

    Output, int I4_POWER, the value of I^J.
*/
{
  int k;
  int value;

  if ( j < 0 )
  {
    if ( i == 1 )
    {
      value = 1;
    }
    else if ( i == 0 )
    {
      printf ( "\n" );
      printf ( "I4_POWER - Fatal error!\n" );
      printf ( "  I^J requested, with I = 0 and J negative.\n" );
      exit ( 1 );
    }
    else
    {
      value = 0;
    }
  }
  else if ( j == 0 )
  {
    if ( i == 0 )
    {
      printf ( "\n" );
      printf ( "I4_POWER - Fatal error!\n" );
      printf ( "  I^J requested, with I = 0 and J = 0.\n" );
      exit ( 1 );
    }
    else
    {
      value = 1;
    }
  }
  else if ( j == 1 )
  {
    value = i;
  }
  else
  {
    value = 1;
    for ( k = 1; k <= j; k++ )
    {
      value = value * i;
    }
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

    12 November 2006

  Author:

    John Burkardt

  Reference:

    Paul Bratley, Bennett Fox, Linus Schrage,
    A Guide to Simulation,
    Springer Verlag, pages 201-202, 1983.

    Pierre L'Ecuyer,
    Random Number Generation,
    in Handbook of Simulation,
    edited by Jerry Banks,
    Wiley Interscience, page 95, 1998.

    Bennett Fox,
    Algorithm 647:
    Implementation and Relative Efficiency of Quasirandom
    Sequence Generators,
    ACM Transactions on Mathematical Software,
    Volume 12, Number 4, pages 362-376, 1986.

    Peter Lewis, Allen Goodman, James Miller
    A Pseudo-Random Number Generator for the System/360,
    IBM Systems Journal,
    Volume 8, pages 136-143, 1969.

  Parameters:

    Input, int A, B, the limits of the interval.

    Input/output, int *SEED, the "seed" value, which should NOT be 0.
    On output, SEED has been updated.

    Output, int I4_UNIFORM, a number between A and B.
*/
{
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

  k = *seed / 127773;

  *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

  if ( *seed < 0 )
  {
    *seed = *seed + 2147483647;
  }

  r = ( float ) ( *seed ) * 4.656612875E-10;
/*
  Scale R to lie between A-0.5 and B+0.5.
*/
  r = ( 1.0 - r ) * ( ( float ) ( i4_min ( a, b ) ) - 0.5 ) 
    +         r   * ( ( float ) ( i4_max ( a, b ) ) + 0.5 );
/*
  Use rounding to convert R to an integer between A and B.
*/
  value = r4_nint ( r );

  value = i4_max ( value, i4_min ( a, b ) );
  value = i4_min ( value, i4_max ( a, b ) );

  return value;
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

int i4vec_search_binary_a ( int n, int a[], int b )

/******************************************************************************/
/*
  Purpose:

    I4VEC_SEARCH_BINARY_A searches an ascending sorted I4VEC for a value.

  Discussion:

    An I4VEC is a vector of I4's.

    Binary search is used.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 September 2008

  Author:

    John Burkardt

  Reference:

    Donald Kreher, Douglas Simpson,
    Algorithm 1.9,
    Combinatorial Algorithms,
    CRC Press, 1998, page 26.

  Parameters:

    Input, int N, the number of elements in the vector.

    Input, int A[N], the array to be searched.  A must
    be sorted in ascending order.

    Input, int B, the value to be searched for.

    Output, int I4VEC_SEARCH_BINARY_A, the result of the search.
    -1, B does not occur in A.
    I, A[I] = B.
*/
{
  int high;
  int index;
  int low;
  int mid;
/*
  Check.
*/
  if ( n <= 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "I4VEC_SEARCH_BINARY_A - Fatal error!\n" );
    fprintf ( stderr, "  The array dimension N is less than 1.\n" );
    exit ( 1 );
  }

  index = -1;

  low = 1;
  high = n;

  while ( low <= high )
  {
    mid = ( low + high ) / 2;

    if ( a[mid-1] == b )
    {
      index = mid;
      break;
    }
    else if ( a[mid-1] < b )
    {
      low = mid + 1;
    }
    else if ( b < a[mid-1] )
    {
      high = mid - 1;
    }
  }
  return index;
}
/******************************************************************************/

float r4_abs ( float x )

/******************************************************************************/
/*
  Purpose:

    R4_ABS returns the absolute value of an R4.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 January 2007

  Author:

    John Burkardt

  Parameters:

    Input, float X, the quantity whose absolute value is desired.

    Output, float R4_ABS, the absolute value of X.
*/
{
  float value;

  if ( 0.0 <= x )
  {
    value = + x;
  }
  else
  {
    value = - x;
  }
  return value;
}
/******************************************************************************/

int r4_nint ( float x )

/******************************************************************************/
/*
  Purpose:

    R4_NINT returns the nearest integer to an R4.

  Example:

        X         R4_NINT

      1.3         1
      1.4         1
      1.5         1 or 2
      1.6         2
      0.0         0
     -0.7        -1
     -1.1        -1
     -1.6        -2

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 May 2006

  Author:

    John Burkardt

  Parameters:

    Input, float X, the value.

    Output, int R4_NINT, the nearest integer to X.
*/
{
  int s;
  int value;

  if ( x < 0.0 )
  {
    s = - 1;
  }
  else
  {
    s = + 1;
  }
  value = s * ( int ) ( fabs ( x ) + 0.5 );

  return value;
}
/******************************************************************************/

float r4_uniform ( float b, float c, int *seed )

/******************************************************************************/
/*
  Purpose:

    R4_UNIFORM returns a scaled pseudorandom R4.

  Discussion:

    The pseudorandom number should be uniformly distributed
    between A and B.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    21 November 2004

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

    Input, float B, C, the limits of the interval.

    Input/output, int *SEED, the "seed" value, which should NOT be 0.
    On output, SEED has been updated.

    Output, float R4_UNIFORM, a number strictly between A and B.
*/
{
  float value;

  if ( *seed == 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R4_UNIFORM - Fatal error!\n" );
    fprintf ( stderr, "  Input value of SEED = 0.\n" );
    exit ( 1 );
  }

  value = b + ( c - b ) * r4_uniform_01 ( seed );

  return value;
}
/******************************************************************************/

float r4_uniform_01 ( int *seed )

/******************************************************************************/
/*
  Purpose:

    R4_UNIFORM_01 returns a real pseudorandom R4.

  Discussion:

    This routine implements the recursion

      seed = 16807 * seed mod ( 2**31 - 1 )
      r4_uniform_01 = seed / ( 2**31 - 1 )

    The integer arithmetic never requires more than 32 bits,
    including a sign bit.

    If the initial seed is 12345, then the first three computations are

      Input     Output      R4_UNIFORM_01
      SEED      SEED

         12345   207482415  0.096616
     207482415  1790989824  0.833995
    1790989824  2035175616  0.947702

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 November 2004

  Author:

    John Burkardt

  Reference:

    Paul Bratley, Bennett Fox, Linus Schrage,
    A Guide to Simulation,
    Springer Verlag, pages 201-202, 1983.

    Pierre L'Ecuyer,
    Random Number Generation,
    in Handbook of Simulation
    edited by Jerry Banks,
    Wiley Interscience, page 95, 1998.

    Bennett Fox,
    Algorithm 647:
    Implementation and Relative Efficiency of Quasirandom
    Sequence Generators,
    ACM Transactions on Mathematical Software,
    Volume 12, Number 4, pages 362-376, 1986.

    Peter Lewis, Allen Goodman, James Miller,
    A Pseudo-Random Number Generator for the System/360,
    IBM Systems Journal,
    Volume 8, pages 136-143, 1969.

  Parameters:

    Input/output, int *SEED, the "seed" value.  Normally, this
    value should not be 0.  On output, SEED has been updated.

    Output, float R4_UNIFORM_01, a new pseudorandom variate, strictly between
    0 and 1.
*/
{
  int k;
  float r;

  k = *seed / 127773;

  *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

  if ( *seed < 0 )
  {
    *seed = *seed + 2147483647;
  }
/*
  Although SEED can be represented exactly as a 32 bit integer,
  it generally cannot be represented exactly as a 32 bit real number!
*/
  r = ( float ) ( *seed ) * 4.656612875E-10;

  return r;
}
/******************************************************************************/

double r8_abs ( double x )

/******************************************************************************/
/*
  Purpose:

    R8_ABS returns the absolute value of an R8.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    07 May 2006

  Author:

    John Burkardt

  Parameters:

    Input, double X, the quantity whose absolute value is desired.

    Output, double R8_ABS, the absolute value of X.
*/
{
  double value;

  if ( 0.0 <= x )
  {
    value = x;
  } 
  else
  {
    value = - x;
  }
  return value;
}
/******************************************************************************/

int r8_is_int ( double r )

/******************************************************************************/
/*
  Purpose:

    R8_IS_INT is 1 if an R8 represents an integer value.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 May 2006

  Author:

    John Burkardt

  Parameters:

    Input, double R, the number to be checked.

    Output, int R8_IS_INT, is 1 if R is an integer value.
*/
{
  if ( ( double ) ( i4_huge ( ) ) < r )
  {
    return 0;
  }
  else if ( r < - ( double ) ( i4_huge ( ) ) ) 
  {
    return 0;
  }
  else if ( r == ( double ) ( ( int ) ( r ) ) )
  {
    return 1;
  }
  else
  {
    return 0;
  }

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

double r8_sign ( double x )

/******************************************************************************/
/*
  Purpose:

    R8_SIGN returns the sign of an R8.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    08 May 2006

  Author:

    John Burkardt

  Parameters:

    Input, double X, the number whose sign is desired.

    Output, double R8_SIGN, the sign of X.
*/
{
  double value;

  if ( x < 0.0 )
  {
    value = -1.0;
  } 
  else
  {
    value = 1.0;
  }
  return value;
}
/******************************************************************************/

double r8_sign2 ( double x, double y )

/******************************************************************************/
/*
  Purpose:

    R8_SIGN2 returns the first argument with the sign of the second.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    10 February 2012

  Author:

    John Burkardt

  Parameters:

    Input, double X, Y, the input arguments.

    Output, double R8_SIGN2, is equal in absolute value to the absolute value 
    of X, and has the sign of Y.
*/
{
  double value;

  if ( 0.0 <= y )
  {
    value = r8_abs ( x );
  } 
  else
  {
    value = - r8_abs ( x );
  }
  return value;
}
/******************************************************************************/

void r8_swap ( double *x, double *y )

/******************************************************************************/
/*
  Purpose:

    R8_SWAP switches two R8's.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    25 March 2009

  Author:

    John Burkardt

  Parameters:

    Input/output, double *X, *Y.  On output, the values of X and
    Y have been interchanged.
*/
{
  double z;

  z = *x;
  *x = *y;
  *y = z;
 
  return;
}
/******************************************************************************/

double r8_uniform ( double b, double c, int *seed )

/******************************************************************************/
/*
  Purpose:

    R8_UNIFORM returns a scaled pseudorandom R8.

  Discussion:

    The pseudorandom number should be uniformly distributed
    between A and B.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    21 November 2004

  Author:

    John Burkardt

  Parameters:

    Input, double B, C, the limits of the interval.

    Input/output, int *SEED, the "seed" value, which should NOT be 0.
    On output, SEED has been updated.

    Output, double R8_UNIFORM, a number strictly between A and B.
*/
{
  double value;

  value = b + ( c - b ) * r8_uniform_01 ( seed );

  return value;
}
/******************************************************************************/

double r8_uniform_01 ( int *seed )

/******************************************************************************/
/*
  Purpose:

    R8_UNIFORM_01 returns a unit pseudorandom R8.

  Discussion:

    This routine implements the recursion

      seed = 16807 * seed mod ( 2^31 - 1 )
      r8_uniform_01 = seed / ( 2^31 - 1 )

    The integer arithmetic never requires more than 32 bits,
    including a sign bit.

    If the initial seed is 12345, then the first three computations are

      Input     Output      R8_UNIFORM_01
      SEED      SEED

         12345   207482415  0.096616
     207482415  1790989824  0.833995
    1790989824  2035175616  0.947702

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 August 2004

  Author:

    John Burkardt

  Reference:

    Paul Bratley, Bennett Fox, Linus Schrage,
    A Guide to Simulation,
    Springer Verlag, pages 201-202, 1983.

    Pierre L'Ecuyer,
    Random Number Generation,
    in Handbook of Simulation
    edited by Jerry Banks,
    Wiley Interscience, page 95, 1998.

    Bennett Fox,
    Algorithm 647:
    Implementation and Relative Efficiency of Quasirandom
    Sequence Generators,
    ACM Transactions on Mathematical Software,
    Volume 12, Number 4, pages 362-376, 1986.

    P A Lewis, A S Goodman, J M Miller,
    A Pseudo-Random Number Generator for the System/360,
    IBM Systems Journal,
    Volume 8, pages 136-143, 1969.

  Parameters:

    Input/output, int *SEED, the "seed" value.  Normally, this
    value should not be 0.  On output, SEED has been updated.

    Output, double R8_UNIFORM_01, a new pseudorandom variate, strictly between
    0 and 1.
*/
{
  int k;
  double r;

  k = *seed / 127773;

  *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

  if ( *seed < 0 )
  {
    *seed = *seed + 2147483647;
  }
/*
  Although SEED can be represented exactly as a 32 bit integer,
  it generally cannot be represented exactly as a 32 bit real number!
*/
  r = ( ( double ) ( *seed ) ) * 4.656612875E-10;

  return r;
}
/******************************************************************************/

double *r83_cr_fa ( int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R83_CR_FA decomposes a real tridiagonal matrix using cyclic reduction.

  Discussion:

    The R83 storage format is used for a tridiagonal matrix.
    The superdiagonal is stored in entries (1,2:N), the diagonal in
    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
    original matrix is "collapsed" vertically into the array.

    Once R83_CR_FA has decomposed a matrix A, then R83_CR_SL may be used to solve
    linear systems A * x = b.

    R83_CR_FA does not employ pivoting.  Hence, the results can be more
    sensitive to ill-conditioning than standard Gauss elimination.  In
    particular, R83_CR_FA will fail if any diagonal element of the matrix
    is zero.  Other matrices may also cause R83_CR_FA to fail.

    R83_CR_FA can be guaranteed to work properly if the matrix is strictly
    diagonally dominant, that is, if the absolute value of the diagonal
    element is strictly greater than the sum of the absolute values of
    the offdiagonal elements, for each equation.

    The algorithm may be illustrated by the following figures:

    The initial matrix is given by:

          D1 U1
          L1 D2 U2
             L2 D3 U3
                L3 D4 U4
                   L4 D U5
                      L5 D6

    Rows and columns are permuted in an odd/even way to yield:

          D1       U1
             D3    L2 U3
                D5    L4 U5
          L1 U2    D2
             L3 U4    D4
                L5       D6

    A block LU decomposition is performed to yield:

          D1      |U1
             D3   |L2 U3
                D5|   L4 U5
          --------+--------
                  |D2'F3
                  |F1 D4'F4
                  |   F2 D6'

    For large systems, this reduction is repeated on the lower right hand
    tridiagonal subsystem until a completely upper triangular system
    is obtained.  The system has now been factored into the product of a
    lower triangular system and an upper triangular one, and the information
    defining this factorization may be used by R83_CR_SL to solve linear
    systems.

  Example:

    Here is how a R83 matrix of order 5 would be stored:

       *  A12 A23 A34 A45
      A11 A22 A33 A44 A55
      A21 A32 A43 A54  *

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    30 May 2009

  Author:

    C version by John Burkardt

  Reference:

    Roger Hockney,
    A fast direct solution of Poisson's equation using Fourier Analysis,
    Journal of the ACM,
    Volume 12, Number 1, pages 95-113, January 1965.

  Parameters:

    Input, int N, the order of the matrix.
    N must be positive.

    Input, double A[3*N], the R83 matrix.

    Output, double R83_CR_FA[3*(2*N+1)], factorization information 
    needed by R83_CR_SL.
*/
{
  double *a_cr;
  int iful;
  int ifulp;
  int ihaf;
  int il;
  int ilp;
  int inc;
  int incr;
  int ipnt;
  int ipntp;
  int j;

  if ( n <= 0 )
  {
    printf ( "\n" );
    printf ( "R83_CR_FA - Fatal error!\n" );
    printf ( "  Nonpositive N = %d\n", n );
    return NULL;
  }

  a_cr = ( double * ) malloc ( 3 * ( 2 * n + 1 ) * sizeof ( double ) );

  if ( n == 1 )
  {
    a_cr[0+0*3] = 0.0;
    a_cr[0+1*3] = 0.0;
    a_cr[0+2*3] = 0.0;
    a_cr[1+0*3] = 0.0;
    a_cr[1+1*3] = 1.0 / a[1+0*3];
    a_cr[1+2*3] = 0.0;
    a_cr[2+0*3] = 0.0;
    a_cr[2+1*3] = 0.0;
    a_cr[2+2*3] = 0.0;

    return a_cr;
  }
/*
  Zero out the workspace entries.
*/
  a_cr[0+0*3] = 0.0;
  for ( j = 1; j <= n-1; j++ )
  {
    a_cr[0+j*3] = a[0+j*3];
  }
  for ( j = n; j <= 2*n; j++ )
  {
    a_cr[0+j*3] = 0.0;
  }

  a_cr[1+0*3] = 0.0;
  for ( j = 1; j <= n; j++ )
  {
    a_cr[1+j*3] = a[1+(j-1)*3];
  }
  for ( j = n+1; j <= 2*n; j++ )
  {
    a_cr[1+j*3] = 0.0;
  }
  a_cr[2+0*3] = 0.0;
  for ( j = 1; j <= n-1; j++ )
  {
    a_cr[2+j*3] = a[2+(j-1)*3];
  }
  for ( j = n; j <= 2*n; j++ )
  {
    a_cr[2+j*3] = 0.0;
  }

  il = n;
  ipntp = 0;

  while ( 1 < il )
  {
    ipnt = ipntp;
    ipntp = ipntp + il;
    if ( ( il % 2 ) == 1 )
    {
      inc = il + 1;
    }
    else
    {
      inc = il;
    }

    incr = inc / 2;
    il = il / 2;
    ihaf = ipntp + incr + 1;
    ifulp = ipnt + inc + 2;

    for ( ilp = incr; 1 <= ilp; ilp-- )
    {
      ifulp = ifulp - 2;
      iful = ifulp - 1;
      ihaf = ihaf - 1;

      a_cr[1+iful*3] = 1.0 / a_cr[1+iful*3];
      a_cr[2+iful*3]  = a_cr[2+iful*3]  * a_cr[1+iful*3];
      a_cr[0+ifulp*3] = a_cr[0+ifulp*3] * a_cr[1+(ifulp+1)*3];
      a_cr[1+ihaf*3]  = a_cr[1+ifulp*3] 
        - a_cr[0+iful*3]  * a_cr[2+iful*3]
        - a_cr[0+ifulp*3] * a_cr[2+ifulp*3];
      a_cr[2+ihaf*3] = -a_cr[2+ifulp*3] * a_cr[2+(ifulp+1)*3];
      a_cr[0+ihaf*3] = -a_cr[0+ifulp*3] * a_cr[0+(ifulp+1)*3];
    }
  }

  a_cr[1+(ipntp+1)*3] = 1.0 / a_cr[1+(ipntp+1)*3];

  return a_cr;
}
/******************************************************************************/

double *r83_cr_sl ( int n, double a_cr[], double b[] )

/******************************************************************************/
/*
  Purpose:

    R83_CR_SL solves a real linear system factored by R83_CR_FA.

  Discussion:

    The matrix A must be tridiagonal.  R83_CR_FA is called to compute the
    LU factors of A.  It does so using a form of cyclic reduction.  If
    the factors computed by R83_CR_FA are passed to R83_CR_SL, then one or 
    many linear systems involving the matrix A may be solved.

    Note that R83_CR_FA does not perform pivoting, and so the solution 
    produced by R83_CR_SL may be less accurate than a solution produced 
    by a standard Gauss algorithm.  However, such problems can be 
    guaranteed not to occur if the matrix A is strictly diagonally 
    dominant, that is, if the absolute value of the diagonal coefficient 
    is greater than the sum of the absolute values of the two off diagonal 
    coefficients, for each row of the matrix.

  Example:

    Here is how a R83 matrix of order 5 would be stored:

       *  A12 A23 A34 A45
      A11 A22 A33 A44 A55
      A21 A32 A43 A54  *

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    30 May 2009

  Author:

    C version by John Burkardt

  Reference:

    Roger Hockney,
    A fast direct solution of Poisson's equation using Fourier Analysis,
    Journal of the ACM,
    Volume 12, Number 1, pages 95-113, January 1965.

  Parameters:

    Input, int N, the order of the matrix.
    N must be positive.

    Input, double A_CR[3*(2*N+1)], factorization information computed by R83_CR_FA.

    Input, double B[N], the right hand side.

    Output, double R83_CR_SL[N], the solution.
*/
{
  int i;
  int iful;
  int ifulm;
  int ihaf;
  int il;
  int ipnt;
  int ipntp;
  int ndiv;
  double *rhs;
  double *x;

  if ( n <= 0 )
  {
    printf ( "\n" );
    printf ( "R83_CR_SL - Fatal error!\n" );
    printf ( "  Nonpositive N = %d\n", n );
    return NULL;
  }

  if ( n == 1 )
  {
    x = ( double * ) malloc ( 1 * sizeof ( double ) );
    x[0] = a_cr[1+1*3] * b[0];
    return x;
  }
/*
  Set up RHS.
*/
  rhs = ( double * ) malloc (  ( 2 * n + 1 ) * sizeof ( double ) );

  rhs[0] = 0.0;
  for ( i = 1; i <= n; i++ )
  {
    rhs[i] = b[i-1];
  }
  for ( i = n+1; i <= 2*n; i++ )
  {
    rhs[i] = 0.0;
  }

  il = n;
  ndiv = 1;
  ipntp = 0;

  while ( 1 < il )
  {
    ipnt = ipntp;
    ipntp = ipntp + il;
    il = il / 2;
    ndiv = ndiv * 2;
    ihaf = ipntp;

    for ( iful = ipnt + 2; iful <= ipntp; iful = iful + 2 )
    {
      ihaf = ihaf + 1;
      rhs[ihaf] = rhs[iful] 
        - a_cr[2+(iful-1)*3] * rhs[iful-1]
        - a_cr[0+iful*3]     * rhs[iful+1];
    }
  }

  rhs[ihaf] = rhs[ihaf] * a_cr[1+ihaf*3];

  ipnt = ipntp;

  while ( 0 < ipnt )
  {
    ipntp = ipnt;
    ndiv = ndiv / 2;
    il = n / ndiv;
    ipnt = ipnt - il;
    ihaf = ipntp;

    for ( ifulm = ipnt + 1; ifulm <= ipntp; ifulm = ifulm + 2 )
    {
      iful = ifulm + 1;
      ihaf = ihaf + 1;
      rhs[iful] = rhs[ihaf];
      rhs[ifulm] = a_cr[1+ifulm*3] * ( 
                               rhs[ifulm] 
        - a_cr[2+(ifulm-1)*3] * rhs[ifulm-1] 
        - a_cr[0+ifulm*3]     * rhs[iful] );
    }
  }

  x = ( double * ) malloc ( n * sizeof ( double ) );

  for ( i = 0; i < n; i++ )
  {
    x[i] = rhs[i+1];
  }

  free ( rhs );

  return x;
}
/******************************************************************************/

double *r83_cr_sls ( int n, double a_cr[], int nb, double b[] )

/******************************************************************************/
/*
  Purpose:

    R83_CR_SLS solves several real linear systems factored by R83_CR_FA.

  Discussion:

    The matrix A must be tridiagonal.  R83_CR_FA is called to compute the
    LU factors of A.  It does so using a form of cyclic reduction.  If
    the factors computed by R83_CR_FA are passed to R83_CR_SLS, then one or 
    many linear systems involving the matrix A may be solved.

    Note that R83_CR_FA does not perform pivoting, and so the solutions
    produced by R83_CR_SLS may be less accurate than a solution produced 
    by a standard Gauss algorithm.  However, such problems can be 
    guaranteed not to occur if the matrix A is strictly diagonally 
    dominant, that is, if the absolute value of the diagonal coefficient 
    is greater than the sum of the absolute values of the two off diagonal 
    coefficients, for each row of the matrix.

  Example:

    Here is how a R83 matrix of order 5 would be stored:

       *  A12 A23 A34 A45
      A11 A22 A33 A44 A55
      A21 A32 A43 A54  *

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    09 May 2010

  Author:

    John Burkardt

  Reference:

    Roger Hockney,
    A fast direct solution of Poisson's equation using Fourier Analysis,
    Journal of the ACM,
    Volume 12, Number 1, pages 95-113, January 1965.

  Parameters:

    Input, int N, the order of the matrix.
    N must be positive.

    Input, double A_CR[3*(2*N+1)], factorization information computed by R83_CR_FA.

    Input, int NB, the number of systems.

    Input, double B[N*NB], the right hand sides.

    Output, double R83_CR_SL[N*NB], the solutions.
*/
{
  int i;
  int iful;
  int ifulm;
  int ihaf;
  int il;
  int ipnt;
  int ipntp;
  int j;
  int ndiv;
  double *rhs;
  double *x;

  if ( n <= 0 )
  {
    printf ( "\n" );
    printf ( "R83_CR_SLS - Fatal error!\n" );
    printf ( "  Nonpositive N = %d\n", n );
    exit ( 1 );
  }

  if ( n == 1 )
  {
    x = ( double * ) malloc ( n * nb * sizeof ( double ) );
    for ( j = 0; j < nb; j++ )
    {
      x[0+j*n] = a_cr[1+1*3] * b[0+j*n];
    }
    return x;
  }
//
//  Set up RHS.
//
  rhs = ( double * ) malloc ( ( 2 * n + 1 ) * nb * sizeof ( double ) );

  for ( j = 0; j < nb; j++ )
  {
    rhs[0+j*n] = 0.0;
    for ( i = 1; i <= n; i++ )
    {
      rhs[i+j*n] = b[i-1+j*n];
    }
    for ( i = n + 1; i <= 2 * n; i++ )
    {
      rhs[i+j*n] = 0.0;
    }
  }

  il = n;
  ndiv = 1;
  ipntp = 0;

  while ( 1 < il )
  {
    ipnt = ipntp;
    ipntp = ipntp + il;
    il = il / 2;
    ndiv = ndiv * 2;

    for ( j = 0; j < nb; j++ )
    {
      ihaf = ipntp;
      for ( iful = ipnt + 2; iful <= ipntp; iful = iful + 2 )
      {
        ihaf = ihaf + 1;
        rhs[ihaf+j*n] = rhs[iful+j*n] 
          - a_cr[2+(iful-1)*3] * rhs[iful-1+j*n]
          - a_cr[0+iful*3]     * rhs[iful+1+j*n];
      }
    }
  }

  for ( j = 0; j < nb; j++ )
  {
    rhs[ihaf+j*n] = rhs[ihaf+j*n] * a_cr[1+ihaf*3];
  }

  ipnt = ipntp;

  while ( 0 < ipnt )
  {
    ipntp = ipnt;
    ndiv = ndiv / 2;
    il = n / ndiv;
    ipnt = ipnt - il;

    for ( j = 0; j < nb; j++ )
    {
      ihaf = ipntp;
      for ( ifulm = ipnt + 1; ifulm <= ipntp; ifulm = ifulm + 2 )
      {
        iful = ifulm + 1;
        ihaf = ihaf + 1;
        rhs[iful+j*n] = rhs[ihaf+j*n];
        rhs[ifulm+j*n] = a_cr[1+ifulm*3] * ( 
                                  rhs[ifulm+j*n] 
          - a_cr[2+(ifulm-1)*3] * rhs[ifulm-1+j*n] 
          - a_cr[0+ifulm*3]     * rhs[iful+j*n] );
      }
    }
  }

  x = ( double * ) malloc ( n * nb * sizeof ( double ) );

  for ( j = 0; j < nb; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      x[i+j*n] = rhs[i+1+j*n];
    }
  }

  free ( rhs );

  return x;
}
/******************************************************************************/

void r83_gs_sl ( int n, double a[], double b[], double x[], int it_max, 
  int job )

/******************************************************************************/
/*
  Purpose:

    R83_GS_SL solves a R83 system using Gauss-Seidel iteration.

  Discussion:

    The R83 storage format is used for a tridiagonal matrix.
    The superdiagonal is stored in entries (1,2:N), the diagonal in
    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
    original matrix is "collapsed" vertically into the array.

    This routine simply applies a given number of steps of the
    iteration to an input approximate solution.  On first call, you can
    simply pass in the zero vector as an approximate solution.  If
    the returned value is not acceptable, you may call again, using
    it as the starting point for additional iterations.

  Example:

    Here is how a R83 matrix of order 5 would be stored:

       *  A12 A23 A34 A45
      A11 A22 A33 A44 A55
      A21 A32 A43 A54  *

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    30 May 2009

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of the matrix.
    N must be at least 2.

    Input, double A[3*N], the R83 matrix.

    Input, double B[N], the right hand side of the linear system.

    Input/output, double X[N], an approximate solution to the system.

    Input, int IT_MAX, the maximum number of iterations to take.

    Input, int JOB, specifies the system to solve.
    0, solve A * x = b.
    nonzero, solve A' * x = b.
*/
{
  int i;
  int it_num;
/*
  No diagonal matrix entry can be zero.
*/
  for ( i = 0; i < n; i++ )
  {
    if ( a[1+i*3] == 0.0 )
    {
      printf ( "\n" );
      printf ( "R83_GS_SL - Fatal error!\n" );
      printf ( "  Zero diagonal entry, index = %d\n", i );
      return;
    }
  }

  if ( job == 0 )
  {
    for ( it_num = 1; it_num <= it_max; it_num++ )
    {
      x[0] =   ( b[0]                   - a[2+0*3] * x[1]     ) / a[1+0*3];
      for ( i = 1; i < n-1; i++ )
      {
        x[i] = ( b[i] - a[0+i*3] * x[i-1] - a[2+i*3] * x[i+1] ) / a[1+i*3];
      }
      x[n-1] =   ( b[n-1] - a[0+(n-1)*3] * x[n-2]             ) / a[1+(n-1)*3];
    }
  }
  else
  {
    for ( it_num = 1; it_num <= it_max; it_num++ )
    {
      x[0] =   ( b[0]                     - a[0+1*3] * x[1]     ) 
           / a[1+0*3];
      for ( i = 1; i < n-1; i++ )
      {
        x[i] = ( b[i] - a[2+(i-1)*3] * x[i-1] - a[0+(i+1)*3] * x[i+1] ) 
             / a[1+i*3];
      }
      x[n-1] =   ( b[n-1] - a[2+(n-2)*3] * x[n-2]                     ) 
             / a[1+(n-1)*3];
   
    }
  }

  return;
}
/******************************************************************************/

double *r83_indicator ( int n )

/******************************************************************************/
/*
  Purpose:

    R83_INDICATOR sets up a R83 indicator matrix.

  Discussion:

    The R83 storage format is used for a tridiagonal matrix.
    The superdiagonal is stored in entries (1,2:N), the diagonal in
    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
    original matrix is "collapsed" vertically into the array.

  Example:

    Here is how a R83 matrix of order 5 would be stored:

       *  A12 A23 A34 A45
      A11 A22 A33 A44 A55
      A21 A32 A43 A54  *

    Here are the values as stored in an indicator matrix:

      00 12 23 34 45
      11 22 33 44 55
      21 32 43 54 00

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    30 May 2009

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of the matrix.
    N must be at least 2.

    Output, double R83_INDICATOR[3*N], the R83 indicator matrix.
*/
{
  double *a;
  int fac;
  int i;
  int j;

  a = ( double * ) malloc ( 3 * n * sizeof ( double ) );

  fac = i4_power ( 10, i4_log_10 ( n ) + 1 );

  a[0+0*3] = 0.0;
  for ( j = 2; j <= n; j++ )
  {
    i = j - 1;
    a[0+(j-1)*3] = ( double ) ( fac * i + j );
  }

  for ( j = 1; j <= n; j++ )
  {
    i = j;
    a[1+(j-1)*3] = ( double ) ( fac * i + j );
  }

  for ( j = 1; j <= n-1; j++ )
  {
    i = j + 1;
    a[2+(j-1)*3] = ( double ) ( fac * i + j );
  }
  a[2+(n-1)*3] = 0.0;

  return a;
}
/******************************************************************************/

void r83_jac_sl ( int n, double a[], double b[], double x[], int it_max, 
  int job )

/******************************************************************************/
/*
  Purpose:

    R83_JAC_SL solves a R83 system using Jacobi iteration.

  Discussion:

    The R83 storage format is used for a tridiagonal matrix.
    The superdiagonal is stored in entries (1,2:N), the diagonal in
    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
    original matrix is "collapsed" vertically into the array.

    This routine simply applies a given number of steps of the
    iteration to an input approximate solution.  On first call, you can
    simply pass in the zero vector as an approximate solution.  If
    the returned value is not acceptable, you may call again, using
    it as the starting point for additional iterations.

  Example:

    Here is how a R83 matrix of order 5 would be stored:

       *  A12 A23 A34 A45
      A11 A22 A33 A44 A55
      A21 A32 A43 A54  *

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    30 May 2009

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of the matrix.
    N must be at least 2.

    Input, double A[3*N], the R83 matrix.

    Input, double B[N], the right hand side of the linear system.

    Input/output, double X[N], an approximate solution to the system.

    Input, int IT_MAX, the maximum number of iterations to take.

    Input, int JOB, specifies the system to solve.
    0, solve A * x = b.
    nonzero, solve A' * x = b.
*/
{
  int i;
  int it_num;
  double *xnew;

  xnew = ( double * ) malloc ( n * sizeof ( double ) );
/*
  No diagonal matrix entry can be zero.
*/
  for ( i = 0; i < n; i++ )
  {
    if ( a[1+i*3] == 0.0 )
    {
      printf ( "\n" );
      printf ( "R83_JAC_SL - Fatal error!\n" );
      printf ( "  Zero diagonal entry, index = %d\n", i );
      return;
    }
  }

  for ( it_num = 1; it_num <= it_max; it_num++ )
  {
/*
  Solve A*x=b:
*/
    if ( job == 0 )
    {
      xnew[0] =   b[0]                           - a[2+0*3] * x[1];
      for ( i = 1; i < n-1; i++ )
      {
        xnew[i] = b[i]   - a[0+i*3]     * x[i-1] - a[2+i*3] * x[i+1];
      }
      xnew[n-1] = b[n-1] - a[0+(n-1)*3] * x[n-2];
    }
/*
  Solve A'*x=b:
*/
    else
    {
      xnew[0] =   b[0]                     - a[0+1*3] * x[1];
      for ( i = 1; i < n-1; i++ )
      {
        xnew[i] = b[i] - a[2+(i-1)*3] * x[i-1] - a[0+(i+1)*3] * x[i+1];
      }
      xnew[n-1] =   b[n-1] - a[2+(n-2)*3] * x[n-2];
    }
/*
  Divide by the diagonal term, and overwrite X.
*/
    for ( i = 0; i < n; i++ )
    {
      xnew[i] = xnew[i] / a[1+i*3];
    }

    for ( i = 0; i < n; i++ )
    {
      x[i] = xnew[i];
    }
  }

  free ( xnew );

  return;
}
/******************************************************************************/

double *r83_mxv ( int n, double a[], double x[] )

/******************************************************************************/
/*
  Purpose:

    R83_MXV multiplies a R83 matrix times a vector.

  Discussion:

    The R83 storage format is used for a tridiagonal matrix.
    The superdiagonal is stored in entries (1,2:N), the diagonal in
    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
    original matrix is "collapsed" vertically into the array.

  Example:

    Here is how a R83 matrix of order 5 would be stored:

       *  A12 A23 A34 A45
      A11 A22 A33 A44 A55
      A21 A32 A43 A54  *

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    30 May 2009

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of the linear system.

    Input, double A[3*N], the R83 matrix.

    Input, double X[N], the vector to be multiplied by A.

    Output, double R83_MXV[N], the product A * x.
*/
{
  double *b;
  int i;

  b = ( double * ) malloc ( n * sizeof ( double ) );

  for ( i = 0; i < n; i++ )
  {
    b[i] =        a[1+i*3] * x[i];
  }
  for ( i = 0; i < n-1; i++ )
  {
    b[i] = b[i] + a[0+(i+1)*3] * x[i+1];
  }
  for ( i = 1; i < n; i++ )
  {
    b[i] = b[i] + a[2+(i-1)*3] * x[i-1];
  }

  return b;
}
/******************************************************************************/

double r83_np_det ( int n, double a_lu[] )

/******************************************************************************/
/*
  Purpose:

    R83_NP_DET: determinant of a tridiagonal system factored by R83_NP_FA.

  Discussion:

    The R83 storage format is used for a tridiagonal matrix.
    The superdiagonal is stored in entries (1,2:N), the diagonal in
    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
    original matrix is "collapsed" vertically into the array.

  Example:

    Here is how a R83 matrix of order 5 would be stored:

       *  A12 A23 A34 A45
      A11 A22 A33 A44 A55
      A21 A32 A43 A54  *

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    30 May 2009

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of the matrix.
    N must be at least 2.

    Input, double A_LU[3*N], the LU factors from R83_NP_FA.

    Output, double R83_NP_DET, the determinant of the matrix.
*/
{
  double det;
  int j;

  det = 1.0;
  for ( j = 0; j < n; j++ )
  {
    det = det * a_lu[1+j*3];
  }

  return det;
}
/******************************************************************************/

int r83_np_fa ( int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R83_NP_FA factors a R83 system without pivoting.

  Discussion:

    The R83 storage format is used for a tridiagonal matrix.
    The superdiagonal is stored in entries (1,2:N), the diagonal in
    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
    original matrix is "collapsed" vertically into the array.

    Because this routine does not use pivoting, it can fail even when
    the matrix is not singular, and it is liable to make larger
    errors.

    R83_NP_FA and R83_NP_SL may be preferable to the corresponding
    LINPACK routine SGTSL for tridiagonal systems, which factors and solves
    in one step, and does not save the factorization.

  Example:

    Here is how a R83 matrix of order 5 would be stored:

       *  A12 A23 A34 A45
      A11 A22 A33 A44 A55
      A21 A32 A43 A54  *

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    30 May 2009

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of the matrix.
    N must be at least 2.

    Input/output, double A[3*N].
    On input, the tridiagonal matrix.  On output, factorization information.

    Output, int R83_NP_FA, singularity flag.
    0, no singularity detected.
    nonzero, the factorization failed on the INFO-th step.
*/
{
  int i;

  for ( i = 1; i <= n-1; i++ )
  {
    if ( a[1+(i-1)*3] == 0.0 )
    {
      printf ( "\n" );
      printf ( "R83_NP_FA - Fatal error!\n" );
      printf ( "  Zero pivot on step %d\n", i );
      return i;
    }
/*
  Store the multiplier in L.
*/
    a[2+(i-1)*3] = a[2+(i-1)*3] / a[1+(i-1)*3];
/*
  Modify the diagonal entry in the next column.
*/
    a[1+i*3] = a[1+i*3] - a[2+(i-1)*3] * a[0+i*3];
  }

  if ( a[1+(n-1)*3] == 0.0 )
  {
    printf ( "\n" );
    printf ( "R83_NP_FA - Fatal error!\n" );
    printf ( "  Zero pivot on step %d\n", n );
    return n;
  }

  return 0;
}
/******************************************************************************/

double *r83_np_fs ( int n, double a[], double b[] )

/******************************************************************************/
/*
  Purpose:

    R83_NP_FS factors and solves a R83 system.

  Discussion:

    The R83 storage format is used for a tridiagonal matrix.
    The superdiagonal is stored in entries (1,2:N), the diagonal in
    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
    original matrix is "collapsed" vertically into the array.

    This algorithm requires that each diagonal entry be nonzero.
    It does not use pivoting, and so can fail on systems that
    are actually nonsingular.

  Example:

    Here is how a R83 matrix of order 5 would be stored:

       *  A12 A23 A34 A45
      A11 A22 A33 A44 A55
      A21 A32 A43 A54  *

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    30 May 2009

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of the linear system.

    Input/output, double A[3*N].
    On input, the nonzero diagonals of the linear system.
    On output, the data in these vectors has been overwritten
    by factorization information.

    Input, double B[N], the right hand side.

    Output, double R83_NP_FS[N], the solution of the linear system.
    This is NULL if there was an error because one of the diagonal
    entries was zero.
*/
{
  int i;
  double *x;
  double xmult;
/*
  Check.
*/
  for ( i = 0; i < n; i++ )
  {
    if ( a[1+i*3] == 0.0 )
    {
      return NULL;
    }
  }

  x = ( double * ) malloc ( n * sizeof ( double ) );

  for ( i = 0; i < n; i++ )
  {
    x[i] = b[i];
  }

  for ( i = 1; i < n; i++ )
  {
    xmult = a[2+(i-1)*3] / a[1+(i-1)*3];
    a[1+i*3] = a[1+i*3] - xmult * a[0+i*3];
    x[i] = x[i] - xmult * x[i-1];
  }

  x[n-1] = x[n-1] / a[1+(n-1)*3];
  for ( i = n-2; 0 <= i; i-- )
  {
    x[i] = ( x[i] - a[0+(i+1)*3] * x[i+1] ) / a[1+i*3];
  }

  return x;
}
/******************************************************************************/

double *r83_np_ml ( int n, double a_lu[], double x[], int job )

/******************************************************************************/
/*
  Purpose:

    R83_NP_ML computes Ax or xA, where A has been factored by R83_NP_FA.

  Discussion:

    The R83 storage format is used for a tridiagonal matrix.
    The superdiagonal is stored in entries (1,2:N), the diagonal in
    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
    original matrix is "collapsed" vertically into the array.

  Example:

    Here is how a R83 matrix of order 5 would be stored:

       *  A12 A23 A34 A45
      A11 A22 A33 A44 A55
      A21 A32 A43 A54  *

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    30 May 2009

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of the matrix.
    N must be at least 2.

    Input, double A_LU[3*N], the LU factors from R83_FA.

    Input, double X[N], the vector to be multiplied by A.

    Output, double B[N], the product.

    Input, int JOB, specifies the product to find.
    0, compute A * x.
    nonzero, compute A' * x.
*/
{
  double *b;
  int i;

  b = ( double * ) malloc ( n * sizeof ( double ) );

  for ( i = 0; i < n; i++ )
  {
    b[i] = x[i];
  }

  if ( job == 0 )
  {
/*
  Compute X := U * X
*/
    for ( i = 1; i <= n; i++ )
    {
      b[i-1] = a_lu[1+(i-1)*3] * b[i-1];

      if ( i < n )
      {
        b[i-1] = b[i-1] + a_lu[0+i*3] * b[i];
      }
    }
/*
  Compute X: = L * X.
*/
    for ( i = n; 2 <= i; i-- )
    {
      b[i-1] = b[i-1] + a_lu[2+(i-2)*3] * b[i-2];
    }
  }
  else
  {
/*
  Compute X: = L' * X.
*/
    for ( i = 1; i <= n-1; i++ )
    {
      b[i-1] = b[i-1] + a_lu[2+(i-1)*3] * b[i];
    }
/*
  Compute X: = U' * X.
*/
    for ( i = n; 1 <= i; i-- )
    {
      b[i-1] = a_lu[1+(i-1)*3] * b[i-1];
      if ( 1 < i )
      {
        b[i-1] = b[i-1] + a_lu[0+(i-1)*3] * b[i-2];
      }
    }
  }

  return b;
}
/******************************************************************************/

double *r83_np_sl ( int n, double a_lu[], double b[], int job )

/******************************************************************************/
/*
  Purpose:

    R83_NP_SL solves a R83 system factored by R83_NP_FA.

  Discussion:

    The R83 storage format is used for a tridiagonal matrix.
    The superdiagonal is stored in entries (1,2:N), the diagonal in
    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
    original matrix is "collapsed" vertically into the array.

  Example:

    Here is how a R83 matrix of order 5 would be stored:

       *  A12 A23 A34 A45
      A11 A22 A33 A44 A55
      A21 A32 A43 A54  *

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    30 May 2009

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of the matrix.
    N must be at least 2.

    Input, double A_LU[3*N], the LU factors from R83_NP_FA.

    Input, double B[N], the right hand side of the linear system.
    On output, B contains the solution of the linear system.

    Input, int JOB, specifies the system to solve.
    0, solve A * x = b.
    nonzero, solve A' * x = b.

    Output, double R83_NP_SL[N], the solution of the linear system.
*/
{
  int i;
  double *x;

  x = ( double * ) malloc ( n * sizeof ( double ) );

  for ( i = 0; i < n; i++ )
  {
    x[i] = b[i];
  }

  if ( job == 0 )
  {
/*
  Solve L * Y = B.
*/
    for ( i = 1; i < n; i++ )
    {
      x[i] = x[i] - a_lu[2+(i-1)*3] * x[i-1];
    }
/*
  Solve U * X = Y.
*/
    for ( i = n; 1 <= i; i-- )
    {
      x[i-1] = x[i-1] / a_lu[1+(i-1)*3];
      if ( 1 < i )
      {
        x[i-2] = x[i-2] - a_lu[0+(i-1)*3] * x[i-1];
      }
    }
  }
  else
  {
/*
  Solve U' * Y = B
*/
    for ( i = 1; i <= n; i++ )
    {
      x[i-1] = x[i-1] / a_lu[1+(i-1)*3];
      if ( i < n )
      {
        x[i] = x[i] - a_lu[0+i*3] * x[i-1];
      }
    }
/*
  Solve L' * X = Y.
*/
    for ( i = n-1; 1 <= i; i-- )
    {
      x[i-1] = x[i-1] - a_lu[2+(i-1)*3] * x[i];
    }
  }

  return x;
}
/******************************************************************************/

void r83_print ( int n, double a[], char *title )

/******************************************************************************/
/*
  Purpose:

    R83_PRINT prints a R83 matrix.

  Discussion:

    The R83 storage format is used for a tridiagonal matrix.
    The superdiagonal is stored in entries (1,2:N), the diagonal in
    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
    original matrix is "collapsed" vertically into the array.

  Example:

    Here is how a R83 matrix of order 5 would be stored:

       *  A12 A23 A34 A45
      A11 A22 A33 A44 A55
      A21 A32 A43 A54  *

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    30 May 2009

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of the matrix.
    N must be positive.

    Input, double A[3*N], the R83 matrix.

    Input, char *TITLE, a title.
*/
{
  r83_print_some ( n, a, 1, 1, n, n, title );

  return;
}
/******************************************************************************/

void r83_print_some ( int n, double a[], int ilo, int jlo, int ihi, int jhi,
  char *title )

/******************************************************************************/
/*
  Purpose:

    R83_PRINT_SOME prints some of a R83 matrix.

  Discussion:

    The R83 storage format is used for a tridiagonal matrix.
    The superdiagonal is stored in entries (1,2:N), the diagonal in
    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
    original matrix is "collapsed" vertically into the array.

  Example:

    Here is how a R83 matrix of order 5 would be stored:

       *  A12 A23 A34 A45
      A11 A22 A33 A44 A55
      A21 A32 A43 A54  *

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    30 May 2009

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of the matrix.
    N must be positive.

    Input, double A[3*N], the R83 matrix.

    Input, int ILO, JLO, IHI, JHI, designate the first row and
    column, and the last row and column, to be printed.

    Input, char *TITLE, a title.
*/
{
# define INCX 5

  int i;
  int i2hi;
  int i2lo;
  int inc;
  int j;
  int j2;
  int j2hi;
  int j2lo;

  printf ( "\n" );
  printf ( "%s\n", title );
/*
  Print the columns of the matrix, in strips of 5.
*/
  for ( j2lo = jlo; j2lo <= jhi; j2lo = j2lo + INCX )
  {
    j2hi = j2lo + INCX - 1;
    j2hi = i4_min ( j2hi, n );
    j2hi = i4_min ( j2hi, jhi );

    inc = j2hi + 1 - j2lo;

    printf ( "\n" );
    printf ( "  Col: " );
    for ( j = j2lo; j <= j2hi; j++ )
    {
      j2 = j + 1 - j2lo;
      printf ( "%7d       ", j );
    }

    printf ( "\n" );
    printf ( "  Row\n" );
    printf ( "  ---\n" );
/*
  Determine the range of the rows in this strip.
*/
    i2lo = i4_max ( ilo, 1 );
    i2lo = i4_max ( i2lo, j2lo - 1 );

    i2hi = i4_min ( ihi, n );
    i2hi = i4_min ( i2hi, j2hi + 1 );

    for ( i = i2lo; i <= i2hi; i++ )
    {
/*
  Print out (up to) 5 entries in row I, that lie in the current strip.
*/
      printf ( "%6d  ", i );

      for ( j2 = 1; j2 <= inc; j2++ )
      {
        j = j2lo - 1 + j2;

        if ( 1 < i-j || 1 < j-i )
        {
          printf ( "              " );
        }
        else if ( j == i+1 )
        {
          printf ( "%12f  ", a[0+(j-1)*3] );
        }
        else if ( j == i )
        {
          printf ( "%12f  ", a[1+(j-1)*3] );
        }
        else if ( j == i-1 )
        {
          printf ( "%12f  ", a[2+(j-1)*3] );
        }

      }
      printf ( "\n" );
    }
  }

  printf ( "\n" );

  return;
# undef INCX
}
/******************************************************************************/

double *r83_random ( int n, int *seed )

/******************************************************************************/
/*
  Purpose:

    R83_RANDOM randomizes a R83 matrix.

  Discussion:

    The R83 storage format is used for a tridiagonal matrix.
    The superdiagonal is stored in entries (1,2:N), the diagonal in
    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
    original matrix is "collapsed" vertically into the array.

  Example:

    Here is how a R83 matrix of order 5 would be stored:

       *  A12 A23 A34 A45
      A11 A22 A33 A44 A55
      A21 A32 A43 A54  *

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    30 May 2009

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of the linear system.

    Input/output, int *SEED, a seed for the random number generator.

    Output, double R83_RANDOM[3*N], the R83 matrix.
*/
{
  double *a;
  int i;
  double *u;
  double *v;
  double *w;

  a = ( double * ) malloc ( 3 * n * sizeof ( double ) );

  u = r8vec_uniform ( n-1, 0.0, 1.0, seed );
  v = r8vec_uniform ( n,   0.0, 1.0, seed );
  w = r8vec_uniform ( n-1, 0.0, 1.0, seed );

  a[0+0*3] = 0.0;
  for ( i = 1; i < n; i++ )
  {
    a[0+i*3] = u[i-1];
  }
   for ( i = 0; i < n; i++ )
  {
    a[1+i*3] = v[i];
  }
  for ( i = 0; i < n-1; i++ )
  {
    a[2+i*3] = w[i];
  }
  a[2+(n-1)*3] = 0.0;

  free ( u );
  free ( v );
  free ( w );

  return a;
}
/******************************************************************************/

double *r83_to_r8ge ( int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R83_TO_R8GE copies a R83 matrix to a R8GE matrix.

  Discussion:

    The R83 storage format is used for a tridiagonal matrix.
    The superdiagonal is stored in entries (1,2:N), the diagonal in
    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
    original matrix is "collapsed" vertically into the array.

  Example:

    Here is how a R83 matrix of order 5 would be stored:

       *  A12 A23 A34 A45
      A11 A22 A33 A44 A55
      A21 A32 A43 A54  *

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    30 May 2009

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of the matrix.
    N must be at least 2.

    Input, double A[3*N], the R83 matrix.

    Output, double R83_TO_R8GE[N*N], the R8GE matrix.
*/
{
  double *b;
  int i;
  int j;

  b = ( double * ) malloc ( n * n * sizeof ( double ) );

  for ( i = 1; i <= n; i++ )
  {
    for ( j = 1; j <= n; j++ )
    {
      if ( j == i-1 )
      {
        b[i-1+(j-1)*n] = a[0+(i-1)*3];
      }
      else if ( i == j )
      {
        b[i-1+(j-1)*n] = a[1+(i-1)*3];
      }
      else if ( j == i+1 )
      {
        b[i-1+(j-1)*n] = a[2+(i-1)*3];
      }
      else
      {
        b[i-1+(j-1)*n] = 0.0;
      }
    }
  }

  return b;
}
/******************************************************************************/

double *r83_vxm ( int n, double a[], double x[] )

/******************************************************************************/
/*
  Purpose:

    R83_VXM multiplies a vector times a R83 matrix.

  Discussion:

    The R83 storage format is used for a tridiagonal matrix.
    The superdiagonal is stored in entries (1,2:N), the diagonal in
    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
    original matrix is "collapsed" vertically into the array.

  Example:

    Here is how a R83 matrix of order 5 would be stored:

       *  A12 A23 A34 A45
      A11 A22 A33 A44 A55
      A21 A32 A43 A54  *

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    30 May 2009

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of the linear system.

    Input, double A[3*N], the R83 matrix.

    Input, double X[N], the vector to be multiplied by A'.

    Output, double R83_VXM[N], the product A' * x.
*/
{
  double *b;
  int i;

  b = ( double * ) malloc ( n * sizeof ( double ) );

  for ( i = 1; i <= n; i++ )
  {
    b[i-1] = a[1+(i-1)*3] * x[i-1];
  }

  for ( i = 1; i <= n-1; i++ )
  {
    b[i-1] = b[i-1] + a[2+(i-1)*3] * x[i];
  }

  for ( i = 2; i <= n; i++ )
  {
    b[i-1] = b[i-1] + a[0+(i-1)*3] * x[i-2];
  }

  return b;
}
/******************************************************************************/

double *r83_zero ( int n )

/******************************************************************************/
/*
  Purpose:

    R83_ZERO zeros a R83 matrix.

  Discussion:

    The R83 storage format is used for a tridiagonal matrix.
    The superdiagonal is stored in entries (1,2:N), the diagonal in
    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
    original matrix is "collapsed" vertically into the array.

  Example:

    Here is how a R83 matrix of order 5 would be stored:

       *  A12 A23 A34 A45
      A11 A22 A33 A44 A55
      A21 A32 A43 A54  *

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    30 May 2009

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of the matrix.
    N must be at least 2.

    Output, double R83_ZERO[3*N], the R83 matrix.
*/
{
  double *a;
  int i;
  int j;

  a = ( double * ) malloc ( 3 * n * sizeof ( double ) );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < 3; i++ )
    {
      a[i+j*3] = 0.0;
    }
  }

  return a;
}
/******************************************************************************/

double *r83np_fs ( int n, double a[], double b[] )

/******************************************************************************/
/*
  Purpose:

    R83NP_FS factors and solves an R83NP system.

  Discussion:

    The R83NP storage format is used for a tridiagonal matrix.
    The subdiagonal   is in entries (0,1:N-1), 
    the diagonal      is in entries (1,0:N-1), 
    the superdiagonal is in entries (2,0:N-2). 

    This algorithm requires that each diagonal entry be nonzero.
    It does not use pivoting, and so can fail on systems that
    are actually nonsingular.

    The "R83NP" format used for this routine is different from the R83 format.
    Here, we insist that the nonzero entries
    for a given row now appear in the corresponding column of the
    packed array.

  Example:

    Here is how a R83 matrix of order 5 would be stored:

       *  A21 A32 A43 A54
      A11 A22 A33 A44 A55
      A12 A23 A34 A45  *

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 May 2009

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of the linear system.

    Input/output, double A[3*N].
    On input, the nonzero diagonals of the linear system.
    On output, the data in these vectors has been overwritten
    by factorization information.

    Input, double B[N], the right hand side.

    Output, double R83NP_FS[N], the solution of the linear system.
*/
{
  int i;
  double *x;
/*
  Check.
*/
  for ( i = 0; i < n; i++ )
  {
    if ( a[1+i*3] == 0.0 )
    {
      printf ( "\n" );
      printf ( "R83NP_FS - Fatal error!\n" );
      printf ( "  A[1+%d*3] = 0.\n", i );
      exit ( 1 );
    }
  }

  x = ( double * ) malloc ( n * sizeof ( double ) );

  for ( i = 0; i < n; i++ )
  {
    x[i] = b[i];
  }

  for ( i = 1; i < n; i++ )
  {
    a[1+i*3] = a[1+i*3] - a[2+(i-1)*3] * a[0+i*3] / a[1+(i-1)*3];
    x[i]     = x[i]     - x[i-1]       * a[0+i*3] / a[1+(i-1)*3];
  }

  x[n-1] = x[n-1] / a[1+(n-1)*3];
  for ( i = n-2; 0 <= i; i-- )
  {
    x[i] = ( x[i] - a[2+i*3] * x[i+1] ) / a[1+i*3];
  }

  return x;
}
/******************************************************************************/

double r83p_det ( int n, double a_lu[], double work4 )

/******************************************************************************/
/*
  Purpose:

    R83P_DET computes the determinant of a matrix factored by R83P_FA.

  Discussion:

    The R83P storage format stores a periodic tridiagonal matrix is stored 
    as a 3 by N array, in which each row corresponds to a diagonal, and 
    column locations are preserved.  The matrix value 
    A(1,N) is stored as the array entry A(1,1), and the matrix value
    A(N,1) is stored as the array entry A(3,N).

  Example:

    Here is how a R83P matrix of order 5 would be stored:

      A51 A12 A23 A34 A45
      A11 A22 A33 A44 A55
      A21 A32 A43 A54 A15

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    13 June 2009

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of the matrix.
    N must be at least 3.

    Input, double A_LU[3*N], the LU factors from R83P_FA.

    Input, double WORK4, factorization information from R83P_FA.

    Output, double R83P_DET, the determinant of the matrix.
*/
{
  double det;
  int i;

  det = work4;
  for ( i = 0; i <= n-2; i++ )
  {
    det = det * a_lu[1+i*3];
  }

  return det;
}

/******************************************************************************/

int r83p_fa ( int n, double a[], double work2[], double work3[], double *work4 )

/******************************************************************************/
/*
  Purpose:

    R83P_FA factors a R83P matrix.

  Discussion:

    The R83P storage format stores a periodic tridiagonal matrix is stored as 
    a 3 by N array, in which each row corresponds to a diagonal, and 
    column locations are preserved.  The matrix value 
    A(1,N) is stored as the array entry A(1,1), and the matrix value
    A(N,1) is stored as the array entry A(3,N).

    Once the matrix has been factored by R83P_FA, R83P_SL may be called
    to solve linear systems involving the matrix.

    The logical matrix has a form which is suggested by this diagram:

      D1 U1          L1
      L2 D2 U2
         L3 R83 U3
            L4 D4 U4
               L5 R85 U5
      U6          L6 D6

    The algorithm treats the matrix as a border banded matrix:

      ( A1  A2 )
      ( A3  A4 )

    where:

      D1 U1          | L1
      L2 D2 U2       |  0
         L3 R83 U3    |  0
            L4 D4 U4 |  0
               L5 R85 | U5
      ---------------+---
      U6  0  0  0 L6 | D6

  Example:

    Here is how a R83P matrix of order 5 would be stored:

      A51 A12 A23 A34 A45
      A11 A22 A33 A44 A55
      A21 A32 A43 A54 A15

  Method:

    The algorithm rewrites the system as:

         X1 + inverse(A1) A2 X2 = inverse(A1) B1

      A3 X1 +             A4 X2 = B2

    The first equation can be "solved" for X1 in terms of X2:

         X1 = - inverse(A1) A2 X2 + inverse(A1) B1

    allowing us to rewrite the second equation for X2 explicitly:

      ( A4 - A3 inverse(A1) A2 ) X2 = B2 - A3 inverse(A1) B1

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    13 June 2009

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of the matrix.
    N must be at least 3.

    Input/output, double A[3*N].
    On input, the periodic tridiagonal matrix.  
    On output, the arrays have been modified to hold information
    defining the border-banded factorization of submatrices A1
    and A3.

    Output, int R83P_FA, singularity flag.
    0, no singularity detected.
    nonzero, the factorization failed on the INFO-th step.

    Output, double WORK2[N-1], WORK3[N-1], *WORK4, factorization information.
*/
{
  int i;
  int info;
  int job;
  double *work1;

  work1 = ( double * ) malloc ( ( n - 2 ) * sizeof ( double ) );
/*
  Compute inverse(A1):
*/
  info = r83_np_fa ( n-1, a );

  if ( info != 0 )
  {
    printf ( "\n" );
    printf ( "R83P_FA - Fatal error!\n" );
    printf ( "  R83_NP_FA returned INFO = %d\n", info );
    printf ( "  Factoring failed for column INFO.\n" );
    printf ( "  The tridiagonal matrix A1 is singular.\n" );
    printf ( "  This algorithm cannot continue!\n" );
    return info;
  }
/*
  WORK2 := inverse(A1) * A2.
*/
  work2[0] = a[2+(n-1)*3];
  for ( i = 1; i < n-2; i++)
  {
    work2[i] = 0.0;
  }
  work2[n-2] = a[0+(n-1)*3];

  job = 0;
  work1 = r83_np_sl ( n-1, a, work2, job );
  for ( i = 0; i < n-1; i++ )
  {
    work2[i] = work1[i];
  }
/*
  WORK3 := inverse ( A1' ) * A3'.
*/
  work3[0] = a[0+0*3];
  for ( i = 1; i < n-2; i++)
  {
    work3[i] = 0.0;
  }
  work3[n-2] = a[2+(n-2)*3];

  job = 1;
  work1 = r83_np_sl ( n-1, a, work3, job );
  for ( i = 0; i < n-1; i++ )
  {
    work3[i] = work1[i];
  }
/*
  A4 := ( A4 - A3 * inverse(A1) * A2 )
*/
  *work4 = a[1+(n-1)*3] - a[0+0*3] * work2[0] - a[2+(n-2)*3] * work2[n-2];

  if ( *work4 == 0.0 )
  {
    printf ( "\n" );
    printf ( "R83P_FA - Fatal error!\n" );
    printf ( "  The factored A4 submatrix is zero.\n" );
    printf ( "  This algorithm cannot continue!\n" );
    return n;
  }

  free ( work1 );

  return 0;
}
/******************************************************************************/

double *r83p_indicator ( int n )

/******************************************************************************/
/*
  Purpose:

    R83P_INDICATOR sets up a R83P indicator matrix.

  Discussion:

    The R83P storage format stores a periodic tridiagonal matrix is stored 
    as a 3 by N array, in which each row corresponds to a diagonal, and 
    column locations are preserved.  The matrix value 
    A(1,N) is stored as the array entry A(1,1), and the matrix value
    A(N,1) is stored as the array entry A(3,N).

  Example:

    Here is how a R83P matrix of order 5 would be stored:

      A51 A12 A23 A34 A45
      A11 A22 A33 A44 A55
      A21 A32 A43 A54 A15

    Here are the values as stored in an indicator matrix:

      51 12 23 34 45
      11 22 33 44 55
      21 32 43 54 15

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    13 June 2009

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of the matrix.
    N must be at least 2.

    Output, double R83P_INDICATOR[3*N], the R83P indicator matrix.
*/
{
  double *a;
  int fac;
  int i;
  int j;

  a = ( double * ) malloc ( 3 * n * sizeof ( double ) );

  fac = i4_power ( 10, i4_log_10 ( n ) + 1 );

  i = n;
  j = 1;
  a[0+(j-1)*3] = ( double ) ( fac * i + j );
  for ( j = 2; j <= n; j++ )
  {
    i = j - 1;
    a[0+(j-1)*3] = ( double ) ( fac * i + j );
  }

  for ( j = 1; j <= n; j++ )
  {
    i = j;
    a[1+(j-1)*3] = ( double ) ( fac * i + j );
  }

  for ( j = 1; j <= n-1; j++ )
  {
    i = j + 1;
    a[2+(j-1)*3] = ( double ) ( fac * i + j );
  }
  i = 1;
  j = n;
  a[2+(j-1)*3] = ( double ) ( fac * i + j );

  return a;
}
/******************************************************************************/

double *r83p_ml ( int n, double a_lu[], double x[], int job )

/******************************************************************************/
/*
  Purpose:

    R83P_ML computes A * x or x * A, where A has been factored by R83P_FA.

  Discussion:

    The R83P storage format stores a periodic tridiagonal matrix is stored 
    as a 3 by N array, in which each row corresponds to a diagonal, and 
    column locations are preserved.  The matrix value 
    A(1,N) is stored as the array entry A(1,1), and the matrix value
    A(N,1) is stored as the array entry A(3,N).

  Example:

    Here is how a R83P matrix of order 5 would be stored:

      A51 A12 A23 A34 A45
      A11 A22 A33 A44 A55
      A21 A32 A43 A54 A15

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    13 June 2009

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of the matrix.
    N must be at least 3.

    Input, double A_LU[3*N], the LU factors from R83P_FA.

    Input, double X[N], the vector to be multiplied by the matrix.

    Input, int JOB, indicates what product should be computed.
    0, compute A * x.
    nonzero, compute A' * x.

    Output, double R83P_ML[N], the result of the multiplication.
*/
{
  double *b;
  double *b_short;
  int i;
/*
  Multiply A(1:N-1,1:N-1) and X(1:N-1).
*/
  b_short = r83_np_ml ( n-1, a_lu, x, job );

  b = ( double * ) malloc ( n * sizeof ( double ) );

  for ( i = 0; i < n-1; i++ )
  {
    b[i] = b_short[i];
  }
  b[n-1] = 0.0;

  free ( b_short );
/*
  Add terms from the border.
*/
  if ( job == 0 )
  {
    b[0] = b[0] + a_lu[2+(n-1)*3] * x[n-1];
    b[n-2] = b[n-2] + a_lu[0+(n-1)*3] * x[n-1];
    b[n-1] = a_lu[0+0*3] * x[0] + a_lu[2+(n-2)*3] * x[n-2] 
      + a_lu[1+(n-1)*3] * x[n-1];
  }
  else
  {
    b[0] = b[0] + a_lu[0+0*3] * x[n-1];
    b[n-2] = b[n-2] + a_lu[2+(n-2)*3] * x[n-1];
    b[n-1] = a_lu[2+(n-1)*3] * x[0] + a_lu[0+(n-1)*3] * x[n-2] 
           + a_lu[1+(n-1)*3] * x[n-1];
  }

  return b;
}
/******************************************************************************/

double *r83p_mxv ( int n, double a[], double x[] )

/******************************************************************************/
/*
  Purpose:

    R83P_MXV multiplies a R83P matrix times a vector.

  Discussion:

    The R83P storage format stores a periodic tridiagonal matrix is stored as 
    a 3 by N array, in which each row corresponds to a diagonal, and 
    column locations are preserved.  The matrix value 
    A(1,N) is stored as the array entry A(1,1), and the matrix value
    A(N,1) is stored as the array entry A(3,N).

  Example:

    Here is how a R83P matrix of order 5 would be stored:

      A51 A12 A23 A34 A45
      A11 A22 A33 A44 A55
      A21 A32 A43 A54 A15

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    13 June 2009

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of the matrix.
    N must be at least 3.

    Input, double A[3*N], the R83P matrix.

    Input, double X[N], the vector to be multiplied by A.

    Output, double R83P_MXV[N], the product A * x.
*/
{
  double *b;
  int i;

  b = ( double * ) malloc ( n * sizeof ( double ) );

  b[0] =   a[2+(n-1)*3] * x[n-1] + a[1+0*3]     * x[0]   + a[0+1*3]     * x[1];

  for ( i = 1; i < n-1; i++ )
  {
    b[i] = a[2+(i-1)*3] * x[i-1] + a[1+i*3]     * x[i]   + a[0+(i+1)*3] * x[i+1];
  }

  b[n-1] = a[2+(n-2)*3] * x[n-2] + a[1+(n-1)*3] * x[n-1] + a[0+0*3]     * x[0];

  return b;
}
/******************************************************************************/

void r83p_print ( int n, double a[], char *title )

/******************************************************************************/
/*
  Purpose:

    R83P_PRINT prints a R83P matrix.

  Discussion:

    The R83P storage format stores a periodic tridiagonal matrix is stored as 
    a 3 by N array, in which each row corresponds to a diagonal, and 
    column locations are preserved.  The matrix value 
    A(1,N) is stored as the array entry A(1,1), and the matrix value
    A(N,1) is stored as the array entry A(3,N).

  Example:

    Here is how a R83P matrix of order 5 would be stored:

      A51 A12 A23 A34 A45
      A11 A22 A33 A44 A55
      A21 A32 A43 A54 A15

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    13 June 2009

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of the matrix.
    N must be positive.

    Input, double A[3*N], the R83P matrix.

    Input, char *TITLE, a title.
*/
{
  r83p_print_some ( n, a, 1, 1, n, n, title );

  return;
}
/******************************************************************************/

void r83p_print_some ( int n, double a[], int ilo, int jlo, int ihi, int jhi,
  char *title )

/******************************************************************************/
/*
  Purpose:

    R83P_PRINT_SOME prints some of a R83P matrix.

  Discussion:

    The R83P storage format stores a periodic tridiagonal matrix is stored as 
    a 3 by N array, in which each row corresponds to a diagonal, and 
    column locations are preserved.  The matrix value 
    A(1,N) is stored as the array entry A(1,1), and the matrix value
    A(N,1) is stored as the array entry A(3,N).

  Example:

    Here is how a R83P matrix of order 5 would be stored:

      A51 A12 A23 A34 A45
      A11 A22 A33 A44 A55
      A21 A32 A43 A54 A15

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    13 June 2009

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of the matrix.
    N must be positive.

    Input, double A[3*N], the R83P matrix.

    Input, int ILO, JLO, IHI, JHI, designate the first row and
    column, and the last row and column, to be printed.

    Input, char *TITLE, a title.
*/
{
# define INCX 5

  int i;
  int i2hi;
  int i2lo;
  int inc;
  int j;
  int j2;
  int j2hi;
  int j2lo;

  printf ( "\n" );
  printf ( "%s\n", title );
/*
  Print the columns of the matrix, in strips of 5.
*/
  for ( j2lo = jlo; j2lo <= jhi; j2lo = j2lo + INCX )
  {
    j2hi = j2lo + INCX - 1;
    j2hi = i4_min ( j2hi, n );
    j2hi = i4_min ( j2hi, jhi );

    inc = j2hi + 1 - j2lo;

    printf ( "\n" );
    printf ( "  Col: " );
    for ( j = j2lo; j <= j2hi; j++ )
    {
      printf ( "%7d       ", j );
    }
    printf ( "\n" );
    printf ( "  Row\n" );
    printf ( "  ---\n" );
/*
  Determine the range of the rows in this strip.
*/
    i2lo = i4_max ( ilo, 1 );

    if ( 1 < i2lo || j2hi < n )
    {
      i2lo = i4_max ( i2lo, j2lo - 1 );
    }

    i2hi = i4_min ( ihi, n );

    if ( i2hi < n || 1 < j2lo )
    {
      i2hi = i4_min ( i2hi, j2hi + 1 );
    }

    for ( i = i2lo; i <= i2hi; i++ )
    {
/*
  Print out (up to) 5 entries in row I, that lie in the current strip.
*/
      printf ( "%4d  ", i );

      for ( j2 = 1; j2 <= inc; j2++ )
      {
        j = j2lo - 1 + j2;

        if ( i == n && j == 1 )
        {
          printf ( "%12f  ", a[0+(j-1)*3] );
        }
        else if ( i == 1 && j == n )
        {
          printf ( "%12f  ", a[2+(j-1)*3] );
        }
        else if ( 1 < i-j || 1 < j-i )
        {
          printf ( "              " );
        }
        else if ( j == i+1 )
        {
          printf ( "%12f  ", a[0+(j-1)*3] );
        }
        else if ( j == i )
        {
          printf ( "%12f  ", a[1+(j-1)*3] );
        }
        else if ( j == i-1 )
        {
          printf ( "%12f  ", a[2+(j-1)*3] );
        }
      }
      printf ( "\n" );
    }
  }

  printf ( "\n" );

  return;
# undef INCX
}
/******************************************************************************/

double *r83p_random ( int n, int *seed )

/******************************************************************************/
/*
  Purpose:

    R83P_RANDOM randomizes a R83P matrix.

  Discussion:

    The R83P storage format stores a periodic tridiagonal matrix is stored as 
    a 3 by N array, in which each row corresponds to a diagonal, and 
    column locations are preserved.  The matrix value 
    A(1,N) is stored as the array entry A(1,1), and the matrix value
    A(N,1) is stored as the array entry A(3,N).

  Example:

    Here is how a R83P matrix of order 5 would be stored:

      A51 A12 A23 A34 A45
      A11 A22 A33 A44 A55
      A21 A32 A43 A54 A15

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    13 June 2009

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of the matrix.
    N must be at least 3.

    Input/output, int *SEED, a seed for the random number generator.

    Output, double R83P_RANDOM[3*N], the R83P matrix.
*/
{
  double *a;
  int i;
  int j;

  a = ( double * ) malloc ( 3 * n * sizeof ( double ) );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < 3; i++ )
    {
      a[i+j*3] = r8_uniform_01 ( seed );
    }
  }
  return a;
}
/******************************************************************************/

double *r83p_sl ( int n, double a_lu[], double b[], int job, double work2[], 
  double work3[], double work4 )

/******************************************************************************/
/*
  Purpose:

    R83P_SL solves a R83P system factored by R83P_FA.

  Discussion:

    The R83P storage format stores a periodic tridiagonal matrix is stored as 
    a 3 by N array, in which each row corresponds to a diagonal, and 
    column locations are preserved.  The matrix value 
    A(1,N) is stored as the array entry A(1,1), and the matrix value
    A(N,1) is stored as the array entry A(3,N).

  Example:

    Here is how a R83P matrix of order 5 would be stored:

      A51 A12 A23 A34 A45
      A11 A22 A33 A44 A55
      A21 A32 A43 A54 A15

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    13 June 2009

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of the matrix.
    N must be at least 3.

    Input, double A_LU[3*N], the LU factors from R83P_FA.

    Input, double B[N], the right hand side of the linear system.

    Input, int JOB, specifies the system to solve.
    0, solve A * x = b.
    nonzero, solve A' * x = b.

    Input, double WORK2(N-1), WORK3(N-1), WORK4, factor data from R83P_FA.

    Output, double R83P_SL[N], the solution to the linear system.
*/
{
  int i;
  double *x;
  double *xnm1;

  x = ( double * ) malloc ( n * sizeof ( double ) );

  for ( i = 0; i < n; i++ )
  {
    x[i] = b[i];
  }

  if ( job == 0 )
  {
/*
  Solve A1 * X1 = B1.
*/
    xnm1 = r83_np_sl ( n-1, a_lu, x, job );
/*
  X2 = B2 - A3 * X1
*/
    for ( i = 0; i < n-1; i++ )
    {
      x[i] = xnm1[i];
    }
    free ( xnm1 );

    x[n-1] = x[n-1] - a_lu[0+0*3] * x[0] - a_lu[2+(n-2)*3] * x[n-2];
/*
  Solve A4 * X2 = X2
*/
    x[n-1] = x[n-1] / work4;
/*
  X1 := X1 - inverse ( A1 ) * A2 * X2.
*/
    for ( i = 0; i < n-1; i++ )
    {
      x[i] = x[i] - work2[i] * x[n-1];
    }
  }
  else
  {
/*
  Solve A1' * X1 = B1.
*/
    xnm1 = r83_np_sl ( n-1, a_lu, x, job );
/*
  X2 := X2 - A2' * B1
*/
    for ( i = 0; i < n-1; i++ )
    {
      x[i] = xnm1[i];
    }
    free ( xnm1 );

    x[n-1] = x[n-1] - a_lu[2+(n-1)*3] * x[0] - a_lu[0+(n-1)*3] * x[n-2];
/*
  Solve A4 * X2 = X2.
*/
    x[n-1] = x[n-1] / work4;
/*
  X1 := X1 - transpose ( inverse ( A1 ) * A3 ) * X2.
*/
    for ( i = 0; i < n-1; i++ )
    {
      x[i] = x[i] - work3[i] * x[n-1];
    }
  }
  return x;
}
/******************************************************************************/

double *r83p_to_r8ge ( int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R83P_TO_R8GE copies a R83P matrix to a R8GE matrix.

  Discussion:

    The R83P storage format stores a periodic tridiagonal matrix is stored as 
    a 3 by N array, in which each row corresponds to a diagonal, and 
    column locations are preserved.  The matrix value 
    A(1,N) is stored as the array entry A(1,1), and the matrix value
    A(N,1) is stored as the array entry A(3,N).

  Example:

    Here is how a R83P matrix of order 5 would be stored:

      A51 A12 A23 A34 A45
      A11 A22 A33 A44 A55
      A21 A32 A43 A54 A15

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    13 June 2009

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of the matrix.
    N must be at least 3.

    Input, double A[3*N], the R83P matrix.

    Output, double R83P_TO_R8GE[N*N], the R8GE matrix.
*/
{
  double *b;
  int i;
  int j;

  b = ( double * ) malloc ( n * n * sizeof ( double ) );

  for ( i = 1; i <= n; i++ )
  {
    for ( j = 1; j <= n; j++ )
    {
      if ( i == j )
      {
        b[i-1+(j-1)*n] = a[1+(j-1)*3];
      }
      else if ( j == i-1 )
      {
        b[i-1+(j-1)*n] = a[2+(j-1)*3];
      }
      else if ( j == i+1 )
      {
        b[i-1+(j-1)*n] = a[0+(j-1)*3];
      }
      else if ( i == 1 && j == n )
      {
        b[i-1+(j-1)*n] = a[2+(j-1)*3];
      }
      else if ( i == n && j == 1 )
      {
        b[i-1+(j-1)*n] = a[0+(j-1)*3];
      }
      else
      {
        b[i-1+(j-1)*n] = 0.0;
      }
    }
  }

  return b;
}
/******************************************************************************/

double *r83p_vxm ( int n, double a[], double x[] )

/******************************************************************************/
/*
  Purpose:

    R83P_VXM multiplies a vector times a R83P matrix.

  Discussion:

    The R83P storage format stores a periodic tridiagonal matrix is stored as 
    a 3 by N array, in which each row corresponds to a diagonal, and 
    column locations are preserved.  The matrix value 
    A(1,N) is stored as the array entry A(1,1), and the matrix value
    A(N,1) is stored as the array entry A(3,N).

  Example:

    Here is how a R83P matrix of order 5 would be stored:

      A51 A12 A23 A34 A45
      A11 A22 A33 A44 A55
      A21 A32 A43 A54 A15

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    13 June 2009

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of the matrix.
    N must be at least 3.

    Input, double A[3*N], the R83P matrix.

    Input, double X, the vector to be multiplied by A.

    Output, double R83P_VXM[N], the product X * A.
*/
{
  double *b;
  int i;

  b = ( double * ) malloc ( n * sizeof ( double ) );

  b[0] = a[0+0*3] * x[n-1] + a[1+0*3] * x[0] + a[2+0*3] * x[1];

  for ( i = 2; i <= n-1; i++ )
  {
    b[i-1] = a[0+(i-1)*3] * x[i-2] + a[1+(i-1)*3] * x[i-1] + a[2+(i-1)*3] * x[i];
  }

  b[n-1] = a[0+(n-1)*3] * x[n-2] + a[1+(n-1)*3] * x[n-1] + a[2+(n-1)*3] * x[0];

  return b;
}
/******************************************************************************/

double *r83p_zero ( int n )

/******************************************************************************/
/*
  Purpose:

    R83P_ZERO zeros a R83P matrix.

  Discussion:

    The R83P storage format stores a periodic tridiagonal matrix is stored as 
    a 3 by N array, in which each row corresponds to a diagonal, and 
    column locations are preserved.  The matrix value 
    A(1,N) is stored as the array entry A(1,1), and the matrix value
    A(N,1) is stored as the array entry A(3,N).

  Example:

    Here is how a R83P matrix of order 5 would be stored:

      A51 A12 A23 A34 A45
      A11 A22 A33 A44 A55
      A21 A32 A43 A54 A15

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    13 June 2009

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of the matrix.
    N must be at least 3.

    Output, double S3P[3*N], the R83P matrix.
*/
{
  double *a;
  int i;
  int j;

  a = ( double * ) malloc ( 3 * n * sizeof ( double ) );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < 3; i++ )
    {
      a[i+j*3] = 0.0;
    }
  }

  return a;
}
/******************************************************************************/

double *r85_indicator ( int n )

/******************************************************************************/
/*
  Purpose:

    R85_INDICATOR sets up a R85 indicator matrix.

  Discussion:

    The R85 storage format represents a pentadiagonal matrix as a 5
    by N array, in which each row corresponds to a diagonal, and
    column locations are preserved.  Thus, the original matrix is
    "collapsed" vertically into the array.

  Example:

    Here is how a R85 matrix of order 6 would be stored:

       *   *  A13 A24 A35 A46
       *  A12 A23 A34 A45 A56
      A11 A22 A33 A44 A55 A66
      A21 A32 A43 A54 A65  *
      A31 A42 A53 A64  *   *

    Here are the values as stored in an indicator matrix:

      00 00 13 24 35 46
      00 12 23 34 45 56
      11 22 33 44 55 66
      21 32 43 54 65 00
      31 42 53 64 00 00

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    05 September 2010

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of the matrix.
    N must be at least 2.

    Output, double R85_INDICATOR[3*N], the R85 indicator matrix.
*/
{
  double *a;
  int fac;
  int i;
  int j;

  a = ( double * ) malloc ( 5 * n * sizeof ( double ) );

  fac = i4_power ( 10, i4_log_10 ( n ) + 1 );

  a[0+0*5] = 0.0;
  a[0+1*5] = 0.0;
  for ( j = 3; j <= n; j++ )
  {
    i = j - 2;
    a[0+(j-1)*5] = ( double ) ( fac * i + j );
  }

  a[1+0*5] = 0.0;
  for ( j = 2; j <= n; j++ )
  {
    i = j - 1;
    a[1+(j-1)*5] = ( double ) ( fac * i + j );
  }

  for ( j = 1; j <= n; j++ )
  {
    i = j;
    a[2+(j-1)*5] = ( double ) ( fac * i + j );
  }

  for ( j = 1; j <= n-1; j++ )
  {
    i = j + 1;
    a[3+(j-1)*5] = ( double ) ( fac * i + j );
  }
  a[3+(n-1)*5] = 0.0;

  for ( j = 1; j <= n-2; j++ )
  {
    i = j + 2;
    a[4+(j-1)*5] = ( double ) ( fac * i + j );
  }
  a[4+(n-2)*5] = 0.0;
  a[4+(n-1)*5] = 0.0;

  return a;
}
/******************************************************************************/

double *r85_np_fs ( int n, double a[], double b[] )

/******************************************************************************/
/*
  Purpose:

    R85_NP_FS factors and solves a R85 system.

  Discussion:

    The R85 storage format represents a pentadiagonal matrix as a 5
    by N array, in which each row corresponds to a diagonal, and
    column locations are preserved.  Thus, the original matrix is
    "collapsed" vertically into the array.

    This algorithm requires that each diagonal entry be nonzero.

    No pivoting is performed, and therefore the algorithm may fail
    in simple cases where the matrix is not singular.

  Example:

    Here is how a R85 matrix of order 6 would be stored:

       *   *  A13 A24 A35 A46
       *  A12 A23 A34 A45 A56
      A11 A22 A33 A44 A55 A66
      A21 A32 A43 A54 A65  *
      A31 A42 A53 A64  *   *

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    05 September 2010

  Author:

    Original FORTRAN77 version by Cheney, Kincaid.
    C version by John Burkardt.

  Reference:

    Ward Cheney, David Kincaid,
    Numerical Mathematics and Computing,
    1985, pages 233-236.

  Parameters:

    Input, int N, the order of the linear system.

    Input/output, double A[5*N],
    On input, the pentadiagonal matrix.
    On output, the data in these vectors has been overwritten
    by factorization information.

    Input/output, double B[N].
    On input, B contains the right hand side of the linear system.
    On output, B has been overwritten by factorization information.

    Output, double R85_NP_FS[N], the solution of the linear system.
*/
{
  int i;
  double *x;
  double xmult;

  for ( i = 0; i < n; i++ )
  {
    if ( a[2+i*5] == 0.0 )
    {
      return NULL;
    }
  }

  x = ( double * ) malloc ( n * sizeof ( double ) );

  for ( i = 2; i <= n-1; i++ )
  {
    xmult = a[1+(i-1)*5] / a[2+(i-2)*5];
    a[2+(i-1)*5] = a[2+(i-1)*5] - xmult * a[3+(i-2)*5];
    a[3+(i-1)*5] = a[3+(i-1)*5] - xmult * a[4+(i-2)*5];

    b[i-1] = b[i-1] - xmult * b[i-2];

    xmult = a[0+i*5] / a[2+(i-2)*5];
    a[1+i*5] = a[1+i*5] - xmult * a[3+(i-2)*5];
    a[2+i*5] = a[2+i*5] - xmult * a[4+(i-2)*5];

    b[i] = b[i] - xmult * b[i-2];
  }

  xmult = a[1+(n-1)*5] / a[2+(n-2)*5];
  a[2+(n-1)*5] = a[2+(n-1)*5] - xmult * a[3+(n-2)*5];

  x[n-1] = ( b[n-1] - xmult * b[n-2] ) / a[2+(n-1)*5];
  x[n-2] = ( b[n-2] - a[3+(n-2)*5] * x[n-1] ) / a[2+(n-2)*5];

  for ( i = n - 2; 1 <= i; i-- )
  {
    x[i-1] = ( b[i-1] - a[3+(i-1)*5] * x[i] - a[4+(i-1)*5] * x[i+1] ) 
      / a[2+(i-1)*5];
  }

  return x;
}
/******************************************************************************/

double *r85_mxv ( int n, double a[], double x[] )

/******************************************************************************/
/*
  Purpose:

    R85_MXV multiplies a R85 matrix times a vector.

  Discussion:

    The R85 storage format represents a pentadiagonal matrix as a 5
    by N array, in which each row corresponds to a diagonal, and
    column locations are preserved.  Thus, the original matrix is
    "collapsed" vertically into the array.

  Example:

    Here is how a R85 matrix of order 6 would be stored:

       *   *  A13 A24 A35 A46
       *  A12 A23 A34 A45 A56
      A11 A22 A33 A44 A55 A66
      A21 A32 A43 A54 A65  *
      A31 A42 A53 A64  *   *

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    05 September 2010

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of the linear system.

    Input, double A[5*N], the pentadiagonal matrix.

    Input, double X[N], the vector to be multiplied by A.

    Output, double R85_MXV[N], the product A * x.
*/
{
  double *b;
  int i;

  b = ( double * ) malloc ( n * sizeof ( double ) );

  for ( i = 0; i < n; i++ )
  {
    b[i] = a[2+i*5] * x[i];
  }
  for ( i = 2; i < n; i++ )
  {
    b[i] = b[i] + a[0+i*5] * x[i-2];
  }
  for ( i = 1; i < n; i++ )
  {
    b[i] = b[i] + a[1+i*5] * x[i-1];
  }

  for ( i = 0; i < n-1; i++ )
  {
    b[i] = b[i] + a[3+i*5] * x[i+1];
  }
  for ( i = 0; i < n-2; i++ )
  {
    b[i] = b[i] + a[4+i*5] * x[i+2];
  }

  return b;
}
/******************************************************************************/

void r85_print ( int n, double a[], char *title )

/******************************************************************************/
/*
  Purpose:

    R85_PRINT prints a R85 matrix.

  Discussion:

    The R85 storage format represents a pentadiagonal matrix as a 5
    by N array, in which each row corresponds to a diagonal, and
    column locations are preserved.  Thus, the original matrix is
    "collapsed" vertically into the array.

  Example:

    Here is how a R85 matrix of order 6 would be stored:

       *   *  A13 A24 A35 A46
       *  A12 A23 A34 A45 A56
      A11 A22 A33 A44 A55 A66
      A21 A32 A43 A54 A65  *
      A31 A42 A53 A64  *   *

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    05 September 2010

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of the matrix.
    N must be positive.

    Input, double A[5*N], the pentadiagonal matrix.

    Input, char *TITLE, a title.
*/
{
  r85_print_some ( n, a, 1, 1, n, n, title );

  return;
}
/******************************************************************************/

void r85_print_some ( int n, double a[], int ilo, int jlo, int ihi, int jhi,
  char *title )

/******************************************************************************/
/*
  Purpose:

    R85_PRINT_SOME prints some of a R85 matrix.

  Discussion:

    The R85 storage format represents a pentadiagonal matrix as a 5
    by N array, in which each row corresponds to a diagonal, and
    column locations are preserved.  Thus, the original matrix is
    "collapsed" vertically into the array.

  Example:

    Here is how a R85 matrix of order 6 would be stored:

       *   *  A13 A24 A35 A46
       *  A12 A23 A34 A45 A56
      A11 A22 A33 A44 A55 A66
      A21 A32 A43 A54 A65  *
      A31 A42 A53 A64  *   *

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    05 September 2010

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of the matrix.
    N must be positive.

    Input, double A[5*N], the pentadiagonal matrix.

    Input, int ILO, JLO, IHI, JHI, designate the first row and
    column, and the last row and column, to be printed.

    Input, char *TITLE, a title.
*/
{
# define INCX 5

  int i;
  int i2hi;
  int i2lo;
  int inc;
  int j;
  int j2;
  int j2hi;
  int j2lo;

  printf ( "\n" );
  printf ( "%s\n", title );
/*
  Print the columns of the matrix, in strips of 5.
*/
  for ( j2lo = jlo; j2lo <= jhi; j2lo = j2lo + INCX )
  {
    j2hi = j2lo + INCX - 1;
    j2hi = i4_min ( j2hi, n );
    j2hi = i4_min ( j2hi, jhi );

    inc = j2hi + 1 - j2lo;

    printf ( "\n" );
    printf ( "  Col:  " );
    for ( j = j2lo; j <= j2hi; j++ )
    {
      printf ( "%7d       ", j );
    }
    printf ( "\n" );
    printf ( "  Row\n" );
    printf ( "  ---\n" );
/*
  Determine the range of the rows in this strip.
*/
    i2lo = i4_max ( ilo, 1 );
    i2lo = i4_max ( i2lo, j2lo - 2 );

    i2hi = i4_min ( ihi, n );
    i2hi = i4_min ( i2hi, j2hi + 2 );

    for ( i = i2lo; i <= i2hi; i++ )
    {
/*
  Print out (up to) 5 entries in row I, that lie in the current strip.
*/
      printf ( "%6d  ", i );

      for ( j2 = 1; j2 <= inc; j2++ )
      {
        j = j2lo - 1 + j2;

        if ( 2 < i-j || 2 < j-i )
        {
          printf ( "            " );
        }
        else if ( j == i+2 )
        {
          printf ( "%10f  ", a[0+(j-1)*5] );
        }
        else if ( j == i+1 )
        {
          printf ( "%10f  ", a[1+(j-1)*5] );
        }
        else if ( j == i )
        {
          printf ( "%10f  ", a[2+(j-1)*5] );
        }
        else if ( j == i-1 )
        {
          printf ( "%10f  ", a[3+(j-1)*5] );
        }
        else if ( j == i-2 )
        {
          printf ( "%10f  ", a[4+(j-1)*5] );
        }
      }
      printf ( "\n" );
    }
    printf ( "\n" );
  }
  printf ( "\n" );

  return;
# undef INCX
}
/******************************************************************************/

double *r85_random ( int n, int *seed )

/******************************************************************************/
/*
  Purpose:

    R85_RANDOM randomizes a R85 matrix.

  Discussion:

    The R85 storage format represents a pentadiagonal matrix as a 5
    by N array, in which each row corresponds to a diagonal, and
    column locations are preserved.  Thus, the original matrix is
    "collapsed" vertically into the array.

  Example:

    Here is how a R85 matrix of order 6 would be stored:

       *   *  A13 A24 A35 A46
       *  A12 A23 A34 A45 A56
      A11 A22 A33 A44 A55 A66
      A21 A32 A43 A54 A65  *
      A31 A42 A53 A64  *   *

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    05 September 2010

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of the linear system.

    Input/output, int *SEED, a seed for the random number generator.

    Output, double R85_RANDOM[5*N], the pentadiagonal matrix.
*/
{
  double *a;
  int i;
  int j;

  a = ( double * ) malloc ( 5  * n * sizeof ( double ) );

  i = 0;
  a[0+0*5]     = 0.0;
  a[0+1*5]     = 0.0;
  for ( j = 2; j < n; j++ )
  {
    a[i+j*5] = r8_uniform_01 ( seed );
  }

  i = 1;
  a[1+0*5]     = 0.0;
  for ( j = 1; j < n; j++ )
  {
    a[i+j*5] = r8_uniform_01 ( seed );
  }

  i = 2;
  for ( j = 0; j < n; j++ )
  {
    a[i+j*5] = r8_uniform_01 ( seed );
  }

  i = 3;
  for ( j = 0; j < n-1; j++ )
  {
    a[i+j*5] = r8_uniform_01 ( seed );
  }
  a[3+(n-1)*5] = 0.0;

  i = 4;
  for ( j = 0; j < n-2; j++ )
  {
    a[i+j*5] = r8_uniform_01 ( seed );
  }
  a[4+(n-2)*5] = 0.0;
  a[4+(n-1)*5] = 0.0;

  return a;
}
/******************************************************************************/

double *r85_to_r8ge ( int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R85_TO_R8GE copies a R85 matrix to a R8GE matrix.

  Discussion:

    The R85 storage format represents a pentadiagonal matrix as a 5
    by N array, in which each row corresponds to a diagonal, and
    column locations are preserved.  Thus, the original matrix is
    "collapsed" vertically into the array.

  Example:

    Here is how a R85 matrix of order 6 would be stored:

       *   *  A13 A24 A35 A46
       *  A12 A23 A34 A45 A56
      A11 A22 A33 A44 A55 A66
      A21 A32 A43 A54 A65  *
      A31 A42 A53 A64  *   *

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    05 September 2010

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of the matrix.
    N must be at least 3.

    Input, double A[5*N], the nonzero diagonals of the matrix.

    Output, double R85_TO_R8GE[N*N], the pentadiagonal matrix.
*/
{
  double *b;
  int i;
  int j;

  b = ( double * ) malloc ( n * n * sizeof ( double ) );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      if ( j == i-2 )
      {
        b[i+j*5] = a[0+i*5];
      }
      else if ( j == i-1 )
      {
        b[i+j*5] = a[1+i*5];
      }
      else if ( i == j )
      {
        b[i+j*5] = a[2+i*5];
      }
      else if ( j == i+1 )
      {
        b[i+j*5] = a[3+i*5];
      }
      else if ( j == i+2 )
      {
        b[i+j*5] = a[4+i*5];
      }
      else
      {
        b[i+j*5] = 0.0;
      }
    }
  }
  return b;
}
/******************************************************************************/

double *r85_vxm ( int n, double a[], double x[] )

/******************************************************************************/
/*
  Purpose:

    R85_VXM multiplies a vector times a R85 matrix.

  Discussion:

    The R85 storage format represents a pentadiagonal matrix as a 5
    by N array, in which each row corresponds to a diagonal, and
    column locations are preserved.  Thus, the original matrix is
    "collapsed" vertically into the array.

  Example:

    Here is how a R85 matrix of order 6 would be stored:

       *   *  A13 A24 A35 A46
       *  A12 A23 A34 A45 A56
      A11 A22 A33 A44 A55 A66
      A21 A32 A43 A54 A65  *
      A31 A42 A53 A64  *   *

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    05 September 2010

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of the linear system.

    Input, double A[5*N], the pentadiagonal matrix.

    Input, double X[N], the vector to be multiplied by A'.

    Output, double R85_VXM[N], the product A' * x.
*/
{
  double *b;
  int j;

  b = ( double * ) malloc ( n * sizeof ( double ) );

  for ( j = 0; j < n; j++ )
  {
    b[j] = a[2+j*5] * x[j];
  }

  for ( j = 1; j < n; j++ )
  {
    b[j] = b[j] + a[3+j*5] * x[j-1];
  }

  for ( j = 2; j < n; j++ )
  {
    b[j] = b[j] + a[4+j*5] * x[j-2];
  }

  for ( j = 0; j < n-1; j++ )
  {
    b[j] = b[j] + a[1+j*5] * x[j+1];
  }

  for ( j = 0; j < n-2; j++ )
  {
    b[j] = b[j] + a[0+j*5] * x[j+2];
  }

  return b;
}
/******************************************************************************/

double *r85_zero ( int n )

/******************************************************************************/
/*
  Purpose:

    R85_ZERO zeros a R85 matrix.

  Discussion:

    The R85 storage format represents a pentadiagonal matrix as a 5
    by N array, in which each row corresponds to a diagonal, and
    column locations are preserved.  Thus, the original matrix is
    "collapsed" vertically into the array.

  Example:

    Here is how a R85 matrix of order 6 would be stored:

       *   *  A13 A24 A35 A46
       *  A12 A23 A34 A45 A56
      A11 A22 A33 A44 A55 A66
      A21 A32 A43 A54 A65  *
      A31 A42 A53 A64  *   *

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    05 September 2010

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of the linear system.

    Output, double R85_ZERO[5*N], the R85 matrix.
*/
{
  double *a;
  int i;
  int j;

  a = ( double * ) malloc ( 5 * n * sizeof ( double ) );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < 5; i++ )
    {
      a[i+j*5] = 0.0;
    }
  }
  return a;
}
/******************************************************************************/

double r8ge_co ( int n, double a[], int pivot[] )

/******************************************************************************/
/*
  Purpose:

    R8GE_CO factors a R8GE matrix and estimates its condition number.

  Discussion:

    The R8GE storage format is used for a "general" M by N matrix.  
    A physical storage space is made for each logical entry.  The two 
    dimensional logical array is mapped to a vector, in which storage is 
    by columns.

    For the system A * X = B, relative perturbations in A and B
    of size EPSILON may cause relative perturbations in X of size
    EPSILON/RCOND.

    If RCOND is so small that the logical expression
      1.0 + rcond == 1.0
    is true, then A may be singular to working precision.  In particular,
    RCOND is zero if exact singularity is detected or the estimate
    underflows.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    10 February 2012

  Author:

    Original FORTRAN77 version by Dongarra, Bunch, Moler, Stewart.
    C version by John Burkardt.

  Reference:

    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
    LINPACK User's Guide,
    SIAM, 1979,
    ISBN13: 978-0-898711-72-1,
    LC: QA214.L56.

  Parameters:

    Input, int N, the order of the matrix A.

    Input/output, double A[N*N].  On input, a matrix to be factored.
    On output, the LU factorization of the matrix.

    Output, int PIVOT[N], the pivot indices.

    Output, double R8GE_CO, an estimate of the reciprocal condition number of A.
*/
{
  double anorm;
  double ek;
  int i;
  int info;
  int j;
  int k;
  int l;
  double rcond;
  double s;
  double sm;
  double t;
  double wk;
  double wkm;
  double ynorm;
  double *z;
/*
  Compute the L1 norm of A.
*/
  anorm = 0.0;
  for ( j = 0; j < n; j++ )
  {
    s = 0.0;
    for ( i = 0; i < n; i++ )
    {
      s = s + r8_abs( a[i+j*n] );
    }
    anorm = r8_max ( anorm, s );
  }
/*
  Compute the LU factorization.
*/
  info = r8ge_fa ( n, a, pivot );

  if ( info != 0 ) 
  {
    rcond = 0.0;
    return rcond;
  }
/*
  RCOND = 1 / ( norm(A) * (estimate of norm(inverse(A))) )

  estimate of norm(inverse(A)) = norm(Z) / norm(Y)

  where
    A * Z = Y
  and
    A' * Y = E

  The components of E are chosen to cause maximum local growth in the
  elements of W, where U'*W = E.  The vectors are frequently rescaled
  to avoid overflow.

  Solve U' * W = E.
*/
  ek = 1.0;
  z = ( double * ) malloc ( n * sizeof ( double ) );
  for ( i = 0; i < n; i++ )
  {
    z[i] = 0.0;
  }

  for ( k = 0; k < n; k++ )
  {
    if ( z[k] != 0.0 ) 
    {
      ek = - r8_sign2 ( ek, z[k] );
    }

    if ( r8_abs ( a[k+k*n] ) < r8_abs ( ek - z[k] ) )
    {
      s = r8_abs ( a[k+k*n] ) / r8_abs ( ek - z[k] );
      for ( i = 0; i < n; i++ )
      {
        z[i] = s * z[i];
      }
      ek = s * ek;
    }

    wk = ek - z[k];
    wkm = -ek - z[k];
    s = r8_abs ( wk );
    sm = r8_abs ( wkm );

    if ( a[k+k*n] != 0.0 )
    {
      wk = wk / a[k+k*n];
      wkm = wkm / a[k+k*n];
    }
    else
    {
      wk = 1.0;
      wkm = 1.0;
    }

    if ( k + 2 <= n )
    {
      for ( j = k+1; j < n; j++ )
      {
        sm = sm + r8_abs ( z[j] + wkm * a[k+j*n] );
        z[j] = z[j] + wk * a[k+j*n];
        s = s + r8_abs ( z[j] );
      }

      if ( s < sm )
      {
        t = wkm - wk;
        wk = wkm;
        for ( j = k+1; j < n; j++ )
        {
          z[j] = z[j] + t * a[k+j*n];
        }
      }
    }
    z[k] = wk;
  }

  s = 0.0;
  for ( i = 0; i < n; i++ )
  {
    s = s + r8_abs ( z[i] );
  }

  for ( i = 0; i < n; i++ )
  {
    z[i] = z[i] / s;
  }
/*
  Solve L' * Y = W
*/
  for ( k = n-1; 0 <= k; k-- )
  {
    for ( i = k+1; i < n; i++ )
    {
      z[k] = z[k] + z[i] * a[i+k*n];
    }
    t = r8_abs ( z[k] );
    if ( 1.0 < t )
    {
      for ( i = 0; i < n; i++ )
      {
        z[i] = z[i] / t;
      }
    }

    l = pivot[k] - 1;

    t    = z[l];
    z[l] = z[k];
    z[k] = t;
  }

  t = 0.0;
  for ( i = 0; i < n; i++ )
  {
    t = t + r8_abs ( z[i] );
  }
  for ( i = 0; i < n; i++ )
  {
    z[i] = z[i] / t;
  }

  ynorm = 1.0;
/*
  Solve L * V = Y.
*/
  for ( k = 0; k < n; k++ )
  {
    l = pivot[k] - 1;

    t    = z[l];
    z[l] = z[k];
    z[k] = t;

    for ( i = k+1; i < n; i++ )
    {
      z[i] = z[i] + t * a[i+k*n];
    }

    t = r8_abs ( z[k] );

    if ( 1.0 < t )
    {
      ynorm = ynorm / t;
      for ( i = 0; i < n; i++ )
      {
        z[i] = z[i] / t;
      }
    }
  }
  s = 0.0;
  for ( i = 0; i < n; i++ )
  {
    s = s + r8_abs ( z[i] );
  }

  for ( i = 0; i < n; i++ )
  {
    z[i] = z[i] / s;
  }
  ynorm = ynorm / s;
/*
  Solve U * Z = V.
*/
  for ( k = n-1; 0 <= k; k-- )
  {
    if ( r8_abs ( a[k+k*n] ) < r8_abs ( z[k] ) )
    {
      s = r8_abs ( a[k+k*n] ) / r8_abs ( z[k] );
      for ( i = 0; i < n; i++ )
      {
        z[i] = s * z[i];
      }
      ynorm = s * ynorm;
    }

    if ( a[k+k*n] != 0.0 )
    {
      z[k] = z[k] / a[k+k*n];
    }
    else
    {
      z[k] = 1.0;
    }

    for ( i = 0; i < k; i++ )
    {
      z[i] = z[i] - a[i+k*n] * z[k];
    }
  }
/*
  Normalize Z in the L1 norm.
*/
  s = 0.0;
  for ( i = 0; i < n; i++ )
  {
    s = s + r8_abs ( z[i] );
  }
  s = 1.0 / s;

  for ( i = 0; i < n; i++ )
  {
    z[i] = s * z[i];
  }
  ynorm = s * ynorm;

  if ( anorm != 0.0 )
  {
    rcond = ynorm / anorm;
  }
  else
  {
    rcond = 0.0;
  }

  free ( z );

  return rcond;
}
/******************************************************************************/

double r8ge_det ( int n, double a_lu[], int pivot[] )

/******************************************************************************/
/*
  Purpose:

    R8GE_DET computes the determinant of a matrix factored by R8GE_FA or R8GE_TRF.

  Discussion:

    The R8GE storage format is used for a "general" M by N matrix.  
    A physical storage space is made for each logical entry.  The two 
    dimensional logical array is mapped to a vector, in which storage is 
    by columns.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    10 February 2012

  Author:

    John Burkardt

  Reference:

    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
    LINPACK User's Guide,
    SIAM, 1979,
    ISBN13: 978-0-898711-72-1,
    LC: QA214.L56.

  Parameters:

    Input, int N, the order of the matrix.
    N must be positive.

    Input, double A_LU[N*N], the LU factors from R8GE_FA or R8GE_TRF.

    Input, int PIVOT[N], as computed by R8GE_FA or R8GE_TRF.

    Output, double R8GE_DET, the determinant of the matrix.
*/
{
  double det;
  int i;

  det = 1.0;

  for ( i = 1; i <= n; i++ )
  {
    det = det * a_lu[i-1+(i-1)*n];
    if ( pivot[i-1] != i )
    {
      det = -det;
    }
  }

  return det;
}
/******************************************************************************/

double *r8ge_dilu ( int m, int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8GE_DILU produces the diagonal incomplete LU factor of a R8GE matrix.

  Discussion:

    The R8GE storage format is used for a "general" M by N matrix.  
    A physical storage space is made for each logical entry.  The two 
    dimensional logical array is mapped to a vector, in which storage is 
    by columns.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    10 February 2012

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows in A.

    Input, int N, the number of columns in A.

    Input, double A[M*N], the R8GE matrix.

    Output, double R8GE_DILU[M], the D-ILU factor.
*/
{
  double *d;
  int i;
  int j;

  d = ( double * ) malloc ( m * sizeof ( double ) );

  for ( i = 0; i < m; i++ )
  {
    if ( i < n ) 
    {
      d[i] = a[i+i*m];
    }
    else
    {
      d[i] = 0.0;
    }
  }

  for ( i = 0; i < m && i < n; i++ )
  {
    d[i] = 1.0 / d[i];
    for ( j = i+1; j < m && j < n; j++ )
    {
      d[j] = d[j] - a[j+i*m] * d[i] * a[i+j*m];
    }
  }

  return d;
}
/******************************************************************************/

int r8ge_fa ( int n, double a[], int pivot[] )

/******************************************************************************/
/*
  Purpose:

    R8GE_FA performs a LINPACK-style PLU factorization of a R8GE matrix.

  Discussion:

    The R8GE storage format is used for a "general" M by N matrix.  
    A physical storage space is made for each logical entry.  The two 
    dimensional logical array is mapped to a vector, in which storage is 
    by columns.

    R8GE_FA is a simplified version of the LINPACK routine SGEFA.

    The two dimensional array is stored by columns in a one dimensional
    array.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    10 February 2012

  Author:

    John Burkardt

  Reference:

    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
    LINPACK User's Guide,
    SIAM, 1979,
    ISBN13: 978-0-898711-72-1,
    LC: QA214.L56.

  Parameters:

    Input, int N, the order of the matrix.
    N must be positive.

    Input/output, double A[N*N], the matrix to be factored.
    On output, A contains an upper triangular matrix and the multipliers
    which were used to obtain it.  The factorization can be written
    A = L * U, where L is a product of permutation and unit lower
    triangular matrices and U is upper triangular.

    Output, int PIVOT[N], a vector of pivot indices.

    Output, int R8GE_FA, singularity flag.
    0, no singularity detected.
    nonzero, the factorization failed on the INFO-th step.
*/
{
  int i;
  int j;
  int k;
  int l;
  double t;

  for ( k = 1; k <= n-1; k++ )
  {
/*
  Find L, the index of the pivot row.
*/
    l = k;

    for ( i = k+1; i <= n; i++ )
    {
      if ( r8_abs ( a[l-1+(k-1)*n] ) < r8_abs ( a[i-1+(k-1)*n] ) )
      {
        l = i;
      }
    }

    pivot[k-1] = l;
/*
  If the pivot index is zero, the algorithm has failed.
*/
    if ( a[l-1+(k-1)*n] == 0.0 )
    {
      fprintf ( stderr, "\n" );
      fprintf ( stderr, "R8GE_FA - Fatal error!\n" );
      fprintf ( stderr, "  Zero pivot on step %d\n", k );
      exit ( 1 );
    }
/*
  Interchange rows L and K if necessary.
*/
    if ( l != k )
    {
      t              = a[l-1+(k-1)*n];
      a[l-1+(k-1)*n] = a[k-1+(k-1)*n];
      a[k-1+(k-1)*n] = t;
    }
/*
  Normalize the values that lie below the pivot entry A(K,K).
*/
    for ( i = k+1; i <= n; i++ )
    {
      a[i-1+(k-1)*n] = -a[i-1+(k-1)*n] / a[k-1+(k-1)*n];
    }
/*
  Row elimination with column indexing.
*/
    for ( j = k+1; j <= n; j++ )
    {
      if ( l != k )
      {
        t              = a[l-1+(j-1)*n];
        a[l-1+(j-1)*n] = a[k-1+(j-1)*n];
        a[k-1+(j-1)*n] = t;
      }

      for ( i = k+1; i <= n; i++ )
      {
        a[i-1+(j-1)*n] = a[i-1+(j-1)*n] + a[i-1+(k-1)*n] * a[k-1+(j-1)*n];
      }

    }

  }

  pivot[n-1] = n;

  if ( a[n-1+(n-1)*n] == 0.0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8GE_FA - Fatal error!\n" );
    fprintf ( stderr, "  Zero pivot on step %d\n", n );
    exit ( 1 );
  }

  return 0;
}
/******************************************************************************/

void r8ge_fs ( int n, double a[], double x[] )

/******************************************************************************/
/*
  Purpose:

    R8GE_FS factors and solves a R8GE system.

  Discussion:

    The R8GE storage format is used for a "general" M by N matrix.  
    A physical storage space is made for each logical entry.  The two 
    dimensional logical array is mapped to a vector, in which storage is 
    by columns.

    The function does not save the LU factors of the matrix, and hence cannot
    be used to efficiently solve multiple linear systems, or even to
    factor A at one time, and solve a single linear system at a later time.

    The function uses partial pivoting, but no pivot vector is required.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    10 February 2012

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of the matrix.
    N must be positive.

    Input/output, double A[N*N].
    On input, A is the coefficient matrix of the linear system.
    On output, A is in unit upper triangular form, and
    represents the U factor of an LU factorization of the
    original coefficient matrix.

    Input/output, double X[N], on input, the right hand side of the linear system.
    On output, the solution of the linear system.
*/
{
  int i;
  int ipiv;
  int j;
  int jcol;
  double piv;
  double t;

  for ( jcol = 1; jcol <= n; jcol++ )
  {
/*
  Find the maximum element in column I.
*/
    piv = r8_abs ( a[jcol-1+(jcol-1)*n] );
    ipiv = jcol;
    for ( i = jcol+1; i <= n; i++ )
    {
      if ( piv < r8_abs ( a[i-1+(jcol-1)*n] ) )
      {
        piv = r8_abs ( a[i-1+(jcol-1)*n] );
        ipiv = i;
      }
    }

    if ( piv == 0.0 )
    {
      fprintf ( stderr, "\n" );
      fprintf ( stderr, "R8GE_FS - Fatal error!\n" );
      fprintf ( stderr, "  Zero pivot on step %d\n", jcol );
      exit ( 1 );
    }
/*
  Switch rows JCOL and IPIV, and X.
*/
    if ( jcol != ipiv )
    {
      for ( j = 1; j <= n; j++ )
      {
        t                 = a[jcol-1+(j-1)*n];
        a[jcol-1+(j-1)*n] = a[ipiv-1+(j-1)*n];
        a[ipiv-1+(j-1)*n] = t;
      }
      t         = x[jcol-1];
      x[jcol-1] = x[ipiv-1];
      x[ipiv-1] = t;
    }
/*
  Scale the pivot row.
*/
    t = a[jcol-1+(jcol-1)*n];
    a[jcol-1+(jcol-1)*n] = 1.0;
    for ( j = jcol+1; j <= n; j++ )
    {
      a[jcol-1+(j-1)*n] = a[jcol-1+(j-1)*n] / t;
    }
    x[jcol-1] = x[jcol-1] / t;
/*
  Use the pivot row to eliminate lower entries in that column.
*/
    for ( i = jcol+1; i <= n; i++ )
    {
      if ( a[i-1+(jcol-1)*n] != 0.0 )
      {
        t = - a[i-1+(jcol-1)*n];
        a[i-1+(jcol-1)*n] = 0.0;
        for ( j = jcol+1; j <= n; j++ )
        {
          a[i-1+(j-1)*n] = a[i-1+(j-1)*n] + t * a[jcol-1+(j-1)*n];
        }
        x[i-1] = x[i-1] + t * x[jcol-1];
      }
    }
  }
/*
  Back solve.
*/
  for ( jcol = n; 2 <= jcol; jcol-- )
  {
    for ( i = 1; i < jcol; i++ )
    {
      x[i-1] = x[i-1] - a[i-1+(jcol-1)*n] * x[jcol-1];
    }
  }

  return;
}
/******************************************************************************/

double *r8ge_fs_new ( int n, double a[], double b[] )

/******************************************************************************/
/*
  Purpose:

    R8GE_FS_NEW factors and solves a R8GE system.

  Discussion:

    The R8GE storage format is used for a "general" M by N matrix.  
    A physical storage space is made for each logical entry.  The two 
    dimensional logical array is mapped to a vector, in which storage is 
    by columns.

    The function does not save the LU factors of the matrix, and hence cannot
    be used to efficiently solve multiple linear systems, or even to
    factor A at one time, and solve a single linear system at a later time.

    The function uses partial pivoting, but no pivot vector is required.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    10 February 2012

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of the matrix.
    N must be positive.

    Input/output, double A[N*N].
    On input, A is the coefficient matrix of the linear system.
    On output, A is in unit upper triangular form, and
    represents the U factor of an LU factorization of the
    original coefficient matrix.

    Input, double B[N], the right hand side of the linear system.

    Output, double R8GE_FS_NEW[N], the solution of the linear system.
*/
{
  int i;
  int ipiv;
  int j;
  int jcol;
  double piv;
  double t;
  double *x;

  x = ( double * ) malloc ( n * sizeof ( double ) );

  for ( i = 0; i < n; i++ )
  {
    x[i] = b[i];
  }

  for ( jcol = 1; jcol <= n; jcol++ )
  {
/*
  Find the maximum element in column I.
*/
    piv = r8_abs ( a[jcol-1+(jcol-1)*n] );
    ipiv = jcol;
    for ( i = jcol+1; i <= n; i++ )
    {
      if ( piv < r8_abs ( a[i-1+(jcol-1)*n] ) )
      {
        piv = r8_abs ( a[i-1+(jcol-1)*n] );
        ipiv = i;
      }
    }

    if ( piv == 0.0 )
    {
      fprintf ( stderr, "\n" );
      fprintf ( stderr, "R8GE_FS_NEW - Fatal error!\n" );
      fprintf ( stderr, "  Zero pivot on step %d\n", jcol );
      exit ( 1 );
    }
/*
  Switch rows JCOL and IPIV, and X.
*/
    if ( jcol != ipiv )
    {
      for ( j = 1; j <= n; j++ )
      {
        t                 = a[jcol-1+(j-1)*n];
        a[jcol-1+(j-1)*n] = a[ipiv-1+(j-1)*n];
        a[ipiv-1+(j-1)*n] = t;
      }
      t         = x[jcol-1];
      x[jcol-1] = x[ipiv-1];
      x[ipiv-1] = t;
    }
/*
  Scale the pivot row.
*/
    t = a[jcol-1+(jcol-1)*n];
    a[jcol-1+(jcol-1)*n] = 1.0;
    for ( j = jcol+1; j <= n; j++ )
    {
      a[jcol-1+(j-1)*n] = a[jcol-1+(j-1)*n] / t;
    }
    x[jcol-1] = x[jcol-1] / t;
/*
  Use the pivot row to eliminate lower entries in that column.
*/
    for ( i = jcol+1; i <= n; i++ )
    {
      if ( a[i-1+(jcol-1)*n] != 0.0 )
      {
        t = - a[i-1+(jcol-1)*n];
        a[i-1+(jcol-1)*n] = 0.0;
        for ( j = jcol+1; j <= n; j++ )
        {
          a[i-1+(j-1)*n] = a[i-1+(j-1)*n] + t * a[jcol-1+(j-1)*n];
        }
        x[i-1] = x[i-1] + t * x[jcol-1];
      }
    }
  }
/*
  Back solve.
*/
  for ( jcol = n; 2 <= jcol; jcol-- )
  {
    for ( i = 1; i < jcol; i++ )
    {
      x[i-1] = x[i-1] - a[i-1+(jcol-1)*n] * x[jcol-1];
    }
  }

  return x;
}
/******************************************************************************/

void r8ge_fss ( int n, double a[], int nb, double x[] )

/******************************************************************************/
/*
  Purpose:

    R8GE_FSS factors and solves a system with multiple right hand sides.

  Discussion:

    This routine uses partial pivoting, but no pivot vector is required.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    28 November 2011

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of the matrix.
    N must be positive.

    Input/output, double A[N*N].
    On input, A is the coefficient matrix of the linear system.
    On output, A is in unit upper triangular form, and
    represents the U factor of an LU factorization of the
    original coefficient matrix.
    Input, int NB, the number of right hand sides.

    Input/output, double X[N*NB], on input, the right hand sides of the
    linear systems.  On output, the solutions of the linear systems.
*/
{
  int i;
  int ipiv;
  int j;
  int jcol;
  double piv;
  double t;

  for ( jcol = 1; jcol <= n; jcol++ )
  {
/*
  Find the maximum element in column I.
*/
    piv = r8_abs ( a[jcol-1+(jcol-1)*n] );
    ipiv = jcol;
    for ( i = jcol+1; i <= n; i++ )
    {
      if ( piv < r8_abs ( a[i-1+(jcol-1)*n] ) )
      {
        piv = r8_abs ( a[i-1+(jcol-1)*n] );
        ipiv = i;
      }
    }

    if ( piv == 0.0 )
    {
      fprintf ( stderr, "\n" );
      fprintf ( stderr, "R8GE_FSS - Fatal error!\n" );
      fprintf ( stderr, "  Zero pivot on step %d\n", jcol );
      exit ( 1 );
    }
/*
  Switch rows JCOL and IPIV, and X.
*/
    if ( jcol != ipiv )
    {
      for ( j = 1; j <= n; j++ )
      {
        t                 = a[jcol-1+(j-1)*n];
        a[jcol-1+(j-1)*n] = a[ipiv-1+(j-1)*n];
        a[ipiv-1+(j-1)*n] = t;
      }
      for ( j = 0; j < nb; j++ )
      {
        t            = x[jcol-1+j*n];
        x[jcol-1+j*n] = x[ipiv-1+j*n];
        x[ipiv-1+j*n] = t;
      }
    }
/*
  Scale the pivot row.
*/
    t = a[jcol-1+(jcol-1)*n];
    a[jcol-1+(jcol-1)*n] = 1.0;
    for ( j = jcol+1; j <= n; j++ )
    {
      a[jcol-1+(j-1)*n] = a[jcol-1+(j-1)*n] / t;
    }
    for ( j = 0; j < nb; j++ )
    {
      x[jcol-1+j*n] = x[jcol-1+j*n] / t;
    }
/*
  Use the pivot row to eliminate lower entries in that column.
*/
    for ( i = jcol+1; i <= n; i++ )
    {
      if ( a[i-1+(jcol-1)*n] != 0.0 )
      {
        t = - a[i-1+(jcol-1)*n];
        a[i-1+(jcol-1)*n] = 0.0;
        for ( j = jcol+1; j <= n; j++ )
        {
          a[i-1+(j-1)*n] = a[i-1+(j-1)*n] + t * a[jcol-1+(j-1)*n];
        }
        for ( j = 0; j < nb; j++ )
        {
          x[i-1+j*n] = x[i-1+j*n] + t * x[jcol-1+j*n];
        }
      }
    }
  }
/*
  Back solve.
*/
  for ( jcol = n; 2 <= jcol; jcol-- )
  {
    for ( i = 1; i < jcol; i++ )
    {
      for ( j = 0; j < nb; j++ )
      {
        x[i-1+j*n] = x[i-1+j*n] - a[i-1+(jcol-1)*n] * x[jcol-1+j*n];
      }
    }
  }

  return;
}
/******************************************************************************/

double *r8ge_fss_new ( int n, double a[], int nb, double b[] )

/******************************************************************************/
/*
  Purpose:

    R8GE_FSS_NEW factors and solves a system with multiple right hand sides.

  Discussion:

    This routine uses partial pivoting, but no pivot vector is required.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    28 November 2011

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of the matrix.
    N must be positive.

    Input/output, double A[N*N].
    On input, A is the coefficient matrix of the linear system.
    On output, A is in unit upper triangular form, and
    represents the U factor of an LU factorization of the
    original coefficient matrix.

    Input, int NB, the number of right hand sides.

    Input, double B[N*NB], the right hand sides of the linear systems.

    Output, double R8GE_FSS_NEW[N*NB], the solutions of the linear systems.
*/
{
  int i;
  int ipiv;
  int j;
  int jcol;
  double piv;
  double t;
  double *x;

  x = ( double * ) malloc ( n * nb * sizeof ( double ) );

  for ( j = 0; j < nb; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      x[i+j*n] = b[i+j*n];
    }
  }
  for ( jcol = 1; jcol <= n; jcol++ )
  {
/*
  Find the maximum element in column I.
*/
    piv = r8_abs ( a[jcol-1+(jcol-1)*n] );
    ipiv = jcol;
    for ( i = jcol+1; i <= n; i++ )
    {
      if ( piv < r8_abs ( a[i-1+(jcol-1)*n] ) )
      {
        piv = r8_abs ( a[i-1+(jcol-1)*n] );
        ipiv = i;
      }
    }

    if ( piv == 0.0 )
    {
      fprintf ( stderr, "\n" );
      fprintf ( stderr, "R8GE_FSS_NEW - Fatal error!\n" );
      fprintf ( stderr, "  Zero pivot on step %d\n", jcol );
      exit ( 1 );
    }
/*
  Switch rows JCOL and IPIV, and X.
*/
    if ( jcol != ipiv )
    {
      for ( j = 1; j <= n; j++ )
      {
        t                 = a[jcol-1+(j-1)*n];
        a[jcol-1+(j-1)*n] = a[ipiv-1+(j-1)*n];
        a[ipiv-1+(j-1)*n] = t;
      }
      for ( j = 0; j < nb; j++ )
      {
        t            = x[jcol-1+j*n];
        x[jcol-1+j*n] = x[ipiv-1+j*n];
        x[ipiv-1+j*n] = t;
      }
    }
/*
  Scale the pivot row.
*/
    t = a[jcol-1+(jcol-1)*n];
    a[jcol-1+(jcol-1)*n] = 1.0;
    for ( j = jcol+1; j <= n; j++ )
    {
      a[jcol-1+(j-1)*n] = a[jcol-1+(j-1)*n] / t;
    }
    for ( j = 0; j < nb; j++ )
    {
      x[jcol-1+j*n] = x[jcol-1+j*n] / t;
    }
/*
  Use the pivot row to eliminate lower entries in that column.
*/
    for ( i = jcol+1; i <= n; i++ )
    {
      if ( a[i-1+(jcol-1)*n] != 0.0 )
      {
        t = - a[i-1+(jcol-1)*n];
        a[i-1+(jcol-1)*n] = 0.0;
        for ( j = jcol+1; j <= n; j++ )
        {
          a[i-1+(j-1)*n] = a[i-1+(j-1)*n] + t * a[jcol-1+(j-1)*n];
        }
        for ( j = 0; j < nb; j++ )
        {
          x[i-1+j*n] = x[i-1+j*n] + t * x[jcol-1+j*n];
        }
      }
    }
  }
/*
  Back solve.
*/
  for ( jcol = n; 2 <= jcol; jcol-- )
  {
    for ( i = 1; i < jcol; i++ )
    {
      for ( j = 0; j < nb; j++ )
      {
        x[i-1+j*n] = x[i-1+j*n] - a[i-1+(jcol-1)*n] * x[jcol-1+j*n];
      }
    }
  }

  return x;
}
/******************************************************************************/

double *r8ge_identity ( int n )

/******************************************************************************/
/*
  Purpose:

    R8GE_IDENTITY sets a R8GE matrix to the identity.

  Discussion:

    The R8GE storage format is used for a "general" M by N matrix.  
    A physical storage space is made for each logical entry.  The two 
    dimensional logical array is mapped to a vector, in which storage is 
    by columns.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    25 February 2012

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of A.

    Output, double R8GE_IDENTITY[N*N], the N by N identity matrix.
*/
{
  double *a;
  int i;
  int j;

  a = ( double * ) malloc ( n * n * sizeof ( double ) );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      if ( i == j )
      {
        a[i+j*n] = 1.0;
      }
      else
      {
        a[i+j*n] = 0.0;
      }
    }
  }

  return a;
}
/******************************************************************************/

void r8ge_ilu ( int m, int n, double a[], double l[], double u[] )

/******************************************************************************/
/*
  Purpose:

    R8GE_ILU produces the incomplete LU factors of a R8GE matrix.

  Discussion:

    The R8GE storage format is used for a "general" M by N matrix.  
    A physical storage space is made for each logical entry.  The two 
    dimensional logical array is mapped to a vector, in which storage is 
    by columns.

    The incomplete LU factors of the M by N matrix A are:

      L, an M by M unit lower triangular matrix,
      U, an M by N upper triangular matrix

    with the property that L and U are computed in the same way as
    the usual LU factors, except that, whenever an off diagonal element
    of the original matrix is zero, then the corresponding value of
    U is forced to be zero.

    This condition means that it is no longer the case that A = L*U.

    On the other hand, L and U will have a simple sparsity structure
    related to that of A.  The incomplete LU factorization is generally
    used as a preconditioner in iterative schemes applied to sparse
    matrices.  It is presented here merely for illustration.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    25 February 2012

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows in A.

    Input, int N, the number of columns in A.

    Input, double A[M*N], the R8GE matrix.

    Output, double L[M*M], the M by M unit lower triangular factor.

    Output, double U[M*N], the M by N upper triangular factor.
*/
{
  int i;
  int j;
  int jhi;
  int k;
/*
  Initialize:

    L := M by M Identity
    U := A
*/
  for ( j = 0; j < m; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      if ( i == j )
      {
        l[i+j*m] = 1.0;
      }
      else
      {
        l[i+j*m] = 0.0;
      }
    }
  }

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      u[i+j*m] = a[i+j*m];
    }
  }

  jhi = i4_min ( m - 1, n );

  for ( j = 0; j < jhi; j++ )
  {
/*
  Zero out the entries in column J, from row J+1 to M.
*/
    for ( i = j+1; i < m; i++ )
    {
      if ( u[i+j*m] != 0.0 )
      {
        l[i+j*m] = u[i+j*m] / u[j+j*m];
        u[i+j*m] = 0.0;

        for ( k = j + 1; k < n; k++ )
        {
          if ( u[i+k*m] != 0.0 )
          {
            u[i+k*m] = u[i+k*m] - l[i+j*m] * u[j+k*m];
          }
        }
      }
    }
  }

  return;
}
/******************************************************************************/

double *r8ge_indicator ( int m, int n )

/******************************************************************************/
/*
  Purpose:

    R8GE_INDICATOR sets up a R8GE indicator matrix.

  Discussion:

    The R8GE storage format is used for a "general" M by N matrix.  
    A physical storage space is made for each logical entry.  The two 
    dimensional logical array is mapped to a vector, in which storage is 
    by columns.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    25 February 2012

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows of the matrix.
    M must be positive.

    Input, int N, the number of columns of the matrix.
    N must be positive.

    Output, double R8GE_INDICATOR[M*N], the R8GE matrix.
*/
{
  double *a;
  int fac;
  int i;
  int j;

  a = ( double * ) malloc ( m * n * sizeof ( double ) );

  fac = i4_power ( 10, i4_log_10 ( n ) + 1 );

  for ( i = 1; i <= m; i++ )
  {
    for ( j = 1; j <= n; j++ )
    {
      a[i-1+(j-1)*m] = ( double ) ( fac * i + j );
    }
  }

  return a;
}
/******************************************************************************/

double *r8ge_inverse ( int n, double a[], int pivot[] )

/******************************************************************************/
/*
  Purpose:

    R8GE_INVERSE computes the inverse of a R8GE matrix factored by R8GE_FA.

  Discussion:

    The R8GE storage format is used for a "general" M by N matrix.  
    A physical storage space is made for each logical entry.  The two 
    dimensional logical array is mapped to a vector, in which storage is 
    by columns.

    R8GE_INVERSE is a simplified standalone version of the LINPACK routine
    SGEDI.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    25 February 2012

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of the matrix A.

    Input, double A[N*N], the factor information computed by R8GE_FA.

    Input, int PIVOT(N), the pivot vector from R8GE_FA.

    Output, double R8GE_INVERSE[N*N], the inverse matrix.
*/
{
  double *b;
  int i;
  int j;
  int k;
  double temp;

  b = ( double * ) malloc ( n * n * sizeof ( double ) );
/*
  Compute Inverse(U).
*/
  for ( k = 1; k <= n; k++ )
  {
    for ( i = 1; i <= k-1; i++ )
    {
      b[i-1+(k-1)*n] = -b[i-1+(k-1)*n] / a[k-1+(k-1)*n];
    }
    b[k-1+(k-1)*n] = 1.0 / a[k-1+(k-1)*n];

    for ( j = k+1; j <= n; j++ )
    {
      b[k-1+(j-1)*n] = 0.0;
      for ( i = 1; i <= k; i++ )
      {
        b[i-1+(j-1)*n] = b[i-1+(j-1)*n] + b[i-1+(k-1)*n] * a[k-1+(j-1)*n];
      }
    }
  }
/*
  Multiply Inverse(U) by Inverse(L).
*/
  for ( k = n-1; 1 <= k; k-- )
  {
    for ( i = k+1; i <= n; i++ )
    {
      b[i-1+(k-1)*n] = 0.0;
    }

    for ( j = k+1; j <= n; j++ )
    {
      for ( i = 1; i <= n; i++ )
      {
        b[i-1+(k-1)*n] = b[i-1+(k-1)*n] + b[i-1+(j-1)*n] * a[j-1+(k-1)*n];
      }
    }

    if ( pivot[k-1] != k )
    {
      for ( i = 1; i <= n; i++ )
      {
        temp = b[i-1+(k-1)*n];
        b[i-1+(k-1)*n] = b[i-1+(pivot[k-1]-1)*n];
        b[i-1+(pivot[k-1]-1)*n] = temp;
      }

    }

  }

  return b;
}
/******************************************************************************/

double *r8ge_ml ( int n, double a_lu[], int pivot[], double x[], int job )

/******************************************************************************/
/*
  Purpose:

    R8GE_ML computes A * x or A' * x, using R8GE_FA factors.

  Discussion:

    The R8GE storage format is used for a "general" M by N matrix.  
    A physical storage space is made for each logical entry.  The two 
    dimensional logical array is mapped to a vector, in which storage is 
    by columns.

    It is assumed that R8GE_FA has overwritten the original matrix
    information by LU factors.  R8GE_ML is able to reconstruct the
    original matrix from the LU factor data.

    R8GE_ML allows the user to check that the solution of a linear
    system is correct, without having to save an unfactored copy
    of the matrix.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    28 February 2012

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of the matrix.
    N must be positive.

    Input, double A_LU[N*N], the LU factors from R8GE_FA.

    Input, int PIVOT[N], the pivot vector computed by R8GE_FA.

    Input, double X[N], the vector to be multiplied.

    Input, int JOB, specifies the operation to be done:
    JOB = 0, compute A * x.
    JOB nonzero, compute A' * X.

    Output, double R8GE_ML[N], the result of the multiplication.
*/
{
  double *b;
  int i;
  int j;
  int k;
  double temp;

  b = ( double * ) malloc ( n * sizeof ( double ) );

  for ( i = 0; i < n; i++ )
  {
    b[i] = x[i];
  }

  if ( job == 0 )
  {
/*
  Y = U * X.
*/
    for ( j = 1; j <= n; j++ )
    {
      for ( i = 1; i <= j-1; i++ )
      {
        b[i-1] = b[i-1] + a_lu[i-1+(j-1)*n] * b[j-1];
      }
      b[j-1] = a_lu[j-1+(j-1)*n] * b[j-1];
    }
/*
  B = PL * Y = PL * U * X = A * x.
*/
    for ( j = n-1; 1 <= j; j-- )
    {
      for ( i = j+1; i <= n; i++ )
      {
        b[i-1] = b[i-1] - a_lu[i-1+(j-1)*n] * b[j-1];
      }

      k = pivot[j-1];

      if ( k != j )
      {
        temp   = b[k-1];
        b[k-1] = b[j-1];
        b[j-1] = temp;
      }
    }
  }
  else
  {
/*
  Y = (PL)' * X:
*/
    for ( j = 1; j <= n-1; j++ )
    {
      k = pivot[j-1];

      if ( k != j )
      {
        temp   = b[k-1];
        b[k-1] = b[j-1];
        b[j-1] = temp;
      }

      temp = 0.0;
      for ( i = j+1; i <= n; i++ )
      {
        temp = temp + b[i-1] * a_lu[i-1+(j-1)*n];
      }
      b[j-1] = b[j-1] - temp;

    }
/*
  B = U' * Y = ( PL * U )' * X = A' * X.
*/
    for ( i = n; 1 <= i; i-- )
    {
      for ( j = i+1; j <= n; j++ )
      {
        b[j-1] = b[j-1] + b[i-1] * a_lu[i-1+(j-1)*n];
      }
      b[i-1] = b[i-1] * a_lu[i-1+(i-1)*n];
    }

  }

  return b;
}
/******************************************************************************/

double *r8ge_mu ( int m, int n, double a_lu[], char trans, int pivot[], 
  double x[] )

/******************************************************************************/
/*
  Purpose:

    R8GE_MU computes A * x or A' * x, using R8GE_TRF factors.

  Discussion:

    The R8GE storage format is used for a "general" M by N matrix.  
    A physical storage space is made for each logical entry.  The two 
    dimensional logical array is mapped to a vector, in which storage is 
    by columns.

    It is assumed that R8GE_TRF has overwritten the original matrix
    information by PLU factors.  R8GE_MU is able to reconstruct the
    original matrix from the PLU factor data.

    R8GE_MU allows the user to check that the solution of a linear
    system is correct, without having to save an unfactored copy
    of the matrix.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    28 February 2012

  Author:

    John Burkardt

  Reference:

    Edward Anderson, Zhaojun Bai, Christian Bischof, Susan Blackford, 
    James Demmel, Jack Dongarra, Jeremy Du Croz, Anne Greenbaum, 
    Sven Hammarling, Alan McKenney, Danny Sorensen,
    LAPACK User's Guide,
    Second Edition,
    SIAM, 1995.

  Parameters:

    Input, int M, the number of rows in the matrix.

    Input, int N, the number of columns in the matrix.

    Input, double A_LU[M*N], the LU factors from R8GE_TRF.

    Input, char TRANS, specifies the form of the system of equations:
    'N':  A * x = b  (No transpose)
    'T':  A'* X = B  (Transpose)
    'C':  A'* X = B  (Conjugate transpose = Transpose)

    Input, int PIVOT[*], the pivot vector computed by R8GE_TRF.

    Input, double X[*], the vector to be multiplied.
    For the untransposed case, X should have N entries.
    For the transposed case, X should have M entries.

    Output, double R8GE_MU[*], the result of the multiplication.
    For the untransposed case, the result should have M entries.
    For the transposed case, the result should have N entries.
*/
{
  double *b;
  int i;
  int j;
  int k;
  int mn_max;
  int npiv;
  double temp;
  double *y;

  npiv = i4_min ( m - 1, n );
  mn_max = i4_max ( m, n );
  y = ( double * ) malloc ( mn_max * sizeof ( double ) );

  if ( trans == 'n' || trans == 'N' )
  {
    b = ( double * ) malloc ( m * sizeof ( double ) );
/*
  Y[MN] = U[MNxN] * X[N].
*/
    for ( i = 0; i < n; i++ )
    {
      y[i] = 0.0;
    }

    for ( j = 1; j <= n; j++ )
    {
      for ( i = 1; i <= i4_min ( j, m ); i++ )
      {
        y[i-1] = y[i-1] + a_lu[i-1+(j-1)*m] * x[j-1];
      }
    }
/*
  Z[M] = L[MxMN] * Y[MN] = L[MxMN] * U[MNxN] * X[N].
*/
    for ( i = 0; i < m; i++ )
    {
      if ( i < n ) 
      {
        b[i] = y[i];
      }
      else
      {
        b[i] = 0.0;
      }
    }

    for ( j = i4_min ( m-1, n ); 1 <= j; j-- )
    {
      for ( i = j+1; i <= m; i++ )
      {
        b[i-1] = b[i-1] + a_lu[i-1+(j-1)*m] * y[j-1];
      }
    }
/*
  B = P * Z = P * L * Y = P * L * U * X = A * x.
*/
    for ( j = npiv; 1 <= j; j-- )
    {
      k = pivot[j-1];

      if ( k != j )
      {
        temp = b[k-1];
        b[k-1] = b[j-1];
        b[j-1] = temp;
      }
    }
  }
  else if ( trans == 't' || trans == 'T' ||
            trans == 'c' || trans == 'C' )
  {
    b = ( double * ) malloc ( n * sizeof ( double ) );
/*
  Y = P' * X:
*/
    for ( i = 1; i <= npiv; i++ )
    {
      k = pivot[i-1];

      if ( k != i )
      {
        temp = x[k-1];
        x[k-1] = x[i-1];
        x[i-1] = temp;
      }
    }

    for ( i = 0; i < n; i++ )
    {
      if ( i < m ) 
      {
        b[i] = x[i];
      }
      else
      {
        b[i] = 0.0;
      }
    }
/*
  Z = L' * Y:
*/
    for ( j = 1; j <= i4_min ( m - 1, n ); j++ )
    {
      for ( i = j+1; i <= m; i++ )
      {
        b[j-1] = b[j-1] + x[i-1] * a_lu[i-1+(j-1)*m];
      }
    }
/*
  B = U' * Z.
*/
    for ( i = m; 1 <= i; i-- )
    {
      for ( j = i+1; j <= n; j++ )
      {
        b[j-1] = b[j-1] + b[i-1] * a_lu[i-1+(j-1)*m];
      }
      if ( i <= n )
      {
        b[i-1] = b[i-1] * a_lu[i-1+(i-1)*m];
      }
    }
/*
  Now restore X.
*/
    for ( i = npiv; 1 <= i; i-- )
    {
      k = pivot[i-1];

      if ( k != i )
      {
        temp = x[k-1];
        x[k-1] = x[i-1];
        x[i-1] = temp;
      }
    }
  }
/*
  Illegal value of TRANS.
*/
  else
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8GE_MU - Fatal error!\n" );
    fprintf ( stderr, "  Illegal value of TRANS = \"%c\"\n", trans );
    exit ( 1 );
  }

  free ( y );

  return b;
}
/******************************************************************************/

double *r8ge_mxm ( int n, double a[], double b[] )

/******************************************************************************/
/*
  Purpose:

    R8GE_MXM multiplies two R8GE matrices.

  Discussion:

    The R8GE storage format is used for a "general" M by N matrix.  
    A physical storage space is made for each logical entry.  The two 
    dimensional logical array is mapped to a vector, in which storage is 
    by columns.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    28 February 2012

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of the matrices.
    N must be positive.

    Input, double A[N*N], B[N*N], the R8GE factor matrices.

    Output, double C[N*N], the R8GE product matrix.
*/
{
  double *c;
  int i;
  int j;
  int k;

  c = ( double * ) malloc ( n * n * sizeof ( double ) );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      c[i+j*n] = 0.0;
      for ( k = 0; k < n; k++ )
      {
        c[i+j*n] = c[i+j*n] + a[i+k*n] * b[k+j*n];
      }
    }
  }

  return c;
}
/******************************************************************************/

double *r8ge_mxv ( int m, int n, double a[], double x[] )

/******************************************************************************/
/*
  Purpose:

    R8GE_MXV multiplies a R8GE matrix times a vector.

  Discussion:

    The R8GE storage format is used for a "general" M by N matrix.  
    A physical storage space is made for each logical entry.  The two 
    dimensional logical array is mapped to a vector, in which storage is 
    by columns.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    28 February 2012

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows of the matrix.
    M must be positive.

    Input, int N, the number of columns of the matrix.
    N must be positive.

    Input, double A[M*N], the R8GE matrix.

    Input, double X[N], the vector to be multiplied by A.

    Output, double R8GE_MXV[M], the product A * x.
*/
{
  double *b;
  int i;
  int j;

  b = ( double * ) malloc ( m * sizeof ( double ) );

  for ( i = 0; i < m; i++ )
  {
    b[i] = 0.0;
    for ( j = 0; j < n; j++ )
    {
      b[i] = b[i] + a[i+j*m] * x[j];
    }
  }

  return b;
}
/******************************************************************************/

double r8ge_np_det ( int n, double a_lu[] )

/******************************************************************************/
/*
  Purpose:

    R8GE_NP_DET computes the determinant of a matrix factored by R8GE_NP_FA.

  Discussion:

    The R8GE storage format is used for a "general" M by N matrix.  
    A physical storage space is made for each logical entry.  The two 
    dimensional logical array is mapped to a vector, in which storage is 
    by columns.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    28 February 2012

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of the matrix.
    N must be positive.

    Input, double A_LU[N*N], the LU factors from R8GE_NP_FA.

    Output, double R8GE_NP_DET, the determinant of the matrix.
*/
{
  double det;
  int i;

  det = 1.0;
  for ( i = 0; i < n; i++ )
  {
    det = det * a_lu[i+i*n];
  }

  return det;
}
/******************************************************************************/

int r8ge_np_fa ( int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8GE_NP_FA factors a R8GE matrix by nonpivoting Gaussian elimination.

  Discussion:

    The R8GE storage format is used for a "general" M by N matrix.  
    A physical storage space is made for each logical entry.  The two 
    dimensional logical array is mapped to a vector, in which storage is 
    by columns.

    R8GE_NP_FA is a version of the LINPACK routine SGEFA, but uses no
    pivoting.  It will fail if the matrix is singular, or if any zero
    pivot is encountered.

    If R8GE_NP_FA successfully factors the matrix, R8GE_NP_SL may be called
    to solve linear systems involving the matrix.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    28 February 2012

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of the matrix.
    N must be positive.

    Input/output, double A[N*N].
    On input, A contains the matrix to be factored.
    On output, A contains information about the factorization,
    which must be passed unchanged to R8GE_NP_SL for solutions.

    Output, int R8GE_NP_FA, singularity flag.
    0, no singularity detected.
    nonzero, the factorization failed on the INFO-th step.
*/
{
  int i;
  int j;
  int k;

  for ( k = 1; k <= n-1; k++ )
  {
    if ( a[k-1+(k-1)*n] == 0.0 )
    {
      return k;
    }

    for ( i = k+1; i <= n; i++ )
    {
      a[i-1+(k-1)*n] = -a[i-1+(k-1)*n] / a[k-1+(k-1)*n];
    }

    for ( j = k+1; j <= n; j++ )
    {
      for ( i = k+1; i <= n; i++ )
      {
        a[i-1+(j-1)*n] = a[i-1+(j-1)*n] + a[i-1+(k-1)*n] * a[k-1+(j-1)*n];
      }
    }
  }

  if ( a[n-1+(n-1)*n] == 0.0 )
  {
    return n;
  }

  return 0;
}
/******************************************************************************/

double *r8ge_np_inverse ( int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8GE_NP_INVERSE computes the inverse of a matrix factored by R8GE_NP_FA.

  Discussion:

    The R8GE storage format is used for a "general" M by N matrix.  
    A physical storage space is made for each logical entry.  The two 
    dimensional logical array is mapped to a vector, in which storage is 
    by columns.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    28 February 2012

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of the matrix A.

    Input, double A[N*N], the factor information computed by R8GE_NP_FA.

    Output, double R8GE_NP_INVERSE[N*N], the inverse matrix.
*/
{
  double *b;
  int i;
  int j;
  int k;
  double temp;
  double *work;

  b = ( double * ) malloc ( n * n * sizeof ( double ) );
  work = ( double * ) malloc ( n * sizeof ( double ) );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      b[i+j*n] = a[i+j*n];
    }
  }
/*
  Compute Inverse(U).
*/
  for ( k = 1; k <= n; k++ )
  {
    b[k-1+(k-1)*n] = 1.0 / b[k-1+(k-1)*n];
    for ( i = 1; i <= k-1; i++ )
    {
      b[i-1+(k-1)*n] = -b[i-1+(k-1)*n] * b[k-1+(k-1)*n];
    }
    for ( j = k+1; j <= n; j++ )
    {
      temp = b[k-1+(j-1)*n];
      b[k-1+(j-1)*n] = 0.0;
      for ( i = 1; i <= k; i++ )
      {
        b[i-1+(j-1)*n] = b[i-1+(j-1)*n] + temp * b[i-1+(k-1)*n];
      }
    }
  }
/*
  Form Inverse(U) * Inverse(L).
*/
  for ( k = n-1; 1 <= k; k-- )
  {
    for ( i = k+1; i <= n; i++ )
    {
      work[i-1] = b[i-1+(k-1)*n];
      b[i-1+(k-1)*n] = 0.0;
    }

    for ( j = k+1; j <= n; j++ )
    {
      for ( i = 1; i <= n; i++ )
      {
        b[i-1+(k-1)*n] = b[i-1+(k-1)*n] + b[i-1+(j-1)*n] * work[j-1];
      }
    }
  }

  free ( work );

  return b;
}
/******************************************************************************/

double *r8ge_np_ml ( int n, double a_lu[], double x[], int job )

/******************************************************************************/
/*
  Purpose:

    R8GE_NP_ML computes A * x or x * A, for a matrix factored by R8GE_NP_FA.

  Discussion:

    The R8GE storage format is used for a "general" M by N matrix.  
    A physical storage space is made for each logical entry.  The two 
    dimensional logical array is mapped to a vector, in which storage is 
    by columns.

    The matrix A is assumed to have been factored by R8GE_NP_FA.

    R8GE_NP_ML allows the user to check that the solution of a linear
    system is correct, without having to save an unfactored copy
    of the matrix.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    28 February 2012

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of the matrix.
    N must be positive.

    Input, double A_LU[N*N], the LU factors from R8GE_NP_FA.

    Input, double X[N], the vector to be multiplied.

    Input, int JOB, determines the multiplication to
    be carried out:
    JOB = 0, compute A * x.
    JOB nonzero, compute A' * X.

    Output, double R8GE_NP_ML[N], the result of the multiplication.
*/
{
  double *b;
  int i;
  int j;
  double t;

  b = ( double * ) malloc ( n * sizeof ( double ) );

  for ( i = 0; i < n; i++ )
  {
    b[i] = x[i];
  }

  if ( job == 0 )
  {
/*
  Compute U * X = Y:
*/
    for ( i = 0; i < n; i++ )
    {
      t = 0.0;
      for ( j = i; j < n; j++ )
      {
        t = t + a_lu[i+j*n] * b[j];
      }
      b[i] = t;
    }
/*
  Compute L * Y = B:
*/
    for ( j = n-2; 0 <= j; j-- )
    {
      for ( i = j+1; i < n; i++ )
      {
        b[i] = b[i] - a_lu[i+j*n] * b[j];
      }
    }
  }
  else
  {
/*
  Compute L' * X = Y:
*/
    for ( j = 0; j < n-1; j++ )
    {
      for ( i = j+1; i < n; i++ )
      {
        b[j] = b[j] - b[i] * a_lu[i+j*n];
      }
    }
/*
  Compute U' * Y = B:
*/
    for ( j = n-1; 0 <= j; j-- )
    {
      b[j] = b[j] * a_lu[j+j*n];
      for ( i = 0; i < j; i++ )
      {
        b[j] = b[j] + b[i] * a_lu[i+j*n];
      }
    }
  }

  return b;
}
/******************************************************************************/

double *r8ge_np_sl ( int n, double a_lu[], double b[], int job )

/******************************************************************************/
/*
  Purpose:

    R8GE_NP_SL solves a system factored by R8GE_NP_FA.

  Discussion:

    The R8GE storage format is used for a "general" M by N matrix.  
    A physical storage space is made for each logical entry.  The two 
    dimensional logical array is mapped to a vector, in which storage is 
    by columns.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    28 February 2012

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of the matrix.
    N must be positive.

    Input, double A_LU[N*N], the LU factors from R8GE_NP_FA.

    Input, double B[N], the right hand side.

    Input, int JOB.
    If JOB is zero, the routine will solve A * x = b.
    If JOB is nonzero, the routine will solve A' * x = b.

    Output, double R8GE_NP_SL[N], the solution.
*/
{
  int i;
  int k;
  double *x;
/*
  Solve A * x = b.
*/
  x = ( double * ) malloc ( n * sizeof ( double ) );
  for ( i = 0; i < n; i++ )
  {
    x[i] = b[i];
  }

  if ( job == 0 )
  {
    for ( k = 0; k < n-1; k++ )
    {
      for ( i = k+1; i < n; i++ )
      {
        x[i] = x[i] + a_lu[i+k*n] * x[k];
      }
    }

    for ( k = n-1; 0 <= k; k-- )
    {
      x[k] = x[k] / a_lu[k+k*n];
      for ( i = 0; i <= k-1; i++ )
      {
        x[i] = x[i] - a_lu[i+k*n] * x[k];
      }
    }
  }
/*
  Solve A' * X = B.
*/
  else
  {
    for ( k = 0; k < n; k++ )
    {
      for ( i = 0; i <= k-1; i++ )
      {
        x[k] = x[k] - x[i] * a_lu[i+k*n];
      }
      x[k] = x[k] / a_lu[k+k*n];
    }

    for ( k = n-2; 0 <= k; k-- )
    {
      for ( i = k+1; i < n; i++ )
      {
        x[k] = x[k] + x[i] * a_lu[i+k*n];
      }
    }

  }

  return x;
}
/******************************************************************************/

int r8ge_np_trf ( int m, int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8GE_NP_TRF computes the LU factorization of a R8GE matrix.

  Discussion:

    The R8GE storage format is used for a "general" M by N matrix.  
    A physical storage space is made for each logical entry.  The two 
    dimensional logical array is mapped to a vector, in which storage is 
    by columns.

    R8GE_NP_TRF is a nonpivoting version of R8GE_TRF, and will fail if
    a zero element is encountered along the diagonal.

    The factorization has the form
      A = L * U
    where L is lower triangular with unit diagonal elements (lower
    trapezoidal if N < M), and U is upper triangular (upper trapezoidal
    if M < N).

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    28 February 2012

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows of the matrix A.  0 <= M.

    Input, int N, the number of columns of the matrix A.  0 <= N.

    Input/output, double A[M*N].
    On entry, the M by N matrix to be factored.
    On exit, the factors L and U from the factorization
    A = L*U; the unit diagonal elements of L are not stored.

    Output, int R8GE_NP_TRF.
    = 0: successful exit
    = -K, the K-th argument had an illegal value
    = K, U(K,K) is exactly zero. The factorization
         has been completed, but the factor U is exactly
         singular, and division by zero will occur if it is used
         to solve a system of equations.
*/
{
  int i;
  int ii;
  int info;
  int j;
/*
  Test the input parameters.
*/
  info = 0;

  if ( m < 0 )
  {
    return (-1);
  }
  else if ( n < 0 )
  {
    return (-2);
  }

  if ( m == 0 || n == 0 )
  {
    return 0;
  }

  for ( j = 1; j <= i4_min ( m, n ); j++ )
  {
/*
  Compute elements J+1:M of the J-th column.
*/
    if ( a[j-1+(j-1)*m] != 0.0 )
    {
      for ( i = j+1; i <= m; i++ )
      {
        a[i-1+(j-1)*m] = a[i-1+(j-1)*m] / a[j-1+(j-1)*m];
      }
    }
    else if ( info == 0 )
    {
      info = j;
    }
/*
  Update the trailing submatrix.
*/
    if ( j < i4_min ( m, n ) )
    {
      for ( ii = j+1; ii <= m; ii++ )
      {
        for ( i = j+1; i <= n; i++ )
        {
          a[ii-1+(i-1)*m] = a[ii-1+(i-1)*m] - a[ii-1+(j-1)*m] * a[j-1+(i-1)*m];
        }
      }
    }
  }

  return info;
}
/******************************************************************************/

double *r8ge_np_trm ( int m, int n, double a[], double x[], int job )

/******************************************************************************/
/*
  Purpose:

    R8GE_NP_TRM computes A * x or A' * x, for a matrix factored by R8GE_NP_TRF.

  Discussion:

    The R8GE storage format is used for a "general" M by N matrix.  
    A physical storage space is made for each logical entry.  The two 
    dimensional logical array is mapped to a vector, in which storage is 
    by columns.

    The matrix A is assumed to have been factored by R8GE_NP_TRF.

    R8GE_NP_TRM allows the user to check that the solution of a linear
    system is correct, without having to save an unfactored copy
    of the matrix.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    28 February 2012

  Author:

    John Burkardt

  Reference:

    Edward Anderson, Zhaojun Bai, Christian Bischof, Susan Blackford, 
    James Demmel, Jack Dongarra, Jeremy Du Croz, Anne Greenbaum, 
    Sven Hammarling, Alan McKenney, Danny Sorensen,
    LAPACK User's Guide,
    Second Edition,
    SIAM, 1995.

  Parameters:

    Input, int M, N, the number of rows and columns in the matrix.
    M and N must be positive.

    Input, double A[M*N], the M by N matrix factors computed by R8GE_NP_TRF.

    Input, double X[*], the vector to be multiplied.
    If JOB is 0, X must have dimension N.
    If JOB is nonzero, X must have dimension M.

    Input, int JOB, determines the multiplication to
    be carried out:
    JOB = 0, compute A * x.
    JOB nonzero, compute A' * X.

    Output, double R8GE_NP_TRM[*], the result of the multiplication.
    If JOB is 0, the output has dimension M.
    If JOB is nonzero, the output has dimension N.
*/
{
  double *b;
  int i;
  int j;
  double temp;

  if ( job == 0 )
  {
    b = ( double * ) malloc ( m * sizeof ( double ) );
    for ( i = 0; i < m; i++ )
    {
      b[i] = 0.0;
    }
/*
  Compute U * X = Y:
*/
    for ( i = 0; i < i4_min ( m, n ); i++ )
    {
      for ( j = i; j < n; j++ )
      {
        b[i] = b[i] + a[i+j*m] * x[j];
      }
    }
/*
  Compute L * Y = B:
*/
    for ( i = i4_min ( m-1, n ); 1 <= i; i-- )
    {
      for ( j = 0; j < i; j++ )
      {
        b[i] = b[i] + a[i+j*m] * b[j];
      }
    }
  }
  else
  {
    b = ( double * ) malloc ( n * sizeof ( double ) );
    for ( i = 0; i < n; i++ )
    {
      b[i] = 0.0;
    }
/*
  Compute L' * X = Y:
*/
    for ( i = 0; i < i4_min ( m, n ); i++ )
    {
      b[i] = x[i];
      for ( j = i+1; j < m; j++ )
      {
        b[i] = b[i] + a[j+i*m] * x[j];
      }
    }
/*
  Compute U' * Y = B:
*/
    for ( i = i4_min ( m, n ) - 1; 0 <= i; i-- )
    {
      temp = 0.0;
      for ( j = 0; j <= i; j++ )
      {
        temp = temp + a[j+i*m] * b[j];
      }
      b[i] = temp;
    }

  }

  return b;
}
/******************************************************************************/

double *r8ge_np_trs ( int n, int nrhs, char trans, double a[], double b[] )

/******************************************************************************/
/*
  Purpose:

    R8GE_NP_TRS solves a system of linear equations factored by R8GE_NP_TRF.

  Discussion:

    The R8GE storage format is used for a "general" M by N matrix.  
    A physical storage space is made for each logical entry.  The two 
    dimensional logical array is mapped to a vector, in which storage is 
    by columns.

    R8GE_NP_TRS is a nonpivoting version of R8GE_TRS.

    R8GE_TRS solves a system of linear equations
      A * x = b  or  A' * X = B
    with a general N by N matrix A using the LU factorization computed
    by R8GE_NP_TRF.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    28 February 2012

  Author:

    John Burkardt

  Reference:

    Edward Anderson, Zhaojun Bai, Christian Bischof, Susan Blackford, 
    James Demmel, Jack Dongarra, Jeremy Du Croz, Anne Greenbaum, 
    Sven Hammarling, Alan McKenney, Danny Sorensen,
    LAPACK User's Guide,
    Second Edition,
    SIAM, 1995.

  Parameters:

    Input, int N, the order of the matrix A.  0 <= N.

    Input, int NRHS, the number of right hand sides.  0 <= NRHS.

    Input, char TRANS, pecifies the form of the system of equations:
    'N':  A * x = b  (No transpose)
    'T':  A'* X = B  (Transpose)
    'C':  A'* X = B  (Conjugate transpose = Transpose)

    Input, double A[N*N], the factors L and U from the factorization
    A = L*U as computed by R8GE_NP_TRF.

    Input, double B[N*NRHS], the right hand side matrix B.

    Output, double R8GE_NP_TRS[N*NRHS], the solution matrix X.
*/
{
  int i;
  int j;
  int k;
  double *x;

  if ( trans != 'n' && trans != 'N' && 
       trans != 't' && trans != 'T' && 
       trans != 'c' && trans != 'C' )
  {
    return NULL;
  }
  if ( n < 0 )
  {
    return NULL;
  }
  if ( nrhs < 0 )
  {
    return NULL;
  }

  if ( n == 0 || nrhs == 0 )
  {
    return NULL;
  }

  x = ( double * ) malloc ( n * nrhs * sizeof ( double ) );

  for ( j = 0; j < nrhs; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      x[i+j*n] = b[i+j*n];
    }
  }

  if ( trans == 'n' || trans == 'N' )
  {
/*
  Solve L * x = b, overwriting b with x.
*/
    for ( k = 0; k < nrhs; k++ )
    {
      for ( j = 1; j <= n-1; j++ )
      {
        for ( i = j+1; i <= n; i++ )
        {
          x[i-1+k*n] = x[i-1+k*n] - a[i-1+(j-1)*n] * x[j-1+k*n];
        }
      }
    }
/*
  Solve U * x = b, overwriting b with x.
*/
    for ( k = 0; k < nrhs; k++ )
    {
      for ( j = n; 1 <= j; j-- )
      {
        x[j-1+k*n] = x[j-1+k*n] / a[j-1+(j-1)*n];
        for ( i = 1; i <= j-1; i++ )
        {
          x[i-1+k*n] = x[i-1+k*n] - a[i-1+(j-1)*n] * x[j-1+k*n];
        }
      }
    }
  }
  else
/*
  Solve U' * x = b, overwriting b with x.
*/
  {
    for ( k = 0; k < nrhs; k++ )
    {
      for ( j = 1; j <= n; j++ )
      {
        x[j-1+k*n] = x[j-1+k*n] / a[j-1+(j-1)*n];
        for ( i = j+1; i <= n; i++ )
        {
          x[i-1+k*n] = x[i-1+k*n] - a[j-1+(i-1)*n] * x[j-1+k*n];
        }
      }
    }
/*
  Solve L' * x = b, overwriting b with x.
*/
    for ( k = 0; k < nrhs; k++ )
    {
      for ( j = n; 2 <= j; j-- )
      {
        for ( i = 1; i <= j-1; i++ )
        {
          x[i-1+k*n] = x[i-1+k*n] - a[j-1+(i-1)*n] * x[j-1+k*n];
        }
      }
    }
  }

  return x;
}
/******************************************************************************/

void r8ge_plu ( int m, int n, double a[], double p[], double l[], double u[] )

/******************************************************************************/
/*
  Purpose:

    R8GE_PLU produces the PLU factors of a R8GE matrix.

  Discussion:

    The R8GE storage format is used for a "general" M by N matrix.  
    A physical storage space is made for each logical entry.  The two 
    dimensional logical array is mapped to a vector, in which storage is 
    by columns.

    The PLU factors of the M by N matrix A are:

      P, an M by M permutation matrix P,
      L, an M by M unit lower triangular matrix,
      U, an M by N upper triangular matrix.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    28 February 2012

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows in A.

    Input, int N, the number of columns in A.

    Input, double A[M,N], the M by N matrix to be factored.

    Output, double P[M*M], the M by M permutation factor.

    Output, double L[M*M], the M by M unit lower triangular factor.

    Output, double U[M*N], the M by N upper triangular factor.
*/
{
  int i;
  int j;
  int k;
  int pivot_row;
  double pivot_value;
  double t;
/*
  Initialize:

    P: = M by M Identity
    L: = M by M Identity
    U: = A
*/

  for ( i = 0; i < m; i++ )
  {
    for ( j = 0; j < m; j++ )
    {
      if ( i == j ) 
      {
        l[i+j*m] = 1.0;
        p[i+j*m] = 1.0;
      }
      else
      {
        l[i+j*m] = 0.0;
        p[i+j*m] = 0.0;
      }
    }
  }

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      u[i+j*m] = a[i+j*m];
    }
  }
/*
  On step J, find the pivot row and the pivot value.
*/
  for ( j = 0; j < i4_min ( m-2, n-1 ); j++ )
  {
    pivot_value = 0.0;
    pivot_row = -1;

    for ( i = j; i < m; i++ )
    {
      if ( pivot_value < r8_abs ( u[i+j*m] ) )
      {
        pivot_value = r8_abs ( u[i+j*m] );
        pivot_row = i;
      }
    }
/*
  If the pivot row is nonzero swap:
  * rows J and PIVOT_ROW of U;
  * rows J and PIVOT_ROW of L and cols J and PIVOT_ROW of L;
  * cols J and PIVOT_ROW of P.
*/
    if ( pivot_row != -1 )
    {
      for ( k = 0; k < n; k++ )
      {
        t                = u[j+k*m];
        u[j+k*m]         = u[pivot_row+k*m];
        u[pivot_row+k*m] = t;
      }

      for ( k = 0; k < m; k++ )
      {
        t                = l[j+k*m];
        l[j+k*m]         = l[pivot_row+k*m];
        l[pivot_row+k*m] = t;
      }
      for ( k = 0; k < m; k++ )
      {
        t                = l[k+j*m];
        l[k+j*m]         = l[k+pivot_row*m];
        l[k+pivot_row*m] = t;
      }

      for ( k = 0; k < m; k++ )
      {
        t                = p[k+j*m];
        p[k+j*m]         = p[k+pivot_row*m];
        p[k+pivot_row*m] = t;
      }
/*
  Zero out the entries in column J, from row J+1 to M.
*/
      for ( i = j+1; i < m; i++ )
      {
        if ( u[i+j*m] != 0.0 )
        {
          l[i+j*m] = u[i+j*m] / u[j+j*m];
          u[i+j*m] = 0.0;
          for ( k = j+1; k < n; k++ )
          {
            u[i+k*m] = u[i+k*m] - l[i+j*m] * u[j+k*m];
          }
        }
      }
    }
  }

  return;
}
/******************************************************************************/

double *r8ge_poly ( int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8GE_POLY computes the characteristic polynomial of a R8GE matrix.

  Discussion:

    The R8GE storage format is used for a "general" M by N matrix.  
    A physical storage space is made for each logical entry.  The two 
    dimensional logical array is mapped to a vector, in which storage is 
    by columns.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    28 February 2012

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of the matrix.
    N must be positive.

    Input, double A[N*N], the R8GE matrix.

    Output, double R8GE_POLY[N+1], the coefficients of the characteristic
    polynomial of A.  P(I) contains the coefficient of X**I.
*/
{
  int i;
  int j;
  int k;
  int order;
  double *p;
  double trace;
  double *work1;
  double *work2;

  p = ( double * ) malloc ( ( n + 1 ) * sizeof ( double ) );
  work2 = ( double * ) malloc ( n * n * sizeof ( double ) );
/*
  Initialize WORK1 to the identity matrix.
*/
  work1 = r8ge_identity ( n );

  p[n] = 1.0;

  for ( order = n-1; 0 <= order; order-- )
  {
/*
  Work2 = A * WORK1.
*/
    for ( j = 0; j < n; j++ )
    {
      for ( i = 0; i < n; i++ )
      {
        work2[i+j*n] = 0.0;
        for ( k = 0; k < n; k++ )
        {
          work2[i+j*n] = work2[i+j*n] + a[i+k*n] * work1[k+j*n];
        }
      }
    }
/*
  Take the trace.
*/
    trace = 0.0;
    for ( i = 0; i < n; i++ )
    {
      trace = trace + work2[i+i*n];
    }
/*
  P(ORDER) = - Trace ( WORK2 ) / ( N - ORDER )
*/
    p[order] = -trace / ( double ) ( n - order );
/*
  WORK1 := WORK2 + P(ORDER) * Identity.
*/
    for ( j = 0; j < n; j++ )
    {
      for ( i = 0; i < n; i++ )
      {
        work1[i+j*n] = work2[i+j*n];
      }
    }
    for ( j = 0; j < n; j++ )
    {
      work1[j+j*n] = work1[j+j*n] + p[order];
    }
  }

  free ( work1 );
  free ( work2 ) ;

  return p;
}
/******************************************************************************/

void r8ge_print ( int m, int n, double a[], char *title )

/******************************************************************************/
/*
  Purpose:

    R8GE_PRINT prints a R8GE matrix.

  Discussion:

    The R8GE storage format is used for a "general" M by N matrix.  
    A physical storage space is made for each logical entry.  The two 
    dimensional logical array is mapped to a vector, in which storage is 
    by columns.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    28 February 2012

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows of the matrix.
    M must be positive.

    Input, int N, the number of columns of the matrix.
    N must be positive.

    Input, double A[M*N], the R8GE matrix.

    Input, char *TITLE, a title.
*/
{
  r8ge_print_some ( m, n, a, 1, 1, m, n, title );

  return;
}
/******************************************************************************/

void r8ge_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi, 
  int jhi, char *title )

/******************************************************************************/
/*
  Purpose:

    R8GE_PRINT_SOME prints some of a R8GE matrix.

  Discussion:

    The R8GE storage format is used for a "general" M by N matrix.  
    A physical storage space is made for each logical entry.  The two 
    dimensional logical array is mapped to a vector, in which storage is 
    by columns.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    28 February 2012

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows of the matrix.
    M must be positive.

    Input, int N, the number of columns of the matrix.
    N must be positive.

    Input, double A[M*N], the R8GE matrix.

    Input, int ILO, JLO, IHI, JHI, designate the first row and
    column, and the last row and column to be printed.

    Input, char *TITLE, a title.
*/
{
# define INCX 5

  int i;
  int i2hi;
  int i2lo;
  int j;
  int j2hi;
  int j2lo;

  printf ( "\n" );
  printf ( "%s\n", title );
/*
  Print the columns of the matrix, in strips of 5.
*/
  for ( j2lo = jlo; j2lo <= jhi; j2lo = j2lo + INCX )
  {
    j2hi = j2lo + INCX - 1;
    j2hi = i4_min ( j2hi, n );
    j2hi = i4_min ( j2hi, jhi );

    printf ( "\n" );
/*
  For each column J in the current range...

  Write the header.
*/
    printf ( "  Col:    " );
    for ( j = j2lo; j <= j2hi; j++ )
    {
      printf ( "%7d       ", j );
    }
    printf ( "\n" );
    printf ( "  Row\n" );
    printf ( "  ---\n" );
/*
  Determine the range of the rows in this strip.
*/
    i2lo = i4_max ( ilo, 1 );
    i2hi = i4_min ( ihi, m );

    for ( i = i2lo; i <= i2hi; i++ )
    {
/*
  Print out (up to) 5 entries in row I, that lie in the current strip.
*/
      printf ( "%5d  ", i );
      for ( j = j2lo; j <= j2hi; j++ )
      {
        printf ( "%12g  ", a[i-1+(j-1)*m] );
      }
      printf ( "\n" );
    }
  }
  return;
# undef INCX
}
/******************************************************************************/

double *r8ge_random ( int m, int n, int *seed )

/******************************************************************************/
/*
  Purpose:

    R8GE_RANDOM randomizes a R8GE matrix.

  Discussion:

    The R8GE storage format is used for a "general" M by N matrix.  
    A physical storage space is made for each logical entry.  The two 
    dimensional logical array is mapped to a vector, in which storage is 
    by columns.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    28 November 2011

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows of the matrix.
    M must be positive.

    Input, int N, the number of columns of the matrix.
    N must be positive.

    Input/output, int *SEED, a seed for the random number generator.

    Output, double R8GE_RANDOM[M*N], the randomized M by N matrix, 
    with entries between 0 and 1.
*/
{
  double *a;
  int i;
  int j;

  a = ( double * ) malloc ( m * n * sizeof ( double ) );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      a[i+j*m] = r8_uniform_01 ( seed );
    }
  }
  return a;
}
/******************************************************************************/

double *r8ge_res ( int m, int n, double a[], double x[], double b[] )

/******************************************************************************/
/*
  Purpose:

    R8GE_RES computes the residual for a R8GE system.

  Discussion:

    The R8GE storage format is used for a "general" M by N matrix.  
    A physical storage space is made for each logical entry.  The two 
    dimensional logical array is mapped to a vector, in which storage is 
    by columns.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 March 2012

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of the matrix.
    N must be positive.

    Input, double A[M*N], the original, UNFACTORED matrix.

    Input, double X[N], an estimate of the solution the linear system.

    Input, double B[M], the right hand side vector.

    Output, double R8GE_RES[M], the residual vector: b - A * x.
*/
{
  double *r;
  int i;
  int j;

  r = ( double * ) malloc ( n * sizeof ( double ) );

  for ( i = 0; i < m; i++ )
  {
    r[i] = b[i];
    for ( j = 0; j < n; j++ )
    {
      r[i] = r[i] - a[i+j*m] * x[j];
    }
  }

  return r;
}
/******************************************************************************/

double *r8ge_sl ( int n, double a_lu[], int pivot[], double b[], int job )

/******************************************************************************/
/*
  Purpose:

    R8GE_SL solves a R8GE system factored by R8GE_FA.

  Discussion:

    The R8GE storage format is used for a "general" M by N matrix.  
    A physical storage space is made for each logical entry.  The two 
    dimensional logical array is mapped to a vector, in which storage is 
    by columns.

    R8GE_SL is a simplified version of the LINPACK routine SGESL.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 March 2012

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of the matrix.
    N must be positive.

    Input, double A_LU[N*N], the LU factors from R8GE_FA.

    Input, int PIVOT[N], the pivot vector from R8GE_FA.

    Input, double B[N], the right hand side vector.

    Input, int JOB, specifies the operation.
    0, solve A * x = b.
    nonzero, solve A' * x = b.

    Output, double R8GE_SL[N], the solution vector.
*/
{
  int i;
  int k;
  int l;
  double t;
  double *x;

  x = ( double * ) malloc ( n * sizeof ( double ) );

  for ( i = 0; i < n; i++ )
  {
    x[i] = b[i];
  }
/*
  Solve A * x = b.
*/
  if ( job == 0 )
  {
/*
  Solve PL * Y = B.
*/
    for ( k = 1; k <= n-1; k++ )
    {
      l = pivot[k-1];

      if ( l != k )
      {
        t      = x[l-1];
        x[l-1] = x[k-1];
        x[k-1] = t;
      }
      for ( i = k+1; i <= n; i++ )
      {
        x[i-1] = x[i-1] + a_lu[i-1+(k-1)*n] * x[k-1];
      }
    }
/*
  Solve U * X = Y.
*/
    for ( k = n; 1 <= k; k-- )
    {
      x[k-1] = x[k-1] / a_lu[k-1+(k-1)*n];
      for ( i = 1; i <= k-1; i++ )
      {
        x[i-1] = x[i-1] - a_lu[i-1+(k-1)*n] * x[k-1];
      }
    }
  }
/*
  Solve A' * X = B.
*/
  else
  {
/*
  Solve U' * Y = B.
*/
    for ( k = 1; k <= n; k++ )
    {
      t = 0.0;
      for ( i = 1; i <= k-1; i++ )
      {
        t = t + x[i-1] * a_lu[i-1+(k-1)*n];
      }
      x[k-1] = ( x[k-1] - t ) / a_lu[k-1+(k-1)*n];
    }
/*
  Solve ( PL )' * X = Y.
*/
    for ( k = n-1; 1 <= k; k-- )
    {
      t = 0.0;
      for ( i = k+1; i <= n; i++ )
      {
        t = t + x[i-1] * a_lu[i-1+(k-1)*n];
      }
      x[k-1] = x[k-1] + t;

      l = pivot[k-1];

      if ( l != k )
      {
        t      = x[l-1];
        x[l-1] = x[k-1];
        x[k-1] = t;
      }

    }

  }

  return x;
}
/******************************************************************************/

double *r8ge_sl_it ( int n, double a[], double a_lu[], int pivot[], double b[], 
  int job, double x[] )

/******************************************************************************/
/*
  Purpose:

    R8GE_SL_IT applies one step of iterative refinement following R8GE_SL.

  Discussion:

    The R8GE storage format is used for a "general" M by N matrix.  
    A physical storage space is made for each logical entry.  The two 
    dimensional logical array is mapped to a vector, in which storage is 
    by columns.

    It is assumed that:

    * the original matrix A has been factored by R8GE_FA;
    * the linear system A * x = b has been solved once by R8GE_SL.

    (Actually, it is not necessary to solve the system once using R8GE_SL.
    You may simply supply the initial estimated solution X = 0.)

    Each time this routine is called, it will compute the residual in
    the linear system, apply one step of iterative refinement, and
    add the computed correction to the current solution.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 March 2012

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of the matrix.
    N must be positive.

    Input, double A[N*N], the original, UNFACTORED R8GE matrix.

    Input, double A_LU[N*N], the LU factors from R8GE_FA.

    Input, int PIVOT[N], the pivot vector from R8GE_FA.

    Input, double B[N], the right hand side vector.

    Input, int JOB, specifies the operation.
    0, solve A*X=B.
    nonzero, solve A'*X=B.

    Input, double X[N], an estimate of the solution of A * x = b.

    Output, double R8GE_SL_IT[N], the solution after one step of 
    iterative refinement.
*/
{
  double *dx;
  int i;
  double *r;
  double *x_new;
/*
  Compute the residual vector.
*/
  r = r8ge_res ( n, n, a, x, b );
/*
  Solve A * dx = r
*/
  dx = r8ge_sl ( n, a_lu, pivot, r, job );
/*
  Add dx to x.
*/
  x_new = ( double * ) malloc ( n * sizeof ( double ) );

  for ( i = 0; i < n; i++ )
  {
    x_new[i] = x[i] + dx[i];
  }

  free ( r );
  free ( dx );

  return x_new;
}
/******************************************************************************/

double *r8ge_to_r8gb ( int m, int n, int ml, int mu, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8GE_TO_R8GB copies a R8GE matrix to a R8GB matrix.

  Discussion:

    The R8GE storage format is used for a "general" M by N matrix.  
    A physical storage space is made for each logical entry.  The two 
    dimensional logical array is mapped to a vector, in which storage is 
    by columns.

    The R8GB storage format is for an M by N banded matrix, with lower
    bandwidth ML and upper bandwidth MU.  Storage includes room for ML
    extra superdiagonals, which may be required to store nonzero entries
    generated during Gaussian elimination.

    It usually doesn't make sense to try to store a general matrix
    in a band matrix format.  You can always do it, but it will take
    more space, unless the general matrix is actually banded.

    The purpose of this routine is to allow a user to set up a
    banded matrix in the easy-to-use general format, and have this
    routine take care of the compression of the data into general
    format.  All the user has to do is specify the bandwidths.

    Note that this routine "believes" what the user says about the
    bandwidth.  It will assume that all entries in the general matrix
    outside of the bandwidth are zero.

    The original M by N matrix is "collapsed" downward, so that diagonals
    become rows of the storage array, while columns are preserved.  The
    collapsed array is logically 2*ML+MU+1 by N.

    LINPACK and LAPACK band storage requires that an extra ML
    superdiagonals be supplied to allow for fillin during Gauss
    elimination.  Even though a band matrix is described as
    having an upper bandwidth of MU, it effectively has an
    upper bandwidth of MU+ML.  This routine will copy nonzero
    values it finds in these extra bands, so that both unfactored
    and factored matrices can be handled.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 March 2012

  Author:

    John Burkardt

  Reference:

    Edward Anderson, Zhaojun Bai, Christian Bischof, Susan Blackford, 
    James Demmel, Jack Dongarra, Jeremy Du Croz, Anne Greenbaum, 
    Sven Hammarling, Alan McKenney, Danny Sorensen,
    LAPACK User's Guide,
    Second Edition,
    SIAM, 1995.

    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
    LINPACK User's Guide,
    SIAM, 1979,
    ISBN13: 978-0-898711-72-1,
    LC: QA214.L56.

  Parameters:

    Input, int M, the number of rows of the matrices.
    M must be positive.

    Input, int N, the number of columns of the matrices.
    N must be positive.

    Input, int ML, MU, the lower and upper bandwidths of A1.
    ML and MU must be nonnegative, and no greater than min(M,N)-1.

    Output, double A[M*N], the R8GE matrix.

    Input, double R8GE_TO_R8GB[(2*ML+MU+1)*N], the R8GB matrix.
*/
{
  double *b;
  int i;
  int j;
  int jhi;
  int jlo;
  int k;

  b = ( double * ) malloc ( (2*ml+mu+1) * n * sizeof ( double ) );
  for ( k = 0; k < (2*ml+mu+1)*n; k++ )
  {
    b[k] = 0.0;
  }

  for ( i = 1; i <= m; i++ )
  {
    jlo = i4_max ( i - ml, 1 );
    jhi = i4_min ( i + mu, n );

    for ( j = jlo; j <= jhi; j++ )
    {
      b[ml+mu+i-j+(j-1)*(2*ml+mu+1)] = a[i-1+(j-1)*m];
    }
  }

  return b;
}
/******************************************************************************/

double *r8ge_to_r8vec ( int m, int n, double *a )

/******************************************************************************/
/*
  Purpose:

    R8GE_TO_R8VEC copies a R8GE matrix to a real vector.

  Discussion:

    The R8GE storage format is used for a "general" M by N matrix.  
    A physical storage space is made for each logical entry.  The two 
    dimensional logical array is mapped to a vector, in which storage is 
    by columns.

    In C++  and FORTRAN, this routine is not really needed.  In MATLAB,
    a data item carries its dimensionality implicitly, and so cannot be
    regarded sometimes as a vector and sometimes as an array.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 March 2012

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns in the array.

    Input, double R8VEC_TO_R8GE[M*N], the array to be copied.

    Output, double X[M*N], the vector.
*/
{
  int i;
  int j;
  int k;
  double *x;

  x = ( double * ) malloc ( m * n * sizeof ( double ) );

  k = 0;
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      x[k] = a[i+j*m];
      k = k + 1;
    }
  }

  return x;
}
/******************************************************************************/

int r8ge_trf ( int m, int n, double a[], int pivot[] )

/******************************************************************************/
/*
  Purpose:

    R8GE_TRF performs a LAPACK-style PLU factorization of a R8GE matrix.

  Discussion:

    The R8GE storage format is used for a "general" M by N matrix.  
    A physical storage space is made for each logical entry.  The two 
    dimensional logical array is mapped to a vector, in which storage is 
    by columns.

    R8GE_TRF is a standalone version of the LAPACK routine SGETRF.

    The factorization uses partial pivoting with row interchanges,
    and has the form
      A = P * L * U
    where P is a permutation matrix, L is lower triangular with unit
    diagonal elements (lower trapezoidal if N < M), and U is upper
    triangular (upper trapezoidal if M < N).

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 March 2012

  Author:

    Original FORTRAN77 version by Anderson, Bai, Bischof, Blackford,
    Demmel, Dongarra, DuCroz, Greenbaum, Hammarling, McKenney, Sorensen.
    C version by John Burkardt

  Reference:

    Edward Anderson, Zhaojun Bai, Christian Bischof, Susan Blackford, 
    James Demmel, Jack Dongarra, Jeremy Du Croz, Anne Greenbaum, 
    Sven Hammarling, Alan McKenney, Danny Sorensen,
    LAPACK User's Guide,
    Second Edition,
    SIAM, 1995.

  Parameters:

    Input, int M, the number of rows of the matrix A.  0 <= M.

    Input, int N, the number of columns of the matrix A.  0 <= N.

    Input/output, double A[M*N].
    On entry, the M by N matrix to be factored.
    On exit, the factors L and U from the factorization
    A = P*L*U; the unit diagonal elements of L are not stored.

    Output, int PIVOT[min(M,N)], the pivot indices.

    Output, int R8GE_TRF.
    = 0: successful exit
    = -K, the K-th argument had an illegal value
    = K: U(K,K) is exactly zero. The factorization
         has been completed, but the factor U is exactly
         singular, and division by zero will occur if it is used
         to solve a system of equations.
*/
{
  int i;
  int ii;
  int info;
  int j;
  int jj;
  int jp;
  double temp;
/*
  Test the input parameters.
*/
  info = 0;

  if ( m < 0 )
  {
    return (-1);
  }
  else if ( n < 0 )
  {
    return (-2);
  }

  if ( m == 0 || n == 0 )
  {
    return 0;
  }

  for ( j = 1; j <= i4_min ( m, n ); j++ )
  {
/*
  Find the pivot.
*/
    temp = r8_abs ( a[j-1+(j-1)*m] );
    jp = j;
    for ( i = j+1; i <= m; i++ )
    {
      if ( temp < r8_abs ( a[i-1+(j-1)*m] ) )
      {
        temp = r8_abs ( a[i-1+(j-1)*m] );
        jp = i;
      }
    }

    pivot[j-1] = jp;
/*
  Apply the interchange to columns 1:N.
  Compute elements J+1:M of the J-th column.
*/
    if ( a[jp-1+(j-1)*m] != 0.0 )
    {
      if ( jp != j )
      {
        for ( jj = 1; jj <= n; jj++ )
        {
          temp             = a[j-1+(jj-1)*m];
          a[j-1+(jj-1)*m]  = a[jp-1+(jj-1)*m];
          a[jp-1+(jj-1)*m] = temp;
        }
      }

      if ( j < m )
      {
        for ( i = j+1; i <= m; i++ )
        {
          a[i-1+(j-1)*m] = a[i-1+(j-1)*m] / a[j-1+(j-1)*m];
        }
      }
    }
    else if ( info == 0 )
    {
      info = j;
    }
/*
  Update the trailing submatrix.
*/
    if ( j < i4_min ( m, n ) )
    {
      for ( ii = j+1; ii <= m; ii++ )
      {
        for ( i = j+1; i <= n; i++ )
        {
          a[ii-1+(i-1)*m] = a[ii-1+(i-1)*m] - a[ii-1+(j-1)*m] * a[j-1+(i-1)*m];
        }
      }
    }
  }

  return info;
}
/******************************************************************************/

double *r8ge_trs ( int n, int nrhs, char trans, double a[], int pivot[], 
  double b[] )

/******************************************************************************/
/*
  Purpose:

    R8GE_TRS solves a system of linear equations factored by R8GE_TRF.

  Discussion:

    The R8GE storage format is used for a "general" M by N matrix.  
    A physical storage space is made for each logical entry.  The two 
    dimensional logical array is mapped to a vector, in which storage is 
    by columns.

    R8GE_TRS is a standalone version of the LAPACK routine SGETRS.

    R8GE_TRS solves a system of linear equations
      A * x = b  or  A' * X = B
    with a general N by N matrix A using the PLU factorization computed
    by R8GE_TRF.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 March 2012

  Author:

    Original FORTRAN77 version by Anderson, Bai, Bischof, Blackford,
    Demmel, Dongarra, DuCroz, Greenbaum, Hammarling, McKenney, Sorensen.
    C version by John Burkardt

  Reference:

    Edward Anderson, Zhaojun Bai, Christian Bischof, Susan Blackford, 
    James Demmel, Jack Dongarra, Jeremy Du Croz, Anne Greenbaum, 
    Sven Hammarling, Alan McKenney, Danny Sorensen,
    LAPACK User's Guide,
    Second Edition,
    SIAM, 1995.

  Parameters:

    Input, int N, the order of the matrix A.  0 <= N.

    Input, int NRHS, the number of right hand sides.  0 <= NRHS.

    Input, char TRANS, specifies the form of the system of equations:
    'N':  A * x = b  (No transpose)
    'T':  A'* X = B  (Transpose)
    'C':  A'* X = B  (Conjugate transpose = Transpose)

    Input, double A[N*N], the factors L and U from the factorization
    A = P*L*U as computed by R8GE_TRF.

    Input, int PIVOT[N], the pivot indices from R8GE_TRF.

    Input, double B[N*NRHS], the right hand side matrix.

    Output, double R8GE_TRS[N*NRHS], the solution matrix X.
*/
{
  int i;
  int j;
  int k;
  double temp;
  double *x;

  if ( trans != 'n' && trans != 'N' && 
       trans != 't' && trans != 'T' && 
       trans != 'c' && trans != 'C' )
  {
    return NULL;
  }

  if ( n < 0 )
  {
    return NULL;
  }

  if ( nrhs < 0 )
  {
    return NULL;
  }

  if ( n == 0 )
  {
    return NULL;
  }
  if ( nrhs == 0 )
  {
    return NULL;
  }

  x = ( double * ) malloc ( n * nrhs * sizeof ( double ) );
  for ( k = 0; k < nrhs; k++ )
  {
    for ( i = 0; i < n; i++ )
    {
      x[i+k*n] = b[i+k*n];
    }
  }

  if ( trans == 'n' || trans == 'N' )
  {
/*
  Apply row interchanges to the right hand sides.
*/
    for ( i = 1; i <= n; i++ )
    {
      if ( pivot[i-1] != i )
      {
        for ( k = 0; k < nrhs; k++ )
        {
          temp                = x[i-1+k*n];
          x[i-1+k*n]          = x[pivot[i-1]-1+k*n];
          x[pivot[i-1]-1+k*n] = temp;
        }
      }
    }
/*
  Solve L * x = b, overwriting b with x.
*/
    for ( k = 0; k < nrhs; k++ )
    {
      for ( j = 1; j <= n-1; j++ )
      {
        for ( i = j+1; i <= n; i++ )
        {
          x[i-1+k*n] = x[i-1+k*n] - a[i-1+(j-1)*n] * x[j-1+k*n];
        }
      }
    }
/*
  Solve U * x = b, overwriting b with x.
*/
    for ( k = 0; k < nrhs; k++ )
    {
      for ( j = n; 1 <= j; j-- )
      {
        x[j-1+k*n] = x[j-1+k*n] / a[j-1+(j-1)*n];
        for ( i = 1; i < j; i++ )
        {
          x[i-1+k*n] = x[i-1+k*n] - a[i-1+(j-1)*n] * x[j-1+k*n];
        }
      }
    }
  }
  else
  {
/*
  Solve U' * x = b, overwriting b with x.
*/
    for ( k = 0; k < nrhs; k++ )
    {
      for ( j = 1; j <= n; j++ )
      {
        x[j-1+k*n] = x[j-1+k*n] / a[j-1+(j-1)*n];
        for ( i = j+1; i <= n; i++ )
        {
          x[i-1+k*n] = x[i-1+k*n] - a[j-1+(i-1)*n] * x[j-1+k*n];
        }
      }
    }
/*
  Solve L' * x = b, overwriting b with x.
*/
    for ( k = 0; k < nrhs; k++ )
    {
      for ( j = n; 2 <= j; j-- )
      {
        for ( i = 1; i < j; i++ )
        {
          x[i-1+k*n] = x[i-1+k*n] - a[j-1+(i-1)*n] * x[j-1+k*n];
        }
      }
    }
/*
  Apply row interchanges to the solution vectors.
*/
    for ( i = n; 1 <= i; i-- )
    {
      if ( pivot[i-1] != i )
      {
        for ( k = 0; k < nrhs; k++ )
        {
          temp                = x[i-1+k*n];
          x[i-1+k*n]          = x[pivot[i-1]-1+k*n];
          x[pivot[i-1]-1+k*n] = temp;
        }
      }
    }
  }

  return x;
}
/******************************************************************************/

double *r8ge_vxm ( int m, int n, double a[], double x[] )

/******************************************************************************/
/*
  Purpose:

    R8GE_VXM multiplies a vector times a R8GE matrix.

  Discussion:

    The R8GE storage format is used for a "general" M by N matrix.  
    A physical storage space is made for each logical entry.  The two 
    dimensional logical array is mapped to a vector, in which storage is 
    by columns.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 March 2012

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows of the matrix.
    M must be positive.

    Input, int N, the number of columns of the matrix.
    N must be positive.

    Input, double A[M*N], the R8GE matrix.

    Input, double X[M], the vector to be multiplied by A.

    Output, double R8GE_VXM[N], the product A' * x.
*/
{
  double *b;
  int i;
  int j;

  b = ( double * ) malloc ( n * sizeof ( double ) );

  for ( i = 0; i < n; i++ )
  {
    b[i] = 0.0;
    for ( j = 0; j < m; j++ )
    {
      b[i] = b[i] + a[j+i*m] * x[j];
    }
  }

  return b;
}
/******************************************************************************/

double *r8ge_zero ( int m, int n )

/******************************************************************************/
/*
  Purpose:

    R8GE_ZERO zeros a R8GE matrix.

  Discussion:

    The R8GE storage format is used for a "general" M by N matrix.  
    A physical storage space is made for each logical entry.  The two 
    dimensional logical array is mapped to a vector, in which storage is 
    by columns.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 March 2012

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows of the matrix.
    M must be positive.

    Input, int N, the number of columns of the matrix.
    N must be positive.

    Output, double R8GE_ZERO[M*N], the M by N matrix.
*/
{
  double *a;
  int i;
  int j;

  a = ( double * ) malloc ( m * n * sizeof ( double ) );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      a[i+j*m] = 0.0;
    }
  }

  return a;
}
/******************************************************************************/

double r8lt_det ( int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8LT_DET computes the determinant of a R8LT matrix.

  Discussion:

    The R8LT storage format is used for an M by N lower triangular 
    matrix A, and allocates storage even for the zero entries.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    07 March 2012

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of the matrix.
    N must be positive.

    Input, double A[N*N], the R8LT matrix.

    Output, double R8LT_DET, the determinant of the matrix.
*/
{
  double det;
  int i;

  det = 1.0;
  for ( i = 0; i < n; i++ )
  {
    det = det * a[i+i*n];
  }

  return det;
}
/******************************************************************************/

double *r8lt_indicator ( int m, int n )

/******************************************************************************/
/*
  Purpose:

    R8LT_INDICATOR sets up a R8LT indicator matrix.

  Discussion:

    The R8LT storage format is used for an M by N lower triangular 
    matrix A, and allocates storage even for the zero entries.  

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    07 March 2012

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns of the matrix.
    M and N must be positive.

    Output, double R8LT_INDICATOR[M*N], the R8LT matrix.
*/
{
  double *a;
  int fac;
  int i;
  int j;

  a = ( double * ) malloc ( m * n * sizeof ( double ) );

  fac = i4_power ( 10, i4_log_10 ( n ) + 1 );

  for ( i = 1; i <= m; i++ )
  {
    for ( j = 1; j <= i4_min ( i, n ); j++ )
    {
      a[i-1+(j-1)*m] = ( double ) ( fac * i + j );
    }
    for ( j = i+1; j <= n; j++ )
    {
      a[i-1+(j-1)*m] = 0.0;
    }
  }

  return a;
}
/******************************************************************************/

double *r8lt_inverse ( int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8LT_INVERSE computes the inverse of a R8LT matrix.

  Discussion:

    The R8LT storage format is used for an M by N lower triangular 
    matrix A, and allocates storage even for the zero entries.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    07 March 2012

  Author:

    John Burkardt

  Reference:

    Albert Nijenhuis, Herbert Wilf,
    Combinatorial Algorithms,
    Second edition,
    Academic Press, 1978,
    ISBN 0-12-519260-6

  Parameters:

    Input, int N, the order of the matrix.

    Input, double A[N*N], the R8LT matrix.

    Output, double R8LT_INVERSE[N*N], the inverse of the matrix.
*/
{
  double *b;
  int i;
  int j;
  int k;
  double t;
/*
  Check.
*/
  for ( i = 0; i < n; i++ )
  {
    if ( a[i+i*n] == 0.0 )
    {
      printf ( "\n" );
      printf ( "R8LT_INVERSE - Fatal error!\n" );
      printf ( "  Zero diagonal element.\n" );
      exit ( 1 );
    }
  }

  b = ( double * ) malloc ( n * n * sizeof ( double ) );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      b[i+j*n] = a[i+j*n];
    }
  }

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      if ( i < j )
      {
        b[i+j*n] = 0.0;
      }
      else if ( i == j )
      {
        b[i+j*n] = 1.0 / b[i+j*n];
      }
      else if ( j < i )
      {
        t = 0.0;
        for ( k = j; k <= i-1; k++ )
        {
          t = t - b[i+k*n] * b[k+j*n];
        }
        b[i+j*n] = t / b[i+i*n];
      }
    }
  }

  return b;
}
/******************************************************************************/

double *r8lt_mxm ( int n, double a[], double b[] )

/******************************************************************************/
/*
  Purpose:

    R8LT_MXM multiplies two R8LT matrices.

  Discussion:

    The R8LT storage format is used for an M by N lower triangular 
    matrix A, and allocates storage even for the zero entries.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    07 March 2012

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of the matrices.
    N must be positive.

    Input, double A[N*N], B[N*N], the R8LT factor matrices.

    Output, double R8LT_MXM[N*N], the R8LT product matrix.
*/
{
  double *c;
  int i;
  int j;
  int k;

  c = ( double * ) malloc ( n * n * sizeof ( double ) );

  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      c[i+j*n] = 0.0;
      for ( k = j; k <= i; k++ )
      {
        c[i+j*n] = c[i+j*n] + a[i+k*n] * b[k+j*n];
      }
    }
  }

  return c;
}
/******************************************************************************/

double *r8lt_mxv ( int m, int n, double a[], double x[] )

/******************************************************************************/
/*
  Purpose:

    R8LT_MXV multiplies a R8LT matrix times a vector.

  Discussion:

    The R8LT storage format is used for an M by N lower triangular 
    matrix A, and allocates storage even for the zero entries.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    07 March 2012

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows of the matrix.
    M must be positive.

    Input, int N, the number of columns of the matrix.
    N must be positive.

    Input, double A[M*N], the R8LT matrix.

    Input, double X[N], the vector to be multiplied by A.

    Output, double R8LT_MXV[M], the product A * x.
*/
{
  double *b;
  int i;
  int j;
  int jmax;

  b = ( double * ) malloc ( m * sizeof ( double ) );

  for ( i = 0; i < m; i++ )
  {
    b[i] = 0.0;
    jmax = i4_min ( i, n-1 );
    for ( j = 0; j <= jmax; j++ )
    {
      b[i] = b[i] + a[i+j*m] * x[j];
    }
  }

  return b;
}
/******************************************************************************/

void r8lt_print ( int m, int n, double a[], char *title )

/******************************************************************************/
/*
  Purpose:

    R8LT_PRINT prints a R8LT matrix.

  Discussion:

    The R8LT storage format is used for an M by N lower triangular 
    matrix A, and allocates storage even for the zero entries.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    07 March 2012

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows of the matrix.
    M must be positive.

    Input, int N, the number of columns of the matrix.
    N must be positive.

    Input, double A[M*N], the R8LT matrix.

    Input, char *TITLE, a title.
*/
{
  r8lt_print_some ( m, n, a, 1, 1, m, n, title );

  return;
}
/******************************************************************************/

void r8lt_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi, 
  int jhi, char *title )

/******************************************************************************/
/*
  Purpose:

    R8LT_PRINT_SOME prints some of a R8LT matrix.

  Discussion:

    The R8LT storage format is used for an M by N lower triangular 
    matrix A, and allocates storage even for the zero entries.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    07 March 2012

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows of the matrix.
    M must be positive.

    Input, int N, the number of columns of the matrix.
    N must be positive.

    Input, double A[M*N], the R8LT matrix.

    Input, int ILO, JLO, IHI, JHI, designate the first row and
    column, and the last row and column to be printed.

    Input, char *TITLE, a title.
*/
{
# define INCX 5

  int i;
  int i2hi;
  int i2lo;
  int j;
  int j2hi;
  int j2lo;

  printf ( "\n" );
  printf ( "%s\n", title );

  if ( ilo < jlo )
  {
    return;
  }
/*
  Print the columns of the matrix, in strips of 5.
*/
  for ( j2lo = jlo; j2lo <= jhi; j2lo = j2lo + INCX )
  {
    j2hi = j2lo + INCX - 1;
    j2hi = i4_min ( j2hi, n );
    j2hi = i4_min ( j2hi, jhi );

    printf ( "\n" );
    printf ( "  Col: " );

    for ( j = j2lo; j <= j2hi; j++ )
    {
      printf ( "%7d       ", j );
    }
    printf ( "\n" );
    printf ( "  Row\n" );
    printf ( "  ---\n" );
/*
  Determine the range of the rows in this strip.
*/
    i2lo = i4_max ( ilo, 1 );
    i2lo = i4_max ( i2lo, j2lo );

    i2hi = i4_min ( ihi, m );

    for ( i = i2lo; i <= i2hi; i++ )
    {
/*
  Print out (up to) 5 entries in row I, that lie in the current strip.
*/
      printf ( "%4d  ", i );

      for ( j = j2lo; j <= j2hi; j++ )
      {
        if ( i < j )
        {
          printf ( "              " );
        }
        else
        {
          printf ( "%12g  ", a[i-1+(j-1)*m] );
        }
      }
      printf ( "\n" );
    }
  }

  return;
# undef INCX
}
/******************************************************************************/

double *r8lt_random ( int m, int n, int *seed )

/******************************************************************************/
/*
  Purpose:

    R8LT_RANDOM randomizes a R8LT matrix.

  Discussion:

    The R8LT storage format is used for an M by N lower triangular 
    matrix A, and allocates storage even for the zero entries.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    07 March 2012

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns of the matrix.
    M and N must be positive.

    Input/output, int SEED, a seed for the random number generator.

    Output, double R8LT_RANDOM[M*N], the R8LT matrix.
*/
{
  double *a;
  int i;
  int j;

  a = ( double * ) malloc ( m * n * sizeof ( double ) );

  for ( i = 0; i < m; i++ )
  {
    for ( j = 0; j <= i4_min ( i, n-1); j++ )
    {
      a[i+j*m] = r8_uniform_01 ( seed );
    }
    for ( j = i+1; j < n; j++ )
    {
      a[i+j*m] = 0.0;
    }
  }

  return a;
}
/******************************************************************************/

double *r8lt_sl ( int n, double a[], double b[], int job )

/******************************************************************************/
/*
  Purpose:

    R8LT_SL solves a R8LT system.

  Discussion:

    The R8LT storage format is used for an M by N lower triangular 
    matrix A, and allocates storage even for the zero entries.

    No factorization of the lower triangular matrix is required.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    07 March 2012

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of the matrix.

    Input, double A[N*N], the R8LT matrix.

    Input, double B[N], the right hand side.

    Input, int JOB, is 0 to solve the untransposed system,
    nonzero to solve the transposed system.

    Output, double R8LT_SL[N], the solution vector.
*/
{
  int i;
  int j;
  double *x;

  x = ( double * ) malloc ( n * sizeof ( double ) );

  for ( i = 0; i < n; i++ )
  {
    x[i] = b[i];
  }

  if ( job == 0 )
  {
    for ( j = 0; j < n; j++ )
    {
      x[j] = x[j] / a[j+j*n];
      for ( i = j+1; i < n; i++ )
      {
        x[i] = x[i] - a[i+j*n] * x[j];
      }
    }
  }
  else
  {
    for ( j = n-1; 0 <= j; j-- )
    {
      x[j] = x[j] / a[j+j*n];
      for ( i = 0; i < j; i++ )
      {
        x[i] = x[i] - a[j+i*n] * x[j];
      }
    }
  }

  return x;
}
/******************************************************************************/

double *r8lt_vxm ( int m, int n, double a[], double x[] )

/******************************************************************************/
/*
  Purpose:

    R8LT_VXM multiplies a vector times a R8LT matrix.

  Discussion:

    The R8LT storage format is used for an M by N lower triangular 
    matrix A, and allocates storage even for the zero entries.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    07 March 2012

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows of the matrix.
    M must be positive.

    Input, int N, the number of columns of the matrix.
    N must be positive.

    Input, double A[M*N], the R8LT matrix.

    Input, double X[M], the vector to be multiplied by A.

    Output, double R8LT_VXM[N], the product A * x.
*/
{
  double *b;
  int i;
  int j;

  b = ( double * ) malloc ( n * sizeof ( double ) );

  for ( j = 0; j < n; j++)
  {
    b[j] = 0.0;
    for ( i = j; i < m; i++ )
    {
      b[j] = b[j] + x[i] * a[i+j*m];
    }
  }

  return b;
}
/******************************************************************************/

double *r8lt_zero ( int m, int n )

/******************************************************************************/
/*
  Purpose:

    R8LT_ZERO zeros a R8LT matrix.

  Discussion:

    The R8LT storage format is used for an M by N lower triangular 
    matrix A, and allocates storage even for the zero entries.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    07 March 2012

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows of the matrix.
    M must be positive.

    Input, int N, the number of columns of the matrix.
    N must be positive.

    Output, double R8LT_ZERO[M*N], the R8LT matrix.
*/
{
  double *a;
  int i;
  int j;

  a = ( double * ) malloc ( m * n * sizeof ( double ) );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      a[i+j*m] = 0.0;
    }
  }

  return a;
}
/******************************************************************************/

void r8mat_print ( int m, int n, double a[], char *title )

/******************************************************************************/
/*
  Purpose:

    R8MAT_PRINT prints an R8MAT.

  Discussion: 							    

    An R8MAT is a doubly dimensioned array of R8's, which
    may be stored as a vector in column-major order.

    Entry A(I,J) is stored as A[I+J*M]

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    28 May 2008

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows in A.

    Input, int N, the number of columns in A.

    Input, double A[M*N], the M by N matrix.

    Input, char *TITLE, a title.
*/
{
  r8mat_print_some ( m, n, a, 1, 1, m, n, title );

  return;
}
/******************************************************************************/

void r8mat_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi, 
  int jhi, char *title )

/******************************************************************************/
/*
  Purpose:

    R8MAT_PRINT_SOME prints some of an R8MAT.

  Discussion: 							    

    An R8MAT is a doubly dimensioned array of R8's, which
    may be stored as a vector in column-major order.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    28 May 2008

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows of the matrix.
    M must be positive.

    Input, int N, the number of columns of the matrix.
    N must be positive.

    Input, double A[M*N], the matrix.

    Input, int ILO, JLO, IHI, JHI, designate the first row and
    column, and the last row and column to be printed.

    Input, char *TITLE, a title.
*/
{
# define INCX 5

  int i;
  int i2hi;
  int i2lo;
  int j;
  int j2hi;
  int j2lo;

  fprintf ( stdout, "\n" );
  fprintf ( stdout, "%s\n", title );
/*
  Print the columns of the matrix, in strips of 5.
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
    fprintf ( stdout, "  Col:  ");
    for ( j = j2lo; j <= j2hi; j++ )
    {
      fprintf ( stdout, "  %7d     ", j );
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
  Print out (up to) 5 entries in row I, that lie in the current strip.
*/
      fprintf ( stdout, "%5d:", i );
      for ( j = j2lo; j <= j2hi; j++ )
      {
        fprintf ( stdout, "  %14f", a[i-1+(j-1)*m] );
      }
      fprintf ( stdout, "\n" );
    }
  }

  return;
# undef INCX
}
/******************************************************************************/

double *r8vec_indicator_new ( int n )

/******************************************************************************/
/*
  Purpose:

    R8VEC_INDICATOR_NEW sets an R8VEC to the indicator vector {1,2,3...}.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    26 August 2008

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of elements of A.

    Output, double R8VEC_INDICATOR_NEW[N], the array.
*/
{
  double *a;
  int i;

  a = ( double * ) malloc ( n * sizeof ( double ) );

  for ( i = 0; i <= n-1; i++ ) 
  {
    a[i] = ( double ) ( i + 1 );
  }

  return a;
}
/******************************************************************************/

void r8vec_print ( int n, double a[], char *title )

/******************************************************************************/
/*
  Purpose:

    R8VEC_PRINT prints an R8VEC.

  Discussion:

    An R8VEC is a vector of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    08 April 2009

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of components of the vector.

    Input, double A[N], the vector to be printed.

    Input, char *TITLE, a title.
*/
{
  int i;

  printf ( "\n" );
  printf ( "%s\n", title );
  printf ( "\n" );

  for ( i = 0; i < n; i++ ) 
  {
    printf ( "  %8d  %14f\n", i, a[i] );
  }

  return;
}
/******************************************************************************/

void r8vec_print_some ( int n, double a[], int i_lo, int i_hi, char *title )

/******************************************************************************/
/*
  Purpose:

    R8VEC_PRINT_SOME prints "some" of an R8VEC.

  Discussion:

    An R8VEC is a vector of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 October 2006

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries of the vector.

    Input, double A[N], the vector to be printed.

    Input, integer I_LO, I_HI, the first and last indices to print.
    The routine expects 1 <= I_LO <= I_HI <= N.

    Input, char *TITLE, a title.
*/
{
  int i;

  printf ( "\n" );
  printf ( "%s\n", title );
  printf ( "\n" );

  for ( i = i4_max ( 1, i_lo ); i <= i4_min ( n, i_hi ); i++ )
  {
    printf ( "  %8d  %14f\n", i, a[i-1] );
  }

  return;
}
/******************************************************************************/

double *r8vec_uniform ( int n, double b, double c, int *seed )

/******************************************************************************/
/*
  Purpose:

    R8VEC_UNIFORM returns a scaled pseudorandom R8VEC.

  Discussion:

    This routine implements the recursion

      seed = 16807 * seed mod ( 2**31 - 1 )
      unif = seed / ( 2**31 - 1 )

    The integer arithmetic never requires more than 32 bits,
    including a sign bit.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    30 January 2005

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

    Input, int N, the number of entries in the vector.

    Input, double B, C, the lower and upper limits of the pseudorandom values.

    Input/output, int *SEED, a seed for the random number generator.

    Output, double R8VEC_UNIFORM_01[N], the vector of pseudorandom values.
*/
{
  int i;
  int i4_huge = 2147483647;
  int k;
  double *r;

  if ( *seed == 0 )
  {
    printf ( "\n" );
    printf ( "R8VEC_UNIFORM - Fatal error!\n" );
    printf ( "  Input value of SEED = 0.\n" );
    exit ( 1 );
  }

  r = malloc ( n * sizeof ( double ) );

  for ( i = 0; i < n; i++ )
  {
    k = *seed / 127773;

    *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

    if ( *seed < 0 )
    {
      *seed = *seed + i4_huge;
    }

    r[i] = b + ( c - b ) * ( double ) ( *seed ) * 4.656612875E-10;
  }

  return r;
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

void sort_heap_external ( int n, int *indx, int *i, int *j, int isgn )

/******************************************************************************/
/*
  Purpose:

    SORT_HEAP_EXTERNAL externally sorts a list of items into ascending order.

  Discussion:

    The actual list is not passed to the routine.  Hence it may
    consist of integers, reals, numbers, names, etc.  The user,
    after each return from the routine, will be asked to compare or
    interchange two items.

    The current version of this code mimics the FORTRAN version,
    so the values of I and J, in particular, are FORTRAN indices.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    05 February 2004

  Author:

    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
    C version by John Burkardt.

  Reference:

    Albert Nijenhuis, Herbert Wilf,
    Combinatorial Algorithms,
    Academic Press, 1978, second edition,
    ISBN 0-12-519260-6.

  Parameters:

    Input, int N, the length of the input list.

    Input/output, int *INDX.
    The user must set INDX to 0 before the first call.
    On return,
      if INDX is greater than 0, the user must interchange
      items I and J and recall the routine.
      If INDX is less than 0, the user is to compare items I
      and J and return in ISGN a negative value if I is to
      precede J, and a positive value otherwise.
      If INDX is 0, the sorting is done.

    Output, int *I, *J.  On return with INDX positive,
    elements I and J of the user's list should be
    interchanged.  On return with INDX negative, elements I
    and J are to be compared by the user.

    Input, int ISGN. On return with INDX negative, the
    user should compare elements I and J of the list.  If
    item I is to precede item J, set ISGN negative,
    otherwise set ISGN positive.
*/
{
  static int i_save = 0;
  static int j_save = 0;
  static int k = 0;
  static int k1 = 0;
  static int n1 = 0;
/*
  INDX = 0: This is the first call.
*/
  if ( *indx == 0 )
  {

    i_save = 0;
    j_save = 0;
    k = n / 2;
    k1 = k;
    n1 = n;
  }
/*
  INDX < 0: The user is returning the results of a comparison.
*/
  else if ( *indx < 0 )
  {
    if ( *indx == -2 ) 
    {
      if ( isgn < 0 ) 
      {
        i_save = i_save + 1;
      }
      j_save = k1;
      k1 = i_save;
      *indx = -1;
      *i = i_save;
      *j = j_save;
      return;
    }

    if ( 0 < isgn ) 
    {
      *indx = 2;
      *i = i_save;
      *j = j_save;
      return;
    }

    if ( k <= 1 ) 
    {
      if ( n1 == 1 ) 
      {
        i_save = 0;
        j_save = 0;
        *indx = 0;
      }
      else 
      {
        i_save = n1;
        j_save = 1;
        n1 = n1 - 1;
        *indx = 1;
      }
      *i = i_save;
      *j = j_save;
      return;
    }
    k = k - 1;
    k1 = k;
  }
/*
  0 < INDX: the user was asked to make an interchange.
*/
  else if ( *indx == 1 ) 
  {
    k1 = k;
  }

  for ( ; ; )
  {

    i_save = 2 * k1;

    if ( i_save == n1 ) 
    {
      j_save = k1;
      k1 = i_save;
      *indx = -1;
      *i = i_save;
      *j = j_save;
      return;
    }
    else if ( i_save <= n1 ) 
    {
      j_save = i_save + 1;
      *indx = -2;
      *i = i_save;
      *j = j_save;
      return;
    }

    if ( k <= 1 ) 
    {
      break;
    }

    k = k - 1;
    k1 = k;
  }

  if ( n1 == 1 ) 
  {
    i_save = 0;
    j_save = 0;
    *indx = 0;
    *i = i_save;
    *j = j_save;
  }
  else 
  {
    i_save = n1;
    j_save = 1;
    n1 = n1 - 1;
    *indx = 1;
    *i = i_save;
    *j = j_save;
  }

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

  printf ( "%s\n", time_buffer );

  return;
# undef TIME_SIZE
}
