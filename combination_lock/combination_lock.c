# include <stdlib.h>
# include <stdio.h>
# include <stdbool.h>
# include <math.h>
# include <time.h>

# include "combination_lock.h"

/******************************************************************************/

int bicycle_lock ( int c )

/******************************************************************************/
/*
  Purpose:

    BICYCLE_LOCK finds the combination on a typical bicycle lock.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    13 May 2012

  Author:

    John Burkardt

  Parameters:

    Input, int C, the combination, a value between 0 and 999.

    Output, int BICYCLE_LOCK, the step on which the combination 
    was found.  A value of -1 means the combination was not found.
*/
{
  int a;
  int step;

  step = -1;

  for ( a = 0; a <= 999; a++ )
  {
    if ( a == c )
    {
      step = a + 1;
      break;
    }
  }
  return step;
}
/******************************************************************************/

int combination_lock ( int m, int n, int c[] )

/******************************************************************************/
/*
  Purpose:

    COMBINATION_LOCK determines the combination of a lock.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    13 May 2012

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of dials.

    Input, int N, the number of symbols on each dial.
    We assume the symbols are the integers 0 to N-1.

    Input, int C[M], the combination.

    Output, int COMBINATION_LOCK, the step on which the combination 
    was found.  A value of -1 means the combination was not found.
*/
{
  int *a;
  int i;
  bool more;
  int step;
/*
  Starting with the guess (0, 0, ... 0),
  generate every possible combination, in order, and try it.
*/
  more = false;
  a = ( int * ) malloc ( m * sizeof ( int ) );
  step = 0;

  while ( 1 )
  {
    combination_next ( m, n, a, &more );

    if ( !more )
    {
      step = -1;
      break;
    }

    step = step + 1;

    if ( i4vec_eq ( m, a, c ) )
    {
      break;
    }
  }

  free ( a );

  return step;
}
/******************************************************************************/

void combination_next ( int m, int base, int a[], bool *more )

/******************************************************************************/
/*
  Purpose:

    COMBINATION_NEXT generates lock combinations in lex order.

  Discussion:

    The vectors are produced in lexical order, starting with
    (0,0,...,0),
    (0,0,...,1),
    ...
    (BASE-1,BASE-1,...,BASE-1).

  Example:

    M = 2,
    BASE = 3

    0   0
    0   1
    0   2
    1   0
    1   1
    1   2
    2   0
    2   1
    2   2

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 May 2012

  Author:

    John Burkardt

  Reference:

    Dennis Stanton, Dennis White,
    Constructive Combinatorics,
    Springer, 1986,
    ISBN: 0387963472,
    LC: QA164.S79.

  Parameters:

    Input, int M, the size of the vectors to be used.

    Input, int BASE, the base to be used.  BASE = 2 will
    give vectors of 0's and 1's, for instance.

    Input/output, int A[M].  The input value of A is
    not important on the first call.  Thereafter, it should simply be the 
    output value from the previous call.  The output value is the next vector
    in the sequence.

    Input/output, bool *MORE.  The input value should be FALSE on the first 
    call, and TRUE on subsequent calls.  The output value will be TRUE as long 
    as the next vector could be computed, and FALSE once there are no more.
*/
{
  int i;

  if ( !(*more) )
  {
    for ( i = 0; i < m; i++ )
    {
      a[i] = 0;
    }
    *more = true;
  }
  else   
  {
    for ( i = m - 1; 0 <= i; i-- )
    {
      a[i] = a[i] + 1;

      if ( a[i] < base )
      {
        return;
      }
      a[i] = 0;
    }
    *more = false;
  }
  return;
}
/******************************************************************************/

int get_seed ( )

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
  Remap SEED from [1,43200] to [1,Huge].
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
      fprintf ( stderr, "\n" );
      fprintf ( stderr, "I4_POWER - Fatal error!\n" );
      fprintf ( stderr, "  I^J requested, with I = 0 and J negative.\n" );
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
      fprintf ( stderr, "\n" );
      fprintf ( stderr, "I4_POWER - Fatal error!\n" );
      fprintf ( stderr, "  I^J requested, with I = 0 and J = 0.\n" );
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

bool i4vec_eq ( int n, int a1[], int a2[] )

/******************************************************************************/
/*
  Purpose:

    I4VEC_EQ is true if two I4VEC's are equal.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 May 2012

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the vectors.

    Input, int A1[N], A2[N], two vectors to compare.

    Output, bool I4VEC_EQ, is TRUE if every pair of elements A1(I) and A2(I) are equal,
    and FALSE otherwise.
*/
{
  int i;
  int value;

  for ( i = 0; i < n; i++ )
  {
    if ( a1[i] != a2[i] )
    {
      value = false;
      return value;
    }
  }
  value = true;

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

float r4_abs ( float x )

/******************************************************************************/
/*
  Purpose:

    R4_ABS returns the absolute value of an R4.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 May 2006

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
    value = x;
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
