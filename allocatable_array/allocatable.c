# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>

int main ( );
void test01 ( int *n, int **a );
int i4_uniform_ab ( int a, int b, int *seed );
void timestamp ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for ALLOCATABLE.

  Discussion:

    ALLOCATABLE demonstrates how a C function can call another function,
    which allocates memory to an array, assigns values to it, and returns
    it in the argument list.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    25 February 2014

  Author:

    John Burkardt
*/
{
  int *a;
  int i;
  int n;

  timestamp ( );
  printf ( "\n" );
  printf ( "ALLOCATABLE:\n" );
  printf ( "  C version\n" );
  printf ( "  Show how a function can declare a pointer to an array,\n" );
  printf ( "  but then call a function to allocate and initialize\n" );
  printf ( "  the array, which is then returned in the argument list.\n" );

  test01 ( &n, &a );

  printf ( "\n" );
  printf ( "  The array size is N = %d\n", n );
  printf ( "\n" );
  printf ( "  The array contents:\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %8d  %8d\n", i, a[i] );
  }
/*
  Free memory.
*/
  free ( a );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "ALLOCATABLE:\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void test01 ( int *n, int **a )

/******************************************************************************/
/*
  Purpose:

    TEST01 allocates and assigns the array.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    25 February 2014

  Author:

    John Burkardt
*/
{
  int i;
  int seed;

  seed = 123456789;
  *n = i4_uniform_ab ( 5, 15, &seed );

  ( *a ) = ( int * ) malloc ( ( *n ) * sizeof ( int ) );

  for ( i = 0; i < *n; i++ )
  {
    ( *a )[i] = i4_uniform_ab ( 0, 100, &seed );
  }

  return;
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
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  printf ( "%s\n", time_buffer );

  return;
# undef TIME_SIZE
}
