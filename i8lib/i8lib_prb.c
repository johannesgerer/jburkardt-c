# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <complex.h>

# include "i8lib.h"

int main ( void );

void test015 ( void );
void test190 ( void );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for I8LIB_PRB.

  Discussion:

    I8LIB_PRB tests the I8LIB library.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    27 June 2010

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "I8LIB_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the I8LIB routines.\n" );

  test015 ( );
  test190 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "I8LIB_PRB\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void test015 ( )

/******************************************************************************/
/*
  Purpose:

    TEST015 tests I8_CHOOSE.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    27 June 2010

  Author:

    John Burkardt
*/
{
  long long int cnk;
  long long int k;
  long long int n;

  fprintf ( stdout, "\n" );
  fprintf ( stdout, "TEST015\n" );
  fprintf ( stdout, "  I8_CHOOSE evaluates C(N,K).\n" );
  fprintf ( stdout, "\n" );
  fprintf ( stdout, "     N     K    CNK\n" );
  fprintf ( stdout, "\n" );

  for ( n = 0; n <= 4; n++ )
  {
    for ( k = 0; k <= n; k++ )
    {
      cnk = i8_choose ( n, k );

      fprintf ( stdout, "%6d  %6d  %6d\n", n, k, cnk );
    }
  }

  return;
}
/******************************************************************************/

void test190 ( )

/******************************************************************************/
/*
  Purpose:

    TEST190 tests I8_XOR.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    27 June 2010

  Author:

    John Burkardt
*/
{
  long long int i;
  long long int i_lo = 0LL;
  long long int i_hi = 100LL;
  long long int j;
  long long int k;
  long long int l;
  long long int seed;
  int test;
  int test_num = 10;

  seed = 123456789LL;

  printf ( "\n" );
  printf ( "TEST190\n" );
  printf ( "  I8_XOR returns the bitwise exclusive OR of\n" );
  printf ( "  two I8's.\n" );
  printf ( "  The operator ^ should generally be used instead.\n" );
  printf ( "\n" );
  printf ( "       I       J  I8_XOR     I^J\n" );
  printf ( "\n" );

  for ( test = 1; test <= test_num; test++ )
  {
    i = i8_uniform ( i_lo, i_hi, &seed );
    j = i8_uniform ( i_lo, i_hi, &seed );
    k = i8_xor ( i, j );
    l = i ^ j;

    printf ( "  %6d  %6d  %6d  %6d\n", i, j, k, l );
  }

  return;
}
