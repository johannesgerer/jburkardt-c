# include <stdlib.h>
# include <stdio.h>

void f ( int n );
int main ( void );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for TEST01.

  Discussion:

    TEST01 calls F, which has a memory "leak".  This memory leak can be
    detected by VALGRID.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    18 May 2011
*/
{
  int n = 10;

  printf ( "\n" );
  printf ( "TEST01:\n" );
  printf ( "  C version\n" );
  printf ( "  A sample code for analysis by VALGRIND.\n" );

  f ( n );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  Normal end of execution.\n" );

  return 0;
}
/******************************************************************************/

void f ( int n )

/******************************************************************************/
/*
  Purpose:

    F computes N+1 entries of the Fibonacci sequence.

  Discussion:

    Unfortunately, F only allocates space for N entries.  Hence, the
    assignment of a value to the N+1 entry causes a memory leak.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    18 May 2011
*/
{
  int i;
  int *x;

  x = ( int * ) malloc ( n * sizeof ( int ) );

  x[0] = 1;
  printf ( "  %2d  %2d\n", 0, x[0] );

  x[1] = 1;
  printf ( "  %2d  %2d\n", 1, x[1] );

  for ( i = 2; i <= n; i++ )
  {
    x[i] = x[i-1] + x[i-2];
    printf ( "  %2d  %2d\n", i, x[i] );
  }

  free ( x );

  return;
}
