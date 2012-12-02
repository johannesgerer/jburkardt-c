# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "asa299.h"

int main ( void );

void test01 ( void );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for ASA299_PRB.

  Discussion:

    ASA299_PRB calls the ASA299 test routines.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    10 January 2007

  Author:

    John Burkardt
*/
{
  timestamp ( );

  printf ( "\n" );
  printf ( "ASA299_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the ASA299 library.\n" );

  test01 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "ASA299_PRB\n" );
  printf ( "  Normal end of execution.\n" );

  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void test01 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST01 tests SIMPLEX_LATTICE_POINT_NEXT.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    10 January 2007

  Author:

    John Burkardt
*/
{
# define N 4

  int i;
  int j;
  int more;
  int t = 4;
  int x[N];

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  SIMPLEX_LATTICE_POINT_NEXT generates lattice points\n" );
  printf ( "  in the simplex\n" );
  printf ( "    0 <= X\n" );
  printf ( "    sum ( X(1:N) ) <= T\n" );
  printf ( "  Here N = %d\n", N );
  printf ( "  and T =  %d\n", t );
  printf ( "\n" );
  printf ( "     Index        X(1)      X(2)      X(3)      X(4)\n" );
  printf ( "\n" );

  more = 0;

  i = 0;

  for ( ; ; )
  {
    simplex_lattice_point_next ( N, t, &more, x );

    i = i + 1;

    printf ( "  %8d  ", i );
    for ( j = 0; j < N; j++ )
    {
      printf ( "  %8d", x[j] );
    }
    printf ( "\n" );

    if ( !more )
    {
      break;
    }

  }
 
  return;
# undef N
}
