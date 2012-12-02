# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "asa159.h"

int main ( void );
void test01 ( void );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for ASA159_PRB.

  Discussion:

    ASA159_PRB tests the routines in ASA159.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 November 2010

  Author:

    John Burkardt
*/
{
  timestamp ( );

  printf ( "\n" );
  printf ( "ASA159_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the ASA159 library.\n" );

  test01 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "ASA159_PRB\n" );
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

    TEST01 tests RCONT2.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 November 2010

  Author:

    John Burkardt
*/
{
# define M 5
# define N 5

  int a[M*N];
  int c[N] = { 2, 2, 2, 2, 1 };
  int i;
  int ierror;
  int key = 0;
  int m = M;
  int n = N;
  int ntest = 10;
  int r[M] = { 3, 2, 2, 1, 1 };
  int seed;

  seed = 123456789;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  RCONT2 constructs a random matrix with\n" );
  printf ( "  given row and column sums.\n" );

  i4vec_print ( m, r, "  The rowsum vector:" );
  i4vec_print ( n, c, "  The columnsum vector:" );

  for ( i = 1; i <= ntest; i++ )
  {
    rcont2 ( m, n, r, c, &key, &seed, a, &ierror );

    if ( ierror != 0 )
    {
      printf ( "\n" );
      printf ( "  RCONT2 returned error flag IERROR = %d\n", ierror  );
      return;
    }

    i4mat_print ( m, n, a, "  The rowcolsum matrix:" );
  }

  return;
}
