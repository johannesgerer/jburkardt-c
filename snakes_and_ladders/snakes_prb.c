# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "snakes.h"

int main ( );
void test01 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for SNAKES_PRB.

  Discussion:

    SNAKES_PRB tests the SNAKES library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 September 2014

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "SNAKES_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the SNAKES library.\n" );

  test01 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "SNAKES_PRB\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void test01 ( )

/******************************************************************************/
/*
  Purpose:

    TEST01 tests SPY_GE for the SNAKES matrix.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 September 2014

  Author:

    John Burkardt
*/
{
  double *a;
  char header[] = "snakes";
  int m;
  int n;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  SNAKES sets up the snakes and ladders matrix.\n" );
  printf ( "  SPY_GE generates a sparsity plot for a matrix stored\n" );
  printf ( "  in general (GE) format.\n" );

  a = snakes ( );
  m = 101;
  n = 101;
  spy_ge ( m, n, a, header );

  free ( a );

  return;
}

