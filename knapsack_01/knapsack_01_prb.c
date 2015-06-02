# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "knapsack_01.h"

int main ( );
void test01 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    KNAPSACK_01_TEST tests the KNAPSACK_01 library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    22 August 2014

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "KNAPSACK_01_TEST\n" );
  printf ( "  C version.\n" );
  printf ( "  Test the KNAPSACK_01 library.\n" );

  test01 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "KNAPSACK_01_TEST\n" );
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

    TEST01 seeks a solution of the 0/1 Knapsack problem.

  Discussion:

    In the 0/1 knapsack problem, a knapsack of capacity C is given,
    as well as N items, with the I-th item of weight W(I).

    A selection is "acceptable" if the total weight is no greater than C.

    It is desired to find an optimal acceptable selection, that is,
    an acceptable selection such that there is no acceptable selection
    of greater weight.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    22 August 2014

  Author:

    John Burkardt
*/
{

  int c;
  int i;
  int n = 6;
  int *s;
  int t;
  int w[6] = {
    16, 17, 23, 24, 39, 40 };

  c = 100;

  printf ( "\n" );
  printf ( "TEST01:\n" );
  printf ( "  Knapsack maximum capacity is %d\n", c );
  printf ( "  Come as close as possible to filling the knapsack.\n" );

  s = knapsack_01 ( n, w, c );

  printf ( "\n" );
  printf ( "   # 0/1  Weight\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %2d  %1d  %4d\n", i, s[i], w[i] );
  }
  t = 0;
  for ( i = 0; i < n; i++ )
  {
    t = t + s[i] * w[i];
  }
  printf ( "\n" );
  printf ( "  Total:   %d\n", t );

  free ( s );

  return;
}
