# include <stdio.h>
# include <stdlib.h>

# include "bellman_ford.h"

int main ( );
void test01 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for BELLMAN_FORD_PRB.

  Discussion:

    BELLMAN_FORD_PRB tests the BELLMAN_FORD library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 November 2014

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "BELLMAN_FORD_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the BELLMAN_FORD library.\n" );

  test01 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "BELLMAN_FORD_PRB\n" );
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

    TEST01 runs a simple test.

  Discussion:

    The correct distances are { 0, -6, -2, 3, 0, 0 }.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 November 2014

  Author:

    John Burkardt
*/
{
  int e[2*10] = 
  {
    1, 0,
    4, 1,
    1, 2,
    2, 4,
    4, 0,
    2, 5,
    5, 0,
    3, 2,
    5, 3,
    3, 0,
  };
  int e_num = 10;
  double e_weight[10] = 
  {
    -3.0,
     6.0,
    -4.0,
    -1.0,
     4.0,
    -2.0,
     2.0,
     8.0,
    -3.0,
     3.0
  };
  int predecessor[6];
  int source = 0;
  int v_num = 6;
  double v_weight[6];
  
  printf ( "\n" );
  printf ( "TEST01:\n" );
  printf ( "  Bellman-Ford shortest path algorithm.\n" );

  printf ( "\n" );
  printf ( "  Number of vertices = %d\n", v_num );
  printf ( "  Number of edges = %d\n", e_num );
  printf ( "  The reference vertex is %d\n", source );
  i4mat_transpose_print ( 2, e_num, e, "  The edge array:" );
  r8vec_print ( e_num, e_weight, "  The edge weights:" );

  bellman_ford ( v_num, e_num, source, e, e_weight, v_weight, predecessor );

  r8vec_print ( v_num, v_weight, "  The shortest distances:" );

  i4vec_print ( v_num, predecessor, "  The vertex predecessor parents for the shortest paths:" );

  return;
}