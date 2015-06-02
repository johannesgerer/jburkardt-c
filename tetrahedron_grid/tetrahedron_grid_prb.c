# include <stdlib.h>
# include <stdio.h>
# include <time.h>

# include "tetrahedron_grid.h"

int main ( void );
void tetrahedron_grid_test01 ( void );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for TETRAHEDRON_GRID_PRB.

  Discussion:

    TETRAHEDRON_GRID_PRB tests the TETRAHEDRON_GRID library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 November 2011

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "TETRAHEDRON_GRID_PRB:\n" );
  printf ( "  C version\n" );
  printf ( "  Test the TETRAHEDRON_GRID library.\n" );

  tetrahedron_grid_test01 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "TETRAHEDRON_GRID_PRB:\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void tetrahedron_grid_test01 ( void )

/******************************************************************************/
/*
  Purpose:

    TETRAHEDRON_GRID_TEST01 tests TETRAHEDRON_GRID.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 November 2011

  Author:

    John Burkardt
*/
{
  char *filename = "tetrahedron_grid_test01.xyz";
  int n;
  int ng;
  double t[3*4] = {
    0.0, 0.0, 0.0, 
    1.0, 0.0, 0.0, 
    0.0, 1.0, 0.0, 
    0.0, 0.0, 1.0 };
  double *tg;

  printf ( "\n" );
  printf ( "TEST01:\n" );
  printf ( "  TETRAHEDRON_GRID can define a tetrahedral grid\n" );
  printf ( "  with N+1 points on a side, based on any tetrahedron.\n" );

  n = 10;
  printf ( "  N = %d\n", n );

  ng = tetrahedron_grid_count ( n );

  r8mat_print ( 3, 4, t, "  Tetrahedron vertices:" );

  tg = tetrahedron_grid ( n, t, ng );

  r83vec_print_part ( ng, tg, 20, "  Part of the grid point array:" );

  r8mat_write ( filename, 3, ng, tg );

  printf ( "\n" );
  printf ( "  Data written to the file \"%s\".\n", filename );

  free ( tg );

  return;
}
