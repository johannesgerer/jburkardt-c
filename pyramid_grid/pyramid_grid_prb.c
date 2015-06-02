# include <stdlib.h>
# include <stdio.h>
# include <string.h>
# include <time.h>

# include "pyramid_grid.h"

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for PYRAMID_GRID_PRB.

  Discussion:

    PYRAMID_GRID_PRB tests the PYRAMID_GRID library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 August 2014

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "PYRAMID_GRID_PRB:\n" );
  printf ( "  C version\n" );
  printf ( "  Test the PYRAMID_GRID library.\n" );

  test01 ( );
  test02 ( );
  test03 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "PYRAMID_GRID_PRB:\n" );
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

    TEST01 tests PYRAMID_GRID_COUNT.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 August 2014

  Author:

    John Burkardt
*/
{
  int n;
  int ng;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  PYRAMID_GRID_COUNT determines the size of a\n" );
  printf ( "  pyramid grid with N+1 points along each edge.\n" );

  printf ( "\n" );
  printf ( "   N    Size\n" );
  printf ( "\n" );
  for ( n = 0; n <= 10; n++ )
  {
    ng = pyramid_grid_count ( n );
    printf ( "  %2d  %6d\n", n, ng );
  }

  return;
}
/******************************************************************************/

void test02 ( )

/******************************************************************************/
/*
  Purpose:

    TEST02 tests PYRAMID_UNIT_GRID.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 August 2014

  Author:

    John Burkardt
*/
{
  int n;
  int ng;
  double *pg;

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  PYRAMID_UNIT_GRID determines a unit pyramid\n" );
  printf ( "  grid with N+1 points along each edge.\n" );

  n = 4;
  r8_print ( n, "  Grid parameter N:" );

  ng = pyramid_grid_count ( n );
  r8_print ( ng, "  Grid size NG:" );

  pg = pyramid_unit_grid ( n, ng );

  r8mat_transpose_print ( 3, ng, pg, "  Pyramid grid points:" );

  free ( pg );

  return;
}
/******************************************************************************/

void test03 ( )

/******************************************************************************/
/*
  Purpose:

    TEST03 tests PYRAMID_UNIT_GRID_PLOT.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 August 2014

  Author:

    John Burkardt
*/
{
  char header[255];
  int n;
  int ng;
  double *pg;

  printf ( "\n" );
  printf ( "TEST03\n" );
  printf ( "  PYRAMID_UNIT_GRID_PLOT plots a unit pyramid\n" );
  printf ( "  grid with N+1 points along each edge.\n" );

  n = 5;
  r8_print ( n, "  Grid parameter N:" );

  ng = pyramid_grid_count ( n );
  r8_print ( ng, "  Grid size NG:" );

  pg = pyramid_unit_grid ( n, ng );

  strcpy ( header, "pyramid_unit" );
  pyramid_unit_grid_plot ( n, ng, pg, header );

  free ( pg );

  return;
}
