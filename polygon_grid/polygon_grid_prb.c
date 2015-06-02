# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <string.h>

# include "polygon_grid.h"

int main ( );
void polygon_grid_count_test ( );
void polygon_grid_display_test ( );
void polygon_grid_points_test01 ( );
void polygon_grid_points_test02 ( );
void polygon_grid_points_test03 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    POLYGON_GRID_TEST tests the POLYGON_GRID library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 May 2015

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "POLYGON_GRID_TEST:\n" );
  printf ( "  C version\n" );
  printf ( "  Test the POLYGON_GRID library.\n" );

  polygon_grid_count_test ( );

  polygon_grid_points_test01 ( );
  polygon_grid_points_test02 ( );
  polygon_grid_points_test03 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "POLYGON_GRID_TEST:\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void polygon_grid_count_test ( )

/******************************************************************************/
/*
  Purpose:

    POLYGON_GRID_COUNT_TEST tests POLYGON_GRID_COUNT.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 May 2015

  Author:

    John Burkardt
*/
{
  int n;
  int ng;
  int nv;

  printf ( "\n" );
  printf ( "POLYGON_GRID_COUNT_TEST:\n" );
  printf ( "  POLYGON_GRID_COUNT counts NG, the number of points in\n" );
  printf ( "  a grid defined with N+1 points on each side of a\n" );
  printf ( "  polygon of NV vertices.\n" );

  for ( nv = 3; nv <= 5; nv++ )
  {
    printf ( "\n" );
    printf ( "  Polygonal vertex count NV = %d\n", nv );
    printf ( "\n" );
    printf ( "   N     NG\n" );
    printf ( "\n" );
    for ( n = 0; n <= 5; n++ )
    {
      ng = polygon_grid_count ( n, nv );
      printf ( "  %2d  %5d\n", n, ng );
    }
  }

  return;
}
/******************************************************************************/

void polygon_grid_points_test01 ( )

/******************************************************************************/
/*
  Purpose:

    POLYGON_GRID_POINTS_TEST01 tests POLYGON_GRID_POINTS

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 May 2015

  Author:

    John Burkardt
*/
{
  char filename[255];
  int n;
  int ng;
  int nv = 3;
  char prefix[255];
  double v[2*3] = {
    0.0, 0.0, 
    1.0, 0.0, 
    0.5, 0.86602540378443860 };
  double *xg;

  printf ( "\n" );
  printf ( "POLYGON_GRID_POINTS_TEST01:\n" );
  printf ( "  POLYGON_GRID_POINTS returns grid points for a polygon\n" );
  printf ( "  of NV vertices, with N+1 points on a side\n" );
  printf ( "\n" );
  printf ( "  For this test, the polygon is a triangle.\n" );

  r8mat_transpose_print ( 2, nv, v, "  Polygon vertices:" );
/*
  Count the grid points.
*/
  n = 5;
  ng = polygon_grid_count ( n, nv );

  printf ( "\n" );
  printf ( "  N = %d\n", n );
  printf ( "  Number of grid points will be NG = %d\n", ng );
/*
  Compute the grid points.
*/
  xg = polygon_grid_points ( n, nv, v, ng );

  r8mat_transpose_print ( 2, ng, xg, "  The grid point array:" );
/*
  Display the points.
*/
  strcpy ( prefix, "triangle" );

  polygon_grid_display ( n, nv, v, ng, xg, prefix );
/*
  Write the points to a file.
*/
  strcpy ( filename, "triangle.xy" );

  r8mat_write ( filename, 2, ng, xg );

  printf ( "\n" );
  printf ( "  Data written to the file '%s'\n", filename );

  free ( xg );

  return;
}
/******************************************************************************/

void polygon_grid_points_test02 ( )

/******************************************************************************/
/*
  Purpose:

    POLYGON_GRID_POINTS_TEST02 tests POLYGON_GRID_POINTS

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 May 2015

  Author:

    John Burkardt
*/
{
  char filename[255];
  int n;
  int ng;
  int nv = 4;
  char prefix[255];
  double v[2*4] = {
    1.0, 1.0, 
    2.0, 0.0, 
    4.0, 3.0, 
    0.0, 5.0 };
  double *xg;

  printf ( "\n" );
  printf ( "POLYGON_GRID_POINTS_TEST02:\n" );
  printf ( "  POLYGON_GRID_POINTS returns grid points for a polygon\n" );
  printf ( "  of NV vertices, with N+1 points on a side\n" );
  printf ( "\n" );
  printf ( "  For this test, the polygon is a convex quadrilateral\n" );
  printf ( "  with sides of varying length.\n" );
/*
  Define the polygon.
*/
  r8mat_transpose_print ( 2, nv, v, "  Polygon vertices:" );
/*
  Count the grid points.
*/
  n = 7;
  ng = polygon_grid_count ( n, nv );

  printf ( "\n" );
  printf ( "  N = %d\n", n );
  printf ( "  Number of grid points will be NG = %d\n", ng );
/*
  Compute the grid points.
*/
  xg = polygon_grid_points ( n, nv, v, ng );

  r8mat_transpose_print ( 2, ng, xg, "  The grid point array:" );
/*
  Display the points.
*/
  strcpy ( prefix, "quad" );

  polygon_grid_display ( n, nv, v, ng, xg, prefix );
/*
  Write the points to a file.
*/
  strcpy ( filename, "quad.xy" );

  r8mat_write ( filename, 2, ng, xg );

  printf ( "\n" );
  printf ( "  Data written to the file '%s'\n", filename );

  free ( xg );

  return;
}
/******************************************************************************/

void polygon_grid_points_test03 ( )

/******************************************************************************/
/*
  Purpose:

    POLYGON_GRID_POINTS_TEST03 tests POLYGON_GRID_POINTS

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 May 2015

  Author:

    John Burkardt
*/
{
  char filename[255];
  int n;
  int ng;
  int nv = 6;
  char prefix[255];
  double v[2*6] = {
    0.0, 0.0, 
    2.0, 0.0, 
    2.0, 1.0, 
    1.0, 1.0, 
    1.0, 2.0, 
    0.0, 2.0 };
  double *xg;

  printf ( "\n" );
  printf ( "POLYGON_GRID_POINTS_TEST03:\n" );
  printf ( "  POLYGON_GRID_POINTS returns grid points for a polygon\n" );
  printf ( "  of NV vertices, with N+1 points on a side\n" );
  printf ( "\n" );
  printf ( "  For this test, the polygon is nonconvex and six sided.\n" );
  printf ( "  Two degenerate triangles are created, and some grid points\n" );
  printf ( "  are generated several times.\n" );
/*
  Define the polygon.
*/
  r8mat_transpose_print ( 2, nv, v, "  Polygon vertices:" );
/*
  Count the grid points.
*/
  n = 5;
  ng = polygon_grid_count ( n, nv );

  printf ( "\n" );
  printf ( "  N = %d\n", n );
  printf ( "  Number of grid points will be NG = %d\n", ng );
/*
  Compute the grid points.
*/
  xg = polygon_grid_points ( n, nv, v, ng );

  r8mat_transpose_print ( 2, ng, xg, "  The grid point array:" );
/*
  Display the points.
*/
  strcpy ( prefix, "ell" );

  polygon_grid_display ( n, nv, v, ng, xg, prefix );
/*
  Write the points to a file.
*/
  strcpy ( filename, "ell.xy" );

  r8mat_write ( filename, 2, ng, xg );

  printf ( "\n" );
  printf ( "  Data written to the file '%s'\n", filename );

  free ( xg );

  return;
}

