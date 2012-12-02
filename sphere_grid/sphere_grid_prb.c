# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "sphere_grid.h"

int main ( void );

void test01 ( void );
void test02 ( void );
void test03 ( void );
void test04 ( void );
void test05 ( void );
void test06 ( void );
void test07 ( void );
void test08 ( void );
void test09 ( void );
void test10 ( void );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for SPHERE_GRID_PRB.

  Discussion:

    SPHERE_GRID_PRB tests routines from the SPHERE_GRID library.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    03 October 2012

  Author:

    John Burkardt
*/
{
  timestamp ( );

  printf ( "\n" );
  printf ( "SPHERE_GRID_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the SPHERE_GRID library.\n" );
 
  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
  test05 ( );
  test06 ( );
  test07 ( );
  test08 ( );
  test09 ( );
  test10 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "SPHERE_GRID_PRB\n" );
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

    TEST01 tests SPHERE_GRID_ICOS_NUM.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    29 August 2010

  Author:

    John Burkardt
*/
{
  int edge_num;
  int factor;
  int factor_log;
  int node_num;
  int triangle_num;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  SPHERE_ICOS_POINT_NUM determines the size\n" );
  printf ( "  (number of vertices, edges and faces) in a grid\n" );
  printf ( "  on a sphere, made by subdividing an initial\n" );
  printf ( "  projected icosahedron.\n" );
  printf ( "\n" );
  printf ( "  N determines the number of subdivisions of each\n" );
  printf ( "  edge of the icosahedral faces.\n" );
  printf ( "\n" );
  printf ( "         N         V         E         F\n" );
  printf ( "  --------  --------  --------  --------\n" );
  printf ( "\n" );

  for ( factor = 1; factor <= 20; factor++ )
  {
    node_num = sphere_icos_point_num ( factor );
    edge_num = sphere_icos_edge_num ( factor );
    triangle_num = sphere_icos_face_num ( factor );
    printf ( "  %8d  %8d  %8d  %8d\n", 
      factor, node_num, edge_num, triangle_num );
  }

  printf ( "\n" );
  printf ( "  Repeat, but using N constrained by doubling:\n" );
  printf ( "\n" );
  printf ( "         N         V         E         F\n" );
  printf ( "  --------  --------  --------  --------\n" );
  printf ( "\n" );

  factor = 1;
  for ( factor_log = 0; factor_log <= 10; factor_log++ )
  {
    node_num = sphere_icos_point_num ( factor );
    edge_num = sphere_icos_edge_num ( factor );
    triangle_num = sphere_icos_face_num ( factor );
    printf ( "  %8d  %8d  %8d  %8d\n", 
      factor, node_num, edge_num, triangle_num );
    factor = factor * 2;
  }

  return;
}
/******************************************************************************/

void test02 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST02 tests SPHERE_ICOS1_POINTS.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    28 August 2010

  Author:

    John Burkardt
*/
{
  int factor;
  char filename[80];
  int node;
  int node_num;
  double *node_xyz;
  FILE *output;

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  SPHERE_GRID_ICOS_NUM \"sizes\" a grid generated\n" );
  printf ( "  on an icosahedron and projected to a sphere.\n" );
  printf ( "  SPHERE_ICOS1_POINTS creates the grid points.\n" );

  factor = 3;

  printf ( "\n" );
  printf ( "  Sizing factor =       %d\n", factor );

  node_num = sphere_icos_point_num ( factor );

  printf ( "\n" );
  printf ( "  Number of vertices =  %d\n", node_num );

  node_xyz = sphere_icos1_points ( factor, node_num );

  r8mat_transpose_print_some ( 3, node_num, node_xyz, 1, 1, 3, 20, 
    "  Initial part of NODE_XYZ array:" );
/*
  Write the nodes to a file.
*/
  if ( 1 )
  {
    sprintf ( filename, "sphere_icos1_points_f%d.xyz", factor );

    output = fopen ( filename, "wt" );
    for ( node = 0; node < node_num; node++ )
    {
      fprintf ( output, "  %f  %f  %f\n", 
        node_xyz[0+node*3], node_xyz[1+node*3], node_xyz[2+node*3] );
    }
    fclose ( output );

    printf ( "\n" );
    printf ( "  Wrote data to \"%s\"\n", filename );
  }

  free ( node_xyz );

  return;
}
/******************************************************************************/

void test03 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST03 tests SPHERE_ICOS2_POINTS.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    29 August 2010

  Author:

    John Burkardt
*/
{
  int factor;
  char filename[80];
  int node;
  int node_num;
  double *node_xyz;
  FILE *output;

  printf ( "\n" );
  printf ( "TEST03\n" );
  printf ( "  SPHERE_GRID_ICOS_NUM \"sizes\" a grid generated\n" );
  printf ( "  on an icosahedron and projected to a sphere.\n" );
  printf ( "  SPHERE_ICOS2_POINTS creates the grid.\n" );

  factor = 3;

  printf ( "\n" );
  printf ( "  Sizing factor FACTOR = %d\n", factor );

  node_num = sphere_icos_point_num ( factor );

  printf ( "\n" );
  printf ( "  Number of nodes =     %d\n", node_num );

  node_xyz = sphere_icos2_points ( factor, node_num );

  r8mat_transpose_print_some ( 3, node_num, node_xyz, 1, 1, 3, 20, 
    "  Initial part of NODE_XYZ array:" );
/*
  Write the nodes to a file.
*/
  if ( 1 )
  {
    sprintf ( filename, "sphere_icos2_points_f%d.xyz", factor );

    output = fopen ( filename, "wt" );
    for ( node = 0; node < node_num; node++ )
    {
      fprintf ( output, "  %f  %f  %f\n", 
        node_xyz[0+node*3], node_xyz[1+node*3], node_xyz[2+node*3] );
    }
    fclose ( output );

    printf ( "\n" );
    printf ( "  Wrote data to \"%s\"\n", filename );
  }

  free ( node_xyz );

  return;
}
/******************************************************************************/

void test04 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST04 tests SPHERE_LL_POINTS.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    08 October 2012

  Author:

    John Burkardt
*/
{
  int lat_num = 3;
  int lon_num = 4;

  double pc[3] = { 0.0, 0.0, 0.0 };
  int i;
  int j;
  int k;
  int node_num;
  double *node_xyz;
  double r = 10.0;

  printf ( "\n" );
  printf ( "TEST04\n" );
  printf ( "  SPHERE_LL_POINTS produces latitude/longitude\n" );
  printf ( "  points on a sphere in 3D.\n" );

  printf ( "\n" );
  printf ( "  Radius = %f\n", r );

  r8vec_print ( 3, pc, "  Center:" );

  printf ( "\n" );
  printf ( "  The number of latitudes =  %d\n", lat_num );
  printf ( "  The number of longitudes = %d\n", lon_num );

  node_num = sphere_ll_point_num ( lat_num, lon_num );

  printf ( "\n" );
  printf ( "  The number of grid points is %d\n", node_num );

  node_xyz = sphere_ll_points ( r, pc, lat_num, lon_num, node_num );

  k = 0;
  printf ( "\n" );
  printf ( "  %8d  %12f  %12f  %12f\n", 
    k, node_xyz[0+k*3], node_xyz[1+k*3], node_xyz[2+k*3] );

  for ( i = 1; i <= lat_num; i++ )
  {
    printf ( "\n" );
    for ( j = 0; j < lon_num; j++ )
    {
      k = k + 1;
      printf ( "  %8d  %12f  %12f  %12f\n", 
        k, node_xyz[0+k*3], node_xyz[1+k*3], node_xyz[2+k*3] );
    }
  }

  printf ( "\n" );
  k = k + 1;
  printf ( "  %8d  %12f  %12f  %12f\n", 
    k, node_xyz[0+k*3], node_xyz[1+k*3], node_xyz[2+k*3] );

  free ( node_xyz );

  return;
}
/******************************************************************************/

void test05 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST05 tests SPHERE_SPIRALPOINTS.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    29 August 2010

  Author:

    John Burkardt
*/
{
  double center_xyz[3] = { 0.0, 0.0, 0.0 };
  char filename[80];
  int node;
  int node_num = 500;
  double *node_xyz;
  FILE *output;
  double r = 1.0;

  printf ( "\n" );
  printf ( "TEST05\n" );
  printf ( "  SPHERE_SPIRALPOINTS produces a spiral of\n" );
  printf ( "  points on an implicit sphere in 3D.\n" );
  printf ( "\n" );
  printf ( "  Radius = %f\n", r );

  r8vec_print ( 3, center_xyz, "  Center:" );

  printf ( "\n" );
  printf ( "  The number of spiral points is %d\n", node_num );

  node_xyz = sphere_spiralpoints ( r, center_xyz, node_num );

  r8mat_transpose_print_some ( 3, node_num, node_xyz, 1, 1, 3, 10,
    "  The spiral points:" );
/*
  Write the nodes to a file.
*/
  if ( 1 )
  {
    sprintf ( filename, "sphere_grid_spiral_n%d.xyz", node_num );

    output = fopen ( filename, "wt" );
    for ( node = 0; node < node_num; node++ )
    {
      fprintf ( output, "  %f  %f  %f\n", 
        node_xyz[0+node*3], node_xyz[1+node*3], node_xyz[2+node*3] );
    }
    fclose ( output );

    printf ( "\n" );
    printf ( "  Wrote data to \"%s\"\n", filename );
  }

  free ( node_xyz );

  return;
}
/******************************************************************************/

void test06 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST06 tests SPHERE_LL_LINES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    08 October 2012

  Author:

    John Burkardt
*/
{
  int lat_num = 3;
  int *line;
  int line_num;
  int long_num = 4;

  printf ( "\n" );
  printf ( "TEST06\n" );
  printf ( "  SPHERE_LL_LINES computes latitude/longitude\n" );
  printf ( "  lines on a sphere in 3D.\n" );
  printf ( "\n" );
  printf ( "  Number of latitudes is  %d\n", lat_num );
  printf ( "  Number of longitudes is %d\n", long_num );

  line_num = sphere_ll_line_num ( lat_num, long_num );

  printf ( "\n" );
  printf ( "  Number of line segments is %d\n", line_num );

  line = sphere_ll_lines ( lat_num, long_num, line_num );

  i4mat_transpose_print ( 2, line_num, line, "  Grid line vertices:" );

  free ( line );

  return;
}
/******************************************************************************/

void test07 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST07 tests SPHERE_GRID_Q4.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    29 August 2010

  Author:

    John Burkardt
*/
{
  int lat_num = 3;
  int long_num = 4;
  int rectangle_num = lat_num * long_num;
  int *rectangle_node;

  printf ( "\n" );
  printf ( "TEST07\n" );
  printf ( "  SPHERE_GRID_Q4 computes a grid\n" );
  printf ( "  of Q4 rectangular elements on a sphere in 3D.\n" );
  printf ( "\n" );
  printf ( "  Number of latitudes is      %d\n", lat_num );
  printf ( "  Number of longitudes is     %d\n", long_num );
  printf ( "  The number of rectangles is %d\n", rectangle_num );

  rectangle_node = sphere_grid_q4 ( lat_num, long_num );

  i4mat_transpose_print ( 4, rectangle_num, rectangle_node, 
    "  Rectangle vertices:" );

  free ( rectangle_node );

  return;
}
/******************************************************************************/

void test08 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST08 tests SPHERE_GRID_T3.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    28 August 2010

  Author:

    John Burkardt
*/
{
  int lat_num = 3;
  int lon_num = 4;

  int triangle_num;
  int *triangle_node;

  printf ( "\n" );
  printf ( "TEST08\n" );
  printf ( "  SPHERE_GRID_T3 computes a grid\n" );
  printf ( "  of T3 triangular elements on a sphere in 3D.\n" );
  printf ( "\n" );
  printf ( "  Number of latitudes is  %d\n", lat_num );
  printf ( "  Number of longitudes is %d\n", lon_num );

  triangle_node = sphere_grid_t3 ( lat_num, lon_num );

  triangle_num = 2 * ( lat_num + 1 ) * lon_num;

  printf ( "\n" );
  printf ( "  The number of triangles is %d\n", triangle_num );

  i4mat_transpose_print ( 3, triangle_num, triangle_node, 
    "  Triangle vertices:" );

  free ( triangle_node );

  return;
}
/******************************************************************************/

void test09 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST09 tests SPHERE_UNIT_SAMPLE.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    30 August 2010

  Author:

    John Burkardt
*/
{
  char filename[80];
  int node;
  int node_num;
  double *node_xyz;
  FILE *output;
  int seed = 123456789;

  printf ( "\n" );
  printf ( "TEST09\n" );
  printf ( "  For the unit sphere in 3 dimensions:\n" );
  printf ( "  SPHERE_UNIT_SAMPLE does a random sampling.\n" );

  node_num = 1000;

  node_xyz = sphere_unit_sample ( node_num, &seed );

  r8mat_transpose_print_some ( 3, node_num, node_xyz, 1, 1, 3, 10, 
    "  The spiral points:" );
/*
  Write the nodes to a file.
*/
  if ( 1 )
  {
    sprintf ( filename, "sphere_sample_n%d.xyz", node_num );

    output = fopen ( filename, "wt" );
    for ( node = 0; node < node_num; node++ )
    {
      fprintf ( output, "  %f  %f  %f\n", 
        node_xyz[0+node*3], node_xyz[1+node*3], node_xyz[2+node*3] );
    }
    fclose ( output );

    printf ( "\n" );
    printf ( "  Wrote data to \"%s\"\n", filename );
  }

  free ( node_xyz );

  return;
}
/******************************************************************************/

void test10 ( )

/******************************************************************************/
/*
  Purpose:

    TEST10 tests SPHERE_CUBED_POINTS.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    03 October 2012

  Author:

    John Burkardt
*/
{
  char filename[80];
  int j;
  int n;
  int ns;
  FILE *output;
  double *xyz;

  printf ( "\n" );
  printf ( "TEST10\n" );
  printf ( "  SPHERE_CUBED_POINTS computes points on a cubed sphere grid.\n" );

  n = 10;
  printf ( "\n" );
  printf ( "  Number of divisions on each face = %d\n", n );

  ns = sphere_cubed_point_num ( n );
  printf ( "  Total number of points = %d\n", ns );

  xyz = sphere_cubed_points ( n, ns );

  r8mat_transpose_print_some ( 3, ns, xyz, 1, 1, 3, 20, "  Initial part of XYZ array:" );
/*
  Write the nodes to a file.
*/
  if ( 1 )
  {
    sprintf ( filename, "sphere_cubed_f%d.xyz", n );

    output = fopen ( filename, "wt" );
    for ( j = 0; j < n; j++ )
    {
      fprintf ( output, "  %f  %f  %f\n", 
        xyz[0+j*3], xyz[1+j*3], xyz[2+j*3] );
    }
    fclose ( output );

    printf ( "\n" );
    printf ( "  Wrote data to \"%s\"\n", filename );
  }

  free ( xyz );

  return;
}
