# include <math.h>
# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <time.h>

# include "sphere_llq_grid.h"

int main ( );
void sphere_llq_grid_point_count_test ( );
void sphere_llq_grid_points_test ( );
void sphere_llq_grid_line_count_test ( );
void sphere_llq_grid_lines_test ( );
void sphere_llq_grid_display_test ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for SPHERE_LLQ_GRID_PRB.

  Discussion:

    SPHERE_LLQ_GRID_PRB tests the SPHERE_LLQ_GRID library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 May 2015

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "SPHERE_LLQ_GRID_TEST\n" );
  printf ( "  C version\n" );
  printf ( "  Test the SPHERE_LLQ_GRID library.\n" );

  sphere_llq_grid_point_count_test ( );
  sphere_llq_grid_points_test ( );
  sphere_llq_grid_line_count_test ( );
  sphere_llq_grid_lines_test ( );
  sphere_llq_grid_display_test ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "SPHERE_LLQ_GRID_TEST\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void sphere_llq_grid_point_count_test ( )

/******************************************************************************/
/*
  Purpose:

    SPHERE_LLQ_GRID_POINT_COUNT_TEST tests SPHERE_LLQ_GRID_POINT_COUNT.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 May 2015

  Author:

    John Burkardt
*/
{
  int lat_num;
  int long_log;
  int long_num;
  int point_num;

  printf ( "\n" );
  printf ( "SPHERE_LLQ_GRID_POINT_COUNT_TEST\n" );
  printf ( "  SPHERE_LLQ_GRID_POINT_COUNT counts the points used for a\n" );
  printf ( "  grid based on quadrilaterals defined by latitude and longitude\n" );
  printf ( "  lines on a sphere in 3D.\n" );
  printf ( "\n" );
  printf ( "     LAT_NUM    LONG_NUM   POINT_NUM\n" );

  for ( lat_num = 1; lat_num <= 17; lat_num = lat_num + 2 );
  {
    printf ( "\n" );
    long_num = 1;
    for ( long_log = 1; long_log <= 4; long_log++ )
    {
      long_num = long_num * 2;
      point_num = sphere_llq_grid_point_count ( lat_num, long_num );
      printf ( "  %8d  %8d  %8d\n", lat_num, long_num, point_num );
    }
  }

  return;
}
/******************************************************************************/

void sphere_llq_grid_points_test ( )

/******************************************************************************/
/*
  Purpose:

    SPHERE_LLQ_GRID_POINTS_TEST tests SPHERE_LLQ_GRID_POINTS.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 May 2015

  Author:

    John Burkardt
*/
{
  int i;
  int j;
  int k;
  int lat_num;
  int long_num;
  int node_num;
  double *node_xyz;
  double pc[3] = { 0.0, 0.0, 0.0 };
  double r;

  lat_num = 3;
  long_num = 4;

  r = 10.0;

  printf ( "\n" );
  printf ( "SPHERE_LLQ_GRID_POINTS_TEST\n" );
  printf ( "  SPHERE_LLQ_POINTS produces latitude/longitude\n" );
  printf ( "  points on a sphere in 3D.\n" );

  printf ( "\n" );
  printf ( "  Radius = %g\n", r );

  r8vec_print ( 3, pc, "  Center:" );

  printf ( "\n" );
  printf ( "  Number of latitudes is  %d\n", lat_num );
  printf ( "  Number of longitudes is %d\n", long_num );

  node_num = sphere_llq_grid_point_count ( lat_num, long_num );

  printf ( "\n" );
  printf ( "  The number of grid points is %d\n", node_num );

  node_xyz = sphere_llq_grid_points ( r, pc, lat_num, long_num, node_num );

  printf ( "\n" );

  k = 0;
  printf ( "  %8d  %14.6g  %14.6g  %14.6g\n", 
    k, node_xyz[0+k*3], node_xyz[1+k*3], node_xyz[2+k*3] );
  k = k + 1;

  printf ( "\n" );

  for ( i = 1; i <= lat_num; i++ )
  {
    printf ( "\n" );
    for ( j = 0; j < long_num; j++ )
    {
      printf ( "  %8d  %14.6g  %14.6g  %14.6g\n", 
        k, node_xyz[0+k*3], node_xyz[1+k*3], node_xyz[2+k*3] );
      k = k + 1;
      printf ( "\n" );
    }
  }

  printf ( "\n" );

  printf ( "  %8d  %14.6g  %14.6g  %14.6g\n", 
    k, node_xyz[0+k*3], node_xyz[1+k*3], node_xyz[2+k*3] );
  k = k + 1;
  printf ( "\n" );

  free ( node_xyz );

  return;
}
/******************************************************************************/

void sphere_llq_grid_line_count_test ( )

/******************************************************************************/
/*
  Purpose:

    SPHERE_LLQ_GRID_LINE_COUNT_TEST tests SPHERE_LLQ_GRID_LINE_COUNT.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 May 2015

  Author:

    John Burkardt
*/
{
  int lat_num;
  int line_num;
  int long_log;
  int long_num;

  lat_num = 3;
  long_num = 4;

  printf ( "\n" );
  printf ( "SPHERE_LLQ_GRID_LINE_COUNT_TEST\n" );
  printf ( "  SPHERE_LLQ_GRID_LINE_COUNT counts the lines used for a\n" );
  printf ( "  grid based on quadrilaterals defined by latitude and longitude\n" );
  printf ( "  lines on a sphere in 3D.\n" );
  printf ( "\n" );
  printf ( "     LAT_NUM    LONG_NUM   LINE_NUM\n" );

  for ( lat_num = 1; lat_num <= 17; lat_num = lat_num + 2 )
  {
    printf ( "\n" );
    long_num = 1;
    for ( long_log = 1; long_log <= 4; long_log++ )
    {
      long_num = long_num * 2;
      line_num = sphere_llq_grid_line_count ( lat_num, long_num );
      printf ( "  %8d  %8d  %8d\n", lat_num, long_num, line_num );
    }
  }

  return;
}
/******************************************************************************/

void sphere_llq_grid_lines_test ( )

/******************************************************************************/
/*
  Purpose:

    SPHERE_LLQ_GRID_LINES_TEST tests SPHERE_LLQ_GRID_LINES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 May 2015

  Author:

    John Burkardt
*/
{
  int lat_num;
  int *line_data;
  int line_num;
  int long_num;

  lat_num = 3;
  long_num = 4;

  printf ( "\n" );
  printf ( "SPHERE_LLQ_GRID_LINES_TEST\n" );
  printf ( "  SPHERE_LLQ_GRID_LINES computes grid lines\n" );
  printf ( "  on a sphere in 3D.\n" );
  printf ( "\n" );
  printf ( "  Number of latitudes is  %d\n", lat_num );
  printf ( "  Number of longitudes is %d\n", long_num );

  line_num = sphere_llq_grid_line_count ( lat_num, long_num );

  printf ( "\n" );
  printf ( "  Number of line segments is %d\n", line_num );

  line_data = sphere_llq_grid_lines ( lat_num, long_num, line_num );

  i4mat_transpose_print ( 2, line_num, line_data, 
    "  Grid line vertices:" );

  free ( line_data );

  return;
}
/******************************************************************************/

void sphere_llq_grid_display_test ( )

/******************************************************************************/
/*
  Purpose:

    SPHERE_LLQ_GRID_DISPLAY_TEST tests SPHERE_LLQ_GRID_DISPLAY.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 May 2015

  Author:

    John Burkardt
*/
{
  int lat_num;
  int *line_data;
  int line_num;
  int long_num;
  int node_num;
  double *node_xyz;
  double pc[3] = { 0.0, 0.0, 0.0 };
  char prefix[255];
  double r;

  lat_num = 10;
  long_num = 12;

  r = 10.0;

  printf ( "\n" );
  printf ( "SPHERE_LLQ_GRID_DISPLAY_TEST\n" );
  printf ( "  SPHERE_LLQ_GRID_DISPLAY displays an LLQ grid on a sphere.\n" );
  printf ( "\n" );
  printf ( "  Number of latitudes is  %d\n", lat_num );
  printf ( "  Number of longitudes is %d\n", long_num );
/*
  Get points.
*/
  node_num = sphere_llq_grid_point_count ( lat_num, long_num );

  printf ( "\n" );
  printf ( "  The number of grid points is %d\n", node_num );

  node_xyz = sphere_llq_grid_points ( r, pc, lat_num, long_num, node_num );
/*
  Get lines.
*/
  line_num = sphere_llq_grid_line_count ( lat_num, long_num );

  printf ( "\n" );
  printf ( "  Number of line segments is %d\n", line_num );

  line_data = sphere_llq_grid_lines ( lat_num, long_num, line_num );

  strcpy ( prefix, "sphere_llq_grid" );

  sphere_llq_grid_display ( node_num, node_xyz, line_num, line_data, prefix );

  free ( line_data );
  free ( node_xyz );

  return;
}
