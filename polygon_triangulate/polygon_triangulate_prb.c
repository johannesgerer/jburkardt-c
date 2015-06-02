# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <string.h>

# include "polygon_triangulate.h"

int main ( );
void test01 ( );
void test02 ( char *prefix );
void test03 ( char *prefix );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for POLYGON_TRIANGULATE_PRB.

  Discussion:

    POLYGON_TRIANGULATE_PRB tests the POLYGON_TRIANGULATE library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 May 2014

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "POLYGON_TRIANGULATE_PRB\n" );
  printf ( "  C version.\n" );
  printf ( "  Test the POLYGON_TRIANGULATE library.\n" );

  test01 ( );

  test02 ( "comb" );
  test02 ( "hand" );
  test02 ( "i18" );

  test03 ( "comb" );
  test03 ( "hand" );
  test03 ( "i18" );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "POLYGON_TRIANGULATE_PRB\n" );
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

    TEST01 tests the "comb_10" polygon.

  Discussion:

    There are N-3 triangles in the triangulation.

    For the first N-2 triangles, the first edge listed is always an
    internal diagonal.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    06 May 2014
*/
{
  int n = 10;
  int *triangles;
  double x[10] = {
    8.0, 7.0, 6.0, 5.0, 4.0, 
    3.0, 2.0, 1.0, 0.0, 4.0 };
  double y[10] = {
    0.0, 10.0,  0.0, 10.0,  0.0, 
   10.0,  0.0, 10.0,  0.0, -2.0 };

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  Triangulate the comb_10 polygon.\n" );

  triangles = polygon_triangulate ( n, x, y );

  i4mat_print ( 3, n - 2, triangles, "  Triangles" );

  free ( triangles );

  return;
}
/******************************************************************************/

void test02 ( char *prefix )

/******************************************************************************/
/*
  Purpose:

    TEST02 triangulates a polygon described in a file.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    04 May 2014

  Author:

    John Burkardt.
*/
{
  int dim_num;
  char element_filename[255];
  int i;
  int n;
  char node_filename[255];
  int triangle_num;
  int *triangles;
  double *x;
  double *xy;
  double *y;
/*
  Create filenames.
*/
  strcpy ( node_filename, prefix );
  strcat ( node_filename, "_nodes.txt" );
  strcpy ( element_filename, prefix );
  strcat ( element_filename, "_elements.txt" );

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  Read polygon coordinates in \"%s\"\n", node_filename );
/*
  Read the node coordinates.
*/
  r8mat_header_read ( node_filename, &dim_num, &n );

  xy = r8mat_data_read ( node_filename, 2, n );
/*
  Get the triangulation.
*/
  x = ( double * ) malloc ( n * sizeof ( double ) );
  y = ( double * ) malloc ( n * sizeof ( double ) );
  for ( i = 0; i < n; i++ )
  {
    x[i] = xy[0+i*2];
    y[i] = xy[1+i*2];
  }
  triangles = polygon_triangulate ( n, x, y );
/*
  Write the triangulation to a file.
*/
  triangle_num = n - 2;
  i4mat_write ( element_filename, 3, triangle_num, triangles );

  printf ( "  Write triangulation to \"%s\"\n", element_filename );
/*
  Free memory.
*/
  free ( triangles );
  free ( x );
  free ( xy );
  free ( y );

  return;
}
/******************************************************************************/

void test03 ( char *prefix )

/******************************************************************************/
/*
  Purpose:

    TEST03 plots a triangulation.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    04 May 2014

  Author:

    John Burkardt.
*/
{
  char command_filename[255];
  FILE *command_unit;
  char diagonal_filename[255];
  FILE *diagonal_unit;
  int dim_num;
  char edge_filename[255];
  FILE *edge_unit;
  int i;
  int i2;
  int j;
  int j2;
  int n;
  int node;
  char node_filename[255];
  char plot_filename[255];
  int triangle_num;
  int *triangles;
  double *x;
  double *xy;
  double *y;
/*
  Create filenames.
*/
  strcpy ( node_filename, prefix );
  strcat ( node_filename, "_nodes.txt" );
  strcpy ( edge_filename, prefix );
  strcat ( edge_filename, "_edges.txt" );
  strcpy ( diagonal_filename, prefix );
  strcat ( diagonal_filename, "_diagonals.txt" );
  strcpy ( command_filename, prefix );
  strcat ( command_filename, "_commands.txt" );
  strcpy ( plot_filename, prefix );
  strcat ( plot_filename, ".png" );

  printf ( "\n" );
  printf ( "TEST03\n" );
  printf ( "  Read node coordinates in \"%s\"\n", node_filename );
/*
  Read the node coordinates.
*/
  r8mat_header_read ( node_filename, &dim_num, &n );

  xy = r8mat_data_read ( node_filename, 2, n );
/*
  Get the triangulation.
*/
  x = ( double * ) malloc ( n * sizeof ( double ) );
  y = ( double * ) malloc ( n * sizeof ( double ) );
  for ( i = 0; i < n; i++ )
  {
    x[i] = xy[0+i*2];
    y[i] = xy[1+i*2];
  }
  triangles = polygon_triangulate ( n, x, y );
/*
  Plot the edges.
*/
  edge_unit = fopen ( edge_filename, "wt" );

  for ( j = 0; j < n + 1; j++ )
  {
    j2 = ( j % n );
    fprintf ( edge_unit, "%g  %g\n", xy[0+j2*2], xy[1+j2*2] );
  }

  fclose ( edge_unit );
/*
  Plot the diagonals.
*/
  diagonal_unit = fopen ( diagonal_filename, "wt" );

  for ( j = 0; j < n - 3; j++ )
  {
    for ( i = 0; i < 2; i++ )
    {
      node = triangles[i+j*3];
      fprintf ( diagonal_unit, "%g  %g\n", xy[0+node*2], xy[1+node*2] );
    }
    fprintf ( diagonal_unit, "\n" );
  }

  fclose ( diagonal_unit );
/*
  Write the GNUPLOT command file.
*/
  command_unit = fopen ( command_filename, "wt" );

  fprintf ( command_unit, "# %s\n", command_filename );
  fprintf ( command_unit, "#\n" );
  fprintf ( command_unit, "# Usage:\n" );
  fprintf ( command_unit, "#  gnuplot < %s\n", command_filename );
  fprintf ( command_unit, "#\n" );
  fprintf ( command_unit, "set term png\n" );
  fprintf ( command_unit, "set output \"%s\"\n", plot_filename ),
  fprintf ( command_unit, "set nokey\n" );
  fprintf ( command_unit, "set size ratio 1\n" );
  fprintf ( command_unit, "set timestamp\n" );
  fprintf ( command_unit, "set xlabel \"<---X--->\"\n" );
  fprintf ( command_unit, "set ylabel \"<---Y--->\n" );
  fprintf ( command_unit, "set title \"Edges (green) and Diagonals (red)\"\n" );
  fprintf ( command_unit, "set grid\n" );
  fprintf ( command_unit, "set style data lines\n" );
  fprintf ( command_unit, 
    "plot \"%s\" using 1:2 lw 3 linecolor rgb \"green\",\\\n", edge_filename );
  fprintf ( command_unit, 
    "     \"%s\" using 1:2 lw 3 linecolor rgb \"red\",\\\n", diagonal_filename );
  fprintf ( command_unit, 
    "     \"%s\" using 1:2 with points pt 7 ps 2 lc rgb \"black\"\n", node_filename );

  fclose ( command_unit );

  printf ( "\n" );
  printf ( "  Write edges to \"%s\"\n", edge_filename );
  printf ( "  Write diagonals to \"%s\"\n", diagonal_filename );
  printf ( "  Write gnuplot commands to \"%s\"\n", command_filename );
/*
  Free memory.
*/
  free ( triangles );
  free ( x );
  free ( xy );
  free ( y );

  return;
}

