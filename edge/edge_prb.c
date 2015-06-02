# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <string.h>

# include "edge.h"

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );
void test035 ( );
void test036 ( );
void test037 ( );
void test04 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    EDGE_PRB tests the EDGE library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    21 September 2014

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "EDGE_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the EDGE library.\n" );

  test01 ( );
  test02 ( );
  test03 ( );
  test035 ( );
  test036 ( );
  test037 ( );
  test04 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "EDGE_PRB\n" );
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

    TEST01 plots functions with jump discontinuities.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    23 September 2014

  Author:

    John Burkardt

  Reference:

    Rick Archibald, Anne Gelb, Jungho Yoon,
    Polynomial fitting for edge detection in irregularly sampled signals 
    and images,
    SIAM Journal on Numerical Analysis,
    Volume 43, Number 1, 2006, pages 259-279.
*/
{
  char command_filename[80];
  FILE *command_unit;
  char data_filename[80];
  FILE *data_unit;
  double *f;
  char header[80];
  int i;
  int n;
  int seed;
  int test;
  int test_num;
  char title[80];
  double *x;
  double x_max;
  double x_min;

  printf ( "\n" );
  printf ( "TEST01:\n" );
  printf ( "  Plot 1D test functions.\n" );

  test_num = 7;

  for ( test = 1; test <= test_num; test++ )
  {
    if ( test == 1 )
    {
      n = 101;
      x_min = -1.0;
      x_max = +1.0;
      x = r8vec_linspace_new ( n, x_min, x_max );
      strcpy ( header, "fx1" );
      f = fx1_vec ( n, x );
      strcpy ( title, "1D Test Function #1" );
    }
    else if ( test == 2 )
    {
      n = 101;
      x_min = -1.0;
      x_max = +1.0;
      x = r8vec_linspace_new ( n, x_min, x_max );
      strcpy ( header, "fx2" );
      f = fx2_vec ( n, x );
      strcpy ( title, "1D Test Function #2" );
    }
    else if ( test == 3 )
    {
      n = 101;
      x_min = -1.0;
      x_max = +1.0;
      x = r8vec_linspace_new ( n, x_min, x_max );
      strcpy ( header, "fx3" );
      f = fx3_vec ( n, x );
      strcpy ( title, "1D Test Function #3" );
    }
    else if ( test == 4 )
    {
      n = 101;
      x_min = 0.0;
      x_max = +1.0;
      x = r8vec_linspace_new ( n, x_min, x_max );
      strcpy ( header, "fx4" );
      f = fx4_vec ( n, x );
      strcpy ( title, "1D Test Function #4" );
    }
    else if ( test == 5 )
    {
      n = 101;
      x_min = -1.0;
      x_max = +1.0;
      x = r8vec_linspace_new ( n, x_min, x_max );
      strcpy ( header, "fx5" );
      f = fx5_vec ( n, x );
      strcpy ( title, "1D Test Function #5" );
    }
    else if ( test == 6 )
    {
      n = 101;
      x_min = 0.0;
      x_max = +1.0;
      x = r8vec_linspace_new ( n, x_min, x_max );
      strcpy ( header, "fx6" );
      f = fx6_vec ( n, x );
      strcpy ( title, "1D Test Function #6" );
    }
    else if ( test == 7 )
    {
      n = 101;
      x_min = 0.0;
      x_max = +1.0;
      x = r8vec_linspace_new ( n, x_min, x_max );
      strcpy ( header, "fx7" );
      f = fx7_vec ( n, x );
      strcpy ( title, "1D Test Function #7" );
    }

    sprintf ( data_filename, "%s_data.txt", header );
    data_unit = fopen ( data_filename, "wt" );
    for ( i = 0; i < n; i++ )
    {
      fprintf ( data_unit, "%g  %g\n", x[i], f[i] );
    }
    fclose ( data_unit );
    printf ( "  Created data file '%s'\n", data_filename );

    sprintf ( command_filename, "%s_commands.txt", header );
    command_unit = fopen ( command_filename, "wt" );
    fprintf ( command_unit, "# %s\n", command_filename );
    fprintf ( command_unit, "#\n" );
    fprintf ( command_unit, "# Usage:\n" );
    fprintf ( command_unit, "#  gnuplot < %s\n", command_filename );
    fprintf ( command_unit, "#\n" );
    fprintf ( command_unit, "set term png\n" );
    fprintf ( command_unit, "set output '%s.png'\n", header );
    fprintf ( command_unit, "set xlabel '<--- X --->'\n" );
    fprintf ( command_unit, "set ylabel '<--- Y --->'\n" );
    fprintf ( command_unit, "set title '%s'\n", title );
    fprintf ( command_unit, "set grid\n" );
    fprintf ( command_unit, "set style data lines\n" );
    fprintf ( command_unit, 
      "plot '%s' using 1:2 with points lt 3 pt 4 linecolor rgb 'blue'\n", 
      data_filename );
    fprintf ( command_unit, "quit\n" );
    fclose ( command_unit );
    printf ( "  Created command file '%s'\n", command_filename );

    free ( f );
    free ( x );
  }

  return;
}
/******************************************************************************/

void test02 ( )

/******************************************************************************/
/*
  Purpose:

    TEST02 plots a function with a jump discontinuity along a circle.

  Discussion:

    This is example 4.1 in the reference.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 February 2014

  Author:

    John Burkardt

  Reference:

    Rick Archibald, Anne Gelb, Jungho Yoon,
    Polynomial fitting for edge detection in irregularly sampled signals 
    and images,
    SIAM Journal on Numerical Analysis,
    Volume 43, Number 1, 2006, pages 259-279.
*/
{
  char command_filename[80];
  FILE *command_unit;
  char data_filename[80];
  FILE *data_unit;
  double fxy;
  char header[80];
  int i;
  int j;
  int n;
  int test;
  int test_num;
  char title[80];
  double *x;
  double x_max;
  double x_min;
  double *y;
  double y_max;
  double y_min;

  printf ( "\n" );
  printf ( "TEST02:\n" );
  printf ( "  Plot 2D test function #1 with jump along circle.\n" );

  strcpy ( header, "fxy1" );
  strcpy ( title, "2D test function #1 with discontinuity along circle" );

  n = 101;
  x_min = -1.0;
  x_max = +1.0;
  x = r8vec_linspace_new ( n, x_min, x_max );
  y_min = -1.0;
  y_max = +1.0;
  y = r8vec_linspace_new ( n, y_min, y_max );

  sprintf ( data_filename, "%s_data.txt", header );
  data_unit = fopen ( data_filename, "wt" );
  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      fxy = fxy1 ( x[i], y[j] );
      fprintf ( data_unit, "%g  %g  %g\n", x[i], y[j], fxy );
    }
    fprintf ( data_unit, "\n" );
  }
  fclose ( data_unit );
  printf ( "  Created data file '%s'\n", data_filename );

  sprintf ( command_filename, "%s_commands.txt", header );
  command_unit = fopen ( command_filename, "wt" );
  fprintf ( command_unit, "# %s\n", command_filename );
  fprintf ( command_unit, "#\n" );
  fprintf ( command_unit, "# Usage:\n" );
  fprintf ( command_unit, "#  gnuplot < %s\n", command_filename );
  fprintf ( command_unit, "#\n" );
  fprintf ( command_unit, "set term png\n" );
  fprintf ( command_unit, "set output '%s.png'\n", header );
  fprintf ( command_unit, "set view 120, 77\n" );
  fprintf ( command_unit, "set hidden3d\n" );
  fprintf ( command_unit, "set timestamp\n" );
  fprintf ( command_unit, "set xlabel '<--- X --->'\n" );
  fprintf ( command_unit, "set ylabel '<--- Y --->'\n" );
  fprintf ( command_unit, "set zlabel '<--- Z --->'\n" );
  fprintf ( command_unit, "set title '%s'\n", title );
  fprintf ( command_unit, "set grid\n" );
  fprintf ( command_unit, "set style data lines\n" );
  fprintf ( command_unit, "splot '%s' with lines\n", data_filename );
  fprintf ( command_unit, "quit\n" );
  fclose ( command_unit );
  printf ( "  Created command file '%s'\n", command_filename );

  free ( x );
  free ( y );

  return;
}
/******************************************************************************/

void test03 ( )

/******************************************************************************/
/*
  Purpose:

    TEST03 plots a function with a jump discontinuity along a circle.

  Discussion:

    This is example 4.2 in the reference.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 February 2014

  Author:

    John Burkardt

  Reference:

    Rick Archibald, Anne Gelb, Jungho Yoon,
    Polynomial fitting for edge detection in irregularly sampled signals 
    and images,
    SIAM Journal on Numerical Analysis,
    Volume 43, Number 1, 2006, pages 259-279.
*/
{
  char command_filename[80];
  FILE *command_unit;
  char data_filename[80];
  FILE *data_unit;
  double fxy;
  char header[80];
  int i;
  int j;
  int n;
  int test;
  int test_num;
  char title[80];
  double *x;
  double x_max;
  double x_min;
  double *y;
  double y_max;
  double y_min;

  printf ( "\n" );
  printf ( "TEST03:\n" );
  printf ( "  Plot 2D test function #2, the Shepp Logan phantom.\n" );

  strcpy ( header, "fxy2" );
  strcpy ( title, "2D test function #2, the Shepp Logan phantom" );

  n = 101;
  x_min = -1.0;
  x_max = +1.0;
  x = r8vec_linspace_new ( n, x_min, x_max );
  y_min = -1.0;
  y_max = +1.0;
  y = r8vec_linspace_new ( n, y_min, y_max );

  sprintf ( data_filename, "%s_data.txt", header );
  data_unit = fopen ( data_filename, "wt" );
  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      fxy = fxy2 ( x[i], y[j] );
      fprintf ( data_unit, "%g  %g  %g\n", x[i], y[j], fxy );
    }
    fprintf ( data_unit, "\n" );
  }
  fclose ( data_unit );
  printf ( "  Created data file '%s'\n", data_filename );

  sprintf ( command_filename, "%s_commands.txt", header );
  command_unit = fopen ( command_filename, "wt" );
  fprintf ( command_unit, "# %s\n", command_filename );
  fprintf ( command_unit, "#\n" );
  fprintf ( command_unit, "# Usage:\n" );
  fprintf ( command_unit, "#  gnuplot < %s\n", command_filename );
  fprintf ( command_unit, "#\n" );
  fprintf ( command_unit, "set term png\n" );
  fprintf ( command_unit, "set output '%s.png'\n", header );
  fprintf ( command_unit, "set view 30, 75\n" );
  fprintf ( command_unit, "set hidden3d\n" );
  fprintf ( command_unit, "set timestamp\n" );
  fprintf ( command_unit, "set xlabel '<--- X --->'\n" );
  fprintf ( command_unit, "set ylabel '<--- Y --->'\n" );
  fprintf ( command_unit, "set zlabel '<--- Z --->'\n" );
  fprintf ( command_unit, "set title '%s'\n", title );
  fprintf ( command_unit, "set grid\n" );
  fprintf ( command_unit, "set style data lines\n" );
  fprintf ( command_unit, "splot '%s' with lines\n", data_filename );
  fprintf ( command_unit, "quit\n" );
  fclose ( command_unit );
  printf ( "  Created command file '%s'\n", command_filename );

  free ( x );
  free ( y );

  return;
}
/******************************************************************************/

void test035 ( )

/******************************************************************************/
/*
  Purpose:

    TEST035 plots a function with a jump discontinuity along a circle.

  Discussion:

    This is example 3.2 in the reference.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 September 2014

  Author:

    John Burkardt

  Reference:

    Rick Archibald, Anne Gelb, Jungho Yoon,
    Determining the locations and discontinuities in the derivatives
    of functions,
    Applied Numerical Mathematics,
    Volume 58, 2008, pages 577-592.
*/
{
  char command_filename[80];
  FILE *command_unit;
  char data_filename[80];
  FILE *data_unit;
  double fxy;
  char header[80];
  int i;
  int j;
  int n;
  int test;
  int test_num;
  char title[80];
  double *x;
  double x_max;
  double x_min;
  double *y;
  double y_max;
  double y_min;

  printf ( "\n" );
  printf ( "TEST035:\n" );
  printf ( "  Plot 2D test function #3, the modified 2D Harten function.\n" );

  strcpy ( header, "fxy3" );
  strcpy ( title, "2D test function #3, the modified 2D Harten function" );

  n = 101;
  x_min = -1.0;
  x_max = +1.0;
  x = r8vec_linspace_new ( n, x_min, x_max );
  y_min = -1.0;
  y_max = +1.0;
  y = r8vec_linspace_new ( n, y_min, y_max );

  sprintf ( data_filename, "%s_data.txt", header );
  data_unit = fopen ( data_filename, "wt" );
  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      fxy = fxy3 ( x[i], y[j] );
      fprintf ( data_unit, "%g  %g  %g\n", x[i], y[j], fxy );
    }
    fprintf ( data_unit, "\n" );
  }
  fclose ( data_unit );
  printf ( "  Created data file '%s'\n", data_filename );

  sprintf ( command_filename, "%s_commands.txt", header );
  command_unit = fopen ( command_filename, "wt" );
  fprintf ( command_unit, "# %s\n", command_filename );
  fprintf ( command_unit, "#\n" );
  fprintf ( command_unit, "# Usage:\n" );
  fprintf ( command_unit, "#  gnuplot < %s\n", command_filename );
  fprintf ( command_unit, "#\n" );
  fprintf ( command_unit, "set term png\n" );
  fprintf ( command_unit, "set output '%s.png'\n", header );
  fprintf ( command_unit, "set view 30, 75\n" );
  fprintf ( command_unit, "set hidden3d\n" );
  fprintf ( command_unit, "set timestamp\n" );
  fprintf ( command_unit, "set xlabel '<--- X --->'\n" );
  fprintf ( command_unit, "set ylabel '<--- Y --->'\n" );
  fprintf ( command_unit, "set zlabel '<--- Z --->'\n" );
  fprintf ( command_unit, "set title '%s'\n", title );
  fprintf ( command_unit, "set grid\n" );
  fprintf ( command_unit, "set style data lines\n" );
  fprintf ( command_unit, "splot '%s' with lines\n", data_filename );
  fprintf ( command_unit, "quit\n" );
  fclose ( command_unit );
  printf ( "  Created command file '%s'\n", command_filename );

  free ( x );
  free ( y );

  return;
}
/******************************************************************************/

void test036 ( )

/******************************************************************************/
/*
  Purpose:

    TEST036 plots a function with a derivative discontinuity.

  Discussion:

    This is example 3.1 in the reference.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    21 September 2014

  Author:

    John Burkardt

  Reference:

    Rick Archibald, Anne Gelb, Jungho Yoon,
    Determining the locations and discontinuities in the derivatives
    of functions,
    Applied Numerical Mathematics,
    Volume 58, 2008, pages 577-592.
*/
{
  char command_filename[80];
  FILE *command_unit;
  char data_filename[80];
  FILE *data_unit;
  double fxy;
  char header[80];
  int i;
  int j;
  int n;
  int test;
  int test_num;
  char title[80];
  double *x;
  double x_max;
  double x_min;
  double *y;
  double y_max;
  double y_min;

  printf ( "\n" );
  printf ( "TEST036:\n" );
  printf ( "  Plot 2D test function #4, the discontinuous medium wave, P(x,t).\n" );

  strcpy ( header, "fxy4" );
  strcpy ( title, "2D test function #4, the discontinuous medium wave, P(x,t)" );

  n = 101;
  x_min = -1.0;
  x_max = 0.0;
  x = r8vec_linspace_new ( n, x_min, x_max );
  y_min = 0.0;
  y_max = 0.1;
  y = r8vec_linspace_new ( n, y_min, y_max );

  sprintf ( data_filename, "%s_data.txt", header );
  data_unit = fopen ( data_filename, "wt" );
  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      fxy = fxy4 ( x[i], y[j] );
      fprintf ( data_unit, "%g  %g  %g\n", x[i], y[j], fxy );
    }
    fprintf ( data_unit, "\n" );
  }
  fclose ( data_unit );
  printf ( "  Created data file '%s'\n", data_filename );

  sprintf ( command_filename, "%s_commands.txt", header );
  command_unit = fopen ( command_filename, "wt" );
  fprintf ( command_unit, "# %s\n", command_filename );
  fprintf ( command_unit, "#\n" );
  fprintf ( command_unit, "# Usage:\n" );
  fprintf ( command_unit, "#  gnuplot < %s\n", command_filename );
  fprintf ( command_unit, "#\n" );
  fprintf ( command_unit, "set term png\n" );
  fprintf ( command_unit, "set output '%s.png'\n", header );
  fprintf ( command_unit, "set view 30, 45\n" );
  fprintf ( command_unit, "set hidden3d\n" );
  fprintf ( command_unit, "set timestamp\n" );
  fprintf ( command_unit, "set xlabel '<--- X --->'\n" );
  fprintf ( command_unit, "set ylabel '<--- Y --->'\n" );
  fprintf ( command_unit, "set zlabel '<--- Z --->'\n" );
  fprintf ( command_unit, "set title '%s'\n", title );
  fprintf ( command_unit, "set grid\n" );
  fprintf ( command_unit, "set style data lines\n" );
  fprintf ( command_unit, "splot '%s' with lines\n", data_filename );
  fprintf ( command_unit, "quit\n" );
  fclose ( command_unit );
  printf ( "  Created command file '%s'\n", command_filename );

  free ( x );
  free ( y );

  return;
}
/******************************************************************************/

void test037 ( )

/******************************************************************************/
/*
  Purpose:

    TEST037 plots a function with a derivative discontinuity.

  Discussion:

    This is example 3.1 in the reference.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    21 September 2014

  Author:

    John Burkardt

  Reference:

    Rick Archibald, Anne Gelb, Jungho Yoon,
    Determining the locations and discontinuities in the derivatives
    of functions,
    Applied Numerical Mathematics,
    Volume 58, 2008, pages 577-592.
*/
{
  char command_filename[80];
  FILE *command_unit;
  char data_filename[80];
  FILE *data_unit;
  double fxy;
  char header[80];
  int i;
  int j;
  int n;
  int test;
  int test_num;
  char title[80];
  double *x;
  double x_max;
  double x_min;
  double *y;
  double y_max;
  double y_min;

  printf ( "\n" );
  printf ( "TEST037:\n" );
  printf ( "  Plot 2D test function #5, the discontinuous medium wave, U(x,t).\n" );

  strcpy ( header, "fxy5" );
  strcpy ( title, "2D test function #5, the discontinuous medium wave, U(x,t)" );

  n = 101;
  x_min = -1.0;
  x_max = 0.0;
  x = r8vec_linspace_new ( n, x_min, x_max );
  y_min = 0.0;
  y_max = 0.1;
  y = r8vec_linspace_new ( n, y_min, y_max );

  sprintf ( data_filename, "%s_data.txt", header );
  data_unit = fopen ( data_filename, "wt" );
  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      fxy = fxy5 ( x[i], y[j] );
      fprintf ( data_unit, "%g  %g  %g\n", x[i], y[j], fxy );
    }
    fprintf ( data_unit, "\n" );
  }
  fclose ( data_unit );
  printf ( "  Created data file '%s'\n", data_filename );

  sprintf ( command_filename, "%s_commands.txt", header );
  command_unit = fopen ( command_filename, "wt" );
  fprintf ( command_unit, "# %s\n", command_filename );
  fprintf ( command_unit, "#\n" );
  fprintf ( command_unit, "# Usage:\n" );
  fprintf ( command_unit, "#  gnuplot < %s\n", command_filename );
  fprintf ( command_unit, "#\n" );
  fprintf ( command_unit, "set term png\n" );
  fprintf ( command_unit, "set output '%s.png'\n", header );
  fprintf ( command_unit, "set view 30, 45\n" );
  fprintf ( command_unit, "set hidden3d\n" );
  fprintf ( command_unit, "set timestamp\n" );
  fprintf ( command_unit, "set xlabel '<--- X --->'\n" );
  fprintf ( command_unit, "set ylabel '<--- Y --->'\n" );
  fprintf ( command_unit, "set zlabel '<--- Z --->'\n" );
  fprintf ( command_unit, "set title '%s'\n", title );
  fprintf ( command_unit, "set grid\n" );
  fprintf ( command_unit, "set style data lines\n" );
  fprintf ( command_unit, "splot '%s' with lines\n", data_filename );
  fprintf ( command_unit, "quit\n" );
  fclose ( command_unit );
  printf ( "  Created command file '%s'\n", command_filename );

  free ( x );
  free ( y );

  return;
}
/******************************************************************************/

void test04 ( )

/******************************************************************************/
/*
  Purpose:

    TEST04 plots slices of a 3D function.

  Discussion:

    Although the slice plots look uninteresting, there is a lot of detail
    hidden in the data in variations that are not obvious at first.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 February 2014

  Author:

    John Burkardt

  Reference:

    Larry Shepp,
    Computerized tomography and nuclear magnetic resonance,
    Journal of Computer Assisted Tomography,
    Volume 4, Number 1, February 1980, pages 94-107.
*/
{
  char command_filename[80];
  FILE *command_unit;
  char data_filename[80];
  FILE *data_unit;
  double fxyz;
  char header[80];
  int i;
  int j;
  int k;
  int n;
  int test;
  int test_num;
  char title[80];
  double *x;
  double x_max;
  double x_min;
  double x_val;
  double *y;
  double y_max;
  double y_min;
  double y_val;
  double *z;
  double z_max;
  double z_min;
  double z_val;

  printf ( "\n" );
  printf ( "TEST04:\n" );
  printf ( "  Plot 3D test function #1, the Shepp Logan 3D phantom.\n" );

  test_num = 3;

  n = 101;
  x_min = -1.5;
  x_max = +1.5;
  x = r8vec_linspace_new ( n, x_min, x_max );
  y_min = -1.5;
  y_max = +1.5;
  y = r8vec_linspace_new ( n, y_min, y_max );
  z_min = -1.5;
  z_max = +1.5;
  z = r8vec_linspace_new ( n, z_min, z_max );

  for ( test = 1; test <= test_num; test++ )
  {
    if ( test == 1 )
    {
      x_val = 0.0;
      strcpy ( title, "Slice X = 0.0" );
      strcpy ( header, "fxyz1_x" );
    }
    else if ( test == 2 )
    {
      y_val = 0.0;
      strcpy ( title, "Slice Y = 0.0" );
      strcpy ( header, "fxyz1_y" );
    }
    else if ( test == 3 )
    {
      z_val = - 0.1;
      strcpy ( title, "Slice Z = - 0.1" );
      strcpy ( header, "fxyz1_z" );
    }

    sprintf ( data_filename, "%s_data.txt", header );
    data_unit = fopen ( data_filename, "wt" );
    for ( i = 0; i < n; i++ )
    {
      for ( j = 0; j < n; j++ )
      {
        if ( test == 1 )
        {
          fxyz = fxyz1 ( x_val, y[j], z[i] );
          fprintf ( data_unit, "%g  %g  %g\n", y[j], z[i], fxyz );
        }
        else if ( test == 2 )
        {
          fxyz = fxyz1 ( x[j], y_val, z[i] );
          fprintf ( data_unit, "%g  %g  %g\n", x[j], z[i], fxyz );
        }
        else if ( test == 3 )
        {
          fxyz = fxyz1 ( x[j], y[i], z_val );
          fprintf ( data_unit, "%g  %g  %g\n", x[j], y[i], fxyz );
        }
      }
      fprintf ( data_unit, "\n" );
    }
    fclose ( data_unit );
    printf ( "  Created data file '%s'\n", data_filename );

    sprintf ( command_filename, "%s_commands.txt", header );
    command_unit = fopen ( command_filename, "wt" );
    fprintf ( command_unit, "# %s\n", command_filename );
    fprintf ( command_unit, "#\n" );
    fprintf ( command_unit, "# Usage:\n" );
    fprintf ( command_unit, "#  gnuplot < %s\n", command_filename );
    fprintf ( command_unit, "#\n" );
    fprintf ( command_unit, "set term png\n" );
    fprintf ( command_unit, "set output '%s.png'\n", header );
    fprintf ( command_unit, "set view 20, 75\n" );
    fprintf ( command_unit, "set hidden3d\n" );
    fprintf ( command_unit, "set timestamp\n" );
    if ( test == 1 )
    {
      fprintf ( command_unit, "set xlabel '<--- Y --->'\n" );
      fprintf ( command_unit, "set ylabel '<--- Z --->'\n" );
      fprintf ( command_unit, "set zlabel '<--- X --->'\n" );
    }
    else if ( test == 2 )
    {
      fprintf ( command_unit, "set xlabel '<--- X --->'\n" );
      fprintf ( command_unit, "set ylabel '<--- Z --->'\n" );
      fprintf ( command_unit, "set zlabel '<--- Y --->'\n" );
    }
    else if ( test == 3 )
    {
      fprintf ( command_unit, "set xlabel '<--- X --->'\n" );
      fprintf ( command_unit, "set ylabel '<--- Y --->'\n" );
      fprintf ( command_unit, "set zlabel '<--- Z --->'\n" );
    }
    fprintf ( command_unit, "set title '%s'\n", title );
    fprintf ( command_unit, "set grid\n" );
    fprintf ( command_unit, "set style data lines\n" );
    fprintf ( command_unit, "splot '%s' with lines\n", data_filename );
    fprintf ( command_unit, "quit\n" );
    fclose ( command_unit );
    printf ( "  Created command file '%s'\n", command_filename );
  }

  free ( x );
  free ( y );
  free ( z );

  return;
}
