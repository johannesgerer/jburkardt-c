# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "hypersphere_properties.h"

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );
void test04 ( );
void test05 ( );
void test06 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for HYPERSPHERE_PROPERTIES_PRB.

  Discussion:

    HYPERSPHERE_PROPERTIES_PRB tests the HYPERSPHERE_PROPERTIES library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 December 2013

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "HYPERSPHERE_PROPERTIES_PRB:\n" );
  printf ( "  C version\n" );
  printf ( "  Test the HYPERSPHERE_PROPERTIES library.\n" );

  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
  test05 ( );
  test06 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "HYPERSPHERE_PROPERTIES_PRB:\n" );
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

    TEST01 tests the coordinate conversion routines.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 December 2013

  Author:

    John Burkardt
*/
{
  double *c;
  double err;
  int m;
  int n;
  double *r;
  int seed;
  int test;
  double *theta;
  double *x;
  double *x2;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  Test the coordinate conversion routines:\n" );
  printf ( "  CARTESIAN_TO_HYPERSPHERE: X       -> R,Theta\n" );
  printf ( "  HYPERSPHERE_TO_CARTESIAN: R,Theta -> X.\n" );
  printf ( "\n" );
  printf ( "  Pick a random X, and compute X2 by converting X\n" );
  printf ( "  to hypersphere and back.  Consider norm of difference.\n" );
  printf ( "\n" );
  printf ( "  M    || X - X2 ||\n" );

  seed = 123456789;

  n = 1;
  r = ( double * ) malloc ( n * sizeof ( double ) );

  for ( m = 1; m <= 5; m++ )
  {
    printf ( "\n" );

    theta = ( double * ) malloc ( ( m - 1 ) * n * sizeof ( double ) );

    for ( test = 1; test <= 5; test++ )
    {
      x = r8mat_uniform_01_new ( m, n, &seed );
      c = r8vec_uniform_01_new ( m, &seed );
      cartesian_to_hypersphere ( m, n, c, x, r, theta );
      x2 = hypersphere_to_cartesian ( m, n, c, r, theta );
      err = r8mat_norm_fro_affine ( m, n, x, x2 );
      printf ( "  %2d  %14.6g\n", m, err );
      free ( c );
      free ( x );
      free ( x2 );
    }
    free ( theta );
  }

  free ( r );

  return;
}
/******************************************************************************/

void test02 ( )

/******************************************************************************/
/*
  Purpose:

    TEST02 tests HYPERSPHERE_01_SURFACE_UNIFORM.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 December 2013

  Author:

    John Burkardt
*/
{
  int m;
  int n;
  int seed;
  int test;
  double *x;

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  HYPERSPHERE_01_SURFACE_UNIFORM samples uniformly from the\n" );
  printf ( "  surface of the unit hypersphere\n" );

  seed = 123456789;

  n = 1;
  for ( m = 1; m <= 5; m++ )
  {
    for ( test = 1; test <= 3; test++ )
    {
      x = hypersphere_01_surface_uniform ( m, n, &seed );
      r8vec_transpose_print ( m, x, "  Random hypersphere point:" );
      free ( x );
    }
  }
  return;
}
/******************************************************************************/

void test03 ( )

/******************************************************************************/
/*
  Purpose:

    TEST03 tests HYPERSPHERE_01_AREA.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 December 2013

  Author:

    John Burkardt
*/
{
  double area;
  double area2;
  int m;
  int n_data;

  printf ( "\n" );
  printf ( "TEST03:\n" );
  printf ( "  HYPERSPHERE_01_AREA evaluates the area of the unit\n" );
  printf ( "  hypersphere in M dimensions.\n" );
  printf ( "\n" );
  printf ( "       M      Exact       Computed\n" );
  printf ( "              Area        Area\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    hypersphere_01_area_values ( &n_data, &m, &area );

    if ( n_data == 0 )
    {
      break;
    }

    area2 = hypersphere_01_area ( m );

    printf ( "  %6d  %10.4f  %10.4f\n", m, area, area2 );
  }

  return;
}
/******************************************************************************/

void test04 ( )

/******************************************************************************/
/*
  Purpose:

    TEST04 tests HYPERSPHERE_01_VOLUME.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 December 2013

  Author:

    John Burkardt
*/
{
  int m;
  int n_data;
  double volume;
  double volume2;

  printf ( "\n" );
  printf ( "TEST04:\n" );
  printf ( "  HYPERSPHERE_01_VOLUME evaluates the area of the unit\n" );
  printf ( "  hypersphere in M dimensions.\n" );
  printf ( "  HYPERSPHERE_01_VOLUME_VALUES returns some test values.\n" );
  printf ( "\n" );
  printf ( "       M      Exact       Computed\n" );
  printf ( "              Volume      Volume\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    hypersphere_01_volume_values ( &n_data, &m, &volume );

    if ( n_data == 0 )
    {
      break;
    }

    volume2 = hypersphere_01_volume ( m );

    printf ( "  %6d  %10.4f  %10.4f\n", m, volume, volume2 );
  }
  return;
}
/******************************************************************************/

void test05 ( )

/******************************************************************************/
/*
  Purpose:

    TEST05 tests HYPERSPHERE_AREA, HYPERSPHERE_VOLUME.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 May 2013

  Author:

    John Burkardt
*/
{
  double area;
  int m;
  double r;
  double volume;

  r = 1.5;

  printf ( "\n" );
  printf ( "TEST05\n" );
  printf ( "  For a hypersphere in M dimensions:\n" );
  printf ( "  HYPERSPHERE_AREA computes the area\n" );
  printf ( "  HYPERSPHERE_VOLUME computes the volume.\n" );
  printf ( "\n" );
  printf ( "  Notice that both quantities eventually decrease\n" );
  printf ( "\n" );
  printf ( "  We use a radius of R = %g\n", r );
  printf ( "\n" );
  printf ( "    M        Area          Volume    Area / Volume \n" );
  printf ( "\n" );

  for ( m = 1; m <= 20; m++ )
  {
    area = hypersphere_area ( m, r );
    volume = hypersphere_volume ( m, r );
    printf ( "  %3d  %14.6g  %14.6g  %14.6g\n", 
      m, area, volume, area / volume );
  }

  return;
}
/******************************************************************************/

void test06 ( )

/******************************************************************************/
/*
  Purpose:

    TEST06 tests the stereographic mapping.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 December 2013

  Author:

    John Burkardt
*/
{
  double err;
  int m;
  int n;
  int seed;
  int test;
  double *x1;
  double *x2;
  double *x3;

  printf ( "\n" );
  printf ( "TEST06\n" );
  printf ( "  Test the stereographic mapping:\n" );
  printf ( "  HYPERSPHERE_STEREOGRAPH maps hypersphere points to the plane.\n" );
  printf ( "  HYPERSPHERE_STEREOGRAPH_INVERSE inverts the mapping.\n" );
  printf ( "\n" );
  printf ( "  Pick a random X1 on the hypersphere.\n" );
  printf ( "  Map it to a point X2 on the plane.\n" );
  printf ( "  Map it back to a point X3 on the hypersphere.\n" );
  printf ( "  Consider norm of difference.\n" );
  printf ( "\n" );
  printf ( "  M    || X1 - X3 ||\n" );

  seed = 123456789;

  n = 1;
  for ( m = 2; m <= 5; m++ )
  {
    printf ( "\n" );
    for ( test = 1; test <= 5; test++ )
    {
      x1 = hypersphere_01_surface_uniform ( m, n, &seed );
      x2 = hypersphere_stereograph ( m, n, x1 );
      x3 = hypersphere_stereograph_inverse ( m, n, x2 );
      err = r8mat_norm_fro_affine ( m, n, x1, x3 );
      printf ( "  %2d  %14.6g\n", m, err );
      free ( x1 );
      free ( x2 );
      free ( x3 );
    }
  }
  return;
}

