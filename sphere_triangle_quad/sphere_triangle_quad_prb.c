# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "sphere_triangle_quad.h"

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );
void test04 ( );
void test05 ( );
void test06 ( );
void polyterm_exponent ( char *s, int e[3] );
double polyterm_value_3d ( double x[] );
/*
  Global data.
*/
int e_save[3];

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for SPHERE_TRIANGLE_QUAD_PRB.

  Discussion:

    SPHERE_TRIANGLE_QUAD_PRB tests the SPHERE_TRIANGLE_QUAD library.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    21 April 2014

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "SPHERE_TRIANGLE_QUAD_PRB\n" );
  printf ( "  C++ version\n" );
  printf ( "  Test the SPHERE_TRIANGLE_QUAD library.\n" );

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
  printf ( "SPHERE_TRIANGLE_QUAD_PRB\n" );
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

    TEST01 tests SPHERE01_TRIANGLE_QUAD_01, 02, 03.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    21 April 2014

  Author:

    John Burkardt
*/
{
  int e[3];
  int i;
  double result_01;
  double result_02;
  double result_03;
  int seed;
  double *v1;
  double *v2;
  double *v3;

  seed = 123456789;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  Approximate an integral on a random spherical triangle.\n" );
  printf ( "\n" );
  printf ( "  QUAD_01 uses centroids of spherical triangles.\n" );
  printf ( "  QUAD_02 uses vertices of spherical triangles.\n" );
  printf ( "  QUAD_03 uses midsides of spherical triangles.\n" );
/*
  Choose three points at random to define a spherical triangle.
*/
  v1 = sphere01_sample ( 1, &seed );
  v2 = sphere01_sample ( 1, &seed );
  v3 = sphere01_sample ( 1, &seed );

  printf ( "\n" );
  printf ( "  Vertices of random spherical triangle:\n" );
  printf ( "\n" );
  r8vec_transpose_print ( 3, v1, "  V1:" );
  r8vec_transpose_print ( 3, v2, "  V2:" );
  r8vec_transpose_print ( 3, v3, "  V3:" );

  printf ( "\n" );
  printf ( "QUAD_01      QUAD_02      QUAD_03\n" );

  for ( i = 1; i <= 17; i++ )
  {
    if ( i == 1 )
    {
     e[0] = 0;
     e[1] = 0;
     e[2] = 0;
    }
    else if ( i == 2 )
    {
     e[0] = 1;
     e[1] = 0;
     e[2] = 0;
    }
    else if ( i == 3 )
    {
     e[0] = 0;
     e[1] = 1;
     e[2] = 0;
    }
    else if ( i == 4 )
    {
     e[0] = 0;
     e[1] = 0;
     e[2] = 1;
    }
    else if ( i == 5 )
    {
     e[0] = 2;
     e[1] = 0;
     e[2] = 0;
    }
    else if ( i == 6 )
    {
     e[0] = 0;
     e[1] = 2;
     e[2] = 2;
    }
    else if ( i == 7 )
    {
     e[0] = 2;
     e[1] = 2;
     e[2] = 2;
    }
    else if ( i == 8 )
    {
     e[0] = 0;
     e[1] = 2;
     e[2] = 4;
    }
    else if ( i == 9 )
    {
     e[0] = 0;
     e[1] = 0;
     e[2] = 6;
    }
    else if ( i == 10 )
    {
     e[0] = 1;
     e[1] = 2;
     e[2] = 4;
    }
    else if ( i == 11 )
    {
     e[0] = 2;
     e[1] = 4;
     e[2] = 2;
    }
    else if ( i == 12 )
    {
     e[0] = 6;
     e[1] = 2;
     e[2] = 0;
    }
    else if ( i == 13 )
    {
     e[0] = 0;
     e[1] = 0;
     e[2] = 8;
    }
    else if ( i == 14 )
    {
     e[0] = 6;
     e[1] = 0;
     e[2] = 4;
    }
    else if ( i == 15 )
    {
     e[0] = 4;
     e[1] = 6;
     e[2] = 2;
    }
    else if ( i == 16 )
    {
     e[0] = 2;
     e[1] = 4;
     e[2] = 8;
    }
    else if ( i == 17 )
    {
     e[0] = 16;
     e[1] = 0;
     e[2] = 0;
    }

    polyterm_exponent ( "SET", e );

    polyterm_exponent ( "PRINT", e );

    result_01 = sphere01_triangle_quad_01 ( v1, v2, v3, polyterm_value_3d );

    result_02 = sphere01_triangle_quad_02 ( v1, v2, v3, polyterm_value_3d );

    result_03 = sphere01_triangle_quad_03 ( v1, v2, v3, polyterm_value_3d );

    printf ( "  %14g  %14g  %14g\n", result_01, result_02, result_03 );

  }

  return;
}
/******************************************************************************/

void test02 ( )

/******************************************************************************/
/*
  Purpose:

    TEST02 tests SPHERE01_TRIANGLE_QUAD_00.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    21 April 2014

  Author:

    John Burkardt
*/
{
  int e[3];
  int i;
  int n_mc1 = 100;
  int n_mc2 = 1000;
  int n_mc3 = 10000;
  double result_01;
  double result_02;
  double result_03;
  int seed;
  double *v1;
  double *v2;
  double *v3;

  seed = 123456789;

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  Approximate the integral of a function on a random spherical triangle.\n" );
  printf ( "\n" );
  printf ( "  SPHERE01_TRIANGLE_QUAD_00 uses the Monte Carlo method.\n" );
  printf ( "  QUAD_MC1 uses a Monte Carlo method with %d points.\n", n_mc1 );
  printf ( "  QUAD_MC2 uses a Monte Carlo method with %d points.\n", n_mc2 );
  printf ( "  QUAD_MC3 uses a Monte Carlo method with %d points.\n", n_mc3 );
/*
  Choose three points at random to define a spherical triangle.
*/
  v1 = sphere01_sample ( 1, &seed );
  v2 = sphere01_sample ( 1, &seed );
  v3 = sphere01_sample ( 1, &seed );

  printf ( "\n" );
  printf ( "  Vertices of random spherical triangle:\n" );
  printf ( "\n" );
  r8vec_transpose_print ( 3, v1, "  V1:" );
  r8vec_transpose_print ( 3, v2, "  V2:" );
  r8vec_transpose_print ( 3, v3, "  V3:" );

  printf ( "\n" );
  printf ( "  QUAD_MC1     QUAD_MC2     QUAD_MC3\n" );

  for ( i = 1; i <= 17; i++ )
  {
    if ( i == 1 )
    {
     e[0] = 0;
     e[1] = 0;
     e[2] = 0;
    }
    else if ( i == 2 )
    {
     e[0] = 1;
     e[1] = 0;
     e[2] = 0;
    }
    else if ( i == 3 )
    {
     e[0] = 0;
     e[1] = 1;
     e[2] = 0;
    }
    else if ( i == 4 )
    {
     e[0] = 0;
     e[1] = 0;
     e[2] = 1;
    }
    else if ( i == 5 )
    {
     e[0] = 2;
     e[1] = 0;
     e[2] = 0;
    }
    else if ( i == 6 )
    {
     e[0] = 0;
     e[1] = 2;
     e[2] = 2;
    }
    else if ( i == 7 )
    {
     e[0] = 2;
     e[1] = 2;
     e[2] = 2;
    }
    else if ( i == 8 )
    {
     e[0] = 0;
     e[1] = 2;
     e[2] = 4;
    }
    else if ( i == 9 )
    {
     e[0] = 0;
     e[1] = 0;
     e[2] = 6;
    }
    else if ( i == 10 )
    {
     e[0] = 1;
     e[1] = 2;
     e[2] = 4;
    }
    else if ( i == 11 )
    {
     e[0] = 2;
     e[1] = 4;
     e[2] = 2;
    }
    else if ( i == 12 )
    {
     e[0] = 6;
     e[1] = 2;
     e[2] = 0;
    }
    else if ( i == 13 )
    {
     e[0] = 0;
     e[1] = 0;
     e[2] = 8;
    }
    else if ( i == 14 )
    {
     e[0] = 6;
     e[1] = 0;
     e[2] = 4;
    }
    else if ( i == 15 )
    {
     e[0] = 4;
     e[1] = 6;
     e[2] = 2;
    }
    else if ( i == 16 )
    {
     e[0] = 2;
     e[1] = 4;
     e[2] = 8;
    }
    else if ( i == 17 )
    {
     e[0] = 16;
     e[1] = 0;
     e[2] = 0;
    }

    polyterm_exponent ( "SET", e );

    polyterm_exponent ( "PRINT", e );

    result_01 = sphere01_triangle_quad_00 ( n_mc1, v1, v2, v3, 
      polyterm_value_3d, &seed );

    result_02 = sphere01_triangle_quad_00 ( n_mc2, v1, v2, v3, 
      polyterm_value_3d, &seed );

    result_03 = sphere01_triangle_quad_00 ( n_mc3, v1, v2, v3, 
      polyterm_value_3d, &seed );

    printf ( "  %14g  %14g  %14g\n", result_01, result_02, result_03 );

  }

  return;
}
/******************************************************************************/

void test03 ( )

/******************************************************************************/
/*
  Purpose:

    TEST03 tests SPHERE01_TRIANGLE_QUAD_ICOS1C.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    21 April 2014

  Author:

    John Burkardt
*/
{
  double best;
  int e[3];
  double error;
  int factor;
  int factor_log;
  int i;
  int node_num;
  double result;
  int seed;
  double *v1;
  double *v2;
  double *v3;

  seed = 123456789;

  printf ( "\n" );
  printf ( "TEST03\n" );
  printf ( "  Approximate the integral of a function on a random spherical triangle.\n" );
  printf ( "\n" );
  printf ( "  SPHERE01_TRIANGLE_QUAD_ICOS1C approximates the\n" );
  printf ( "  integral of a function over a spherical triangle on\n" );
  printf ( "  the surface of the unit sphere using a centroid rule.\n" );
  printf ( "\n" );
  printf ( "  We do not have an exact result, so we compare each\n" );
  printf ( "  estimate to the final one.\n" );
/*
  Choose three points at random to define a spherical triangle.
*/
  v1 = sphere01_sample ( 1, &seed );
  v2 = sphere01_sample ( 1, &seed );
  v3 = sphere01_sample ( 1, &seed );

  printf ( "\n" );
  printf ( "  Vertices of random spherical triangle:\n" );
  printf ( "\n" );
  r8vec_transpose_print ( 3, v1, "  V1:" );
  r8vec_transpose_print ( 3, v2, "  V2:" );
  r8vec_transpose_print ( 3, v3, "  V3:" );

  printf ( "\n" );
  printf ( "FACTOR   N   RESULT\n" );

  for ( i = 1; i <= 17; i++ )
  {
    if ( i == 1 )
    {
     e[0] = 0;
     e[1] = 0;
     e[2] = 0;
    }
    else if ( i == 2 )
    {
     e[0] = 1;
     e[1] = 0;
     e[2] = 0;
    }
    else if ( i == 3 )
    {
     e[0] = 0;
     e[1] = 1;
     e[2] = 0;
    }
    else if ( i == 4 )
    {
     e[0] = 0;
     e[1] = 0;
     e[2] = 1;
    }
    else if ( i == 5 )
    {
     e[0] = 2;
     e[1] = 0;
     e[2] = 0;
    }
    else if ( i == 6 )
    {
     e[0] = 0;
     e[1] = 2;
     e[2] = 2;
    }
    else if ( i == 7 )
    {
     e[0] = 2;
     e[1] = 2;
     e[2] = 2;
    }
    else if ( i == 8 )
    {
     e[0] = 0;
     e[1] = 2;
     e[2] = 4;
    }
    else if ( i == 9 )
    {
     e[0] = 0;
     e[1] = 0;
     e[2] = 6;
    }
    else if ( i == 10 )
    {
     e[0] = 1;
     e[1] = 2;
     e[2] = 4;
    }
    else if ( i == 11 )
    {
     e[0] = 2;
     e[1] = 4;
     e[2] = 2;
    }
    else if ( i == 12 )
    {
     e[0] = 6;
     e[1] = 2;
     e[2] = 0;
    }
    else if ( i == 13 )
    {
     e[0] = 0;
     e[1] = 0;
     e[2] = 8;
    }
    else if ( i == 14 )
    {
     e[0] = 6;
     e[1] = 0;
     e[2] = 4;
    }
    else if ( i == 15 )
    {
     e[0] = 4;
     e[1] = 6;
     e[2] = 2;
    }
    else if ( i == 16 )
    {
     e[0] = 2;
     e[1] = 4;
     e[2] = 8;
    }
    else if ( i == 17 )
    {
     e[0] = 16;
     e[1] = 0;
     e[2] = 0;
    }

    polyterm_exponent ( "SET", e );

    polyterm_exponent ( "PRINT", e );

    factor = i4_power ( 2, 9 );

    best = sphere01_triangle_quad_icos1c ( v1, v2, v3, factor, 
      polyterm_value_3d, &node_num  );

    factor = 1;

    for ( factor_log = 0; factor_log <= 9; factor_log++ )
    {
      result = sphere01_triangle_quad_icos1c ( v1, v2, v3, factor, 
        polyterm_value_3d, &node_num );

      error = r8_abs ( result - best );

      printf ( "  %4d  %8d  %16g  %10g\n", factor, node_num, result, error );

      factor = factor * 2;
    }
  }

  return;
}
/******************************************************************************/

void test04 ( )

/******************************************************************************/
/*
  Purpose:

    TEST04 tests SPHERE01_TRIANGLE_QUAD_ICOS1M.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    21 April 2014

  Author:

    John Burkardt
*/
{
  double best;
  int e[3];
  double error;
  int factor;
  int factor_log;
  int i;
  int node_num;
  double result;
  int seed;
  double *v1;
  double *v2;
  double *v3;

  seed = 123456789;

  printf ( "\n" );
  printf ( "TEST04\n" );
  printf ( "  Approximate the integral of a function on a random spherical triangle.\n" );
  printf ( "\n" );
  printf ( "  SPHERE01_TRIANGLE_QUAD_ICOS1M approximates the\n" );
  printf ( "  integral of a function over a spherical triangle on\n" );
  printf ( "  the surface of the unit sphere using a midside rule.\n" );
  printf ( "\n" );
  printf ( "  We do not have an exact result, so we compare each\n" );
  printf ( "  estimate to the final one.\n" );
/*
  Choose three points at random to define a spherical triangle.
*/
  v1 = sphere01_sample ( 1, &seed );
  v2 = sphere01_sample ( 1, &seed );
  v3 = sphere01_sample ( 1, &seed );

  printf ( "\n" );
  printf ( "  Vertices of random spherical triangle:\n" );
  printf ( "\n" );
  r8vec_transpose_print ( 3, v1, "  V1:" );
  r8vec_transpose_print ( 3, v2, "  V2:" );
  r8vec_transpose_print ( 3, v3, "  V3:" );

  printf ( "\n" );
  printf ( "FACTOR   N   RESULT\n" );

  for ( i = 1; i <= 17; i++ )
  {
    if ( i == 1 )
    {
     e[0] = 0;
     e[1] = 0;
     e[2] = 0;
    }
    else if ( i == 2 )
    {
     e[0] = 1;
     e[1] = 0;
     e[2] = 0;
    }
    else if ( i == 3 )
    {
     e[0] = 0;
     e[1] = 1;
     e[2] = 0;
    }
    else if ( i == 4 )
    {
     e[0] = 0;
     e[1] = 0;
     e[2] = 1;
    }
    else if ( i == 5 )
    {
     e[0] = 2;
     e[1] = 0;
     e[2] = 0;
    }
    else if ( i == 6 )
    {
     e[0] = 0;
     e[1] = 2;
     e[2] = 2;
    }
    else if ( i == 7 )
    {
     e[0] = 2;
     e[1] = 2;
     e[2] = 2;
    }
    else if ( i == 8 )
    {
     e[0] = 0;
     e[1] = 2;
     e[2] = 4;
    }
    else if ( i == 9 )
    {
     e[0] = 0;
     e[1] = 0;
     e[2] = 6;
    }
    else if ( i == 10 )
    {
     e[0] = 1;
     e[1] = 2;
     e[2] = 4;
    }
    else if ( i == 11 )
    {
     e[0] = 2;
     e[1] = 4;
     e[2] = 2;
    }
    else if ( i == 12 )
    {
     e[0] = 6;
     e[1] = 2;
     e[2] = 0;
    }
    else if ( i == 13 )
    {
     e[0] = 0;
     e[1] = 0;
     e[2] = 8;
    }
    else if ( i == 14 )
    {
     e[0] = 6;
     e[1] = 0;
     e[2] = 4;
    }
    else if ( i == 15 )
    {
     e[0] = 4;
     e[1] = 6;
     e[2] = 2;
    }
    else if ( i == 16 )
    {
     e[0] = 2;
     e[1] = 4;
     e[2] = 8;
    }
    else if ( i == 17 )
    {
     e[0] = 16;
     e[1] = 0;
     e[2] = 0;
    }

    polyterm_exponent ( "SET", e );

    polyterm_exponent ( "PRINT", e );

    factor = i4_power ( 2, 9 );

    best = sphere01_triangle_quad_icos1m ( v1, v2, v3, factor, 
      polyterm_value_3d, &node_num  );

    factor = 1;

    for ( factor_log = 0; factor_log <= 9; factor_log++ )
    {
      result = sphere01_triangle_quad_icos1m ( v1, v2, v3, factor, 
        polyterm_value_3d, &node_num );

      error = r8_abs ( result - best );

      printf ( "  %4d  %8d  %16g  %10e\n", factor, node_num, result, error );

      factor = factor * 2;
    }
  }

  return;
}
/******************************************************************************/

void test05 ( )

/******************************************************************************/
/*
  Purpose:

    TEST05 tests SPHERE01_TRIANGLE_QUAD_ICOS1V.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    21 April 2014

  Author:

    John Burkardt
*/
{
  double best;
  int e[3];
  double error;
  int factor;
  int factor_log;
  int i;
  int node_num;
  double result;
  int seed;
  double *v1;
  double *v2;
  double *v3;

  seed = 123456789;

  printf ( "\n" );
  printf ( "TEST05\n" );
  printf ( "  Approximate an integral on a random spherical triangle.\n" );
  printf ( "\n" );
  printf ( "  SPHERE01_TRIANGLE_QUAD_ICOS1V approximates the\n" );
  printf ( "  integral of a function over a spherical triangle on\n" );
  printf ( "  the surface of the unit sphere using a vertex rule.\n" );
  printf ( "\n" );
  printf ( "  We do not have an exact result, so we compare each\n" );
  printf ( "  estimate to the final one.\n" );
/*
  Choose three points at random to define a spherical triangle.
*/
  v1 = sphere01_sample ( 1, &seed );
  v2 = sphere01_sample ( 1, &seed );
  v3 = sphere01_sample ( 1, &seed );

  printf ( "\n" );
  printf ( "  Vertices of random spherical triangle:\n" );
  printf ( "\n" );
  r8vec_transpose_print ( 3, v1, "  V1:" );
  r8vec_transpose_print ( 3, v2, "  V2:" );
  r8vec_transpose_print ( 3, v3, "  V3:" );

  printf ( "\n" );
  printf ( "FACTOR   N   RESULT\n" );

  for ( i = 1; i <= 17; i++ )
  {
    if ( i == 1 )
    {
     e[0] = 0;
     e[1] = 0;
     e[2] = 0;
    }
    else if ( i == 2 )
    {
     e[0] = 1;
     e[1] = 0;
     e[2] = 0;
    }
    else if ( i == 3 )
    {
     e[0] = 0;
     e[1] = 1;
     e[2] = 0;
    }
    else if ( i == 4 )
    {
     e[0] = 0;
     e[1] = 0;
     e[2] = 1;
    }
    else if ( i == 5 )
    {
     e[0] = 2;
     e[1] = 0;
     e[2] = 0;
    }
    else if ( i == 6 )
    {
     e[0] = 0;
     e[1] = 2;
     e[2] = 2;
    }
    else if ( i == 7 )
    {
     e[0] = 2;
     e[1] = 2;
     e[2] = 2;
    }
    else if ( i == 8 )
    {
     e[0] = 0;
     e[1] = 2;
     e[2] = 4;
    }
    else if ( i == 9 )
    {
     e[0] = 0;
     e[1] = 0;
     e[2] = 6;
    }
    else if ( i == 10 )
    {
     e[0] = 1;
     e[1] = 2;
     e[2] = 4;
    }
    else if ( i == 11 )
    {
     e[0] = 2;
     e[1] = 4;
     e[2] = 2;
    }
    else if ( i == 12 )
    {
     e[0] = 6;
     e[1] = 2;
     e[2] = 0;
    }
    else if ( i == 13 )
    {
     e[0] = 0;
     e[1] = 0;
     e[2] = 8;
    }
    else if ( i == 14 )
    {
     e[0] = 6;
     e[1] = 0;
     e[2] = 4;
    }
    else if ( i == 15 )
    {
     e[0] = 4;
     e[1] = 6;
     e[2] = 2;
    }
    else if ( i == 16 )
    {
     e[0] = 2;
     e[1] = 4;
     e[2] = 8;
    }
    else if ( i == 17 )
    {
     e[0] = 16;
     e[1] = 0;
     e[2] = 0;
    }

    polyterm_exponent ( "SET", e );

    polyterm_exponent ( "PRINT", e );

    factor = i4_power ( 2, 9 );

    best = sphere01_triangle_quad_icos1v ( v1, v2, v3, factor, 
      polyterm_value_3d, &node_num  );

    factor = 1;

    for ( factor_log = 0; factor_log <= 9; factor_log++ )
    {
      result = sphere01_triangle_quad_icos1v ( v1, v2, v3, factor, 
        polyterm_value_3d, &node_num );

      error = r8_abs ( result - best );

      printf ( "  %4d  %8d  %16g  %10e\n", factor, node_num, result, error );

      factor = factor * 2;
    }
  }

  return;
}
/******************************************************************************/

void test06 ( )

/******************************************************************************/
/*
  Purpose:

    TEST06 tests SPHERE01_TRIANGLE_QUAD_ICOS2V.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    21 April 2014

  Author:

    John Burkardt
*/
{
  double best;
  int e[3];
  double error;
  int factor;
  int factor_log;
  int i;
  int node_num;
  double result;
  int seed;
  double *v1;
  double *v2;
  double *v3;

  seed = 123456789;

  printf ( "\n" );
  printf ( "TEST06\n" );
  printf ( "  Approximate the integral of a function on a random spherical triangle.\n" );
  printf ( "\n" );
  printf ( "  SPHERE01_TRIANGLE_QUAD_ICOS2V approximates the\n" );
  printf ( "  integral of a function over a spherical triangle on\n" );
  printf ( "  the surface of the unit sphere using a vertex rule.\n" );
  printf ( "\n" );
  printf ( "  We do not have an exact result, so we compare each\n" );
  printf ( "  estimate to the final one.\n" );
/*
  Choose three points at random to define a spherical triangle.
*/
  v1 = sphere01_sample ( 1, &seed );
  v2 = sphere01_sample ( 1, &seed );
  v3 = sphere01_sample ( 1, &seed );

  printf ( "\n" );
  printf ( "  Vertices of random spherical triangle:\n" );
  printf ( "\n" );
  r8vec_transpose_print ( 3, v1, "  V1:" );
  r8vec_transpose_print ( 3, v2, "  V2:" );
  r8vec_transpose_print ( 3, v3, "  V3:" );

  printf ( "\n" );
  printf ( "FACTOR   N   RESULT\n" );

  for ( i = 1; i <= 17; i++ )
  {
    if ( i == 1 )
    {
     e[0] = 0;
     e[1] = 0;
     e[2] = 0;
    }
    else if ( i == 2 )
    {
     e[0] = 1;
     e[1] = 0;
     e[2] = 0;
    }
    else if ( i == 3 )
    {
     e[0] = 0;
     e[1] = 1;
     e[2] = 0;
    }
    else if ( i == 4 )
    {
     e[0] = 0;
     e[1] = 0;
     e[2] = 1;
    }
    else if ( i == 5 )
    {
     e[0] = 2;
     e[1] = 0;
     e[2] = 0;
    }
    else if ( i == 6 )
    {
     e[0] = 0;
     e[1] = 2;
     e[2] = 2;
    }
    else if ( i == 7 )
    {
     e[0] = 2;
     e[1] = 2;
     e[2] = 2;
    }
    else if ( i == 8 )
    {
     e[0] = 0;
     e[1] = 2;
     e[2] = 4;
    }
    else if ( i == 9 )
    {
     e[0] = 0;
     e[1] = 0;
     e[2] = 6;
    }
    else if ( i == 10 )
    {
     e[0] = 1;
     e[1] = 2;
     e[2] = 4;
    }
    else if ( i == 11 )
    {
     e[0] = 2;
     e[1] = 4;
     e[2] = 2;
    }
    else if ( i == 12 )
    {
     e[0] = 6;
     e[1] = 2;
     e[2] = 0;
    }
    else if ( i == 13 )
    {
     e[0] = 0;
     e[1] = 0;
     e[2] = 8;
    }
    else if ( i == 14 )
    {
     e[0] = 6;
     e[1] = 0;
     e[2] = 4;
    }
    else if ( i == 15 )
    {
     e[0] = 4;
     e[1] = 6;
     e[2] = 2;
    }
    else if ( i == 16 )
    {
     e[0] = 2;
     e[1] = 4;
     e[2] = 8;
    }
    else if ( i == 17 )
    {
     e[0] = 16;
     e[1] = 0;
     e[2] = 0;
    }

    polyterm_exponent ( "SET", e );

    polyterm_exponent ( "PRINT", e );

    factor = i4_power ( 2, 9 );

    best = sphere01_triangle_quad_icos2v ( v1, v2, v3, factor, 
      polyterm_value_3d, &node_num  );

    factor = 1;

    for ( factor_log = 0; factor_log <= 9; factor_log++ )
    {
      result = sphere01_triangle_quad_icos2v ( v1, v2, v3, factor, 
        polyterm_value_3d, &node_num );

      error = r8_abs ( result - best );

      printf ( "  %4d  %8d  %16g  %10e\n", factor, node_num, result, error );

      factor = factor * 2;
    }
  }

  return;
}
/******************************************************************************/

void polyterm_exponent ( char *action, int e[3] )

/******************************************************************************/
/*
  Purpose:

    POLYTERM_EXPONENT gets or sets the exponents for the polynomial term.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    21 April 2014

  Author:

    John Burkardt

  Parameters:

    Input, char *ACTION.
    'GET' asks the routine to return the current values in E.
    'SET' asks the routine to set the current values to E.

    Input/output, int E[3], storage used to set or get values.
*/
{
  int i;

  if ( action[0] == 'G' )
  {
    for ( i = 0; i < 3; i++ )
    {
      e[i] = e_save[i];
    }
  }
  else if ( action[0] == 'P' )
  {
    printf ( "\n" );

    if ( e_save[0] == 0 && e_save[1] == 0 && e_save[2] == 0 )
    {
      printf ( "P(X,Y,Z) = 1\n" );
    }
    else
    {
      printf ( "P(X,Y,Z) = " );

      if ( e_save[0] == 0 )
      {
      }
      else if ( e_save[0] == 1 )
      {
        printf ( " X" );
      }
      else
      {
        printf ( " X^%d", e_save[0] );
      }
      if ( e_save[1] == 0 )
      {
      }
      else if ( e_save[1] == 1 )
      {
        printf ( " Y" );
      }
      else
      {
        printf ( " Y^%d", e_save[1] );
      }
      if ( e_save[2] == 0 )
      {
      }
      else if ( e_save[2] == 1 )
      {
        printf ( " Z" );
      }
      else
      {
        printf ( " Z^%d", e_save[2] );
      }
      printf ( "\n" );
    }
  } 
  else if ( action[0] == 'S' )
  {
    for ( i = 0; i < 3; i++ )
    {
      e_save[i] = e[i];
    }
  }

  return;
}
/******************************************************************************/

double polyterm_value_3d ( double x[] )

/******************************************************************************/
/*
  Purpose:

    POLYTERM_VALUE_3D evaluates a single polynomial term in 3D.

  Discussion:

    The polynomial term has the form:

      F(X) = X(1)^E(1) * X(2)^E(2) * X(3)^E(3)

    The exponents E(1:3) are set by calling POLYTERM_EXPONENT.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    21 April 2014

  Author:

    John Burkardt

  Parameters:

    Input, double X[3], the points where the polynomial term 
    is to be evaluated.

    Output, double POLYTERM_VALUE_3D, the value of the polynomial term.
*/
{
  int e[3];
  int i;
  double value;

  polyterm_exponent ( "GET", e );

  value = 1.0;
  for ( i = 0; i < 3; i++ )
  {
    if ( e[i] != 0 )
    {
      value = value * pow ( x[i], e[i] );
    }
  }
  
  return value;
}
