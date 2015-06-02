# include <stdlib.h>
# include <stdio.h>

# include "line_cvt_lloyd.h"

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );
void test04 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    LINE_CVT_LLOYD_PRB tests the line_cvt_lloyd library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 July 2014

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "LINE_CVT_LLOYD_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the LINE_CVT_LLOYD library.\n" );

  test01 ( );
  test02 ( );
/*
  Repeat, using sorted initial points.
*/
  test03 ( );
  test04 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "LINE_CVT_LLOYD_PRB\n" );
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

    LINE_CVT_LLOYD_TEST01 tests the unconstrained computation.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    29 July 2014

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double h;
  char header[] = "test01";
  int it_num;
  int n = 25;
  int seed;
  double *x;

  printf ( "\n" );
  printf ( "LINE_CVT_LLOYD_TEST01:\n" );
  printf ( "  Test the unconstrained computation.\n" );

  a = 0.0;
  b = 1.0;
  it_num = 200;
  seed = 123456789;
  x = r8vec_uniform_ab_new ( n, a, b, &seed );

  printf ( "\n" );
  printf ( "  Use %d points in the interval [%f,%f]\n", n, a, b );
  printf ( "  Number of iterations to take is %d\n", it_num );
  printf ( "  Call this calculation '%s'\n", header );
  h = ( b - a ) / ( double ) ( n - 1 );
  printf ( "  Expect a uniform spacing of %g\n", h );

  r8vec_print ( n, x, "  Initial generators:" );

  line_cvt_lloyd ( n, a, b, it_num, header, x );

  r8vec_print ( n, x, "  Final generators:" );

  free ( x );

  return;
}
/******************************************************************************/

void test02 ( )

/******************************************************************************/
/*
  Purpose:

    LINE_CVT_LLOYD_TEST02 tests the constrained computation.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    29 July 2014

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double h;
  char header[] = "test02";
  int it_num;
  int n = 25;
  int seed;
  double *x;

  printf ( "\n" );
  printf ( "LINE_CVT_LLOYD_TEST02:\n" );
  printf ( "  Test the constrained computation.\n" );

  a = 0.0;
  b = 1.0;
  it_num = 200;
  seed = 123456789;
  x = r8vec_uniform_ab_new ( n, a, b, &seed );

  printf ( "\n" );
  printf ( "  Use %d points in the interval [%f,%f]\n", n, a, b );
  printf ( "  Number of iterations to take is %d\n", it_num );
  printf ( "  Call this calculation '%s'\n", header );
  h = ( b - a ) / ( double ) ( n );
  printf ( "  Expect a uniform spacing of %g\n", h );

  r8vec_print ( n, x, "  Initial generators:" );

  line_ccvt_lloyd ( n, a, b, it_num, header, x );

  r8vec_print ( n, x, "  Final generators:" );

  free ( x );

  return;
}
/******************************************************************************/

void test03 ( )

/******************************************************************************/
/*
  Purpose:

    LINE_CVT_LLOYD_TEST03 tests the unconstrained computation.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    29 July 2014

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double h;
  char header[] = "test03";
  int it_num;
  int n = 25;
  int seed;
  double *x;

  printf ( "\n" );
  printf ( "LINE_CVT_LLOYD_TEST03:\n" );
  printf ( "  Test the unconstrained computation.\n" );
  printf ( "  SORT the random initial values before use.\n" );

  a = 0.0;
  b = 1.0;
  it_num = 200;
  seed = 123456789;
  x = r8vec_uniform_ab_new ( n, a, b, &seed );
  r8vec_sort_insert_a ( n, x );

  printf ( "\n" );
  printf ( "  Use %d points in the interval [%f,%f]\n", n, a, b );
  printf ( "  Number of iterations to take is %d\n", it_num );
  printf ( "  Call this calculation '%s'\n", header );
  h = ( b - a ) / ( double ) ( n - 1 );
  printf ( "  Expect a uniform spacing of %g\n", h );

  r8vec_print ( n, x, "  Initial generators:" );

  line_cvt_lloyd ( n, a, b, it_num, header, x );

  r8vec_print ( n, x, "  Final generators:" );

  free ( x );

  return;
}
/******************************************************************************/

void test04 ( )

/******************************************************************************/
/*
  Purpose:

    LINE_CVT_LLOYD_TEST04 tests the constrained computation.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    29 July 2014

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double h;
  char header[] = "test04";
  int it_num;
  int n = 25;
  int seed;
  double *x;

  printf ( "\n" );
  printf ( "LINE_CVT_LLOYD_TEST04:\n" );
  printf ( "  Test the constrained computation.\n" );
  printf ( "  SORT the initial points before use.\n" );

  a = 0.0;
  b = 1.0;
  it_num = 200;
  seed = 123456789;
  x = r8vec_uniform_ab_new ( n, a, b, &seed );
  r8vec_sort_insert_a ( n, x );

  printf ( "\n" );
  printf ( "  Use %d points in the interval [%f,%f]\n", n, a, b );
  printf ( "  Number of iterations to take is %d\n", it_num );
  printf ( "  Call this calculation '%s'\n", header );
  h = ( b - a ) / ( double ) ( n );
  printf ( "  Expect a uniform spacing of %g\n", h );

  r8vec_print ( n, x, "  Initial generators:" );

  line_ccvt_lloyd ( n, a, b, it_num, header, x );

  r8vec_print ( n, x, "  Final generators:" );

  free ( x );

  return;
}
