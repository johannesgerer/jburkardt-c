# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "wavelet.h"

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );
void test04 ( );
void test05 ( );
void test06 ( );
void test07 ( );
void test08 ( );
void test09 ( );
void test10 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for WAVELET_PRB.

  Discussion:

    WAVELET_PRB tests the WAVELET library.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    09 May 2012

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "WAVELET_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the WAVELET library.\n" );

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
  printf ( "WAVELET_PRB\n" );
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

    TEST01 tests DAUB2_TRANSFORM and DAUB2_TRANSFORM_INVERSE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    28 April 2012

  Author:

    John Burkardt
*/
{
  double a_first;
  double a_last;
  int i;
  int n;
  int seed;
  double *u;
  double *v;
  double *w;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  DAUB2_TRANSFORM computes the DAUB2 transform of a vector.\n" );
  printf ( "  DAUB2_TRANSFORM_INVERSE inverts it.\n" );
/*
  Random data.
*/
  n = 16;
  seed = 123456789;

  u = r8vec_uniform_01_new ( n, &seed );

  v = daub2_transform ( n, u );

  w = daub2_transform_inverse ( n, v );

  printf ( "\n" );
  printf ( "   i      U          D2(U)    D2inv(D2(U))\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %2d  %10g  %10g  %10g\n", i, u[i], v[i], w[i] );
  }

  free ( u );
  free ( v );
  free ( w );
/*
  Constant signal.
*/
  n = 8;
  u = ( double * ) malloc ( n * sizeof ( double ) );
  for ( i = 0; i < n; i++ )
  {
    u[i] = 1.0;
  }

  v = daub2_transform ( n, u );

  w = daub2_transform_inverse ( n, v );

  printf ( "\n" );
  printf ( "   i      U          D2(U)    D2inv(D2(U))\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %2d  %10g  %10g  %10g\n", i, u[i], v[i], w[i] );
  }

  free ( u );
  free ( v );
  free ( w );
/*
  Linear signal.
*/
  n = 16;
  a_first = 1.0;
  a_last = ( double ) ( n );
  u = r8vec_linspace_new ( n, a_first, a_last );

  v = daub2_transform ( n, u );

  w = daub2_transform_inverse ( n, v );

  printf ( "\n" );
  printf ( "   i      U          D2(U)    D2inv(D2(U))\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %2d  %10g  %10g  %10g\n", i, u[i], v[i], w[i] );
  }

  free ( u );
  free ( v );
  free ( w );
/*
  Quadratic data.
*/
  n = 8;
  u = ( double * ) malloc ( n * sizeof ( double ) );
  for ( i = 0; i < n; i++ )
  {
    u[i] = pow ( ( double ) ( i - 5 ), 2 );
  }

  v = daub2_transform ( n, u );

  w = daub2_transform_inverse ( n, v );

  printf ( "\n" );
  printf ( "   i      U          D2(U)    D2inv(D2(U))\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %2d  %10g  %10g  %10g\n", i, u[i], v[i], w[i] );
  }

  free ( u );
  free ( v );
  free ( w );

  return;
}
/******************************************************************************/

void test02 ( )

/******************************************************************************/
/*
  Purpose:

    TEST02 tests DAUB4_TRANSFORM and DAUB4_TRANSFORM_INVERSE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    28 April 2012

  Author:

    John Burkardt
*/
{
  double a_first;
  double a_last;
  int i;
  int n;
  int seed;
  double *u;
  double *v;
  double *w;

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  DAUB4_TRANSFORM computes the DAUB4 transform of a vector.\n" );
  printf ( "  DAUB4_TRANSFORM_INVERSE inverts it.\n" );
/*
  Random data.
*/
  n = 16;
  seed = 123456789;

  u = r8vec_uniform_01_new ( n, &seed );

  v = daub4_transform ( n, u );

  w = daub4_transform_inverse ( n, v );

  printf ( "\n" );
  printf ( "   i      U          D4(U)    D4inv(D4(U))\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %2d  %10g  %10g  %10g\n", i, u[i], v[i], w[i] );
  }

  free ( u );
  free ( v );
  free ( w );
/*
  Constant signal.
*/
  n = 8;
  u = ( double * ) malloc ( n * sizeof ( double ) );
  for ( i = 0; i < n; i++ )
  {
    u[i] = 1.0;
  }

  v = daub4_transform ( n, u );

  w = daub4_transform_inverse ( n, v );

  printf ( "\n" );
  printf ( "   i      U          D4(U)    D4inv(D4(U))\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %2d  %10g  %10g  %10g\n", i, u[i], v[i], w[i] );
  }

  free ( u );
  free ( v );
  free ( w );
/*
  Linear signal.
*/
  n = 16;
  a_first = 1.0;
  a_last = ( double ) ( n );
  u = r8vec_linspace_new ( n, a_first, a_last );

  v = daub4_transform ( n, u );

  w = daub4_transform_inverse ( n, v );

  printf ( "\n" );
  printf ( "   i      U          D4(U)    D4inv(D4(U))\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %2d  %10g  %10g  %10g\n", i, u[i], v[i], w[i] );
  }

  free ( u );
  free ( v );
  free ( w );
/*
  Quadratic data.
*/
  n = 8;
  u = ( double * ) malloc ( n * sizeof ( double ) );
  for ( i = 0; i < n; i++ )
  {
    u[i] = pow ( ( double ) ( i - 5 ), 2 );
  }

  v = daub4_transform ( n, u );

  w = daub4_transform_inverse ( n, v );

  printf ( "\n" );
  printf ( "   i      U          D4(U)    D4inv(D4(U))\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %2d  %10g  %10g  %10g\n", i, u[i], v[i], w[i] );
  }

  free ( u );
  free ( v );
  free ( w );

  return;
}
/******************************************************************************/

void test03 ( )

/******************************************************************************/
/*
  Purpose:

    TEST03 tests DAUB6_TRANSFORM and DAUB6_TRANSFORM_INVERSE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 April 2012

  Author:

    John Burkardt
*/
{
  double a_first;
  double a_last;
  int i;
  int n;
  int seed;
  double *u;
  double *v;
  double *w;

  printf ( "\n" );
  printf ( "TEST03\n" );
  printf ( "  DAUB6_TRANSFORM computes the DAUB6 transform of a vector.\n" );
  printf ( "  DAUB6_TRANSFORM_INVERSE inverts it.\n" );
/*
  Random data.
*/
  n = 16;
  seed = 123456789;

  u = r8vec_uniform_01_new ( n, &seed );

  v = daub6_transform ( n, u );

  w = daub6_transform_inverse ( n, v );

  printf ( "\n" );
  printf ( "   i      U          D6(U)    D6inv(D6(U))\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %2d  %10g  %10g  %10g\n", i, u[i], v[i], w[i] );
  }

  free ( u );
  free ( v );
  free ( w );
/*
  Constant signal.
*/
  n = 8;
  u = ( double * ) malloc ( n * sizeof ( double ) );
  for ( i = 0; i < n; i++ )
  {
    u[i] = 1.0;
  }

  v = daub6_transform ( n, u );

  w = daub6_transform_inverse ( n, v );

  printf ( "\n" );
  printf ( "   i      U          D6(U)    D6inv(D6(U))\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %2d  %10g  %10g  %10g\n", i, u[i], v[i], w[i] );
  }

  free ( u );
  free ( v );
  free ( w );
/*
  Linear signal.
*/
  n = 16;
  a_first = 1.0;
  a_last = ( double ) ( n );
  u = r8vec_linspace_new ( n, a_first, a_last );

  v = daub6_transform ( n, u );

  w = daub6_transform_inverse ( n, v );

  printf ( "\n" );
  printf ( "   i      U          D6(U)    D6inv(D6(U))\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %2d  %10g  %10g  %10g\n", i, u[i], v[i], w[i] );
  }

  free ( u );
  free ( v );
  free ( w );
/*
  Quadratic data.
*/
  n = 8;
  u = ( double * ) malloc ( n * sizeof ( double ) );
  for ( i = 0; i < n; i++ )
  {
    u[i] = pow ( ( double ) ( i - 5 ), 2 );
  }

  v = daub6_transform ( n, u );

  w = daub6_transform_inverse ( n, v );

  printf ( "\n" );
  printf ( "   i      U          D6(U)    D6inv(D6(U))\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %2d  %10g  %10g  %10g\n", i, u[i], v[i], w[i] );
  }

  free ( u );
  free ( v );
  free ( w );

  return;
}
/******************************************************************************/

void test04 ( )

/******************************************************************************/
/*
  Purpose:

    TEST04 tests DAUB8_TRANSFORM and DAUB8_TRANSFORM_INVERSE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 April 2012

  Author:

    John Burkardt
*/
{
  double a_first;
  double a_last;
  int i;
  int n;
  int seed;
  double *u;
  double *v;
  double *w;

  printf ( "\n" );
  printf ( "TEST04\n" );
  printf ( "  DAUB8_TRANSFORM computes the DAUB8 transform of a vector.\n" );
  printf ( "  DAUB8_TRANSFORM_INVERSE inverts it.\n" );
/*
  Random data.
*/
  n = 16;
  seed = 123456789;

  u = r8vec_uniform_01_new ( n, &seed );

  v = daub8_transform ( n, u );

  w = daub8_transform_inverse ( n, v );

  printf ( "\n" );
  printf ( "   i      U          D8(U)    D8inv(D8(U))\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %2d  %10g  %10g  %10g\n", i, u[i], v[i], w[i] );
  }

  free ( u );
  free ( v );
  free ( w );
/*
  Constant signal.
*/
  n = 8;
  u = ( double * ) malloc ( n * sizeof ( double ) );
  for ( i = 0; i < n; i++ )
  {
    u[i] = 1.0;
  }

  v = daub8_transform ( n, u );

  w = daub8_transform_inverse ( n, v );

  printf ( "\n" );
  printf ( "   i      U          D8(U)    D8inv(D8(U))\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %2d  %10g  %10g  %10g\n", i, u[i], v[i], w[i] );
  }

  free ( u );
  free ( v );
  free ( w );
/*
  Linear signal.
*/
  n = 16;
  a_first = 1.0;
  a_last = ( double ) ( n );
  u = r8vec_linspace_new ( n, a_first, a_last );

  v = daub8_transform ( n, u );

  w = daub8_transform_inverse ( n, v );

  printf ( "\n" );
  printf ( "   i      U          D8(U)    D8inv(D8(U))\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %2d  %10g  %10g  %10g\n", i, u[i], v[i], w[i] );
  }

  free ( u );
  free ( v );
  free ( w );
/*
  Quadratic data.
*/
  n = 8;
  u = ( double * ) malloc ( n * sizeof ( double ) );
  for ( i = 0; i < n; i++ )
  {
    u[i] = pow ( ( double ) ( i - 5 ), 2 );
  }

  v = daub8_transform ( n, u );

  w = daub8_transform_inverse ( n, v );

  printf ( "\n" );
  printf ( "   i      U          D8(U)    D8inv(D8(U))\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %2d  %10g  %10g  %10g\n", i, u[i], v[i], w[i] );
  }

  free ( u );
  free ( v );
  free ( w );

  return;
}
/******************************************************************************/

void test05 ( )

/******************************************************************************/
/*
  Purpose:

    TEST05 tests DAUB10_TRANSFORM and DAUB10_TRANSFORM_INVERSE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    04 May 2012

  Author:

    John Burkardt
*/
{
  double a_first;
  double a_last;
  int i;
  int n;
  int seed;
  double *u;
  double *v;
  double *w;

  printf ( "\n" );
  printf ( "TEST05\n" );
  printf ( "  DAUB10_TRANSFORM computes the DAUB10 transform of a vector.\n" );
  printf ( "  DAUB10_TRANSFORM_INVERSE inverts it.\n" );
/*
  Random data.
*/
  n = 16;
  seed = 123456789;

  u = r8vec_uniform_01_new ( n, &seed );

  v = daub10_transform ( n, u );

  w = daub10_transform_inverse ( n, v );

  printf ( "\n" );
  printf ( "   i      U          D10(U)    D10inv(D10(U))\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %2d  %10g  %10g  %10g\n", i, u[i], v[i], w[i] );
  }

  free ( u );
  free ( v );
  free ( w );
/*
  Constant signal.
*/
  n = 8;
  u = ( double * ) malloc ( n * sizeof ( double ) );
  for ( i = 0; i < n; i++ )
  {
    u[i] = 1.0;
  }

  v = daub10_transform ( n, u );

  w = daub10_transform_inverse ( n, v );

  printf ( "\n" );
  printf ( "   i      U          D10(U)    D10inv(D10(U))\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %2d  %10g  %10g  %10g\n", i, u[i], v[i], w[i] );
  }

  free ( u );
  free ( v );
  free ( w );
/*
  Linear signal.
*/
  n = 16;
  a_first = 1.0;
  a_last = ( double ) ( n );
  u = r8vec_linspace_new ( n, a_first, a_last );

  v = daub10_transform ( n, u );

  w = daub10_transform_inverse ( n, v );

  printf ( "\n" );
  printf ( "   i      U          D10(U)    D10inv(D10(U))\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %2d  %10g  %10g  %10g\n", i, u[i], v[i], w[i] );
  }

  free ( u );
  free ( v );
  free ( w );
/*
  Quadratic data.
*/
  n = 8;
  u = ( double * ) malloc ( n * sizeof ( double ) );
  for ( i = 0; i < n; i++ )
  {
    u[i] = pow ( ( double ) ( i - 5 ), 2 );
  }

  v = daub10_transform ( n, u );

  w = daub10_transform_inverse ( n, v );

  printf ( "\n" );
  printf ( "   i      U          D10(U)    D10inv(D10(U))\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %2d  %10g  %10g  %10g\n", i, u[i], v[i], w[i] );
  }

  free ( u );
  free ( v );
  free ( w );

  return;
}
/******************************************************************************/

void test06 ( )

/******************************************************************************/
/*
  Purpose:

    TEST06 tests DAUB12_TRANSFORM and DAUB12_TRANSFORM_INVERSE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 May 2012

  Author:

    John Burkardt
*/
{
  double a_first;
  double a_last;
  int i;
  int n;
  int seed;
  double *u;
  double *v;
  double *w;

  printf ( "\n" );
  printf ( "TEST06\n" );
  printf ( "  DAUB12_TRANSFORM computes the DAUB12 transform of a vector.\n" );
  printf ( "  DAUB12_TRANSFORM_INVERSE inverts it.\n" );
/*
  Random data.
*/
  n = 16;
  seed = 123456789;

  u = r8vec_uniform_01_new ( n, &seed );

  v = daub12_transform ( n, u );

  w = daub12_transform_inverse ( n, v );

  printf ( "\n" );
  printf ( "   i      U          D12(U)    D12inv(D12(U))\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %2d  %10g  %10g  %10g\n", i, u[i], v[i], w[i] );
  }

  free ( u );
  free ( v );
  free ( w );
/*
  Constant signal.
*/
  n = 8;
  u = ( double * ) malloc ( n * sizeof ( double ) );
  for ( i = 0; i < n; i++ )
  {
    u[i] = 1.0;
  }

  v = daub12_transform ( n, u );

  w = daub12_transform_inverse ( n, v );

  printf ( "\n" );
  printf ( "   i      U          D12(U)    D12inv(D12(U))\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %2d  %10g  %10g  %10g\n", i, u[i], v[i], w[i] );
  }

  free ( u );
  free ( v );
  free ( w );
/*
  Linear signal.
*/
  n = 16;
  a_first = 1.0;
  a_last = ( double ) ( n );
  u = r8vec_linspace_new ( n, a_first, a_last );

  v = daub12_transform ( n, u );

  w = daub12_transform_inverse ( n, v );

  printf ( "\n" );
  printf ( "   i      U          D12(U)    D12inv(D12(U))\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %2d  %10g  %10g  %10g\n", i, u[i], v[i], w[i] );
  }

  free ( u );
  free ( v );
  free ( w );
/*
  Quadratic data.
*/
  n = 8;
  u = ( double * ) malloc ( n * sizeof ( double ) );
  for ( i = 0; i < n; i++ )
  {
    u[i] = pow ( ( double ) ( i - 5 ), 2 );
  }

  v = daub12_transform ( n, u );

  w = daub12_transform_inverse ( n, v );

  printf ( "\n" );
  printf ( "   i      U          D12(U)    D12inv(D12(U))\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %2d  %10g  %10g  %10g\n", i, u[i], v[i], w[i] );
  }

  free ( u );
  free ( v );
  free ( w );

  return;
}
/******************************************************************************/

void test07 ( )

/******************************************************************************/
/*
  Purpose:

    TEST07 tests DAUB14_TRANSFORM and DAUB14_TRANSFORM_INVERSE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 May 2012

  Author:

    John Burkardt
*/
{
  double a_first;
  double a_last;
  int i;
  int n;
  int seed;
  double *u;
  double *v;
  double *w;

  printf ( "\n" );
  printf ( "TEST07\n" );
  printf ( "  DAUB14_TRANSFORM computes the DAUB14 transform of a vector.\n" );
  printf ( "  DAUB14_TRANSFORM_INVERSE inverts it.\n" );
/*
  Random data.
*/
  n = 16;
  seed = 123456789;

  u = r8vec_uniform_01_new ( n, &seed );

  v = daub14_transform ( n, u );

  w = daub14_transform_inverse ( n, v );

  printf ( "\n" );
  printf ( "   i      U          D14(U)    D14inv(D14(U))\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %2d  %10g  %10g  %10g\n", i, u[i], v[i], w[i] );
  }

  free ( u );
  free ( v );
  free ( w );
/*
  Constant signal.
*/
  n = 8;
  u = ( double * ) malloc ( n * sizeof ( double ) );
  for ( i = 0; i < n; i++ )
  {
    u[i] = 1.0;
  }

  v = daub14_transform ( n, u );

  w = daub14_transform_inverse ( n, v );

  printf ( "\n" );
  printf ( "   i      U          D14(U)    D14inv(D14(U))\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %2d  %10g  %10g  %10g\n", i, u[i], v[i], w[i] );
  }

  free ( u );
  free ( v );
  free ( w );
/*
  Linear signal.
*/
  n = 16;
  a_first = 1.0;
  a_last = ( double ) ( n );
  u = r8vec_linspace_new ( n, a_first, a_last );

  v = daub14_transform ( n, u );

  w = daub14_transform_inverse ( n, v );

  printf ( "\n" );
  printf ( "   i      U          D14(U)    D14inv(D14(U))\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %2d  %10g  %10g  %10g\n", i, u[i], v[i], w[i] );
  }

  free ( u );
  free ( v );
  free ( w );
/*
  Quadratic data.
*/
  n = 8;
  u = ( double * ) malloc ( n * sizeof ( double ) );
  for ( i = 0; i < n; i++ )
  {
    u[i] = pow ( ( double ) ( i - 5 ), 2 );
  }

  v = daub14_transform ( n, u );

  w = daub14_transform_inverse ( n, v );

  printf ( "\n" );
  printf ( "   i      U          D14(U)    D14inv(D14(U))\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %2d  %10g  %10g  %10g\n", i, u[i], v[i], w[i] );
  }

  free ( u );
  free ( v );
  free ( w );

  return;
}
/******************************************************************************/

void test08 ( )

/******************************************************************************/
/*
  Purpose:

    TEST08 tests DAUB16_TRANSFORM and DAUB16_TRANSFORM_INVERSE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 May 2012

  Author:

    John Burkardt
*/
{
  double a_first;
  double a_last;
  int i;
  int n;
  int seed;
  double *u;
  double *v;
  double *w;

  printf ( "\n" );
  printf ( "TEST08\n" );
  printf ( "  DAUB16_TRANSFORM computes the DAUB16 transform of a vector.\n" );
  printf ( "  DAUB16_TRANSFORM_INVERSE inverts it.\n" );
/*
  Random data.
*/
  n = 16;
  seed = 123456789;

  u = r8vec_uniform_01_new ( n, &seed );

  v = daub16_transform ( n, u );

  w = daub16_transform_inverse ( n, v );

  printf ( "\n" );
  printf ( "   i      U          D16(U)    D16inv(D16(U))\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %2d  %10g  %10g  %10g\n", i, u[i], v[i], w[i] );
  }

  free ( u );
  free ( v );
  free ( w );
/*
  Constant signal.
*/
  n = 8;
  u = ( double * ) malloc ( n * sizeof ( double ) );
  for ( i = 0; i < n; i++ )
  {
    u[i] = 1.0;
  }

  v = daub16_transform ( n, u );

  w = daub16_transform_inverse ( n, v );

  printf ( "\n" );
  printf ( "   i      U          D16(U)    D16inv(D16(U))\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %2d  %10g  %10g  %10g\n", i, u[i], v[i], w[i] );
  }

  free ( u );
  free ( v );
  free ( w );
/*
  Linear signal.
*/
  n = 16;
  a_first = 1.0;
  a_last = ( double ) ( n );
  u = r8vec_linspace_new ( n, a_first, a_last );

  v = daub16_transform ( n, u );

  w = daub16_transform_inverse ( n, v );

  printf ( "\n" );
  printf ( "   i      U          D16(U)    D16inv(D16(U))\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %2d  %10g  %10g  %10g\n", i, u[i], v[i], w[i] );
  }

  free ( u );
  free ( v );
  free ( w );
/*
  Quadratic data.
*/
  n = 8;
  u = ( double * ) malloc ( n * sizeof ( double ) );
  for ( i = 0; i < n; i++ )
  {
    u[i] = pow ( ( double ) ( i - 5 ), 2 );
  }

  v = daub16_transform ( n, u );

  w = daub16_transform_inverse ( n, v );

  printf ( "\n" );
  printf ( "   i      U          D16(U)    D16inv(D16(U))\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %2d  %10g  %10g  %10g\n", i, u[i], v[i], w[i] );
  }

  free ( u );
  free ( v );
  free ( w );

  return;
}
/******************************************************************************/

void test09 ( )

/******************************************************************************/
/*
  Purpose:

    TEST09 tests DAUB18_TRANSFORM and DAUB18_TRANSFORM_INVERSE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 May 2012

  Author:

    John Burkardt
*/
{
  double a_first;
  double a_last;
  int i;
  int n;
  int seed;
  double *u;
  double *v;
  double *w;

  printf ( "\n" );
  printf ( "TEST09\n" );
  printf ( "  DAUB18_TRANSFORM computes the DAUB18 transform of a vector.\n" );
  printf ( "  DAUB18_TRANSFORM_INVERSE inverts it.\n" );
/*
  Random data.
*/
  n = 16;
  seed = 123456789;

  u = r8vec_uniform_01_new ( n, &seed );

  v = daub18_transform ( n, u );

  w = daub18_transform_inverse ( n, v );

  printf ( "\n" );
  printf ( "   i      U          D18(U)    D18inv(D18(U))\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %2d  %10g  %10g  %10g\n", i, u[i], v[i], w[i] );
  }

  free ( u );
  free ( v );
  free ( w );
/*
  Constant signal.
*/
  n = 8;
  u = ( double * ) malloc ( n * sizeof ( double ) );
  for ( i = 0; i < n; i++ )
  {
    u[i] = 1.0;
  }

  v = daub18_transform ( n, u );

  w = daub18_transform_inverse ( n, v );

  printf ( "\n" );
  printf ( "   i      U          D18(U)    D18inv(D18(U))\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %2d  %10g  %10g  %10g\n", i, u[i], v[i], w[i] );
  }

  free ( u );
  free ( v );
  free ( w );
/*
  Linear signal.
*/
  n = 16;
  a_first = 1.0;
  a_last = ( double ) ( n );
  u = r8vec_linspace_new ( n, a_first, a_last );

  v = daub18_transform ( n, u );

  w = daub18_transform_inverse ( n, v );

  printf ( "\n" );
  printf ( "   i      U          D18(U)    D18inv(D18(U))\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %2d  %10g  %10g  %10g\n", i, u[i], v[i], w[i] );
  }

  free ( u );
  free ( v );
  free ( w );
/*
  Quadratic data.
*/
  n = 8;
  u = ( double * ) malloc ( n * sizeof ( double ) );
  for ( i = 0; i < n; i++ )
  {
    u[i] = pow ( ( double ) ( i - 5 ), 2 );
  }

  v = daub18_transform ( n, u );

  w = daub18_transform_inverse ( n, v );

  printf ( "\n" );
  printf ( "   i      U          D18(U)    D18inv(D18(U))\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %2d  %10g  %10g  %10g\n", i, u[i], v[i], w[i] );
  }

  free ( u );
  free ( v );
  free ( w );

  return;
}
/******************************************************************************/

void test10 ( )

/******************************************************************************/
/*
  Purpose:

    TEST10 tests DAUB20_TRANSFORM and DAUB20_TRANSFORM_INVERSE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 May 2012

  Author:

    John Burkardt
*/
{
  double a_first;
  double a_last;
  int i;
  int n;
  int seed;
  double *u;
  double *v;
  double *w;

  printf ( "\n" );
  printf ( "TEST10\n" );
  printf ( "  DAUB20_TRANSFORM computes the DAUB20 transform of a vector.\n" );
  printf ( "  DAUB20_TRANSFORM_INVERSE inverts it.\n" );
/*
  Random data.
*/
  n = 16;
  seed = 123456789;

  u = r8vec_uniform_01_new ( n, &seed );

  v = daub20_transform ( n, u );

  w = daub20_transform_inverse ( n, v );

  printf ( "\n" );
  printf ( "   i      U          D20(U)    D20inv(D20(U))\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %2d  %10g  %10g  %10g\n", i, u[i], v[i], w[i] );
  }

  free ( u );
  free ( v );
  free ( w );
/*
  Constant signal.
*/
  n = 8;
  u = ( double * ) malloc ( n * sizeof ( double ) );
  for ( i = 0; i < n; i++ )
  {
    u[i] = 1.0;
  }

  v = daub20_transform ( n, u );

  w = daub20_transform_inverse ( n, v );

  printf ( "\n" );
  printf ( "   i      U          D20(U)    D20inv(D20(U))\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %2d  %10g  %10g  %10g\n", i, u[i], v[i], w[i] );
  }

  free ( u );
  free ( v );
  free ( w );
/*
  Linear signal.
*/
  n = 16;
  a_first = 1.0;
  a_last = ( double ) ( n );
  u = r8vec_linspace_new ( n, a_first, a_last );

  v = daub20_transform ( n, u );

  w = daub20_transform_inverse ( n, v );

  printf ( "\n" );
  printf ( "   i      U          D20(U)    D20inv(D20(U))\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %2d  %10g  %10g  %10g\n", i, u[i], v[i], w[i] );
  }

  free ( u );
  free ( v );
  free ( w );
/*
  Quadratic data.
*/
  n = 8;
  u = ( double * ) malloc ( n * sizeof ( double ) );
  for ( i = 0; i < n; i++ )
  {
    u[i] = pow ( ( double ) ( i - 5 ), 2 );
  }

  v = daub20_transform ( n, u );

  w = daub20_transform_inverse ( n, v );

  printf ( "\n" );
  printf ( "   i      U          D20(U)    D20inv(D20(U))\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %2d  %10g  %10g  %10g\n", i, u[i], v[i], w[i] );
  }

  free ( u );
  free ( v );
  free ( w );

  return;
}
