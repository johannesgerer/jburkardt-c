# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <string.h>

# include "differ.h"

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );
void test04 ( );
void test05 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for DIFFER_PRB.

  Discussion:

    DIFFER_PRB tests the DIFFER library.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    10 November 2013

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "DIFFER_PRB:\n" );
  printf ( "  C version\n" );
  printf ( "  Test the DIFFER library.\n" );

  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
  test05 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "DIFFER_PRB:\n" );
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

    TEST01 tests DIFFER_MATRIX.

  Discussion:

    DIFFER_MATRIX computes a modified Vandermonde matrix A1.

    The solution of a system A1 * X1 = B is related to the solution
    of the system A2 * X2 = B, where A2 is the standard Vandermonde
    matrix, simply by X2(I) = X1(I) * A(I,1).

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    03 November 2013

  Author:

    John Burkardt
*/
{
  double *a;
  double *b;
  int i;
  int info;
  int job;
  int n = 4;
  double stencil[4] = { 2.5, 3.3, -1.3, 0.5 };
  double x[4] = { 1.0, 2.0, 3.0, 4.0 };
  double *x1;
  double *x2;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  Demonstrate that the DIFFER matrix is 'really'\n" );
  printf ( "  a Vandermonde matrix.\n" );

  a = differ_matrix ( n, stencil );
  r8mat_print ( n, n, a, "  Stencil matrix:" );
  b = r8mat_mv_new ( n, n, a, x );
/*
  Set up and solve the DIFFER system.
*/
  a = differ_matrix ( n, stencil );
  x1 = r8mat_fs_new ( n, a, b );

  r8vec_print ( n, x1, "  Solution of DIFFER system:" );
/*
  R8VM_SL_NEW solves the related Vandermonde system.
  A simple transformation gives us the solution to the DIFFER system.
*/
  job = 0;
  x2 = r8vm_sl_new ( n, stencil, b, job, &info );

  if ( info != 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "TEST01 - Warning!\n" );
    fprintf ( stderr, "  VANDERMONDE system is singular.\n" );
    exit ( 1 );
  }

  r8vec_print ( n, x2, "  Solution of VANDERMONDE system:" );

  for ( i = 0; i < n; i++ )
  {
    x2[i] = x2[i] / stencil[i];
  }
  r8vec_print ( n, x2, "  Transformed solution of VANDERMONDE system:" );

  free ( a );
  free ( b );
  free ( x1 );
  free ( x2 );

  return;
}
/******************************************************************************/

void test02 ( )

/******************************************************************************/
/*
  Purpose:

    TEST02 tests DIFFER_INVERSE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    03 November 2013

  Author:

    John Burkardt
*/
{
  double *a;
  double *b;
  double err;
  int n;
  int n_max = 8;
  int seed;
  int test;
  double *x;

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  DIFFER_INVERSE returns the inverse of a DIFFER matrix;\n" );
  printf ( "\n" );
  printf ( "   N    Inverse error\n" );

  seed = 123456789;

  for ( n = 2; n <= n_max; n++ )
  {
    printf ( "\n" );

    for ( test = 1; test <= 5; test++ )
    {
      x = r8vec_uniform_01_new ( n, &seed );
      a = differ_matrix ( n, x );
      b = differ_inverse ( n, x );
      err = inverse_error ( n, a, b );
      printf ( "  %2d  %14.6g\n", n, err );
      free ( a );
      free ( b );
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

    TEST03 tests DIFFER_MATRIX.

  Discussion:

    Reproduce a specific example.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    28 October 2013

  Author:

    John Burkardt
*/
{
  double *a;
  double *b;
  double *c;
  double df;
  double dfdx;
  double dx;
  int i;
  int n = 4;
  int order;
  double stencil[4] = { -3.0, -2.0, -1.0, 1.0 };
  double x0;

  printf ( "\n" );
  printf ( "TEST03\n" );
  printf ( "  Reproduce a specific example.\n" );
/*
  Compute the coefficients for a specific stencil.
*/
  b = ( double * ) malloc ( n * sizeof ( double ) );
  for ( i = 0; i < n; i++ )
  {
    b[i] = 0.0;
  }
  order = 1;
  b[order-1] = 1.0;
  a = differ_matrix ( n, stencil );
  c = r8mat_fs_new ( n, a, b );

  r8vec_print ( n, c, "  Solution of DIFFER system:" );
/*
  Use the coefficients C to estimate the first derivative of EXP(X)
  at X0, using a spacing of DX = 0.1.
*/
  x0 = 1.3;
  dx = 0.1;
  df = 0.0;
  for ( i = 0; i < n; i++ )
  {
    df = df + c[i] * ( exp ( x0 + stencil[i] * dx ) - exp ( x0 ) );
  }
  dfdx = df / dx;

  printf ( "\n" );
  printf ( "  DFDX =         %g\n", dfdx );
  printf ( "  d exp(x) /dx = %g\n", exp ( x0 ) );

  free ( a );
  free ( b );
  free ( c );

  return;
}
/******************************************************************************/

void test04 ( )

/******************************************************************************/
/*
  Purpose:

    TEST04 tests DIFFER_FORWARD, DIFFER_BACKWARD, DIFFER_CENTRAL.

  Discussion:

    Evaluate the coefficients for uniformly spaced finite difference
    approximations of derivatives.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    10 November 2013

  Author:

    John Burkardt
*/
{
  double *c;
  double h;
  char label[80];
  int n;
  int o;
  int p;
  double *x;

  printf ( "\n" );
  printf ( "TEST04\n" );
  printf ( "  DIFFER_FORWARD,\n" );
  printf ( "  DIFFER_BACKWARD, and\n" );
  printf ( "  DIFFER_CENTRAL produce coefficients for difference\n" );
  printf ( "  approximations of the O-th derivative,\n" );
  printf ( "  with error of order H^P, for a uniform spacing of H.\n" );

  h = 1.0;
  printf ( "\n" );
  printf ( "  Use a spacing of H = %g for all examples.\n", h );
/*
  Forward difference approximation to the third derivative with error of O(h).
*/
  o = 3;
  p = 1;
  n = o + p;
  c = ( double * ) malloc ( n * sizeof ( double ) );
  x = ( double * ) malloc ( n * sizeof ( double ) );
  differ_forward ( h, o, p, c, x );
  sprintf ( label, "  Forward difference coefficients, O = %d P = %d", o, p );
  r8vec2_print ( n, x, c, label );
  free ( c );
  free ( x );
/*
  Backward difference approximation to the third derivative with error of O(h).
*/
  o = 3;
  p = 1;
  n = o + p;
  c = ( double * ) malloc ( n * sizeof ( double ) );
  x = ( double * ) malloc ( n * sizeof ( double ) );
  differ_backward ( h, o, p, c, x );
  sprintf ( label, "  Backward difference coefficients, O = %d P = %d", o, p );
  r8vec2_print ( n, x, c, label );
  free ( c );
  free ( x );
/*
  Central difference approximation to the third derivative with error of O(h^2).
*/
  o = 3;
  p = 2;
  n = o + p;
  c = ( double * ) malloc ( n * sizeof ( double ) );
  x = ( double * ) malloc ( n * sizeof ( double ) );
  differ_central ( h, o, p, c, x );
  sprintf ( label, "  Central difference coefficients, O = %d P = %d", o, p );
  r8vec2_print ( n, x, c, label );
  free ( c );
  free ( x );
/*
  Central difference approximation to the third derivative with error of O(h^4).
*/
  o = 3;
  p = 4;
  n = o + p;
  c = ( double * ) malloc ( n * sizeof ( double ) );
  x = ( double * ) malloc ( n * sizeof ( double ) );
  differ_central ( h, o, p, c, x );
  sprintf ( label, "  Central difference coefficients, O = %d P = %d", o, p );
  r8vec2_print ( n, x, c, label );
  free ( c );
  free ( x );
/*
  Forward difference approximation to the fourth derivative with error of O(h).
*/
  o = 4;
  p = 1;
  n = o + p;
  c = ( double * ) malloc ( n * sizeof ( double ) );
  x = ( double * ) malloc ( n * sizeof ( double ) );
  differ_forward ( h, o, p, c, x );
  sprintf ( label, "  Forward difference coefficients, O = %d P = %d", o, p );
  r8vec2_print ( n, x, c, label );
  free ( c );
  free ( x );
/*
  Backward difference approximation to the fourth derivative with error of O(h).
*/
  o = 4;
  p = 1;
  n = o + p;
  c = ( double * ) malloc ( n * sizeof ( double ) );
  x = ( double * ) malloc ( n * sizeof ( double ) );
  differ_backward ( h, o, p, c, x );
  sprintf ( label, "  Backward difference coefficients, O = %d P = %d", o, p );
  r8vec2_print ( n, x, c, label );
  free ( c );
  free ( x );
/*
   Central difference approximation to the fourth derivative with error of O(h^3).
*/
  o = 4;
  p = 3;
  n = o + p;
  c = ( double * ) malloc ( n * sizeof ( double ) );
  x = ( double * ) malloc ( n * sizeof ( double ) );
  differ_central ( h, o, p, c, x );
  sprintf ( label, "  Central difference coefficients, O = %d P = %d", o, p );
  r8vec2_print ( n, x, c, label );
  free ( c );
  free ( x );

  return;
}
/******************************************************************************/

void test05 ( )

/******************************************************************************/
/*
  Purpose:

    TEST05 tests DIFFER_STENCIL.

  Discussion:

    Evaluate the coefficients for uniformly spaced finite difference
    approximations of derivatives.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    10 November 2013

  Author:

    John Burkardt
*/
{
  double *c;
  double h;
  int i;
  char label[80];
  int n;
  int o;
  int p;
  double *x;
  double x0;

  printf ( "\n" );
  printf ( "TEST05\n" );
  printf ( "  DIFFER_STENCIL produces coefficients for difference\n" );
  printf ( "  approximations of the O-th derivative,\n" );
  printf ( "  using arbitrarily spaced data, with maximum spacing H\n" );
  printf ( "  with error of order H^P.\n" );
/*
  Let X0 = 1.0.
*/
  x0 = 0.0;
  h = 1.0;
  printf ( "\n" );
  printf ( "  Use a spacing of H = %g for all examples.\n", h );
/*
  Forward difference approximation to the third derivative with error of O(h).
*/
  o = 3;
  p = 1;
  n = o + p;
  c = ( double * ) malloc ( n * sizeof ( double ) );
  x = ( double * ) malloc ( n * sizeof ( double ) );
  for ( i = 0; i < n; i++ )
  {
    x[i] = ( double ) ( i ) * h;
  }
  differ_stencil ( x0, o, p, x, c );
  sprintf ( label, "  Forward difference coefficients, O = %d P = %d", o, p );
  r8vec2_print ( n, x, c, label );
  free ( c );
  free ( x );
/*
  Backward difference approximation to the third derivative with error of O(h).
*/
  o = 3;
  p = 1;
  n = o + p;
  c = ( double * ) malloc ( n * sizeof ( double ) );
  x = ( double * ) malloc ( n * sizeof ( double ) );
  for ( i = 0; i < n; i++ )
  {
    x[i] = ( double ) ( i + 1 - n ) * h;
  }
  differ_stencil ( x0, o, p, x, c );
  sprintf ( label, "  Backward difference coefficients, O = %d P = %d", o, p );
  r8vec2_print ( n, x, c, label );
  free ( c );
  free ( x );
/*
  Central difference approximation to the third derivative with error of O(h^2).
*/
  o = 3;
  p = 2;
  n = o + p;
  c = ( double * ) malloc ( n * sizeof ( double ) );
  x = ( double * ) malloc ( n * sizeof ( double ) );
  for ( i = 0; i < n; i++ )
  {
    x[i] = ( double ) ( - n + 1 + 2 * i ) * h / 2.0;
  }
  differ_stencil ( x0, o, p, x, c );
  sprintf ( label, "  Central difference coefficients, O = %d P = %d", o, p );
  r8vec2_print ( n, x, c, label );
  free ( c );
  free ( x );
/*
  Central difference approximation to the third derivative with error of O(h^4).
*/
  o = 3;
  p = 4;
  n = o + p;
  c = ( double * ) malloc ( n * sizeof ( double ) );
  x = ( double * ) malloc ( n * sizeof ( double ) );
  for ( i = 0; i < n; i++ )
  {
    x[i] = ( double ) ( - n + 1 + 2 * i ) * h / 2.0;
  }
  differ_stencil ( x0, o, p, x, c );
  sprintf ( label, "  Central difference coefficients, O = %d P = %d", o, p );
  r8vec2_print ( n, x, c, label );
  free ( c );
  free ( x );
/*
  Forward difference approximation to the fourth derivative with error of O(h).
*/
  o = 4;
  p = 1;
  n = o + p;
  c = ( double * ) malloc ( n * sizeof ( double ) );
  x = ( double * ) malloc ( n * sizeof ( double ) );
  for ( i = 0; i < n; i++ )
  {
    x[i] = ( double ) ( i ) * h;
  }
  differ_stencil ( x0, o, p, x, c );
  sprintf ( label, "  Forward difference coefficients, O = %d P = %d", o, p );
  r8vec2_print ( n, x, c, label );
  free ( c );
  free ( x );
/*
  Backward difference approximation to the fourth derivative with error of O(h).
*/
  o = 4;
  p = 1;
  n = o + p;
  c = ( double * ) malloc ( n * sizeof ( double ) );
  x = ( double * ) malloc ( n * sizeof ( double ) );
  for ( i = 0; i < n; i++ )
  {
    x[i] = ( double ) ( i + 1 - n ) * h;
  }
  differ_stencil ( x0, o, p, x, c );
  sprintf ( label, "  Backward difference coefficients, O = %d P = %d", o, p );
  r8vec2_print ( n, x, c, label );
  free ( c );
  free ( x );
/*
   Central difference approximation to the fourth derivative with error of O(h^3).
*/
  o = 4;
  p = 3;
  n = o + p;
  c = ( double * ) malloc ( n * sizeof ( double ) );
  x = ( double * ) malloc ( n * sizeof ( double ) );
  for ( i = 0; i < n; i++ )
  {
    x[i] = ( double ) ( - n + 1 + 2 * i ) * h / 2.0;
  }
  differ_stencil ( x0, o, p, x, c );
  sprintf ( label, "  Central difference coefficients, O = %d P = %d", o, p );
  r8vec2_print ( n, x, c, label );
  free ( c );
  free ( x );

  return;
}
