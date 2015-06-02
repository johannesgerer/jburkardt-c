# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "haar.h"

int main ( );
void test01 ( );
void test02 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for HAAR_PRB.

  Discussion:

    HAAR_PRB tests the HAAR library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    18 March 2011

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "HAAR_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the HAAR library.\n" );

  test01 ( );
  test02 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "HAAR_PRB\n" );
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

    TEST01 tests HAAR_1D.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    06 March 2014

  Author:

    John Burkardt
*/
{
  double err;
  int i;
  int n;
  int seed;
  double *u;
  double *v;
  double *w;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  HAAR_1D computes the Haar transform of a vector.\n" );
/*
  Random data.
*/
  n = 16;
  seed = 123456789;
  u = r8vec_uniform_01_new ( n, &seed );
  v = r8vec_copy_new ( n, u );

  haar_1d ( n, v );

  w = r8vec_copy_new ( n, v );
  haar_1d_inverse ( n, w );

  printf ( "\n" );
  printf ( "   i      U(i)        H(U)(i)  Hinv(H(U))(i)\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %2d  %10f  %10f  %10f\n", i, u[i], v[i], w[i] );
  }

  free ( u );
  free ( v );
  free ( w );
/*
  Constant signal.
*/
  n = 8;
  u = r8vec_ones_new ( n );
  v = r8vec_copy_new ( n, u );

  haar_1d ( n, v );

  w = r8vec_copy_new ( n, v );
  haar_1d_inverse ( n, w );

  printf ( "\n" );
  printf ( "   i      U(i)        H(U)(i)  Hinv(H(U))(i)\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %2d  %10f  %10f  %10f\n", i, u[i], v[i], w[i] );
  }
  free ( u );
  free ( v );
  free ( w );
/*
  Linear signal.
*/
  n = 16;
  u = r8vec_linspace_new ( n, 1.0, ( double ) n );
  v = r8vec_copy_new ( n, u );

  haar_1d ( n, v );

  w = r8vec_copy_new ( n, v );
  haar_1d_inverse ( n, w );

  printf ( "\n" );
  printf ( "   i      U(i)        H(U)(i)  Hinv(H(U))(i)\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %2d  %10f  %10f  %10f\n", i, u[i], v[i], w[i] );
  }

  free ( u );
  free ( v );
  free ( w );
/*
  Quadratic data.
*/
  n = 8;
  u = ( double * ) malloc ( n * sizeof ( double ) );
  u[0] = 25.0;
  u[1] = 16.0;
  u[2] = 9.0;
  u[3] = 4.0;
  u[4] = 1.0;
  u[5] = 0.0;
  u[6] = 1.0;
  u[7] = 4.0;
  v = r8vec_copy_new ( n, u );

  haar_1d ( n, v );

  w = r8vec_copy_new ( n, v );
  haar_1d_inverse ( n, w );

  printf ( "\n" );
  printf ( "   i      U(i)        H(U)(i)  Hinv(H(U))(i)\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %2d  %10f  %10f  %10f\n", i, u[i], v[i], w[i] );
  }

  free ( u );
  free ( v );
  free ( w );
/*
  N not a power of 2.
*/
  n = 99;
  seed = 123456789;
  u = r8vec_uniform_01_new ( n, &seed );

  v = r8vec_copy_new ( n, u );
  haar_1d ( n, v );

  w = r8vec_copy_new ( n, v );
  haar_1d_inverse ( n, w );

  err = r8vec_diff_norm ( n, u, w );

  printf ( "\n" );
  printf ( "  For N = %d, ||u-haar_1d_inverse(haar_1d(u))|| = %g\n", n, err );

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

    TEST02 tests HAAR_2D and HAAR_2D_INVERSE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    06 March 2014

  Author:

    John Burkardt
*/
{
  double err;
  int i;
  int j;
  int m = 16;
  int n = 4;
  int seed;
  double *u;
  double *v;
  double *w;

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  HAAR_2D computes the Haar transform of an array.\n" );
  printf ( "  HAAR_2D_INVERSE inverts the transform.\n" );
/*
  Demonstrate successful inversion.
*/
  seed = 123456789;
  u = r8mat_uniform_01_new ( m, n, &seed );

  r8mat_print ( m, n, u, "  Input array U:" );

  v = r8mat_copy_new ( m, n, u );

  haar_2d ( m, n, v );

  r8mat_print ( m, n, v, "  Transformed array V:" );

  w = r8mat_copy_new ( m, n, v );

  haar_2d_inverse ( m, n, w );

  r8mat_print ( m, n, w, "  Recovered array W:" );

  free ( u );
  free ( v );
  free ( w );
/*
  M, N not powers of 2.
*/
  m = 37;
  n = 53;
  seed = 123456789;
  u = r8mat_uniform_01_new ( m, n, &seed );

  v = r8mat_copy_new ( m, n, u );
  haar_2d ( m, n, v );

  w = r8mat_copy_new ( m, n, v );
  haar_2d_inverse ( m, n, w );

  err = r8mat_dif_fro ( m, n, u, w );

  printf ( "\n" );
  printf ( "  M = %d, N = %d, ||haar_2d_inverse(haar_2d(u))-u|| = %g\n", m, n, err );

  free ( u );
  free ( v );
  free ( w );

  return;
}
