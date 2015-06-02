# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <string.h>

# include "test_interp_nd.h"
# include "r8lib.h"

int main ( );
void test01 ( );
void test02 ( int m, int n );
void test03 ( int m, int n );
void test04 ( int m, int n );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for TEST_INTERP_ND_PRB.

  Discussion:

    TEST_INTERP_ND_PRB tests the TEST_INTERP_ND library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 August 2012

  Author:

    John Burkardt
*/
{
  int m;
  int n;

  timestamp (  );
  printf ( "\n" );
  printf ( "TEST_INTERP_ND_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the TEST_INTERP_ND library.\n" );
  printf ( "  The R8LIB library is also needed.\n" );

  test01 ( );

  n = 10;
  for ( m = 2; m <= 4; m++ )
  {
    test02 ( m, n );
  }

  n = 2;
  for ( m = 2; m <= 4; m++ )
  {
    test03 ( m, n );
  }

  n = 10000;
  m = 4;
  test04 ( m, n );
//
//  Terminate.
//
  printf ( "\n" );
  printf ( "TEST_INTERP_ND_PRB\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
//****************************************************************************80

void test01 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST01 prints the title of each test function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 August 2012
//
//  Author:
//
//    John Burkardt
//
{
  int prob;
  int prob_num;
  char *title;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  P00_TITLE returns the problem title.\n" );

  prob_num = p00_prob_num ( );
  printf ( "\n" );
  printf ( "  There are a total of %d problems.\n", prob_num );
  printf ( "\n" );

  for ( prob = 1; prob <= prob_num; prob++ )
  {
    title = p00_title ( prob );
    printf ( "  %d  \"%s\"\n", prob, title );
    free ( title );
  }

  return;
}
//****************************************************************************80

void test02 ( int m, int n )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 samples each function in M dimensions, at N points.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 August 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Input, int N, the number of evaluation points.
//
{
  double *c;
  double *f;
  int i;
  int j;
  int prob;
  int prob_num;
  int seed;
  double *w;
  double *x;

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  P00_F evaluates the function.\n" );
  printf ( "  Here, we use spatial dimension M = %d\n", m );
  printf ( "  The number of points is N = %d\n", n );

  seed = 123456789;
  x = r8mat_uniform_01_new ( m, n, &seed );

  prob_num = p00_prob_num ( );

  for ( prob = 1; prob <= prob_num; prob++ )
  {
    printf ( "\n" );
    printf ( "  Problem  %d\n", prob );

    c = p00_c ( prob, m, &seed );
    r8vec_print ( m, c, "  C parameters:" );

    w = p00_w ( prob, m, &seed );
    r8vec_print ( m, w, "  W parameters:" );

    printf ( "\n" );
    printf ( "      F(X)              X(1)      X(2) ...\n" );
    printf ( "\n" );

    f = p00_f ( prob, m, c, w, n, x );

    for ( j = 0; j < n; j++ )
    {
      printf ( "  %14g", f[j] );
      for ( i = 0; i < m; i++ )
      {
        printf ( "  %10f", x[i+j*m] );
      }
      printf ( "\n" );
    }
    free ( c );
    free ( f );
    free ( w );
  }
  free ( x );

  return;
}
//****************************************************************************80

void test03 ( int m, int n )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 samples each derivative component in M dimensions, at N points.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 August 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Input, int N, the number of evaluation points.
//
{
  double *c;
  double *d;
  double *dmat;
  double *f;
  int i;
  int id;
  int j;
  int prob;
  int prob_num;
  int seed;
  double *w;
  double *x;

  printf ( "\n" );
  printf ( "TEST03\n" );
  printf ( "  P00_D evaluates derivative components.\n" );
  printf ( "  Here, we use spatial dimension M = %d\n", m );
  printf ( "  The number of points is N = %d\n", n );

  seed = 123456789;
  x = r8mat_uniform_01_new ( m, n, &seed );

  prob_num = p00_prob_num ( );

  for ( prob = 1; prob <= prob_num; prob++ )
  {
    printf ( "\n" );
    printf ( "  Problem  %d\n", prob );

    c = p00_c ( prob, m, &seed );
    r8vec_print ( m, c, "  C parameters:" );

    w = p00_w ( prob, m, &seed );
    r8vec_print ( m, w, "  W parameters:" );

    printf ( "\n" );
    printf ( "                         X(1)      X(2) ...\n" );
    printf ( "      F(X)            dFdX(1)   dFdX(2) ...\n" );

    f = p00_f ( prob, m, c, w, n, x );

    dmat = ( double * ) malloc ( m * n * sizeof ( double ) );

    for ( id = 0; id < m; id++ )
    {
      d = p00_d ( prob, m, id, c, w, n, x );
      for ( j = 0; j < n; j++ )
      {
        dmat[id+j*m] = d[j];
      }
      free ( d );
    }

    for ( j = 0; j < n; j++ )
    {
      printf ( "\n" );
      printf ( "                " );
      for ( i = 0; i < m; i++ )
      {
        printf ( "  %10f", x[i+j*m] );
      }
      printf ( "\n" );
      printf ( "  %14g", f[j] );
      for ( i = 0; i < m; i++ )
      {
        printf ( "  %10g", dmat[i+j*m] );
      }
      printf ( "\n" );
    }
    free ( c );
    free ( dmat );
    free ( f );
    free ( w );
  }

  free ( x );

  return;
}
/******************************************************************************/

void test04 ( int m, int n )

/******************************************************************************/
/*
  Purpose:

    TEST04 estimates integrals in M dimensions, using N points.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 August 2012

  Author:

    John Burkardt

  Parameters:

    Input, int M, the spatial dimension.

    Input, int N, the number of evaluation points.
*/
{
  double *c;
  double *f;
  int j;
  int prob;
  int prob_num;
  double q1;
  double q2;
  int seed;
  double *w;
  double *x;

  printf ( "\n" );
  printf ( "TEST04\n" );
  printf ( "P00_Q returns the integral of F over [0,1]^m.\n" );
  printf ( "  Here, we use spatial dimension M = %d\n", m );
  printf ( "  The number of sample points is N = %d\n", n );

  seed = 123456789;
  x = r8mat_uniform_01_new ( m, n, &seed );

  prob_num = p00_prob_num ( );

  for ( prob = 1; prob <= prob_num; prob++ )
  {
    printf ( "\n" );
    printf ( "  Problem  %d\n", prob );

    c = p00_c ( prob, m, &seed );
    r8vec_print ( m, c, "  C parameters:" );

    w = p00_w ( prob, m, &seed );
    r8vec_print ( m, w, "  W parameters:" );

    printf ( "\n" );
    printf ( "      Exact Integral     Q\n" );
    printf ( "\n" );

    q1 = p00_q ( prob, m, c, w );

    f = p00_f ( prob, m, c, w, n, x );
    q2 = r8vec_sum ( n, f ) / ( double ) ( n );

    printf ( "  %14g  %14g\n", q1, q2 );

    free ( c );
    free ( f );
    free ( w );
  }

  free ( x );

  return;
}

