# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "bernstein.h"

int main ( void );
void test01 ( void );
void test02 ( void );
void test03 ( void );
void test04 ( void );
void test05 ( void );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for BERNSTEIN_PRB.

  Discussion:

    BERNSTEIN_PRB calls the BERNSTEIN test routines.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 February 2012

  Author:

    John Burkardt
*/
{
  timestamp ( );

  printf ( "\n" );
  printf ( "BERNSTEIN_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the BERNSTEIN library.\n" );

  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
  test05 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "BERNSTEIN_PRB\n" );
  printf ( "  Normal end of execution.\n" );

  printf ( " \n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void test01 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST01 tests BERNSTEIN_POLY and BERNSTEIN_POLY_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 February 2012

  Author:

    John Burkardt
*/
{
  double b;
  double *bvec;
  int k;
  int n;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST01:\n" );
  printf ( "  BERNSTEIN_POLY evaluates the Bernstein polynomials\n" );
  printf ( "  based on the interval [0,1].\n" );
  printf ( "  BERNSTEIN_POLY_VALUES returns some exact values.\n" );
  printf ( "\n" );
  printf ( "     N     K     X       Exact         BP01(N,K)(X)\n" );
  printf ( "\n" );

  n_data = 0;

  while ( 1 )
  {
    bernstein_poly_values ( &n_data, &n, &k, &x, &b );

    if ( n_data == 0 )
    {
      break;
    }

    bvec = bernstein_poly ( n, x );

    printf ( "  %4d  %4d  %7f  %14g  %14g\n", n, k, x, b, bvec[k] );

    free ( bvec );
  }

  return;
}
/******************************************************************************/

void test02 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST02 tests BPAB.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 February 2012

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double *bern;
  int k;
  int n = 10;
  double x;

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  BPAB evaluates Bernstein polynomials over an\n" );
  printf ( "  arbitrary interval [A,B].\n" );
  printf ( "\n" );
  printf ( "  Here, we demonstrate that \n" );
  printf ( "    BPAB(N,K,A1,B1)(X1) = BPAB(N,K,A2,B2)(X2)\n" );
  printf ( "  provided only that\n" );
  printf ( "    (X1-A1)/(B1-A1) = (X2-A2)/(B2-A2).\n" );

  x = 0.3;
  a = 0.0;
  b = 1.0;
  bern = bpab ( n, a, b, x );
 
  printf ( "\n" );
  printf ( "     N     K     A        B        X       BPAB(N,K,A,B)(X)\n" );
  printf ( "\n" );
  for ( k = 0; k <= n; k++ )
  {
    printf ( "  %4d  %4d  %7f  %7f  %7f  %14g\n", n, k, a, b, x, bern[k] );
  }

  free ( bern );
 
  x = 1.3;
  a = 1.0;
  b = 2.0;
  bern = bpab ( n, a, b, x );
 
  printf ( "\n" );
  printf ( "     N     K     A        B        X       BPAB(N,K,A,B)(X)\n" );
  printf ( "\n" ); 
  for ( k = 0; k <= n; k++ )
  {
    printf ( "  %4d  %4d  %7f  %7f  %7f  %14g\n", n, k, a, b, x, bern[k] );
  }

  free ( bern );

  x = 2.6;
  a = 2.0;
  b = 4.0;
  bern = bpab ( n, a, b, x );
 
  printf ( "\n" );
  printf ( "     N     K     A        B        X       BPAB(N,K,A,B)(X)\n" );
  printf ( "\n" );
 
  for ( k = 0; k <= n; k++ )
  {
    printf ( "  %4d  %4d  %7f  %7f  %7f  %14g\n", n, k, a, b, x, bern[k] );
  }

  free ( bern );

  return;
}
/******************************************************************************/

void test03 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST03 tests the Partition-of-Unity property.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 February 2012

  Author:

    John Burkardt
*/
{
  double *bvec;
  int n;
  int n_data;
  int seed;
  double x;

  printf ( "\n" );
  printf ( "TEST03:\n" );
  printf ( "  BERNSTEIN_POLY evaluates the Bernstein polynomials\n" );
  printf ( "  based on the interval [0,1].\n" );
  printf ( "\n" );
  printf ( "  Here we test the partition of unity property.\n" );
  printf ( "\n" );
  printf ( "     N     X          Sum ( 0 <= K <= N ) BP01(N,K)(X)\n" );
  printf ( "\n" );

  seed = 123456789;

  for ( n = 0; n <= 10; n++ )
  {
    x = r8_uniform_01 ( &seed );

    bvec = bernstein_poly ( n, x );

    printf ( "  %4d  %7f  %14g\n", n, x, r8vec_sum ( n + 1, bvec ) );

    free ( bvec );
  }
  return;
}
/******************************************************************************/

void test04 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST04 tests BPAB_APPROX.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 February 2012

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double error_max;
  int i;
  int maxdata = 20;
  int ndata;
  int nsample;
  int nval = 501;
  double *xdata;
  double *xval;
  double *ydata;
  double *yval;

  printf ( "\n" );
  printf ( "TEST04\n" );
  printf ( "  BPAB_APPROX evaluates the Bernstein polynomial\n" );
  printf ( "  approximant to a function F(X).\n" );

  a = 1.0;
  b = 3.0;

  printf ( "\n" );
  printf ( "     N      Max Error\n" );
  printf ( "\n" );

  for ( ndata = 0; ndata <= maxdata; ndata++ )
  {
/*
  Generate data values.
*/
    xdata = ( double * ) malloc ( ( ndata + 1 ) * sizeof ( double ) );
    ydata = ( double * ) malloc ( ( ndata + 1 ) * sizeof ( double ) );
    for ( i = 0; i <= ndata; i++)
    {
      if ( ndata == 0 )
      {
        xdata[i] = 0.5 * ( a + b );
      }
      else
      {
        xdata[i] = ( ( double ) ( ndata - i ) * a   
                   + ( double ) (         i ) * b ) 
                   / ( double ) ( ndata     );
      }
      ydata[i] = sin ( xdata[i] );
    }
/*
  Compare the true function and the approximant.
*/
    xval = r8vec_linspace_new ( nval, a, b );

    error_max = 0.0;

    yval = bpab_approx ( ndata, a, b, ydata, nval, xval );

    error_max = 0.0;
    for ( i = 0; i < nval; i++ )
    {
      error_max = r8_max ( error_max, r8_abs ( yval[i] - sin ( xval[i] ) ) );
    }
    printf ( "  %4d  %14g\n", ndata, error_max );

    free ( xdata );
    free ( xval );
    free ( ydata );
    free ( yval );
  }
  return;
}
/******************************************************************************/

void test05 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST05 tests BERNSTEIN_MATRIX and BERNSTEIN_MATRIX_INVERSE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 February 2012

  Author:

    John Burkardt
*/
{
  double *a;
  double a_norm_frobenius;
  double *b;
  double b_norm_frobenius;
  double *c;
  double error_norm_frobenius;
  int n;

  printf ( "\n" );
  printf ( "TEST05\n" );
  printf ( "  BERNSTEIN_MATRIX returns a matrix A which transforms a\n" );
  printf ( "  polynomial coefficient vector from the power basis to\n" );
  printf ( "  the Bernstein basis.\n" );
  printf ( "  BERNSTEIN_MATRIX_INVERSE computes the inverse B.\n" );
  printf ( "\n" );
  printf ( "     N     ||A||            ||B||      ||I-A*B||\n" );
  printf ( "\n" );

  for ( n = 5; n <= 15; n++ )
  {
    a = bernstein_matrix ( n );
    a_norm_frobenius = r8mat_norm_fro ( n, n, a );

    b = bernstein_matrix_inverse ( n );
    b_norm_frobenius = r8mat_norm_fro ( n, n, b );

    c = r8mat_mm_new ( n, n, n, a, b );
    error_norm_frobenius = r8mat_is_identity ( n, c );

    printf ( "  %4d  %14g  %14g  %14g\n", 
      n, a_norm_frobenius, b_norm_frobenius, error_norm_frobenius );

    free ( a );
    free ( b );
    free ( c );
  }
  return;
}

