# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "jacobi_polynomial.h"

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

    MAIN is the main program for JACOBI_POLYNOMIAL_PRB.

  Discussion:

    JACOBI_POLYNOMIAL_PRB tests the JACOBI_POLYNOMIAL library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 August 2013

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "JACOBI_POLYNOMIAL_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the JACOBI_POLYNOMIAL library.\n" );

  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "JACOBI_POLYNOMIAL_PRB\n" );
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

    TEST01 tests J_POLYNOMIAL.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 August 2013

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double e;
  double fx1;
  double fx2;
  double *fx2_vec;
  int m;
  int n;
  int n_data;
  double x;
  double x_vec[1];

  printf ( "\n" );
  printf ( "TEST01:\n" );
  printf ( "  J_POLYNOMIAL_VALUES stores values of\n" );
  printf ( "  the Jacobi polynomials.\n" );
  printf ( "  J_POLYNOMIAL evaluates the polynomial.\n" );
  printf ( "\n" );
  printf ( "                                    Tabulated                 Computed\n" );
  printf ( "     N     A     B        X           J(N,A,B,X)                    J(N,A,B,X)                     Error\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    j_polynomial_values ( &n_data, &n, &a, &b, &x, &fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    m = 1;
    x_vec[0] = x;
    fx2_vec = j_polynomial ( m, n, a, b, x_vec );
    fx2 = fx2_vec[0+n*1];
    e = fx1 - fx2;

    printf ( "  %4d  %6g  %6g  %6g  %24g  %24g  %8g\n",
      n, a, b, x, fx1, fx2, e );

    free ( fx2_vec );
  }
  return;
}
/******************************************************************************/

void test02 ( )

/******************************************************************************/
/*
  Purpose:

    TEST02 tests J_POLYNOMIAL_ZEROS.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 August 2013

  Author:

    John Burkardt
*/
{
  double a;
  double a_test[3] = { 0.5, 1.0, 2.0 };
  double b;
  double b_test[3] = { 0.5, 1.5, 0.5 };
  int degree;
  double *hz;
  int test;
  int test_num = 3;
  char title[80];
  double *z;

  printf ( "\n" );
  printf ( "TEST02:\n" );
  printf ( "  J_POLYNOMIAL_ZEROS computes the zeros of J(n,a,b,x);\n" );
  printf ( "  Check by calling J_POLYNOMIAL there.\n" );

  for ( test = 0; test < test_num; test++ )
  {
    a = a_test[test];
    b = b_test[test];

    for ( degree = 1; degree <= 5; degree++ )
    {
      z = j_polynomial_zeros ( degree, a, b );
      sprintf ( title, "Zeros for J(%d,%f,%f)", degree, a, b );
      r8vec_print ( degree, z, title );

      hz = j_polynomial ( degree, degree, a, b, z );
      sprintf ( title, "Evaluate J(%d,%f,%f)", degree, a, b );
      r8vec_print ( degree, hz + degree * degree, title );

      free ( hz );
      free ( z );
    }
  }
  return;
}
/******************************************************************************/

void test03 ( )

/******************************************************************************/
/*
  Purpose:

    TEST03 tests J_QUADRATURE_RULE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 August 2013

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  int i;
  int j;
  double *ji;
  double *jj;
  int k;
  int n;
  double q;
  double q_exact;
  double *w;
  double *x;

  printf ( "\n" );
  printf ( "TEST03:\n" );
  printf ( "  J_QUADRATURE_RULE computes the quadrature rule\n" );
  printf ( "  associated with J(n,a,b,x);\n" );

  n = 7;
  a = 1.0;
  b = 2.5;

  x = ( double * ) malloc ( n * sizeof ( double ) );
  w = ( double * ) malloc ( n * sizeof ( double ) );

  j_quadrature_rule ( n, a, b, x, w );

  r8vec2_print ( n, x, w, "      X            W" );

  printf ( "\n" );
  printf ( "  Use the quadrature rule to estimate:\n" );
  printf ( "\n" );
  printf ( "    Q = Integral (-1<x<+1) J(i,a,b,x) J(j,a,b,x) (1-x)^a (1+x)^b dx\n" );
  printf ( "\n" );
  printf ( "   I   J      Q_Estimate         Q_Exact\n" );
  printf ( "\n" );

  for ( i = 0; i <= 5; i++ )
  {
    ji = j_polynomial ( n, i, a, b, x );
    for ( j = i; j <= 5; j++ )
    {
      jj = j_polynomial ( n, j, a, b, x );
      q = 0.0;
      for ( k = 0; k < n; k++ )
      {
        q = q + w[k] * ji[k+i*n] * jj[k+j*n];
      }
      q_exact = j_double_product_integral ( i, j, a, b );
      printf ( "  %2d  %2d  %14g  %14g\n", i, j, q, q_exact );
      free ( jj );
    }
    free ( ji );
  }
  free ( x );
  free ( w );

  return;
}
/******************************************************************************/

void test04 ( )

/******************************************************************************/
/*
  Purpose:

    TEST04 tests J_DOUBLE_PRODUCT_INTEGRAL.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 August 2013

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  int i;
  int j;
  double q;

  printf ( "\n" );
  printf ( "TEST04:\n" );
  printf ( "  J_DOUBLE_PRODUCT_INTEGRAL returns the weighted integral of\n" );
  printf ( "  J(i,a,b,x) * J(j,a,b,x);\n" );

  a = 1.0;
  b = 2.5;

  printf ( "\n" );
  printf ( "    Q = Integral (-1<x<+1) J(i,a,b,x) J(j,a,b,x) (1-x)^a (1+x)^b dx\n" );
  printf ( "\n" );
  printf ( "   I   J      Q\n" );
  printf ( "\n" );

  for ( i = 0; i <= 5; i++ )
  {
    for ( j = i; j <= 5; j++ )
    {
      q = j_double_product_integral ( i, j, a, b );
      printf ( "  %2d  %2d  %14g\n", i, j, q );
    }
  }
  return;
}

