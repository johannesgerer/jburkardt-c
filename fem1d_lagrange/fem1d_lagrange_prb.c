# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "fem1d_lagrange.h"

int main ( );
void legendre_set_test ( );
void lagrange_value_test ( );
void lagrange_derivative_test ( );
void fem1d_lagrange_stiffness_test ( int x_num, int q_num );
double f ( double x );
double *exact ( int x_num, double x[] );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for FEM1D_LAGRANGE_PRB.

  Discussion:

    FEM1D_LAGRANGE_PRB tests FEM1D_LAGRANGE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 November 2014

  Author:

    John Burkardt
*/
{
  int q_num;
  int x_num;

  timestamp ( );
  printf ( "\n" );
  printf ( "FEM1D_LAGRANGE_PRB\n" );
  printf ( "  C version.\n" );
  printf ( "  Test the FEM1D_LAGRANGE library.\n" );

  legendre_set_test ( );
  lagrange_value_test ( );
  lagrange_derivative_test ( );

  x_num = 11;
  q_num = 5;
  fem1d_lagrange_stiffness_test ( x_num, q_num );

  x_num = 11;
  q_num = 10;
  fem1d_lagrange_stiffness_test ( x_num, q_num );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "FEM1D_LAGRANGE_PRB\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void legendre_set_test ( )

/******************************************************************************/
/*
  Purpose:

    LEGENDRE_SET_TEST tests LEGENDRE_SET.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 November 2014

  Author:

    John Burkardt
*/
{
  double e1;
  double e2;
  double e3;
  int i;
  int n;
  double *w;
  double *x;

  printf ( "\n" );
  printf ( "LEGENDRE_SET_TEST\n" );
  printf ( "  LEGENDRE_SET returns points and weights of\n" );
  printf ( "  Gauss-Legendre quadrature rules.\n" );
  printf ( "\n" );
  printf ( "   N               1             X^4           Runge\n" );
  printf ( "\n" );

  for ( n = 1; n <= 10; n++ )
  {
    x = ( double * ) malloc ( n * sizeof ( double ) );
    w = ( double * ) malloc ( n * sizeof ( double ) );

    legendre_set ( n, x, w );
    e1 = 0.0;
    e2 = 0.0;
    e3 = 0.0;
    for ( i = 0; i < n; i++ )
    {
      e1 = e1 + w[i];
      e2 = e2 + w[i] * pow ( x[i], 4 );
      e3 = e3 + w[i] / ( 1.0 + 25.0 * x[i] * x[i] );
    }
    printf ( "  %2d  %14.6g  %14.6g  %14.6g\n", n, e1, e2, e3 );
    free ( w );
    free ( x );
  }

  return;
}
/******************************************************************************/

void lagrange_value_test ( )

/******************************************************************************/
/*
  Purpose:

    LAGRANGE_VALUE_TEST tests LAGRANGE_VALUE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 November 2014

  Author:

    John Burkardt
*/
{
  int i;
  int j;
  double *li;
  int nd;
  int ni;
  double *xd;
  double xhi;
  double *xi;
  double xlo;

  printf ( "\n" );
  printf ( "LAGRANGE_VALUE_TEST\n" );
  printf ( "  LAGRANGE_VALUE evaluates the Lagrange basis polynomials.\n" );

  nd = 5;
  xlo = 0.0;
  xhi = ( double ) ( nd - 1 ) ;
  xd = r8vec_linspace_new ( nd, xlo, xhi );

  r8vec_print ( nd, xd, "  Lagrange basis points:" );
/*
  Evaluate the polynomials.
*/
  printf ( "\n" );
  printf ( "   I      X          L1(X)       L2(X)       L3(X)" );
  printf ( "       L4(X)       L5(X)\n" );
  printf ( "\n" );
 
  ni = 2 * nd - 1;
  xi = r8vec_linspace_new ( ni, xlo, xhi );

  li = lagrange_value ( nd, xd, ni, xi );

  for ( i = 0; i < ni; i++ )
  {
    printf ( "  %2d  %10.4f", i, xi[i] );
    for ( j = 0; j < nd; j++ )
    {
      printf ( "  %10.4f", li[i+j*ni] );
    }
    printf ( "\n" );
  }

  free ( li );
  free ( xd );
  free ( xi );

  return;
}
/******************************************************************************/

void lagrange_derivative_test ( )

/******************************************************************************/
/*
  Purpose:

    LAGRANGE_DERIVATIVE_TEST tests LAGRANGE_DERIVATIVE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 November 2014

  Author:

    John Burkardt
*/
{
  int i;
  int j;
  double *lpi;
  int nd;
  int ni;
  double *xd;
  double xhi;
  double *xi;
  double xlo;

  printf ( "\n" );
  printf ( "LAGRANGE_DERIVATIVE_TEST\n" );
  printf ( "  LAGRANGE_DERIVATIVE evaluates the Lagrange basis derivative.\n" );

  nd = 5;
  xlo = 0.0;
  xhi = ( double ) ( nd - 1 );
  xd = r8vec_linspace_new ( nd, xlo, xhi );

  r8vec_print ( nd, xd, "  Lagrange basis points:" );
/*
  Evaluate the polynomials.
*/
  printf ( "\n" );
  printf ( "   I      X         L1'(X)      L2'(X)      L3'(X)" );
  printf ( "      L4'(X)      L5'(X)\n" );
  printf ( "\n" );
 
  ni = 2 * nd - 1;
  xi = r8vec_linspace_new ( ni, xlo, xhi );
  lpi = lagrange_derivative ( nd, xd, ni, xi );

  for ( i = 0; i < ni; i++ )
  {
    printf ( "  %2d  %10.4f", i, xi[i] );
    for ( j = 0; j < nd; j++ )
    {
      printf ( "  %10.4f", lpi[i+j*ni] );
    }
    printf ( "\n" );
  }

  free ( lpi );
  free ( xd );
  free ( xi );

  return;
}
/******************************************************************************/

void fem1d_lagrange_stiffness_test ( int x_num, int q_num )

/******************************************************************************/
/*
  Purpose:

    FEM1D_LAGRANGE_STIFFNESS_TEST tests FEM1D_LAGRANGE_STIFFNESS.

  Discussion:

    The results are very sensitive to the quadrature rule.

    In particular, if X_NUM points are used, the mass matrix will
    involve integrals of polynomials of degree 2*(X_NUM-1), so the
    quadrature rule should use at least Q_NUM = X_NUM - 1 points.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 November 2014

  Author:

    John Burkardt

  Parameters:

    Input, int X_NUM, the number of nodes.

    Input, int Q_NUM, the number of quadrature points.
*/
{
  double *a;
  double *b;
  int i;
  int info;
  int j;
  double *k;
  double *m;
  double *u;
  double *u_e;
  double *x;
  double x_hi;
  double x_lo;

  printf ( "\n" );
  printf ( "FEM1D_LAGRANGE_STIFFNESS_TEST\n" );
  printf ( "  FEM1D_LAGRANGE_STIFFNESS computes the stiffness matrix,\n" );
  printf ( "  the mass matrix, and right hand side vector for a\n" );
  printf ( "  finite element problem using Lagrange interpolation\n" );
  printf ( "  basis polynomials.\n" );
  printf ( "\n" );
  printf ( "  Solving:\n" );
  printf ( "    -u''+u=x on 0 < x < 1\n" );
  printf ( "    u(0) = u(1) = 0\n" );
  printf ( "  Exact solution:\n" );
  printf ( "    u(x) = x - sinh(x)/sinh(1)\n" );
  printf ( "\n" );
  printf ( "  Number of mesh points = %d\n", x_num );
  printf ( "  Number of quadrature points = %d\n", q_num );

  x_lo = 0.0;
  x_hi = 1.0;
  x = r8vec_linspace_new ( x_num, x_lo, x_hi );

  a = ( double * ) malloc ( x_num * x_num * sizeof ( double ) );
  m = ( double * ) malloc ( x_num * x_num * sizeof ( double ) );
  b = ( double * ) malloc ( x_num * sizeof ( double ) );

  fem1d_lagrange_stiffness ( x_num, x, q_num, f, a, m, b );

  k = ( double * ) malloc ( x_num * x_num * sizeof ( double ) );

  for ( j = 0; j < x_num; j++ )
  {
    for ( i = 0; i < x_num; i++ )
    {
      k[i+j*x_num] = a[i+j*x_num] + m[i+j*x_num];
    }
  }
  for ( j = 0; j < x_num; j++ )
  {
    k[0+j*x_num] = 0.0;
  }
  k[0+0*x_num] = 1.0;
  b[0] = 0.0;

  for ( j = 0; j < x_num; j++ )
  {
    k[x_num-1+j*x_num] = 0.0;
  }
  k[x_num-1+(x_num-1)*x_num] = 1.0;
  b[x_num-1] = 0.0;

  u = r8mat_fs_new ( x_num, k, b );

  u_e = exact ( x_num, x );

  printf ( "\n" );
  printf ( "   I      X             U              U(exact)         Error\n" );
  printf ( "\n" );

  for ( i = 0; i < x_num; i++ )
  {
    printf ( "  %2d  %8.4f  %14.6g  %14.6g  %14.6g\n",
      i, x[i], u[i], u_e[i], fabs ( u[i] - u_e[i] ) );
  }

  free ( a );
  free ( b );
  free ( k );
  free ( m );
  free ( u );
  free ( u_e );
  free ( x );

  return;
}
/******************************************************************************/

double f ( double x )

/******************************************************************************/
/*
  Purpose:

    F evaluates the right hand side function.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 November 2014

  Author:

    John Burkardt

  Parameters:

    Input, double X, the evaluation point.

    Output, double F, the value of the right hand side at X.
*/
{
  double value;

  value = x;

  return value;
}
/******************************************************************************/

double *exact ( int x_num, double x[] )

/******************************************************************************/
/*
  Purpose:

    EXACT returns the exact solution.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 November 2014

  Author:

    John Burkardt

  Parameters:

    Input, int X_NUM, the number of nodes.

    Input, double X[X_NUM], the nodes.

    Output, double UE[X_NUM], the exact solution at the nodes.
*/
{
  double *ue;
  int x_i;

  ue = ( double * ) malloc ( x_num * sizeof ( double ) );

  for ( x_i = 0; x_i < x_num; x_i++ )
  {
    ue[x_i] = x[x_i] - sinh ( x[x_i] ) / sinh ( 1.0 );
  }
  return ue;
}

