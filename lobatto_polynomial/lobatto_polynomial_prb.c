# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <string.h>

# include "lobatto_polynomial.h"

int main ( );
void lobatto_polynomial_value_test ( );
void lobatto_polynomial_derivative_test ( );
void lobatto_polynomial_plot_test ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    LOBATTO_POLYNOMIAL_TEST tests the LOBATTO_POLYNOMIAL library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 November 2014

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "LOBATTO_POLYNOMIAL_TEST:\n" );
  printf ( "  C version.\n" );
  printf ( "  Test the LOBATTO_POLYNOMIAL library.\n" );

  lobatto_polynomial_value_test ( );
  lobatto_polynomial_derivative_test ( );
  lobatto_polynomial_plot_test ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "LOBATTO_POLYNOMIAL_TEST\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void lobatto_polynomial_value_test ( )

/******************************************************************************/
/*
  Purpose:

    LOBATTO_POLYNOMIAL_VALUE_TEST tests LOBATTO_POLYNOMIAL_VALUE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 November 2014

  Author:

    John Burkardt
*/
{
  double e;
  double fx1;
  double fx2;
  double *l;
  int m;
  int n;
  int n_data;
  double x;
  double xvec[1];

  m = 1;

  printf ( "\n" );
  printf ( "LOBATTO_POLYNOMIAL_VALUE_TEST:\n" );
  printf ( "  LOBATTO_POLYNOMIAL_VALUES stores values of\n" );
  printf ( "  the completed Lobatto polynomial L(n,x).\n" );
  printf ( "  LOBATTO_POLYNOMIAL_VALUE evaluates the completed Lobatto polynomial.\n" );
  printf ( "\n" );
  printf ( "                                       Tabulated                 Computed\n" );
  printf ( "     N        X                        L(N,X)                    L(N,X)        Error\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    lobatto_polynomial_values ( &n_data, &n, &x, &fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    xvec[0] = x;

    l = lobatto_polynomial_value ( m, n, xvec );

    fx2 = l[0+(n-1)*m];

    e = fx1 - fx2;

    printf ( "  %4d  %12.4f  %24.16g  %24.16g  %8.1g\n", n, x, fx1, fx2, e );

    free ( l );
  }

  return;
}
/******************************************************************************/

void lobatto_polynomial_derivative_test ( )

/******************************************************************************/
/*
  Purpose:

    LOBATTO_POLYNOMIAL_DERIVATIVE_TEST tests LOBATTO_POLYNOMIAL_DERIVATIVE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 November 2014

  Author:

    John Burkardt
*/
{
  double e;
  double fx1;
  double fx2;
  double *lp;
  int m;
  int n;
  int n_data;
  double x;
  double xvec[1];

  m = 1;

  printf ( "\n" );
  printf ( "LOBATTO_POLYNOMIAL_DERIVATIVE_TEST:\n" );
  printf ( "  LOBATTO_POLYNOMIAL_DERIVATIVES stores derivatives of\n" );
  printf ( "  the completed Lobatto polynomial L(n,x).\n" );
  printf ( "  LOBATTO_POLYNOMIAL_DERIVATIVE evaluates the completed Lobatto polynomial.\n" );
  printf ( "\n" );
  printf ( "                                       Tabulated                 Computed\n" );
  printf ( "     N        X                        L''(N,X)                   L''(N,X)       Error\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    lobatto_polynomial_derivatives ( &n_data, &n, &x, &fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    xvec[0] = x;
    lp = lobatto_polynomial_derivative ( m, n, xvec );
    fx2 = lp[0+(n-1)*m];

    e = fx1 - fx2;

    printf ( "  %4d  %12.4f  %24.16g  %24.16g  %8.1g\n", n, x, fx1, fx2, e );

    free ( lp );
  }

  return;
}
/******************************************************************************/

void lobatto_polynomial_plot_test ( )

/******************************************************************************/
/*
  Purpose:

    LOBATTO_POLYNOMIAL_PLOT_TEST tests LOBATTO_POLYNOMIAL_PLOT.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 November 2014

  Author:

    John Burkardt
*/
{
  int ndx[7] = { 1, 2, 3, 4, 5, 6, 7 };
  int ndx_num = 7;
  char prefix[] = "test";

  printf ( "\n" );
  printf ( "LOBATTO_POLYNOMIAL_PLOT_TEST:\n" );
  printf ( "  LOBATTO_POLYNOMIAL_PLOT plots Lobatto polynomials.\n" );

  lobatto_polynomial_plot ( ndx_num, ndx, prefix );

  return;
}
