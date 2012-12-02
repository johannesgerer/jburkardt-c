# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "fem1d_heat_steady.h"

int main ( void );
void fem1d_heat_steady_test01 ( void );
double k1 ( double x );
double f1 ( double x );
double exact1 ( double x );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    FEM1D_HEAT_STEADY_PRB tests the routines in FEM1D_HEAT_STEADY.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 August 2010

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "FEM1D_BVP_LINEAR_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the FEM1D_HEAT_STEADY library.\n" );

  fem1d_heat_steady_test01 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "FEM1D_HEAT_STEADY_PRB\n" );
  printf ( "  Normal end of execution.\n" );

  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void fem1d_heat_steady_test01 ( void )

/******************************************************************************/
/*
  Purpose:

    FEM1D_HEAT_STEADY_TEST01 carries out test case #1.

  Discussion:

    Use K1, F1, EXACT1.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 April 2011

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  int i;
  int n = 11;
  double *u;
  double ua;
  double ub;
  double uexact;
  double *x;

  printf ( "\n" );
  printf ( "FEM1D_HEAT_STEADY_TEST01\n" );
  printf ( "  K1(X)  = 1.0\n" );
  printf ( "  F1(X)  = X * ( X + 3 ) * exp ( X )\n" );
  printf ( "  U1(X)  = X * ( 1 - X ) * exp ( X )\n" );
/*
  Geometry definitions.
*/
  a = 0.0;
  b = 1.0;
  ua = 0.0;
  ub = 0.0;
  x = r8vec_even_new ( n, a, b );

  printf ( "\n" );
  printf ( "  Number of nodes = %d\n", n );
  printf ( "  Left endpoint A = %f\n", a );
  printf ( "  Right endpoint B = %f\n", b );
  printf ( "  Prescribed U(A) = %f\n", ua );
  printf ( "  Prescribed U(B) = %f\n", ub );

  u = fem1d_heat_steady ( n, a, b, ua, ub, k1, f1, x );

  printf ( "\n" );
  printf ( "     I         X          U                Uexact      Error\n" );
  printf ( "\n" );

  for ( i = 0; i < n; i++ )
  {
    uexact = exact1 ( x[i] );
    printf ( "  %4d  %8f  %14f  %14f  %14g\n",
      i, x[i], u[i], uexact, r8_abs ( u[i] - uexact ) );
  }

  free ( u );
  free ( x );

  return;
}
/******************************************************************************/

double k1 ( double x )

/******************************************************************************/
/*
  Purpose:

    K1 evaluates K function #1.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 August 2010

  Author:

    John Burkardt

  Parameters:

    Input, double X, the evaluation point.

    Output, double K1, the value of K(X).
*/
{
  double value;

  value = 1.0;

  return value;
}
/******************************************************************************/

double f1 ( double x )

/******************************************************************************/
/*
  Purpose:

    F1 evaluates right hand side function #1.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 August 2010

  Author:

    John Burkardt

  Parameters:

    Input, double X, the evaluation point.

    Output, double F1, the value of F(X).
*/
{
  double value;

  value = x * ( x + 3.0 ) * exp ( x );

  return value;
}
/******************************************************************************/

double exact1 ( double x )

/******************************************************************************/
/*
  Purpose:

    EXACT1 evaluates exact solution #1.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 August 2010

  Author:

    John Burkardt

  Parameters:

    Input, double X, the evaluation point.

    Output, double EXACT1, the value of U(X).
*/
{
  double value;

  value = x * ( 1.0 - x ) * exp ( x );

  return value;
}
