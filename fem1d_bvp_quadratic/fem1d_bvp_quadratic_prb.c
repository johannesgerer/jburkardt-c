# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "fem1d_bvp_quadratic.h"

int main ( );
void test01 ( );
double a1 ( double x );
double c1 ( double x );
double exact1 ( double x );
double exact_ux1 ( double x );
double f1 ( double x );

void test02 ( );
double a2 ( double x );
double c2 ( double x );
double exact2 ( double x );
double exact_ux2 ( double x );
double f2 ( double x );

void test03 ( );
double a3 ( double x );
double c3 ( double x );
double exact3 ( double x );
double exact_ux3 ( double x );
double f3 ( double x );

void test04 ( );
double a4 ( double x );
double c4 ( double x );
double exact4 ( double x );
double exact_ux4 ( double x );
double f4 ( double x );

void test05 ( );
double a5 ( double x );
double c5 ( double x );
double exact5 ( double x );
double exact_ux5 ( double x );
double f5 ( double x );

void test06 ( );
double a6 ( double x );
double c6 ( double x );
double exact6 ( double x );
double exact_ux6 ( double x );
double f6 ( double x );

void test07 ( );
double a7 ( double x );
double c7 ( double x );
double exact7 ( double x );
double exact_ux7 ( double x );
double f7 ( double x );

void test08 ( );
double a8 ( double x );
double c8 ( double x );
double exact8 ( double x );
double exact_ux8 ( double x );
double f8 ( double x );

void test09 ( );
double a9 ( double x );
double c9 ( double x );
double exact9 ( double x );
double exact_ux9 ( double x );
double f9 ( double x );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for FEM1D_BVP_QUADRATIC_PRB.

  Discussion:

    FEM1D_BVP_QUADRATIC_PRB tests the FEM1D_BVP_QUADRATIC library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    18 June 2014

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "FEM1D_BVP_QUADRATIC_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the FEM1D_BVP_QUADRATIC library.\n" );

  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
  test05 ( );
  test06 ( );
  test07 ( );
  test08 ( );
  test09 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "FEM1D_BVP_QUADRATIC_PRB\n" );
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

    TEST01 carries out test case #1.

  Discussion:

    Use A1, C1, F1, EXACT1, EXACT_UX1.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 June 2014

  Author:

    John Burkardt

  Reference:

    Dianne O'Leary,
    Scientific Computing with Case Studies,
    SIAM, 2008,
    ISBN13: 978-0-898716-66-5,
    LC: QA401.O44.
*/
{
  int i;
  int n = 11;
  double e1;
  double e2;
  double h1s;
  double *u;
  double uexact;
  double *x;
  double x_first;
  double x_last;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  Solve -( A(x) U'(x) )' + C(x) U(x) = F(x)\n" );
  printf ( "  for 0 < x < 1, with U(0) = U(1) = 0.\n" );
  printf ( "  A1(X)  = 1.0\n" );
  printf ( "  C1(X)  = 0.0\n" );
  printf ( "  F1(X)  = X * ( X + 3 ) * exp ( X )\n" );
  printf ( "  U1(X)  = X * ( 1 - X ) * exp ( X )\n" );
  printf ( "\n" );
  printf ( "  Number of nodes = %d\n", n );
/*
  Geometry definitions.
*/
  x_first = 0.0;
  x_last = 1.0;
  x = r8vec_even_new ( n, x_first, x_last );

  u = fem1d_bvp_quadratic ( n, a1, c1, f1, x );

  printf ( "\n" );
  printf ( "     I    X         U         Uexact    Error\n" );
  printf ( "\n" );

  for ( i = 0; i < n; i++ )
  {
    uexact = exact1 ( x[i] );
    printf ( "  %4d  %8f  %14f  %14f %14e\n", 
      i, x[i], u[i], uexact, fabs ( u[i] - uexact ) );
  }

  e1 = l1_error ( n, x, u, exact1 );
  e2 = l2_error_quadratic ( n, x, u, exact1 );
  h1s = h1s_error_quadratic ( n, x, u, exact_ux1 );
  printf ( "\n" );
  printf ( "  l1 norm of error  = %g\n", e1 );
  printf ( "  L2 norm of error  = %g\n", e2 );
  printf ( "  Seminorm of error = %g\n", h1s );

  free ( u );
  free ( x );

  return;
}
/******************************************************************************/

double a1 ( double x )

/******************************************************************************/
/*
  Purpose:

    A1 evaluates A function #1.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 August 2010

  Author:

    John Burkardt

  Parameters:

    Input, double X, the evaluation point.

    Output, double A1, the value of A(X).
*/
{
  double value;

  value = 1.0;

  return value;
}
/******************************************************************************/

double c1 ( double x )

/******************************************************************************/
/*
  Purpose:

    C1 evaluates C function #1.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 August 2010

  Author:

    John Burkardt

  Parameters:

    Input, double X, the evaluation point.

    Output, double C1, the value of C(X).
*/
{
  double value;

  value = 0.0;

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

    14 June 2014

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
/******************************************************************************/

double exact_ux1 ( double x )

/******************************************************************************/
/*
  Purpose:

    EXACT_UX1 evaluates the derivative of exact solution #1.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 June 2014

  Author:

    John Burkardt

  Parameters:

    Input, double X, the evaluation point.

    Output, double EXACT_UX1, the value of dUdX(X).
*/
{
  double value;

  value = ( 1.0 - x - x * x ) * exp ( x );

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

void test02 ( )

/******************************************************************************/
/*
  Purpose:

    TEST02 carries out test case #2.

  Discussion:

    Use A2, C2, F2, EXACT2, EXACT_UX2.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 June 2014

  Author:

    John Burkardt

  Reference:

    Dianne O'Leary,
    Scientific Computing with Case Studies,
    SIAM, 2008,
    ISBN13: 978-0-898716-66-5,
    LC: QA401.O44.
*/
{
  int i;
  int n = 11;
  double e1;
  double e2;
  double h1s;
  double *u;
  double uexact;
  double *x;
  double x_first;
  double x_last;

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  Solve -( A(x) U'(x) )' + C(x) U(x) = F(x)\n" );
  printf ( "  for 0 < x < 1, with U(0) = U(1) = 0.\n" );
  printf ( "  A2(X)  = 1.0\n" );
  printf ( "  C2(X)  = 2.0\n" );
  printf ( "  F2(X)  = X * ( 5 - X ) * exp ( X )\n" );
  printf ( "  U2(X)  = X * ( 1 - X ) * exp ( X )\n" );
  printf ( "\n" );
  printf ( "  Number of nodes = %d\n", n );
/*
  Geometry definitions.
*/
  x_first = 0.0;
  x_last = 1.0;
  x = r8vec_even_new ( n, x_first, x_last );

  u = fem1d_bvp_quadratic ( n, a2, c2, f2, x );

  printf ( "\n" );
  printf ( "     I    X         U         Uexact    Error\n" );
  printf ( "\n" );

  for ( i = 0; i < n; i++ )
  {
    uexact = exact2 ( x[i] );
    printf ( "  %4d  %8f  %14f  %14f %14e\n", 
      i, x[i], u[i], uexact, fabs ( u[i] - uexact ) );
  }

  e1 = l1_error ( n, x, u, exact2 );
  e2 = l2_error_quadratic ( n, x, u, exact2 );
  h1s = h1s_error_quadratic ( n, x, u, exact_ux2 );
  printf ( "\n" );
  printf ( "  l1 norm of error  = %g\n", e1 );
  printf ( "  L2 norm of error  = %g\n", e2 );
  printf ( "  Seminorm of error = %g\n", h1s );

  free ( u );
  free ( x );

  return;
}
/******************************************************************************/

double a2 ( double x )

/******************************************************************************/
/*
  Purpose:

    A2 evaluates A function #2.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 June 2014

  Author:

    John Burkardt

  Parameters:

    Input, double X, the evaluation point.

    Output, double A2, the value of A(X).
*/
{
  double value;

  value = 1.0;

  return value;
}
/******************************************************************************/

double c2 ( double x )

/******************************************************************************/
/*
  Purpose:

    C2 evaluates C function #2.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 August 2010

  Author:

    John Burkardt

  Parameters:

    Input, double X, the evaluation point.

    Output, double C2, the value of C(X).
*/
{
  double value;

  value = 2.0;

  return value;
}
/******************************************************************************/

double exact2 ( double x )

/******************************************************************************/
/*
  Purpose:

    EXACT2 evaluates exact solution #2.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 June 2014

  Author:

    John Burkardt

  Parameters:

    Input, double X, the evaluation point.

    Output, double EXACT2, the value of U(X).
*/
{
  double value;

  value = x * ( 1.0 - x ) * exp ( x );

  return value;
}
/******************************************************************************/

double exact_ux2 ( double x )

/******************************************************************************/
/*
  Purpose:

    EXACT_UX2 evaluates the derivative of exact solution #2.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 June 2014

  Author:

    John Burkardt

  Parameters:

    Input, double X, the evaluation point.

    Output, double EXACT_UX2, the value of dUdX(X).
*/
{
  double value;

  value = ( 1.0 - x - x * x ) * exp ( x );

  return value;
}
/******************************************************************************/

double f2 ( double x )

/******************************************************************************/
/*
  Purpose:

    F2 evaluates right hand side function #2.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 August 2010

  Author:

    John Burkardt

  Parameters:

    Input, double X, the evaluation point.

    Output, double F2, the value of F(X).
*/
{
  double value;

  value = x * ( 5.0 - x ) * exp ( x );

  return value;
}
/******************************************************************************/

void test03 ( )

/******************************************************************************/
/*
  Purpose:

    TEST03 carries out test case #3.

  Discussion:

    Use A3, C3, F3, EXACT3, EXACT_UX3.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 June 2014

  Author:

    John Burkardt

  Reference:

    Dianne O'Leary,
    Scientific Computing with Case Studies,
    SIAM, 2008,
    ISBN13: 978-0-898716-66-5,
    LC: QA401.O44.
*/
{
  int i;
  int n = 11;
  double e1;
  double e2;
  double h1s;
  double *u;
  double uexact;
  double *x;
  double x_first;
  double x_last;

  printf ( "\n" );
  printf ( "TEST03\n" );
  printf ( "  Solve -( A(x) U'(x) )' + C(x) U(x) = F(x)\n" );
  printf ( "  for 0 < x < 1, with U(0) = U(1) = 0.\n" );
  printf ( "  A3(X)  = 1.0\n" );
  printf ( "  C3(X)  = 2.0 * X\n" );
  printf ( "  F3(X)  = - X * ( 2 * X * X - 3 * X - 3 ) * exp ( X )\n" );
  printf ( "  U3(X)  = X * ( 1 - X ) * exp ( X )\n" );
  printf ( "\n" );
  printf ( "  Number of nodes = %d\n", n );
/*
  Geometry definitions.
*/
  x_first = 0.0;
  x_last = 1.0;
  x = r8vec_even_new ( n, x_first, x_last );

  u = fem1d_bvp_quadratic ( n, a3, c3, f3, x );

  printf ( "\n" );
  printf ( "     I    X         U         Uexact    Error\n" );
  printf ( "\n" );

  for ( i = 0; i < n; i++ )
  {
    uexact = exact3 ( x[i] );
    printf ( "  %4d  %8f  %14f  %14f %14e\n", 
      i, x[i], u[i], uexact, fabs ( u[i] - uexact ) );
  }

  e1 = l1_error ( n, x, u, exact3 );
  e2 = l2_error_quadratic ( n, x, u, exact3 );
  h1s = h1s_error_quadratic ( n, x, u, exact_ux3 );
  printf ( "\n" );
  printf ( "  l1 norm of error  = %g\n", e1 );
  printf ( "  L2 norm of error  = %g\n", e2 );
  printf ( "  Seminorm of error = %g\n", h1s );

  free ( u );
  free ( x );

  return;
}
/******************************************************************************/

double a3 ( double x )

/******************************************************************************/
/*
  Purpose:

    A3 evaluates A function #3.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 June 2014

  Author:

    John Burkardt

  Parameters:

    Input, double X, the evaluation point.

    Output, double A3, the value of A(X).
*/
{
  double value;

  value = 1.0;

  return value;
}
/******************************************************************************/

double c3 ( double x )

/******************************************************************************/
/*
  Purpose:

    C3 evaluates C function #3.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 August 2010

  Author:

    John Burkardt

  Parameters:

    Input, double X, the evaluation point.

    Output, double C3, the value of C(X).
*/
{
  double value;

  value = 2.0 * x;

  return value;
}
/******************************************************************************/

double exact3 ( double x )

/******************************************************************************/
/*
  Purpose:

    EXACT3 evaluates exact solution #3.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 June 2014

  Author:

    John Burkardt

  Parameters:

    Input, double X, the evaluation point.

    Output, double EXACT3, the value of U(X).
*/
{
  double value;

  value = x * ( 1.0 - x ) * exp ( x );

  return value;
}
/******************************************************************************/

double exact_ux3 ( double x )

/******************************************************************************/
/*
  Purpose:

    EXACT_UX3 evaluates the derivative of exact solution #3.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 June 2014

  Author:

    John Burkardt

  Parameters:

    Input, double X, the evaluation point.

    Output, double EXACT_UX3, the value of dUdX(X).
*/
{
  double value;

  value = ( 1.0 - x - x * x ) * exp ( x );

  return value;
}
/******************************************************************************/

double f3 ( double x )

/******************************************************************************/
/*
  Purpose:

    F3 evaluates right hand side function #3.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 August 2010

  Author:

    John Burkardt

  Parameters:

    Input, double X, the evaluation point.

    Output, double F3, the value of F(X).
*/
{
  double value;

  value = - x * ( 2.0 * x * x - 3.0 * x - 3.0 ) * exp ( x );

  return value;
}
/******************************************************************************/

void test04 ( )

/******************************************************************************/
/*
  Purpose:

    TEST04 carries out test case #4.

  Discussion:

    Use A4, C4, F4, EXACT4, EXACT_UX4.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 June 2014

  Author:

    John Burkardt

  Reference:

    Dianne O'Leary,
    Scientific Computing with Case Studies,
    SIAM, 2008,
    ISBN13: 978-0-898716-66-5,
    LC: QA401.O44.
*/
{
  int i;
  int n = 11;
  double e1;
  double e2;
  double h1s;
  double *u;
  double uexact;
  double *x;
  double x_first;
  double x_last;

  printf ( "\n" );
  printf ( "TEST04\n" );
  printf ( "  Solve -( A(x) U'(x) )' + C(x) U(x) = F(x)\n" );
  printf ( "  for 0 < x < 1, with U(0) = U(1) = 0.\n" );
  printf ( "  A4(X)  = 1.0 + X * X\n" );
  printf ( "  C4(X)  = 0.0\n" );
  printf ( "  F4(X)  = ( X + 3 X^2 + 5 X^3 + X^4 ) * exp ( X )\n" );
  printf ( "  U4(X)  = X * ( 1 - X ) * exp ( X )\n" );
  printf ( "\n" );
  printf ( "  Number of nodes = %d\n", n );
/*
  Geometry definitions.
*/
  x_first = 0.0;
  x_last = 1.0;
  x = r8vec_even_new ( n, x_first, x_last );

  u = fem1d_bvp_quadratic ( n, a4, c4, f4, x );

  printf ( "\n" );
  printf ( "     I    X         U         Uexact    Error\n" );
  printf ( "\n" );

  for ( i = 0; i < n; i++ )
  {
    uexact = exact4 ( x[i] );
    printf ( "  %4d  %8f  %14f  %14f %14e\n", 
      i, x[i], u[i], uexact, fabs ( u[i] - uexact ) );
  }

  e1 = l1_error ( n, x, u, exact4 );
  e2 = l2_error_quadratic ( n, x, u, exact4 );
  h1s = h1s_error_quadratic ( n, x, u, exact_ux4 );
  printf ( "\n" );
  printf ( "  l1 norm of error  = %g\n", e1 );
  printf ( "  L2 norm of error  = %g\n", e2 );
  printf ( "  Seminorm of error = %g\n", h1s );

  free ( u );
  free ( x );

  return;
}
/******************************************************************************/

double a4 ( double x )

/******************************************************************************/
/*
  Purpose:

    A4 evaluates A function #4.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 June 2014

  Author:

    John Burkardt

  Parameters:

    Input, double X, the evaluation point.

    Output, double A4, the value of A(X).
*/
{
  double value;

  value = 1.0 + x * x;

  return value;
}
/******************************************************************************/

double c4 ( double x )

/******************************************************************************/
/*
  Purpose:

    C4 evaluates C function #4.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 June 2014

  Author:

    John Burkardt

  Parameters:

    Input, double X, the evaluation point.

    Output, double C4, the value of C(X).
*/
{
  double value;

  value = 0.0;

  return value;
}
/******************************************************************************/

double exact4 ( double x )

/******************************************************************************/
/*
  Purpose:

    EXACT4 evaluates exact solution #4.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 June 2014

  Author:

    John Burkardt

  Parameters:

    Input, double X, the evaluation point.

    Output, double EXACT4, the value of U(X).
*/
{
  double value;

  value = x * ( 1.0 - x ) * exp ( x );

  return value;
}
/******************************************************************************/

double exact_ux4 ( double x )

/******************************************************************************/
/*
  Purpose:

    EXACT_UX4 evaluates the derivative of exact solution #4.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 June 2014

  Author:

    John Burkardt

  Parameters:

    Input, double X, the evaluation point.

    Output, double EXACT_UX4, the value of dUdX(X).
*/
{
  double value;

  value = ( 1.0 - x - x * x ) * exp ( x );

  return value;
}
/******************************************************************************/

double f4 ( double x )

/******************************************************************************/
/*
  Purpose:

    F4 evaluates right hand side function #4.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 August 2010

  Author:

    John Burkardt

  Parameters:

    Input, double X, the evaluation point.

    Output, double F4, the value of F(X).
*/
{
  double value;

  value = ( x + 3.0 * x * x + 5.0 * x * x * x + x * x * x * x ) * exp ( x );

  return value;
}
/******************************************************************************/

void test05 ( )

/******************************************************************************/
/*
  Purpose:

    TEST05 carries out test case #5.

  Discussion:

    Use A5, C5, F5, EXACT5, EXACT_UX5.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 June 2014

  Author:

    John Burkardt

  Reference:

    Dianne O'Leary,
    Scientific Computing with Case Studies,
    SIAM, 2008,
    ISBN13: 978-0-898716-66-5,
    LC: QA401.O44.
*/
{
  int i;
  int n = 11;
  double e1;
  double e2;
  double h1s;
  double *u;
  double uexact;
  double *x;
  double x_first;
  double x_last;

  printf ( "\n" );
  printf ( "TEST05\n" );
  printf ( "  Solve -( A(x) U'(x) )' + C(x) U(x) = F(x)\n" );
  printf ( "  for 0 < x < 1, with U(0) = U(1) = 0.\n" );
  printf ( "  A5(X)  = 1.0 + X * X for X <= 1/3\n" );
  printf ( "         = 7/9 + X     for      1/3 < X\n" );
  printf ( "  C5(X)  = 0.0\n" );
  printf ( "  F5(X)  = ( X + 3 X^2 + 5 X^3 + X^4 ) * exp ( X )\n" );
  printf ( "                       for X <= 1/3\n" );
  printf ( "         = ( - 1 + 10/3 X + 43/9 X^2 + X^3 ) .* exp ( X )\n" );
  printf ( "                       for      1/3 <= X\n" );
  printf ( "  U5(X)  = X * ( 1 - X ) * exp ( X )\n" );
  printf ( "\n" );
  printf ( "  Number of nodes = %d\n", n );
/*
  Geometry definitions.
*/
  x_first = 0.0;
  x_last = 1.0;
  x = r8vec_even_new ( n, x_first, x_last );

  u = fem1d_bvp_quadratic ( n, a5, c5, f5, x );

  printf ( "\n" );
  printf ( "     I    X         U         Uexact    Error\n" );
  printf ( "\n" );

  for ( i = 0; i < n; i++ )
  {
    uexact = exact5 ( x[i] );
    printf ( "  %4d  %8f  %14f  %14f %14e\n", 
      i, x[i], u[i], uexact, fabs ( u[i] - uexact ) );
  }

  e1 = l1_error ( n, x, u, exact5 );
  e2 = l2_error_quadratic ( n, x, u, exact5 );
  h1s = h1s_error_quadratic ( n, x, u, exact_ux5 );
  printf ( "\n" );
  printf ( "  l1 norm of error  = %g\n", e1 );
  printf ( "  L2 norm of error  = %g\n", e2 );
  printf ( "  Seminorm of error = %g\n", h1s );

  free ( u );
  free ( x );

  return;
}
/******************************************************************************/

double a5 ( double x )

/******************************************************************************/
/*
  Purpose:

    A5 evaluates A function #5.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 June 2014

  Author:

    John Burkardt

  Parameters:

    Input, double X, the evaluation point.

    Output, double A5, the value of A(X).
*/
{
  double value;

  if ( x <= 1.0 / 3.0 )
  {
    value = 1.0 + x * x;
  }
  else
  {
    value = x + 7.0 / 9.0;
  }

  return value;
}
/******************************************************************************/

double c5 ( double x )

/******************************************************************************/
/*
  Purpose:

    C5 evaluates C function #5.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 June 2014

  Author:

    John Burkardt

  Parameters:

    Input, double X, the evaluation point.

    Output, double C5, the value of C(X).
*/
{
  double value;

  value = 0.0;

  return value;
}
/******************************************************************************/

double exact5 ( double x )

/******************************************************************************/
/*
  Purpose:

    EXACT5 evaluates exact solution #5.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 June 2014

  Author:

    John Burkardt

  Parameters:

    Input, double X, the evaluation point.

    Output, double EXACT5, the value of U(X).
*/
{
  double value;

  value = x * ( 1.0 - x ) * exp ( x );

  return value;
}
/******************************************************************************/

double exact_ux5 ( double x )

/******************************************************************************/
/*
  Purpose:

    EXACT_UX5 evaluates the derivative of exact solution #5.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 June 2014

  Author:

    John Burkardt

  Parameters:

    Input, double X, the evaluation point.

    Output, double EXACT_UX5, the value of dUdX(X).
*/
{
  double value;

  value = ( 1.0 - x - x * x ) * exp ( x );

  return value;
}
/******************************************************************************/

double f5 ( double x )

/******************************************************************************/
/*
  Purpose:

    F5 evaluates right hand side function #5.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 August 2010

  Author:

    John Burkardt

  Parameters:

    Input, double X, the evaluation point.

    Output, double F5, the value of F(X).
*/
{
  double value;

  if ( x <= 1.0 / 3.0 )
  {
    value = ( x + 3.0 * x * x + 5.0 * x * x * x + x * x * x * x ) * exp ( x );
  }
  else
  {
    value = ( - 1.0 + ( 10.0 / 3.0 ) * x 
      + ( 43.0 / 9.0 ) * x * x + x * x * x ) * exp ( x );
  }

  return value;
}
/******************************************************************************/

void test06 ( )

/******************************************************************************/
/*
  Purpose:

    TEST06 does an error analysis.

  Discussion:

    Use A6, C6, F6, EXACT6, EXACT_UX6.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 June 2014

  Author:

    John Burkardt
*/
{
  int i;
  int n;
  double e1;
  double e2;
  double h1s;
  double *u;
  double uexact;
  double *x;
  double x_first;
  double x_last;

  printf ( "\n" );
  printf ( "TEST06\n" );
  printf ( "  Solve -( A(x) U'(x) )' + C(x) U(x) = F(x)\n" );
  printf ( "  for 0 < x < 1, with U(0) = U(1) = 0.\n" );
  printf ( "  A6(X)  = 1.0\n" );
  printf ( "  C6(X)  = 0.0\n" );
  printf ( "  F6(X)  = pi*pi*sin(pi*X)\n" );
  printf ( "  U6(X)  = sin(pi*x)\n" );
  printf ( "\n" );
  printf ( "  Compute l1norm, L2norm and seminorm of error for various N.\n" );
  printf ( "\n" );
  printf ( "     N        l1 error       L2 error      Seminorm error\n" );
  printf ( "\n" );

  n = 11;
  for ( i = 0; i <= 4; i++ )
  {
/*
  Geometry definitions.
*/
    x_first = 0.0;
    x_last = 1.0;
    x = r8vec_even_new ( n, x_first, x_last );

    u = fem1d_bvp_quadratic ( n, a6, c6, f6, x );

    e1 = l1_error ( n, x, u, exact6 );
    e2 = l2_error_quadratic ( n, x, u, exact6 );
    h1s = h1s_error_quadratic ( n, x, u, exact_ux6 );

    printf ( "  %4d  %14f  %14f  %14f\n", 
      n, e1, e2, h1s );

    n = 2 * ( n - 1 ) + 1;

    free ( u );
    free ( x );
  }

  return;
}
/******************************************************************************/

double a6 ( double x )

/******************************************************************************/
/*
  Purpose:

    A6 evaluates A function #6.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 June 2014

  Author:

    John Burkardt

  Parameters:

    Input, double X, the evaluation point.

    Output, double A6, the value of A(X).
*/
{
  double value;

  value = 1.0;

  return value;
}
/******************************************************************************/

double c6 ( double x )

/******************************************************************************/
/*
  Purpose:

    C6 evaluates C function #6.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 June 2014

  Author:

    John Burkardt

  Parameters:

    Input, double X, the evaluation point.

    Output, double C6, the value of C(X).
*/
{
  double value;

  value = 0.0;

  return value;
}
/******************************************************************************/

double exact6 ( double x )

/******************************************************************************/
/*
  Purpose:

    EXACT6 returns exact solution #6.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 June 2014

  Author:

    John Burkardt

  Parameters:

    Input, double X, the evaluation point.

    Output, double EXACT6, the value of U(X).
*/
{
  const double pi = 3.141592653589793;
  double value;

  value = sin ( pi * x );

  return value;
}
/******************************************************************************/

double exact_ux6 ( double x )

/******************************************************************************/
/*
  Purpose:

    EXACT_UX6 returns the derivative of exact solution #6.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 June 2014

  Author:

    John Burkardt

  Parameters:

    Input, double X, the evaluation point.

    Output, double EXACT_UX6, the value of U(X).
*/
{
  const double pi = 3.141592653589793;
  double value;

  value = pi * cos ( pi * x );

  return value;
}
/******************************************************************************/

double f6 ( double x )

/******************************************************************************/
/*
  Purpose:

    F6 evaluates right hand side function #6.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 February 2012

  Author:

    John Burkardt

  Parameters:

    Input, double X, the evaluation point.

    Output, double F6, the value of F(X).
*/
{
  static double pi = 3.141592653589793;
  double value;

  value = pi * pi * sin ( pi * x );

  return value;
}
/******************************************************************************/

void test07 ( )

/******************************************************************************/
/*
  Purpose:

    TEST07 does an error analysis.

  Discussion:

    Use A7, C7, F7, EXACT7, EXACT_UX7.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 June 2014

  Author:

    John Burkardt

  Reference:

    Eric Becker, Graham Carey, John Oden,
    Finite Elements, An Introduction, Volume I,
    Prentice-Hall, 1981, page 123-124,
    ISBN: 0133170578,
    LC: TA347.F5.B4.
*/
{
  int i;
  int n;
  double e1;
  double e2;
  double h1s;
  double *u;
  double uexact;
  double *x;
  double x_first;
  double x_last;

  printf ( "\n" );
  printf ( "TEST07\n" );
  printf ( "  Solve -( A(x) U'(x) )' + C(x) U(x) = F(x)\n" );
  printf ( "  for 0 < x < 1, with U(0) = U(1) = 0.\n" );
  printf ( "  Becker/Carey/Oden Example\n" );
  printf ( "\n" );
  printf ( "  Compute l1 norm, L2 norm and seminorm of error for various N.\n" );
  printf ( "\n" );
  printf ( "     N        l1 error      L2 error      Seminorm error\n" );
  printf ( "\n" );

  n = 11;
  for ( i = 0; i <= 4; i++ )
  {
/*
  Geometry definitions.
*/
    x_first = 0.0;
    x_last = 1.0;
    x = r8vec_even_new ( n, x_first, x_last );

    u = fem1d_bvp_quadratic ( n, a7, c7, f7, x );

    e1 = l1_error ( n, x, u, exact7 );
    e2 = l2_error_quadratic ( n, x, u, exact7 );
    h1s = h1s_error_quadratic ( n, x, u, exact_ux7 );

    printf ( "  %4d  %14f  %14f  %14f\n", 
      n, e1, e2, h1s );

    n = 2 * ( n - 1 ) + 1;

    free ( u );
    free ( x );
  }

  return;
}
/******************************************************************************/

double a7 ( double x )

/******************************************************************************/
/*
  Purpose:

    A7 evaluates A function #7.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 June 2014

  Author:

    John Burkardt

  Parameters:

    Input, double X, the evaluation point.

    Output, double A7, the value of A(X).
*/
{
  double alpha;
  double value;
  double x0;

  alpha = 30.0;
  x0 = 1.0 / 3.0;
  value = 1.0 / alpha + alpha * pow ( x - x0, 2 );

  return value;
}
/******************************************************************************/

double c7 ( double x )

/******************************************************************************/
/*
  Purpose:

    C7 evaluates C function #7.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 June 2014

  Author:

    John Burkardt

  Parameters:

    Input, double X, the evaluation point.

    Output, double C7, the value of C(X).
*/
{
  double value;

  value = 0.0;

  return value;
}
/******************************************************************************/

double exact7 ( double x )

/******************************************************************************/
/*
  Purpose:

    EXACT7 returns exact solution #7.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 June 2014

  Author:

    John Burkardt

  Parameters:

    Input, double X, the evaluation point.

    Output, double EXACT7, the value of U(X).
*/
{
  double alpha;
  double value;
  double x0;

  alpha = 30.0;
  x0 = 1.0 / 3.0;
  value = ( 1.0 - x ) 
    * ( atan ( alpha * ( x - x0 ) ) + atan ( alpha * x0 ) );

  return value;
}
/******************************************************************************/

double exact_ux7 ( double x )

/******************************************************************************/
/*
  Purpose:

    EXACT_UX7 returns the derivative of exact solution #7.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 June 2014

  Author:

    John Burkardt

  Parameters:

    Input, double X, the evaluation point.

    Output, double EXACT_UX7, the value of U(X).
*/
{
  double alpha;
  double value;
  double x0;

  alpha = 30.0;
  x0 = 1.0 / 3.0;
  value = - atan ( alpha * ( x - x0 ) ) - atan ( alpha * x0 ) 
    + ( 1.0 - x ) * alpha / ( 1.0 + alpha * alpha * pow ( x - x0, 2 ) );

  return value;
}
/******************************************************************************/

double f7 ( double x )

/******************************************************************************/
/*
  Purpose:

    F7 evaluates right hand side function #7.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 June 2014

  Author:

    John Burkardt

  Parameters:

    Input, double X, the evaluation point.

    Output, double F7, the value of F(X).
*/
{
  double alpha;
  double value;
  double x0;

  alpha = 30.0;
  x0 = 1.0 / 3.0;
  value = 2.0 * ( 1.0 + alpha * ( x - x0 ) * 
    ( atan ( alpha * ( x - x0 ) ) + atan ( alpha * x0 ) ) );

  return value;
}
/******************************************************************************/

void test08 ( )

/******************************************************************************/
/*
  Purpose:

    TEST08 carries out test case #8.

  Discussion:

    Use A8, C8, F8, EXACT8, EXACT_UX8.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 June 2014

  Author:

    John Burkardt

  Reference:

    Dianne O'Leary,
    Scientific Computing with Case Studies,
    SIAM, 2008,
    ISBN13: 978-0-898716-66-5,
    LC: QA401.O44.
*/
{
  int i;
  int n = 11;
  double e1;
  double e2;
  double h1s;
  double *u;
  double uexact;
  double *x;
  double x_first;
  double x_last;

  printf ( "\n" );
  printf ( "TEST08\n" );
  printf ( "  Solve -( A(x) U'(x) )' + C(x) U(x) = F(x)\n" );
  printf ( "  for 0 < x < 1, with U(0) = U(1) = 0.\n" );
  printf ( "  A8(X) = 1.0\n" );
  printf ( "  C8(X) = 0.0\n" );
  printf ( "  F8(X) = X * ( X + 3 ) * exp ( X ),   X <= 2/3\n" );
  printf ( "        = 2 * exp ( 2/3),                   2/3 < X\n" );
  printf ( "  U8(X) = X * ( 1 - X ) * exp ( X ),   X <= 2/3\n" );
  printf ( "        = X * ( 1 - X ) * exp ( 2/3 ),      2/3 < X\n" );
  printf ( "\n" );
  printf ( "  Number of nodes = %d\n", n );
/*
  Geometry definitions.
*/
  x_first = 0.0;
  x_last = 1.0;
  x = r8vec_even_new ( n, x_first, x_last );

  u = fem1d_bvp_quadratic ( n, a8, c8, f8, x );

  printf ( "\n" );
  printf ( "     I    X         U         Uexact    Error\n" );
  printf ( "\n" );

  for ( i = 0; i < n; i++ )
  {
    uexact = exact8 ( x[i] );
    printf ( "  %4d  %8f  %14f  %14f %14e\n", 
      i, x[i], u[i], uexact, fabs ( u[i] - uexact ) );
  }

  e1 = l1_error ( n, x, u, exact8 );
  e2 = l2_error_quadratic ( n, x, u, exact8 );
  h1s = h1s_error_quadratic ( n, x, u, exact_ux8 );
  printf ( "\n" );
  printf ( "  l1 norm of error  = %g\n", e1 );
  printf ( "  L2 norm of error  = %g\n", e2 );
  printf ( "  Seminorm of error = %g\n", h1s );

  free ( u );
  free ( x );

  return;
}
/******************************************************************************/

double a8 ( double x )

/******************************************************************************/
/*
  Purpose:

    A8 evaluates A function #8.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 June 2014

  Author:

    John Burkardt

  Parameters:

    Input, double X, the evaluation point.

    Output, double A8, the value of A(X).
*/
{
  double value;

  value = 1.0;

  return value;
}
/******************************************************************************/

double c8 ( double x )

/******************************************************************************/
/*
  Purpose:

    C8 evaluates C function #8.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 June 2014

  Author:

    John Burkardt

  Parameters:

    Input, double X, the evaluation point.

    Output, double C8, the value of C(X).
*/
{
  double value;

  value = 0.0;

  return value;
}
/******************************************************************************/

double exact8 ( double x )

/******************************************************************************/
/*
  Purpose:

    EXACT8 evaluates exact solution #8.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 June 2014

  Author:

    John Burkardt

  Parameters:

    Input, double X, the evaluation point.

    Output, double EXACT8, the value of U(X).
*/
{
  double value;

  if ( x <= 2.0 / 3.0 )
  {
    value = x * ( 1.0 - x ) * exp ( x );
  }
  else
  {
    value = x * ( 1.0 - x ) * exp ( 2.0 / 3.0 );
  }

  return value;
}
/******************************************************************************/

double exact_ux8 ( double x )

/******************************************************************************/
/*
  Purpose:

    EXACT_UX8 evaluates the derivative of exact solution #8.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 June 2014

  Author:

    John Burkardt

  Parameters:

    Input, double X, the evaluation point.

    Output, double EXACT_UX8, the value of dUdX(X).
*/
{
  double value;

  if ( x <= 2.0 / 3.0 )
  {
    value = ( 1.0 - x - x * x ) * exp ( x );
  }
  else
  {
    value = ( 1.0 - 2.0 * x ) * exp ( 2.0 / 3.0 );
  }

  return value;
}
/******************************************************************************/

double f8 ( double x )

/******************************************************************************/
/*
  Purpose:

    F8 evaluates right hand side function #8.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 June 2014

  Author:

    John Burkardt

  Parameters:

    Input, double X, the evaluation point.

    Output, double F8, the value of F(X).
*/
{
  double value;

  if ( x <= 2.0 / 3.0 )
  {
    value = x * ( x + 3.0 ) * exp ( x );
  }
  else
  {
    value = 2.0 * exp ( 2.0 / 3.0 );
  }

  return value;
}
/******************************************************************************/

void test09 ( )

/******************************************************************************/
/*
  Purpose:

    TEST09 carries out test case #9.

  Discussion:

    Use A9, C9, F9, EXACT9, EXACT_UX9.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 June 2014

  Author:

    John Burkardt

  Reference:

    Dianne O'Leary,
    Scientific Computing with Case Studies,
    SIAM, 2008,
    ISBN13: 978-0-898716-66-5,
    LC: QA401.O44.
*/
{
  int i;
  int n = 11;
  double e1;
  double e2;
  double h1s;
  double *u;
  double uexact;
  double *x;
  double x_first;
  double x_last;

  printf ( "\n" );
  printf ( "TEST09\n" );
  printf ( "  Solve -( A(x) U'(x) )' + C(x) U(x) = F(x)\n" );
  printf ( "  for 0 < x < 1, with U(0) = U(1) = 0.\n" );
  printf ( "  A9(X) = 1.0\n" );
  printf ( "  C9(X) = 0.0\n" );
  printf ( "  F9(X) = X * ( X + 3 ) * exp ( X ),   X <= 2/3\n" );
  printf ( "        = 2 * exp ( 2/3),                   2/3 < X\n" );
  printf ( "  U9(X) = X * ( 1 - X ) * exp ( X ),   X <= 2/3\n" );
  printf ( "        = X * ( 1 - X ),                    2/3 < X\n" );
  printf ( "\n" );
  printf ( "  Number of nodes = %d\n", n );
/*
  Geometry definitions.
*/
  x_first = 0.0;
  x_last = 1.0;
  x = r8vec_even_new ( n, x_first, x_last );

  u = fem1d_bvp_quadratic ( n, a9, c9, f9, x );

  printf ( "\n" );
  printf ( "     I    X         U         Uexact    Error\n" );
  printf ( "\n" );

  for ( i = 0; i < n; i++ )
  {
    uexact = exact9 ( x[i] );
    printf ( "  %4d  %8f  %14f  %14f %14e\n", 
      i, x[i], u[i], uexact, fabs ( u[i] - uexact ) );
  }

  e1 = l1_error ( n, x, u, exact9 );
  e2 = l2_error_quadratic ( n, x, u, exact9 );
  h1s = h1s_error_quadratic ( n, x, u, exact_ux9 );
  printf ( "\n" );
  printf ( "  l1 norm of error  = %g\n", e1 );
  printf ( "  L2 norm of error  = %g\n", e2 );
  printf ( "  Seminorm of error = %g\n", h1s );

  free ( u );
  free ( x );

  return;
}
/******************************************************************************/

double a9 ( double x )

/******************************************************************************/
/*
  Purpose:

    A9 evaluates A function #9.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 June 2014

  Author:

    John Burkardt

  Parameters:

    Input, double X, the evaluation point.

    Output, double A9, the value of A(X).
*/
{
  double value;

  value = 1.0;

  return value;
}
/******************************************************************************/

double c9 ( double x )

/******************************************************************************/
/*
  Purpose:

    C9 evaluates C function #9.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 June 2014

  Author:

    John Burkardt

  Parameters:

    Input, double X, the evaluation point.

    Output, double C9, the value of C(X).
*/
{
  double value;

  value = 0.0;

  return value;
}
/******************************************************************************/

double exact9 ( double x )

/******************************************************************************/
/*
  Purpose:

    EXACT9 evaluates exact solution #9.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 June 2014

  Author:

    John Burkardt

  Parameters:

    Input, double X, the evaluation point.

    Output, double EXACT9, the value of U(X).
*/
{
  double value;

  if ( x <= 2.0 / 3.0 )
  {
    value = x * ( 1.0 - x ) * exp ( x );
  }
  else
  {
    value = x * ( 1.0 - x );
  }

  return value;
}
/******************************************************************************/

double exact_ux9 ( double x )

/******************************************************************************/
/*
  Purpose:

    EXACT_UX9 evaluates the derivative of exact solution #9.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 June 2014

  Author:

    John Burkardt

  Parameters:

    Input, double X, the evaluation point.

    Output, double EXACT_UX9, the value of dUdX(X).
*/
{
  double value;

  if ( x <= 2.0 / 3.0 )
  {
    value = ( 1.0 - x - x * x ) * exp ( x );
  }
  else
  {
    value = 1.0 - 2.0 * x;
  }

  return value;
}
/******************************************************************************/

double f9 ( double x )

/******************************************************************************/
/*
  Purpose:

    F9 evaluates right hand side function #9.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 June 2014

  Author:

    John Burkardt

  Parameters:

    Input, double X, the evaluation point.

    Output, double F9, the value of F(X).
*/
{
  double value;

  if ( x <= 2.0 / 3.0 )
  {
    value = x * ( x + 3.0 ) * exp ( x );
  }
  else
  {
    value = 2.0;
  }

  return value;
}
