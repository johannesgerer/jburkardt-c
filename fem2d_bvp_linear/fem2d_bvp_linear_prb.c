# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "fem2d_bvp_linear.h"

int main ( );

void test01 ( );
double a1 ( double x, double y );
double c1 ( double x, double y );
double exact1 ( double x, double y );
double exact_ux1 ( double x, double y );
double exact_uy1 ( double x, double y );
double f1 ( double x, double y );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for FEM2D_BVP_LINEAR_PRB.

  Discussion:

    FEM2D_BVP_LINEAR_PRB tests the FEM2D_BVP_LINEAR library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 June 2014

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "FEM2D_BVP_LINEAR_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the FEM2D_BVP_LINEAR library.\n" );

  test01 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "FEM2D_BVP_LINEAR_PRB\n" );
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

    Use A1, C1, F1, EXACT1, EXACT_UX1, EXACT_UX2.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 June 2014

  Author:

    John Burkardt
*/
{
  double e1;
  double e2;
  double h1s;
  int i;
  int j;
  int mn;
  int nx = 3;
  int ny = 3;
  double *u;
  double uexact;
  double *x;
  double x_first;
  double x_last;
  double *y;
  double y_first;
  double y_last;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  Solve - del ( A del U ) + C U = F \n" );
  printf ( "  on the unit square with zero boundary conditions.\n" );
  printf ( "  A1(X,Y) = 1.0\n" );
  printf ( "  C1(X,Y) = 0.0\n" );
  printf ( "  F1(X,Y) = 2*X*(1-X)+2*Y*(1-Y)\n" );
  printf ( "  U1(X,Y) = X * ( 1 - X ) * Y * ( 1 - Y )\n" );
  printf ( "\n" );
  printf ( "  Number of X grid values NX = %d\n", nx );
  printf ( "  Number of Y grid values NY = %d\n", ny );
  mn = nx * ny;
/*
  Geometry definitions.
*/
  x_first = 0.0;
  x_last = 1.0;
  x = r8vec_even_new ( nx, x_first, x_last );

  y_first = 0.0;
  y_last = 1.0;
  y = r8vec_even_new ( ny, y_first, y_last );

  u = fem2d_bvp_linear ( nx, ny, a1, c1, f1, x, y );

  printf ( "\n" );
  printf ( "     I     J    X         Y         U         Uexact    Error\n" );
  printf ( "\n" );

  for ( j = 0; j < ny; j++ )
  {
    for ( i = 0; i < nx; i++ )
    {
      uexact = exact1 ( x[i], y[j] );
      printf ( "  %4d  %4d  %8f  %8f  %14f  %14f %14e\n", 
        i, j, x[i], y[j], u[i+j*nx], uexact, fabs ( u[i+j*nx] - uexact ) );
    }
  }

  e1 = fem2d_l1_error ( nx, ny, x, y, u, exact1 );
  e2 = fem2d_l2_error_linear ( nx, ny, x, y, u, exact1 );
  h1s = fem2d_h1s_error_linear ( nx, ny, x, y, u, exact_ux1, exact_uy1 );
  printf ( "\n" );
  printf ( "  l1 norm of error  = %g\n", e1 );
  printf ( "  L2 norm of error  = %g\n", e2 );
  printf ( "  Seminorm of error = %g\n", h1s );

  free ( u );
  free ( x );
  free ( y );

  return;
}
/******************************************************************************/

double a1 ( double x, double y )

/******************************************************************************/
/*
  Purpose:

    A1 evaluates A function #1.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 June 2014

  Author:

    John Burkardt

  Parameters:

    Input, double X, Y, the evaluation point.

    Output, double A1, the value of A(X,Y).
*/
{
  double value;

  value = 1.0;

  return value;
}
/******************************************************************************/

double c1 ( double x, double y )

/******************************************************************************/
/*
  Purpose:

    C1 evaluates C function #1.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 June 2014

  Author:

    John Burkardt

  Parameters:

    Input, double X, Y, the evaluation point.

    Output, double C1, the value of C(X,Y).
*/
{
  double value;

  value = 0.0;

  return value;
}
/******************************************************************************/

double exact1 ( double x, double y )

/******************************************************************************/
/*
  Purpose:

    EXACT1 evaluates exact solution #1.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 June 2014

  Author:

    John Burkardt

  Parameters:

    Input, double X, Y, the evaluation point.

    Output, double EXACT1, the value of U(X,Y).
*/
{
  double value;

  value = x * ( 1.0 - x ) * y * ( 1.0 - y );

  return value;
}
/******************************************************************************/

double exact_ux1 ( double x, double y )

/******************************************************************************/
/*
  Purpose:

    EXACT_UX1 evaluates the derivative of exact solution #1.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 June 2014

  Author:

    John Burkardt

  Parameters:

    Input, double X, Y, the evaluation point.

    Output, double EXACT_UX1, the value of dUdX(X,Y).
*/
{
  double value;

  value = ( 1.0 - 2.0 * x ) * ( y - y * y );

  return value;
}
/******************************************************************************/

double exact_uy1 ( double x, double y )

/******************************************************************************/
/*
  Purpose:

    EXACT_UY1 evaluates the derivative of exact solution #1.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 June 2014

  Author:

    John Burkardt

  Parameters:

    Input, double X, Y, the evaluation point.

    Output, double EXACT_UY1, the value of dUdY(X,Y).
*/
{
  double value;

  value = ( x - x * x ) * ( 1.0 - 2.0 * y );

  return value;
}
/******************************************************************************/

double f1 ( double x, double y )

/******************************************************************************/
/*
  Purpose:

    F1 evaluates right hand side function #1.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 June 2014

  Author:

    John Burkardt

  Parameters:

    Input, double X, Y, the evaluation point.

    Output, double F1, the value of F(X,Y).
*/
{
  double value;

  value = 2.0 * x * ( 1.0 - x ) 
        + 2.0 * y * ( 1.0 - y );

  return value;
}

