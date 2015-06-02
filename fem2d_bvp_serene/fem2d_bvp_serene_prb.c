# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "fem2d_bvp_serene.h"

int main ( );

void test01 ( );
double a1 ( double x, double y );
double c1 ( double x, double y );
double exact1 ( double x, double y );
double exact_ux1 ( double x, double y );
double exact_uy1 ( double x, double y );
double f1 ( double x, double y );
void test02 ( );
void test03 ( );
double a3 ( double x, double y );
double c3 ( double x, double y );
double exact3 ( double x, double y );
double exact_ux3 ( double x, double y );
double exact_uy3 ( double x, double y );
double f3 ( double x, double y );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for FEM2D_BVP_SERENE_PRB.

  Discussion:

    FEM2D_BVP_SERENE_PRB tests the FEM2D_BVP_SERENE library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 July 2014

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "FEM2D_BVP_SERENE_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the FEM2D_BVP_SERENE library.\n" );

  test01 ( );
  test02 ( );
  test03 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "FEM2D_BVP_SERENE_PRB\n" );
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

    Use A1, C1, F1, EXACT1, EXACT_UX1, EXACT_UX1.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 July 2014

  Author:

    John Burkardt
*/
{
  double e1;
  double e2;
  double h1s;
  int i;
  int inc;
  int j;
  int k;
  int mn;
  int nx = 5;
  int ny = 5;
  int show11;
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
  x = r8vec_linspace_new ( nx, x_first, x_last );

  y_first = 0.0;
  y_last = 1.0;
  y = r8vec_linspace_new ( ny, y_first, y_last );

  show11 = 0;
  u = fem2d_bvp_serene ( nx, ny, a1, c1, f1, x, y, show11 );

  if ( nx * ny <= 25 )
  {
    printf ( "\n" );
    printf ( "     I     J    X         Y               U               Uexact     Error\n" );
    printf ( "\n" );

    k = 0;

    for ( j = 0; j < ny; j++ )
    {
      if ( ( j % 2 ) == 0 )
      {
        inc = 1;
      }
      else
      {
        inc = 2;
      }
      for ( i = 0; i < nx; i = i + inc )
      {
        uexact = exact1 ( x[i], y[j] );
        printf ( "  %4d  %4d  %8f  %8f  %14f  %14f %14e\n", 
          i, j, x[i], y[j], u[k], uexact, fabs ( u[k] - uexact ) );
        k = k + 1;
      }
    }
  }

  e1 = fem2d_l1_error_serene ( nx, ny, x, y, u, exact1 );
  e2 = fem2d_l2_error_serene ( nx, ny, x, y, u, exact1 );
  h1s = fem2d_h1s_error_serene ( nx, ny, x, y, u, exact_ux1, exact_uy1 );

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

    07 July 2014

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

    07 July 2014

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

    07 July 2014

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

    07 July 2014

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

    07 July 2014

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

    07 July 2014

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
/******************************************************************************/

void test02 ( )

/******************************************************************************/
/*
  Purpose:

    TEST02 checks the basis functions.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 July 2014

  Author:

    John Burkardt
*/
{
  int i;
  int j;
  int seed;
  double *v;
  double *vx;
  double *vy;
  double xe;
  double xq;
  double xw;
  double xx[8] = { 2.0, 1.0, 0.0, 0.0, 0.0, 1.0, 2.0, 2.0 };
  double yn;
  double yq;
  double ys;
  double yy[8] = { 5.0, 5.0, 5.0, 4.0, 3.0, 3.0, 3.0, 4.0 };

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  Basis function checks.\n" );
/*
  Check that V is identity matrix at nodes.
*/
  printf ( "\n" );
  printf ( "  The matrix Aij = V(j)(X(i),Y(i)) should be the identity.\n" );
  printf ( "\n" );

  xw = 0.0;
  ys = 3.0;
  xe = 2.0;
  yn = 5.0;

  for ( i = 0; i < 8; i++ )
  {
    xq = xx[i];
    yq = yy[i];
    v = basis_serene ( xq, yq, xw, ys, xe, yn, xx, yy );
    for ( j = 0; j < 8; j++ )
    {
      printf ( "  %8.4f", v[j] );
    }
    printf ( "\n" );
    free ( v );
  }
/*
  Check that VX and VY sum to zero anywhere.
*/
  printf ( "\n" );
  printf ( "  The vectors dVdX(1:8)(X,Y) and dVdY(1:8)(X,Y)\n" );
  printf ( "  should both sum to zero for any (X,Y).\n" );

  seed = 123456789;
  xq = 2.0 * r8_uniform_01 ( &seed );
  yq = 3.0 + 2.0 * r8_uniform_01 ( &seed );
  xw = 0.0;
  ys = 3.0;
  xe = 2.0;
  yn = 5.0;

  vx = basis_dx_serene ( xq, yq, xw, ys, xe, yn, xx, yy );
  vy = basis_dy_serene ( xq, yq, xw, ys, xe, yn, xx, yy );

  printf ( "\n" );
  printf ( "  Random evaluation point is (%g,%g)\n", xq, yq );
  printf ( "\n" );
  printf ( "              dVdX        dVdY\n" );
  printf ( "\n" );
  for ( i = 0; i < 8; i++ )
  {
    printf ( "  %4d  %10.4g  %10.4g\n", i, vx[i], vy[i] );
  }
  printf ( "\n" );
  printf ( "  Sum:  %10.4g  %10.4g\n", r8vec_sum ( 8, vx ), r8vec_sum ( 8, vy ) );

  free ( vx );
  free ( vy );

  return;
}
/******************************************************************************/

void test03 ( )

/******************************************************************************/
/*
  Purpose:

    TEST03 carries out test case #3.

  Discussion:

    Use A3, C3, F3, EXACT3, EXACT_UX3, EXACT_UX3.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 July 2014

  Author:

    John Burkardt
*/
{
  double *amat;
  double e1;
  double e2;
  double h1s;
  int i;
  int inc;
  int j;
  int k;
  int mn;
  int nx = 5;
  int ny = 5;
  double scale;
  int show11;
  double *u;
  double uexact;
  double *x;
  double x_first;
  double x_last;
  double *y;
  double y_first;
  double y_last;

  printf ( "\n" );
  printf ( "TEST03\n" );
  printf ( "  Solve - del ( A del U ) + C U = F \n" );
  printf ( "  on the unit square with zero boundary conditions.\n" );
  printf ( "  A1(X,Y) = 0.0\n" );
  printf ( "  C1(X,Y) = 1.0\n" );
  printf ( "  F1(X,Y) = X * ( 1 - X ) * Y * ( 1 - Y )\n" );
  printf ( "  U1(X,Y) = X * ( 1 - X ) * Y * ( 1 - Y )\n" );
  printf ( "\n" );
  printf ( "  This example is contrived so that the system matrix\n" );
  printf ( "  is the WATHEN matrix.\n" );
  printf ( "\n" );
  printf ( "  Number of X grid values NX = %d\n", nx );
  printf ( "  Number of Y grid values NY = %d\n", ny );
  mn = nx * ny;
/*
  Geometry definitions.
*/
  x_first = 0.0;
  x_last = 1.0;
  x = r8vec_linspace_new ( nx, x_first, x_last );

  y_first = 0.0;
  y_last = 1.0;
  y = r8vec_linspace_new ( ny, y_first, y_last );

  show11 = 1;
  u = fem2d_bvp_serene ( nx, ny, a3, c3, f3, x, y, show11 );

  if ( nx * ny <= 25 )
  {
    printf ( "\n" );
    printf ( "     I     J    X         Y               U               Uexact     Error\n" );
    printf ( "\n" );

    k = 0;

    for ( j = 0; j < ny; j++ )
    {
      if ( ( j % 2 ) == 0 )
      {
        inc = 1;
      }
      else
      {
        inc = 2;
      }
      for ( i = 0; i < nx; i = i + inc )
      {
        uexact = exact3 ( x[i], y[j] );
        printf ( "  %4d  %4d  %8f  %8f  %14f  %14f %14e\n", 
          i, j, x[i], y[j], u[k], uexact, fabs ( u[k] - uexact ) );
        k = k + 1;
      }
    }
  }

  e1 = fem2d_l1_error_serene ( nx, ny, x, y, u, exact3 );
  e2 = fem2d_l2_error_serene ( nx, ny, x, y, u, exact3 );
  h1s = fem2d_h1s_error_serene ( nx, ny, x, y, u, exact_ux3, exact_uy3 );

  printf ( "\n" );
  printf ( "  l1 norm of error  = %g\n", e1 );
  printf ( "  L2 norm of error  = %g\n", e2 );
  printf ( "  Seminorm of error = %g\n", h1s );

  free ( u );
  free ( x );
  free ( y );
/*
  Pull out the Wathen matrix from MATLAB.
  It will have been multiplied by a random scale factor.
  While my numbering scheme is
    3  2  1
    4     8
    5  6  7
  the numbering scheme used here is 
    1  2  3
    4     5
    6  7  8
*/    
  amat = wathen ( 1, 1, 8 );
 
  scale = 0.5 * amat[0+2*8];
  for ( j = 0; j < 8; j++ )
  {
    for ( i = 0; i < 8; i++ )
    {
      amat[i+j*8] = amat[i+j*8] / scale;
    }
  }
 
  r8mat_print ( 8, 8, amat, "  Wathen matrix:" );

  free ( amat );

  return;
}
/******************************************************************************/

double a3 ( double x, double y )

/******************************************************************************/
/*
  Purpose:

    A3 evaluates A function #3.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 July 2014

  Author:

    John Burkardt

  Parameters:

    Input, double X, Y, the evaluation point.

    Output, double A3, the value of A(X,Y).
*/
{
  double value;

  value = 0.0;

  return value;
}
/******************************************************************************/

double c3 ( double x, double y )

/******************************************************************************/
/*
  Purpose:

    C3 evaluates C function #3.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 July 2014

  Author:

    John Burkardt

  Parameters:

    Input, double X, Y, the evaluation point.

    Output, double C3, the value of C(X,Y).
*/
{
  double value;

  value = 1.0;

  return value;
}
/******************************************************************************/

double exact3 ( double x, double y )

/******************************************************************************/
/*
  Purpose:

    EXACT3 evaluates exact solution #3.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 July 2014

  Author:

    John Burkardt

  Parameters:

    Input, double X, Y, the evaluation point.

    Output, double EXACT3, the value of U(X,Y).
*/
{
  double value;

  value = x * ( 1.0 - x ) * y * ( 1.0 - y );

  return value;
}
/******************************************************************************/

double exact_ux3 ( double x, double y )

/******************************************************************************/
/*
  Purpose:

    EXACT_UX3 evaluates the derivative of exact solution #3.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 July 2014

  Author:

    John Burkardt

  Parameters:

    Input, double X, Y, the evaluation point.

    Output, double EXACT_UX3, the value of dUdX(X,Y).
*/
{
  double value;

  value = ( 1.0 - 2.0 * x ) * ( y - y * y );

  return value;
}
/******************************************************************************/

double exact_uy3 ( double x, double y )

/******************************************************************************/
/*
  Purpose:

    EXACT_UY3 evaluates the derivative of exact solution #3.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 July 2014

  Author:

    John Burkardt

  Parameters:

    Input, double X, Y, the evaluation point.

    Output, double EXACT_UY3, the value of dUdY(X,Y).
*/
{
  double value;

  value = ( x - x * x ) * ( 1.0 - 2.0 * y );

  return value;
}
/******************************************************************************/

double f3 ( double x, double y )

/******************************************************************************/
/*
  Purpose:

    F3 evaluates right hand side function #3.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 July 2014

  Author:

    John Burkardt

  Parameters:

    Input, double X, Y, the evaluation point.

    Output, double F3, the value of F(X,Y).
*/
{
  double value;

  value = x * ( 1.0 - x ) * y * ( 1.0 - y );

  return value;
}
