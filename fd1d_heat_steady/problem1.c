# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "fd1d_heat_steady.h"

int main ( void );
double k1 ( double x );
double f1 ( double x );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for problem 1.

 Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 May 2009

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  int header;
  int n = 11;
  double *u;
  char *u_file = "problem1_values.txt";
  double ua;
  double ub;
  double *x;
  char *x_file = "problem1_nodes.txt";

  timestamp ( );
  printf ( "\n" );
  printf ( "PROBLEM1:\n" );
  printf ( "  C version\n" );
  printf ( "  A test problem for FD1D_HEAT_STEADY.\n" );

  a = 0.0;
  b = 1.0;

  x = r8vec_even ( n, a, b );

  ua = 0.0;
  ub = 1.0;

  u = fd1d_heat_steady ( n, a, b, ua, ub, k1, f1, x );

  r8mat_write ( x_file, 1, n, x );

  printf ( "\n" );
  printf ( "  X data written to \"%s\".\n", x_file );

  r8mat_write ( u_file, 1, n, u );

  printf ( "  U data written to \"%s\".\n", u_file );

  free ( u );
  free ( x );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "PROBLEM1:\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

double k1 ( double x )

/******************************************************************************/
/*
  Purpose:

    K1 evaluates the heat transfer coefficient K(X).

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 May 2009

  Author:

    John Burkardt

  Parameters:

    Input, double X, the position.

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

    F1 evaluates the right hand side of the steady state heat equation.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 May 2009

  Author:

    John Burkardt

  Parameters:

    Input, double X, the position.

    Output, double F1, the value of F(X).
*/
{
  double value;

  value = 0.0;

  return value;
}


