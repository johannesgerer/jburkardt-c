# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "fd1d_heat_steady.h"

int main ( void );
double k2 ( double x );
double f2 ( double x );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for problem 2.

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
  char *u_file = "problem2_values.txt";
  double ua;
  double ub;
  double *x;
  char *x_file = "problem2_nodes.txt";

  timestamp ( );
  printf ( "\n" );
  printf ( "PROBLEM2:\n" );
  printf ( "  C version\n" );
  printf ( "  A test problem for FD1D_HEAT_STEADY.\n" );
  printf ( "  Low K, then high K, then moderate K.\n" );

  a = 0.0;
  b = 1.0;

  x = r8vec_even ( n, a, b );

  ua = 0.0;
  ub = 1.0;

  u = fd1d_heat_steady ( n, a, b, ua, ub, k2, f2, x );

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
  printf ( "PROBLEM2\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

double k2 ( double x )

/******************************************************************************/
/*
  Purpose:

    K2 evaluates the heat transfer coefficient K(X).

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 May 2009

  Author:

    John Burkardt

  Parameters:

    Input, double X, the position.

    Output, double K2, the value of K(X).
*/
{
  double value;

  if ( x < 0.5 )
  {
    value = 0.25;
  }
  else if ( x < 0.75 )
  {
    value = 4.0;
  }
  else
  {
    value = 1.0;
  }

  return value;
}
/******************************************************************************/

double f2 ( double x )

/******************************************************************************/
/*
  Purpose:

    F2 evaluates the right hand side of the steady state heat equation.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 May 2009

  Author:

    John Burkardt

  Parameters:

    Input, double X, the position.

    Output, double F2, the value of F(X).
*/
{
  double value;

  value = 0.0;

  return value;
}


