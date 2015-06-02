# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <string.h>

# include "zero_rc.h"

int main ( );
void test01 ( );
void test02 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for ZERO_RC_PRB.

  Discussion:

    ZERO_RC_PRB tests the ZERO_RC library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    25 January 2013

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "ZERO_RC_PRB:\n" );
  printf ( "  C version\n" );
  printf ( "  Test the ZERO_RC library.\n" );

  test01 ( );
  test02 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "ZERO_RC_PRB:\n" );
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

    TEST01 tests ROOT_RC.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    22 January 2013

  Author:

    John Burkardt
*/
{
  double ferr;
  double fx;
  int i;
  int it;
  int it_max;
  double q[9];
  double x;
  double xerr;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  Test ROOT_RC, which searches for an\n" );
  printf ( "  approximate solution of F(X) = 0.\n" );
  printf ( "\n" );
  printf ( "       X              XERR            FX              FERR\n" );
  printf ( "\n" );
/*
  Initialization.
*/
  it = 0;
  it_max = 30;
  for ( i = 0; i < 9; i++ )
  {
    q[i] = 0.0;
  }
  x = - 2.1;
/*
  Each call takes one more step of improvement.
*/
  for ( ; ; )
  {
    fx = cos ( x ) - x;

    if ( it == 0 )
    {
      printf ( "  %14.6g                  %14.6g\n", x, fx );
    }
    else
    {
      printf ( "  %14.6g  %14.6g  %14.6g  %14.6g\n", x, xerr, fx, ferr );
    }

    x = root_rc ( x, fx, &ferr, &xerr, q );

    if ( ferr < 1.0E-08 )
    {
      printf ( "\n" );
      printf ( "  Uncertainty in F(X) less than tolerance\n" );
      break;
    }

    if ( xerr < 1.0E-08 )
    {
      printf ( "\n" );
      printf ( "  Width of X interal less than tolerance.\n" );
      break;
    }

    if ( it_max < it )
    {
      printf ( "\n" );
      printf ( "  Too many iterations!'\n" );
      break;
    }
    it = it + 1;     
  }

  return;
}
/******************************************************************************/

void test02 ( )

/******************************************************************************/
/*
  Purpose:

    TEST02 tests ROOTS_RC.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    25 January 2013

  Author:

    John Burkardt
*/
{
  double ferr;
  double *fx;
  int i;
  int it;
  int it_max = 30;
  int j;
  int n = 4;
  double *q;
  double *x;
  double *xnew;

  fx = ( double * ) malloc ( n * sizeof ( double ) );
  q = ( double * ) malloc ( ( 2 * n + 2 ) * ( n + 2 ) * sizeof ( double ) );
  x = ( double * ) malloc ( n * sizeof ( double ) );
  xnew = ( double * ) malloc ( n * sizeof ( double ) );

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  Test ROOTS_RC, which seeks a solution of\n" );
  printf ( "  the N-dimensional nonlinear system F(X) = 0.\n" );
  printf ( "\n" );
  printf ( "       FERR          X\n" );
  printf ( "\n" );
/*
  Initialization.
*/
  for ( j = 0; j < n + 2; j++ )
  {
    for ( i = 0; i < 2 * n + 2; i++ )
    {
      q[i+j*(2*n+2)] = 0.0;
    }
  }

  xnew[0] = 1.2;
  for ( i = 1; i < n; i++ )
  {
    xnew[i] = 1.0;
  }

  it = 0;

  for ( ; ; )
  {
    for ( i = 0; i < n; i++ )
    {
      x[i] = xnew[i];
    }

    fx[0] = 1.0 - x[0];
    for ( i = 1; i < n; i++ )
    {
      fx[i] = 10.0 * ( x[i] - x[i-1] * x[i-1] );
    }

    if ( it == 0 )
    {
      printf ( "                " );
    }
    else
    {
      printf ( "  %14.6g", ferr );
    }
    for ( i = 0; i < n; i++ )
    {
      printf ( "  %14.6g", x[i] );
    }
    printf ( "\n" );

    roots_rc ( n, x, fx, &ferr, xnew, q );

    if ( ferr < 1.0E-07 )
    {
      printf ( "\n" );
      printf ( "  Sum of |f(x)| less than tolerance.\n" );
      break;
    }

    if ( it_max < it )
    {
      printf ( "\n" );
      printf ( "  Too many iterations!\n" );
      break;
    }
    it = it + 1;
  }

  free ( fx );
  free ( q );
  free ( x );
  free ( xnew );

  return;
}
