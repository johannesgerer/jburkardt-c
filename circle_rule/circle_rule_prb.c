# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "circle_rule.h"

int main ( );
void test01 ( int nt );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for CIRCLE_RULE_PRB.

  Discussion:

    CIRCLE_RULE_PRB tests the CIRCLE_RULE library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    06 April 2014

  Author:

    John Burkardt
*/
{
  int nt;

  timestamp ( );
  printf ( "\n" );
  printf ( "CIRCLE_RULE:\n" );
  printf ( "  C version\n" );
  printf ( "  Test the CIRCLE_RULE library.\n" );

  nt = 8;
  test01 ( nt );

  nt = 32;
  test01 ( nt );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "CIRCLE_RULE:\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void test01 ( int nt )

/******************************************************************************/
/*
  Purpose:

    TEST01 tests CIRCLE_RULE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 April 2014

  Author:

    John Burkardt
*/
{
  int e[2];
  int e1;
  int e2;
  double exact;
  int i;
  double q;
  double r8_pi = 3.141592653589793;
  double *t;
  double *w;
  double x;
  double y;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  CIRCLE_RULE can compute a rule Q(f) for the unit circle\n" );
  printf ( "  using NT equally spaced angles.\n" );
  printf ( "  Estimate integrals I(f) where f = x^e(1) * y^e(2)\n" );
  printf ( "  using %d points.\n", nt );
/*
  Compute the quadrature rule.
*/
  w = ( double * ) malloc ( nt * sizeof ( double ) );
  t = ( double * ) malloc ( nt * sizeof ( double ) );

  circle_rule ( nt, w, t );
/*
  Apply it to integrands.
*/
  printf ( "\n" );
  printf ( "  E(1)  E(2)    I(f)            Q(f)\n" );
  printf ( "\n" );
/*
  Specify a monomial.
*/
  for ( e1 = 0; e1 <= 6; e1 = e1 + 2 )
  {
    e[0] = e1;

    for ( e2 = e1; e2 <= 6; e2 = e2 + 2 )
    {
      e[1] = e2;

      q = 0.0;
      for ( i = 0; i < nt; i++ )
      {
        x = cos ( t[i] );
        y = sin ( t[i] );
        q = q + w[i] * pow ( x, e[0] ) * pow ( y, e[1] );
      }

      q = 2.0 * r8_pi * q;

      exact = circle01_monomial_integral ( e );

      printf ( "   %2d  %2d  %14.6g  %14.6g\n", e[0], e[1], exact, q );
    }
  }

  free ( t );
  free ( w );

  return;
}

