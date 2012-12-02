# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>

# include "bisection_integer.h"

int main ( );
void test01 ( );
int f01 ( int x );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for BISECTION_INTEGER_PRB.

  Discussion:

    BISECTION_INTEGER_PRB tests the BISECTION_INTEGER library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    23 August 2012

  Author:

    John Burkardt
*/
{
  timestamp ( );

  printf ( "\n" );
  printf ( "BISECTION_INTEGER_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the BISECTION_INTEGER library.\n" );

  test01 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "BISECTION_INTEGER_PRB\n" );
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

    TEST01 tests BISECTION_INTEGER;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    23 August 2012

  Author:

    John Burkardt
*/
{
  int a;
  int b;
  int c;
  int fc;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  BISECTION_INTEGER attempts to locate an integer root C\n" );
  printf ( "  of an equation F(C) = 0.\n" );
  printf ( "  The user supplies a change of sign interval [A,B].\n" );
  printf ( "  The function considered here has two real roots\n" );
  printf ( "  as well as an integer root, so the algorithm can\n" );
  printf ( "  fail depending on how the change of sign interval is chosen.\n" );

  a = 4;
  b = 100;

  printf ( "\n" );
  printf ( "  The initial change of sign interval is:\n" );
  printf ( "  F(%d) = %d\n", a, f01 ( a ) );
  printf ( "  F(%d) = %d\n", b, f01 ( b ) );

  bisection_integer ( f01, &a, &b, &c, &fc );

  if ( fc == 0 )
  {
    printf ( "\n" );
    printf ( "  An exact root was found at C = %d\n", c );
  }
  else
  {
    printf ( "\n" );
    printf ( "  An exact root was NOT found.\n" );
    printf ( "  The change of sign interval is now:\n" );
    printf ( "  F(%d) = %d\n", a, f01 ( a ) );
    printf ( "  F(%d) = %d\n", b, f01 ( b ) );
  }

  a = -10;
  b = 15;

  printf ( "\n" );
  printf ( "  The initial change of sign interval is:\n" );
  printf ( "  F(%d) = %d\n", a, f01 ( a ) );
  printf ( "  F(%d) = %d\n", b, f01 ( b ) );

  bisection_integer ( f01, &a, &b, &c, &fc );

  if ( fc == 0 )
  {
    printf ( "\n" );
    printf ( "  An exact root was found at C = %d\n", c );
  }
  else
  {
    printf ( "\n" );
    printf ( "  An exact root was NOT found.\n" );
    printf ( "  The change of sign interval is now:\n" );
    printf ( "  F(%d) = %d\n", a, f01 ( a ) );
    printf ( "  F(%d) = %d\n", b, f01 ( b ) );
  }

  return;
}
/******************************************************************************/

int f01 ( int n )

/******************************************************************************/
/*
  Purpose:

    F01 is a test function.

  Discussion:

    The polynomial has roots 1/2, 7/2, and 10.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    23 August 2012

  Author:

    John Burkardt

  Parameters:

    Input, int N, the argument.

    Output, int F01, the function value.
*/
{
  int value;

  value = ( 2 * n - 7 ) * ( 2 * n - 1 ) * ( n - 10 );

  return value;
}
