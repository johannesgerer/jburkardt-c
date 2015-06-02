# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "asa091.h"

int main ( );
void test01 ( );
void test02 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for ASA091_PRB.

  Discussion:

    ASA091_PRB tests the ASA091 library.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    03 November 2010

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "ASA091_PRB:\n" );
  printf ( "  C version\n" );
  printf ( "  Test the ASA091 library.\n" );

  test01 ( );
  test02 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "ASA091_PRB:\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void test01 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST01 makes a single simple calculation with PPCHI2.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    03 November 2010

  Author:

    John Burkardt
*/
{
  double g;
  int ifault;
  double p;
  double v;
  double value;
  double value_correct = 0.4;

  p = 0.017523;
  v = 4.0;

  printf ( "\n" );
  printf ( "TEST01:\n" );
  printf ( "  Perform a simple sample calculation using\n" );
  printf ( "  PPCHI2 to invert the Chi-Squared CDF.\n" );

  g = lgamma ( v / 2.0 );

  printf ( "\n" );
  printf ( "  P =                  %24.16f\n", p );
  printf ( "  V =                  %24.16f\n", v );
  printf ( "  G Log(Gamma(V/2)) =  %24.16f\n", g );

  value = ppchi2 ( p, v, g, &ifault );

  printf ( "  VALUE =              %24.16f\n", value );
  printf ( "  VALUE (correct) =    %24.16f\n", value_correct );

  printf ( "\n" );
  printf ( "  Error flag IFAULT = %d\n", ifault );

  return;
}
/******************************************************************************/

void test02 ( )

/******************************************************************************/
/*
  Purpose:

    TEST02 compares PPCHI2 against tabulated values.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    03 November 2010

  Author:

    John Burkardt
*/
{
  int a;
  double fx;
  double g;
  int ifault;
  int n_data;
  double v;
  double x;
  double x2;

  printf ( "\n" );
  printf ( "TEST02:\n" );
  printf ( "  PPCHI2 computes percentage points of the Chi-Square CDF.\n" );
  printf ( "  Compare to tabulated values.\n" );
  printf ( "\n" );
  printf ( "         N        CDF       X                        " );
  printf ( " X2                    DIFF\n" );
  printf ( "                           (tabulated)               " );
  printf ( "(PPCHI2)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    chi_square_cdf_values ( &n_data, &a, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    v = ( double ) ( a );

    g = lgamma ( v / 2.0 );

    x2 = ppchi2 ( fx, v, g, &ifault );

    printf ( "  %10.4f  %10.4f  %24.16f  %24.16f  %10.4e\n",
      v, fx, x, x2, fabs ( x - x2 ) );
  }

  return;
}
