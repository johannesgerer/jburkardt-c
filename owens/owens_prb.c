# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "owens.h"

int main ( void );
void test01 ( void );
void test02 ( void );
void test03 ( void );
void test04 ( void );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for OWENS_PRB.

  Discussion:

    OWENS_PRB tests the OWENS library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    23 November 2010

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "OWENS_PRB:\n" );
  printf ( "  C version\n" );
  printf ( "  Test the OWENS library.\n" );

  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "OWENS_PRB:\n" );
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

    TEST01 demonstrates the use of T.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    21 November 2010

  Author:

    John Burkardt
*/
{
  double a;
  double h;
  int n_data;
  double t1;
  double t2;

  printf ( "\n" );
  printf ( "TEST01:\n" );
  printf ( "  T computes Owen's T function.\n" );
  printf ( "  Compare to tabulated values.\n" );
  printf ( "\n" );
  printf ( "             H             A      " );
  printf ( "    T                         T  \n" );
  printf ( "                                  " );
  printf ( "    (Tabulated)               (TFN)               DIFF\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    owen_values ( &n_data, &h, &a, &t1 );

    if ( n_data == 0 )
    {
      break;
    }

    t2 = t ( h, a );

    printf ( "  %12.4f  %12.4f  %24.16f  %24.16f  %10.4g\n",
      h, a, t1, t2, r8_abs ( t1 - t2 ) );
  }

  return;
}
/******************************************************************************/

void test02 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST02 demonstrates the use of BIVNOR.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 April 2012

  Author:

    John Burkardt
*/
{
  double fxy1;
  double fxy2;
  int n_data;
  double r;
  double x;
  double y;

  printf ( "\n" );
  printf ( "TEST02:\n" );
  printf ( "  BIVNOR computes the bivariate normal probability.\n" );
  printf ( "  Compare to tabulated values.\n" );
  printf ( "\n" );
  printf ( "          X               Y               " );
  printf ( "R           P                         P                       DIFF\n" );
  printf ( "                                          " );
  printf ( "           (Tabulated)               (BIVNOR)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    bivariate_normal_cdf_values ( &n_data, &x, &y, &r, &fxy1 );

    if ( n_data == 0 )
    {
      break;
    }

    fxy2 = bivnor ( - x, - y, r );

    printf ( "  %12.4f  %12.4f  %12.4f  %24.16f  %24.16f  %10.4g\n",
      x, y, r, fxy1, fxy2, r8_abs ( fxy1 - fxy2 ) );
  }

  return;
}
/******************************************************************************/

void test03 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST03 demonstrates the use of ZNORM1.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    21 November 2010

  Author:

    John Burkardt
*/
{
  double fx1;
  double fx2;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST03:\n" );
  printf ( "  ZNORM1 computes the normal CDF starting at 0.\n" );
  printf ( "  Compare to tabulated values.\n" );
  printf ( "\n" );
  printf ( "          X           P                         P                       DIFF\n" );
  printf ( "                     (Tabulated)               (ZNORM1)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    normal_01_cdf_values ( &n_data, &x, &fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx1 = fx1 - 0.5;

    fx2 = znorm1 ( x );

    printf ( "  %12.4f  %24.16f  %24.16f  %10.4g\n",
      x, fx1, fx2, r8_abs ( fx1 - fx2 ) );
  }

  return;
}
/******************************************************************************/

void test04 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST04 demonstrates the use of ZNORM2.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    23 November 2010

  Author:

    John Burkardt
*/
{
  double fx1;
  double fx2;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST04:\n" );
  printf ( "  ZNORM2 computes the complementary normal CDF.\n" );
  printf ( "  Compare to tabulated values.\n" );
  printf ( "\n" );
  printf ( "          X           P                         P                       DIFF\n" );
  printf ( "                     (Tabulated)               (ZNORM2)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    normal_01_cdf_values ( &n_data, &x, &fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx1 = 1.0 - fx1;

    fx2 = znorm2 ( x );

    printf ( "  %12.4f  %24.16f  %24.16f  %10.4g\n",
      x, fx1, fx2, r8_abs ( fx1 - fx2 ) );
  }

  return;
}

