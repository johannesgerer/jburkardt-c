# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "asa245.h"

int main ( void );
void test01 ( void );
void test02 ( void );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for ASA245_PRB.

  Discussion:

    ASA245_PRB calls the ASA245 routines.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    20 November 2010

  Author:

    John Burkardt
*/
{
  timestamp ( );

  printf ( "\n" );
  printf ( "ASA245_PRB:\n" );
  printf ( "  C version\n" );
  printf ( "  Test the ASA245 library.\n" );

  test01 ( );
  test02 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "ASA245_PRB:\n" );
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

    TEST01 demonstrates the use of ALNGAM.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    20 November 2010

  Author:

    John Burkardt
*/
{
  double fx;
  double fx2;
  int ifault;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST01:\n" );
  printf ( "  ALNGAM computes the logarithm of the \n" );
  printf ( "  Gamma function.  We compare the result\n" );
  printf ( "  to tabulated values.\n" );
  printf ( "\n" );
  printf ( "          X                     " );
  printf ( "FX                        FX2\n" );
  printf ( "                                " );
  printf ( "(Tabulated)               (ALNGAM)                DIFF\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    gamma_log_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = alngam ( x, &ifault );

    printf ( "  %24.16f  %24.16f  %24.16f  %10.4e\n",
      x, fx, fx2, fabs ( fx - fx2 ) );
  }

  return;
}
/******************************************************************************/

void test02 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST02 demonstrates the use of LNGAMMA.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    20 November 2010

  Author:

    John Burkardt
*/
{
  double fx;
  double fx2;
  int ier;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST02:\n" );
  printf ( "  LNGAMMA computes the logarithm of the \n" );
  printf ( "  Gamma function.  We compare the result\n" );
  printf ( "  to tabulated values.\n" );
  printf ( "\n" );
  printf ( "          X                     " );
  printf ( "FX                        FX2\n" );
  printf ( "                                " );
  printf ( "(Tabulated)               (LNGAMMA)                DIFF\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    gamma_log_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = lngamma ( x, &ier );

    printf ( "  %24.16f  %24.16f  %24.16f  %10.4e\n",
      x, fx, fx2, fabs ( fx - fx2 ) );
  }

  return;
}
