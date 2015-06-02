# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>

# include "asa241.h"

int main ( );
void test01 ( );
void test02 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for ASA241_PRB.

  Discussion:

    ASA241_PRB tests the ASA241 library.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    19 March 2010

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "ASA241_PRB:\n" );
  printf ( "  C version\n" );
  printf ( "  Test the ASA241 library.\n" );

  test01 ( );
  test02 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "ASA241_PRB:\n" );
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

    TEST01 tests R4_NORMAL_01_CDF_INVERSE, NORMAL_01_CDF_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    19 March 2010

  Author:

    John Burkardt
*/
{
  double fx;
  float fx2;
  int n_data;
  double x;
  float x2;

  printf ( "\n" );
  printf ( "TEST01:\n" );
  printf ( "  Let FX = NormalCDF ( X ).\n" );
  printf ( "\n" );
  printf ( "  NORMAL_01_CDF_VALUES returns some values of ( X, FX ).\n" );
  printf ( "\n" );
  printf ( "  R4_NORMAL_01_CDF_INVERSE takes the value of FX, and computes an\n" );
  printf ( "    estimate X2, of the corresponding input argument,\n" );
  printf ( "    accurate to about 7 decimal places.\n" );
  printf ( "\n" );
  printf ( "          FX                        X                        X2          DIFF\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    normal_01_cdf_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = ( float ) ( fx );
    x2 = r4_normal_01_cdf_inverse ( fx2 );

    printf ( "  %24.16f  %24.16f  %24.16f  %10.6g\n", fx, x, x2, fabs ( x - x2 ) );
  }

  return;
}
/******************************************************************************/

void test02 ( )

/******************************************************************************/
/*
  Purpose:

    TEST02 tests R8_NORMAL_01_CDF_INVERSE, NORMAL_01_CDF_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    19 March 2010

  Author:

    John Burkardt
*/
{
  double fx;
  double fx2;
  int n_data;
  double x;
  double x2;

  printf ( "\n" );
  printf ( "TEST02:\n" );
  printf ( "  Let FX = NormalCDF ( X ).\n" );
  printf ( "\n" );
  printf ( "  NORMAL_01_CDF_VALUES returns some values of ( X, FX ).\n" );
  printf ( "\n" );
  printf ( "  R8_NORMAL_01_CDF_INVERSE takes the value of FX, and computes an\n" );
  printf ( "    estimate X2, of the corresponding input argument,\n" );
  printf ( "    accurate to about 16 decimal places.\n" );
  printf ( "\n" );
  printf ( "          FX                        X                        X2          DIFF\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    normal_01_cdf_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    x2 = r8_normal_01_cdf_inverse ( fx );

    printf ( "  %24.16f  %24.16f  %24.16f  %10.6g\n", fx, x, x2, fabs ( x - x2 ) );
  }

  return;
}
