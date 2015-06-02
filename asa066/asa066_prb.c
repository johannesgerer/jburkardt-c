# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>

# include "asa066.h"

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for ASA066_PRB.

  Discussion:

    ASA066_PRB tests the ASA066 library.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    01 November 2010

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "ASA066_PRB:\n" );
  printf ( "  C version\n" );
  printf ( "  Test the ASA066 library.\n" );

  test01 ( );
  test02 ( );
  test03 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "ASA066_PRB:\n" );
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

    TEST01 compares ALNORM against tabulated values.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    01 November 2010

  Author:

    John Burkardt
*/
{
  double fx;
  double fx2;
  int n_data;
  int upper = 0;
  double x;

  printf ( "\n" );
  printf ( "TEST01:\n" );
  printf ( "  Compare tabulated values of the normal\n" );
  printf ( "  Cumulative Density Function against values\n" );
  printf ( "  computed by ALNORM.\n" );
  printf ( "\n" );
  printf ( "         X        CDF                       CDF" );
  printf ( "                    DIFF\n" );
  printf ( "               (tabulated)                 (ALNORM)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    normal_01_cdf_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = alnorm ( x, upper );

    printf ( "  %10.4f  %24.16f  %24.16f  %20.4e\n",
      x, fx, fx2, r8_abs ( fx - fx2 ) );
  }

  return;
}
/******************************************************************************/

void test02 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST02 compares NORMP against tabulated values.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    01 November 2010

  Author:

    John Burkardt
*/
{
  double fx;
  double fx2;
  int n_data;
  double p;
  double pdf;
  double q;
  double x;

  printf ( "\n" );
  printf ( "TEST02:\n" );
  printf ( "  Compare tabulated values of the normal\n" );
  printf ( "  Cumulative Density Function against values\n" );
  printf ( "  computed by NORMP.\n" );
  printf ( "\n" );
  printf ( "         X        CDF                       CDF" );
  printf ( "                    DIFF\n" );
  printf ( "               (tabulated)                 (NORMP)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    normal_01_cdf_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    normp ( x, &p, &q, &pdf );
    fx2 = p;

    printf ( "  %10.4f  %24.16f  %24.16f  %20.4e\n",
      x, fx, fx2, r8_abs ( fx - fx2 ) );
  }

  return;
}
/******************************************************************************/

void test03 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST03 compares NPROB against tabulated values.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    01 November 2010

  Author:

    John Burkardt
*/
{
  double fx;
  double fx2;
  int n_data;
  double p;
  double pdf;
  double q;
  double x;

  printf ( "\n" );
  printf ( "TEST03\n" );
  printf ( "  Compare tabulated values of the normal\n" );
  printf ( "  Cumulative Density Function against values\n" );
  printf ( "  computed by NPROBP.\n" );
  printf ( "\n" );
  printf ( "         X        CDF                       CDF" );
  printf ( "                    DIFF\n" );
  printf ( "               (tabulated)                 (NPROB)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    normal_01_cdf_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    nprob ( x, &p, &q, &pdf );
    fx2 = p;

    printf ( "  %10.4f  %24.16f  %24.16f  %20.4e\n",
      x, fx, fx2, r8_abs ( fx - fx2 ) );
  }

  return;
}
