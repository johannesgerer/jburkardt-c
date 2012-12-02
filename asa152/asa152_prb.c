# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "asa152.h"

int main ( void );
void test01 ( void );
void test02 ( void );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for ASA152_PRB.

  Discussion:

    ASA152_PRB calls the ASA152 routines.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 November 2010

  Author:

    John Burkardt
*/
{
  timestamp ( );

  printf ( "\n" );
  printf ( "ASA152_PRB:\n" );
  printf ( "  C version\n" );
  printf ( "  Test the ASA152 library.\n" );

  test01 ( );
  test02 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "ASA152_PRB:\n" );
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

    TEST01 demonstrates CHYPER for cumulative probabilities.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 November 2010

  Author:

    John Burkardt
*/
{
  double fx;
  double fx2;
  int ifault;
  int n_data;
  int point;
  int pop;
  int sam;
  int suc;
  int x;

  printf ( "\n" );
  printf ( "TEST01:\n" );
  printf ( "  CHYPER computes cumulative probablities\n" );
  printf ( "  of the hypergeometric PDF.\n" );
  printf ( "  Compare to tabulated values.\n" );
  printf ( "\n" );
  printf ( "   SAM   SUC   POP     X    " );
  printf ( "  CDF                       CDF                     DIFF\n" );
  printf ( "                            " );
  printf ( " (tabulated)               (CHYPER)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    hypergeometric_cdf_values ( &n_data, &sam, &suc, &pop, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    point = 0;

    fx2 = chyper ( point, sam, x, pop, suc, &ifault );

    printf ( "  %4d  %4d  %4d  %4d  %24.16f  %24.16f  %10.4e\n",
      sam, suc, pop, x, fx, fx2, r8_abs ( fx - fx2 ) );
  }

  return;
}
/******************************************************************************/

void test02 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST02 demonstrates CHYPER for point probabilities.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 November 2010

  Author:

    John Burkardt
*/
{
  double fx;
  double fx2;
  int ifault;
  int n_data;
  int point;
  int pop;
  int sam;
  int suc;
  int x;

  printf ( "\n" );
  printf ( "TEST02:\n" );
  printf ( "  CHYPER computes point probablities\n" );
  printf ( "  of the hypergeometric PDF.\n" );
  printf ( "  Compare to tabulated values.\n" );
  printf ( "\n" );
  printf ( "   SAM   SUC   POP     X    " );
  printf ( "  PDF                       PDF                     DIFF\n" );
  printf ( "                            " );
  printf ( " (tabulated)               (CHYPER)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    hypergeometric_pdf_values ( &n_data, &sam, &suc, &pop, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    point = 1;

    fx2 = chyper ( point, sam, x, pop, suc, &ifault );

    printf ( "  %4d  %4d  %4d  %4d  %24.16f  %24.16f  %10.4e\n",
      sam, suc, pop, x, fx, fx2, r8_abs ( fx - fx2 ) );
  }

  return;
}
