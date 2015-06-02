# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "asa076.h"

int main ( );
void test01 ( );
void test02 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for ASA076_PRB.

  Discussion:

    ASA076_PRB tests the ASA076 library.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    02 November 2010

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "ASA076_PRB:\n" );
  printf ( "  C version\n" );
  printf ( "  Test the ASA076 library.\n" );

  test01 ( );
  test02 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "ASA076_PRB:\n" );
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

    TEST01 demonstrates the use of TFN.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    02 November 2010

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
  printf ( "  TFN computes Owen's T function.\n" );
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

    t2 = tfn ( h, a );

    printf ( "  %12.4f  %12.4f  %24.16g  %24.16g  %10.4e\n",
      h, a, t1, t2, fabs ( t1 - t2 ) );
  }

  return;
}
/******************************************************************************/

void test02 ( )

/******************************************************************************/
/*
  Purpose:

    TEST02 demonstrates the use of THA.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    02 November 2010

  Author:

    John Burkardt
*/
{
  double a;
  double h;
  int n_data;
  double one = 1.0;
  double t1;
  double t2;

  printf ( "\n" );
  printf ( "TEST02:\n" );
  printf ( "  THA computes Owen's T function.\n" );
  printf ( "  Compare to tabulated values.\n" );
  printf ( "\n" );
  printf ( "             H             A      " );
  printf ( "    T                         T  \n" );
  printf ( "                                  " );
  printf ( "    (Tabulated)               (THA)               DIFF\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    owen_values ( &n_data, &h, &a, &t1 );

    if ( n_data == 0 )
    {
      break;
    }

    t2 = tha ( h, one, a, one );

    printf ( "  %12.4f  %12.4f  %24.16g  %24.16g  %10.4e\n",
      h, a, t1, t2, fabs ( t1 - t2 ) );
  }

  return;
}
