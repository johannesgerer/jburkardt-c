# include <stdlib.h>
# include <stdio.h>

# include "machar.h"

int main ( void );
void test01 ( void );
void test02 ( void );

/**********************************************************************/

int main ( void )

/**********************************************************************/
/*
  Purpose:

    MAIN is the main program for MACHAR_PRB.

  Discussion:

    MACHAR_PRB runs the MACHAR tests.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    30 November 2006

  Author:

    John Burkardt
*/
{
  timestamp ( );

  printf ( "\n" );
  printf ( "MACHAR_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the MACHAR library.\n" );

  test01 ( );
  test02 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "MACHAR_PRB\n" );
  printf ( "  Normal end of execution.\n" );

  printf ( "\n" );
  timestamp ( );

  return 0;
}
/**********************************************************************/

void test01 ( void )

/**********************************************************************/
/*
  Purpose:

    TEST01 tests R4_MACHAR.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    30 November 2006

  Author:

    John Burkardt
*/
{
  float eps;
  float epsneg;
  long int ibeta;
  long int iexp;
  long int irnd;
  long int it;
  long int machep;
  long int maxexp;
  long int minexp;
  long int negep;
  long int ngrd;
  float xmax;
  float xmin;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  Test R4_MACHAR, which computes single\n" );
  printf ( "  precision machine constants.\n" );

  r4_machar ( &ibeta, &it, &irnd, &ngrd, &machep, &negep, &iexp,
    &minexp, &maxexp, &eps, &epsneg, &xmin, &xmax );

  printf ( "\n" );
  printf ( "  IBETA =  %d\n", ibeta );
  printf ( "  IT =     %d\n", it );
  printf ( "  IRND =   %d\n", irnd );
  printf ( "  NGRD =   %d\n", ngrd );
  printf ( "  MACHEP = %d\n", machep );
  printf ( "  NEGEP =  %d\n", negep );
  printf ( "  IEXP =   %d\n", iexp );
  printf ( "  MINEXP = %d\n", minexp );
  printf ( "  MAXEXP = %d\n", maxexp );
  printf ( "  EPS =    %e\n", eps );
  printf ( "  EPSNEG = %e\n", epsneg );
  printf ( "  XMIN =   %e\n", xmin );
  printf ( "  XMAX =   %e\n", xmax );

  return;
}
/**********************************************************************/

void test02 ( void )

/**********************************************************************/
/*
  Purpose:

    TEST02 tests R8_MACHAR.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    30 November 2006

  Author:

    John Burkardt
*/
{
  double eps;
  double epsneg;
  long int ibeta;
  long int iexp;
  long int irnd;
  long int it;
  long int machep;
  long int maxexp;
  long int minexp;
  long int negep;
  long int ngrd;
  double xmax;
  double xmin;

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  Test R8_MACHAR, which computes double\n" );
  printf ( "  precision machine constants.\n" );

  r8_machar ( &ibeta, &it, &irnd, &ngrd, &machep, &negep, &iexp,
    &minexp, &maxexp, &eps, &epsneg, &xmin, &xmax );

  printf ( "\n" );
  printf ( "  IBETA =  %d\n", ibeta );
  printf ( "  IT =     %d\n", it );
  printf ( "  IRND =   %d\n", irnd );
  printf ( "  NGRD =   %d\n", ngrd );
  printf ( "  MACHEP = %d\n", machep );
  printf ( "  NEGEP =  %d\n", negep );
  printf ( "  IEXP =   %d\n", iexp );
  printf ( "  MINEXP = %d\n", minexp );
  printf ( "  MAXEXP = %d\n", maxexp );
  printf ( "  EPS =    %e\n", eps );
  printf ( "  EPSNEG = %e\n", epsneg );
  printf ( "  XMIN =   %e\n", xmin );
  printf ( "  XMAX =   %e\n", xmax );

  return;
}
