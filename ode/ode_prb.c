# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "ode.h"

int main ( void );
void test01 ( void );
void test02 ( void );
void f01 ( double t, double y[], double yp[] );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for ODE_PRB.

  Discussion:

    ODE_PRB tests the ODE library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 February 2012

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "ODE_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the ODE library.\n" );

  test01 ( );
  test02 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "ODE_PRB\n" );
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

    TEST01 tests ODE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 February 2012

  Author:

    John Burkardt
*/
{
  double abserr;
  int i;
  int iflag;
  int iwork[5];
  int neqn = 2;
  double pi = 3.141592653589793;
  double relerr;
  int step_num = 12;
  double t;
  double tout;
  double *work;
  double *y;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  ODE solves a system of ordinary differential\n" );
  printf ( "  equations.\n" );
  printf ( "\n" );
  printf ( "      T           Y(1)         Y(2)\n" );
  printf ( "\n" );

  abserr = 0.00001;
  relerr = 0.00001;

  iflag = 1;

  t = 0.0;
  y = ( double * ) malloc ( neqn * sizeof ( double ) );
  y[0] = 1.0;
  y[1] = 0.0;

  printf ( "  %8g  %14g  %14g\n", t, y[0], y[1] );

  work = ( double * ) malloc ( ( 100 + 21 * neqn ) * sizeof ( double ) );

  for ( i = 1; i <= step_num; i++ )
  {
    tout = ( double ) ( i ) * 2.0 * pi / ( double ) ( step_num );

    ode ( f01, neqn, y, &t, tout, relerr, abserr, &iflag, work, iwork );

    if ( iflag != 2 )
    {
      printf ( "\n" );
      printf ( "TEST01 - Fatal error!\n" );
      printf ( "  ODE returned IFLAG = %d\n", iflag );
      break;
    }
    printf ( "  %8g  %14g  %14g\n", t, y[0], y[1] );
  }

  free ( work );
  free ( y );

  return;
}
/******************************************************************************/

void test02 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST02 tests ODE by integrating in the NEGATIVE time direction.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 February 2012

  Author:

    John Burkardt
*/
{
  double abserr;
  int i;
  int iflag;
  int iwork[5];
  int neqn = 2;
  double pi = 3.141592653589793;
  double relerr;
  int step_num = 12;
  double t;
  double tout;
  double *work;
  double *y;

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  ODE solves a system of ordinary differential\n" );
  printf ( "  equations.\n" );
  printf ( "\n" );
  printf ( "  In this example, we integrate in the negative\n" );
  printf ( "  time direction.\n" );
  printf ( "\n" );
  printf ( "      T           Y(1)         Y(2)\n" );
  printf ( "\n" );

  abserr = 0.00001;
  relerr = 0.00001;

  iflag = 1;

  t = 0.0;
  y = ( double * ) malloc ( neqn * sizeof ( double ) );
  y[0] = 1.0;
  y[1] = 0.0;

  printf ( "  %8g  %14g  %14g\n", t, y[0], y[1] );

  work = ( double * ) malloc ( ( 100 + 21 * neqn ) * sizeof ( double ) );

  for ( i = 1; i <= step_num; i++ )
  {
    tout = - ( double ) ( i ) * 2.0 * pi / ( double ) ( step_num );

    ode ( f01, neqn, y, &t, tout, relerr, abserr, &iflag, work, iwork );

    if ( iflag != 2 )
    {
      printf ( "\n" );
      printf ( "TEST02 - Fatal error!\n" );
      printf ( "  ODE returned IFLAG = %d\n", iflag );
      break;
    }
    printf ( "  %8g  %14g  %14g\n", t, y[0], y[1] );
  }

  free ( work );
  free ( y );

  return;
}
/******************************************************************************/

void f01 ( double t, double y[], double yp[] )

/******************************************************************************/
/*
  Purpose:

    F01 supplies the right hand side of the ODE for problem 1.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 February 2012

  Author:

    John Burkardt

  Parameters:

    Input, double T, the time.

    Input, double Y[], the dependent variable.

    Output, double YP[], the value of the derivative.
*/
{
  yp[0] =   y[1];
  yp[1] = - y[0];

  return;
}
