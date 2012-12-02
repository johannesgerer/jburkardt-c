# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "rkf45.h"

int main ( main );

void test07 ( main );
void r8_f_7 ( double t, double y[], double yp[] );
double r8_exact_7 ( double t );
/*
  Variables declared here can be accessed by any other functions that
  are declared in this file.

  One way to remind ourselves of the special status of these variables
  is to use CAPITAL LETTERS for their names.

  To be safe, we can also assign them initial values.
*/
  double ALPHA = 0.0;
  double BETA = 0.0;
  double GAMMA = 0.0;

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for RKF45_PRB2.

  Discussion:

    RKF45_PRB2 tests the RKF45 ODE integrator.

    TEST07, in particular, demonstrates one way to define parameters
    that are accessible to the derivative routine.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    20 November 2007

  Author:

    John Burkardt
*/
{
  timestamp ( );

  printf ( "\n" );
  printf ( "RKF45_PRB2\n" );
  printf ( "  C version\n" );
  printf ( "  Test the RKF45 library.\n" );
/*
  Pick some values for the parameters, call the test.
*/
  ALPHA = 1.0;
  BETA = 2.0;
  GAMMA = 3.0;

  test07 ( );

  ALPHA = 1.0;
  BETA = 0.5;
  GAMMA = 0.0;
/*
  Pick some values for the parameters, call the test again.
*/
  test07 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "RKF45_PRB2\n" );
  printf ( "  Normal end of execution.\n" );

  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void test07 ( )

/******************************************************************************/
/*
  Purpose:

    TEST07 solves a scalar ODE with parameters ALPHA, BETA, and GAMMA.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    20 November 2007

  Author:

    John Burkardt
*/
{
# define NEQN 1

  double abserr;
  int    flag;
  int    i_step;
  int    n_step;
  double relerr;
  double t;
  double t_out;
  double t_start;
  double t_stop;
  double y[NEQN];
  double yp[NEQN];

  printf ( "\n" );
  printf ( "TEST07\n" );
  printf ( "  Solve a scalar equation using R8_RKF:\n" );
  printf ( "\n" );
  printf ( "  There are three parameters, ALPHA, BETA, GAMMA\n" );
  printf ( "\n" );
  printf ( "  Y' = ALPHA * BETA * cos ( BETA * T + GAMMA )\n" );
  printf ( "\n" );
  printf ( "  ALPHA = %f\n", ALPHA );
  printf ( "  BETA  = %f\n", BETA );
  printf ( "  GAMMA = %f\n", GAMMA );

  abserr = sqrt ( r8_epsilon ( ) );
  relerr = sqrt ( r8_epsilon ( ) );

  flag = 1;

  t_start = 0.0;
  t_stop = 20.0;

  n_step = 5;

  t = 0.0;
  t_out = 0.0;
  y[0] = r8_exact_7 ( t );
  r8_f_7 ( t, y, yp );

  printf ( "\n" );
  printf ( "FLAG             T          Y         Y'          Y_Exact         Error\n" );
  printf ( "\n" );

  printf ( "%4d  %12f  %12f  %12f  %12f  %12f\n",
    flag, t, y[0], yp[0], r8_exact_7 ( t ), y[0] - r8_exact_7 ( t ) );

  for ( i_step = 1; i_step <= n_step; i_step++ )
  {
    t = ( ( double ) ( n_step - i_step + 1 ) * t_start  
        + ( double ) (          i_step - 1 ) * t_stop ) 
        / ( double ) ( n_step              );

    t_out = ( ( double ) ( n_step - i_step ) * t_start  
            + ( double ) (          i_step ) * t_stop )  
            / ( double ) ( n_step          );

    flag = r8_rkf45 ( r8_f_7, NEQN, y, yp, &t, t_out, &relerr, abserr, flag );

    printf ( "%4d  %12f  %12f  %12f  %12f  %12f\n",
      flag, t, y[0], yp[0], r8_exact_7 ( t ), y[0] - r8_exact_7 ( t ) );
  }

  return;
# undef NEQN
}
/******************************************************************************/

void r8_f_7 ( double t, double y[], double yp[] )

/******************************************************************************/
/*
  Purpose:

    R8_F_7 evaluates the derivative for the ODE for test #7.

  Discussion:

    Problem parameters ALPHA, BETA and GAMMA are in "global" memory.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    20 November 2007

  Author:

    John Burkardt

  Parameters:

    Input, double T, the value of the independent variable.

    Input, double Y[NEQN], the value of the dependent variable.

    Output, double YP[NEQN], the value of the derivative dY(1:NEQN)/dT.
*/
{
  yp[0] = ALPHA * BETA * cos ( BETA * t + GAMMA );

  return;
}
/******************************************************************************/

double r8_exact_7 ( double t )

/******************************************************************************/
/*
  Purpose:

    R8_EXACT_7 evaluates the exact solution of the ODE for test #7.

  Discussion:

    Problem parameters ALPHA, BETA and GAMMA are in "global" memory.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    20 November 2007

  Author:

    John Burkardt

  Parameters:

    Input, double T, the value of the independent variable.

    Output, double R8_EXACT_7, the exact solution.
*/
{
  double value;

  value = ALPHA * sin ( BETA * t + GAMMA );

  return value;
}
