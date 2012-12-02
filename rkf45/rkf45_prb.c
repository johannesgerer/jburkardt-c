# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "rkf45.h"

int main ( void );

void test01 ( void );
void test02 ( void );
void test03 ( void );
void test04 ( void );
void test05 ( void );
void test06 ( void );
void r4_f1 ( float t, float y[], float yp[] );
float r4_y1x ( float t );
void r4_f2 ( float t, float y[], float yp[] );
void r8_f1 ( double t, double y[], double yp[] );
double r8_y1x ( double t );
void r8_f2 ( double t, double y[], double yp[] );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for RKF45_PRB.

  Discussion:

    RKF45_PRB tests the RKF45 ODE integrator.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 June 2006

  Author:

    John Burkardt
*/
{
  timestamp ( );

  printf ( "\n" );
  printf ( "RKF45_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the RKF45 library.\n" );

  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
  test05 ( );
  test06 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "RKF45_PRB\n" );
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

    TEST01 solves a scalar ODE.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 June 2006

  Author:

    John Burkardt
*/
{
# define NEQN 1

  float abserr;
  int flag;
  int i_step;
  int n_step;
  float relerr;
  float t;
  float t_out;
  float t_start;
  float t_stop;
  float y[NEQN];
  float yp[NEQN];

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  Solve a scalar equation using R4_RKF:\n" );
  printf ( "\n" );
  printf ( "  Y' = 0.25 * Y * ( 1 - Y / 20 )\n" );
  printf ( "\n" );

  abserr = sqrt ( r4_epsilon ( ) );
  relerr = sqrt ( r4_epsilon ( ) );

  flag = 1;

  t_start = 0.0;
  t_stop = 20.0;

  n_step = 5;

  t = 0.0;
  t_out = 0.0;
  y[0] = 1.0;
  r4_f1 ( t, y, yp );

  printf ( "\n" );
  printf ( "FLAG             T          Y         Y'          Y_Exact         Error\n" );
  printf ( "\n" );

  printf ( "%4d  %12f  %12f  %12f  %12f  %12f\n",
    flag, t, y[0], yp[0], r4_y1x ( t ), y[0] - r4_y1x ( t ) );

  for ( i_step = 1; i_step <= n_step; i_step++ )
  {
    t = ( ( float ) ( n_step - i_step + 1 ) * t_start  
        + ( float ) (          i_step - 1 ) * t_stop ) 
        / ( float ) ( n_step              );

    t_out = ( ( float ) ( n_step - i_step ) * t_start  
            + ( float ) (          i_step ) * t_stop )  
            / ( float ) ( n_step          );

    flag = r4_rkf45 ( r4_f1, NEQN, y, yp, &t, t_out, &relerr, abserr, flag );

    printf ( "%4d  %12f  %12f  %12f  %12f  %12f\n",
      flag, t, y[0], yp[0], r4_y1x ( t ), y[0] - r4_y1x ( t ) );
  }

  return;
# undef NEQN
}
/******************************************************************************/

void test02 ( )

/******************************************************************************/
/*
  Purpose:

    TEST02 solves a vector ODE.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 June 2006

  Author:

    John Burkardt
*/
{
# define NEQN 2

  float abserr;
  int flag;
  int i_step;
  int n_step;
  float relerr;
  float t;
  float t_out;
  float t_start;
  float t_stop;
  float y[NEQN];
  float yp[NEQN];

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  Solve a vector equation using R4_RKF:\n" );
  printf ( "\n" );
  printf ( "  Y'(1) =  Y(2)\n" );
  printf ( "  Y'(2) = -Y(1)\n" );
  printf ( "\n" );
  printf ( "\n" );
  printf ( "  This system is equivalent to the following\n" );
  printf ( "  second order system:\n" );
  printf ( "\n" );
  printf ( "  Z\" = - Z.\n" );

  abserr = sqrt ( r4_epsilon ( ) );
  relerr = sqrt ( r4_epsilon ( ) );

  flag = 1;

  t_start = 0.0;
  t_stop = 2.0 * 3.14159265;

  n_step = 12;

  t = 0.0;
  t_out = 0.0;

  y[0] = 1.0;
  y[1] = 0.0;

  printf ( "\n" );
  printf ( "FLAG             T          Y(1)       Y(2)\n" );
  printf ( "\n" );

  printf ( "%4d  %12f  %12f  %12f\n", flag, t, y[0], y[1] );

  for ( i_step = 1; i_step <= n_step; i_step++ )
  {
    t = ( ( float ) ( n_step - i_step + 1 ) * t_start 
        + ( float ) (          i_step - 1 ) * t_stop ) 
        / ( float ) ( n_step              );

    t_out = ( ( float ) ( n_step - i_step ) * t_start 
            + ( float ) (	   i_step ) * t_stop ) 
            / ( float ) ( n_step );

    flag = r4_rkf45 ( r4_f2, NEQN, y, yp, &t, t_out, &relerr, abserr, flag );

    printf ( "%4d  %12f  %12f  %12f\n", flag, t, y[0], y[1] );
  }

  return;
# undef NEQN
}
/******************************************************************************/

void test03 ( )

/******************************************************************************/
/*
  Purpose:

    TEST03 solves a scalar ODE using single step mode.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 June 2006

  Author:

    John Burkardt
*/
{
# define NEQN 1

  float abserr;
  int flag;
  int i_step;
  int n_step;
  float relerr;
  float t;
  float t_out;
  float t_start;
  float t_stop;
  float y[NEQN];
  float yp[NEQN];

  printf ( "\n" );
  printf ( "TEST03\n" );
  printf ( "  Solve a scalar equation using R4_RKF:\n" );
  printf ( "\n" );
  printf ( "  Y' = 0.25 * Y * ( 1 - Y / 20 )\n" );
  printf ( "\n" );
  printf ( "  This routine uses the SINGLE STEP mode.\n" );
  printf ( "\n" );

  abserr = sqrt ( r4_epsilon ( ) );
  relerr = sqrt ( r4_epsilon ( ) );

  flag = -1;

  t_start = 0.0;
  t_stop = 20.0;

  n_step = 5;

  t = 0.0;
  t_out = 0.0;
  y[0] = 1.0;
  r4_f1 ( t, y, yp );

  printf ( "\n" );
  printf ( "FLAG             T          Y         Y'        Y_Exact         Error\n" );
  printf ( "\n" );

  printf ( "%4d  %12f  %12f  %12f  %12f  %12f\n",
    flag, t, y[0], yp[0], r4_y1x ( t ), y[0] - r4_y1x ( t ) );

  for ( i_step = 1; i_step <= n_step; i_step++ )
  {
    t = ( ( float ) ( n_step - i_step + 1 ) * t_start  
        + ( float ) (          i_step - 1 ) * t_stop ) 
        / ( float ) ( n_step              );

    t_out = ( ( float ) ( n_step - i_step ) * t_start  
            + ( float ) (          i_step ) * t_stop )  
            / ( float ) ( n_step          );

    while ( flag < 0 )
    {
      flag = r4_rkf45 ( r4_f1, NEQN, y, yp, &t, t_out, &relerr, abserr, flag );

      printf ( "%4d  %12f  %12f  %12f  %12f  %12f\n",
        flag, t, y[0], yp[0], r4_y1x ( t ), y[0] - r4_y1x ( t ) );
    }
    flag = -2;
  }

  return;
# undef NEQN
}
/******************************************************************************/

void test04 ( )

/******************************************************************************/
/*
  Purpose:

    TEST04 solves a scalar ODE.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 June 2006

  Author:

    John Burkardt
*/
{
# define NEQN 1

  double abserr;
  int flag;
  int i_step;
  int n_step;
  double relerr;
  double t;
  double t_out;
  double t_start;
  double t_stop;
  double y[NEQN];
  double yp[NEQN];

  printf ( "\n" );
  printf ( "TEST04\n" );
  printf ( "  Solve a scalar equation using R8_RKF:\n" );
  printf ( "\n" );
  printf ( "  Y' = 0.25 * Y * ( 1 - Y / 20 )\n" );
  printf ( "\n" );

  abserr = sqrt ( r8_epsilon ( ) );
  relerr = sqrt ( r8_epsilon ( ) );

  flag = 1;

  t_start = 0.0;
  t_stop = 20.0;

  n_step = 5;

  t = 0.0;
  t_out = 0.0;
  y[0] = 1.0;
  r8_f1 ( t, y, yp );

  printf ( "\n" );
  printf ( "FLAG             T          Y         Y'          Y_Exact         Error\n" );
  printf ( "\n" );

  printf ( "%4d  %12f  %12f  %12f  %12f  %12f\n",
    flag, t, y[0], yp[0], r8_y1x ( t ), y[0] - r8_y1x ( t ) );

  for ( i_step = 1; i_step <= n_step; i_step++ )
  {
    t = ( ( double ) ( n_step - i_step + 1 ) * t_start  
        + ( double ) (          i_step - 1 ) * t_stop ) 
        / ( double ) ( n_step              );

    t_out = ( ( double ) ( n_step - i_step ) * t_start  
            + ( double ) (          i_step ) * t_stop )  
            / ( double ) ( n_step          );

    flag = r8_rkf45 ( r8_f1, NEQN, y, yp, &t, t_out, &relerr, abserr, flag );

    printf ( "%4d  %12f  %12f  %12f  %12f  %12f\n",
      flag, t, y[0], yp[0], r4_y1x ( t ), y[0] - r4_y1x ( t ) );
  }

  return;
# undef NEQN
}
/******************************************************************************/

void test05 ( )

/******************************************************************************/
/*
  Purpose:

    TEST05 solves a vector ODE.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 June 2006

  Author:

    John Burkardt
*/
{
# define NEQN 2

  double abserr;
  int flag;
  int i_step;
  int n_step;
  double relerr;
  double t;
  double t_out;
  double t_start;
  double t_stop;
  double y[NEQN];
  double yp[NEQN];

  printf ( "\n" );
  printf ( "TEST05\n" );
  printf ( "  Solve a vector equation using R8_RKF:\n" );
  printf ( "\n" );
  printf ( "  Y'(1) =  Y(2)\n" );
  printf ( "  Y'(2) = -Y(1)\n" );
  printf ( "\n" );

  abserr = sqrt ( r8_epsilon ( ) );
  relerr = sqrt ( r8_epsilon ( ) );

  flag = 1;

  t_start = 0.0;
  t_stop = 2.0 * 3.14159265;

  n_step = 12;

  t = 0.0;
  t_out = 0.0;

  y[0] = 1.0;
  y[1] = 0.0;
  r8_f2 ( t, y,  yp );

  printf ( "\n" );
  printf ( "FLAG             T          Y(1)       Y(2)\n" );
  printf ( "\n" );

  printf ( "%4d  %12f  %12f  %12f\n", flag, t, y[0], y[1] );

  for ( i_step = 1; i_step <= n_step; i_step++ )
  {
    t = ( ( double ) ( n_step - i_step + 1 ) * t_start 
        + ( double ) (          i_step - 1 ) * t_stop ) 
        / ( double ) ( n_step              );

    t_out = ( ( double ) ( n_step - i_step ) * t_start 
            + ( double ) (	   i_step ) * t_stop ) 
            / ( double ) ( n_step );

    flag = r8_rkf45 ( r8_f2, NEQN, y, yp, &t, t_out, &relerr, abserr, flag );

    printf ( "%4d  %12f  %12f  %12f\n", flag, t, y[0], y[1] );
  }

  return;
# undef NEQN
}
/******************************************************************************/

void test06 ( )

/******************************************************************************/
/*
  Purpose:

    TEST06 solves a scalar ODE using single step mode.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 June 2006

  Author:

    John Burkardt
*/
{
# define NEQN 1

  double abserr;
  int flag;
  int i_step;
  int n_step;
  double relerr;
  double t;
  double t_out;
  double t_start;
  double t_stop;
  double y[NEQN];
  double yp[NEQN];

  printf ( "\n" );
  printf ( "TEST06\n" );
  printf ( "  Solve a scalar equation using R8_RKF:\n" );
  printf ( "\n" );
  printf ( "  Y' = 0.25 * Y * ( 1 - Y / 20 )\n" );
  printf ( "\n" );
  printf ( "  This routine uses the SINGLE STEP mode.\n" );
  printf ( "\n" );

  abserr = sqrt ( r8_epsilon ( ) );
  relerr = sqrt ( r8_epsilon ( ) );

  flag = -1;

  t_start = 0.0;
  t_stop = 20.0;

  n_step = 5;

  t = 0.0;
  t_out = 0.0;
  y[0] = 1.0;
  r8_f1 ( t, y, yp );

  printf ( "\n" );
  printf ( "FLAG             T          Y         Y'        Y_Exact         Error\n" );
  printf ( "\n" );

  printf ( "%4d  %12f  %12f  %12f  %12f  %12f\n",
    flag, t, y[0], yp[0], r4_y1x ( t ), y[0] - r4_y1x ( t ) );

  for ( i_step = 1; i_step <= n_step; i_step++ )
  {
    t = ( ( double ) ( n_step - i_step + 1 ) * t_start  
        + ( double ) (          i_step - 1 ) * t_stop ) 
        / ( double ) ( n_step              );

    t_out = ( ( double ) ( n_step - i_step ) * t_start  
            + ( double ) (          i_step ) * t_stop )  
            / ( double ) ( n_step          );

    while ( flag < 0 )
    {
      flag = r8_rkf45 ( r8_f1, NEQN, y, yp, &t, t_out, &relerr, abserr, flag );

      printf ( "%4d  %12f  %12f  %12f  %12f  %12f\n",
        flag, t, y[0], yp[0], r4_y1x ( t ), y[0] - r4_y1x ( t ) );
    }
    flag = -2;
  }

  return;
# undef NEQN
}
/******************************************************************************/

void r4_f1 ( float t, float y[], float yp[] )

/******************************************************************************/
/*
  Purpose:

    R4_F1 evaluates the derivative for the ODE.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    26 March 2004

  Author:

    John Burkardt

  Parameters:

    Input, float T, the value of the independent variable.

    Input, float Y[NEQN], the value of the dependent variable.

    Output, float YP[NEQN], the value of the derivative dY(1:NEQN)/dT.
*/
{
  yp[0] = 0.25 * y[0] * ( 1.0 - y[0] / 20.0 );

  return;
}
/******************************************************************************/

float r4_y1x ( float t )

/******************************************************************************/
/*
  Purpose:

    R4_Y1X evaluates the exact solution of the ODE.

  Modified:

    26 March 2004

  Author:

    John Burkardt

  Parameters:

    Input, float T, the value of the independent variable.

    Output, float R4_Y1X, the exact solution.
*/
{
  float value;

  value = 20.0 / ( 1.0 + 19.0 * exp ( -0.25 * t ) );

  return value;
}
/******************************************************************************/

void r4_f2 ( float t, float y[], float yp[] )

/******************************************************************************/
/*
  Purpose:

    R4_F2 evaluates the derivative for the ODE.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    26 March 2004

  Author:

    John Burkardt

  Parameters:

    Input, float T, the value of the independent variable.

    Input, float Y(NEQN), the value of the dependent variable.

    Output float YP(NEQN), the value of the derivative dY(1:NEQN)/dT.
*/
{
  yp[0] =  y[1];
  yp[1] = -y[0];

  return;
}
/******************************************************************************/

void r8_f1 ( double t, double y[], double yp[] )

/******************************************************************************/
/*
  Purpose:

    R8_F1 evaluates the derivative for the ODE.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    26 March 2004

  Author:

    John Burkardt

  Parameters:

    Input, double T, the value of the independent variable.

    Input, double Y[NEQN], the value of the dependent variable.

    Output, double YP[NEQN], the value of the derivative dY(1:NEQN)/dT.
*/
{
  yp[0] = 0.25 * y[0] * ( 1.0 - y[0] / 20.0 );

  return;
}
/******************************************************************************/

double r8_y1x ( double t )

/******************************************************************************/
/*
  Purpose:

    R8_Y1X evaluates the exact solution of the ODE.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    26 March 2004

  Author:

    John Burkardt

  Parameters:

    Input, double T, the value of the independent variable.

    Output, double R8_Y1X, the exact solution.
*/
{
  double value;

  value = 20.0 / ( 1.0 + 19.0 * exp ( -0.25 * t ) );

  return value;
}
/******************************************************************************/

void r8_f2 ( double t, double y[], double yp[] )

/******************************************************************************/
/*
  Purpose:

    R8_F2 evaluates the derivative for the ODE.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    26 March 2004

  Author:

    John Burkardt

  Parameters:

    Input, double T, the value of the independent variable.

    Input, double Y(NEQN), the value of the dependent variable.

    Output double YP(NEQN), the value of the derivative dY(1:NEQN)/dT.
*/
{
  yp[0] =  y[1];
  yp[1] = -y[0];

  return;
}
