# include <stdlib.h>
# include <stdio.h>
# include <time.h>
# include <math.h>

# include "fd1d_wave.h"

int main ( void );
void fd1d_wave_test01 ( void );
double u_x1_01 ( double t );
double u_x2_01 ( double t );
double *u_t1_01 ( int x_num, double x_vec[] );
double *ut_t1_01 ( int x_num, double x_vec[] );
void fd1d_wave_test02 ( void );
double u_x1_02 ( double t );
double u_x2_02 ( double t );
double *u_t1_02 ( int x_num, double x_vec[] );
double *ut_t1_02 ( int x_num, double x_vec[] );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    FD1D_WAVE_PRB tests the FD1D finite difference wave computation.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    24 January 2012

  Author:

    John Burkardt
*/
{
  timestamp ( );

  printf ( "\n" );
  printf ( "FD1D_WAVE_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the FD1D_WAVE library.\n" );

  fd1d_wave_test01 ( );
  fd1d_wave_test02 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "FD1D_WAVE_PRB\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void fd1d_wave_test01 ( void )

/******************************************************************************/
/*
  Purpose:

    FD1D_WAVE_TEST_01 tests the FD1D finite difference wave computation.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    24 January 2012

  Author:

    John Burkardt
*/
{
  double alpha;
  double c;
  int i;
  int j;
  double t;
  double t_delta;
  int t_num = 41;
  double *t_vec;
  double t1;
  double t2;
  double *u;
  double *u1;
  double *u2;
  double *u3;
  int x_num = 16;
  double *x_vec;
  double x1;
  double x2;

  printf ( "\n" );
  printf ( "FD1D_WAVE_TEST01\n" );
  printf ( "  Try the \"shark\" wave.\n" );

  x1 = 0.0;
  x2 = 1.5;
  x_vec = r8vec_linspace_new ( x_num, x1, x2 );

  t1 = 0.0;
  t2 = 4.0;
  t_vec = r8vec_linspace_new ( t_num, t1, t2 );
  t_delta = ( t2 - t1 ) / ( double ) ( t_num - 1 );

  c = 1.0;
  alpha = fd1d_wave_alpha ( x_num, x1, x2, t_num, t1, t2, c );
/*
  Load the initial condition.
*/
  u = ( double * ) malloc ( t_num * x_num * sizeof ( double ) );

  u1 = u_t1_01 ( x_num, x_vec );

  for ( j = 0; j < x_num; j++ )
  {
    u[0+j*t_num] = u1[j];
  }
/*
  Take the first step.
*/
  t = t_vec[1];
  u2 = fd1d_wave_start ( x_num, x_vec, t, t_delta, alpha, u_x1_01, u_x2_01, 
    ut_t1_01, u1 );

  for ( j = 0; j < x_num; j++ )
  {
    u[1+j*t_num] = u2[j];
  }
/*
  Take all the other steps.
*/
  for ( i = 2; i < t_num; i++ )
  {
    t = t_vec[i];
    u3 = fd1d_wave_step ( x_num, t, alpha, u_x1_01, u_x2_01, u1, u2 );
    for ( j = 0; j < x_num; j++ )
    {
      u[i+j*t_num] = u3[j];
      u1[j] = u2[j];
      u2[j] = u3[j];
    }
    free ( u3 );
  }
/*
  Write the solution to a file.
*/
  r8mat_write ( "test01_plot.txt", t_num, x_num, u );

  printf ( "\n" );
  printf ( "  Plot data written to \"test01_plot.txt\".\n" );

  free ( t_vec );
  free ( u );
  free ( u1 );
  free ( u2 );
  free ( x_vec );

  return;
}
/******************************************************************************/

double u_x1_01 ( double t )

/******************************************************************************/
/*
  Purpose:

    U_X1_01 evaluates U at the boundary X1.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    24 January 2012

  Author:

    John Burkardt

  Parameters:

    Input, double T, the time.

    Output, double U, the value of U(T,X1).
*/
{
  int nd = 6;
  double td[6] = { 0.0, 0.10, 0.20, 0.30, 0.40, 0.50 };
  double tv[1];
  double u;
  double ud[6] = { 0.0, 2.0, 10.0, 8.0, 5.0, 0.0 };
  double *uv;

  tv[0] = t;

  uv = piecewise_linear ( nd, td, ud, 1, tv );

  u = uv[0];

  free ( uv );

  return u;
}
/******************************************************************************/

double u_x2_01 ( double t )

/******************************************************************************/
/*
  Purpose:

    U_X2_01 evaluates U at the boundary X2.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    24 January 2012

  Author:

    John Burkardt

  Parameters:

    Input, double T, the time.

    Output, double U, the value of U(T,X2).
*/
{
  double u;

  u = 0.0;

  return u;
}
/******************************************************************************/

double *u_t1_01 ( int x_num, double x_vec[] )

/******************************************************************************/
/*
  Purpose:

    U_T1_01 evaluates U at the initial time T1.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    24 January 2012

  Author:

    John Burkardt

  Parameters:

    Input, int X_NUM, the number of nodes.

    Input, double X_VEC[X_NUM], the coordinates of the nodes.

    Output, double U_T1_01[X_NUM], the value of U at the initial time. 
*/
{
  int j;
  double *u;

  u = ( double * ) malloc ( x_num * sizeof ( double ) );

  for ( j = 0; j < x_num; j++ )
  {
    u[j] = 0.0;
  }

  return u;
}
/******************************************************************************/

double *ut_t1_01 ( int x_num, double x_vec[] )

/******************************************************************************/
/*
  Purpose:

    UT_T1_01 evaluates dUdT at the initial time T1.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    24 January 2012

  Author:

    John Burkardt

  Parameters:

    Input, int X_NUM, the number of nodes.

    Input, double X_VEC[X_NUM], the coordinates of the nodes.

    Output, double UT_T1_01[X_NUM], the value of dUdT at the initial time. 
*/
{
  int j;
  double *ut;

  ut = ( double * ) malloc ( x_num * sizeof ( double ) );

  for ( j = 0; j < x_num; j++ )
  {
    ut[j] = 0.0;
  }

  return ut;
}
/******************************************************************************/

void fd1d_wave_test02 ( void )

/******************************************************************************/
/*
  Purpose:

    FD1D_WAVE_TEST_02 tests the FD1D finite difference wave computation.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    24 January 2012

  Author:

    John Burkardt
*/
{
  double alpha;
  double c;
  int i;
  int j;
  double t;
  double t_delta;
  double *t_vec;
  int t_num = 41;
  double t1;
  double t2;
  double *u;
  double *u1;
  double *u2;
  double *u3;
  int x_num = 16;
  double *x_vec;
  double x1;
  double x2;

  printf ( "\n" );
  printf ( "FD1D_WAVE_TEST02\n" );
  printf ( "  Try a sine curve.\n" );

  x1 = 0.0;
  x2 = 1.5;
  x_vec = r8vec_linspace_new ( x_num, x1, x2 );

  t1 = 0.0;
  t2 = 4.0;
  t_vec = r8vec_linspace_new ( t_num, t1, t2 );
  t_delta = ( t2 - t1 ) / ( double ) ( t_num - 1 );
/*
  Changing T2 to 4.5 is enough to push the algorithm into instability.

  t2 = 4.5;
*/
  c = 1.0;
  alpha = fd1d_wave_alpha ( x_num, x1, x2, t_num, t1, t2, c );
/*
  Load the initial condition.
*/
  u = ( double * ) malloc ( t_num * x_num * sizeof ( double ) );

  u1 = u_t1_02 ( x_num, x_vec );

  for ( j = 0; j < x_num; j++ )
  {
    u[0+j*t_num] = u1[j];
  }
/*
  Take the first step.
*/
  t = t_vec[1];
  u2 = fd1d_wave_start ( x_num, x_vec, t, t_delta, alpha, u_x1_02, u_x2_02, 
    ut_t1_02, u1 );

  for ( j = 0; j < x_num; j++ )
  {
    u[1+j*t_num] = u2[j];
  }
/*
  Take all the other steps.
*/
  for ( i = 2; i < t_num; i++ )
  {
    t = t_vec[i];
    u3 = fd1d_wave_step ( x_num, t, alpha, u_x1_02, u_x2_02, u1, u2 );
    for ( j = 0; j < x_num; j++ )
    {
      u[i+j*t_num] = u3[j];
      u1[j] = u2[j];
      u2[j] = u3[j];
    }
    free ( u3 );
  }
/*
  Write the solution to a file.
*/
  r8mat_write ( "test02_plot.txt", t_num, x_num, u );

  printf ( "\n" );
  printf ( "  Plot data written to \"test02_plot.txt\".\n" );

  free ( t_vec );
  free ( u );
  free ( u1 );
  free ( u2 );
  free ( x_vec );

  return;
}
/******************************************************************************/

double u_x1_02 ( double t )

/******************************************************************************/
/*
  Purpose:

    U_X1_02 evaluates U at the boundary X1.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    24 January 2012

  Author:

    John Burkardt

  Parameters:

    Input, double T, the time.

    Output, double U_X1_02, the value of U(T,X1).
*/
{
  double u;

  u = 0.0;

  return u;
}
/******************************************************************************/

double u_x2_02 ( double t )

/******************************************************************************/
/*
  Purpose:

    U_X2_02 evaluates U at the boundary X2.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    24 January 2012

  Author:

    John Burkardt

  Parameters:

    Input, double T, the time.

    Output, double U, the value of U(T,X2).
*/
{
  double u;

  u = 0.0;

  return u;
}
/******************************************************************************/

double *u_t1_02 ( int x_num, double x_vec[] )

/******************************************************************************/
/*
  Purpose:

    U_T1_02 evaluates U at the initial time T1.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    24 January 2012

  Author:

    John Burkardt

  Parameters:

    Input, int X_NUM, the number of nodes.

    Input, double X_VEC[X_NUM], the spatial node coordinates.

    Output, double U[X_NUM], the value of U at the initial time,
    and every node.
*/
{
  int j;
  double pi = 3.141592653589793;
  double *u;

  u = ( double * ) malloc ( x_num * sizeof ( double ) );

  for ( j = 0; j < x_num; j++ )
  {
    u[j] = sin ( 2.0 * pi * x_vec[j] );
  }

  return u;
}
/******************************************************************************/

double *ut_t1_02 ( int x_num, double x_vec[] )

/******************************************************************************/
/*
  Purpose:

    UT_T1_02 evaluates dUdT at the initial time T1.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    24 January 2012

  Author:

    John Burkardt

  Parameters:

    Input, int X_NUM, the number of spatial intervals.

    Input, double X_VEC[X_NUM], the spatial node coordinates.

    Output, double UT_T1_02[X_NUM], the value of dUdT at the initial time,
    and every node.
*/
{
  int j;
  double *ut;

  ut = ( double * ) malloc ( x_num * sizeof ( double ) );

  for ( j = 0; j < x_num; j++ )
  {
    ut[j] = 0.0;
  }

  return ut;
}
