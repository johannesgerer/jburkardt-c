# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <string.h>

# include "stochastic_heat2d.h"

int main ( );
void test01 ( );
void test02 ( );
void boundary ( int nx, int ny, double x[], double y[], int n, double a[], 
  double rhs[] );
double test01_f ( double x, double y );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for STOCHASTIC_HEAT2D_PRB.

  Discussion:

    STOCHASTIC_HEAT2D_PRB tests the STOCHASTIC_HEAT2D library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    03 September 2013

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "STOCHASTIC_HEAT2D_PRB:\n" );
  printf ( "  C version\n" );
  printf ( "  Test the STOCHASTIC_HEAT2D library.\n" );

  test01 ( );
  test02 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "STOCHASTIC_HEAT2D_PRB:\n" );
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

    TEST01 plots a sample solution of a 2D stochastic diffusivity equation.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    03 September 2013

  Author:

    John Burkardt
*/
{
  char command_filename[] = "solution_commands.txt";
  FILE *command_unit;
  char data_filename[] = "solution_data.txt";
  FILE *data_unit;
  int i;
  int j;
  int nx;
  int ny;
  double *omega;
  int seed;
  double *umat;
  double u_mean;
  double *xmat;
  double xmax;
  double xmin;
  double *xvec;
  double *ymat;
  double ymax;
  double ymin;
  double *yvec;

  printf ( "\n" );
  printf ( "TEST01:\n" );
  printf ( "  Consider the steady heat equation in the unit square,\n" );
  printf ( "  with 0 Dirichlet boundary conditions, \n" );
  printf ( "  and a heat source term F that is a Gaussian centered at (0.60,0.80).\n" );
  printf ( "\n" );
  printf ( "  Model the diffusivity coefficient as spatially varying,\n" );
  printf ( "  with a stochastic dependence on parameters OMEGA(1:4),\n" );
  printf ( "  as described in Babuska, Nobile, Tempone (BNT).\n" );
  printf ( "\n" );
  printf ( "  Compute and display the solution U for a given choice\n" );
  printf ( "  of the parameters OMEGA.\n" );
/*
  Create the X and Y coordinate vectors.
*/
  nx = 21;
  xmin = 0.0;
  xmax = 1.0;
  xvec = r8vec_linspace_new ( nx, xmin, xmax );

  ny = 21;
  ymin = 0.0;
  ymax = 1.0;
  yvec = r8vec_linspace_new ( ny, ymin, ymax );
/*
  Create the X and Y coordinate matrices.
*/
  xmat = ( double * ) malloc ( nx * ny * sizeof ( double ) );
  ymat = ( double * ) malloc ( nx * ny * sizeof ( double ) );
  r8vec_mesh_2d ( nx, ny, xvec, yvec, xmat, ymat );
/*
  Sample OMEGA:
*/
  seed = 123456789;
  omega = r8vec_normal_01_new ( 4, &seed );
  for ( i = 0; i < 4; i++ )
  {
    omega[i] = 2.0 * omega[i];
  }

  r8vec_print ( 4, omega, "  Sampled OMEGA values:" );
/*
  Solve the finite difference approximation to the steady 2D heat equation
  for this set of OMEGA values.
*/
  umat = stochastic_heat2d ( omega, nx, ny, xvec, yvec, test01_f );
/*
  Create a data file.
*/
  data_unit = fopen ( data_filename, "wt" );
  for ( j = 0; j < ny; j++ )
  {
    for ( i = 0; i < nx; i++ )
    {
      fprintf ( data_unit, "%g  %g  %g\n", 
        xmat[i+j*nx], ymat[i+j*nx], umat[i+j*nx] );
    }
    fprintf ( data_unit, "\n" );
  }
  fclose ( data_unit );
  printf ( "\n" );
  printf ( "  Created graphics data file '%s'\n", data_filename );
/*
  Create the command file.
*/
  command_unit = fopen ( command_filename, "wt" );

  fprintf ( command_unit, "# %s\n", command_filename );
  fprintf ( command_unit, "#\n" );
  fprintf ( command_unit, "# Usage:\n" );
  fprintf ( command_unit, "#  gnuplot < %s\n", command_filename );
  fprintf ( command_unit, "#\n" );
  fprintf ( command_unit, "set term png\n" );
  fprintf ( command_unit, "set output 'solution.png'\n" );
  fprintf ( command_unit, "set xlabel '<---X--->'\n" );
  fprintf ( command_unit, "set ylabel '<---Y--->'\n" );
  fprintf ( command_unit, "set zlabel '<---U(X,Y)--->'\n" );
  fprintf ( command_unit, "set title 'Sample Solution'\n" );
  fprintf ( command_unit, "set contour\n" );
  fprintf ( command_unit, "set timestamp\n" );
  fprintf ( command_unit, "set cntrparam levels 10\n" );
  fprintf ( command_unit, "set view 75, 75\n" );
  fprintf ( command_unit, "unset key\n" );
  fprintf ( command_unit, "splot '%s'\n", data_filename );

  fclose ( command_unit );

  printf ( "  Created graphics command file '%s'\n", command_filename );
/*
  Report the average value of U.
*/
  u_mean = r8mat_mean ( nx, ny, umat );

  printf ( "\n" );
  printf ( "  Mean value of U is %g\n", u_mean );
/*
  Free memory.
*/
  free ( omega );
  free ( umat );
  free ( xmat );
  free ( xvec );
  free ( ymat );
  free ( yvec );

  return;
}
/******************************************************************************/

void test02 ( )

/******************************************************************************/
/*
  Purpose:

    TEST02 looks at mean temperature as a function of OMEGA(1) and OMEGA(2).

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    03 September 2013

  Author:

    John Burkardt
*/
{
  char command_filename[] = "umean_commands.txt";
  FILE *command_unit;
  char data_filename[] = "umean_data.txt";
  FILE *data_unit;
  int i;
  int j;
  int nx;
  int ny;
  double omega[4];
  double *omega1_mat;
  double omega1_max;
  double omega1_min;
  int omega1_num;
  double *omega1_vec;
  double *omega2_mat;
  double omega2_max;
  double omega2_min;
  int omega2_num;
  double *omega2_vec;
  double *umat;
  double *u_mean_mat;
  double u_mean_max;
  double *xmat;
  double xmax;
  double xmin;
  double *xvec;
  double *ymat;
  double ymax;
  double ymin;
  double *yvec;

  printf ( "\n" );
  printf ( "TEST02:\n" );
  printf ( "  Fix OMEGA(3) = 4, OMEGA(4) = 0, and\n" );
  printf ( "  examine dependence of average temperature on OMEGA(1) and OMEGA(2)\n" );
  printf ( "  over the range [-10,+10].\n" );
/*
  Create the X and Y coordinate vectors.
*/
  nx = 21;
  xmin = 0.0;
  xmax = 1.0;
  xvec = r8vec_linspace_new ( nx, xmin, xmax );

  ny = 21;
  ymin = 0.0;
  ymax = 1.0;
  yvec = r8vec_linspace_new ( ny, ymin, ymax );
/*
  Create the X and Y coordinate matrices.
*/
  xmat = ( double * ) malloc ( nx * ny * sizeof ( double ) );
  ymat = ( double * ) malloc ( nx * ny * sizeof ( double ) );
  r8vec_mesh_2d ( nx, ny, xvec, yvec, xmat, ymat );
/*
  Create OMEGA1 and OMEGA2 vectors.
*/
  omega1_num = 21;
  omega1_min = -10.0;
  omega1_max = +10.0;
  omega1_vec = r8vec_linspace_new ( omega1_num, omega1_min, omega1_max );

  omega2_num = 21;
  omega2_min = -10.0;
  omega2_max = +10.0;
  omega2_vec = r8vec_linspace_new ( omega2_num, omega2_min, omega2_max );
/*
  Create the OMEGA1 and OMEGA2 coordinate matrices.
*/
  omega1_mat = ( double * ) malloc ( omega1_num * omega2_num * sizeof ( double ) );
  omega2_mat = ( double * ) malloc ( omega1_num * omega2_num * sizeof ( double ) );
  r8vec_mesh_2d ( omega1_num, omega2_num, omega1_vec, omega2_vec, omega1_mat, omega2_mat );
/*
  Set OMEGA(3) and OMEGA(4).
*/
  omega[2] = 4.0;
  omega[3] = 0.0;

  printf ( "\n" );
  printf ( "  Omega(3) fixed at %g\n", omega[2] );
  printf ( "  Omega(4) fixed at %g\n", omega[3] );
/*
  Solve the finite difference approximation to the steady 2D heat equation,
  and save the mean value of the solution, which is a slightly biased
  estimate of the heat integral over the unit square.
*/
  u_mean_mat = ( double * ) malloc ( omega1_num * omega2_num * sizeof ( double ) );

  for ( j = 0; j < omega2_num; j++ )
  {
    omega[1] = omega2_vec[j];
    for ( i = 0; i < omega1_num; i++ )
    {
      omega[0] = omega1_vec[i];
      umat = stochastic_heat2d ( omega, nx, ny, xvec, yvec, test01_f );
      u_mean_mat[i+j*omega1_num] = r8mat_mean ( nx, ny, umat );
      free ( umat );
    }
  }
/*
  Create a data file.
*/
  data_unit = fopen ( data_filename, "wt" );
  for ( j = 0; j < ny; j++ )
  {
    for ( i = 0; i < nx; i++ )
    {
      fprintf ( data_unit, "%g  %g  %g\n",
        omega1_mat[i+j*omega1_num], omega2_mat[i+j*omega1_num], u_mean_mat[i+j*omega1_num] );
    }
    fprintf ( data_unit, "\n" );
  }
  fclose ( data_unit );
  printf ( "\n" );
  printf ( "  Created graphics data file '%s'\n", data_filename );
/*
  Create the command file.
*/
  command_unit = fopen ( command_filename, "wt" );

  fprintf ( command_unit, "# %s\n", command_filename );
  fprintf ( command_unit, "#\n" );
  fprintf ( command_unit, "# Usage:\n" );
  fprintf ( command_unit, "#  gnuplot < %s\n", command_filename );
  fprintf ( command_unit, "#\n" );
  fprintf ( command_unit, "set term png\n" );
  fprintf ( command_unit, "set output 'umean.png'\n" );
  fprintf ( command_unit, "set xlabel '<---OMEGA1--->'\n" );
  fprintf ( command_unit, "set ylabel '<---OMEGA2--->'\n" );
  fprintf ( command_unit, "set zlabel '<---U_MEAN(OMEGA1,OMEGA2)--->'\n" );
  fprintf ( command_unit, "set title 'Solution Mean as Function of Omega1, Omega2'\n" );
  fprintf ( command_unit, "set contour\n" );
  fprintf ( command_unit, "set timestamp\n" );
  fprintf ( command_unit, "set cntrparam levels 10\n" );
  fprintf ( command_unit, "set view 75, 75\n" );
  fprintf ( command_unit, "unset key\n" );
  fprintf ( command_unit, "splot '%s'\n", data_filename );

  fclose ( command_unit );

  printf ( "  Created graphics command file '%s'\n", command_filename );
/*
  Print the maximum value of the mean.
*/
  u_mean_max = r8mat_max ( omega1_num, omega2_num, u_mean_mat );

  printf ( "\n" );
  printf ( "  U_Mean_Max = %g\n", u_mean_max );
/*
  Free memory.
*/
  free ( omega1_mat );
  free ( omega1_vec );
  free ( omega2_mat );
  free ( omega2_vec );
  free ( u_mean_mat );
  free ( xmat );
  free ( xvec );
  free ( ymat );
  free ( yvec );

  return;
}
/******************************************************************************/

void boundary ( int nx, int ny, double x[], double y[], int n, double a[], 
  double rhs[] )

/******************************************************************************/
/*
  Purpose:

    BOUNDARY sets up the matrix and right hand side at boundary nodes.

  Discussion:

    For this simple problem, the boundary conditions specify that the solution
    is 100 on the left side, and insulated on the right, top and bottom.

    Nodes are assigned a single index K, which increases as:

    (NY-1)*NX+1  (NY-1)*NX+2  ...  NY * NX
           ....         ....  ...    .....
           NX+1         NX+2  ...   2 * NX
              1            2  ...       NX

    The index K of a node on the lower boundary satisfies:
      1 <= K <= NX
    The index K of a node on the upper boundary satisfies:
      (NY-1)*NX+1 <= K <= NY * NX
    The index K of a node on the left boundary satisfies:
      mod ( K, NX ) = 1
    The index K of a node on the right boundary satisfies:
      mod ( K, NX ) = 0

    If we number rows from bottom I = 1 to top I = NY
    and columns from left J = 1 to right J = NX, then the relationship
    between the single index K and the row and column indices I and J is:
      K = ( I - 1 ) * NX + J
    and
      J = 1 + mod ( K - 1, NX )
      I = 1 + ( K - J ) / NX
      
  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    03 September 2013

  Author:

    John Burkardt

  Parameters:

    Input, int NX, NY, the number of grid points in X and Y.

    Input, double X[NX], Y[NY], the coordinates of grid lines.

    Input, int N, the number of nodes.

    Input/output, double A[N*N].  On input, the system matrix, with the 
    entries for the interior nodes filled in.  On output, the entries for
    the boundary nodes have been set as well.

    Input, double RHS[N], on input, the system right hand side, 
    with the entries for the interior nodes filled in.  On output, the entries for
    the boundary nodes have been set as well.
*/
{
  int i;
  int j;
  int kc;
/*
  Left boundary.
*/
  j = 0;
  for ( i = 1; i < ny - 1; i++ )
  {
    kc = i * nx + j;
    a[kc+kc*n] = a[kc+kc*n] + 1.0;
    rhs[kc] = 0.0;
  }
/*
  Right boundary.
*/
  j = nx - 1;
  for ( i = 1; i < ny - 1; i++ )
  {
    kc = i * nx + j;
    a[kc+kc*n] = a[kc+kc*n] + 1.0;
    rhs[kc] = 0.0;
  }
/*
  Lower boundary.
*/
  i = 0;
  for ( j = 0; j < nx; j++ )
  {
    kc = i * nx + j;
    a[kc+kc*n] = a[kc+kc*n] + 1.0;
    rhs[kc] = 0.0;
  }
/*
  Upper boundary.
*/
  i = ny - 1;
  for ( j = 0; j < nx; j++ )
  {
    kc = i * nx + j;
    a[kc+kc*n] = a[kc+kc*n] + 1.0;
    rhs[kc] = 0.0;
  }

  return;
}
/******************************************************************************/

double test01_f ( double x, double y )

/******************************************************************************/
/*
  Purpose:

    TEST01_F evaluates the heat source term.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    03 September 2013

  Author:

    John Burkardt

  Parameters:

    Input, double X, Y, the evaluation point.

    Output, double TEST01_F, the value of the heat source term at (X,Y).
*/
{
  double arg;
  double v;
  double value;

  v = 0.05;
  arg = ( pow ( x - 0.60, 2 ) + pow ( y - 0.80, 2 ) ) / pow ( v, 2 );
  value = 2000.0 * exp ( - arg );

  return value;
}
