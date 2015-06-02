# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>

# include "fd2d_heat_steady.h"

int main ( );
void test01 ( );
double d ( double x, double y );
double f ( double x, double y );
void boundary ( int nx, int ny, double x[], double y[], int n, double a[], 
  double rhs[] );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for FD2D_HEAT_STEADY_PRB.

  Discussion:

    FD2D_HEAT_STEADY_PRB tests the FD2D_HEAT_STEADY library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    29 August 2013

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "FD2D_HEAT_STEADY_PRB:\n" );
  printf ( "  C version\n" );
  printf ( "  Test the FD2D_HEAT_STEADY library.\n" );

  test01 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "FD2D_HEAT_STEADY_PRB:\n" );
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

    TEST01 computes the solution for a steady state heat equation problem.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 August 2013

  Author:

    John Burkardt
*/
{
  char command_filename[] = "test01_commands.txt";
  FILE *command_unit;
  char data_filename[] = "test01_data.txt";
  FILE *data_unit;
  int i;
  int j;
  int nx;
  int ny;
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
/*
  Specify the spatial grid.
*/
  nx = 41;
  xvec = r8vec_linspace_new ( nx, 0.0, 2.0 );

  ny = 21;
  yvec = r8vec_linspace_new ( ny, 0.0, 1.0 );

  xmat = ( double * ) malloc ( nx * ny * sizeof ( double ) );
  ymat = ( double * ) malloc ( nx * ny * sizeof ( double ) );
  r8vec_mesh_2d ( nx, ny, xvec, yvec, xmat, ymat );
/*
  Solve the finite difference approximation to the steady 2D heat equation.
*/
  umat = fd2d_heat_steady ( nx, ny, xvec, yvec, d, f );
/*
  Create a data file.
*/
  data_unit = fopen ( data_filename, "wt" );
  for ( j = 0; j < ny; j++ )
  {
    for ( i = 0; i < nx; i++ )
    {
      fprintf ( data_unit, "  %14g  %14g  %14g\n", 
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
  fprintf ( command_unit, "set output 'test01.png'\n" );
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
  u_mean = 0.0;
  for ( j = 0; j < ny; j++ )
  {
    for ( i = 0; i < nx; i++ )
    {
      u_mean = u_mean + umat[i+j*nx];
    }
  }
  u_mean = u_mean / ( double ) ( nx * ny );

  printf ( "\n" );
  printf ( "  Mean value of U is %g\n", u_mean );
/*
  Free memory.
*/
  free ( umat );
  free ( xmat );
  free ( xvec );
  free ( ymat );
  free ( yvec );

  return;
}
/******************************************************************************/

double d ( double x, double y )

/******************************************************************************/
/*
  Purpose:

    D evaluates the heat conductivity coefficient.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    23 July 2013

  Author:

    John Burkardt

  Parameters:

    Input, double X, Y, the evaluation point.

    Output, double D, the value of the heat conductivity at (X,Y).
*/
{
  double value;

  value = 1.0;

  return value;
}
/******************************************************************************/

double f ( double x, double y )

/******************************************************************************/
/*
  Purpose:

    F evaluates the heat source term.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    23 July 2013

  Author:

    John Burkardt

  Parameters:

    Input, double X, Y, the evaluation point.

    Output, double F, the value of the heat source term at (X,Y).
*/
{
  double value;

  value = 0.0;

  return value;
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

    29 August 2013

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
    rhs[kc] = 10.0;
  }
/*
  Right boundary.
*/
  j = nx - 1;
  for ( i = 1; i < ny - 1; i++ )
  {
    kc = i * nx + j;
    a[kc+kc*n] = a[kc+kc*n] + 1.0;
    rhs[kc] = 100.0;
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

