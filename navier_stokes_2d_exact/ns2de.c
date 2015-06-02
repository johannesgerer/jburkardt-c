# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <string.h>
# include <time.h>

# include "ns2de.h"

/******************************************************************************/

void grid_2d ( int x_num, double x_lo, double x_hi, int y_num, double y_lo, 
  double y_hi, double x[], double y[] )

/******************************************************************************/
/*
  Purpose:

    GRID_2D returns a regular 2D grid.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 January 2015

  Author:

    John Burkardt

  Parameters:

    Input, int X_NUM, the number of X values to use.

    Input, double X_LO, X_HI, the range of X values.

    Input, int Y_NUM, the number of Y values to use.

    Input, double Y_LO, Y_HI, the range of Y values.

    Output, double X[X_NUM*Y_NUM], Y[X_NUM*Y_NUM], 
    the coordinates of the grid.
*/
{
  int i;
  int j;
  double xi;
  double yj;

  if ( x_num == 1 )
  {
    for ( j = 0; j < y_num; j++ )
    {
      for ( i = 0; i < x_num; i++ )
      {
        x[i+j*x_num] = ( x_lo + x_hi ) / 2.0;
      }
    }
  }
  else
  {
    for ( i = 0; i < x_num; i++ )
    {
      xi = ( ( double ) ( x_num - i - 1 ) * x_lo   
           + ( double ) (         i     ) * x_hi ) 
           / ( double ) ( x_num     - 1 );
      for ( j = 0; j < y_num; j++ )
      {
        x[i+j*x_num] = xi;
      }
    }
  }

  if ( y_num == 1 )
  {
    for ( j = 0; j < y_num; j++ )
    {
      for ( i = 0; i < x_num; i++ )
      {
        y[i+j*x_num] = ( y_lo + y_hi ) / 2.0;
      }
    }
  }
  else
  {
    for ( j = 0; j < y_num; j++ )
    {
      yj = ( ( double ) ( y_num - j - 1 ) * y_lo   
           + ( double ) (         j     ) * y_hi ) 
           / ( double ) ( y_num     - 1 );
      for ( i = 0; i < x_num; i++ )
      {
        y[i+j*x_num] = yj;
      }
    }
  }

  return;
}
/******************************************************************************/

void ns2de_gnuplot ( char *header, int n, double x[], double y[], double u[], 
  double v[], double s )

/******************************************************************************/
/*
  Purpose:

    NS2DE_GNUPLOT writes the velocity field to files for GNUPLOT.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    21 January 2015

  Author:

    John Burkardt

  Parameters:

    Input, char *HEADER, a header to be used to name the files.

    Input, int N, the number of evaluation points.

    Input, double X[N], Y[N], the coordinates of the evaluation points.

    Input, double U[N], V[N], the velocity components.

    Input, double S, a scale factor for the velocity vectors.
*/
{
  char command_filename[255];
  FILE *command_unit;
  char data_filename[255];
  FILE *data_unit;
  int i;
  char plot_filename[255];
/*
  Write the data file.
*/
  strcpy ( data_filename, header );
  strcat ( data_filename, "_data.txt" );

  data_unit = fopen ( data_filename, "wt" );

  for ( i = 0; i < n; i++ )
  {
    fprintf ( data_unit, "  %g  %g  %g  %g\n", x[i], y[i], s * u[i], s * v[i] );
  }

  fclose ( data_unit );

  printf ( "\n" );
  printf ( "  Data written to '%s'\n", data_filename );
/*
  Write the command file.
*/
  strcpy ( command_filename, header );
  strcat ( command_filename, "_commands.txt" );

  strcpy ( plot_filename, header );
  strcat ( plot_filename, ".png" );

  command_unit = fopen ( command_filename, "wt" );

  fprintf ( command_unit, "#  %s\n", command_filename );
  fprintf ( command_unit, "#\n" );
  fprintf ( command_unit, "set term png\n" );
  fprintf ( command_unit, "set output '%s'\n", plot_filename );
  fprintf ( command_unit, "#\n" );
  fprintf ( command_unit, "#  Add titles and labels.\n" );
  fprintf ( command_unit, "#\n" );
  fprintf ( command_unit, "set xlabel '<--- X --->'\n" );
  fprintf ( command_unit, "set ylabel '<--- Y --->'\n" );
  fprintf ( command_unit, "set title 'Navier-Stokes velocity field'\n" );
  fprintf ( command_unit, "unset key\n" );
  fprintf ( command_unit, "#\n" );
  fprintf ( command_unit, "#  Add grid lines.\n" );
  fprintf ( command_unit, "#\n" );
  fprintf ( command_unit, "set grid\n" );
  fprintf ( command_unit, "set size ratio -1\n" );
  fprintf ( command_unit, "#\n" );
  fprintf ( command_unit, "#  Timestamp the plot.\n" );
  fprintf ( command_unit, "#\n" );
  fprintf ( command_unit, "set timestamp\n" );
  fprintf ( command_unit, 
    "plot '%s' using 1:2:3:4 with vectors \\\n", data_filename );
  fprintf ( command_unit, "  head filled lt 2 linecolor rgb 'blue'\n" );
  fprintf ( command_unit, "quit\n" );

  fclose ( command_unit );

  printf ( "  Commands written to '%s'\n", command_filename );

  return;
}
/******************************************************************************/

double r8vec_amax ( int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_AMAX returns the maximum absolute value in an R8VEC.

  Discussion:

    An R8VEC is a vector of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 May 2009

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the array.

    Input, double A[N], the array.

    Output, double AMAX, the value of the entry
    of largest magnitude.
*/
{
  double amax;
  int i;

  amax = 0.0;
  for ( i = 0; i < n; i++ )
  {
    if ( amax < fabs ( a[i] ) )
    {
      amax = fabs ( a[i] );
    }
  }
  return amax;
}
/******************************************************************************/

double r8vec_amin ( int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_AMIN returns the minimum absolute value in an R8VEC.

  Discussion:

    An R8VEC is a vector of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 May 2009

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the array.

    Input, double A[N], the array.

    Output, double R8VEC_AMIN, the value of the entry
    of smallest magnitude.
*/
{
  double amin;
  int i;
  const double r8_huge = 1.79769313486231571E+308;

  amin = r8_huge;
  for ( i = 0; i < n; i++ )
  {
    if ( fabs ( a[i] ) < amin )
    {
      amin = fabs ( a[i] );
    }
  }

  return amin;
}
/******************************************************************************/

void r8vec_linspace ( int n, double a, double b, double x[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_LINSPACE creates a vector of linearly spaced values.

  Discussion:

    An R8VEC is a vector of R8's.

    4 points evenly spaced between 0 and 12 will yield 0, 4, 8, 12.
 
    In other words, the interval is divided into N-1 even subintervals,
    and the endpoints of intervals are used as the points.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    10 April 2014

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the vector.

    Input, double A, B, the first and last entries.

    Output, double X[N], a vector of linearly spaced data.
*/
{
  int i;

  if ( n == 1 )
  {
    x[0] = ( a + b ) / 2.0;
  }
  else
  {
    for ( i = 0; i < n; i++ )
    {
      x[i] = ( ( double ) ( n - 1 - i ) * a 
             + ( double ) (         i ) * b ) 
             / ( double ) ( n - 1     );
    }
  }
  return;
}
/******************************************************************************/

double r8vec_max ( int n, double r8vec[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_MAX returns the value of the maximum element in a R8VEC.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 May 2006

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the array.

    Input, double R8VEC[N], a pointer to the first entry of the array.

    Output, double R8VEC_MAX, the value of the maximum element.  This
    is set to 0.0 if N <= 0.
*/
{
  int i;
  double value;

  if ( n <= 0 )
  {
    value = 0.0;
    return value;
  }

  value = r8vec[0];

  for ( i = 1; i < n; i++ )
  {
    if ( value < r8vec[i] )
    {
      value = r8vec[i];
    }
  }
  return value;
}
/******************************************************************************/

double r8vec_min ( int n, double r8vec[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_MIN returns the value of the minimum element in a R8VEC.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 May 2006

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the array.

    Input, double R8VEC[N], the array to be checked.

    Output, double R8VEC_MIN, the value of the minimum element.
*/
{
  int i;
  double value;

  value = r8vec[0];

  for ( i = 1; i < n; i++ )
  {
    if ( r8vec[i] < value )
    {
      value = r8vec[i];
    }
  }
  return value;
}
/******************************************************************************/

double r8vec_norm_l2 ( int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_NORM_L2 returns the L2 norm of an R8VEC.

  Discussion:

    The vector L2 norm is defined as:

      R8VEC_NORM_L2 = sqrt ( sum ( 1 <= I <= N ) A(I)^2 ).

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 March 2003

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in A.

    Input, double A[N], the vector whose L2 norm is desired.

    Output, double R8VEC_NORM_L2, the L2 norm of A.
*/
{
  int i;
  double v;

  v = 0.0;

  for ( i = 0; i < n; i++ )
  {
    v = v + a[i] * a[i];
  }
  v = sqrt ( v );

  return v;
}
/******************************************************************************/

double *r8vec_uniform_ab_new ( int n, double a, double b, int *seed )

/******************************************************************************/
/*
  Purpose:

    R8VEC_UNIFORM_AB_NEW returns a scaled pseudorandom R8VEC.

  Discussion:

    Each dimension ranges from A to B.

    This routine implements the recursion

      seed = 16807 * seed mod ( 2^31 - 1 )
      unif = seed / ( 2^31 - 1 )

    The integer arithmetic never requires more than 32 bits,
    including a sign bit.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    30 January 2005

  Author:

    John Burkardt

  Reference:

    Paul Bratley, Bennett Fox, Linus Schrage,
    A Guide to Simulation,
    Second Edition,
    Springer, 1987,
    ISBN: 0387964673,
    LC: QA76.9.C65.B73.

    Bennett Fox,
    Algorithm 647:
    Implementation and Relative Efficiency of Quasirandom
    Sequence Generators,
    ACM Transactions on Mathematical Software,
    Volume 12, Number 4, December 1986, pages 362-376.

    Pierre L'Ecuyer,
    Random Number Generation,
    in Handbook of Simulation,
    edited by Jerry Banks,
    Wiley, 1998,
    ISBN: 0471134031,
    LC: T57.62.H37.

    Peter Lewis, Allen Goodman, James Miller,
    A Pseudo-Random Number Generator for the System/360,
    IBM Systems Journal,
    Volume 8, Number 2, 1969, pages 136-143.

  Parameters:

    Input, int N, the number of entries in the vector.

    Input, double A, B, the lower and upper limits of the pseudorandom values.

    Input/output, int *SEED, a seed for the random number generator.

    Output, double R8VEC_UNIFORM_AB_NEW[N], the vector of pseudorandom values.
*/
{
  int i;
  const int i4_huge = 2147483647;
  int k;
  double *r;

  if ( *seed == 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8VEC_UNIFORM_AB_NEW - Fatal error!\n" );
    fprintf ( stderr, "  Input value of SEED = 0.\n" );
    exit ( 1 );
  }

  r = ( double * ) malloc ( n * sizeof ( double ) );

  for ( i = 0; i < n; i++ )
  {
    k = *seed / 127773;

    *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

    if ( *seed < 0 )
    {
      *seed = *seed + i4_huge;
    }

    r[i] = a + ( b - a ) * ( double ) ( *seed ) * 4.656612875E-10;
  }

  return r;
}
/******************************************************************************/

void resid_lucas ( double nu, double rho, int n, double x[], double y[], 
  double t, double ur[], double vr[], double pr[] )

/******************************************************************************/
/*
  Purpose:

    RESID_LUCAS returns Lucas Bystricky residuals.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    06 March 2015

  Author:

    John Burkardt

  Parameters:

    Input, double NU, the kinematic viscosity.

    Input, double RHO, the density.

    Input, int N, the number of evaluation points.

    Input, double X[N], Y[N], the coordinates of the points.

    Input, double T, the time coordinate or coordinates.

    Output, double UR[N], VR[N], PR[N], the residuals in the U, 
    V and P equations.
*/
{
  double dpdx;
  double dpdy;
  double dudt;
  double dudx;
  double dudxx;
  double dudy;
  double dudyy;
  double dvdt;
  double dvdx;
  double dvdxx;
  double dvdy;
  double dvdyy;
  double *f;
  double *g;
  double *h;
  int i;
  double p;
  const double r8_pi = 3.141592653589793;
  double u;
  double v;
/*
  Get the right hand sides.
*/
  f = ( double * ) malloc ( n * sizeof ( double ) );
  g = ( double * ) malloc ( n * sizeof ( double ) );
  h = ( double * ) malloc ( n * sizeof ( double ) );

  rhs_lucas ( nu, rho, n, x, y, t, f, g, h );
/*
  Form the functions and derivatives of the left hand side.
*/
  for ( i = 0; i < n; i++ )
  {
    u = - cos ( r8_pi * x[i] ) / r8_pi;
    dudt = 0.0;
    dudx = sin ( r8_pi * x[i] );
    dudxx = r8_pi * cos ( r8_pi * x[i] );
    dudy = 0.0;
    dudyy = 0.0;

    v = - y[i] * sin ( r8_pi * x[i] );
    dvdt = 0.0;
    dvdx = - r8_pi * y[i] * cos ( r8_pi * x[i] );
    dvdxx = + r8_pi * r8_pi * y[i] * sin ( r8_pi * x[i] );
    dvdy = - sin ( r8_pi * x[i] );
    dvdyy = 0.0;

    p = 0.0;
    dpdx = 0.0;
    dpdy = 0.0;
/*
  Evaluate the residuals.
*/
    ur[i] = dudt + u * dudx + v * dudy 
      + ( 1.0 / rho ) * dpdx - nu * ( dudxx + dudyy ) - f[i];

    vr[i] = dvdt + u * dvdx + v * dvdy 
      + ( 1.0 / rho ) * dpdy - nu * ( dvdxx + dvdyy ) - g[i];

    pr[i] = dudx + dvdy - h[i];
  }
/*
  Free memory.
*/
  free ( f );
  free ( g );
  free ( h );

  return;
}
/******************************************************************************/

void resid_spiral ( double nu, double rho, int n, double x[], double y[], double t, 
  double ur[], double vr[], double pr[] )

/******************************************************************************/
/*
  Purpose:

    RHS_SPIRAL evaluates the residual of the spiral flow problem.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 January 2015

  Author:

    John Burkardt

  Reference:

    Maxim Olshanskii, Leo Rebholz,
    Application of barycenter refined meshes in linear elasticity
    and incompressible fluid dynamics,
    ETNA: Electronic Transactions in Numerical Analysis, 
    Volume 38, pages 258-274, 2011.

  Parameters:

    Input, double NU, the kinematic viscosity.

    Input, double RHO, the fluid density.

    Input, int N, the number of nodes.

    Input, double X[N], Y[N], the coordinates of nodes.

    Input, double T, the current time.

    Output, double UR[N], VR[N], PR[N], the residuals sides.
*/
{
  double *dpdx;
  double *dpdy;
  double *dudt;
  double *dudx;
  double *dudxx;
  double *dudy;
  double *dudyy;
  double *dvdt;
  double *dvdx;
  double *dvdxx;
  double *dvdy;
  double *dvdyy;
  double *f;
  double *g;
  double *h;
  int i;
  double *p;
  double *u;
  double *v;

  dpdx = ( double * ) malloc ( n * sizeof ( double ) );
  dpdy = ( double * ) malloc ( n * sizeof ( double ) );
  dudt = ( double * ) malloc ( n * sizeof ( double ) );
  dudx = ( double * ) malloc ( n * sizeof ( double ) );
  dudxx = ( double * ) malloc ( n * sizeof ( double ) );
  dudy = ( double * ) malloc ( n * sizeof ( double ) );
  dudyy = ( double * ) malloc ( n * sizeof ( double ) );
  dvdt = ( double * ) malloc ( n * sizeof ( double ) );
  dvdx = ( double * ) malloc ( n * sizeof ( double ) );
  dvdxx = ( double * ) malloc ( n * sizeof ( double ) );
  dvdy = ( double * ) malloc ( n * sizeof ( double ) );
  dvdyy = ( double * ) malloc ( n * sizeof ( double ) );
  f = ( double * ) malloc ( n * sizeof ( double ) );
  g = ( double * ) malloc ( n * sizeof ( double ) );
  h = ( double * ) malloc ( n * sizeof ( double ) );
  p = ( double * ) malloc ( n * sizeof ( double ) );
  u = ( double * ) malloc ( n * sizeof ( double ) );
  v = ( double * ) malloc ( n * sizeof ( double ) );
/*
  Get the right hand side functions.
*/
  rhs_spiral ( nu, rho, n, x, y, t, f, g, h );
/*
  Form the functions and derivatives for the left hand side.
*/
  for ( i = 0; i < n; i++ )
  {
    u[i] = ( 1.0 + nu * t ) * 2.0 
      * ( pow ( x[i], 4 ) - 2.0 * pow ( x[i], 3 ) + pow ( x[i], 2 ) )
      * ( 2.0 * pow ( y[i], 3 ) - 3.0 * pow ( y[i], 2 ) + y[i] );

    dudt[i] = nu * 2.0 
      * ( pow ( x[i], 4 ) - 2.0 * pow ( x[i], 3 ) + pow ( x[i], 2 ) )
      * ( 2.0 * pow ( y[i], 3 ) - 3.0 * pow ( y[i], 2 ) + y[i] );

    dudx[i] = ( 1.0 + nu * t ) * 2.0 
      * ( 4.0 * pow ( x[i], 3 ) - 6.0 * pow ( x[i], 2 ) + 2.0 * x[i] )
      * ( 2.0 * pow ( y[i], 3 ) - 3.0 * pow ( y[i], 2 ) + y[i] );

    dudxx[i] = ( 1.0 + nu * t ) * 2.0 
      * ( 12.0 * pow ( x[i], 2 ) - 12.0 * x[i] + 2.0 )
      * ( 2.0 * pow ( y[i], 3 ) - 3.0 * pow ( y[i], 2 ) + y[i] );

    dudy[i] = ( 1.0 + nu * t ) * 2.0 
      * ( pow ( x[i], 4 ) - 2.0 * pow ( x[i], 3 ) + pow ( x[i], 2 ) )
      * ( 6.0 * pow ( y[i], 2 ) - 6.0 * y[i] + 1.0 );

    dudyy[i] = ( 1.0 + nu * t ) * 2.0 
      * ( pow ( x[i], 4 ) - 2.0 * pow ( x[i], 3 ) + pow ( x[i], 2 ) )
      * ( 12.0 * y[i] - 6.0 );

    v[i] = - ( 1.0 + nu * t ) * 2.0 
      * ( 2.0 * pow ( x[i], 3 ) - 3.0 * pow ( x[i], 2 ) + x[i] )  
      * ( pow ( y[i], 4 ) - 2.0 * pow ( y[i], 3 ) + pow ( y[i], 2 ) );

    dvdt[i] = - nu * 2.0 
      * ( 2.0 * pow ( x[i], 3 ) - 3.0 * pow ( x[i], 2 ) + x[i] )  
      * ( pow ( y[i], 4 ) - 2.0 * pow ( y[i], 3 ) + pow ( y[i], 2 ) );

    dvdx[i] = - ( 1.0 + nu * t ) * 2.0 
      * ( 6.0 * pow ( x[i], 2 ) - 6.0 * x[i] + 1.0 )  
      * ( pow ( y[i], 4 ) - 2.0 * pow ( y[i], 3 ) + pow ( y[i], 2 ) );

    dvdxx[i] = - ( 1.0 + nu * t ) * 2.0 
      * ( 12.0 * x[i] - 6.0 )  
      * ( pow ( y[i], 4 ) - 2.0 * pow ( y[i], 3 ) + pow ( y[i], 2 ) );

    dvdy[i] = - ( 1.0 + nu * t ) * 2.0 
      * ( 2.0 * pow ( x[i], 3 ) - 3.0 * pow ( x[i], 2 ) + x[i] )  
      * ( 4.0 * pow ( y[i], 3 ) - 6.0 * pow ( y[i], 2 ) + 2.0 * y[i] );

    dvdyy[i] = - ( 1.0 + nu * t ) * 2.0 
      * ( 2.0 * pow ( x[i], 3 ) - 3.0 * pow ( x[i], 2 ) + x[i] )  
      * ( 12.0 * pow ( y[i], 2 ) - 12.0 * y[i] + 2.0 );

    p[i] = rho * y[i];
    dpdx[i] = 0.0;
    dpdy[i] = rho;
/*
  Evaluate the residuals.
*/
    ur[i] = dudt[i] - nu * ( dudxx[i] + dudyy[i] ) 
      + u[i] * dudx[i] + v[i] * dudy[i] + dpdx[i] / rho - f[i];

    vr[i] = dvdt[i] - nu * ( dvdxx[i] + dvdyy[i] ) 
      + u[i] * dvdx[i] + v[i] * dvdy[i] + dpdy[i] / rho - g[i];

    pr[i] = dudx[i] + dvdy[i] - h[i];
  }

  free ( dpdx );
  free ( dpdy );
  free ( dudt );
  free ( dudx );
  free ( dudxx );
  free ( dudy );
  free ( dudyy );
  free ( dvdt );
  free ( dvdx );
  free ( dvdxx );
  free ( dvdy );
  free ( dvdyy );
  free ( f );
  free ( g );
  free ( h );
  free ( p );
  free ( u );
  free ( v );

  return;
}
/******************************************************************************/

void resid_taylor ( double nu, double rho, int n, double x[], double y[], 
  double t, double ur[], double vr[], double pr[] )

/******************************************************************************/
/*
  Purpose:

    RESID_TAYLOR returns residuals of the Taylor exact Navier Stokes solution.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 January 2015

  Author:

    John Burkardt

  Reference:

    Geoffrey Taylor,
    On the decay of vortices in a viscous fluid,
    Philosophical Magazine,
    Volume 46, 1923, pages 671-674.

    Geoffrey Taylor, A E Green,
    Mechanism for the production of small eddies from large ones,
    Proceedings of the Royal Society of London, 
    Series A, Volume 158, 1937, pages 499-521.

  Parameters:

    Input, double NU, the kinematic viscosity.

    Input, double RHO, the density.

    Input, int N, the number of evaluation points.

    Input, double X[N], Y[N], the coordinates of the points.

    Input, double T, the time coordinate or coordinates.

    Output, double UR[N], VR[N], PR[N], the residuals in the U, 
    V and P equations.
*/
{
  double dpdx;
  double dpdy;
  double dudt;
  double dudx;
  double dudxx;
  double dudy;
  double dudyy;
  double dvdt;
  double dvdx;
  double dvdxx;
  double dvdy;
  double dvdyy;
  double *f;
  double *g;
  double *h;
  int i;
  double p;
  const double r8_pi = 3.141592653589793;
  double u;
  double v;

/*
  Get the right hand sides.
*/
  f = ( double * ) malloc ( n * sizeof ( double ) );
  g = ( double * ) malloc ( n * sizeof ( double ) );
  h = ( double * ) malloc ( n * sizeof ( double ) );

  rhs_taylor ( nu, rho, n, x, y, t, f, g, h );
/*
  Form the functions and derivatives of the left hand side.
*/
  for ( i = 0; i < n; i++ )
  {
    u  =  -                          
        cos ( r8_pi * x[i] ) * sin ( r8_pi * y[i] );
    dudx =                       r8_pi 
      * sin ( r8_pi * x[i] ) * sin ( r8_pi * y[i] );
    dudxx =              r8_pi * r8_pi 
      * cos ( r8_pi * x[i] ) * sin ( r8_pi * y[i] );
    dudy =  -                    r8_pi 
      * cos ( r8_pi * x[i] ) * cos ( r8_pi * y[i] );
    dudyy =              r8_pi * r8_pi 
      * cos ( r8_pi * x[i] ) * sin ( r8_pi * y[i] );
    dudt =  + 2.0 * nu * r8_pi * r8_pi 
      * cos ( r8_pi * x[i] ) * sin ( r8_pi * y[i] );

    v  =                             
        sin ( r8_pi * x[i] ) * cos ( r8_pi * y[i] );
    dvdx =                       r8_pi 
      * cos ( r8_pi * x[i] ) * cos ( r8_pi * y[i] );
    dvdxx = -            r8_pi * r8_pi 
      * sin ( r8_pi * x[i] ) * cos ( r8_pi * y[i] );
    dvdy =  -                    r8_pi 
      * sin ( r8_pi * x[i] ) * sin ( r8_pi * y[i] );
    dvdyy = -            r8_pi * r8_pi 
      * sin ( r8_pi * x[i] ) * cos ( r8_pi * y[i] );
    dvdt =  - 2.0 * nu * r8_pi * r8_pi 
      * sin ( r8_pi * x[i] ) * cos ( r8_pi * y[i] );

    p =   - 0.25 * rho * 
      ( cos ( 2.0 * r8_pi * x[i] ) + cos ( 2.0 * r8_pi * y[i] ) );
    dpdx =  + 0.5  * rho * r8_pi * sin ( 2.0 * r8_pi * x[i] );
    dpdy =  + 0.5  * rho * r8_pi * sin ( 2.0 * r8_pi * y[i] );
/*
  Time scaling.
*/
    u     = u     * exp ( - 2.0 * r8_pi * r8_pi * nu * t );
    dudx  = dudx  * exp ( - 2.0 * r8_pi * r8_pi * nu * t );
    dudxx = dudxx * exp ( - 2.0 * r8_pi * r8_pi * nu * t );
    dudy  = dudy  * exp ( - 2.0 * r8_pi * r8_pi * nu * t );
    dudyy = dudyy * exp ( - 2.0 * r8_pi * r8_pi * nu * t );
    dudt  = dudt  * exp ( - 2.0 * r8_pi * r8_pi * nu * t );

    v     = v     * exp ( - 2.0 * r8_pi * r8_pi * nu * t );
    dvdx  = dvdx  * exp ( - 2.0 * r8_pi * r8_pi * nu * t );
    dvdxx = dvdxx * exp ( - 2.0 * r8_pi * r8_pi * nu * t );
    dvdy  = dvdy  * exp ( - 2.0 * r8_pi * r8_pi * nu * t );
    dvdyy = dvdyy * exp ( - 2.0 * r8_pi * r8_pi * nu * t );
    dvdt  = dvdt  * exp ( - 2.0 * r8_pi * r8_pi * nu * t );

    p =     p     * exp ( - 4.0 * r8_pi * r8_pi * nu * t );
    dpdx =  dpdx  * exp ( - 4.0 * r8_pi * r8_pi * nu * t );
    dpdy =  dpdy  * exp ( - 4.0 * r8_pi * r8_pi * nu * t );
/*
  Evaluate the residuals.
*/
    ur[i] = dudt + u * dudx + v * dudy 
      + ( 1.0 / rho ) * dpdx - nu * ( dudxx + dudyy ) - f[i];

    vr[i] = dvdt + u * dvdx + v * dvdy 
      + ( 1.0 / rho ) * dpdy - nu * ( dvdxx + dvdyy ) - g[i];

    pr[i] = dudx + dvdy - h[i];
  }
/*
  Free memory.
*/
  free ( f );
  free ( g );
  free ( h );

  return;
}
/******************************************************************************/

void rhs_lucas ( double nu, double rho, int n, double x[], double y[], double t, 
  double f[], double g[], double h[] )

/******************************************************************************/
/*
  Purpose:

    RHS_LUCAS evaluates the right hand side of Lucas Bystricky's problem.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    06 March 2015

  Author:

    John Burkardt

  Parameters:

    Input, double NU, the kinematic viscosity.

    Input, double RHO, the fluid density.

    Input, int N, the number of nodes.

    Input, double X[N], Y[N], the coordinates of nodes.

    Input, double T, the current time.

    Output, double F[N], G[N], H[N], the right hand sides.
*/
{
  double *dpdx;
  double *dpdy;
  double *dudt;
  double *dudx;
  double *dudxx;
  double *dudy;
  double *dudyy;
  double *dvdt;
  double *dvdx;
  double *dvdxx;
  double *dvdy;
  double *dvdyy;
  int i;
  double *p;
  const double r8_pi = 3.141592653589793;
  double *u;
  double *v;

  dpdx = ( double * ) malloc ( n * sizeof ( double ) );
  dpdy = ( double * ) malloc ( n * sizeof ( double ) );
  dudt = ( double * ) malloc ( n * sizeof ( double ) );
  dudx = ( double * ) malloc ( n * sizeof ( double ) );
  dudxx = ( double * ) malloc ( n * sizeof ( double ) );
  dudy = ( double * ) malloc ( n * sizeof ( double ) );
  dudyy = ( double * ) malloc ( n * sizeof ( double ) );
  dvdt = ( double * ) malloc ( n * sizeof ( double ) );
  dvdx = ( double * ) malloc ( n * sizeof ( double ) );
  dvdxx = ( double * ) malloc ( n * sizeof ( double ) );
  dvdy = ( double * ) malloc ( n * sizeof ( double ) );
  dvdyy = ( double * ) malloc ( n * sizeof ( double ) );
  p = ( double * ) malloc ( n * sizeof ( double ) );
  u = ( double * ) malloc ( n * sizeof ( double ) );
  v = ( double * ) malloc ( n * sizeof ( double ) );

  for ( i = 0; i < n; i++ )
  {
    u[i] = - cos ( r8_pi * x[i] ) / r8_pi;
    dudt[i] = 0.0;
    dudx[i] = sin ( r8_pi * x[i] );
    dudxx[i] = r8_pi * cos ( r8_pi * x[i] );
    dudy[i] = 0.0;
    dudyy[i] = 0.0;

    v[i] = - y[i] * sin ( r8_pi * x[i] );
    dvdt[i] = 0.0;
    dvdx[i] = - r8_pi * y[i] * cos ( r8_pi * x[i] );
    dvdxx[i] = + r8_pi * r8_pi * y[i] * sin ( r8_pi * x[i] );
    dvdy[i] = - sin ( r8_pi * x[i] );
    dvdyy[i] = 0.0;

    p[i] = 0.0;
    dpdx[i] = 0.0;
    dpdy[i] = 0.0;

    f[i] = dudt[i] - nu * ( dudxx[i] + dudyy[i] ) 
      + u[i] * dudx[i] + v[i] * dudy[i] + dpdx[i] / rho;

    g[i] = dvdt[i] - nu * ( dvdxx[i] + dvdyy[i] ) 
      + u[i] * dvdx[i] + v[i] * dvdy[i] + dpdy[i] / rho;

    h[i] = dudx[i] + dvdy[i];
  }

  free ( dpdx );
  free ( dpdy );
  free ( dudt );
  free ( dudx );
  free ( dudxx );
  free ( dudy );
  free ( dudyy );
  free ( dvdt );
  free ( dvdx );
  free ( dvdxx );
  free ( dvdy );
  free ( dvdyy );
  free ( p );
  free ( u );
  free ( v );

  return;
}
/******************************************************************************/

void rhs_spiral ( double nu, double rho, int n, double x[], double y[], double t, 
  double f[], double g[], double h[] )

/******************************************************************************/
/*
  Purpose:

    RHS_SPIRAL evaluates the right hand side of the spiral flow problem.

  Discussion:

    The right hand side is artificially determined by the requirement
    that the specified values of U, V and P satisfy the discretized
    Navier Stokes and continuity equations.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 January 2015

  Author:

    John Burkardt

  Reference:

    Maxim Olshanskii, Leo Rebholz,
    Application of barycenter refined meshes in linear elasticity
    and incompressible fluid dynamics,
    ETNA: Electronic Transactions in Numerical Analysis, 
    Volume 38, pages 258-274, 2011.

  Parameters:

    Input, double NU, the kinematic viscosity.

    Input, double RHO, the fluid density.

    Input, int N, the number of nodes.

    Input, double X[N], Y[N], the coordinates of nodes.

    Input, double T, the current time.

    Output, double F[N], G[N], H[N], the right hand sides.
*/
{
  double *dpdx;
  double *dpdy;
  double *dudt;
  double *dudx;
  double *dudxx;
  double *dudy;
  double *dudyy;
  double *dvdt;
  double *dvdx;
  double *dvdxx;
  double *dvdy;
  double *dvdyy;
  int i;
  double *p;
  double *u;
  double *v;

  dpdx = ( double * ) malloc ( n * sizeof ( double ) );
  dpdy = ( double * ) malloc ( n * sizeof ( double ) );
  dudt = ( double * ) malloc ( n * sizeof ( double ) );
  dudx = ( double * ) malloc ( n * sizeof ( double ) );
  dudxx = ( double * ) malloc ( n * sizeof ( double ) );
  dudy = ( double * ) malloc ( n * sizeof ( double ) );
  dudyy = ( double * ) malloc ( n * sizeof ( double ) );
  dvdt = ( double * ) malloc ( n * sizeof ( double ) );
  dvdx = ( double * ) malloc ( n * sizeof ( double ) );
  dvdxx = ( double * ) malloc ( n * sizeof ( double ) );
  dvdy = ( double * ) malloc ( n * sizeof ( double ) );
  dvdyy = ( double * ) malloc ( n * sizeof ( double ) );
  p = ( double * ) malloc ( n * sizeof ( double ) );
  u = ( double * ) malloc ( n * sizeof ( double ) );
  v = ( double * ) malloc ( n * sizeof ( double ) );

  for ( i = 0; i < n; i++ )
  {
    u[i] = ( 1.0 + nu * t ) * 2.0 
      * ( pow ( x[i], 4 ) - 2.0 * pow ( x[i], 3 ) + pow ( x[i], 2 ) )
      * ( 2.0 * pow ( y[i], 3 ) - 3.0 * pow ( y[i], 2 ) + y[i] );

    dudt[i] = nu * 2.0 
      * ( pow ( x[i], 4 ) - 2.0 * pow ( x[i], 3 ) + pow ( x[i], 2 ) )
      * ( 2.0 * pow ( y[i], 3 ) - 3.0 * pow ( y[i], 2 ) + y[i] );

    dudx[i] = ( 1.0 + nu * t ) * 2.0 
      * ( 4.0 * pow ( x[i], 3 ) - 6.0 * pow ( x[i], 2 ) + 2.0 * x[i] )
      * ( 2.0 * pow ( y[i], 3 ) - 3.0 * pow ( y[i], 2 ) + y[i] );

    dudxx[i] = ( 1.0 + nu * t ) * 2.0 
      * ( 12.0 * pow ( x[i], 2 ) - 12.0 * x[i] + 2.0 )
      * ( 2.0 * pow ( y[i], 3 ) - 3.0 * pow ( y[i], 2 ) + y[i] );

    dudy[i] = ( 1.0 + nu * t ) * 2.0 
      * ( pow ( x[i], 4 ) - 2.0 * pow ( x[i], 3 ) + pow ( x[i], 2 ) )
      * ( 6.0 * pow ( y[i], 2 ) - 6.0 * y[i] + 1.0 );

    dudyy[i] = ( 1.0 + nu * t ) * 2.0 
      * ( pow ( x[i], 4 ) - 2.0 * pow ( x[i], 3 ) + pow ( x[i], 2 ) )
      * ( 12.0 * y[i] - 6.0 );

    v[i] = - ( 1.0 + nu * t ) * 2.0 
      * ( 2.0 * pow ( x[i], 3 ) - 3.0 * pow ( x[i], 2 ) + x[i] )  
      * ( pow ( y[i], 4 ) - 2.0 * pow ( y[i], 3 ) + pow ( y[i], 2 ) );

    dvdt[i] = - nu * 2.0 
      * ( 2.0 * pow ( x[i], 3 ) - 3.0 * pow ( x[i], 2 ) + x[i] )  
      * ( pow ( y[i], 4 ) - 2.0 * pow ( y[i], 3 ) + pow ( y[i], 2 ) );

    dvdx[i] = - ( 1.0 + nu * t ) * 2.0 
      * ( 6.0 * pow ( x[i], 2 ) - 6.0 * x[i] + 1.0 )  
      * ( pow ( y[i], 4 ) - 2.0 * pow ( y[i], 3 ) + pow ( y[i], 2 ) );

    dvdxx[i] = - ( 1.0 + nu * t ) * 2.0 
      * ( 12.0 * x[i] - 6.0 )  
      * ( pow ( y[i], 4 ) - 2.0 * pow ( y[i], 3 ) + pow ( y[i], 2 ) );

    dvdy[i] = - ( 1.0 + nu * t ) * 2.0 
      * ( 2.0 * pow ( x[i], 3 ) - 3.0 * pow ( x[i], 2 ) + x[i] )  
      * ( 4.0 * pow ( y[i], 3 ) - 6.0 * pow ( y[i], 2 ) + 2.0 * y[i] );

    dvdyy[i] = - ( 1.0 + nu * t ) * 2.0 
      * ( 2.0 * pow ( x[i], 3 ) - 3.0 * pow ( x[i], 2 ) + x[i] )  
      * ( 12.0 * pow ( y[i], 2 ) - 12.0 * y[i] + 2.0 );

    p[i] = rho * y[i];
    dpdx[i] = 0.0;
    dpdy[i] = rho;

    f[i] = dudt[i] - nu * ( dudxx[i] + dudyy[i] ) 
      + u[i] * dudx[i] + v[i] * dudy[i] + dpdx[i] / rho;

    g[i] = dvdt[i] - nu * ( dvdxx[i] + dvdyy[i] ) 
      + u[i] * dvdx[i] + v[i] * dvdy[i] + dpdy[i] / rho;

    h[i] = dudx[i] + dvdy[i];
  }

  free ( dpdx );
  free ( dpdy );
  free ( dudt );
  free ( dudx );
  free ( dudxx );
  free ( dudy );
  free ( dudyy );
  free ( dvdt );
  free ( dvdx );
  free ( dvdxx );
  free ( dvdy );
  free ( dvdyy );
  free ( p );
  free ( u );
  free ( v );

  return;
}
/******************************************************************************/

void rhs_taylor ( double nu, double rho, int n, double x[], double y[], 
  double t, double f[], double g[], double h[] )

/******************************************************************************/
/*
  Purpose:

    RHS_TAYLOR returns the right hand sides of the Taylor vortex equations.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 January 2015

  Author:

    John Burkardt

  Reference:

    Geoffrey Taylor,
    On the decay of vortices in a viscous fluid,
    Philosophical Magazine,
    Volume 46, 1923, pages 671-674.

    Geoffrey Taylor, A E Green,
    Mechanism for the production of small eddies from large ones,
    Proceedings of the Royal Society of London, 
    Series A, Volume 158, 1937, pages 499-521.

  Parameters:

    Input, double NU, the kinematic viscosity.

    Input, double RHO, the density.

    Input, int N, the number of evaluation points.

    Input, double X[N], Y[N], the coordinates of the points.

    Input, double T, the time coordinate or coordinates.

    Output, double F[N], G[N], H[N], the residuals in the U, 
    V and P equations.
*/
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    f[i] = 0.0;
    g[i] = 0.0;
    h[i] = 0.0;
  }

  return;
}
/******************************************************************************/

void timestamp ( )

/******************************************************************************/
/*
  Purpose:

    TIMESTAMP prints the current YMDHMS date as a time stamp.

  Example:

    31 May 2001 09:45:54 AM

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    24 September 2003

  Author:

    John Burkardt

  Parameters:

    None
*/
{
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct tm *tm;
  size_t len;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  len = strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  fprintf ( stdout, "%s\n", time_buffer );

  return;
# undef TIME_SIZE
}
/******************************************************************************/

void uvp_lucas ( double nu, double rho, int n, double x[], double y[], 
  double t, double u[], double v[], double p[] )

/******************************************************************************/
/*
  Purpose:

    UVP_LUCAS evaluates Lucas Bystricky's exact Navier Stokes solution.

  Discussion:

    There is no time dependence.

    The pressure is 0.

    The preferred domain is the unit square.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    06 March 2015

  Author:

    John Burkardt

  Parameters:

    Input, double NU, the kinematic viscosity.

    Input, double RHO, the density.

    Input, int N, the number of evaluation points.

    Input, double X[N], Y[N], the coordinates of the points.

    Input, double T, the time coordinate or coordinates.

    Output, double U[N], V[N], P[N], the velocity components and
    pressure at each of the points.
*/
{
  int i;
  const double r8_pi = 3.141592653589793;

  for ( i = 0; i < n; i++ )
  {
    u[i] = - cos ( r8_pi * x[i] ) / r8_pi;
    v[i] = - y[i] *  sin ( r8_pi * x[i] );
    p[i] = - 0.0;
  }

  return;
}
/******************************************************************************/

void uvp_spiral ( double nu, double rho, int n, double x[], double y[], 
  double t, double u[], double v[], double p[] )

/******************************************************************************/
/*
  Purpose:

    UVP_SPIRAL returns velocity and pressure for the spiral flow.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    25 February 2011

  Author:

    John Burkardt

  Reference:

    Maxim Olshanskii, Leo Rebholz,
    Application of barycenter refined meshes in linear elasticity
    and incompressible fluid dynamics,
    ETNA: Electronic Transactions in Numerical Analysis, 
    Volume 38, pages 258-274, 2011.

  Parameters:

    Input, double NU, the kinematic viscosity.

    Input, double RHO, the fluid density.

    Input, int N, the number of nodes.

    Input, double X[N], Y[N], the coordinates of nodes.

    Input, double T, the current time.

    Output, double U[N], V[N], the X and Y velocity.

    Output, double P[N], the pressure.
*/
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    u[i] = ( 1.0 + nu * t ) * 2.0 
      * pow ( x[i], 2 ) * pow ( x[i] - 1.0, 2 )
      * y[i] * ( 2.0 * y[i] - 1.0 ) * ( y[i] - 1.0 );

    v[i] = - ( 1.0 + nu * t ) * 2.0 
      * x[i] * ( 2.0 * x[i] - 1.0 ) * ( x[i] - 1.0 )  
      * pow ( y[i], 2 ) * pow ( y[i] - 1.0, 2 );

    p[i] = y[i];
  }

  return;
}
/******************************************************************************/

void uvp_taylor ( double nu, double rho, int n, double x[], double y[], 
  double t, double u[], double v[], double p[] )

/******************************************************************************/
/*
  Purpose:

    UVP_TAYLOR evaluates the Taylor exact Navier Stokes solution.

  Discussion:

    This flow is known as a Taylor-Green vortex.

    The given velocity and pressure fields are exact solutions for the 2D 
    incompressible time-dependent Navier Stokes equations over any region.

    To define a typical problem, one chooses a bounded spatial region 
    and a starting time, and then imposes boundary and initial conditions
    by referencing the exact solution appropriately.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 January 2015

  Author:

    John Burkardt

  Reference:

    Geoffrey Taylor,
    On the decay of vortices in a viscous fluid,
    Philosophical Magazine,
    Volume 46, 1923, pages 671-674.

    Geoffrey Taylor, A E Green,
    Mechanism for the production of small eddies from large ones,
    Proceedings of the Royal Society of London, 
    Series A, Volume 158, 1937, pages 499-521.

  Parameters:

    Input, double NU, the kinematic viscosity.

    Input, double RHO, the density.

    Input, int N, the number of evaluation points.

    Input, double X[N], Y[N], the coordinates of the points.

    Input, double T, the time coordinate or coordinates.

    Output, double U[N], V[N], P[N], the velocity components and
    pressure at each of the points.
*/
{
  int i;
  const double r8_pi = 3.141592653589793;

  for ( i = 0; i < n; i++ )
  {
    u[i] = - cos ( r8_pi * x[i] ) * sin ( r8_pi * y[i] );
    v[i] =   sin ( r8_pi * x[i] ) * cos ( r8_pi * y[i] );
    p[i] = - 0.25 * rho 
      * ( cos ( 2.0 * r8_pi * x[i] ) + cos ( 2.0 * r8_pi * y[i] ) );

    u[i] = u[i] * exp ( - 2.0 * r8_pi * r8_pi * nu * t );
    v[i] = v[i] * exp ( - 2.0 * r8_pi * r8_pi * nu * t );
    p[i] = p[i] * exp ( - 4.0 * r8_pi * r8_pi * nu * t );
  }

  return;
}

