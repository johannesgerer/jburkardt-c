# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <string.h>
# include <time.h>

# include "s2de.h"

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

void resid_stokes1 ( int n, double x[], double y[], double ur[], double vr[], 
  double pr[] )

/******************************************************************************/
/*
  Purpose:

    RESID_STOKES1 returns residuals of the exact Stokes solution #1.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    23 January 2015

  Author:

    John Burkardt

  Reference:

    Junping Wang, Yanqiu Wang, Xiu Ye,
    A robust numerical method for Stokes equations based on divergence-free
    H(div) finite element methods,
    SIAM Journal on Scientific Computing,
    Volume 31, Number 4, 2009, pages 2784-2802.

  Parameters:

    Input, int N, the number of evaluation points.

    Input, double X[N], Y[N], the coordinates of the points.

    Output, double UR[N], VR[N], PR[N], the residuals in the U, 
    V and P equations.
*/
{
  double *f;
  double *g;
  double *h;
  int i;
  double p;
  double px;
  double py;
  double u;
  double ux;
  double uxx;
  double uy;
  double uyy;
  double v;
  double vx;
  double vxx;
  double vy;
  double vyy;
/*
  Get the right hand sides.
*/
  f = ( double * ) malloc ( n * sizeof ( double ) );
  g = ( double * ) malloc ( n * sizeof ( double ) );
  h = ( double * ) malloc ( n * sizeof ( double ) );

  rhs_stokes1 ( n, x, y, f, g, h );
/*
  Form the functions and derivatives.
*/
  for ( i = 0; i < n; i++ )
  {
    u = - 2.0 
          * pow ( x[i], 2 ) * pow ( x[i] - 1.0, 2 ) 
          * y[i] * ( y[i] - 1.0 ) * ( 2.0 * y[i] - 1.0 );

    ux = - 2.0 
          * ( 4.0 * pow ( x[i], 3 ) - 6.0 * pow ( x[i], 2 ) 
          + 2.0 * x[i] ) 
          * y[i] * ( y[i] - 1.0 ) * ( 2.0 * y[i] - 1.0 );

    uxx = - 2.0 
          * ( 12.0 * pow ( x[i], 2 ) - 12.0 * x[i] + 2.0 ) 
          * ( 2.0 * pow ( y[i], 3 ) - 3.0 * pow ( y[i], 2 ) + y[i] );

    uy = - 2.0  
          * pow ( x[i], 2 ) * pow ( x[i] - 1.0, 2 )  
          * ( 6.0 * pow ( y[i], 2 ) - 3.0 * y[i] + 1.0 );

    uyy = - 2.0 
          * ( pow ( x[i], 4 ) - 2.0 * pow ( x[i], 3 ) + pow ( x[i], 2 ) ) 
          * ( 12.0 * y[i] - 6.0 );

    v =   2.0 
          * x[i] * ( x[i] - 1.0 ) * ( 2.0 * x[i] - 1.0 ) 
          * pow ( y[i], 2 ) * pow ( y[i] - 1.0, 2 );

    vx =   2.0 
          * ( 6.0 * pow ( x[i], 2 ) - 6.0 * x[i] + 1.0 ) 
          * pow ( y[i], 2 ) * pow ( y[i] - 1.0, 2 );

    vxx =   2.0 
          * ( 12.0 * x[i] - 6.0 ) 
          * pow ( y[i], 2 ) * pow ( y[i] - 1.0, 2 );

    vy =   2.0 
          * x[i] * ( x[i] - 1.0 ) * ( 2.0 * x[i] - 1.0 ) 
          * ( 4.0 * pow ( y[i], 3 ) - 6.0 * pow ( y[i], 2 )  
          + 2.0 * y[i] );

    vyy =   2.0 
          * x[i] * ( x[i] - 1.0 ) * ( 2.0 * x[i] - 1.0 ) 
          * ( 12.0 * pow ( y[i], 2 ) - 12.0 * y[i] + 2.0 );

    p = 0.0;
    px = 0.0;
    py = 0.0;

    ur[i] = px - ( uxx + uyy ) - f[i];
    vr[i] = py - ( vxx + vyy ) - g[i];
    pr[i] = ux + vy - h[i];
  }
/*
  Deallocate memory.
*/
  free ( f );
  free ( g );
  free ( h );

  return;
}
/******************************************************************************/

void resid_stokes2 ( int n, double x[], double y[], double ur[], double vr[], 
  double pr[] )

/******************************************************************************/
/*
  Purpose:

    RESID_STOKES2 returns residuals of the exact Stokes solution #2.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    24 January 2015

  Author:

    John Burkardt

  Reference:

    Junping Wang, Yanqiu Wang, Xiu Ye,
    A robust numerical method for Stokes equations based on divergence-free
    H(div) finite element methods,
    SIAM Journal on Scientific Computing,
    Volume 31, Number 4, 2009, pages 2784-2802.

  Parameters:

    Input, int N, the number of evaluation points.

    Input, double X[N], Y[N], the coordinates of the points.

    Output, double UR[N], VR[N], PR[N], the residuals in the U, 
    V and P equations.
*/
{
  double *f;
  double *g;
  double *h;
  int i;
  double p;
  double px;
  double py;
  const double r8_pi = 3.141592653589793;
  double u;
  double ux;
  double uxx;
  double uy;
  double uyy;
  double v;
  double vx;
  double vxx;
  double vy;
  double vyy;
/*
  Get the right hand sides.
*/
  f = ( double * ) malloc ( n * sizeof ( double ) );
  g = ( double * ) malloc ( n * sizeof ( double ) );
  h = ( double * ) malloc ( n * sizeof ( double ) );

  rhs_stokes2 ( n, x, y, f, g, h );
/*
  Form the functions and derivatives.
*/
  for ( i = 0; i < n; i++ )
  {
    u =   2.0 
          * sin ( 2.0 * r8_pi * x[i] ) 
          * cos ( 2.0 * r8_pi * y[i] );

    ux =   4.0 * r8_pi 
          * cos ( 2.0 * r8_pi * x[i] ) 
          * cos ( 2.0 * r8_pi * y[i] );

    uxx = - 8.0 * pow ( r8_pi, 2 )
          * sin ( 2.0 * r8_pi * x[i] ) 
          * cos ( 2.0 * r8_pi * y[i] );

    uy = - 4.0 * r8_pi 
          * sin ( 2.0 * r8_pi * x[i] ) 
          * sin ( 2.0 * r8_pi * y[i] );

    uyy = - 8.0 * pow ( r8_pi, 2 )
          * sin ( 2.0 * r8_pi * x[i] ) 
          * cos ( 2.0 * r8_pi * y[i] );

    v = - 2.0 
          * cos ( 2.0 * r8_pi * x[i] ) 
          * sin ( 2.0 * r8_pi * y[i] );

    vx =   4.0 * r8_pi 
          * sin ( 2.0 * r8_pi * x[i] ) 
          * sin ( 2.0 * r8_pi * y[i] );

    vxx =   8.0 * pow ( r8_pi, 2 )
          * cos ( 2.0 * r8_pi * x[i] ) 
          * sin ( 2.0 * r8_pi * y[i] );

    vy = - 4.0 * r8_pi 
          * cos ( 2.0 * r8_pi * x[i] ) 
          * cos ( 2.0 * r8_pi * y[i] );

    vyy =   8.0 * pow ( r8_pi, 2 )
          * cos ( 2.0 * r8_pi * x[i] ) 
          * sin ( 2.0 * r8_pi * y[i] );

    p = pow ( x[i], 2 ) + pow ( y[i], 2 );

    px = 2.0 * x[i];
    py = 2.0 * y[i];

    ur[i] = px - ( uxx + uyy ) - f[i];
    vr[i] = py - ( vxx + vyy ) - g[i];
    pr[i] = ux + vy - h[i];
  }
/*
  Deallocate memory.
*/
  free ( f );
  free ( g );
  free ( h );

  return;
}
/******************************************************************************/

void resid_stokes3 ( int n, double x[], double y[], double ur[], double vr[], 
  double pr[] )

/******************************************************************************/
/*
  Purpose:

    RESID_STOKES3 returns residuals of the exact Stokes solution #3.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 February 2015

  Author:

    John Burkardt

  Reference:

    Howard Elman, Alison Ramage, David Silvester,
    Finite Elements and Fast Iterative Solvers with
    Applications in Incompressible Fluid Dynamics,
    Oxford, 2005,
    ISBN: 978-0198528678,
    LC: QA911.E39.

  Parameters:

    Input, int N, the number of evaluation points.

    Input, double X[N], Y[N], the coordinates of the points.

    Output, double UR[N], VR[N], PR[N], the residuals in the U, 
    V and P equations.
*/
{
  double *f;
  double *g;
  double *h;
  int i;
  double p;
  double px;
  double py;
  double u;
  double ux;
  double uxx;
  double uy;
  double uyy;
  double v;
  double vx;
  double vxx;
  double vy;
  double vyy;
/*
  Get the right hand sides.
*/
  f = ( double * ) malloc ( n * sizeof ( double ) );
  g = ( double * ) malloc ( n * sizeof ( double ) );
  h = ( double * ) malloc ( n * sizeof ( double ) );

  rhs_stokes3 ( n, x, y, f, g, h );
/*
  Form the functions and derivatives.
*/
  for ( i = 0; i < n; i++ )
  {
    u =   20.0 * x[i] * pow ( y[i], 3 );
    ux = 20.0 * pow ( y[i], 3 );
    uxx = 0.0;
    uy = 60.0 * x[i] * pow ( y[i], 2 );
    uyy = 120.0 * x[i] * y[i];

    v = 5.0 * ( pow ( x[i], 4 )  - pow ( y[i], 4 ) );
    vx = 20.0 * pow ( x[i], 3 );
    vxx = 60.0 * pow ( x[i], 2 );
    vy = - 20.0 * pow ( y[i], 3 );
    vyy = - 60.0 * pow ( y[i], 2 );

    p =   60.0 * pow ( x[i], 2 ) * y[i] - 20.0 * pow ( y[i], 3 ) + 10.0;
    px = 120.0 * x[i] * y[i];
    py =  60.0 * pow ( x[i], 2 ) - 60.0 * pow ( y[i], 2 );

    ur[i] = px - ( uxx + uyy ) - f[i];
    vr[i] = py - ( vxx + vyy ) - g[i];
    pr[i] = ux + vy - h[i];
  }
/*
  Deallocate memory.
*/
  free ( f );
  free ( g );
  free ( h );

  return;
}
/******************************************************************************/

void rhs_stokes1 ( int n, double x[], double y[], double f[], double g[], 
  double h[] )

/******************************************************************************/
/*
  Purpose:

    RHS_STOKES1 returns the right hand sides of the exact Stokes solution #1.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    24 January 2015

  Author:

    John Burkardt

  Reference:

    Junping Wang, Yanqiu Wang, Xiu Ye,
    A robust numerical method for Stokes equations based on divergence-free
    H(div) finite element methods,
    SIAM Journal on Scientific Computing,
    Volume 31, Number 4, 2009, pages 2784-2802.

  Parameters:

    Input, int N, the number of evaluation points.

    Input, double X[N], Y[N], the coordinates of the points.

    Output, double F[N], G[N], H[N], the right hand sides in the U,
    V and P equations.
*/
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    f[i] = + 2.0 
          * ( 12.0 * pow ( x[i], 2 ) - 12.0 * x[i] + 2.0 ) 
          * ( 2.0 * pow ( y[i], 3 ) - 3.0 * pow ( y[i], 2 ) + y[i] ) 
          + 2.0 
          * ( pow ( x[i], 4 ) - 2.0 * pow ( x[i], 3 ) + pow ( x[i], 2 ) ) 
          * ( 12.0 * y[i] - 6.0 );

    g[i] = - 2.0 
          * ( 12.0 * x[i] - 6.0 ) 
          * ( pow ( y[i], 4 ) - 2.0 * pow ( y[i], 3 ) + pow ( y[i], 2 ) ) 
          - 2.0 
          * ( 2.0 * pow ( x[i], 3 ) - 3.0 * pow ( x[i], 2 ) + x[i] ) 
          * ( 12.0 * pow ( y[i], 2 ) - 12.0 * y[i] + 2.0 );

    h[i] = 0.0;
  }

  return;
}
/******************************************************************************/

void rhs_stokes2 ( int n, double x[], double y[], double f[], double g[], 
  double h[] )

/******************************************************************************/
/*
  Purpose:

    RHS_STOKES2 returns the right hand sides of the exact Stokes solution #2.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    24 January 2015

  Author:

    John Burkardt

  Reference:

    Junping Wang, Yanqiu Wang, Xiu Ye,
    A robust numerical method for Stokes equations based on divergence-free
    H(div) finite element methods,
    SIAM Journal on Scientific Computing,
    Volume 31, Number 4, 2009, pages 2784-2802.

  Parameters:

    Input, int N, the number of evaluation points.

    Input, double X[N], Y[N], the coordinates of the points.

    Output, double F[N], G[N], H[N], the right hand sides in the U,
    V and P equations.
*/
{
  int i;
  double p;
  double px;
  double py;
  const double r8_pi = 3.141592653589793;
  double u;
  double ux;
  double uxx;
  double uy;
  double uyy;
  double v;
  double vx;
  double vxx;
  double vy;
  double vyy;

  for ( i = 0; i < n; i++ )
  {
    u =   2.0 
          * sin ( 2.0 * r8_pi * x[i] ) 
          * cos ( 2.0 * r8_pi * y[i] );

    ux =   4.0 * r8_pi 
          * cos ( 2.0 * r8_pi * x[i] ) 
          * cos ( 2.0 * r8_pi * y[i] );

    uxx = - 8.0 * pow ( r8_pi, 2 )
          * sin ( 2.0 * r8_pi * x[i] ) 
          * cos ( 2.0 * r8_pi * y[i] );

    uy = - 4.0 * r8_pi 
          * sin ( 2.0 * r8_pi * x[i] ) 
          * sin ( 2.0 * r8_pi * y[i] );

    uyy = - 8.0 * pow ( r8_pi, 2 )
          * sin ( 2.0 * r8_pi * x[i] ) 
          * cos ( 2.0 * r8_pi * y[i] );

    v = - 2.0 
          * cos ( 2.0 * r8_pi * x[i] ) 
          * sin ( 2.0 * r8_pi * y[i] );

    vx =   4.0 * r8_pi 
          * sin ( 2.0 * r8_pi * x[i] ) 
          * sin ( 2.0 * r8_pi * y[i] );

    vxx =   8.0 * pow ( r8_pi, 2 )
          * cos ( 2.0 * r8_pi * x[i] ) 
          * sin ( 2.0 * r8_pi * y[i] );

    vy = - 4.0 * r8_pi 
          * cos ( 2.0 * r8_pi * x[i] ) 
          * cos ( 2.0 * r8_pi * y[i] );

    vyy =   8.0 * pow ( r8_pi, 2 )
          * cos ( 2.0 * r8_pi * x[i] ) 
          * sin ( 2.0 * r8_pi * y[i] );

    p = pow ( x[i], 2 ) + pow ( y[i], 2 );

    px = 2.0 * x[i];
    py = 2.0 * y[i];

    f[i] = px - ( uxx + uyy );
    g[i] = py - ( vxx + vyy );
    h[i] = ux + vy;
  }

  return;
}
/******************************************************************************/

void rhs_stokes3 ( int n, double x[], double y[], double f[], double g[], 
  double h[] )

/******************************************************************************/
/*
  Purpose:

    RHS_STOKES3 returns the right hand sides of the exact Stokes solution #3.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 February 2015

  Author:

    John Burkardt

  Reference:

    Howard Elman, Alison Ramage, David Silvester,
    Finite Elements and Fast Iterative Solvers with
    Applications in Incompressible Fluid Dynamics,
    Oxford, 2005,
    ISBN: 978-0198528678,
    LC: QA911.E39.

  Parameters:

    Input, int N, the number of evaluation points.

    Input, double X[N], Y[N], the coordinates of the points.

    Output, double F[N], G[N], H[N], the right hand sides in the U,
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

void stokes_gnuplot ( char *header, int n, double x[], double y[], double u[], 
  double v[], double s )

/******************************************************************************/
/*
  Purpose:

    STOKES_GNUPLOT writes the Stokes velocity field to files for GNUPLOT.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    24 January 2015

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
  fprintf ( command_unit, "set title 'Stokes velocity field'\n" );
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

void uvp_stokes1 ( int n, double x[], double y[], double u[], double v[], 
  double p[] )

/******************************************************************************/
/*
  Purpose:

    UVP_STOKES1 evaluates the exact Stokes solution #1.

  Discussion:

    The solution is defined over the unit square [0,1]x[0,1].

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    24 January 2015

  Author:

    John Burkardt

  Reference:

    Junping Wang, Yanqiu Wang, Xiu Ye,
    A robust numerical method for Stokes equations based on divergence-free
    H(div) finite element methods,
    SIAM Journal on Scientific Computing,
    Volume 31, Number 4, 2009, pages 2784-2802.

  Parameters:

    Input, int N, the number of evaluation points.

    Input, double X[N], Y[N], the coordinates of the points.

    Output, double U[N], V[N], P[N], the velocity components and
    pressure at each of the points.
*/
{
  int i;

  for ( i = 0; i < n; i++ )
  {

    u[i] = - 2.0 
          * pow ( x[i], 2 ) * pow ( x[i] - 1.0, 2 )
          * y[i] * ( y[i] - 1.0 ) * ( 2.0 * y[i] - 1.0 );

    v[i] =   2.0 
          * x[i] * ( x[i] - 1.0 ) * ( 2.0 * x[i] - 1.0 ) 
          * pow ( y[i], 2 ) * pow ( y[i] - 1.0, 2 );

    p[i] = 0.0;
  }

  return;
}
/******************************************************************************/

void uvp_stokes2 ( int n, double x[], double y[], double u[], double v[], 
  double p[] )

/******************************************************************************/
/*
  Purpose:

    UVP_STOKES2 evaluates the exact Stokes solution #2.

  Discussion:

    The solution is defined over the unit square [0,1]x[0,1].

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    24 January 2015

  Author:

    John Burkardt

  Reference:

    Junping Wang, Yanqiu Wang, Xiu Ye,
    A robust numerical method for Stokes equations based on divergence-free
    H(div) finite element methods,
    SIAM Journal on Scientific Computing,
    Volume 31, Number 4, 2009, pages 2784-2802.

  Parameters:

    Input, int N, the number of evaluation points.

    Input, double X[N], Y[N], the coordinates of the points.

    Output, double U[N], V[N], P[N], the velocity components and
    pressure at each of the points.
*/
{
  int i;
  const double r8_pi = 3.141592653589793;

  for ( i = 0; i < n; i++ )
  {
    u[i] =   2.0 
          * sin ( 2.0 * r8_pi * x[i] ) 
          * cos ( 2.0 * r8_pi * y[i] );


    v[i] = - 2.0 
          * cos ( 2.0 * r8_pi * x[i] ) 
          * sin ( 2.0 * r8_pi * y[i] );

    p[i] = pow ( x[i], 2 ) + pow ( y[i], 2 );
  }

  return;
}
/******************************************************************************/

void uvp_stokes3 ( int n, double x[], double y[], double u[], double v[], 
  double p[] )

/******************************************************************************/
/*
  Purpose:

    UVP_STOKES3 evaluates the exact Stokes solution #3.

  Discussion:

    The solution is defined over the unit square [-1,+1]x[-1,+1].

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 February 2015

  Author:

    John Burkardt

  Reference:

    Howard Elman, Alison Ramage, David Silvester,
    Finite Elements and Fast Iterative Solvers with
    Applications in Incompressible Fluid Dynamics,
    Oxford, 2005,
    ISBN: 978-0198528678,
    LC: QA911.E39.

  Parameters:

    Input, int N, the number of evaluation points.

    Input, double X[N], Y[N], the coordinates of the points.

    Output, double U[N], V[N], P[N], the velocity components and
    pressure at each of the points.
*/
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    u[i] =   20.0 * x[i] * pow ( y[i], 3 );
    v[i] =    5.0 * ( pow ( x[i], 4 )  - pow ( y[i], 4 ) );
    p[i] =   60.0 * pow ( x[i], 2 ) * y[i] - 2.0 * pow ( y[i], 3 ) + 10.0;
  }

  return;
}
