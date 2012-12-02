# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>

int main ( int argc, char **argv );
double potential ( double a, double x );
double r8_abs ( double x );
double r8_uniform_01 ( int *seed );
void timestamp ( void );

/******************************************************************************/

int main ( int argc, char **argv )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for FEYNMAN_KAC_1D.

  Discussion:

    This program is derived from section 2.5, exercise 2.2 of Petersen 
    and Arbenz.

    The problem is to determine the solution U(X) of the following 
    partial differential equation:

      (1/2) Laplacian U - V(X) * U = 0,

    inside the domain D:
 
      D = { X | (X/A)^2 <= 1 }
   
    with the boundary condition U(boundary(D)) = 1.

    V(X) is the potential function:

      V = 2 * ( (X/A^2)^2 + 1/A^2.

    The analytic solution of this problem is already known:

      U(X) = exp ( (X/A)^2 - 1 ).

    Our method is via the Feynman-Kac Formula.

    The idea is to start from any x in D, and
    compute x+Wx(t) where 1D Brownian motion
    Wx is updated each step by sqrt(h)*z,
    each z is an independent approximately Gaussian 
    random variable with zero mean and variance 1. 

    Each x1(t) is advanced until x1(t) exits the domain D.  

    Upon its first exit from D, the sample path x1 is stopped and a 
    new sample path at x is started until N such paths are completed.
 
    The Feynman-Kac formula gives the solution here as

      U(X) = (1/N) sum(1 <= I <= N) Y(tau_i),

    where

      Y(tau) = exp( -int(s=0..tau) v(x1(s)) ds),

    and tau = first exit time for path x1. 

    The integration procedure is a second order weak accurate method:

      X(t+h)  = x1(t) + sqrt ( h ) * z

    Here Z is an approximately normal univariate Gaussian. 

    An Euler predictor approximates Y at the end of the step

      Y_e     = (1 - h*v(X(t)) * Y(t), 

    A trapezoidal rule completes the step:

      Y(t+h)  = Y(t) - (h/2)*[v(X(t+h))*Y_e + v(X(t))*Y(t)].

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    29 May 2012

  Author:

    Original C 3D version by Wesley Petersen.
    C 1D version by John Burkardt.

  Reference:

    Peter Arbenz, Wesley Petersen,
    Introduction to Parallel Computing:
    A Practical Guide with Examples in C,
    Oxford, 2004,
    ISBN: 0-19-851577-4,
    LC: QA76.59.P47.
*/
{
  double a = 2.0;
  double chk;
  double dx;
  double err;
  double h = 0.0001;
  int i;
  int it;
  int j;
  int k;
  int n1;
  int n = 10000;
  int n_int;
  int ni;
  double rth;
  int seed = 123456789;
  int steps;
  int steps_ave;
  double sum;
  double test;
  double us;
  double vh;
  double vs;
  double x;
  double x1;
  double w;
  double w_exact;
  double we;
  double wt;

  timestamp ( );

  printf ( "\n" );
  printf ( "FEYNMAN_KAC_1D:\n" );
  printf ( "  C version.\n" );
  printf ( "\n" );
  printf ( "  Program parameters:\n" );
  printf ( "\n" );
  printf ( "  The calculation takes place inside an interval.\n" );
  printf ( "  The solution will be estimated at points\n" );
  printf ( "  on a regular spaced grid within the interval.\n" );
  printf ( "  Each solution will be estimated by computing %d trajectories\n", n );
  printf ( "  from the point to the boundary.\n" );
  printf ( "\n" );
  printf ( "    (X/A)^2 = 1\n" );
  printf ( "\n" );
  printf ( "  The interval parameter A is:\n" );
  printf ( "\n" );
  printf ( "    A = %g\n", a );
  printf ( "\n" );
  printf ( "  Path stepsize H = %g\n", h );
/*
  Choose the spacing so we have about ni points on or in the interval.
*/
  ni = 21;

  printf ( "\n" );
  printf ( "  X coordinate discretized by %d points\n", ni + 2 );
/*
  RTH is the scaled stepsize.
*/
  rth = sqrt ( h );

  err = 0.0;
/*
  Loop over the points.
*/
  printf ( "\n" );
  printf ( "     I     K       X           W exact" );
  printf ( "      W Approx        Error      Ave Steps  Test\n" );
  printf ( "\n" );

  k = 0;
  n_int = 0;

  for ( i = 0; i <= ni + 1; i++ )
  {
    x = ( ( double ) ( ni - i     ) * ( - a ) 
        + ( double ) (      i - 1 ) *     a ) 
        / ( double ) ( ni     - 1 );

    k = k + 1;

    test = a * a - x * x;

    if ( test < 0.0 )
    {
      w_exact = 1.0;
      wt = 1.0;
      steps_ave = 0;
      printf ( "  %4d  %4d  %12g  %12g  %12g  %12g  %8d  %8g\n",
        i, k, x, w_exact, wt, r8_abs ( w_exact - wt ), steps_ave, test );
      continue;
    }

    n_int = n_int + 1;
/*
  Compute the exact solution at this point (x,y,z).
*/
    w_exact = exp ( pow ( x / a, 2 ) - 1.0 );
/*
  Now try to estimate the solution at this point.
*/
    wt = 0.0;
    steps = 0;

    for ( it = 1; it <= n; it++ )
    {

      x1 = x;
/* 
  W = exp(-int(s=0..t) v(X)ds) 
*/
      w = 1.0;
/*
  CHK is < 1.0 while the point is inside the interval.
*/
      chk = 0.0;

      while ( chk < 1.0 )
      {
/*
  Determine DX.
*/
        us = r8_uniform_01 ( &seed ) - 0.5;
        if ( us < 0.0 )
        {
          dx = - rth;
        }
        else
        {
          dx = + rth;
        }

        vs = potential ( a, x1 );
/*
  Move to the new point.
*/
        x1 = x1 + dx;

        steps = steps + 1;

        vh = potential ( a, x1 );

        we = ( 1.0 - h * vs ) * w;
        w = w - 0.5 * h * ( vh * we + vs * w );

        chk = pow ( x1 / a, 2 );
      }
      wt = wt + w;
    }
/*
   WT is the average of the sum of the different trials.
*/
    wt = wt / ( double ) ( n );
    steps_ave = steps / n;
/*
  Add error in WT to the running L2 error in the solution.
*/
    err = err + pow ( w_exact - wt, 2 );

    printf ( "  %4d  %4d  %12g  %12g  %12g  %12g  %8d  %8g\n",
      i, k, x, w_exact, wt, r8_abs ( w_exact - wt ), steps_ave, test );
  }
/*
  Compute the RMS error for all the points.
*/
  err = sqrt ( err / ( double ) ( n_int ) );

  printf ( "\n" );
  printf ( "  RMS absolute error in solution = %g\n", err );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "FEYNMAN_KAC_1D:\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );

  timestamp ( );

  return;
}
/******************************************************************************/

double potential ( double a, double x )

/******************************************************************************/
/*
  Purpose:

    POTENTIAL evaluates the potential function V(X).

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    19 February 2008

  Author:

    John Burkardt

  Parameters:

    Input, double A, the parameter that defines the region.

    Input, double X, the coordinate of the point.

    Output, double POTENTIAL, the value of the potential function.
*/
{
  double value;

  value = 2.0 * pow ( x / a / a, 2 ) + 1.0 / a / a;

  return value;
}
/******************************************************************************/

double r8_abs ( double x )

/******************************************************************************/
/*
  Purpose:

    R8_ABS returns the absolute value of an R8.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    07 May 2006

  Author:

    John Burkardt

  Parameters:

    Input, double X, the quantity whose absolute value is desired.

    Output, double R8_ABS, the absolute value of X.
*/
{
  double value;

  if ( 0.0 <= x )
  {
    value = x;
  } 
  else
  {
    value = - x;
  }
  return value;
}
/******************************************************************************/

double r8_uniform_01 ( int *seed )

/******************************************************************************/
/*
  Purpose:

    R8_UNIFORM_01 returns a unit pseudorandom R8.

  Discussion:

    This routine implements the recursion

      seed = 16807 * seed mod ( 2**31 - 1 )
      r8_uniform_01 = seed / ( 2**31 - 1 )

    The integer arithmetic never requires more than 32 bits,
    including a sign bit.

    If the initial seed is 12345, then the first three computations are

      Input     Output      R8_UNIFORM_01
      SEED      SEED

         12345   207482415  0.096616
     207482415  1790989824  0.833995
    1790989824  2035175616  0.947702

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 August 2004

  Author:

    John Burkardt

  Reference:

    Paul Bratley, Bennett Fox, Linus Schrage,
    A Guide to Simulation,
    Springer Verlag, pages 201-202, 1983.

    Pierre L'Ecuyer,
    Random Number Generation,
    in Handbook of Simulation
    edited by Jerry Banks,
    Wiley Interscience, page 95, 1998.

    Bennett Fox,
    Algorithm 647:
    Implementation and Relative Efficiency of Quasirandom
    Sequence Generators,
    ACM Transactions on Mathematical Software,
    Volume 12, Number 4, pages 362-376, 1986.

    Peter Lewis, Allen Goodman, James Miller,
    A Pseudo-Random Number Generator for the System/360,
    IBM Systems Journal,
    Volume 8, pages 136-143, 1969.

  Parameters:

    Input/output, int *SEED, the "seed" value.  Normally, this
    value should not be 0.  On output, SEED has been updated.

    Output, double R8_UNIFORM_01, a new pseudorandom variate, strictly between
    0 and 1.
*/
{
  int k;
  double r;

  k = *seed / 127773;

  *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

  if ( *seed < 0 )
  {
    *seed = *seed + 2147483647;
  }
/*
  Although SEED can be represented exactly as a 32 bit integer,
  it generally cannot be represented exactly as a 32 bit real number!
*/
  r = ( double ) ( *seed ) * 4.656612875E-10;

  return r;
}
/******************************************************************************/

void timestamp ( void )

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

  printf ( "%s\n", time_buffer );

  return;
# undef TIME_SIZE
}
