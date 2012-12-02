# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>

int main ( int argc, char **argv );
int i4_ceiling ( double x );
int i4_min ( int i1, int i2 );
double potential ( double a, double b, double c, double x, double y, double z );
double r8_abs ( double x );
double r8_uniform_01 ( int *seed );
void timestamp ( void );

/******************************************************************************/

int main ( int arc, char **argv )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for FEYNMAN_KAC_3D.

  Discussion:

    This program is derived from section 2.5, exercise 2.2 of Petersen and Arbenz.

    The problem is to determine the solution U(X,Y,Z) of the following 
    partial differential equation:

      (1/2) Laplacian U - V(X,Y,Z) * U = 0,

    inside the elliptic domain D:
 
      D = { (X,Y,Z) | (X/A)^2+(Y/B)^2+(Z/C)^2 <= 1 }
   
    with the boundary condition U(boundary(D)) = 1.

    The V(X,Y,Z) is the potential function:

      V = 2 * ( (X/A^2)^2 + (Y/B^2)^2 + (Z/C^2)^2 ) + 1/A^2 + 1/B^2 + 1/C^2.

    The analytic solution of this problem is already known:

      U(X,Y,Z) = exp ( (X/A)^2 + (Y/B)^2 + (Z/C)^2 - 1 ).

    Our method is via the Feynman-Kac Formula.

    The idea is to start from any (x,y,z) in D, and
    compute (x+Wx(t),y+Wy(t),z+Wz(t)) where 3-D Brownian motion
    (Wx,Wy,Wz) is updated each step by sqrt(h)*(z1,z2,z3),
    each z1,z2,z3 are independent approximately Gaussian 
    random variables with zero mean and variance 1. 

    Each (x1(t),x2(t),x3(t)) is advanced until (x1,x2,x3) exits 
    the domain D.  

    Upon its first exit from D, the sample path (x1,x2,x3) is stopped and a 
    new sample path at (x,y,z) is started until N such paths are completed.
 
    The Feynman-Kac formula gives the solution here as

      U(X,Y,Z) = (1/N) sum(1 <= I <= N) Y(tau_i),

    where

      Y(tau) = exp( -int(s=0..tau) v(x1(s),x2(s),x3(s)) ds),

    and tau = first exit time for path (x1,x2,x3). 

    The integration procedure is a second order weak accurate method:

      X(t+h)  = [ x1(t) + sqrt ( h ) * z1 ]
                [ x2(t) + sqrt ( h ) * z2 ]
                [ x3(t) + sqrt ( h ) * z3 ]

    Here Z1, Z2, and Z3 are approximately normal univariate Gaussians. 

    An Euler predictor approximates Y at the end of the step

      Y_e     = (1 - h*v(X(t)) * Y(t), 

    A trapezoidal rule completes the step:

      Y(t+h)  = Y(t) - (h/2)*[v(X(t+h))*Y_e + v(X(t))*Y(t)].

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    29 May 2012

  Author:

    Original C version by Wesley Petersen.
    C++ version by John Burkardt.

  Reference:

    Peter Arbenz, Wesley Petersen,
    Introduction to Parallel Computing:
    A Practical Guide with Examples in C,
    Oxford, 2004,
    ISBN: 0-19-851577-4,
    LC: QA76.59.P47.
*/
{
  double a = 3.0;
  double b = 2.0;
  double c = 1.0;
  double chk;
  int dim = 3;
  double dx;
  double dy;
  double dz;
  double err;
  double h = 0.001;
  int i;
  int j;
  int k;
  int N = 10000;
  int n_inside;
  int ni;
  int nj;
  int nk;
  double rth;
  int seed = 123456789;
  int steps;
  int steps_ave;
  int trial;
  double us;
  double ut;
  double vh;
  double vs;
  double x;
  double x1;
  double x2;
  double x3;
  double y;
  double w;
  double w_exact;
  double we;
  double wt;
  double z;
  double z1;
  double z2;
  double z3;

  printf ( "\n" );
  timestamp ( );

  printf ( "\n" );
  printf ( "FEYNMAN_KAC_3D:\n" );
  printf ( "  C version.\n" );
  printf ( "\n" );
  printf ( "  Program parameters:\n" );
  printf ( "\n" );
  printf ( "  The calculation takes place inside a 3D ellipsoid.\n" );
  printf ( "  A rectangular grid of points will be defined.\n" );
  printf ( "  The solution will be estimated for those grid points\n" );
  printf ( "  that lie inside the ellipsoid.\n" );
  printf ( "\n" );
  printf ( "  Each solution will be estimated by computing %d trajectories\n", N );
  printf ( "  from the point to the boundary.\n" );
  printf ( "\n" );
  printf ( "    (X/A)^2 + (Y/B)^2 + (Z/C)^2 = 1\n" );
  printf ( "\n" );
  printf ( "  The ellipsoid parameters A, B, C are set to:\n" );
  printf ( "\n" );
  printf ( "    A = %f\n", a );
  printf ( "    B = %f\n", b );
  printf ( "    C = %f\n", c );
  printf ( "  Stepsize H = %6.4f\n", h );
/*
  RTH is the scaled stepsize.
*/
  rth = sqrt ( ( double ) dim * h );
/*
  Choose the spacing so we have about 10 points in the shortest direction.
*/
  if ( a == i4_min ( i4_min ( a, b ), c ) )
  {
    ni = 11;
    nj = 1 + i4_ceiling ( b / a ) * ( ni - 1 );
    nk = 1 + i4_ceiling ( c / a ) * ( ni - 1 );
  }
  else if ( b == i4_min ( i4_min ( a, b ), c ) )
  {
    nj = 11;
    ni = 1 + i4_ceiling ( a / b ) * ( nj - 1 );
    nk = 1 + i4_ceiling ( c / b ) * ( nj - 1 );
  }
  else
  {
    nk = 11;
    ni = 1 + i4_ceiling ( a / c ) * ( nk - 1 );
    nj = 1 + i4_ceiling ( b / c ) * ( nk - 1 );
  }

  printf ( "\n" );
  printf ( "  X coordinate marked by %d points\n", ni );
  printf ( "  Y coordinate marked by %d points\n", nj );
  printf ( "  Z coordinate marked by %d points\n", nk );
/*
  Loop over grid points.
*/
  err = 0.0;
  n_inside = 0;

  printf ( "\n" );
  printf ( "     X        Y        Z       W approx     W exact    Error    Ave Steps\n" );
  printf ( "\n" );

  for ( i = 1; i <= ni; i++ )
  {
    x = ( ( double ) ( ni - i     ) * ( - a )
        + ( double ) (      i - 1 ) *     a )
        / ( double ) ( ni     - 1 );

    for ( j = 1; j <= nj; j++ )
    {
      y = ( ( double ) ( nj - j     ) * ( - b )
          + ( double ) (      j - 1 ) *     b ) 
          / ( double ) ( nj     - 1 );

      for ( k = 1; k <= nk; k++ )
      {
        z = ( ( double ) ( nk - k     ) * ( - c )
            + ( double ) (      k - 1 ) *     c ) 
            / ( double ) ( nk     - 1 );

        chk = pow ( x / a, 2 ) + pow ( y / b, 2 ) + pow ( z / c, 2 );

        if ( 1.0 < chk )
        {
          w_exact = 1.0;
          wt = 1.0;
          steps_ave = 0;
          printf ( "  %7.4f  %7.4f  %7.4f  %10.4e  %10.4e  %10.4e  %8d\n",
            x, y, z, wt, w_exact, r8_abs ( w_exact - wt ), steps_ave );
          continue;
        }

        n_inside = n_inside + 1;
/*
  Compute the exact solution at this point (x,y,z).
*/
        w_exact = exp ( pow ( x / a, 2 )
                      + pow ( y / b, 2 )
                      + pow ( z / c, 2 ) - 1.0 );
/*
  Now try to estimate the solution at this point.
*/
        wt = 0.0;
        steps = 0;
/* 
  Do N paths 
*/
        for ( trial = 0; trial < N; trial++ )
        {
          x1  = x;
          x2  = y;
          x3  = z;
/* 
  W = exp(-int(s=0..t) v(X)ds) 
*/
          w = 1.0;  
/*
  CHK is < 1.0 while the point is inside the ellipsoid.
*/
          chk = 0.0;
          while ( chk < 1.0 )
          {
/*
  Determine DX, DY, DZ.
*/
            ut = r8_uniform_01 ( &seed );
            if ( ut < 1.0 / 3.0 )
            {
              us = r8_uniform_01 ( &seed ) - 0.5;
              if ( us < 0.0)
              {
                dx = - rth;
              } 
              else
              {
                dx = rth;
              }
            } 
            else
            {
              dx = 0.0;
            }

            ut = r8_uniform_01 ( &seed );
            if ( ut < 1.0 / 3.0 )
            {
              us = r8_uniform_01 ( &seed ) - 0.5;
              if ( us < 0.0 )
              { 
                dy = - rth;
              }
              else
              {
                dy = rth;
              }
            }
            else
            {
              dy = 0.0;
            }

            ut = r8_uniform_01 ( &seed );
            if ( ut < 1.0 / 3.0 )
            {
              us = r8_uniform_01 ( &seed ) - 0.5;
              if ( us < 0.0 )
              {
                dz = - rth;
              }
              else
              {
                dz = rth;
              }
            } 
            else
            {
              dz = 0.0;
            }

            vs = potential ( a, b, c, x1, x2, x3 );
/*
  Move to the new point.
*/
            x1 = x1 + dx;
            x2 = x2 + dy;
            x3 = x3 + dz;

            steps = steps + 1;

            vh = potential ( a, b, c, x1, x2, x3 );

            we = ( 1.0 - h * vs ) * w;
            w = w - 0.5 * h * ( vh * we + vs * w ); 

            chk = pow ( x1 / a, 2 ) 
                + pow ( x2 / b, 2 )
                + pow ( x3 / c, 2 );
          }
          wt = wt + w;
        }
/*
  WT is the average of the sum of the different trials.
*/
        wt = wt / ( double ) ( N ); 
        steps_ave = steps / ( double ) ( N );
/*
  Add error in WT to the running L2 error in the solution.
*/
        err = err + pow ( w_exact - wt, 2 );

        printf ( "  %7.4f  %7.4f  %7.4f  %10.4e  %10.4e  %10.4e  %8d\n",
          x, y, z, wt, w_exact, r8_abs ( w_exact - wt ), steps_ave );
      }
    }
  }
/*
  Compute the RMS error for all the points.
*/
  err = sqrt ( err / ( double ) ( n_inside ) );

  printf ( "\n" );
  printf ( "  RMS absolute error in solution = %e\n", err );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "FEYNMAN_KAC_3D:\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

int i4_ceiling ( double x )

/******************************************************************************/
/*
  Purpose:

    I4_CEILING rounds an R8 up to the nearest I4.

  Discussion:

    The "ceiling" of X is the value of X rounded towards plus infinity.

  Example:

    X        I4_CEILING(X)

   -1.1      -1
   -1.0      -1
   -0.9       0
   -0.1       0
    0.0       0
    0.1       1
    0.9       1
    1.0       1
    1.1       2
    2.9       3
    3.0       3
    3.14159   4

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    10 November 2011

  Author:

    John Burkardt

  Parameters:

    Input, double X, the number whose ceiling is desired.

    Output, int I4_CEILING, the ceiling of X.
*/
{
  int value;

  value = ( int ) x;

  if ( value < x )
  {
    value = value + 1;
  }

  return value;
}
/******************************************************************************/

int i4_min ( int i1, int i2 )

/******************************************************************************/
/*
  Purpose:

    I4_MIN returns the smaller of two I4's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    29 August 2006

  Author:

    John Burkardt

  Parameters:

    Input, int I1, I2, two integers to be compared.

    Output, int I4_MIN, the smaller of I1 and I2.
*/
{
  int value;

  if ( i1 < i2 )
  {
    value = i1;
  }
  else
  {
    value = i2;
  }
  return value;
}
/******************************************************************************/

double potential ( double a, double b, double c, double x, double y, double z )

/******************************************************************************/
/*
  Purpose:

    POTENTIAL evaluates the potential function V(X,Y,Z).

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    19 February 2008

  Author:

    John Burkardt

  Parameters:

    Input, double A, B, C, the parameters that define the ellipse.

    Input, double X, Y, Z, the coordinates of the point.

    Output, double POTENTIAL, the value of the potential function at (X,Y,Z).
*/
{
  double value;

  value = 2.0 * ( pow ( x / a / a, 2 )
                + pow ( y / b / b, 2 )
                + pow ( z / c / c, 2 ) )
              + 1.0 / a / a
              + 1.0 / b / b 
              + 1.0 / c / c;

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
