# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <string.h>

# include "disk_rule.h"

/******************************************************************************/

void disk_rule ( int nr, int nt, double w[], double r[], double t[] )

/******************************************************************************/
/*
  Purpose:

    DISK_RULE computes a quadrature rule for the unit disk.

  Discussion:

    The unit disk is the region:

      x * x + y * y <= 1.

    The integral I(f) is then approximated by

      Q(f) = pi * sum ( 1 <= j <= NT ) sum ( 1 <= i <= NR ) 
        W(i) * F ( R(i) * cos(T(j)), R(i) * sin(T(j)) ).

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 March 2014

  Author:

    John Burkardt

  Parameters:

    Input, int NR, the number of points in the radial rule.

    Input, int NT, the number of angles to use.

    Output, double W[NR], the weights for the disk rule.

    Output, double R[NR], T[NT], the (R,Theta) points for the rule.
*/
{
  int ir;
  int it;
  double r8_pi = 3.141592653589793;
  double *wr;
  double *xr;
/*
  Request a Legendre rule for [-1,+1].
*/
  xr = ( double * ) malloc ( nr * sizeof ( double ) );
  wr = ( double * ) malloc ( nr * sizeof ( double ) );

  legendre_ek_compute ( nr, xr, wr );
/*
  Shift the rule to [0,1].
*/
  for ( ir = 0; ir < nr; ir++ )
  {
    xr[ir] = ( xr[ir] + 1.0 ) / 2.0;
    wr[ir] = wr[ir] / 2.0;
  }
/*
  Compute the disk rule.
*/
  for ( it = 0; it < nt; it++ )
  {
    t[it] = 2.0 * r8_pi * ( double ) ( it  ) / ( double ) ( nt );
  }

  for ( ir = 0; ir < nr; ir++ )
  {
    w[ir] = wr[ir] / ( double ) ( nt );
    r[ir] = sqrt ( xr[ir] );
  }

  free ( wr );
  free ( xr );

  return;
}
/******************************************************************************/

double disk01_monomial_integral ( int e[2] )

/******************************************************************************/
/*
  Purpose:

    DISK01_MONOMIAL_INTEGRAL returns monomial integrals in the unit disk in 2D.

  Discussion:

    The integration region is 

      X^2 + Y^2 <= 1.

    The monomial is F(X,Y) = X^E(1) * Y^E(2).

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    03 January 2014

  Author:

    John Burkardt

  Reference:

    Philip Davis, Philip Rabinowitz,
    Methods of Numerical Integration,
    Second Edition,
    Academic Press, 1984, page 263.

  Parameters:

    Input, int E[2], the exponents of X and Y in the 
    monomial.  Each exponent must be nonnegative.

    Output, double DISK01_MONOMIAL_INTEGRAL, the integral.
*/
{
  double arg;
  int i;
  double integral;
  const double r = 1.0;
  const double r8_pi = 3.141592653589793;
  double s;

  if ( e[0] < 0 || e[1] < 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "DISK01_MONOMIAL_INTEGRAL - Fatal error!\n" );
    fprintf ( stderr, "  All exponents must be nonnegative.\n" );
    fprintf ( stderr, "  E[0] = %d\n", e[0] );
    fprintf ( stderr, "  E[1] = %d\n", e[1] );
    exit ( 1 );
  }

  if ( ( e[0] % 2 ) == 1 || ( e[1] % 2 ) == 1 )
  {
    integral = 0.0;
  }
  else
  {
    integral = 2.0;

    for ( i = 0; i < 2; i++ )
    {
      arg = 0.5 * ( double ) ( e[i] + 1 );
      integral = integral * r8_gamma ( arg );
    }
    arg = 0.5 * ( double ) ( e[0] + e[1] + 2 );
    integral = integral / r8_gamma ( arg );
  }
/*
  Adjust the surface integral to get the volume integral.
*/
  s = e[0] + e[1] + 2;
  integral = integral * pow ( r, s ) / ( double ) ( s );

  return integral;
}
/******************************************************************************/

void imtqlx ( int n, double d[], double e[], double z[] )

/******************************************************************************/
/*
  Purpose:

    IMTQLX diagonalizes a symmetric tridiagonal matrix.

  Discussion:

    This routine is a slightly modified version of the EISPACK routine to 
    perform the implicit QL algorithm on a symmetric tridiagonal matrix. 

    The authors thank the authors of EISPACK for permission to use this
    routine. 

    It has been modified to produce the product Q' * Z, where Z is an input 
    vector and Q is the orthogonal matrix diagonalizing the input matrix.  
    The changes consist (essentially) of applying the orthogonal transformations
    directly to Z as they are generated.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 January 2010

  Author:

    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
    C version by John Burkardt.

  Reference:

    Sylvan Elhay, Jaroslav Kautsky,
    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
    Interpolatory Quadrature,
    ACM Transactions on Mathematical Software,
    Volume 13, Number 4, December 1987, pages 399-415.

    Roger Martin, James Wilkinson,
    The Implicit QL Algorithm,
    Numerische Mathematik,
    Volume 12, Number 5, December 1968, pages 377-383.

  Parameters:

    Input, int N, the order of the matrix.

    Input/output, double D(N), the diagonal entries of the matrix.
    On output, the information in D has been overwritten.

    Input/output, double E(N), the subdiagonal entries of the 
    matrix, in entries E(1) through E(N-1).  On output, the information in
    E has been overwritten.

    Input/output, double Z(N).  On input, a vector.  On output,
    the value of Q' * Z, where Q is the matrix that diagonalizes the
    input symmetric tridiagonal matrix.
*/
{
  double b;
  double c;
  double f;
  double g;
  int i;
  int ii;
  int itn = 30;
  int j;
  int k;
  int l;
  int m;
  int mml;
  double p;
  double prec;
  double r;
  double s;

  prec = r8_epsilon ( );

  if ( n == 1 )
  {
    return;
  }

  e[n-1] = 0.0;

  for ( l = 1; l <= n; l++ )
  {
    j = 0;
    for ( ; ; )
    {
      for ( m = l; m <= n; m++ )
      {
        if ( m == n )
        {
          break;
        }

        if ( r8_abs ( e[m-1] ) <= prec * ( r8_abs ( d[m-1] ) + r8_abs ( d[m] ) ) )
        {
          break;
        }
      }
      p = d[l-1];
      if ( m == l )
      {
        break;
      }
      if ( itn <= j )
      {
        printf ( "\n" );
        printf ( "IMTQLX - Fatal error!\n" );
        printf ( "  Iteration limit exceeded\n" );
        exit ( 1 );
      }
      j = j + 1;
      g = ( d[l] - p ) / ( 2.0 * e[l-1] );
      r =  sqrt ( g * g + 1.0 );
      g = d[m-1] - p + e[l-1] / ( g + r8_abs ( r ) * r8_sign ( g ) );
      s = 1.0;
      c = 1.0;
      p = 0.0;
      mml = m - l;

      for ( ii = 1; ii <= mml; ii++ )
      {
        i = m - ii;
        f = s * e[i-1];
        b = c * e[i-1];

        if ( r8_abs ( g ) <= r8_abs ( f ) )
        {
          c = g / f;
          r =  sqrt ( c * c + 1.0 );
          e[i] = f * r;
          s = 1.0 / r;
          c = c * s;
        }
        else
        {
          s = f / g;
          r =  sqrt ( s * s + 1.0 );
          e[i] = g * r;
          c = 1.0 / r;
          s = s * c;
        }
        g = d[i] - p;
        r = ( d[i-1] - g ) * s + 2.0 * c * b;
        p = s * r;
        d[i] = g + p;
        g = c * r - b;
        f = z[i];
        z[i] = s * z[i-1] + c * f;
        z[i-1] = c * z[i-1] - s * f;
      }
      d[l-1] = d[l-1] - p;
      e[l-1] = g;
      e[m-1] = 0.0;
    }
  }
/*
  Sorting.
*/
  for ( ii = 2; ii <= m; ii++ )
  {
    i = ii - 1;
    k = i;
    p = d[i-1];

    for ( j = ii; j <= n; j++ )
    {
      if ( d[j-1] < p )
      {
         k = j;
         p = d[j-1];
      }
    }

    if ( k != i )
    {
      d[k-1] = d[i-1];
      d[i-1] = p;
      p = z[i-1];
      z[i-1] = z[k-1];
      z[k-1] = p;
    }
  }
  return;
}
/******************************************************************************/

void legendre_ek_compute ( int n, double x[], double w[] )

/******************************************************************************/
/*
  Purpose:

    LEGENDRE_EK_COMPUTE: Legendre quadrature rule by the Elhay-Kautsky method.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 March 2014

  Author:

    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
    C version by John Burkardt.

  Reference:

    Sylvan Elhay, Jaroslav Kautsky,
    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
    Interpolatory Quadrature,
    ACM Transactions on Mathematical Software,
    Volume 13, Number 4, December 1987, pages 399-415.

  Parameters:

    Input, int N, the order.

    Output, double X[N], the abscissas.

    Output, double W[N], the weights.
*/
{
  double *bj;
  int i;
  double ip1;
  double zemu;
/*
  Define the zero-th moment.
*/
  zemu = 2.0;
/*
  Define the Jacobi matrix.
*/
  bj = ( double * ) malloc ( n * sizeof ( double ) );

  for ( i = 0; i < n; i++ )
  {
    ip1 = ( double ) ( i + 1 );
    bj[i] = sqrt ( ip1 * ip1 / ( 4.0 * ip1 * ip1 - 1.0 ) );
  }

  for ( i = 0; i < n; i++ )
  {
    x[i] = 0.0;
  }

  w[0] = sqrt ( zemu );
  for ( i = 1; i < n; i++ )
  {
    w[i] = 0.0;
  }
/*
  Diagonalize the Jacobi matrix.
*/
  imtqlx ( n, x, bj, w );

  for ( i = 0; i < n; i++ )
  {
    w[i] = w[i] * w[i];
  }

  free ( bj );

  return;
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
    value = + x;
  }
  else
  {
    value = - x;
  }
  return value;
}
/******************************************************************************/

double r8_epsilon ( )

/******************************************************************************/
/*
  Purpose:

    R8_EPSILON returns the R8 round off unit.

  Discussion:

    R8_EPSILON is a number R which is a power of 2 with the property that,
    to the precision of the computer's arithmetic,
      1 < 1 + R
    but
      1 = ( 1 + R / 2 )

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 September 2012

  Author:

    John Burkardt

  Parameters:

    Output, double R8_EPSILON, the R8 round-off unit.
*/
{
  const double value = 2.220446049250313E-016;

  return value;
}
/******************************************************************************/

double r8_gamma ( double x )

/******************************************************************************/
/*
  Purpose:

    R8_GAMMA evaluates Gamma(X) for a real argument.

  Discussion:

    The C math library includes the GAMMA ( X ) function which should generally
    be used instead of this function.

    This routine calculates the gamma function for a real argument X.

    Computation is based on an algorithm outlined in reference 1.
    The program uses rational functions that approximate the gamma
    function to at least 20 significant decimal digits.  Coefficients
    for the approximation over the interval (1,2) are unpublished.
    Those for the approximation for 12 <= X are from reference 2.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 January 2010

  Author:

    Original FORTRAN77 version by William Cody, Laura Stoltz.
    C version by John Burkardt.

  Reference:

    William Cody,
    An Overview of Software Development for Special Functions,
    in Numerical Analysis Dundee, 1975,
    edited by GA Watson,
    Lecture Notes in Mathematics 506,
    Springer, 1976.

    John Hart, Ward Cheney, Charles Lawson, Hans Maehly,
    Charles Mesztenyi, John Rice, Henry Thatcher,
    Christoph Witzgall,
    Computer Approximations,
    Wiley, 1968,
    LC: QA297.C64.

  Parameters:

    Input, double X, the argument of the function.

    Output, double R8_GAMMA, the value of the function.
*/
{
  double c[7] = {
   -1.910444077728E-03, 
    8.4171387781295E-04, 
   -5.952379913043012E-04, 
    7.93650793500350248E-04, 
   -2.777777777777681622553E-03, 
    8.333333333333333331554247E-02, 
    5.7083835261E-03 };
  double eps = 2.22E-16;
  double fact;
  int i;
  int n;
  double p[8] = {
  -1.71618513886549492533811E+00,
   2.47656508055759199108314E+01, 
  -3.79804256470945635097577E+02,
   6.29331155312818442661052E+02, 
   8.66966202790413211295064E+02,
  -3.14512729688483675254357E+04, 
  -3.61444134186911729807069E+04,
   6.64561438202405440627855E+04 };
  int parity;
  const double pi = 3.1415926535897932384626434;
  double q[8] = {
  -3.08402300119738975254353E+01,
   3.15350626979604161529144E+02, 
  -1.01515636749021914166146E+03,
  -3.10777167157231109440444E+03, 
   2.25381184209801510330112E+04,
   4.75584627752788110767815E+03, 
  -1.34659959864969306392456E+05,
  -1.15132259675553483497211E+05 };
  double res;
  const double sqrtpi = 0.9189385332046727417803297;
  double sum;
  double value;
  double xbig = 171.624;
  double xden;
  double xinf = 1.79E+308;
  double xminin = 2.23E-308;
  double xnum;
  double y;
  double y1;
  double ysq;
  double z;

  parity = 0;
  fact = 1.0;
  n = 0;
  y = x;
/*
  Argument is negative.
*/
  if ( y <= 0.0 )
  {
    y = - x;
    y1 = ( double ) ( int ) ( y );
    res = y - y1;

    if ( res != 0.0 )
    {
      if ( y1 != ( double ) ( int ) ( y1 * 0.5 ) * 2.0 )
      {
        parity = 1;
      }

      fact = - pi / sin ( pi * res );
      y = y + 1.0;
    }
    else
    {
      res = xinf;
      value = res;
      return value;
    }
  }
/*
  Argument is positive.
*/
  if ( y < eps )
  {
/*
  Argument < EPS.
*/
    if ( xminin <= y )
    {
      res = 1.0 / y;
    }
    else
    {
      res = xinf;
      value = res;
      return value;
    }
  }
  else if ( y < 12.0 )
  {
    y1 = y;
/*
  0.0 < argument < 1.0.
*/
    if ( y < 1.0 )
    {
      z = y;
      y = y + 1.0;
    }
/*
  1.0 < argument < 12.0.
  Reduce argument if necessary.
*/
    else
    {
      n = ( int ) ( y ) - 1;
      y = y - ( double ) ( n );
      z = y - 1.0;
    }
/*
  Evaluate approximation for 1.0 < argument < 2.0.
*/
    xnum = 0.0;
    xden = 1.0;
    for ( i = 0; i < 8; i++ )
    {
      xnum = ( xnum + p[i] ) * z;
      xden = xden * z + q[i];
    }
    res = xnum / xden + 1.0;
/*
  Adjust result for case  0.0 < argument < 1.0.
*/
    if ( y1 < y )
    {
      res = res / y1;
    }
/*
  Adjust result for case 2.0 < argument < 12.0.
*/
    else if ( y < y1 )
    {
      for ( i = 1; i <= n; i++ )
      {
        res = res * y;
        y = y + 1.0;
      }
    }
  }
  else
  {
/*
  Evaluate for 12.0 <= argument.
*/
    if ( y <= xbig )
    {
      ysq = y * y;
      sum = c[6];
      for ( i = 0; i < 6; i++ )
      {
        sum = sum / ysq + c[i];
      }
      sum = sum / y - y + sqrtpi;
      sum = sum + ( y - 0.5 ) * log ( y );
      res = exp ( sum );
    }
    else
    {
      res = xinf;
      value = res;
      return value;
    }
  }
/*
  Final adjustments and return.
*/
  if ( parity )
  {
    res = - res;
  }

  if ( fact != 1.0 )
  {
    res = fact / res;
  }

  value = res;

  return value;
}
/******************************************************************************/

double r8_sign ( double x )

/******************************************************************************/
/*
  Purpose:

    R8_SIGN returns the sign of an R8.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 May 2006

  Author:

    John Burkardt

  Parameters:

    Input, double X, the number whose sign is desired.

    Output, double R8_SIGN, the sign of X.
*/
{
  double value;

  if ( x < 0.0 )
  {
    value = - 1.0;
  }
  else
  {
    value = + 1.0;
  }
  return value;
}
/******************************************************************************/

void r8vec_print ( int n, double a[], char *title )

/******************************************************************************/
/*
  Purpose:

    R8VEC_PRINT prints an R8VEC.

  Discussion:

    An R8VEC is a vector of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 April 2009

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of components of the vector.

    Input, double A[N], the vector to be printed.

    Input, char *TITLE, a title.
*/
{
  int i;

  fprintf ( stdout, "\n" );
  fprintf ( stdout, "%s\n", title );
  fprintf ( stdout, "\n" );
  for ( i = 0; i < n; i++ )
  {
    fprintf ( stdout, "  %8d: %14g\n", i, a[i] );
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
