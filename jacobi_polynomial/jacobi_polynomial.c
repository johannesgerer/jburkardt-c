# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <string.h>

# include "jacobi_polynomial.h"

/******************************************************************************/

int i4_max ( int i1, int i2 )

/******************************************************************************/
/*
  Purpose:

    I4_MAX returns the maximum of two I4's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    29 August 2006

  Author:

    John Burkardt

  Parameters:

    Input, int I1, I2, are two integers to be compared.

    Output, int I4_MAX, the larger of I1 and I2.
*/
{
  int value;

  if ( i2 < i1 )
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
    The changes consist (essentialy) of applying the orthogonal transformations
    directly to Z as they are generated.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 January 2010

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

double j_double_product_integral ( int i, int j, double a, double b )

/******************************************************************************/
/*
  Purpose:

    J_DOUBLE_PRODUCT_INTEGRAL: integral of J(i,x)*J(j,x)*(1-x)^a*(1+x)^b.

  Discussion:

    VALUE = integral ( -1 <= x <= +1 ) J(i,x)*J(j,x)*(1-x)^a*(1+x)^b dx

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 August 2013

  Author:

    John Burkardt

  Parameters:

    Input, int I, J, the polynomial indices.

    Input, double A, B, the parameters.
    -1 < A, B.

    Output, double VALUE, the value of the integral.
*/
{
  double i_r8;
  double value;

  if ( i != j )
  {
    value = 0.0;
  }
  else
  {
    i_r8 = ( double ) ( i );

    value = pow ( 2, a + b + 1.0 )
      / ( 2.0 * i_r8 + a + b + 1.0 )
      * tgamma ( i_r8 + a + 1.0 )
      * tgamma ( i_r8 + b + 1.0 )
      / r8_factorial ( i )
     / tgamma ( i_r8 + a + b + 1.0 );
  }
  return value;
}
/******************************************************************************/

double j_integral ( int n )

/******************************************************************************/
/*
  Purpose:

    J_INTEGRAL evaluates a monomial integral associated with J(n,a,b,x).

  Discussion:

    The integral:

      integral ( -1 <= x < +1 ) x^n dx

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 August 2013

  Author:

    John Burkardt

  Parameters:

    Input, int N, the exponent.
    0 <= N.

    Output, double J_INTEGRAL, the value of the integral.
*/
{
  double value;

  if ( ( n % 2 ) == 1 )
  {
    value = 0.0;
  }
  else
  {
    value = 2.0 / ( double ) ( n + 1 );
  }

  return value;
}
/******************************************************************************/

double *j_polynomial ( int m, int n, double alpha, double beta, double x[] )

/******************************************************************************/
/*
  Purpose:

    JACOBI_POLY evaluates the Jacobi polynomial J(n,a,b,x).

  Differential equation:

    (1-X*X) Y'' + (BETA-ALPHA-(ALPHA+BETA+2) X) Y' + N (N+ALPHA+BETA+1) Y = 0

  Recursion:

    P(0,ALPHA,BETA,X) = 1,

    P(1,ALPHA,BETA,X) = ( (2+ALPHA+BETA)*X + (ALPHA-BETA) ) / 2

    P(N,ALPHA,BETA,X)  =
      (
        (2*N+ALPHA+BETA-1)
        * ((ALPHA^2-BETA**2)+(2*N+ALPHA+BETA)*(2*N+ALPHA+BETA-2)*X)
        * P(N-1,ALPHA,BETA,X)
        -2*(N-1+ALPHA)*(N-1+BETA)*(2*N+ALPHA+BETA) * P(N-2,ALPHA,BETA,X)
      ) / 2*N*(N+ALPHA+BETA)*(2*N-2+ALPHA+BETA)

  Restrictions:

    -1 < ALPHA
    -1 < BETA

  Norm:

    Integral ( -1 <= X <= 1 ) ( 1 - X )^ALPHA * ( 1 + X )^BETA
      * P(N,ALPHA,BETA,X)^2 dX
    = 2^(ALPHA+BETA+1) * Gamma ( N + ALPHA + 1 ) * Gamma ( N + BETA + 1 ) /
      ( 2 * N + ALPHA + BETA ) * N! * Gamma ( N + ALPHA + BETA + 1 )

  Special values:

    P(N,ALPHA,BETA,1) = (N+ALPHA)!/(N!*ALPHA!) for integer ALPHA.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 August 2013

  Author:

    John Burkardt

  Reference:

    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.

  Parameters:

    Input, int M, the number of evaluation points.

    Input, int N, the highest order polynomial to compute.  Note
    that polynomials 0 through N will be computed.

    Input, double ALPHA, one of the parameters defining the Jacobi
    polynomials, ALPHA must be greater than -1.

    Input, double BETA, the second parameter defining the Jacobi
    polynomials, BETA must be greater than -1.

    Input, double X[M], the evaluation points.

    Output, double J_POLYNOMIAL[M*(N+1)], the values.
*/
{
  double c1;
  double c2;
  double c3;
  double c4;
  int i;
  int j;
  double *v;

  if ( alpha <= -1.0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "J_POLYNOMIAL - Fatal error!\n" );
    fprintf ( stderr, "  Illegal input value of ALPHA = %g\n", alpha );
    fprintf ( stderr, "  But ALPHA must be greater than -1.\n" );
    exit ( 1 );
  }

  if ( beta <= -1.0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "J_POLYNOMIAL - Fatal error!\n" );
    fprintf ( stderr, "  Illegal input value of BETA = %g\n", beta );
    fprintf ( stderr, "  But BETA must be greater than -1.\n" );
    exit ( 1 );
  }

  if ( n < 0 )
  {
    return NULL;
  }

  v = ( double * ) malloc ( m*(n+1) * sizeof ( double ) );

  for ( i = 0; i < m; i++ )
  {
    v[i+0*m] = 1.0;
  }

  if ( n == 0 )
  {
    return v;
  }

  for ( i = 0; i < m; i++ )
  {
    v[i+1*m] = ( 1.0 + 0.5 * ( alpha + beta ) ) * x[i]
      + 0.5 * ( alpha - beta );
  }

  for ( i = 0; i < m; i++ )
  {
    for ( j = 2; j <= n; j++ )
    {
      c1 = 2.0 * ( double ) ( j ) * ( ( double ) ( j ) + alpha + beta )
        * ( ( double ) ( 2 * j - 2 ) + alpha + beta );

      c2 = ( ( double ) ( 2 * j - 1 ) + alpha + beta )
        * ( ( double ) ( 2 * j ) + alpha + beta )
        * ( ( double ) ( 2 * j - 2 ) + alpha + beta );

      c3 = ( ( double ) ( 2 * j - 1 ) + alpha + beta )
        * ( alpha + beta ) * ( alpha - beta );

      c4 = - ( double ) ( 2 ) * ( ( double ) ( j - 1 ) + alpha )
        * ( ( double ) ( j - 1 ) + beta )
        * ( ( double ) ( 2 * j ) + alpha + beta );

      v[i+j*m] = ( ( c3 + c2 * x[i] ) * v[i+(j-1)*m] + c4 * v[i+(j-2)*m] ) / c1;
    }
  }

  return v;
}
/******************************************************************************/

void j_polynomial_values ( int *n_data, int *n, double *a, double *b, double *x,
  double *fx )

/******************************************************************************/
/*
  Purpose:

    J_POLYNOMIAL_VALUES returns some values of the Jacobi polynomial.

  Discussion:

    In Mathematica, the function can be evaluated by:

      JacobiP[ n, a, b, x ]

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 August 2013

  Author:

    John Burkardt

  Reference:

    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.

    Stephen Wolfram,
    The Mathematica Book,
    Fourth Edition,
    Cambridge University Press, 1999,
    ISBN: 0-521-64314-7,
    LC: QA76.95.W65.

  Parameters:

    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
    first call.  On each call, the routine increments N_DATA by 1, and
    returns the corresponding data; when there is no more data, the
    output value of N_DATA will be 0 again.

    Output, int *N, the degree of the polynomial.

    Output, double *A, *B, parameters of the function.

    Output, double *X, the argument of the function.

    Output, double *FX, the value of the function.
*/
{
# define N_MAX 26

  static double a_vec[N_MAX] = {
     0.0, 0.0, 0.0, 0,
     0.0, 0.0, 1.0, 2,
     3.0, 4.0, 5.0, 0,
     0.0, 0.0, 0.0, 0,
     0.0, 0.0, 0.0, 0,
     0.0, 0.0, 0.0, 0,
     0.0, 0.0 };

  static double b_vec[N_MAX] = {
    1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 2.0,
    3.0, 4.0, 5.0, 1.0,
    1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0,
    1.0, 1.0 };

  static double fx_vec[N_MAX] = {
      0.1000000000000000E+01,
      0.2500000000000000E+00,
     -0.3750000000000000E+00,
     -0.4843750000000000E+00,
     -0.1328125000000000E+00,
      0.2753906250000000E+00,
     -0.1640625000000000E+00,
     -0.1174804687500000E+01,
     -0.2361328125000000E+01,
     -0.2616210937500000E+01,
      0.1171875000000000E+00,
      0.4218750000000000E+00,
      0.5048828125000000E+00,
      0.5097656250000000E+00,
      0.4306640625000000E+00,
     -0.6000000000000000E+01,
      0.3862000000000000E-01,
      0.8118400000000000E+00,
      0.3666000000000000E-01,
     -0.4851200000000000E+00,
     -0.3125000000000000E+00,
      0.1891200000000000E+00,
      0.4023400000000000E+00,
      0.1216000000000000E-01,
     -0.4396200000000000E+00,
      0.1000000000000000E+01 };

  static int n_vec[N_MAX] = {
     0, 1, 2, 3,
     4, 5, 5, 5,
     5, 5, 5, 5,
     5, 5, 5, 5,
     5, 5, 5, 5,
     5, 5, 5, 5,
     5, 5 };

  static double x_vec[N_MAX] = {
      0.5E+00,
      0.5E+00,
      0.5E+00,
      0.5E+00,
      0.5E+00,
      0.5E+00,
      0.5E+00,
      0.5E+00,
      0.5E+00,
      0.5E+00,
      0.5E+00,
      0.5E+00,
      0.5E+00,
      0.5E+00,
      0.5E+00,
     -1.0E+00,
     -0.8E+00,
     -0.6E+00,
     -0.4E+00,
     -0.2E+00,
      0.0E+00,
      0.2E+00,
      0.4E+00,
      0.6E+00,
      0.8E+00,
      1.0E+00 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  *n_data = *n_data + 1;

  if ( N_MAX < *n_data )
  {
    *n_data = 0;
    *n = 0;
    *a = 0.0;
    *b = 0.0;
    *x = 0.0;
    *fx = 0.0;
  }
  else
  {
    *n = n_vec[*n_data-1];
    *a = a_vec[*n_data-1];
    *b = b_vec[*n_data-1];
    *x = x_vec[*n_data-1];
    *fx = fx_vec[*n_data-1];
  }

  return;
# undef N_MAX
}
/******************************************************************************/

double *j_polynomial_zeros ( int n, double alpha, double beta )

/******************************************************************************/
/*
  Purpose:

    J_POLYNOMIAL_ZEROS: zeros of Jacobi polynomial J(n,a,b,x).

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 August 2013

  Author:

    John Burkardt.

  Reference:

    Sylvan Elhay, Jaroslav Kautsky,
    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
    Interpolatory Quadrature,
    ACM Transactions on Mathematical Software,
    Volume 13, Number 4, December 1987, pages 399-415.

  Parameters:

    Input, int, N, the order.

    Input, double, ALPHA, BETA, the parameters.
    -1 < ALPHA, BETA.

    Output, double J_POLYNOMIAL_ZEROS[N], the zeros.
*/
{
  double a2b2;
  double ab;
  double abi;
  double *bj;
  int i;
  double i_r8;
  double *w;
  double *x;
  double zemu;

  ab = alpha + beta;
  abi = 2.0 + ab;
/*
  Define the zero-th moment.
*/
  zemu = pow ( 2.0, ab + 1.0 ) * tgamma ( alpha + 1.0 )
    * tgamma ( beta + 1.0 ) / tgamma ( abi );
/*
  Define the Jacobi matrix.
*/
  x = ( double * ) malloc ( n * sizeof ( double ) );
  x[0] = ( beta - alpha ) / abi;
  for ( i = 1; i < n; i++ )
  {
    x[i] = 0.0;
  }

  bj = ( double * ) malloc ( n * sizeof ( double ) );

  bj[0] = 4.0 * ( 1.0 + alpha ) * ( 1.0 + beta )
    / ( ( abi + 1.0 ) * abi * abi );
  for ( i = 1; i < n; i++ )
  {
    bj[i] = 0.0;
  }

  a2b2 = beta * beta - alpha * alpha;

  for ( i = 1; i < n; i++ )
  {
    i_r8 = ( double ) ( i + 1 );
    abi = 2.0 * i_r8 + ab;
    x[i] = a2b2 / ( ( abi - 2.0 ) * abi );
    abi = abi * abi;
    bj[i] = 4.0 * i_r8 * ( i_r8 + alpha ) * ( i_r8 + beta )
      * ( i_r8 + ab ) / ( ( abi - 1.0 ) * abi );
  }

  for ( i = 0; i < n; i++ )
  {
    bj[i] = sqrt ( bj[i] );
  }

  w = ( double * ) malloc ( n * sizeof ( double ) );

  w[0] = sqrt ( zemu );
  for ( i = 1; i < n; i++ )
  {
    w[i] = 0.0;
  }
/*
  Diagonalize the Jacobi matrix.
*/
  imtqlx ( n, x, bj, w );

  free ( bj );
  free ( w );

  return x;
}
/******************************************************************************/

void j_quadrature_rule ( int n, double alpha, double beta, double x[],
  double w[] )

/******************************************************************************/
/*
  Purpose:

    J_QUADRATURE_RULE: Gauss-Jacobi quadrature based on J(n,a,b,x).

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 August 2013

  Author:

    John Burkardt.

  Reference:

    Sylvan Elhay, Jaroslav Kautsky,
    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
    Interpolatory Quadrature,
    ACM Transactions on Mathematical Software,
    Volume 13, Number 4, December 1987, pages 399-415.

  Parameters:

    Input, int, N, the order.

    Input, double, ALPHA, BETA, the parameters.
    -1 < ALPHA, BETA.

    Output, double X[N], the abscissas.

    Output, double W[N], the weights.
*/
{
  double a2b2;
  double ab;
  double abi;
  double *bj;
  int i;
  double i_r8;
  double zemu;

  ab = alpha + beta;
  abi = 2.0 + ab;
/*
  Define the zero-th moment.
*/
  zemu = pow ( 2.0, ab + 1.0 ) * tgamma ( alpha + 1.0 )
    * tgamma ( beta + 1.0 ) / tgamma ( abi );
/*
  Define the Jacobi matrix.
*/
  x[0] = ( beta - alpha ) / abi;
  for ( i = 1; i < n; i++ )
  {
    x[i] = 0.0;
  }

  bj = ( double * ) malloc ( n * sizeof ( double ) );

  bj[0] = 4.0 * ( 1.0 + alpha ) * ( 1.0 + beta )
    / ( ( abi + 1.0 ) * abi * abi );
  for ( i = 1; i < n; i++ )
  {
    bj[i] = 0.0;
  }

  a2b2 = beta * beta - alpha * alpha;

  for ( i = 1; i < n; i++ )
  {
    i_r8 = ( double ) ( i + 1 );
    abi = 2.0 * i_r8 + ab;
    x[i] = a2b2 / ( ( abi - 2.0 ) * abi );
    abi = abi * abi;
    bj[i] = 4.0 * i_r8 * ( i_r8 + alpha ) * ( i_r8 + beta )
      * ( i_r8 + ab ) / ( ( abi - 1.0 ) * abi );
  }

  for ( i = 0; i < n; i++ )
  {
    bj[i] = sqrt ( bj[i] );
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

double r8_choose ( int n, int k )

/******************************************************************************/
/*
  Purpose:

    R8_CHOOSE computes the binomial coefficient C(N,K) as an R8.

  Discussion:

    The value is calculated in such a way as to avoid overflow and
    roundoff.  The calculation is done in R8 arithmetic.

    The formula used is:

      C(N,K) = N! / ( K! * (N-K)! )

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 July 2008

  Author:

    John Burkardt

  Reference:

    ML Wolfson, HV Wright,
    Algorithm 160:
    Combinatorial of M Things Taken N at a Time,
    Communications of the ACM,
    Volume 6, Number 4, April 1963, page 161.

  Parameters:

    Input, int N, K, the values of N and K.

    Output, double R8_CHOOSE, the number of combinations of N
    things taken K at a time.
*/
{
  int i;
  int mn;
  int mx;
  double value;

  if ( k < n - k )
  {
    mn = k;
    mx = n - k;
  }
  else
  {
    mn = n - k;
    mx = k;
  }

  if ( mn < 0 )
  {
    value = 0.0;
  }
  else if ( mn == 0 )
  {
    value = 1.0;
  }
  else
  {
    value = ( double ) ( mx + 1 );

    for ( i = 2; i <= mn; i++ )
    {
      value = ( value * ( double ) ( mx + i ) ) / ( double ) i;
    }
  }

  return value;
}
/******************************************************************************/

double r8_epsilon ( void )

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

double r8_factorial ( int n )

/******************************************************************************/
/*
  Purpose:

    R8_FACTORIAL computes the factorial of N.

  Discussion:

    factorial ( N ) = product ( 1 <= I <= N ) I

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 June 2008

  Author:

    John Burkardt

  Parameters:

    Input, int N, the argument of the factorial function.
    If N is less than 1, the function value is returned as 1.

    Output, double R8_FACTORIAL, the factorial of N.
*/
{
  int i;
  double value;

  value = 1.0;

  for ( i = 1; i <= n; i++ )
  {
    value = value * ( double ) ( i );
  }

  return value;
}
/******************************************************************************/

double r8_max ( double x, double y )

/******************************************************************************/
/*
  Purpose:

    R8_MAX returns the maximum of two R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 May 2006

  Author:

    John Burkardt

  Parameters:

    Input, double X, Y, the quantities to compare.

    Output, double R8_MAX, the maximum of X and Y.
*/
{
  double value;

  if ( y < x )
  {
    value = x;
  }
  else
  {
    value = y;
  }
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

void r8mat_print ( int m, int n, double a[], char *title )

/******************************************************************************/
/*
  Purpose:

    R8MAT_PRINT prints an R8MAT.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

    Entry A(I,J) is stored as A[I+J*M]

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    28 May 2008

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows in A.

    Input, int N, the number of columns in A.

    Input, double A[M*N], the M by N matrix.

    Input, char *TITLE, a title.
*/
{
  r8mat_print_some ( m, n, a, 1, 1, m, n, title );

  return;
}
/******************************************************************************/

void r8mat_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi,
  int jhi, char *title )

/******************************************************************************/
/*
  Purpose:

    R8MAT_PRINT_SOME prints some of an R8MAT.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 June 2013

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows of the matrix.
    M must be positive.

    Input, int N, the number of columns of the matrix.
    N must be positive.

    Input, double A[M*N], the matrix.

    Input, int ILO, JLO, IHI, JHI, designate the first row and
    column, and the last row and column to be printed.

    Input, char *TITLE, a title.
*/
{
# define INCX 5

  int i;
  int i2hi;
  int i2lo;
  int j;
  int j2hi;
  int j2lo;

  fprintf ( stdout, "\n" );
  fprintf ( stdout, "%s\n", title );

  if ( m <= 0 || n <= 0 )
  {
    fprintf ( stdout, "\n" );
    fprintf ( stdout, "  (None)\n" );
    return;
  }
/*
  Print the columns of the matrix, in strips of 5.
*/
  for ( j2lo = jlo; j2lo <= jhi; j2lo = j2lo + INCX )
  {
    j2hi = j2lo + INCX - 1;
    if ( n < j2hi )
    {
      j2hi = n;
    }
    if ( jhi < j2hi )
    {
      j2hi = jhi;
    }

    fprintf ( stdout, "\n" );
/*
  For each column J in the current range...

  Write the header.
*/
    fprintf ( stdout, "  Col:  ");
    for ( j = j2lo; j <= j2hi; j++ )
    {
      fprintf ( stdout, "  %7d     ", j - 1 );
    }
    fprintf ( stdout, "\n" );
    fprintf ( stdout, "  Row\n" );
    fprintf ( stdout, "\n" );
/*
  Determine the range of the rows in this strip.
*/
    if ( 1 < ilo )
    {
      i2lo = ilo;
    }
    else
    {
      i2lo = 1;
    }
    if ( m < ihi )
    {
      i2hi = m;
    }
    else
    {
      i2hi = ihi;
    }

    for ( i = i2lo; i <= i2hi; i++ )
    {
/*
  Print out (up to) 5 entries in row I, that lie in the current strip.
*/
      fprintf ( stdout, "%5d:", i - 1 );
      for ( j = j2lo; j <= j2hi; j++ )
      {
        fprintf ( stdout, "  %14f", a[i-1+(j-1)*m] );
      }
      fprintf ( stdout, "\n" );
    }
  }

  return;
# undef INCX
}
/******************************************************************************/

double r8vec_dot_product ( int n, double a1[], double a2[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_DOT_PRODUCT computes the dot product of a pair of R8VEC's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 July 2007

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the vectors.

    Input, double A1[N], A2[N], the two vectors to be considered.

    Output, double R8VEC_DOT_PRODUCT, the dot product of the vectors.
*/
{
  int i;
  double value;

  value = 0.0;
  for ( i = 0; i < n; i++ )
  {
    value = value + a1[i] * a2[i];
  }
  return value;
}
/******************************************************************************/

double *r8vec_linspace_new ( int n, double a, double b )

/******************************************************************************/
/*
  Purpose:

    R8VEC_LINSPACE_NEW creates a vector of linearly spaced values.

  Discussion:

    An R8VEC is a vector of R8's.

    4 points evenly spaced between 0 and 12 will yield 0, 4, 8, 12.
 
    In other words, the interval is divided into N-1 even subintervals,
    and the endpoints of intervals are used as the points.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    29 March 2011

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the vector.

    Input, double A, B, the first and last entries.

    Output, double R8VEC_LINSPACE_NEW[N], a vector of linearly spaced data.
*/
{
  int i;
  double *x;

  x = ( double * ) malloc ( n * sizeof ( double ) );

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
  return x;
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
    fprintf ( stdout, "  %8d: %14f\n", i, a[i] );
  }

  return;
}
/******************************************************************************/

void r8vec2_print ( int n, double a1[], double a2[], char *title )

/******************************************************************************/
/*
  Purpose:

    R8VEC2_PRINT prints an R8VEC2.

  Discussion:

    An R8VEC2 is a dataset consisting of N pairs of real values, stored
    as two separate vectors A1 and A2.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 March 2009

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of components of the vector.

    Input, double A1[N], double A2[N], the vectors to be printed.

    Input, char *TITLE, a title.
*/
{
  int i;

  fprintf ( stdout, "\n" );
  fprintf ( stdout, "%s\n", title );
  fprintf ( stdout, "\n" );
  for ( i = 0; i < n; i++ )
  {
    fprintf ( stdout, "  %4d: %14f  %14f\n", i, a1[i], a2[i] );
  }

  return;
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

  fprintf ( stdout, "%s\n", time_buffer );

  return;
# undef TIME_SIZE
}
