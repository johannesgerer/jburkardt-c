# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <string.h>

# include "quadmom.h"
# include "toms655.h"

/******************************************************************************/

void jacobi_eigenvalue ( int n, double a[], int it_max, double v[], 
  double d[], int *it_num, int *rot_num )

/******************************************************************************/
/*
  Purpose:

    JACOBI_EIGENVALUE carries out the Jacobi eigenvalue iteration.

  Discussion:

    This function computes the eigenvalues and eigenvectors of a
    real symmetric matrix, using Rutishauser's modfications of the classical
    Jacobi rotation method with threshold pivoting. 

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 September 2013

  Author:

    C version by John Burkardt

  Parameters:

    Input, int N, the order of the matrix.

    Input, double A[N*N], the matrix, which must be square, real,
    and symmetric.

    Input, int IT_MAX, the maximum number of iterations.

    Output, double V[N*N], the matrix of eigenvectors.

    Output, double D[N], the eigenvalues, in descending order.

    Output, int *IT_NUM, the total number of iterations.

    Output, int *ROT_NUM, the total number of rotations.
*/
{
  double *bw;
  double c;
  double g;
  double gapq;
  double h;
  int i;
  int j;
  int k;
  int l;
  int m;
  int p;
  int q;
  double s;
  double t;
  double tau;
  double term;
  double termp;
  double termq;
  double theta;
  double thresh;
  double w;
  double *zw;

  r8mat_identity ( n, v );

  r8mat_diag_get_vector ( n, a, d );

  bw = ( double * ) malloc ( n * sizeof ( double ) );
  zw = ( double * ) malloc ( n * sizeof ( double ) );

  for ( i = 0; i < n; i++ )
  {
    bw[i] = d[i];
    zw[i] = 0.0;
  }
  *it_num = 0;
  *rot_num = 0;

  while ( *it_num < it_max )
  {
    *it_num = *it_num + 1;
/*
  The convergence threshold is based on the size of the elements in
  the strict upper triangle of the matrix.
*/
    thresh = 0.0;
    for ( j = 0; j < n; j++ )
    {
      for ( i = 0; i < j; i++ )
      {
        thresh = thresh + a[i+j*n] * a[i+j*n];
      }
    }

    thresh = sqrt ( thresh ) / ( double ) ( 4 * n );

    if ( thresh == 0.0 )
    {
      break;
    }

    for ( p = 0; p < n; p++ )
    {
      for ( q = p + 1; q < n; q++ )
      {
        gapq = 10.0 * fabs ( a[p+q*n] );
        termp = gapq + fabs ( d[p] );
        termq = gapq + fabs ( d[q] );
/*
  Annihilate tiny offdiagonal elements.
*/
        if ( 4 < *it_num &&
             termp == fabs ( d[p] ) &&
             termq == fabs ( d[q] ) )
        {
          a[p+q*n] = 0.0;
        }
/*
  Otherwise, apply a rotation.
*/
        else if ( thresh <= fabs ( a[p+q*n] ) )
        {
          h = d[q] - d[p];
          term = fabs ( h ) + gapq;

          if ( term == fabs ( h ) )
          {
            t = a[p+q*n] / h;
          }
          else
          {
            theta = 0.5 * h / a[p+q*n];
            t = 1.0 / ( fabs ( theta ) + sqrt ( 1.0 + theta * theta ) );
            if ( theta < 0.0 )
            {
              t = - t;
            }
          }
          c = 1.0 / sqrt ( 1.0 + t * t );
          s = t * c;
          tau = s / ( 1.0 + c );
          h = t * a[p+q*n];
/*
  Accumulate corrections to diagonal elements.
*/
          zw[p] = zw[p] - h;                 
          zw[q] = zw[q] + h;
          d[p] = d[p] - h;
          d[q] = d[q] + h;

          a[p+q*n] = 0.0;
/*
  Rotate, using information from the upper triangle of A only.
*/
          for ( j = 0; j < p; j++ )
          {
            g = a[j+p*n];
            h = a[j+q*n];
            a[j+p*n] = g - s * ( h + g * tau );
            a[j+q*n] = h + s * ( g - h * tau );
          }

          for ( j = p + 1; j < q; j++ )
          {
            g = a[p+j*n];
            h = a[j+q*n];
            a[p+j*n] = g - s * ( h + g * tau );
            a[j+q*n] = h + s * ( g - h * tau );
          }

          for ( j = q + 1; j < n; j++ )
          {
            g = a[p+j*n];
            h = a[q+j*n];
            a[p+j*n] = g - s * ( h + g * tau );
            a[q+j*n] = h + s * ( g - h * tau );
          }
/*
  Accumulate information in the eigenvector matrix.
*/
          for ( j = 0; j < n; j++ )
          {
            g = v[j+p*n];
            h = v[j+q*n];
            v[j+p*n] = g - s * ( h + g * tau );
            v[j+q*n] = h + s * ( g - h * tau );
          }
          *rot_num = *rot_num + 1;
        }
      }
    }

    for ( i = 0; i < n; i++ )
    {
      bw[i] = bw[i] + zw[i];
      d[i] = bw[i];
      zw[i] = 0.0;
    }
  }
/*
  Restore upper triangle of input matrix.
*/
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < j; i++ )
    {
      a[i+j*n] = a[j+i*n];
    }
  }
/*
  Ascending sort the eigenvalues and eigenvectors.
*/
  for ( k = 0; k < n - 1; k++ )
  {
    m = k;
    for ( l = k + 1; l < n; l++ )
    {
      if ( d[l] < d[m] )
      {
        m = l;
      }
    }

    if ( m != k )
    {
      t    = d[m];
      d[m] = d[k];
      d[k] = t;
      for ( i = 0; i < n; i++ )
      {
        w        = v[i+m*n];
        v[i+m*n] = v[i+k*n];
        v[i+k*n] = w;
      }
    }
  }

  free ( bw );
  free ( zw );

  return;
}
/******************************************************************************/

void moment_method ( int n, double moment[], double x[], double w[] )

/******************************************************************************/
/*
  Purpose:

    MOMENT_METHOD computes a quadrature rule by the method of moments.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    18 September 2013

  Author:

    John Burkardt

  Reference:

    Gene Golub, John Welsch,
    Calculation of Gaussian Quadrature Rules,
    Mathematics of Computation,
    Volume 23, Number 106, April 1969, pages 221-230.

  Parameters:

    Input, int N, the order of the quadrature rule.

    Input, double MOMENT[2*N+1], moments 0 through 2*N.

    Output, double X[N], W[N], the points and weights of the quadrature rule.
*/
{
  double *alpha;
  double *beta;
  int debug;
  int flag;
  double *h;
  int i;
  int it_max;
  int it_num;
  int j;
  double *jacobi;
  double *r;
  int rot_num;
  double *v;

  debug = 0;

  if ( debug )
  {
    r8vec_print ( 2 * n + 1, moment, "  Moments:" );
  }

/*
  Define the N+1 by N+1 Hankel matrix H(I,J) = moment(I+J).
*/
  h = ( double * ) malloc ( ( n + 1 ) * ( n + 1 ) * sizeof ( double ) );

  for ( i = 0; i <= n; i++ )
  {
    for ( j = 0; j <= n; j++ )
    {
      h[i+j*(n+1)] = moment[i+j];
    }
  }

  if ( debug )
  {
    r8mat_print ( n + 1, n + 1, h, "  Hankel matrix:" );
  }
/*
  Compute R, the upper triangular Cholesky factor of H.
*/
  r = r8mat_cholesky_factor_upper ( n + 1, h, &flag );

  if ( flag != 0 )
  {
    printf ( "\n" );
    printf ( "MOMENT_METHOD - Fatal error!\n" );
    printf ( "  R8MAT_CHOLESKY_FACTOR_UPPER returned FLAG = %d\n", flag );
    exit ( 1 );
  }

  if ( debug )
  {
    r8mat_print ( n + 1, n + 1, r, "  Cholesky factor:" );
  }
/*
  Compute ALPHA and BETA from R, using Golub and Welsch's formula.
*/
  alpha = ( double * ) malloc ( n * sizeof ( double ) );

  alpha[0] = r[0+1*(n+1)] / r[0+0*(n+1)];
  for ( i = 1; i < n; i++ )
  {
    alpha[i] = r[i+(i+1)*(n+1)] / r[i+i*(n+1)] - r[i-1+i*(n+1)] / r[i-1+(i-1)*(n+1)];
  }
  beta = ( double * ) malloc ( ( n - 1 ) * sizeof ( double ) );

  for ( i = 0; i < n - 1; i++ )
  {
    beta[i] = r[i+1+(i+1)*(n+1)] / r[i+i*(n+1)];
  }
/*
  Compute the points and weights from the moments.
*/
  jacobi = ( double * ) malloc ( n * n * sizeof ( double ) );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      jacobi[i+j*n] = 0.0;
    }
  }

  for ( i = 0; i < n; i++ )
  {
    jacobi[i+i*n] = alpha[i];
  }

  for ( i = 0; i < n - 1; i++ )
  {
    jacobi[i+(i+1)*n] = beta[i];
    jacobi[i+1+i*n] = beta[i];
  }

  if ( debug )
  {
    r8mat_print ( n, n, jacobi, "  The Jacobi matrix:" );
  }
/*
  Get the eigendecomposition of the Jacobi matrix.
*/
  it_max = 100;
  v = ( double * ) malloc ( n * n * sizeof ( double ) );

  jacobi_eigenvalue ( n, jacobi, it_max, v, x, &it_num, &rot_num );

  if ( debug )
  {
    r8mat_print ( n, n, v, "  Eigenvector" );
  }

  for ( i = 0; i < n; i++ )
  {
    w[i] = moment[0] * pow ( v[0+i*n], 2 );
  }

  free ( alpha );
  free ( beta );
  free ( h );
  free ( jacobi );
  free ( r );
  free ( v );

  return;
}
/******************************************************************************/

double *moments_laguerre ( int m )

/******************************************************************************/
/*
  Purpose:

    MOMENTS_LAGUERRE returns moments of the Laguerre distribution.

  Discussion:

    pdf(x) = exp ( -x )
    mu(k) = integral ( 0 <= x < +oo ) x^k pdf(x) dx

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 September 2013

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of moments desired.

    Output, double W(0:M-1), the weighted integrals of X^0 
    through X^(M-1).
*/
{
  int k;
  double *w;

  w = ( double * ) malloc ( m * sizeof ( double ) );

  for ( k = 0; k < m; k++ )
  {
    w[k] = r8_factorial ( k );
  }

  return w;
}
/******************************************************************************/

double *moments_legendre ( int m, double a, double b )

/******************************************************************************/
/*
  Purpose:

    MOMENTS_LEGENDRE returns moments of the Legendre weight on [A,B].

  Discussion:

    mu(k) = integral ( a <= x <= b ) x^k dx

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 September 2013

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of moments desired.

    Input, double A, B, the left and right endpoints 
    of the interval.

    Output, double W(0:M-1), the weighted integrals of X^0 
    through X^(M-1).
*/
{
  double ak;
  double bk;
  int k;
  double *w;

  w = ( double * ) malloc ( m * sizeof ( double ) );

  bk = 1.0;
  ak = 1.0;
  for ( k = 0; k < m; k++ )
  {
    bk = bk * b;
    ak = ak * a;
    w[k] = ( bk - ak ) / ( double ) ( k + 1 );
  }

  return w;
}
/******************************************************************************/

double *moments_normal_01 ( int m )

/******************************************************************************/
/*
  Purpose:

    MOMENTS_NORMAL_01 returns moments of the standard Normal distribution.

  Discussion:

    pdf(x) = exp ( -x^2/2 ) / sqrt ( pi * 2 )
    mu(k) = integral ( -oo < x < +oo ) x^k pdf(x) dx

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 September 2013

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of moments desired.

    Output, double W(0:M-1), the weighted integrals of X^0 
    through X^(M-1).
*/
{
  int k;
  double *w;

  w = ( double * ) malloc ( m * sizeof ( double ) );

  w[0] = 1.0;

  for ( k = 2; k < m; k = k + 2 )
  {
    w[k] = r8_factorial2 ( k - 1 );
  }

  for ( k = 1; k < m; k = k + 2 )
  {
    w[k] = 0.0;
  }

  return w;
}
/******************************************************************************/

double *moments_normal ( int m, double mu, double sigma )

/******************************************************************************/
/*
  Purpose:

    MOMENTS_NORMAL returns moments of the standard Normal distribution.

  Discussion:

    pdf(x) = exp ( -((x-mu)/sigma)^2/2 ) / sigma / sqrt ( pi * 2 )
    mu(k) = integral ( -oo < x < +oo ) x^k pdf(x) dx

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 September 2013

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of moments desired.

    Input, double MU, SIGMA, the mean and standard deviation.

    Output, double W(0:M-1), the weighted integrals of X^0 
    through X^(M-1).
*/
{
  int j;
  int j_hi;
  int k;
  double t;
  double *w;

  w = ( double * ) malloc ( m * sizeof ( double ) );

  for ( k = 0; k < m; k++ )
  {
    t = 0.0;
    j_hi = k / 2;
    for ( j = 0; j <= j_hi; j++ )
    {
      t = t + r8_choose ( k, 2 * j ) * r8_factorial2 ( 2 * j - 1 ) 
        * pow ( sigma, 2 * j ) * pow ( mu, k - 2 * j );
    }
    w[k] = t;
  }

  return w;
}
/******************************************************************************/

double *moments_truncated_normal_ab ( int m, double mu, double sigma,
  double a, double b )

/******************************************************************************/
/*
  Purpose:

    MOMENTS_TRUNCATED_NORMAL_AB: moments of the truncated Normal distribution.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 September 2013

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of moments desired.

    Input, double MU, SIGMA, the mean and standard deviation.

    Input, double A, B, the lower and upper truncation limits.

    Output, double W(0:M-1), the weighted integrals of X^0 
    through X^(M-1).
*/
{
  int order;
  double *w;

  w = ( double * ) malloc ( m * sizeof ( double ) );

  for ( order = 0; order < m; order++ )
  {
    w[order] = truncated_normal_ab_moment ( order, mu, sigma, a, b );
  }

  return w;
}
/******************************************************************************/

double *moments_truncated_normal_a ( int m, double mu, double sigma,
  double a )

/******************************************************************************/
/*
  Purpose:

    MOMENTS_TRUNCATED_NORMAL_A: moments of the lower truncated Normal.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 September 2013

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of moments desired.

    Input, double MU, SIGMA, the mean and standard deviation.

    Input, double A, the lower truncation limit.

    Output, double W(0:M-1), the weighted integrals of X^0 
    through X^(M-1).
*/
{
  int order;
  double *w;

  w = ( double * ) malloc ( m * sizeof ( double ) );

  for ( order = 0; order < m; order++ )
  {
    w[order] = truncated_normal_a_moment ( order, mu, sigma, a );
  }

  return w;
}
/******************************************************************************/

double *moments_truncated_normal_b ( int m, double mu, double sigma,
  double b )

/******************************************************************************/
/*
  Purpose:

    MOMENTS_TRUNCATED_NORMAL_B: moments of the upper truncated Normal.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 September 2013

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of moments desired.

    Input, double MU, SIGMA, the mean and standard deviation.

    Input, double B, the upper truncation limit.

    Output, double W(0:M-1), the weighted integrals of X^0 
    through X^(M-1).
*/
{
  int order;
  double *w;

  w = ( double * ) malloc ( m * sizeof ( double ) );

  for ( order = 0; order < m; order++ )
  {
    w[order] = truncated_normal_b_moment ( order, mu, sigma, b );
  }

  return w;
}
/******************************************************************************/

double normal_01_cdf ( double x )

/******************************************************************************/
/*
  Purpose:

    NORMAL_01_CDF evaluates the Normal 01 CDF.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    10 February 1999

  Author:

    John Burkardt

  Reference:

    A G Adams,
    Areas Under the Normal Curve,
    Algorithm 39,
    Computer j.,
    Volume 12, pages 197-198, 1969.

  Parameters:

    Input, double X, the argument of the CDF.

    Output, double CDF, the value of the CDF.
*/
{
  double a1 = 0.398942280444;
  double a2 = 0.399903438504;
  double a3 = 5.75885480458;
  double a4 = 29.8213557808;
  double a5 = 2.62433121679;
  double a6 = 48.6959930692;
  double a7 = 5.92885724438;
  double b0 = 0.398942280385;
  double b1 = 3.8052E-08;
  double b2 = 1.00000615302;
  double b3 = 3.98064794E-04;
  double b4 = 1.98615381364;
  double b5 = 0.151679116635;
  double b6 = 5.29330324926;
  double b7 = 4.8385912808;
  double b8 = 15.1508972451;
  double b9 = 0.742380924027;
  double b10 = 30.789933034;
  double b11 = 3.99019417011;
  double cdf;
  double q;
  double y;
/*
  |X| <= 1.28.
*/
  if ( fabs ( x ) <= 1.28 )
  {
    y = 0.5 * x * x;

    q = 0.5 - fabs ( x ) * ( a1 - a2 * y / ( y + a3 - a4 / ( y + a5
      + a6 / ( y + a7 ) ) ) );
/*
  1.28 < |X| <= 12.7
*/
  }
  else if ( fabs ( x ) <= 12.7 )
  {
    y = 0.5 * x * x;

    q = exp ( - y ) * b0 / ( fabs ( x ) - b1
      + b2 / ( fabs ( x ) + b3
      + b4 / ( fabs ( x ) - b5
      + b6 / ( fabs ( x ) + b7
      - b8 / ( fabs ( x ) + b9
      + b10 / ( fabs ( x ) + b11 ) ) ) ) ) );
/*
  12.7 < |X|
*/
  }
  else
  {
    q = 0.0;
  }
/*
  Take account of negative X.
*/
  if ( x < 0.0 )
  {
    cdf = q;
  }
  else
  {
    cdf = 1.0 - q;
  }

  return cdf;
}
/******************************************************************************/

double normal_01_pdf ( double x )

/******************************************************************************/
/*
  Purpose:

    NORMAL_01_PDF evaluates the Normal 01 PDF.

  Discussion:

    The Normal 01 PDF is also called the "Standard Normal" PDF, or
    the Normal PDF with 0 mean and variance 1.

    PDF(X) = exp ( - 0.5 * X^2 ) / sqrt ( 2 * PI )

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    18 September 2004

  Author:

    John Burkardt

  Parameters:

    Input, double X, the argument of the PDF.

    Output, double PDF, the value of the PDF.
*/
{
  double pdf;
  const double pi = 3.14159265358979323;

  pdf = exp ( -0.5 * x * x ) / sqrt ( 2.0 * pi );

  return pdf;
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

double r8_factorial2 ( int n )

/******************************************************************************/
/*
  Purpose:

    R8_FACTORIAL2 computes the double factorial function.

  Discussion:

    FACTORIAL2( N ) = Product ( N * (N-2) * (N-4) * ... * 2 )  (N even)
                    = Product ( N * (N-2) * (N-4) * ... * 1 )  (N odd)

     N Value
    -- -----
     0     1
     1     1
     2     2
     3     3
     4     8
     5    15
     6    48
     7   105
     8   384
     9   945
    10  3840

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 July 2008

  Author:

    John Burkardt

  Parameters:

    Input, int N, the argument of the double factorial
    function.  If N is less than 1, R8_FACTORIAL2 is returned as 1.0.

    Output, double R8_FACTORIAL2, the value of Factorial2(N).
*/
{
  int n_copy;
  double value;

  value = 1.0;

  if ( n < 1 )
  {
    return value;
  }

  n_copy = n;

  while ( 1 < n_copy )
  {
    value = value * ( double ) n_copy;
    n_copy = n_copy - 2;
  }

  return value;
}
/******************************************************************************/

double r8_mop ( int i )

/******************************************************************************/
/*
  Purpose:

    R8_MOP returns the I-th power of -1 as an R8 value.

  Discussion:

    An R8 is an double value.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 July 2008

  Author:

    John Burkardt

  Parameters:

    Input, int I, the power of -1.

    Output, double R8_MOP, the I-th power of -1.
*/
{
  double value;

  if ( ( i % 2 ) == 0 )
  {
    value = + 1.0;
  }
  else
  {
    value = - 1.0;
  }

  return value;
}
/******************************************************************************/

double *r8mat_cholesky_factor_upper ( int n, double a[], int *flag )

/******************************************************************************/
/*
  Purpose:

    R8MAT_CHOLESKY_FACTOR_UPPER: upper Cholesky factor of a symmetric R8MAT.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

    The matrix must be symmetric and positive semidefinite.

    For a positive semidefinite symmetric matrix A, the Cholesky factorization
    is an upper triangular matrix R such that:

      A = R' * R

    Note that the usual Cholesky factor is a LOWER triangular matrix L
    such that

      A = L * L'

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    03 August 2013

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of rows and columns of the matrix A.

    Input, double A[N*N], the N by N matrix.

    Output, int *FLAG, an error flag.
    0, no error was detected.
    1, the matrix was not positive definite.  A NULL factor was returned.

    Output, double R8MAT_CHOLESKY_FACTOR_UPPER[N*N], the N by N upper triangular
    "Choresky" factor.
*/
{
  double *c;
  int i;
  int j;
  int k;
  double sum2;

  *flag = 0;

  c = r8mat_copy_new ( n, n, a );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < j; i++ )
    {
      c[j+i*n] = 0.0;
    }
    for ( i = j; i < n; i++ )
    {
      sum2 = c[i+j*n];
      for ( k = 0; k < j; k++ )
      {
        sum2 = sum2 - c[k+j*n] * c[k+i*n];
      }
      if ( i == j )
      {
        if ( sum2 <= 0.0 )
        {
          *flag = 1;
          return NULL;
        }
        c[j+i*n] = sqrt ( sum2 );
      }
      else
      {
        if ( c[j+j*n] != 0.0 )
        {
          c[j+i*n] = sum2 / c[j+j*n];
        }
        else
        {
          c[j+i*n] = 0.0;
        }
      }
    }
  }

  return c;
}
/******************************************************************************/

double *r8mat_copy_new ( int m, int n, double a1[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_COPY_NEW copies one R8MAT to a "new" R8MAT.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 July 2008

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns.

    Input, double A1[M*N], the matrix to be copied.

    Output, double R8MAT_COPY_NEW[M*N], the copy of A1.
*/
{
  double *a2;
  int i;
  int j;

  a2 = ( double * ) malloc ( m * n * sizeof ( double ) );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      a2[i+j*m] = a1[i+j*m];
    }
  }

  return a2;
}
/******************************************************************************/

void r8mat_diag_get_vector ( int n, double a[], double v[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_DIAG_GET_VECTOR gets the value of the diagonal of an R8MAT.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 July 2013

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of rows and columns of the matrix.

    Input, double A[N*N], the N by N matrix.

    Output, double V[N], the diagonal entries
    of the matrix.
*/
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    v[i] = a[i+i*n];
  }

  return;
}
/******************************************************************************/

void r8mat_identity  ( int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_IDENTITY sets an R8MAT to the identity matrix.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    06 September 2005

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of A.

    Output, double A[N*N], the N by N identity matrix.
*/
{
  int i;
  int j;
  int k;

  k = 0;
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      if ( i == j )
      {
        a[k] = 1.0;
      }
      else
      {
        a[k] = 0.0;
      }
      k = k + 1;
    }
  }

  return;
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
        fprintf ( stdout, "  %14g", a[i-1+(j-1)*m] );
      }
      fprintf ( stdout, "\n" );
    }
  }

  return;
# undef INCX
}
/******************************************************************************/

void r8vec_print_dupe ( int n, double a[], char *title )

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
    fprintf ( stdout, "  %4d: %14f  %14g\n", i, a1[i], a2[i] );
  }

  return;
}
/******************************************************************************/

void timestamp_dupe ( void )

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

double truncated_normal_ab_moment ( int order, double mu, double s, double a,
  double b )

/******************************************************************************/
/*
  Purpose:

    TRUNCATED_NORMAL_AB_MOMENT: moments of the truncated Normal PDF.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 September 2013

  Author:

    John Burkardt

  Reference:

    Phoebus Dhrymes,
    Moments of Truncated Normal Distributions,
    May 2005.

  Parameters:

    Input, int ORDER, the order of the moment.
    0 <= ORDER.

    Input, double MU, S, the mean and standard deviation of the
    parent Normal distribution.
    0.0 < S.

    Input, double A, B, the lower and upper truncation limits.
    A < B.

    Output, double TRUNCATED_NORMAL_AB_MOMENT, the moment of the PDF.
*/
{
  double a_h;
  double a_cdf;
  double a_pdf;
  double b_h;
  double b_cdf;
  double b_pdf;
  double ir;
  double irm1;
  double irm2;
  double moment;
  int r;

  if ( order < 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "TRUNCATED_NORMAL_AB_MOMENT - Fatal error!\n" );
    fprintf ( stderr, "  ORDER < 0.\n" );
    exit ( 1 );
  }

  if ( s <= 0.0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "TRUNCATED_NORMAL_AB_MOMENT - Fatal error!\n" );
    fprintf ( stderr, "  S <= 0.\n" );
    exit ( 1 );
  }

  if ( b <= a )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "TRUNCATED_NORMAL_AB_MOMENT - Fatal error!\n" );
    fprintf ( stderr, "  B <= A.\n" );
    exit ( 1 );
  }

  a_h = ( a - mu ) / s;
  a_pdf = normal_01_pdf ( a_h );
  a_cdf = normal_01_cdf ( a_h );

  if ( a_cdf == 0.0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "TRUNCATED_NORMAL_AB_MOMENT - Fatal error!\n" );
    fprintf ( stderr, "  PDF/CDF ratio fails, A_CDF too small.\n" );
    fprintf ( stderr, "  A_PDF = %g\n", a_pdf );
    fprintf ( stderr, "  A_CDF = %g\n", a_cdf );
    exit ( 1 );
  }

  b_h = ( b - mu ) / s;
  b_pdf = normal_01_pdf ( b_h );
  b_cdf = normal_01_cdf ( b_h );

  if ( b_cdf == 0.0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "TRUNCATED_NORMAL_AB_MOMENT - Fatal error!\n" );
    fprintf ( stderr, "  PDF/CDF ratio fails, B_CDF too small.\n" );
    fprintf ( stderr, "  B_PDF = %g\n", b_pdf );
    fprintf ( stderr, "  B_CDF = %g\n", b_cdf );
    exit ( 1 );
  }

  moment = 0.0;
  irm2 = 0.0;
  irm1 = 0.0;

  for ( r = 0; r <= order; r++ )
  {
    if ( r == 0 )
    {
      ir = 1.0;
    }
    else if ( r == 1 )
    {
      ir = - ( b_pdf - a_pdf ) / ( b_cdf - a_cdf );
    }
    else
    {
      ir = ( double ) ( r - 1 ) * irm2 
        - ( pow ( b_h, r - 1 ) * b_pdf - pow ( a_h, r - 1 ) * a_pdf )
        / ( b_cdf - a_cdf );
    }

    moment = moment + r8_choose ( order, r ) * pow ( mu, order - r ) 
      * pow ( s, r ) * ir;

    irm2 = irm1;
    irm1 = ir;
  }

  return moment;
}
/******************************************************************************/

double truncated_normal_a_moment ( int order, double mu, double s, double a )

/******************************************************************************/
/*
  Purpose:

    TRUNCATED_NORMAL_A_MOMENT: moments of the lower truncated Normal PDF.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 September 2013

  Author:

    John Burkardt

  Reference:

    Phoebus Dhrymes,
    Moments of Truncated Normal Distributions,
    May 2005.

  Parameters:

    Input, int ORDER, the order of the moment.
    0 <= ORDER.

    Input, double MU, S, the mean and standard deviation of the
    parent Normal distribution.

    Input, double A, the lower truncation limit.

    Output, double TRUNCATED_NORMAL_A_MOMENT, the moment of the PDF.
*/
{
  double moment;

  moment = r8_mop ( order ) 
    * truncated_normal_b_moment ( order, - mu, s, - a );

  return moment;
}
/******************************************************************************/

double truncated_normal_b_moment ( int order, double mu, double s, double b )

/******************************************************************************/
/*
  Purpose:

    TRUNCATED_NORMAL_B_MOMENT: moments of the upper truncated Normal PDF.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 September 2013

  Author:

    John Burkardt

  Reference:

    Phoebus Dhrymes,
    Moments of Truncated Normal Distributions,
    May 2005.

  Parameters:

    Input, int ORDER, the order of the moment.
    0 <= ORDER.

    Input, double MU, S, the mean and standard deviation of the
    parent Normal distribution.

    Input, double B, the upper truncation limit.

    Output, double TRUNCATED_NORMAL_B_MOMENT, the moment of the PDF.
*/
{
  double f;
  double h;
  double h_cdf;
  double h_pdf;
  double ir;
  double irm1;
  double irm2;
  double moment;
  int r;

  if ( order < 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "TRUNCATED_NORMAL_B_MOMENT - Fatal error!\n" );
    fprintf ( stderr, "  ORDER < 0.\n" );
    exit ( 1 );
  }

  h = ( b - mu ) / s;
  h_pdf = normal_01_pdf ( h );
  h_cdf = normal_01_cdf ( h );

  if ( h_cdf == 0.0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "TRUNCATED_NORMAL_B_MOMENT - Fatal error!\n" );
    fprintf ( stderr, "  CDF((B-MU)/S) = 0.\n" );
    exit ( 1 );
  }

  f = h_pdf / h_cdf;

  moment = 0.0;
  irm2 = 0.0;
  irm1 = 0.0;

  for ( r = 0; r <= order; r++ )
  {
    if ( r == 0 )
    {
      ir = 1.0;
    }
    else if ( r == 1 )
    {
      ir = - f;
    }
    else
    {
      ir = - pow ( h, r - 1 ) * f + ( double ) ( r - 1 ) * irm2;
    }

    moment = moment + r8_choose ( order, r ) * pow ( mu, order - r ) 
      * pow ( s, r ) * ir;

    irm2 = irm1;
    irm1 = ir;
  }

  return moment;
}
