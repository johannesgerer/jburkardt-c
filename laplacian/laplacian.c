# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <string.h>

# include "laplacian.h"

/******************************************************************************/

double cholesky_upper_error ( int n, double a[], double c[] )

/******************************************************************************/
/*
  Purpose:

    CHOLESKY_UPPER_ERROR determines the error in an upper Cholesky factor.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 October 2013

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of the matrix.

    Input, double A[N*N], the matrix.

    Input, double C[N*N], the upper triangular Cholesky factor.

    Output, double CHOLESKY_UPPER_ERROR, the Frobenius norm
    of the difference matrix A - C' * C.
*/
{
  double *ctc;
  double *d;
  double value;

  ctc = r8mat_mtm_new ( n, n, n, c, c );

  d = r8mat_sub_new ( n, n, a, ctc );
 
  value = r8mat_norm_fro ( n, n, d );

  free ( ctc );
  free ( d );

  return value;
}
/******************************************************************************/

double eigen_error ( int n, int k, double a[], double x[], double lambda[] )

/******************************************************************************/
/*
  Purpose:

    EIGEN_ERROR determines the error in a (right) eigensystem.

  Discussion:

    An R8MAT is a matrix of double values.

    This routine computes the Frobenius norm of

      A * X - X * LAMBDA

    where

      A is an N by N matrix,
      X is an N by K matrix (each of K columns is an eigenvector)
      LAMBDA is a K by K diagonal matrix of eigenvalues.

    This routine assumes that A, X and LAMBDA are all real!

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 October 2013

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of the matrix.

    Input, int K, the number of eigenvectors.
    K is usually 1 or N.

    Input, double A[N*N], the matrix.

    Input, double X[N*K]), the K eigenvectors.

    Input, double LAMBDA[K], the K eigenvalues.

    Output, double EIGEN_ERROR, the Frobenius norm
    of the difference matrix A * X - X * LAMBDA, which would be exactly zero
    if X and LAMBDA were exact eigenvectors and eigenvalues of A.
*/
{
  double *c;
  int i;
  int j;
  double value;

  c = r8mat_mm_new ( n, n, k, a, x );

  for ( j = 0; j < k; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      c[i+n*j] = c[i+n*j] - lambda[j] * x[i+n*j];
    }
  }

  value = r8mat_norm_fro ( n, k, c );

  free ( c );

  return value;
}
/******************************************************************************/

double inverse_error ( int n, double a[], double b[] )

/******************************************************************************/
/*
  Purpose:

    INVERSE_ERROR determines the error in an inverse matrix.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 October 2013

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of the matrix.

    Input, double A[N*N], the matrix.

    Input, double B[N*N], the inverse.

    Output, double ERROR_FROBENIUS, the Frobenius norm
    of (A*B-I) + (B*A-I).
*/
{
  double *c;
  int j;
  double value;

  c = r8mat_mm_new ( n, n, n, a, b );

  for ( j = 0; j < n; j++ )
  {
    c[j+j*n] = c[j+j*n] - 1.0;
  }

  value = r8mat_norm_fro ( n, n, c );

  free ( c );

  c = r8mat_mm_new ( n, n, n, b, a );

  for ( j = 0; j < n; j++ )
  {
    c[j+j*n] = c[j+j*n] - 1.0;
  }

  value = value + r8mat_norm_fro ( n, n, c );

  free ( c );

  return value;
}
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

double *l1dd_apply ( int n, double h, double u[] )

/******************************************************************************/
/*
  Purpose:

    L1DD_APPLY applies the 1D DD Laplacian to a vector.

  Discussion:

    The N grid points are assumed to be evenly spaced by H.

    For N = 5, the discrete Laplacian with Dirichlet boundary conditions
    at both ends of [0,6] is applied to a vector of 7 values, with a spacing
    of H = 6/(N+1) = 1 at the points X:

      0  1  2  3  4  5  6

    and has the matrix form L:

       2 -1  0  0  0
      -1  2 -1  0  0
       0 -1  2 -1  0
       0  0 -1  2 -1
       0  0  0 -1  2

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 October 2013

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of points.
    N must be at least 3.

    Input, double H, the spacing between points.

    Input, double U[N], the value at each point.

    Output, double L1DD_APPLY[N], the Laplacian evaluated at each point.
*/
{
  int i;
  double *lu;

  if ( n < 3 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "L1DD_APPLY - Fatal error!\n" );
    fprintf ( stderr, "  N < 3.\n" );
    exit ( 1 );
  }

  lu = ( double * ) malloc ( n * sizeof ( double ) );

  i = 0;
  lu[i] = ( 2.0 * u[i] - u[i+1] ) / h / h;
  for ( i = 1; i < n - 1; i++ )
  {
    lu[i] = ( - u[i-1] + 2.0 * u[i] - u[i+1] ) / h / h;
  }
  i = n - 1;
  lu[i] = ( - u[i-1] + 2.0 * u[i] ) / h / h;

  return lu;
}
/******************************************************************************/

double *l1dd_cholesky ( int n, double h )

/******************************************************************************/
/*
  Purpose:

    L1DD_CHOLESKY computes the Cholesky factor of the 1D DD Laplacian.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 October 2013

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of points.
    N must be at least 3.

    Input, double H, the spacing between points.

    Output, double L1DD_CHOLESKY[N*N], the Cholesky factor.
*/
{
  double *c;
  int i;
  int j;

  if ( n < 3 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "L1DD_CHOLESKY - Fatal error!\n" );
    fprintf ( stderr, "  N < 3.\n" );
    exit ( 1 );
  }

  c = ( double * ) malloc ( n * n * sizeof ( double ) );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      c[i+j*n] = 0.0;
    }
  }

  for ( i = 0; i < n; i++ )
  {
    c[i+i*n] = sqrt ( ( double ) ( i + 2 ) ) / sqrt ( ( double ) ( i + 1 ) );
  }

  for ( i = 0; i < n - 1; i++ )
  {
    c[i+(i+1)*n] = - sqrt ( ( double ) ( i + 1 ) ) 
      / sqrt ( ( double ) ( i + 2 ) );
  }

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      c[i+j*n] = c[i+j*n] / h;
    }
  }

  return c;
}
/******************************************************************************/

void l1dd_eigen ( int n, double h, double v[], double lambda[] )

/******************************************************************************/
/*
  Purpose:

    L1DD_EIGEN returns eigeninformation for the 1D DD Laplacian.

  Discussion:

    The grid points are assumed to be evenly spaced by H.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 October 2013

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of points.

    Input, double H, the spacing between points.

    Output, double V[N*N], the eigenvectors.

    Output, double LAMBDA[N], the eigenvalues.
*/
{
  int i;
  double i_r8;
  int j;
  double j_r8;
  double n_r8;
  const double pi = 3.141592653589793;
  double s;
  double theta;

  n_r8 = ( double ) ( n );

  for ( j = 0; j < n; j++ )
  {
    j_r8 = ( double ) ( j + 1 );
    theta = 0.5 * pi * j_r8 / ( n_r8 + 1.0 );
    lambda[j] = pow ( 2.0 * sin ( theta ) / h, 2 );
    for ( i = 0; i < n; i++ )
    {
      i_r8 = ( double ) ( i + 1 );
      theta = pi * i_r8 * j_r8 / ( n_r8 + 1.0 );
      v[i+j*n] = sqrt ( 2.0 / ( n_r8 + 1.0 ) ) * sin ( theta );
    }
  }

  return;
}
/******************************************************************************/

double *l1dd ( int n, double h )

/******************************************************************************/
/*
  Purpose:

    L1DD stores the 1D DD Laplacian as a full matrix.

  Discussion:

    The N grid points are assumed to be evenly spaced by H.

    For N = 5, the discrete Laplacian with Dirichlet boundary conditions
    at both ends of [0,6] has the matrix form L:

       2 -1  0  0  0
      -1  2 -1  0  0
       0 -1  2 -1  0
       0  0 -1  2 -1
       0  0  0 -1  2

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 October 2013

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of points.
    N must be at least 3.

    Input, double H, the spacing between points.

    Output, double L1DD[N*N], the Laplacian matrix.
*/
{
  int i;
  int j;
  double *l;

  if ( n < 3 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "L1DD - Fatal error!\n" );
    fprintf ( stderr, "  N < 3.\n" );
    exit ( 1 );
  }

  l = ( double * ) malloc ( n * n * sizeof ( double ) );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      l[i+j*n] = 0.0;
    }
  }

  i = 0;
  l[i+i*n]     =  2.0 / h / h;
  l[i+(i+1)*n] = -1.0 / h / h;

  for ( i = 1; i < n - 1; i++ )
  {
    l[i+(i-1)*n] = -1.0 / h / h;
    l[i+i*n] =      2.0 / h / h;
    l[i+(i+1)*n] = -1.0 / h / h;
  }

  i = n - 1;
  l[i+(i-1)*n] = -1.0 / h / h;
  l[i+i*n] =      2.0 / h / h;

  return l;
}
/******************************************************************************/

double *l1dd_inverse ( int n, double h )

/******************************************************************************/
/*
  Purpose:

    L1DD_INVERSE stores the inverse of the 1D DD Laplacian.

  Discussion:

    The N grid points are assumed to be evenly spaced by H.

    For N = 5, the discrete Laplacian with Dirichlet boundary conditions
    at both ends of [0,6] has the matrix form L:

       2 -1  0  0  0
      -1  2 -1  0  0
       0 -1  2 -1  0
       0  0 -1  2 -1
       0  0  0 -1  2

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 October 2013

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of points.
    N must be at least 3.

    Input, double H, the spacing between points.

    Output, double L1DD_INVERSE[N*N], the inverse of the Laplacian matrix.
*/
{
  int i;
  int j;
  double *l;

  if ( n < 3 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "L1DD_INVERSE - Fatal error!\n" );
    fprintf ( stderr, "  N < 3.\n" );
    exit ( 1 );
  }

  l = ( double * ) malloc ( n * n * sizeof ( double ) );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      l[i+j*n] = ( double ) ( i4_min ( i + 1, j + 1 ) * ( n - i4_max ( i, j ) ) ) 
             * h * h / ( double ) ( n + 1 );
    }
  }

  return l;
}
/******************************************************************************/

void l1dd_lu ( int n, double h, double l[], double u[] )

/******************************************************************************/
/*
  Purpose:

    L1DD_LU computes the LU factors of the 1D DD Laplacian.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 November 2013

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of points.
    N must be at least 3.

    Input, double H, the spacing between points.

    Output, double L[N*N], U[N*N], the LU factors.
*/
{
  int i;
  double i_r8;
  int j;

  if ( n < 3 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "L1DD_LU - Fatal error!\n" );
    fprintf ( stderr, "  N < 3.\n" );
    exit ( 1 );
  }

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      l[i+j*n] = 0.0;
    }
  }

  for ( i = 0; i < n; i++ )
  {
    l[i+i*n] = 1.0;
  }

  for ( i = 1; i < n; i++ )
  {
    i_r8 = ( double ) ( i + 1 );
    l[i+(i-1)*n] = - ( i_r8 - 1.0 ) / i_r8;
  }

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      l[i+j*n] = l[i+j*n] / h;
    }
  }

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      u[i+j*n] = 0.0;
    }
  }

  for ( i = 0; i < n; i++ )
  {
    i_r8 = ( double ) ( i + 1 );
    u[i+i*n] = ( i_r8 + 1.0 ) / i_r8;
  }

  for ( i = 0; i < n - 1; i++ )
  {
    u[i+(i+1)*n] = - 1.0;
  }

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      u[i+j*n] = u[i+j*n] / h;
    }
  }

  return;
}
/******************************************************************************/

double *l1dn_apply ( int n, double h, double u[] )

/******************************************************************************/
/*
  Purpose:

    L1DN_APPLY applies the 1D DN Laplacian to a vector.

  Discussion:

    The N grid points are assumed to be evenly spaced by H.

    For N = 5, the discrete Laplacian with left Dirichlet and right
    Neumann condition on [0,6] has the matrix form L:

       2 -1  0  0  0
      -1  2 -1  0  0
       0 -1  2 -1  0
       0  0 -1  2 -1
       0  0  0 -1  1

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 October 2013

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of points.
    N must be at least 3.

    Input, double H, the spacing between points.

    Input, double U[N], the value at each point.

    Output, double L1DN_APPLY[N], the Laplacian evaluated at each point.
*/
{
  int i;
  double *lu;

  if ( n < 3 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "L1DN_APPLY - Fatal error!\n" );
    fprintf ( stderr, "  N < 3.\n" );
    exit ( 1 );
  }

  lu = ( double * ) malloc ( n * sizeof ( double ) );

  i = 0;
  lu[i] = ( 2.0 * u[i] - u[i+1] ) / h / h;
  for ( i = 1; i < n - 1; i++ )
  {
    lu[i] = ( - u[i-1] + 2.0 * u[i] - u[i+1] ) / h / h;
  }
  i = n - 1;
  lu[i] = ( - u[i-1] + u[i] ) / h / h;

  return lu;
}
/******************************************************************************/

double *l1dn_cholesky ( int n, double h )

/******************************************************************************/
/*
  Purpose:

    L1DN_CHOLESKY computes the Cholesky factor of the 1D DN Laplacian.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 October 2013

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of points.
    N must be at least 3.

    Input, double H, the spacing between points.

    Output, double L1DN_CHOLESKY[N*N], the Cholesky factor.
*/
{
  double *c;
  int i;
  double i_r8;
  int j;

  if ( n < 3 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "L1DN_CHOLESKY - Fatal error!\n" );
    fprintf ( stderr, "  N < 3.\n" );
    exit ( 1 );
  }

  c = ( double * ) malloc ( n * n * sizeof ( double ) );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      c[i+j*n] = 0.0;
    }
  }

  for ( i = 0; i < n - 1; i++ )
  {
    i_r8 = ( double ) ( i + 1 );
    c[i+i*n]   =   sqrt ( i_r8 + 1.0 ) / sqrt ( i_r8 );
    c[i+(i+1)*n] = - sqrt ( i_r8 ) / sqrt ( i_r8 + 1.0 );
  }
  c[n-1+(n-1)*n] = 1.0 / sqrt ( ( double ) ( n ) );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      c[i+j*n] = c[i+j*n] / h;
    }
  }

  return c;
}
/******************************************************************************/

void l1dn_eigen ( int n, double h, double v[], double lambda[] )

/******************************************************************************/
/*
  Purpose:

    L1DN_EIGEN returns eigeninformation for the 1D DN Laplacian.

  Discussion:

    The grid points are assumed to be evenly spaced by H.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 October 2013

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of points.

    Input, double H, the spacing between points.

    Output, double V[N*N], the eigenvectors.

    Output, double LAMBDA[N], the eigenvalues.
*/
{
  int i;
  double i_r8;
  int j;
  double j_r8;
  double n_r8;
  const double pi = 3.141592653589793;
  double s;
  double theta;

  n_r8 = ( double ) ( n );

  for ( j = 0; j < n; j++ )
  {
    j_r8 = ( double ) ( j + 1 );
    theta = pi * ( j_r8 - 0.5 ) / ( 2.0 * n_r8 + 1.0 );
    lambda[j] = pow ( 2.0 * sin ( theta ) / h, 2 );
    for ( i = 0; i < n; i++ )
    {
      i_r8 = ( double ) ( i + 1 );
      theta = pi * i_r8 * ( 2.0 * j_r8 - 1.0 ) / 
        ( 2.0 * n_r8 + 1.0 );
      v[i+j*n] = sqrt ( 2.0 / ( n_r8 + 0.5 ) ) * sin ( theta );
    }
  }

  return;
}
/******************************************************************************/

double *l1dn ( int n, double h )

/******************************************************************************/
/*
  Purpose:

    L1DN stores the 1D DN Laplacian as a full matrix.

  Discussion:

    The N grid points are assumed to be evenly spaced by H.

    For N = 5, the discrete Laplacian with left Dirichlet and right
    Neumann condition on [0,6] has the matrix form L:

       2 -1  0  0  0
      -1  2 -1  0  0
       0 -1  2 -1  0
       0  0 -1  2 -1
       0  0  0 -1  1

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 October 2013

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of points.
    N must be at least 3.

    Input, double H, the spacing between points.

    Output, double L1DN[N*N], the Laplacian matrix.
*/
{
  int i;
  int j;
  double *l;

  if ( n < 3 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "L1DN - Fatal error!\n" );
    fprintf ( stderr, "  N < 3.\n" );
    exit ( 1 );
  }

  l = ( double * ) malloc ( n * n * sizeof ( double ) );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      l[i+j*n] = 0.0;
    }
  }

  i = 0;
  l[i+i*n]     =  2.0 / h / h;
  l[i+(i+1)*n] = -1.0 / h / h;

  for ( i = 1; i < n - 1; i++ )
  {
    l[i+(i-1)*n] = -1.0 / h / h;
    l[i+i*n] =      2.0 / h / h;
    l[i+(i+1)*n] = -1.0 / h / h;
  }

  i = n - 1;
  l[i+(i-1)*n] = -1.0 / h / h;
  l[i+i*n] =      1.0 / h / h;

  return l;
}
/******************************************************************************/

double *l1dn_inverse ( int n, double h )

/******************************************************************************/
/*
  Purpose:

    L1DN_INVERSE stores the inverse of the 1D DN Laplacian.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 October 2013

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of points.
    N must be at least 3.

    Input, double H, the spacing between points.

    Output, double L1DN_INVERSE[N*N], the inverse of the Laplacian matrix.
*/
{
  int i;
  int j;
  double *l;

  if ( n < 3 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "L1DN_INVERSE - Fatal error!\n" );
    fprintf ( stderr, "  N < 3.\n" );
    exit ( 1 );
  }

  l = ( double * ) malloc ( n * n * sizeof ( double ) );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      l[i+j*n] = ( double ) ( i4_min ( i + 1, j + 1 ) ) * h * h;
    }
  }

  return l;
}
/******************************************************************************/

void l1dn_lu ( int n, double h, double l[], double u[] )

/******************************************************************************/
/*
  Purpose:

    L1DD_LU computes the LU factors of the 1D DN Laplacian.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 November 2013

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of points.
    N must be at least 3.

    Input, double H, the spacing between points.

    Output, double L[N*N], U[N*N], the LU factors.
*/
{
  int i;
  double i_r8;
  int j;

  if ( n < 3 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "L1DN_LU - Fatal error!\n" );
    fprintf ( stderr, "  N < 3.\n" );
    exit ( 1 );
  }

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      l[i+j*n] = 0.0;
    }
  }

  for ( i = 0; i < n; i++ )
  {
    l[i+i*n] = 1.0;
  }

  for ( i = 1; i < n; i++ )
  {
    i_r8 = ( double ) ( i + 1 );
    l[i+(i-1)*n] = - ( i_r8 - 1.0 ) / i_r8;
  }

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      l[i+j*n] = l[i+j*n] / h;
    }
  }

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      u[i+j*n] = 0.0;
    }
  }

  for ( i = 0; i < n - 1; i++ )
  {
    i_r8 = ( double ) ( i + 1 );
    u[i+i*n] = ( i_r8 + 1.0 ) / i_r8;
  }
  i = n - 1;
  i_r8 = ( double ) ( i + 1 );
  u[i+i*n] = 1.0 / i_r8;

  for ( i = 0; i < n - 1; i++ )
  {
    u[i+(i+1)*n] = - 1.0;
  }

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      u[i+j*n] = u[i+j*n] / h;
    }
  }

  return;
}
/******************************************************************************/

double *l1nd_apply ( int n, double h, double u[] )

/******************************************************************************/
/*
  Purpose:

    L1ND_APPLY applies the 1D ND Laplacian to a vector.

  Discussion:

    The N grid points are assumed to be evenly spaced by H.

    For N = 5, the discrete Laplacian with left Neumann and right Dirichlet
    boundary conditions on [0,6] has the matrix form L:

       1 -1  0  0  0
      -1  2 -1  0  0
       0 -1  2 -1  0
       0  0 -1  2 -1
       0  0  0 -1  2

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 October 2013

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of points.
    N must be at least 3.

    Input, double H, the spacing between points.

    Input, double U[N], the value at each point.

    Output, double L1ND_APPLY[N], the Laplacian evaluated at each point.
*/
{
  int i;
  double *lu;

  if ( n < 3 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "L1ND_APPLY - Fatal error!\n" );
    fprintf ( stderr, "  N < 3.\n" );
    exit ( 1 );
  }

  lu = ( double * ) malloc ( n * sizeof ( double ) );

  i = 0;
  lu[i] = ( u[i] - u[i+1] ) / h / h;
  for ( i = 1; i < n - 1; i++ )
  {
    lu[i] = ( - u[i-1] + 2.0 * u[i] - u[i+1] ) / h / h;
  }
  i = n - 1;
  lu[i] = ( - u[i-1] + 2.0 * u[i] ) / h / h;

  return lu;
}
/******************************************************************************/

double *l1nd_cholesky ( int n, double h )

/******************************************************************************/
/*
  Purpose:

    L1ND_CHOLESKY computes the Cholesky factor of the 1D ND Laplacian.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 October 2013

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of points.
    N must be at least 3.

    Input, double H, the spacing between points.

    Output, double L1ND_CHOLESKY[N*N], the Cholesky factor.
*/
{
  double *c;
  int i;
  double i_r8;
  int j;

  if ( n < 3 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "L1DN_CHOLESKY - Fatal error!\n" );
    fprintf ( stderr, "  N < 3.\n" );
    exit ( 1 );
  }

  c = ( double * ) malloc ( n * n * sizeof ( double ) );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      c[i+j*n] = 0.0;
    }
  }

  for ( i = 0; i < n; i++ )
  {
    c[i+i*n] = 1.0;
  }

  for ( i = 0; i < n - 1; i++ )
  {
    c[i+(i+1)*n] = - 1.0;
  }

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      c[i+j*n] = c[i+j*n] / h;
    }
  }

  return c;
}
/******************************************************************************/

void l1nd_eigen ( int n, double h, double v[], double lambda[] )

/******************************************************************************/
/*
  Purpose:

    L1ND_EIGEN returns eigeninformation for the 1D ND Laplacian.

  Discussion:

    The grid points are assumed to be evenly spaced by H.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 October 2013

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of points.

    Input, double H, the spacing between points.

    Output, double V[N*N], the eigenvectors.

    Output, double LAMBDA[N], the eigenvalues.
*/
{
  int i;
  double i_r8;
  int j;
  double j_r8;
  double n_r8;
  const double pi = 3.141592653589793;
  double s;
  double theta;

  n_r8 = ( double ) ( n );

  for ( j = 0; j < n; j++ )
  {
    j_r8 = ( double ) ( j + 1 );
    theta = pi * ( j_r8 - 0.5 ) / ( 2.0 * n_r8 + 1.0 );
    lambda[j] = pow ( 2.0 * sin ( theta ) / h, 2 );
    for ( i = 0; i < n; i++ )
    {
      i_r8 = ( double ) ( i + 1 );
      theta = pi * ( i_r8 - 0.5 ) * ( 2.0 * j_r8 - 1.0 ) / 
        ( 2.0 * n_r8 + 1.0 );
      v[i+j*n] = sqrt ( 2.0 / ( n_r8 + 0.5 ) ) * cos ( theta );
    }
  }

  return;
}
/******************************************************************************/

double *l1nd ( int n, double h )

/******************************************************************************/
/*
  Purpose:

    L1ND stores the 1D ND Laplacian as a full matrix.

  Discussion:

    The N grid points are assumed to be evenly spaced by H.

    For N = 5, the discrete Laplacian with left Neumann and right Dirichlet
    boundary conditions on [0,6] has the matrix form L:

       1 -1  0  0  0
      -1  2 -1  0  0
       0 -1  2 -1  0
       0  0 -1  2 -1
       0  0  0 -1  2

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    29 October 2013

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of points.
    N must be at least 3.

    Input, double H, the spacing between points.

    Output, double L1ND[N*N], the Laplacian matrix.
*/
{
  int i;
  int j;
  double *l;

  if ( n < 3 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "L1ND - Fatal error!\n" );
    fprintf ( stderr, "  N < 3.\n" );
    exit ( 1 );
  }

  l = ( double * ) malloc ( n * n * sizeof ( double ) );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      l[i+j*n] = 0.0;
    }
  }

  i = 0;
  l[i+i*n]     =  1.0 / h / h;
  l[i+(i+1)*n] = -1.0 / h / h;

  for ( i = 1; i < n - 1; i++ )
  {
    l[i+(i-1)*n] = -1.0 / h / h;
    l[i+i*n] =      2.0 / h / h;
    l[i+(i+1)*n] = -1.0 / h / h;
  }

  i = n - 1;
  l[i+(i-1)*n] = -1.0 / h / h;
  l[i+i*n] =      2.0 / h / h;

  return l;
}
/******************************************************************************/

double *l1nd_inverse ( int n, double h )

/******************************************************************************/
/*
  Purpose:

    L1ND_INVERSE stores the inverse of the 1D ND Laplacian.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 October 2013

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of points.
    N must be at least 3.

    Input, double H, the spacing between points.

    Output, double L1ND_INVERSE[N*N], the inverse of the Laplacian matrix.
*/
{
  int i;
  int j;
  double *l;

  if ( n < 3 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "L1ND_INVERSE - Fatal error!\n" );
    fprintf ( stderr, "  N < 3.\n" );
    exit ( 1 );
  }

  l = ( double * ) malloc ( n * n * sizeof ( double ) );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      l[i+j*n] = ( double ) ( n - i4_max ( i, j ) ) * h * h;
    }
  }

  return l;
}
/******************************************************************************/

void l1nd_lu ( int n, double h, double l[], double u[] )

/******************************************************************************/
/*
  Purpose:

    L1ND_LU computes the LU factors of the 1D ND Laplacian.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 November 2013

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of points.
    N must be at least 3.

    Input, double H, the spacing between points.

    Output, double L[N*N], U[N*N], the LU factors.
*/
{
  int i;
  int j;

  if ( n < 3 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "L1ND_LU - Fatal error!\n" );
    fprintf ( stderr, "  N < 3.\n" );
    exit ( 1 );
  }

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      l[i+j*n] = 0.0;
    }
  }

  for ( i = 0; i < n; i++ )
  {
    l[i+i*n] = 1.0;
  }

  for ( i = 1; i < n; i++ )
  {
    l[i+(i-1)*n] = - 1.0;
  }

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      l[i+j*n] = l[i+j*n] / h;
    }
  }

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      u[i+j*n] = 0.0;
    }
  }

  for ( i = 0; i < n; i++ )
  {
    u[i+i*n] = 1.0;
  }

  for ( i = 0; i < n - 1; i++ )
  {
    u[i+(i+1)*n] = - 1.0;
  }

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      u[i+j*n] = u[i+j*n] / h;
    }
  }

  return;
}
/******************************************************************************/

double *l1nn_apply ( int n, double h, double u[] )

/******************************************************************************/
/*
  Purpose:

    L1NN_APPLY applies the 1D NN Laplacian to a vector.

  Discussion:

    The N grid points are assumed to be evenly spaced by H.

    For N = 5, the discrete Laplacian with left Neumann and right Neumann
    boundary conditions on [0,6] has the matrix form L:

       1 -1  0  0  0
      -1  2 -1  0  0
       0 -1  2 -1  0
       0  0 -1  2 -1
       0  0  0 -1  1

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 October 2013

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of points.
    N must be at least 3.

    Input, double H, the spacing between points.

    Input, double U[N], the value at each point.

    Output, double L1NN_APPLY[N], the Laplacian evaluated at each point.
*/
{
  int i;
  double *lu;

  if ( n < 3 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "L1NN_APPLY - Fatal error!\n" );
    fprintf ( stderr, "  N < 3.\n" );
    exit ( 1 );
  }

  lu = ( double * ) malloc ( n * sizeof ( double ) );

  i = 0;
  lu[i] = ( u[i] - u[i+1] ) / h / h;
  for ( i = 1; i < n - 1; i++ )
  {
    lu[i] = ( - u[i-1] + 2.0 * u[i] - u[i+1] ) / h / h;
  }
  i = n - 1;
  lu[i] = ( - u[i-1] +  u[i] ) / h / h;

  return lu;
}
/******************************************************************************/

double *l1nn_cholesky ( int n, double h )

/******************************************************************************/
/*
  Purpose:

    L1NN_CHOLESKY computes the Cholesky factor of the 1D NN Laplacian.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 October 2013

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of points.
    N must be at least 3.

    Input, double H, the spacing between points.

    Output, double L1NN_CHOLESKY[N*N], the Cholesky factor.
*/
{
  double *c;
  int i;
  int j;

  if ( n < 3 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "L1NN_CHOLESKY - Fatal error!\n" );
    fprintf ( stderr, "  N < 3.\n" );
    exit ( 1 );
  }

  c = ( double * ) malloc ( n * n * sizeof ( double ) );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      c[i+j*n] = 0.0;
    }
  }

  for ( i = 0; i < n - 1; i++ )
  {
    c[i+i*n]     = + 1.0;
    c[i+(i+1)*n] = - 1.0;
  }

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      c[i+j*n] = c[i+j*n] / h;
    }
  }

  return c;
}
/******************************************************************************/

void l1nn_eigen ( int n, double h, double v[], double lambda[] )

/******************************************************************************/
/*
  Purpose:

    L1NN_EIGEN returns eigeninformation for the 1D NN Laplacian.

  Discussion:

    The grid points are assumed to be evenly spaced by H.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 October 2013

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of points.

    Input, double H, the spacing between points.

    Output, double V[N*N], the eigenvectors.

    Output, double LAMBDA[N], the eigenvalues.
*/
{
  int i;
  double i_r8;
  int j;
  double j_r8;
  double n_r8;
  const double pi = 3.141592653589793;
  double s;
  double theta;

  n_r8 = ( double ) ( n );

  for ( j = 0; j < n; j++ )
  {
    j_r8 = ( double ) ( j + 1 );
    theta = pi * ( j_r8 - 1.0 ) / ( 2.0 * n_r8 );
    lambda[j] = pow ( 2.0 * sin ( theta ) / h, 2 );
    if ( j == 0 )
    {
      for ( i = 0; i < n; i++ )
      {
        v[i+j*n] = sqrt ( n_r8 );
      }
    }
    else
    {
      for ( i = 0; i < n; i++ )
      {
        i_r8 = ( double ) ( i + 1 );
        theta = pi * ( i_r8 - 0.5 ) * ( j_r8 - 1.0 ) / n_r8;
        v[i+j*n] = sqrt ( 2.0 / n_r8 ) * cos ( theta );
      }
    }
  }

  return;
}
/******************************************************************************/

double *l1nn ( int n, double h )

/******************************************************************************/
/*
  Purpose:

    L1NN stores the 1D NN Laplacian as a full matrix.

  Discussion:

    The N grid points are assumed to be evenly spaced by H.

    For N = 5, the discrete Laplacian with Neumann boundary conditions
    at both ends of [0,6] has the matrix form L:

       1 -1  0  0  0
      -1  2 -1  0  0
       0 -1  2 -1  0
       0  0 -1  2 -1
       0  0  0 -1  1

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 October 2013

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of points.
    N must be at least 3.

    Input, double H, the spacing between points.

    Output, double L1NN[N*N], the Laplacian matrix.
*/
{
  int i;
  int j;
  double *l;

  if ( n < 3 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "L1NN - Fatal error!\n" );
    fprintf ( stderr, "  N < 3.\n" );
    exit ( 1 );
  }

  l = ( double * ) malloc ( n * n * sizeof ( double ) );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      l[i+j*n] = 0.0;
    }
  }

  i = 0;
  l[i+i*n]     =  1.0 / h / h;
  l[i+(i+1)*n] = -1.0 / h / h;

  for ( i = 1; i < n - 1; i++ )
  {
    l[i+(i-1)*n] = -1.0 / h / h;
    l[i+i*n] =      2.0 / h / h;
    l[i+(i+1)*n] = -1.0 / h / h;
  }

  i = n - 1;
  l[i+(i-1)*n] = -1.0 / h / h;
  l[i+i*n] =      1.0 / h / h;

  return l;
}
/******************************************************************************/

void l1nn_lu ( int n, double h, double l[], double u[] )

/******************************************************************************/
/*
  Purpose:

    L1NN_LU computes the LU factors of the 1D NN Laplacian.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 November 2013

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of points.
    N must be at least 3.

    Input, double H, the spacing between points.

    Output, double L[N*N], U[N*N], the LU factors.
*/
{
  int i;
  int j;

  if ( n < 3 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "L1NN_LU - Fatal error!\n" );
    fprintf ( stderr, "  N < 3.\n" );
    exit ( 1 );
  }

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      l[i+j*n] = 0.0;
    }
  }

  for ( i = 0; i < n; i++ )
  {
    l[i+i*n] = 1.0;
  }

  for ( i = 1; i < n; i++ )
  {
    l[i+(i-1)*n] = - 1.0;
  }

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      l[i+j*n] = l[i+j*n] / h;
    }
  }

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      u[i+j*n] = 0.0;
    }
  }

  for ( i = 0; i < n - 1; i++ )
  {
    u[i+i*n] = 1.0;
  }
  i = n - 1;
  u[i+i*n] = 0.0;

  for ( i = 0; i < n - 1; i++ )
  {
    u[i+(i+1)*n] = - 1.0;
  }

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      u[i+j*n] = u[i+j*n] / h;
    }
  }

  return;
}
/******************************************************************************/

double *l1pp_apply ( int n, double h, double u[] )

/******************************************************************************/
/*
  Purpose:

    L1PP_APPLY applies the 1D PP Laplacian to a vector.

  Discussion:

    The N grid points are assumed to be evenly spaced by H.

    For N = 5, the discrete Laplacian with periodic boundary conditions
    on [0,6] has the matrix form L:

       2 -1  0  0 -1
      -1  2 -1  0  0
       0 -1  2 -1  0
       0  0 -1  2 -1
      -1  0  0 -1  2

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 October 2013

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of points.
    N must be at least 3.

    Input, double H, the spacing between points.

    Input, double U[N], the value at each point.

    Output, double L1PP_APPLY[N], the Laplacian evaluated at each point.
*/
{
  int i;
  double *lu;

  if ( n < 3 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "L1PP_APPLY - Fatal error!\n" );
    fprintf ( stderr, "  N < 3.\n" );
    exit ( 1 );
  }

  lu = ( double * ) malloc ( n * sizeof ( double ) );

  i = 0;
  lu[i] = ( - u[n-1] + 2.0 * u[i] - u[i+1] ) / h / h;
  for ( i = 1; i < n - 1; i++ )
  {
    lu[i] = ( - u[i-1] + 2.0 * u[i] - u[i+1] ) / h / h;
  }
  i = n - 1;
  lu[i] = ( - u[i-1] + 2.0 * u[i] - u[0] ) / h / h;

  return lu;
}
/******************************************************************************/

double *l1pp_cholesky ( int n, double h )

/******************************************************************************/
/*
  Purpose:

    L1PP_CHOLESKY computes the Cholesky factor of the 1D PP Laplacian.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 October 2013

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of points.
    N must be at least 3.

    Input, double H, the spacing between points.

    Output, double L1PP_CHOLESKY[N*N], the Cholesky factor.
*/
{
  double *c;
  int i;
  double i_r8;
  int j;

  if ( n < 3 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "L1PP_CHOLESKY - Fatal error!\n" );
    fprintf ( stderr, "  N < 3.\n" );
    exit ( 1 );
  }

  c = ( double * ) malloc ( n * n * sizeof ( double ) );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      c[i+j*n] = 0.0;
    }
  }

  for ( i = 0; i < n - 1; i++ )
  {
    i_r8 = ( double ) ( i + 1 );
    c[i+i*n] = sqrt ( i_r8 + 1.0 ) / sqrt ( i_r8 );
  }

  for ( i = 0; i < n - 2; i++ )
  {
    i_r8 = ( double ) ( i + 1 );
    c[i+(i+1)*n] = - i_r8 / ( i_r8 + 1.0 ) * sqrt ( i_r8 + 1.0 ) 
      / sqrt ( i_r8 );
  }

  for ( i = 0; i < n - 2; i++ )
  {
    i_r8 = ( double ) ( i + 1 );
    c[i+(n-1)*n] = - 1.0 / ( i_r8 + 1.0 ) * sqrt ( i_r8 + 1.0 )  
      / sqrt ( i_r8 );
  }

  i = n - 2;
  i_r8 = ( double ) ( i + 1 );
  c[i+(n-1)*n] = - ( double ) ( n ) / ( i_r8 + 1.0 ) 
    * sqrt ( i_r8 + 1.0 ) / sqrt ( i_r8 );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      c[i+j*n] = c[i+j*n] / h;
    }
  }

  return c;
}
/******************************************************************************/

void l1pp_eigen ( int n, double h, double v[], double lambda[] )

/******************************************************************************/
/*
  Purpose:

    L1PP_EIGEN returns eigeninformation for the 1D PP Laplacian.

  Discussion:

    The grid points are assumed to be evenly spaced by H.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 October 2013

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of points.

    Input, double H, the spacing between points.

    Output, double V[N*N], the eigenvectors.

    Output, double LAMBDA[N], the eigenvalues.
*/
{
  int i;
  double i_r8;
  int j;
  double j_r8;
  double n_r8;
  const double pi = 3.141592653589793;
  double s;
  double theta;

  n_r8 = ( double ) ( n );

  for ( j = 0; j < n; j++ )
  {
    j_r8 = ( double ) ( j + 1 );
    if ( ( j % 2 ) == 0 )
    {
      theta = pi * ( j_r8 - 1.0 ) / ( 2.0 * n_r8 );
    }
    else
    {
      theta = pi *   j_r8         / ( 2.0 * n_r8 );
    }
    lambda[j] = pow ( 2.0 * sin ( theta ) / h, 2 );

    if ( ( j % 2 ) == 0 )
    {
      if ( j == 0 )
      {
        for ( i = 0; i < n; i++ )
        {
          v[i+j*n] = 1.0 / sqrt ( n_r8 );
        }
      }
      else
      {
        for ( i = 0; i < n; i++ )
        {
          i_r8 = ( double ) ( i + 1 );
          theta = pi * ( i_r8 - 0.5 ) * ( j_r8 - 1.0 ) /  n_r8;
          v[i+j*n] = sqrt ( 2.0 / n_r8 ) * cos ( theta );
        }
      }
    }
    else
    {
      if ( j == n - 1 )
      {
        s = - 1.0 / sqrt ( n_r8 );
        for ( i = 0; i < n; i++ )
        {
          v[i+j*n] = s;
          s = - s;
        }
      }
      else
      {
        for ( i = 0; i < n; i++ )
        {
          i_r8 = ( double ) ( i + 1 );
          theta = pi * ( i_r8 - 0.5 ) * j_r8 / n_r8;
          v[i+j*n] = sqrt ( 2.0 / n_r8 ) * sin ( theta );
        }
      }
    }

  }

  return;
}
/******************************************************************************/

double *l1pp ( int n, double h )

/******************************************************************************/
/*
  Purpose:

    L1PP stores the 1D PP Laplacian as a full matrix.

  Discussion:

    The N grid points are assumed to be evenly spaced by H.

    For N = 5, the discrete Laplacian with periodic boundary conditions
    has the matrix form L:

       2 -1  0  0 -1
      -1  2 -1  0  0
       0 -1  2 -1  0
       0  0 -1  2 -1
      -1  0  0 -1  2

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 October 2013

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of points.
    N must be at least 3.

    Input, double H, the spacing between points.

    Output, double L1PP[N*N], the Laplacian matrix.
*/
{
  int i;
  int j;
  double *l;

  if ( n < 3 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "L1PP - Fatal error!\n" );
    fprintf ( stderr, "  N < 3.\n" );
    exit ( 1 );
  }

  l = ( double * ) malloc ( n * n * sizeof ( double ) );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      l[i+j*n] = 0.0;
    }
  }

  i = 0;
  l[i+i*n]       =  2.0 / h / h;
  l[i+(i+1)*n]   = -1.0 / h / h;
  l[i+(n-1)*n]   = -1.0 / h / h;

  for ( i = 1; i < n - 1; i++ )
  {
    l[i+(i-1)*n] = -1.0 / h / h;
    l[i+i*n] =      2.0 / h / h;
    l[i+(i+1)*n] = -1.0 / h / h;
  }

  i = n - 1;
  l[i+0*n] =     -1.0 / h / h;
  l[i+(i-1)*n] = -1.0 / h / h;
  l[i+i*n] =      2.0 / h / h;

  return l;
}
/******************************************************************************/

void l1pp_lu ( int n, double h, double l[], double u[] )

/******************************************************************************/
/*
  Discussion:

    L1PP_LU computes the LU factors of the 1D PP Laplacian.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 November 2013

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of points.
    N must be at least 3.

    Input, double H, the spacing between points.

    Output, double L[N*N], U[N*N], the LU factors.
*/
{
  int i;
  double i_r8;
  int j;

  if ( n < 3 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "L1PP_LU - Fatal error!\n" );
    fprintf ( stderr, "  N < 3.\n" );
    exit ( 1 );
  }

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      l[i+j*n] = 0.0;
    }
  }

  for ( i = 0; i < n; i++ )
  {
    l[i+i*n] = 1.0;
  }

  for ( i = 1; i < n - 1; i++ )
  {
    i_r8 = ( double ) ( i + 1 );
    l[i+(i-1)*n] =   - ( i_r8 - 1.0 ) / i_r8;
    l[n-1+(i-1)*n] =          - 1.0   / i_r8;
  }
  i = n - 1;
  l[i+(i-1)*n] = -1.0;

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      l[i+j*n] = l[i+j*n] / h;
    }
  }

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      u[i+j*n] = 0.0;
    }
  }

  for ( i = 0; i < n - 2; i++ )
  {
    i_r8 = ( double ) ( i + 1 );
    u[i+i*n] = ( i_r8 + 1.0 ) / i_r8;
    u[i+(i+1)*n] = - 1.0;
    u[i+(n-1)*n] =   - 1.0 / i_r8;
  }

  i = n - 2;
  i_r8 = ( double ) ( i + 1 );
  u[i+i*n] =       ( i_r8 + 1.0 ) / i_r8;
  u[i+(i+1)*n] = - ( i_r8 + 1.0 ) / i_r8;

  i = n - 1;
  u[i+i*n] = 0.0;

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      u[i+j*n] = u[i+j*n] / h;
    }
  }

  return;
}
/******************************************************************************/

double lu_error ( int n, double a[], double l[], double u[] )

/******************************************************************************/
/*
  Purpose:

    LU_ERROR determines the error in an LU factorization.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 November 2013

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of the matrix.

    Input, double A[N*N], the matrix.

    Input, double L[N*N], U[N*N], the LU factorization.

    Output, double LU_ERROR, the Frobenius norm
    of the difference matrix A - L * U.
*/
{
  double *d;
  double *lu;
  double value;

  lu = r8mat_mm_new ( n, n, n, l, u );

  d = r8mat_sub_new ( n, n, a, lu );
 
  value = r8mat_norm_fro ( n, n, d );

  free ( d );
  free ( lu );

  return value;
}
/******************************************************************************/

double *r8mat_mm_new ( int n1, int n2, int n3, double a[], double b[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_MM_NEW multiplies two matrices.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

    For this routine, the result is returned as the function value.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 April 2009

  Author:

    John Burkardt

  Parameters:

    Input, int N1, N2, N3, the order of the matrices.

    Input, double A[N1*N2], double B[N2*N3], the matrices to multiply.

    Output, double R8MAT_MM_NEW[N1*N3], the product matrix C = A * B.
*/
{
  double *c;
  int i;
  int j;
  int k;

  c = ( double * ) malloc ( n1 * n3 * sizeof ( double ) );

  for ( i = 0; i < n1; i++ )
  {
    for ( j = 0; j < n3; j++ )
    {
      c[i+j*n1] = 0.0;
      for ( k = 0; k < n2; k++ )
      {
        c[i+j*n1] = c[i+j*n1] + a[i+k*n1] * b[k+j*n2];
      }
    }
  }

  return c;
}
/******************************************************************************/

double *r8mat_mtm_new ( int n1, int n2, int n3, double a[], double b[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_MTM_NEW computes C = A' * B.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

    For this routine, the result is returned as the function value.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 September 2012

  Author:

    John Burkardt

  Parameters:

    Input, int N1, N2, N3, the order of the matrices.

    Input, double A[N2*N1], double B[N2*N3], the matrices to multiply.

    Output, double R8MAT_MTM_NEW[N1*N3], the product matrix C = A' * B.
*/
{
  double *c;
  int i;
  int j;
  int k;

  c = ( double * ) malloc ( n1 * n3 * sizeof ( double ) );

  for ( i = 0; i < n1; i++ )
  {
    for ( j = 0; j < n3; j++ )
    {
      c[i+j*n1] = 0.0;
      for ( k = 0; k < n2; k++ )
      {
        c[i+j*n1] = c[i+j*n1] + a[k+i*n2] * b[k+j*n2];
      }
    }
  }

  return c;
}
/******************************************************************************/

double r8mat_norm_fro ( int m, int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_NORM_FRO returns the Frobenius norm of an R8MAT.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

    The Frobenius norm is defined as

      R8MAT_NORM_FRO = sqrt (
        sum ( 1 <= I <= M ) sum ( 1 <= j <= N ) A(I,J)^2 )
    The matrix Frobenius norm is not derived from a vector norm, but
    is compatible with the vector L2 norm, so that:

      r8vec_norm_l2 ( A * x ) <= r8mat_norm_fro ( A ) * r8vec_norm_l2 ( x ).

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 July 2008

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows in A.

    Input, int N, the number of columns in A.

    Input, double A[M*N], the matrix whose Frobenius
    norm is desired.

    Output, double R8MAT_NORM_FRO, the Frobenius norm of A.
*/
{
  int i;
  int j;
  double value;

  value = 0.0;
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      value = value + pow ( a[i+j*m], 2 );
    }
  }
  value = sqrt ( value );

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
        fprintf ( stdout, "  %14g", a[i-1+(j-1)*m] );
      }
      fprintf ( stdout, "\n" );
    }
  }

  return;
# undef INCX
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

double *r8mat_sub_new ( int m, int n, double a[], double b[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_SUB_NEW computes C = A - B.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

    For this routine, the result is returned as the function value.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 October 2013

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the order of the matrices.

    Input, double A[M*N], double B[M*N], the matrices.

    Output, double R8MAT_SUB_NEW[M*N], the value of A-B.
*/
{
  double *c;
  int i;
  int j;

  c = ( double * ) malloc ( m * n * sizeof ( double ) );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      c[i+j*m] = a[i+j*m] - b[i+j*m];
    }
  }

  return c;
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

