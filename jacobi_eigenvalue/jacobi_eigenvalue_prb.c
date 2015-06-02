# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <string.h>

# include "jacobi_eigenvalue.h"

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for JACOBI_EIGENVALUE_PRB.

  Discussion:

    JACOBI_EIGENVALUE_PRB tests the JACOBI_EIGENVALUE library.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    15 July 2013

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "JACOBI_EIGENVALUE_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the JACOBI_EIGENVALUE library.\n" );

  test01 ( );
  test02 ( );
  test03 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "JACOBI_EIGENVALUE_PRB\n" );
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

    TEST01 uses a 4x4 test matrix.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    15 July 2013

  Author:

    John Burkardt
*/
{
# define N 4

  double a[N*N] = {
      4.0,  -30.0,    60.0,   -35.0, 
    -30.0,  300.0,  -675.0,   420.0, 
     60.0, -675.0,  1620.0, -1050.0, 
    -35.0,  420.0, -1050.0,   700.0 };
  double d[N];
  double error_frobenius;
  int it_max;
  int it_num;
  int n = N;
  int rot_num;
  double v[N*N];

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  For a symmetric matrix A,\n" );
  printf ( "  JACOBI_EIGENVALUE computes the eigenvalues D\n" );
  printf ( "  and eigenvectors V so that A * V = D * V.\n" );

  r8mat_print ( n, n, a, "  Input matrix A:" );

  it_max = 100;

  jacobi_eigenvalue ( n, a, it_max, v, d, &it_num, &rot_num );

  printf ( "\n" );
  printf ( "  Number of iterations = %d\n", it_num );
  printf ( "  Number of rotations  = %d\n", rot_num );

  r8vec_print ( n, d, "  Eigenvalues D:" );

  r8mat_print ( n, n, v, "  Eigenvector matrix V:" );
/*
  Compute eigentest.
*/
  error_frobenius = r8mat_is_eigen_right ( n, n, a, v, d );
  printf ( "\n" );
  printf ( "  Frobenius norm error in eigensystem A*V-D*V = %g\n",
    error_frobenius );

  return;
# undef N
}
/******************************************************************************/

void test02 ( )

/******************************************************************************/
/*
  Purpose:

    TEST02 uses a 4x4 test matrix.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    15 July 2013

  Author:

    John Burkardt
*/
{
# define N 4

  double a[N*N] = {
    4.0, 0.0, 0.0, 0.0, 
    0.0, 1.0, 0.0, 0.0, 
    0.0, 0.0, 3.0, 0.0, 
    0.0, 0.0, 0.0, 2.0 };
  double d[N];
  double error_frobenius;
  int it_max;
  int it_num;
  int n = N;
  int rot_num;
  double v[N*N];

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  For a symmetric matrix A,\n" );
  printf ( "  JACOBI_EIGENVALUE computes the eigenvalues D\n" );
  printf ( "  and eigenvectors V so that A * V = D * V.\n" );
  printf ( "\n" );
  printf ( "As a sanity check, input a diagonal matrix.\n" );

  r8mat_print ( n, n, a, "  Input matrix A:" );

  it_max = 100;

  jacobi_eigenvalue ( n, a, it_max, v, d, &it_num, &rot_num );

  printf ( "\n" );
  printf ( "  Number of iterations = %d\n", it_num );
  printf ( "  Number of rotations  = %d\n", rot_num );

  r8vec_print ( n, d, "  Eigenvalues D:" );

  r8mat_print ( n, n, v, "  Eigenvector matrix V:" );
/*
  Compute eigentest.
*/
  error_frobenius = r8mat_is_eigen_right ( n, n, a, v, d );
  printf ( "\n" );
  printf ( "  Frobenius norm error in eigensystem A*V-D*V = %g\n",
    error_frobenius );

  return;
# undef N
}
/******************************************************************************/

void test03 ( )

/******************************************************************************/
/*
  Purpose:

    TEST03 uses a 5x5 test matrix.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    15 July 2013

  Author:

    John Burkardt
*/
{
# define N 5

  double a[N*N];
  double d[N];
  double error_frobenius;
  int i;
  int it_max;
  int it_num;
  int j;
  int n = N;
  int rot_num;
  double v[N*N];

  printf ( "\n" );
  printf ( "TEST03\n" );
  printf ( "  For a symmetric matrix A,\n" );
  printf ( "  JACOBI_EIGENVALUE computes the eigenvalues D\n" );
  printf ( "  and eigenvectors V so that A * V = D * V.\n" );
  printf ( "\n" );
  printf ( "  Use the discretized second derivative matrix.\n" );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      if ( i == j )
      {
        a[i+j*n] = -2.0;
      }
      else if ( i == j + 1 || i == j - 1 )
      {
        a[i+j*n] = 1.0;
      }
      else
      {
        a[i+j*n] = 0.0;
      }
    }
  }

  r8mat_print ( n, n, a, "  Input matrix A:" );

  it_max = 100;

  jacobi_eigenvalue ( n, a, it_max, v, d, &it_num, &rot_num );

  printf ( "\n" );
  printf ( "  Number of iterations = %d\n", it_num );
  printf ( "  Number of rotations  = %d\n", rot_num );

  r8vec_print ( n, d, "  Eigenvalues D:" );

  r8mat_print ( n, n, v, "  Eigenvector matrix V:" );
/*
  Compute eigentest.
*/
  error_frobenius = r8mat_is_eigen_right ( n, n, a, v, d );
  printf ( "\n" );
  printf ( "  Frobenius norm error in eigensystem A*V-D*V = %g\n",
    error_frobenius );

  return;
# undef N
}
