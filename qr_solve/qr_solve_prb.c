# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "qr_solve.h"
# include "test_ls.h"
# include "r8lib.h"

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );
void test04 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for QR_SOLVE_PRB.

  Discussion:

    QR_SOLVE_PRB tests the QR_SOLVE library.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    04 October 2012

  Author:

    John Burkardt
*/
{
  timestamp ( );

  printf ( "\n" );
  printf ( "QR_SOLVE_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the QR_SOLVE library.\n" );
  printf ( "  QR_SOLVE needs the R8LIB library.\n" );
  printf ( "  This test also needs the TEST_LS library.\n" );
 
  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "QR_SOLVE_PRB\n" );
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

    TEST01 tests NORMAL_SOLVE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 September 2012

  Author:

    John Burkardt
*/
{
  double *a;
  double *b;
  double b_norm;
  int flag;
  int i;
  int m;
  int n;
  int prob;
  int prob_num;
  double *r1;
  double r1_norm;
  double *r2;
  double r2_norm;
  double x_diff_norm;
  double *x1;
  double x1_norm;
  double *x2;
  double x2_norm;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  NORMAL_SOLVE is a function with a simple interface which\n" );
  printf ( "  solves a linear system A*x = b in the least squares sense.\n" );
  printf ( "  Compare a tabulated solution X1 to the NORMAL_SOLVE result X2.\n" );
  printf ( "\n" );
  printf ( "  NORMAL_SOLVE cannot be applied when N < M,\n" );
  printf ( "  or if the matrix does not have full column rank.\n" );

  prob_num = p00_prob_num ( );

  printf ( "\n" );
  printf ( "  Number of problems = %d\n", prob_num );
  printf ( "\n" );
  printf ( "  Index     M     N     ||B||         ||X1 - X2||   ||X1||       ||X2||        ||R1||        ||R2||\n" );
  printf ( "\n" );

  for ( prob = 1; prob <= prob_num; prob++ )
  {
/*
  Get problem size.
*/
    m = p00_m ( prob );
    n = p00_n ( prob );
/*
  Retrieve problem data.
*/
    a = p00_a ( prob, m, n );
    b = p00_b ( prob, m );
    x1 = p00_x ( prob, n );

    b_norm = r8vec_norm ( m, b );
    x1_norm = r8vec_norm ( n, x1 );
    r1 = r8mat_mv_new ( m, n, a, x1 );
    for ( i = 0; i < m; i++ )
    {
      r1[i] = r1[i] - b[i];
    }
    r1_norm = r8vec_norm ( m, r1 );
/*
  Use NORMAL_SOLVE on the problem.
*/
    x2 = normal_solve ( m, n, a, b, &flag );

    if ( flag != 0 )
    {
      printf ( "  %5d  %4d  %4d  %12g  ------------  %12g   ------------  %12g  ------------\n",
      prob, m, n, b_norm, x1_norm, r1_norm );
    }
    else
    {
      x2_norm = r8vec_norm ( n, x2 );
      r2 = r8mat_mv_new ( m, n, a, x2 );
      for ( i = 0; i < m; i++ )
      {
        r2[i] = r2[i] - b[i];
      }
      r2_norm = r8vec_norm ( m, r2 );
/*
  Compare tabulated and computed solutions.
*/
      x_diff_norm = r8vec_norm_affine ( n, x1, x2 );
/*
  Report results for this problem.
*/
      printf ( "  %5d  %4d  %4d  %12g  %12g  %12g  %12g  %12g  %12g\n",
        prob, m, n, b_norm, x_diff_norm, x1_norm, x2_norm, r1_norm, r2_norm );

      free ( r2 );
      free ( x2 );
    }
/*
  Deallocate memory.
*/
    free ( a );
    free ( b );
    free ( r1 );
    free ( x1 );
  }
  return;
}
/******************************************************************************/

void test02 ( )

/******************************************************************************/
/*
  Purpose:

    TEST02 tests QR_SOLVE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 September 2012

  Author:

    John Burkardt
*/
{
  double *a;
  double *b;
  double b_norm;
  int i;
  int m;
  int n;
  int prob;
  int prob_num;
  double *r1;
  double r1_norm;
  double *r2;
  double r2_norm;
  double x_diff_norm;
  double *x1;
  double x1_norm;
  double *x2;
  double x2_norm;

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  QR_SOLVE is a function with a simple interface which\n" );
  printf ( "  solves a linear system A*x = b in the least squares sense.\n" );
  printf ( "  Compare a tabulated solution X1 to the QR_SOLVE result X2.\n" );

  prob_num = p00_prob_num ( );

  printf ( "\n" );
  printf ( "  Number of problems = %d\n", prob_num );
  printf ( "\n" );
  printf ( "  Index     M     N     ||B||         ||X1 - X2||   ||X1||       ||X2||        ||R1||        ||R2||\n" );
  printf ( "\n" );

  for ( prob = 1; prob <= prob_num; prob++ )
  {
/*
  Get problem size.
*/
    m = p00_m ( prob );
    n = p00_n ( prob );
/*
  Retrieve problem data.
*/
    a = p00_a ( prob, m, n );
    b = p00_b ( prob, m );
    x1 = p00_x ( prob, n );

    b_norm = r8vec_norm ( m, b );
    x1_norm = r8vec_norm ( n, x1 );
    r1 = r8mat_mv_new ( m, n, a, x1 );
    for ( i = 0; i < m; i++ )
    {
      r1[i] = r1[i] - b[i];
    }
    r1_norm = r8vec_norm ( m, r1 );
/*
  Use QR_SOLVE on the problem.
*/
    x2 = qr_solve ( m, n, a, b );

    x2_norm = r8vec_norm ( n, x2 );
    r2 = r8mat_mv_new ( m, n, a, x2 );
    for ( i = 0; i < m; i++ )
    {
      r2[i] = r2[i] - b[i];
    }
    r2_norm = r8vec_norm ( m, r2 );
/*
  Compare tabulated and computed solutions.
*/
    x_diff_norm = r8vec_norm_affine ( n, x1, x2 );
/*
  Report results for this problem.
*/
    printf ( "  %5d  %4d  %4d  %12g  %12g  %12g  %12g  %12g  %12g\n",
      prob, m, n, b_norm, x_diff_norm, x1_norm, x2_norm, r1_norm, r2_norm );
/*
  Deallocate memory.
*/
    free ( a );
    free ( b );
    free ( r1 );
    free ( r2 );
    free ( x1 );
    free ( x2 );
  }
  return;
}
/******************************************************************************/

void test03 ( )

/******************************************************************************/
/*
  Purpose:

    TEST03 tests SVD_SOLVE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 September 2012

  Author:

    John Burkardt
*/
{
  double *a;
  double *b;
  double b_norm;
  int i;
  int m;
  int n;
  int prob;
  int prob_num;
  double *r1;
  double r1_norm;
  double *r2;
  double r2_norm;
  double x_diff_norm;
  double *x1;
  double x1_norm;
  double *x2;
  double x2_norm;

  printf ( "\n" );
  printf ( "TEST03\n" );
  printf ( "  SVD_SOLVE is a function with a simple interface which\n" );
  printf ( "  solves a linear system A*x = b in the least squares sense.\n" );
  printf ( "  Compare a tabulated solution X1 to the SVD_SOLVE result X2.\n" );

  prob_num = p00_prob_num ( );

  printf ( "\n" );
  printf ( "  Number of problems = %d\n", prob_num );
  printf ( "\n" );
  printf ( "  Index     M     N     ||B||         ||X1 - X2||   ||X1||       ||X2||        ||R1||        ||R2||\n" );
  printf ( "\n" );

  for ( prob = 1; prob <= prob_num; prob++ )
  {
/*
  Get problem size.
*/
    m = p00_m ( prob );
    n = p00_n ( prob );
/*
  Retrieve problem data.
*/
    a = p00_a ( prob, m, n );
    b = p00_b ( prob, m );
    x1 = p00_x ( prob, n );

    b_norm = r8vec_norm ( m, b );
    x1_norm = r8vec_norm ( n, x1 );
    r1 = r8mat_mv_new ( m, n, a, x1 );
    for ( i = 0; i < m; i++ )
    {
      r1[i] = r1[i] - b[i];
    }
    r1_norm = r8vec_norm ( m, r1 );
/*
  Use SVD_SOLVE on the problem.
*/
    x2 = svd_solve ( m, n, a, b );

    x2_norm = r8vec_norm ( n, x2 );
    r2 = r8mat_mv_new ( m, n, a, x2 );
    for ( i = 0; i < m; i++ )
    {
      r2[i] = r2[i] - b[i];
    }
    r2_norm = r8vec_norm ( m, r2 );
/*
  Compare tabulated and computed solutions.
*/
    x_diff_norm = r8vec_norm_affine ( n, x1, x2 );
/*
  Report results for this problem.
*/
    printf ( "  %5d  %4d  %4d  %12g  %12g  %12g  %12g  %12g  %12g\n",
      prob, m, n, b_norm, x_diff_norm, x1_norm, x2_norm, r1_norm, r2_norm );
/*
  Deallocate memory.
*/
    free ( a );
    free ( b );
    free ( r1 );
    free ( r2 );
    free ( x1 );
    free ( x2 );
  }
  return;
}
/******************************************************************************/

void test04 ( )

/******************************************************************************/
/*
  Purpose:

    TEST04 tests DQRLS.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 September 2012

  Author:

    John Burkardt
*/
{
  double *a;
  double b[5] = { 1.0, 2.3, 4.6, 3.1, 1.2 };
  int i;
  int ind;
  int itask;
  int j;
  int *jpvt;
  int kr;
  int m = 5;
  int n = 3;
  double *qraux;
  double tol;
  double *x;

  a = ( double * ) malloc ( m * n * sizeof ( double ) );
  jpvt = ( int * ) malloc ( n * sizeof ( int ) );
  qraux = ( double * ) malloc ( n * sizeof ( double ) );
  x = ( double * ) malloc ( n * sizeof ( double ) );
/*
  Set up least-squares problem
  quadratic model, equally-spaced points
*/
  printf ( "\n" );
  printf ( "TEST04\n" );
  printf ( "  DQRLS solves a linear system A*x = b in the least squares sense.\n" );

  for ( i = 0; i < m; i++ )
  {
    a[i+0*m] = 1.0;
    for ( j = 1; j < n; j++ )
    {
      a[i+j*m] = a[i+(j-1)*m] * ( double ) ( i + 1 );
    }
  }

  tol = 1.0E-06;

  r8mat_print ( m, n, a, "  Coefficient matrix A:" );

  r8vec_print ( m, b, "  Right hand side b:" );
/*
  Solve least-squares problem
*/
  itask = 1;
  ind = dqrls ( a, m, m, n, tol, &kr, b, x, b, jpvt, qraux, itask );
/*
  Print results
*/
  printf ( "\n" );
  printf ( "  Error code = %d\n", ind );
  printf ( "  Estimated matrix rank = %d\n", kr );

  r8vec_print ( n, x, "  Least squares solution x:" );

  r8vec_print ( m, b, "  Residuals A*x-b" );

  free ( a );
  free ( jpvt );
  free ( qraux );
  free ( x );

  return;
}
