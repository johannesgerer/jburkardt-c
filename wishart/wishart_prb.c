# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <string.h>

# include "wishart.h"
# include "pdflib.h"
# include "rnglib.h"

int main ( );
void wishart_test01 ( );
void wishart_test02 ( );
void wishart_test03 ( );
void wishart_test04 ( );
void wishart_test05 ( );
void wishart_test06 ( );
void wishart_test07 ( );
void wishart_test08 ( );
void wishart_test09 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for WISHART_PRB.

  Discussion:

    WISHART_PRB tests the WISHART library.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 October 2013

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "WISHART_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the WISHART library.\n" );

  wishart_test01 ( );
  wishart_test02 ( );
  wishart_test03 ( );
  wishart_test04 ( );
  wishart_test05 ( );
  wishart_test06 ( );
  wishart_test07 ( );
  wishart_test08 ( );
  wishart_test09 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "WISHART_PRB\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void wishart_test01 ( )

/******************************************************************************/
/*
  Purpose:

    WISHART_TEST01 demonstrates the unit Wishart sampling function.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 July 2013

  Author:

    John Burkardt
*/
{
  int df;
  int it_max;
  int it_num;
  double *lambda;
  int n;
  int rot_num;
  double *v;
  double *w;
/*
  Initialize the RNGLIB library.
*/
  initialize ( );

  printf ( "\n" );
  printf ( "WISHART_TEST01:\n" );
  printf ( "  We can compute sample unit Wishart matrices by:\n" );
  printf ( "    W = wishart_unit_sample ( n, df );\n" );
/*
  Set the parameters and call.
*/
  n = 5;
  df = 8;
  w = wishart_unit_sample ( n, df );
  r8mat_print ( n, n, w, "  wishart_unit_sample ( 5, 8 ):" );
  free ( w );
/*
  Calling again yields a new matrix.
*/
  w = wishart_unit_sample ( n, df );
  r8mat_print ( n, n, w, "  wishart_unit_sample ( 5, 8 ):" );
  free ( w );
/*
  Reduce DF
*/
  n = 5;
  df = 5;
  w = wishart_unit_sample ( n, df );
  r8mat_print ( n, n, w, "  wishart_unit_sample ( 5, 5 ):" );
  free ( w );
/*
  Try a smaller matrix.
*/
  n = 3;
  df = 5;
  w = wishart_unit_sample ( n, df );
  r8mat_print ( n, n, w, "  wishart_unit_sample ( 3, 5 ):" );
/*
  What is the eigendecomposition of the matrix?
*/
  it_max = 50;
  v = ( double * ) malloc ( n * n * sizeof ( double ) );
  lambda = ( double * ) malloc ( n * sizeof ( double ) );

  jacobi_eigenvalue ( n, w, it_max, v, lambda, &it_num, &rot_num );
  r8mat_print ( n, n, v, "  Eigenvectors of previous matrix:" );
  r8vec_print ( n, lambda, "  Eigenvalues of previous matrix:" );
/*
  Free memory.
*/
  free ( lambda );
  free ( v );
  free ( w );

  return;
}
/******************************************************************************/

void wishart_test02 ( )

/******************************************************************************/
/*
  Purpose:

    WISHART_TEST02 demonstrates the unit Bartlett sampling function.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 July 2013

  Author:

    John Burkardt
*/
{
  int df;
  int it_max;
  int it_num;
  double *lambda;
  int n;
  int rot_num;
  double *t;
  double *v;
  double *w;
/*
   Initialize the RNGLIB library.
*/
  initialize ( );

  printf ( "\n" );
  printf ( "WISHART_TEST02:\n" );
  printf ( "  We can compute sample unit Bartlett matrices by:\n" );
  printf ( "    T = bartlett_unit_sample ( n, df );\n" );
/*
   Set the parameters and call.
*/
  n = 5;
  df = 8;
  t = bartlett_unit_sample ( n, df );
  r8mat_print ( n, n, t, "  bartlett_unit_sample ( 5, 8 ):" );
  free ( t );
/*
   Calling again yields a new matrix.
*/
  t = bartlett_unit_sample ( n, df );
  r8mat_print ( n, n, t, "  bartlett_unit_sample ( 5, 8 ):" );
  free ( t );
/*
   Reduce DF.
*/
  n = 5;
  df = 5;
  t = bartlett_unit_sample ( n, df );
  r8mat_print ( n, n, t, "  bartlett_unit_sample ( 5, 5 ):" );
  free ( t );
/*
   Try a smaller matrix.
*/
  n = 3;
  df = 5;
  t = bartlett_unit_sample ( n, df );
  r8mat_print ( n, n, t, "  bartlett_unit_sample ( 3, 5 ):" );
/*
   What is the eigendecomposition of the matrix T' * T?
*/
  w = r8mat_mtm_new ( n, n, n, t, t );

  it_max = 50;
  v = ( double * ) malloc ( n * n * sizeof ( double ) );
  lambda = ( double * ) malloc ( n * sizeof ( double ) );

  jacobi_eigenvalue ( n, w, it_max, v, lambda, &it_num, &rot_num );
  r8mat_print ( n, n, v, "  Eigenvectors of previous matrix:" );
  r8vec_print ( n, lambda, "  Eigenvalues of previous matrix:" );
/*
  Free memory.
*/
  free ( lambda );
  free ( t );
  free ( v );
  free ( w );

  return;
}
/******************************************************************************/

void wishart_test03 ( )

/******************************************************************************/
/*
  Purpose:

    WISHART_TEST03 compares the unit Wishart and Bartlett sample matrices.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 July 2013

  Author:

    John Burkardt
*/
{
  int df;
  double diff;
  int n;
  double *t;
  double *tt;
  double *w;
/*
   Initialize the RNGLIB library.
   Normally, we would do this just once, here at the beginning.
   In this example, however, we really want to do it just before
   we call each of the sampling routines, so that they both access
   the same set of random numbers...
*/
  initialize ( );

  printf ( "\n" );
  printf ( "WISHART_TEST03:\n" );
  printf ( "  Verify that, if using the same set of random numbers,\n" );
  printf ( "    W = T' * T,\n" );
  printf ( "  where\n" );
  printf ( "    W = wishart_unit_sample ( n, df );\n" );
  printf ( "    T = bartlett_unit_sample ( n, df );\n" );
/*
   Set the parameters.
*/
  n = 5;
  df = 8;
/*
   Initialize the random number package and compute W.
*/
  initialize ( );
  w = wishart_unit_sample ( n, df );
/*
   Initialize the random number package again, and compute T.
*/
  initialize ( );
  t = bartlett_unit_sample ( n, df );
/*
   Compute T' * T.
*/
  tt = r8mat_mtm_new ( n, n, n, t, t );
/*
   Compare T'T to W.
*/
  diff = r8mat_norm_fro_affine ( n, n, w, tt );
  printf ( "\n" );
  printf ( "  Frobenius norm of error is %g\n", diff );
/*
  Free memory.
*/
  free ( t );
  free ( tt );
  free ( w );

  return;
}
/******************************************************************************/

void wishart_test04 ( )

/******************************************************************************/
/*
  Purpose:

    WISHART_TEST04 demonstrates the Wishart sampling function.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 July 2013

  Author:

    John Burkardt
*/
{
  int df;
  int i;
  int it_max;
  int it_num;
  int j;
  double *lambda;
  int n;
  int rot_num;
/*
  Note that R is an upper triangular matrix,
  whose entries here are listed in column major order.
*/
  double r[3*3] = { 
    5.0, 0.0, 0.0,
    1.0, 4.0, 0.0,
    3.0, 2.0, 6.0 };
  double *sigma;
  double sigma_diag[5] = { 1.0, 2.0, 3.0, 4.0, 5.0 };
  double *v;
  double *w;
/*
   Initialize the RNGLIB library.
*/
  initialize ( );

  printf ( "\n" );
  printf ( "WISHART_TEST04:\n" );
  printf ( "  We can compute sample Wishart matrices by:\n" );
  printf ( "    W = wishart_sample ( n, df, sigma );\n" );
/*
   Set the parameters and call.
*/
  n = 5;
  df = 8;
  sigma = r8mat_identity_new ( n );
  w = wishart_sample ( n, df, sigma );
  r8mat_print ( n, n, w, "  wishart_sample ( 5, 8, Identity ):" );
  free ( w );
/*
   Calling again yields a new matrix.
*/
  w = wishart_sample ( n, df, sigma );
  r8mat_print ( n, n, w, "  wishart_sample ( 5, 8, Identity ):" );
  free ( sigma );
  free ( w );
/*
   Try a diagonal matrix.
*/
  sigma = r8mat_diagonal_new ( n, sigma_diag );
  w = wishart_sample ( n, df, sigma );
  r8mat_print ( n, n, w, "  wishart_sample ( 5, 8, diag(1,2,3,4,5) ):" );
  free ( sigma );
  free ( w );
/*
   Try a smaller matrix.  Sigma must be positive definite symmetric.
*/
  n = 3;
  df = 3;
  sigma = r8mat_mtm_new ( n, n, n, r, r );
  r8mat_print ( n, n, sigma, "  Set covariance SIGMA:" );
  w = wishart_sample ( n, df, sigma );
  r8mat_print ( n, n, w, "  wishart_sample ( 3, 3, sigma ):" );
/*
   What is the eigendecomposition of this matrix?
*/
  it_max = 50;
  v = ( double * ) malloc ( n * n * sizeof ( double ) );
  lambda = ( double * ) malloc ( n * sizeof ( double ) );

  jacobi_eigenvalue ( n, w, it_max, v, lambda, &it_num, &rot_num );
  r8mat_print ( n, n, v, "  Eigenvectors of previous matrix:" );
  r8vec_print ( n, lambda, "  Eigenvalues of previous matrix:" );
/*
  Free memory.
*/
  free ( lambda );
  free ( sigma );
  free ( v );
  free ( w );

  return;
}
/******************************************************************************/

void wishart_test05 ( )

/******************************************************************************/
/*
  Purpose:

    WISHART_TEST05 demonstrates the Bartlett sampling function.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 July 2013

  Author:

    John Burkardt
*/
{
  int df;
  int it_max;
  int it_num;
  double *lambda;
  int n;
/*
  Note that R is an upper triangular matrix,
  whose entries here are listed in column major order.
*/
  double r[3*3] = { 
    5.0, 0.0, 0.0,
    1.0, 4.0, 0.0,
    3.0, 2.0, 6.0 };
  int rot_num;
  double *sigma;
  double sigma_diag[5] = { 1.0, 2.0, 3.0, 4.0, 5.0 };
  double *t;
  double *v;
  double *w;
/*
   Initialize the RNGLIB library.
*/
  initialize ( );

  printf ( "\n" );
  printf ( "WISHART_TEST05:\n" );
  printf ( "  We can compute sample Bartlett matrices by:\n" );
  printf ( "    T = bartlett_sample ( n, df, sigma );\n" );
/*
   Set the parameters and call.
*/
  n = 5;
  df = 8;
  sigma = r8mat_identity_new ( n );
  t = bartlett_sample ( n, df, sigma );
  r8mat_print ( n, n, t, "  bartlett_sample ( 5, 8, Identity ):" );
  free ( t );
/*
   Calling again yields a new matrix.
*/
  t = bartlett_sample ( n, df, sigma );
  r8mat_print ( n, n, t, "  bartlett_sample ( 5, 8, Identity ):" );
  free ( sigma );
  free ( t );
/*
   Try a diagonal matrix.
*/
  sigma = r8mat_diagonal_new ( n, sigma_diag );
  t = bartlett_sample ( n, df, sigma );
  r8mat_print ( n, n, t, "  bartlett_sample ( 5, 8, diag(1,2,3,4,5) ):" );
  free ( sigma );
  free ( t );
/*
   Try a smaller matrix.
*/
  n = 3;
  df = 3;
  sigma = r8mat_mtm_new ( n, n, n, r, r );
  r8mat_print ( n, n, sigma, "  Set covariance SIGMA:" );
  t = bartlett_sample ( n, df, sigma );
  r8mat_print ( n, n, t, "  bartlett_sample ( 3, 3, sigma ):" );
/*
   What is the eigendecomposition of T' * T?
*/
  w = r8mat_mtm_new ( n, n, n, t, t );
  it_max = 50;
  v = ( double * ) malloc ( n * n * sizeof ( double ) );
  lambda = ( double * ) malloc ( n * sizeof ( double ) );

  jacobi_eigenvalue ( n, w, it_max, v, lambda, &it_num, &rot_num );
  r8mat_print ( n, n, v, "  Eigenvectors of previous matrix:" );
  r8vec_print ( n, lambda, "  Eigenvalues of previous matrix:" );
/*
  Free memory.
*/
  free ( lambda );
  free ( sigma );
  free ( t );
  free ( v );
  free ( w );

  return;
}
/******************************************************************************/

void wishart_test06 ( )

/******************************************************************************/
/*
  Purpose:

    WISHART_TEST06 compares the Wishart and Bartlett sample matrices.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 July 2013

  Author:

    John Burkardt
*/
{
  int df;
  double diff;
  int n;
/*
  Note that R is an upper triangular matrix,
  whose entries here are listed in column major order.
*/
  double r[3*3] = { 
    5.0, 0.0, 0.0,
    1.0, 4.0, 0.0,
    3.0, 2.0, 6.0 };
  double *sigma;
  double *t;
  double *tt;
  double *w;
/*
   Initialize the RNGLIB library.
   Normally, we would do this just once, here at the beginning.
   In this example, however, we really want to do it just before
   we call each of the sampling routines, so that they both access
   the same set of random numbers...
*/
  initialize ( );

  printf ( "\n" );
  printf ( "WISHART_TEST06:\n" );
  printf ( "  Verify that, if using the same set of random numbers,\n" );
  printf ( "    W = T'' * T,\n" );
  printf ( "  where\n" );
  printf ( "    W = wishart_sample ( n, df, sigma );\n" );
  printf ( "    T = bartlett_sample ( n, df, sigma );\n" );
/*
   Set the parameters.
*/
  n = 3;
  df = 5;
  sigma = r8mat_mtm_new ( n, n, n, r, r );
  r8mat_print ( n, n, sigma, "  Covariance SIGMA:" );
/*
   Initialize the random number package and compute W.
*/
  initialize ( );
  w = wishart_sample ( n, df, sigma );
/*
   Initialize the random number package again, and compute T.
*/
  initialize ( );
  t = bartlett_sample ( n, df, sigma );
/*
   Compute T' * T.
*/
  tt = r8mat_mtm_new ( n, n, n, t, t );
/*
   Compare T'T to W.
*/
  diff = r8mat_norm_fro_affine ( n, n, w, tt );
  printf ( "\n" );
  printf ( "  Frobenius norm of error is %g\n", diff );
/*
  Free memory.
*/
  free ( sigma );
  free ( t );
  free ( tt );
  free ( w );

  return;
}
/******************************************************************************/

void wishart_test07 ( )

/******************************************************************************/
/*
  Purpose:

    WISHART_TEST07 demonstrates a property of the Wishart distribution.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 July 2013

  Author:

    John Burkardt
*/
{
  int df;
  double diff;
  double divisor;
  int i;
  int n;
/*
  Note that R is an upper triangular matrix,
  whose entries here are listed in column major order.
*/
  double r[3*3] = { 
    5.0, 0.0, 0.0,
    1.0, 4.0, 0.0,
    3.0, 2.0, 6.0 };
  int sample_num;
  double *sigma;
  double *w;
  double *w_average;
/*
   Initialize the RNGLIB library.
*/
  initialize ( );

  printf ( "\n" );
  printf ( "WISHART_TEST07:\n" );
  printf ( "  For given values of N, DF, SIGMA, the random\n" );
  printf ( "  matrices from the Wishart distribution:\n" );
  printf ( "    W = wishart_sample ( n, df, sigma );\n" );
  printf ( "  should have mean DF * SIGMA.\n" );
/*
   Set the parameters.
*/
  n = 3;
  printf ( "  Fix N = %d\n", n );
  df = 5;
  printf ( "  Fix DF = %d\n", df );
  sigma = r8mat_mtm_new ( n, n, n, r, r );
  r8mat_print ( n, n, sigma, "  Fix covariance SIGMA:" );
/*
   Sample many times and average.
*/
  sample_num = 1000;
  w_average = r8mat_zero_new ( n, n );
  for ( i = 1; i <= sample_num; i++ )
  {
    w = wishart_sample ( n, df, sigma );
    r8mat_add ( n, n, w, w_average );
    free ( w );
  }
  divisor = ( double ) sample_num;
  r8mat_divide ( n, n, divisor, w_average );
/*
   Compare SIGMA and W_SAMPLE / DF.
*/
  divisor = ( double ) df;
  r8mat_divide ( n, n, divisor, w_average );

  r8mat_print ( n, n, w_average, "  W_Average / DF: " );

  diff = r8mat_norm_fro_affine ( n, n, sigma, w_average );
  printf ( "\n" );
  printf ( "  Frobenius norm of SIGMA-W_average/DF = %g\n", diff );
/*
  Free memory.
*/
  free ( sigma );
  free ( w_average );

  return;
}
/******************************************************************************/

void wishart_test08 ( )

/******************************************************************************/
/*
  Purpose:

    WISHART_TEST08 samples the unit Wishart and unit Wishart inverse matrices.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 October 2013

  Author:

    John Burkardt
*/
{
  int df;
  double diff;
  int i;
  double *ident;
  int j;
  double *m;
  int n;
  double *w;
  double *wm;
/*
   Initialize the RNGLIB library.
   Normally, we would do this just once, here at the beginning.
   In this example, however, we really want to do it just before
   we call each of the sampling routines, so that they both access
   the same set of random numbers...
*/
  initialize ( );

  printf ( "\n" );
  printf ( "WISHART_TEST08:\n" );
  printf ( "  Verify that, if using the same set of random numbers,\n" );
  printf ( "    inverse(W) = M,\n" );
  printf ( "  where\n" );
  printf ( "    W = wishart_unit_sample ( n, df );\n" );
  printf ( "    M = wishart_unit_sample_inverse ( n, df );\n" );
/*
   Set the parameters.
*/
  n = 5;
  df = 8;
/*
   Initialize the random number package and compute W.
*/
  initialize ( );
  w = wishart_unit_sample ( n, df );
/*
   Initialize the random number package again, and compute M.
*/
  initialize ( );
  m = wishart_unit_sample_inverse ( n, df );
/*
   Compute W * M.
*/
  wm = r8mat_mm_new ( n, n, n, w, m );
/*
   Compare W*M to I.
*/
  ident = r8mat_identity_new ( n );
  diff = r8mat_norm_fro_affine ( n, n, wm, ident );
  printf ( "\n" );
  printf ( "  Frobenius norm of error is %g\n", diff );
/*
  Free memory.
*/
  free ( ident );
  free ( m );
  free ( w );
  free ( wm );

  return;
}
/******************************************************************************/

void wishart_test09 ( )

/******************************************************************************/
/*
  Purpose:

    WISHART_TEST09 samples the Wishart and Wishart inverse matrices.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 October 2013

  Author:

    John Burkardt
*/
{
  int df;
  double diff;
  int i;
  double *ident;
  int j;
  double *m;
  int n;
/*
  Note that R is an upper triangular matrix,
  whose entries here are listed in column major order.
*/
  double r[5*5] = { 
    3.0, 0.0, 0.0, 0.0, 0.0, 
    1.0, 7.0, 0.0, 0.0, 0.0, 
    1.0, 1.0, 5.0, 0.0, 0.0, 
    1.0, 2.0, 1.0, 4.0, 0.0, 
    1.0, 3.0, 3.0, 2.0, 6.0 };
  double *sigma;
  double *w;
  double *wm;
/*
   Initialize the RNGLIB library.
   Normally, we would do this just once, here at the beginning.
   In this example, however, we really want to do it just before
   we call each of the sampling routines, so that they both access
   the same set of random numbers...
*/
  initialize ( );

  printf ( "\n" );
  printf ( "WISHART_TEST09:\n" );
  printf ( "  Verify that, if using the same set of random numbers,\n" );
  printf ( "    inverse(W) = M,\n" );
  printf ( "  where\n" );
  printf ( "    W = wishart_sample ( n, df, sigma );\n" );
  printf ( "    M = wishart_sample_inverse ( n, df, sigma );\n" );
/*
   Set the parameters.
*/
  n = 5;
  df = 8;
  sigma = r8mat_mtm_new ( n, n, n, r, r );
/*
   Initialize the random number package and compute W.
*/
  initialize ( );
  w = wishart_sample ( n, df, sigma );
/*
   Initialize the random number package again, and compute M.
*/
  initialize ( );
  m = wishart_sample_inverse ( n, df, sigma );
/*
   Compute W * M.
*/
  wm = r8mat_mm_new ( n, n, n, w, m );
/*
   Compare W*M to I.
*/
  ident = r8mat_identity_new ( n );
  diff = r8mat_norm_fro_affine ( n, n, wm, ident );
  printf ( "\n" );
  printf ( "  Frobenius norm of error is %g\n", diff );
/*
  Free memory.
*/
  free ( ident );
  free ( m );
  free ( sigma );
  free ( w );
  free ( wm );

  return;
}
