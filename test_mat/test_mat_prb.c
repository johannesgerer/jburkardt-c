# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <complex.h>
# include <string.h>

# include "test_mat.h"

int main ( void );
void test_cond ( void );
void test_determinant ( void );
void test_eigen ( void );
void test_inverse ( void );
void test_null ( void );
void test_plu ( void );
void test_solution ( void );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for TEST_MAT_PRB.

  Discussion:

    TEST_MAT_PRB calls the TEST_MAT test routines.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    07 April 2012

  Author:

    John Burkardt
*/
{
  timestamp ( );

  printf ( "\n" );
  printf ( "TEST_MAT_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the TEST_MAT library.\n" );

  test_cond ( );
  test_determinant ( );
  test_eigen ( );
  test_inverse ( );
  test_null ( );
  test_plu ( );
  test_solution ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "TEST_MAT_PRB\n" );
  printf ( "  Normal end of execution.\n" );

  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void test_cond ( void )

/******************************************************************************/
/*
  Purpose:

    TEST_COND tests the condition number computations.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    08 April 2012

  Author:

    John Burkardt
*/
{
  double alpha;
  double beta;
  double cond;
  int n;
  int seed;
  char title[21];

  printf ( "\n" );
  printf ( "TEST_COND\n" );
  printf ( "  Compute the condition number of an example of each\n" );
  printf ( "  test matrix\n" );
  printf ( "\n" );
  printf ( "  Matrix title             N      COND\n" );
  printf ( "\n" );
/*
  AEGERTER matrix.
*/
  strcpy ( title, "AEGERTER            " );
  n = 5;
  cond = aegerter_condition ( n );
  printf ( "  %20s  %4d  %14.6g\n", title, n, cond );
/*
  BAB matrix.
*/
  strcpy ( title, "BAB                 " );
  n = 5;
  seed = 123456789;
  alpha = r8_nint ( 50.0 * r8_uniform_01 ( &seed ) - 25.0 ) / 5.0;
  beta = r8_nint ( 50.0 * r8_uniform_01 ( &seed ) - 25.0 ) / 5.0;
  cond = bab_condition ( n, alpha, beta );
  printf ( "  %20s  %4d  %14.6g\n", title, n, cond );
/*
  BODEWIG matrix.
*/
  strcpy ( title, "BODEWIG             " );
  n = 4;
  cond = bodewig_condition ( );
  printf ( "  %20s  %4d  %14.6g\n", title, n, cond );
/*
  COMBIN matrix.
*/
  strcpy ( title, "COMBIN              " );
  n = 5;
  seed = 123456789;
  alpha = r8_uniform_01 ( &seed );
  alpha = r8_nint ( 50.0 * alpha ) / 5.0;
  beta = r8_uniform_01 ( &seed );
  beta = r8_nint ( 50.0 * beta ) / 5.0;
  cond = combin_condition ( alpha, beta, n );
  printf ( "  %20s  %4d  %14.6g\n", title, n, cond );
/*
  CONEX3 matrix.
*/
  strcpy ( title, "CONEX3              " );
  n = 5;
  cond = conex3_condition ( n );
  printf ( "  %20s  %4d  %14.6g\n", title, n, cond );
/*
  RUTIS5 matrix.
*/
  strcpy ( title, "RUTIS5              " );
  n = 4;
  cond = rutis5_condition ( );
  printf ( "  %20s  %4d  %14.6g\n", title, n, cond );
/*
  SUMMATION matrix.
*/
  strcpy ( title, "SUMMATION           " );
  n = 5;
  cond = summation_condition ( n );
  printf ( "  %20s  %4d  %14.6g\n", title, n, cond );
/*
  TRI_UPPER matrix.
*/
  strcpy ( title, "TRI_UPPER           " );
  n = 5;
  seed = 123456789;
  alpha = r8_uniform_01 ( &seed );
  alpha = r8_nint ( 50.0 * alpha - 25.0 ) / 5.0;
  cond = tri_upper_condition ( alpha, n );
  printf ( "  %20s  %4d  %14.6g\n", title, n, cond );
/*
  WILK03 matrix.
*/
  strcpy ( title, "WILK03              " );
  n = 3;
  cond = wilk03_condition ( );
  printf ( "  %20s  %4d  %14.6g\n", title, n, cond );
/*
  WILSON matrix.
*/
  strcpy ( title, "WILSON              " );
  n = 4;
  cond = wilson_condition ( );
  printf ( "  %20s  %4d  %14.6g\n", title, n, cond );

  return;
}
/******************************************************************************/

void test_determinant ( void )

/******************************************************************************/
/*
  Purpose:

    TEST_DETERMINANT tests the determinant computations.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 July 2011

  Author:

    John Burkardt
*/
{
  double *a;
  double alpha;
  double b;
  double beta;
  int col_num;
  double d1;
  double d2;
  double d3;
  double d4;
  double d5;
  double da;
  double determ1;
  double determ2;
  double di;
  double gamma;
  int i;
  int ii;
  int jj;
  int k;
  double *l;
  int m;
  int n;
  double norm_frobenius;
  double *p;
  double perturb;
  int *pivot;
  double prob;
  int rank;
  int row_num;
  int seed;
  int seed_save;
  char title[21];
  double *u;
  double *v1;
  double *v2;
  double *v3;
  double *w;
  double *x;
  int x_n;
  double *y;
  int y_n;
  double y_sum;
  double *z;
  
  printf ( "\n" );
  printf ( "TEST_DETERMINANT\n" );
  printf ( "  Compute the determinants of an example of each\n" );
  printf ( "  test matrix; compare with the determinant routine,\n" );
  printf ( "  if available.  Print the matrix Frobenius norm\n" );
  printf ( "  for an estimate of magnitude.\n" );
  printf ( "\n" );
  printf ( "  Matrix title             N          Determ          Determ          ||A||\n" );
  printf ( "\n" );
/*
  AEGERTER matrix.
*/
  strcpy ( title, "AEGERTER            " );
  n = 5;
  a = aegerter ( n );
  determ1 = aegerter_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  %20s  %4d  %14g  %14g  %14f\n", 
    title, n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  ANTICIRCULANT matrix.
*/
  strcpy ( title, "ANTICIRCULANT       " );
  n = 3;
  seed = 123456789;
  x = r8vec_uniform_01_new ( n, &seed );
  for ( i = 0; i < n; i++ )
  {
    x[i] = ( r8_nint ( 50.0 * x[i] - 25.0 ) ) / 5.0;
  }
  a = anticirculant ( n, n, x );
  determ1 = anticirculant_determinant ( n, x );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  %20s  %4d  %14g  %14g  %14f\n", 
    title, n, determ1, determ2, norm_frobenius );
  free ( a );
  free ( x );
/*
  ANTICIRCULANT matrix.
*/
  strcpy ( title, "ANTICIRCULANT       " );
  n = 4;
  seed = 123456789;
  x = r8vec_uniform_01_new ( n, &seed );
  for ( i = 0; i < n; i++ )
  {
    x[i] = ( r8_nint ( 50.0 * x[i] - 25.0 ) ) / 5.0;
  }
  a = anticirculant ( n, n, x );
  determ1 = anticirculant_determinant ( n, x );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  %20s  %4d  %14g  %14g  %14f\n", 
    title, n, determ1, determ2, norm_frobenius );
  free ( a );
  free ( x );
/*
  ANTICIRCULANT matrix.
*/
  strcpy ( title, "ANTICIRCULANT       " );
  n = 5;
  seed = 123456789;
  x = r8vec_uniform_01_new ( n, &seed );
  for ( i = 0; i < n; i++ )
  {
    x[i] = ( r8_nint ( 50.0 * x[i] - 25.0 ) ) / 5.0;
  }
  a = anticirculant ( n, n, x );
  determ1 = anticirculant_determinant ( n, x );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  %20s  %4d  %14g  %14g  %14f\n", 
    title, n, determ1, determ2, norm_frobenius );
  free ( a );
  free ( x );
/*
  ANTIHADAMARD matrix.
*/
  strcpy ( title, "ANTIHADAMARD        " );
  n = 5;
  a = antihadamard ( n );
  determ1 = antihadamard_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  %20s  %4d  %14g  %14g  %14f\n", 
    title, n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  ANTISYMM_RANDOM matrix.
*/
  strcpy ( title, "ANTISYMM_RANDOM     " );
  n = 5;
  seed = 123456789;
  a = antisymm_random ( n, &seed );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  %20s  %4d                  %14g  %14f\n", 
    title, n,          determ2, norm_frobenius );
  free ( a );
/*
  ANTISYMM_RANDOM matrix.
*/
  strcpy ( title, "ANTISYMM_RANDOM     " );
  n = 6;
  seed = 123456789;
  a = antisymm_random ( n, &seed );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  %20s  %4d                  %14g  %14f\n", 
    title, n,          determ2, norm_frobenius );
  free ( a );
/*
  BAB matrix.
*/
  strcpy ( title, "BAB                 " );
  n = 5;
  seed = 123456789;
  alpha = r8_nint ( 50.0 * r8_uniform_01 ( &seed ) - 25.0 ) / 5.0;
  beta = r8_nint ( 50.0 * r8_uniform_01 ( &seed ) - 25.0 ) / 5.0;
  a = bab ( n, alpha, beta );
  determ1 = bab_determinant ( n, alpha, beta );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  %20s  %4d  %14g  %14g  %14f\n", 
    title, n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  BIMARKOV_RANDOM matrix.
*/
  strcpy ( title, "BIMARKOV_RANDOM     " );
  n = 5;
  seed = 123456789;
  a = bimarkov_random ( n, &seed );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  %20s  %4d                  %14g  %14f\n", 
    title, n,          determ2, norm_frobenius );
  free ( a );
/*
  BIS matrix.
*/
  strcpy ( title, "BIS                 " );
  n = 5;
  seed = 123456789;
  alpha = r8_nint ( 50.0 * r8_uniform_01 ( &seed ) ) / 5.0;
  beta = r8_nint ( 50.0 * r8_uniform_01 ( &seed ) ) / 5.0;
  a = bis ( alpha, beta, n, n );
  determ1 = bis_determinant ( alpha, beta, n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  %20s  %4d  %14g  %14g  %14f\n", 
    title, n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  BODEWIG matrix.
*/
  strcpy ( title, "BODEWIG             " );
  n = 4;
  a = bodewig ( );
  determ1 = bodewig_determinant ( );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  %20s  %4d  %14g  %14g  %14f\n", 
    title, n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  BOOTHROYD matrix.
*/
  strcpy ( title, "BOOTHROYD           " );
  n = 5;
  a = boothroyd ( n );
  determ1 = boothroyd_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  %20s  %4d  %14g  %14g  %14f\n", 
    title, n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  BORDERBAND matrix.
*/
  strcpy ( title, "BORDERBAND          " );
  n = 5;
  a = borderband ( n );
  determ1 = borderband_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  %20s  %4d  %14g  %14g  %14f\n", 
    title, n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  CARRY matrix.
*/
  n = 5;
  seed = 123456789;
  k = i4_uniform ( 2, 20, &seed );
  a = carry ( k, n );
  determ1 = carry_determinant ( k, n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  CARRY               %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  CAUCHY matrix.
*/
  n = 5;
  seed = 123456789;
  x = r8vec_uniform_01_new ( n, &seed );
  y = r8vec_uniform_01_new ( n, &seed );
  a = cauchy ( n, x, y );
  determ1 = cauchy_determinant ( n, x, y );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  CAUCHY              %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
  free ( x );
  free ( y );
/*
  CHEBY_DIFF1 matrix.
*/
  n = 5;
  a = cheby_diff1 ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  CHEBY_DIFF1         %4d                  %14g  %14f\n", 
    n,          determ2, norm_frobenius );
  free ( a );
/*
  CHEBY_DIFF1 matrix.
*/
  n = 6;
  a = cheby_diff1 ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  CHEBY_DIFF1         %4d                  %14g  %14f\n", 
    n,          determ2, norm_frobenius );
  free ( a );
/*
  CHEBY_T matrix.
*/
  n = 5;
  a = cheby_t ( n );
  determ1 = cheby_t_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  CHEBY_T             %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  CHEBY_U matrix.
*/
  n = 5;
  a = cheby_u ( n );
  determ1 = cheby_u_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  CHEBY_U             %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  CHEBY_VAN1 matrix.
*/
  n = 5;
  seed = 123456789;
  x = r8vec_uniform_01_new ( n, &seed );
  for ( i = 0; i < n; i++ )
  {
    x[i] = r8_nint ( 10.0 * x[i] - 5.0 );
  }
  a = cheby_van1 ( n, x );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  CHEBY_VAN1          %4d                  %14g  %14f\n", 
    n,          determ2, norm_frobenius );
  free ( a );
  free ( x );
/*
  CHEBY_VAN2 matrix.
*/
  for ( n = 2; n <= 10; n++ )
  {
    a = cheby_van2 ( n );
    determ1 = cheby_van2_determinant ( n );
    determ2 = r8mat_determinant ( n, a );
    norm_frobenius = r8mat_norm_fro ( n, n, a );
    printf ( "  CHEBY_VAN2          %4d  %14g  %14g  %14f\n", 
      n, determ1, determ2, norm_frobenius );
    free ( a );
  }
/*
  CHEBY_VAN3 matrix.
*/
  n = 5;
  a = cheby_van3 ( n );
  determ1 = cheby_van3_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  CHEBY_VAN3          %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  CHOW matrix.
*/
  n = 5;
  seed = 123456789;
  alpha = r8_uniform_01 ( &seed );
  alpha = r8_nint ( 50.0 * alpha ) / 5.0;
  beta = r8_uniform_01 ( &seed );
  beta = r8_nint ( 50.0 * beta ) / 5.0;
  a = chow ( alpha, beta, n, n );
  determ1 = chow_determinant ( alpha, beta, n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  CHOW                %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  CIRCULANT matrix.
*/
  n = 5;
  seed = 123456789;
  x = r8vec_uniform_01_new ( n, &seed );
  for ( i = 0; i < n; i++ )
  {
    x[i] = r8_nint ( 50.0 * x[i] - 25.0 ) / 5.0;
  }
  a = circulant ( n, n, x );
  determ1 = circulant_determinant ( n, x );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  CIRCULANT           %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
  free ( x );
/*
  CIRCULANT2 matrix.
*/
  n = 3;
  a = circulant2 ( n );
  determ1 = circulant2_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  CIRCULANT2          %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  CIRCULANT2 matrix.
*/
  n = 4;
  a = circulant2 ( n );
  determ1 = circulant2_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  CIRCULANT2          %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  CIRCULANT2 matrix.
*/
  n = 5;
  a = circulant2 ( n );
  determ1 = circulant2_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  CIRCULANT2          %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  CLEMENT1 matrix.
*/
  n = 5;
  a = clement1 ( n );
  determ1 = clement1_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  CLEMENT1            %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  CLEMENT1 matrix.
*/
  n = 6;
  a = clement1 ( n );
  determ1 = clement1_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  CLEMENT1            %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  CLEMENT2 matrix.
*/
  n = 5;
  a = clement2 ( n );
  determ1 = clement2_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  CLEMENT2            %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  CLEMENT2 matrix.
*/
  n = 6;
  a = clement2 ( n );
  determ1 = clement2_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  CLEMENT2            %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  CLEMENT3 matrix.
*/
  n = 5;
  seed = 123456789;
  x = r8vec_uniform_01_new ( n - 1, &seed );
  for ( i = 0; i < n - 1; i++ )
  {
    x[i] = r8_nint ( 10.0 * x[i] - 5.0 );
  }
  y = r8vec_uniform_01_new ( n - 1, &seed );
  for ( i = 0; i < n - 1; i++ )
  {
    y[i] = r8_nint ( 10.0 * y[i] - 5.0 );
  }
  a = clement3 ( n, x, y );
  determ1 = clement3_determinant ( n, x, y );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  CLEMENT3            %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
  free ( x );
  free ( y );
/*
  CLEMENT3 matrix.
*/
  n = 6;
  seed = 123456789;
  x = r8vec_uniform_01_new ( n - 1, &seed );
  for ( i = 0; i < n - 1; i++ )
  {
    x[i] = r8_nint ( 10.0 * x[i] - 5.0 );
  }
  y = r8vec_uniform_01_new ( n - 1, &seed );
  for ( i = 0; i < n - 1; i++ )
  {
    y[i] = r8_nint ( 10.0 * y[i] - 5.0 );
  }
  a = clement3 ( n, x, y );
  determ1 = clement3_determinant ( n, x, y );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  CLEMENT3            %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
  free ( x );
  free ( y );
/*
  COMBIN matrix.
*/
  n = 5;
  seed = 123456789;
  alpha = r8_uniform_01 ( &seed );
  alpha = r8_nint ( 50.0 * alpha ) / 5.0;
  beta = r8_uniform_01 ( &seed );
  beta = r8_nint ( 50.0 * beta ) / 5.0;
  a = combin ( alpha, beta, n );
  determ1 = combin_determinant ( alpha, beta, n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  COMBIN              %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  COMPANION matrix.
*/
  n = 5;
  seed = 123456789;
  x = r8vec_uniform_01_new ( n, &seed );
  for ( i = 0; i < n; i++ )
  {
    x[i] = r8_nint ( 10.0 * x[i] - 5.0 );
  }
  a = companion ( n, x );
  determ1 = companion_determinant ( n, x );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  COMPANION           %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
  free ( x );
/*
  COMPLEX_I matrix.
*/
  n = 2;
  a = complex_i ( );
  determ1 = complex_i_determinant ( );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  COMPLEX_I           %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  CONEX1 matrix.
*/
  n = 4;
  seed = 123456789;
  alpha = r8_uniform_01 ( &seed );
  alpha = r8_nint ( 50.0 * alpha ) / 5.0;
  a = conex1 ( alpha );
  determ1 = conex1_determinant ( alpha );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  CONEX1              %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  CONEX2 matrix.
*/
  n = 3;
  seed = 123456789;
  alpha = r8_uniform_01 ( &seed );
  alpha = r8_nint ( 50.0 * alpha ) / 5.0;
  a = conex2 ( alpha );
  determ1 = conex2_determinant ( alpha );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  CONEX2              %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  CONEX3 matrix.
*/
  n = 5;
  a = conex3 ( n );
  determ1 = conex3_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  CONEX3              %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  CONFERENCE matrix.
*/
  n = 6;
  a = conference ( n );
  determ1 = conference_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  CONFERENCE          %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  CREATION matrix.
*/
  n = 5;
  a = creation ( n, n );
  determ1 = creation_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  CREATION            %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  DAUB2 matrix.
*/
  n = 4;
  a = daub2 ( n );
  determ1 = daub2_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  DAUB2               %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  DAUB4 matrix.
*/
  n = 8;
  a = daub4 ( n );
  determ1 = daub4_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  DAUB4               %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  DAUB6 matrix.
*/
  n = 12;
  a = daub6 ( n );
  determ1 = daub6_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  DAUB6               %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  DAUB8 matrix.
*/
  n = 16;
  a = daub8 ( n );
  determ1 = daub8_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  DAUB8               %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  DAUB10 matrix.
*/
  n = 20;
  a = daub10 ( n );
  determ1 = daub10_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  DAUB10              %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  DAUB12 matrix.
*/
  n = 24;
  a = daub12 ( n );
  determ1 = daub12_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  DAUB12              %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  DIAGONAL matrix.
*/
  n = 5;
  seed = 123456789;
  x = r8vec_uniform_01_new ( n, &seed );
  for ( i = 0; i < n; i++ )
  {
    x[i] = r8_nint ( 10.0 * x[i] - 5.0 );
  }
  a = diagonal ( n, n, x );
  determ1 = diagonal_determinant ( n, x );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  DIAGONAL            %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
  free ( x );
/*
  DIF1 matrix.
*/
  n = 5;
  a = dif1 ( n, n );
  determ1 = dif1_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  DIF1                %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  DIF1CYCLIC matrix.
*/
  n = 5;
  a = dif1cyclic ( n );
  determ1 = dif1cyclic_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  DIF1CYCLIC          %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  DIF2 matrix.
*/
  n = 5;
  a = dif2 ( n, n );
  determ1 = dif2_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  DIF2                %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  DIF2CYCLIC matrix.
*/
  n = 5;
  a = dif2cyclic ( n );
  determ1 = dif2cyclic_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  DIF2CYCLIC          %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  DORR matrix.
*/
  n = 5;
  seed = 123456789;
  alpha = r8_uniform_01 ( &seed );
  alpha = r8_nint ( 50.0 * alpha ) / 5.0;
  a = dorr ( alpha, n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  DORR                %4d                  %14g  %14f\n", 
    n,          determ2, norm_frobenius );
  free ( a );
/*
  DOWNSHIFT matrix.
*/
  n = 5;
  a = downshift ( n );
  determ1 = downshift_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  DOWNSHIFT           %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  EBERLEIN matrix.
*/
  n = 5;
  seed = 123456789;
  alpha = r8_uniform_01 ( &seed );
  alpha = r8_nint ( 50.0 * alpha ) / 5.0;
  a = eberlein ( alpha, n );
  determ1 = eberlein_determinant ( alpha, n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  EBERLEIN            %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  EULERIAN matrix.
*/
  n = 5;
  a = eulerian ( n, n );
  determ1 = eulerian_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  EULERIAN            %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  EXCHANGE matrix.
*/
  n = 5;
  a = exchange ( n, n );
  determ1 = exchange_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  EXCHANGE            %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  FIBONACCI1 matrix.
*/
  n = 5;
  seed = 123456789;
  alpha = r8_uniform_01 ( &seed );
  alpha = r8_nint ( 50.0 * alpha ) / 5.0;
  beta = r8_uniform_01 ( &seed );
  beta = r8_nint ( 50.0 * beta ) / 5.0;
  a = fibonacci1 ( n, alpha, beta );
  determ1 = fibonacci1_determinant ( n, alpha, beta );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  FIBONACCI1          %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  FIBONACCI2 matrix.
*/
  n = 5;
  a = fibonacci2 ( n );
  determ1 = fibonacci2_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  FIBONACCI2          %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  FIBONACCI3 matrix.
*/
  n = 5;
  a = fibonacci3 ( n );
  determ1 = fibonacci3_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  FIBONACCI3          %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  FIEDLER matrix.
*/
  n = 7;
  seed = 123456789;
  x = r8vec_uniform_01_new ( n, &seed );
  for ( i = 0; i < n; i++ )
  {
    x[i] = r8_nint ( 50.0 * x[i] - 25.0 ) / 5.0;
  }
  a = fiedler ( n, n, x );
  determ1 = fiedler_determinant ( n, x );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  FIEDLER             %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
  free ( x );
/*
  FORSYTHE matrix.
*/
  n = 7;
  seed = 123456789;
  alpha = r8_uniform_01 ( &seed );
  alpha = r8_nint ( 50.0 * alpha ) / 5.0;
  beta = r8_uniform_01 ( &seed );
  beta = r8_nint ( 50.0 * beta ) / 5.0;
  a = forsythe ( alpha, beta, n );
  determ1 = forsythe_determinant ( alpha, beta, n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  FORSYTHE            %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  FOURIER_COSINE matrix.
*/
  n = 5;
  a = fourier_cosine ( n );
  determ1 = fourier_cosine_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  FOURIER_COSINE      %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  FOURIER_SINE matrix.
*/
  n = 5;
  a = fourier_sine ( n );
  determ1 = fourier_sine_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  FOURIER_SINE        %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  FRANK matrix.
*/
  n = 5;
  a = frank ( n );
  determ1 = frank_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  FRANK               %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  GEAR matrix.
*/
  for ( n = 4; n <= 8; n++ )
  {
    seed = 123456789;
    ii = i4_uniform ( -n, n, &seed );
    jj = i4_uniform ( -n, n, &seed );
    a = gear ( ii, jj, n );
    determ1 = gear_determinant ( ii, jj, n );
    determ2 = r8mat_determinant ( n, a );
    norm_frobenius = r8mat_norm_fro ( n, n, a );
    printf ( "  GEAR                %4d  %14g  %14g  %14f\n", 
      n, determ1, determ2, norm_frobenius );
    free ( a );
  }
/*
  GFPP matrix.
*/
  n = 5;
  seed = 123456789;
  alpha = r8_uniform_01 ( &seed );
  alpha = r8_nint ( 50.0 * alpha ) / 5.0;
  a = gfpp ( n, alpha );
  determ1 = gfpp_determinant ( n, alpha );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  GFPP                %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  GIVENS matrix.
*/
  n = 5;
  a = givens ( n, n );
  determ1 = givens_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  GIVENS              %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  GK316 matrix.
*/
  n = 5;
  a = gk316 ( n );
  determ1 = gk316_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  GK316               %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  GK323 matrix.
*/
  n = 5;
  a = gk323 ( n, n );
  determ1 = gk323_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  GK323               %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  GK324 matrix.
*/
  n = 5;
  seed = 123456789;
  x = r8vec_uniform_01_new ( n, &seed );
  for ( i = 0; i < n; i++ )
  {
    x[i] = r8_nint ( 50.0 * x[i] - 25.0 ) / 5.0;
  }
  a = gk324 ( n, n, x );
  determ1 = gk324_determinant ( n, x );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  GK324               %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
  free ( x );
/*
  GRCAR matrix.
*/
  n = 5;
  seed = 123456789;
  k = i4_uniform ( 1, n - 1, &seed );
  a = grcar ( n, n, k );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  GRCAR               %4d                  %14g  %14f\n", 
    n,          determ2, norm_frobenius );
  free ( a );
/*
  HADAMARD matrix.
*/
  n = 5;
  a = hadamard ( n, n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  HADAMARD            %4d                  %14g  %14f\n", 
    n,          determ2, norm_frobenius );
  free ( a );
/*
  HANKEL matrix.
*/
  n = 5;
  seed = 123456789;
  x = r8vec_uniform_01_new ( 2 * n - 1, &seed );
  for ( i = 0; i < 2 * n - 1; i++ )
  {
    x[i] = r8_nint ( 50.0 * x[i] - 25.0 ) / 5.0;
  }
  a = hankel ( n, x );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  HANKEL              %4d                  %14g  %14f\n", 
    n,          determ2, norm_frobenius );
  free ( a );
  free ( x );
/*
  HANOWA matrix.
*/
  n = 6;
  seed = 123456789;
  alpha = r8_uniform_01 ( &seed );
  alpha = r8_nint ( 50.0 * alpha ) / 5.0;
  a = hanowa ( alpha, n );
  determ1 = hanowa_determinant ( alpha, n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  HANOWA              %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  HARMAN matrix.
*/
  n = 8;
  a = harman ( );
  determ1 = harman_determinant ( );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  HARMAN              %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  HARTLEY matrix.
*/
  for ( n = 5; n <= 8; n++ )
  {
    a = hartley ( n );
    determ1 = hartley_determinant ( n );
    determ2 = r8mat_determinant ( n, a );
    norm_frobenius = r8mat_norm_fro ( n, n, a );
    printf ( "  HARTLEY             %4d  %14g  %14g  %14f\n", 
      n, determ1, determ2, norm_frobenius );
    free ( a );
  }
/*
  HELMERT matrix.
*/
  n = 5;
  a = helmert ( n );
  determ1 = helmert_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  HELMERT             %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  HELMERT2 matrix.
*/
  n = 5;
  seed = 123456789;
  x = r8vec_uniform_01_new ( n, &seed );
  for ( i = 0; i < n; i++ )
  {
    x[i] = r8_nint ( 50.0 * x[i] - 25.0 ) / 5.0;
  }
  a = helmert2 ( n, x );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  HELMERT2            %4d                  %14g  %14f\n", 
    n,          determ2, norm_frobenius );
  free ( a );
  free ( x );
/*
  HERMITE matrix.
*/
  n = 5;
  a = hermite ( n );
  determ1 = hermite_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  HERMITE             %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  HERNDON matrix.
*/
  n = 5;
  a = herndon ( n );
  determ1 = herndon_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  HERNDON             %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  HILBERT matrix.
*/
  n = 5;
  a = hilbert ( n, n );
  determ1 = hilbert_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  HILBERT             %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  HOUSEHOLDER matrix.
*/
  n = 5;
  seed = 123456789;
  x = r8vec_uniform_01_new ( n, &seed );
  for ( i = 0; i < n; i++ )
  {
    x[i] = r8_nint ( 50.0 * x[i] - 25.0 ) / 5.0;
  }
  a = householder ( n, x );
  determ1 = householder_determinant ( n, x );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  HOUSEHOLDER         %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
  free ( x );
/*
  IDEM_RANDOM matrix.
*/
  n = 5;
  seed = 123456789;
  rank = i4_uniform ( 0, n, &seed );
  a = idem_random ( n, rank, &seed );
  determ1 = idem_random_determinant ( n, rank );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  IDEM_RANDOM         %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  IDENTITY matrix.
*/
  n = 5;
  a = identity ( n, n );
  determ1 = identity_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  IDENTITY            %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  IJFACT1 matrix.
*/
  n = 5;
  a = ijfact1 ( n );
  determ1 = ijfact1_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  IJFACT1             %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  IJFACT2 matrix.
*/
  n = 5;
  a = ijfact2 ( n );
  determ1 = ijfact2_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  IJFACT2             %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  ILL3 matrix.
*/
  n = 3;
  a = ill3 ( );
  determ1 = ill3_determinant ( );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  ILL3                %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  INTEGRATION matrix.
*/
  n = 6;
  seed = 123456789;
  alpha = r8_uniform_01 ( &seed );
  alpha = r8_nint ( 50.0 * alpha ) / 5.0;
  a = integration ( alpha, n );
  determ1 = integration_determinant ( alpha, n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  INTEGRATION         %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  INVOL matrix.
*/
  n = 5;
  a = invol ( n );
  determ1 = invol_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  INVOL               %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  INVOL_RANDOM matrix.
*/
  n = 5;
  seed = 123456789;
  rank = i4_uniform ( 0, n, &seed );
  a = invol_random ( n, rank, &seed );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  INVOL_RANDOM        %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  JACOBI matrix.
*/
  n = 5;
  a = jacobi ( n, n );
  determ1 = jacobi_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  JACOBI              %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  JACOBI matrix.
*/
  n = 6;
  a = jacobi ( n, n );
  determ1 = jacobi_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  JACOBI              %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  JORDAN matrix.
*/
  n = 6;
  seed = 123456789;
  alpha = r8_uniform_01 ( &seed );
  alpha = r8_nint ( 50.0 * alpha ) / 5.0;
  a = jordan ( alpha, n, n );
  determ1 = jordan_determinant ( alpha, n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  JORDAN              %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  KAHAN matrix.
*/
  n = 5;
  seed = 123456789;
  alpha = r8_uniform_01 ( &seed );
  alpha = r8_nint ( 50.0 * alpha ) / 5.0;
  a = kahan ( alpha, n, n );
  determ1 = kahan_determinant ( alpha, n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  KAHAN               %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  KERSHAW matrix.
*/
  n = 4;
  a = kershaw ( );
  determ1 = kershaw_determinant ( );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  KERSHAW             %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  KERSHAWTRI matrix.
*/
  n = 5;
  x_n = ( n + 1 ) / 2;
  seed = 123456789;
  x = r8vec_uniform_01_new ( x_n, &seed );
  for ( i = 0; i < x_n; i++ )
  {
    x[i] = r8_nint ( 50.0 * x[i] - 25.0 ) / 5.0;
  }
  a = kershawtri ( n, x );
  determ1 = kershawtri_determinant ( n, x );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  KERSHAWTRI          %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
  free ( x );
/*
  KMS matrix.
*/
  n = 5;
  seed = 123456789;
  alpha = r8_uniform_01 ( &seed );
  alpha = r8_nint ( 50.0 * alpha ) / 5.0;
  a = kms ( alpha, n, n );
  determ1 = kms_determinant ( alpha, n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  KMS                 %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  LAGUERRE matrix.
*/
  n = 5;
  a = laguerre ( n );
  determ1 = laguerre_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  LAGUERRE            %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  LEHMER matrix.
*/
  n = 5;
  a = lehmer ( n, n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  LEHMER              %4d                  %14g  %14f\n", 
    n,          determ2, norm_frobenius );
  free ( a );
/*
  LESLIE matrix.
*/
  n = 4;
  b = 0.025;
  di = 0.010;
  da = 0.100;
  a = leslie ( b, di, da );
  determ1 = leslie_determinant ( b, di, da );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  LESLIE              %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  LESP matrix.
*/
  n = 5;
  a = lesp ( n, n );
  determ1 = lesp_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  LESP                %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  LIETZKE matrix.
*/
  n = 5;
  a = lietzke ( n );
  determ1 = lietzke_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  LIETZKE             %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  LIGHTS_OUT matrix.
*/
  if ( 0 )
  {
    row_num = 5;
    col_num = 5;
    n = row_num * col_num;
/*
    a = lights_out ( row_num, col_num, n );
*/
    determ2 = r8mat_determinant ( n, a );
    norm_frobenius = r8mat_norm_fro ( n, n, a );
    printf ( "  LIGHTS_OUT          %4d                  %14g  %14f\n", 
      n,          determ2, norm_frobenius );
    }
  else
  {
    printf ( "  LIGHTS_OUT          -----Not ready----\n" );
  }
/*
  LINE_ADJ matrix.
*/
  n = 5;
  a = line_adj ( n );
  determ1 = line_adj_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  LINE_ADJ            %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  LINE_LOOP_ADJ matrix.
*/
  n = 5;
  a = line_loop_adj ( n );
  determ1 = line_loop_adj_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  LINE_LOOP_ADJ       %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  LOEWNER matrix.
*/
  n = 5;
  seed = 123456789;
  w = r8vec_uniform_01_new ( n, &seed );
  x = r8vec_uniform_01_new ( n, &seed );
  y = r8vec_uniform_01_new ( n, &seed );
  z = r8vec_uniform_01_new ( n, &seed );
  a = loewner ( w, x, y, z, n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  LOEWNER             %4d                  %14g  %14f\n", 
    n,          determ2, norm_frobenius );
  free ( a );
  free ( w );
  free ( x );
  free ( y );
  free ( z );
/*
  LOTKIN matrix.
*/
  n = 5;
  a = lotkin ( n, n );
  determ1 = lotkin_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  LOTKIN              %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  MARKOV_RANDOM matrix.
*/
  n = 5;
  seed = 123456789;
  a = markov_random ( n, &seed );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  MARKOV_RANDOM       %4d                  %14g  %14f\n", 
    n,          determ2, norm_frobenius );
  free ( a );
/*
  MAXIJ matrix.
*/
  n = 5;
  a = maxij ( n, n );
  determ1 = maxij_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  MAXIJ               %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  MILNES matrix.
*/
  n = 5;
  seed = 123456789;
  x = r8vec_uniform_01_new ( n, &seed );
  for ( i = 0; i < n; i++ )
  {
    x[i] = r8_nint ( 50.0 * x[i] - 25.0 ) / 5.0;
  }
  a = milnes ( n, n, x );
  determ1 = milnes_determinant ( n, x );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  MILNES              %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
  free ( x );
/*
  MINIJ matrix.
*/
  n = 5;
  a = minij ( n, n );
  determ1 = minij_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  MINIJ               %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  MOLER1 matrix.
*/
  n = 5;
  seed = 123456789;
  alpha = r8_uniform_01 ( &seed );
  alpha = r8_nint ( 50.0 * alpha ) / 5.0;
  a = moler1 ( alpha, n, n );
  determ1 = moler1_determinant ( alpha, n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  MOLER1              %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  MOLER2 matrix.
*/
  n = 5;
  a = moler2 ( );
  determ1 = moler2_determinant ( );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  MOLER2              %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  MOLER3 matrix.
*/
  n = 5;
  a = moler3 ( n, n );
  determ1 = moler3_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  MOLER3              %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  NEUMANN matrix.
*/
  row_num = 5;
  col_num = 5;
  n = row_num * col_num;
  a = neumann ( row_num, col_num );
  determ1 = neumann_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  NEUMANN             %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  ONE matrix.
*/
  n = 5;
  a = one ( n, n );
  determ1 = one_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  ONE                 %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  ORTEGA matrix.
*/
  n = 5;
  seed = 123456789;
  v1 = r8vec_uniform_01_new ( n, &seed );
  for ( i = 0; i < n; i++ )
  {
    v1[i] = r8_nint ( 50.0 * v1[i] - 25.0 ) / 5.0;
  }
  v2 = r8vec_uniform_01_new ( n, &seed );
  for ( i = 0; i < n; i++ )
  {
    v2[i] = r8_nint ( 50.0 * v2[i] - 25.0 ) / 5.0;
  }
  v3 = r8vec_uniform_01_new ( n, &seed );
  for ( i = 0; i < n; i++ )
  {
    v3[i] = r8_nint ( 50.0 * v3[i] - 25.0 ) / 5.0;
  }
  a = ortega ( n, v1, v2, v3 );
  determ1 = ortega_determinant ( n, v1, v2, v3 );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  ORTEGA              %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
  free ( v1 );
  free ( v2 );
  free ( v3 );
/*
  ORTH_RANDOM matrix.
*/
  n = 5;
  seed = 123456789;
  a = orth_random ( n, &seed );
  determ1 = orth_random_determinant ( n, &seed );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  ORTH_RANDOM         %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  ORTH_SYMM matrix.
*/
  n = 5;
  a = orth_symm ( n );
  determ1 = orth_symm_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  ORTH_SYMM           %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  OTO matrix.
*/
  n = 5;
  a = oto ( n, n );
  determ1 = oto_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  OTO                 %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  PARTER matrix.
*/
  n = 5;
  a = parter ( n, n );
  determ1 = parter_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  PARTER              %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  PASCAL1 matrix.
*/
  n = 5;
  a = pascal1 ( n );
  determ1 = pascal1_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  PASCAL1             %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  PASCAL2 matrix.
*/
  n = 5;
  a = pascal2 ( n );
  determ1 = pascal2_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  PASCAL2             %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  PASCAL3 matrix.
*/
  n = 5;
  seed = 123456789;
  alpha = r8_uniform_01 ( &seed );
  alpha = r8_nint ( 50.0 * alpha - 25.0 ) / 5.0;
  a = pascal3 ( n, alpha );
  determ1 = pascal3_determinant ( n, alpha );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  PASCAL3             %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  PDS_RANDOM matrix.
*/
  n = 5;
  seed = 123456789;
  seed_save = seed;
  a = pds_random ( n, &seed );
  seed = seed_save;
  determ1 = pds_random_determinant ( n, &seed );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  PDS_RANDOM          %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  PEI matrix.
*/
  n = 5;
  seed = 123456789;
  alpha = r8_uniform_01 ( &seed );
  alpha = r8_nint ( 50.0 * alpha ) / 5.0;
  a = pei ( alpha, n );
  determ1 = pei_determinant ( alpha, n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  PEI                 %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  PERMUTATION_RANDOM matrix.
*/
  n = 5;
  seed_save = 123456789;
  seed = seed_save;
  a = permutation_random ( n, &seed );
  seed = seed_save;
  determ1 = permutation_random_determinant ( n, &seed );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  PERMUTATION_RANDOM  %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  PLU matrix.
*/
  n = 5;
  l = ( double * ) malloc ( n * n * sizeof ( double ) );
  p = ( double * ) malloc ( n * n * sizeof ( double ) );
  pivot = ( int * ) malloc ( n * sizeof ( int ) );
  u = ( double * ) malloc ( n * n * sizeof ( double ) );
  for ( i = 0; i < n; i++ )
  {
    pivot[i] = i + 1;
  }
  a = plu ( n, pivot, p, l, u );
  determ1 = plu_determinant ( n, p, l, u );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  PLU                 %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
  free ( l );
  free ( p );
  free ( pivot );
  free ( u );
/*
  POISSON matrix.
*/
  row_num = 5;
  col_num = 5;
  n = row_num * col_num;
  a = poisson ( row_num, col_num, n );
  determ1 = poisson_determinant ( row_num, col_num, n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  POISSON             %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  PROLATE matrix.
*/
  n = 5;
  seed = 123456789;
  alpha = r8_uniform_01 ( &seed );
  alpha = r8_nint ( 50.0 * alpha ) / 5.0;
  a = prolate ( alpha, n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  PROLATE             %4d                  %14g  %14f\n", 
    n,          determ2, norm_frobenius );
  free ( a );
/*
  RECTANGLE_ADJ matrix.
*/
  if ( 0 )
  {
    row_num = 5;
    col_num = 5;
    n = row_num * col_num;
/*  a = rectangle_adj ( row_num, col_num, n ); */
    determ1 = rectangle_adj_determinant ( row_num, col_num );
    determ2 = r8mat_determinant ( n, a );
    norm_frobenius = r8mat_norm_fro ( n, n, a );
    printf ( "  RECTANGLE_ADJ       %4d  %14g  %14g  %14f\n", 
      n, determ1, determ2, norm_frobenius );
    free ( a );
  }
  else
  {
    printf ( "  RECTANGLE_ADJ       -----Not ready-----\n" );
  }
/*
  REDHEFFER matrix.
*/
  n = 5;
  a = redheffer ( n );
  determ1 = redheffer_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  REDHEFFER           %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  REF_RANDOM matrix.
*/
  n = 5;
  prob = 0.65;
  seed_save = 123456789;
  seed = seed_save;
  a = ref_random ( n, n, prob, &seed );
  seed = seed_save;
  determ1 = ref_random_determinant ( n, prob, &seed );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  REF_RANDOM          %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  REF_RANDOM matrix.
*/
  n = 5;
  prob = 0.85;
  seed_save = 123456789;
  seed = seed_save;
  a = ref_random ( n, n, prob, &seed );
  seed = seed_save;
  determ1 = ref_random_determinant ( n, prob, &seed );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  REF_RANDOM          %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  RIEMANN matrix.
*/
  n = 5;
  a = riemann ( n, n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  RIEMANN             %4d                  %14g  %14f\n", 
    n,          determ2, norm_frobenius );
  free ( a );
/*
  RING_ADJ matrix.
*/
  for ( n = 1; n <= 8; n++ )
  {
    a = ring_adj ( n );
    determ1 = ring_adj_determinant ( n );
    determ2 = r8mat_determinant ( n, a );
    norm_frobenius = r8mat_norm_fro ( n, n, a );
    printf ( "  RING_ADJ            %4d  %14g  %14g  %14f\n", 
      n, determ1, determ2, norm_frobenius );
    free ( a );
  }
/*
  RIS matrix.
*/
  n = 5;
  a = ris ( n );
  determ1 = ris_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  RIS                 %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  RODMAN matrix.
*/
  n = 5;
  seed = 123456789;
  alpha = r8_uniform_01 ( &seed );
  alpha = r8_nint ( 50.0 * alpha - 25.0 ) / 5.0;
  a = rodman ( alpha, n, n );
  determ1 = rodman_determinant ( alpha, n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  RODMAN              %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  ROSSER1 matrix.
*/
  n = 8;
  a = rosser1 ( );
  determ1 = rosser1_determinant ( );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  ROSSER1             %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  ROUTH matrix.
*/
  n = 5;
  seed = 123456789;
  x = r8vec_uniform_01_new ( n, &seed );
  for ( i = 0; i < n; i++ )
  {
    x[i] = r8_nint ( 50.0 * x[i] ) / 5.0;
  }
  a = routh ( n, x );
  determ1 = routh_determinant ( n, x );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  ROUTH               %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
  free ( x );
/*
  RUTIS1 matrix.
*/
  n = 4;
  a = rutis1 ( );
  determ1 = rutis1_determinant ( );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  RUTIS1              %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  RUTIS2 matrix.
*/
  n = 4;
  a = rutis2 ( );
  determ1 = rutis2_determinant ( );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  RUTIS2              %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  RUTIS3 matrix.
*/
  n = 4;
  a = rutis3 ( );
  determ1 = rutis3_determinant ( );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  RUTIS3              %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  RUTIS4 matrix.
*/
  n = 4;
  a = rutis4 ( n );
  determ1 = rutis4_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  RUTIS4              %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  RUTIS5 matrix.
*/
  n = 4;
  a = rutis5 ( );
  determ1 = rutis5_determinant ( );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  RUTIS5              %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  SCHUR_BLOCK matrix.
*/
  n = 5;
  x_n = ( n + 1 ) / 2;
  y_n = n / 2;
  seed = 123456789;
  x = r8vec_uniform_01_new ( x_n, &seed );
  for ( i = 0; i < x_n; i++ )
  {
    x[i] = r8_nint ( 50.0 * x[i] - 25.0 ) / 5.0;
  }
  y = r8vec_uniform_01_new ( y_n, &seed );
  for ( i = 0; i < y_n; i++ )
  {
    y[i] = r8_nint ( 50.0 * y[i] - 25.0 ) / 5.0;
  }
  a = schur_block ( n, x, y );
  determ1 = schur_block_determinant ( n, x, y );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  SCHUR_BLOCK         %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
  free ( x );
  free ( y );
/*
  SKEW_CIRCULANT matrix.
*/
  n = 5;
  seed = 123456789;
  x = r8vec_uniform_01_new ( n, &seed );
  for ( i = 0; i < n; i++ )
  {
    x[i] = r8_nint ( 50.0 * x[i] - 25.0 ) / 5.0;
  }
  a = skew_circulant ( n, n, x );
  determ1 = skew_circulant_determinant ( n, x );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  SKEW_CIRCULANT      %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
  free ( x );
/*
  SPLINE matrix.
*/
  n = 5;
  seed = 123456789;
  x = r8vec_uniform_01_new ( n, &seed );
  for ( i = 0; i < n; i++ )
  {
    x[i] = r8_nint ( 50.0 * x[i] - 25.0 ) / 5.0;
  }
  a = spline ( n, x );
  determ1 = spline_determinant ( n, x );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  SPLINE              %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
  free ( x );
/*
  STIRLING matrix.
*/
  n = 5;
  a = stirling ( n, n );
  determ1 = stirling_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  STIRLING            %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  STRIPE matrix.
*/
  n = 5;
  a = stripe ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  STRIPE              %4d                  %14g  %14f\n", 
    n,          determ2, norm_frobenius );
  free ( a );
/*
  SUMMATION matrix.
*/
  n = 5;
  a = summation ( n, n );
  determ1 = summation_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  SUMMATION           %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  SWEET1 matrix.
*/
  n = 6;
  perturb = 0.0;
  a = sweet1 ( perturb );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  SWEET1              %4d                  %14g  %14f\n", 
    n,          determ2, norm_frobenius );
  free ( a );
/*
  SWEET2 matrix.
*/
  n = 6;
  perturb = 0.0;
  a = sweet2 ( perturb );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  SWEET2              %4d                  %14g  %14f\n", 
    n,          determ2, norm_frobenius );
  free ( a );
/*
  SWEET3 matrix.
*/
  n = 6;
  perturb = 0.0;
  a = sweet3 ( perturb );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  SWEET3              %4d                  %14g  %14f\n", 
    n,          determ2, norm_frobenius );
  free ( a );
/*
  SWEET4 matrix.
*/
  n = 13;
  perturb = 0.0;
  a = sweet4 ( perturb );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  SWEET4              %4d                  %14g  %14f\n", 
    n,          determ2, norm_frobenius );
  free ( a );
/*
  SYLVESTER matrix.
*/
  n = 5;
  x_n = 3 + 1;
  y_n = 2 + 1;
  seed = 123456789;
  x = r8vec_uniform_01_new ( x_n, &seed );
  for ( i = 0; i < x_n; i++ )
  {
    x[i] = r8_nint ( 50.0 * x[i] - 25.0 ) / 5.0;
  }
  y = r8vec_uniform_01_new ( y_n, &seed );
  for ( i = 0; i < y_n; i++ )
  {
    y[i] = r8_nint ( 50.0 * y[i] - 25.0 ) / 5.0;
  }
  a = sylvester ( n, x_n - 1, x, y_n - 1, y );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  SYLVESTER           %4d                  %14g  %14f\n", 
    n,          determ2, norm_frobenius );
  free ( a );
  free ( x );
  free ( y );
/*
  SYMM_RANDOM matrix.
*/
  n = 5;
  seed = 123456789;
  x = r8vec_uniform_01_new ( n, &seed );
  for ( i = 0; i < n; i++ )
  {
    x[i] = r8_nint ( 50.0 * x[i] - 25.0 ) / 5.0;
  }
  a = symm_random ( n, x, &seed );
  determ1 = symm_random_determinant ( n, x );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  SYMM_RANDOM         %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
  free ( x );
/*
  TOEPLITZ matrix.
*/
  n = 5;
  seed = 123456789;
  x = r8vec_uniform_01_new ( 2 * n - 1, &seed );
  for ( i = 0; i < 2 * n - 1; i++ )
  {
    x[i] = r8_nint ( 50.0 * x[i] - 25.0 ) / 5.0;
  }
  a = toeplitz ( n, x );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  TOEPLITZ            %4d                  %14g  %14f\n", 
    n,          determ2, norm_frobenius );
  free ( a );
  free ( x );
/*
  TOEPLITZ_5DIAG matrix.
*/
  n = 5;
  seed = 123456789;
  d1 = r8_uniform_01 ( &seed );
  d1 = r8_nint ( 50.0 * d1 - 25.0 ) / 5.0;
  d2 = r8_uniform_01 ( &seed );
  d2 = r8_nint ( 50.0 * d2 - 25.0 ) / 5.0;
  d3 = r8_uniform_01 ( &seed );
  d3 = r8_nint ( 50.0 * d3 - 25.0 ) / 5.0;
  d4 = r8_uniform_01 ( &seed );
  d4 = r8_nint ( 50.0 * d4 - 25.0 ) / 5.0;
  d5 = r8_uniform_01 ( &seed );
  d5 = r8_nint ( 50.0 * d5 - 25.0 ) / 5.0;
  a = toeplitz_5diag ( n, d1, d2, d3, d4, d5 );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  TOEPLITZ_5DIAG      %4d                  %14g  %14f\n", 
    n,          determ2, norm_frobenius );
  free ( a );
/*
  TOEPLITZ_5S matrix.
*/
  row_num = 5;
  col_num = 5;
  n = row_num * col_num;
  seed = 123456789;
  alpha = r8_uniform_01 ( &seed );
  alpha = r8_nint ( 50.0 * alpha - 25.0 ) / 5.0;
  beta = r8_uniform_01 ( &seed );
  beta = r8_nint ( 50.0 * beta - 25.0 ) / 5.0;
  gamma = r8_uniform_01 ( &seed );
  gamma = r8_nint ( 50.0 * gamma - 25.0 ) / 5.0;
  a = toeplitz_5s ( row_num, col_num, alpha, beta, gamma, n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  TOEPLITZ_5DS        %4d                  %14g  %14f\n", 
    n,          determ2, norm_frobenius );
  free ( a );
/*
  TOEPLITZ_PDS matrix.
*/
  m = 3;
  n = 5;
  seed = 123456789;
  x = r8vec_uniform_01_new ( m, &seed );
  y = r8vec_uniform_01_new ( m, &seed );
  y_sum = r8vec_sum ( m, y );
  for ( i = 0; i < m; i++ )
  {
    y[i] = y[i] / y_sum;
  }
  a = toeplitz_pds ( m, n, x, y );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  TOEPLITZ_PDS        %4d                  %14g  %14f\n", 
    n,          determ2, norm_frobenius );
  free ( a );
  free ( x );
  free ( y );
/*
  TOURNAMENT_RANDOM matrix.
*/
  n = 5;
  seed_save = 123456789;
  seed = seed_save;
  a = tournament_random ( n, &seed );
  seed = seed_save;
  determ1 = tournament_random_determinant ( n, &seed );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  TOURNAMENT_RANDOM   %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  TRANSITION_RANDOM matrix.
*/
  n = 5;
  seed = 123456789;
  a = transition_random ( n, &seed );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  TRANSITION_RANDOM   %4d                  %14g  %14f\n", 
    n,          determ2, norm_frobenius );
  free ( a );
/*
  TRENCH matrix.
*/
  n = 5;
  seed = 123456789;
  alpha = r8_uniform_01 ( &seed );
  alpha = r8_nint ( 50.0 * alpha - 25.0 ) / 5.0;
  a = trench ( alpha, n, n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  TRENCH              %4d                  %14g  %14f\n", 
    n,          determ2, norm_frobenius );
  free ( a );
/*
  TRI_UPPER matrix.
*/
  n = 5;
  seed = 123456789;
  alpha = r8_uniform_01 ( &seed );
  alpha = r8_nint ( 50.0 * alpha - 25.0 ) / 5.0;
  a = tri_upper ( alpha, n );
  determ1 = tri_upper_determinant ( alpha, n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  TRI_UPPER           %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  TRIS matrix.
*/
  n = 5;
  seed = 123456789;
  alpha = r8_uniform_01 ( &seed );
  alpha = r8_nint ( 50.0 * alpha - 25.0 ) / 5.0;
  beta = r8_uniform_01 ( &seed );
  beta = r8_nint ( 50.0 * beta - 25.0 ) / 5.0;
  gamma = r8_uniform_01 ( &seed );
  gamma = r8_nint ( 50.0 * gamma - 25.0 ) / 5.0;
  a = tris ( n, n, alpha, beta, gamma );
  determ1 = tris_determinant ( n, alpha, beta, gamma );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  TRIS                %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  TRIV matrix.
*/
  n = 5;
  seed = 123456789;
  x = r8vec_uniform_01_new ( n - 1, &seed );
  for ( i = 0; i < n - 1; i++ )
  {
    x[i] = r8_nint ( 50.0 * x[i] - 25.0 ) / 5.0;
  }
  y = r8vec_uniform_01_new ( n, &seed );
  for ( i = 0; i < n; i++ )
  {
    y[i] = r8_nint ( 50.0 * y[i] - 25.0 ) / 5.0;
  }
  z = r8vec_uniform_01_new ( n - 1, &seed );
  for ( i = 0; i < n - 1; i++ )
  {
    z[i] = r8_nint ( 50.0 * z[i] - 25.0 ) / 5.0;
  }
  a = triv ( n, x, y, z );
  determ1 = triv_determinant ( n, x, y, z );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  TRIV                %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
  free ( x );
  free ( y );
  free ( z );
/*
  TRIW matrix.
*/
  n = 5;
  seed = 123456789;
  k = i4_uniform ( 0, n - 1, &seed );
  alpha = r8_uniform_01 ( &seed );
  alpha = r8_nint ( 50.0 * alpha - 25.0 ) / 5.0;
  a = triw ( alpha, k, n );
  determ1 = triw_determinant ( alpha, k, n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  TRIW                %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  UPSHIFT matrix.
*/
  n = 5;
  a = upshift ( n );
  determ1 = upshift_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  UPSHIFT             %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  VAND1 matrix.
*/
  n = 5;
  seed = 123456789;
  x = r8vec_uniform_01_new ( n, &seed );
  for ( i = 0; i < n; i++ )
  {
    x[i] = r8_nint ( 50.0 * x[i] - 25.0 ) / 5.0;
  }
  a = vand1 ( n, x );
  determ1 = vand1_determinant ( n, x );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  VAND1               %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
  free ( x );
/*
  VAND2 matrix.
*/
  n = 5;
  seed = 123456789;
  x = r8vec_uniform_01_new ( n, &seed );
  for ( i = 0; i < n; i++ )
  {
    x[i] = r8_nint ( 50.0 * x[i] - 25.0 ) / 5.0;
  }
  a = vand2 ( n, x );
  determ1 = vand2_determinant ( n, x );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  VAND2               %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
  free ( x );
/*
  WATHEN matrix.
*/
  if ( 0 )
  {
  row_num = 5;
  col_num = 5;
  n = wathen_order ( row_num, col_num );
  a = wathen ( row_num, col_num, n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  WATHEN              %4d                  %14g  %14f\n", 
    n,          determ2, norm_frobenius );
  free ( a );
  }
  else
  {
  printf ( "  WATHEN             -----Not ready-----\n" );
  }
/*
  WILK03 matrix.
*/
  n = 3;
  a = wilk03 ( );
  determ1 = wilk03_determinant ( );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  WILK03              %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  WILK04 matrix.
*/
  n = 4;
  a = wilk04 ( );
  determ1 = wilk04_determinant ( );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  WILK04              %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  WILK05 matrix.
*/
  n = 5;
  a = wilk05 ( );
  determ1 = wilk05_determinant ( );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  WILK05              %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  WILK12 matrix.
*/
  n = 12;
  a = wilk12 ( );
  determ1 = wilk12_determinant ( );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  WILK12              %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  WILK20 matrix.
*/
  n = 20;
  seed = 123456789;
  alpha = r8_uniform_01 ( &seed );
  alpha = r8_nint ( 50 * alpha - 25.0 ) / 5.0;
  a = wilk20 ( alpha );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  WILK20              %4d                  %14g  %14f\n", 
    n,          determ2, norm_frobenius );
  free ( a );
/*
  WILK21 matrix.
*/
  n = 21;
  a = wilk21 ( n );
  determ1 = wilk21_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  WILK21              %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  WILSON matrix.
*/
  n = 4;
  a = wilson ( );
  determ1 = wilson_determinant ( );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  WILSON              %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  ZERO matrix.
*/
  n = 5;
  a = zero ( n, n );
  determ1 = zero_determinant ( n );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  ZERO                %4d  %14g  %14g  %14f\n", 
    n, determ1, determ2, norm_frobenius );
  free ( a );
/*
  ZIELKE matrix.
*/
  n = 5;
  seed = 123456789;
  d1 = r8_uniform_01 ( &seed );
  d1 = r8_nint ( 50.0 * d1 - 25.0 ) / 5.0;
  d2 = r8_uniform_01 ( &seed );
  d2 = r8_nint ( 50.0 * d2 - 25.0 ) / 5.0;
  d3 = r8_uniform_01 ( &seed );
  d3 = r8_nint ( 50.0 * d3 - 25.0 ) / 5.0;
  a = zielke ( n, d1, d2, d3 );
  determ2 = r8mat_determinant ( n, a );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  ZIELKE              %4d                  %14g  %14f\n", 
    n,          determ2, norm_frobenius );
  free ( a );

  return;
}
/******************************************************************************/

void test_eigen ( void )

/******************************************************************************/
/*
  Purpose:

    TEST_EIGEN tests the eigenvalue computations.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    03 July 2011

  Author:

    John Burkardt
*/
{
  double *a;
  double alpha;
  double beta;
  double error_frobenius;
  double gamma;
  int i;
  int i1;
  int k;
  double *lambda;
  int n;
  double norm_frobenius;
  int rank;
  int seed;
  int seed_save;
  char title[21];
  double *v1;
  double *v2;
  double *v3;
  double *x;

  printf ( "\n" );
  printf ( "TEST_EIGEN\n" );
  printf ( "  Compute the Frobenius norm of the eigenvalue error:\n" );
  printf ( "    A * X - X * LAMBDA\n" );
  printf ( "  given a set of K eigenvectors X and eigenvalues LAMBDA.\n" );
  printf ( "\n" );
  printf ( "  Matrix title             N     K      ||A||          ||(A-Lambda*I)*X||\n" );
  printf ( "\n" );
//
//  BODEWIG matrix.
//
  n = 4;
  k = 4;
  a = bodewig ( );
  lambda = bodewig_eigenvalues ( );
  x = bodewig_right ( );
  error_frobenius = r8mat_is_eigen_right ( n, k, a, x, lambda );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  BODEWIG               %4d  %4d  %14g  %14g\n",
    n, k, norm_frobenius, error_frobenius );
  free ( a );
  free ( lambda );
  free ( x );
/*
  CARRY matrix.
*/
  n = 5;
  k = 5;
  seed = 123456789;
  i1 = i4_uniform ( 2, 20, &seed );
  a = carry ( i1, n );
  lambda = carry_eigenvalues ( i1, n );
  x = carry_right ( n );
  error_frobenius = r8mat_is_eigen_right ( n, k, a, x, lambda );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  CARRY                 %4d  %4d  %14g  %14g\n",
    n, k, norm_frobenius, error_frobenius );
  free ( a );
  free ( lambda );
  free ( x );
/*
  CHOW matrix.
*/
  n = 5;
  k = 5;
  seed = 123456789;
  alpha = r8_uniform_01 ( &seed );
  alpha = r8_nint ( 50.0 * alpha ) / 5.0;
  beta = r8_uniform_01 ( &seed );
  beta = r8_nint ( 50.0 * beta ) / 5.0;
  a = chow ( alpha, beta, n, n );
  lambda = chow_eigenvalues ( alpha, beta, n );
  x = chow_right ( alpha, beta, n );
  error_frobenius = r8mat_is_eigen_right ( n, k, a, x, lambda );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  CHOW                  %4d  %4d  %14g  %14g\n",
    n, k, norm_frobenius, error_frobenius );
  free ( a );
  free ( lambda );
  free ( x );
/*
  COMBIN matrix.
*/
  n = 5;
  k = 5;
  seed = 123456789;
  alpha = r8_uniform_01 ( &seed );
  alpha = r8_nint ( 50.0 * alpha ) / 5.0;
  beta = r8_uniform_01 ( &seed );
  beta = r8_nint ( 50.0 * beta ) / 5.0;
  a = combin ( alpha, beta, n );
  lambda = combin_eigenvalues ( alpha, beta, n );
  x = combin_right ( alpha, beta, n );
  error_frobenius = r8mat_is_eigen_right ( n, k, a, x, lambda );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  COMBIN                %4d  %4d  %14g  %14g\n",
    n, k, norm_frobenius, error_frobenius );
  free ( a );
  free ( lambda );
  free ( x );
/*
  DIF2 matrix.
*/
  n = 5;
  k = 5;
  a = dif2 ( n, n );
  lambda = dif2_eigenvalues ( n );
  x = dif2_right ( n );
  error_frobenius = r8mat_is_eigen_right ( n, k, a, x, lambda );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  DIF2                  %4d  %4d  %14g  %14g\n",
    n, k, norm_frobenius, error_frobenius );
  free ( a );
  free ( lambda );
  free ( x );
/*
  EXCHANGE matrix.
*/
  n = 5;
  k = 5;
  a = exchange ( n, n );
  lambda = exchange_eigenvalues ( n );
  x = exchange_right ( n );
  error_frobenius = r8mat_is_eigen_right ( n, k, a, x, lambda );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  EXCHANGE              %4d  %4d  %14g  %14g\n",
    n, k, norm_frobenius, error_frobenius );
  free ( a );
  free ( lambda );
  free ( x );
/*
  IDEM_RANDOM matrix.
*/
  n = 5;
  k = 5;
  rank = 3;
  seed_save = 987654321;
  seed = seed_save;
  a = idem_random ( n, rank, &seed );
  lambda = idem_random_eigenvalues ( n, rank );
  seed = seed_save;
  x = idem_random_right ( n, rank, &seed );
  error_frobenius = r8mat_is_eigen_right ( n, k, a, x, lambda );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  IDEM_RANDOM           %4d  %4d  %14g  %14g\n",
    n, k, norm_frobenius, error_frobenius );
  free ( a );
  free ( lambda );
  free ( x );
/*
  IDENTITY matrix.
*/
  n = 5;
  k = 5;
  a = identity ( n, n );
  lambda = identity_eigenvalues ( n );
  x = identity_right ( n );
  error_frobenius = r8mat_is_eigen_right ( n, k, a, x, lambda );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  IDENTITY              %4d  %4d  %14g  %14g\n",
    n, k, norm_frobenius, error_frobenius );
  free ( a );
  free ( lambda );
  free ( x );
/*
  ILL3 matrix.
*/
  n = 3;
  k = 3;
  a = ill3 ( );
  lambda = ill3_eigenvalues ( );
  x = ill3_right ( );
  error_frobenius = r8mat_is_eigen_right ( n, k, a, x, lambda );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  ILL3                  %4d  %4d  %14g  %14g\n",
    n, k, norm_frobenius, error_frobenius );
  free ( a );
  free ( lambda );
  free ( x );
/*
  KMS matrix.
*/
  n = 5;
  k = 5;
  seed = 123456789;
  alpha = r8_uniform_01 ( &seed );
  a = kms ( alpha, n, n );
  lambda = kms_eigenvalues ( alpha, n );
  x = kms_right ( alpha, n );
  error_frobenius = r8mat_is_eigen_right ( n, k, a, x, lambda );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  KMS                   %4d  %4d  %14g  %14g\n",
    n, k, norm_frobenius, error_frobenius );
  free ( a );
  free ( lambda );
  free ( x );
/*
  ONE matrix.
*/
  n = 5;
  k = 5;
  a = one ( n, n );
  lambda = one_eigenvalues ( n );
  x = one_right ( n );
  error_frobenius = r8mat_is_eigen_right ( n, k, a, x, lambda );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  ONE                   %4d  %4d  %14g  %14g\n",
    n, k, norm_frobenius, error_frobenius );
  free ( a );
  free ( lambda );
  free ( x );
/*
  ORTEGA matrix.
*/
  n = 5;
  k = 5;
  seed = 123456789;
  v1 = r8vec_uniform_01_new ( n, &seed );
  for ( i = 0; i < n; i++ )
  {
    v1[i] = r8_nint ( 50.0 * v1[i] - 25.0 ) / 5.0;
  }
  v2 = r8vec_uniform_01_new ( n, &seed );
  for ( i = 0; i < n; i++ )
  {
    v2[i] = r8_nint ( 50.0 * v2[i] - 25.0 ) / 5.0;
  }
  v3 = r8vec_uniform_01_new ( n, &seed );
  for ( i = 0; i < n; i++ )
  {
    v3[i] = r8_nint ( 50.0 * v3[i] - 25.0 ) / 5.0;
  }
  a = ortega ( n, v1, v2, v3 );
  lambda = ortega_eigenvalues ( n, v1, v2, v3 );
  x = ortega_right ( n, v1, v2, v3 );
  error_frobenius = r8mat_is_eigen_right ( n, k, a, x, lambda );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  ORTEGA                %4d  %4d  %14g  %14g\n",
    n, k, norm_frobenius, error_frobenius );
  free ( a );
  free ( lambda );
  free ( v1 );
  free ( v2 );
  free ( v3 );
  free ( x );
/*
  OTO matrix.
*/
  n = 5;
  k = 5;
  a = oto ( n, n );
  lambda = oto_eigenvalues ( n );
  x = oto_right ( n );
  error_frobenius = r8mat_is_eigen_right ( n, k, a, x, lambda );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  OTO                   %4d  %4d  %14g  %14g\n",
    n, k, norm_frobenius, error_frobenius );
  free ( a );
  free ( lambda );
  free ( x );
/*
  PDS_RANDOM matrix.
*/
  n = 5;
  k = 5;
  seed_save = 123456789;
  seed = seed_save;
  a = pds_random ( n, &seed );
  seed = seed_save;
  lambda = pds_random_eigenvalues ( n, &seed );
  seed = seed_save;
  x = pds_random_right ( n, &seed );
  error_frobenius = r8mat_is_eigen_right ( n, k, a, x, lambda );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  PDS_RANDOM            %4d  %4d  %14g  %14g\n",
    n, k, norm_frobenius, error_frobenius );
  free ( a );
  free ( lambda );
  free ( x );
/*
  PEI matrix.
*/
  n = 5;
  k = 5;
  seed = 123456789;
  alpha = r8_uniform_01 ( &seed );
  alpha = r8_nint ( 50.0 * alpha ) / 5.0;
  a = pei ( alpha, n );
  lambda = pei_eigenvalues ( alpha, n );
  x = pei_right ( alpha, n );
  error_frobenius = r8mat_is_eigen_right ( n, k, a, x, lambda );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  PEI                   %4d  %4d  %14g  %14g\n",
    n, k, norm_frobenius, error_frobenius );
  free ( a );
  free ( lambda );
  free ( x );
/*
  RODMAN matrix.
*/
  n = 5;
  k = 5;
  seed = 123456789;
  alpha = r8_uniform_01 ( &seed );
  alpha = r8_nint ( 50.0 * alpha ) / 5.0;
  a = rodman ( alpha, n, n );
  lambda = rodman_eigenvalues ( alpha, n );
  x = rodman_right ( alpha, n );
  error_frobenius = r8mat_is_eigen_right ( n, k, a, x, lambda );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  RODMAN                %4d  %4d  %14g  %14g\n",
    n, k, norm_frobenius, error_frobenius );
  free ( a );
  free ( lambda );
  free ( x );
/*
  ROSSER1 matrix.
*/
  n = 8;
  k = 8;
  a = rosser1 ( );
  lambda = rosser1_eigenvalues ( );
  x = rosser1_right ( );
  error_frobenius = r8mat_is_eigen_right ( n, k, a, x, lambda );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  ROSSER1               %4d  %4d  %14g  %14g\n",
    n, k, norm_frobenius, error_frobenius );
  free ( a );
  free ( lambda );
  free ( x );
/*
  RUTIS1 matrix.
*/
  n = 4;
  k = 4;
  a = rutis1 ( );
  lambda = rutis1_eigenvalues ( );
  x = rutis1_right ( );
  error_frobenius = r8mat_is_eigen_right ( n, k, a, x, lambda );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  RUTIS1                %4d  %4d  %14g  %14g\n",
    n, k, norm_frobenius, error_frobenius );
  free ( a );
  free ( lambda );
  free ( x );
/*
  RUTIS2 matrix.
*/
  n = 4;
  k = 4;
  a = rutis2 ( );
  lambda = rutis2_eigenvalues ( );
  x = rutis2_right ( );
  error_frobenius = r8mat_is_eigen_right ( n, k, a, x, lambda );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  RUTIS2                %4d  %4d  %14g  %14g\n",
    n, k, norm_frobenius, error_frobenius );
  free ( a );
  free ( lambda );
  free ( x );
/*
  RUTIS5 matrix.
*/
  n = 4;
  k = 4;
  a = rutis5 ( );
  lambda = rutis5_eigenvalues ( );
  x = rutis5_right ( );
  error_frobenius = r8mat_is_eigen_right ( n, k, a, x, lambda );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  RUTIS5                %4d  %4d  %14g  %14g\n",
    n, k, norm_frobenius, error_frobenius );
  free ( a );
  free ( lambda );
  free ( x );
/*
  WILK12 matrix.
*/
  n = 12;
  k = 12;
  a = wilk12 ( );
  lambda = wilk12_eigenvalues ( );
  x = wilk12_right ( );
  error_frobenius = r8mat_is_eigen_right ( n, k, a, x, lambda );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  WILK12                %4d  %4d  %14g  %14g\n",
    n, k, norm_frobenius, error_frobenius );
  free ( a );
  free ( lambda );
  free ( x );
/*
  WILSON matrix.
*/
  n = 4;
  k = 4;
  a = wilson ( );
  lambda = wilson_eigenvalues ( );
  x = wilson_right ( );
  error_frobenius = r8mat_is_eigen_right ( n, k, a, x, lambda );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  WILSON                %4d  %4d  %14g  %14g\n",
    n, k, norm_frobenius, error_frobenius );
  free ( a );
  free ( lambda );
  free ( x );
/*
  ZERO matrix.
*/
  n = 5;
  k = 5;
  a = zero ( n, n );
  lambda = zero_eigenvalues ( n );
  x = zero_right ( n );
  error_frobenius = r8mat_is_eigen_right ( n, k, a, x, lambda );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  ZERO                  %4d  %4d  %14g  %14g\n",
    n, k, norm_frobenius, error_frobenius );
  free ( a );
  free ( lambda );
  free ( x );

  return;
}
/******************************************************************************/

void test_inverse ( void )

/******************************************************************************/
/*
  Purpose:

    TEST_INVERSE tests the inverse computations.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 July 2011

  Author:

    John Burkardt
*/
{
  double *a;
  double alpha;
  double *b;;
  double beta;
  double *c;
  double error_ab;
  double error_ac;;
  double gamma;
  int i;
  int ii;;
  int jj;
  int k;
  double *l;
  int n;
  double norma_frobenius;
  double normc_frobenius;
  double *p;
  int *pivot;
  int seed;
  int seed_save;
  char title[21];
  double *u;
  double *v1;
  double *v2;
  double *v3;
  double *w;
  double *x;
  int x_n;
  double *y;
  int y_n;
  double *z;

  printf ( "\n" );
  printf ( "TEST_INVERSE\n" );
  printf ( "  A = a test matrix of order N;\n" );
  printf ( "  B = inverse as computed by a routine.\n" );
  printf ( "  C = inverse as computed by R8MAT_INVERSE.\n" );
  printf ( "\n" );
  printf ( "  ||I-AB|| = Frobenius norm of I-A*B.\n" );
  printf ( "  ||I-AC|| = Frobenius norm of I-A*C.\n" );
  printf ( "  ||I-AB|| = Frobenius norm of I-A*B.\n" );
  printf ( "\n" );
  printf ( "  Matrix title             N        " );
  printf ( "   ||A||          ||C||      ||I-AC||        ||I-AB||\n" );
  printf ( "\n" );
/*
  AEGERTER matrix.
*/
  strcpy ( title, "AEGERTER            " );
  n = 5;
  a = aegerter ( n );
  b = aegerter_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  printf ( "  %20s  %4d  %14g  %14g  %14g  %14g\n",
    title, n, norma_frobenius, normc_frobenius, error_ac, error_ab );
  free ( a );
  free ( b );
  free ( c );
/*
  BAB matrix.
*/
  n = 5;
  seed = 123456789;
  alpha = r8_nint ( 50.0 * r8_uniform_01 ( &seed ) - 25.0 ) / 5.0;
  beta = r8_nint ( 50.0 * r8_uniform_01 ( &seed ) - 25.0 ) / 5.0;
  a = bab ( n, alpha, beta );
  b = bab_inverse ( n, alpha, beta );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  printf ( "  BAB                   %4d  %14g  %14g  %14g  %14g\n",
    n, norma_frobenius, normc_frobenius, error_ac, error_ab );
  free ( a );
  free ( b );
  free ( c );
/*
  BERNSTEIN matrix.
*/
  n = 5;
  a = bernstein ( n );
  b = bernstein_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  printf ( "  BERNSTEIN             %4d  %14g  %14g  %14g  %14g\n",
    n, norma_frobenius, normc_frobenius, error_ac, error_ab );
  free ( a );
  free ( b );
  free ( c );
/*
  BIS matrix.
*/
  n = 5;
  seed = 123456789;
  alpha = r8_nint ( 50.0 * r8_uniform_01 ( &seed ) - 25.0 ) / 5.0;
  beta = r8_nint ( 50.0 * r8_uniform_01 ( &seed ) - 25.0 ) / 5.0;
  a = bis ( alpha, beta, n, n );
  b = bis_inverse ( alpha, beta, n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  printf ( "  BIS                   %4d  %14g  %14g  %14g  %14g\n",
    n, norma_frobenius, normc_frobenius, error_ac, error_ab );
  free ( a );
  free ( b );
  free ( c );
/*
  BODEWIG matrix.
*/
  n = 4;
  a = bodewig ( );
  b = bodewig_inverse ( );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  printf ( "  BODEWIG               %4d  %14g  %14g  %14g  %14g\n",
    n, norma_frobenius, normc_frobenius, error_ac, error_ab );
  free ( a );
  free ( b );
  free ( c );
/*
  BOOTHROYD matrix.
*/
  n = 5;
  a = boothroyd ( n );
  b = boothroyd_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  printf ( "  BOOTHROYD             %4d  %14g  %14g  %14g  %14g\n",
    n, norma_frobenius, normc_frobenius, error_ac, error_ab );
  free ( a );
  free ( b );
  free ( c );
/*
  BORDERBAND matrix.
*/
  n = 5;
  a = borderband ( n );
  b = borderband_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  printf ( "  BORDERBAND            %4d  %14g  %14g  %14g  %14g\n",
    n, norma_frobenius, normc_frobenius, error_ac, error_ab );
  free ( a );
  free ( b );
  free ( c );
/*
  CARRY matrix.
*/
  n = 5;
  seed = 123456789;
  k = i4_uniform ( 2, 20, &seed );
  a = carry ( k, n );
  b = carry_inverse ( k, n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  printf ( "  CARRY                 %4d  %14g  %14g  %14g  %14g\n",
    n, norma_frobenius, normc_frobenius, error_ac, error_ab );
  free ( a );
  free ( b );
  free ( c );
/*
  CAUCHY matrix.
*/
  n = 5;
  seed = 123456789;
  x = r8vec_uniform_01_new ( n, &seed );
  y = r8vec_uniform_01_new ( n, &seed );
  a = cauchy ( n, x, y );
  b = cauchy_inverse ( n, x, y );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  printf ( "  CAUCHY                %4d  %14g  %14g  %14g  %14g\n",
    n, norma_frobenius, normc_frobenius, error_ac, error_ab );
  free ( a );
  free ( b );
  free ( c );
  free ( x );
  free ( y );
/*
  CHEBY_T matrix.
*/
  n = 5;
  a = cheby_t ( n );
  b = cheby_t_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  printf ( "  CHEBY_T               %4d  %14g  %14g  %14g  %14g\n",
    n, norma_frobenius, normc_frobenius, error_ac, error_ab );
  free ( a );
  free ( b );
  free ( c );
/*
  CHEBY_U matrix.
*/
  n = 5;
  a = cheby_u ( n );
  b = cheby_u_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  printf ( "  CHEBY_U               %4d  %14g  %14g  %14g  %14g\n",
    n, norma_frobenius, normc_frobenius, error_ac, error_ab );
  free ( a );
  free ( b );
  free ( c );
/*
  CHEBY_VAN2 matrix.
*/
  n = 5;
  a = cheby_van2 ( n );
  b = cheby_van2_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  printf ( "  CHEBY_VAN2            %4d  %14g  %14g  %14g  %14g\n",
    n, norma_frobenius, normc_frobenius, error_ac, error_ab );
  free ( a );
  free ( b );
  free ( c );
/*
  CHEBY_VAN3 matrix.
*/
  n = 5;
  a = cheby_van3 ( n );
  b = cheby_van3_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  printf ( "  CHEBY_VAN3            %4d  %14g  %14g  %14g  %14g\n",
    n, norma_frobenius, normc_frobenius, error_ac, error_ab );
  free ( a );
  free ( b );
  free ( c );
/*
  CHOW matrix.
*/
  n = 5;
  seed = 123456789;
  alpha = r8_uniform_01 ( &seed );
  alpha = r8_nint ( 50.0 * alpha ) / 5.0;
  beta = r8_uniform_01 ( &seed );
  beta = r8_nint ( 50.0 * beta ) / 5.0;
  a = chow ( alpha, beta, n, n );
  b = chow_inverse ( alpha, beta, n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  printf ( "  CHOW                  %4d  %14g  %14g  %14g  %14g\n",
    n, norma_frobenius, normc_frobenius, error_ac, error_ab );
  free ( a );
  free ( b );
  free ( c );
/*
  CIRCULANT matrix.
*/
  n = 5;
  seed = 123456789;
  x = r8vec_uniform_01_new ( n, &seed );
  for ( i = 0; i < n; i++ )
  {
    x[i] = r8_nint ( 50.0 * x[i] - 25.0 ) / 5.0;
  }
  a = circulant ( n, n, x );
  b = circulant_inverse ( n, x );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  printf ( "  CIRCULANT             %4d  %14g  %14g  %14g  %14g\n",
    n, norma_frobenius, normc_frobenius, error_ac, error_ab );
  free ( a );
  free ( b );
  free ( c );
  free ( x );
/*
  CIRCULANT2 matrix.
*/
  if ( 1 )
  {
  n = 5;
  a = circulant2 ( n );
  b = circulant2_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  printf ( "  CIRCULANT2            %4d  %14g  %14g  %14g  %14g\n",
    n, norma_frobenius, normc_frobenius, error_ac, error_ab );
  free ( a );
  free ( b );
  free ( c );
  }
  else
  {
    printf ( "  CIRCULANT2 ---- Not ready!\n" );
  }
/*
  CLEMENT1 matrix.
*/
  n = 6;
  a = clement1 ( n );
  b = clement1_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  printf ( "  CLEMENT1              %4d  %14g  %14g  %14g  %14g\n",
    n, norma_frobenius, normc_frobenius, error_ac, error_ab );
  free ( a );
  free ( b );
  free ( c );
/*
  CLEMENT2 matrix.
*/
  n = 6;
  a = clement2 ( n );
  b = clement2_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  printf ( "  CLEMENT2              %4d  %14g  %14g  %14g  %14g\n",
    n, norma_frobenius, normc_frobenius, error_ac, error_ab );
  free ( a );
  free ( b );
  free ( c );
/*
  CLEMENT3.
*/
  n = 6;
  seed = 123456789;
  x = r8vec_uniform_01_new ( n - 1, &seed );
  for ( i = 0; i < n - 1; i++ )
  {
    x[i] = r8_nint ( 50.0 * x[i] - 25.0 ) / 5.0;
  }
  y = r8vec_uniform_01_new ( n - 1, &seed );
  for ( i = 0; i < n - 1; i++ )
  {
    y[i] = r8_nint ( 50.0 * y[i] - 25.0 ) / 5.0;
  }
  a = clement3 ( n, x, y );
  b = clement3_inverse ( n, x, y );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  printf ( "  CLEMENT3              %4d  %14g  %14g  %14g  %14g\n",
    n, norma_frobenius, normc_frobenius, error_ac, error_ab );
  free ( a );
  free ( b );
  free ( c );
  free ( x );
  free ( y );
/*
  COMBIN matrix.
*/
  n = 5;
  seed = 123456789;
  alpha = r8_uniform_01 ( &seed );
  alpha = r8_nint ( 50.0 * alpha ) / 5.0;
  beta = r8_uniform_01 ( &seed );
  beta = r8_nint ( 50.0 * beta ) / 5.0;
  a = combin ( alpha, beta, n );
  b = combin_inverse ( alpha, beta, n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  printf ( "  COMBIN                %4d  %14g  %14g  %14g  %14g\n",
    n, norma_frobenius, normc_frobenius, error_ac, error_ab );
  free ( a );
  free ( b );
  free ( c );
/*
  COMPANION.
*/
  n = 5;
  seed = 123456789;
  x = r8vec_uniform_01_new ( n, &seed );
  for ( i = 0; i < n; i++ )
  {
    x[i] = r8_nint ( 10.0 * x[i] - 5.0 );
  }
  a = companion ( n, x );
  b = companion_inverse ( n, x );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  printf ( "  COMPANION             %4d  %14g  %14g  %14g  %14g\n",
    n, norma_frobenius, normc_frobenius, error_ac, error_ab );
  free ( a );
  free ( b );
  free ( c );
  free ( x );
/*
  COMPLEX_I
*/
  n = 2;
  a = complex_i ( );
  b = complex_i_inverse ( );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  printf ( "  COMPLEX_I             %4d  %14g  %14g  %14g  %14g\n",
    n, norma_frobenius, normc_frobenius, error_ac, error_ab );
  free ( a );
  free ( b );
  free ( c );
/*
  CONEX1 matrix.
*/
  n = 4;
  seed = 123456789;
  alpha = r8_uniform_01 ( &seed );
  alpha = r8_nint ( 50.0 * alpha ) / 5.0;
  a = conex1 ( alpha );
  b = conex1_inverse ( alpha );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  printf ( "  CONEX1                %4d  %14g  %14g  %14g  %14g\n",
    n, norma_frobenius, normc_frobenius, error_ac, error_ab );
  free ( a );
  free ( b );
  free ( c );
/*
  CONEX2 matrix.
*/
  n = 3;
  seed = 123456789;
  alpha = r8_uniform_01 ( &seed );
  alpha = r8_nint ( 50.0 * alpha ) / 5.0;
  a = conex2 ( alpha );
  b = conex2_inverse ( alpha );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  printf ( "  CONEX2                %4d  %14g  %14g  %14g  %14g\n",
    n, norma_frobenius, normc_frobenius, error_ac, error_ab );
  free ( a );
  free ( b );
  free ( c );
/*
  CONEX3 matrix.
*/
  n = 5;
  a = conex3 ( n );
  b = conex3_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  printf ( "  CONEX3                %4d  %14g  %14g  %14g  %14g\n",
    n, norma_frobenius, normc_frobenius, error_ac, error_ab );
  free ( a );
  free ( b );
  free ( c );
/*
  CONFERENCE matrix.
*/
  n = 6;
  a = conference ( n );
  b = conference_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  printf ( "  CONFERENCE            %4d  %14g  %14g  %14g  %14g\n",
    n, norma_frobenius, normc_frobenius, error_ac, error_ab );
  free ( a );
  free ( b );
  free ( c );
/*
  DAUB2 matrix.
*/
  n = 4;
  a = daub2 ( n );
  b = daub2_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  printf ( "  DAUB2                 %4d  %14g  %14g  %14g  %14g\n",
    n, norma_frobenius, normc_frobenius, error_ac, error_ab );
  free ( a );
  free ( b );
  free ( c );
/*
  DAUB4 matrix.
*/
  n = 8;
  a = daub4 ( n );
  b = daub4_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  printf ( "  DAUB4                 %4d  %14g  %14g  %14g  %14g\n",
    n, norma_frobenius, normc_frobenius, error_ac, error_ab );
  free ( a );
  free ( b );
  free ( c );
/*
  DAUB6 matrix.
*/
  n = 12;
  a = daub6 ( n );
  b = daub6_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  printf ( "  DAUB6                 %4d  %14g  %14g  %14g  %14g\n",
    n, norma_frobenius, normc_frobenius, error_ac, error_ab );
  free ( a );
  free ( b );
  free ( c );
/*
  DAUB8 matrix.
*/
  n = 16;
  a = daub8 ( n );
  b = daub8_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  printf ( "  DAUB8                 %4d  %14g  %14g  %14g  %14g\n",
    n, norma_frobenius, normc_frobenius, error_ac, error_ab );
  free ( a );
  free ( b );
  free ( c );
/*
  DAUB10 matrix.
*/
  n = 20;
  a = daub10 ( n );
  b = daub10_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  printf ( "  DAUB10                %4d  %14g  %14g  %14g  %14g\n",
    n, norma_frobenius, normc_frobenius, error_ac, error_ab );
  free ( a );
  free ( b );
  free ( c );
/*
  DAUB12 matrix.
*/
  n = 24;
  a = daub12 ( n );
  b = daub12_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  printf ( "  DAUB12                %4d  %14g  %14g  %14g  %14g\n",
    n, norma_frobenius, normc_frobenius, error_ac, error_ab );
  free ( a );
  free ( b );
  free ( c );
/*
  DIAGONAL.
*/
  n = 5;
  seed = 123456789;
  x = r8vec_uniform_01_new ( n, &seed );
  for ( i = 0; i < n; i++ )
  {
    x[i] = r8_nint ( 50.0 * x[i] - 25.0 ) / 5.0;
  }
  a = diagonal ( n, n, x );
  b = diagonal_inverse ( n, x );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  printf ( "  DIAGONAL              %4d  %14g  %14g  %14g  %14g\n",
    n, norma_frobenius, normc_frobenius, error_ac, error_ab );
  free ( a );
  free ( b );
  free ( c );
  free ( x );
/*
  DIF2 matrix.
*/
  n = 5;
  a = dif2 ( n, n );
  b = dif2_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  printf ( "  DIF2                  %4d  %14g  %14g  %14g  %14g\n",
    n, norma_frobenius, normc_frobenius, error_ac, error_ab );
  free ( a );
  free ( b );
  free ( c );
/*
  DOWNSHIFT matrix.
*/
  n = 5;
  a = downshift ( n );
  b = downshift_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  printf ( "  DOWNSHIFT             %4d  %14g  %14g  %14g  %14g\n",
    n, norma_frobenius, normc_frobenius, error_ac, error_ab );
  free ( a );
  free ( b );
  free ( c );
/*
  DRMAC
*/

/*
  EULERIAN matrix.
*/
  n = 5;
  a = eulerian ( n, n );
  b = eulerian_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  printf ( "  EULERIAN              %4d  %14g  %14g  %14g  %14g\n",
    n, norma_frobenius, normc_frobenius, error_ac, error_ab );
  free ( a );
  free ( b );
  free ( c );
/*
  EXCHANGE matrix.
*/
  n = 5;
  a = exchange ( n, n );
  b = exchange_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  printf ( "  EXCHANGE              %4d  %14g  %14g  %14g  %14g\n",
    n, norma_frobenius, normc_frobenius, error_ac, error_ab );
  free ( a );
  free ( b );
  free ( c );
/*
  FIBONACCI2 matrix.
*/
  n = 5;
  a = fibonacci2 ( n );
  b = fibonacci2_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  printf ( "  FIBONACCI2            %4d  %14g  %14g  %14g  %14g\n",
    n, norma_frobenius, normc_frobenius, error_ac, error_ab );
  free ( a );
  free ( b );
  free ( c );
/*
  FIBONACCI3 matrix.
*/
  n = 5;
  a = fibonacci3 ( n );
  b = fibonacci3_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  printf ( "  FIBONACCI3            %4d  %14g  %14g  %14g  %14g\n",
    n, norma_frobenius, normc_frobenius, error_ac, error_ab );
  free ( a );
  free ( b );
  free ( c );
/*
  FIEDLER.
  The FIEDLER_INVERSE routine assumes the X vector is sorted.
*/
  n = 7;
  seed = 123456789;
  x = r8vec_uniform_01_new ( n, &seed );
  for ( i = 0; i < n; i++ )
  {
    x[i] = r8_nint ( 50.0 * x[i] - 25.0 ) / 5.0;
  }
  r8vec_sort_bubble_a ( n, x );
  a = fiedler ( n, n, x );
  b = fiedler_inverse ( n, x );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  printf ( "  FIEDLER               %4d  %14g  %14g  %14g  %14g\n",
    n, norma_frobenius, normc_frobenius, error_ac, error_ab );
  free ( a );
  free ( b );
  free ( c );
  free ( x );
/*
  FORSYTHE matrix.
*/
  n = 5;
  seed = 123456789;
  alpha = r8_uniform_01 ( &seed );
  alpha = r8_nint ( 50.0 * alpha ) / 5.0;
  beta = r8_uniform_01 ( &seed );
  beta = r8_nint ( 50.0 * beta ) / 5.0;
  a = forsythe ( alpha, beta, n );
  b = forsythe_inverse ( alpha, beta, n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  printf ( "  FORSYTHE              %4d  %14g  %14g  %14g  %14g\n",
    n, norma_frobenius, normc_frobenius, error_ac, error_ab );
  free ( a );
  free ( b );
  free ( c );
/*
  FOURIER_COSINE matrix.
*/
  n = 5;
  a = fourier_cosine ( n );
  b = fourier_cosine_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  printf ( "  FOURIER_COSINE        %4d  %14g  %14g  %14g  %14g\n",
    n, norma_frobenius, normc_frobenius, error_ac, error_ab );
  free ( a );
  free ( b );
  free ( c );
/*
  FOURIER_SINE matrix.
*/
  n = 5;
  a = fourier_sine ( n );
  b = fourier_sine_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  printf ( "  FOURIER_SINE          %4d  %14g  %14g  %14g  %14g\n",
    n, norma_frobenius, normc_frobenius, error_ac, error_ab );
  free ( a );
  free ( b );
  free ( c );
/*
  FRANK matrix.
*/
  n = 5;
  a = frank ( n );
  b = frank_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  printf ( "  FRANK                 %4d  %14g  %14g  %14g  %14g\n",
    n, norma_frobenius, normc_frobenius, error_ac, error_ab );
  free ( a );
  free ( b );
  free ( c );
/*
  GFPP matrix.
*/
  n = 5;
  seed = 123456789;
  alpha = r8_uniform_01 ( &seed );
  a = gfpp ( n, alpha );
  b = gfpp_inverse ( n, alpha );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  printf ( "  GFPP                  %4d  %14g  %14g  %14g  %14g\n",
    n, norma_frobenius, normc_frobenius, error_ac, error_ab );
  free ( a );
  free ( b );
  free ( c );
/*
  GIVENS matrix.
*/
  n = 5;
  a = givens ( n, n );
  b = givens_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  printf ( "  GIVENS                %4d  %14g  %14g  %14g  %14g\n",
    n, norma_frobenius, normc_frobenius, error_ac, error_ab );
  free ( a );
  free ( b );
  free ( c );
/*
  GK316 matrix.
*/
  n = 5;
  a = gk316 ( n );
  b = gk316_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  printf ( "  GK316                 %4d  %14g  %14g  %14g  %14g\n",
    n, norma_frobenius, normc_frobenius, error_ac, error_ab );
  free ( a );
  free ( b );
  free ( c );
/*
  GK323 matrix.
*/
  n = 5;
  a = gk323 ( n, n );
  b = gk323_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  printf ( "  GK323                 %4d  %14g  %14g  %14g  %14g\n",
    n, norma_frobenius, normc_frobenius, error_ac, error_ab );
  free ( a );
  free ( b );
  free ( c );
/*
  GK324 matrix.
*/
  n = 5;
  seed = 123456789;
  x = r8vec_uniform_01_new ( n, &seed );
  for ( i = 0; i < n; i++ )
  {
    x[i] = r8_nint ( 50.0 * x[i] - 25.0 ) / 5.0;
  }
  a = gk324 ( n, n, x );
  b = gk324_inverse ( n, x );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  printf ( "  GK324                 %4d  %14g  %14g  %14g  %14g\n",
    n, norma_frobenius, normc_frobenius, error_ac, error_ab );
  free ( a );
  free ( b );
  free ( c );
  free ( x );
/*
  HANOWA matrix.
*/
  n = 8;
  seed = 123456789;
  alpha = r8_uniform_01 ( &seed );
  a = hanowa ( alpha, n );
  b = hanowa_inverse ( alpha, n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  printf ( "  HANOWA                %4d  %14g  %14g  %14g  %14g\n",
    n, norma_frobenius, normc_frobenius, error_ac, error_ab );
  free ( a );
  free ( b );
  free ( c );
/*
  HARMAN matrix.
*/
  n = 8;
  a = harman (  );
  b = harman_inverse (  );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  printf ( "  HARMAN                %4d  %14g  %14g  %14g  %14g\n",
    n, norma_frobenius, normc_frobenius, error_ac, error_ab );
  free ( a );
  free ( b );
  free ( c );
/*
  HARTLEY matrix.
*/
  n = 5;
  a = hartley ( n );
  b = hartley_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  printf ( "  HARTLEY               %4d  %14g  %14g  %14g  %14g\n",
    n, norma_frobenius, normc_frobenius, error_ac, error_ab );
  free ( a );
  free ( b );
  free ( c );
/*
  HELMERT matrix.
*/
  n = 5;
  a = helmert ( n );
  b = helmert_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  printf ( "  HELMERT               %4d  %14g  %14g  %14g  %14g\n",
    n, norma_frobenius, normc_frobenius, error_ac, error_ab );
  free ( a );
  free ( b );
  free ( c );
/*
  HELMERT2 matrix.
*/
  n = 5;
  seed = 123456789;
  x = r8vec_uniform_01_new ( n, &seed );
  for ( i = 0; i < n; i++ )
  {
    x[i] = r8_nint ( 50.0 * x[i] - 25.0 ) / 5.0;
  }
  a = helmert2 ( n, x );
  b = helmert2_inverse ( n, x );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  printf ( "  HELMERT2              %4d  %14g  %14g  %14g  %14g\n",
    n, norma_frobenius, normc_frobenius, error_ac, error_ab );
  free ( a );
  free ( b );
  free ( c );
  free ( x );
/*
  HERMITE matrix.
*/
  n = 5;
  a = hermite ( n );
  b = hermite_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  printf ( "  HERMITE               %4d  %14g  %14g  %14g  %14g\n",
    n, norma_frobenius, normc_frobenius, error_ac, error_ab );
  free ( a );
  free ( b );
  free ( c );
/*
  HERNDON matrix.
*/
  n = 5;
  a = herndon ( n );
  b = herndon_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  printf ( "  HERNDON               %4d  %14g  %14g  %14g  %14g\n",
    n, norma_frobenius, normc_frobenius, error_ac, error_ab );
  free ( a );
  free ( b );
  free ( c );
/*
  HILBERT matrix.
*/
  n = 5;
  a = hilbert ( n, n );
  b = hilbert_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  printf ( "  HILBERT               %4d  %14g  %14g  %14g  %14g\n",
    n, norma_frobenius, normc_frobenius, error_ac, error_ab );
  free ( a );
  free ( b );
  free ( c );
/*
  HOUSEHOLDER matrix.
*/
  n = 5;
  seed = 123456789;
  x = r8vec_uniform_01_new ( n, &seed );
  for ( i = 0; i < n; i++ )
  {
    x[i] = r8_nint ( 50.0 * x[i] - 25.0 ) / 5.0;
  }
  a = householder ( n, x );
  b = householder_inverse ( n, x );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  printf ( "  HOUSEHOLDER           %4d  %14g  %14g  %14g  %14g\n",
    n, norma_frobenius, normc_frobenius, error_ac, error_ab );
  free ( a );
  free ( b );
  free ( c );
  free ( x );
/*
  IDENTITY matrix.
*/
  n = 5;
  a = identity ( n, n );
  b = identity_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  printf ( "  IDENTITY              %4d  %14g  %14g  %14g  %14g\n",
    n, norma_frobenius, normc_frobenius, error_ac, error_ab );
  free ( a );
  free ( b );
  free ( c );
/*
  ILL3 matrix.
*/
  n = 3;
  a = ill3 ( );
  b = ill3_inverse ( );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  printf ( "  ILL3                  %4d  %14g  %14g  %14g  %14g\n",
    n, norma_frobenius, normc_frobenius, error_ac, error_ab );
  free ( a );
  free ( b );
  free ( c );
/*
  INTEGRATION matrix.
*/
  n = 6;
  seed = 123456789;
  alpha = r8_uniform_01 ( &seed );
  alpha = r8_nint ( 50.0 * alpha ) / 5.0;
  a = integration ( alpha, n );
  b = integration_inverse ( alpha, n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  printf ( "  INTEGRATION           %4d  %14g  %14g  %14g  %14g\n",
    n, norma_frobenius, normc_frobenius, error_ac, error_ab );
  free ( a );
  free ( b );
  free ( c );
/*
  INVOL matrix.
*/
  n = 5;
  a = invol ( n );
  b = invol_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  printf ( "  INVOL                 %4d  %14g  %14g  %14g  %14g\n",
    n, norma_frobenius, normc_frobenius, error_ac, error_ab );
  free ( a );
  free ( b );
  free ( c );
/*
  JORDAN matrix.
*/
  n = 6;
  seed = 123456789;
  alpha = r8_uniform_01 ( &seed );
  alpha = r8_nint ( 50.0 * alpha ) / 5.0;
  a = jordan ( alpha, n, n );
  b = jordan_inverse ( alpha, n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  printf ( "  JORDAN                %4d  %14g  %14g  %14g  %14g\n",
    n, norma_frobenius, normc_frobenius, error_ac, error_ab );
  free ( a );
  free ( b );
  free ( c );
/*
  KAHAN matrix.
*/
  n = 5;
  seed = 123456789;
  alpha = r8_uniform_01 ( &seed );
  alpha = r8_nint ( 50.0 * alpha ) / 5.0;
  a = kahan ( alpha, n, n );
  b = kahan_inverse ( alpha, n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  printf ( "  KAHAN                 %4d  %14g  %14g  %14g  %14g\n",
    n, norma_frobenius, normc_frobenius, error_ac, error_ab );
  free ( a );
  free ( b );
  free ( c );
/*
  KERSHAW matrix.
*/
  n = 4;
  a = kershaw ( );
  b = kershaw_inverse ( );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  printf ( "  KERSHAW               %4d  %14g  %14g  %14g  %14g\n",
    n, norma_frobenius, normc_frobenius, error_ac, error_ab );
  free ( a );
  free ( b );
  free ( c );
/*
  KERSHAWTRI matrix.
*/
  n = 5;
  x_n = ( n + 1 ) / 2;
  seed = 123456789;
  x = r8vec_uniform_01_new ( x_n, &seed );
  for ( i = 0; i < x_n; i++ )
  {
    x[i] = r8_nint ( 50.0 * x[i] - 25.0 ) / 5.0;
  }
  a = kershawtri ( n, x );
  b = kershawtri_inverse ( n, x );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  printf ( "  KERSHAWTRI            %4d  %14g  %14g  %14g  %14g\n",
    n, norma_frobenius, normc_frobenius, error_ac, error_ab );
  free ( a );
  free ( b );
  free ( c );
  free ( x );
/*
  KMS matrix.
*/
  n = 5;
  seed = 123456789;
  alpha = r8_uniform_01 ( &seed );
  alpha = r8_nint ( 50.0 * alpha ) / 5.0;
  a = kms ( alpha, n, n );
  b = kms_inverse ( alpha, n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  printf ( "  KMS                   %4d  %14g  %14g  %14g  %14g\n",
    n, norma_frobenius, normc_frobenius, error_ac, error_ab );
  free ( a );
  free ( b );
  free ( c );
/*
  LAGUERRE matrix.
*/
  n = 5;
  a = laguerre ( n );
  b = laguerre_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  printf ( "  LAGUERRE              %4d  %14g  %14g  %14g  %14g\n",
    n, norma_frobenius, normc_frobenius, error_ac, error_ab );
  free ( a );
  free ( b );
  free ( c );
/*
  LEGENDRE matrix.
*/
  n = 5;
  a = legendre ( n );
  b = legendre_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  printf ( "  LEGENDRE              %4d  %14g  %14g  %14g  %14g\n",
    n, norma_frobenius, normc_frobenius, error_ac, error_ab );
  free ( a );
  free ( b );
  free ( c );
/*
  LEHMER matrix.
*/
  n = 5;
  a = lehmer ( n, n );
  b = lehmer_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  printf ( "  LEHMER                %4d  %14g  %14g  %14g  %14g\n",
    n, norma_frobenius, normc_frobenius, error_ac, error_ab );
  free ( a );
  free ( b );
  free ( c );
/*
  LIETZKE matrix.
*/
  n = 5;
  a = lietzke ( n );
  b = lietzke_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  printf ( "  LIETZKE               %4d  %14g  %14g  %14g  %14g\n",
    n, norma_frobenius, normc_frobenius, error_ac, error_ab );
  free ( a );
  free ( b );
  free ( c );
/*
  LOTKIN matrix.
*/
  n = 5;
  a = lotkin ( n, n );
  b = lotkin_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  printf ( "  LOTKIN                %4d  %14g  %14g  %14g  %14g\n",
    n, norma_frobenius, normc_frobenius, error_ac, error_ab );
  free ( a );
  free ( b );
  free ( c );
/*
  MAXIJ matrix.
*/
  n = 5;
  a = maxij ( n, n );
  b = maxij_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  printf ( "  JORDAN                %4d  %14g  %14g  %14g  %14g\n",
    n, norma_frobenius, normc_frobenius, error_ac, error_ab );
  free ( a );
  free ( b );
  free ( c );
/*
  MILNES matrix.
*/
  n = 5;
  seed = 123456789;
  x = r8vec_uniform_01_new ( n, &seed );
  for ( i = 0; i < n; i++ )
  {
    x[i] = r8_nint ( 50.0 * x[i] - 25.0 ) / 5.0;
  }
  a = milnes ( n, n, x );
  b = milnes_inverse ( n, x );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  printf ( "  MILNES                %4d  %14g  %14g  %14g  %14g\n",
    n, norma_frobenius, normc_frobenius, error_ac, error_ab );
  free ( a );
  free ( b );
  free ( c );
  free ( x );
/*
  MINIJ matrix.
*/
  n = 5;
  a = minij ( n, n );
  b = minij_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  printf ( "  MINIJ                 %4d  %14g  %14g  %14g  %14g\n",
    n, norma_frobenius, normc_frobenius, error_ac, error_ab );
  free ( a );
  free ( b );
  free ( c );
/*
  MOLER1 matrix.
*/
  n = 5;
  seed = 123456789;
  alpha = r8_uniform_01 ( &seed );
  alpha = r8_nint ( 50.0 * alpha ) / 5.0;
  a = moler1 ( alpha, n, n );
  b = moler1_inverse ( alpha, n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  printf ( "  MOLER1                %4d  %14g  %14g  %14g  %14g\n",
    n, norma_frobenius, normc_frobenius, error_ac, error_ab );
  free ( a );
  free ( b );
  free ( c );
/*
  MOLER3 matrix.
*/
  n = 5;
  a = moler3 ( n, n );
  b = moler3_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  printf ( "  MOLER3                %4d  %14g  %14g  %14g  %14g\n",
    n, norma_frobenius, normc_frobenius, error_ac, error_ab );
  free ( a );
  free ( b );
  free ( c );
/*
  ORTEGA matrix.
*/
  n = 5;
  seed = 123456789;
  v1 = r8vec_uniform_01_new ( n, &seed );
  for ( i = 0; i < n; i++ )
  {
    v1[i] = r8_nint ( 50.0 * v1[i] - 25.0  ) / 5.0;
  }
  v2 = r8vec_uniform_01_new ( n, &seed );
  for ( i = 0; i < n; i++ )
  {
    v2[i] = r8_nint ( 50.0 * v2[i] - 25.0  ) / 5.0;
  }
  v3 = r8vec_uniform_01_new ( n, &seed );
  for ( i = 0; i < n; i++ )
  {
    v3[i] = r8_nint ( 50.0 * v3[i] - 25.0  ) / 5.0;
  }
  a = ortega ( n, v1, v2, v3 );
  b = ortega_inverse ( n, v1, v2, v3 );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  printf ( "  ORTEGA                %4d  %14g  %14g  %14g  %14g\n",
    n, norma_frobenius, normc_frobenius, error_ac, error_ab );
  free ( a );
  free ( b );
  free ( c );
  free ( v1 );
  free ( v2 );
  free ( v3 );
/*
  ORTH_SYMM matrix.
*/
  n = 5;
  a = orth_symm ( n );
  b = orth_symm_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  printf ( "  ORTH_SYMM             %4d  %14g  %14g  %14g  %14g\n",
    n, norma_frobenius, normc_frobenius, error_ac, error_ab );
  free ( a );
  free ( b );
  free ( c );
/*
  OTO matrix.
*/
  n = 5;
  a = oto ( n, n );
  b = oto_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  printf ( "  OTO                   %4d  %14g  %14g  %14g  %14g\n",
    n, norma_frobenius, normc_frobenius, error_ac, error_ab );
  free ( a );
  free ( b );
  free ( c );
/*
  PARTER matrix.
*/
  n = 5;
  a = parter ( n, n );
  b = parter_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  printf ( "  PARTER                %4d  %14g  %14g  %14g  %14g\n",
    n, norma_frobenius, normc_frobenius, error_ac, error_ab );
  free ( a );
  free ( b );
  free ( c );
/*
  PASCAL1 matrix.
*/
  n = 5;
  a = pascal1 ( n );
  b = pascal1_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  printf ( "  PASCAL1               %4d  %14g  %14g  %14g  %14g\n",
    n, norma_frobenius, normc_frobenius, error_ac, error_ab );
  free ( a );
  free ( b );
  free ( c );
/*
  PASCAL2 matrix.
*/
  n = 5;
  a = pascal2 ( n );
  b = pascal2_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  printf ( "  PASCAL2               %4d  %14g  %14g  %14g  %14g\n",
    n, norma_frobenius, normc_frobenius, error_ac, error_ab );
  free ( a );
  free ( b );
  free ( c );
/*
  PASCAL3 matrix.
*/
  n = 5;
  seed = 123456789;
  alpha = r8_uniform_01 ( &seed );
  alpha = r8_nint ( 50.0 * alpha - 25.0 ) / 5.0;
  a = pascal3 ( n, alpha );
  b = pascal3_inverse ( n, alpha );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  printf ( "  PASCAL3               %4d  %14g  %14g  %14g  %14g\n",
    n, norma_frobenius, normc_frobenius, error_ac, error_ab );
  free ( a );
  free ( b );
  free ( c );
/*
  PDS_RANDOM matrix.
*/
  n = 5;
  seed_save = 123456789;
  seed = seed_save;
  a = pds_random ( n, &seed );
  seed = seed_save;
  b = pds_random_inverse ( n, &seed );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  printf ( "  PDS_RANDOM            %4d  %14g  %14g  %14g  %14g\n",
    n, norma_frobenius, normc_frobenius, error_ac, error_ab );
  free ( a );
  free ( b );
  free ( c );
/*
  PEI matrix.
*/
  n = 5;
  seed = 123456789;
  alpha = r8_uniform_01 ( &seed );
  alpha = r8_nint ( 50.0 * alpha ) / 5.0;
  a = pei ( alpha, n );
  b = pei_inverse ( alpha, n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  printf ( "  PEI                   %4d  %14g  %14g  %14g  %14g\n",
    n, norma_frobenius, normc_frobenius, error_ac, error_ab );
  free ( a );
  free ( b );
  free ( c );
/*
  PERMUTATION_RANDOM matrix.
*/
  n = 5;
  seed = 123456789;
  seed_save = seed;
  a = permutation_random ( n, &seed );
  seed = seed_save;
  b = permutation_random_inverse ( n, &seed );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  printf ( "  PERMUTATION_RANDOM    %4d  %14g  %14g  %14g  %14g\n",
    n, norma_frobenius, normc_frobenius, error_ac, error_ab );
  free ( a );
  free ( b );
  free ( c );
/*
  PLU matrix.
*/
  n = 5;
  pivot = ( int * ) malloc ( n * sizeof ( int ) );
  p = ( double * ) malloc ( n * n * sizeof ( double ) );
  l = ( double * ) malloc ( n * n * sizeof ( double ) );
  u = ( double * ) malloc ( n * n * sizeof ( double ) );
  for ( i = 0; i < n; i++ )
  {
    pivot[i] = i + 1;
  }
  a = plu ( n, pivot, p, l, u );
  b = plu_inverse ( n, p, l, u );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  printf ( "  PLU                   %4d  %14g  %14g  %14g  %14g\n",
    n, norma_frobenius, normc_frobenius, error_ac, error_ab );
  free ( a );
  free ( b );
  free ( c );
  free ( l );
  free ( p );
  free ( pivot );
  free ( u );
/*
  RIS matrix.
*/
  n = 5;
  a = ris ( n );
  b = ris_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  printf ( "  RIS                   %4d  %14g  %14g  %14g  %14g\n",
    n, norma_frobenius, normc_frobenius, error_ac, error_ab );
  free ( a );
  free ( b );
  free ( c );
/*
  RODMAN matrix.
*/
  n = 5;
  seed = 123456789;
  alpha = r8_uniform_01 ( &seed );
  alpha = r8_nint ( 50.0 * alpha - 25.0 ) / 5.0;
  a = rodman ( alpha, n, n );
  b = rodman_inverse ( alpha, n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  printf ( "  RODMAN                %4d  %14g  %14g  %14g  %14g\n",
    n, norma_frobenius, normc_frobenius, error_ac, error_ab );
  free ( a );
  free ( b );
  free ( c );
/*
  RUTIS1 matrix.
*/
  n = 4;
  a = rutis1 ( );
  b = rutis1_inverse ( );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  printf ( "  RUTIS1                %4d  %14g  %14g  %14g  %14g\n",
    n, norma_frobenius, normc_frobenius, error_ac, error_ab );
  free ( a );
  free ( b );
  free ( c );
/*
  RUTIS2 matrix.
*/
  n = 4;
  a = rutis2 ( );
  b = rutis2_inverse ( );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  printf ( "  RUTIS2                %4d  %14g  %14g  %14g  %14g\n",
    n, norma_frobenius, normc_frobenius, error_ac, error_ab );
  free ( a );
  free ( b );
  free ( c );
/*
  RUTIS3 matrix.
*/
  n = 4;
  a = rutis3 ( );
  b = rutis3_inverse ( );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  printf ( "  RUTIS3                %4d  %14g  %14g  %14g  %14g\n",
    n, norma_frobenius, normc_frobenius, error_ac, error_ab );
  free ( a );
  free ( b );
  free ( c );
/*
  RUTIS4 matrix.
*/
  n = 5;
  a = rutis4 ( n );
  b = rutis4_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  printf ( "  RUTIS4                %4d  %14g  %14g  %14g  %14g\n",
    n, norma_frobenius, normc_frobenius, error_ac, error_ab );
  free ( a );
  free ( b );
  free ( c );
/*
  RUTIS5 matrix.
*/
  n = 4;
  a = rutis5 ( );
  b = rutis5_inverse ( );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  printf ( "  RUTIS5                %4d  %14g  %14g  %14g  %14g\n",
    n, norma_frobenius, normc_frobenius, error_ac, error_ab );
  free ( a );
  free ( b );
  free ( c );
/*
  SCHUR_BLOCK matrix.
*/
  n = 5;
  x_n = ( n + 1 ) / 2;
  y_n = n / 2;
  seed = 123456789;
  x = r8vec_uniform_01_new ( x_n, &seed );
  for ( i = 0; i < x_n; i++ )
  {
    x[i] = r8_nint ( 50.0 * x[i] - 25.0 ) / 5.0;
  }
  y = r8vec_uniform_01_new ( y_n, &seed );
  for ( i = 0; i < y_n; i++ )
  {
    y[i] = r8_nint ( 50.0 * y[i] - 25.0 ) / 5.0;
  }
  a = schur_block ( n, x, y );
  b = schur_block_inverse ( n, x, y );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  printf ( "  SCHUR_BLOCK           %4d  %14g  %14g  %14g  %14g\n",
    n, norma_frobenius, normc_frobenius, error_ac, error_ab );
  free ( a );
  free ( b );
  free ( c );
  free ( x );
  free ( y );
/*
  SPLINE matrix.
*/
  n = 5;
  seed = 123456789;
  x = r8vec_uniform_01_new ( n - 1, &seed );
  for ( i = 0; i < n - 1; i++ )
  {
    x[i] = r8_nint ( 50.0 * x[i] - 25.0 ) / 5.0;
  }
  a = spline ( n, x );
  b = spline_inverse ( n, x );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  printf ( "  SPLINE                %4d  %14g  %14g  %14g  %14g\n",
    n, norma_frobenius, normc_frobenius, error_ac, error_ab );
  free ( a );
  free ( b );
  free ( c );
  free ( x );
/*
  STIRLING matrix.
*/
  n = 5;
  a = stirling ( n, n );
  b = stirling_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  printf ( "  STIRLING              %4d  %14g  %14g  %14g  %14g\n",
    n, norma_frobenius, normc_frobenius, error_ac, error_ab );
  free ( a );
  free ( b );
  free ( c );
/*
  SUMMATION matrix.
*/
  n = 5;
  a = summation ( n, n );
  b = summation_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  printf ( "  SUMMATION             %4d  %14g  %14g  %14g  %14g\n",
    n, norma_frobenius, normc_frobenius, error_ac, error_ab );
  free ( a );
  free ( b );
  free ( c );
/*
  TRI_UPPER matrix.
*/
  n = 5;
  seed = 123456789;
  alpha = r8_uniform_01 ( &seed );
  alpha = r8_nint ( 50.0 * alpha - 25.0 ) / 5.0;
  a = tri_upper ( alpha, n );
  b = tri_upper_inverse ( alpha, n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  printf ( "  TRI_UPPER             %4d  %14g  %14g  %14g  %14g\n",
    n, norma_frobenius, normc_frobenius, error_ac, error_ab );
  free ( a );
  free ( b );
  free ( c );
/*
  TRIS matrix.
*/
  n = 5;
  seed = 123456789;
  alpha = r8_uniform_01 ( &seed );
  alpha = r8_nint ( 50.0 * alpha - 25.0 ) / 5.0;
  beta = r8_uniform_01 ( &seed );
  beta = r8_nint ( 50.0 * beta - 25.0 ) / 5.0;
  gamma = r8_uniform_01 ( &seed );
  gamma = r8_nint ( 50.0 * gamma - 25.0 ) / 5.0;
  a = tris ( n, n, alpha, beta, gamma );
  b = tris_inverse ( n, alpha, beta, gamma );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  printf ( "  TRIS                  %4d  %14g  %14g  %14g  %14g\n",
    n, norma_frobenius, normc_frobenius, error_ac, error_ab );
  free ( a );
  free ( b );
  free ( c );
/*
  TRIV matrix.
*/
  n = 5;
  seed = 123456789;
  x = r8vec_uniform_01_new ( n - 1, &seed );
  for ( i = 0; i < n - 1; i++ )
  {
    x[i] = r8_nint ( 50.0 * x[i] - 25.0 ) / 5.0;
  }
  y = r8vec_uniform_01_new ( n, &seed );
  for ( i = 0; i < n; i++ )
  {
    y[i] = r8_nint ( 50.0 * y[i] - 25.0 ) / 5.0;
  }
  z = r8vec_uniform_01_new ( n - 1, &seed );
  for ( i = 0; i < n - 1; i++ )
  {
    z[i] = r8_nint ( 50.0 * z[i] - 25.0 ) / 5.0;
  }
  a = triv ( n, x, y, z );
  b = triv_inverse ( n, x, y, z );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  printf ( "  TRIV                  %4d  %14g  %14g  %14g  %14g\n",
    n, norma_frobenius, normc_frobenius, error_ac, error_ab );
  free ( a );
  free ( b );
  free ( c );
  free ( x );
  free ( y );
  free ( z );
/*
  TRIW matrix.
*/
  n = 5;
  seed = 123456789;
  k = i4_uniform ( 0, n - 1, &seed );
  alpha = r8_uniform_01 ( &seed );
  alpha = r8_nint ( 50.0 * alpha - 25.0 ) / 5.0;
  a = triw ( alpha, k, n );
  b = triw_inverse ( alpha, k, n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  printf ( "  TRIW                  %4d  %14g  %14g  %14g  %14g\n",
    n, norma_frobenius, normc_frobenius, error_ac, error_ab );
  free ( a );
  free ( b );
  free ( c );
/*
  UPSHIFT matrix.
*/
  n = 5;
  a = upshift ( n );
  b = upshift_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  printf ( "  UPSHIFT               %4d  %14g  %14g  %14g  %14g\n",
    n, norma_frobenius, normc_frobenius, error_ac, error_ab );
  free ( a );
  free ( b );
  free ( c );
/*
  VAND1 matrix.
*/
  n = 5;
  seed = 123456789;
  x = r8vec_uniform_01_new ( n, &seed );
  for ( i = 0; i < n; i++ )
  {
    x[i] = r8_nint ( 50.0 * x[i] - 25.0 ) / 5.0;
  }
  a = vand1 ( n, x );
  b = vand1_inverse ( n, x );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  printf ( "  VAND1                 %4d  %14g  %14g  %14g  %14g\n",
    n, norma_frobenius, normc_frobenius, error_ac, error_ab );
  free ( a );
  free ( b );
  free ( c );
  free ( x );
/*
  VAND2 matrix.
*/
  n = 5;
  seed = 123456789;
  x = r8vec_uniform_01_new ( n, &seed );
  for ( i = 0; i < n; i++ )
  {
    x[i] = r8_nint ( 50.0 * x[i] - 25.0 ) / 5.0;
  }
  a = vand2 ( n, x );
  b = vand2_inverse ( n, x );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  printf ( "  VAND2                 %4d  %14g  %14g  %14g  %14g\n",
    n, norma_frobenius, normc_frobenius, error_ac, error_ab );
  free ( a );
  free ( b );
  free ( c );
/*
  WILK03 matrix.
*/
  n = 3;
  a = wilk03 ( );
  b = wilk03_inverse ( );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  printf ( "  WILK03                %4d  %14g  %14g  %14g  %14g\n",
    n, norma_frobenius, normc_frobenius, error_ac, error_ab );
  free ( a );
  free ( b );
  free ( c );
/*
  WILK04 matrix.
*/
  n = 4;
  a = wilk04 ( );
  b = wilk04_inverse ( );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  printf ( "  WILK04                %4d  %14g  %14g  %14g  %14g\n",
    n, norma_frobenius, normc_frobenius, error_ac, error_ab );
  free ( a );
  free ( b );
  free ( c );
/*
  WILK05 matrix.
*/
  n = 5;
  a = wilk05 ( );
  b = wilk05_inverse ( );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  printf ( "  WILK05                %4d  %14g  %14g  %14g  %14g\n",
    n, norma_frobenius, normc_frobenius, error_ac, error_ab );
  free ( a );
  free ( b );
  free ( c );
/*
  WILK21 matrix.
*/
  n = 21;
  a = wilk21 ( n );
  b = wilk21_inverse ( n );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  printf ( "  WILK21                %4d  %14g  %14g  %14g  %14g\n",
    n, norma_frobenius, normc_frobenius, error_ac, error_ab );
  free ( a );
  free ( b );
  free ( c );
/*
  WILSON matrix.
*/
  n = 4;
  a = wilson ( );
  b = wilson_inverse ( );
  c = r8mat_inverse ( n, a );
  error_ab = r8mat_is_inverse ( n, a, b );
  error_ac = r8mat_is_inverse ( n, a, c );
  norma_frobenius = r8mat_norm_fro ( n, n, a );
  normc_frobenius = r8mat_norm_fro ( n, n, c );
  printf ( "  WILSON                %4d  %14g  %14g  %14g  %14g\n",
    n, norma_frobenius, normc_frobenius, error_ac, error_ab );
  free ( a );
  free ( b );
  free ( c );

  return;
}
/******************************************************************************/

void test_null ( void )

/******************************************************************************/
/*
  Purpose:

    TEST_NULL tests the null vectors.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 June 2011

  Author:

    John Burkardt
*/
{
  double *a;
  double *at;
  double alpha;
  int col_num;
  double error_l2;
  double f1;
  double f2;
  int m;
  int mt;
  int n;
  int nt;
  double norm_a_frobenius;
  double norm_x_l2;
  int row_num;
  int seed;
  char title[21];
  double *x;

  printf ( "\n" );
  printf ( "TEST_NULL\n" );
  printf ( "  A = a test matrix of order M by N\n" );
  printf ( "  x = an N vector, candidate for a null vector.\n" );
  printf ( "\n" );
  printf ( "  ||A|| = Frobenius norm of A.\n" );
  printf ( "  ||x|| = L2 norm of x.\n" );
  printf ( "  ||A*x||/||x|| = L2 norm of A*x over L2 norm of x.\n" );
  printf ( "\n" );
  printf ( "  Matrix title	           M     N      " );
  printf ( "||A||            ||x||        ||A*x||/||x||\n" );
  printf ( "\n" );
/*
  ARCHIMEDES matrix.
*/
  m = 7;
  n = 8;
  a = archimedes ( );
  x = archimedes_null ( );
  error_l2 = r8mat_is_null_vector ( m, n, a, x );
  norm_a_frobenius = r8mat_norm_fro ( m, n, a );
  norm_x_l2 = r8vec_norm_l2 ( n, x );
  printf ( "  ARCHIMEDES            %4d  %4d  %14g  %14g  %14g\n",
    m, n, norm_a_frobenius, norm_x_l2, error_l2 );
  free ( a );
  free ( x );
/*
  CHEBY_DIFF1 matrix.
*/
  m = 5;
  n = 5;
  a = cheby_diff1 ( n );
  x = cheby_diff1_null ( n );
  error_l2 = r8mat_is_null_vector ( m, n, a, x );
  norm_a_frobenius = r8mat_norm_fro ( m, n, a );
  norm_x_l2 = r8vec_norm_l2 ( n, x );
  printf ( "  CHEBY_DIFF1           %4d  %4d  %14g  %14g  %14g\n",
    m, n, norm_a_frobenius, norm_x_l2, error_l2 );
  free ( a );
  free ( x );
/*
  CREATION matrix.
*/
  m = 5;
  n = 5;
  a = creation ( m, n );
  x = creation_null ( n );
  error_l2 = r8mat_is_null_vector ( m, n, a, x );
  norm_a_frobenius = r8mat_norm_fro ( m, n, a );
  norm_x_l2 = r8vec_norm_l2 ( n, x );
  printf ( "  CREATION              %4d  %4d  %14g  %14g  %14g\n",
    m, n, norm_a_frobenius, norm_x_l2, error_l2 );
  free ( a );
  free ( x );
/*
  DIF1 matrix.
  Only has null vectors for N odd.
*/
  if ( 1 )
  {
    m = 5;
    n = 5;
    a = dif1 ( m, n );
    x = dif1_null ( n );
    error_l2 = r8mat_is_null_vector ( m, n, a, x );
    norm_a_frobenius = r8mat_norm_fro ( m, n, a );
    norm_x_l2 = r8vec_norm_l2 ( n, x );
    printf ( "  DIF1                  %4d  %4d  %14g  %14g  %14g\n",
      m, n, norm_a_frobenius, norm_x_l2, error_l2 );
    free ( a );
    free ( x );
  }
/*
  DIF1CYCLIC matrix.
*/
  m = 5;
  n = 5;
  a = dif1cyclic ( n );
  x = dif1cyclic_null ( n );
  error_l2 = r8mat_is_null_vector ( m, n, a, x );
  norm_a_frobenius = r8mat_norm_fro ( m, n, a );
  norm_x_l2 = r8vec_norm_l2 ( n, x );
  printf ( "  DIF1CYCLIC            %4d  %4d  %14g  %14g  %14g\n",
    m, n, norm_a_frobenius, norm_x_l2, error_l2 );
  free ( a );
  free ( x );
/*
  DIF2CYCLIC matrix.
*/
  m = 5;
  n = 5;
  a = dif2cyclic ( n );
  x = dif2cyclic_null ( n );
  error_l2 = r8mat_is_null_vector ( m, n, a, x );
  norm_a_frobenius = r8mat_norm_fro ( m, n, a );
  norm_x_l2 = r8vec_norm_l2 ( n, x );
  printf ( "  DIF2CYCLIC            %4d  %4d  %14g  %14g  %14g\n",
    m, n, norm_a_frobenius, norm_x_l2, error_l2 );
  free ( a );
  free ( x );
/*
  EBERLEIN matrix.
  We have a LEFT null vector.
*/
  m = 5;
  n = 5;
  seed = 123456789;
  alpha = r8_uniform_01 ( &seed );
  alpha = r8_nint ( 50.0 * alpha - 25.0 ) / 5.0;
  a = eberlein ( alpha, n );
  mt = n;
  nt = m;
  at = r8mat_transpose_new ( m, n, a );
  x = eberlein_null_left ( n );
  error_l2 = r8mat_is_null_vector ( mt, nt, at, x );
  norm_a_frobenius = r8mat_norm_fro ( mt, nt, at );
  norm_x_l2 = r8vec_norm_l2 ( nt, x );
   printf ( "  EBERLEIN (left)       %4d  %4d  %14g  %14g  %14g\n",
    m, n, norm_a_frobenius, norm_x_l2, error_l2 );
  free ( a );
  free ( at );
  free ( x );
/*
  FIBONACCI1 matrix.
*/
  m = 5;
  n = 5;
  seed = 123456789;
  f1 = r8_uniform_01 ( &seed );
  f1 = r8_nint ( 50.0 * f1 - 25.0 ) / 5.0;
  f2 = r8_uniform_01 ( &seed );
  f2 = r8_nint ( 50.0 * f2 - 25.0 ) / 5.0;
  a = fibonacci1 ( n, f1, f2 );
  x = fibonacci1_null ( n, f1, f2 );
  error_l2 = r8mat_is_null_vector ( m, n, a, x );
  norm_a_frobenius = r8mat_norm_fro ( m, n, a );
  norm_x_l2 = r8vec_norm_l2 ( n, x );
  printf ( "  FIBONACCI1            %4d  %4d  %14g  %14g  %14g\n",
    m, n, norm_a_frobenius, norm_x_l2, error_l2 );
  free ( a );
  free ( x );
/*
  LAUCHLI matrix.
  We have a LEFT null vector of a RECTANGULAR matrix.
*/
  m = 6;
  n = m - 1;
  seed = 123456789;
  alpha = r8_uniform_01 ( &seed );
  alpha = r8_nint ( 50.0 * alpha - 25.0 ) / 5.0;
  a = lauchli ( alpha, m, n );
  mt = n;
  nt = m;
  at = r8mat_transpose_new ( m, n, a );
  x = lauchli_null_left ( alpha, m, n );
  error_l2 = r8mat_is_null_vector ( mt, nt, at, x );
  norm_a_frobenius = r8mat_norm_fro ( mt, nt, at );
  norm_x_l2 = r8vec_norm_l2 ( nt, x );
  printf ( "  LAUCHLI (left)        %4d  %4d  %14g  %14g  %14g\n",
    m, n, norm_a_frobenius, norm_x_l2, error_l2 );
  free ( a );
  free ( at );
  free ( x );
/*
  LINE_ADJ matrix.
*/
  m = 7;
  n = 7;
  a = line_adj ( n );
  x = line_adj_null ( n );
  error_l2 = r8mat_is_null_vector ( m, n, a, x );
  norm_a_frobenius = r8mat_norm_fro ( m, n, a );
  norm_x_l2 = r8vec_norm_l2 ( n, x );
  printf ( "  LINE_ADJ              %4d  %4d  %14g  %14g  %14g\n",
    m, n, norm_a_frobenius, norm_x_l2, error_l2 );
  free ( a );
  free ( x );
/*
  MOLER2 matrix.
*/
  m = 5;
  n = 5;
  a = moler2 ( );
  x = moler2_null ( );
  error_l2 = r8mat_is_null_vector ( m, n, a, x );
  norm_a_frobenius = r8mat_norm_fro ( m, n, a );
  norm_x_l2 = r8vec_norm_l2 ( n, x );
  printf ( "  MOLER2                %4d  %4d  %14g  %14g  %14g\n",
    m, n, norm_a_frobenius, norm_x_l2, error_l2 );
  free ( a );
  free ( x );
/*
  NEUMANN matrix.
*/
  row_num = 5;
  col_num = 5;
  m = row_num * col_num;
  n = row_num * col_num;
  a = neumann ( row_num, col_num );
  x = neumann_null ( row_num, col_num );
  error_l2 = r8mat_is_null_vector ( m, n, a, x );
  norm_a_frobenius = r8mat_norm_fro ( m, n, a );
  norm_x_l2 = r8vec_norm_l2 ( n, x );
  printf ( "  NEUMANN               %4d  %4d  %14g  %14g  %14g\n",
    m, n, norm_a_frobenius, norm_x_l2, error_l2 );
  free ( a );
  free ( x );
/*
  ONE matrix.
*/
  m = 5;
  n = 5;
  a = one ( n, n );
  x = one_null ( n );
  error_l2 = r8mat_is_null_vector ( m, n, a, x );
  norm_a_frobenius = r8mat_norm_fro ( m, n, a );
  norm_x_l2 = r8vec_norm_l2 ( n, x );
  printf ( "  ONE                   %4d  %4d  %14g  %14g  %14g\n",
    m, n, norm_a_frobenius, norm_x_l2, error_l2 );
  free ( a );
  free ( x );
/*
  RING_ADJ matrix.
  N must be a multiple of 4 for there to be a null vector.
*/
  m = 12;
  n = 12;
  a = ring_adj ( n );
  x = ring_adj_null ( n );
  error_l2 = r8mat_is_null_vector ( m, n, a, x );
  norm_a_frobenius = r8mat_norm_fro ( m, n, a );
  norm_x_l2 = r8vec_norm_l2 ( n, x );
  printf ( "  RING_ADJ              %4d  %4d  %14g  %14g  %14g\n",
    m, n, norm_a_frobenius, norm_x_l2, error_l2 );
  free ( a );
  free ( x );
/*
  ROSSER1 matrix.
*/
  m = 8;
  n = 8;
  a = rosser1 ( );
  x = rosser1_null ( );
  error_l2 = r8mat_is_null_vector ( m, n, a, x );
  norm_a_frobenius = r8mat_norm_fro ( m, n, a );
  norm_x_l2 = r8vec_norm_l2 ( n, x );
  printf ( "  ROSSER1               %4d  %4d  %14g  %14g  %14g\n",
    m, n, norm_a_frobenius, norm_x_l2, error_l2 );
  free ( a );
  free ( x );
/*
  ZERO matrix.
*/
  m = 5;
  n = 5;
  a = zero ( m, n );
  x = zero_null ( n );
  error_l2 = r8mat_is_null_vector ( m, n, a, x );
  norm_a_frobenius = r8mat_norm_fro ( m, n, a );
  norm_x_l2 = r8vec_norm_l2 ( n, x );
  printf ( "  ZERO                  %4d  %4d  %14g  %14g  %14g\n",
    m, n, norm_a_frobenius, norm_x_l2, error_l2 );
  free ( a );
  free ( x );

  return;
}
/******************************************************************************/

void test_plu ( void )

/******************************************************************************/
/*
  Purpose:

    TEST_PLU tests the PLU factors.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 June 2011

  Author:

    John Burkardt
*/
{
  double *a;
  double alpha;
  double error_frobenius;
  double *l;
  int m;
  int n;
  double norm_a_frobenius;
  double *p;
  int seed;
  char title[21];
  double *u;

  printf ( "\n" );
  printf ( "TEST_PLU\n" );
  printf ( "  A = a test matrix of order M by N\n" );
  printf ( "  P, L, U are the PLU factors.\n" );
  printf ( "\n" );
  printf ( "  ||A|| = Frobenius norm of A.\n" );
  printf ( "  ||A-PLU|| = Frobenius norm of A-P*L*U.\n" );
  printf ( "\n" );
  printf ( "  Matrix title	           M     N      " );
  printf ( "||A||            ||A-PLU||\n" );
  printf ( "\n" );
/*
  BODEWIG matrix.
*/
  m = 4;
  n = 4;
  p = ( double * ) malloc ( m * m * sizeof ( double ) );
  l = ( double * ) malloc ( m * m * sizeof ( double ) );
  u = ( double * ) malloc ( m * n * sizeof ( double ) );
  a = bodewig ( );
  bodewig_plu ( p, l, u );
  error_frobenius = r8mat_is_plu ( m, n, a, p, l, u );
  norm_a_frobenius = r8mat_norm_fro ( m, n, a );
  printf ( "  BODEWIG               %4d  %4d  %14g  %14g\n",
    m, n, norm_a_frobenius, error_frobenius );
  free ( a );
  free ( l );
  free ( p );
  free ( u );
/*
  BORDERBAND matrix.
*/
  m = 5;
  n = 5;
  p = ( double * ) malloc ( m * m * sizeof ( double ) );
  l = ( double * ) malloc ( m * m * sizeof ( double ) );
  u = ( double * ) malloc ( m * n * sizeof ( double ) );
  a = borderband ( n );
  borderband_plu ( n, p, l, u );
  error_frobenius = r8mat_is_plu ( m, n, a, p, l, u );
  norm_a_frobenius = r8mat_norm_fro ( m, n, a );
  printf ( "  BORDERBAND            %4d  %4d  %14g  %14g\n",
    m, n, norm_a_frobenius, error_frobenius );
  free ( a );
  free ( l );
  free ( p );
  free ( u );
/*
  DIF2 matrix.
*/
  m = 5;
  n = 5;
  p = ( double * ) malloc ( m * m * sizeof ( double ) );
  l = ( double * ) malloc ( m * m * sizeof ( double ) );
  u = ( double * ) malloc ( m * n * sizeof ( double ) );
  a = dif2 ( m, n );
  dif2_plu ( n, p, l, u );
  error_frobenius = r8mat_is_plu ( m, n, a, p, l, u );
  norm_a_frobenius = r8mat_norm_fro ( m, n, a );
  printf ( "  DIF2                  %4d  %4d  %14g  %14g\n",
    m, n, norm_a_frobenius, error_frobenius );
  free ( a );
  free ( l );
  free ( p );
  free ( u );
/*
  GFPP matrix.
*/
  m = 5;
  n = 5;
  seed = 123456789;
  alpha = r8_uniform_01 ( &seed );
  alpha = r8_nint ( 50.0 * alpha ) / 5.0;
  p = ( double * ) malloc ( m * m * sizeof ( double ) );
  l = ( double * ) malloc ( m * m * sizeof ( double ) );
  u = ( double * ) malloc ( m * n * sizeof ( double ) );
  a = gfpp ( n, alpha );
  gfpp_plu ( n, alpha, p, l, u );
  error_frobenius = r8mat_is_plu ( m, n, a, p, l, u );
  norm_a_frobenius = r8mat_norm_fro ( m, n, a );
  printf ( "  GFPP                  %4d  %4d  %14g  %14g\n",
    m, n, norm_a_frobenius, error_frobenius );
  free ( a );
  free ( l );
  free ( p );
  free ( u );
/*
  GIVENS matrix.
*/
  m = 5;
  n = 5;
  p = ( double * ) malloc ( m * m * sizeof ( double ) );
  l = ( double * ) malloc ( m * m * sizeof ( double ) );
  u = ( double * ) malloc ( m * n * sizeof ( double ) );
  a = givens ( n, n );
  givens_plu ( n, p, l, u );
  error_frobenius = r8mat_is_plu ( m, n, a, p, l, u );
  norm_a_frobenius = r8mat_norm_fro ( m, n, a );
  printf ( "  GIVENS                %4d  %4d  %14g  %14g\n",
    m, n, norm_a_frobenius, error_frobenius );
  free ( a );
  free ( l );
  free ( p );
  free ( u );
/*
  KMS matrix.
*/
  m = 5;
  n = 5;
  seed = 123456789;
  alpha = r8_uniform_01 ( &seed );
  alpha = r8_nint ( 50.0 * alpha ) / 5.0;
  p = ( double * ) malloc ( m * m * sizeof ( double ) );
  l = ( double * ) malloc ( m * m * sizeof ( double ) );
  u = ( double * ) malloc ( m * n * sizeof ( double ) );
  a = kms ( alpha, m, n );
  kms_plu ( alpha, n, p, l, u );
  error_frobenius = r8mat_is_plu ( m, n, a, p, l, u );
  norm_a_frobenius = r8mat_norm_fro ( m, n, a );
  printf ( "  KMS                   %4d  %4d  %14g  %14g\n",
    m, n, norm_a_frobenius, error_frobenius );
  free ( a );
  free ( l );
  free ( p );
  free ( u );
/*
  MAXIJ matrix.
*/
  m = 5;
  n = 5;
  p = ( double * ) malloc ( m * m * sizeof ( double ) );
  l = ( double * ) malloc ( m * m * sizeof ( double ) );
  u = ( double * ) malloc ( m * n * sizeof ( double ) );
  a = maxij ( n, n );
  maxij_plu ( n, p, l, u );
  error_frobenius = r8mat_is_plu ( m, n, a, p, l, u );
  norm_a_frobenius = r8mat_norm_fro ( m, n, a );
  printf ( "  MAXIJ                 %4d  %4d  %14g  %14g\n",
    m, n, norm_a_frobenius, error_frobenius );
  free ( a );
  free ( l );
  free ( p );
  free ( u );
/*
  MINIJ matrix.
*/
  m = 5;
  n = 5;
  p = ( double * ) malloc ( m * m * sizeof ( double ) );
  l = ( double * ) malloc ( m * m * sizeof ( double ) );
  u = ( double * ) malloc ( m * n * sizeof ( double ) );
  a = minij ( m, n );
  minij_plu ( n, p, l, u );
  error_frobenius = r8mat_is_plu ( m, n, a, p, l, u );
  norm_a_frobenius = r8mat_norm_fro ( m, n, a );
  printf ( "  MINIJ                 %4d  %4d  %14g  %14g\n",
    m, n, norm_a_frobenius, error_frobenius );
  free ( a );
  free ( l );
  free ( p );
  free ( u );
/*
  MOLER1 matrix.
*/
  m = 5;
  n = 5;
  seed = 123456789;
  alpha = r8_uniform_01 ( &seed );
  alpha = r8_nint ( 50.0 * alpha ) / 5.0;
  p = ( double * ) malloc ( m * m * sizeof ( double ) );
  l = ( double * ) malloc ( m * m * sizeof ( double ) );
  u = ( double * ) malloc ( m * n * sizeof ( double ) );
  a = moler1 ( alpha, m, n );
  moler1_plu ( alpha, n, p, l, u );
  error_frobenius = r8mat_is_plu ( m, n, a, p, l, u );
  norm_a_frobenius = r8mat_norm_fro ( m, n, a );
  printf ( "  MOLER1                %4d  %4d  %14g  %14g\n",
    m, n, norm_a_frobenius, error_frobenius );
  free ( a );
  free ( l );
  free ( p );
  free ( u );
/*
  MOLER3 matrix.
*/
  m = 5;
  n = 5;
  p = ( double * ) malloc ( m * m * sizeof ( double ) );
  l = ( double * ) malloc ( m * m * sizeof ( double ) );
  u = ( double * ) malloc ( m * n * sizeof ( double ) );
  a = moler3 ( m, n );
  moler3_plu ( n, p, l, u );
  error_frobenius = r8mat_is_plu ( m, n, a, p, l, u );
  norm_a_frobenius = r8mat_norm_fro ( m, n, a );
  printf ( "  MOLER3                %4d  %4d  %14g  %14g\n",
    m, n, norm_a_frobenius, error_frobenius );
  free ( a );
  free ( l );
  free ( p );
  free ( u );
/*
  OTO matrix.
*/
  m = 5;
  n = 5;
  p = ( double * ) malloc ( m * m * sizeof ( double ) );
  l = ( double * ) malloc ( m * m * sizeof ( double ) );
  u = ( double * ) malloc ( m * n * sizeof ( double ) );
  a = oto ( m, n );
  oto_plu ( n, p, l, u );
  error_frobenius = r8mat_is_plu ( m, n, a, p, l, u );
  norm_a_frobenius = r8mat_norm_fro ( m, n, a );
  printf ( "  OTO                   %4d  %4d  %14g  %14g\n",
    m, n, norm_a_frobenius, error_frobenius );
  free ( a );
  free ( l );
  free ( p );
  free ( u );
/*
  PASCAL2 matrix.
*/
  m = 5;
  n = 5;
  p = ( double * ) malloc ( m * m * sizeof ( double ) );
  l = ( double * ) malloc ( m * m * sizeof ( double ) );
  u = ( double * ) malloc ( m * n * sizeof ( double ) );
  a = pascal2 ( n );
  pascal2_plu ( n, p, l, u );
  error_frobenius = r8mat_is_plu ( m, n, a, p, l, u );
  norm_a_frobenius = r8mat_norm_fro ( m, n, a );
  printf ( "  PASCAL2               %4d  %4d  %14g  %14g\n",
    m, n, norm_a_frobenius, error_frobenius );
  free ( a );
  free ( l );
  free ( p );
  free ( u );
/*
  WILSON matrix.
*/
  m = 4;
  n = 4;
  p = ( double * ) malloc ( m * m * sizeof ( double ) );
  l = ( double * ) malloc ( m * m * sizeof ( double ) );
  u = ( double * ) malloc ( m * n * sizeof ( double ) );
  a = wilson ( );
  wilson_plu ( p, l, u );
  error_frobenius = r8mat_is_plu ( m, n, a, p, l, u );
  norm_a_frobenius = r8mat_norm_fro ( m, n, a );
  printf ( "  WILSON                %4d  %4d  %14g  %14g\n",
    m, n, norm_a_frobenius, error_frobenius );
  free ( a );
  free ( l );
  free ( p );
  free ( u );

  return;
}
/******************************************************************************/

void test_solution ( )

/******************************************************************************/
/*
  Purpose:

    TEST_SOLUTION tests the linear solution computations.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 June 2011

  Author:

    John Burkardt
*/
{
  double *a;
  double alpha;
  double *b;
  double beta;
  double error_frobenius;
  double gamma;
  int i1;
  int k;
  int m;
  int n;
  int ncol;
  double norm_frobenius;
  int nrow;
  int seed;
  int seed_save;
  char title[21];
  double *x;

  printf ( "\n" );
  printf ( "TEST_SOLUTION\n" );
  printf ( "  Compute the Frobenius norm of the solution error:\n" );
  printf ( "    A * X - B\n" );
  printf ( "  given MxN matrix A, NxK solution X, MxK right hand side B.\n" );
  printf ( "\n" );
  printf ( "  Matrix title             M     N     K      ||A||         ||A*X-B||\n" );
  printf ( "\n" );
/*
  BODEWIG matrix.
*/
  m = 4;
  n = 4;
  k = 1;
  a = bodewig ( );
  b = bodewig_rhs ( );
  x = bodewig_solution ( );
  error_frobenius = r8mat_is_solution ( m, n, k, a, x, b );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  BODEWIG               %4d  %4d  %4d  %14g  %14g\n",
    m, n, k, norm_frobenius, error_frobenius );
  free ( a );
  free ( b );
  free ( x );
/*
  DIF2 matrix.
*/
  m = 10;
  n = 10;
  k = 2;
  a = dif2 ( m, n );
  b = dif2_rhs ( m, k );
  x = dif2_solution ( n, k );
  error_frobenius = r8mat_is_solution ( m, n, k, a, x, b );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  DIF2                  %4d  %4d  %4d  %14g  %14g\n",
    m, n, k, norm_frobenius, error_frobenius );
  free ( a );
  free ( b );
  free ( x );
/*
  FRANK matrix.
*/
  m = 10;
  n = 10;
  k = 2;
  a = frank ( n );
  b = frank_rhs ( m, k );
  x = frank_solution ( n, k );
  error_frobenius = r8mat_is_solution ( m, n, k, a, x, b );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  FRANK                 %4d  %4d  %4d  %14g  %14g\n",
    m, n, k, norm_frobenius, error_frobenius );
  free ( a );
  free ( b );
  free ( x );
/*
  POISSON matrix.
*/
  nrow = 4;
  ncol = 5;
  m = nrow * ncol;
  n = nrow * ncol;
  k = 1;
  a = poisson ( nrow, ncol, n );
  b = poisson_rhs ( nrow, ncol, n );
  x = poisson_solution ( nrow, ncol, n );
  error_frobenius = r8mat_is_solution ( m, n, k, a, x, b );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  POISSON               %4d  %4d  %4d  %14g  %14g\n",
    m, n, k, norm_frobenius, error_frobenius );
  free ( a );
  free ( b );
  free ( x );
/*
  WILK03 matrix.
*/
  m = 3;
  n = 3;
  k = 1;
  a = wilk03 ( );
  b = wilk03_rhs ( );
  x = wilk03_solution ( );
  error_frobenius = r8mat_is_solution ( m, n, k, a, x, b );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  WILK03                %4d  %4d  %4d  %14g  %14g\n",
    m, n, k, norm_frobenius, error_frobenius );
  free ( a );
  free ( b );
  free ( x );
/*
  WILK04 matrix.
*/
  m = 4;
  n = 4;
  k = 1;
  a = wilk04 ( );
  b = wilk04_rhs ( );
  x = wilk04_solution ( );
  error_frobenius = r8mat_is_solution ( m, n, k, a, x, b );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  WILK04                %4d  %4d  %4d  %14g  %14g\n",
    m, n, k, norm_frobenius, error_frobenius );
  free ( a );
  free ( b );
  free ( x );
/*
  WILSON matrix.
*/
  m = 4;
  n = 4;
  k = 1;
  a = wilson ( );
  b = wilson_rhs ( );
  x = wilson_solution ( );
  error_frobenius = r8mat_is_solution ( m, n, k, a, x, b );
  norm_frobenius = r8mat_norm_fro ( n, n, a );
  printf ( "  WILSON                %4d  %4d  %4d  %14g  %14g\n",
    m, n, k, norm_frobenius, error_frobenius );
  free ( a );
  free ( b );
  free ( x );

  return;
}
