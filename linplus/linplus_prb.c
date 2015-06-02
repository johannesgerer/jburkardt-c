# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <complex.h>

# include "linplus.h"

int main ( void );
void test02 ( void );
void test03 ( void );
void test035 ( void );
void test037 ( void );
void test04 ( void );
void test05 ( void );
void test06 ( void );
void test07 ( void );
void test08 ( void );
void test09 ( void );
void test10 ( void );
void test105 ( void );
void test11 ( void );
void test115 ( void );
void test12 ( void );
void test13 ( void );
void test14 ( void );
void test15 ( void );
void test151 ( void );
void test152 ( void );
void test153 ( void );
void test154 ( void );
void test155 ( void );
void test1565 ( void );
void test1566 ( void );
void test157 ( void );
void test16 ( void );
void test17 ( void );
void test18 ( void );
void test19 ( void );
void test193 ( void );
void test195 ( void );
void test197 ( void );
void test20 ( void );
void test21 ( void );
void test22 ( void );
void test23 ( void );
void test235 ( void );
void test24 ( void );
void test25 ( void );
void test26 ( void );
void test265 ( void );
void test2655 ( void );
void test27 ( void );
void test275 ( void );
void test28 ( void );
void test285 ( void );
void test29 ( void );
void test295 ( void );
void test30 ( void );
void test31 ( void );
void test315 ( void );
void test317 ( void );
void test32 ( void );
void test33 ( void );
void test34 ( void );
void test345 ( void );
void test35 ( void );
void test36 ( void );
void test37 ( void );
void test38 ( void );
void test385 ( void );
void test39 ( void );
void test40 ( void );
void test41 ( void );
void test42 ( void );
void test422 ( void );
void test423 ( void );
void test425 ( void );
void test426 ( void );
void test428 ( void );
void test43 ( void );
void test44 ( void );
void test443 ( void );
void test445 ( void );
void test45 ( void );
void test46 ( void );
void test47 ( void );
void test48 ( void );
void test485 ( void );
void test49 ( void );
void test50 ( void );
void test505 ( void );
void test51 ( void );
void test515 ( void );
void test517 ( void );
void test52 ( void );
void test525 ( void );
void test527 ( void );
void test53 ( void );
void test534 ( void );
void test535 ( void );
void test54 ( void );
void test55 ( void );
void test555 ( void );
void test56 ( void );
void test57 ( void );
void test5705 ( void );
void test571 ( void );
void test572 ( void );
void test5722 ( void );
void test5724 ( void );
void test5725 ( void );
void test573 ( void );
void test574 ( void );
void test5745 ( void );
void test575 ( void );
void test577 ( void );
void test58 ( void );
void test581 ( void );
void test583 ( void );
void test585 ( void );
void test587 ( void );
void test589 ( void );
void test59 ( void );
void test60 ( void );
void test605 ( void );
void test61 ( void );
void test62 ( void );
void test63 ( void );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for LINPLUS_PRB.

  Discussion:

    LINPLUS_PRB tests the LINPLUS library.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    20 April 2013

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "LINPLUS_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the LINPLUS library.\n" );

  test02 ( );
  test03 ( );
  test035 ( );
  test037 ( );
  test04 ( );
  test05 ( );
  test06 ( );
  test07 ( );
  test08 ( );
  test09 ( );

  test10 ( );
  test105 ( );
  test11 ( );
  test115 ( );
  test12 ( );
  test13 ( );
  test14 ( );
  test15 ( );
  test151 ( );
  test152 ( );
  test153 ( );
  test154 ( );
  test155 ( );
  test1565 ( );
  test1566 ( );
  test157 ( );
  test16 ( );
  test17 ( );
  test18 ( );
  test19 ( );
  test193 ( );
  test195 ( );
  test197 ( );

  test20 ( );
  test21 ( );
  test22 ( );
  test23 ( );
  test235 ( );
  test24 ( );
  test25 ( );
  test26 ( );
  test265 ( );
  test2655 ( );
  test27 ( );
  test275 ( );
  test28 ( );
  test285 ( );
  test29 ( );
  test295 ( );

  test30 ( );
  test31 ( );
  test315 ( );
  test317 ( );
  test32 ( );
  test33 ( );
  test34 ( );
  test345 ( );
  test35 ( );
  test36 ( );
  test37 ( );
  test38 ( );
  test385 ( );
  test39 ( );

  test40 ( );
  test41 ( );
  test42 ( );
  test422 ( );
  test423 ( );
  test425 ( );
  test426 ( );
  test428 ( );
  test43 ( );
  test44 ( );
  test443 ( );
  test445 ( );
  test45 ( );
  test46 ( );
  test47 ( );
  test48 ( );
  test485 ( );
  test49 ( );

  test50 ( );
  test505 ( );
  test51 ( );
  test515 ( );
  test517 ( );
  test52 ( );
  test525 ( );
  test527 ( );
  test53 ( );
  test534 ( );
  test535 ( );
  test54 ( );
  test55 ( );
  test555 ( );
  test56 ( );
  test57 ( );
  test5705 ( );
  test571 ( );
  test572 ( );
  test5722 ( );
  test5724 ( );
  test5725 ( );
  test573 ( );
  test574 ( );
  test5745 ( );
  test575 ( );
  test577 ( );
  test58 ( );
  test581 ( );
  test583 ( );
  test585 ( );
  test587 ( );
  test589 ( );
  test59 ( );

  test60 ( );
  test605 ( );
  test61 ( );
  test62 ( );
  test63 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "LINPLUS_PRB\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void test02 ( )

/******************************************************************************/
/*
  Purpose:

    TEST02 tests R83_CR_FA, R83_CR_SL.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    31 May 2009

  Author:

    John Burkardt
*/
{
# define N 5

  double a[3*N];
  double *a_cr;
  double b[N];
  int debug = 0;
  int i;
  int j;
  int n = N;
  double *x;

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  R83_CR_FA factors a real tridiagonal matrix;\n" );
  printf ( "  R83_CR_SL solves a factored system.\n" );
  printf ( "\n" );
  printf ( "  Matrix order N = %d\n", n );
  printf ( "  Demonstrate multiple system solution method.\n" );
/*
  Set the matrix values.
*/
  a[0+0*3] = 0.0;
  for ( j = 1; j < n; j++ )
  {
    a[0+j*3] = -1.0;
  }
  for ( j = 0; j < n; j++ )
  {
    a[1+j*3] = 2.0;
  }
  for ( j = 0; j < n - 1; j++ )
  {
    a[2+j*3] = -1.0;
  }
  a[2+(n-1)*3] = 0.0;

  if ( debug )
  {
    r83_print ( n, a, "  Input matrix:" );
  }
/*
  Factor the matrix once.
*/
  a_cr = r83_cr_fa ( n, a );

  if ( debug )
  {
    r83_print ( 2 * n + 1, a_cr, "  Cyclic reduction factor information:" );
  }
/*
  Set up the linear systems.
*/
  for ( j = 1; j <= 2; j++ )
  {
    printf ( "\n" );
    printf ( "  Solve linear system number #%d.\n", j );

    if ( j == 1 )
    {
      for ( i = 0; i < n - 1; i++ )
      {
        b[i] = 0.0;
      }
      b[n-1] = ( double ) ( n + 1 );
    }
    else
    {
      b[0] = 1.0;
      for ( i = 1; i < n - 1; i++ )
      {
        b[i] = 0.0;
      }
      b[n-1] = 1.0;
    }
/*
  Solve the linear system.
*/
    x = r83_cr_sl ( n, a_cr, b );

    r8vec_print_some ( n, x, 1, 10, "  Solution:" );

    free ( x );
  }

  free ( a_cr );

  return;
# undef N
}
/******************************************************************************/

void test03 ( )

/******************************************************************************/
/*
  Purpose:

    TEST03 tests R83_CR_FA, R83_CR_SL.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    31 May 2009

  Author:

    John Burkardt
*/
{
# define N 10

  double a[3*N];
  double *a_cr;
  double *b;
  int i;
  int j;
  int n = N;
  double *x;

  printf ( "\n" );
  printf ( "TEST03:\n" );
  printf ( "  For a real tridiagonal matrix,\n" );
  printf ( "  using CYCLIC REDUCTION,\n" );
  printf ( "  R83_CR_FA factors;\n" );
  printf ( "  R83_CR_SL solves.\n" );
  printf ( "\n" );
  printf ( "  Matrix order N = %d\n", n );
  printf ( "  The matrix is NOT symmetric.\n" );
/*
  Set the matrix values.
*/
  a[0+0*3] = 0.0;
  for ( j = 1; j < n; j++ )
  {
    a[0+j*3] = ( double ) ( j + 1 );
  }

  for ( j = 0; j < n; j++ )
  {
    a[1+j*3] = 4.0 * ( double ) ( j + 1 );
  }

  for ( j = 0; j < n - 1; j++ )
  {
    a[2+j*3] = ( double ) ( j + 1 );
  }
  a[2+(n-1)*3] = 0.0;

  r83_print ( n, a, "  The matrix:" );
/*
  Set the desired solution.
*/
  x = r8vec_indicator_new ( n );
/*
  Compute the corresponding right hand side.
*/
  b = r83_mxv ( n, a, x );
/*
  Factor the matrix.
*/
  a_cr = r83_cr_fa ( n, a );
/*
  Solve the linear system.
*/
  x = r83_cr_sl ( n, a_cr, b );

  r8vec_print_some ( n, x, 1, 10, "  Solution:" );

  free ( a_cr );
  free ( b );
  free ( x );

  return;
# undef N
}
/******************************************************************************/

void test035 ( )

/******************************************************************************/
/*
  Purpose:

    TEST035 tests R83_GS_SL.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    25 February 2013

  Author:

    John Burkardt
*/
{
  double *a;
  double *b;
  int i;
  int j;
  int job;
  int maxit = 1000;
  int n = 100;
  double *x;

  printf ( "\n" );
  printf ( "TEST035\n" );
  printf ( "  For a real tridiagonal system,\n" );
  printf ( "  R83_GS_SL solves a linear system using Gauss-Seidel iteration.\n" );
  printf ( "\n" );
  printf ( "  Matrix order N = %d\n", n );
  printf ( "  Iterations per call = %d\n", maxit );
  printf ( "\n" );
/*
  Set the matrix values.
*/
  a = ( double * ) malloc ( 3 * n * sizeof ( double ) );

  a[0+0*3] = 0.0;
  for ( j = 1; j < n; j++ )
  {
    a[0+j*3] = -1.0;
  }
  for ( j = 0; j < n; j++ )
  {
    a[1+j*3] = 2.0;
  }
  for ( j = 0; j < n-1; j++ )
  {
    a[2+j*3] = -1.0;
  }
  a[2+(n-1)*3] = 0.0;

  for ( job = 0; job <= 1; job++ )
  {
    if ( job == 0 )
    {
      printf ( "\n" );
      printf ( "  Solving A * x = b.\n" );
      printf ( "\n" );
    }
    else
    {
      printf ( "\n" );
      printf ( "  Solving A' * x = b.\n" );
      printf ( "\n" );
    }
/*
  Set the desired solution.
*/
    x = r8vec_indicator_new ( n );
/*
  Compute the corresponding right hand side.
*/
    if ( job == 0 )
    {
      b = r83_mxv ( n, a, x );
    }
    else
    {
      b = r83_vxm ( n, a, x );
    }
/*
  Set the starting solution.
*/
    for ( i = 0; i < n; i++ )
    {
      x[i] = 0.0;
    }
/*
  Solve the linear system.
*/
    for ( i = 1; i <= 3; i++ )
    {
      r83_gs_sl ( n, a, b, x, maxit, job );

      r8vec_print_some ( n, x, 1, 10, "  Current estimated solution:" );
    }

    free ( b );
    free ( x );
  }

  free ( a );

  return;
}
/******************************************************************************/

void test037 ( )

/******************************************************************************/
/*
  Purpose:

    TEST037 tests R83_INDICATOR.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    15 February 2013

  Author:

    John Burkardt
*/
{
  double *a;
  int n = 10;

  printf ( "\n" );
  printf ( "TEST037\n" );
  printf ( "  R83_INDICATOR sets up an S3 indicator matrix.\n" );
  printf ( "\n" );
  printf ( "  Matrix order N = %d\n", n );

  a = r83_indicator ( n );

  r83_print ( n, a, "  The S3 indicator matrix:" );

  free ( a );

  return;
}
/******************************************************************************/

void test04 ( )

/******************************************************************************/
/*
  Purpose:

    TEST04 tests R83_JAC_SL.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    03 March 2013

  Author:

    John Burkardt
*/
{
# define N 100

  double a[3*N];
  double *b;
  int i;
  int j;
  int job;
  int maxit = 1000;
  double *x;

  printf ( "\n" );
  printf ( "TEST04\n" );
  printf ( "  For a real tridiagonal system,\n" );
  printf ( "  R83_JAC_SL solves a linear system using Jacobi iteration.\n" );
  printf ( "\n" );
  printf ( "  Matrix order N = %d\n", N );
  printf ( "  Iterations per call = %d\n", maxit );
  printf ( "\n" );
/*
  Set the matrix values.
*/
  a[0+0*3] = 0.0;
  for ( j = 1; j < N; j++ )
  {
    a[0+j*3] = -1.0;
  }
  for ( j = 0; j < N; j++ )
  {
    a[1+j*3] = 2.0;
  }
  for ( j = 0; j < N-1; j++ )
  {
    a[2+j*3] = -1.0;
  }
  a[2+(N-1)*3] = 0.0;

  for ( job = 0; job <= 1; job++ )
  {
    if ( job == 0 )
    {
      printf ( "\n" );
      printf ( "  Solving A * x = b.\n" );
      printf ( "\n" );
    }
    else
    {
      printf ( "\n" );
      printf ( "  Solving A' * x = b.\n" );
      printf ( "\n" );
    }
/*
  Set the desired solution.
*/
    x = r8vec_indicator_new ( N );
/*
  Compute the corresponding right hand side.
*/
    if ( job == 0 )
    {
      b = r83_mxv ( N, a, x );
    }
    else
    {
      b = r83_vxm ( N, a, x );
    }

    r8vec_print_some ( N, b, 1, 10, "  The right hand side:" );
/*
  Set the starting solution.
*/
    for ( i = 0; i < N; i++ )
    {
      x[i] = 0.0;
    }
/*
  Solve the linear system.
*/
    for ( i = 1; i <= 3; i++ )
    {
      r83_jac_sl ( N, a, b, x, maxit, job );

      r8vec_print_some ( N, x, 1, 10, "  Current estimated solution:" );
    }

    free ( b );
    free ( x );
  }

  return;
# undef N
}
/******************************************************************************/

void test05 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST05 tests R83_NP_DET, R83_NP_FA.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    03 March 2013

  Author:

    John Burkardt
*/
{
# define N 10

  double *a;
  double *b;
  double det;
  int info;
  int seed = 123456789;

  printf ( "\n" );
  printf ( "TEST05\n" );
  printf ( "  For a tridiagonal matrix that can be factored\n" );
  printf ( "  with no pivoting,\n" );
  printf ( "  R83_NP_FA factors,\n" );
  printf ( "  R83_NP_DET computes the determinant.\n" );
  printf ( "\n" );
  printf ( "  Matrix order N = %d\n", N );
/*
  Set the matrix.
*/
  a = r83_random ( N, &seed );

  r83_print ( N, a, "  The S3 matrix:" );
/*
  Copy the matrix into general storage.
*/
  b = r83_to_r8ge ( N, a );
/*
  Factor the matrix.
*/
  info = r83_np_fa ( N, a );

  r83_print ( N, a, "  The factored S3 matrix:" );

  if ( info != 0 )
  {
    printf ( "\n" );
    printf ( "TEST05 - Warning!\n" );
    printf ( "  R83_NP_FA returns INFO = %d\n", info );
  }
/*
  Compute the determinant.
*/
  det = r83_np_det ( N, a );

  printf ( "\n" );
  printf ( "  R83_NP_DET computes determinant =  %g\n", det );
/*
  Factor the matrix in R8GE storage.
*/
  info = r8ge_np_fa ( N, b );

  if ( info != 0 )
  {
    printf ( "\n" );
    printf ( "TEST05 - Warning!\n" );
    printf ( "  R8GE_NP_FA returns INFO = %d\n", info );
  }
/*
  Compute the determinant of the matrix in R8GE storage.
*/
  det = r8ge_np_det ( N, b );

  printf ( "  R8GE_NP_DET computes determinant = %g\n", det );

  free ( a );
  free ( b );

  return;
# undef N
}
/******************************************************************************/

void test06 ( )

/******************************************************************************/
/*
 TEST06 tests R83_NP_FA, R83_NP_SL.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    05 March 2013

  Author:

    John Burkardt
*/
{
# define N 10

  double *a;
  double *b;
  int i;
  int info;
  int job;
  int seed = 123456789;
  double *x;

  printf ( "\n" );
  printf ( "TEST06\n" );
  printf ( "  For a tridiagonal matrix that can be factored\n" );
  printf ( "  with no pivoting,\n" );
  printf ( "  R83_NP_FA factors\n" );
  printf ( "  R83_NP_SL solves a factored system.\n" );
  printf ( "\n" );
  printf ( "  Matrix order N = %d\n", N );
/*
  Set the matrix.
*/
  a = r83_random ( N, &seed );

  r83_print ( N, a, "  The tridiagonal matrix:" );
/*
  Set the desired solution.
*/
  x = r8vec_indicator_new ( N );
/*
  Compute the corresponding right hand side.
*/
  b = r83_mxv ( N, a, x );
  free ( x );
/*
  Factor the matrix.
*/
  info = r83_np_fa ( N, a );

  if ( info != 0 )
  {
    printf ( "\n" );
    printf ( "TEST06 - Fatal error!\n" );
    printf ( "  The test matrix is singular.\n" );
    return;
  }
/*
  Solve the linear system.
*/
  job = 0;
  x = r83_np_sl ( N, a, b, job );
 
  r8vec_print ( N, x, "  Solution to A*x=b:" );
/*
  Set the desired solution
*/
  free ( x );
  x = r8vec_indicator_new ( N );
/*
  Compute the corresponding right hand side, using the factored matrix.
*/
  job = 1;
  free ( b );
  b = r83_np_ml ( N, a, x, job );
/*
  Solve the linear system.
*/
  job = 1;
  free ( x );
  x = r83_np_sl ( N, a, b, job );
 
  r8vec_print ( N, x, "  Solution to A'*x=b:" );
 
  free ( a );
  free ( b );
  free ( x );

  return;
# undef N
}
/******************************************************************************/

void test07 ( )

/******************************************************************************/
/*
  Purpose:

    TEST07 tests R83_NP_FS.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    04 March 2013

  Author:

    John Burkardt
*/
{
# define N 10

  double *a;
  double *b;
  int i;
  int seed = 123456789;
  double *x;

  printf ( "\n" );
  printf ( "TEST07\n" );
  printf ( "  R83_NP_FS factors and solves a tridiagonal\n" );
  printf ( "  linear system.\n" );
  printf ( "\n" );
  printf ( "  Matrix order N = %d\n", N );
/*
  Set the matrix elements.
*/
  a = r83_random ( N, &seed );
/*
  Set the desired solution.
*/
  x = r8vec_indicator_new ( N );
/*
  Compute b = A * x.
*/
  b = r83_mxv ( N, a, x );
/*
  Wipe out the solution.
*/
  free ( x );
/*
  Solve the system.
*/
  x = r83_np_fs ( N, a, b );

  r8vec_print ( N, x, "  Solution to A*x=b:" );

  free ( a );
  free ( b );
  free ( x );

  return;
# undef N
}
/******************************************************************************/

void test08 ( )

/******************************************************************************/
/*
  Purpose:

    TEST08 tests R83_NP_ML.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    04 March 2013

  Author:

    John Burkardt
*/
{
# define N 10

  double *a;
  double *b;
  double *b2;
  int info;
  int i;
  int job;
  int seed = 123456789;
  double *x;

  printf ( "\n" );
  printf ( "TEST08\n" );
  printf ( "  R83_NP_ML computes A*x or A'*x\n" );
  printf ( "  where A has been factored by R83_FA.\n" );
  printf ( "\n" );
  printf ( "  Matrix order N = %d\n", N );

  for ( job = 0; job <= 1; job++ )
  {
/*
  Set the matrix.
*/
    a = r83_random ( N, &seed );
/*
  Set the desired solution.
*/
    x = r8vec_indicator_new ( N );
/*
  Compute the corresponding right hand side.
*/
    if ( job == 0 )
    {
      b = r83_mxv ( N, a, x );
    }
    else
    {
      b = r83_vxm ( N, a, x );
    }
/*
  Factor the matrix.
*/
    info = r83_np_fa ( N, a );

    if ( info != 0 )
    {
      printf ( "\n" );
      printf ( "TEST08 - Fatal error!\n" );
      printf ( "  R83_NP_FA declares the matrix is singular!\n" );
      printf ( "  The value of INFO is %d\n", info );
      return;
    }
/*
  Now multiply factored matrix times solution to get right hand side again.
*/
    b2 = r83_np_ml ( N, a, x, job );

    if ( job == 0 )
    {
      r8vec2_print_some ( N, b, b2, 10, "  A*x and PLU*x:" );
    }
    else
    {
      r8vec2_print_some ( N, b, b2, 10, "  A'*x and (PLU)'*x" );
    }

    free ( a );
    free ( b );
    free ( b2 );
    free ( x );
  }

  return;
# undef N
}
/******************************************************************************/

void test09 ( )

/******************************************************************************/
/*
  Purpose:

    TEST09 tests R83P_DET.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    05 March 2013

  Author:

    John Burkardt
*/
{
# define N 12

  double *a;
  double *b;
  double det;
  int info;
  int pivot[N];
  int seed = 123456789;
  double work2[N-1];
  double work3[N-1];
  double work4;

  printf ( "\n" );
  printf ( "TEST09\n" );
  printf ( "  R83P_DET, determinant of a tridiagonal\n" );
  printf ( "  periodic matrix.\n" );
  printf ( "\n" );
  printf ( "  Matrix order N = %d\n", N );
/*
  Set the matrix.
*/
  a = r83p_random ( N, &seed );

  r83p_print ( N, a, "  The periodic tridiagonal matrix:" );
/*
  Copy the matrix into a general array.
*/
  b = r83p_to_r8ge ( N, a );
/*
  Factor the matrix.
*/
  info = r83p_fa ( N, a, work2, work3, &work4 );
/*
  Compute the determinant.
*/
  det = r83p_det ( N, a, work4 );

  printf ( "\n" );
  printf ( "  R83P_DET computes the determinant = %g\n", det );
/*
  Factor the general matrix.
*/
  info = r8ge_fa ( N, b, pivot );
/*
  Compute the determinant.
*/
  det = r8ge_det ( N, b, pivot );

  printf ( "  R8GE_DET computes the determinant = %g\n", det );

  free ( a );
  free ( b );

  return;
# undef N
}
/******************************************************************************/

void test10 ( )

/******************************************************************************/
/*
  Purpose:

    TEST10 tests R83P_FA, R83P_SL.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 March 2013

  Author:

    John Burkardt
*/
{
# define N 10

  double *a;
  double *b;
  int i;
  int info;
  int job;
  int seed = 123456789;
  double work2[N-1];
  double work3[N-1];
  double work4;
  double *x;

  printf ( "\n" );
  printf ( "TEST10\n" );
  printf ( "  R83P_FA factors a tridiagonal periodic system.\n" );
  printf ( "  R83P_SL solves a tridiagonal periodic system.\n" );
  printf ( "\n" );
  printf ( "  Matrix order N = %d\n", N );
/*
  Set the matrix.
*/
  a = r83p_random ( N, &seed );
/*
  Factor the matrix.
*/
  info = r83p_fa ( N, a, work2, work3, &work4 );

  if ( info != 0 )
  {
    printf ( "\n" );
    printf ( "  TEST10 - Fatal error!\n" );
    printf ( "  R83P_FA returns INFO = %d\n", info );
    return;
  }

  for ( job = 0; job <= 1; job++ )
  {
/*
  Set the desired solution.
*/
    x = r8vec_indicator_new ( N );
/*
  Compute the corresponding right hand side.
*/
    b = r83p_ml ( N, a, x, job );
/*
  Solve the linear system.
*/
    free ( x );

    x = r83p_sl ( N, a, b, job, work2, work3, work4 );

    if ( job == 0 )
    {
      r8vec_print ( N, x, "  Solution to A*x=b:" );
    }
    else
    {
      r8vec_print ( N, x, "  Solution to A'*x=b:" );
    }

    free ( x );
    free ( b );
  }
 
  free ( a );

  return;
# undef N
}
/******************************************************************************/

void test105 ( )

/******************************************************************************/
/*
  Purpose:

    TEST105 tests R83P_INDICATOR.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    07 March 2013

  Author:

    John Burkardt
*/
{
# define N 10

  double *a;

  printf ( "\n" );
  printf ( "TEST105\n" );
  printf ( "  R83P_INDICATOR sets up an S3P indicator matrix.\n" );
  printf ( "\n" );
  printf ( "  Matrix order N = %d\n", N );

  a = r83p_indicator ( N );

  r83p_print ( N, a, "  The S3P indicator matrix:" );

  free ( a );

  return;
# undef N
}
/******************************************************************************/

void test11 ( )

/******************************************************************************/
/*
  Purpose:

    TEST11 tests R83P_ML.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    07 March 2013

  Author:

    John Burkardt
*/
{
# define N 10

  double *a;
  double *b;
  double *b2;
  int info;
  int i;
  int job;
  int seed = 123456789;
  double work2[N-1];
  double work3[N-1];
  double work4;
  double *x;

  printf ( "\n" );
  printf ( "TEST11\n" );
  printf ( "  R83P_ML computes A*x or A'*X\n" );
  printf ( "  where A has been factored by R83P_FA.\n" );
  printf ( "\n" );
  printf ( "  Matrix order N = %d\n", N );

  for ( job = 0; job <= 1; job++ )
  {
/*
  Set the matrix.
*/
    a = r83p_random ( N, &seed );
/*
  Set the desired solution.
*/
    x = r8vec_indicator_new ( N );
/*
  Compute the corresponding right hand side.
*/
    if ( job == 0 )
    {
      b = r83p_mxv ( N, a, x );
    }
    else
    {
      b = r83p_vxm ( N, a, x );
    }
/*
  Factor the matrix.
*/
    info = r83p_fa ( N, a, work2, work3, &work4 );

    if ( info != 0 )
    {
      printf ( "\n" );
      printf ( "TEST11 - Fatal error!\n" );
      printf ( "  R83P_FA declares the matrix is singular!\n" );
      printf ( "  The value of INFO is %d\n", info );
      return;
    }
/*
  Now multiply factored matrix times solution to get right hand side again.
*/
    b2 = r83p_ml ( N, a, x, job );

    if ( job == 0 )
    {
      r8vec2_print_some ( N, b, b2, 10, "  A*x and PLU*x" );
    }
    else
    {
      r8vec2_print_some ( N, b, b2, 10, "  A'*x and (PLU)'*x" );
    }

    free ( a );
    free ( x );
    free ( b );
    free ( b2 );
  }

  return;
# undef N
}
/******************************************************************************/

void test115 ( )

/******************************************************************************/
/*
  Purpose:

    TEST115 tests R85_INDICATOR.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    08 March 2013

  Author:

    John Burkardt
*/
{
# define N 10

  double *a;

  printf ( "\n" );
  printf ( "TEST115\n" );
  printf ( "  R85_INDICATOR sets up an R85 indicator matrix.\n" );
  printf ( "\n" );
  printf ( "  Matrix order N = %d\n", N );

  a = r85_indicator ( N );

  r85_print ( N, a, "  The R85 indicator matrix:" );

  free ( a );

  return;
# undef N
}
/******************************************************************************/

void test12 ( )

/******************************************************************************/
/*
  Purpose:

    TEST12 tests R85_NP_FS.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    08 March 2013

  Author:

    John Burkardt
*/
{
# define N 10

  double *a;
  double *b;
  int i;
  int seed = 123456789;
  double *x;

  printf ( "\n" );
  printf ( "TEST12\n" );
  printf ( "  R85_NP_FS factors and solves a pentadiagonal\n" );
  printf ( "  linear system, with no pivoting.\n" );
  printf ( "\n" );
  printf ( "  Matrix order N = %d\n", N );
/*
  Set the matrix to a random value.
*/
  a = r85_random ( N, &seed );

  r85_print ( N, a, "  The pentadiagonal matrix:" );
/*
  Set the desired solution.
*/
  x = r8vec_indicator_new ( N );
/*
  Compute b = A * x.
*/
  b = r85_mxv ( N, a, x );

  r8vec_print ( N, b, "  Right hand side:" );
/*
  Wipe out the solution.
*/
  for ( i = 0; i < N; i++ )
  {
    x[i] = 0.0;
  }
/*
  Solve the system.
*/
  free ( x );
  x = r85_np_fs ( N, a, b );

  r8vec_print ( N, x, "  Solution to A*x=b:" );

  free ( a );
  free ( b );
  free ( x );

  return;
# undef N
}
/******************************************************************************/

void test13 ( )

/******************************************************************************/
/*
  Purpose:

    TEST13 tests R8BB_FA, R8BB_PRINT, R8BB_RANDOM, R8BB_SL.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    08 March 2013

  Author:

    John Burkardt
*/
{
# define N1 8
# define N2 2
# define ML 1
# define MU 1
# define N N1+N2

  double *a;
  double *b;
  int i;
  int info;
  int pivot[N];
  int seed = 123456789;
  double *x;

  printf ( "\n" );
  printf ( "TEST13\n" );
  printf ( "  For a border banded matrix:\n" );
  printf ( "  R8BB_FA factors;\n" );
  printf ( "  R8BB_PRINT prints;\n" );
  printf ( "  R8BB_RANDOM randomizes;\n" );
  printf ( "  R8BB_SL solves.\n" );
  printf ( "\n" );
  printf ( "  Matrix order N     = %d\n", N );
  printf ( "  Matrix suborder N1 = %d\n", N1 );
  printf ( "  Matrix suborder N2 = %d\n", N2 );
  printf ( "  Lower bandwidth ML = %d\n", ML );
  printf ( "  Upper bandwidth MU = %d\n", MU );
/*
  Set the matrix.
*/
  a = r8bb_random ( N1, N2, ML, MU, &seed );

  r8bb_print ( N1, N2, ML, MU, a, "  The border-banded matrix:" );
/*
  Set the desired solution.
*/
  x = r8vec_indicator_new ( N );
/*
  Compute the corresponding right hand side.
*/
  b = r8bb_mxv ( N1, N2, ML, MU, a, x );

  r8vec_print ( N, b, "  The right hand side vector:" );
/*
  Factor the matrix.
*/
  info = r8bb_fa ( N1, N2, ML, MU, a, pivot );

  r8bb_print ( N1, N2, ML, MU, a, "  The FACTORED border-banded matrix:" );

  if ( info != 0 )
  {
    printf ( "\n" );
    printf ( "TEST13 - Fatal error!\n" );
    printf ( "  R8BB_FA claims the matrix is singular.\n" );
    printf ( "  The value of INFO is %d\n", info );
    return;
  }
/*
  Solve the system.
*/
  free ( x );
  x = r8bb_sl ( N1, N2, ML, MU, a, pivot, b );

  r8vec_print ( N, x, "  Solution to A*x=b:" );

  free ( a );
  free ( b );
  free ( x );

  return;
# undef ML
# undef MU
# undef N
# undef N1
# undef N2
}
/******************************************************************************/

void test14 ( )

/******************************************************************************/
/*
  Purpose:

    TEST14 tests R8BB_FA, R8BB_SL.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    08 March 2013

  Author:

    John Burkardt
*/
{
# define N1 0
# define N2 10
# define N N1+N2
# define ML 0
# define MU 0

  double *a;
  double *b;
  int i;
  int info;
  int pivot[N];
  int seed = 123456789;
  double *x;

  printf ( "\n" );
  printf ( "TEST14\n" );
  printf ( "  For a border banded matrix:\n" );
  printf ( "  R8BB_FA factors;\n" );
  printf ( "  R8BB_SL solves.\n" );
  printf ( "\n" );
  printf ( "  Matrix order N     = %d\n", N );
  printf ( "  Matrix suborder N1 = %d\n", N1 );
  printf ( "  Matrix suborder N2 = %d\n", N2 );
  printf ( "  Lower bandwidth ML = %d\n", ML );
  printf ( "  Upper bandwidth MU = %d\n", MU );
/*
  Set the matrix.
*/
  a = r8bb_random ( N1, N2, ML, MU, &seed );
/*
  Set the desired solution.
*/
  x = r8vec_indicator_new ( N );
/*
  Compute the corresponding right hand side.
*/
  b = r8bb_mxv ( N1, N2, ML, MU, a, x );
/*
  Factor the matrix.
*/
  info = r8bb_fa ( N1, N2, ML, MU, a, pivot );

  if ( info != 0 )
  {
    printf ( "\n" );
    printf ( "TEST14 - Fatal error!\n" );
    printf ( "  R8BB_FA claims the matrix is singular.\n" );
    printf ( "  The value of INFO is %d\n", info );
    return;
  }
/*
  Solve the system.
*/
  free ( x );
  x = r8bb_sl ( N1, N2, ML, MU, a, pivot, b );

  r8vec_print ( N, x, "  Solution to A*x=b:" );
 
  free ( a );
  free ( b ); 
  free ( x );

  return;
# undef ML
# undef MU
# undef N
# undef N1
# undef N2
}
/******************************************************************************/

void test15 ( )

/******************************************************************************/
/*
  Purpose:

    TEST15 tests R8BB_FA, R8BB_SL.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    09 March 2013

  Author:

    John Burkardt
*/
{
# define N1 10
# define N2 0
# define N N1+N2
# define ML 1
# define MU 1

  double *a;
  double *b;
  int i;
  int info;
  int pivot[N];
  int seed = 123456789;
  double *x;

  printf ( "\n" );
  printf ( "TEST15\n" );
  printf ( "  For a border banded matrix:\n" );
  printf ( "  R8BB_FA factors;\n" );
  printf ( "  R8BB_SL solves.\n" );
  printf ( "\n" );
  printf ( "  Matrix order N     = %d\n", N );
  printf ( "  Matrix suborder N1 = %d\n", N1 );
  printf ( "  Matrix suborder N2 = %d\n", N2 );
  printf ( "  Lower bandwidth ML = %d\n", ML );
  printf ( "  Upper bandwidth MU = %d\n", MU );
/*
  Set the matrix.
*/
  a = r8bb_random ( N1, N2, ML, MU, &seed );
/*
  Set the desired solution.
*/
  x = r8vec_indicator_new ( N );
/*
  Compute the corresponding right hand side.
*/
  b = r8bb_mxv ( N1, N2, ML, MU, a, x );
/*
  Factor the matrix
*/
  info = r8bb_fa ( N1, N2, ML, MU, a, pivot );
 
  if ( info != 0 )
  {
    printf ( "\n" );
    printf ( "TEST15 - Fatal error!\n" );
    printf ( "  R8BB_FA claims the matrix is singular.\n" );
    printf ( "  The value of INFO is %d\n", info );
    return;
  }
/*
  Solve the system.
*/
  free ( x );
  x = r8bb_sl ( N1, N2, ML, MU, a, pivot, b );

  r8vec_print ( N, x, "  Solution to A*x=b:" );
 
  free ( a );
  free ( b );
  free ( x );

  return;
# undef ML
# undef MU
# undef N
# undef N1
# undef N2
}
/******************************************************************************/

void test151 ( )

/******************************************************************************/
/*
  Purpose:

    TEST151 tests R8BB_INDICATOR.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    09 March 2013

  Author:

    John Burkardt
*/
{
# define N1 6
# define N2 2
# define ML 1
# define MU 1

  double *a;

  printf ( "\n" );
  printf ( "TEST151\n" );
  printf ( "  R8BB_INDICATOR sets up an R8BB indicator matrix.\n" );
  printf ( "\n" );
  printf ( "  Matrix order N     = %d\n", N1 + N2 );
  printf ( "  Matrix suborder N1 = %d\n", N1 );
  printf ( "  Matrix suborder N2 = %d\n", N2 );
  printf ( "  Lower bandwidth ML = %d\n", ML );
  printf ( "  Upper bandwidth MU = %d\n",  MU );

  a = r8bb_indicator ( N1, N2, ML, MU );

  r8bb_print ( N1, N2, ML, MU, a, "  The R8BB indicator matrix:" );

  free ( a );

  return;
# undef ML
# undef MU
# undef N1
# undef N2
}
/******************************************************************************/

void test152 ( )

/******************************************************************************/
/*
  Purpose:

    TEST152 tests R8BLT_INDICATOR.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    10 March 2013

  Author:

    John Burkardt
*/
{
  double *a;
  int n = 6;
  int ml = 2;

  printf ( "\n" );
  printf ( "TEST152\n" );
  printf ( "  R8BLT_INDICATOR sets up an R8BLT indicator matrix.\n" );
  printf ( "\n" );
  printf ( "  Matrix order N = %d\n", n );
  printf ( "  Lower bandwidth ML = %d\n", ml );

  a = r8blt_indicator ( n, ml );

  r8blt_print ( n, ml, a, "  The R8BLT indicator matrix:" );

  free ( a );

  return;
}
/******************************************************************************/

void test153 ( )

/******************************************************************************/
/*
  Purpose:

    TEST153 tests R8BLT_MXV, R8BLT_PRINT, R8BLT_RANDOM, R8BLT_SL, R8BLT_VXM.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    10 March 2013

  Author:

    John Burkardt
*/
{
# define ML 3
# define N 10

  double *a;
  double *b;
  int i;
  int j;
  int job;
  int seed = 123456789;
  double *x;

  printf ( "\n" );
  printf ( "TEST153\n" );
  printf ( "  For a band matrix in lower triangular storage,\n" );
  printf ( "  R8BLT_RANDOM sets a random value;\n" );
  printf ( "  R8BLT_SL solves systems;\n" );
  printf ( "  R8BLT_MXV computes A*x;\n" );
  printf ( "  R8BLT_VXM computes A'*x;\n" );
  printf ( "\n" );
  printf ( "  Matrix order N = %d\n", N );
  printf ( "  Lower bandwidth ML = %d\n", ML );

  a = r8blt_random ( N, ML, &seed );

  r8blt_print ( N, ML, a, "  The R8BLT matrix:" );

  for ( job = 0; job <= 1; job++ )
  {
/*
  Set the desired solution.
*/
    x = r8vec_indicator_new ( N );
/*
  Compute the corresponding right hand side.
*/
    if ( job == 0 )
    {
      b = r8blt_mxv ( N, ML, a, x );
    }
    else
    {
      b = r8blt_vxm ( N, ML, a, x );
    }

    r8vec_print ( N, b, "  The right hand side:" );
/*
  Solve the linear system.
*/
    free ( x );
    x = r8blt_sl ( N, ML, a, b, job );
 
    if ( job == 0 )
    {
      r8vec_print ( N, x, "  Solution to A*x=b:" );
    }
    else
    {
      r8vec_print ( N, x, "  Solution to A'*x=b:" );
    }

    free ( b );
    free ( x );
  }

  free ( a );

  return;
# undef ML
# undef N
}
/******************************************************************************/

void test154 ( )

/******************************************************************************/
/*
  Purpose:

    TEST154 tests R8BTO_INDICATOR.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 March 2013

  Author:

    John Burkardt
*/
{
# define L 3
# define M 2

  double *a;

  printf ( "\n" );
  printf ( "TEST154\n" );
  printf ( "  For a real block Toeplitz matrix,\n" );
  printf ( "  R8BTO_INDICATOR sets up an indicator matrix.\n" );
  printf ( "\n" );
  printf ( "  Block order M =  %d\n", M );
  printf ( "  Block number L = %d\n", L );
  printf ( "  Matrix order N = %d\n", M * L );

  a = r8bto_indicator ( M, L );

  r8bto_print ( M, L, a, "  The block Toeplitz matrix:" );

  free ( a );

  return;
# undef L
# undef M
}
/******************************************************************************/

void test155 ( )

/******************************************************************************/
/*
  Purpose:

    TEST155 tests R8BTO_MXV, R8BTO_VXM.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 March 2013

  Author:

    John Burkardt
*/
{
# define L 3
# define M 2
# define N ( M * L )

  double a[M*M*(2*L-1)] = {
    1.0, 5.0, 2.0, 5.0, 
    3.0, 6.0, 4.0, 6.0, 
    5.0, 7.0, 6.0, 7.0, 
    7.0, 8.0, 8.0, 8.0, 
    9.0, 9.0, 0.0, 9.0 };
  double *b;
  double *x;

  printf ( "\n" );
  printf ( "TEST155\n" );
  printf ( "  For a real block Toeplitz matrix,\n" );
  printf ( "  R8BTO_MXV computes A * x.\n" );
  printf ( "  R8BTO_VXM computes A'* x.\n" );
  printf ( "\n" );
  printf ( "  Block order M =  %d\n", M );
  printf ( "  Block number L = %d\n", L );
  printf ( "  Matrix order N = %d\n", N );

  r8bto_print ( M, L, a, "  The block Toeplitz matrix:" );

  x = r8ge_indicator ( M, L );

  r8ge_print ( M, L, x, "  The 'vector' x:" );

  b = r8bto_mxv ( M, L, a, x );

  r8ge_print ( M, L, b, "  The product A*x:" );

  free ( b );
  b = r8bto_vxm ( M, L, a, x );

  r8ge_print ( M, L, b, "  The product A'*x:" );

  free ( b );
  free ( x );

  return;
# undef L
# undef M
# undef N
}
/******************************************************************************/

void test1565 ( )

/******************************************************************************/
/*
  Purpose:

    TEST1565 tests R8BUT_INDICATOR.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 March 2013

  Author:

    John Burkardt
*/
{
# define N 6
# define MU 2

  double *a;

  printf ( "\n" );
  printf ( "TEST1565\n" );
  printf ( "  R8BUT_INDICATOR sets up an SBUT indicator matrix.\n" );
  printf ( "\n" );
  printf ( "  Matrix order N     = %d\n", N );
  printf ( "  Upper bandwidth MU = %d\n", MU );
/*
  Set the matrix.
*/
  a = r8but_indicator ( N, MU );

  r8but_print ( N, MU, a, "  The R8BUT indicator matrix:" );

  free ( a );

  return;
# undef N
# undef MU
}
/******************************************************************************/

void test1566 ( )

/******************************************************************************/
/*
  Purpose:

    TEST1566 tests R8BUT_MXV, R8BUT_PRINT, R8BUT_RANDOM, R8BUT_SL, R8BUT_VXM.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 March 2013

  Author:

    John Burkardt
*/
{
# define MU 3
# define N 10

  double *a;
  double *b;
  int i;
  int j;
  int job;
  int seed = 123456789;
  double *x;

  printf ( "\n" );
  printf ( "TEST1566\n" );
  printf ( "  For a band matrix in upper triangular storage,\n" );
  printf ( "  R8BUT_RANDOM sets a random value;\n" );
  printf ( "  R8BUT_SL solves systems;\n" );
  printf ( "  R8BUT_MXV computes matrix-vector products;\n" );
  printf ( "  R8BUT_VXM computes vector-matrix products;\n" );
  printf ( "\n" );
  printf ( "  Matrix order N =     %d\n", N );
  printf ( "  Upper bandwidth MU = %d\n", MU );

  a = r8but_random ( N, MU, &seed );

  r8but_print ( N, MU, a, "  The R8BUT matrix:" );

  for ( job = 0; job <= 1; job++ )
  {
/*
  Set the desired solution.
*/
    x = r8vec_indicator_new ( N );
/*
  Compute the corresponding right hand side.
*/
    if ( job == 0 )
    {
      b = r8but_mxv ( N, MU, a, x );
    }
    else
    {
      b = r8but_vxm ( N, MU, a, x );
    }

    r8vec_print ( N, b, "  The right hand side:" );
/*
  Solve the linear system.
*/
    free ( x );
    x = r8but_sl ( N, MU, a, b, job );
 
    if ( job == 0 )
    {
      r8vec_print ( N, x, "  Solution to the untransposed system:" );
    }
    else
    {
      r8vec_print ( N, x, "  Solution to the transposed system:" );
    }
    free ( b );
    free ( x );

  }

  free ( a );

  return;
# undef MU
# undef N
}
/******************************************************************************/

void test157 ( )

/******************************************************************************/
/*
  Purpose:

    TEST157 tests R8CB_INDICATOR.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    13 March 2013

  Author:

    John Burkardt
*/
{
# define M 8
# define N 10
# define ML 2
# define MU 3

  double *a;

  printf ( "\n" );
  printf ( "TEST157\n" );
  printf ( "  R8CB_INDICATOR sets up an R8CB indicator matrix.\n" );
  printf ( "\n" );
  printf ( "  Matrix rows M =    %d\n", M );
  printf ( "  Matrix columns N = %d\n", N );
  printf ( "  Lower bandwidth ML = %d\n", ML );
  printf ( "  Upper bandwidth MU = %d\n", MU );

  a = r8cb_indicator ( M, N, ML, MU );

  r8cb_print ( M, N, ML, MU, a, "  The R8CB indicator matrix:" );

  free ( a );

  return;
# undef M
# undef N
# undef ML
# undef MU
}
/******************************************************************************/

void test16 ( )

/******************************************************************************/
/*
  Purpose:

    TEST16 tests R8CB_NP_FA, R8CB_DET.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    13 March 2013

  Author:

    John Burkardt
*/
{
# define N 10
# define ML 2
# define MU 3

  double *a;
  double *a2;
  double det;
  int info;
  int pivot[N];
  int seed = 123456789;

  printf ( "\n" );
  printf ( "TEST16\n" );
  printf ( "  For a compact band matrix, no pivoting:\n" );
  printf ( "  R8CB_NP_FA factors;\n" );
  printf ( "  R8CB_DET computes the determinant;\n" );
  printf ( "\n" );
  printf ( "  Matrix order N     = %d\n", N );
  printf ( "  Lower bandwidth ML = %d\n", ML );
  printf ( "  Upper bandwidth MU = %d\n", MU );
/*
  Set the matrix.
*/
  a = r8cb_random ( N, ML, MU, &seed );

  r8cb_print ( N, N, ML, MU, a, "  The compact band matrix:" );
/*
  Copy the matrix into a general array.
*/
  a2 = r8cb_to_r8ge ( N, ML, MU, a );
/*
  Factor the matrix.
*/
  info = r8cb_np_fa ( N, ML, MU, a );
/*
  Compute the determinant.
*/
  det = r8cb_det ( N, ML, MU, a );

  printf ( "\n" );
  printf ( "  R8CB_DET computes the determinant = %g\n", det );
/*
  Factor the general matrix.
*/
  info = r8ge_fa ( N, a2, pivot );
/*
  Compute the determinant.
*/
  det = r8ge_det ( N, a2, pivot );

  printf ( "  R8GE_DET computes the determinant = %g\n", det );

  free ( a );
  free ( a2 );

  return;
# undef N
# undef ML
# undef MU
}
/******************************************************************************/

void test17 ( )

/******************************************************************************/
/*
  Purpose:

    TEST17 tests R8CB_NP_FA, R8CB_NP_SL.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 March 2013

  Author:

    John Burkardt
*/
{
# define N 10
# define ML 1
# define MU 2

  double *a;
  double *b;
  int i;
  int info;
  int job;
  int seed = 123456789;
  double *x;

  printf ( "\n" );
  printf ( "TEST17\n" );
  printf ( "  For a compact band matrix, no pivoting:\n" );
  printf ( "  R8CB_NP_FA factors;\n" );
  printf ( "  R8CB_NP_SL solves.\n" );
  printf ( "\n" );
  printf ( "  Matrix order N     = %d\n", N );
  printf ( "  Lower bandwidth ML = %d\n", ML );
  printf ( "  Upper bandwidth MU = %d\n", MU );

  for ( job = 0; job <= 1; job++ )
  {
/*
  Set the matrix.
*/
    a = r8cb_random ( N, ML, MU, &seed );
/*
  Set the desired solution.
*/
    x = r8vec_indicator_new ( N );
/*
  Compute the right hand side.
*/
    if ( job == 0 )
    {
      b = r8cb_mxv ( N, ML, MU, a, x );
    }
    else
    {
      b = r8cb_vxm ( N, ML, MU, a, x );
    }
/*
  Factor the matrix.
*/
    info = r8cb_np_fa ( N, ML, MU, a );

    if ( info != 0 )
    {
      printf ( "\n" );
      printf ( "TEST17 - Fatal error!\n" );
      printf ( "  R8CB_NP_FA claims the matrix is singular.\n" );
      printf ( "  The value of info is %d\n", info );
      return;
    }
/*
  Solve the system.
*/
    free ( x );

    x = r8cb_np_sl ( N, ML, MU, a, b, job );

    if ( job == 0 )
    {
      r8vec_print ( N, x, "  Solution to A*x=b:" );
    }
    else
    {
      r8vec_print ( N, x, "  Solution to A'x=b:" );
    }

    free ( a );
    free ( b );
    free ( x );
  }

  return;
# undef N
# undef ML
# undef MU
}
/******************************************************************************/

void test18 ( )

/******************************************************************************/
/*
  Purpose:

    TEST18 tests R8CB_ML.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 March 2013

  Author:

    John Burkardt
*/
{
# define N 10
# define ML 1
# define MU 2

  double *a;
  double *b;
  double *b2;
  int i;
  int info;
  int job;
  int seed = 123456789;
  double *x;

  printf ( "\n" );
  printf ( "TEST18\n" );
  printf ( "  For a compact band matrix:\n" );
  printf ( "  R8CB_ML computes A*x or A'*X\n" );
  printf ( "  where A has been factored by R8CB_FA.\n" );
  printf ( "\n" );
  printf ( "  Matrix order N     = %d\n", N );
  printf ( "  Lower bandwidth ML = %d\n", ML );
  printf ( "  Upper bandwidth MU = %d\n", MU );

  for ( job = 0; job <= 1; job++ )
  {
/*
  Set the matrix.
*/
    a = r8cb_random ( N, ML, MU, &seed );
/*
  Set the desired solution.
*/
    x = r8vec_indicator_new ( N );
/*
  Compute the corresponding right hand side.
*/
    if ( job == 0 )
    {
      b = r8cb_mxv ( N, ML, MU, a, x );
    }
    else
    {
      b = r8cb_vxm ( N, ML, MU, a, x );
    }
/*
  Factor the matrix.
*/
    info = r8cb_np_fa ( N, ML, MU, a );

    if ( info != 0 )
    {
      printf ( "\n" );
      printf ( "TEST18 - Fatal error!\n" );
      printf ( "  R8CB_FA declares the matrix is singular!\n" );
      printf ( "  The value of INFO is %d\n", info );
      return;
    }
/*
  Now multiply factored matrix times solution to get right hand side again.
*/
    b2 = r8cb_ml ( N, ML, MU, a, x, job );

    if ( job == 0 )
    {
      r8vec2_print_some ( N, b, b2, 10, "  A*x and PLU*x" );
    }
    else
    {
      r8vec2_print_some ( N, b, b2, 10, "  A'*x and (PLU)'*x" );
    }

    free ( a );
    free ( b );
    free ( b2 );
    free ( x );
  }
  return;
# undef N
# undef ML
# undef MU
}
/******************************************************************************/

void test19 ( )

/******************************************************************************/
/*
  Purpose:

    TEST19 tests R8CBB_FA, R8CBB_PRINT, R8CBB_RANDOM, R8CBB_SL.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 March 2013

  Author:

    John Burkardt
*/
{
# define N1 8
# define N2 2
# define ML 1
# define MU 1

  double *a;
  double *b;
  int i;
  int info;
  int seed = 123456789;
  double *x;

  printf ( "\n" );
  printf ( "TEST19\n" );
  printf ( "  For a compressed border banded matrix:\n" );
  printf ( "  R8CBB_RANDOM randomly generates;\n" );
  printf ( "  R8CBB_PRINT prints;\n" );
  printf ( "  R8CBB_FA factors (no pivoting);\n" );
  printf ( "  R8CBB_SL solves.\n" );
  printf ( "\n" );
  printf ( "  Matrix order N     = %d\n", N1 + N2 );
  printf ( "  Matrix suborder N1 = %d\n", N1 );
  printf ( "  Matrix suborder N2 = %d\n", N2 );
  printf ( "  Lower bandwidth ML = %d\n", ML );
  printf ( "  Upper bandwidth MU = %d\n", MU );
/*
  Set the matrix.
*/
  a = r8cbb_random ( N1, N2, ML, MU, &seed );

  r8cbb_print ( N1, N2, ML, MU, a, "  The R8CBB matrix:" );
/*
  Set the desired solution.
*/
  x = r8vec_indicator_new ( N1 + N2 );
/*
  Compute the corresponding right hand side.
*/
  b = r8cbb_mxv ( N1, N2, ML, MU, a, x );
/*
  Factor the matrix
*/
  info = r8cbb_fa ( N1, N2, ML, MU, a );

  r8cbb_print ( N1, N2, ML, MU, a, "  The factored R8CBB matrix:" );

  if ( info != 0 )
  {
    printf ( "\n" );
    printf ( "  Fatal error!\n" );
    printf ( "  R8CBB_FA claims the matrix is singular.\n" );
    printf ( "  The value of INFO is %d\n", info );
    return;
  }
/*
  Solve the system.
*/
  r8vec_print ( N1 + N2, b, "  The right hand side vector b:" );

  free ( x );
  x = r8cbb_sl ( N1, N2, ML, MU, a, b );

  r8vec_print ( N1 + N2, x, "  Solution to A*x=b:" );

  free ( a );
  free ( b );
  free ( x );

  return;
# undef ML
# undef MU
# undef N1
# undef N2
}
/******************************************************************************/

void test193 ( )

/******************************************************************************/
/*
  Purpose:

    TEST193 tests R8CBB_INDICATOR.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    18 March 2013

  Author:

    John Burkardt
*/
{
# define N1 6
# define N2 2
# define ML 1
# define MU 1

  double *a;

  printf ( "\n" );
  printf ( "TEST193\n" );
  printf ( "  R8CBB_INDICATOR sets up an R8CBB indicator matrix.\n" );
  printf ( "\n" );
  printf ( "  Matrix order N     = %d\n", N1 + N2 );
  printf ( "  Matrix suborder N1 = %d\n", N1 );
  printf ( "  Matrix suborder N2 = %d\n", N2 );
  printf ( "  Lower bandwidth ML = %d\n", ML );
  printf ( "  Upper bandwidth MU = %d\n", MU );

  a = r8cbb_indicator ( N1, N2, ML, MU );

  r8cbb_print ( N1, N2, ML, MU, a, "  The R8CBB indicator matrix:" );

  free ( a );

  return;
# undef ML
# undef MU
# undef N1
# undef N2
}
/******************************************************************************/

void test195 ( )

/******************************************************************************/
/*
  Purpose:

    TEST195 tests R8CI_EVAL.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    18 March 2013

  Author:

    John Burkardt
*/
{
# define N 5

  double *a;
  double complex *lambda;
  int seed = 123456789;

  printf ( "\n" );
  printf ( "TEST195\n" );
  printf ( "  R8CI_EVAL finds the eigenvalues of \n" );
  printf ( "  a real circulant system.\n" );
  printf ( "\n" );
  printf ( "  Matrix order N = %d\n", N );
/*
  Set the matrix.
*/
  a = r8ci_random ( N, &seed );

  r8ci_print ( N, a, "  The circulant matrix:" );

  lambda = r8ci_eval ( N, a );

  c8vec_print ( N, lambda, "  The eigenvalues:" );

  free ( a );
  free ( lambda );

  return;
# undef N
}
/******************************************************************************/

void test197 ( )

/******************************************************************************/
/*
  Purpose:

    TEST197 tests R8CI_INDICATOR.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    18 March 2013

  Author:

    John Burkardt
*/
{
# define N 5

  double *a;

  printf ( "\n" );
  printf ( "TEST197\n" );
  printf ( "  R8CI_INDICATOR sets up an R8CI indicator matrix.\n" );
  printf ( "\n" );
  printf ( "  Matrix order N = %d\n", N );

  a = r8ci_indicator ( N );

  r8ci_print ( N, a, "  The circulant matrix:" );

  free ( a );

  return;
# undef N
}
/******************************************************************************/

void test20 ( )

/******************************************************************************/
/*
  Purpose:

    TEST20 tests R8CI_SL.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    19 March 2013

  Author:

    John Burkardt
*/
{
  double *a;
  double *b;
  int i;
  int job;
  int n = 10;
  int seed = 123456789;
  double *x;

  printf ( "\n" );
  printf ( "TEST20\n" );
  printf ( "  R8CI_SL solves a circulant system.\n" );
  printf ( "\n" );
  printf ( "  Matrix order N = %d\n", n );
/*
  Set the matrix.
*/
  a = r8ci_random ( n, &seed );

  r8ci_print ( n, a, "  The circulant matrix:" );

  for ( job = 0; job <= 1; job++ )
  {
/*
  Set the desired solution.
*/
    x = r8vec_indicator_new ( n );
/*
  Compute the corresponding right hand side.
*/
    if ( job == 0 )
    {
      b = r8ci_mxv ( n, a, x );
    }
    else
    {
      b = r8ci_vxm ( n, a, x );
    }
/*
  Solve the linear system.
*/
    free ( x );
    x = r8ci_sl ( n, a, b, job );

    if ( job == 0 )
    {
      r8vec_print ( n, x, "  Solution to A*x=b:" );
    }
    else
    {
      r8vec_print ( n, x, "  Solution to A'*x=b:" );
    }

    free ( b );
    free ( x );

  }
 
  free ( a );

  return;
}
/******************************************************************************/

void test21 ( )

/******************************************************************************/
/*
  Purpose:

    TEST21 tests R8GB_DET.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    20 March 2013

  Author:

    John Burkardt
*/
{
# define M 10
# define ML 3
# define MU 2
# define N 10

  double *a;
  double *a2;
  double det;
  int info;
  int pivot[N];
  int seed = 123456789;

  printf ( "\n" );
  printf ( "TEST21\n" );
  printf ( "  For a general banded matrix,\n" );
  printf ( "  R8GB_DET computes the determinant.\n" );
  printf ( "\n" );
  printf ( "  Matrix rows M =    %d\n", M );
  printf ( "  Matrix columns N = %d\n", N );
  printf ( "  Lower bandwidth ML  = %d\n", ML );
  printf ( "  Upper bandwidth MU  = %d\n", MU );
/*
  Set the matrix.
*/
  a = r8gb_random ( M, N, ML, MU, &seed );
/*
  Copy the matrix into a general array.
*/
  a2 = r8gb_to_r8ge ( M, N, ML, MU, a );
/*
  Print the matrix just to show what it looks like.
*/
  r8ge_print ( M, N, a, "  The banded matrix:" );
/*
  Factor the matrix.
*/
  info = r8gb_fa ( N, ML, MU, a, pivot );
/*
  Compute the determinant.
*/
  det = r8gb_det ( N, ML, MU, a, pivot );

  printf ( "\n" );
  printf ( "  R8GB_DET computes the determinant = %g\n", det );
/*
  Factor the general matrix.
*/
  info = r8ge_fa ( N, a2, pivot );
/*
  Compute the determinant.
*/
  det = r8ge_det ( N, a2, pivot );

  printf ( "  R8GE_DET computes the determinant = %g\n", det );

  free ( a );
  free ( a2 );

  return;
# undef M
# undef ML
# undef MU
# undef N
}
/******************************************************************************/

void test22 ( )

/******************************************************************************/
/*
  Purpose:

    TEST22 tests R8GB_FA, R8GB_SL.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    21 March 2013

  Author:

    John Burkardt
*/
{
# define M 5
# define ML 1
# define MU 2
# define N 5

  double *a;
  double *b;
  int i;
  int info;
  int job;
  int pivot[N];
  int seed = 123456789;
  double *x;

  printf ( "\n" );
  printf ( "TEST22\n" );
  printf ( "  For a general banded matrix,\n" );
  printf ( "  R8GB_FA computes the PLU factors.\n" );
  printf ( "  R8GB_SL solves a factored linear system.\n" );
  printf ( "\n" );
  printf ( "  Matrix rows M =    %d\n", M );
  printf ( "  Matrix columns N = %d\n", N );
  printf ( "  Lower bandwidth ML = %d\n", ML );
  printf ( "  Upper bandwidth MU = %d\n", MU );
/*
  Set the matrix.
*/
  a = r8gb_random ( M, N, ML, MU, &seed );

  r8gb_print ( M, N, ML, MU, a, "  The banded matrix:" );
/*
  Set the desired solution.
*/
  x = r8vec_indicator_new ( N );
/*
  Compute the corresponding right hand side.
*/
  b = r8gb_mxv ( M, N, ML, MU, a, x );
/*
  Factor the matrix.
*/
  info = r8gb_fa ( N, ML, MU, a, pivot );

  if ( info != 0 )
  {
    printf ( "\n" );
    printf ( "  Fatal error!\n" );
    printf ( "  R8GB_FA declares the matrix is singular!\n" );
    printf ( "  The value of INFO is %d\n", info );
    return;
  }
/*
  Solve the linear system.
*/
  job = 0;
  free ( x );
  x = r8gb_sl ( N, ML, MU, a, pivot, b, job );
 
  r8vec_print ( N, x, "  Solution:" );
/*
  Set the desired solution.
*/
  free ( x );
  x = r8vec_indicator_new ( N );
/*
  Compute the corresponding right hand side.
*/
  job = 1;
  free ( b );
  b = r8gb_ml ( N, ML, MU, a, pivot, x, job );
  r8vec_print ( N, b, "  Right hand side of transposed system:" );
/*
  Solve the linear system.
*/
  job = 1;
  free ( x );
  x = r8gb_sl ( N, ML, MU, a, pivot, b, job );
 
  r8vec_print ( N, x, "  Solution to transposed system:" );

  free ( a );
  free ( b );
  free ( x );

  return;
# undef M
# undef ML
# undef MU
# undef N
}
/******************************************************************************/

void test23 ( )

/******************************************************************************/
/*
  Purpose:

    TEST23 tests R8GB_FA, R8GB_TRF.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    25 March 2013

  Author:

    John Burkardt
*/
{
# define M 5
# define ML 1
# define MU 1
# define N 5

  double *a;
  int info;
  int pivot[N];
  int seed;

  printf ( "\n" );
  printf ( "TEST23\n" );
  printf ( "  For a general banded matrix,\n" );
  printf ( "  R8GB_FA factors, using LINPACK conventions;\n" );
  printf ( "  R8GB_TRF factors, using LAPACK conventions;\n" );
  printf ( "\n" );
  printf ( "  Matrix rows M =    %d\n", M );
  printf ( "  Matrix columns N = %d\n", N );
  printf ( "  Lower bandwidth ML = %d\n", ML );
  printf ( "  Upper bandwidth MU = %d\n", MU );
/*
  Set the matrix.
*/
  seed = 123456789;
  a = r8gb_random ( M, N, ML, MU, &seed );
/*
  Factor the matrix.
*/
  info = r8gb_fa ( N, ML, MU, a, pivot );

  r8gb_print ( M, N, ML, MU, a, "  The R8GB_FA factors:" );
/*
  Set the matrix.
*/
  free ( a );
  seed = 123456789;
  a = r8gb_random ( M, N, ML, MU, &seed );
/*
  Factor the matrix.
*/
  info = r8gb_trf ( M, N, ML, MU, a, pivot );

  r8gb_print ( M, N, ML, MU, a, "  The R8GB_TRF factors:" );

  free ( a );
  
  return;
# undef M
# undef ML
# undef MU
# undef N
}
/******************************************************************************/

void test235 ( )

/******************************************************************************/
/*
  Purpose:

    TEST235 tests R8GB_INDICATOR.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    26 March 2013

  Author:

    John Burkardt
*/
{
# define M 10
# define ML 3
# define MU 2
# define N 8

  double *a;

  printf ( "\n" );
  printf ( "TEST235\n" );
  printf ( "  For a general banded matrix,\n" );
  printf ( "  R8GB_INDICATOR sets up an indicator matrix.\n" );
  printf ( "\n" );
  printf ( "  Matrix rows M =    %d\n", M );
  printf ( "  Matrix columns N = %d\n", N );
  printf ( "  Lower bandwidth ML = %d\n", ML );
  printf ( "  Upper bandwidth MU = %d\n", MU );

  a = r8gb_indicator ( M, N, ML, MU );

  r8ge_print ( 2*ML+MU+1, N, a, "  The banded matrix in R8GE format:" );

  r8gb_print ( M, N, ML, MU, a, "  The R8GB indicator matrix:" );

  free ( a );

  return;
# undef M
# undef ML
# undef MU
# undef N
}
/******************************************************************************/

void test24 ( )

/******************************************************************************/
/*
  Purpose:

    TEST24 tests R8GB_ML.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    27 March 2013

  Author:

    John Burkardt
*/
{
# define M 10
# define ML 1
# define MU 2
# define N 10

  double *a;
  double *b;
  double *b2;
  int i;
  int info;
  int job;
  int pivot[N];
  int seed = 123456789;
  double *x;

  printf ( "\n" );
  printf ( "TEST24\n" );
  printf ( "  For a general banded matrix,\n" );
  printf ( "  R8GB_ML computes A*x or A'*X\n" );
  printf ( "  where A has been factored by R8GB_FA.\n" );
  printf ( "\n" );
  printf ( "  Matrix rows M =    %d\n", M );
  printf ( "  Matrix columns N = %d\n", N );
  printf ( "  Lower bandwidth ML = %d\n", ML );
  printf ( "  Upper bandwidth MU = %d\n", MU );

  for ( job = 0; job <= 1; job++ )
  {
/*
  Set the matrix.
*/
    a = r8gb_random ( M, N, ML, MU, &seed );
/*
  Set the desired solution.
*/
    x = r8vec_indicator_new ( N );
/*
  Compute the corresponding right hand side.
*/
    if ( job == 0 )
    {
      b = r8gb_mxv ( M, N, ML, MU, a, x );
    }
    else
    {
      b = r8gb_vxm ( M, N, ML, MU, a, x );
    }
/*
  Factor the matrix.
*/
    info = r8gb_fa ( N, ML, MU, a, pivot );

    if ( info != 0 )
    {
      printf ( "\n" );
      printf ( "TEST24 - Fatal error!\n" );
      printf ( "  R8GB_FA declares the matrix is singular!\n" );
      printf ( "  The value of INFO is %d\n", info );
      return;
    }
/*
  Now multiply factored matrix times solution to get right hand side again.
*/
    b2 = r8gb_ml ( N, ML, MU, a, pivot, x, job );

    if ( job == 0 )
    {
      r8vec2_print_some ( N, b, b2, 10, "  A*x and PLU*x" );
    }
    else
    {
      r8vec2_print_some ( N, b, b2, 10, "  A'*x and (PLU)'*x" );
    }

    free ( a );
    free ( b );
    free ( b2 );
    free ( x );
  }

  return;
# undef M
# undef ML
# undef MU
# undef N
}
/******************************************************************************/

void test25 ( )

/******************************************************************************/
/*
  Purpose:

    TEST25 tests R8GB_PRINT.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    28 March 2013

  Author:

    John Burkardt
*/
{
# define M 8
# define ML 1
# define MU 3
# define N 10

  double *a;

  printf ( "\n" );
  printf ( "TEST25\n" );
  printf ( "  For a general banded matrix,\n" );
  printf ( "  R8GB_PRINT prints the matrix.\n" );
  printf ( "\n" );
  printf ( "  Matrix rows M =    %d\n", M );
  printf ( "  Matrix columns N = %d\n", N );
  printf ( "  Lower bandwidth ML = %d\n", ML );
  printf ( "  Upper bandwidth MU = %d\n", MU );

  a = r8gb_indicator ( M, N, ML, MU );

  r8gb_print ( M, N, ML, MU, a, "  The banded matrix:" );

  free ( a );

  return;
# undef M
# undef ML
# undef MU
# undef N
}
/******************************************************************************/

void test26 ( )

/******************************************************************************/
/*
  Purpose:

    TEST26 tests R8GB_NZ_NUM.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    29 March 2013

  Author:

    John Burkardt
*/
{
# define M 10
# define N 10
# define ML 1
# define MU 2

  double *a;
  int diag;
  int i;
  int j;
  int nz_num;
  int seed = 123456789;

  printf ( "\n" );
  printf ( "TEST26\n" );
  printf ( "  For a general banded matrix,\n" );
  printf ( "  R8GB_NZ_NUM counts the nonzero entries.\n" );
  printf ( "\n" );
  printf ( "  Matrix rows M =    %d\n", M );
  printf ( "  Matrix columns N = %d\n", N );
  printf ( "  Lower bandwidth ML = %d\n", ML );
  printf ( "  Upper bandwidth MU = %d\n", MU );
/*
  Set the matrix.
*/
  a = r8gb_random ( M, N, ML, MU, &seed );
/*
  Make some zero entries.
*/
  for ( j = 0; j < N; j++ )
  {
    for ( diag = 0; diag < 2*ML+MU+1; diag++ )
    {
      if ( a[diag+j*(2*ML+MU+1)] < 0.3 )
      {
        a[diag+j*(2*ML+MU+1)] = 0.0;
      }
    }
  }

  r8gb_print ( M, N, ML, MU, a, "  The R8GB matrix:" );

  nz_num = r8gb_nz_num ( M, N, ML, MU, a );

  printf ( "\n" );
  printf ( "  Nonzero entries = %d\n", nz_num );

  free ( a );

  return;
# undef M
# undef N
# undef ML
# undef MU
}
/******************************************************************************/

void test265 ( )

/******************************************************************************/
/*
  Purpose:

    TEST265 tests R8GB_TO_R8GE and R8GE_TO_R8GB.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    30 March 2013

  Author:

    John Burkardt
*/
{
# define M 5
# define ML 2
# define MU 1
# define N 8

  double *a;
  double *b;
  double *c;

  printf ( "\n" );
  printf ( "TEST265\n" );
  printf ( "  R8GB_TO_R8GE copies a R8GB matrix to a R8GE matrix.\n" );
  printf ( "  R8GE_TO_R8GB copies a R8GE matrix to a R8GB matrix.\n" );
  printf ( "\n" );
  printf ( "  Matrix rows M =    %d\n", M );
  printf ( "  Matrix columns N = %d\n", N );
  printf ( "  Lower bandwidth ML = %d\n", ML );
  printf ( "  Upper bandwidth MU = %d\n", MU );

  a = r8gb_indicator ( M, N, ML, MU );

  r8gb_print ( M, N, ML, MU, a, "  The R8GB matrix:" );

  b = r8gb_to_r8ge ( M, N, ML, MU, a );

  r8ge_print ( M, N, b, "  The R8GE matrix:" );

  c = r8ge_to_r8gb ( M, N, ML, MU, b );

  r8gb_print ( M, N, ML, MU, c, "  The recovered R8GB matrix:" );

  free ( a );
  free ( b );
  free ( c );

  return;
# undef M
# undef ML
# undef MU
# undef N
}
/******************************************************************************/

void test2655 ( )

/******************************************************************************/
/*
  Purpose:

    TEST2655 tests R8GB_TO_R8S3.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    31 March 2013

  Author:

    John Burkardt
*/
{
# define M 5
# define ML 2
# define MU 1
# define N 8

  double *a;
  double *b;
  int *col;
  int isym;
  int nz_num;
  int *row;

  printf ( "\n" );
  printf ( "TEST2655\n" );
  printf ( "  R8GB_TO_R8S3 copies a R8GB matrix to a R8S3 matrix.\n" );
  printf ( "\n" );
  printf ( "  Matrix rows M =    %d\n", M );
  printf ( "  Matrix columns N = %d\n", N );
  printf ( "  Lower bandwidth ML = %d\n", ML );
  printf ( "  Upper bandwidth MU = %d\n", MU );

  a = r8gb_indicator ( M, N, ML, MU );

  r8gb_print ( M, N, ML, MU, a, "  The R8GB matrix:" );

  nz_num = r8gb_nz_num ( M, N, ML, MU, a );

  printf ( "  Nonzeros NZ_NUM =    %d\n", nz_num );

  row = ( int * ) malloc ( nz_num * sizeof ( int ) );
  col = ( int * ) malloc ( nz_num * sizeof ( int ) );
  b = ( double * ) malloc ( nz_num * sizeof ( double ) );

  r8gb_to_r8s3 ( M, N, ML, MU, a, nz_num, &isym, row, col, b );

  r8s3_print ( M, N, nz_num, isym, row, col, b, "  The R8S3 matrix:" );

  free ( a );
  free ( b );
  free ( col );
  free ( row );

  return;
# undef M
# undef ML
# undef MU
# undef N
}
/******************************************************************************/

void test27 ( )

/******************************************************************************/
/*
  Purpose:

    TEST27 tests R8GB_TRF, R8GB_TRS.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    31 March 2013

  Author:

    John Burkardt
*/
{
# define M 10
# define N 10
# define ML 1
# define MU 2
# define NRHS 1

  double *a;
  double *b;
  int i;
  int info;
  int job;
  int pivot[N];
  int seed = 123456789;
  double *x;

  printf ( "\n" );
  printf ( "TEST27\n" );
  printf ( "  For a general banded matrix,\n" );
  printf ( "  R8GB_TRF computes the PLU factors.\n" );
  printf ( "  R8GB_TRS solves a factored linear system.\n" );
  printf ( "\n" );
  printf ( "  Matrix rows M =    %d\n", M );
  printf ( "  Matrix columns N = %d\n", N );
  printf ( "  Lower bandwidth ML = %d\n", ML );
  printf ( "  Upper bandwidth MU = %d\n", MU );
/*
  Set the matrix.
*/
  a = r8gb_random ( M, N, ML, MU, &seed );
/*
  Set the desired solution.
*/
  x = r8vec_indicator_new ( N );
/*
  Compute the corresponding right hand side.
*/
  b = r8gb_mxv ( M, N, ML, MU, a, x );
/*
  Factor the matrix.
*/
  info = r8gb_trf ( M, N, ML, MU, a, pivot );

  if ( info != 0 )
  {
    printf ( "\n" );
    printf ( "TEST27 - Fatal error!\n" );
    printf ( "  R8GB_TRF declares the matrix is singular!\n" );
    printf ( "  The value of INFO is %d\n", info );
    return;
  }
/*
  Solve the linear system.
*/
  free ( x );
  x = r8gb_trs ( N, ML, MU, NRHS, 'N', a, pivot, b );

  r8vec_print ( N, x, "  Solution:" );
/*
  Set the desired solution.
*/
  free ( x );
  x = r8vec_indicator_new ( N );
/*
  Compute the corresponding right hand side.
*/
  job = 1;
  free ( b );
  b = r8gb_mu ( N, ML, MU, a, pivot, x, job );
/*
  Solve the linear system.
*/
  free ( x );
  x = r8gb_trs ( N, ML, MU, NRHS, 'T', a, pivot, b );

  r8vec_print ( N, x, "  Solution to transposed system:" );

  free ( a );
  free ( b );
  free ( x );

  return;
# undef M
# undef N
# undef ML
# undef MU
# undef NRHS
}
/******************************************************************************/

void test275 ( )

/******************************************************************************/
/*
  Purpose:

    TEST275 tests R8GD_INDICATOR.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    31 March 2013

  Author:

    John Burkardt
*/
{
# define N 10
# define NDIAG 4

  double *a;
  int offset[NDIAG] = { -2, 0, 1, N-1 };

  printf ( "\n" );
  printf ( "TEST275\n" );
  printf ( "  For a general diagonal matrix:\n" );
  printf ( "  R8GD_INDICATOR sets up an indicator matrix;\n" );
  printf ( "\n" );
  printf ( "  Matrix order N            = %d\n", N );
  printf ( "  Number of diagonals NDIAG = %d\n", NDIAG );

  i4vec_print ( NDIAG, offset, "  The offset vector:" );

  a = r8gd_indicator ( N, NDIAG, offset );

  r8gd_print ( N, NDIAG, offset, a, "  The R8GD indicator matrix:" );

  free ( a );

  return;
# undef N
# undef NDIAG
}
/******************************************************************************/

void test28 ( )

/******************************************************************************/
/*
  Purpose:

    TEST28 tests R8GD_MXV. R8GD_PRINT, R8GD_RANDOM, R8GD_VXM.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    01 April 2013

  Author:

    John Burkardt
*/
{
# define N 10
# define NDIAG 4

  double *a;
  double *b;
  int i;
  int j;
  int offset[NDIAG] = { -2, 0, 1, N-1 };
  int seed = 123456789;
  double *x;

  printf ( "\n" );
  printf ( "TEST28\n" );
  printf ( "  For a general diagonal matrix:\n" );
  printf ( "  R8GD_MXV computes A * x;\n" );
  printf ( "  R8GD_PRINT prints it;\n" );
  printf ( "  R8GD_RANDOM randomly generates one;\n" );
  printf ( "  R8GD_VXM computes A'*x;\n" );
  printf ( "\n" );
  printf ( "  Matrix order N            = %d\n", N );
  printf ( "  Number of diagonals NDIAG = %d\n", NDIAG );

  i4vec_print ( NDIAG, offset, "  The offset vector:" );
/*
  Set the matrix.
*/
  a = r8gd_random ( N, NDIAG, offset, &seed );

  r8ge_print ( N, NDIAG, a, "  The raw matrix: " );

  r8gd_print ( N, NDIAG, offset, a, "  The general diagonal matrix:" );
/*
  Set the desired solution.
*/
  x = r8vec_indicator_new ( N );
/*
  Compute the corresponding right hand side.
*/
  b = r8gd_mxv ( N, NDIAG, offset, a, x );

  r8vec_print ( N, b, "  A * x:" );
/*
  Compute the corresponding right hand side.
*/
  free ( b );
  b = r8gd_vxm ( N, NDIAG, offset, a, x );

  r8vec_print ( N, b, "  A' * x:" );

  free ( a );
  free ( b );
  free ( x );

  return;
# undef N
# undef NDIAG
}
/******************************************************************************/

void test285 ( )

/******************************************************************************/
/*
  Purpose:

    TEST285 tests R8GE_CO.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    02 April 2013

  Author:

    John Burkardt
*/
{
# define N 4

  double a[N*N];
  double *a_inverse;
  double a_inverse_norm_l1;
  double a_lu[N*N];
  double a_norm_l1;
  double cond_l1;
  int i;
  int info;
  int j;
  int pivot[N];
  double rcond;
  double row_sum;
  double x = 2.0;
  double y = 3.0;

  printf ( "\n" );
  printf ( "TEST285\n" );
  printf ( "  For a matrix in general storage,\n" );
  printf ( "  R8GE_CO estimates the condition number.\n" );
  printf ( "\n" );
  printf ( "  Matrix order N = %d\n",N );
/*
  Set the matrix.
*/
  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      if ( i == j )
      {
        a[i+j*N] = x + y;
      }
      else
      {
        a[i+j*N] = y;
      }
    }
  }

  a_norm_l1 = 0.0;
  for ( j = 0; j < N; j++ )
  {
    row_sum = 0.0;
    for ( i = 0; i < N; i++ )
    {
      row_sum = row_sum + r8_abs ( a[i+j*N] );
    }
    a_norm_l1 = r8_max ( a_norm_l1, row_sum );
  }

  for ( j = 0; j < N; j++ )
  {
    for ( i = 0; i < N; i++ )
    {
      a_lu[i+j*N] = a[i+j*N];
    }
  }

  info = r8ge_fa ( N, a_lu, pivot );

  a_inverse = r8ge_inverse ( N, a_lu, pivot );

  a_inverse_norm_l1 = 0.0;
  for ( j = 0; j < N; j++ )
  {
    row_sum = 0.0;
    for ( i = 0; i < N; i++ )
    {
      row_sum = row_sum + r8_abs ( a_inverse[i+j*N] );
    }
    a_inverse_norm_l1 = r8_max ( a_inverse_norm_l1, row_sum );
  }

  cond_l1 = a_norm_l1 * a_inverse_norm_l1;

  printf ( "\n" );
  printf ( "  The L1 condition number is %g\n", cond_l1 );
/*
  Factor the matrix.
*/
  rcond = r8ge_co ( N, a, pivot );

  printf ( "\n" );
  printf ( "  The R8GE_CO estimate is     %g\n", 1.0 / rcond );

  free ( a_inverse );

  return;
# undef N
}
/******************************************************************************/

void test29 ( )

/******************************************************************************/
/*
  Purpose:

    TEST29 tests R8GE_DET.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    03 April 2013

  Author:

    John Burkardt
*/
{
# define N 4

  double a[N*N];
  double det;
  double exact;
  int i;
  int info;
  int j;
  int pivot[N];
  double x = 2.0;
  double y = 3.0;

  printf ( "\n" );
  printf ( "TEST29\n" );
  printf ( "  R8GE_DET, determinant of a general matrix.\n" );
  printf ( "\n" );
  printf ( "  Matrix order N = %d\n", N );
/*
  Set the matrix.
*/
  for ( j = 0; j < N; j++ )
  {
    for ( i = 0; i < N; i++ )
    {
      if ( i == j )
      {
        a[i+j*N] = x + y;
      }
      else
      {
        a[i+j*N] = y;
      }
    }
  }
/*
  Factor the matrix.
*/
  info = r8ge_fa ( N, a, pivot );
/*
  Compute the determinant.
*/
  det = r8ge_det ( N, a, pivot );
  exact = pow ( x, N - 1 ) * ( x + ( ( double ) N ) * y );

  printf ( "\n" );
  printf ( "  R8GE_DET computes the determinant = %g\n", det );
  printf ( "  Correct determinant =               %g\n", exact );

  return;
# undef N
}
/******************************************************************************/

void test295 ( )

/******************************************************************************/
/*
  Purpose:

    TEST295 tests R8GE_DILU.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    04 April 2013

  Author:

    John Burkardt
*/
{
# define NCOL 3
# define NROW 3
# define N NROW * NCOL
# define M N

  double a[M*N];
  double *d;
  int i;
  int j;
  int k;

  printf ( "\n" );
  printf ( "TEST295\n" );
  printf ( "  For a matrix in general storage,\n" );
  printf ( "  R8GE_DILU returns the DILU factors.\n" );
  printf ( "\n" );
  printf ( "  Matrix rows M =    %d\n", M );
  printf ( "  Matrix columns N = %d\n", N );

  for ( i = 1; i <= NROW * NCOL; i++ )
  {
    for ( j = 1; j <= NROW * NCOL; j++ )
    {
      if ( i == j )
      {
        a[i-1+(j-1)*M] = 4.0;
      }
      else if ( 
        i == j + 1 ||
        i == j - 1 ||
        i == j + NROW ||
        i == j - NROW )
      {
        a[i-1+(j-1)*M] = -1.0;
      }
      else
      {
        a[i-1+(j-1)*M] = 0.0;
      }
    }
  }

  r8ge_print ( M, N, a, "  Matrix A:" );
/*
  Compute the incomplete LU factorization.
*/
  d = r8ge_dilu ( M, N, a );

  r8vec_print ( M, d, "  DILU factor of A:" );

  free ( d );

  return;
# undef M
# undef N
# undef NCOL
# undef NROW
}
/******************************************************************************/

void test30 ( )

/******************************************************************************/
/*
  Purpose:

    TEST30 tests R8GE_FA, R8GE_SL;

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    04 April 2013

  Author:

    John Burkardt
*/
{
# define N 5

  double *a;
  double *b;
  int i;
  int info;
  int job;
  int pivot[N];
  int seed = 123456789;
  double *x;

  printf ( "\n" );
  printf ( "TEST30\n" );
  printf ( "  R8GE_FA factors a general linear system,\n" );
  printf ( "  R8GE_SL solves a factored system.\n" );
  printf ( "\n" );
  printf ( "  Matrix order N = %d\n", N );
/*
  Set the matrix.
*/
  a = r8ge_random ( N, N, &seed );

  r8ge_print ( N, N, a, "  Random matrix A:" );
/*
  Set the desired solution.
*/
  x = r8vec_indicator_new ( N );
/*
  Compute the corresponding right hand side.
*/
  b = r8ge_mxv ( N, N, a, x );
/*/
  Factor the matrix.
*/
  info = r8ge_fa ( N, a, pivot );

  if ( info != 0 )
  {
    printf ( "\n" );
    printf ( "  Fatal error!\n" );
    printf ( "  R8GE_FA declares the matrix is singular!\n" );
    printf ( "  The value of INFO is %d\n", info );
    return;
  }
/*
  Solve the linear system.
*/
  job = 0;
  free ( x );
  x = r8ge_sl_new ( N, a, pivot, b, job );

  r8vec_print ( N, x, "  Solution:" );
/*
  Set the desired solution.
*/
  for ( i = 0; i < N; i++ )
  {
    x[i] = 1.0;
  }
/*
  Compute the corresponding right hand side.
*/
  job = 0;
  free ( b );
  b = r8ge_ml ( N, a, pivot, x, job );
/*
  Solve the system
*/
  job = 0;
  free ( x );
  x = r8ge_sl_new ( N, a, pivot, b, job );

  r8vec_print ( N, x, "  Solution:" );
/*
  Set the desired solution.
*/
  free ( x );
  x = r8vec_indicator_new ( N );
/*
  Compute the corresponding right hand side.
*/
  job = 1;
  free ( b );
  b = r8ge_ml ( N, a, pivot, x, job );
/*
  Solve the system
*/
  job = 1;
  free ( x );
  x = r8ge_sl_new ( N, a, pivot, b, job );

  r8vec_print ( N, x, "  Solution of transposed system:" );

  free ( a );
  free ( b );
  free ( x );

  return;
# undef N
}
/******************************************************************************/

void test31 ( )

/******************************************************************************/
/*
  Purpose:

    TEST31 tests R8GE_FA, R8GE_SL.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    05 April 2013

  Author:

    John Burkardt
*/
{
# define N 5

  double *a;
  double *b;
  int i;
  int info;
  int j;
  int job;
  int pivot[N];
  int seed = 123456789;
  double *x;

  printf ( "\n" );
  printf ( "TEST31\n" );
  printf ( "  For a matrix in general storage,\n" );
  printf ( "  R8GE_FA computes the LU factors,\n" );
  printf ( "  R8GE_SL solves a factored system.\n" );
  printf ( "\n" );
  printf ( "  Matrix order N = %d\n", N );
/*
  Set the matrix.
*/
  a = r8ge_random ( N, N, &seed );

  r8ge_print ( N, N, a, "  The matrix:" );
/*
  Set the desired solution.
*/
  x = r8vec_indicator_new ( N );
/*
  Compute the corresponding right hand side.
*/
  b = r8ge_mxv ( N, N, a, x );
/*
  Factor the matrix.
*/
  info = r8ge_fa ( N, a, pivot );
 
  if ( info != 0 )
  {
    printf ( "\n" );
    printf ( "TEST31 - Fatal error!\n" );
    printf ( "  R8GE_FA declares the matrix is singular!\n" );
    printf ( "  The value of INFO is %d\n", info );
    return;
  }
/*
  Display the gory details.
*/
  r8mat_print ( N, N, a, "  The compressed LU factors:" );

  i4vec_print ( N, pivot, "  The pivot vector P:" );
/*
  Solve the linear system.
*/
  job = 0;
  free ( x );
  x = r8ge_sl_new ( N, a, pivot, b, job );

  r8vec_print ( N, x, "  Solution:" );

  free ( a );
  free ( b );
  free ( x );

  return;
# undef N
}
/******************************************************************************/

void test315 ( )

/******************************************************************************/
/*
  Purpose:

    TEST315 tests R8GE_ILU.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    05 April 2013

  Author:

    John Burkardt
*/
{
# define NCOL 3
# define NROW 3
# define N NROW * NCOL
# define M N

  double a[M*N];
  int i;
  int j;
  int k;
  double l[M*M];
  double lu[M*N];
  double u[M*N];

  printf ( "\n" );
  printf ( "TEST315\n" );
  printf ( "  For a matrix in general storage,\n" );
  printf ( "  R8GE_ILU returns the ILU factors.\n" );
  printf ( "\n" );
  printf ( "  Matrix rows M =    %d\n", M );
  printf ( "  Matrix columns N = %d\n", N );

  for ( i = 0; i < NROW * NCOL; i++ )
  {
    for ( j = 0; j < NROW * NCOL; j++ )
    {
      if ( i == j )
      {
        a[i+j*M] = 4.0;
      }
      else if ( 
        i == j + 1 | 
        i == j - 1 | 
        i == j + NROW | 
        i == j - NROW 
      )
      {
        a[i+j*M] = -1.0;
      }
      else
      {
        a[i+j*M] = 0.0;
      }
    }
  }

  r8ge_print ( M, N, a, "  Matrix A:" );
/*
  Compute the incomplete LU factorization.
*/
  r8ge_ilu ( M, N, a, l, u );

  r8ge_print ( M, M, l, "  Factor L:" );

  r8ge_print ( M, N, u, "  Factor U:" );

  for ( j = 0; j < N; j++ )
  {
    for ( i = 0; i < M; i++ )
    {
      lu[i+j*M] = 0.0;
      for ( k = 0; k < M; k++ )
      {
        lu[i+j*M] = lu[i+j*M] + l[i+k*M] * u[k+j*M];
      }
    }
  }

  r8ge_print ( M, N, lu, "  Product L*U:" );

  return;
# undef M
# undef N
# undef NCOL
# undef NROW
}
/******************************************************************************/

void test317 ( )

/******************************************************************************/
/*
  Purpose:

    TEST317 tests R8GE_INDICATOR.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    05 April 2013

  Author:

    John Burkardt
*/
{
# define M 7
# define N 5

  double *a;

  printf ( "\n" );
  printf ( "TEST317\n" );
  printf ( "  For a matrix in general storage,\n" );
  printf ( "  R8GE_INDICATOR sets up an indicator matrix.\n" );
  printf ( "\n" );
  printf ( "  Matrix rows M =    %d\n", M );
  printf ( "  Matrix columns N = %d\n", N );

  a = r8ge_indicator ( M, N );

  r8ge_print ( M, N, a, "  The R8GE indicator matrix:" );

  free ( a );

  return;
# undef M
# undef N
}
/******************************************************************************/

void test32 ( )

/******************************************************************************/
/*
  Purpose:

    TEST32 tests R8GE_NP_FA, R8GE_NP_SL.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    05 April 2013

  Author:

    John Burkardt
*/
{
# define N 10

  double *a;
  double *b;
  int i;
  int info;
  int job;
  int seed = 123456789;
  double *x;

  printf ( "\n" );
  printf ( "TEST32\n" );
  printf ( "  For a matrix in general storage,\n" );
  printf ( "  R8GE_NP_FA computes the LU factors without pivoting,\n" );
  printf ( "  R8GE_NP_SL solves factored systems.\n" );
  printf ( "\n" );
  printf ( "  Matrix order N = %d\n", N );
/*
  Set the matrix.
*/
  a = r8ge_random ( N, N, &seed );
/*
  Set the desired solution.
*/
  x = ( double * ) malloc ( N * sizeof ( double ) );

  for ( i = 0; i < N; i++ )
  {
    x[i] = 1.0;
  }
/*
  Compute the corresponding right hand side.
*/
  b = r8ge_mxv ( N, N, a, x );
/*
  Factor the matrix.
*/
  info = r8ge_np_fa ( N, a );

  if ( info != 0 )
  {
    printf ( "\n" );
    printf ( "TEST32 - Fatal error!\n" );
    printf ( "  R8GE_NP_FA declares the matrix is singular!\n" );
    printf ( "  The value of INFO is %d\n", info );
    return;
  }
/*
  Solve the linear system.
*/
  job = 0;
  free ( x );
  x = r8ge_np_sl ( N, a, b, job );
 
  r8vec_print_some ( N, x, 1, 10, "  Solution:" );
/*
  Set the desired solution.
*/
  free ( x );
  x = r8vec_indicator_new ( N );
/*
  Compute the corresponding right hand side.
*/
  job = 0;
  free ( b );
  b = r8ge_np_ml ( N, a, x, job );
/*
  Solve the system
*/
  job = 0;
  free ( x );
  x = r8ge_np_sl ( N, a, b, job );

  r8vec_print_some ( N, x, 1, 10, "  Solution:" );
/*
  Set the desired solution.
*/
  free ( x );
  x = r8vec_indicator_new ( N );
/*
  Compute the corresponding right hand side.
*/
  job = 1;
  free ( b );
  b = r8ge_np_ml ( N, a, x, job );
/*
  Solve the system
*/
  job = 1;
  free ( x );
  x = r8ge_np_sl ( N, a, b, job );

  r8vec_print ( N, x, "  Solution of transposed system:" );

  free ( a );
  free ( b );
  free ( x );

  return;
# undef N
}
/******************************************************************************/

void test33 ( )

/******************************************************************************/
/*
  Purpose:

    TEST33 tests R8GE_NP_FA, R8GE_NP_INVERSE.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 April 2013

  Author:

    John Burkardt
*/
{
# define N 5

  double *a;
  double *a_lu;
  double *b;
  double *c;
  int i;
  int info;
  int seed = 123456789;
  int j;

  printf ( "\n" );
  printf ( "TEST33\n" );
  printf ( "  For a matrix in general storage,\n" );
  printf ( "  R8GE_NP_FA computes LU factors without pivoting,\n" );
  printf ( "  R8GE_NP_INVERSE computes the inverse.\n" );
  printf ( "\n" );
  printf ( "  Matrix order N = %d\n", N );
/*
  Set the matrix.
*/
  a = r8ge_random ( N, N, &seed );

  r8ge_print ( N, N, a, "  The random matrix:" );
/*
  Factor and invert the matrix.
*/
  a_lu = ( double * ) malloc ( N * N * sizeof ( double ) );

  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      a_lu[i+j*N] = a[i+j*N];
    }
  }

  info = r8ge_np_fa ( N, a_lu );

  if ( info != 0 )
  {
    printf ( "\n" );
    printf ( "TEST33 - Fatal error!\n" );
    printf ( "  R8GE_NP_FA declares the matrix is singular!\n" );
    printf ( "  The value of INFO is %d\n", info );
    return;
  }

  b = r8ge_np_inverse ( N, a_lu );

  r8ge_print ( N, N, b, "  The inverse matrix:" );
/*
  Compute A * B = C.
*/
  c = r8ge_mxm ( N, a, b );

  r8ge_print ( N, N, c, "  The product:" );

  free ( a );
  free ( a_lu );
  free ( b );
  free ( c );

  return;
# undef N
}
/******************************************************************************/

void test34 ( )

/******************************************************************************/
/*
  Purpose:

    TEST34 tests R8GE_FS_NEW.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    07 April 2013

  Author:

    John Burkardt
*/
{
# define N 10

  double *a;
  double *b;
  int i;
  int info;
  int seed = 123456789;
  double *x;

  printf ( "\n" );
  printf ( "TEST34\n" );
  printf ( "  For a matrix in general storage,\n" );
  printf ( "  R8GE_FS_NEW factors and solves a linear system.\n" );
  printf ( "\n" );
  printf ( "  Matrix order N = %d\n", N );
/*
  Set the matrix.
*/
  a = r8ge_random ( N, N, &seed );
/*
  Set the desired solution.
*/
  x = r8vec_indicator_new ( N );
/*
  Compute the corresponding right hand side.
*/
  b = r8ge_mxv ( N, N, a, x );
/*
  Factor and solve the system.
*/
  free ( x );
  x = r8ge_fs_new ( N, a, b );
  
  if ( x == NULL )
  {
    printf ( "\n" );
    printf ( "TEST34 - Fatal error!\n" );
    printf ( "  R8GE_FS_NEW reports the matrix is singular.\n" );
    return;
  }
  r8vec_print ( N, x, "  Solution:" );

  free ( a );
  free ( b );
  free ( x );

  return;
# undef N
}
/******************************************************************************/

void test345 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST345 tests R8GE_FSS_NEW.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    28 November 2011

  Author:

    John Burkardt
*/
{
# define N 10
# define NB 3

  double *a;
  double *b;
  int i;
  int info;
  int j;
  int k;
  int n = N;
  int nb = NB;
  int seed = 123456789;
  double *x;

  printf ( "\n" );
  printf ( "TEST345\n" );
  printf ( "  For a matrix in general storage,\n" );
  printf ( "  R8GE_FSS_NEW factors and solves multiple linear systems.\n" );
  printf ( "\n" );
  printf ( "  Matrix order N = %d\n", n );
  printf ( "  Number of systems NB = %d\n", nb );
/*
  Set the matrix.
*/
  a = r8ge_random ( n, n, &seed );
/*
  Set the desired solutions.
*/
  b = ( double * ) malloc ( n * nb * sizeof ( double ) );

  x = ( double * ) malloc ( n * sizeof ( double ) );

  for ( i = 0; i < n; i++ )
  {
    x[i] = 1.0;
  }
  k = 0;
  for ( i = 0; i < n; i++ )
  {
    b[i+k*n] = 0.0;
    for ( j = 0; j < n; j++ )
    {
      b[i+k*n] = b[i+k*n] + a[i+j*n] * x[j];
    }
  }
  for ( i = 0; i < n; i++ )
  {
    x[i] = ( double ) ( i + 1 );
  }
  k = 1;
  for ( i = 0; i < n; i++ )
  {
    b[i+k*n] = 0.0;
    for ( j = 0; j < n; j++ )
    {
      b[i+k*n] = b[i+k*n] + a[i+j*n] * x[j];
    }
  }
  for ( i = 0; i < n; i++ )
  {
    x[i] = ( i % 3 ) + 1;
  }
  k = 2;
  for ( i = 0; i < n; i++ )
  {
    b[i+k*n] = 0.0;
    for ( j = 0; j < n; j++ )
    {
      b[i+k*n] = b[i+k*n] + a[i+j*n] * x[j];
    }
  }
/*
  Factor and solve the system.
*/
  free ( x );

  x = r8ge_fss_new ( n, a, nb, b );
  
  r8mat_print ( n, nb, x, "  Solutions:" );

  free ( a );
  free ( b );
  free ( x );

  return;
# undef N
# undef NB
}
/******************************************************************************/

void test35 ( )

/******************************************************************************/
/*
  Purpose:

    TEST35 tests R8GE_INVERSE.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    08 April 2013

  Author:

    John Burkardt
*/
{
# define N 4

  double a[N*N];
  double a_lu[N*N];
  double *a_inverse;
  double *c;
  int i;
  int info;
  int j;
  int pivot[N];
  double x = 2.0;
  double y = 3.0;

  printf ( "\n" );
  printf ( "TEST35\n" );
  printf ( "  R8GE_INVERSE inverts a general matrix.\n" );
  printf ( "\n" );
  printf ( "  Matrix order N = %d\n", N );
/*
  Set the matrix.
*/
  for ( j = 0; j < N; j++ )
  {
    for ( i = 0; i < N; i++ )
    {
      if ( i == j )
      {
        a[i+j*N] = x + y;
      }
      else
      {
        a[i+j*N] = y;
      }
    }
  }

  r8ge_print ( N, N, a, "  Matrix A:" );
/*
  Factor and invert the matrix.
*/
  for ( j = 0; j < N; j++ )
  {
    for ( i = 0; i < N; i++ )
    {
      a_lu[i+j*N] = a[i+j*N];
    }
  }

  info = r8ge_fa ( N, a_lu, pivot );

  if ( info != 0 )
  {
    printf ( "\n" );
    printf ( "TEST35 - Fatal error!\n" );
    printf ( "  R8GE_FA reports the matrix is singular.\n" );
    return;
  }
  a_inverse = r8ge_inverse ( N, a_lu, pivot );

  r8ge_print ( N, N, a_inverse, "  Inverse matrix B:" );
/*
  Check.
*/
  c = r8ge_mxm ( N, a, a_inverse );

  r8ge_print ( N, N, c, "  Product matrix:" );

  free ( a_inverse );
  free ( c );

  return;
# undef N
}
/******************************************************************************/

void test36 ( )

/******************************************************************************/
/*
  Purpose:

    TEST36 tests R8GE_ML.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    09 April 2013

  Author:

    John Burkardt
*/
{
# define N 10

  double *a;
  double *b;
  double *b2;
  int info;
  int i;
  int job;
  int pivot[N];
  int seed = 123456789;
  double *x;

  printf ( "\n" );
  printf ( "TEST36\n" );
  printf ( "  R8GE_ML computes A*x or A'*X\n" );
  printf ( "  where A has been factored by R8GE_FA.\n" );
  printf ( "\n" );
  printf ( "  Matrix order N = %d\n", N );

  for ( job = 0; job <= 1; job++ )
  {
/*
  Set the matrix.
*/
    a = r8ge_random ( N, N, &seed );
/*
  Set the desired solution.
*/
    x = r8vec_indicator_new ( N );
/*
  Compute the corresponding right hand side.
*/
    if ( job == 0 )
    {
      b = r8ge_mxv ( N, N, a, x );
    }
    else
    {
      b = r8ge_vxm ( N, N, a, x );
    }
/*
  Factor the matrix.
*/
    info = r8ge_fa ( N, a, pivot );

    if ( info != 0 )
    {
      printf ( "\n" );
      printf ( "  Fatal error!\n" );
      printf ( "  R8GE_FA declares the matrix is singular!\n" );
      printf ( "  The value of INFO is %d\n", info );
      continue;
    }
/*
  Now multiply factored matrix times solution to get right hand side again.
*/
    b2 = r8ge_ml ( N, a, pivot, x, job );

    if ( job == 0 )
    {
      r8vec2_print_some ( N, b, b2, 10, "  A*x and PLU*x\n" );
    }
    else
    {
      r8vec2_print_some ( N, b, b2, 10, "  A'*x and (PLU)'*x\n" );
    }

    free ( a );
    free ( b );
    free ( b2 );
    free ( x );
  }

  return;
# undef N
}
/******************************************************************************/

void test37 ( )

/******************************************************************************/
/*
  Purpose:

    TEST37 tests R8GE_NP_ML.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    09 April 2013

  Author:

    John Burkardt
*/
{
# define N 10

  double *a;
  double *b;
  double *b2;
  int info;
  int i;
  int job;
  int seed = 123456789;
  double *x;

  printf ( "\n" );
  printf ( "TEST37\n" );
  printf ( "  For a matrix in general storage,\n" ); 
  printf ( "  R8GE_NP_ML computes A*x or A'*X\n" );
  printf ( "  where A has been factored by R8GE_NP_FA.\n" );
  printf ( "\n" );
  printf ( "  Matrix order N = %d\n", N );

  for ( job = 0; job <= 1; job++ )
  {
/*
  Set the matrix.
*/
    a = r8ge_random ( N, N, &seed );
/*
  Set the desired solution.
*/
    x = r8vec_indicator_new ( N );
/*
  Compute the corresponding right hand side.
*/
    if ( job == 0 )
    {
      b = r8ge_mxv ( N, N, a, x );
    }
    else
    {
      b = r8ge_vxm ( N, N, a, x );
    }
/*
  Factor the matrix.
*/
    info = r8ge_np_fa ( N, a );

    if ( info != 0 )
    {
      printf ( "\n" );
      printf ( "TEST37 - Fatal error!\n" );
      printf ( "  R8GE_NP_FA declares the matrix is singular!\n" );
      printf ( "  The value of INFO is %d\n", info );
      continue;
    }
/*
  Now multiply factored matrix times solution to get right hand side again.
*/
    b2 = r8ge_np_ml ( N, a, x, job );

    if ( job == 0 )
    {
      r8vec2_print_some ( N, b, b2, 10, "  A*x and PLU*x" );
    }
    else
    {
      r8vec2_print_some ( N, b, b2, 10, "  A'*x and (PLU)'*x" );
    }

    free ( a );
    free ( b );
    free ( b2 );
    free ( x );
  }

  return;
# undef N
}
/******************************************************************************/

void test38 ( )

/******************************************************************************/
/*
  Purpose:

    TEST38 tests R8GE_MU.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    10 April 2013

  Author:

    John Burkardt
*/
{
# define M 5
# define N 3

  double *amn;
  double *anm;
  double *bm;
  double *bn;
  double *cm;
  double *cn;
  int info;
  int i;
  int job;
  int pivot[M+N];
  int seed = 123456789;
  char trans;
  double *xm;
  double *xn;

  printf ( "\n" );
  printf ( "TEST38\n" );
  printf ( "  For a matrix in general storage,\n" );
  printf ( "  R8GE_MU computes A*x or A'*X\n" );
  printf ( "  where A has been factored by R8GE_TRF.\n" );
/*
  First test.
  A is 5 x 3, and we compute A * x.
*/
  printf ( "\n" );
  printf ( "  Matrix rows M =    %d\n", M );
  printf ( "  Matrix columns N = %d\n", N );

  job = 0;
  trans = 'N';

  amn = r8ge_random ( M, N, &seed );
  xn = r8vec_indicator_new ( N );
  cm = r8ge_mxv ( M, N, amn, xn );

  info = r8ge_trf ( M, N, amn, pivot );

  if ( info != 0 )
  {
    printf ( "\n" );
    printf ( "TEST38 - Fatal error!\n" );
    printf ( "  R8GE_TRF declares the matrix is singular!\n" );
    printf ( "  The value of INFO is %d\n", info );
    return;
  }

  bm = r8ge_mu ( M, N, amn, trans, pivot, xn );

  r8vec2_print_some ( M, cm, bm, 10, "  A*x and PLU*x" );

  free ( amn );
  free ( bm );
  free ( cm );
  free ( xn );
/*
  Second test.
  A is 5 x 3, and we compute A' * x.
*/
  printf ( "\n" );
  printf ( "  Matrix rows M =    %d\n", M );
  printf ( "  Matrix columns N = %d\n", N );

  job = 1;
  trans = 'T';

  amn = r8ge_random ( M, N, &seed );

  xm = r8vec_indicator_new ( M );

  cn = r8ge_vxm ( M, N, amn, xm );

  info = r8ge_trf ( M, N, amn, pivot );

  if ( info != 0 )
  {
    printf ( "\n" );
    printf ( "TEST38 - Fatal error!\n" );
    printf ( "  R8GE_TRF declares the matrix is singular!\n" );
    printf ( "  The value of INFO is %d\n", info );
    return;
  }

  bn = r8ge_mu ( M, N, amn, trans, pivot, xm );

  r8vec2_print_some ( N, cn, bn, 10, "  A'*x and (PLU)'*x" );

  free ( amn );
  free ( bn );
  free ( cn );
  free ( xm );
/*
  Third test.
  A is 3 x 5, and we compute A * x.
*/
  printf ( "\n" );
  printf ( "  Matrix rows M =    %d\n", N );
  printf ( "  Matrix columns N = %d\n", M );

  job = 0;
  trans = 'N';

  anm = r8ge_random ( N, M, &seed );

  xm = r8vec_indicator_new ( M );
  cn = r8ge_mxv ( N, M, anm, xm );
 
  info = r8ge_trf ( N, M, anm, pivot );

  if ( info != 0 )
  {
    printf ( "\n" );
    printf ( "TEST38 - Fatal error!\n" );
    printf ( "  R8GE_TRF declares the matrix is singular!\n" );
    printf ( "  The value of INFO is %d\n", info );
    return;
  }

  bn = r8ge_mu ( N, M, anm, trans, pivot, xm );

  r8vec2_print_some ( N, cn, bn, 10, "  A*x and PLU*x" );

  free ( anm );
  free ( bn );
  free ( cn );
  free ( xm );
/*
  Fourth test.
  A is 3 x 5, and we compute A' * x.
*/
  printf ( "\n" );
  printf ( "  Matrix rows M =    %d\n", N );
  printf ( "  Matrix columns N = %d\n", M );

  job = 1;
  trans = 'T';

  anm = r8ge_random ( N, M, &seed );

  xn = r8vec_indicator_new ( N );
  cm = r8ge_vxm ( N, M, anm, xn );

  info = r8ge_trf ( N, M, anm, pivot );

  if ( info != 0 )
  {
    printf ( "\n" );
    printf ( "TEST38 - Fatal error!\n" );
    printf ( "  R8GE_TRF declares the matrix is singular!\n" );
    printf ( "  The value of INFO is %d\n", info );
    return;
  }

  bm = r8ge_mu ( N, M, anm, trans, pivot, xn );

  r8vec2_print_some ( M, cm, bm, 10, "  A'*x and (PLU)'*x" );

  free ( anm );
  free ( bm );
  free ( cm );
  free ( xn );

  return;
# undef M
# undef N
}
/******************************************************************************/

void test385 ( )

/******************************************************************************/
//
//  Purpose:
/*
    TEST385 tests R8GE_PLU.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    10 April 2013

  Author:

    John Burkardt
*/
{
# define M 5
# define N 4

  double *a;
  int i;
  int j;
  int k;
  double l[M*M];
  double lu[M*N];
  double p[M*M];
  double plu[M*N];
  int seed;
  double u[M*N];

  printf ( "\n" );
  printf ( "TEST385\n" );
  printf ( "  For a matrix in general storage,\n" );
  printf ( "  R8GE_PLU returns the PLU factors of a matrix.\n" );
  printf ( "\n" );
  printf ( "  Matrix rows M =    %d\n", M );
  printf ( "  Matrix columns N = %d\n", N );

  seed = 123456789;

  a = r8ge_random ( M, N, &seed );

  r8ge_print ( M, N, a, "  Matrix A:" );
/*
  Compute the PLU factors.
*/
  r8ge_plu ( M, N, a, p, l, u );

  r8ge_print ( M, M, p, "  Factor P:" );

  r8ge_print ( M, M, l, "  Factor L:" );

  r8ge_print ( M, N, u, "  Factor U:" );

  for ( i = 0; i < M; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      lu[i+j*M] = 0.0;
      for ( k = 0; k < M; k++ )
      {
        lu[i+j*M] = lu[i+j*M] + l[i+k*M] * u[k+j*M];
      }
    }
  }

  for ( i = 0; i < M; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      plu[i+j*M] = 0.0;
      for ( k = 0; k < M; k++ )
      {
        plu[i+j*M] = plu[i+j*M] + p[i+k*M] * lu[k+j*M];
      }
    }
  }
  r8ge_print ( M, N, plu, "  Product P*L*U:" );

  free ( a );

  return;
# undef M
# undef N
}
/******************************************************************************/

void test39 ( )

/******************************************************************************/
/*
  Purpose:

    TEST39 tests R8GE_POLY.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    10 April 2013

  Author:

    John Burkardt
*/
{
# define N 12

  double a[N*N];
  int i;
  int j;
  double *p;
  double p_true[N+1] = {   
         1.0,    -23.0,    231.0,  -1330.0,    4845.0, 
    -11628.0,  18564.0, -19448.0,  12870.0,   -5005.0, 
      1001.0,    -78.0,      1.0 };

  printf ( "\n" );
  printf ( "TEST39\n" );
  printf ( "  For a matrix in general storage,\n" );
  printf ( "  R8GE_POLY computes the characteristic polynomial.\n" );
  printf ( "\n" );
  printf ( "  Matrix order N = %d\n", N );
/*
  Set the matrix.
*/
  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      a[i+j*N] = ( double ) i4_min ( i + 1, j + 1 );
    }
  }
/*
  Get the characteristic polynomial.
*/
  p = r8ge_poly ( N, a );
/*
  Compare.
*/
  r8vec2_print_some ( N+1, p, p_true, 10, "I, P(I), True P(I)" );

  free ( p );

  return;
# undef N
}
/******************************************************************************/

void test40 ( )

/******************************************************************************/
/*
  Purpose:

    TEST40 tests R8GE_SL_IT.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 April 2013

  Author:

    John Burkardt
*/
{
# define N 6

  double *a;
  double a_lu[N*N];
  double b[N];
  int i;
  int info;
  int j;
  int job;
  int pivot[N];
  double *r;
  double *x;
  double *x_new;

  printf ( "\n" );
  printf ( "TEST40\n" );
  printf ( "  For a matrix in general storage,\n" );
  printf ( "  R8GE_SL_IT applies one step of iterative \n" ); 
  printf ( "  refinement to an R8GE_SL solution.\n" );
  printf ( "\n" );
  printf ( "  Matrix order N = %d\n", N );
/*
  Set the coefficient matrix.
*/
  a = hilbert_inverse ( N );
/*
  Set the right hand side b.
*/
  for ( i = 0; i < N-1; i++ )
  {
    b[i] = 0.0;
  }
  b[N-1] = 1.0;
/*
  It is necessary to keep both an unfactored and factored copy
  of the coefficient matrix.
*/
  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      a_lu[i+j*N] = a[i+j*N];
    }
  }
/*
  Compute the factored coefficient matrix.
*/
  info = r8ge_fa ( N, a_lu, pivot );

  if ( info != 0 )
  {
    printf ( "\n" );
    printf ( "TEST40 - Fatal error!\n" );
    printf ( "  R8GE_FA declares the matrix is singular!\n" );
    printf ( "  The value of INFO is %d\n", info );
    return;
  }
/*
  Solve the system.
  (Careful!  R8GE_SL overwrites the right hand side with the solution!)
*/
  job = 0;
  x = r8ge_sl_new ( N, a_lu, pivot, b, job );
/*
  Compute and print the residual.
*/
  r = r8ge_res ( N, N, a, x, b );

  r8vec2_print_some ( N, x, r, 10, "  i, x, b-A*x" );
/*
  Take a few steps of iterative refinement.
*/
  for ( j = 1; j <= 5; j++ )
  {
    printf ( "\n" );
    printf ( "Iterative refinement step %d\n", j );
    printf ( "\n" );
/*
  Improve the solution.
*/
    job = 0;
    x_new = r8ge_sl_it ( N, a, a_lu, pivot, b, job, x );
/*
  Compute and print the residual.
*/
    free ( r );

    r = r8ge_res ( N, N, a, x_new, b );

    r8vec2_print_some ( N, x_new, r, 10, "  i, x, b-A*x" );

    for ( i = 0; i < N; i++ )
    {
      x[i] = x_new[i];
    }
    free ( x_new );
  }

  free ( a );
  free ( r );
  free ( x );

  return;
# undef N
}
/******************************************************************************/

void test41 ( )

/******************************************************************************/
/*
  Purpose:

    TEST41 tests R8GE_TRF, R8GE_TRS.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 September 2006

  Author:

    John Burkardt
*/
{
# define N 5
# define M 5
# define NRHS 1

  double a[N*N];
  double b[N*NRHS];
  int i;
  int info;
  int j;
  int pivot[N];
  double *x;

  printf ( "\n" );
  printf ( "TEST41\n" );
  printf ( "  For a matrix in general storage,\n" );
  printf ( "  R8GE_TRF computes the LU factors,\n" );
  printf ( "  R8GE_TRS solves a factored system.\n" );
  printf ( "\n" );
  printf ( "  Matrix rows M =    %d\n", M );
  printf ( "  Matrix columns N = %d\n", N );

  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      if ( i == j )
      {
        a[i+j*N] = 2.0;
      }
      else if ( i == j - 1 )
      {
        a[i+j*N] = - 1.0;
      }
      else if ( i == j + 1 )
      {
        a[i+j*N] = - 1.0;
      }
      else
      {
        a[i+j*N] = 0.0;
      }
    }
  }

  info = r8ge_trf ( M, N, a, pivot );

  if ( info != 0 )
  {
    printf ( "\n" );
    printf ( "TEST41 - Fatal error!\n" );
    printf ( "  R8GE_TRF declares the matrix is singular!\n" );
    printf ( "  The value of INFO is %d\n", info );
    return;
  }

  for ( i = 0; i < N-1; i++ )
  {
    b[i+0*N] = 0.0;
  }
  b[N-1+0*N] = ( double ) ( N + 1 );

  x = r8ge_trs ( N, NRHS, 'N', a, pivot, b );

  if ( x == NULL )
  {
    printf ( "\n" );
    printf ( "TEST41 - Fatal error!\n" );
    printf ( "  R8GE_TRS returned an error condition!\n" );
    return;
  }
  r8vec_print ( N, x, "  Solution:" );

  for ( i = 0; i < N-1; i++ )
  {
    b[i+0*N] = 0.0;
  }
  b[N-1+0*N] = ( double ) ( N + 1 );

  free ( x );

  x = r8ge_trs ( N, NRHS, 'T', a, pivot, b );

  if ( x == NULL )
  {
    printf ( "\n" );
    printf ( "TEST41 - Fatal error!\n" );
    printf ( "  R8GE_TRS returned an error condition!\n" );
    return;
  }

  r8vec_print ( N, x, "  Solution to transposed system:" );

  free ( x );

  return;
# undef M
# undef N
# undef NRHS
}
/******************************************************************************/

void test42 ( )

/******************************************************************************/
/*
  Purpose:

    TEST42 tests R8GE_NP_TRF, R8GE_NP_TRM, R8GE_NP_TRS.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 April 2013

  Author:

    John Burkardt
*/
{
# define M 10
# define N 10
# define NRHS 1

  double *a;
  double *b;
  int i;
  int info;
  int job;
  int seed = 123456789;
  double *x;

  printf ( "\n" );
  printf ( "TEST42\n" );
  printf ( "  For a matrix in general storage,\n" );
  printf ( "  R8GE_NP_TRF factors without pivoting,\n" );
  printf ( "  R8GE_NP_TRS solves factored systems.\n" );
  printf ( "  R8GE_NP_TRM computes A*X for factored A.\n" );
  printf ( "\n" );
  printf ( "  Matrix rows M =    %d\n", M );
  printf ( "  Matrix columns N = %d\n", N );
/*
  Set the matrix.
*/
  a = r8ge_random ( M, N, &seed );
/*
  Set the desired solution.
*/
  x = ( double * ) malloc ( N * sizeof ( double ) );
  for ( i = 0; i < N; i++ )
  {
    x[i] = 1.0;
  }
/*
  Compute the corresponding right hand side.
*/
  b = r8ge_mxv ( M, N, a, x );
/*
  Factor the matrix.
*/
  info = r8ge_np_trf ( M, N, a );

  if ( info != 0 )
  {
    printf ( "\n" );
    printf ( "TEST42 - Fatal error!\n" );
    printf ( "  R8GE_NP_TRF declares the matrix is singular!\n" );
    printf ( "  The value of INFO is %d\n", info );
    return;
  }
/*
  Solve the linear system.
*/
  free ( x );
  x = r8ge_np_trs ( N, NRHS, 'N', a, b );
 
  if ( info != 0 )
  {
    printf ( "\n" );
    printf ( "TEST42 - Fatal error!\n" );
    printf ( "  R8GE_TRS returned an error condition!\n" );
    printf ( "  The value of INFO is %d\n", info );
    return;
  }
  r8vec_print ( N, x, "  Solution:" );
/*
  Set the desired solution.
*/
  free ( x );
  x = r8vec_indicator_new ( N );
/*
  Compute the corresponding right hand side.
*/
  job = 0;
  free ( b );
  b = r8ge_np_trm ( M, N, a, x, job );
/*
  Solve the system
*/
  x = r8ge_np_trs ( N, NRHS, 'N', a, b );

  if ( info != 0 )
  {
    printf ( "\n" );
    printf ( "TEST42 - Fatal error!\n" );
    printf ( "  R8GE_TRS returned an error condition!\n" );
    printf ( "  The value of INFO is %d\n", info );
    return;
  }

  r8vec_print ( N, x, "  Solution:" );
/*
  Set the desired solution.
*/
  free ( x );
  x = r8vec_indicator_new ( N );
/*
  Compute the corresponding right hand side.
*/
  job = 1;
  free ( b );
  b = r8ge_np_trm ( M, N, a, x, job );
/*
  Solve the system.
*/
  free ( x );
  x = r8ge_np_trs ( N, NRHS, 'T', a, b );

  if ( info != 0 )
  {
    printf ( "\n" );
    printf ( "TEST42 - Fatal error!\n" );
    printf ( "  R8GE_TRS returned an error condition!\n" );
    printf ( "  The value of INFO is %d\n", info );
    return;
  }

  r8vec_print ( N, x, "  Solution of transposed system:" );

  free ( a );
  free ( b );
  free ( x );

  return;
# undef M
# undef N
# undef NRHS
}
/******************************************************************************/

void test422 ( )

/******************************************************************************/
/*
  Purpose:

    TEST422 tests R8CC_GET, R8CC_IJK, R8CC_INC, R8CC_KIJ, R8CC_SET.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 April 2013

  Author:

    John Burkardt
*/
{
# define M 5
# define N 5
# define NZ_NUM 12

  double *a;
  int colptr[N+1] = { 1, 4, 6, 8, 10, 13 };
  int i;
  int j;
  int k;
  int rowind[NZ_NUM] = {
    1, 2, 4, 1, 2, 3, 5, 4, 5, 1, 2, 5 };
  int seed = 123456789;
  int test;
  int test_num = 20;
  double value;

  printf ( "\n" );
  printf ( "TEST422\n" );
  printf ( "  For a matrix in the R8CC format,\n" );
  printf ( "  (double precision Harwell-Boeing Unsymmetric Assembled)\n" );
  printf ( "  R8CC_GET gets an entry;\n" );
  printf ( "  R8CC_IJK gets K from (I,J)\n" );
  printf ( "  R8CC_INC increments an entry;\n" );
  printf ( "  R8CC_KIJ gets (I,J) from K;\n" );
  printf ( "  R8CC_SET sets an entry;\n" );
  printf ( "\n" );
  printf ( "  Matrix rows M    = %d\n", M );
  printf ( "  Matrix columns N = %d\n", N );
  printf ( "  Nonzeros NZ_NUM  = %d\n", NZ_NUM );

  i4vec_print ( N+1, colptr, "  The COLPTR vector:" );

  i4vec_print ( NZ_NUM, rowind, "  The ROWIND vector:" );
/*
  Initialize the matrix to random values.
*/
  a = r8cc_random ( M, N, NZ_NUM, colptr, rowind, &seed );

  r8cc_print ( M, N, NZ_NUM, colptr, rowind, a, 
    "  The initial R8CC matrix:" );

  printf ( "\n" );
  printf ( "  R8CC_IJK locates some (I,J) entries.\n" );
  printf ( "\n" );
  printf ( "         I         J         K\n" );
  printf ( "\n" );

  for ( test = 1; test <= test_num; test++ )
  {
    i = i4_uniform ( 1, M, &seed );
    j = i4_uniform ( 1, N, &seed );
    k = r8cc_ijk ( M, N, NZ_NUM, colptr, rowind, i, j );
    printf ( "  %8d  %8d  %8d\n", i, j, k );
  }

  printf ( "\n" );
  printf ( "  R8CC_KIJ locates some K entries.\n" );
  printf ( "\n" );
  printf ( "         K         I         J\n" );
  printf ( "\n" );

  for ( test = 1; test <= test_num; test++ )
  {
    k = i4_uniform ( 1, NZ_NUM, &seed );
    r8cc_kij ( M, N, NZ_NUM, colptr, rowind, k, &i, &j );
    printf ( "  %8d  %8d  %8d\n", k, i, j );
  }

  printf ( "\n" );
  printf ( "  R8CC_SET sets 10 entries at random.\n" );
  printf ( "\n" );
  printf ( "         I         J         K      NEW_VALUE\n" );
  printf ( "\n" );

  for ( test = 1; test <= 10; test++ )
  {
    k = i4_uniform ( 1, NZ_NUM, &seed );
    r8cc_kij ( M, N, NZ_NUM, colptr, rowind, k, &i, &j );
    value = 100.0 + ( double ) ( test );
    r8cc_set ( M, N, NZ_NUM, colptr, rowind, a, i, j, value );
    printf ( "  %8d  %8d  %8d  %14g\n", i, j, k, value );
  }

  printf ( "\n" );
  printf ( "  R8CC_INC increments 10 entries at random.\n" );
  printf ( "\n" );
  printf ( "         I         J         K       NEW_VALUE\n" );
  printf ( "\n" );

  for ( test = 1; test <= 10; test++ )
  {
    k = i4_uniform ( 1, NZ_NUM, &seed );
    r8cc_kij ( M, N, NZ_NUM, colptr, rowind, k, &i, &j );
    value = 20.0 + ( double ) ( test );
    r8cc_inc ( M, N, NZ_NUM, colptr, rowind, a, i, j, value );
    value = r8cc_get ( M, N, NZ_NUM, colptr, rowind, a, i, j );
    printf ( "  %8d  %8d  %8d  %14g\n", i, j, k, value );
  }

  printf ( "\n" );
  printf ( "  R8CC_GET retrieves 10 entries.\n" );
  printf ( "\n" );
  printf ( "         I         J         K      VALUE\n" );
  printf ( "\n" );

  for ( test = 1; test <= 10; test++ )
  {
    k = i4_uniform ( 1, NZ_NUM, &seed );
    r8cc_kij ( M, N, NZ_NUM, colptr, rowind, k, &i, &j );
    value = r8cc_get ( M, N, NZ_NUM, colptr, rowind, a, i, j );
    printf ( "  %8d  %8d  %8d  %14g\n", i, j, k, value );
  }

  r8cc_print ( M, N, NZ_NUM, colptr, rowind, a, "  The final R8CC matrix:" );

  free ( a );

  return;
# undef M
# undef N
# undef NZ_NUM
}
/******************************************************************************/

void test423 ( )

/******************************************************************************/
/*
  Purpose:

    TEST423 tests R8CC_INDICATOR.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 April 2013

  Author:

    John Burkardt
*/
{
# define M 5
# define N 5
# define NZ_NUM 12

  double *a;
  int colptr[N+1] = { 1, 4, 6, 8, 10, 13 };
  int rowind[NZ_NUM] = { 1, 2, 4, 1, 2, 3, 5, 4, 5, 1, 2, 5 };

  printf ( "\n" );
  printf ( "TEST423\n" );
  printf ( "  R8CC_INDICATOR sets up an SHBUA indicator matrix;\n" );
  printf ( "\n" );
  printf ( "  Matrix rows M    = %d\n", M );
  printf ( "  Matrix columns N = %d\n", N );
  printf ( "  Nonzeros NZ_NUM  = %d\n", NZ_NUM );

  a = r8cc_indicator ( M, N, NZ_NUM, colptr, rowind );

  r8cc_print ( M, N, NZ_NUM, colptr, rowind, a, "  The R8CC indicator matrix:" );

  free ( a );

  return;
# undef M
# undef N
# undef NZ_NUM
}
/******************************************************************************/

void test425 ( )

/******************************************************************************/
/*
  Purpose:

    TEST425 tests R8CC_MXV, R8CC_VXM.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 April 2013

  Author:

    John Burkardt
*/
{
# define M 5
# define N 5
# define NZ_NUM 12

  double *a;
  double *bm;
  double *bn;
  double *c;
  int colptr[N+1] = { 1, 4, 6, 8, 10, 13 };
  int i;
  int rowind[NZ_NUM] = { 1, 2, 4, 1, 2, 3, 5, 4, 5, 1, 2, 5 };
  int seed = 123456789;
  double xn[N];
  double xm[M];

  printf ( "\n" );
  printf ( "TEST425\n" );
  printf ( "  R8CC_MXV multiplies an SHBUA matrix by a vector;\n" );
  printf ( "  R8CC_VXM multiplies a vector by an SHBUA matrix;\n" );
  printf ( "\n" );
  printf ( "  Matrix rows M    = %d\n", M );
  printf ( "  Matrix columns N = %d\n", N );
  printf ( "  Nonzeros NZ_NUM  = %d\n", NZ_NUM );
/*
  Set the matrix.
*/
  a = r8cc_random ( M, N, NZ_NUM, colptr, rowind, &seed );
/*
  Make an R8GE copy.
*/
  c = r8cc_to_r8ge ( M, N, NZ_NUM, colptr, rowind, a );
/*
  Print the R8GE copy.
*/
  r8ge_print ( N, N, c, "  The R8CC matrix, in R8GE form:" );
/*
  Compute A * xn = bm.
*/
  xn[0] = 1.0;
  for ( i = 1; i < N-1; i++ )
  {
    xn[i] = 0.0;
  }
  xn[N-1] = -1.0;

  r8vec_print ( N, xn, "  The vector xn:" );

  bm = r8cc_mxv ( M, N, NZ_NUM, colptr, rowind, a, xn );

  r8vec_print ( M, bm, "  The product A * xn:" );
/*
  Compute xm * A = bn.
*/
  xm[0] = 1.0;
  for ( i = 1; i < M-1; i++ )
  {
    xm[i] = 0.0;
  }
  xm[M-1] = -1.0;

  bn = r8cc_vxm ( M, N, NZ_NUM, colptr, rowind, a, xm );

  r8vec_print ( N, bn, "  The product A' * xm:" );

  free ( a );
  free ( bm );
  free ( bn );
  free ( c );

  return;
# undef M
# undef N
# undef NZ_NUM
}
/******************************************************************************/

void test426 ( )

/******************************************************************************/
/*
  Purpose:

    TEST426 tests R8CC_PRINT.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 April 2013

  Author:

    John Burkardt
*/
{
# define M 5
# define N 5
# define NZ_NUM 12

  double *a;
  int colptr[N+1] = { 1, 4, 6, 8, 10, 13 };
  int rowind[NZ_NUM] = { 1, 2, 4, 1, 2, 3, 5, 4, 5, 1, 2, 5 };
  int seed = 123456789;

  printf ( "\n" );
  printf ( "TEST426\n" );
  printf ( "  R8CC_PRINT prints an R8CC matrix.\n" );
  printf ( "\n" );
  printf ( "  Matrix rows M    = %d\n", M );
  printf ( "  Matrix columns N = %d\n", N );
  printf ( "  Nonzeros NZ_NUM  = %d\n", NZ_NUM );
/*
  Set the matrix.
*/
  a = r8cc_random ( M, N, NZ_NUM, colptr, rowind, &seed );
/*
  Print the matrix.
*/
  r8cc_print ( M, N, NZ_NUM, colptr, rowind, a, "  The R8CC matrix:" );

  free ( a );

  return;
# undef M
# undef N
# undef NZ_NUM
}
/******************************************************************************/

void test428 ( )

/******************************************************************************/
/*
  Purpose:

    TEST428 tests R8LT_INDICATOR;

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    13 April 2013

  Author:

    John Burkardt
*/
{
# define M 6
# define N 5

  double *a;

  printf ( "\n" );
  printf ( "TEST428\n" );
  printf ( "  For a matrix in lower triangular storage,\n" );
  printf ( "  R8LT_INDICATOR sets up an indicator matrix;\n" );
  printf ( "\n" );
  printf ( "  Matrix rows M =    %d\n", M );
  printf ( "  Matrix columns N = %d\n", N );

  a = r8lt_indicator ( M, N );

  r8lt_print ( M, N, a, "  The R8LT indicator matrix:" );

  free ( a );

  return;
# undef M
# undef N
}
/******************************************************************************/

void test43 ( )

/******************************************************************************/
/*
  Purpose:

    TEST43 tests R8LT_SL, R8LT_MXV, R8LT_VXM.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    13 April 2013

  Author:

    John Burkardt
*/
{
# define N 10

  double a[N*N];
  double *b;
  int i;
  int j;
  int job;
  double *x;

  printf ( "\n" );
  printf ( "TEST43\n" );
  printf ( "  For a matrix in lower triangular storage,\n" );
  printf ( "  R8LT_SL solves systems;\n" );
  printf ( "  R8LT_MXV computes matrix-vector products;\n" );
  printf ( "  R8LT_VXM computes vector-matrix products;\n" );
  printf ( "\n" );
  printf ( "  Matrix order N = %d\n", N );

  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      if ( j <= i )
      {
        a[i+j*N] = ( double ) ( j + 1 );
      }
      else
      {
        a[i+j*N] = 0.0;
      }
    }
  }

  r8lt_print ( N, N, a, "  The lower triangular matrix:" );

  for ( job = 0; job <= 1; job++ )
  {
/*
  Set the desired solution.
*/
    x = r8vec_indicator_new ( N );
/*
  Compute the corresponding right hand side.
*/
    if ( job == 0 )
    {
      b = r8lt_mxv ( N, N, a, x );
    }
    else
    {
      b = r8lt_vxm ( N, N, a, x );
    }
/*
  Solve the linear system.
*/
    free ( x );
    x = r8lt_sl ( N, a, b, job );
 
    if ( job == 0 )
    {
      r8vec_print ( N, x, "  Solution:" );
    }
    else
    {
      r8vec_print ( N, x, "  Solution to the transposed system:" );
    }
    free ( b );
    free ( x );
  }

  return;
# undef N
}
/******************************************************************************/

void test44 ( )

/******************************************************************************/
/*
  Purpose:

    TEST44 tests R8LT_DET, R8LT_INVERSE, R8LT_MXM, R8LT_RANDOM.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    13 April 2013

  Author:

    John Burkardt
*/
{
# define N 5

  double *a;
  double *b;
  double *c;
  double det;
  int i;
  int j;
  int seed = 123456789;

  printf ( "\n" );
  printf ( "TEST44\n" );
  printf ( "  For a matrix in lower triangular storage,\n" );
  printf ( "  R8LT_DET computes the determinant.\n" );
  printf ( "  R8LT_INVERSE computes the inverse.\n" );
  printf ( "  R8LT_MXM computes matrix products.\n" );
  printf ( "  R8LT_RANDOM sets a random value.\n" );
  printf ( "\n" );
  printf ( "  Matrix order N = %d\n", N );

  a = r8lt_random ( N, N, &seed );

  r8lt_print ( N, N, a, "  Matrix A:" );
/*
  Compute the determinant.
*/
  det = r8lt_det ( N, a );

  printf ( "\n" );
  printf ( "  Determinant is %g\n", det );
/*
  Compute the inverse matrix.
*/
  b = r8lt_inverse ( N, a );

  r8lt_print ( N, N, b, "  Inverse matrix B:" );
/*
  Check.
*/
  c = r8lt_mxm ( N, a, b );

  r8lt_print ( N, N, c, "  Product A * B:" );

  free ( a );
  free ( b );
  free ( c );

  return;
# undef N
}
/******************************************************************************/

void test443 ( )

/******************************************************************************/
/*
  Purpose:

    TEST443 tests R8NCF_INDICATOR.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    14 April 2013

  Author:

   John Burkardt
*/
{
# define M 7
# define N 5
# define NZ_NUM 15

  double *a;
  int rowcol[2*NZ_NUM] = {
    1, 1, 
    2, 2,
    3, 3,
    4, 4, 
    5, 5, 
    2, 1, 
    5, 1, 
    1, 2, 
    5, 2,
    1, 4, 
    2, 4, 
    3, 4,
    4, 5, 
    4, 6, 
    1, 7 };

  printf ( "\n" );
  printf ( "TEST443\n" );
  printf ( "  R8NCF_INDICATOR sets up a R8NCF indicator matrix;\n" );
  printf ( "\n" );
  printf ( "  Matrix rows M =    %d\n", M );
  printf ( "  Matrix columns N = %d\n", N );
  printf ( "  Matrix nonzeros =  %d\n", NZ_NUM );

  a = r8ncf_indicator ( M, N, NZ_NUM, rowcol );

  r8ncf_print ( M, N, NZ_NUM, rowcol, a, "  The R8NCF indicator matrix:" );

  free ( a );

  return;
# undef M
# undef N
# undef NZ_NUM
}
/******************************************************************************/

void test445 ( )

/******************************************************************************/
/*
  Purpose:

    TEST445 tests R8PBL_INDICATOR.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    14 April 2013

  Author:

    John Burkardt
*/
{
# define N 9
# define MU 3

  double *a;
  int i;
  int j;

  printf ( "\n" );
  printf ( "TEST445\n" );
  printf ( "  R8PBL_INDICATOR sets up a R8PBL indicator matrix.\n" );
  printf ( "\n" );
  printf ( "  Matrix order N = %d\n", N );
  printf ( "  Bandwidth MU = %d\n", MU );

  a = r8pbl_indicator ( N, MU );

  r8pbl_print ( N, MU, a, "  The R8PBL indicator matrix:" );

  free ( a );

  return;
# undef MU
# undef N
}
/******************************************************************************/

void test45 ( )

/******************************************************************************/
/*
  Purpose:

    TEST45 tests R8PBU_CG.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    14 April 2013

  Author:

    John Burkardt
*/
{
# define N 50
# define MU 1

  double a[(MU+1)*N];
  double *b;
  int i;
  int j;
  double *r;
  double res_max;
  double *x;
  double *x_init;

  printf ( "\n" );
  printf ( "TEST45\n" );
  printf ( "  R8PBU_CG applies the conjugate gradient method\n" );
  printf ( "  to a symmetric positive definite banded\n" );
  printf ( "  linear system.\n" );
  printf ( "\n" );
  printf ( "  Matrix order N = %d\n", N );
  printf ( "  Upper bandwidth MU = %d\n", MU );
/*
  Set the matrix values.
*/
  a[0+0*(MU+1)] = 0.0;
  for ( j = 1; j < N; j++ )
  {
    a[0+j*(MU+1)] = -1.0;
  }
  for ( j = 0; j < N; j++ )
  {
    a[1+j*(MU+1)] = 2.0;
  }

  r8pbu_print_some ( N, MU, a, 1, 1, 10, 10, "The symmetric banded matrix:" );
/*
  Set the desired solution.
*/
  x = r8vec_indicator_new ( N );
/*
  Compute the right hand side.
*/
  b = r8pbu_mxv ( N, MU, a, x );
/*
  Set the approximate solution.
*/
  x_init = ( double * ) malloc ( N * sizeof ( double ) );
  for ( i = 0; i < N; i++ )
  {
    x_init[i] = 1.0;
  }
/*
  Call the conjugate gradient method.
*/
  free ( x );

  x = r8pbu_cg ( N, MU, a, b, x_init );

  r8vec_print_some ( N, x, 1, 10, "  Solution:" );
/*
  Compute the residual, A*x-b
*/
  r = r8pbu_mxv ( N, MU, a, x );

  res_max = 0.0;
  for ( i = 0; i < N; i++ )
  { 
    res_max = r8_max ( res_max, r8_abs ( r[i] - b[i] ) );
  }

  printf ( "\n" );
  printf ( "  Maximum residual = %d\n", res_max );

  free ( r );
  free ( x );
  free ( x_init );
 
  return;
# undef MU
# undef N
}
/******************************************************************************/

void test46 ( )

/******************************************************************************/
/*
  Purpose:

    TEST46 tests R8PBU_DET.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    14 April 2013

  Author:

    John Burkardt
*/
{
# define N 10
# define MU 3

  double *a;
  double *a2;
  double *a_lu;
  double det;
  int info;
  int pivot[N];
  int seed = 123456789;

  printf ( " \n" );  
  printf ( "TEST46\n" );
  printf ( "  R8PBU_DET, determinant of a positive definite\n" );
  printf ( "  symmetric banded matrix.\n" );
  printf ( "\n" );
  printf ( "  Matrix order N = %d\n", N );
  printf ( "  Upper bandwidth MU = %d\n", MU );
/*
  Set the matrix.
*/
  a = r8pbu_random ( N, MU, &seed );

  r8pbu_print ( N, MU, a, "  The R8PBU matrix:" );
/*
  Copy the matrix into a general array.
*/
  a2 = r8pbu_to_r8ge ( N, MU, a );
/*
  Factor the matrix.
*/
  a_lu = r8pbu_fa ( N, MU, a );
  r8pbu_print ( N, MU, a_lu, "  The factored R8PBU matrix:" );
/*
  Compute the determinant.
*/
  det = r8pbu_det ( N, MU, a_lu );

  printf ( "\n" );
  printf ( "  R8PBU_DET computes the determinant = %g\n", det );
/*
  Factor the general matrix.
*/
  info = r8ge_fa ( N, a2, pivot );
/*
  Compute the determinant.
*/
  det = r8ge_det ( N, a2, pivot );

  printf ( "  R8GE_DET computes the determinant =  %g\n", det );

  free ( a );
  free ( a2 );
  free ( a_lu );

  return;
# undef MU
# undef N
}
/******************************************************************************/

void test47 ( )

/******************************************************************************/
/*
  Purpose:

    TEST47 tests R8PBU_FA, R8PBU_RANDOM, R8PBU_SL.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    14 April 2013

  Author:

    John Burkardt
*/
{
# define N 50
# define MU 1

  double *a;
  double *a_lu;
  double *b;
  int i;
  int info;
  int j;
  int seed = 123456789;
  double *x;

  printf ( "\n" );
  printf ( "TEST47\n" );
  printf ( "  For a banded positive definite symmetric matrix,\n" );
  printf ( "  R8PBU_FA computes the LU factors.\n" );
  printf ( "  R8PBU_RANDOM sets a random value.\n" );
  printf ( "  R8PBU_SL solves a linear system.\n" );
  printf ( "\n" );
  printf ( "  Matrix order N = %d\n", N );
  printf ( "  Upper bandwidth MU = %d\n", MU );
  printf ( "\n" );
/*
  Set the matrix values.
*/
  a = r8pbu_random ( N, MU, &seed );
/*
  Set the desired solution.
*/
  x = r8vec_indicator_new ( N );
/*
  Compute the right hand side.
*/
  b = r8pbu_mxv ( N, MU, a, x );
/*
  Factor the matrix.
*/
  a_lu = r8pbu_fa ( N, MU, a );
/*
  Solve the linear system.
*/
  free ( x );
  x = r8pbu_sl ( N, MU, a_lu, b );

  r8vec_print_some ( N, x, 1, 10, "  Solution:" );

  free ( a );
  free ( a_lu );
  free ( b ); 
  free ( x );

  return;
# undef MU
# undef N
}
/******************************************************************************/

void test48 ( )

/******************************************************************************/
/*
  Purpose:

    TEST48 tests R8PBU_ML.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    14 April 2013

  Author:

    John Burkardt
*/
{
# define N 10
# define MU 3

  double *a;
  double *a_lu;
  double *b;
  double *b2;
  int i;
  int info;
  int seed = 123456789;
  double *x;

  printf ( "\n" );
  printf ( "TEST48\n" );
  printf ( "  R8PBU_ML computes A*x \n" );
  printf ( "  where A has been factored by R8PBU_FA.\n" );
  printf ( "\n" );
  printf ( "  Matrix order N = %d\n", N );
  printf ( "  Upper bandwidth MU = %d\n", MU );
/*
  Set the matrix.
*/
  a = r8pbu_random ( N, MU, &seed );
/*
  Set the desired solution.
*/
  x = r8vec_indicator_new ( N );
/*
  Compute the corresponding right hand side.
*/
  b = r8pbu_mxv ( N, MU, a, x );
/*
  Factor the matrix.
*/
  a_lu = r8pbu_fa ( N, MU, a );
/*
  Now multiply the factored matrix times solution to get right hand side again.
*/
  b2 = r8pbu_ml ( N, MU, a_lu, x );

  r8vec2_print_some ( N, b, b2, 10, "  A*x and PLU*x" );

  free ( a );
  free ( a_lu );
  free ( b );
  free ( b2 );
  free ( x );

  return;
# undef MU
# undef N
}
/******************************************************************************/

void test485 ( )

/******************************************************************************/
/*
  Purpose:

    TEST485 tests R8PBU_INDICATOR.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    14 April 2013

  Author:

    John Burkardt
*/
{
# define N 9
# define MU 3

  double *a;

  printf ( " \n" );
  printf ( "TEST485\n" );
  printf ( "  R8PBU_INDICATOR sets up an SPBU indicator matrix.\n" );
  printf ( "\n" );
  printf ( "  Matrix order N = %d\n", N );
  printf ( "  Bandwidth MU = %d\n", MU );

  a = r8pbu_indicator ( N, MU );

  r8pbu_print ( N, MU, a, "  The R8PBU indicator matrix:" );

  free ( a );

  return;
# undef MU
# undef N
}
/******************************************************************************/

void test49 ( )

/******************************************************************************/
/*
  Purpose:

    TEST49 tests R8PBU_SOR.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    14 April 2013

  Author:

    John Burkardt
*/
{
# define N 50
# define MU 1

  double a[(MU+1)*N];
  double *b;
  double *b2;
  double eps;
  double err;
  int i;
  int itchk;
  int itmax;
  int j;
  int k;
  double omega;
  double pi = 3.14159265;
  double t;
  double *x;
  double x_init[N];

  printf ( "\n" );
  printf ( "TEST49\n" );
  printf ( "  R8PBU_SOR, SOR routine for iterative\n" );
  printf ( "  solution of A*x=b.\n" );
  printf ( "\n" );
  printf ( "  Matrix order N = %d\n", N );
  printf ( "  Upper bandwidth MU = %d\n", MU );

  for ( k = 1; k <= 3; k++ )
  {
    if ( k == 1 )
    {
      omega = 0.25;
    }
    else if ( k == 2 )
    {
      omega = 0.75;
    }
    else
    {
      omega = 1.00;
    }
/*
  Set matrix values.
*/
    a[0+0*(MU+1)] = 0.0;
    for ( j = 1; j < N; j++ )
    {
      a[0+j*(MU+1)] = -1.0;
    }

    for ( j = 0; j < N; j++ )
    {
      a[1+j*(MU+1)] = 2.0;
    }
/*
  Set the desired solution.
*/
    x = ( double * ) malloc ( N * sizeof ( double ) );
    for ( i = 0; i < N; i++ )
    {
      t = pi * ( double ) ( i ) / ( double ) ( N - 1 );
      x[i] = sin ( t );
    }
/*
  Compute the right hand side.
*/
    b = r8pbu_mxv ( N, MU, a, x );
/*
  Set the initial solution estimate.
*/
    for ( i = 0; i < N; i++ )
    {
      x_init[i] = 1.0;
    }
 
    itchk = 1;
    itmax = 8000;
    eps = 0.0001;

    free ( x );
    x = r8pbu_sor ( N, MU, a, b, eps, itchk, itmax, omega, x_init );
/*
  Compute residual, A*x-b
*/
    b2 = r8pbu_mxv ( N, MU, a, x );
 
    err = 0.0;
    for ( i = 0; i < N; i++ )
    {
      err = r8_max ( err, r8_abs ( b2[i] - b[i] ) );
    }
 
    printf ( "\n" );
    printf ( "SOR iteration.\n" );
    printf ( "\n" );
    printf ( "  Relaxation factor OMEGA = %g\n", omega );

    r8vec_print_some ( N, x, 1, 10, "  Solution:" );

    printf ( "\n" );
    printf ( "  Maximum error = %g\n", err );

    free ( b );
    free ( b2 );
    free ( x );
  }
 
  return;
# undef MU
# undef N
}
/******************************************************************************/

void test50 ( )

/******************************************************************************/
/*
  Purpose:

    TEST50 tests R8PO_FA, R8PO_SL.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    14 April 2013

  Author:

    John Burkardt
*/
{
# define N 10

  double a[N*N];
  double *b;
  int i;
  int info;
  int j;
  double *r;
  double *x;

  printf ( "\n" );
  printf ( "TEST50\n" );
  printf ( "  R8PO_FA factors a positive definite symmetric\n" );
  printf ( "  linear system,\n" );
  printf ( "  R8PO_SL solves a factored system.\n" );
  printf ( "\n" );
  printf ( "  Matrix order N = %d\n", N );

  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      a[i+j*N] = ( double ) ( i4_min ( i + 1, j + 1 ) );
    }
  }
/*
  Set the desired solution.
*/
  x = r8vec_indicator_new ( N );
/*
  Compute the corresponding right hand side.
*/
  b = r8po_mxv ( N, a, x );
/*
  Factor the matrix.
*/
  r = r8po_fa ( N, a );

  if ( r == NULL )
  {
    printf ( "\n" );
    printf ( "  Fatal error!\n" );
    printf ( "  R8PO_FA declares the matrix is singular!\n" );
    return;
  }
/*
  Solve the linear system.
*/
  free ( x );
  x = r8po_sl ( N, r, b );
 
  r8vec_print ( N, x, "  Solution:" );
/*
  Set the desired solution.
*/
  free ( x );
  x = ( double * ) malloc ( N * sizeof ( double ) );

  for ( i = 0; i < N; i++ )
  {
    x[i] = 1.0;
  }
/*
  Compute the corresponding right hand side, using the factored matrix.
*/
  free ( b );
  b = r8po_ml ( N, r, x );
/*
  Solve the linear system.
*/
  free ( x );
  x = r8po_sl ( N, r, b );
 
  r8vec_print ( N, x, "  Solution:" );

  free ( b );
  free ( r );
  free ( x );

  return;
# undef N
}
/******************************************************************************/

void test505 ( )

/******************************************************************************/
/*
  Purpose:

    TEST505 tests R8PO_FA;

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    15 April 2013

  Author:

    John Burkardt
*/
{
# define N 5

  double a[N*N];
  int i;
  int j;
  int k;
  double *r;

  printf ( "\n" );
  printf ( "TEST505\n" );
  printf ( "  R8PO_FA factors a positive definite symmetric\n" );
  printf ( "  linear system,\n" );
  printf ( "\n" );
  printf ( "  Matrix order N = %d\n", N );

  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      a[i+j*N] = ( double ) ( i4_min ( i + 1, j + 1 ) );
    }
  }

  r8po_print ( N, a, "  The matrix A:" );
/*
  Factor the matrix.
*/
  r = r8po_fa ( N, a );

  if ( r == NULL )
  {
    printf ( "\n" );
    printf ( "  Fatal error!\n" );
    printf ( "  R8PO_FA declares the matrix is singular!\n" );
    return;
  }
  r8ut_print ( N, N, r, "  The factor R (an R8UT matrix):" );
/*
  Compute the product R' * R.
*/
  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      a[i+j*N] = 0.0;
      for ( k = 0; k < N; k++ )
      {
        a[i+j*N] = a[i+j*N] + r[k+i*N] * r[k+j*N];
      }
    }
  }

  r8ge_print ( N, N, a, "  The product R' * R:" );

  free ( r );

  return;
# undef N
}
/******************************************************************************/

void test51 ( )

/******************************************************************************/
/*
  Purpose:

    TEST51 tests R8PO_DET, R8PO_INVERSE.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    15 April 2013

  Author:

    John Burkardt
*/
{
# define N 4

  double a[N*N];
  double *a_inv;
  double *c;
  double det;
  int i;
  int info;
  int j;
  int k;
  double *r;

  printf ( "\n" );
  printf ( "TEST51\n" );
  printf ( "  For a symmetric positive definite matrix\n" );
  printf ( "  factored by R8PO_FA,\n" );
  printf ( "  R8PO_DET computes the determinant;\n" );
  printf ( "  R8PO_INVERSE computes the inverse.\n" );
  printf ( "\n" );
  printf ( "  Matrix order N = %d\n", N );
/*
  Set the matrix.
*/
  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      a[i+j*N] = ( double ) ( i4_min ( i + 1, j + 1 ) );
    }
  }

  r8po_print ( N, a, "  Matrix A:" );
/*
  Factor the matrix.
*/
  r = r8po_fa ( N, a );
/*
  Compute the determinant.
*/
  det = r8po_det ( N, r );
 
  printf ( "\n" );
  printf ( "  Matrix determinant = %g\n", det );
/*
  Compute the inverse.
*/
  a_inv = r8po_inverse ( N, r );

  r8po_print ( N, a_inv, "  Inverse matrix A_INV:" );
/*
  Check.
*/
  c = r8po_mxm ( N, a, a_inv );

  r8po_print ( N, c, "  Product A * A_INV:" );

  free ( a_inv );
  free ( c );
  free ( r );

  return;
# undef N
}
/******************************************************************************/

void test515 ( )

/******************************************************************************/
/*
  Purpose:

    TEST515 tests R8PO_FA, R8PO_SL, R8PO_ML, R8PO_MXV.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    15 April 2013

  Author:

    John Burkardt
*/
{
# define N 5

  double *a;
  double *b;
  int info;
  int seed = 123456789;
  double *r;
  double *x;

  printf ( "\n" );
  printf ( "TEST515\n" );
  printf ( "  For a positive definite symmetric matrix,\n" );
  printf ( "  R8PO_FA computes the Cholesky factor.\n" );
  printf ( "  R8PO_SL solves a factored linear system.\n" );
  printf ( "  R8PO_MXV multiplies unfactored A * x\n" );
  printf ( "  R8PO_ML multiplies factored A * x\n" );
  printf ( "\n" );
  printf ( "  Matrix order N = %d\n", N );
/*
  Set the matrix.
*/
  a = r8po_random ( N, &seed );

  r8po_print ( N, a, "  The matrix A:" );
/*
  Compute the desired right hand side.
*/
  x = r8vec_indicator_new ( N );

  b = r8po_mxv ( N, a, x );

  r8vec_print ( N, b, "  Right hand side, computed by R8PO_MXV" );
/*
  Factor the matrix.
*/
  r = r8po_fa ( N, a );
/*
  Solve the linear system.
*/
  free ( x );
  x = r8po_sl ( N, r, b );

  r8vec_print ( N, x, "  Solution (should be 1,2,3...)" );
/*
  Recompute the desired right hand side.
  Since A has been factored, we have to call R8PO_ML now.
*/
  free ( x );
  x = r8vec_indicator_new ( N );

  free ( b );
  b = r8po_ml ( N, a, x );

  r8vec_print ( N, b, "  Right hand side, computed by R8PO_ML" );
/*
  Solve the linear system.
*/
  free ( x );
  x = r8po_sl ( N, a, b );

  r8vec_print ( N, x, "  Solution (should be 1,2,3...)" );

  free ( a );
  free ( b );
  free ( r );
  free ( x );

  return;
}
/******************************************************************************/

void test517 ( )

/******************************************************************************/
/*
  Purpose:

    TEST517 tests R8PO_INDICATOR.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 April 2013

  Author:

    John Burkardt
*/
{
# define N 5

  double *a;

  printf ( "\n" );
  printf ( "TEST517\n" );
  printf ( "  R8PO_INDICATOR sets up an R8PO indicator matrix.\n" );
  printf ( "\n" );
  printf ( "  Matrix order N = %d\n", N );
/*
  Set the matrix.
*/
  a = r8po_indicator ( N );

  r8po_print ( N, a, "  The R8PO indicator matrix:" );
 
  free ( a );

  return;
# undef N
}
/******************************************************************************/

void test52 ( )

/******************************************************************************/
/*
  Purpose:

    TEST52 tests R8PO_RANDOM, R8PO_TO_R8GE.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 April 2013

  Author:

    John Burkardt
*/
{
# define N 5

  double *a;
  double *b;
  int seed = 123456789;

  printf ( "\n" );
  printf ( "TEST52\n" );
  printf ( "  R8PO_RANDOM computes a random positive definite\n" );
  printf ( "  symmetric matrix.\n" );
  printf ( "  R8PO_TO_R8GE converts an R8PO matrix to R8GE format.\n" );
  printf ( "\n" );
  printf ( "  Matrix order N = %d\n", N );
/*
  Set the matrix.
*/
  a = r8po_random ( N, &seed );

  r8po_print ( N, a, "  The random R8PO matrix:" );
 
  r8ge_print ( N, N, a, "  The random R8PO matrix (printed by R8GE_PRINT):" );

  b = r8po_to_r8ge ( N, a );

  r8ge_print ( N, N, b, "  The random R8GE matrix (printed by R8GE_PRINT):" );

  free ( a );
  free ( b );

  return;
# undef N
}
/******************************************************************************/

void test525 ( )

/******************************************************************************/
/*
  Purpose:

    TEST525 tests R8PP_FA, R8PP_SL.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 April 2013

  Author:

    John Burkardt
*/
{
# define N 5

  double *a;
  double *b;
  int i;
  int info;
  int j;
  double *r;
  int seed = 123456789;
  double *x;

  printf ( "\n" );
  printf ( "TEST525\n" );
  printf ( "  R8PP_FA factors an R8PP system,\n" );
  printf ( "  R8PP_SL solves an R8PP system.\n" );
  printf ( "\n" );
  printf ( "  Matrix order N = %d\n", N );
/*
  Set the matrix.
*/
  a = r8pp_random ( N, &seed );

  r8pp_print ( N, a, "  The R8PP matrix:" );
/*
  Set the desired solution.
*/
  x = r8vec_indicator_new ( N );

  r8vec_print ( N, x, "  The desired solution:" );
/*
  Compute the corresponding right hand side.
*/
  b = r8pp_mxv ( N, a, x );

  r8vec_print ( N, b, "  The right hand side:" );
/*
  Factor the matrix.
*/
  r = r8pp_fa ( N, a );

  if ( r == NULL )
  {
    printf ( "\n" );
    printf ( "  Fatal error!\n" );
    printf ( "  R8PP_FA declares the matrix is singular!\n" );
    return;
  }
  printf ( "\n" );
  printf ( "  The R8PP matrix has been factored.\n" );
/*
  Solve the linear system.
*/
  free ( x );
  x = r8pp_sl ( N, r, b );
 
  r8vec_print ( N, x, "  Solution:" );

  free ( a );
  free ( b );
  free ( r );
  free ( x );

  return;
# undef N
}
/******************************************************************************/

void test527 ( )

/******************************************************************************/
/*
  Purpose:

    TEST527 tests R8PP_INDICATOR.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 April 2013

  Author:

    John Burkardt
*/
{
# define N 5

  double *a;

  printf ( "\n" );
  printf ( "TEST527\n" );
  printf ( "  R8PP_INDICATOR sets up an R8PP indicator matrix.\n" );
  printf ( "\n" );
  printf ( "  Matrix order N = %d\n", N );

  a = r8pp_indicator ( N );

  r8pp_print ( N, a, "  The R8PP indicator matrix:" );
 
  free ( a );

  return;
# undef N
}
/******************************************************************************/

void test53 ( )

/******************************************************************************/
/*
/  Purpose:

    TEST53 tests R8PP_RANDOM.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 April 2013

  Author:

    John Burkardt
*/
{
# define N 5

  double *a;
  double *b;
  int seed = 123456789;

  printf ( "\n" );
  printf ( "TEST53\n" );
  printf ( "  R8PP_RANDOM, compute a random positive definite\n" );
  printf ( "  symmetric packed matrix.\n" );
  printf ( "\n" );
  printf ( "  Matrix order N = %d\n", N );
/*
  Set the matrix.
*/
  a = r8pp_random ( N, &seed );

  r8pp_print ( N, a, "  The matrix (printed by R8PP_PRINT):" );
 
  b = r8pp_to_r8ge ( N, a );

  r8ge_print ( N, N, b, "  The random R8PP matrix (printed by R8GE_PRINT):" );

  free ( a );
  free ( b );

  return;
# undef N
}
/******************************************************************************/

void test534 ( )

/******************************************************************************/
/*
  Purpose:

    TEST534 tests R8S3_WRITE.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    18 April 2013

  Author:

    John Burkardt
*/
{
# define N 100
# define NZ_NUM ( 3 * N - 2 )

  double a[NZ_NUM];
  int col[NZ_NUM];
  int i;
  int isym;
  int j;
  int k;
  char *output_file = "r8s3_matrix.txt";
  int row[NZ_NUM];

  printf ( "\n" );
  printf ( "TEST534\n" );
  printf ( "  For a R8S3 matrix,\n" );
  printf ( "  R8S3_WRITE writes the matrix to a file.\n" );
  printf ( "\n" );
  printf ( "  Matrix order N =         %d\n", N );
  printf ( "  Matrix nonzeros NZ_NUM = %d\n", NZ_NUM );

  isym = 0;
/*
  Set the matrix values.
*/
  k = 0;
  for ( i = 1; i <= N; i++ )
  {

    j = i;
    row[k] = i;
    col[k] = j;
    a[k] = ( double ) ( 100 * i + j );
    k = k + 1;
  }

  for ( i = 2; i <= N; i++ )
  {
    j = i - 1;
    row[k] = i;
    col[k] = j;
    a[k] = ( double ) ( 100 * i + j );
    k = k + 1;
  }

  for ( i = 1; i <= N-1; i++ )
  {
    j = i + 1;
    row[k] = i;
    col[k] = j;
    a[k] = ( double ) ( 100 * i + j );
    k = k + 1;
  }

  r8s3_print_some ( N, N, NZ_NUM, isym, row, col, a, 1, 1, 
    10, 10, "  Initial 10x10 block of R8S3 matrix:" );

  r8s3_write ( N, NZ_NUM, isym, row, col, a, output_file );

  printf ( "  R8S3_WRITE wrote the matrix data to \"%s\".\n", output_file );

  return;
# undef N
# undef NZ_NUM
}
/******************************************************************************/

void test535 ( )

/******************************************************************************/
/*
  Purpose:

    TEST535 tests R8S3_READ, R8S3_READ_SIZE..

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    18 April 2013

  Author:

    John Burkardt
*/
{
  double *a;
  int *col;
  int isym;
  int n;
  int nz_num;
  char *input_file = "r8s3_matrix.txt";
  int *row;

  printf ( "\n" );
  printf ( "TEST535\n" );
  printf ( "  For a R8S3 matrix,\n" );
  printf ( "  R8S3_READ reads a matrix from a file.\n" );
  printf ( "  R8S3_READ_SIZE reads the sizes of the matrix from a file.\n" );

  r8s3_read_size ( input_file, &n, &nz_num );

  printf ( "\n" );
  printf ( "  R8S3_READ_SIZE reports matrix size data:\n" );
  printf ( "\n" );

  printf ( "\n" );
  printf ( "  Matrix order N =         %d\n", n );
  printf ( "  Matrix nonzeros NZ_NUM = %d\n", nz_num );

  row = ( int * ) malloc ( nz_num * sizeof ( int ) );
  col = ( int * ) malloc ( nz_num * sizeof ( int ) );
  a = ( double * ) malloc ( nz_num * sizeof ( double ) );

  r8s3_read ( input_file, n, nz_num, row, col, a );

  isym = 0;

  r8s3_print_some ( n, n, nz_num, isym, row, col, a, 1, 1, 
    10, 10, "  Initial 10x10 block of R8S3 matrix:" );

  printf ( "\n" );
  printf ( "  Deleting the matrix data file \"%s\"\n", input_file );

  file_delete ( input_file );

  free ( a );
  free ( col );
  free ( row );

  return;
}
/******************************************************************************/

void test54 ( )

/******************************************************************************/
/*
  Purpose:

    TEST54 tests R8SD_CG.

  Discussion:

    NX and NY are the number of grid points in the X and Y directions.
    N is the number of unknowns.
    NDIAG is the number of nonzero diagonals we will store.  We only
    store the main diagonal, and the superdiagonals.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    18 April 2013

  Author:

    John Burkardt
*/
{
# define NDIAG 3
# define NX 10
# define NY 10
# define N ( NX * NY )

  double a[N*NDIAG];
  double *b;
  double *b2;
  double err;
  int i;
  int j;
  int k;
  int offset[NDIAG] = { 0, 1, NX };
  double *x;
  double x_init[N];

  printf ( "\n" );
  printf ( "TEST54\n" );
  printf ( "  R8SD_CG applies the conjugate gradient method\n" );
  printf ( "  to a symmetric positive definite linear\n" );
  printf ( "  system stored by diagonals.\n" );
  printf ( "\n" );
  printf ( "  Matrix order N = %d\n", N );
  printf ( "  Matrix diagonals NDIAG = %d\n", NDIAG );
  printf ( "\n" );
/*
  OFFSET tells us where the nonzero diagonals are.  It does this
  by recording how "high" or to the right the diagonals are from
  the main diagonal.


  Now we compute the numbers that go into the diagonals.  For this
  problem, we could simply store a column of 4's, and two columns of
  -1's, but I wanted to go through the motions of thinking about the
  value of each entry.  "K" counts the row of the original matrix
  that we are working on.
*/
  k = 0;

  for ( j = 1; j <= NY; j++ )
  {
    for ( i = 1; i <= NX; i++ )
    {
/*
  Central
*/
      a[k+0*N] = 4.0;
/*
  East ( = West )
*/
      if ( i == NX )
      {
        a[k+1*N] = 0.0;
      }
      else
      {
        a[k+1*N] = -1.0;
      }
/*
  North ( = South )
*/
      if ( j == NY )
      {
        a[k+2*N] = 0.0;
      }
      else
      {
        a[k+2*N] = -1.0;
      }
      k = k + 1;
    }
  }
/*
  Print some of the matrix.
*/
  r8sd_print_some ( N, NDIAG, offset, a, 1, 1, 10, 10, 
    "  First 10 rows and columns of matrix." );
/*
  Set the desired solution.
*/
  x = ( double * ) malloc ( N * sizeof ( double ) );
  k = 0;
  for ( j = 1; j <= NY; j++ )
  {
    for ( i = 1; i <= NX; i++ )
    {
      x[k] = ( double ) ( 10 * i + j );
      k = k + 1;
    }
  }
/*
  Compute the corresponding right hand side.
*/
  b = r8sd_mxv ( N, NDIAG, offset, a, x );

  r8vec_print_some ( N, b, 1, 10, "  Right hand side:" );
/*
  Set X to zero so no one accuses us of cheating.
*/
  for ( i = 0; i < N; i++ )
  {
    x_init[i] = 0.0;
  }
/*
  Call the conjugate gradient method.
*/
  free ( x );
  x = r8sd_cg ( N, NDIAG, offset, a, b, x_init );
/*
  Compute the residual, A*x-b
*/
  b2 = r8sd_mxv ( N, NDIAG, offset, a, x );
 
  err = 0.0;
  for ( i = 0; i < N; i++ )
  {
    err = r8_max ( err, r8_abs ( b2[i] - b[i] ) );
  }
 
  r8vec_print_some ( N, x, 1, 10, "  Solution:" );

  printf ( "\n" );
  printf ( "  Maximum residual = %g\n", err );
/*
  Note that if we're not satisfied with the solution, we can
  call again, using the computed X as our starting estimate.

  Call the conjugate gradient method AGAIN.
*/
  for ( i = 0; i < N; i++ )
  {
    x_init[i] = x[i];
  }
  free ( x );

  x = r8sd_cg ( N, NDIAG, offset, a, b, x_init );
/*
  Compute the residual, A*x-b
*/
  free ( b2 );
  b2 = r8sd_mxv ( N, NDIAG, offset, a, x );
 
  err = 0.0;
  for ( i = 0; i < N; i++ )
  {
    err = r8_max ( err, r8_abs ( b2[i] - b[i] ) );
  }
 
  r8vec_print_some ( N, x, 1, 10, "  Second attempt at solution:" );

  printf ( "\n" );
  printf ( "  Maximum residual of second attempt = %g\n", err );

  free ( b );
  free ( b2 );
  free ( x );

  return;
# undef N
# undef NDIAG
# undef NX
# undef NY
}
/******************************************************************************/

void test55 ( )

/******************************************************************************/
/*
  Purpose:

    TEST55 tests R8SD_CG.

  Discussion:

    This is a sample demonstration of how to compute some eigenvalues
    and corresponding eigenvectors of a matrix.  The matrix is the
    discretized Laplacian operator, which can be stored by diagonals,
    and handled by the conjugate gradient method.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    20 April 2013

  Author:

    John Burkardt
*/
{
# define MAXVEC 3
# define NDIAG 3
# define NX 10
# define NY 10
# define N ( NX * NY )
# define PI 3.141592653589

  double a[N*NDIAG];
  double del;
  double dot;
  double eval;
  int i;
  int iter;
  int i4vec;
  int j;
  int k;
  double lambda;
  double lambda_old;
  double lamvec[MAXVEC];
  double norm;
  int nvec;
  int offset[NDIAG] = { 0, 1, NX };
  double vec[N*MAXVEC];
  double x[N];
  double *xnew;

  printf ( "\n" );
  printf ( "TEST55\n" );
  printf ( "  R8SD_CG is used for linear equation solving\n" );
  printf ( "  in a demonstration of inverse iteration to\n" );
  printf ( "  compute eigenvalues and eigenvectors of a\n" );
  printf ( "  symmetric matrix stored by diagonals.\n" );
  printf ( "\n" );
  printf ( "  Matrix order N = %d\n", N );
  printf ( "  Matrix diagonals NDIAG = %d\n", NDIAG );
  printf ( "\n" );
  printf ( "  Here are 25 of the smallest eigenvalues:\n" );
  printf ( "\n" );
  printf ( "  I, J, eigenvalue(I,J):\n" );
  printf ( "\n" );

  for ( i = 1; i <= i4_min ( 5, NX ); i++ )
  {
    for ( j = 1; j <= i4_min ( 5, NY ); j++ )
    {
      eval = 4.0 
        - 2.0 * cos ( ( double ) ( i ) * PI / ( double ) ( NX + 1 ) ) 
        - 2.0 * cos ( ( double ) ( j ) * PI / ( double ) ( NY + 1 ) );
      printf ( "  %6d  %6d  %12g\n", i, j, eval );
    }
  }
/*
  Now we compute the numbers that go into the diagonals.  For this
  problem, we could simply store a column of 4's, and two columns of
  -1's, but I wanted to go through the motions of thinking about the
  value of each entry.  "K" counts the row of the original matrix
  that we are working on.
*/
  k = 0;
  for ( j = 1; j <= NY; j++ )
  {
    for ( i = 1; i <= NX; i++ )
    {
/*
  Central
*/
      a[k+0*N] = 4.0;
/*
  East ( = West )
*/
      if ( i == NX )
      {
        a[k+1*N] = 0.0;
      }
      else
      {
        a[k+1*N] = -1.0;
      }
/*
  North ( = South )
*/
      if ( j == NY )
      {
        a[k+2*N] = 0.0;
      }
      else
      {
        a[k+2*N] = -1.0;
      }
      k = k + 1;
    }
  }

  nvec = 0;
/*
  Set the starting eigenvector and eigenvalue estimates.
*/
  for ( ; ; )
  {
    printf ( "\n" );

    lambda = 0.0;

    k = 0;
    for ( j = 1; j <= NY; j++ )
    {
      for ( i = 1; i <= NX; i++ )
      {
        x[k] = 1.0;
        k = k + 1;
      }
    }
/*
  Remove any components of previous eigenvectors.
*/
    for ( i4vec = 0; i4vec < nvec; i4vec++ )
    {
      dot = 0.0;
      for ( i = 0; i < N; i++ )
      {
        dot = dot + x[i] * vec[i+i4vec*N];
      }
      for ( i = 0; i < N; i++ )
      {
        x[i] = x[i] - dot * vec[i+i4vec*N];
      }
    }

    xnew = ( double * ) malloc ( N * sizeof ( double ) );
    for ( i = 0; i < N; i++ )
    {
      xnew[i] = x[i];
    }
/*
  Iterate
*/
    for ( iter = 1; iter <= 40; iter++ )
    {
      norm = 0.0;
      for ( i = 0; i < N; i++ )
      {
        norm = norm + xnew[i] * xnew[i];
      }
      norm = sqrt ( norm );

      for ( i = 0; i < N; i++ )
      {
        xnew[i] = xnew[i] / norm;
      }

      lambda_old = lambda;
      lambda = 1.0 / norm;
/*
  Check for convergence.
*/
      if ( 1 < iter )
      {
        del = r8_abs ( lambda - lambda_old );

        if ( del < 0.000001 )
        {
          printf ( "\n" );
          printf ( "Lambda estimate = %g\n", lambda );
          printf ( "Converged on step %d\n", iter );
          break;
        }
      }
/*
  Call the conjugate gradient method, solving A * XNEW = X.
*/
      for ( i = 0; i < N; i++ )
      {
        x[i] = xnew[i];
      }
      free ( xnew );

      xnew = r8sd_cg ( N, NDIAG, offset, a, x, x );

      for ( i4vec = 0; i4vec < nvec; i4vec++ )
      {
        dot = 0.0;
        for ( i = 0; i < N; i++ )
        {
          dot = dot + xnew[i] * vec[i+i4vec*N];
        }
        for ( i = 0; i < N; i++ )
        {
          xnew[i] = xnew[i] - dot * vec[i+i4vec*N];
        }
      }

    }

    lamvec[nvec] = lambda;
    for ( i = 0; i < N; i++ )
    {
      vec[i+nvec*N] = xnew[i];
    }
    nvec = nvec + 1;

    if ( MAXVEC <= nvec )
    {
      break;
    }

  }

  free ( xnew );

  return;
# undef MAXVEC
# undef N
# undef NDIAG
# undef NX
# undef NY
# undef PI
}
/******************************************************************************/

void test555 ( )

/******************************************************************************/
/*
  Purpose:

    TEST555 tests R8SD_INDICATOR.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    20 April 2013

  Author:

    John Burkardt
*/
{
# define N 10
# define NDIAG 3

  double *a;
  int offset[NDIAG] = { 0, 1, 3 };

  printf ( "\n" );
  printf ( "TEST555\n" );
  printf ( "  R8SD_INDICATOR sets up an R8SD indicator matrix.\n" );
  printf ( "\n" );
  printf ( "  Matrix order N = %d\n", N );
  printf ( "  Matrix diagonals NDIAG = %d\n", NDIAG );

  a = r8sd_indicator ( N, NDIAG, offset );

  r8sd_print ( N, NDIAG, offset, a, "  The R8SD indicator matrix:" );

  free ( a );

  return;
# undef N
# undef NDIAG
}
/******************************************************************************/

void test56 ( )

/******************************************************************************/
/*
  Purpose:

    TEST56 tests R8SM_ML.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    20 April 2013

  Author:

    John Burkardt
*/
{
# define M 7
# define N 7

  double a[M*N];
  double *b;
  double *b2;
  int info;
  int i;
  int job;
  int pivot[N];
  int seed = 123456789;
  double u[M];
  double v[N];
  double *x;

  printf ( "\n" );
  printf ( "TEST56\n" );
  printf ( "  R8SM_ML computes A*x or A'*X\n" );
  printf ( "  where A is a Sherman Morrison matrix.\n" );
  printf ( "\n" );
  printf ( "  Matrix rows M =    %d\n", M );
  printf ( "  Matrix columns N = %d\n", N );

  for ( job = 0; job <= 1; job++ )
  {
/*
  Set the matrix.
*/
    r8sm_random ( M, N, a, u, v, &seed );

    r8sm_print ( M, N, a, u, v, "  The Sherman Morrison matrix:" );
/*
  Set the desired solution.
*/
    x = r8vec_indicator_new ( N );
/*
  Compute the corresponding right hand side.
*/
    if ( job == 0 )
    {
      b = r8sm_mxv ( M, N, a, u, v, x );
    }
    else
    {
      b = r8sm_vxm ( M, N, a, u, v, x );
    }
/*
  Factor the matrix.
*/
    info = r8ge_fa ( N, a, pivot );

    if ( info != 0 )
    {
      printf ( "\n" );
      printf ( "  Fatal error!\n" );
      printf ( "  R8GE_FA declares the matrix is singular!\n" );
      printf ( "  The value of INFO is %d\n", info );
      return;
    }
/*
  Now multiply factored matrix times solution to get right hand side again.
*/
    b2 = r8sm_ml ( N, a, u, v, pivot, x, job );

    if ( job == 0 )
    {
      r8vec2_print_some ( N, b, b2, 10, "  A*x and PLU*x" );
    }
    else
    {
      r8vec2_print_some ( N, b, b2, 10, "  A'*x and (PLU)'*x" );
    }

    free ( b );
    free ( b2 );
    free ( x );
  }

  return;
# undef M
# undef N
}
/******************************************************************************/

void test57 ( )

/******************************************************************************/
/*
  Purpose:

    TEST57 tests R8SM_SL.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    19 April 2013

  Author:

    John Burkardt
*/
{
# define M 5
# define N 5

  double a[M*N];
  double *b;
  int i;
  int ierror;
  int info;
  int job;
  int pivot[N];
  int seed = 123456789;
  double u[M];
  double v[N];
  double *x;

  printf ( "\n" );
  printf ( "TEST57\n" );
  printf ( "  R8SM_SL implements the Sherman-Morrison method \n" );
  printf ( "  for solving a perturbed linear system.\n" );
  printf ( "\n" );
  printf ( "  Matrix rows M =    %d\n", M );
  printf ( "  Matrix columns N = %d\n", N );

  for ( job = 1; 0 <= job; job-- )
  {
/*
  Set the matrix.
*/
    r8sm_random ( M, N, a, u, v, &seed );

    r8sm_print ( M, N, a, u, v, "  The Sherman-Morrison matrix A:" );
/*
  Set the desired solution.
*/
    x = r8vec_indicator_new ( N );
/*
  Compute the corresponding right hand side.
*/
    if ( job == 0 )
    {
      b = r8sm_mxv ( M, N, a, u, v, x );
    }
    else
    {
      b = r8sm_vxm ( M, N, a, u, v, x );
    }

    r8vec_print ( N, b, "  The right hand side vector B:" );
/*
  Factor the matrix.
*/
    info = r8ge_fa ( N, a, pivot );

    if ( info != 0 )
    {
      printf ( "\n" );
      printf ( "  Fatal error!\n" );
      printf ( "  R8GE_FA declares the matrix is singular!\n" );
      printf ( "  The value of INFO is %d\n", info );
      continue;
    }
/*
  Solve the linear system.
*/
    b = r8sm_sl ( N, a, u, v, b, pivot, job );
 
    if ( job == 0 )
    {
      r8vec_print ( N, b, "  Solution to A * X = B:" );
    }
    else
    {
      r8vec_print ( N, b, "  Solution to A' * X = B:" );
    }

    free ( b );
  }

  return;
# undef M
# undef N
}
/******************************************************************************/

void test5705 ( )

/******************************************************************************/
/*
  Purpose:

    TEST5705 tests R8SP_IJ_TO_K.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    20 April 2013

  Author:

    John Burkardt
*/
{
# define NZ_NUM 10

  int check;
  int col[NZ_NUM] = { 2, 5, 1, 5, 1, 2, 3, 4, 4, 1 };
  int i;
  int j;
  int k;
  int m = 7;
  int n = 5;
  int nz_num = NZ_NUM;
  int row[NZ_NUM] = { 1, 1, 2, 2, 4, 4, 4, 5, 6, 7 };

  printf ( "\n" );
  printf ( "TEST5705\n" );
  printf ( "  R8SP_IJ_TO_K returns the R8SP index of (I,J).\n" );
  printf ( "\n" );
  printf ( "  Matrix rows M =    %d\n", m );
  printf ( "  Matrix columns N = %d\n", n );
  printf ( "  Matrix nonzeros =  %d\n", nz_num );

  check = r8sp_check ( m, n, nz_num, row, col );

  if ( !check )
  {
    printf ( "\n" );
    printf ( "R8SP_CHECK - Error!\n" );
    printf ( "  The matrix is not in the proper sorted format.\n" );
    return;
  }

  printf ( "\n" );
  printf ( "         I         J         K\n" );
  printf ( "\n" );
  for ( i = 1; i <= m; i++ )
  {
    for ( j = 1; j <= n; j++ )
    {
      k = r8sp_ij_to_k ( nz_num, row, col, i, j );

      printf ( "  %8d  %8d  %8d\n", i, j, k );
    }
  }

  return;
# undef NZ_NUM
}
/******************************************************************************/

void test571 ( )

/******************************************************************************/
/*
  Purpose:

    TEST571 tests R8SP_INDICATOR.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    20 April 2013

  Author:

    John Burkardt
*/
{
# define M 7
# define N 5
# define NZ_NUM 10

  double *a;
  int col[NZ_NUM] = { 2, 5, 1, 5, 1, 2, 3, 4, 4, 1 };
  int row[NZ_NUM] = { 1, 1, 2, 2, 4, 4, 4, 5, 6, 7 };

  printf ( "\n" );
  printf ( "TEST571\n" );
  printf ( "  R8SP_INDICATOR sets up an R8SP indicator matrix;\n" );
  printf ( "\n" );
  printf ( "  Matrix rows M =          %d\n", M );
  printf ( "  Matrix columns N =       %d\n", N );
  printf ( "  Matrix nonzeros NZ_NUM = %d\n", NZ_NUM );

  a = r8sp_indicator ( M, N, NZ_NUM, row, col );

  r8sp_print ( M, N, NZ_NUM, row, col, a, "  The R8SP indicator matrix:" );

  free ( a );

  return;
# undef M
# undef N
# undef NZ_NUM
}
/******************************************************************************/

void test572 ( )

/******************************************************************************/
/*
  Purpose:

    TEST572 tests R8SP_MXV, R8SP_VXM.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    20 April 2013

  Author:

    John Burkardt
*/
{
# define M 7
# define N 5
# define NZ_NUM 10

  double *a;
  double *b;
  double *c;
  int col[NZ_NUM] = { 2, 5, 1, 5, 1, 2, 3, 4, 4, 1 };
  int i;
  int row[NZ_NUM] = { 1, 1, 2, 2, 4, 4, 4, 5, 6, 7 };
  int seed = 123456789;
  double *x;

  printf ( "\n" );
  printf ( "TEST572\n" );
  printf ( "  R8SP_MXV multiplies an R8SP matrix by a vector;\n" );
  printf ( "  R8SP_VXM multiplies a vector by an R8SP matrix;\n" );
  printf ( "\n" );
  printf ( "  Matrix rows M =          %d\n", M );
  printf ( "  Matrix columns N =       %d\n", N );
  printf ( "  Matrix nonzeros NZ_NUM = %d\n", NZ_NUM );
/*
  Set the matrix.
*/
  a = r8sp_random ( M, N, NZ_NUM, row, col, &seed );
/*
  Make an R8GE copy.
*/
  c = r8sp_to_r8ge ( M, N, NZ_NUM, row, col, a );
/*
  Print the R8GE copy.
*/
  r8ge_print ( M, N, c, "  The R8SP matrix, in R8GE form:" );

  x = ( double * ) malloc ( N * sizeof ( double ) );

  x[0] = 1.0;
  for ( i = 1; i < N-1; i++ )
  {
    x[i] = 0.0;
  }
  x[N-1] = -1.0;

  r8vec_print ( N, x, "  The vector x:" );

  b = r8sp_mxv ( M, N, NZ_NUM, row, col, a, x );

  r8vec_print ( M, b, "  The product A * x:" );

  free ( x );
  x = ( double * ) malloc ( M * sizeof ( double ) );

  x[0] = 1.0;
  for ( i = 1; i < M-1; i++ )
  {
    x[i] = 0.0;
  }
  x[M-1] = -1.0;

  r8vec_print ( M, x, "  The vector x:" );

  free ( b );

  b = r8sp_vxm ( M, N, NZ_NUM, row, col, a, x );

  r8vec_print ( N, b, "  The product A' * x:" );

  free ( a );
  free ( b );
  free ( c );
  free ( x );

  return;
# undef M
# undef N
# undef NZ_NUM
}
/******************************************************************************/

void test5722 ( )

/******************************************************************************/
/*
  Purpose:

    TEST5722 tests R8SP_PRINT.

  Discussion:

    Because MATLAB seems to allow a R8SP matrix to store the same index
    several times, presumably with the matrix entry being the SUM of
    these occurrences, I modified R8SP_PRINT to handle this situation
    (I hope).

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    20 April 2013

  Author:

    John Burkardt
*/
{
# define M 5
# define N 7
# define NZ_NUM 12

  double a[NZ_NUM] = {
    21.0,  51.0, 12.0, 52.0, 14.0, 
    24.0,  34.0, 45.0, 46.0, 17.0, 
   100.0, 200.0 };
  int col[NZ_NUM] = { 1, 1, 2, 2, 4, 4, 4, 5, 6, 7, 2, 4 };
  int row[NZ_NUM] = { 2, 5, 1, 5, 1, 2, 3, 4, 4, 1, 1, 3 };

  printf ( "\n" );
  printf ( "TEST5722\n" );
  printf ( "  R8SP_PRINT prints a R8SP matrix;\n" );
  printf ( "  In this example, we have listed several matrix\n" );
  printf ( "  locations TWICE.  R8SP_PRINT should compute the\n" );
  printf ( "  sum of these values.\n" );
  printf ( "\n" );
  printf ( "  In particular, we want A(1,2) = 112 and A(3,4) = 234;\n" );
  printf ( "\n" );
  printf ( "  Matrix rows M =          %d\n", M );
  printf ( "  Matrix columns N =       %d\n", N );
  printf ( "  Matrix nonzeros NZ_NUM = %d\n", NZ_NUM );

  r8sp_print ( M, N, NZ_NUM, row, col, a, "  The R8SP matrix:" );

  return;
# undef M
# undef N
# undef NZ_NUM
}
/******************************************************************************/

void test5724 ( )

/******************************************************************************/
/*
  Purpose:

    TEST5724 tests R8SP_WRITE.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    20 April 2013

  Author:

    John Burkardt
*/
{
# define M 100
# define N 100
# define NZ_NUM ( 3 * N - 2 )

  double a[NZ_NUM];
  int col[NZ_NUM];
  int i;
  int j;
  int k;
  char output_file[] = "r8sp_matrix.txt";
  int row[NZ_NUM];

  printf ( "\n" );
  printf ( "TEST5724\n" );
  printf ( "  For a R8SP matrix,\n" );
  printf ( "  R8SP_WRITE writes the matrix to a file.\n" );
  printf ( "\n" );
  printf ( "  Matrix number of rows M =    %d\n", M );
  printf ( "  Matrix number of columns N = %d\n", N );
  printf ( "  Matrix nonzeros NZ_NUM =     %d\n", NZ_NUM );
/*
  Set the matrix values.
*/
  k = 0;
  for ( i = 1; i <= N; i++ )
  {

    j = i;
    row[k] = i;
    col[k] = j;
    a[k] = ( double ) ( 100 * i + j );
    k = k + 1;
  }

  for ( i = 2; i <= N; i++ )
  {
    j = i - 1;
    row[k] = i;
    col[k] = j;
    a[k] = ( double ) ( 100 * i + j );
    k = k + 1;
  }

  for ( i = 1; i <= N-1; i++ )
  {
    j = i + 1;
    row[k] = i;
    col[k] = j;
    a[k] = ( double ) ( 100 * i + j );
    k = k + 1;
  }

  r8sp_print_some ( M, N, NZ_NUM, row, col, a, 1, 1, 
    10, 10, "  Initial 10x10 block of R8SP matrix:" );

  r8sp_write ( M, N, NZ_NUM, row, col, a, output_file );

  printf ( "  R8SP_WRITE wrote the matrix data to \"%s\".\n", output_file );

  return;
# undef M
# undef N
# undef NZ_NUM
}
/******************************************************************************/

void test5725 ( )

/******************************************************************************/
/*
  Purpose:

    TEST5725 tests R8SP_READ, R8SP_READ_SIZE..

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    20 April 2013

  Author:

    John Burkardt
*/
{
  double *a;
  int *col;
  int m;
  int n;
  int nz_num;
  char input_file[] = "r8sp_matrix.txt";
  int *row;

  printf ( "\n" );
  printf ( "TEST5725\n" );
  printf ( "  For a R8SP matrix,\n" );
  printf ( "  R8SP_READ reads a matrix from a file.\n" );
  printf ( "  R8SP_READ_SIZE reads the sizes of the matrix from a file.\n" );

  r8sp_read_size ( input_file, &m, &n, &nz_num );

  printf ( "\n" );
  printf ( "  R8SP_READ_SIZE reports matrix size data:\n" );
  printf ( "\n" );

  printf ( "\n" );
  printf ( "  Matrix number of rows M =    %d\n", m );
  printf ( "  Matrix number of columns N = %d\n",  n );
  printf ( "  Matrix nonzeros NZ_NUM =     %d\n", nz_num );

  row = ( int * ) malloc ( nz_num * sizeof ( int ) );
  col = ( int * ) malloc ( nz_num * sizeof ( int ) );
  a = ( double * ) malloc ( nz_num * sizeof ( double ) );

  r8sp_read ( input_file, m, n, nz_num, row, col, a );

  r8sp_print_some ( m, n, nz_num, row, col, a, 1, 1, 
    10, 10, "  Initial 10x10 block of R8SP matrix:" );

  printf ( "\n" );
  printf ( "  Deleting the matrix data file \"%s\"\n", input_file );

  file_delete ( input_file );

  free ( a );
  free ( col );
  free ( row );

  return;
}
/******************************************************************************/

void test573 ( )

/******************************************************************************/
/*
  Purpose:

    TEST573 tests R8SR_INDICATOR.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    19 April 2013

  Author:

    John Burkardt
*/
{
# define N 5
# define NZ 7

  int col[NZ] = { 2, 5, 5, 1, 1, 2, 3 };
  double diag[N];
  double off[NZ];
  int row[N+1] = { 1, 3, 4, 5, 6, 8 };

  printf ( "\n" );
  printf ( "TEST573\n" );
  printf ( "  R8SR_INDICATOR sets up an R8SR indicator matrix;\n" );
  printf ( "\n" );
  printf ( "  Matrix order N = %d\n", N );

  r8sr_indicator ( N, NZ, row, col, diag, off );

  r8sr_print ( N, NZ, row, col, diag, off, "  The R8SR indicator matrix:" );

  return;
# undef N
# undef NZ
}
/******************************************************************************/

void test574 ( )

/******************************************************************************/
/*
  Purpose:

    TEST574 tests R8SR_MXV, R8SR_VXM.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 April 2013

  Author:

    John Burkardt
*/
{
# define N 5
# define NZ 7

  double *b;
  double *c;
  int col[NZ] = { 2, 5, 5, 1, 1, 2, 3 };
  double diag[N];
  int i;
  double off[NZ];
  int row[N+1] = { 1, 3, 4, 5, 6, 8 };
  int seed = 123456789;
  double x[N];

  printf ( "\n" );
  printf ( "TEST574\n" );
  printf ( "  R8SR_MXV multiplies an R8SR matrix by a vector;\n" );
  printf ( "  R8SR_VXM multiplies a vector by an R8SR matrix;\n" );
  printf ( "\n" );
  printf ( "  Matrix order N = %d\n", N );
/*
  Set the matrix.
*/
  r8sr_random ( N, NZ, row, col, diag, off, &seed );
/*
  Make an R8GE copy.
*/
  c = r8sr_to_r8ge ( N, NZ, row, col, diag, off );
/*
  Print the R8GE copy.
*/
  r8ge_print ( N, N, c, "  The R8SR matrix, in R8GE form:" );

  x[0] = 1.0;
  for ( i = 1; i < N-1; i++ )
  {
    x[i] = 0.0;
  }
  x[N-1] = -1.0;

  r8vec_print ( N, x, "  The vector x:" );

  b = r8sr_mxv ( N, NZ, row, col, diag, off, x );

  r8vec_print ( N, b, "  The product A * x:" );

  free ( b );
  b = r8sr_vxm ( N, NZ, row, col, diag, off, x );

  r8vec_print ( N, b, "  The product A' * x:" );

  free ( b );
  free ( c );

  return;
# undef N
# undef NZ
}
/******************************************************************************/

void test5745 ( )

/******************************************************************************/
/*
  Purpose:

    TEST5745 tests R8SR_PRINT.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 April 2013

  Author:

    John Burkardt
*/
{
# define N 5
# define NZ 7

  int col[NZ] = { 2, 5, 5, 1, 1, 2, 3 };
  double diag[N];
  double off[NZ];
  int row[N+1] = { 1, 3, 4, 5, 6, 8 };
  int seed = 123456789;

  printf ( "\n" );
  printf ( "TEST5745\n" );
  printf ( "  R8SR_PRINT prints an R8SR matrix.\n" );
  printf ( "\n" );
  printf ( "  Matrix order N = %d\n", N );
/*
  Set the matrix.
*/
  r8sr_random ( N, NZ, row, col, diag, off, &seed );
/*
  Print the matrix.
*/
  r8sr_print ( N, NZ, row, col, diag, off, "  The R8SR matrix:" );

  return;
# undef N
# undef NZ
}
/******************************************************************************/

void test575 ( )

/******************************************************************************/
/*
  Purpose:

    TEST575 tests R8SR_RANDOM.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 April 2013

  Author:

    John Burkardt
*/
{
# define N 5
# define NZ 7

  double *b;
  int col[NZ] = { 2, 5, 5, 1, 1, 2, 3 };
  double diag[N];
  double off[NZ];
  int row[N+1] = { 1, 3, 4, 5, 6, 8 };
  int seed = 123456789;

  printf ( "\n" );
  printf ( "TEST575\n" );
  printf ( "  R8SR_RANDOM randomizes an R8SR matrix\n" );
  printf ( "\n" );
  printf ( "  Matrix order N = %d\n", N );
/*
  Set the matrix.
*/
  r8sr_random ( N, NZ, row, col, diag, off, &seed );
/*
  Make an R8GE copy.
*/
  b = r8sr_to_r8ge ( N, NZ, row, col, diag, off );
/*
  Print the R8GE copy.
*/
  r8ge_print ( N, N, b, "  The R8SR matrix, in R8GE form:" );

  free ( b );

  return;
# undef N
# undef NZ
}
/******************************************************************************/

void test577 ( )

/******************************************************************************/
/*
  Purpose:

    TEST577 tests R8SS_INDICATOR.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 April 2013

  Author:

    John Burkardt
*/
{
# define N 9

  double *a;
  int diag[N];
  int na;

  printf ( "\n" );
  printf ( "TEST577\n" );
  printf ( "  For a symmetric skyline storage matrix,\n" );
  printf ( "  R8SS_INDICATOR sets up an indicator matrix.\n" );
  printf ( "\n" );
  printf ( "  Matrix order N = %d\n", N );

  a = r8ss_indicator ( N, &na, diag );

  r8ss_print ( N, na, diag, a, "  The R8SS indicator matrix:" );

  free ( a );

  return;  
# undef N
}
/******************************************************************************/

void test58 ( )

/******************************************************************************/
/*
  Purpose:

    TEST58 tests R8SS_MXV, R8SS_PRINT.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 April 2013

  Author:

    John Burkardt
*/
{
# define N 9

  double a[(N*(N+1))/2];
  double *a2;
  double *b;
  double *b2;
  int diag[N];
  int i;
  int ij;
  int ilo;
  int j;
  int na;
  int seed = 123456789;
  double *x;

  printf ( "\n" );
  printf ( "TEST58\n" );
  printf ( "  For a symmetric skyline storage matrix,\n" );
  printf ( "  R8SS_MXV computes A*x,\n" );
  printf ( "  R8SS_PRINT prints it.\n" );
  printf ( "\n" );
  printf ( "  Matrix order N = %d\n", N );
/*
  Set the matrix.
*/
  r8ss_random ( N, &na, diag, a, &seed );

  printf ( "\n" );
  printf ( "  Number of nonzero entries stored is %d\n", na );
  printf ( "\n" );
  printf ( "  Diagonal storage indices:\n" );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    printf ( "%6d  %g\n", i + 1, diag[i] );
  }
/*
  Replace the random entries by marker values.
*/
  ij = 0;
  for ( j = 1; j <= N; j++ )
  {
    if ( j == 1 )
    {
      ilo = 1;
    }
    else
    {
      ilo = diag[j-2] - diag[j-1] + j + 1;
    }

    for ( i = ilo; i <= j; i++ )
    {
      a[ij] = ( double ) ( 10 * i + j );
      ij = ij + 1;
    }
  }
  r8ss_print ( N, na, diag, a, "  The R8SS matrix:" );
/*
  Copy the matrix into a general matrix.
*/
  a2 = r8ss_to_r8ge ( N, na, diag, a );
/*
  Set the vector X.
*/
  x = r8vec_indicator_new ( N );
/*
  Compute the product.
*/
  b = r8ss_mxv ( N, na, diag, a, x );
/*
  Compute the product using the general matrix.
*/
  b2 = r8ge_mxv ( N, N, a2, x );
/*
  Compare the results.
*/
  r8vec2_print_some ( N, b, b2, 10, "  R8SS_MXV verse R8GE_MXV" );

  free ( a2 );
  free ( b );
  free ( b2 );
  free ( x );

  return;  
# undef N
}
/******************************************************************************/

void test581 ( )

/******************************************************************************/
/*
  Purpose:

    TEST581 tests R8STO_INDICATOR.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 April 2013

  Author:

    John Burkardt
*/
{
# define N 4

  double *a;

  printf ( "\n" );
  printf ( "TEST581\n" );
  printf ( "  R8STO_INDICATOR sets up an R8STO indicator matrix.\n" );
  printf ( "\n" );
  printf ( "  Matrix order N = %d\n", N );

  a = r8sto_indicator ( N );

  r8sto_print ( N, a, "  The R8STO indicator matrix:" );

  free ( a );

  return;
# undef N
}
/******************************************************************************/

void test583 ( )

/******************************************************************************/
/*
  Purpose:

    TEST583 tests R8STO_INVERSE.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 April 2013

  Author:

    John Burkardt
*/
{
# define N 3

  double a[N] = { 4.0, 2.0, 0.8 };
  double *a2;
  double *b;
  double *c;

  printf ( "\n" );
  printf ( "TEST583\n" );
  printf ( "  R8STO_INVERSE computes the inverse of a positive\n" );
  printf ( "  definite symmetric Toeplitz matrix.\n" );
  printf ( "\n" );
  printf ( "  Matrix order N = %d\n", N );

  r8sto_print ( N, a, "  The symmetric Toeplitz matrix A:" );

  b = r8sto_inverse ( N, a );

  r8ge_print ( N, N, b, "  The inverse matrix B:" );

  a2 = r8sto_to_r8ge ( N, a );

  c = r8ge_mxm ( N, a2, b );

  r8ge_print ( N, N, c, "  The product C = A * B:" );

  free ( a2 );
  free ( b );
  free ( c );

  return;
# undef N
}
/******************************************************************************/

void test585 ( )

/******************************************************************************/
/*
  Purpose:

    TEST585 tests R8STO_MXV, R8STO_YW_SL.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

  09 April 2013

  Author:

    John Burkardt
*/
{
# define N 3

  double a[N];
  double *b;
  int i;
  int job;
  double r[N+1] = { 1.0, 0.5, 0.2, 0.1 };
  double *x;

  for ( i = 0; i < N; i++ )
  {
    a[i] = r[i];
  }

  printf ( "\n" );
  printf ( " TEST585\n" );
  printf ( "  R8STO_YW_SL solves the Yule-Walker equations for a\n" );
  printf ( "  symmetric Toeplitz system.\n" );
  printf ( "\n" );
  printf ( "  Matrix order N = %d\n", N );

  r8sto_print ( N, a, "  The symmetric Toeplitz matrix:" );

  b = ( double * ) malloc ( N * sizeof ( double ) );
  for ( i = 0; i < N; i++ )
  {
    b[i] = -r[i+1];
  }

  r8vec_print ( N, b, "  The right hand side, B:" );

  for ( i = 0; i < N; i++ )
  {
    b[i] = -b[i];
  }

  x = r8sto_yw_sl ( N, b );

  r8vec_print ( N, x, "  The computed solution, X:" );

  free ( b );
  b = r8sto_mxv ( N, a, x );

  r8vec_print ( N, b, "  The product A*x:" );

  free ( b );
  free ( x );

  return;
# undef N
}
/******************************************************************************/

void test587 ( )

/******************************************************************************/
/*
  Purpose:

    TEST587 tests R8STO_SL.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    08 April 2013

  Author:

    John Burkardt
*/
{
# define N 3

  double a[N] = { 1.0, 0.5, 0.2 };
  double *ax;
  double b[N] = { 4.0, -1.0, 3.0 };
  double *x;

  printf ( "\n" );
  printf ( "TEST587\n" );
  printf ( "  R8STO_SL solves a positive definite symmetric Toeplitz system.\n" );
  printf ( "\n" );
  printf ( "  Matrix order N = %d\n", N );

  r8sto_print ( N, a, "  The symmetric Toeplitz matrix A:" );

  r8vec_print ( N, b, "  The right hand side vector b:" );

  x = r8sto_sl ( N, a, b );

  r8vec_print ( N, x, "  The computed solution x:" );

  ax = r8sto_mxv ( N, a, x );

  r8vec_print ( N, ax, "  The product vector A * x:" );

  free ( ax );
  free ( x );

  return;
# undef N
}
/******************************************************************************/

void test589 ( )

/******************************************************************************/
/*
  Purpose:

    TEST589 tests R8TO_INDICATOR.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    07 April 2013

  Author:

    John Burkardt
*/
{
# define N 10

  double *a;

  printf ( "\n" );
  printf ( "TEST589\n" );
  printf ( "  R8TO_INDICATOR sets up an R8TO indicator matrix.\n" );
  printf ( "\n" );
  printf ( "  Matrix order N = %d\n", N );

  a = r8to_indicator ( N );

  r8to_print ( N, a, "  The R8TO indicator matrix:" );

  free ( a );

  return;
# undef N
}
/******************************************************************************/

void test59 ( )

/******************************************************************************/
/*
  Purpose:

    TEST59 tests R8TO_SL.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 April 2013

  Author:

    John Burkardt
*/
{
# define N 10

  double *a;
  double *b;
  int job;
  int seed = 123456789;
  double *x;

  printf ( "\n" );
  printf ( "TEST59\n" );
  printf ( "  R8TO_SL solves a Toeplitz system.\n" );
  printf ( "\n" );
  printf ( "  Matrix order N = %d\n", N );
/*
  Set the matrix.
*/
  a = r8to_random ( N, &seed );

  r8to_print ( N, a, "  The Toeplitz matrix:" );

  for ( job = 0; job <= 1; job++ )
  {
//
//  Set the desired solution.
//
    x = r8vec_indicator_new ( N );
//
//  Compute the corresponding right hand side.
//
    if ( job == 0 )
    {
      b = r8to_mxv ( N, a, x );
    }
    else
    {
      b = r8to_vxm ( N, a, x );
    }
//
//  Solve the linear system.
//
    free ( x );
    x = r8to_sl ( N, a, b, job );

    if ( job == 0 )
    {
      r8vec_print_some ( N, x, 1, 10, "  Solution:" );
    }
    else
    {
      r8vec_print_some ( N, x, 1, 10, "  Solution to transposed system:" );
    }
    free ( b );
    free ( x );
  }
  free ( a );

  return;
# undef N
}
/******************************************************************************/

void test60 ( )

/******************************************************************************/
/*
  Purpose:

    TEST60 tests R8UT_INVERSE, R8UT_DET, R8UT_MXM, R8UT_RANDOM.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    03 April 2013

  Author:

    John Burkardt
*/
{
# define N 5

  double *a;
  double *b;
  double *c;
  double det;
  int i;
  int j;
  int seed = 123456789;

  printf ( "\n" );
  printf ( "TEST60\n" );
  printf ( "  For an upper triangular matrix,\n" );
  printf ( "  R8UT_INVERSE computes the inverse.\n" );
  printf ( "  R8UT_DET computes the determinant.\n" );
  printf ( "  R8UT_MXM computes matrix products.\n" );
  printf ( "\n" );
  printf ( "  Matrix order N = %d\n", N );

  a = r8ut_random ( N, N, &seed );

  r8ut_print ( N, N, a, "  The matrix A:" );
/*
  Compute the determinant.
*/
  det = r8ut_det ( N, a );

  printf ( "\n" );
  printf ( "  Determinant is %g\n", det );
/*
  Compute the inverse matrix B.
*/
  b = r8ut_inverse ( N, a );

  r8ut_print ( N, N, b, "  The inverse matrix B:" );
/*
  Check
*/
  c = r8ut_mxm ( N, a, b );

  r8ut_print ( N, N, c, "  The product A * B:" );

  free ( a );
  free ( b );
  free ( c );

  return;
# undef N
}
/******************************************************************************/

void test605 ( )

/******************************************************************************/
/*
  Purpose:

    TEST605 tests R8UT_INDICATOR.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    02 April 2013

  Author:

    John Burkardt
*/
{
# define M 8
# define N 5

  double *a;

  printf ( "\n" );
  printf ( "TEST605\n" );
  printf ( "  For an upper triangular matrix,\n" );
  printf ( "  R8UT_INDICATOR sets up an indicator matrix.\n" );
  printf ( "\n" );
  printf ( "  Matrix rows M =    %d\n", M );
  printf ( "  Matrix columns N = %d\n", N );

  a = r8ut_indicator ( M, N );

  r8ut_print ( M, N, a, "  The R8UT indicator matrix:" );

  free ( a );

  return;
# undef M
# undef N
}
/******************************************************************************/

void test61 ( )

/******************************************************************************/
/*
  Purpose:

    TEST61 tests R8UT_SL, R8UT_MXV, R8UT_VXM.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    01 April 2013

  Author:

    John Burkardt
*/
{
# define N 10

  double a[N*N];
  double *b;
  int i;
  int j;
  int job;
  double *x;

  printf ( "\n" );
  printf ( "TEST61\n" );
  printf ( "  For an upper triangular matrix,\n" );
  printf ( "  R8UT_SL solves systems;\n" );
  printf ( "  R8UT_MXV computes matrix-vector products.\n" );
  printf ( "  R8UT_VXM computes vector-matrix products.\n" );
  printf ( "\n" );
  printf ( "  Matrix order N = %d \n", N );

  for ( i = 1; i <= N; i++ )
  {
    for ( j = 1; j <= N; j++ )
    {
      if ( i <= j )
      {
        a[i-1+(j-1)*N] = ( double ) j;
      }
      else
      {
        a[i-1+(j-1)*N] = 0.0;
      }
    }
  }

  r8ut_print ( N, N, a, "  The upper triangular matrix:" );

  for ( job = 0; job <= 1; job++ )
  {
/*
  Set the desired solution.
*/
    x = r8vec_indicator_new ( N );
/*
  Compute the corresponding right hand side.
*/
    if ( job == 0 )
    {
      b = r8ut_mxv ( N, N, a, x );
    }
    else
    {
      b = r8ut_vxm ( N, N, a, x );
    }
/*
  Solve the linear system.
*/
    free ( x );
    x = r8ut_sl ( N, a, b, job );

    if ( job == 0 )
    {
      r8vec_print ( N, x, "  Solution:" );
    }
    else
    {
      r8vec_print ( N, x, "  Solution to transposed system:" );
    }
    free ( b );
    free ( x );
  }

  return;
# undef N
}
/******************************************************************************/

void test62 ( )

/******************************************************************************/
/*
  Purpose:

    TEST62 tests R8VM_DET.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    31 March 2013

  Author:

    John Burkardt
*/
{
# define N 10

  double *a;
  double *a2;
  double det;
  int info;
  int pivot[N];
  int seed = 123456789;

  printf ( "\n" );
  printf ( "TEST62\n" );
  printf ( "  R8VM_DET, determinant of a Vandermonde matrix.\n" );
  printf ( "\n" );
  printf ( "  Matrix order N = %d\n", N );
/*
  Set the matrix.
*/
  a = r8vm_random ( N, N, &seed );

  r8vm_print ( N, N, a, "  The Vandermonde matrix:" );
/*
  Copy the matrix into a general array.
*/
  a2 = r8vm_to_r8ge ( N, N, a );
/*
  Compute the determinant.
*/
  det = r8vm_det ( N, a );

  printf ( "\n" );
  printf ( "  R8VM_DET computes the determinant = %g\n", det );
/*
  Factor the general matrix.
*/
  info = r8ge_fa ( N, a2, pivot );
/*
  Compute the determinant.
*/
  det = r8ge_det ( N, a2, pivot );

  printf ( "  R8GE_DET computes the determinant = %g\n", det );

  free ( a );
  free ( a2 );

  return;
# undef N
}
/******************************************************************************/

void test63 ( )

/******************************************************************************/
/*
  Purpose:

    TEST63 tests R8VM_SL_NEW.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    08 March 2013

  Author:

    John Burkardt
*/
{
# define N 5

  double *a;
  double *b;
  int i;
  int info;
  int job;
  int seed = 123456789;
  double *x;

  printf ( "\n" );
  printf ( "TEST63\n" );
  printf ( "  R8VM_SL_NEW solves a Vandermonde system.\n" );
  printf ( "\n" );
  printf ( "  Matrix order N = %d\n", N );
/*
  Set the matrix.
*/
  a = r8vm_random ( N, N, &seed );

  for ( job = 0; job <= 1; job++ )
  {
/*
  Set the desired solution.
*/
    x = r8vec_indicator_new ( N );
/*
  Compute the corresponding right hand side.
*/
    if ( job == 0 )
    {
      b = r8vm_mxv ( N, N, a, x );
    }
    else
    {
      b = r8vm_vxm ( N, N, a, x );
    }
/*
  Solve the linear system.
*/
    free ( x );
    x = r8vm_sl_new ( N, a, b, job, &info );

    if ( job == 0 )
    {
      r8vec_print_some ( N, x, 1, 10, "  Solution:" );
    }
    else
    {
      r8vec_print_some ( N, x, 1, 10, "  Solution to transposed system:" );
    }

    free ( b );
    free ( x );
  }

  free ( a );

  return;
# undef N
}
