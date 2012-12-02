# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>

# include "mgmres.h"

int main ( void );
void test01 ( void );
void test02 ( void );
void test03 ( void );
void test04 ( void );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    MAIN runs the quick checks for the MGMRES code.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 July 2007

  Author:

    John Burkardt
*/
{
  timestamp ( );

  printf ( "\n" );
  printf ( "MGMRES_PRB:\n" );
  printf ( "  C version\n" );
  printf ( "  Test the MGMRES library.\n" );

  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "MGMRES_PRB:\n" );
  printf ( "  Normal end of execution.\n" );

  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void test01 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST01 tests MGMRES_ST on the simple -1,2-1 matrix.

  Discussion:

    This is a very weak test, since the matrix has such a simple
    structure, is diagonally dominant (though not strictly), 
    and is symmetric.

    To make the matrix bigger, simply increase the value of N.

    Note that MGMRES_ST expects the matrix to be stored using the
    sparse triplet storage format.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 July 2007

  Author:

    John Burkardt
*/
{
# define N 20
# define NZ_NUM 3 * N - 2

  double a[NZ_NUM];
  int i;
  int ia[NZ_NUM];
  int itr_max;
  int j;
  int ja[NZ_NUM];
  int k;
  int mr;
  int n = N;
  int nz_num = NZ_NUM;
  double rhs[N];
  int test;
  double tol_abs;
  double tol_rel;
  double x_error;
  double x_estimate[N];
  double x_exact[N];

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  Test MGMRES_ST on the simple -1,2-1 matrix.\n" );
/*
  Set the matrix.
  Note that we use zero based index valuesin IA and JA.
*/
  k = 0;

  for ( i = 0; i < n; i++ )
  {
    if ( 0 < i )
    {
      ia[k] = i;
      ja[k] = i-1;
      a[k] = -1.0;
      k = k + 1;
    }

    ia[k] = i;
    ja[k] = i;
    a[k] = 2.0;
    k = k + 1;

    if ( i < n - 1 )
    {
      ia[k] = i;
      ja[k] = i+1;
      a[k] = -1.0;
      k = k + 1;
    }

  }
/*
  Set the right hand side:
*/
  for ( i = 0; i < n-1; i++ )
  {
    rhs[i] = 0.0;
  }
  rhs[n-1] = ( double ) ( n + 1 );
/*
  Set the exact solution.
*/
  for ( i = 0; i < n; i++ )
  {
    x_exact[i] = ( double ) ( i + 1 );
  }

  for ( test = 1; test <= 3; test++ )
  {
/*
  Set the initial solution estimate.
*/
    for ( i = 0; i < n; i++ )
    {
      x_estimate[i] = 0.0;
    }
    x_error = 0.0;
    for ( i = 0; i < n; i++ )
    {
      x_error = x_error + pow ( x_exact[i] - x_estimate[i], 2 );
    }
    x_error = sqrt ( x_error );

    if ( test == 1 )
    {
      itr_max = 1;
      mr = 20;
    }
    else if ( test == 2 )
    {
      itr_max = 2;
      mr = 10;
    }
    else if ( test == 3 )
    {
      itr_max = 5;
      mr = 4;
    }
    tol_abs = 1.0E-08;
    tol_rel = 1.0E-08;

    printf ( "\n" );
    printf ( "  Test %d\n", test );
    printf ( "  Matrix order N = %d\n", n );
    printf ( "  Inner iteration limit = %d\n", mr );
    printf ( "  Outer iteration limit = %d\n", itr_max );
    printf ( "  Initial X_ERROR = %g\n", x_error );

    mgmres_st ( n, nz_num, ia, ja, a, x_estimate, rhs, itr_max, mr, tol_abs,
      tol_rel );

    x_error = 0.0;
    for ( i = 0; i < n; i++ )
    {
      x_error = x_error + pow ( x_exact[i] - x_estimate[i], 2 );
    }
    x_error = sqrt ( x_error );

    printf ( "  Final X_ERROR = %g\n", x_error );
  }

  return;
# undef N
# undef NZ_NUM
}
/******************************************************************************/

void test02 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST02 tests MGMRES_ST on a 9 by 9 matrix.

  Discussion:

    Note that MGMRES_ST expects the matrix to be stored using the
    sparse triplet storage format.

    A = 
      2  0  0 -1  0  0  0  0  0
      0  2 -1  0  0  0  0  0  0
      0 -1  2  0  0  0  0  0  0
     -1  0  0  2 -1  0  0  0  0
      0  0  0 -1  2 -1  0  0  0
      0  0  0  0 -1  2 -1  0  0
      0  0  0  0  0 -1  2 -1  0
      0  0  0  0  0  0 -1  2 -1
      0  0  0  0  0  0  0 -1  2

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 July 2007

  Author:

    John Burkardt
*/
{
# define N 9
# define NZ_NUM 23

  double a[NZ_NUM] = {
    2.0, -1.0,
    2.0, -1.0,
    -1.0, 2.0,
    -1.0, 2.0, -1.0,
    -1.0, 2.0, -1.0,
    -1.0, 2.0, -1.0,
    -1.0, 2.0, -1.0,
    -1.0, 2.0, -1.0,
    -1.0, 2.0 };
  int i;
  int ia[NZ_NUM] = {
    0, 0,
    1, 1,
    2, 2,
    3, 3, 3,
    4, 4, 4,
    5, 5, 5,
    6, 6, 6,
    7, 7, 7,
    8, 8 };
  int itr_max;
  int j;
  int ja[NZ_NUM] = {
    0, 3,
    1, 2,
    1, 2,
    0, 3, 4,
    3, 4, 5,
    4, 5, 6,
    5, 6, 7,
    6, 7, 8,
    7, 8 };
  int k;
  int mr;
  int n = N;
  int nz_num = NZ_NUM;
  double rhs[N] = {
    1.0,
    1.0,
    1.0,
    1.0,
    1.0,
    1.0,
    1.0,
    1.0,
    1.0 };
  int seed = 123456789;
  int test;
  double tol_abs;
  double tol_rel;
  double x_error;
  double *x_estimate;
  double x_exact[N] = {
    3.5,
    1.0,
    1.0,
    6.0,
    7.5,
    8.0,
    7.5,
    6.0,
    3.5 };

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  Test MGMRES_ST on matrix that is not quite the -1,2,-1 matrix,\n" );
  printf ( "  of order N = %d\n", n );

  for ( test = 1; test <= 2; test++ )
  {
    if ( test == 1 )
    {
      printf ( "\n" );
      printf ( "  First try, use zero initial vector:\n" );
      printf ( "\n" );

      x_estimate = ( double * ) malloc ( ( size_t ) ( n * sizeof ( double ) ) );

      for ( i = 0; i < n; i++ )
      {
        x_estimate[i] = 0.0;
      }
    }
    else
    {
      printf ( "\n" );
      printf ( "  Second try, use random initial vector:\n" );
      printf ( "\n" );

      x_estimate = r8vec_uniform_01 ( n, &seed );
    }
/*
  Set the initial solution estimate.
*/
    x_error = 0.0;
    for ( i = 0; i < n; i++ )
    {
      x_error = x_error + pow ( x_exact[i] - x_estimate[i], 2 );
    }
    x_error = sqrt ( x_error );

    printf ( "  Before calling the solver, X_ERROR = %g\n", x_error );

    itr_max = 20;
    mr = n - 1;
    tol_abs = 1.0E-08;
    tol_rel = 1.0E-08;

    mgmres_st ( n, nz_num, ia, ja, a, x_estimate, rhs, itr_max, mr, tol_abs,
      tol_rel );

    x_error = 0.0;
    for ( i = 0; i < n; i++ )
    {
      x_error = x_error + pow ( x_exact[i] - x_estimate[i], 2 );
    }
    x_error = sqrt ( x_error );

    printf ( "  After calling the solver, X_ERROR = %g\n", x_error );

    printf ( "\n" );
    printf ( "  Final solution estimate:\n" );
    printf ( "\n" );
    for ( i = 0; i < n; i++ )
    {
      printf ( "  %8d  %12f\n", i, x_estimate[i] );
    }

    free ( ( char * ) x_estimate );
  }

  return;
# undef N
# undef NZ_NUM
}
/******************************************************************************/

void test03 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST03 tests PMGMRES_ILU_CR on the simple -1,2-1 matrix.

  Discussion:

    This is a very weak test, since the matrix has such a simple
    structure, is diagonally dominant (though not strictly), 
    and is symmetric.

    To make the matrix bigger, simply increase the value of N.

    Note that PGMRES_ILU_CR expects the matrix to be stored using the
    sparse compressed row format.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 July 2007

  Author:

    John Burkardt
*/
{
# define N 20
# define NZ_NUM 3 * N - 2

  double a[NZ_NUM];
  int i;
  int ia[N+1];
  int itr_max;
  int j;
  int ja[NZ_NUM];
  int k;
  int mr;
  int n = N;
  int nz_num = NZ_NUM;
  double rhs[N];
  int test;
  double tol_abs;
  double tol_rel;
  double x_error;
  double x_estimate[N];
  double x_exact[N];

  printf ( "\n" );
  printf ( "TEST03\n" );
  printf ( "  Test PMGMRES_ILU_CR on the simple -1,2-1 matrix.\n" );
/*
  Set the matrix.
  Note that we use zero based index valuesin IA and JA.
*/
  k = 0;
  ia[0] = 0;

  printf ( "\n" );
  printf ( "  ia[%d] = %d\n", 0, ia[0] );
  for ( i = 0; i < n; i++ )
  {
    ia[i+1] = ia[i];
    if ( 0 < i )
    {
      ia[i+1] = ia[i+1] + 1;
      ja[k] = i-1;
      a[k] = -1.0;
      k = k + 1;
    }

    ia[i+1] = ia[i+1] + 1;
    ja[k] = i;
    a[k] = 2.0;
    k = k + 1;

    if ( i < N-1 )
    {
      ia[i+1] = ia[i+1] + 1;
      ja[k] = i+1;
      a[k] = -1.0;
      k = k + 1;
    }
    printf ( "  ia[%d] = %d\n", i+1, ia[i+1] );
  }
/*
  Set the right hand side:
*/
  for ( i = 0; i < n-1; i++ )
  {
    rhs[i] = 0.0;
  }
  rhs[n-1] = ( double ) ( n + 1 );
/*
  Set the exact solution.
*/
  for ( i = 0; i < n; i++ )
  {
    x_exact[i] = ( double ) ( i + 1 );
  }

  for ( test = 1; test <= 3; test++ )
  {
/*
  Set the initial solution estimate.
*/
    for ( i = 0; i < n; i++ )
    {
      x_estimate[i] = 0.0;
    }
    x_error = 0.0;
    for ( i = 0; i < n; i++ )
    {
      x_error = x_error + pow ( x_exact[i] - x_estimate[i], 2 );
    }
    x_error = sqrt ( x_error );

    if ( test == 1 )
    {
      itr_max = 1;
      mr = 20;
    }
    else if ( test == 2 )
    {
      itr_max = 2;
      mr = 10;
    }
    else if ( test == 3 )
    {
      itr_max = 5;
      mr = 4;
    }
    tol_abs = 1.0E-08;
    tol_rel = 1.0E-08;

    printf ( "\n" );
    printf ( "  Test %d\n", test );
    printf ( "  Matrix order N = %d\n", n );
    printf ( "  Inner iteration limit = %d\n", mr );
    printf ( "  Outer iteration limit = %d\n", itr_max );
    printf ( "  Initial X_ERROR = %g\n", x_error );

    pmgmres_ilu_cr ( n, nz_num, ia, ja, a, x_estimate, rhs, itr_max, mr, 
      tol_abs, tol_rel );

    x_error = 0.0;
    for ( i = 0; i < n; i++ )
    {
      x_error = x_error + pow ( x_exact[i] - x_estimate[i], 2 );
    }
    x_error = sqrt ( x_error );

    printf ( "  Final X_ERROR = %g\n", x_error );
  }

  return;
# undef N
# undef NZ_NUM
}
/******************************************************************************/

void test04 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST04 tests PMGMRES_ILU_CR on a simple 5 by 5 matrix.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    25 July 2007

  Author:

    John Burkardt
*/
{
# define N 5
# define NZ_NUM 9

  double a[NZ_NUM] = { 
     1.0, 2.0, 1.0, 
     2.0,
     3.0, 3.0,
     4.0, 
     1.0, 5.0 };
  int i;
  int ia[N+1] = { 0, 3, 4, 6, 7, 9 };
  int itr_max;
  int j;
  int ja[NZ_NUM] = {
    0, 3, 4, 
    1, 
    0, 2, 
    3, 
    1, 4 };
  int k;
  int mr;
  int n = N;
  int nz_num = NZ_NUM;
  double rhs[N] = { 14.0, 4.0, 12.0, 16.0, 27.0 };
  int test;
  double tol_abs;
  double tol_rel;
  double x_error;
  double x_estimate[N];
  double x_exact[N] = { 1.0, 2.0, 3.0, 4.0, 5.0 };

  printf ( "\n" );
  printf ( "TEST04\n" );
  printf ( "  Test PMGMRES_ILU_CR on a simple 5 x 5 matrix.\n" );

  printf ( "\n" );
  for ( i = 0; i <= n; i++ )
  {
    printf ( "  ia[%d] = %d\n", i, ia[i] );
  }
  for ( test = 1; test <= 3; test++ )
  {
/*
  Set the initial solution estimate.
*/
    for ( i = 0; i < n; i++ )
    {
      x_estimate[i] = 0.0;
    }
    x_error = 0.0;
    for ( i = 0; i < n; i++ )
    {
      x_error = x_error + pow ( x_exact[i] - x_estimate[i], 2 );
    }
    x_error = sqrt ( x_error );

    if ( test == 1 )
    {
      itr_max = 1;
      mr = 20;
    }
    else if ( test == 2 )
    {
      itr_max = 2;
      mr = 10;
    }
    else if ( test == 3 )
    {
      itr_max = 5;
      mr = 4;
    }
    tol_abs = 1.0E-08;
    tol_rel = 1.0E-08;

    printf ( "\n" );
    printf ( "  Test %d\n", test );
    printf ( "  Matrix order N = %d\n", n );
    printf ( "  Inner iteration limit = %d\n", mr );
    printf ( "  Outer iteration limit = %d\n", itr_max );
    printf ( "  Initial X_ERROR = %g\n", x_error );

    pmgmres_ilu_cr ( n, nz_num, ia, ja, a, x_estimate, rhs, itr_max, mr, 
      tol_abs, tol_rel );

    x_error = 0.0;
    for ( i = 0; i < n; i++ )
    {
      x_error = x_error + pow ( x_exact[i] - x_estimate[i], 2 );
    }
    x_error = sqrt ( x_error );

    printf ( "  Final X_ERROR = %g\n", x_error );
  }

  return;
# undef N
# undef NZ_NUM
}
