# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "linplus.h"

int main ( void );
void test02 ( void );
void test03 ( void );
void test345 ( void );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    LINPLUS_PRB tests routines from the LINPLUS library.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    28 November 2011

  Author:

    John Burkardt
*/
{
  timestamp ( );

  printf ( "\n" );
  printf ( "LINPLUS_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test routines in the LINPLUS library.\n" );

  test02 ( );
  test03 ( );
  test345 ( );
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