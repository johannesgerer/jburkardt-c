# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "cyclic_reduction.h"

int main ( void );
void test01 ( void );
void test02 ( void );
void test03 ( void );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for CYCLIC_REDUCTION_PRB.

  Discussion:

    CYCLIC_REDUCTION_PRB tests the CYCLIC_REDUCTION library.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    10 May 2010

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "CYCLIC_REDUCTION_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the CYCLIC_REDUCTION library.\n" );

  test02 ( );
  test03 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "CYCLIC_REDUCTION_PRB\n" );
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

    TEST02 tests R83_CR_FA, R83_CR_SLS.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    10 May 2010

  Author:

    John Burkardt
*/
{
  double *a;
  double *a_cr;
  double *b;
  int debug = 1;
  int i;
  int j;
  int n = 5;
  int nb = 2;
  double *x;

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  R83_CR_FA factors a real tridiagonal matrix;\n" );
  printf ( "  R83_CR_SLS solves 1 or more systems.\n" );
  printf ( "\n" );
  printf ( "  Matrix order N = %d\n", n );
  printf ( "  Demonstrate multiple system solution method.\n" );
/*
  Set the matrix.
*/
  a = ( double * ) malloc ( 3 * n * sizeof ( double ) );

  a[0+0*3] = 0.0;
  for ( j = 1; j < n; j++ )
  {
    a[0+j*3] = - 1.0;
  }
  for ( j = 0; j < n; j++ )
  {
    a[1+j*3] = 2.0;
  }
  for ( j = 0; j < n - 1; j++ )
  {
    a[2+j*3] = - 1.0;
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
  Solve 2 systems simultaneously.
*/
  b = ( double * ) malloc ( n * nb * sizeof ( double ) );

  for ( i = 0; i < n - 1; i++ )
  {
    b[i+0*n] = 0.0;
  }
  b[n-1+0*n] = ( double ) ( n + 1 );

  b[0+1*n] = 1.0;
  for ( i = 1; i < n - 1; i++ )
  {
    b[i+1*n] = 0.0;
  }
  b[n-1+1*n] = 1.0;
/*
  Solve the linear systems.
*/
  x = r83_cr_sls ( n, a_cr, nb, b );

  r8mat_print_some ( n, nb, x, 1, 1, 10, nb, "  Solutions:" );

  free ( a );
  free ( a_cr );
  free ( b );
  free ( x );

  return;
}
/******************************************************************************/

void test03 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST03 tests R83_CR_FA, R83_CR_SL.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    07 May 2010

  Author:

    John Burkardt
*/
{
  double *a;
  double *a_cr;
  double *b;
  int debug = 0;
  int i;
  int j;
  int n = 10;
  double *x;

  printf ( "\n" );
  printf ( "TEST03\n" );
  printf ( "  For a real tridiagonal matrix,\n" );
  printf ( "  using CYCLIC REDUCTION,\n" );
  printf ( "  R83_CR_FA factors;\n" );
  printf ( "  R83_CR_SL solves a system.\n" );
  printf ( "\n" );
  printf ( "  Matrix order N = %d\n", n );
  printf ( "  The matrix is NOT symmetric.\n" );
/*
  Set the matrix values.
*/
  a = ( double * ) malloc ( 3 * n * sizeof ( double ) );

  a[0+0*3] = 0.0;
  for ( j = 2; j <= n; j++ )
  {
    a[0+(j-1)*3] = ( double ) ( j );
  }
  for ( j = 1; j <= n; j++ )
  {
    a[1+(j-1)*3] = 4.0 * ( double ) ( j );
  }
  for ( j = 1; j <= n - 1; j++ )
  {
    a[2+(j-1)*3] = ( double ) ( j );
  }
  a[2+(n-1)*3] = 0.0;

  if ( debug )
  {
    r83_print ( n, a, "  The matrix:" );
  }
/*
  Set the desired solution.
*/
  x = r8vec_indicator_new ( n );
/*
  Compute the corresponding right hand side.
*/
  b = r83_mxv_new ( n, a, x );

  if ( debug )
  {
    r8vec_print  ( n, b, "  The right hand side:" );
  }
  free ( x );
/*
  Factor the matrix.
*/
  a_cr = r83_cr_fa ( n, a );

  if ( debug )
  {
    r83_print ( 2 * n + 1, a_cr, "  The factor information:" );
  }
/*
  Solve the linear system.
*/
  x = r83_cr_sl ( n, a_cr, b );

  r8vec_print  ( n, x, "  The solution:" );

  free ( a );
  free ( a_cr );
  free ( b );
  free ( x );

  return;
}
