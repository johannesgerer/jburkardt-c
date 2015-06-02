# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "laplacian.h"

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );
void test04 ( );
void test05 ( );
void test06 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for LAPLACIAN_PRB.

  Discussion:

    LAPLACIAN_PRB tests the LAPLACIAN library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 November 2013

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "LAPLACIAN_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the LAPLACIAN library.\n" );

  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
  test05 ( );
  test06 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "LAPLACIAN_PRB\n" );
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

    TEST01 tests L1DD and similar routines.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 October 2013

  Author:

    John Burkardt
*/
{
  double h;
  double *l;
  int n = 5;
  int test;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  A full-storage matrix is returned by:\n" );
  printf ( "  L1DD: Dirichlet/Dirichlet BC;\n" );
  printf ( "  L1DN: Dirichlet/Neumann BC;\n" );
  printf ( "  L1ND: Neumann/Dirichlet BC;\n" );
  printf ( "  L1NN: Neumann/Neumann BC;\n" );
  printf ( "  L1PP: Periodic BC;\n" );

  for ( test = 1; test <= 2; test++ )
  {
    if ( test == 1 )
    {
      h = 1.0;
    }
    else
    {
      h = 1.0 / ( double ) ( n + 1 );
    }

    printf ( "\n" );
    printf ( "  Using spacing H = %g\n", h );

    l = l1dd ( n, h );
    r8mat_print ( n, n, l, "  L1DD:" );
    free ( l );

    l = l1dn ( n, h );
    r8mat_print ( n, n, l, "  L1DN:" );
    free ( l );

    l = l1nd ( n, h );
    r8mat_print ( n, n, l, "  L1ND:" );
    free ( l );

    l = l1nn ( n, h );
    r8mat_print ( n, n, l, "  L1NN:" );
    free ( l );

    l = l1pp ( n, h );
    r8mat_print ( n, n, l, "  L1PP:" );
    free ( l );
  }

  return;
}
/******************************************************************************/

void test02 ( )

/******************************************************************************/
/*
  Purpose:

    TEST02 tests L1DD_APPLY and similar functions.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 October 2013

  Author:

    John Burkardt
*/
{
  double h;
  int i;
  double *lu;
  int n = 9;
  double *u;
  double *x;

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  The Laplacian L is applied to data U by:\n" );
  printf ( "  L1DD_APPLY for Dirichlet/Dirichlet BC;\n" );
  printf ( "  L1DN_APPLY for Dirichlet/Neumann BC;\n" );
  printf ( "  L1ND_APPLY for Neumann/Dirichlet BC;\n" );
  printf ( "  L1NN_APPLY for Neumann/Neumann BC;\n" );
  printf ( "  L1PP_APPLY for Periodic BC;\n" );

  x = ( double * ) malloc ( n * sizeof ( double ) );
  for ( i = 0; i < n; i++ )
  {
    x[i] = ( double ) ( i + 1 ) / ( double ) ( n + 1 );
  }
  h = 1.0 / ( double ) ( n + 1 );

  printf ( "\n" );
  printf ( " Using spacing H = %g\n", h );

  u = ( double * ) malloc ( n * sizeof ( double ) );
  for ( i = 0; i < n; i++ )
  {
    u[i] = x[i] * ( 1.0 - x[i] );
  }

  r8vec_print ( n, u, "  Vector U:" );

  lu = l1dd_apply ( n, h, u );
  r8vec_print ( n, lu, "  L1DD(U):" );
  free ( lu );

  lu = l1dn_apply ( n, h, u );
  r8vec_print ( n, lu, "  L1DN(U):" );
  free ( lu );

  lu = l1nd_apply ( n, h, u );
  r8vec_print ( n, lu, "  L1ND(U):" );
  free ( lu );

  lu = l1nn_apply ( n, h, u );
  r8vec_print ( n, lu, "  L1NN(U):" );
  free ( lu );

  lu = l1pp_apply ( n, h, u );
  r8vec_print ( n, lu, "  L1PP(U):" );
  free ( lu );
  
  free ( u );
  free ( x );

  return;
}
/******************************************************************************/

void test03 ( )

/******************************************************************************/
/*
  Purpose:

    TEST03 tests L1DD_EIGEN and similar routines.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 October 2013

  Author:

    John Burkardt
*/
{
  double *a;
  double err;
  double h;
  int i;
  double *lambda;
  int n = 5;
  int test;
  double *v;

  printf ( "\n" );
  printf ( "TEST03\n" );
  printf ( "  Compute eigen information for the Laplacian:\n" );
  printf ( "  L1DD_EIGEN for Dirichlet/Dirichlet BC;\n" );
  printf ( "  L1DN_EIGEN for Dirichlet/Neumann BC;\n" );
  printf ( "  L1ND_EIGEN for Neumann/Dirichlet BC;\n" );
  printf ( "  L1NN_EIGEN for Neumann/Neumann BC;\n" );
  printf ( "  L1PP_EIGEN for Periodic BC;\n" );

  v = ( double * ) malloc ( n * n * sizeof ( double ) );
  lambda = ( double * ) malloc ( n * sizeof ( double ) );

  for ( test = 1; test <= 2; test++ )
  {
    if ( test == 1 )
    {
      h = 1.0;
    }
    else
    {
      h = 1.0 / ( double ) ( n + 1 );
    }

    printf ( "\n" );
    printf ( "  Using spacing H = %g\n", h );

    a = l1dd ( n, h );
    l1dd_eigen ( n, h, v, lambda );
    r8vec_print ( n, lambda, "  L1DD Eigenvalues:" );
    r8mat_print ( n, n, v, "  L1DD Eigenvectors:" );
    err = eigen_error ( n, n, a, v, lambda );
    printf ( "\n" );
    printf ( "  L1DD eigenerror = %g\n", err );
    free ( a );

    a = l1dn ( n, h );
    l1dn_eigen ( n, h, v, lambda );
    r8vec_print ( n, lambda, "  L1DN Eigenvalues:" );
    r8mat_print ( n, n, v, "  L1DN Eigenvectors:" );
    err = eigen_error ( n, n, a, v, lambda );
    printf ( "\n" );
    printf ( "  L1DN eigenerror = %g\n", err );
    free ( a );

    a = l1nd ( n, h );
    l1nd_eigen ( n, h, v, lambda );
    r8vec_print ( n, lambda, "  L1ND Eigenvalues:" );
    r8mat_print ( n, n, v, "  L1ND Eigenvectors:" );
    err = eigen_error ( n, n, a, v, lambda );
    printf ( "\n" );
    printf ( "  L1ND eigenerror = %g\n", err );
    free ( a );

    a = l1nn ( n, h );
    l1nn_eigen ( n, h, v, lambda );
    r8vec_print ( n, lambda, "  L1NN Eigenvalues:" );
    r8mat_print ( n, n, v, "  L1NN Eigenvectors:" );
    err = eigen_error ( n, n, a, v, lambda );
    printf ( "\n" );
    printf ( "  L1NN eigenerror = %g\n", err );
    free ( a );

    a = l1pp ( n, h );
    l1pp_eigen ( n, h, v, lambda );
    r8vec_print ( n, lambda, "  L1PP Eigenvalues:" );
    r8mat_print ( n, n, v, "  L1PP Eigenvectors:" );
    err = eigen_error ( n, n, a, v, lambda );
    printf ( "\n" );
    printf ( "  L1PP eigenerror = %g\n", err );
    free ( a );
  }

  free ( lambda );
  free ( v );

  return;
}
/******************************************************************************/

void test04 ( )

/******************************************************************************/
/*
  Purpose:

    TEST04 tests L1DD_INVERSE and similar routines.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 October 2013

  Author:

    John Burkardt
*/
{
  double err;
  double h;
  double *l;
  double *linv;
  int n = 5;
  int test;

  printf ( "\n" );
  printf ( "TEST04\n" );
  printf ( "  The inverse of a full-storage matrix is returned by:\n" );
  printf ( "  L1DD_INVERSE: Dirichlet/Dirichlet BC;\n" );
  printf ( "  L1DN_INVERSE: Dirichlet/Neumann BC;\n" );
  printf ( "  L1ND_INVERSE: Neumann/Dirichlet BC;\n" );

  for ( test = 1; test <= 2; test++ )
  {
    if ( test == 1 )
    {
      h = 1.0;
    }
    else
    {
      h = 1.0 / ( double ) ( n + 1 );
    }

    printf ( "\n" );
    printf ( "  Using spacing H = %g\n", h );

    l = l1dd ( n, h );
    r8mat_print ( n, n, l, "  L1DD:" );
    linv = l1dd_inverse ( n, h );
    r8mat_print ( n, n, linv, "  L1DD_INVERSE:" );
    err = inverse_error ( n, l, linv );
    printf ( "\n" );
    printf ( "  L1DD inverse error = %g\n", err );
    free ( l );
    free ( linv );

    l = l1dn ( n, h );
    r8mat_print ( n, n, l, "  L1DN:" );
    linv = l1dn_inverse ( n, h );
    r8mat_print ( n, n, linv, "  L1DN_INVERSE:" );
    err = inverse_error ( n, l, linv );
    printf ( "\n" );
    printf ( "  L1DN inverse error = %g\n", err );
    free ( l );
    free ( linv );

    l = l1nd ( n, h );
    r8mat_print ( n, n, l, "  L1ND:" );
    linv = l1nd_inverse ( n, h );
    r8mat_print ( n, n, linv, "  L1ND_INVERSE:" );
    err = inverse_error ( n, l, linv );
    printf ( "\n" );
    printf ( "  L1ND inverse error = %g\n", err );
    free ( l );
    free ( linv );
  }

  return;
}
/******************************************************************************/

void test05 ( )

/******************************************************************************/
/*
  Purpose:

    TEST05 tests L1DD_CHOLESKY and similar routines.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 October 2013

  Author:

    John Burkardt
*/
{
  double *a;
  double *c;
  double err;
  double h;
  int i;
  int n = 5;
  int test;

  printf ( "\n" );
  printf ( "TEST05\n" );
  printf ( "  Compute upper Cholesky factors for the Laplacian:\n" );
  printf ( "  L1DD_CHOLESKY for Dirichlet/Dirichlet BC;\n" );
  printf ( "  L1DN_CHOLESKY for Dirichlet/Neumann BC;\n" );
  printf ( "  L1ND_CHOLESKY for Neumann/Dirichlet BC;\n" );
  printf ( "  L1NN_CHOLESKY for Neumann/Neumann BC;\n" );
  printf ( "  L1PP_CHOLESKY for Periodic BC;\n" );

  for ( test = 1; test <= 2; test++ )
  {
    if ( test == 1 )
    {
      h = 1.0;
    }
    else
    {
      h = 1.0 / ( double ) ( n + 1 );
    }

    printf ( "\n" );
    printf ( "  Using spacing H = %g\n", h );

    a = l1dd ( n, h );
    c = l1dd_cholesky ( n, h );
    r8mat_print ( n, n, c, "  L1DD Cholesky factor:" );
    err = cholesky_upper_error ( n, a, c );
    printf ( "\n" );
    printf ( "  L1DD Cholesky error = %g\n", err );
    free ( a );
    free ( c );

    a = l1dn ( n, h );
    c = l1dn_cholesky ( n, h );
    r8mat_print ( n, n, c, "  L1DN Cholesky factor:" );
    err = cholesky_upper_error ( n, a, c );
    printf ( "\n" );
    printf ( "  L1DN Cholesky error = %g\n", err );
    free ( a );
    free ( c );

    a = l1nd ( n, h );
    c = l1nd_cholesky ( n, h );
    r8mat_print ( n, n, c, "  L1ND Cholesky factor:" );
    err = cholesky_upper_error ( n, a, c );
    printf ( "\n" );
    printf ( "  L1ND Cholesky error = %g\n", err );
    free ( a );
    free ( c );

    a = l1nn ( n, h );
    c = l1nn_cholesky ( n, h );
    r8mat_print ( n, n, c, "  L1NN Cholesky factor:" );
    err = cholesky_upper_error ( n, a, c );
    printf ( "\n" );
    printf ( "  L1NN Cholesky error = %g\n", err );
    free ( a );
    free ( c );

    a = l1pp ( n, h );
    c = l1pp_cholesky ( n, h );
    r8mat_print ( n, n, c, "  L1PP Cholesky factor:" );
    err = cholesky_upper_error ( n, a, c );
    printf ( "\n" );
    printf ( "  L1PP Cholesky error = %g\n", err );
    free ( a );
    free ( c );
  }

  return;
}
/******************************************************************************/

void test06 ( )

/******************************************************************************/
/*
  Purpose:

    TEST06 tests L1DD_LU and similar routines.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 November 2013

  Author:

    John Burkardt
*/
{
  double *a;
  double err;
  double h;
  int i;
  double *l;
  int n = 5;
  int test;
  double *u;

  printf ( "\n" );
  printf ( "TEST06\n" );
  printf ( "  Compute LU factors for the Laplacian:\n" );
  printf ( "  L1DD_LU for Dirichlet/Dirichlet BC;\n" );
  printf ( "  L1DN_LU for Dirichlet/Neumann BC;\n" );
  printf ( "  L1ND_LU for Neumann/Dirichlet BC;\n" );
  printf ( "  L1NN_LU for Neumann/Neumann BC;\n" );
  printf ( "  L1PP_LU for Periodic BC;\n" );

  l = ( double * ) malloc ( n * n * sizeof ( double ) );
  u = ( double * ) malloc ( n * n * sizeof ( double ) );

  for ( test = 1; test <= 2; test++ )
  {
    if ( test == 1 )
    {
      h = 1.0;
    }
    else
    {
      h = 1.0 / ( double ) ( n + 1 );
    }

    printf ( "\n" );
    printf ( "  Using spacing H = %g\n", h );

    a = l1dd ( n, h );
    l1dd_lu ( n, h, l, u );
    r8mat_print ( n, n, l, "  L1DD L factor:" );
    r8mat_print ( n, n, u, "  L1DD U factor:" );
    err = lu_error ( n, a, l, u );
    printf ( "\n" );
    printf ( "  L1DD LU error = %g\n", err );
    free ( a );

    a = l1dn ( n, h );
    l1dn_lu ( n, h, l, u );
    r8mat_print ( n, n, l, "  L1DN L factor:" );
    r8mat_print ( n, n, u, "  L1DN U factor:" );
    err = lu_error ( n, a, l, u );
    printf ( "\n" );
    printf ( "  L1DN LU error = %g\n", err );
    free ( a );

    a = l1nd ( n, h );
    l1nd_lu ( n, h, l, u );
    r8mat_print ( n, n, l, "  L1ND L factor:" );
    r8mat_print ( n, n, u, "  L1ND U factor:" );
    err = lu_error ( n, a, l, u );
    printf ( "\n" );
    printf ( "  L1ND LU error = %g\n", err );
    free ( a );

    a = l1nn ( n, h );
    l1nn_lu ( n, h, l, u );
    r8mat_print ( n, n, l, "  L1NN L factor:" );
    r8mat_print ( n, n, u, "  L1NN U factor:" );
    err = lu_error ( n, a, l, u );
    printf ( "\n" );
    printf ( "  L1NN LU error = %g\n", err );
    free ( a );

    a = l1pp ( n, h );
    l1pp_lu ( n, h, l, u );
    r8mat_print ( n, n, l, "  L1PP L factor:" );
    r8mat_print ( n, n, u, "  L1PP U factor:" );
    err = lu_error ( n, a, l, u );
    printf ( "\n" );
    printf ( "  L1PP LU error = %g\n", err );
    free ( a );
  }

  free ( l );
  free ( u );

  return;
}

