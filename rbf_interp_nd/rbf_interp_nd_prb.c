# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "rbf_interp_nd.h"
# include "r8lib.h"

int main ( void );
void rbf_interp_nd_test01 ( void );
void rbf_interp_nd_test02 ( void );
void rbf_interp_nd_test03 ( void );
void rbf_interp_nd_test04 ( void );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    RBF_INTERP_ND_TEST tests RBF_INTERP_ND.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 July 2012

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "RBF_INTERP_ND_TEST:\n" );
  printf ( "  C version\n" );
  printf ( "  Test the RBF_INTERP_ND library.\n" );
  printf ( "  The R8LIB library is also required.\n" );

  rbf_interp_nd_test01 ( );
  rbf_interp_nd_test02 ( );
  rbf_interp_nd_test03 ( );
  rbf_interp_nd_test04 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "RBF_INTERP_ND_TEST:\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void rbf_interp_nd_test01 ( void )

/******************************************************************************/
/*
  Purpose:

    RBF_INTERP_ND_TEST01 tests RBF_WEIGHTS and RBF_INTERP_ND with PHI1.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 July 2012

  Author:

    John Burkardt
*/
{
  double a;
  double app_error;
  double b;
  int debug = 0;
  double *fd;
  double *fe;
  double *fi;
  int i;
  double int_error;
  int j;
  int m = 2;
  int n1d = 5;
  int nd;
  int ni;
  double r0;
  int seed;
  double *w;
  double *x1d;
  double *xd;
  double *xi;

  printf ( "\n" );
  printf ( "RBF_INTERP_ND_TEST01:\n" );
  printf ( "  RBF_WEIGHT computes weights for RBF interpolation.\n" );
  printf ( "  RBF_INTERP_ND evaluates the RBF interpolant.\n" );
  printf ( "  Use the multiquadratic basis function PHI1(R).\n" );
  a = 0.0;
  b = 2.0;
  r0 = ( b - a ) / ( double ) ( n1d );
  printf ( "  Scale factor R0 = %g\n", r0 );

  x1d = r8vec_linspace_new ( n1d, a, b );
  nd = i4_power ( n1d, m );
  xd = ( double * ) malloc ( m * nd * sizeof ( double ) );

  for ( i = 0; i < m; i++ )
  {
    r8vec_direct_product ( i, n1d, x1d, m, nd, xd );
  }

  if ( debug )
  {
    r8mat_transpose_print ( m, nd, xd, "  The product points:" );
  }

  fd = ( double * ) malloc ( nd * sizeof ( double ) );
  for ( j = 0; j < nd; j++ )
  {
    fd[j] = xd[0+j*m] * xd[1+j*m] * exp ( - xd[0+j*m] * xd[1+j*m] );
  }
  if ( debug )
  {
    r8vec_print ( nd, fd, "  Function data:" );
  }

  w = rbf_weight ( m, nd, xd, r0, phi1, fd );
  if ( debug )
  {
    r8vec_print ( nd, w, "  Weight vector:" );
  }
/*
  #1: Interpolation test.  Does interpolant match function at interpolation points?
*/
  ni = nd;
  xi = r8mat_copy_new ( m, ni, xd );

  fi = rbf_interp_nd ( m, nd, xd, r0, phi1, w, ni, xi );

  int_error = r8vec_norm_affine ( nd, fd, fi ) / ( double ) ( nd );

  printf ( "\n" );
  printf ( "  L2 interpolation error averaged per interpolant node = %g\n", int_error );

  free ( fi );
  free ( xi );
/*
  #2: Approximation test.  Estimate the integral (f-interp(f))^2.
*/
  ni = 1000;
  seed = 123456789;

  xi = r8mat_uniform_ab_new ( m, ni, a, b, &seed );

  fi = rbf_interp_nd ( m, nd, xd, r0, phi1, w, ni, xi );

  fe = ( double * ) malloc ( ni * sizeof ( double ) );
  for ( j = 0; j < ni; j++ )
  {
    fe[j] = xi[0+j*m] * xi[1+j*m] * exp ( - xi[0+j*m] * xi[1+j*m] );
  }

  app_error = pow ( b - a, m ) * r8vec_norm_affine ( ni, fi, fe ) / ( double ) ( ni );

  printf ( "  L2 approximation error averaged per 1000 samples = %g\n", app_error );

  free ( fd );
  free ( fe );
  free ( fi );
  free ( w );
  free ( x1d );
  free ( xd );
  free ( xi );

  return;
}
/******************************************************************************/

void rbf_interp_nd_test02 ( void )

/******************************************************************************/
/*
  Purpose:

    RBF_INTERP_ND_TEST02 tests RBF_WEIGHTS and RBF_INTERP_ND with PHI2.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 July 2012

  Author:

    John Burkardt
*/
{
  double a;
  double app_error;
  double b;
  int debug = 0;
  double *fd;
  double *fe;
  double *fi;
  int i;
  double int_error;
  int j;
  int m = 2;
  int n1d = 5;
  int nd;
  int ni;
  double r0;
  int seed;
  double *w;
  double *x1d;
  double *xd;
  double *xi;

  printf ( "\n" );
  printf ( "RBF_INTERP_ND_TEST02:\n" );
  printf ( "  RBF_WEIGHT computes weights for RBF interpolation.\n" );
  printf ( "  RBF_INTERP_ND evaluates the RBF interpolant.\n" );
  printf ( "  Use the inverse multiquadratic basis function PHI2(R).\n" );
  a = 0.0;
  b = 2.0;
  r0 = ( b - a ) / ( double ) ( n1d );
  printf ( "  Scale factor R0 = %g\n", r0 );

  x1d = r8vec_linspace_new ( n1d, a, b );
  nd = i4_power ( n1d, m );
  xd = ( double * ) malloc ( m * nd * sizeof ( double ) );

  for ( i = 0; i < m; i++ )
  {
    r8vec_direct_product ( i, n1d, x1d, m, nd, xd );
  }

  if ( debug )
  {
    r8mat_transpose_print ( m, nd, xd, "  The product points:" );
  }

  fd = ( double * ) malloc ( nd * sizeof ( double ) );
  for ( j = 0; j < nd; j++ )
  {
    fd[j] = xd[0+j*m] * xd[1+j*m] * exp ( - xd[0+j*m] * xd[1+j*m] );
  }
  if ( debug )
  {
    r8vec_print ( nd, fd, "  Function data:" );
  }

  w = rbf_weight ( m, nd, xd, r0, phi2, fd );

  if ( debug )
  {
    r8vec_print ( nd, w, "  Weight vector:" );
  }
/*
  #1: Interpolation test.  Does interpolant match function at interpolation points?
*/
  ni = nd;
  xi = r8mat_copy_new ( m, ni, xd );

  fi = rbf_interp_nd ( m, nd, xd, r0, phi2, w, ni, xi );

  int_error = r8vec_norm_affine ( nd, fd, fi ) / ( double ) ( nd );

  printf ( "\n" );
  printf ( "  L2 interpolation error averaged per interpolant node = %g\n", int_error );

  free ( fi );
  free ( xi );
/*
  #2: Approximation test.  Estimate the integral (f-interp(f))^2.
*/
  ni = 1000;
  seed = 123456789;

  xi = r8mat_uniform_ab_new ( m, ni, a, b, &seed );

  fi = rbf_interp_nd ( m, nd, xd, r0, phi2, w, ni, xi );

  fe = ( double * ) malloc ( ni * sizeof ( double ) );
  for ( j = 0; j < ni; j++ )
  {
    fe[j] = xi[0+j*m] * xi[1+j*m] * exp ( - xi[0+j*m] * xi[1+j*m] );
  }

  app_error = pow ( b - a, m ) * r8vec_norm_affine ( ni, fi, fe ) / ( double ) ( ni );

  printf ( "  L2 approximation error averaged per 1000 samples = %g\n", app_error );

  free ( fd );
  free ( fe );
  free ( fi );
  free ( w );
  free ( x1d );
  free ( xd );
  free ( xi );

  return;
}
/******************************************************************************/

void rbf_interp_nd_test03 ( void )

/******************************************************************************/
/*
  Purpose:

    RBF_INTERP_ND_TEST03 tests RBF_WEIGHTS and RBF_INTERP_ND with PHI3.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 July 2012

  Author:

    John Burkardt
*/
{
  double a;
  double app_error;
  double b;
  int debug = 0;
  double *fd;
  double *fe;
  double *fi;
  int i;
  double int_error;
  int j;
  int m = 2;
  int n1d = 5;
  int nd;
  int ni;
  double r0;
  int seed;
  double *w;
  double *x1d;
  double *xd;
  double *xi;

  printf ( "\n" );
  printf ( "RBF_INTERP_ND_TEST03:\n" );
  printf ( "  RBF_WEIGHT computes weights for RBF interpolation.\n" );
  printf ( "  RBF_INTERP_ND evaluates the RBF interpolant.\n" );
  printf ( "  Use the thin-plate spline basis function PHI3(R).\n" );
  a = 0.0;
  b = 2.0;
  r0 = ( b - a ) / ( double ) ( n1d );
  printf ( "  Scale factor R0 = %g\n", r0 );

  x1d = r8vec_linspace_new ( n1d, a, b );
  nd = i4_power ( n1d, m );
  xd = ( double * ) malloc ( m * nd * sizeof ( double ) );

  for ( i = 0; i < m; i++ )
  {
    r8vec_direct_product ( i, n1d, x1d, m, nd, xd );
  }

  if ( debug )
  {
    r8mat_transpose_print ( m, nd, xd, "  The product points:" );
  }

  fd = ( double * ) malloc ( nd * sizeof ( double ) );
  for ( j = 0; j < nd; j++ )
  {
    fd[j] = xd[0+j*m] * xd[1+j*m] * exp ( - xd[0+j*m] * xd[1+j*m] );
  }
  if ( debug )
  {
    r8vec_print ( nd, fd, "  Function data:" );
  }

  w = rbf_weight ( m, nd, xd, r0, phi3, fd );

  if ( debug )
  {
    r8vec_print ( nd, w, "  Weight vector:" );
  }
/*
  #1: Interpolation test.  Does interpolant match function at interpolation points?
*/
  ni = nd;
  xi = r8mat_copy_new ( m, ni, xd );

  fi = rbf_interp_nd ( m, nd, xd, r0, phi3, w, ni, xi );

  int_error = r8vec_norm_affine ( nd, fd, fi ) / ( double ) ( nd );

  printf ( "\n" );
  printf ( "  L2 interpolation error averaged per interpolant node = %g\n", int_error );

  free ( fi );
  free ( xi );
/*
  #2: Approximation test.  Estimate the integral (f-interp(f))^2.
*/
  ni = 1000;
  seed = 123456789;

  xi = r8mat_uniform_ab_new ( m, ni, a, b, &seed );

  fi = rbf_interp_nd ( m, nd, xd, r0, phi3, w, ni, xi );

  fe = ( double * ) malloc ( ni * sizeof ( double ) );
  for ( j = 0; j < ni; j++ )
  {
    fe[j] = xi[0+j*m] * xi[1+j*m] * exp ( - xi[0+j*m] * xi[1+j*m] );
  }

  app_error = pow ( b - a, m ) * r8vec_norm_affine ( ni, fi, fe ) / ( double ) ( ni );

  printf ( "  L2 approximation error averaged per 1000 samples = %g\n", app_error );

  free ( fd );
  free ( fe );
  free ( fi );
  free ( w );
  free ( x1d );
  free ( xd );
  free ( xi );

  return;
}
/******************************************************************************/

void rbf_interp_nd_test04 ( void )

/******************************************************************************/
/*
  Purpose:

    RBF_INTERP_ND_TEST04 tests RBF_WEIGHTS and RBF_INTERP_ND with PHI4.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 July 2012

  Author:

    John Burkardt
*/
{
  double a;
  double app_error;
  double b;
  int debug = 0;
  double *fd;
  double *fe;
  double *fi;
  int i;
  double int_error;
  int j;
  int m = 2;
  int n1d = 5;
  int nd;
  int ni;
  double r0;
  int seed;
  double *w;
  double *x1d;
  double *xd;
  double *xi;

  printf ( "\n" );
  printf ( "RBF_INTERP_ND_TEST04:\n" );
  printf ( "  RBF_WEIGHT computes weights for RBF interpolation.\n" );
  printf ( "  RBF_INTERP_ND evaluates the RBF interpolant.\n" );
  printf ( "  Use the gaussian basis function PHI4(R).\n" );
  a = 0.0;
  b = 2.0;
  r0 = ( b - a ) / ( double ) ( n1d );
  printf ( "  Scale factor R0 = %g\n", r0 );

  x1d = r8vec_linspace_new ( n1d, a, b );
  nd = i4_power ( n1d, m );
  xd = ( double * ) malloc ( m * nd * sizeof ( double ) );

  for ( i = 0; i < m; i++ )
  {
    r8vec_direct_product ( i, n1d, x1d, m, nd, xd );
  }

  if ( debug )
  {
    r8mat_transpose_print ( m, nd, xd, "  The product points:" );
  }
  r0 = ( b - a ) / ( double ) ( n1d );

  fd = ( double * ) malloc ( nd * sizeof ( double ) );
  for ( j = 0; j < nd; j++ )
  {
    fd[j] = xd[0+j*m] * xd[1+j*m] * exp ( - xd[0+j*m] * xd[1+j*m] );
  }

  if ( debug )
  {
    r8vec_print ( nd, fd, "  Function data:" );
  }

  w = rbf_weight ( m, nd, xd, r0, phi4, fd );

  if ( debug )
  {
    r8vec_print ( nd, w, "  Weight vector:" );
  }
/*
  #1: Interpolation test.  Does interpolant match function at interpolation points?
*/
  ni = nd;
  xi = r8mat_copy_new ( m, ni, xd );

  fi = rbf_interp_nd ( m, nd, xd, r0, phi4, w, ni, xi );

  int_error = r8vec_norm_affine ( nd, fd, fi ) / ( double ) ( nd );

  printf ( "\n" );
  printf ( "  L2 interpolation error averaged per interpolant node = %g\n", int_error );

  free ( fi );
  free ( xi );
/*
  #2: Approximation test.  Estimate the integral (f-interp(f))^2.
*/
  ni = 1000;
  seed = 123456789;

  xi = r8mat_uniform_ab_new ( m, ni, a, b, &seed );

  fi = rbf_interp_nd ( m, nd, xd, r0, phi4, w, ni, xi );

  fe = ( double * ) malloc ( ni * sizeof ( double ) );
  for ( j = 0; j < ni; j++ )
  {
    fe[j] = xi[0+j*m] * xi[1+j*m] * exp ( - xi[0+j*m] * xi[1+j*m] );
  }

  app_error = pow ( b - a, m ) * r8vec_norm_affine ( ni, fi, fe ) / ( double ) ( ni );

  printf ( "  L2 approximation error averaged per 1000 samples = %g\n", app_error );

  free ( fd );
  free ( fe );
  free ( fi );
  free ( w );
  free ( x1d );
  free ( xd );
  free ( xi );

  return;
}
