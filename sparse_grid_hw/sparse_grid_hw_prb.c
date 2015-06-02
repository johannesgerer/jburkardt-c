# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <string.h>

# include "sparse_grid_hw.h"

int main ( );
void ccl_test ( );
void ccl_sparse_test ( );
void ccs_test ( );
void ccs_sparse_test ( );
void cce_test ( );
void cce_sparse_test ( );
void get_seq_test ( );
void gqn_test ( );
void gqn_sparse_test ( );
void gqn2_sparse_test ( );
void gqu_test ( );
void gqu_sparse_test ( );
void kpn_test ( );
void kpn_sparse_test ( );
void kpu_test ( );
void kpu_sparse_test ( );
void nwspgr_size_test ( );
void nwspgr_time_test ( );
void nwspgr_test ( );
void order_report ( );
void symmetric_sparse_size_test ( );
void tensor_product_test ( );
void tensor_product_cell_test ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for SPARSE_GRID_HW_PRB.

  Discussion:

    SPARSE_GRID_HW_PRB tests the SPARSE_GRID_HW library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 February 2014

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "SPARSE_GRID_HW_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the SPARSE_GRID_HW library.\n" );

  ccl_test ( );
  ccl_sparse_test ( );

  ccs_test ( );  
  ccs_sparse_test ( );

  cce_test ( );
  cce_sparse_test ( );

  get_seq_test ( );

  gqn_test ( );
  gqn_sparse_test ( );
  gqn2_sparse_test ( );

  gqu_test ( );
  gqu_sparse_test ( );

  kpn_test ( );
  kpn_sparse_test ( );

  kpu_test ( );
  kpu_sparse_test ( );

  nwspgr_size_test ( );
  nwspgr_time_test ( );
  nwspgr_test ( );

  order_report ( );

  symmetric_sparse_size_test ( );

  tensor_product_test ( );
  tensor_product_cell_test ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "SPARSE_GRID_HW_PRB\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void ccl_test ( )

/******************************************************************************/
/*
  Purpose:

    CCL_TEST uses CCL_ORDER + CC for 1D quadrature over [0,1].

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 February 2014

  Author:

    John Burkardt
*/
{
  int d;
  double e;
  double exact;
  double *fx;
  int l;
  int n;
  double q;
  double *w;
  double *x;

  printf ( "\n" );
  printf ( "CCL_TEST:\n" );
  printf ( "  CCL_ORDER + CC\n" );
  printf ( "  Clenshaw Curtis Linear (CCL) 1D quadrature:\n" );
  printf ( "\n" );
  printf ( "   Level   Nodes    Estimate  Error\n" );
  printf ( "\n" );

  d = 1;
  exact = fu_integral ( d );

  for ( l = 1; l <= 5; l++ )
  {
    n = ccl_order ( l );

    x = ( double * ) malloc ( n * sizeof ( double ) );
    w = ( double * ) malloc ( n * sizeof ( double ) );

    cc ( n, x, w );

    fx = fu_value ( d, n, x );

    q = r8vec_dot_product ( n, w, fx );

    e = r8_abs ( q - exact ) / exact;

    printf ( "  %2d    %6d  %14.6g  %14.6g\n", l, n, q, e );

    free ( fx );
    free ( w );
    free ( x );
  }
  return;
}
/******************************************************************************/

void ccl_sparse_test ( )

/******************************************************************************/
/*
  Purpose:

    CCL_SPARSE_TEST uses CCL_ORDER + CC to build a sparse grid.

  Discussion:

    CCL is the Clenshaw Curtis Linear growth rule.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 February 2014

  Author:

    John Burkardt

  Local parameters:

    Local, int D, the spatial dimension.

    Local, int MAXK, the maximum level to check.
*/
{
  int d;
  double error_mc;
  double error_sg;
  double estimate;
  double *fx;
  int k;
  int maxk;
  int n;
  int n2;
  int r;
  double *s;
  int s_num;
  int seed;
  double trueval;
  double *w;
  double *x;

  d = 10;
  maxk = 7;

  trueval = fu_integral ( d );

  printf ( "\n" );
  printf ( "CCL_SPARSE_TEST:\n" );
  printf ( "  CCL_ORDER + CC\n" );
  printf ( "  Sparse Clenshaw Curtis Linear sparse grid.\n" );
  printf ( "\n" );
  printf ( "   D  Level   Nodes    SG error    MC error\n" );
  printf ( "\n" );

  for ( k = 2; k <= maxk; k++ )
  {
/*
  Compute sparse grid estimate.
*/
    n = nwspgr_size ( ccl_order, d, k );

    x = ( double * ) malloc ( d * n * sizeof ( double ) );
    w = ( double * ) malloc ( n * sizeof ( double ) );

    nwspgr ( cc, ccl_order, d, k, n, &n2, x, w );

    fx = fu_value ( d, n2, x );
    estimate = r8vec_dot_product ( n2, w, fx );
    error_sg = r8_abs ( ( estimate - trueval ) / trueval );
    free ( fx );
    free ( w );
    free ( x );
/*
  Compute 1000 Monte Carlo estimates with same number of points, and average.
*/
    s_num = 1000;
    s = ( double * ) malloc ( s_num * sizeof ( double ) );
    seed = 123456789;
    for ( r = 0; r < 1000; r++ )
    {
      x = r8mat_uniform_01_new ( d, n2, &seed );
      fx = fu_value ( d, n2, x );
      s[r] = r8vec_sum ( n2, fx ) / ( double ) ( n2 );
      free ( fx );
      free ( x );
    }

    error_mc = 0.0;
    for ( r = 0; r < s_num; r++ )
    {
      error_mc = error_mc + pow ( s[r] - trueval, 2 );
    }
    error_mc = sqrt ( error_mc / ( double ) ( s_num ) ) / trueval;

    printf ( "  %2d     %2d  %6d  %10.5g  %10.5g\n", 
      d, k, n2, error_sg, error_mc );

    free ( s );
  }

  return;
}
/******************************************************************************/

void ccs_test ( )

/******************************************************************************/
/*
  Purpose:

    CCS_TEST uses CCS_ORDER + CC for 1D quadrature over [0,1].

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 December 2012

  Author:

    John Burkardt
*/
{
  int d;
  double e;
  double exact;
  double *fx;
  int l;
  int n;
  double q;
  double *w;
  double *x;

  printf ( "\n" );
  printf ( "CCS_TEST:\n" );
  printf ( "  CCS_ORDER + CC\n" );
  printf ( "  Clenshaw Curtis Slow 1D quadrature:\n" );
  printf ( "\n" );
  printf ( "   Level   Nodes    Estimate  Error\n" );
  printf ( "\n" );

  d = 1;
  exact = fu_integral ( d );

  for ( l = 1; l <= 5; l++ )
  {
    n = ccs_order ( l );

    x = ( double * ) malloc ( n * sizeof ( double ) );
    w = ( double * ) malloc ( n * sizeof ( double ) );

    cc ( n, x, w );

    fx = fu_value ( d, n, x );

    q = r8vec_dot_product ( n, w, fx );

    e = r8_abs ( q - exact ) / exact;

    printf ( "  %2d    %6d  %14.6g  %14.6g\n", l, n, q, e );
    free ( fx );
    free ( w );
    free ( x );
  }
  return;
}
/******************************************************************************/

void ccs_sparse_test ( )

/******************************************************************************/
/*
  Purpose:

    CCS_SPARSE_TEST uses CCS_ORDER + CC to build a sparse grid.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 December 2012

  Author:

    John Burkardt

  Local parameters:

    Local, int D, the spatial dimension.

    Local, int MAXK, the maximum level to check.
*/
{
  int d;
  double error_mc;
  double error_sg;
  double estimate;
  double *fx;
  int k;
  int maxk;
  int n;
  int n2;
  int r;
  double *s;
  int s_num;
  int seed;
  double trueval;
  double *w;
  double *x;

  d = 10;
  maxk = 7;

  trueval = fu_integral ( d );

  printf ( "\n" );
  printf ( "CCS_SPARSE_TEST:\n" );
  printf ( "  CCS_ORDER + CC\n" );
  printf ( "  Sparse Clenshaw Curtis Slow sparse grid.\n" );
  printf ( "\n" );
  printf ( "   D  Level   Nodes    SG error    MC error\n" );
  printf ( "\n" );

  for ( k = 2; k <= maxk; k++ )
  {
/*
  Compute sparse grid estimate.
*/
    n = nwspgr_size ( ccs_order, d, k );

    x = ( double * ) malloc ( d * n * sizeof ( double ) );
    w = ( double * ) malloc ( n * sizeof ( double ) );

    nwspgr ( cc, ccs_order, d, k, n, &n2, x, w );

    fx = fu_value ( d, n2, x );
    estimate = r8vec_dot_product ( n2, w, fx );
    error_sg = r8_abs ( ( estimate - trueval ) / trueval );
    free ( fx );
    free ( w );
    free ( x );
/*
  Compute 1000 Monte Carlo estimates with same number of points, and average.
*/
    s_num = 1000;
    s = ( double * ) malloc ( s_num * sizeof ( double ) );
    seed = 123456789;
    for ( r = 0; r < 1000; r++ )
    {
      x = r8mat_uniform_01_new ( d, n2, &seed );
      fx = fu_value ( d, n2, x );
      s[r] = r8vec_sum ( n2, fx ) / ( double ) ( n2 );
      free ( fx );
      free ( x );
    }

    error_mc = 0.0;
    for ( r = 0; r < s_num; r++ )
    {
      error_mc = error_mc + pow ( s[r] - trueval, 2 );
    }
    error_mc = sqrt ( error_mc / ( double ) ( s_num ) ) / trueval;

    printf ( "  %2d     %2d  %6d  %10.5g  %10.5g\n", 
      d, k, n2, error_sg, error_mc );

    free ( s );
  }

  return;
}
/******************************************************************************/

void cce_test ( )

/******************************************************************************/
/*
  Purpose:

    CCE_TEST uses CCE_ORDER + CC for 1D quadrature over [0,1].

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 December 2012

  Author:

    John Burkardt
*/
{
  int d;
  double e;
  double exact;
  double *fx;
  int l;
  int n;
  double q;
  double *w;
  double *x;

  printf ( "\n" );
  printf ( "CCE_TEST:\n" );
  printf ( "  CCE_ORDER + CC\n" );
  printf ( "  Clenshaw Curtis Exponential 1D quadrature:\n" );
  printf ( "\n" );
  printf ( "   Level   Nodes    Estimate  Error\n" );
  printf ( "\n" );

  d = 1;
  exact = fu_integral ( d );

  for ( l = 1; l <= 5; l++ )
  {
    n = cce_order ( l );

    x = ( double * ) malloc ( n * sizeof ( double ) );
    w = ( double * ) malloc ( n * sizeof ( double ) );

    cc ( n, x, w );

    fx = fu_value ( d, n, x );

    q = r8vec_dot_product ( n, w, fx );

    e = r8_abs ( q - exact ) / exact;

    printf ( "  %2d    %6d  %14.6g  %14.6g\n", l, n, q, e );
    free ( fx );
    free ( w );
    free ( x );
  }
  return;
}
/******************************************************************************/

void cce_sparse_test ( )

/******************************************************************************/
/*
  Purpose:

    CCE_SPARSE_TEST uses the CCE + CC to build a sparse grid.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 December 2012

  Author:

    John Burkardt

  Local parameters:

    Local, int D, the spatial dimension.

    Local, int MAXK, the maximum level to check.
*/
{
  int d;
  double error_mc;
  double error_sg;
  double estimate;
  double *fx;
  int k;
  int maxk;
  int n;
  int n2;
  int r;
  double *s;
  int s_num;
  int seed;
  double trueval;
  double *w;
  double *x;

  d = 10;
  maxk = 7;

  trueval = fu_integral ( d );

  printf ( "\n" );
  printf ( "CCE_SPARSE_TEST:\n" );
  printf ( "  CCE_ORDER + CC\n" );
  printf ( "  Sparse Clenshaw Curtis Exponential sparse grid.\n" );
  printf ( "\n" );
  printf ( "   D  Level   Nodes    SG error    MC error\n" );
  printf ( "\n" );

  for ( k = 2; k <= maxk; k++ )
  {
/*
  Compute sparse grid estimate.
*/
    n = nwspgr_size ( cce_order, d, k );

    x = ( double * ) malloc ( d * n * sizeof ( double ) );
    w = ( double * ) malloc ( n * sizeof ( double ) );

    nwspgr ( cc, cce_order, d, k, n, &n2, x, w );

    fx = fu_value ( d, n2, x );
    estimate = r8vec_dot_product ( n2, w, fx );
    error_sg = r8_abs ( ( estimate - trueval ) / trueval );
    free ( fx );
    free ( w );
    free ( x );
/*
  Compute 1000 Monte Carlo estimates with same number of points, and average.
*/
    s_num = 1000;
    s = ( double * ) malloc ( s_num * sizeof ( double ) );
    seed = 123456789;
    for ( r = 0; r < 1000; r++ )
    {
      x = r8mat_uniform_01_new ( d, n2, &seed );
      fx = fu_value ( d, n2, x );
      s[r] = r8vec_sum ( n2, fx ) / ( double ) ( n2 );
      free ( fx );
      free ( x );
    }

    error_mc = 0.0;
    for ( r = 0; r < s_num; r++ )
    {
      error_mc = error_mc + pow ( s[r] - trueval, 2 );
    }
    error_mc = sqrt ( error_mc / ( double ) ( s_num ) ) / trueval;

    printf ( "  %2d     %2d  %6d  %10.5g  %10.5g\n", 
      d, k, n2, error_sg, error_mc );

    free ( s );
  }

  return;
}
/******************************************************************************/

void get_seq_test ( )

/******************************************************************************/
/*
  Purpose:

    GET_SEQ_TEST tests GET_SEQ.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 December 2012

  Author:

    John Burkardt
*/
{
  int d;

  int *fs;
  int norm;
  int seq_num;

  printf ( "\n" );
  printf ( "GET_SEQ_TEST\n" );
  printf ( "  GET_SEQ returns all D-dimensional vectors that sum to NORM.\n" );

  d = 3;
  norm = 6;

  printf ( "\n" );
  printf ( "  D = %d\n", d );
  printf ( "  NORM = %d\n", norm );

  seq_num = num_seq ( norm - d, d );

  fs = get_seq ( d, norm, seq_num );

  i4mat_print ( seq_num, d, fs, "  The compositions" );

  free ( fs );

  return;
}
/******************************************************************************/

void gqn_test ( )

/******************************************************************************/
/*
  Purpose:

    GQN_TEST uses the GQN function for 1D quadrature over (-oo,+oo).

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 December 2012

  Author:

    John Burkardt
*/
{
  int d;
  double e;
  double exact;
  double *fx;
  int i;
  int l;
  int n;
  double q;
  double *w;
  double *x;

  printf ( "\n" );
  printf ( "\n" );
  printf ( "GQN_TEST:\n" );
  printf ( "  Gauss-Hermite quadrature over (-oo,+oo):\n" );
  printf ( "\n" );
  printf ( "   Level   Nodes    Estimate  Error\n" );
  printf ( "\n" );

  d = 1;
  exact = fn_integral ( d );

  for ( l = 1; l <= 5; l++ )
  {
    n = l;
    x = ( double * ) malloc ( n * sizeof ( double ) );
    w = ( double * ) malloc ( n * sizeof ( double ) );

    gqn ( n, x, w );

    fx = fn_value ( d, n, x );

    q = r8vec_dot_product ( n, w, fx );

    e = r8_abs ( q - exact ) / exact;

    printf ( "  %2d    %6d  %14.6g  %14.6g\n", l, n, q, e );

    free ( fx );
    free ( w );
    free ( x );

  }
  return;
}
/******************************************************************************/

void gqn_sparse_test ( )

/******************************************************************************/
/*
  Purpose:

    GQN_SPARSE_TEST uses the GQN function to build a sparse grid.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 December 2012

  Author:

    John Burkardt

  Local parameters:

    Local, int D, the spatial dimension.

    Local, int MAXK, the maximum level to check.
*/
{
  int d;
  double error_mc;
  double error_sg;
  double estimate;
  double *fx;
  int k;
  int maxk;
  int n;
  int n2;
  int r;
  double *s;
  int s_num;
  int seed;
  double trueval;
  double *w;
  double *x;

  d = 10;
  maxk = 7;

  trueval = fn_integral ( d );

  printf ( "\n" );
  printf ( "GQN_SPARSE_TEST:\n" );
  printf ( "  GQN sparse grid:\n" );
  printf ( "  Sparse Gaussian quadrature with Hermite weight over (-oo,+oo).\n" );
  printf ( "\n" );
  printf ( "   D  Level   Nodes    SG error    MC error\n" );
  printf ( "\n" );

  for ( k = 2; k <= maxk; k++ )
  {
/*
  Compute sparse grid estimate.
*/
    n = nwspgr_size ( gqn_order, d, k );
    x = ( double * ) malloc ( d * n * sizeof ( double ) );
    w = ( double * ) malloc ( n * sizeof ( double ) );
    nwspgr ( gqn, gqn_order, d, k, n, &n2, x, w );
    fx = fn_value ( d, n2, x );
    estimate = r8vec_dot_product ( n2, w, fx );

    error_sg = r8_abs ( ( estimate - trueval ) / trueval );

    free ( fx );
    free ( w );
    free ( x );
/*
  Compute 1000 Monte Carlo estimates with same number of points, and average.
*/
    s_num = 1000;
    s = ( double * ) malloc ( s_num * sizeof ( double ) );
    seed = 123456789;
    for ( r = 0; r < 1000; r++ )
    {
      x = r8mat_normal_01_new ( d, n2, &seed );
      fx = fn_value ( d, n2, x );
      s[r] = r8vec_sum ( n2, fx ) / ( double ) ( n2 );
      free ( fx );
      free ( x );
    }

    error_mc = 0.0;
    for ( r = 0; r < s_num; r++ )
    {
      error_mc = error_mc + pow ( s[r] - trueval, 2 );
    }
    error_mc = sqrt ( error_mc / ( double ) ( s_num ) ) / trueval;

    printf ( "  %2d     %2d  %6d  %10.5g  %10.5g\n", d, k, n2, error_sg, error_mc );

    free ( s );
  }
  return;
}
/******************************************************************************/

void gqn2_sparse_test ( )

/******************************************************************************/
/*
  Purpose:

    GQN2_SPARSE_TEST uses the GQN and GQN2_ORDER functions.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 February 2014

  Author:

    John Burkardt

  Local parameters:

    Local, int D, the spatial dimension.

    Local, int MAXK, the maximum level to check.
*/
{
  int d;
  int j;
  int k;
  int maxk;
  int n;
  int n2;
  double *w;
  double *x;

  d = 2;
  maxk = 4;

  printf ( "\n" );
  printf ( "GQN2_SPARSE_TEST:\n" );
  printf ( "  GQN sparse grid:\n" );
  printf ( "  Gauss-Hermite sparse grids over (-oo,+oo).\n" );
  printf ( "  Use GQN2_ORDER, the growth rule N = 2 * L - 1.\n" );

  for ( k = 2; k <= maxk; k++ )
  {
    printf ( "\n" );
    printf ( "     J      W                X               Y\n" );
    printf ( "\n" );

    n = nwspgr_size ( gqn2_order, d, k );

    x = ( double * ) malloc ( d * n * sizeof ( double ) );
    w = ( double * ) malloc ( n * sizeof ( double ) );
    nwspgr ( gqn, gqn2_order, d, k, n, &n2, x, w );

    for ( j = 0; j < n2; j++ )
    {
      printf ( "  %4d  %14.6g  %14.6g  %14.6g\n", j, w[j], x[0+j*d], x[1+j*d] );
    }
    free ( w );
    free ( x );
  }
  return;
}
/******************************************************************************/

void gqu_test ( )

/******************************************************************************/
/*
  Purpose:

    GQU_TEST uses the GQU function for 1D quadrature over [0,1].

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    10 December 2012

  Author:

    John Burkardt
*/
{
  int d;
  double e;
  double exact;
  double *fx;
  int i;
  int l;
  int n;
  double q;
  double *w;
  double *x;

  printf ( "\n" );
  printf ( "GQU_TEST:\n" );
  printf ( "  Gauss-Legendre quadrature over [0,1]:\n" );
  printf ( "\n" );
  printf ( "   Level   Nodes    Estimate  Error\n" );
  printf ( "\n" );

  d = 1;
  exact = fu_integral ( d );

  for ( l = 1; l <= 5; l++ )
  {
    n = l;
    x = ( double * ) malloc ( n * sizeof ( double ) );
    w = ( double * ) malloc ( n * sizeof ( double ) );

    gqu ( n, x, w );

    fx = fu_value ( d, n, x );

    q = r8vec_dot_product ( n, w, fx );

    e = r8_abs ( q - exact ) / exact;

    printf ( "  %2d    %4d  %14.6g  %14.6g\n", l, n, q, e );

    free ( fx );
    free ( w );
    free ( x );
  }

  return;
}
/******************************************************************************/

void gqu_sparse_test ( )

/******************************************************************************/
/*
  Purpose:

    GQU_SPARSE_TEST uses the GQU function to build a sparse grid.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 December 2012

  Author:

    John Burkardt

  Local parameters:

    Local, int D, the spatial dimension.

    Local, int MAXK, the maximum level to check.
*/
{
  int d;
  double error_mc;
  double error_sg;
  double estimate;
  double *fx;
  int k;
  int maxk;
  int n;
  int n2;
  int r;
  double *s;
  int s_num;
  int seed;
  double trueval;
  double *w;
  double *x;

  d = 10;
  maxk = 7;

  trueval = fu_integral ( d );

  printf ( "\n" );
  printf ( "GQU_SPARSE_TEST:\n" );
  printf ( "  GQU sparse grid:\n" );
  printf ( "  Sparse Gaussian unweighted quadrature over [0,1].\n" );
  printf ( "\n" );
  printf ( "   D  Level   Nodes    SG error    MC error\n" );
  printf ( "\n" );

  for ( k = 2; k <= maxk; k++ )
  {
/*
  Compute sparse grid estimate.
*/
    n = nwspgr_size ( gqu_order, d, k );
    x = ( double * ) malloc ( d * n * sizeof ( double ) );
    w = ( double * ) malloc ( n * sizeof ( double ) );
    nwspgr ( gqu, gqu_order, d, k, n, &n2, x, w );
    fx = fu_value ( d, n2, x );
    estimate = r8vec_dot_product ( n2, w, fx );
    error_sg = r8_abs ( ( estimate - trueval ) / trueval );
    free ( fx );
    free ( w );
    free ( x );
/*
  Compute 1000 Monte Carlo estimates with same number of points, and average.
*/
    s_num = 1000;
    s = ( double * ) malloc ( s_num * sizeof ( double ) );
    seed = 123456789;
    for ( r = 0; r < 1000; r++ )
    {
      x = r8mat_uniform_01_new ( d, n2, &seed );
      fx = fu_value ( d, n2, x );
      s[r] = r8vec_sum ( n2, fx ) / ( double ) ( n2 );
      free ( fx );
      free ( x );
    }

    error_mc = 0.0;
    for ( r = 0; r < s_num; r++ )
    {
      error_mc = error_mc + pow ( s[r] - trueval, 2 );
    }
    error_mc = sqrt ( error_mc / ( double ) ( s_num ) ) / trueval;

    printf ( "  %2d     %2d  %6d  %10.5g  %10.5g\n", d, k, n2, error_sg, error_mc );

    free ( s );
  }

  return;
}
/******************************************************************************/

void kpn_test ( )

/******************************************************************************/
/*
  Purpose:

    KPN_TEST uses the KPN function for 1D quadrature over (-oo,+oo).

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    10 December 2012

  Author:

    John Burkardt
*/
{
  int d;
  double e;
  double exact;
  double *fx;
  int i;
  int l;
  int n;
  double q;
  double *w;
  double *x;

  printf ( "\n" );
  printf ( "KPN_TEST:\n" );
  printf ( "  Kronrod-Patterson-Hermite quadrature over (-oo,+oo):\n" );
  printf ( "\n" );
  printf ( "   Level   Nodes    Estimate  Error\n" );
  printf ( "\n" );

  d = 1;
  exact = fn_integral ( d );

  for ( l = 1; l <= 5; l++ )
  {
    n = kpn_order ( l );
    x = ( double * ) malloc ( n * sizeof ( double ) );
    w = ( double * ) malloc ( n * sizeof ( double ) );
    kpn ( n, x, w );

    fx = fn_value ( d, n, x );

    q = r8vec_dot_product ( n, w, fx );

    e = r8_abs ( q - exact ) / exact;

    printf ( "  %2d    %4d  %14.6g  %14.6g\n", l, n, q, e );

    free ( fx );
    free ( w );
    free ( x );
  }
  return;
}
/******************************************************************************/

void kpn_sparse_test ( )

/******************************************************************************/
/*
  Purpose:

    KPN_SPARSE_TEST uses the KPN function to build a sparse grid.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 December 2012

  Author:

    John Burkardt

  Local parameters:

    Local, int D, the spatial dimension.

    Local, int MAXK, the maximum level to check.
*/
{
  int d;
  double error_mc;
  double error_sg;
  double estimate;
  double *fx;
  int k;
  int maxk;
  int n;
  int n2;
  int r;
  double *s;
  int s_num;
  int seed;
  double trueval;
  double *w;
  double *x;

  d = 10;
  maxk = 7;

  trueval = fn_integral ( d );

  printf ( "\n" );
  printf ( "KPN_SPARSE_TEST:\n" );
  printf ( "  KPN sparse grid:\n" );
  printf ( "  Sparse Kronrod-Patterson quadrature with Hermite weight over (-oo,+oo).\n" );
  printf ( "\n" );
  printf ( "   D  Level   Nodes    SG error    MC error\n" );
  printf ( "\n" );

  for ( k = 2; k <= maxk; k++ )
  {
/*
  Compute sparse grid estimate.
*/
    n = nwspgr_size ( kpn_order, d, k );
    x = ( double * ) malloc ( d * n * sizeof ( double ) );
    w = ( double * ) malloc ( n * sizeof ( double ) );
    nwspgr ( kpn, kpn_order, d, k, n, &n2, x, w );
    fx = fn_value ( d, n2, x );
    estimate = r8vec_dot_product ( n2, w, fx );

    error_sg = r8_abs ( ( estimate - trueval ) / trueval );

    free ( fx );
    free ( w );
    free ( x );
/*
  Compute 1000 Monte Carlo estimates with same number of points, and average.
*/
    s_num = 1000;
    s = ( double * ) malloc ( s_num * sizeof ( double ) );
    seed = 123456789;
    for ( r = 0; r < 1000; r++ )
    {
      x = r8mat_normal_01_new ( d, n2, &seed );
      fx = fn_value ( d, n2, x );
      s[r] = r8vec_sum ( n2, fx ) / ( double ) ( n2 );
      free ( fx );
      free ( x );
    }

    error_mc = 0.0;
    for ( r = 0; r < s_num; r++ )
    {
      error_mc = error_mc + pow ( s[r] - trueval, 2 );
    }
    error_mc = sqrt ( error_mc / ( double ) ( s_num ) ) / trueval;

    printf ( "  %2d     %2d  %6d  %10.5g  %10.5g\n", d, k, n2, error_sg, error_mc );

    free ( s );
  }
  return;
}
/******************************************************************************/

void kpu_test ( )

/******************************************************************************/
/*
  Purpose:

    KPU_TEST uses the KPU function for 1D quadrature over [0,1].

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 December 2012

  Author:

    John Burkardt
*/
{
  int d;
  double e;
  double exact;
  double *fx;
  int i;
  int l;
  int n;
  double q;
  double *w;
  double *x;

  printf ( "\n" );
  printf ( "KPU_TEST:\n" );
  printf ( "  Kronrod-Patterson quadrature over [0,1]:\n" );
  printf ( "\n" );
  printf ( "   Level   Nodes    Estimate  Error\n" );
  printf ( "\n" );

  d = 1;
  exact = fu_integral ( d );

  for ( l = 1; l <= 5; l++ )
  {
    n = kpu_order ( l );
    x = ( double * ) malloc ( n * sizeof ( double ) );
    w = ( double * ) malloc ( n * sizeof ( double ) );
    kpu ( n, x, w );

    fx = fu_value ( d, n, x );

    q = r8vec_dot_product ( n, w, fx );

    e = r8_abs ( q - exact ) / exact;

    printf ( "  %2d    %4d  %14.6g  %14.6g\n", l, n, q, e );

    free ( fx );
    free ( w );
    free ( x );
  }
  return;
}
/******************************************************************************/

void kpu_sparse_test ( )

/******************************************************************************/
/*
  Purpose:

    KPU_SPARSE_TEST uses the KPU function to build a sparse grid.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 December 2012

  Author:

    John Burkardt

  Local parameters:

    Local, int D, the spatial dimension.

    Local, int MAXK, the maximum level to check.
*/
{
  int d;
  double error_mc;
  double error_sg;
  double estimate;
  double *fx;
  int k;
  int maxk;
  int n;
  int n2;
  int r;
  double *s;
  int s_num;
  int seed;
  double trueval;
  double *w;
  double *x;

  d = 10;
  maxk = 7;

  trueval = fu_integral ( d );

  printf ( "\n" );
  printf ( "KPU_SPARSE_TEST:\n" );
  printf ( "  KPU sparse grid:\n" );
  printf ( "  Sparse Kronrod-Patterson unweighted quadrature over [0,1].\n" );
  printf ( "\n" );
  printf ( "   D  Level   Nodes    SG error    MC error\n" );
  printf ( "\n" );

  for ( k = 2; k <= maxk; k++ )
  {
/*
  Compute sparse grid estimate.
*/
    n = nwspgr_size ( kpu_order, d, k );
    x = ( double * ) malloc ( d * n * sizeof ( double ) );
    w = ( double * ) malloc ( n * sizeof ( double ) );
    nwspgr ( kpu, kpu_order, d, k, n, &n2, x, w );
    fx = fu_value ( d, n2, x );
    estimate = r8vec_dot_product ( n2, w, fx );
    error_sg = r8_abs ( ( estimate - trueval ) / trueval );
    free ( fx );
    free ( w );
    free ( x );
/*
  Compute 1000 Monte Carlo estimates with same number of points, and average.
*/
    s_num = 1000;
    s = ( double * ) malloc ( s_num * sizeof ( double ) );
    seed = 123456789;
    for ( r = 0; r < 1000; r++ )
    {
      x = r8mat_uniform_01_new ( d, n2, &seed );
      fx = fu_value ( d, n2, x );
      s[r] = r8vec_sum ( n2, fx ) / ( double ) ( n2 );
      free ( fx );
      free ( x );
    }

    error_mc = 0.0;
    for ( r = 0; r < s_num; r++ )
    {
      error_mc = error_mc + pow ( s[r] - trueval, 2 );
    }
    error_mc = sqrt ( error_mc / ( double ) ( s_num ) ) / trueval;

    printf ( "  %2d     %2d  %6d  %10.5g  %10.5g\n", d, k, n2, error_sg, error_mc );

    free ( s );
  }

  return;
}
/******************************************************************************/

void nwspgr_size_test ( )

/******************************************************************************/
/*
  Purpose:

    NWSPGR_SIZE_TEST tests NWSPGR_SIZE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 January 2013

  Author:

    John Burkardt.
*/
{
  int dim;
  int k;
  int r_size;

  printf ( "\n" );
  printf ( "NWSPGR_SIZE_TEST:\n" );
  printf ( "  NWSPGR_SIZE returns the size of a sparse grid, based on either:\n" );
  printf ( "  one of the built-in 1D rules, or a family of 1D rules\n" );
  printf ( "  supplied by the user.\n" );

  dim = 2;
  k = 3;
  printf ( "\n" );
  printf ( "  Kronrod-Patterson, [0,1], Dim %d, Level %d, Symmetric\n", dim, k );
  printf ( "\n" );
  r_size = nwspgr_size ( kpu_order, dim, k );
  printf ( "  Full          %d\n", r_size );

  dim = 2;
  k = 3;
  printf ( "\n" );
  printf ( "  Kronrod-Patterson, (-oo,+oo), Dim %d, Level %d, Symmetric\n", dim, k );
  printf ( "\n" );
  r_size = nwspgr_size ( kpn_order, dim, k );
  printf ( "  Full          %d\n", r_size );

  dim = 2;
  k = 3;
  printf ( "\n" );
  printf ( "  Gauss-Legendre, [0,1], Dim %d, Level %d, Symmetric\n", dim, k );
  printf ( "\n" );
  r_size = nwspgr_size ( gqu_order, dim, k );
  printf ( "  Full          %d\n", r_size );

  dim = 2;
  k = 3;
  printf ( "\n" );
  printf ( "  Gauss Hermite, (-oo,+oo), [0,1], Dim %d, Level %d, Symmetric\n", dim, k );
  printf ( "\n" );
  r_size = nwspgr_size ( gqn_order, dim, k );
  printf ( "  Full          %d\n", r_size );

  dim = 2;
  k = 3;
  printf ( "\n" );
  printf ( "  Clenshaw Curtis Exponential, [-1,+1], [0,1], Dim %d, Level %d, Unsymmetric\n", dim, k );
  printf ( "\n" );
  r_size = nwspgr_size ( cce_order, dim, k );
  printf ( "  Full          %d\n", r_size );
/*
  Do a table.
*/
  printf ( "\n" );
  printf ( "  Dimension / Level table for Clenshaw Curtis Exponential\n" );
  printf ( "\n" );
  printf ( " Dim: " );
  for ( dim = 1; dim <= 10; dim++ )
  {
    printf ( "  %6d", dim );
  }
  printf ( "\n" );
  printf ( "Level\n" );
  for ( k = 1; k <= 5; k++ )
  {
    printf ( "  %2d  ", k );
    for ( dim = 1; dim <= 10; dim++ )
    {
      r_size = nwspgr_size ( cce_order, dim, k );
      printf ( "  %6d", r_size );
    }
    printf ( "\n" );
  }

  return;
}
/******************************************************************************/

void nwspgr_time_test ( )

/******************************************************************************/
/*
  Purpose:

    NWSPGR_TIME_TEST times NWSPGR.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 January 2013

  Author:

    John Burkardt.
*/
{
  int dim;
  int k;
  double *nodes;
  int r_size;
  int s_size;
  double t1;
  double t2;
  double *weights;

  printf ( "\n" );
  printf ( "  This function measures the time in seconds required by NWSPGR\n" );
  printf ( "  to compute a sparse grid, based on either:\n" );
  printf ( "  one of the built-in 1D rules, or a family of 1D rules\n" );
  printf ( "  supplied by the user.\n" );

  dim = 20;
  k = 5;
  printf ( "\n" );
  printf ( "  Kronrod-Patterson, [0,1], Dim %d, Level %d, Symmetric\n", dim, k );
  printf ( "\n" );
  r_size = nwspgr_size ( kpu_order, dim, k );
  nodes = ( double * ) malloc ( dim * r_size * sizeof ( double ) );
  weights = ( double * ) malloc ( r_size * sizeof ( double ) );
  t1 = cpu_time ( );
  nwspgr ( kpu, kpu_order, dim, k, r_size, &s_size, nodes, weights );
  t2 = cpu_time ( );
  free ( nodes );
  free ( weights );
  printf ( "  Full          %g\n", t2 - t1 );

  dim = 20;
  k = 5;
  printf ( "\n" );
  printf ( "  Kronrod-Patterson, (-oo,+oo), Dim %d, Level %d, Symmetric\n", dim, k );
  printf ( "\n" );
  r_size = nwspgr_size ( kpn_order, dim, k );
   nodes = ( double * ) malloc ( dim * r_size * sizeof ( double ) );
  weights = ( double * ) malloc ( r_size * sizeof ( double ) );
  t1 = cpu_time ( );
  nwspgr ( kpn, kpn_order, dim, k, r_size, &s_size, nodes, weights );
  t2 = cpu_time ( );
  free ( nodes );
  free ( weights );
  printf ( "  Full          %g\n", t2 - t1 );

  dim = 20;
  k = 5;
  printf ( "\n" );
  printf ( "  Gauss-Legendre, [0,1], Dim %d, Level %d, Symmetric\n", dim, k );
  printf ( "\n" );
  r_size = nwspgr_size ( gqu_order, dim, k );
  nodes = ( double * ) malloc ( dim * r_size * sizeof ( double ) );
  weights = ( double * ) malloc ( r_size * sizeof ( double ) );
  t1 = cpu_time ( );
  nwspgr ( gqu, gqu_order, dim, k, r_size, &s_size, nodes, weights );
  t2 = cpu_time ( );
  free ( nodes );
  free ( weights );
  printf ( "  Full          %g\n", t2 - t1 );

  dim = 20;
  k = 5;
  printf ( "\n" );
  printf ( "  Gauss Hermite, (-oo,+oo), [0,1], Dim %d, Level %d, Symmetric\n", dim, k );
  printf ( "\n" );
  r_size = nwspgr_size ( gqn_order, dim, k );
  nodes = ( double * ) malloc ( dim * r_size * sizeof ( double ) );
  weights = ( double * ) malloc ( r_size * sizeof ( double ) );
  t1 = cpu_time ( );
  nwspgr ( gqn, gqn_order, dim, k, r_size, &s_size, nodes, weights );
  t2 = cpu_time ( );
  free ( nodes );
  free ( weights );
  printf ( "  Full          %g\n", t2 - t1 );

  dim = 20;
  k = 5;
  printf ( "\n" );
  printf ( "  Clenshaw Curtis Exponential, [-1,+1], [0,1], Dim %d, Level %d, Unsymmetric\n", dim, k );
  printf ( "\n" );
  r_size = nwspgr_size ( cce_order, dim, k );
  nodes = ( double * ) malloc ( dim * r_size * sizeof ( double ) );
  weights = ( double * ) malloc ( r_size * sizeof ( double ) );
  t1 = cpu_time ( );
  nwspgr ( cc, cce_order, dim, k, r_size, &s_size, nodes, weights );
  t2 = cpu_time ( );
  free ( nodes );
  free ( weights );
  printf ( "  Full          %g\n", t2 - t1 );
/*
  Do a table.
*/
  printf ( "\n" );
  printf ( "  Dimension / Level table for Clenshaw Curtis Exponential\n" );
  printf ( "\n" );
  printf ( " Dim: " );
  for ( dim = 1; dim <= 10; dim++ )
  {
    printf ( "  %6d", dim );
  }
  printf ( "\n" );
  printf ( "Level\n" );
  for ( k = 1; k <= 5; k++ )
  {
    printf ( "  %2d  ", k );
    for ( dim = 1; dim <= 10; dim++ )
    {
      r_size = nwspgr_size ( cce_order, dim, k );
      nodes = ( double * ) malloc ( dim * r_size * sizeof ( double ) );
      weights = ( double * ) malloc ( r_size * sizeof ( double ) );
      t1 = cpu_time ( );
      nwspgr ( cc, cce_order, dim, k, r_size, &s_size, nodes, weights );
      t2 = cpu_time ( );
      free ( nodes );
      free ( weights );
      printf ( "  %10.3g", t2 - t1 );
    }
    printf ( "\n" );
  }
  printf ( "\n" );
  printf ( " Dim: " );
  for ( dim = 11; dim <= 20; dim++ )
  {
    printf ( "  %6d", dim );
  }
  printf ( "\n" );
  printf ( "Level\n" );
  for ( k = 1; k <= 5; k++ )
  {
    printf ( "  %2d  ", k );
    for ( dim = 11; dim <= 20; dim++ )
    {
      r_size = nwspgr_size ( cce_order, dim, k );
      nodes = ( double * ) malloc ( dim * r_size * sizeof ( double ) );
      weights = ( double * ) malloc ( r_size * sizeof ( double ) );
      t1 = cpu_time ( );
      nwspgr ( cc, cce_order, dim, k, r_size, &s_size, nodes, weights );
      t2 = cpu_time ( );
      free ( nodes );
      free ( weights );
      printf ( "  %10.3g", t2 - t1 );
    }
    printf ( "\n" );
  }
  return;
}
/******************************************************************************/

void nwspgr_test ( )

/******************************************************************************/
/*
  Purpose:

    NWSPGR_TEST tests NWSPGR.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 December 2012

  Author:

    John Burkardt.
*/
{
  int dim;
  int k;
  double *nodes;
  int r_size;
  int s_size;
  double *weights;

  printf ( "\n" );
  printf ( "NWSPGR_TEST:\n" );
  printf ( "  NWSPGR generates a sparse grid, based on either:\n" );
  printf ( "  one of the built-in 1D rules, or a family of 1D rules\n" );
  printf ( "  supplied by the user.\n" );

  dim = 2;
  k = 3;
  r_size = nwspgr_size ( kpu_order, dim, k );
  nodes = ( double * ) malloc ( dim * r_size * sizeof ( double ) );
  weights = ( double * ) malloc ( r_size * sizeof ( double ) );
  nwspgr ( kpu, kpu_order, dim, k, r_size, &s_size, nodes, weights );
  quad_rule_print ( dim, s_size, nodes, weights, "  Kronrod-Patterson, [0,1], Dim 2, Level 3" );
  free ( nodes );
  free ( weights );

  dim = 2;
  k = 3;
  r_size = nwspgr_size ( kpn_order, dim, k );
  nodes = ( double * ) malloc ( dim * r_size * sizeof ( double ) );
  weights = ( double * ) malloc ( r_size * sizeof ( double ) );
  nwspgr ( kpn, kpn_order, dim, k, r_size, &s_size, nodes, weights );
  quad_rule_print ( dim, s_size, nodes, weights, "  Kronrod-Patterson, (-oo,+oo), Dim 2, Level 3" );
  free ( nodes );
  free ( weights );

  dim = 2;
  k = 3;
  r_size = nwspgr_size ( gqu_order, dim, k );
  nodes = ( double * ) malloc ( dim * r_size * sizeof ( double ) );
  weights = ( double * ) malloc ( r_size * sizeof ( double ) );
  nwspgr ( gqu, gqu_order, dim, k, r_size, &s_size, nodes, weights );
  quad_rule_print ( dim, s_size, nodes, weights, "  Gauss-Legendre, [0,1], Dim 2, Level 3" );
  free ( nodes );
  free ( weights );

  dim = 2;
  k = 3;
  r_size = nwspgr_size ( gqn_order, dim, k );
  nodes = ( double * ) malloc ( dim * r_size * sizeof ( double ) );
  weights = ( double * ) malloc ( r_size * sizeof ( double ) );
  nwspgr ( gqn, gqn_order, dim, k, r_size, &s_size, nodes, weights );
  quad_rule_print ( dim, s_size, nodes, weights, 
    "  Gauss Hermite, (-oo,+oo), Dim 2, Level 3" );
  free ( nodes );
  free ( weights );

  dim = 2;
  k = 3;
  r_size = nwspgr_size ( cce_order, dim, k );
  nodes = ( double * ) malloc ( dim * r_size * sizeof ( double ) );
  weights = ( double * ) malloc ( r_size * sizeof ( double ) );
  nwspgr ( cc, cce_order, dim, k, r_size, &s_size, nodes, weights );
  quad_rule_print ( dim, s_size, nodes, weights, "  Clenshaw Curtis Exponential, [-1,+1], Dim 2, Level 3" );
  free ( nodes );
  free ( weights );

  return;
}
/******************************************************************************/

void order_report ( )

/******************************************************************************/
/*
  Purpose:

    ORDER_REPORT reports on the order of each family of rules.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    10 December 2012

  Author:

    John Burkardt
*/
{
  int ap;
  int k;
  int kpn_order[5] = {
    1, 3, 9, 19, 35 };
  int l;
  int o;
  int rp;

  printf ( "\n" );
  printf ( "ORDER_REPORT\n" );
  printf ( "  For each family of rules, report:\n" );
  printf ( "  L,  the level index,\n" );
  printf ( "  RP, the required polynomial precision,\n" );
  printf ( "  AP, the actual polynomial precision,\n" );
  printf ( "  O,  the rule order (number of points).\n" );

  printf ( "\n" );
  printf ( "  GQN family\n" );
  printf ( "  Gauss quadrature, exponential weight, (-oo,+oo)\n" );
  printf ( "\n" );
  printf ( "   L  RP  AP   O\n" );
  printf ( "\n" );

  for ( l = 1; l <= 25; l++ )
  {
    rp = 2 * l - 1;
    o = l;
    ap = 2 * o - 1;
    printf ( "  %2d  %2d  %2d  %2d\n", l, rp, ap, o );
  }

  printf ( "\n" );
  printf ( "  GQU family\n" );
  printf ( "  Gauss quadrature, unit weight, [0,1]\n" );
  printf ( "\n" );
  printf ( "   L  RP  AP   O\n" );
  printf ( "\n" );

  for ( l = 1; l <= 25; l++ )
  {
    rp = 2 * l - 1;
    o = l;
    ap = 2 * o - 1;
    printf ( "  %2d  %2d  %2d  %2d\n", l, rp, ap, o );
  }

  printf ( "\n" );
  printf ( "  KPN family\n" );
  printf ( "  Gauss-Kronrod-Patterson quadrature, exponential weight, (-oo,+oo)\n" );
  printf ( "\n" );
  printf ( "   L  RP  AP   O\n" );
  printf ( "\n" );

  k = 1;
  o = 1;
  ap = 1;

  for ( l = 1; l <= 25; l++ )
  {
    rp = 2 * l - 1;

    while ( ap < rp )
    {
      if ( k == 5 )
      {
        printf ( "\n" );
        printf ( "  No higher order rule is available!\n" );
        break;
      }
/*
  Can we use a simple rule?
*/
      if ( rp < kpn_order[k] )
      {
        o = rp;
        ap = rp;
      }
/*
  Otherwise, move to next higher rule.
*/
      else
      {
        k = k + 1;
        ap = 2 * kpn_order[k-1] - kpn_order[k-2];
        o = kpn_order[k-1];
      }
    }
    printf ( "  %2d  %2d  %2d  %2d\n", l, rp, ap, o );
  }

  printf ( "\n" );
  printf ( "  KPU family\n" );
  printf ( "  Gauss-Kronrod-Patterson quadrature, unit weight, [0,1]\n" );
  printf ( "\n" );
  printf ( "   L  RP  AP   O\n" );
  printf ( "\n" );

  for ( l = 1; l <= 25; l++ )
  {
    rp = 2 * l - 1;
    o = 1;
    ap = 1;
    while ( ap < rp )
    {
      o = 2 * ( o + 1 ) - 1;
      ap = ( 3 * o + 1 ) / 2;
    }
    printf ( "  %2d  %2d  %2d  %2d\n", l, rp, ap, o );
  }

  return;
}
/******************************************************************************/

void symmetric_sparse_size_test ( )

/******************************************************************************/
/*
  Purpose:

    SYMMETRIC_SPARSE_SIZE_TEST tests SYMMETRIC_SPARSE_SIZE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 December 2012

  Author:

    Original MATLAB version by Florian Heiss, Viktor Winschel.
    C version by John Burkardt.

  Local parameters:

    Local, int D, the spatial dimension.

    Local, int MAXK, the maximum level to check.
*/
{
  int test_num = 3;

  int dim;
  int dim_test[3] = { 5, 5, 3 };
  double nodes1[6*5] = {
   0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 
   0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 
   0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 
   0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 
   0.0, 1.0, 0.0, 0.0, 0.0, 0.0 };
  double nodes2[21*5] = {
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.73205, 
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.73205, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.73205, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 
    0.0, 0.0, 0.0, 1.0, 1.0, 1.73205, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 
    0.0, 1.0, 1.73205, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0 };
  double nodes3[23*3] = {
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.741964, 1.0, 
    1.0, 1.0, 1.0, 1.0, 1.0, 1.73205, 1.73205, 1.73205, 2.33441, 
    0.0, 0.0, 0.0, 0.0, 0.0, 0.741964, 1.0, 1.0, 1.0, 1.73205, 1.73205, 2.33441, 
    0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.73205, 0.0, 0.0, 1.0, 0.0, 
    0.0, 0.741964, 1.0, 1.73205, 2.33441, 0.0, 0.0, 1.0, 1.73205, 0.0, 1.0, 0.0, 
    0.0, 0.0, 1.0, 1.73205, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0 };
  int r;
  int r_test[3] = { 6, 21, 23 };
  int r2;
  int test;
  double x0;

  printf ( "\n" );
  printf ( "SYMMETRIC_SPARSE_SIZE_TEST\n" );
  printf ( "  Given a symmetric sparse grid rule represented only by\n" );
  printf ( "  the points with positive values, determine the total number\n" );
  printf ( "  of points in the grid.\n" );
  printf ( "\n" );
  printf ( "  For dimension DIM, we report\n" );
  printf ( "  R, the number of points in the positive orthant, and\n" );
  printf ( "  R2, the total number of points.\n" );
  printf ( "\n" );
  printf ( "       DIM         R        R2\n" );
  printf ( "\n" );

  x0 = 0.0;

  for ( test = 0; test < test_num; test++ )
  {
    r = r_test[test];
    dim = dim_test[test];
    if ( test == 0 )
    {
      r2 = symmetric_sparse_size ( r, dim, nodes1, x0 );
    }
    else if ( test == 1 )
    {
      r2 = symmetric_sparse_size ( r, dim, nodes2, x0 );
    }
    else if ( test == 2 )
    {
      r2 = symmetric_sparse_size ( r, dim, nodes3, x0 );
    }
    printf ( "  %8d  %8d  %8d\n", dim, r, r2 );
  }

  return;
}
/******************************************************************************/

void tensor_product_test ( )

/******************************************************************************/
/*
  Purpose:

    TENSOR_PRODUCT_TEST tests TENSOR_PRODUCT.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    10 December 2012

  Author:

    John Burkardt
*/
{
  int d;
  int i;
  int i1;
  int i2;
  int j;
  int n;
  int n1d;
  int order1 = 2;
  int order2 = 3;
  int order3 = 2;
  int *order1d;
  double *w1d;
  double *wnd;
  double w1_1d[2] = { 1.0, 1.0 };
  double w2_1d[3] = { 0.25, 0.50, 0.25 };
  double w3_1d[2] = { 2.50, 2.50 };
  double x1_1d[2] = { -1.0, +1.0 };
  double x2_1d[3] = { 2.0, 2.5, 3.0 };
  double x3_1d[2] = { 10.0, 15.0 };
  double *x1d;
  double *xnd;

  printf ( "\n" );
  printf ( "TENSOR_PRODUCT_TEST:\n" );
  printf ( "  Given a sequence of 1D quadrature rules, construct the\n" );
  printf ( "  tensor product rule.\n" );
/*
  1D rule.
*/
  d = 1;
  order1d = ( int * ) malloc ( d * sizeof ( int ) );

  order1d[0] = order1;

  n1d = i4vec_sum ( d, order1d );
  x1d = ( double * ) malloc ( n1d * sizeof ( double ) );
  w1d = ( double * ) malloc ( n1d * sizeof ( double ) );

  n = i4vec_product ( d, order1d );
  xnd = ( double * ) malloc ( d * n * sizeof ( double ) );
  wnd = ( double * ) malloc ( n * sizeof ( double ) );

  j = 0;
  for ( i = 0; i < order1; i++ )
  {
    x1d[j] = x1_1d[i];
    w1d[j] = w1_1d[i];
    j = j + 1;
  }
  tensor_product ( d, order1d, n1d, x1d, w1d, n, xnd, wnd );
  
  quad_rule_print ( d, n, xnd, wnd, "  A 1D rule over [-1,+1]:" );

  free ( order1d );
  free ( w1d );
  free ( wnd );
  free ( x1d );
  free ( xnd );
/*
  2D rule.
*/
  d = 2;
  order1d = ( int * ) malloc ( d * sizeof ( int ) );

  order1d[0] = order1;
  order1d[1] = order2;

  n1d = i4vec_sum ( d, order1d );
  x1d = ( double * ) malloc ( n1d * sizeof ( double ) );
  w1d = ( double * ) malloc ( n1d * sizeof ( double ) );

  n = i4vec_product ( d, order1d );
  xnd = ( double * ) malloc ( d * n * sizeof ( double ) );
  wnd = ( double * ) malloc ( n * sizeof ( double ) );

  j = 0;
  for ( i = 0; i < order1; i++ )
  {
    x1d[j] = x1_1d[i];
    w1d[j] = w1_1d[i];
    j = j + 1;
  }
  for ( i = 0; i < order2; i++ )
  {
    x1d[j] = x2_1d[i];
    w1d[j] = w2_1d[i];
    j = j + 1;
  }

  tensor_product ( d, order1d, n1d, x1d, w1d, n, xnd, wnd );
  
  quad_rule_print ( d, n, xnd, wnd, "  A 2D rule over [-1,+1] x [2.0,3.0]:" );

  free ( order1d );
  free ( w1d );
  free ( wnd );
  free ( x1d );
  free ( xnd );
/*
  3D rule.
*/
  d = 3;
  order1d = ( int * ) malloc ( d * sizeof ( int ) );

  order1d[0] = order1;
  order1d[1] = order2;
  order1d[2] = order3;

  n1d = i4vec_sum ( d, order1d );
  x1d = ( double * ) malloc ( n1d * sizeof ( double ) );
  w1d = ( double * ) malloc ( n1d * sizeof ( double ) );

  n = i4vec_product ( d, order1d );
  xnd = ( double * ) malloc ( d * n * sizeof ( double ) );
  wnd = ( double * ) malloc ( n * sizeof ( double ) );

  j = 0;
  for ( i = 0; i < order1; i++ )
  {
    x1d[j] = x1_1d[i];
    w1d[j] = w1_1d[i];
    j = j + 1;
  }
  for ( i = 0; i < order2; i++ )
  {
    x1d[j] = x2_1d[i];
    w1d[j] = w2_1d[i];
    j = j + 1;
  }
  for ( i = 0; i < order3; i++ )
  {
    x1d[j] = x3_1d[i];
    w1d[j] = w3_1d[i];
    j = j + 1;
  }

  tensor_product ( d, order1d, n1d, x1d, w1d, n, xnd, wnd );

  quad_rule_print ( d, n, xnd, wnd, 
    "  A 3D rule over [-1,+1] x [2.0,3.0] x [10.0,15.0]:" );

  free ( order1d );
  free ( w1d );
  free ( wnd );
  free ( x1d );
  free ( xnd );

  return;
}
/******************************************************************************/

void tensor_product_cell_test ( )

/******************************************************************************/
/*
  Purpose:

    TENSOR_PRODUCT_CELL_TEST tests TENSOR_PRODUCT_CELL.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    10 December 2012

  Author:

    John Burkardt
*/
{
  int d;
  int i1;
  int i2;
  int n1d;
  int nc;
  int np;
  int nr[3] = { 2, 3, 2 };
  int order1 = 2;
  int order2 = 3;
  int order3 = 2;
  int *order1d;
  int *roff;
  double *w1d;
  double *wc;
  double *wp;
  double w1_1d[2] = { 1.0, 1.0 };
  double w2_1d[3] = { 0.25, 0.50, 0.25 };
  double w3_1d[2] = { 2.50, 2.50 };
  double x1_1d[2] = { -1.0, +1.0 };
  double x2_1d[3] = { 2.0, 2.5, 3.0 };
  double x3_1d[2] = { 10.0, 15.0 };
  double *x1d;
  double *xc;
  double *xp;

  printf ( "\n" );
  printf ( "TENSOR_PRODUCT_TEST_CELL:\n" );
  printf ( "  Given a set of 1D quadrature rules stored in a cell array,\n" );
  printf ( "  construct the tensor product rule.\n" );
/*
  We can construct ROFF once and for all.
*/
  roff = r8cvv_offset ( 3, nr );
/*
  1D rule.
*/
  d = 1;
  nc = i4vec_sum ( d, nr );
  xc = ( double * ) malloc ( nc * sizeof ( double ) );
  r8cvv_rset ( nc, xc, d, roff, 0, x1_1d );
  wc = ( double * ) malloc ( nc * sizeof ( double ) );
  r8cvv_rset ( nc, wc, d, roff, 0, w1_1d );
  np = i4vec_product ( d, nr );
  xp = ( double * ) malloc ( d * np * sizeof ( double ) );
  wp = ( double * ) malloc ( np * sizeof ( double ) );
  tensor_product_cell ( nc, xc, wc, d, nr, roff, np, xp, wp );
  quad_rule_print ( d, np, xp, wp, "  A 1D rule over [-1,+1]:" );

  free ( wc );
  free ( wp );
  free ( xc );
  free ( xp );
/*
  2D rule.
*/
  d = 2;
  nc = i4vec_sum ( d, nr );
  xc = ( double * ) malloc ( nc * sizeof ( double ) );
  r8cvv_rset ( nc, xc, d, roff, 0, x1_1d );
  r8cvv_rset ( nc, xc, d, roff, 1, x2_1d );
  wc = ( double * ) malloc ( nc * sizeof ( double ) );
  r8cvv_rset ( nc, wc, d, roff, 0, w1_1d );
  r8cvv_rset ( nc, wc, d, roff, 1, w2_1d );
  np = i4vec_product ( d, nr );
  xp = ( double * ) malloc ( d * np * sizeof ( double ) );
  wp = ( double * ) malloc ( np * sizeof ( double ) );

  tensor_product_cell ( nc, xc, wc, d, nr, roff, np, xp, wp );

  quad_rule_print ( d, np, xp, wp, "  A 1D rule over [-1,+1]:" );

  free ( wc );
  free ( wp );
  free ( xc );
  free ( xp );
/*
  3D rule.
*/
  d = 3;
  nc = i4vec_sum ( d, nr );
  xc = ( double * ) malloc ( nc * sizeof ( double ) );
  r8cvv_rset ( nc, xc, d, roff, 0, x1_1d );
  r8cvv_rset ( nc, xc, d, roff, 1, x2_1d );
  r8cvv_rset ( nc, xc, d, roff, 2, x3_1d );
  wc = ( double * ) malloc ( nc * sizeof ( double ) );
  r8cvv_rset ( nc, wc, d, roff, 0, w1_1d );
  r8cvv_rset ( nc, wc, d, roff, 1, w2_1d );
  r8cvv_rset ( nc, wc, d, roff, 2, w3_1d );
  np = i4vec_product ( d, nr );
  xp = ( double * ) malloc ( d * np * sizeof ( double ) );
  wp = ( double * ) malloc ( np * sizeof ( double ) );

  tensor_product_cell ( nc, xc, wc, d, nr, roff, np, xp, wp );

  quad_rule_print ( d, np, xp, wp, "  A 1D rule over [-1,+1]:" );

  free ( wc );
  free ( wp );
  free ( xc );
  free ( xp );

  free ( roff );

  return;
}
