# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "spline.h"

int main ( void );
void test001 ( void );
void test002 ( void );
void test003 ( void );
void test004 ( void );
void test005 ( void );
void test006 ( void );

void test01 ( void );
void test02 ( void );
void test03 ( void );
void test04 ( void );
void test05 ( void );
void test06 ( void );
void test07 ( void );
void test08 ( void );
void test09 ( void );

void test10 ( void );
void test11 ( void );
void test115 ( void );
void test116 ( void );
void test12 ( void );
void test125 ( void );
void test126 ( void );
void test127 ( void );
void test13 ( void );
void test14 ( void );
void test145 ( void );
void test15 ( void );
void test16 ( void );
void test17 ( void );
void test18 ( void );
void test19 ( void );

void test20 ( void );
void test205 ( void );
void test21 ( void );
void test215 ( void );
void test22 ( void );
void test225 ( void );
void test23 ( void );
void test235 ( void );
void test24 ( void );

double frunge ( double x );
double fprunge ( double x );
double fpprunge ( double x );
double fcube ( double x );
double fpcube ( double x );
double fppcube ( double x );
void parabola_formula ( double x, double *y, double *yp, double *ypp );

/******************************************************************************/

int main ( )

/******************************************************************************/
//
//  Purpose:
//
//    MAIN is the main program for SPLINE_PRB.
//
//  Discussion:
//
//    SPLINE_PRB tests the SPLINE routines.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 February 2012
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );

  printf ( "\n" );
  printf ( "SPLINE_PRB\n" );
  printf ( "  C version:\n" );
  printf ( "  Test the SPLINE library.\n" );

  test001 ( );
  test002 ( );
  test003 ( );
  test004 ( );
  test005 ( );
  test006 ( );

  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
  test05 ( );
  test06 ( );
  test07 ( );
  test08 ( );
  test09 ( );

  test10 ( );
  test11 ( );
  test115 ( );
  test116 ( );
  test12 ( );
  test125 ( );
  test126 ( );
  test127 ( );
  test13 ( );
  test14 ( );
  test145 ( );
  test15 ( );
  test16 ( );
  test17 ( );
  test18 ( );
  test19 ( );

  test20 ( );
  test205 ( );
  test21 ( );
  test215 ( );
  test22 ( );
  test225 ( );
  test23 ( );
  test235 ( );
  test24 ( );
//
//  Terminate.
//
  printf ( "\n" );
  printf ( "SPLINE_PRB\n" );
  printf ( "  Normal end of execution.\n" );

  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void test001 ( )

/******************************************************************************/
//
//  Purpose:
//
//    TEST001 tests PARABOLA_VAL2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 February 2012
//
//  Author:
//
//    John Burkardt
//
# define NDIM 1
# define NDATA 5
{
  int i;
  int left;
  double xdata[NDATA];
  double xval;
  double ydata[NDIM*NDATA];
  double yval[NDIM];
  double zdata[NDATA];
  double zval[NDIM];

  printf ( "\n" );
  printf ( "TEST001\n" );
  printf ( "  PARABOLA_VAL2 evaluates parabolas through\n" );
  printf ( "    3 points in a table\n" );
  printf ( "\n" );
  printf ( "  Our data tables will actually be parabolas:\n" );
  printf ( "    A: 2*x**2 + 3 * x + 1.\n" );
  printf ( "    B: 4*x**2 - 2 * x + 5.\n" );
  printf ( "\n" );

  for ( i = 0; i < NDATA; i++ )
  {
    xval = 2.0 * ( double ) ( i + 1 );
    xdata[i] = xval;
    ydata[0+i*NDIM] = 2.0 * xval * xval + 3.0 * xval + 1.0;
    zdata[i] = 4.0 * xval * xval - 2.0 * xval + 5.0;
    printf ( "%6d  %10g  %10g  %10g\n", i, xdata[i], ydata[i], zdata[i] );
  }

  printf ( "\n" );
  printf ( "  Interpolated data:\n" );
  printf ( "\n" );
  printf ( "  LEFT, X, Y1, Y2\n" );
  printf ( "\n" );

  for ( i = 1; i <= 5; i++ )
  {
    xval = ( double ) ( 2 * i - 1 );
    left = i;
    if ( NDATA - 2 < left )
    {
      left = NDATA - 2;
    }
    if ( left < 1 )
    {
      left = 1;
    }
    parabola_val2 ( NDIM, NDATA, xdata, ydata, left, xval, yval );
    parabola_val2 ( NDIM, NDATA, xdata, zdata, left, xval, zval );

    printf ( "  %6d  %10g  %10g  %10g\n", left, xval, yval[0], zval[0] );
  }

  return;
# undef NDATA
# undef NDIM
}
/******************************************************************************/

void test002 ( )

/******************************************************************************/
//
//  Purpose:
//
//    TEST002 tests R8VEC_BRACKET.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 February 2012
//
//  Author:
//
//    John Burkardt
//
# define N 10
# define NTEST 6
{
  int i;
  int itest;
  int left;
  int right;
  double x[N];
  double xtest[NTEST] = { -10.0, 1.0, 4.5, 5.0, 10.0, 12.0 };
  double xval;

  printf ( "\n" );
  printf ( "TEST002\n" );
  printf ( "  R8VEC_BRACKET finds a pair of entries in a\n" );
  printf ( "    sorted real array which bracket a value.\n" );

  for ( i = 1; i <= N; i++ )
  {
    x[i-1] = ( double ) i;
  }
  x[5] = x[4];

  r8vec_print ( N, x, "  Sorted array:" );

  for ( itest = 0; itest < NTEST; itest++ )
  {
    xval = xtest[itest];

    printf ( "\n" );
    printf ( "  Search for XVAL = %g\n", xval );

    r8vec_bracket ( N, x, xval, &left, &right );

    printf ( "  X[%d] = %g\n", left - 1, x[left-1] );

    printf ( "  X[%d] = %g\n", right - 1, x[right-1] );

  }

  return;

# undef N
# undef NTEST
}
/******************************************************************************/

void test003 ( )

/******************************************************************************/
//
//  Purpose:
//
//    TEST003 tests R8VEC_BRACKET3.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 February 2012
//
//  Author:
//
//    John Burkardt
//
# define N 10
# define NTEST 6
{
  int i;
  int itest;
  int left;
  double x[N];
  double xtest[NTEST] = { -10.0, 1.0, 4.5, 5.0, 10.0, 12.0 };
  double xval;

  printf ( "\n" );
  printf ( "TEST003\n" );
  printf ( "  R8VEC_BRACKET3 finds a pair of entries in a\n" );
  printf ( "    sorted real array which bracket a value.\n" );

  for ( i = 0; i < N; i++ )
  {
    x[i] = ( double ) ( i + 1 );
  }
  x[5] = x[4];

  r8vec_print ( N, x, "  Sorted array:" );

  left = ( N + 1 ) / 2;

  for ( itest = 0; itest < NTEST; itest++ )
  {
    xval = xtest[itest];

    printf ( "\n" );
    printf ( "  Search for XVAL = %g\n", xval );

    printf ( "  Starting guess for interval is = %d\n", left );

    r8vec_bracket3 ( N, x, xval, &left );

    printf ( "  Nearest interval:\n" );
    printf ( "   X[%d]= %g\n", left-1, x[left-1] );
    printf ( "   X[%d]= %g\n", left, x[left] );

  }

  return;
# undef N
# undef NTEST
}
/******************************************************************************/

void test004 ( )

/******************************************************************************/
//
//  Purpose:
//
//    TEST004 tests R8VEC_ORDER_TYPE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 February 2012
//
//  Author:
//
//    John Burkardt
//
# define N 4
# define NTEST 6
{
  int itest;
  int j;
  int order;
  double x[N];

  printf ( "\n" );
  printf ( "TEST004\n" );
  printf ( "  R8VEC_ORDER_TYPE classifies a real vector as\n" );
  printf ( "  -1: no order\n" );
  printf ( "   0: all equal;\n" );
  printf ( "   1: ascending;\n" );
  printf ( "   2: strictly ascending;\n" );
  printf ( "   3: descending;\n" );
  printf ( "   4: strictly descending.\n" );
  printf ( "\n" );

  for ( itest = 1; itest <= NTEST; itest++ )
  {
    if ( itest == 1 )
    {
      x[0] = 1.0;
      x[1] = 3.0;
      x[2] = 2.0;
      x[3] = 4.0;
    }
    else if ( itest == 2 )
    {
      x[0] = 2.0;
      x[1] = 2.0;
      x[2] = 2.0;
      x[3] = 2.0;
    }
    else if ( itest == 3 )
    {
      x[0] = 1.0;
      x[1] = 2.0;
      x[2] = 2.0;
      x[3] = 4.0;
    }
    else if ( itest == 4 )
    {
      x[0] = 1.0;
      x[1] = 2.0;
      x[2] = 3.0;
      x[3] = 4.0;
    }
    else if ( itest == 5 )
    {
      x[0] = 4.0;
      x[1] = 4.0;
      x[2] = 3.0;
      x[3] = 1.0;
    }
    else if ( itest == 6 )
    {
      x[0] = 9.0;
      x[1] = 7.0;
      x[2] = 3.0;
      x[3] = 0.0;
    }

    order = r8vec_order_type ( N, x );

    printf ( "\n" );
    printf ( "  Vector of order type %d:\n", order );
    printf ( "\n" );

    for ( j = 0; j < N; j++ )
    {
      printf ( "  %6d  %10g\n", j, x[j] );
    }

  }

  return;
# undef N
# undef NTEST
}
/******************************************************************************/

void test005 ( )

/******************************************************************************/
//
//  Purpose:
//
//    TEST005 tests D3_NP_FS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 February 2012
//
//  Author:
//
//    John Burkardt
//
# define N 10
{
  double *a;
  double *b;
  int i;
  int seed;
  double *x;

  printf ( "\n" );
  printf ( "TEST005\n" );
  printf ( "  D3_NP_FS factors and solves a tridiagonal\n" );
  printf ( "    linear system.\n" );
//
//  Set the matrix.
//
  seed = 123456789;
  a = d3_uniform ( N, &seed );
//
//  Set the desired solution.
//
  x = r8vec_indicator_new ( N );
//
//  Compute b = A * x.
//
  b = d3_mxv ( N, a, x );
//
//  Wipe out the solution.
//  Solve the system.
//
  free ( x );

  x = d3_np_fs ( N, a, b );
//
//  Print the solution.
//
  r8vec_print ( N, x, "  Computed solution:" );

  free ( a );
  free ( b );
  free ( x );

  return;
# undef N
}
/******************************************************************************/

void test006 ( )

/******************************************************************************/
//
//  Purpose:
//
//    TEST006 tests DATA_TO_DIF and DIF_VAL.
//
//  Discussion:
//
//    This test demonstrates how divided difference approximation
//    improves with N.
//
//    Evaluate these polynomials at 2.5.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 February 2012
//
//  Author:
//
//    John Burkardt
//
{
# define MAXTAB 8

  double diftab[MAXTAB];
  double error;
  int j;
  int ntab;
  double true_value;
  double xtab[MAXTAB];
  double xval;
  double ytab[MAXTAB];
  double yval;

  printf ( "\n" );
  printf ( "TEST006\n" );
  printf ( "  Approximate Y = EXP(X) using orders 1 to %d.\n", MAXTAB );

  printf ( "\n" );
  printf ( "  Original data:\n" );
  printf ( "\n" );
  printf ( "       X          Y\n" );
  printf ( "\n" );
  for ( j = 0; j < MAXTAB; j++ )
  {
    xtab[j] = ( double ) j;
    ytab[j] = exp ( xtab[j] );
    printf ( "  %12g  %12g\n", xtab[j], ytab[j] );
  }

  xval = 2.5;
  true_value = exp ( xval );
  printf ( "\n" );
  printf ( "  Evaluate at X = %g where EXP(X) = %g\n", xval, true_value );
  printf ( "\n" );
  printf ( "  Order  Approximate Y     Error\n" );
  printf ( "\n" );

  for ( ntab = 1; ntab <= MAXTAB; ntab++ )
  {

    for ( j = 0; j < ntab; j++ )
    {
      xtab[j] = ( double ) j;
      ytab[j] = exp ( xtab[j] );
    }

    data_to_dif ( ntab, xtab, ytab, diftab );

    yval = dif_val ( ntab, xtab, diftab, xval );

    error = yval - true_value;

    printf ( "  %6d  %10g  %10g\n", ntab, yval, error );
  }

  return;
# undef MAXTAB
}
/******************************************************************************/

void test01 ( )

/******************************************************************************/
//
//  Purpose:
//
//    TEST01 tests BASIS_FUNCTION_B_VAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 February 2012
//
//  Author:
//
//    John Burkardt
//
{
# define NDATA 5

  int i;
  int j;
  int jhi;
  char mark;
  int nsample = 4;
  double tdata[NDATA] = { 0.0, 1.0, 4.0, 6.0, 10.0 };
  double thi;
  double tlo;
  double tval;
  double yval;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  BASIS_FUNCTION_B_VAL evaluates the \n" );
  printf ( "    B spline basis function.\n" );
  printf ( "\n" );
  printf ( "           T            B(T)\n" );
  printf ( "\n" );

  for ( i = 0; i <= NDATA; i++ )
  {
    if ( i == 0 )
    {
      tlo = tdata[0] - 0.5 * ( tdata[1] - tdata[0] );
      thi = tdata[0];
    }
    else if ( i < NDATA )
    {
      tlo = tdata[i-1];
      thi = tdata[i];
    }
    else if ( NDATA <= i )
    {
      tlo = tdata[NDATA-1];
      thi = tdata[NDATA-1] + 0.5 * ( tdata[NDATA-1] - tdata[NDATA-2] );
    }

    if ( i < NDATA )
    {
      jhi = nsample - 1;
    }
    else
    {
      jhi = nsample;
    }

    for ( j = 0; j <= jhi; j++ )
    {

      tval = ( ( double ) ( nsample - j ) * tlo
             + ( double ) (           j ) * thi )
             / ( double ) ( nsample     );

      yval = basis_function_b_val ( tdata, tval );

      if ( 0 < i & j == 0 )
      {
        mark = '*';
      }
      else
      {
        mark = ' ';
      }
      printf ( "  %c  %12g  %12g\n", mark, tval, yval );
    }
  }

  return;
# undef NDATA
}
/******************************************************************************/

void test02 ( )

/******************************************************************************/
//
//  Purpose:
//
//    TEST02 tests BASIS_FUNCTION_BETA_VAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 February 2012
//
//  Author:
//
//    John Burkardt
//
{
# define NDATA 5

  double beta1;
  double beta2;
  int i;
  int j;
  int jhi;
  char mark;
  int nsample = 4;
  double tdata[NDATA] = { 0.0, 1.0, 4.0, 6.0, 10.0 };
  double thi;
  double tlo;
  double tval;
  double yval;

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  BASIS_FUNCTION_BETA_VAL evaluates the \n" );
  printf ( "    Beta spline basis function.\n" );

  beta1 = 1.0;
  beta2 = 0.0;

  printf ( "\n" );
  printf ( "  BETA1 = %g\n", beta1 );
  printf ( "  BETA2 = %g\n", beta2 );
  printf ( "\n" );
  printf ( "              T           B(T)\n" );
  printf ( "\n" );

  for ( i = 0; i <= NDATA; i++ )
  {
    if ( i == 0 )
    {
      tlo = tdata[0] - 0.5 * ( tdata[1] - tdata[0] );
      thi = tdata[0];
    }
    else if ( i < NDATA )
    {
      tlo = tdata[i-1];
      thi = tdata[i];
    }
    else if ( NDATA <= i )
    {
      tlo = tdata[NDATA-1];
      thi = tdata[NDATA-1] + 0.5 * ( tdata[NDATA-1] - tdata[NDATA-2] );
    }

    if ( i < NDATA )
    {
      jhi = nsample - 1;
    }
    else
    {
      jhi = nsample;
    }

    for ( j = 0; j <= jhi; j++ )
    {
      tval = ( ( double ) ( nsample - j ) * tlo
             + ( double ) (           j ) * thi )
             / ( double ) ( nsample     );

      yval = basis_function_beta_val ( beta1, beta2, tdata, tval );

      if ( 0 < i & j == 0 )
      {
        mark = '*';
      }
      else
      {
        mark = ' ';
      }
      printf ( "  %c  %12g  %12g\n", mark, tval, yval );;
    }
  }

  beta1 = 1.0;
  beta2 = 100.0;

  printf ( "\n" );
  printf ( "  BETA1 = %g\n", beta1 );
  printf ( "  BETA2 = %g\n", beta2 );
  printf ( "\n" );
  printf ( "              T           B(T)\n" );
  printf ( "\n" );

  for ( i = 0; i <= NDATA; i++ )
  {
    if ( i == 0 )
    {
      tlo = tdata[0] - 0.5 * ( tdata[1] - tdata[0] );
      thi = tdata[0];
    }
    else if ( i < NDATA )
    {
      tlo = tdata[i-1];
      thi = tdata[i];
    }
    else if ( NDATA <= i )
    {
      tlo = tdata[NDATA-1];
      thi = tdata[NDATA-1] + 0.5 * ( tdata[NDATA-1] - tdata[NDATA-2] );
    }

    if ( i < NDATA )
    {
      jhi = nsample - 1;
    }
    else
    {
      jhi = nsample;
    }

    for ( j = 0; j <= jhi; j++ )
    {
      tval = ( ( double ) ( nsample - j ) * tlo
             + ( double ) (           j ) * thi )
             / ( double ) ( nsample     );

      yval = basis_function_beta_val ( beta1, beta2, tdata, tval );

      if ( 0 < i & j == 0 )
      {
        mark = '*';
      }
      else
      {
        mark = ' ';
      }

      printf ( "  %c  %12g  %12g\n", mark, tval, yval );
    }
  }

  beta1 = 100.0;
  beta2 = 0.0;

  printf ( "\n" );
  printf ( "  BETA1 = %g\n", beta1 );
  printf ( "  BETA2 = %g\n", beta2 );
  printf ( "\n" );
  printf ( "              T           B(T)\n" );
  printf ( "\n" );

  for ( i = 0; i <= NDATA; i++ )
  {
    if ( i == 0 )
    {
      tlo = tdata[0] - 0.5 * ( tdata[1] - tdata[0] );
      thi = tdata[0];
    }
    else if ( i < NDATA )
    {
      tlo = tdata[i-1];
      thi = tdata[i];
    }
    else if ( NDATA <= i )
    {
      tlo = tdata[NDATA-1];
      thi = tdata[NDATA-1] + 0.5 * ( tdata[NDATA-1] - tdata[NDATA-2] );
    }

    if ( i < NDATA )
    {
      jhi = nsample - 1;
    }
    else
    {
      jhi = nsample;
    }

    for ( j = 0; j <= jhi; j++ )
    {
      tval = ( ( double ) ( nsample - j ) * tlo
             + ( double ) (           j ) * thi )
             / ( double ) ( nsample     );

      yval = basis_function_beta_val ( beta1, beta2, tdata, tval );

      if ( 0 < i & j == 0 )
      {
        mark = '*';
      }
      else
      {
        mark = ' ';
      }
      printf ( "  %c  %12g  %12g\n", mark, tval, yval );
    }
  }

  return;
# undef NDATA
}
/******************************************************************************/

void test03 ( )

/******************************************************************************/
//
//  Purpose:
//
//    TEST03 tests BASIS_MATRIX_B_UNI and BASIS_MATRIX_TMP.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 February 2012
//
//  Author:
//
//    John Burkardt
//
{
# define N 4
# define NDATA 4

  int i;
  int j;
  int jhi;
  int left;
  char mark;
  double *mbasis;
  int nsample = 4;
  double tdata[NDATA] = { -1.0, 0.0, 1.0, 2.0 };
  double thi;
  double tlo;
  double tval;
  double ydata[NDATA] = { 4.0, 7.0, 12.0, 19.0 };
  double yval;

  printf ( "\n" );
  printf ( "TEST03\n" );
  printf ( "  BASIS_MATRIX_B_UNI sets up the basis matrix\n" );
  printf ( "    for the uniform B spline.\n" );

  mbasis = basis_matrix_b_uni ( );

  printf ( "\n" );
  printf ( "    TDATA, YDATA\n" );
  printf ( "\n" );
  for ( i = 0; i < NDATA; i++ )
  {
      printf ( "  %12g  %12g\n", tdata[i], ydata[i] );
  }

  printf ( "\n" );
  printf ( "              T      Spline(T)\n" );
  printf ( "\n" );

  left = 2;

  for ( i = 0; i <= NDATA; i++ )
  {
    if ( i == 0 )
    {
      tlo = tdata[0] - 0.5 * ( tdata[1] - tdata[0] );
      thi = tdata[0];
    }
    else if ( i < NDATA )
    {
      tlo = tdata[i-1];
      thi = tdata[i];
    }
    else if ( NDATA <= i )
    {
      tlo = tdata[NDATA-1];
      thi = tdata[NDATA-1] + 0.5 * ( tdata[NDATA-1] - tdata[NDATA-2] );
    }

    if ( i < NDATA )
    {
      jhi = nsample - 1;
    }
    else
    {
      jhi = nsample;
    }

    for ( j = 0; j <= jhi; j++ )
    {
      tval = ( ( double ) ( nsample - j ) * tlo
             + ( double ) (           j ) * thi )
             / ( double ) ( nsample     );

      yval = basis_matrix_tmp ( left, N, mbasis, NDATA, tdata, ydata, tval );

      if ( 0 < i & j == 0 )
      {
        mark = '*';
      }
      else
      {
        mark = ' ';
      }
      printf ( "  %c  %12g  %12g\n", mark, tval, yval );
    }
  }

  free ( mbasis );

  return;
# undef N
# undef NDATA
}
/******************************************************************************/

void test04 ( )

/******************************************************************************/
//
//  Purpose:
//
//    TEST04 tests BASIS_MATRIX_BETA_UNI and BASIS_MATRIX_TMP.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 February 2012
//
//  Author:
//
//    John Burkardt
//
{
# define N 4
# define NDATA 4

  double beta1;
  double beta2;
  int i;
  int j;
  int jhi;
  int left;
  char mark;
  double *mbasis;
  int nsample = 4;
  double tdata[NDATA] = { -1.0, 0.0, 1.0, 2.0 };
  double thi;
  double tlo;
  double tval;
  double ydata[NDATA] = { 4.0, 7.0, 12.0, 19.0 };
  double yval;

  printf ( "\n" );
  printf ( "TEST04\n" );
  printf ( "  BASIS_MATRIX_BETA_UNI sets up the basis matrix\n" );
  printf ( "    for the uniform beta spline.\n" );
//
//  First test
//
  beta1 = 1.0;
  beta2 = 0.0;

  printf ( "\n" );
  printf ( "  BETA1 = %g\n", beta1 );
  printf ( "  BETA2 = %g\n", beta2 );

  mbasis = basis_matrix_beta_uni ( beta1, beta2 );

  printf ( "\n" );
  printf ( "    TDATA, YDATA\n" );
  printf ( "\n" );
  for ( i = 0; i < NDATA; i++ )
  {
    printf ( "  %12g  %12g\n", tdata[i], ydata[i] );
  }

  left = 2;

  printf ( "\n" );
  printf ( "              T      Spline(T)\n" );
  printf ( "\n" );

  for ( i = 0; i <= NDATA; i++ )
  {
    if ( i == 0 )
    {
      tlo = tdata[0] - 0.5 * ( tdata[1] - tdata[0] );
      thi = tdata[0];
    }
    else if ( i < NDATA )
    {
      tlo = tdata[i-1];
      thi = tdata[i];
    }
    else if ( NDATA <= i )
    {
      tlo = tdata[NDATA-1];
      thi = tdata[NDATA-1] + 0.5 * ( tdata[NDATA-1] - tdata[NDATA-2] );
    }

    if ( i < NDATA )
    {
      jhi = nsample - 1;
    }
    else
    {
      jhi = nsample;
    }

    for ( j = 0; j <= jhi; j++ )
    {
      tval = ( ( double ) ( nsample - j ) * tlo
             + ( double ) (           j ) * thi )
             / ( double ) ( nsample     );

      yval = basis_matrix_tmp ( left, N, mbasis, NDATA, tdata, ydata, tval );

      if ( 0 < i & j == 0 )
      {
        mark = '*';
      }
      else
      {
        mark = ' ';
      }
      printf ( "  %c  %12g  %12g\n", mark, tval, yval );
    }
  }
//
//  Second test
//
  beta1 = 1.0;
  beta2 = 100.0;

  printf ( "\n" );
  printf ( "  BETA1 = %g\n", beta1 );
  printf ( "  BETA2 = %g\n", beta2 );

  free ( mbasis );
  mbasis = basis_matrix_beta_uni ( beta1, beta2 );

  printf ( "\n" );
  printf ( "    TDATA, YDATA\n" );
  printf ( "\n" );
  for ( i = 0; i < NDATA; i++ )
  {
    printf ( "  %12g  %12g\n", tdata[i], ydata[i] );
  }

  left = 2;

  printf ( "\n" );
  printf ( "              T      Spline(T)\n" );
  printf ( "\n" );

  for ( i = 0; i <= NDATA; i++ )
  {
    if ( i == 0 )
    {
      tlo = tdata[0] - 0.5 * ( tdata[1] - tdata[0] );
      thi = tdata[0];
    }
    else if ( i < NDATA )
    {
      tlo = tdata[i-1];
      thi = tdata[i];
    }
    else if ( NDATA <= i )
    {
      tlo = tdata[NDATA-1];
      thi = tdata[NDATA-1] + 0.5 * ( tdata[NDATA-1] - tdata[NDATA-2] );
    }

    if ( i < NDATA )
    {
      jhi = nsample - 1;
    }
    else
    {
      jhi = nsample;
    }

    for ( j = 0; j <= jhi; j++ )
    {
      tval = ( ( double ) ( nsample - j ) * tlo
             + ( double ) (           j ) * thi )
             / ( double ) ( nsample     );

      yval = basis_matrix_tmp ( left, N, mbasis, NDATA, tdata, ydata, tval );

      if ( 0 < i & j == 0 )
      {
        mark = '*';
      }
      else
      {
        mark = ' ';
      }
      printf ( "  %c  %12g  %12g\n", mark, tval, yval );
    }
  }
//
//  Third test
//
  beta1 = 100.0;
  beta2 = 0.0;

  printf ( "\n" );
  printf ( "  BETA1 = %g\n", beta1 );
  printf ( "  BETA2 = %g\n", beta2 );

  free ( mbasis );
  mbasis = basis_matrix_beta_uni ( beta1, beta2 );

  printf ( "\n" );
  printf ( "    TDATA, YDATA\n" );
  printf ( "\n" );
  for ( i = 0; i < NDATA; i++ )
  {
    printf ( "  %12g  %12g\n", tdata[i], ydata[i] );
  }

  left = 2;

  printf ( "\n" );
  printf ( "              T      Spline(T)\n" );
  printf ( "\n" );

  for ( i = 0; i <= NDATA; i++ )
  {
    if ( i == 0 )
    {
      tlo = tdata[0] - 0.5 * ( tdata[1] - tdata[0] );
      thi = tdata[0];
    }
    else if ( i < NDATA )
    {
      tlo = tdata[i-1];
      thi = tdata[i];
    }
    else if ( NDATA <= i )
    {
      tlo = tdata[NDATA-1];
      thi = tdata[NDATA-1] + 0.5 * ( tdata[NDATA-1] - tdata[NDATA-2] );
    }

    if ( i < NDATA )
    {
      jhi = nsample - 1;
    }
    else
    {
      jhi = nsample;
    }

    for ( j = 0; j <= jhi; j++ )
    {
      tval = ( ( double ) ( nsample - j ) * tlo
             + ( double ) (           j ) * thi )
             / ( double ) ( nsample     );

      yval = basis_matrix_tmp ( left, N, mbasis, NDATA, tdata, ydata, tval );

      if ( 0 < i & j == 0 )
      {
        mark = '*';
      }
      else
      {
        mark = ' ';
      }
      printf ( "  %c  %12g  %12g\n", mark, tval, yval );
    }
  }

  free ( mbasis );

  return;
# undef N
# undef NDATA
}
/******************************************************************************/

void test05 ( )

/******************************************************************************/
//
//  Purpose:
//
//    TEST05 tests BASIS_MATRIX_BEZIER and BASIS_MATRIX_TMP.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 4
# define NDATA 4

  int i;
  int j;
  int jhi;
  int left;
  char mark;
  double *mbasis;
  int nsample = 4;
  double tdata[NDATA] = { 0.0, 0.0, 1.0, 1.0 };
  double thi;
  double tlo;
  double tval;
  double ydata[NDATA] = { 7.0,  8.3333333,   10.0, 12.0 };
  double yval;

  printf ( "\n" );
  printf ( "TEST05\n" );
  printf ( "  BASIS_MATRIX_BEZIER sets up the basis\n" );
  printf ( "    matrix for the uniform Bezier spline.\n" );

  mbasis = basis_matrix_bezier ( );

  printf ( "\n" );
  printf ( "    TDATA, YDATA\n" );
  printf ( "\n" );
  for ( i = 0; i < NDATA; i++ )
  {
    printf ( "  %12g  %12g\n", tdata[i], ydata[i] );
  }

  left = 2;

  printf ( "\n" );
  printf ( "              T      Spline(T)\n" );
  printf ( "\n" );


  for ( i = 0; i <= NDATA; i++ )
  {
    if ( i == 0 )
    {
      tlo = tdata[0] - 0.5 * ( tdata[1] - tdata[0] );
      thi = tdata[0];
    }
    else if ( i < NDATA )
    {
      tlo = tdata[i-1];
      thi = tdata[i];
    }
    else if ( NDATA <= i )
    {
      tlo = tdata[NDATA-1];
      thi = tdata[NDATA-1] + 0.5 * ( tdata[NDATA-1] - tdata[NDATA-2] );
    }

    if ( i < NDATA )
    {
      jhi = nsample - 1;
    }
    else
    {
      jhi = nsample;
    }

    for ( j = 0; j <= jhi; j++ )
    {
      tval = ( ( double ) ( nsample - j ) * tlo
             + ( double ) (           j ) * thi )
             / ( double ) ( nsample     );

      yval = basis_matrix_tmp ( left, N, mbasis, NDATA, tdata, ydata, tval );

      if ( 0 < i & j == 0 )
      {
        mark = '*';
      }
      else
      {
        mark = ' ';
      }
      printf ( "  %c  %12g  %12g\n", mark, tval, yval );
    }
  }

  free ( mbasis );

  return;
# undef N
# undef NDATA
}
/******************************************************************************/

void test06 ( )

/******************************************************************************/
//
//  Purpose:
//
//    TEST06 tests BASIS_MATRIX_HERMITE and BASIS_MATRIX_TMP.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 4
# define NDATA 4

  int i;
  int j;
  int jhi;
  int left;
  char mark;
  double *mbasis;
  int nsample = 4;
  double tdata[NDATA] = { 0.0, 0.0, 1.0, 1.0 };
  double thi;
  double tlo;
  double tval;
  double ydata[NDATA] = { 7.0, 12.0, 4.0, 6.0 };
  double yval;

  printf ( "\n" );
  printf ( "TEST06\n" );
  printf ( "  BASIS_MATRIX_HERMITE sets up the basis matrix\n" );
  printf ( "    for the Hermite spline.\n" );

  mbasis = basis_matrix_hermite ( );

  printf ( "\n" );
  printf ( "    TDATA, YDATA\n" );
  printf ( "\n" );
  for ( i = 0; i < NDATA; i++ )
  {
    printf ( "  %12g  %12g\n", tdata[i], ydata[i] );
  }

  left = 2;

  printf ( "\n" );
  printf ( "              T      Spline(T)\n" );
  printf ( "\n" );


  for ( i = 0; i <= NDATA; i++ )
  {
    if ( i == 0 )
    {
      tlo = tdata[0] - 0.5 * ( tdata[1] - tdata[0] );
      thi = tdata[0];
    }
    else if ( i < NDATA )
    {
      tlo = tdata[i-1];
      thi = tdata[i];
    }
    else if ( NDATA <= i )
    {
      tlo = tdata[NDATA-1];
      thi = tdata[NDATA-1] + 0.5 * ( tdata[NDATA-1] - tdata[NDATA-2] );
    }

    if ( i < NDATA )
    {
      jhi = nsample - 1;
    }
    else
    {
      jhi = nsample;
    }

    for ( j = 0; j <= jhi; j++ )
    {
      tval = ( ( double ) ( nsample - j ) * tlo
             + ( double ) (           j ) * thi )
             / ( double ) ( nsample     );

      yval = basis_matrix_tmp ( left, N, mbasis, NDATA, tdata, ydata, tval );

      if ( 0 < i & j == 0 )
      {
        mark = '*';
      }
      else
      {
        mark = ' ';
      }
      printf ( "  %c  %12g  %12g\n", mark, tval, yval );
    }
  }

  free ( mbasis );

  return;
# undef N
# undef NDATA
}
/******************************************************************************/

void test07 ( )

/******************************************************************************/
//
//  Purpose:
//
//    TEST07 tests BASIS_MATRIX_OVERHAUSER_UNI and BASIS_MATRIX_TMP.
//
//  Discussion:
//
//    YDATA(1:NDATA) = ( TDATA(1:NDATA) + 2 )**2 + 3
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 4
# define NDATA 4

  int i;
  int j;
  int jhi;
  int left;
  char mark;
  double *mbasis;
  int nsample = 4;
  double tdata[NDATA] = { -1.0, 0.0, 1.0, 2.0 };
  double thi;
  double tlo;
  double tval;
  double ydata[NDATA] = { 4.0, 7.0, 12.0, 19.0 };
  double yval;

  printf ( "\n" );
  printf ( "TEST07\n" );
  printf ( "  BASIS_MATRIX_OVERHAUSER_UNI sets up the basis\n" );
  printf ( "    matrix for the uniform Overhauser spline.\n" );

  mbasis = basis_matrix_overhauser_uni ( );

  printf ( "\n" );
  printf ( "    TDATA, YDATA\n" );
  printf ( "\n" );
  for ( i = 0; i < NDATA; i++ )
  {
    printf ( "  %12g  %12g\n", tdata[i], ydata[i] );
  }

  left = 2;

  printf ( "\n" );
  printf ( "              T      Spline(T)\n" );
  printf ( "\n" );

  for ( i = 0; i <= NDATA; i++ )
  {
    if ( i == 0 )
    {
      tlo = tdata[0] - 0.5 * ( tdata[1] - tdata[0] );
      thi = tdata[0];
    }
    else if ( i < NDATA )
    {
      tlo = tdata[i-1];
      thi = tdata[i];
    }
    else if ( NDATA <= i )
    {
      tlo = tdata[NDATA-1];
      thi = tdata[NDATA-1] + 0.5 * ( tdata[NDATA-1] - tdata[NDATA-2] );
    }

    if ( i < NDATA )
    {
      jhi = nsample - 1;
    }
    else
    {
      jhi = nsample;
    }

    for ( j = 0; j <= jhi; j++ )
    {
      tval = ( ( double ) ( nsample - j ) * tlo
             + ( double ) (           j ) * thi )
             / ( double ) ( nsample     );

      yval = basis_matrix_tmp ( left, N, mbasis, NDATA, tdata, ydata, tval );

      if ( 0 < i & j == 0 )
      {
        mark = '*';
      }
      else
      {
        mark = ' ';
      }
      printf ( "  %c  %12g  %12g\n", mark, tval, yval );
    }
  }

  free ( mbasis );

  return;
# undef N
# undef NDATA
}
/******************************************************************************/

void test08 ( )

/******************************************************************************/
//
//  Purpose:
//
//    TEST08 tests BASIS_MATRIX_OVERHAUSER_NONUNI and BASIS_MATRIX_TMP.
//
//  Discussion:
//
//    YDATA(1:NDATA) = ( TDATA(1:NDATA) - 2 )**2 + 3
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 4
# define NDATA 4

  double alpha;
  double beta;
  int i;
  int j;
  int jhi;
  int left;
  char mark;
  double *mbasis;
  int nsample = 4;
  double tdata[NDATA];
  double thi;
  double tlo;
  double tval;
  double ydata[NDATA];
  double yval;

  printf ( "\n" );
  printf ( "TEST08\n" );
  printf ( "  BASIS_MATRIX_OVERHAUSER_NONUNI sets up the\n" );
  printf ( "    basis matrix for the nonuniform Overhauser\n" );
  printf ( "    spline.\n" );

  tdata[0] = 0.0;
  tdata[1] = 1.0;
  tdata[2] = 2.0;
  tdata[3] = 3.0;

  alpha = ( tdata[2] - tdata[1] ) / ( tdata[2] - tdata[0] );
  beta =  ( tdata[2] - tdata[1] ) / ( tdata[3] - tdata[1] );

  printf ( "\n" );
  printf ( "  ALPHA = %g\n", alpha );
  printf ( "  BETA  = %g\n", beta );

  mbasis = basis_matrix_overhauser_nonuni ( alpha, beta );

  for ( i = 0; i < NDATA; i++ )
  {
    ydata[i] = pow ( ( tdata[i] - 2.0 ), 2 ) + 3.0;
  }

  printf ( "\n" );
  printf ( "    TDATA, YDATA\n" );
  printf ( "\n" );
  for ( i = 0; i < NDATA; i++ )
  {
    printf ( "  %12g  %12g\n", tdata[i], ydata[i] );
  }

  left = 2;

  printf ( "\n" );
  printf ( "              T      Spline(T)\n" );
  printf ( "\n" );

  for ( i = 0; i <= NDATA; i++ )
  {
    if ( i == 0 )
    {
      tlo = tdata[0] - 0.5 * ( tdata[1] - tdata[0] );
      thi = tdata[0];
    }
    else if ( i < NDATA )
    {
      tlo = tdata[i-1];
      thi = tdata[i];
    }
    else if ( NDATA <= i )
    {
      tlo = tdata[NDATA-1];
      thi = tdata[NDATA-1] + 0.5 * ( tdata[NDATA-1] - tdata[NDATA-2] );
    }

    if ( i < NDATA )
    {
      jhi = nsample - 1;
    }
    else
    {
      jhi = nsample;
    }

    for ( j = 0; j <= jhi; j++ )
    {
      tval = ( ( double ) ( nsample - j ) * tlo
             + ( double ) (           j ) * thi )
             / ( double ) ( nsample     );

      yval = basis_matrix_tmp ( left, N, mbasis, NDATA, tdata, ydata, tval );

      if ( 0 < i & j == 0 )
      {
        mark = '*';
      }
      else
      {
        mark = ' ';
      }
      printf ( "  %c  %12g  %12g\n", mark, tval, yval );
    }
  }

  tdata[0] = 0.0;
  tdata[1] = 1.0;
  tdata[2] = 2.0;
  tdata[3] = 5.0;

  alpha = ( tdata[2] - tdata[1] ) / ( tdata[2] - tdata[0] );
  beta =  ( tdata[2] - tdata[1] ) / ( tdata[3] - tdata[1] );

  printf ( "\n" );
  printf ( "  ALPHA = %g\n", alpha );
  printf ( "  BETA  = %g\n", beta );

  mbasis = basis_matrix_overhauser_nonuni ( alpha, beta );

  for ( i = 0; i < NDATA; i++ )
  {
    ydata[i] = pow ( ( tdata[i] - 2.0 ), 2 ) + 3.0;
  }

  printf ( "\n" );
  printf ( "    TDATA, YDATA\n" );
  printf ( "\n" );
  for ( i = 0; i < NDATA; i++ )
  {
    printf ( "  %12g  %12g\n", tdata[i], ydata[i] );
  }

  left = 2;

  printf ( "\n" );
  printf ( "              T      Spline(T)\n" );
  printf ( "\n" );

  for ( i = 0; i <= NDATA; i++ )
  {
    if ( i == 0 )
    {
      tlo = tdata[0] - 0.5 * ( tdata[1] - tdata[0] );
      thi = tdata[0];
    }
    else if ( i < NDATA )
    {
      tlo = tdata[i-1];
      thi = tdata[i];
    }
    else if ( NDATA <= i )
    {
      tlo = tdata[NDATA-1];
      thi = tdata[NDATA-1] + 0.5 * ( tdata[NDATA-1] - tdata[NDATA-2] );
    }

    if ( i < NDATA )
    {
      jhi = nsample - 1;
    }
    else
    {
      jhi = nsample;
    }

    for ( j = 0; j <= jhi; j++ )
    {
      tval = ( ( double ) ( nsample - j ) * tlo
             + ( double ) (           j ) * thi )
             / ( double ) ( nsample     );

      yval = basis_matrix_tmp ( left, N, mbasis, NDATA, tdata, ydata, tval );

      if ( 0 < i & j == 0 )
      {
        mark = '*';
      }
      else
      {
        mark = ' ';
      }
      printf ( "  %c  %12g  %12g\n", mark, tval, yval );
    }
  }

  tdata[0] = 0.0;
  tdata[1] = 3.0;
  tdata[2] = 4.0;
  tdata[3] = 5.0;

  alpha = ( tdata[2] - tdata[1] ) / ( tdata[2] - tdata[0] );
  beta =  ( tdata[2] - tdata[1] ) / ( tdata[3] - tdata[1] );

  printf ( "\n" );
  printf ( "  ALPHA = %g\n", alpha );
  printf ( "  BETA  = %g\n", beta );

  free ( mbasis );
  mbasis = basis_matrix_overhauser_nonuni ( alpha, beta );

  for ( i = 0; i < NDATA; i++ )
  {
    ydata[i] = pow ( ( tdata[i] - 2.0 ), 2 ) + 3.0;
  }

  printf ( "\n" );
  printf ( "    TDATA, YDATA\n" );
  printf ( "\n" );
  for ( i = 0; i < NDATA; i++ )
  {
    printf ( "  %12g  %12g\n", tdata[i], ydata[i] );
  }

  left = 2;

  printf ( "\n" );
  printf ( "              T      Spline(T)\n" );
  printf ( "\n" );

  for ( i = 0; i <= NDATA; i++ )
  {
    if ( i == 0 )
    {
      tlo = tdata[0] - 0.5 * ( tdata[1] - tdata[0] );
      thi = tdata[0];
    }
    else if ( i < NDATA )
    {
      tlo = tdata[i-1];
      thi = tdata[i];
    }
    else if ( NDATA <= i )
    {
      tlo = tdata[NDATA-1];
      thi = tdata[NDATA-1] + 0.5 * ( tdata[NDATA-1] - tdata[NDATA-2] );
    }

    if ( i < NDATA )
    {
      jhi = nsample - 1;
    }
    else
    {
      jhi = nsample;
    }

    for ( j = 0; j <= jhi; j++ )
    {
      tval = ( ( double ) ( nsample - j ) * tlo
             + ( double ) (           j ) * thi )
             / ( double ) ( nsample     );

      yval = basis_matrix_tmp ( left, N, mbasis, NDATA, tdata, ydata, tval );

      if ( 0 < i & j == 0 )
      {
        mark = '*';
      }
      else
      {
        mark = ' ';
      }
      printf ( "  %c  %12g  %12g\n", mark, tval, yval );
    }
  }

  free ( mbasis );

  return;
# undef N
# undef NDATA
}
/******************************************************************************/

void test09 ( )

/******************************************************************************/
//
//  Purpose:
//
//    TEST09 tests BASIS_MATRIX_OVERHAUSER_NONUNI and BASIS_MATRIX_TMP.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 4
# define NDATA 4

  double alpha;
  double beta;
  int i;
  int j;
  int jhi;
  int left;
  char mark;
  double *mbasis;
  int nsample = 4;
  double tdata[NDATA];
  double thi;
  double tlo;
  double tval;
  double ydata[NDATA];
  double yval;

  printf ( "\n" );
  printf ( "TEST09\n" );
  printf ( "  BASIS_MATRIX_OVERHAUSER_NONUNI sets up the\n" );
  printf ( "    basis matrix for the nonuniform Overhauser \n" );
  printf ( "    spline.\n" );
  printf ( "\n" );
  printf ( "  First test that the nonuniform code can match\n" );
  printf ( "  the uniform code.  Compare these results with\n" );
  printf ( "  the uniform output.\n" );
  printf ( "\n" );

  tdata[0] = -1.0;
  tdata[1] =  0.0;
  tdata[2] =  1.0;
  tdata[3] =  2.0;

  alpha = ( tdata[2] - tdata[1] ) / ( tdata[2] - tdata[0] );
  beta =  ( tdata[2] - tdata[1] ) / ( tdata[3] - tdata[1] );

  printf ( "\n" );
  printf ( "  ALPHA = %g\n", alpha );
  printf ( "  BETA  = %g\n", beta );

  mbasis = basis_matrix_overhauser_nonuni ( alpha, beta );

  for ( i = 0; i < NDATA; i++ )
  {
    ydata[i] = pow ( ( tdata[i] + 2.0 ), 2 ) + 3.0;
  }

  printf ( "\n" );
  printf ( "    TDATA, YDATA\n" );
  printf ( "\n" );
  for ( i = 0; i < NDATA; i++ )
  {
    printf ( "  %12g  %12g\n", tdata[i], ydata[i] );
  }

  left = 2;

  printf ( "\n" );
  printf ( "              T      Spline(T)\n" );
  printf ( "\n" );

  for ( i = 0; i <= NDATA; i++ )
  {
    if ( i == 0 )
    {
      tlo = tdata[0] - 0.5 * ( tdata[1] - tdata[0] );
      thi = tdata[0];
    }
    else if ( i < NDATA )
    {
      tlo = tdata[i-1];
      thi = tdata[i];
    }
    else if ( NDATA <= i )
    {
      tlo = tdata[NDATA-1];
      thi = tdata[NDATA-1] + 0.5 * ( tdata[NDATA-1] - tdata[NDATA-2] );
    }

    if ( i < NDATA )
    {
      jhi = nsample - 1;
    }
    else
    {
      jhi = nsample;
    }

    for ( j = 0; j <= jhi; j++ )
    {
      tval = ( ( double ) ( nsample - j ) * tlo
             + ( double ) (           j ) * thi )
             / ( double ) ( nsample     );

      yval = basis_matrix_tmp ( left, N, mbasis, NDATA, tdata, ydata, tval );

      if ( 0 < i & j == 0 )
      {
        mark = '*';
      }
      else
      {
        mark = ' ';
      }
      printf ( "  %c  %12g  %12g\n", mark, tval, yval );
    }
  }

  printf ( "\n" );
  printf ( "  Now test that the nonuniform code on a\n" );
  printf ( "  nonuniform grid.\n" );
  printf ( "\n" );

  tdata[0] = -4.0;
  tdata[1] = -3.0;
  tdata[2] = -1.0;
  tdata[3] =  2.0;

  alpha = ( tdata[2] - tdata[1] ) / ( tdata[2] - tdata[0] );
  beta =  ( tdata[2] - tdata[1] ) / ( tdata[3] - tdata[1] );

  printf ( "\n" );
  printf ( "  ALPHA = %g\n", alpha );
  printf ( "  BETA  = %g\n", beta );

  free ( mbasis );
  mbasis = basis_matrix_overhauser_nonuni ( alpha, beta );

  for ( i = 0; i < NDATA; i++ )
  {
    ydata[i] = pow ( ( tdata[i] + 2.0 ), 2 ) + 3.0;
  }

  printf ( "\n" );
  printf ( "    TDATA, YDATA\n" );
  printf ( "\n" );
  for ( i = 0; i < NDATA; i++ )
  {
    printf ( "  %12g  %12g\n", tdata[i], ydata[i] );
  }

  left = 2;

  printf ( "\n" );
  printf ( "              T      Spline(T)\n" );
  printf ( "\n" );

  for ( i = 0; i <= NDATA; i++ )
  {
    if ( i == 0 )
    {
      tlo = tdata[0] - 0.5 * ( tdata[1] - tdata[0] );
      thi = tdata[0];
    }
    else if ( i < NDATA )
    {
      tlo = tdata[i-1];
      thi = tdata[i];
    }
    else if ( NDATA <= i )
    {
      tlo = tdata[NDATA-1];
      thi = tdata[NDATA-1] + 0.5 * ( tdata[NDATA-1] - tdata[NDATA-2] );
    }

    if ( i < NDATA )
    {
      jhi = nsample - 1;
    }
    else
    {
      jhi = nsample;
    }

    for ( j = 0; j <= jhi; j++ )
    {
      tval = ( ( double ) ( nsample - j ) * tlo
             + ( double ) (           j ) * thi )
             / ( double ) ( nsample     );

      yval = basis_matrix_tmp ( left, N, mbasis, NDATA, tdata, ydata, tval );

      if ( 0 < i & j == 0 )
      {
        mark = '*';
      }
      else
      {
        mark = ' ';
      }
      printf ( "  %c  %12g  %12g\n", mark, tval, yval );
    }
  }

  free ( mbasis );

  return;
# undef N
# undef NDATA
}
/******************************************************************************/

void test10 ( )

/******************************************************************************/
//
//  Purpose:
//
//    TEST10 tests BC_VAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 2

  int i;
  int nsample = 101;
  double t;
  double xcon[N+1] = { 0.0, 0.75, 1.0 };
  double xval;
  double ycon[N+1] = { 1.0, 0.0,  1.0 };
  double yval;

  printf ( "\n" );
  printf ( "TEST10\n" );
  printf ( "  BC_VAL evaluates a general Bezier function.\n" );
//
//  One point on the curve should be about (0.75, 0.536).
//
  printf ( "\n" );
  printf ( "        T             X(T)           Y(T)\n" );
  printf ( "\n" );

  for ( i = 1; i <= nsample; i++ )
  {
    t = ( double ) ( i - 1 ) / ( double ) ( nsample - 1 );

    bc_val ( N, t, xcon, ycon, &xval, &yval );

    printf ( "  %12g  %12g  %12g\n", t, xval, yval );
  }

  printf ( "\n" );
  printf ( "  The point ( 0.75, 0.536 ) should be on the curve.\n" );

  return;
# undef N
}
/******************************************************************************/

void test11 ( )

/******************************************************************************/
//
//  Purpose:
//
//    TEST11 tests BEZ_VAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 2

  double a = 0.0;
  double b = 1.0;
  double bval;
  int i;
  int nsample = 21;
  double x;
  double y[N+1] = { 1.0, 0.0, 1.0 };

  printf ( "\n" );
  printf ( "TEST11\n" );
  printf ( "  BEZ_VAL evaluates a Bezier function.\n" );
//
//  One point on the curve should be (0.75, 20/32).
//
  printf ( "\n" );
  printf ( "        I             X           B\n" );
  printf ( "\n" );

  for ( i = 1; i <= nsample; i++ )
  {
    x = ( ( double ) ( nsample - i     ) * a
        + ( double ) (           i - 1 ) * b )
        / ( double ) ( nsample     - 1 );

    bval = bez_val ( N, x, a, b, y );

    printf ( "  %6d  %12g  %12g\n", i, x, bval );
  }

  printf ( "\n" );
  printf ( "  When X = 0.75\n" );
  printf ( "  BEZ_VAL(X) should be 0.625\n" );

  return;
# undef N
}
/******************************************************************************/

void test115 ( )

/******************************************************************************/
//
//  Purpose:
//
//    TEST115 tests BP01.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N_MAX 3

  double a = 0.0;
  double b = 1.0;
  double *bern;
  int i;
  int j;
  int n;
  int nsample = 11;
  double x;

  printf ( "\n" );
  printf ( "TEST115\n" );
  printf ( "  BP01 evaluates the Bernstein basis polynomials\n" );
  printf ( "  for the interval [0,1].\n" );

  for ( n = 0; n <= N_MAX; n++ )
  {
    printf ( "\n" );
    printf ( "  Degree N = %d\n", n );
    printf ( "\n" );
    printf ( "   X         BERN(N,0,X)  BERN(N,1,X)  BERN(N,2,X)  BERN(N,3,X)\n" );
    printf ( "\n" );

    for ( i = 1; i <= nsample; i++ )
    {
      x = ( ( double ) ( nsample - i     ) * a
          + ( double ) (           i - 1 ) * b )
          / ( double ) ( nsample     - 1 );

      bern = bp01 ( n, x );

      printf ( "  %10g  ", x );
      for ( j = 0; j <= n; j++ )
      {
        printf ( "%12g  ", bern[j] );
      }
      printf ( "\n" );

      free ( bern );
    }
  }

  return;
# undef N_MAX
}
/******************************************************************************/

void test116 ( )

/******************************************************************************/
//
//  Purpose:
//
//    TEST116 tests BPAB.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N_MAX 3

  double a = 1.0;
  double b = 3.0;
  double *bern;
  int i;
  int j;
  int n;
  int nsample = 11;
  double x;

  printf ( "\n" );
  printf ( "TEST116\n" );
  printf ( "  BPAB evaluates the Bernstein basis polynomials\n" );
  printf ( "  for the interval [A,B].\n" );

  for ( n = 0; n <= N_MAX; n++ )
  {
    printf ( "\n" );
    printf ( "  Degree N = %d\n", n );
    printf ( "\n" );
    printf ( "   X         BERN(N,0,X)  BERN(N,1,X)  BERN(N,2,X)  BERN(N,3,X)\n" );
    printf ( "\n" );

    for ( i = 1; i <= nsample; i++ )
    {
      x = ( ( double ) ( nsample - i     ) * a
          + ( double ) (           i - 1 ) * b )
          / ( double ) ( nsample     - 1 );

      bern = bpab ( n, a, b, x );

      printf ( "  %10g  ", x );
      for ( j = 0; j <= n; j++ )
      {
        printf ( "%12g  ", bern[j] );
      }
      printf ( "\n" );

      free ( bern );
    }
  }

  return;
# undef N_MAX
}
/******************************************************************************/

void test12 ( )

/******************************************************************************/
//
//  Purpose:
//
//    TEST12 tests BPAB_APPROX.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define MAXDATA 10

  double a;
  double b;
  int i;
  int ndata;
  int nsample;
  double xdata[MAXDATA+1];
  double xval;
  double ydata[MAXDATA+1];
  double yval;

  printf ( "\n" );
  printf ( "TEST12\n" );
  printf ( "  BPAB_APPROX evaluates the Bernstein polynomial\n" );
  printf ( "  approximant to a function F(X).\n" );

  a = 1.0;
  b = 3.0;

  for ( ndata = 0; ndata <= 9; ndata = ndata + 3 )
  {
    for ( i = 0; i <= ndata; i++ )
    {
      if ( ndata == 0 )
      {
        xdata[i] = 0.5 * ( a + b );
      }
      else
      {
        xdata[i] = ( ( double ) ( ndata - i ) * a
                   + ( double ) (         i ) * b )
                   / ( double ) ( ndata );
      }

      ydata[i] = sin ( xdata[i] );

    }

    printf ( "\n" );
    printf ( "    XDATA    YDATA\n" );
    printf ( "\n" );
    for ( i = 0; i <= ndata; i++ )
    {
      printf ( "  %12g  %12g\n", xdata[i], ydata[i] );
    }

    printf ( "\n" );
    printf ( "  Bernstein approximant of degree N = %d\n", ndata );
    printf ( "\n" );
    printf ( "    X      F(X)     BERN(X)    ERROR\n" );
    printf ( "\n" );

    nsample = 2 * ndata + 1;

    for ( i = 1; i<= nsample; i++ )
    {
      if ( nsample == 1 )
      {
        xval = 0.5 * ( a + b );
      }
      else
      {
        xval = ( ( double ) ( nsample - i     ) * a
               + ( double ) (           i - 1 ) * b )
               / ( double ) ( nsample     - 1 );
      }

      yval = bpab_approx ( ndata, a, b, ydata, xval );

      printf ( "  %12g  %12g  %12g  %12g\n", xval, sin ( xval ), yval, yval - sin ( xval ) );
    }
  }

  return;
# undef MAXDATA
}
/******************************************************************************/

void test125 ( )

/******************************************************************************/
//
//  Purpose:
//
//    TEST125 tests LEAST_SET_OLD and LEAST_VAL_OLD.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define MAXDEG 6
# define NTAB 21

  double b[MAXDEG];
  double c[MAXDEG+1];
  double d[MAXDEG-1];
  double eps;
  double error;
  int i;
  int ierror;
  int j;
  int jhi;
  int ndeg;
  double ptab[NTAB];
  double xtab[NTAB];
  double xval;
  double ytab[NTAB];
  double ytrue;
  double yval;

  printf ( "\n" );
  printf ( "TEST125\n" );
  printf ( "  LEAST_SET_OLD sets a least squares polynomial,\n" );
  printf ( "  LEAST_VAL_OLD evaluates it.\n" );

  for ( i = 0; i < NTAB; i++ )
  {
    xtab[i] = ( ( double ) ( NTAB - i - 1 ) * ( -1.0 )
              + ( double ) (        i     ) * ( +1.0 ) )
              / ( double ) ( NTAB     - 1 );
    ytab[i] = ( double ) ( ( int ) ( exp ( xtab[i] ) * 100.0 + 0.5 ) )
      / 100.0;
  }

  printf ( "\n" );
  printf ( "  The data to be interpolated:\n" );
  printf ( "\n" );
  printf ( "  Number of data values = %d\n", NTAB );
  printf ( "\n" );
  printf ( "       X             Y\n" );
  printf ( "\n" );

    printf ( "\n" );
    printf ( "    XTAB    YTAB\n" );
    printf ( "\n" );
    for ( i = 0; i < NTAB; i++ )
    {
      printf ( "  %12g  %12g\n", xtab[i], ytab[i] );
    }

  for ( ndeg = 1; ndeg <= MAXDEG; ndeg++ )
  {
    printf ( "\n" );
    printf ( "  Use a polynomial of degree: %d\n", ndeg );
    printf ( "\n" );

    least_set_old ( NTAB, xtab, ytab, ndeg, ptab, b, c, d, &eps, &ierror );

    printf ( "\n" );
    printf ( "  Total approximation error = %g\n", eps );
    printf ( "\n" );
    printf ( "       X            F(X)          P(X)          Error\n" );
    printf ( "\n" );

    for ( i = 1; i <= NTAB; i++ )
    {
      if ( i < NTAB )
      {
        jhi = 2;
      }
      else
      {
        jhi = 0;
      }

      for ( j = 0; j <= jhi; j++ )
      {
        if ( i < NTAB )
        {
          xval = ( ( double ) ( 3 - j ) * xtab[i-1]
                 + ( double ) (     j ) * xtab[i]   )
                 / ( double ) ( 3     );
        }
        else
        {
          xval = xtab[i-1];
        }

        yval = least_val_old ( xval, ndeg, b, c, d );

        ytrue = ( double ) ( ( int ) ( exp ( xval ) * 100.0 + 0.5 ) )
          / 100.0;

        error = yval - ytrue;

        printf ( "  %12g  %12g  %12g  %12g\n", xval, yval, ytrue, error );
      }
    }
  }

  return;
# undef MAXDEG
}
/******************************************************************************/

void test126 ( )

/******************************************************************************/
//
//  Purpose:
//
//    TEST126 tests LEAST_SET and LEAST_VAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define POINT_NUM 21
# define NTERMS 4

  double b[NTERMS];
  double c[NTERMS];
  double d[NTERMS];
  double f[POINT_NUM];
  double fp[POINT_NUM];
  int i;
  int nterms2;
  double px;
  double w[POINT_NUM];
  double x[POINT_NUM];

  for ( i = 0; i < POINT_NUM; i++ )
  {
    w[i] = 1.0;
    x[i] = -1.0 + ( double ) ( i ) / 10.0;
    f[i] = x[i] * x[i] - x[i] - 6.0;
    fp[i] = 2.0 * x[i] - 1.0;
  }

  least_set ( POINT_NUM, x, f, w, NTERMS, b, c, d );

  printf ( "\n" );
  printf ( "TEST126\n" );
  printf ( "  LEAST_SET sets a least squares polynomial,\n" );
  printf ( "  LEAST_VAL evaluates it.\n" );
  printf ( "\n" );
  printf ( "  X, F(X), P(X), Error\n" );
  printf ( "\n" );
  for ( nterms2 = 1; nterms2 <= NTERMS; nterms2++ )
  {
    printf ( "\n" );
    printf ( "  Using polynomial order = %d\n", nterms2 );
    printf ( "\n" );
    for ( i = 0; i < POINT_NUM; i++ )
    {
      px = least_val ( nterms2, b, c, d, x[i] );
      printf ( "  %12g  %12g  %12g  %12g\n", x[i], f[i], px, px - f[i] );
    }
  }

  return;
# undef POINT_NUM
# undef NTERMS
}
/******************************************************************************/

void test127 ( )

/******************************************************************************/
//
//  Purpose:
//
//    TEST127 tests LEAST_SET and LEAST_VAL2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define POINT_NUM 21
# define NTERMS 4

  double b[NTERMS];
  double c[NTERMS];
  double d[NTERMS];
  double f[POINT_NUM];
  double fp[POINT_NUM];
  int i;
  int nterms2;
  double px;
  double pxp;
  double w[POINT_NUM];
  double x[POINT_NUM];

  for ( i = 0; i < POINT_NUM; i++ )
  {
    w[i] = 1.0;
    x[i] = -1.0 + ( double ) ( i ) / 10.0;
    f[i] = x[i] * x[i] - x[i] - 6.0;
    fp[i] = 2.0 * x[i] - 1.0;
  }

  least_set ( POINT_NUM, x, f, w, NTERMS, b, c, d );

  printf ( "\n" );
  printf ( "TEST127\n" );
  printf ( "  LEAST_SET sets a least squares polynomial,\n" );
  printf ( "  LEAST_VAL2 evaluates it.\n" );
  printf ( "\n" );
  printf ( "  X, F(X), P(X), FP(X), PP(X)\n" );
  printf ( "\n" );

  for ( nterms2 = 1; nterms2 <= NTERMS; nterms2++ )
  {
    printf ( "\n" );
    printf ( "  Using polynomial order = %d\n", nterms2 );
    printf ( "\n" );
    for ( i = 0; i < POINT_NUM; i++ )
    {
      least_val2 ( nterms2, b, c, d, x[i], &px, &pxp );
      printf ( "  %12g  %12g  %12g  %12g  %12g\n", x[i], f[i], px, fp[i], pxp );
    }
  }

  return;
# undef POINT_NUM
# undef NTERMS
}
/******************************************************************************/

void test13 ( )

/******************************************************************************/
//
//  Purpose:
//
//    TEST13 tests SPLINE_B_VAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define NDATA 11

  int i;
  int j;
  int jhi;
  char mark;
  int nsample = 4;
  double pi = 3.141592653589793;
  double tdata[NDATA];
  double thi;
  double tlo;
  double tval;
  double ydata[NDATA];
  double yval;

  printf ( "\n" );
  printf ( "TEST13\n" );
  printf ( "  SPLINE_B_VAL evaluates the\n" );
  printf ( "    B spline.\n" );
  printf ( "\n" );
  printf ( "  TDATA   YDATA\n" );
  printf ( "\n" );

  for ( i = 0; i < NDATA; i++ )
  {
    tdata[i] = ( double ) i;
    ydata[i] = sin ( 2.0 * pi * tdata[i] / ( double ) ( NDATA - 1) );
    printf ( "  %12g  %12g\n", tdata[i], ydata[i] );
  }

  printf ( "\n" );
  printf ( "    T, Spline(T)\n" );
  printf ( "\n" );

  for ( i = 0; i <= NDATA; i++ )
  {
    if ( i == 0 )
    {
      tlo = tdata[1] - 0.5 * ( tdata[1] - tdata[0] );
      thi = tdata[0];
    }
    else if ( i < NDATA )
    {
      tlo = tdata[i-1];
      thi = tdata[i];
    }
    else if ( NDATA <= i )
    {
      tlo = tdata[NDATA-1];
      thi = tdata[NDATA-1] + 0.5 * ( tdata[NDATA-1] - tdata[NDATA-2] );
    }

    if ( i < NDATA )
    {
      jhi = nsample - 1;
    }
    else
    {
      jhi = nsample;
    }

    for ( j = 0; j <= jhi; j++ )
    {
      tval = ( ( double ) ( nsample - j ) * tlo
             + ( double ) (           j ) * thi )
             / ( double ) ( nsample     );

      yval = spline_b_val ( NDATA, tdata, ydata, tval );

      if ( 0 < i && j == 0 )
      {
        mark = '*';
      }
      else
      {
        mark = ' ';
      }
      printf ( "  %c  %12g  %12g\n", mark, tval, yval );

    }

  }

  return;
# undef NDATA
}
/******************************************************************************/

void test14 ( )

/******************************************************************************/
//
//  Purpose:
//
//    TEST14 tests SPLINE_BETA_VAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define NDATA 11

  double beta1;
  double beta2;
  int i;
  int j;
  int jhi;
  char mark;
  int nsample = 4;
  double pi = 3.141592653589793;
  double tdata[NDATA];
  double thi;
  double tlo;
  double tval;
  double ydata[NDATA];
  double yval;

  printf ( "\n" );
  printf ( "TEST14\n" );
  printf ( "  SPLINE_BETA_VAL evaluates the BETA spline.\n" );
  printf ( "\n" );
  printf ( "       TDATA         YDATA\n" );
  printf ( "\n" );
  for ( i = 0; i < NDATA; i++ )
  {
    tdata[i] = ( double ) i;
    ydata[i] = sin ( 2.0 * pi * tdata[i] / ( double ) ( NDATA - 1) );
    printf ( "  %12g  %12g\n", tdata[i], ydata[i] );
  }

  beta1 = 1.0;
  beta2 = 0.0;

  printf ( "\n" );
  printf ( "  BETA1 = %g\n", beta1 );
  printf ( "  BETA2 = %g\n", beta2 );
  printf ( "\n" );
  printf ( "    T, Spline(T)\n" );
  printf ( "\n" );

  for ( i = 0; i <= NDATA; i++ )
  {

    if ( i == 0 )
    {
      tlo = tdata[0] - 0.5 * ( tdata[1] - tdata[0] );
      thi = tdata[0];
    }
    else if ( i < NDATA )
    {
      tlo = tdata[i-1];
      thi = tdata[i];
    }
    else if ( NDATA <= i )
    {
      tlo = tdata[NDATA-1];
      thi = tdata[NDATA-1] + 0.5 * ( tdata[NDATA-1] - tdata[NDATA-2] );
    }

    if ( i < NDATA )
    {
      jhi = nsample - 1;
    }
    else
    {
      jhi = nsample;
    }

    for ( j = 0; j <= jhi; j++ )
    {
      tval = ( ( double ) ( nsample - j ) * tlo
             + ( double ) (           j ) * thi )
             / ( double ) ( nsample     );

      yval = spline_beta_val ( beta1, beta2, NDATA, tdata, ydata, tval );

      if ( 0 < i && j == 0 )
      {
        mark = '*';
      }
      else
      {
        mark = ' ';
      }
      printf ( "  %c  %12g  %12g\n", mark, tval, yval );
    }
  }

  beta1 = 1.0;
  beta2 = 100.0;

  printf ( "\n" );
   printf ( "  BETA1 = %g\n", beta1 );
  printf ( "  BETA2 = %g\n", beta2 );
  printf ( "\n" );
  printf ( "    T, Spline(T)\n" );
  printf ( "\n" );

  for ( i = 0; i <= NDATA; i++ )
  {
    if ( i == 0 )
    {
      tlo = tdata[0] - 0.5 * ( tdata[1] - tdata[0] );
      thi = tdata[0];
    }
    else if ( i < NDATA )
    {
      tlo = tdata[i-1];
      thi = tdata[i];
    }
    else if ( NDATA <= i )
    {
      tlo = tdata[NDATA-1];
      thi = tdata[NDATA-1] + 0.5 * ( tdata[NDATA-1] - tdata[NDATA-2] );
    }

    if ( i < NDATA )
    {
      jhi = nsample - 1;
    }
    else
    {
      jhi = nsample;
    }

    for ( j = 0; j <= jhi; j++ )
    {
      tval = ( ( double ) ( nsample - j ) * tlo
             + ( double ) (           j ) * thi )
             / ( double ) ( nsample     );

      yval = spline_beta_val ( beta1, beta2, NDATA, tdata, ydata, tval );

      if ( 0 < i && j == 0 )
      {
        mark = '*';
      }
      else
      {
        mark = ' ';
      }
      printf ( "  %c  %12g  %12g\n", mark, tval, yval );
    }

  }

  beta1 = 100.0;
  beta2 = 0.0;

  printf ( "\n" );
  printf ( "  BETA1 = %g\n", beta1 );
  printf ( "  BETA2 = %g\n", beta2 );
  printf ( "\n" );
  printf ( "    T, Spline(T)\n" );
  printf ( "\n" );

  for ( i = 0; i <= NDATA; i++ )
  {
    if ( i == 0 )
    {
      tlo = tdata[0] - 0.5 * ( tdata[1] - tdata[0] );
      thi = tdata[0];
    }
    else if ( i < NDATA )
    {
      tlo = tdata[i-1];
      thi = tdata[i];
    }
    else if ( NDATA <= i )
    {
      tlo = tdata[NDATA-1];
      thi = tdata[NDATA-1] + 0.5 * ( tdata[NDATA-1] - tdata[NDATA-2] );
    }

    if ( i < NDATA )
    {
      jhi = nsample - 1;
    }
    else
    {
      jhi = nsample;
    }

    for ( j = 0; j <= jhi; j++ )
    {
      tval = ( ( double ) ( nsample - j ) * tlo
             + ( double ) (           j ) * thi )
             / ( double ) ( nsample     );

      yval = spline_beta_val ( beta1, beta2, NDATA, tdata, ydata, tval );

      if ( 0 < i && j == 0 )
      {
        mark = '*';
      }
      else
      {
        mark = ' ';
      }
      printf ( "  %c  %12g  %12g\n", mark, tval, yval );
    }
  }

  return;
# undef NDATA
}
/******************************************************************************/

void test145 ( )

/******************************************************************************/
//
//  Purpose:
//
//    TEST145 tests SPLINE_CONSTANT_VAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define NDATA 12
# define N_TEST 20

  double fval;
  int i;
  int j;
  int seed = 123456789;
  double *tdata;
  double thi;
  double tlo;
  double *t_test;
  double tval;
  double ydata[NDATA];
  double yval;

  printf ( "\n" );
  printf ( "TEST145\n" );
  printf ( "  SPLINE_CONSTANT_VAL evaluates a piecewise\n" );
  printf ( "  constant spline.\n" );
  printf ( "\n" );
  printf ( "  Runge's function, evenly spaced knots.\n" );
//
//  Set the data.
//
  tlo = -1.0;
  thi = +1.0;

  tdata = r8vec_even_new ( NDATA-1, tlo, thi );

  for ( i = 0; i < NDATA; i++ )
  {
    if ( i == 0 )
    {
      tval = tdata[0];
    }
    else if ( i < NDATA - 1 )
    {
      tval = 0.5 * ( tdata[i-1] + tdata[i] );
    }
    else if ( i == NDATA - 1 )
    {
      tval = tdata[i-1];
    }

    ydata[i] = frunge ( tval );

  }

  printf ( "\n" );
  printf ( "  The data to be interpolated:\n" );
  printf ( "\n" );
  printf ( "  Number of data values = %d\n", NDATA );
  printf ( "\n" );
  printf ( "       T             Y\n" );
  printf ( "\n" );

  for ( i = 0; i < NDATA; i++ )
  {
    printf ( "  *              %12g\n", ydata[i] );
    if ( i < NDATA -1 )
    {
      printf ( "  *%12g\n", tdata[i] );
    }
  }
//
//  Sample the spline.
//
  t_test = r8vec_uniform_new ( N_TEST, tlo-1.0, thi+1.0, &seed );

  r8vec_sort_bubble_a ( N_TEST, t_test );

  printf ( "\n" );
  printf ( "     T     Y(interp)    Y(exact)\n" );
  printf ( "\n" );

  j = 0;
  printf ( "  *              %12g\n", ydata[j] );
  j = j + 1;

  for ( i = 0; i < N_TEST; i++ )
  {
    tval = t_test[i];

    yval = spline_constant_val ( NDATA, tdata, ydata, tval );

    if ( j <= NDATA - 1 )
    {
      while ( tdata[j-1] <= tval )
      {
        fval = frunge ( tdata[j-1] );
        printf ( "  *%12g                %12g\n", tdata[j-1], fval );
        printf ( "  *              %12g\n", ydata[j] );
        j = j + 1;
        if ( NDATA <= j )
        {
          break;
        }
      }
    }

    fval = frunge ( tval );

    printf ( "   %12g  %12g  %12g\n", tval, yval, fval );
  }

  free ( tdata );
  free ( t_test );

  return;
# undef NDATA
# undef N_TEST
}
/******************************************************************************/

void test15 ( )

/******************************************************************************/
//
//  Purpose:
//
//    TEST15 tests SPLINE_CUBIC_SET and SPLINE_CUBIC_VAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 11

  int i;
  int ibcbeg;
  int ibcend;
  int j;
  int jhi;
  int k;
  double t[N];
  double tval;
  double y[N];
  double ybcbeg;
  double ybcend;
  double *ypp;
  double yppval;
  double ypval;
  double yval;
//
//  Set up the data.
//
  printf ( "\n" );
  printf ( "TEST15\n" );
  printf ( "  SPLINE_CUBIC_SET sets up a cubic spline;\n" );
  printf ( "  SPLINE_CUBIC_VAL evaluates it.\n" );
  printf ( "\n" );
  printf ( "  Runge's function, evenly spaced knots.\n" );
  printf ( "\n" );
  printf ( "     I     T         Y\n" );
  printf ( "\n" );

  for ( i = 0; i < N; i++ )
  {
    t[i] =  ( ( double ) ( N - i     ) * (-1.0)
            + ( double ) (     i - 1 ) * (+1.0) )
            / ( double ) ( N     - 1 );
    y[i] =  frunge ( t[i] );
    printf ( "  %6d  %12g  %12g\n", i, t[i], y[i] );
  }
//
//  Try all three types of boundary condition.
//
  for ( k = 0; k < 3; k++ )
  {
    if ( k == 0 )
    {
      ibcbeg = 0;
      ybcbeg = 0.0;

      ibcend = 0;
      ybcend = 0.0;

      printf ( "\n" );
      printf ( "  Boundary condition 0 at both ends:\n" );
      printf ( "  Spline is quadratic in boundary intervals.\n" );
    }
    else if ( k == 1 )
    {
      ibcbeg = 1;
      ybcbeg = fprunge ( t[0] );

      ibcend = 1;
      ybcend = fprunge ( t[N-1] );

      printf ( "\n" );
      printf ( "  Boundary condition 1 at both ends:\n" );
      printf ( "  Y'(left) =  %g\n", ybcbeg );
      printf ( "  Y'(right) = %g\n", ybcend );

    }
    else if ( k == 2 )
    {
      ibcbeg = 2;
      ybcbeg = fpprunge ( t[0] );

      ibcend = 2;
      ybcend = fpprunge ( t[N-1] );

      printf ( "\n" );
      printf ( "  Boundary condition 2 at both ends:\n" );
      printf ( "  YP''(left) =  %g\n", ybcbeg );
      printf ( "  YP''(right) = %g\n", ybcend );
    }
    else if ( k == 3 )
    {
      ibcbeg = 2;
      ybcbeg = 0.0;

      ibcend = 2;
      ybcend = 0.0;

      printf ( "\n" );
      printf ( "  Natural spline:\n" );
      printf ( "  YP''(left) =  %g\n", ybcbeg );
      printf ( "  YP''(right) = %g\n", ybcend );
    }

    ypp = spline_cubic_set ( N, t, y, ibcbeg, ybcbeg, ibcend, ybcend );

    printf ( "\n" );
    printf ( "  SPLINE''(T), F''(T):\n" );
    printf ( "\n" );
    for ( i = 0; i < N; i++ )
    {
      printf ( "%10g  %10g\n", ypp[i], fpprunge ( t[i] ) );
    }

    printf ( "\n" );
    printf ( "  T, SPLINE(T), F(T)\n" );
    printf ( "\n" );

    for ( i = 0; i <= N; i++ )
    {
      if ( i == 0 )
      {
        jhi = 1;
      }
      else if ( i < N )
      {
        jhi = 2;
      }
      else
      {
        jhi = 2;
      }

      for ( j = 1; j <= jhi; j++ )
      {
        if ( i == 0 )
        {
          tval = t[0] - 1.0;
        }
        else if ( i < N )
        {
          tval = (
              ( double ) ( jhi - j + 1 ) * t[i-1]
            + ( double ) (       j - 1 ) * t[i] )
            / ( double ) ( jhi         );
        }
        else
        {
          if ( j == 1 )
          {
            tval = t[N-1];
          }
          else
          {
            tval = t[N-1] + 1.0;
          }
        }

        yval = spline_cubic_val ( N, t, y, ypp, tval, &ypval, &yppval );

        printf ( "%10g  %10g  %10g\n", tval, yval, frunge ( tval ) );
      }
    }
    free ( ypp );
  }

  return;
# undef N
}
/******************************************************************************/

void test16 ( )

/******************************************************************************/
//
//  Purpose:
//
//    TEST16 tests SPLINE_CUBIC_SET and SPLINE_CUBIC_VAL2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 11

  int i;
  int ibcbeg;
  int ibcend;
  int j;
  int jhi;
  int k;
  int left;
  int left_in;
  double t[N];
  double tval;
  double y[N];
  double ybcbeg;
  double ybcend;
  double *ypp;
  double yppval;
  double ypval;
  double yval;
//
//  Set up the data.
//
  printf ( "\n" );
  printf ( "TEST16\n" );
  printf ( "  SPLINE_CUBIC_SET sets up a cubic spline;\n" );
  printf ( "  SPLINE_CUBIC_VAL2 evaluates it.\n" );
  printf ( "\n" );
  printf ( "  Runge's function, evenly spaced knots.\n" );
  printf ( "\n" );
  printf ( "     I      T       Y\n" );
  printf ( "\n" );

  for ( i = 0; i < N; i++ )
  {
    t[i] =  ( ( double ) ( N - i     ) * (-1.0)
            + ( double ) (     i - 1 ) * (+1.0) )
            / ( double ) ( N     - 1 );
    y[i] =  frunge ( t[i] );
    printf ( "  %6d  %12g  %12g\n", i, t[i], y[i] );
  }
//
//  Try all three types of boundary condition.
//
  for ( k = 0; k < 3; k++ )
  {
    if ( k == 0 )
    {
      ibcbeg = 0;
      ybcbeg = 0.0;

      ibcend = 0;
      ybcend = 0.0;

      printf ( "\n" );
      printf ( "  Boundary condition 0 at both ends:\n" );
      printf ( "  Spline is quadratic in boundary intervals.\n" );
    }
    else if ( k == 1 )
    {
      ibcbeg = 1;
      ybcbeg = fprunge ( t[0] );

      ibcend = 1;
      ybcend = fprunge ( t[N-1] );

      printf ( "\n" );
      printf ( "  Boundary condition 1 at both ends:\n" );
      printf ( "  Y'(left) =  %g\n", ybcbeg );
      printf ( "  Y'(right) = %g\n", ybcend );
    }
    else if ( k == 2 )
    {
      ibcbeg = 2;
      ybcbeg = fpprunge ( t[0] );

      ibcend = 2;
      ybcend = fpprunge ( t[N-1] );

      printf ( "\n" );
      printf ( "  Boundary condition 2 at both ends:\n" );
      printf ( "  YP''(left) =  %g\n", ybcbeg );
      printf ( "  YP''(right) = %g\n", ybcend );
    }

    ypp = spline_cubic_set ( N, t, y, ibcbeg, ybcbeg, ibcend, ybcend );

    printf ( "\n" );
    printf ( "  SPLINE''(T), F''(T):\n" );
    printf ( "\n" );
    for ( i = 0; i < N; i++ )
    {
      printf ( "%12g  %12g\n", ypp[i], fpprunge(t[i]) );
    }

    left = 0;

    printf ( "\n" );
    printf ( "  T, SPLINE(T), F(T), LEFT_IN, LEFT_OUT\n" );
    printf ( "\n" );

    for ( i = 0; i <= N; i++ )
    {
      if ( i == 0 )
      {
        jhi = 1;
      }
      else if ( i < N )
      {
        jhi = 2;
      }
      else
      {
        jhi = 2;
      }

      for ( j = 1; j <= jhi; j++ )
      {
        if ( i == 0 )
        {
          tval = t[0] - 1.0;
        }
        else if ( i < N )
        {
          tval = ( ( double ) ( jhi - j + 1 ) * t[i-1]
                 + ( double ) (       j - 1 ) * t[i] )
                 / ( double ) ( jhi         );
        }
        else
        {
          if ( j == 1 )
          {
            tval = t[N-1];
          }
          else
          {
            tval = t[N-1] + 1.0;
          }
        }

        left_in = left;
        spline_cubic_val2 ( N, t, tval, &left, y, ypp, &yval, &ypval,
          &yppval );

        printf ( "%12g  %12g  %12g  %6d  %6d\n", tval, yval, frunge ( tval ), left_in, left );
      }
    }
    free ( ypp );
  }

  return;

# undef N
}
/******************************************************************************/

void test17 ( )

/******************************************************************************/
//
//  Purpose:
//
//    TEST17 tests SPLINE_CUBIC_SET and SPLINE_CUBIC_VAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 11

  int i;
  int ibcbeg;
  int ibcend;
  int j;
  int jhi;
  int k;
  double t[N];
  double tval;
  double y[N];
  double ybcbeg;
  double ybcend;
  double *ypp;
  double yppval;
  double ypval;
  double yval;
//
//  Set up the data.
//
  printf ( "\n" );
  printf ( "TEST17\n" );
  printf ( "  SPLINE_CUBIC_SET sets up a cubic spline;\n" );
  printf ( "  SPLINE_CUBIC_VAL evaluates it.\n" );
  printf ( "\n" );
  printf ( "  Cubic function, unevenly spaced knots.\n" );
  printf ( "\n" );
  printf ( "  T, Y\n" );
  printf ( "\n" );

  for ( i = 0; i < N; i++ )
  {
    t[i] =  ( double ) ( i ) / ( double ) ( N - 1 );
    t[i] = t[i] * t[i];
    y[i] =  fcube ( t[i] );
    printf ( "  %6d  %12g  %12g\n", i, t[i], y[i] );
  }
//
//  Try all three types of boundary condition.
//
  for ( k = 0; k < 3; k++ )
  {
    if ( k == 0 )
    {
      ibcbeg = 0;
      ybcbeg = 0.0;

      ibcend = 0;
      ybcend = 0.0;

      printf ( "\n" );
      printf ( "  Boundary condition 0 at both ends:\n" );
      printf ( "  Spline is quadratic in boundary intervals.\n" );

    }
    else if ( k == 1 )
    {
      ibcbeg = 1;
      ybcbeg = fpcube ( t[0] );

      ibcend = 1;
      ybcend = fpcube ( t[N-1] );

      printf ( "\n" );
      printf ( "  Boundary condition 1 at both ends:\n" );
      printf ( "  Y'(left) =  %g\n", ybcbeg );
      printf ( "  Y'(right) = %g\n", ybcend );

    }
    else if ( k == 2 )
    {
      ibcbeg = 2;
      ybcbeg = fppcube ( t[0] );

      ibcend = 2;
      ybcend = fppcube ( t[N-1] );

      printf ( "\n" );
      printf ( "  Boundary condition 2 at both ends:\n" );
      printf ( "  YP''(left) =  %g\n", ybcbeg );
      printf ( "  YP''(right) = %g\n", ybcend );

    }

    ypp = spline_cubic_set ( N, t, y, ibcbeg, ybcbeg, ibcend, ybcend );

    printf ( "\n" );
    printf ( "  SPLINE''(T), F''(T):\n" );
    printf ( "\n" );
    for ( i = 0; i < N; i++ )
    {
      printf ( "%12g  %12g\n", ypp[i], fppcube(t[i]) );
    }

    printf ( "\n" );
    printf ( "  T, SPLINE(T), F(T)\n" );
    printf ( "\n" );

    for ( i = 0; i <= N; i++ )
    {
      if ( i == 0 )
      {
        jhi = 1;
      }
      else if ( i < N )
      {
        jhi = 2;
      }
      else
      {
        jhi = 2;
      }

      for ( j = 1; j <= jhi; j++ )
      {
        if ( i == 0 )
        {
          tval = t[0] - 1.0;
        }
        else if ( i < N )
        {
          tval = ( ( double ) ( jhi - j + 1 ) * t[i-1]
                 + ( double ) (       j - 1 ) * t[i] )
                 / ( double ) ( jhi         );
        }
        else
        {
          if ( j == 1 )
          {
            tval = t[N-1];
          }
          else
          {
            tval = t[N-1] + 1.0;
          }
        }

        yval = spline_cubic_val ( N, t, y, ypp, tval, &ypval, &yppval );

        printf ( "%12g  %12g  %12g\n", tval, yval, fcube ( tval ) );
      }
    }
    free ( ypp );
  }

  return;
# undef N
}
/******************************************************************************/

void test18 ( )

/******************************************************************************/
//
//  Purpose:
//
//    TEST18 tests SPLINE_CUBIC_SET and SPLINE_CUBIC_VAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 11

  int i;
  int ibcbeg;
  int ibcend;
  int j;
  int jhi;
  int k;
  double t[N];
  double tval;
  double y[N];
  double ybcbeg;
  double ybcend;
  double *ypp;
  double yppval;
  double ypval;
  double yval;
//
//  Set up the data.
//
  printf ( "\n" );
  printf ( "TEST18\n" );
  printf ( "  SPLINE_CUBIC_SET sets up a cubic spline;\n" );
  printf ( "  SPLINE_CUBIC_VAL evaluates it.\n" );
  printf ( "\n" );
  printf ( "  Cubic function, evenly spaced knots.\n" );
  printf ( "\n" );
  printf ( "        T            Y\n" );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    t[i] =  ( double ) ( i ) / ( double ) ( N - 1 );
    y[i] =  fcube ( t[i] );
    printf ( "  %12g  %12g\n", t[i], y[i] );
  }
//
//  Try all three types of boundary condition.
//
  for ( k = 0; k < 3; k++ )
  {
    if ( k == 0 )
    {
      ibcbeg = 0;
      ybcbeg = 0.0;

      ibcend = 0;
      ybcend = 0.0;

      printf ( "\n" );
      printf ( "  Boundary condition 0 at both ends:\n" );
      printf ( "  Spline is quadratic in boundary intervals.\n" );

    }
    else if ( k == 1 )
    {
      ibcbeg = 1;
      ybcbeg = fpcube ( t[0] );

      ibcend = 1;
      ybcend = fpcube ( t[N-1] );

      printf ( "\n" );
      printf ( "  Boundary condition 1 at both ends:\n" );
      printf ( "  Y'(left) =  %g\n", ybcbeg );
      printf ( "  Y'(right) = %g\n", ybcend );

    }
    else if ( k == 2 )
    {
      ibcbeg = 2;
      ybcbeg = fppcube ( t[0] );

      ibcend = 2;
      ybcend = fppcube ( t[N-1] );

      printf ( "\n" );
      printf ( "  Boundary condition 2 at both ends:\n" );
      printf ( "  YP''(left) =  %g\n", ybcbeg );
      printf ( "  YP''(right) = %g\n", ybcend );

    }

    ypp = spline_cubic_set ( N, t, y, ibcbeg, ybcbeg, ibcend, ybcend );

    printf ( "\n" );
    printf ( "     SPLINE''(T)       F''(T):\n" );
    printf ( "\n" );
    for ( i = 0; i < N; i++ )
    {
      printf ( "%12g  %12g\n", ypp[i], fppcube(t[i]) );
    }

    printf ( "\n" );
    printf ( "        T   SPLINE(T)     F(T)\n" );
    printf ( "\n" );

    for ( i = 0; i <= N; i++ )
    {
      if ( i == 0 )
      {
        jhi = 1;
      }
      else if ( i < N )
      {
        jhi = 2;
      }
      else
      {
        jhi = 2;
      }

      for ( j = 1; j <= jhi; j++ )
      {
        if ( i == 0 )
        {
          tval = t[0] - 1.0;
        }
        else if ( i < N )
        {
          tval = ( ( double ) ( jhi - j + 1 ) * t[i-1]
                 + ( double ) (       j - 1 ) * t[i] )
                 / ( double ) ( jhi         );
        }
        else
        {
          if ( j == 1 )
          {
            tval = t[N-1];
          }
          else
          {
            tval = t[N-1] + 1.0;
          }
        }

        yval = spline_cubic_val ( N, t, y, ypp, tval, &ypval, &yppval );

        printf ( "\n" );
        printf ( "%12g  %12g  %12g\n", tval, yval, fcube ( tval ) );
        printf ( "            %12g  %12g\n", ypval, fpcube ( tval ) );
        printf ( "            %12g  %12g\n", yppval, fppcube ( tval ) );
      }
    }
    free ( ypp );
  }

  return;
# undef N
}
/******************************************************************************/

void test19 ( )

/******************************************************************************/
//
//  Purpose:
//
//    TEST19 tests SPLINE_CUBIC_SET and SPLINE_CUBIC_VAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 2

  int i;
  int ibcbeg;
  int ibcend;
  int j;
  int jhi;
  int k1;
  int k2;
  double t[N];
  double tval;
  double y[N];
  double ybcbeg;
  double ybcend;
  double *ypp;
  double yppval;
  double ypval;
  double yval;
//
//  Set up the data.
//
  printf ( "\n" );
  printf ( "TEST19\n" );
  printf ( "  SPLINE_CUBIC_SET sets up a cubic spline;\n" );
  printf ( "  SPLINE_CUBIC_VAL evaluates it.\n" );
  printf ( "\n" );
  printf ( "  Cubic function, evenly spaced knots.\n" );
  printf ( "  ONLY TWO KNOTS!\n" );
  printf ( "\n" );
  printf ( "  The data to be interpolated:\n" );
  printf ( "\n" );
  printf ( "  Number of data values = %d\n", N );
  printf ( "\n" );
  printf ( "           T             Y\n" );
  printf ( "\n" );

  for ( i = 0; i < N; i++ )
  {
    t[i] =  ( double ) ( i ) / ( double ) ( N - 1 );
    y[i] =  fcube ( t[i] );
    printf ( "  %12g  %12g\n", t[i], y[i] );
  }
//
//  Try all nine pairs of boundary condition.
//
  for ( k1 = 0; k1 < 3; k1++ )
  {
    if ( k1 == 0 )
    {
      ibcbeg = 0;
      ybcbeg = 0.0;

      printf ( "\n" );
      printf ( "  Boundary condition 0 at left end.\n" );
    }
    else if ( k1 == 1 )
    {
      ibcbeg = 1;
      ybcbeg = fpcube ( t[0] );

      printf ( "\n" );
      printf ( "  Boundary condition 1 at left end.\n" );
      printf ( "  Y'(left) =  %g\n", ybcbeg );
    }
    else if ( k1 == 2 )
    {
      ibcbeg = 2;
      ybcbeg = fppcube ( t[0] );

      printf ( "\n" );
      printf ( "  Boundary condition 2 at left end.\n" );
      printf ( "  YP''(left) =  %g\n", ybcbeg );
    }

    for ( k2 = 0; k2 < 3; k2++ )
    {
      if ( k2 == 0 )
      {
        ibcend = 0;
        ybcend = 0.0;

        printf ( "  Boundary condition 0 at right end.\n" );
      }
      else if ( k2 == 1 )
      {
        ibcend = 1;
        ybcend = fpcube ( t[N-1] );

        printf ( "\n" );
        printf ( "  Boundary condition 1 at right end.\n" );
        printf ( "  Y'(right) = %g\n", ybcend );
      }
      else if ( k2 == 2 )
      {
        ibcend = 2;
        ybcend = fppcube ( t[N-1] );

        printf ( "\n" );
        printf ( "  Boundary condition 2 at right end.\n" );
        printf ( "  YP''(right) = %g\n", ybcend );

      }

      ypp = spline_cubic_set ( N, t, y, ibcbeg, ybcbeg, ibcend, ybcend );

      printf ( "\n" );
      printf ( "  SPLINE''(T)        F''(T)\n" );
      printf ( "\n" );
      for ( i = 0; i < N; i++ )
      {
        printf ( "12g  %12g\n", ypp[i], fppcube(t[i]) );
      }

      printf ( "\n" );
      printf ( "           T    SPLINE(T)         F(T)\n" );
      printf ( "\n" );

      for ( i = 0; i <= N; i++ )
      {

        if ( i == 0 )
        {
          jhi = 1;
        }
        else if ( i < N )
        {
          jhi = 2;
        }
        else
        {
          jhi = 2;
        }

        for ( j = 1; j <= jhi; j++ )
        {

          if ( i == 0 )
          {
            tval = t[0] - 1.0;
          }
          else if ( i < N )
          {
            tval = ( ( double ) ( jhi - j + 1 ) * t[i-1]
                   + ( double ) (       j - 1 ) * t[i] )
                   / ( double ) ( jhi         );
          }
          else
          {
            if ( j == 1 )
            {
              tval = t[N-1];
            }
            else
            {
              tval = t[N-1] + 1.0;
            }
          }

          yval = spline_cubic_val ( N, t, y, ypp, tval, &ypval, &yppval );

          printf ( "\n" );
          printf ( "%12g  %12g  %12g\n", tval, yval, fcube ( tval ) );
          printf ( "            %12g  %12g\n", ypval, fpcube ( tval ) );
          printf ( "            %12g  %12g\n", yppval, fppcube ( tval ) );
        }
      }
    }
    free ( ypp );
  }

  return;
# undef N
}
/******************************************************************************/

void test20 ( )

/******************************************************************************/
//
//  Purpose:
//
//    TEST20 tests SPLINE_HERMITE_SET and SPLINE_HERMITE_VAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define NDATA 4

  double *c;
  double fpval;
  double fval;
  int i;
  int j;
  int jhi;
  char mark;
  double pi = 3.141592653589793;
  double tdata[NDATA];
  double tval;
  double ydata[NDATA];
  double ypdata[NDATA];
  double ypval;
  double yval;

  printf ( "\n" );
  printf ( "TEST20\n" );
  printf ( "  SPLINE_HERMITE_SET sets up a Hermite spline;\n" );
  printf ( "  SPLINE_HERMITE_VAL evaluates it.\n" );
//
//  Set the data.
//
  for ( i = 0; i < NDATA; i++ )
  {
    tdata[i] = 0.5 * ( double ) i * pi / ( double ) ( NDATA - 1 );
    ydata[i] = sin ( tdata[i] );
    ypdata[i] = cos ( tdata[i] );
  }

  printf ( "\n" );
  printf ( "  Data\n" );
  printf ( "\n" );
  printf ( "     TDATA(I)     YDATA[I]     Y'DATA[I]\n" );
  printf ( "\n" );

  for ( i = 0; i < NDATA; i++ )
  {
    printf ( "  %12g  %12g  %12g\n", tdata[i], ydata[i], ypdata[i] );
  }
//
//  Set up the spline.
//
  c = spline_hermite_set ( NDATA, tdata, ydata, ypdata );
//
//  Now evaluate the spline all over the place.
//
  printf ( "\n" );
  printf ( "              T     Y(hermite)     Y(exact)   Y'(hermite)     Y'(exact)\n" );
  printf ( "\n" );

  for ( i = 0; i < NDATA; i++ )
  {
    if ( i == NDATA-1 )
    {
      jhi = 0;
    }
    else
    {
      jhi = 2;
    }

    for ( j = 0; j <= jhi; j++ )
    {
      tval = 0.5 * ( double ) ( 3 * i + j ) * pi
        / ( double ) ( 3 * ( NDATA - 1 ) );

      fval = sin ( tval );
      fpval = cos ( tval );

      spline_hermite_val ( NDATA, tdata, c, tval, &yval, &ypval );

      if ( j == 0 )
      {
        mark = '*';
      }
      else
      {
        mark = ' ';
      }
      printf ( "  %c  %12g  %12g  %12g  %12g  %12g\n", mark, tval, yval, fval, ypval, fpval );
    }

  }

  free ( c );

  return;
# undef NDATA
}
/******************************************************************************/

void test205 ( )

/******************************************************************************/
//
//  Purpose:
//
//    TEST205 tests SPLINE_LINEAR_INT and SPLINE_LINEAR_INTSET.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 4

  double a;
  double b;
  double data_x[N];
  double data_y[N];
  int i;
  double int_x[N+1] = { 0.0, 1.0, 4.0, 5.0, 10.0 };
  double int_v[N] = { 10.0, 2.0, 8.0, 27.5 };
  double value;

  printf ( "\n" );
  printf ( "TEST205\n" );
  printf ( "  SPLINE_LINEAR_INTSET is given some interval endpoints,\n" );
  printf ( "  and a value associated with each interval.\n" );
  printf ( "\n" );
  printf ( "  It determines a linear spline, with breakpoints\n" );
  printf ( "  at the centers of each interval, whose integral\n" );
  printf ( "  over each interval is equal to the given value.\n" );

  r8vec_print ( N+1, int_x, "  The interval end points:" );
  r8vec_print ( N, int_v, "  The desired interval integral values:" );

  spline_linear_intset ( N, int_x, int_v, data_x, data_y );

  r8vec_print ( N, data_x, "  The spline break points:" );
  r8vec_print ( N, data_y, "  The spline data values: " );

  printf ( "\n" );
  printf ( "  As a check, call SPLINE_LINEAR_INT to compute\n" );
  printf ( "  the integral of the spline over each interval,\n" );
  printf ( "  and compare to the desired value.\n" );
  printf ( "\n" );
  printf ( "       A         B       Desired      Computed\n" );
  printf ( "\n" );

  for ( i = 1; i <= N; i++ )
  {
    a = int_x[i-1];
    b = int_x[i];
    value = spline_linear_int ( N, data_x, data_y, a, b );
    printf ( "%8g  %8g  %12g  %12g\n", a, b, int_v[i-1], value );
  }

  return;
# undef N
}
/******************************************************************************/

void test21 ( )

/******************************************************************************/
//
//  Purpose:
//
//    TEST21 tests SPLINE_LINEAR_VAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 January 2007
//
//  Author:
//
//    John Burkardt
//
{
#define N 11

  double fval;
  int i;
  int j;
  int jhi;
  double t[N];
  double tval;
  double y[N];
  double ypval;
  double yval;

  printf ( "\n" );
  printf ( "TEST21\n" );
  printf ( "  SPLINE_LINEAR_VAL evaluates a linear spline.\n" );
  printf ( "\n" );
  printf ( "  Runge's function, evenly spaced knots.\n" );

  for ( i = 0; i < N; i++ )
  {
    t[i] =  ( ( double ) ( N - i     ) * (-1.0)
            + ( double ) (     i - 1 ) * (+1.0) )
            / ( double ) ( N     - 1 );
    y[i] =  frunge ( t[i] );
  }

  printf ( "\n" );
  printf ( "  The data to be interpolated:\n" );
  printf ( "\n" );
  printf ( "  Number of data values = %d\n", N );
  printf ( "\n" );
  printf ( "        T             Y\n" );
  printf ( "\n" );

  for ( i = 0; i < N; i++ )
  {
    printf ( "  %12g  %12g\n", t[i], y[i] );
  }

  printf ( "\n" );
  printf ( "  Interpolation:\n" );
  printf ( "\n" );
  printf ( "       T             Y            Yexact\n" );
  printf ( "\n" );
  for ( i = 0; i <= N; i++ )
  {

    if ( i == 0 )
    {
      jhi = 1;
    }
    else if ( i < N )
    {
      jhi = 2;
    }
    else
    {
      jhi = 2;
    }

    for ( j = 1; j <= jhi; j++ )
    {
      if ( i == 0 )
      {
        tval = t[0] - 1.0;
      }
      else if ( i < N )
      {
        tval = ( ( double ) ( jhi - j + 1 ) * t[i-1]
               + ( double ) (       j - 1 ) * t[i] )
               / ( double ) ( jhi         );
      }
      else
      {
        if ( j == 1 )
        {
          tval = t[N-1];
        }
        else
        {
          tval = t[N-1] + 1.0;
        }
      }

      spline_linear_val ( N, t, y, tval, &yval, &ypval );

      fval = frunge ( tval );

      printf ( "%12g  %12g  %12g\n", tval, yval, fval );

    }

  }

  return;
#undef N
}
/******************************************************************************/

void test215 ( )

/******************************************************************************/
//
//  Purpose:
//
//    TEST215 tests SPLINE_LINEAR_INT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 3

  double a;
  double b;
  int i;
  double int_val;
  double t[N] = { 2.0, 4.5,  7.5 };
  double y[N] = { 3.0, 3.75, 5.5 };

  printf ( "\n" );
  printf ( "TEST215\n" );
  printf ( "  SPLINE_LINEAR_INT computes the integral of a linear spline.\n" );
  printf ( "\n" );

  printf ( "\n" );
  printf ( "  The data to be interpolated:\n" );
  printf ( "\n" );
  printf ( "  Number of data values = %d\n", N );
  printf ( "\n" );
  printf ( "	  T		Y\n" );
  printf ( "\n" );

  for ( i = 0; i < N; i++ )
  {
    printf ( "  %12g  %12g\n", t[i], y[i] );
  }

  printf ( "\n" );
  printf ( "    A             B           Integral\n" );
  printf ( "\n" );

  for ( i = 1; i <= 5; i++ )
  {
    if ( i == 1 )
    {
      a = 0.0;
      b = 4.0;
    }
    else if ( i == 2 )
    {
      a = 4.0;
      b = 5.0;
    }
    else if ( i == 3 )
    {
      a = 5.0;
      b = 10.0;
    }
    else if ( i == 4 )
    {
      a = 0.0;
      b = 10.0;
    }
    else
    {
      a = 10.0;
      b = 0.0;
    }

    int_val = spline_linear_int ( N, t, y, a, b );

    printf ( "  %12g  %12g  %12g\n", a, b, int_val );
  }

  return;
# undef N
}
/******************************************************************************/

void test22 ( )

/******************************************************************************/
//
//  Purpose:
//
//    TEST22 tests SPLINE_OVERHAUSER_UNI_VAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define NDATA 11

  int i;
  int j;
  int jhi;
  char mark;
  int nsample = 4;
  double pi = 3.141592653589793;
  double tdata[NDATA];
  double thi;
  double tlo;
  double tval;
  double ydata[NDATA];
  double yval;

  printf ( "\n" );
  printf ( "TEST22\n" );
  printf ( "  SPLINE_OVERHAUSER_UNI_VAL evaluates the\n" );
  printf ( "    uniform Overhauser spline.\n" );

  for ( i = 0; i < NDATA; i++ )
  {
    tdata[i] = ( double ) ( i );
    ydata[i] = sin ( 2.0 * pi * tdata[i] / ( double ) ( NDATA - 1) );
  }

  printf ( "\n" );
  printf ( "  The data to be interpolated:\n" );
  printf ( "\n" );
  printf ( "  Number of data values = %d\n", NDATA );
  printf ( "\n" );
  printf ( "       T             Y\n" );
  printf ( "\n" );
  for ( i = 0; i < NDATA; i++ )
  {
    printf ( "  %12g  %12g\n", tdata[i], ydata[i] );
  }

  printf ( "\n" );
  printf ( "    T, Spline(T)\n" );
  printf ( "\n" );

  for ( i = 0; i <= NDATA; i++ )
  {
    if ( i == 0 )
    {
      tlo = tdata[0] - 0.5 * ( tdata[1] - tdata[0] );
      thi = tdata[0];
    }
    else if ( i < NDATA )
    {
      tlo = tdata[i-1];
      thi = tdata[i];
    }
    else if ( NDATA <= i )
    {
      tlo = tdata[NDATA-1];
      thi = tdata[NDATA-1] + 0.5 * ( tdata[NDATA-1] - tdata[NDATA-2] );
    }

    if ( i < NDATA )
    {
      jhi = nsample - 1;
    }
    else
    {
      jhi = nsample;
    }

    for ( j = 0; j <= jhi; j++ )
    {
      tval = ( ( double ) ( nsample - j ) * tlo
             + ( double ) (           j ) * thi )
             / ( double ) ( nsample     );

      yval = spline_overhauser_uni_val ( NDATA, tdata, ydata, tval );

      if ( 0 < i & j == 0 )
      {
        mark = '*';
      }
      else
      {
        mark = ' ';
      }
      printf ( "  %c  %12g  %12g\n", mark, tval, yval );

    }

  }

  return;
# undef NDATA
}
/******************************************************************************/

void test225 ( )

/******************************************************************************/
//
//  Purpose:
//
//    TEST225 tests SPLINE_OVERHAUSER_NONUNI_VAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 February 2012
//
//  Author:
//
//    John Burkardt
//
{
# define NDATA 11

  int i;
  int j;
  int jhi;
  char mark;
  int nsample = 4;
  double pi = 3.141592653589793;
  double tdata[NDATA];
  double thi;
  double tlo;
  double tval;
  double ydata[NDATA];
  double yval;

  printf ( "\n" );
  printf ( "TEST225\n" );
  printf ( "  SPLINE_OVERHAUSER_NONUNI_VAL evaluates the\n" );
  printf ( "    nonuniform Overhauser spline.\n" );
  printf ( "\n" );
  printf ( "  In this draft of a test, we simply repeat the work\n" );
  printf ( "  for the uniform test.\n" );

  for ( i = 0; i < NDATA; i++ )
  {
    tdata[i] = ( double ) ( i );
    ydata[i] = sin ( 2.0 * pi * tdata[i] / ( double ) ( NDATA - 1) );
  }

  printf ( "\n" );
  printf ( "  The data to be interpolated:\n" );
  printf ( "\n" );
  printf ( "  Number of data values = %d\n", NDATA );
  printf ( "\n" );
  printf ( "       T             Y\n" );
  printf ( "\n" );
  for ( i = 0; i < NDATA; i++ )
  {
    printf ( "  %12g  %12g\n", tdata[i], ydata[i] );
  }

  printf ( "\n" );
  printf ( "    T, Spline(T)\n" );
  printf ( "\n" );

  for ( i = 0; i <= NDATA; i++ )
  {
    if ( i == 0 )
    {
      tlo = tdata[0] - 0.5 * ( tdata[1] - tdata[0] );
      thi = tdata[0];
    }
    else if ( i < NDATA )
    {
      tlo = tdata[i-1];
      thi = tdata[i];
    }
    else if ( NDATA <= i )
    {
      tlo = tdata[NDATA-1];
      thi = tdata[NDATA-1] + 0.5 * ( tdata[NDATA-1] - tdata[NDATA-2] );
    }

    if ( i < NDATA )
    {
      jhi = nsample - 1;
    }
    else
    {
      jhi = nsample;
    }

    for ( j = 0; j <= jhi; j++ )
    {
      tval = ( ( double ) ( nsample - j ) * tlo
             + ( double ) (           j ) * thi )
             / ( double ) ( nsample     );

      yval = spline_overhauser_nonuni_val ( NDATA, tdata, ydata, tval );

      if ( 0 < i & j == 0 )
      {
        mark = '*';
      }
      else
      {
        mark = ' ';
      }
      printf ( "  %c  %12g  %12g\n", mark, tval, yval );

    }

  }

  return;
# undef NDATA
}
/******************************************************************************/

void test23 ( )

/******************************************************************************/
//
//  Purpose:
//
//    TEST23 tests SPLINE_OVERHAUSER_VAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 February 2012
//
//  Author:
//
//    John Burkardt
//
{
# define NDATA 4
# define NDIM 1

  int i;
  double tdata[NDATA];
  double tval;
  double ydata[NDIM*NDATA];
  double yval[NDIM];
  double zdata[NDATA];
  double zval[NDIM];

  printf ( "\n" );
  printf ( "TEST23\n" );
  printf ( "  SPLINE_OVERHAUSER_VAL evaluates the\n" );
  printf ( "    Overhauser spline.\n" );
//
//  Set the data.
//
  tdata[0] = 1.0;
  ydata[0+0*NDIM] =   0.0;
  zdata[0+0*NDIM] =   0.0;

  tdata[1] = 2.0;
  ydata[0+1*NDIM] =   1.0;
  zdata[0+1*NDIM] =   1.0;

  tdata[2] = 3.0;
  ydata[0+2*NDIM] =   2.0;
  zdata[0+2*NDIM] = - 1.0;

  tdata[3] = 4.0;
  ydata[0+3*NDIM] =   3.0;
  zdata[0+3*NDIM] =   0.0;

  printf ( "\n" );
  printf ( "  Data\n" );
  printf ( "  TDATA[I], YDATA[I], ZDATA[I]\n" );
  printf ( "\n" );
  for ( i = 0; i < NDATA; i++ )
  {
    printf ( "%12g  %12g  %12g\n",  tdata[i], ydata[0+i*NDIM], zdata[0+i*NDIM] );
  }
//
//  Now evaluate the spline all over the place.
//
  printf ( "\n" );
  printf ( "  T, Spline value\n" );
  printf ( "\n" );

  for ( i = 0; i <= 6 * NDATA + 3; i++ )
  {
    tval = ( ( double ) i ) / 6.0;
    spline_overhauser_val ( NDIM, NDATA, tdata, ydata, tval, yval );
    spline_overhauser_val ( NDIM, NDATA, tdata, zdata, tval, zval );

    printf ( "%12g  %12g  %12g\n", tval, yval[0], zval[0] );
  }

  return;
# undef NDATA
# undef NDIM
}
/******************************************************************************/

void test235 ( )

/******************************************************************************/
//
//  Purpose:
//
//    TEST235 tests SPLINE_PCHIP_SET and SPLINE_PCHIP_VAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 February 2012
//
//  Author:
//
//    John Burkardt
//
{
# define N 21
# define NE 101

  double d[N];
  double diff;
  double f[N];
  double fd[NE];
  double fe[NE];
  int i;
  double x[N];
  double xe[NE];

  printf ( "\n" );
  printf ( "TEST235\n" );
  printf ( "  SPLINE_PCHIP_SET sets up a piecewise cubic\n" );
  printf ( "    Hermite interpolant.\n" );
  printf ( "  SPLINE_PCHIP_VAL evaluates the interpolant.\n" );
  printf ( "\n" );
//
//  Compute Runge's function at N points in [-1,1].
//
  for ( i = 0; i < N; i++ )
  {
    x[i] = -1.0 + ( double ) ( i ) / 10.0;
    f[i] = frunge ( x[i] );
  }
//
//  SPLINE_PCHIP_SET takes the data in X and F, and constructs a table in D
//  that defines the interpolant.
//
  spline_pchip_set ( N, x, f, d );
//
//  Evaluate the interpolant and derivative at NE points from -1 to 0.
//
  for ( i = 0; i < NE; i++ )
  {
    xe[i] = -1.0 + ( double ) ( i ) / ( double ) ( NE - 1 );
  }

  spline_pchip_val ( N, x, f, d, NE, xe, fe );
//
//  Print the table of X, F(exact) and F(interpolated)
//
  for ( i = 0; i < NE; i++ )
  {
    diff = fe[i] - frunge ( xe[i] );

    printf ( "  %8g  %10g  %10g  %10g\n", xe[i], frunge ( xe[i] ), fe[i], diff );
  }

  return;
# undef N
# undef NE
}
/******************************************************************************/

void test24 ( )

/******************************************************************************/
//
//  Purpose:
//
//    TEST24 tests SPLINE_QUADRATIC_VAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 February 2012
//
//  Author:
//
//    John Burkardt
//
{
# define N 11

  double fval;
  int i;
  int j;
  int jhi;
  double t[N];
  double tval;
  double y[N];
  double ypval;
  double yval;

  printf ( "\n" );
  printf ( "TEST24\n" );
  printf ( "  SPLINE_QUADRATIC_VAL evaluates a\n" );
  printf ( "    quadratic spline.\n" );
  printf ( "\n" );
  printf ( "  Runge''s function, evenly spaced knots.\n" );

  for ( i = 0; i < N; i++ )
  {
    t[i] =  ( ( double ) ( N - i     ) * (-1.0)
            + ( double ) (     i - 1 ) * (+1.0) )
            / ( double ) ( N     - 1 );
    y[i] =  frunge ( t[i] );
  }

  printf ( "\n" );
  printf ( "  The data to be interpolated:\n" );
  printf ( "\n" );
  printf ( "  Number of data values = %d\n", N );
  printf ( "\n" );
  printf ( "       T             Y\n" );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    printf ( "  %12g  %12g\n", t[i], y[i] );
  }

  printf ( "\n" );
  printf ( "  Interpolated values\n" );
  printf ( "\n" );
  printf ( "       T             Y           Y(exact)\n" );
  printf ( "\n" );

  for ( i = 0; i <= N; i++ )
  {
    if ( i == 0 )
    {
      jhi = 1;
    }
    else if ( i < N )
    {
      jhi = 2;
    }
    else
    {
      jhi = 2;
    }

    for ( j = 1; j <= jhi; j++ )
    {
      if ( i == 0 )
      {
        tval = t[0] - 1.0;
      }
      else if ( i < N )
      {
        tval = ( ( double ) ( jhi - j + 1 ) * t[i-1]
               + ( double ) (       j - 1 ) * t[i] )
               / ( double ) ( jhi         );
      }
      else
      {
        if ( j == 1 )
        {
          tval = t[N-1];
        }
        else
        {
          tval = t[N-1] + 1.0;
        }
      }

      spline_quadratic_val ( N, t, y, tval, &yval, &ypval );

      fval = frunge ( tval );

      printf ( "%12g  %12g  %12g\n", tval, yval, fval );
    }
  }

  return;
# undef N
}
/******************************************************************************/

void parabola_formula ( double x, double *y, double *yp, double *ypp )

/******************************************************************************/
//
//  Purpose:
//
//    PARABOLA_FORMULA evaluates a parabola for us.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 February 2012
//
//  Author:
//
//    John Burkardt
//
{
  *y = 2.0 * x * x + 3.0 * x + 1.0;
  *yp = 4.0 * x + 3.0;
  *ypp = 4.0;

  return;
}
/******************************************************************************/

double frunge ( double x )

/******************************************************************************/
//
//  Purpose:
//
//    FRUNGE sets the Runge data values.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 February 2012
//
//  Author:
//
//    John Burkardt
//
{
  double fx;

  fx = 1.0 / ( 1.0 + 25.0 * x * x );

  return fx;
}
/******************************************************************************/

double fprunge ( double x )

/******************************************************************************/
//
//  Purpose:
//
//    FPRUNGE sets the Runge derivative values at the endpoints.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 February 2012
//
//  Author:
//
//    John Burkardt
//
{
  double bot;
  double fx;

  bot = 1.0 + 25.0 * x * x;
  fx = -50.0 * x / ( bot * bot );

  return fx;
}
/******************************************************************************/

double fpprunge ( double x )

/******************************************************************************/
//
//  Purpose:
//
//    FPPRUNGE sets the Runge second derivative values at the endpoints.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 January 2007
//
//  Author:
//
//    John Burkardt
//
{
  double bot;
  double fx;

  bot = 1.0 + 25.0 * x * x;
  fx = ( -50.0 + 3750.0 * x * x ) / ( bot * bot * bot );

  return fx;
}
/******************************************************************************/

double fcube ( double x )

/******************************************************************************/
//
//  Purpose:
//
//    FCUBE evaluates a cubic function.
//
//  Discussion:
//
//    Y(X) = ( ( X + 2 ) * X + 3 ) * X + 4
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 February 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the point at which the function is evaluated.
//
//    Output, double FCUBE, the value of the function.
//
{
  double fx;

  fx = ( ( x + 2.0 ) * x + 3.0 ) * x + 4.0;

  return fx;
}
/******************************************************************************/

double fpcube ( double x )

/******************************************************************************/
//
//  Purpose:
//
//    FPCUBE evaluates the derivative of a cubic function.
//
//  Discussion:
//
//    Y(X) = ( ( X + 2 ) * X + 3 ) * X + 4
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 February 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the point at which the function is evaluated.
//
//    Output, double FPCUBE, the value of the derivative of the function.
//
{
  double fx;

  fx = ( 3.0 * x + 4.0 ) * x + 3.0;

  return fx;
}
/******************************************************************************/

double fppcube ( double x )

/******************************************************************************/
//
//  Purpose:
//
//    FPPCUBE evaluates the second derivative of a cubic function.
//
//  Discussion:
//
//    Y(X) = ( ( X + 2 ) * X + 3 ) * X + 4
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 February 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the point at which the function is evaluated.
//
//    Output, double FPPCUBE, the value of the second derivative of the function.
//
{
  double fx;

  fx = 6.0 * x + 4.0;

  return fx;
}
