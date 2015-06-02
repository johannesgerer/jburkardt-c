# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <string.h>

# include "interp.h"

int main ( );
void test01 ( );
void test02 ( );
void test03 ( int data_num );
void test04 ( int data_num );
double *f_runge ( int m, int n, double x[] );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for INTERP_PRB.

  Discussion:

    INTERP_PRB tests the INTERP library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 March 2014

  Author:

    John Burkardt
*/
{
  int data_num;

  timestamp ( );
  printf ( "\n" );
  printf ( "INTERP_PRB\n" );
  printf ( "  C version:\n" );
  printf ( "  Test the INTERP library.\n" );

  test01 ( );

  test02 ( );

  data_num = 6;
  test03 ( data_num );

  data_num = 11;
  test03 ( data_num );

  data_num = 6;
  test04 ( data_num );

  data_num = 11;
  test04 ( data_num );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "INTERP_PRB\n" );
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

    TEST01 tests INTERP_NEAREST on 1-dimensional data.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 March 2014

  Author:

    John Burkardt
*/
{
  int after;
  int before;
  int data_num = 11;
  int fat;
  int i;
  int interp;
  int interp_num;
  int j;
  int m = 1;
  double p;
  double *p_data;
  double *p_interp;
  double *p_value;
  double t;
  double *t_data;
  double *t_interp;
  double t_max;
  double t_min;

  printf ( " \n" );
  printf ( "TEST01\n" );
  printf ( "  INTERP_NEAREST evaluates a nearest-neighbor interpolant.\n" );
  printf ( " \n" );
  printf ( "  In this example, the function we are interpolating is\n" );
  printf ( "  Runge''s function, with Chebyshev knots.\n" );

  t_min = -1.0;
  t_max = +1.0;

  t_data = cc_abscissas_ab ( t_min, t_max, data_num );

  p_data = f_runge ( m, data_num, t_data );

  printf ( "\n" );
  printf ( "  The data to be interpolated:\n" );
  printf ( "\n" );
  printf ( "  Spatial dimension =     %d\n", m );
  printf ( "  Number of data values = %d\n", data_num );
  printf ( "\n" );
  printf ( "       T_data        P_data\n" );
  printf ( "\n" );
  for ( j = 0; j < data_num; j++ )
  {
    printf ( "  %14.6g  %14.6g\n", t_data[j], p_data[0+j*m] );
  }
/*
  Our interpolation values will include the original T values, plus
  3 new values in between each pair of original values.
*/
  before = 4;
  fat = 3;
  after = 2;

  interp_num = before + 1 + ( data_num - 1 ) * ( fat + 1 ) + after;

  t_interp = r8vec_expand_linear2 ( data_num, t_data, before, fat, after );

  p_interp = interp_nearest ( m, data_num, t_data, p_data, interp_num, t_interp );

  p_value = f_runge ( m, interp_num, t_interp );

  printf ( "\n" );
  printf ( "  Interpolation:\n" );
  printf ( "\n" );
  printf ( "    T_interp      P_interp        P_exact        Error\n" );
  printf ( "\n" );

  for ( j = 0; j < interp_num; j++ )
  {
    printf ( "  %10.4f  %14.6g  %14.6g  %10.2g\n",
      t_interp[j], p_interp[0+j*m], p_value[0+j*m], 
      p_interp[0+j*m] - p_value[0+j*m] ); 
  }

  free ( p_data );
  free ( p_interp );
  free ( p_value );
  free ( t_data );
  free ( t_interp );

  return;
}
/******************************************************************************/

void test02 ( )

/******************************************************************************/
/*
  Purpose:

    TEST02 tests INTERP_LINEAR on 1-dimensional data.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 March 2014

  Author:

    John Burkardt
*/
{
  int after;
  int before;
  int data_num = 11;
  int fat;
  int i;
  int interp;
  int interp_num;
  int j;
  int m = 1;
  double p;
  double *p_data;
  double *p_interp;
  double *p_value;
  double t;
  double *t_data;
  double *t_interp;
  double t_max;
  double t_min;

  printf ( " \n" );
  printf ( "TEST02\n" );
  printf ( "  INTERP_LINEAR evaluates a piecewise linear spline.\n" );
  printf ( " \n" );
  printf ( "  In this example, the function we are interpolating is\n" );
  printf ( "  Runge's function, with evenly spaced knots.\n" );

  t_min = -1.0;
  t_max = +1.0;

  t_data = ncc_abscissas_ab ( t_min, t_max, data_num );

  p_data = f_runge ( m, data_num, t_data );

  printf ( "\n" );
  printf ( "  The data to be interpolated:\n" );
  printf ( "\n" );
  printf ( "  Spatial dimension =     %d\n", m );
  printf ( "  Number of data values = %d\n", data_num );
  printf ( "\n" );
  printf ( "       T_data        P_data\n" );
  printf ( "\n" );
  for ( j = 0; j < data_num; j++ )
  {
    printf ( "  %14.6g  %14.6g\n", t_data[j], p_data[0+j*m] );
  }
/*
  Our interpolation values will include the original T values, plus
  3 new values in between each pair of original values.
*/
  before = 4;
  fat = 3;
  after = 2;

  interp_num = before + 1 + ( data_num - 1 ) * ( fat + 1 ) + after;

  t_interp = r8vec_expand_linear2 ( data_num, t_data, before, fat, after );

  p_interp = interp_linear ( m, data_num, t_data, p_data, interp_num, 
    t_interp );

  p_value = f_runge ( m, interp_num, t_interp );

  printf ( "\n" );
  printf ( "  Interpolation:\n" );
  printf ( "\n" );
  printf ( "    T_interp      P_interp        P_exact        Error\n" );
  printf ( "\n" );

  for ( j = 0; j < interp_num; j++ )
  {
    printf ( "  %10.4f  %14.6g  %14.6g  %10.2g\n",
      t_interp[j], p_interp[0+j*m], p_value[0+j*m], 
      p_interp[0+j*m] - p_value[0+j*m] ); 
  }

  free ( p_data );
  free ( p_interp );
  free ( p_value );
  free ( t_data );
  free ( t_interp );

  return;
}
/******************************************************************************/

void test03 ( int data_num )

/******************************************************************************/
/*
  Purpose:

    TEST03 tests INTERP_LAGRANGE on 1-dimensional data, equally spaced data.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 March 2014

  Author:

    John Burkardt

  Parameters:

    Input, int DATA_NUM, the number of data values.
*/
{
  int after;
  int before;
  int fat;
  int i;
  int interp;
  int interp_num;
  int j;
  int m = 1;
  double p;
  double *p_data;
  double *p_interp;
  double *p_value;
  double t;
  double *t_data;
  double *t_interp;
  double t_max;
  double t_min;

  printf ( " \n" );
  printf ( "TEST03\n" );
  printf ( "  INTERP_LAGRANGE evaluates a polynomial interpolant.\n" );
  printf ( "  In this example, the function we are interpolating is\n" );
  printf ( "  Runge's function, with evenly spaced knots.\n" );

  t_min = -1.0;
  t_max = +1.0;

  t_data = ncc_abscissas_ab ( t_min, t_max, data_num );

  p_data = f_runge ( m, data_num, t_data );

  printf ( "\n" );
  printf ( "  The data to be interpolated:\n" );
  printf ( "\n" );
  printf ( "  Spatial dimension =     %d\n", m );
  printf ( "  Number of data values = %d\n", data_num );
  printf ( "\n" );
  printf ( "       T_data        P_data\n" );
  printf ( "\n" );
  for ( j = 0; j < data_num; j++ )
  {
    printf ( "  %14.6g  %14.6g\n", t_data[j], p_data[0+j*m] );
  }
/*
  Our interpolation values will include the original T values, plus
  3 new values in between each pair of original values.
*/
  before = 4;
  fat = 3;
  after = 2;

  interp_num = before + 1 + ( data_num - 1 ) * ( fat + 1 ) + after;

  t_interp = r8vec_expand_linear2 ( data_num, t_data, before, fat, after );

  p_interp = interp_lagrange ( m, data_num, t_data, p_data, interp_num, 
    t_interp );

  p_value = f_runge ( m, interp_num, t_interp );

  printf ( "\n" );
  printf ( "  Interpolation:\n" );
  printf ( "\n" );
  printf ( "    T_interp      P_interp        P_exact        Error\n" );
  printf ( "\n" );

  for ( j = 0; j < interp_num; j++ )
  {
    printf ( "  %10.4f  %14.6g  %14.6g  %10.2g\n",
      t_interp[j], p_interp[0+j*m], p_value[0+j*m], 
      p_interp[0+j*m] - p_value[0+j*m] ); 
  }

  free ( p_data );
  free ( p_interp );
  free ( p_value );
  free ( t_data );
  free ( t_interp );

  return;
}
/******************************************************************************/

void test04 ( int data_num )

/******************************************************************************/
/*
  Purpose:

    TEST04 tests INTERP_LAGRANGE on 1-dimensional data, Clenshaw-Curtis data.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 March 2014

  Author:

    John Burkardt

  Parameters:

    Input, int DATA_NUM, the number of data values.
*/
{
  int after;
  int before;
  int fat;
  int i;
  int interp;
  int interp_num;
  int j;
  int m = 1;
  double p;
  double *p_data;
  double *p_interp;
  double *p_value;
  double t;
  double *t_data;
  double *t_interp;
  double t_max;
  double t_min;

  printf ( " \n" );
  printf ( "TEST04\n" );
  printf ( "  INTERP_LAGRANGE evaluates a polynomial interpolant.\n" );
  printf ( " \n" );
  printf ( "  In this example, the function we are interpolating is\n" );
  printf ( "  Runge's function, with Clenshaw Curtis knots.\n" );

  t_min = -1.0;
  t_max = +1.0;

  t_data = cc_abscissas_ab ( t_min, t_max, data_num );

  p_data = f_runge ( m, data_num, t_data );

  printf ( "\n" );
  printf ( "  The data to be interpolated:\n" );
  printf ( "\n" );
  printf ( "  Spatial dimension =     %d\n", m );
  printf ( "  Number of data values = %d\n", data_num );
  printf ( "\n" );
  printf ( "       T_data        P_data\n" );
  printf ( "\n" );
  for ( j = 0; j < data_num; j++ )
  {
    printf ( "  %14.6g  %14.6g\n", t_data[j], p_data[0+j*m] );
  }
/*
  Our interpolation values will include the original T values, plus
  3 new values in between each pair of original values.
*/
  before = 4;
  fat = 3;
  after = 2;

  interp_num = before + 1 + ( data_num - 1 ) * ( fat + 1 ) + after;

  t_interp = r8vec_expand_linear2 ( data_num, t_data, before, fat, after );

  p_interp = interp_lagrange ( m, data_num, t_data, p_data, interp_num, 
    t_interp );

  p_value = f_runge ( m, interp_num, t_interp );

  printf ( "\n" );
  printf ( "  Interpolation:\n" );
  printf ( "\n" );
  printf ( "    T_interp      P_interp        P_exact        Error\n" );
  printf ( "\n" );

  for ( j = 0; j < interp_num; j++ )
  {
    printf ( "  %10.4f  %14.6g  %14.6g  %10.2g\n",
      t_interp[j], p_interp[0+j*m], p_value[0+j*m], 
      p_interp[0+j*m] - p_value[0+j*m] ); 
  }

  free ( p_data );
  free ( p_interp );
  free ( p_value );
  free ( t_data );
  free ( t_interp );

  return;
}
/******************************************************************************/

double *f_runge ( int m, int n, double x[] )

/******************************************************************************/
/*
  Purpose:

    F_RUNGE evaluates the Runge function.

  Discussion:

    Interpolation of the Runge function at evenly spaced points in [-1,1]
    is a common test.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 March 2014

  Author:

    John Burkardt

  Parameters:

    Input, int M, the spatial dimension.

    Input, int N, the number of evaluation points.

    Input, double X[M*N], the evaluation points.

    Output, double F_RUNGE[N], the function values.
*/
{
  double *f;
  int i;
  int j;
  double t;

  f = ( double * ) malloc ( n * sizeof ( double ) );

  for ( j = 0; j < n; j++ )
  {
    t = 0.0;
    for ( i = 0; i < m; i++ )
    {
      t = t + pow ( x[i+j*m], 2 );
    }
    f[j] = 1.0 / ( 1.0 + 25.0 * t );
  }

  return f;
}
