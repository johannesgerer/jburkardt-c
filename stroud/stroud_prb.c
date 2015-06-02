# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <string.h>

# include "stroud.h"
/*
  GLOBAL VARIABLES USED TO SELECT INTEGRATION FUNCTION.
*/
int function_1d_index = 0;
int function_2d_index = 0;
int function_3d_index = 0;
int function_nd_index = 0;

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );
void test04 ( );
void test045 ( );
void test05 ( );
void test052 ( );
void test054 ( );
void test07 ( );
void test08 ( );
void test085 ( );
void test09 ( );
void test10 ( );
void test11 ( );
void test12 ( );
void test13 ( );
void test14 ( );
void test15 ( );
void test16 ( );
void test163 ( );
void cn_geg_test ( int n, double alpha, int expon[] );
void test165 ( );
void cn_jac_test ( int n, double alpha, double beta, int expon[] );
void test167 ( );
void cn_leg_test ( int n, int expon[] );
void test17 ( );
void test18 ( );
void test19 ( );
void test20 ( );
void test205 ( );
void test207 ( );
void en_r2_test ( int n, int expon[] );
void test2075 ( );
void epn_glg_test ( int n, int expon[], double alpha );
void test208 ( );
void epn_lag_test ( int n, int expon[] );
void test21 ( );
void test215 ( );
void test22 ( );
void test23 ( );
void test24 ( );
void test25 ( );
void test255 ( );
void test26 ( );
void test27 ( );
void test28 ( );
void test29 ( );
void test30 ( );
void test31 ( );
void test32 ( );
void test322 ( );
void test324 ( );
void test326 ( );
void test33 ( );
void test335 ( );
void test34 ( );
void test345 ( );
void test35 ( );
void test36 ( );
void test37 ( );
void test38 ( );
void test39 ( );
void test40 ( );
void test41 ( );
void test42 ( );
void test425 ( );
void test43 ( );
void test44 ( );
void test45 ( );
void test46 ( );
void test47 ( );
void test48 ( );
void test49 ( );
double fu18 ( double x );
double fl18 ( double x );
double fu28 ( double x, double y );
double fl28 ( double x, double y );

double function_1d ( double x );
void function_1d_name ( char *name );
int function_1d_num ( );

double function_2d ( double x, double y );
void function_2d_name ( char *name );
int function_2d_num ( );

double function_3d ( double x, double y, double z );
void function_3d_name ( char *name );
int function_3d_num ( );

double function_nd ( int n, double x[] );
void function_nd_name ( char *name );
int function_nd_num ( );

double f_1_2d ( double x, double y );
double f_x_2d ( double x, double y );
double f_r_2d ( double x, double y );
double mono_000_3d ( int n, double x[] );
double mono_111_3d ( int n, double x[] );
double mono_202_3d ( int n, double x[] );
double mono_422_3d ( int n, double x[] );
double *setsim ( int n );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for STROUD_PRB.

  Discussion:

    STROUD_PRB tests the STROUD library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 September 2010

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "STROUD_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the STROUD library.\n" );

  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
  test045 ( );
  test05 ( );
  test052 ( );
  test054 ( );
  test07 ( );
  test08 ( );
  test085 ( );
  test09 ( );

  test10 ( );
  test11 ( );
  test12 ( );
  test13 ( );
  test14 ( );
  test15 ( );
  test16 ( );
  test163 ( );
  test165 ( );
  test167 ( );
  test17 ( );
  test18 ( );
  test19 ( );

  test20 ( );
  test205 ( );
  test207 ( );
  test2075 ( );
  test208 ( );
  test21 ( );
  test215 ( );
  test22 ( );
  test23 ( );
  test24 ( );
  test25 ( );
  test255 ( );
  test26 ( );
  test27 ( );
  test28 ( );
  test29 ( );

  test30 ( );
  test31 ( );
  test32 ( );
  test322 ( );
  test324 ( );
  test326 ( );
  test33 ( );
  test335 ( );
  test34 ( );
  test345 ( );
  test35 ( );
  test36 ( );
  test37 ( );
  test38 ( );
  test39 ( );

  test40 ( );
  test41 ( );
  test42 ( );
  test425 ( );
  test43 ( );
  test44 ( );
  test45 ( );
  test46 ( );
  test47 ( );
  test48 ( );
  test49 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "STROUD_PRB\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/*****************************************************************************80*/

void test01 ( )

/*****************************************************************************80

  Purpose:

    TEST01 tests BALL_F1_ND, BALL_F3_ND.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 March 2008

  Author:

    John Burkardt
*/
{
  double *center;
  int i;
  int n;
  int n_max = 3;
  char name[8];
  int num;
  double result1;
  double result2;
  double r;

  r = 2.0;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  For integrals in a ball in ND:\n" );
  printf ( "  BALL_F1_ND approximates the integral;\n" );
  printf ( "  BALL_F3_ND approximates the integral.\n" );

  for ( n = 2; n <= n_max; n++ )
  {
    center = ( double * ) malloc ( n * sizeof ( double ) );
    for ( i = 0; i < n; i++ )
    {
      center[i] = ( double ) ( i + 1 );
    }
    printf ( "\n" );
    printf ( "  Spatial dimension N = %d\n", n );
    printf ( "  Ball center:\n" );
    for ( i = 0; i < n; i++ )
    {
      printf ( "  %10f", center[i] );
    }
    printf ( "\n" );
    printf ( "  Ball radius = %f\n", r );
    printf ( "  Ball volume = %f\n", ball_volume_nd ( n, r ) );
    printf ( "\n" );
    printf ( "    Rule:      F1          F3\n" );
    printf ( "    F(X)\n" );
    printf ( "\n" );

    num = function_nd_num ( );

    for ( i = 1; i <= num; i++ )
    {
      function_nd_index = i;
      function_nd_name ( name );

      result1 = ball_f1_nd ( function_nd, n, center, r );
      result2 = ball_f3_nd ( function_nd, n, center, r );
      printf ( "  %8s  %14f  %14f\n", name, result1, result2 );
    }
    free ( center );
  }
  return;
}
/*****************************************************************************80*/

void test02 ( )

/*****************************************************************************80

  Purpose:

    TEST02 tests BALL_MONOMIAL_ND.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 April 2008

  Author:

    John Burkardt
*/
{
  double center[3] = { 0.0, 0.0, 0.0 };
  int dim_num = 3;
  int p[3];
  double result1;
  double result2;
  double r = 2.0;
  char string[11];
  int test;
  int test_num = 4;

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  For the integral of a monomial in a ball in ND:\n" );
  printf ( "  BALL_MONOMIAL_ND approximates the integral.\n" );
  printf ( "  BALL_F1_ND, which can handle general integrands,\n" );
  printf ( "    will be used for comparison.\n" );
  printf ( "\n" );
  printf ( "  Spatial dimension N = %d\n", dim_num );
  printf ( "  Ball radius = %f\n", r );
  printf ( "  Ball volume = %f\n", ball_volume_nd ( dim_num, r ) );
  printf ( "\n" );
  printf ( "    Rule:     MONOMIAL    F1\n" );
  printf ( "    F(X)\n" );
  printf ( "\n" );

  for ( test = 1; test <= test_num; test++ )
  {
    if ( test == 1 )
    {
      strcpy ( string, "         1" );
      p[0] = 0;
      p[1] = 0;
      p[2] = 0;
      result2 = ball_f1_nd ( mono_000_3d, dim_num, center, r );
    }
    else if ( test == 2 )
    {
      strcpy ( string, "       xyz" );
      p[0] = 1;
      p[1] = 1;
      p[2] = 1;
      result2 = ball_f1_nd ( mono_111_3d, dim_num, center, r );
    }
    else if ( test == 3 )
    {
      strcpy ( string, "   x^2 z^2" );
      p[0] = 2;
      p[1] = 0;
      p[2] = 2;
      result2 = ball_f1_nd ( mono_202_3d, dim_num, center, r );
    }
    else if ( test == 4 )
    {
      strcpy ( string, " x^4y^2z^2" );
      p[0] = 4;
      p[1] = 2;
      p[2] = 2;
      result2 = ball_f1_nd ( mono_422_3d, dim_num, center, r );
    }

    result1 = ball_monomial_nd ( dim_num, p, r );

    printf ( "  %11s  %14f  %14f\n", string, result1, result2 );
  }

  return;
}
/*****************************************************************************80*/

void test03 ( )

/*****************************************************************************80

  Purpose:

    TEST03 tests BALL_UNIT_**_3D.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 April 2008

  Author:

    John Burkardt
*/
{
  int i;
  char name[8];
  int num;
  double result1;
  double result2;
  double result3;

  printf ( "\n" );
  printf ( "TEST03\n" );
  printf ( "  For integrals in the unit ball in 3D:\n" );
  printf ( "  BALL_UNIT_07_3D uses a formula of degree 7;\n" );
  printf ( "  BALL_UNIT_14_3D uses a formula of degree 14;\n" );
  printf ( "  BALL_UNIT_15_3D uses a formula of degree 15.\n" );
  printf ( "\n" );
  printf ( "  Unit ball volume = %f\n", ball_unit_volume_nd ( 3 ) );
  printf ( "\n" );
  printf ( "    Rule:      #7             #14           #15\n" );
  printf ( "    F(X)\n" );
  printf ( "\n" );

  num = function_3d_num ( );

  for ( i = 1; i <= num; i++ )
  {
    function_3d_index = i;
    function_3d_name ( name );

    result1 = ball_unit_07_3d ( function_3d );
    result2 = ball_unit_14_3d ( function_3d );
    result3 = ball_unit_15_3d ( function_3d );

    printf ( "  %8s  %14f  %14f  %14f\n", name, result1, result2, result3 );
  }
  return;
}
/*****************************************************************************80*/

void test04 ( )

/*****************************************************************************80

  Purpose:

    TEST04 tests BALL_UNIT_F1_ND, BALL_UNIT_F3_ND.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 April 2008

  Author:

    John Burkardt
*/
{
  int i;
  int n;
  int n_max = 3;
  char name[8];
  int num;
  double result1;
  double result2;

  printf ( "\n" );
  printf ( "TEST04\n" );
  printf ( "  For integrals inside the unit ball in ND:\n" );
  printf ( "  BALL_UNIT_F1_ND approximates the integral;\n" );
  printf ( "  BALL_UNIT_F3_ND approximates the integral.\n" );
  printf ( "\n" );

  for ( n = 2; n <= n_max; n++ )
  {
    printf ( "\n" );
    printf ( "  Spatial dimension N = %d\n", n );
    printf ( "  Unit ball volume = %f\n", ball_unit_volume_nd ( n ) );
    printf ( "\n" );
    printf ( "\n" );
    printf ( "    Rule:      F1          F3\n" );
    printf ( "    F(X)\n" );
    printf ( "\n" );

    num = function_nd_num ( );

    for ( i = 1; i <= num; i++ )
    {
      function_nd_index = i;
      function_nd_name ( name );

      result1 = ball_unit_f1_nd ( function_nd, n );
      result2 = ball_unit_f3_nd ( function_nd, n );

      printf ( "  %8s  %14f  %14f\n", name, result1, result2 );
    }
  }
  return;
}
/*****************************************************************************80*/

void test045 ( )

/*****************************************************************************80

  Purpose:

    TEST045 tests BALL_UNIT_VOLUME_3D.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    04 March 2008

  Author:

    John Burkardt
*/
{
  int dim_num = 3;

  printf ( "\n" );
  printf ( "TEST045\n" );
  printf ( "  In 3 dimensions:\n" );
  printf ( "  BALL_UNIT_VOLUME_3D gets the volume of the unit ball.\n" );
  printf ( "  BALL_UNIT_VOLUME_ND will be called for comparison.\n" );
  printf ( "\n" );
  printf ( "    N    Volume    Method\n" );
  printf ( "\n" );

  printf ( "  %3d  %14f  BALL_UNIT_VOLUME_3D\n",
    dim_num, ball_unit_volume_3d ( ) );

  printf ( "  %3d  %14f  BALL_UNIT_VOLUME_ND\n",
    dim_num, ball_unit_volume_nd ( dim_num ) );

  return;
}
/*****************************************************************************80*/

void test05 ( )

/*****************************************************************************80

  Purpose:

    TEST05 tests BALL_UNIT_VOLUME_ND.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    04 March 2008

  Author:

    John Burkardt
*/
{
  int dim_num;

  printf ( "\n" );
  printf ( "TEST05\n" );
  printf ( "  BALL_UNIT_VOLUME_ND computes the volume\n" );
  printf ( "    of the unit ball in ND.\n" );
  printf ( "\n" );
  printf ( "    N    Volume\n" );
  printf ( "\n" );

  for ( dim_num = 2; dim_num <= 10; dim_num++ )
  {
    printf ( "  %3d  %14f\n", dim_num, ball_unit_volume_nd ( dim_num ) );
  }
  return;
}
/*****************************************************************************80*/

void test052 ( )

/*****************************************************************************80

  Purpose:

    TEST052 tests BALL_VOLUME_3D.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    04 March 2008

  Author:

    John Burkardt
*/
{
  int i;
  int n = 3;
  double r;

  printf ( "\n" );
  printf ( "TEST052\n" );
  printf ( "  In 3 dimensions:\n" );
  printf ( "  BALL_VOLUME_3D computes the volume of a unit ball.\n" );
  printf ( "  BALL_VOLUME_ND will be called for comparison.\n" );
  printf ( "\n" );
  printf ( "    N    R      Volume    Method\n" );
  printf ( "\n" );

  r = 1.0;

  for ( i = 1; i <= 3; i++ )
  {
    printf ( "  %3d  %14f  %14f  BALL_VOLUME_3D\n", n, r, ball_volume_3d ( r ) );

    printf ( "  %3d  %14f  %14f  BALL_VOLUME_ND\n", n, r, ball_volume_nd ( n, r ) );

    r = r * 2.0;
  }
  return;
}
/*****************************************************************************80*/

void test054 ( )

/*****************************************************************************80

  Purpose:

    TEST054 tests BALL_VOLUME_ND.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    04 March 2008

  Author:

    John Burkardt
*/
{
  int i;
  int n;
  double r;

  printf ( "\n" );
  printf ( "TEST054\n" );
  printf ( "  BALL_UNIT_VOLUME_ND computes the volume of\n" );
  printf ( "    the unit ball in N dimensions.\n" );
  printf ( "\n" );
  printf ( "    N        R      Volume\n" );
  printf ( "\n" );

  for ( n = 2; n <= 10; n++ )
  {
    r = 0.5;
    for ( i = 1; i <= 3; i++ )
    {
      printf ( "  %3d  %14f  %14f\n", n, r,  ball_volume_nd ( n, r ) );
      r = r * 2.0;
    }
  }
  return;
}
/*****************************************************************************80*/

void test07 ( )

/*****************************************************************************80

  Purpose:

    TEST07 tests CIRCLE_ANNULUS.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 April 2008

  Author:

    John Burkardt
*/
{
  double area;
  double center[2];
  double center_test[2*2] = {
    0.0, 0.0,
    0.0, 0.0 };
  int dim;
  int dim_num = 2;
  int i;
  int j;
  char name[8];
  int num;
  int nr;
  double radius1;
  double radius1_test[2] = { 0.0, 1.0 };
  double radius2;
  double radius2_test[2] = { 1.0, 2.0 };
  double result;
  int test_num = 2;

  printf ( "\n" );
  printf ( "TEST07\n" );
  printf ( "  CIRCLE_ANNULUS estimates integrals\n" );
  printf ( "    in a circular annulus.\n" );
  printf ( "\n" );
  printf ( "        F       CENTER         Radius1   Radius2   NR  Result\n" );
  printf ( "\n" );

  for ( i = 0; i < test_num; i++ )
  {
    for ( dim = 0; dim < dim_num; dim++ )
    {
      center[dim] = center_test[dim+i*dim_num];
    }
    radius1 = radius1_test[i];
    radius2 = radius2_test[i];

    area = circle_annulus_area_2d ( radius1, radius2 );

    printf ( "\n" );
    printf ( "     Area  %8f  %8f  %8f  %8f  %10f\n",
      center[0], center[1], radius1, radius2, area );

    num = function_2d_num ( );

    for ( j = 1; j <= num; j++ )
    {
      function_2d_index = j;

      for ( nr = 1; nr <= 4; nr++ )
      {
        result = circle_annulus ( function_2d, center, radius1, radius2, nr );

        function_2d_name ( name );

        printf ( "  %8s  %10f  %10f  %10f  %10f  %2d  %10f\n",
          name, center[0], center[1], radius1, radius2, nr, result );
      }
    }
  }
  return;
}
/*****************************************************************************80*/

void test08 ( )

/*****************************************************************************80

  Purpose:

    TEST08 tests CIRCLE_ANNULUS, CIRCLE_RT_SET, CIRCLE_RT_SUM.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 April 2008

  Author:

    John Burkardt
*/
{
  double area;
  double center[2];
  double center_test[2*3] = {
    0.0, 0.0,
    0.0, 0.0,
    0.0, 0.0 };
  int dim;
  int dim_num = 2;
  int i;
  int j;
  char name[8];
  int nc;
  int num;
  int nr;
  int nr2;
  int nt;
  double ra[5];
  double radius1;
  double radius1_test[3] = { 0.0, 1.0, 1.0 };
  double radius2;
  double radius2_test[3] = { 1.0, 2.0, 3.0 };
  double result1;
  double result2;
  double result3;
  double rw[5];
  int rule;
  double ta[20];
  int test_num = 3;
  double tw[20];
  double zw;

  printf ( "\n" );
  printf ( "TEST08\n" );
  printf ( "  CIRCLE_ANNULUS estimates integrals in a\n" );
  printf ( "    circular annulus.\n" );
  printf ( "  CIRCLE_RT_SET sets up a rule for a circle;\n" );
  printf ( "  CIRCLE_RT_SUM applies the rule.\n" );
  printf ( "\n" );
  printf ( "  RESULT1 = CIRCLE_ANNULUS result.\n" );
  printf ( "  RESULT2 = Difference of two CIRCLE_RT_SUM results.\n" );
  printf ( "\n" );
  printf ( "        F      CENTER       Radius1   Radius2   Result1 Result2\n" );
  printf ( "\n" );

  for ( i = 0; i < test_num; i++ )
  {
    for ( dim = 0; dim < dim_num; dim++ )
    {
      center[dim] = center_test[dim+i*dim_num];
    }
    radius1 = radius1_test[i];
    radius2 = radius2_test[i];

    area = circle_annulus_area_2d ( radius1, radius2 );

    strcpy ( name, "   Area" );
    printf ( "\n" );
    printf ( "  %8s  %9f  %9f  %9f  %9f  %9f\n",
     name, center[0], center[1], radius1, radius2, area );

    rule = 9;
    circle_rt_size ( rule, &nr2, &nt, &nc );
    circle_rt_set ( rule, nr2, nt, nc, ra, rw, ta, tw, &zw );

    num = function_2d_num ( );

    for ( j = 1; j <= num; j++ )
    {
      function_2d_index = j;
      function_2d_name ( name );

      nr = 5;
      result1 = circle_annulus ( function_2d, center, radius1, radius2, nr );

      result2 = circle_rt_sum ( function_2d, center, radius1, nr2, ra, rw, nt,
        ta, tw, zw );

      result3 = circle_rt_sum ( function_2d, center, radius2, nr2, ra, rw, nt,
        ta, tw, zw );

     printf ( "  %8s  %9f  %9f  %9f  %9f  %9f  %9f\n",
       name, center[0], center[1], radius1, radius2, result1, result3 - result2 );
    }
  }
  return;
}
/*****************************************************************************80*/

void test085 ( )

/*****************************************************************************80

  Purpose:

    TEST085 tests CIRCLE_ANNULUS_AREA_2D.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 April 2008

  Author:

    John Burkardt
*/
{
  double area;
  double center[2];
  double center_test[2*3] = {
    0.0, 0.0,
    1.0, 0.0,
    3.0, 4.0 };
  int dim;
  int dim_num = 2;
  int i;
  int ntest = 3;
  double radius1;
  double radius1_test[3] = { 0.0, 1.0, 1.0 };
  double radius2;
  double radius2_test[3] = { 1.0, 2.0, 3.0 };

  printf ( "\n" );
  printf ( "TEST085\n" );
  printf ( "  CIRCLE_ANNULUS_AREA_2D computes the area of a\n" );
  printf ( "    circular annulus.\n" );
  printf ( "\n" );
  printf ( "      CENTER       Radius1   Radius2   Area\n" );
  printf ( "\n" );

  for ( i = 0; i < ntest; i++ )
  {
    for ( dim = 0; dim < dim_num; dim++ )
    {
      center[dim] = center_test[dim+i*dim_num];
    }

    radius1 = radius1_test[i];
    radius2 = radius2_test[i];

    area = circle_annulus_area_2d ( radius1, radius2 );

    printf ( "\n" );
    printf ( "  %9f  %9f  %9f  %9f  %9f\n",
      center[0], center[1], radius1, radius2, area );
  }

  return;
}
/*****************************************************************************80*/

void test09 ( )

/*****************************************************************************80

  Purpose:

    TEST09 tests CIRCLE_ANNULUS_SECTOR, CIRCLE_RT_SET, CIRCLE_RT_SUM.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 April 2008

  Author:

    John Burkardt
*/
{
  double as1;
  double as2;
  double as3;
  double as4;
  double center[2];
  int dim_num = 2;
  int j;
  int nc;
  char name[8];
  int num;
  int nr;
  int nr2;
  int nt;
  double pi = 3.141592653589793;
  double ra[5];
  double radius;
  double radius1a;
  double radius2a;
  double radius1b;
  double radius2b;
  double radius1c;
  double radius2c;
  double radius1d;
  double radius2d;
  double result1;
  double result2;
  int rule;
  double rw[5];
  double ta[20];
  int test_num = 4;
  double theta1a;
  double theta2a;
  double theta1b;
  double theta2b;
  double theta1c;
  double theta2c;
  double theta1d;
  double theta2d;
  double tw[20];
  double zw;

  nr = 5;

  rule = 9;
  circle_rt_size ( rule, &nr2, &nt, &nc );
  circle_rt_set ( rule, nr2, nt, nc, ra, rw, ta, tw, &zw );

  printf ( "\n" );
  printf ( "TEST09\n" );
  printf ( "  CIRCLE_ANNULUS_SECTOR estimates an integral in a\n" );
  printf ( "    circular annulus sector.\n" );
  printf ( "  CIRCLE_RT_SET sets an integration rule in a circle.\n" );
  printf ( "  CIRCLE_RT_SUM uses an integration rule in a circle.\n" );
  printf ( "\n" );
  printf ( "  To test CIRCLE_ANNULUS_SECTOR, we estimate an integral\n" );
  printf ( "  over 4 annular sectors that make up the unit circle,\n" );
  printf ( "  and add to get RESULT1.\n" );
  printf ( "\n" );
  printf ( "  We will also estimate the integral over the unit circle\n" );
  printf ( "  using CIRCLE_RT_SET and CIRCLE_RT_SUM to get RESULT2.\n" );
  printf ( "\n" );
  printf ( "  We will then compare RESULT1 and RESULT2.\n" );
  printf ( "\n" );
  printf ( "  CIRCLE_ANNULUS_SECTOR computations will use NR = %d\n", nr );
  printf ( "  CIRCLE_RT_SET/CIRCLE_RT_SUM will use rule %d\n", rule );
  printf ( "\n" );
  printf ( "  RESULT1 is the sum of Annulus Sector calculations.\n" );
  printf ( "  RESULT2 is for CIRCLE_RT_SET/CIRCLE_RT_SUM.\n" );
  printf ( "\n" );

  center[0] = 0.0;
  center[1] = 0.0;
  radius = 1.0;

  radius1a = 0.0;
  radius2a = 0.25;
  theta1a = 0.0;
  theta2a = 0.5 * pi;

  radius1b = 0.0;
  radius2b = 0.25;
  theta1b = 0.5 * pi;
  theta2b = 2.0 * pi;

  radius1c = 0.25;
  radius2c = 1.0;
  theta1c = 0.0;
  theta2c = 0.25 * pi;

  radius1d = 0.25;
  radius2d = 1.0;
  theta1d = 0.25 * pi;
  theta2d = 2.0 * pi;

  printf ( "\n" );
  printf ( "       F  Result1  Result2\n" );
  printf ( "\n" );

  num = function_2d_num ( );

  for ( j = 1; j <= num; j++ )
  {
    function_2d_index = j;

    as1 = circle_annulus_sector ( function_2d, center, radius1a, radius2a, theta1a,
      theta2a, nr );

    as2 = circle_annulus_sector ( function_2d, center, radius1b, radius2b, theta1b,
      theta2b, nr );

    as3 = circle_annulus_sector ( function_2d, center, radius1c, radius2c, theta1c,
      theta2c, nr );

    as4 = circle_annulus_sector ( function_2d, center, radius1d, radius2d, theta1d,
      theta2d, nr );

    result1 = as1 + as2 + as3 + as4;

    result2 = circle_rt_sum ( function_2d, center, radius, nr2, ra, rw, nt,
      ta, tw, zw );

    function_2d_name ( name );

    printf ( "  %8s  %14f  %14f\n", name, result1, result2 );
  }

  return;
}
/*****************************************************************************80*/

void test10 ( )

/*****************************************************************************80

  Purpose:

    TEST10 tests CIRCLE_CUM.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 April 2008

  Author:

    John Burkardt
*/
{
  double center[2] = { 0.0, 0.0 };
  int dim_num = 2;
  int i;
  int j;
  char name[8];
  int num;
  int order;
  double r = 3.0;

  printf ( "\n" );
  printf ( "TEST10\n" );
  printf ( "  CIRCLE_CUM approximates an integral over a circle.\n" );
  printf ( "\n" );
  printf ( "  We use radius R = %f\n", r );
  printf ( "  and center:\n" );
  printf ( "  CENTER = ( %f, %f).\n", center[0], center[1] );
  printf ( "\n" );

  printf ( "\n" );
  printf ( "    Order:      2             4              8            16\n" );
  printf ( "  F(X)\n" );
  printf ( "\n" );

  num = function_2d_num ( );

  for ( i = 1; i <= num; i++ )
  {
    function_2d_index = i;
    function_2d_name ( name );

    printf ( "  %8s", name );

    for ( j = 1; j <= 4; j++ )
    {
      order = i4_power ( 2, j );

      printf ( "%14f", circle_cum ( function_2d, center, r, order ) );
    }
    printf ( "\n" );
  }

  return;
}
/*****************************************************************************80*/

void test11 ( )

/*****************************************************************************80

  Purpose:

    TEST11 tests LENS_HALF_AREA_2D, CIRCLE_SECTOR_AREA_2D, CIRCLE_TRIANGLE_AREA_2D.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 April 2008

  Author:

    John Burkardt
*/
{
  double area1;
  double area2;
  double area3;
  int i;
  double pi = 3.141592653589793;
  double r;
  double theta1;
  double theta2;

  r = 1.0;

  printf ( "\n" );
  printf ( "TEST11\n" );
  printf ( "  LENS_HALF_AREA_2D computes the area of a\n" );
  printf ( "    circular half lens, defined by joining the endpoints\n" );
  printf ( "    of a circular arc.\n" );
  printf ( "  CIRCLE_SECTOR_AREA_2D computes the area of a\n" );
  printf ( "    circular sector, defined by joining the endpoints\n" );
  printf ( "    of a circular arc to the center.\n" );
  printf ( "  CIRCLE_TRIANGLE_AREA_2D computes the signed area of a\n" );
  printf ( "    triangle, defined by joining the endpoints\n" );
  printf ( "    of a circular arc and the center.\n" );
  printf ( "\n" );
  printf ( "      R       Theta1 Theta2        Sector       Triangle     Half Lens\n" );
  printf ( "\n" );

  for ( i = 0; i <= 12; i++ )
  {
    theta1 = 0.0;
    theta2 = ( double ) ( i ) * 2.0 * pi / 12.0;

    area1 = circle_sector_area_2d ( r, theta1, theta2 );

    area2 = circle_triangle_area_2d ( r, theta1, theta2 );

    area3 = lens_half_area_2d ( r, theta1, theta2 );

    printf ( "  %9f  %9f  %14f  %14f  %14f  %14f\n",
      r, theta1, theta2, area1, area2, area3 );
  }

  return;
}
/*****************************************************************************80*/

void test12 ( )

/*****************************************************************************80

  Purpose:

    TEST12 tests LENS_HALF_AREA_2D, LENS_HALF_H_AREA_2D, LENS_HALF_W_AREA_2D.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    06 April 2008

  Author:

    John Burkardt
*/
{
  double area1;
  double area2;
  double area3;
  double pi = 3.141592653589793;
  double h;
  int i;
  double r;
  double theta1;
  double theta2;
  double w;

  r = 50.0;

  printf ( "\n" );
  printf ( "TEST12\n" );
  printf ( "  For the area of a circular half lens,\n" );
  printf ( "  LENS_HALF_AREA_2D uses two angles;\n" );
  printf ( "  LENS_HALF_H_AREA_2D works from the height;\n" );
  printf ( "  LENS_HALF_W_AREA_2D works from the width.\n" );
  printf ( "\n" );
  printf ( "  The circle has radius R = %f\n", r );
  printf ( "\n" );
  printf ( "  THETA1 THETA2  H     W  Area(THETA) Area(H)  Area(W)\n" );
  printf ( "\n" );

  for ( i = 0; i <= 12; i++ )
  {
    theta1 = 0.0;
    theta2 = ( double ) ( i ) * 2.0 * pi / 12.0;
    w = 2.0 * r * sin ( 0.5 * ( theta2 - theta1 ) );
    h = r * ( 1.0 - cos ( 0.5 * ( theta2 - theta1 ) ) );

    area1 = lens_half_area_2d ( r, theta1, theta2 );

    area2 = lens_half_h_area_2d ( r, h );

    area3 = lens_half_w_area_2d ( r, w );

    printf ( "  %6f  %6f  %6f  %6f  %10f  %10f  %10f\n",
      theta1, theta2, h, w, area1, area2, area3 );
  }

  return;
}
/*****************************************************************************80*/

void test13 ( )

/*****************************************************************************80

  Purpose:

    TEST13 tests CIRCLE_SECTOR, CIRCLE_SECTOR_AREA_2D.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 April 2008

  Author:

    John Burkardt
*/
{
  double area;
  double center[2];
  double center_test[2*4] = {
    0.0, 0.0,
    0.0, 0.0,
    0.0, 0.0,
    0.0, 0.0 };
  int dim;
  int dim_num = 2;
  int i;
  int j;
  char name[8];
  int num;
  int nr;
  int nrhi = 5;
  int nrlo = 1;
  double pi = 3.141592653589793;
  double radius;
  double radius_test[4] = { 1.0, 2.0, 4.0, 8.0 };
  double result;
  int test_num = 4;
  double theta1;
  double theta1_test[4] = { 0.0, 0.0, 0.0, 0.0 };
  double theta2;
  double theta2_test[4] = { 2.0, 1.0, 0.5, 0.25 };

  printf ( "\n" );
  printf ( "TEST13\n" );
  printf ( "  CIRCLE_SECTOR_AREA_2D computes the area\n" );
  printf ( "    of a circular sector.\n" );
  printf ( "  CIRCLE_SECTOR estimates an integral \n" );
  printf ( "    in a circular sector.\n" );
  printf ( "\n" );
  printf ( "  The user can specify NR, the number of radial values\n" );
  printf ( "  used to approximated the integral.\n" );
  printf ( "\n" );
  printf ( "  In this test, computations will use values of NR\n" );
  printf ( "  from %d\n", nrlo );
  printf ( "  to   %d\n", nrhi );
  printf ( "\n" );

  for ( i = 0; i < test_num; i++ )
  {
    for ( dim = 0; dim < dim_num; dim++ )
    {
      center[dim] = center_test[dim+i*dim_num];
    }
    radius = radius_test[i];
    theta1 = theta1_test[i] * pi;
    theta2 = theta2_test[i] * pi;

    area = circle_sector_area_2d ( radius, theta1, theta2 );

    printf ( "\n" );
    printf ( "       CENTER      RADIUS  THETA1  THETA2  Area\n" );
    printf ( "\n" );
    printf ( "  %6f  %6f  %6f  %6f  %6f  %6f\n",
      center[0], center[1], radius, theta1, theta2, area );
    printf ( "\n" );
    printf ( "       F   " );
    for ( nr = nrlo; nr <= nrhi; nr++ )
    {
      printf ( "      %2d      ", nr );
    }
    printf ( "\n" );
    printf ( "\n" );

    num = function_2d_num ( );

    for ( j = 1; j <= num; j++ )
    {
      function_2d_index = j;

      function_2d_name ( name );

      printf ( "  %8s", name );

      for ( nr = nrlo; nr <= nrhi; nr++ )
      {
        result = circle_sector ( function_2d, center, radius, theta1, theta2,
          nr );
        printf ( "%14f", result );
      }
      printf ( "\n" );
    }
  }

  return;
}
/*****************************************************************************80*/

void test14 ( )

/*****************************************************************************80

  Purpose:

    TEST14 tests CIRCLE_SECTOR, CIRCLE_RT_SET, CIRCLE_RT_SUM.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 April 2008

  Author:

    John Burkardt
*/
{
  double area1;
  double area2;
  double area3;
  double center[2];
  double center_test[2*4] = {
    0.0, 0.0,
    0.0, 0.0,
    0.0, 0.0,
    0.0, 0.0 };
  int dim;
  int dim_num = 2;
  int i;
  int j;
  char name[8];
  int nc;
  int num;
  int nr;
  int nr2;
  int nt;
  double pi = 3.141592653589793;
  double ra[5];
  double radius;
  double radius_test[4] = { 1.0, 2.0, 4.0, 8.0 };
  double result1;
  double result2;
  double resulta;
  double resultb;
  int rule;
  double rw[5];
  double ta[20];
  int test_num = 4;
  double theta1;
  double theta1_test[4] = { 0.0, 0.0, 0.0, 0.0 };
  double theta2;
  double theta2_test[4] = { 2.0, 1.0, 0.5, 0.25 };
  double theta3;
  double tw[20];
  double zw;

  nr = 5;

  rule = 9;
  circle_rt_size ( rule, &nr2, &nt, &nc );
  circle_rt_set ( rule, nr2, nt, nc, ra, rw, ta, tw, &zw );

  printf ( "\n" );
  printf ( "TEST14\n" );
  printf ( "  CIRCLE_SECTOR estimates integrals in a circular sector.\n" );
  printf ( "  CIRCLE_RT_SET sets an integration rule in a circle.\n" );
  printf ( "  CIRCLE_RT_SUM uses an integration rule in a circle.\n" );
  printf ( "\n" );
  printf ( "  To test CIRCLE_SECTOR, we estimate an integral over\n" );
  printf ( "  a sector, and over its complement and add the results\n" );
  printf ( "  to get RESULT1.\n" );
  printf ( "\n" );
  printf ( "  We also estimate the integral over the whole circle\n" );
  printf ( "  using CIRCLE_RT_SET and CIRCLE_RT_SUM to get RESULT2.\n" );
  printf ( "\n" );
  printf ( "  We will then compare RESULT1 and RESULT2.\n" );
  printf ( "\n" );
  printf ( "  CIRCLE_SECTOR computations will use NR = %d\n", nr );
  printf ( "  CIRCLE_RT_SET/CIRCLE_RT_SUM will use rule %d\n", rule );
  printf ( "\n" );
  printf ( "  'Sector1' and 'Sector2' are the CIRCLE_SECTOR\n" );
  printf ( "  computations\n" );
  printf ( "  for the sector and its complement.\n" );
  printf ( "  'Sum' is the sum of Sector1 and Sector2.\n" );
  printf ( "  'Circle' is the computation for \n" );
  printf ( "  CIRCLE_RT_SET + CIRCLE_RT_SUM.\n" );
  printf ( "\n" );

  for ( i = 0; i < test_num; i++ )
  {
    for ( dim = 0; dim < dim_num; dim++ )
    {
      center[dim] = center_test[dim+i*dim_num];
    }
    radius = radius_test[i];

    theta1 = theta1_test[i] * pi;
    theta2 = theta2_test[i] * pi;
    theta3 = theta2 + 2.0 * pi - ( theta2 - theta1 );

    area1 = circle_sector_area_2d ( radius, theta1, theta2 );
    area2 = circle_sector_area_2d ( radius, theta2, theta3 );
    area3 = circle_area_2d ( radius );

    printf ( "\n" );
    printf ( "      CENTER       RADIUS   THETA1   THETA2   Area1   Area2  Circle\n" );
    printf ( "\n" );
    printf ( "  %7f  %7f  %7f  %7f  %7f  %7f  %7f  %7f\n",
      center[0], center[1], radius, theta1, theta2, area1, area2, area3 );
    printf ( "\n" );
    printf ( "       F   Sector1       Sector2         Sum         Circle\n" );
    printf ( "\n" );

    num = function_2d_num ( );

    for ( j = 1; j <= num; j++ )
    {
      function_2d_index = j;
      function_2d_name ( name );

      resulta = circle_sector ( function_2d, center, radius, theta1,
        theta2, nr );

      resultb = circle_sector ( function_2d, center, radius, theta2,
        theta3, nr );

      result1 = resulta + resultb;

      result2 = circle_rt_sum ( function_2d, center, radius, nr2, ra, rw,
        nt, ta, tw, zw );

      printf ( "  %8s  %14f  %14f  %14f  %14f\n",
        name, resulta, resultb, result1, result2 );
    }
  }
  return;
}
/*****************************************************************************80*/

void test15 ( )

/*****************************************************************************80

  Purpose:

    TEST15 tests CIRCLE_RT_SET and CIRCLE_RT_SUM.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 April 2008

  Author:

    John Burkardt
*/
{
  double center[2] = { 1.0, 1.0 };
  int dim_num = 2;
  int i;
  int ihi;
  int ilo;
  char name[8];
  int nc;
  int num;
  int nr;
  int nt;
  double pi = 3.141592653589793;
  double r = 1.0;
  double *ra;
  double result;
  double *rw;
  int rule;
  int rule_max = 9;
  double *ta;
  double *tw;
  double zw;

  printf ( "\n" );
  printf ( "TEST15\n" );
  printf ( "  For R, Theta product rules on the unit circle,\n" );
  printf ( "  CIRCLE_RT_SET sets a rule.\n" );
  printf ( "  CIRCLE_RT_SUM uses the rule in an arbitrary circle.\n" );
  printf ( "\n" );
  printf ( "  We use a radius %f\n", r );
  printf ( "  and center:\n" );
  printf ( "  CENTER = ( %f, %f )\n", center[0], center[1] );
  printf ( "\n" );

  for ( ilo = 1; ilo <= rule_max; ilo = ilo + 5 )
  {
    ihi = i4_min ( ilo +  4, rule_max );

    printf ( "\n" );
    printf ( "  Rule:  " );
    for ( rule = ilo; rule <= ihi; rule++ )
    {
      printf ( "       %7d", rule );
    }
    printf ( "\n" );
    printf ( "Function\n" );

    num = function_2d_num ( );

    for ( i = 1; i <= num; i++ )
    {
      function_2d_index = i;
      function_2d_name ( name );
      printf ( "  %8s", name );

      for ( rule = ilo; rule <= ihi; rule++ )
      {
        circle_rt_size ( rule, &nr, &nt, &nc );

        ra = ( double * ) malloc ( nr * sizeof ( double ) );
        rw = ( double * ) malloc ( nr * sizeof ( double ) );
        ta = ( double * ) malloc ( nt * sizeof ( double ) );
        tw = ( double * ) malloc ( nt * sizeof ( double ) );

        circle_rt_set ( rule, nr, nt, nc, ra, rw, ta, tw, &zw );

        result = circle_rt_sum ( function_2d, center, r, nr, ra, rw, nt, ta,
          tw, zw );
        printf ( "%14f", result );

        free ( ra );
        free ( rw );
        free ( ta );
        free ( tw );
      }
      printf ( "\n" );
    }
  }
  return;
}
/*****************************************************************************80*/

void test16 ( )

/*****************************************************************************80

  Purpose:

    TEST16 tests CIRCLE_XY_SET and CIRCLE_XY_SUM.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 March 2008

  Author:

    John Burkardt
*/
{
  double center[2] = { 1.0, 1.0 };
  int dim_num = 2;
  int i;
  int ihi;
  int ilo;
  char name[8];
  int num;
  int order;
  double r;
  double result;
  int rule;
  int rule_max = 13;
  double *weight;
  double *xtab;
  double *ytab;

  r = 1.0;
  printf ( "\n" );
  printf ( "TEST16\n" );
  printf ( "  CIRCLE_XY_SET sets a quadrature rule\n" );
  printf ( "    for the unit circle.\n" );
  printf ( "  CIRCLE_XY_SUM evaluates the quadrature rule\n" );
  printf ( "    in an arbitrary circle.\n" );
  printf ( "\n" );
  printf ( "  We use a radius %f\n", r );
  printf ( "  and center:\n" );
  printf ( "  CENTER = (%f, %f).\n", center[0], center[1] );
  printf ( "\n" );

  for ( ilo = 1; ilo <= rule_max; ilo = ilo + 5 )
  {
    ihi = i4_min ( ilo +  4, rule_max );

    printf ( "\n" );
    printf ( "  Rule:  " );
    for ( rule = ilo; rule <= ihi; rule++ )
    {
      printf ( "       %7d", rule );
    }
    printf ( "\n" );
    printf ( " 'Function\n" );

    num = function_2d_num ( );

    for ( i = 1; i <= num; i++ )
    {
      function_2d_index = i;
      function_2d_name ( name );

      printf ( "  %8s", name );

      for ( rule = ilo; rule <= ihi; rule++ )
      {
        order = circle_xy_size ( rule );

        xtab = ( double * ) malloc ( order * sizeof ( double ) );
        ytab = ( double * ) malloc ( order * sizeof ( double ) );
        weight = ( double * ) malloc ( order * sizeof ( double ) );

        circle_xy_set ( rule, order, xtab, ytab, weight );

        result = circle_xy_sum ( function_2d, center, r, order, xtab, ytab,
          weight );

        printf ( "%14f", result );

        free ( weight );
        free ( xtab );
        free ( ytab );
      }
      printf ( "\n" );
    }
  }
  return;
}
/*****************************************************************************80*/

void test163 ( )

/*****************************************************************************80

  Purpose:

    TEST163 tests the rules for CN with Gegenbauer weight on monomials.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 January 2010

  Author:

    John Burkardt
*/
{
# define TEST_NUM 5

  double alpha;
  double alpha_test[TEST_NUM] = { -0.5, 0.0, 0.5, 1.0, 1.5 };
  int *expon;
  int i;
  int n;
  int test;

  printf ( "\n" );
  printf ( "TEST163\n" );
  printf ( "  Demonstrate the use of quadrature rules for the region\n" );
  printf ( "  CN_GEG, that is, the hypercube [-1,+1]^N, with the\n" );
  printf ( "  weight W(ALPHA;X) = product ( 1 <= I <= N )\n" );
  printf ( "    (1-X(I)^2)^ALPHA\n" );
  printf ( "\n" );
  printf ( "  We use the formulas to integrate various monomials of\n" );
  printf ( "  the form X(1)^E(1) * X(2)^E(2) * ... X(N)^E(N)\n" );
  printf ( "  and compare to the exact integral.\n" );
  printf ( "\n" );
  printf ( "  The precision of each formula is known, and we only use\n" );
  printf ( "  a formula if its precision indicates it should be able to\n" );
  printf ( "  produce an exact result.\n" );

  for ( n = 1; n <= 6; n++ )
  {
    expon = ( int * ) malloc ( n * sizeof ( int ) );

    for ( test = 0; test < TEST_NUM; test++ )
    {
      alpha = alpha_test[test];

      i4vec_zero ( n, expon );
      cn_geg_test ( n, alpha, expon );
    }

    for ( test = 0; test < TEST_NUM; test++ )
    {
      alpha = alpha_test[test];

      i4vec_zero ( n, expon );
      expon[n-1] = 1;
      cn_geg_test ( n, alpha, expon );
    }

    if ( 2 <= n )
    {
      for ( test = 0; test < TEST_NUM; test++ )
      {
        alpha = alpha_test[test];

        i4vec_zero ( n, expon );
        expon[0] = 1;
        expon[1] = 1;
        cn_geg_test ( n, alpha, expon );
      }
    }

    for ( test = 0; test < TEST_NUM; test++ )
    {
      alpha = alpha_test[test];

      i4vec_zero ( n, expon );
      expon[0] = 2;
      cn_geg_test ( n, alpha, expon );

    }
    free ( expon );
  }

  return;
# undef TEST_NUM
}
/*****************************************************************************80*/

void cn_geg_test ( int n, double alpha, int expon[] )

/*****************************************************************************80

  Purpose:

    CN_GEG_TEST tests the rules for CN with Gegenbauer weight on a monomial.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 March 2010

  Author:

    John Burkardt
*/
{
  double c1;
  int d;
  double delta0;
  double err;
  double exact;
  double gamma0;
  int i;
  int o;
  int option;
  int p;
  double pi = 3.141592653589793;
  double quad;
  double *v;
  double volume_1d;
  double *w;
  double *x;

  printf ( "\n" );
  printf ( "  N = %d\n", n );
  printf ( "  ALPHA = %f\n", alpha );
  printf ( "  EXPON = " );
  for ( i = 0; i < n; i++ )
  {
    printf ( "%4d", expon[i] );
  }
  printf ( "\n" );
  d = i4vec_sum ( n, expon );
  printf ( "  Degree = %d\n", d );
  printf ( "\n" );

  exact = cn_geg_monomial_integral ( n, alpha, expon );

  p = 0;

  if ( d <= p )
  {
    o = cn_geg_00_1_size ( n, alpha );
    x = ( double * ) malloc ( n * o * sizeof ( double ) );
    w = ( double * ) malloc ( o * sizeof ( double ) );
    cn_geg_00_1 ( n, alpha, o, x, w );
    v = monomial_value ( n, o, x, expon );
    quad = r8vec_dot_product ( o, w, v );
    err = r8_abs ( quad - exact );
    printf ( "  CN_GEG_00_1:     %6d  %14f  %14f\n", o, quad, err );
    free ( v );
    free ( w );
    free ( x );
  }

  p = 1;

  if ( d <= p )
  {
    o = cn_geg_01_1_size ( n, alpha );
    x = ( double * ) malloc ( n * o * sizeof ( double ) );
    w = ( double * ) malloc ( o * sizeof ( double ) );
    cn_geg_01_1 ( n, alpha, o, x, w );
    v = monomial_value ( n, o, x, expon );
    quad = r8vec_dot_product ( o, w, v );
    err = r8_abs ( quad - exact );
    printf ( "  CN_GEG_01_1:     %6d  %14f  %14f\n", o, quad, err );
    free ( v );
    free ( w );
    free ( x );
  }

  p = 2;

  if ( d <= p )
  {
    o = cn_geg_02_xiu_size ( n, alpha );
    x = ( double * ) malloc ( n * o * sizeof ( double ) );
    w = ( double * ) malloc ( o * sizeof ( double ) );
    cn_geg_02_xiu ( n, alpha, o, x, w );
    v = monomial_value ( n, o, x, expon );
    quad = r8vec_dot_product ( o, w, v );
    err = r8_abs ( quad - exact );
    printf ( "  CN_GEG_02_XIU:   %6d  %14f  %14f\n", o, quad, err );
    free ( v );
    free ( w );
    free ( x );

    o = gw_02_xiu_size ( n );
    gamma0 = 1.0;
    delta0 = 0.0;
    c1 = 1.0 / ( 2.0 * alpha + 3.0 );
    volume_1d = sqrt ( pi ) * r8_gamma ( alpha + 1.0 )
      / r8_gamma ( alpha + 1.5 );
    x = ( double * ) malloc ( n * o * sizeof ( double ) );
    w = ( double * ) malloc ( o * sizeof ( double ) );
    gw_02_xiu ( n, o, gamma0, delta0, c1, volume_1d, x, w );
    v = monomial_value ( n, o, x, expon );
    quad = r8vec_dot_product ( o, w, v );
    err = r8_abs ( quad - exact );
    printf ( "  GW_02_XIU:       %6d  %14f  %14f\n", o, quad, err );
    free ( v );
    free ( w );
    free ( x );
  }

  p = 3;

  if ( d <= p )
  {
    o = cn_geg_03_xiu_size ( n, alpha );
    x = ( double * ) malloc ( n * o * sizeof ( double ) );
    w = ( double * ) malloc ( o * sizeof ( double ) );
    cn_geg_03_xiu ( n, alpha, o, x, w );
    v = monomial_value ( n, o, x, expon );
    quad = r8vec_dot_product ( o, w, v );
    err = r8_abs ( quad - exact );
    printf ( "  CN_GEG_03_XIU:   %6d  %14f  %14f\n", o, quad, err );
    free ( v );
    free ( w );
    free ( x );
  }

  printf ( "  EXACT                    %14f\n", exact );

  return;
}
/*****************************************************************************80*/

void test165 ( )

/*****************************************************************************80

  Purpose:

    TEST165 tests the rules for CN with Jacobi weight on monomials.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 January 2010

  Author:

    John Burkardt
*/
{
# define TEST_NUM 4

  double alpha;
  double alpha_test[TEST_NUM] = { 0.0, 1.0, 0.0, 0.5 };
  double beta;
  double beta_test[TEST_NUM] = { 0.0, 0.0, 2.0, 1.5 };
  int *expon;
  int i;
  int n;
  int test;

  printf ( "\n" );
  printf ( "TEST165\n" );
  printf ( "  Demonstrate the use of quadrature rules for the region\n" );
  printf ( "  CN_JAC, that is, the hypercube [-1,+1]^N, with the\n" );
  printf ( "  weight W(ALPHA,BETA;X) = product ( 1 <= I <= N )\n" );
  printf ( "    (1-X(I))^ALPHA (1+X(I))^BETA\n" );
  printf ( "\n" );
  printf ( "  We use the formulas to integrate various monomials of\n" );
  printf ( "  the form X(1)^E(1) * X(2)^E(2) * ... X(N)^E(N)\n" );
  printf ( "  and compare to the exact integral.\n" );
  printf ( "\n" );
  printf ( "  The precision of each formula is known, and we only use\n" );
  printf ( "  a formula if its precision indicates it should be able to\n" );
  printf ( "  produce an exact result.\n" );

  for ( n = 1; n <= 6; n++ )
  {
    expon = ( int * ) malloc ( n * sizeof ( int ) );

    for ( test = 0; test < TEST_NUM; test++ )
    {
      alpha = alpha_test[test];
      beta  = beta_test[test];

      i4vec_zero ( n, expon );
      cn_jac_test ( n, alpha, beta, expon );
    }

    for ( test = 0; test < TEST_NUM; test++ )
    {
      alpha = alpha_test[test];
      beta  = beta_test[test];

      i4vec_zero ( n, expon );
      expon[n-1] = 1;
      cn_jac_test ( n, alpha, beta, expon );
    }

    if ( 2 <= n )
    {
      for ( test = 0; test < TEST_NUM; test++ )
      {
        alpha = alpha_test[test];
        beta  = beta_test[test];

        i4vec_zero ( n, expon );
        expon[0] = 1;
        expon[1] = 1;
        cn_jac_test ( n, alpha, beta, expon );
      }
    }

    for ( test = 0; test < TEST_NUM; test++ )
    {
      alpha = alpha_test[test];
      beta  = beta_test[test];

      i4vec_zero ( n, expon );
      expon[0] = 2;
      cn_jac_test ( n, alpha, beta, expon );

    }
    free ( expon );
  }

  return;
# undef TEST_NUM
}
/*****************************************************************************80*/

void cn_jac_test ( int n, double alpha, double beta, int expon[] )

/*****************************************************************************80

  Purpose:

    CN_JAC_TEST tests the rules for CN with Jacobi weight on a monomial.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 March 2010

  Author:

    John Burkardt
*/
{
  double c1;
  int d;
  double delta0;
  double err;
  double exact;
  double gamma0;
  int i;
  int o;
  int option;
  int p;
  double quad;
  double *v;
  double volume_1d;
  double *w;
  double *x;

  printf ( "\n" );
  printf ( "  N = %d\n", n );
  printf ( "  ALPHA = %f\n", alpha );
  printf ( "  BETA =  %f\n", beta );
  printf ( "  EXPON = " );
  for ( i = 0; i < n; i++ )
  {
    printf ( "%4d", expon[i] );
  }
  printf ( "\n" );
  d = i4vec_sum ( n, expon );
  printf ( "  Degree = %d\n", d );
  printf ( "\n" );

  exact = cn_jac_monomial_integral ( n, alpha, beta, expon );

  p = 0;

  if ( d <= p )
  {
    o = cn_jac_00_1_size ( n, alpha, beta );
    x = ( double * ) malloc ( n * o * sizeof ( double ) );
    w = ( double * ) malloc ( o * sizeof ( double ) );
    cn_jac_00_1 ( n, alpha, beta, o, x, w );
    v = monomial_value ( n, o, x, expon );
    quad = r8vec_dot_product ( o, w, v );
    err = r8_abs ( quad - exact );
    printf ( "  CN_JAC_00_1:     %6d  %14f  %14f\n", o, quad, err );
    free ( v );
    free ( w );
    free ( x );
  }

  p = 1;

  if ( d <= p )
  {
    o = cn_jac_01_1_size ( n, alpha, beta );
    x = ( double * ) malloc ( n * o * sizeof ( double ) );
    w = ( double * ) malloc ( o * sizeof ( double ) );
    cn_jac_01_1 ( n, alpha, beta, o, x, w );
    v = monomial_value ( n, o, x, expon );
    quad = r8vec_dot_product ( o, w, v );
    err = r8_abs ( quad - exact );
    printf ( "  CN_JAC_01_1:     %6d  %14f  %14f\n", o, quad, err );
    free ( v );
    free ( w );
    free ( x );
  }

  p = 2;

  if ( d <= p )
  {
    o = cn_jac_02_xiu_size ( n, alpha, beta );
    x = ( double * ) malloc ( n * o * sizeof ( double ) );
    w = ( double * ) malloc ( o * sizeof ( double ) );
    cn_jac_02_xiu ( n, alpha, beta, o, x, w );
    v = monomial_value ( n, o, x, expon );
    quad = r8vec_dot_product ( o, w, v );
    err = r8_abs ( quad - exact );
    printf ( "  CN_JAC_02_XIU:   %6d  %14f  %14f\n", o, quad, err );
    free ( v );
    free ( w );
    free ( x );

    o = gw_02_xiu_size ( n );
    gamma0 = ( alpha + beta + 2.0 ) / 2.0;
    delta0 = ( alpha - beta ) / 2.0;
    c1 = 2.0 * ( alpha + 1.0 ) * ( beta + 1.0 ) / ( alpha + beta + 3.0 )
      / ( alpha + beta + 2.0 );
    volume_1d = pow ( 2.0, alpha + beta + 1.0 ) * r8_gamma ( alpha + 1.0 )
      * r8_gamma ( beta + 1.0 ) / ( alpha + beta + 1.0 ) / r8_gamma ( alpha + beta + 1.0 );
    x = ( double * ) malloc ( n * o * sizeof ( double ) );
    w = ( double * ) malloc ( o * sizeof ( double ) );
    gw_02_xiu ( n, o, gamma0, delta0, c1, volume_1d, x, w );
    v = monomial_value ( n, o, x, expon );
    quad = r8vec_dot_product ( o, w, v );
    err = r8_abs ( quad - exact );
    printf ( "  GW_02_XIU:       %6d  %14f  %14f\n", o, quad, err );
    free ( v );
    free ( w );
    free ( x );
  }
  printf ( "  EXACT                    %14f\n", exact );

  return;
}
/*****************************************************************************80*/

void test167 ( )

/*****************************************************************************80

  Purpose:

    TEST167 tests the rules for CN with Legendre weight on monomials.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    03 March 2010

  Author:

    John Burkardt
*/
{
  int *expon;
  int i;
  int n;
  int test;

  printf ( "\n" );
  printf ( "TEST167\n" );
  printf ( "  Demonstrate the use of quadrature rules for the region\n" );
  printf ( "  CN_LEG, that is, the hypercube [-1,+1]^N, with the\n" );
  printf ( "  Legendre weight W(X) = 1.\n" );
  printf ( "\n" );
  printf ( "  We use the formulas to integrate various monomials of\n" );
  printf ( "  the form X(1)^E(1) * X(2)^E(2) * ... X(N)^E(N)\n" );
  printf ( "  and compare to the exact integral.\n" );
  printf ( "\n" );
  printf ( "  The precision of each formula is known, and we only use\n" );
  printf ( "  a formula if its precision indicates it should be able to\n" );
  printf ( "  produce an exact result.\n" );

  for ( n = 1; n <= 6; n++ )
  {
    expon = ( int * ) malloc ( n * sizeof ( int ) );

    i4vec_zero ( n, expon );
    cn_leg_test ( n, expon );

    i4vec_zero ( n, expon );
    expon[n-1] = 1;
    cn_leg_test ( n, expon );

    if ( 2 <= n )
    {
      i4vec_zero ( n, expon );
      expon[0] = 1;
      expon[1] = 1;
      cn_leg_test ( n, expon );
    }

    i4vec_zero ( n, expon );
    expon[0] = 2;
    cn_leg_test ( n, expon );

    i4vec_zero ( n, expon );
    expon[0] = 3;
    cn_leg_test ( n, expon );

    i4vec_zero ( n, expon );
    expon[n-1] = 4;
    cn_leg_test ( n, expon );

    if ( 2 <= n )
    {
      i4vec_zero ( n, expon );
      expon[0] = 3;
      expon[1] = 2;
      cn_leg_test ( n, expon );
    }

    free ( expon );
  }

  return;
}
/*****************************************************************************80*/

void cn_leg_test ( int n, int expon[] )

/*****************************************************************************80

  Purpose:

    CN_LEG_TEST tests the rules for CN with Legendre weight on a monomial.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    03 March 2010

  Author:

    John Burkardt
*/
{
  double c1;
  int d;
  double delta0;
  double err;
  double exact;
  double gamma0;
  int i;
  int o;
  int option;
  int p;
  double quad;
  double *v;
  double volume_1d;
  double *w;
  double *x;

  printf ( "\n" );
  printf ( "  N = %d\n", n );
  printf ( "  EXPON = " );
  for ( i = 0; i < n; i++ )
  {
    printf ( "%4d", expon[i] );
  }
  printf ( "\n" );
  d = i4vec_sum ( n, expon );
  printf ( "  Degree = %d\n", d );
  printf ( "\n" );

  exact = cn_leg_monomial_integral ( n, expon );

  p = 1;

  if ( d <= p )
  {
    o = cn_leg_01_1_size ( n );
    x = ( double * ) malloc ( n * o * sizeof ( double ) );
    w = ( double * ) malloc ( o * sizeof ( double ) );
    cn_leg_01_1 ( n, o, x, w );
    v = monomial_value ( n, o, x, expon );
    quad = r8vec_dot_product ( o, w, v );
    err = r8_abs ( quad - exact );
    printf ( "  CN_LEG_01_1:     %6d  %14f  %14f\n", o, quad, err );
    free ( v );
    free ( w );
    free ( x );
  }

  p = 2;

  if ( d <= p )
  {
    o = cn_leg_02_xiu_size ( n );
    x = ( double * ) malloc ( n * o * sizeof ( double ) );
    w = ( double * ) malloc ( o * sizeof ( double ) );
    cn_leg_02_xiu ( n, o, x, w );
    v = monomial_value ( n, o, x, expon );
    quad = r8vec_dot_product ( o, w, v );
    err = r8_abs ( quad - exact );
    printf ( "  CN_LEG_02_XIU:   %6d  %14f  %14f\n", o, quad, err );
    free ( v );
    free ( w );
    free ( x );

    o = gw_02_xiu_size ( n );
    gamma0 = 1.0;
    delta0 = 0.0;
    c1 = 1.0 / 3.0;
    volume_1d = 2.0;
    x = ( double * ) malloc ( n * o * sizeof ( double ) );
    w = ( double * ) malloc ( o * sizeof ( double ) );
    gw_02_xiu ( n, o, gamma0, delta0, c1, volume_1d, x, w );
    v = monomial_value ( n, o, x, expon );
    quad = r8vec_dot_product ( o, w, v );
    err = r8_abs ( quad - exact );
    printf ( "  GW_02_XIU:       %6d  %14f  %14f\n", o, quad, err );
    free ( v );
    free ( w );
    free ( x );
  }

  p = 3;

  if ( d <= p )
  {
    o = cn_leg_03_1_size ( n );
    x = ( double * ) malloc ( n * o * sizeof ( double ) );
    w = ( double * ) malloc ( o * sizeof ( double ) );
    cn_leg_03_1 ( n, o, x, w );
    v = monomial_value ( n, o, x, expon );
    quad = r8vec_dot_product ( o, w, v );
    err = r8_abs ( quad - exact );
    printf ( "  CN_LEG_03_1:     %6d  %14f  %14f\n", o, quad, err );
    free ( v );
    free ( w );
    free ( x );

    o = cn_leg_03_xiu_size ( n );
    x = ( double * ) malloc ( n * o * sizeof ( double ) );
    w = ( double * ) malloc ( o * sizeof ( double ) );
    cn_leg_03_xiu ( n, o, x, w );
    v = monomial_value ( n, o, x, expon );
    quad = r8vec_dot_product ( o, w, v );
    err = r8_abs ( quad - exact );
    printf ( "  CN_LEG_03_XIU:   %6d  %14f  %14f\n", o, quad, err );
    free ( v );
    free ( w );
    free ( x );
  }

  p = 5;

  if ( d <= p )
  {
    if ( 4 <= n && n <= 6 )
    {
      o = cn_leg_05_1_size ( n );
      x = ( double * ) malloc ( n * o * sizeof ( double ) );
      w = ( double * ) malloc ( o * sizeof ( double ) );
      option = 1;
      cn_leg_05_1 ( n, option, o, x, w );
      v = monomial_value ( n, o, x, expon );
      quad = r8vec_dot_product ( o, w, v );
      err = r8_abs ( quad - exact );
      printf ( "  CN_LEG_05_1(1):  %6d  %14f  %14f\n", o, quad, err );
      free ( v );
      free ( w );
      free ( x );
    }
    if ( 4 <= n && n <= 5 )
    {
      o = cn_leg_05_1_size ( n );
      x = ( double * ) malloc ( n * o * sizeof ( double ) );
      w = ( double * ) malloc ( o * sizeof ( double ) );
      option = 2;
      cn_leg_05_1 ( n, option, o, x, w );
      v = monomial_value ( n, o, x, expon );
      quad = r8vec_dot_product ( o, w, v );
      err = r8_abs ( quad - exact );
      printf ( "  CN_LEG_05_1(2):  %6d  %14f  %14f\n", o, quad, err );
      free ( v );
      free ( w );
      free ( x );
    }

    if ( 2 <= n )
    {
      o = cn_leg_05_2_size ( n );
      x = ( double * ) malloc ( n * o * sizeof ( double ) );
      w = ( double * ) malloc ( o * sizeof ( double ) );
      cn_leg_05_2 ( n, o, x, w );
      v = monomial_value ( n, o, x, expon );
      quad = r8vec_dot_product ( o, w, v );
      err = r8_abs ( quad - exact );
      printf ( "  CN_LEG_05_2:     %6d  %14f  %14f\n", o, quad, err );
      free ( v );
      free ( w );
      free ( x );
    }
  }

  printf ( "  EXACT                    %14f\n", exact );

  return;
}
/*****************************************************************************80*/

void test17 ( )

/*****************************************************************************80

  Purpose:

    TEST17 tests CONE_UNIT_3D, CONE_VOLUME_3D.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 April 2008

  Author:

    John Burkardt
*/
{
  double h;
  int i;
  char name[8];
  int num;
  double r = 1.0;
  double result;

  h = 1.0;

  printf ( "\n" );
  printf ( "TEST17\n" );
  printf ( "  CONE_UNIT_3D approximates integrals in a unit cone.\n" );
  printf ( "\n" );
  printf ( "  Volume = %f\n", cone_volume_3d ( r, h ) );
  printf ( "\n" );
  printf ( "    F(X)    CONE_3D\n" );
  printf ( "\n" );

  num = function_3d_num ( );

  for ( i = 1; i <= num; i++ )
  {
    function_3d_index = i;
    function_3d_name ( name );

    result = cone_unit_3d ( function_3d );

    printf ( "  %8s  %14f\n", name, result );
  }
  return;
}
/*****************************************************************************80*/

void test18 ( )

/*****************************************************************************80

  Purpose:

    TEST18 tests CUBE_SHELL_ND.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 April 2008

  Author:

    John Burkardt
*/
{
  int i;
  int n;
  int n_max = 4;
  char name[8];
  int num;
  double r1;
  double r1_test[2] = { 0.0, 1.0 };
  double r2;
  double r2_test[2] = { 1.0, 2.0 };
  double result;
  int test;
  int test_num = 2;

  printf ( "\n" );
  printf ( "TEST18\n" );
  printf ( "  CUBE_SHELL_ND approximates integrals in a\n" );
  printf ( "    cubical shell in ND.\n" );

  for ( test = 0; test < test_num; test++ )
  {
    r1 = r1_test[test];
    r2 = r2_test[test];

    printf ( "\n" );
    printf ( "  Inner radius = %f\n", r1 );
    printf ( "  Outer radius = %f\n", r2 );
    printf ( "\n" );

    for ( n = 2; n <= n_max; n++ )
    {
      printf ( "\n" );
      printf ( "  Spatial dimension N = %d\n", n );
      printf ( "  Volume = %f\n", cube_shell_volume_nd ( n, r1, r2 ) );
      printf ( "\n" );
      printf ( "    F(X)      CUBE_SHELL_ND\n" );
      printf ( "\n" );

      num = function_nd_num ( );

      for ( i = 1; i <= num; i++ )
      {
        function_nd_index = i;
        function_nd_name ( name );

        result = cube_shell_nd ( function_nd, n, r1, r2 );

        printf ( "  %8s  %14f\n", name, result );
      }
    }
  }
  return;
}
/*****************************************************************************80*/

void test19 ( )

/*****************************************************************************80

  Purpose:

    TEST19 tests CUBE_UNIT_3D, QMULT_3D, RECTANGLE_3D.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 March 2008

  Author:

    John Burkardt
*/
{
  double a1;
  double a[3];
  double b1;
  double b[3];
  int i;
  int n = 3;
  char name[8];
  int num;
  double result1;
  double result2;
  double result3;

  a1 = -1.0;
  b1 = +1.0;

  a[0]= -1.0;
  a[1] = -1.0;
  a[2] = -1.0;
  b[0] = 1.0;
  b[1] = 1.0;
  b[2] = 1.0;

  printf ( "\n" );
  printf ( "TEST19\n" );
  printf ( "  CUBE_UNIT_3D approximates integrals\n" );
  printf ( "    in the unit cube in 3D.\n" );
  printf ( "  QMULT_3D approximates triple integrals.\n" );
  printf ( "  RECTANGLE_3D approximates integrals\n" );
  printf ( "    in a rectangular block.\n" );
  printf ( "\n" );
  printf ( "    F(X)      CUBE_UNIT_3D    QMULT_3D        RECTANGLE_3D\n" );
  printf ( "\n" );

  num = function_3d_num ( );

  for ( i = 1; i <= num; i++ )
  {
    function_3d_index = i;
    function_3d_name ( name );

    result1 = cube_unit_3d ( function_3d );
    result2 = qmult_3d ( function_3d, a1, b1, fu18, fl18, fu28, fl28 );
    result3 = rectangle_3d ( function_3d, a, b );

    printf ( "  %8s  %14f  %14f  %14f\n", name, result1, result2, result3 );
  }
  return;
}
/*****************************************************************************80*/

void test20 ( )

/*****************************************************************************80

  Purpose:

    TEST20 tests CUBE_UNIT_ND.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 April 2008

  Author:

    John Burkardt
*/
{
  int i;
  int i_test;
  int k;
  int k2;
  int khi;
  int klo;
  int k_test[2] = { 10, 5 };
  int max_k = 10;
  int max_test = 2;
  int n;
  int n_test[2] = { 2, 3 };
  char name[8];
  int num;
  double qa[10];
  double qb[10];

  printf ( "\n" );
  printf ( "TEST20\n" );
  printf ( "  CUBE_UNIT_ND approximates integrals inside \n" );
  printf ( "    the unit cube in ND.\n" );

  for ( i_test = 0; i_test < max_test; i_test++ )
  {
    n = n_test[i_test];
    k = k_test[i_test];
    printf ( "\n" );
    printf ( "\n" );
    printf ( "  Spatial dimension N = %d\n", n );
    printf ( "  Value of K = %d\n", k );
    printf ( "\n" );
    printf ( "    F(X)    CUBE_UNIT_ND\n" );
    printf ( "\n" );

    num = function_nd_num ( );

    for ( i = 1; i <= num; i++ )
    {
      function_nd_index = i;
      function_nd_name ( name );

      cube_unit_nd ( function_nd, qa, qb, n, k );

      for ( klo = 0; klo <= k - 1; klo = klo + 5 )
      {
        khi = i4_min ( klo + 4, k-1 );
        if ( klo == 0 )
        {
          printf ( "  %8s", name );
        }
        else
        {
          printf ( "          " );
        }
        for ( k2 = klo; k2 <= khi; k2++ )
        {
          printf ( "%14f", qa[k2] );
        }
        printf ( "\n" );
      }
      for ( klo = 0; klo <= k - 1; klo = klo + 5 )
      {
        khi = i4_min ( klo + 4, k - 1 );
        printf ( "         " );
        for ( k2 = klo; k2 <= khi; k2++ )
        {
          printf ( "%14f", qb[k2] );
        }
        printf ( "\n" );
      }
    }
  }
  return;
}
/*****************************************************************************80*/

void test205 ( )

/*****************************************************************************80

  Purpose:

    TEST205 tests ELLIPSE_AREA_2D, ELLIPSE_CIRCUMFERENCE_2D, ELLIPSE_ECCENTRICITY_2D.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    06 April 2008

  Author:

    John Burkardt
*/
{
  double area;
  double e;
  int i;
  double p;
  double r1;
  double r2;
  int seed = 123456789;

  printf ( "\n" );
  printf ( "TEST205\n" );
  printf ( "  ELLIPSE_AREA_2D returns the area of an ellipse.\n" );
  printf ( "  ELLIPSE_ECCENTRICITY_2D returns the\n" );
  printf ( "    eccentricity of an ellipse.\n" );
  printf ( "  ELLIPSE_CIRCUMFERENCE_2D returns the \n" );
  printf ( "    circumference of an ellipse.\n" );
  printf ( "\n" );
  printf ( "        R1        R2        E         Circum    Area\n" );
  printf ( "\n" );

  for ( i = 1; i <= 5; i++ )
  {
    if ( i == 1 )
    {
      r1 = 25.0;
      r2 = 20.0;
    }
    else
    {
      r1 = r8_uniform_01 ( &seed );
      r2 = r8_uniform_01 ( &seed );
    }

    e = ellipse_eccentricity_2d ( r1, r2 );
    p = ellipse_circumference_2d ( r1, r2 ) ;
    area = ellipse_area_2d ( r1, r2 );

    printf ( "  %10f  %10f  %10f  %10f  %10f\n", r1, r2, e, p, area );
  }
  printf ( "\n" );
  printf ( "  (For the first example, \n" );
  printf ( "  the eccentricity should be 0.6,\n" );
  printf ( "  the circumference should be about 141.8).\n" );

  return;
}
/*****************************************************************************80*/

void test207 ( )

/*****************************************************************************80

  Purpose:

    TEST207 tests the Stroud EN_R2 rules on monomials.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    23 April 2013

  Author:

    John Burkardt
*/
{
  int *expon;
  int expon_dim;
  int i;
  int n;

  printf ( "\n" );
  printf ( "TEST207\n" );
  printf ( "  Demonstrate the use of Stroud rules for the region\n" );
  printf ( "  EN_R2, that is, all of N-dimensional space, with the\n" );
  printf ( "  weight function W(X) = exp ( - X1^2 - X2^2 ... -XN^2 )\n" );
  printf ( "\n" );
  printf ( "  We use the formulas to integrate various monomials of\n" );
  printf ( "  the form X1^EXPON1 * X2^EXPON2 * ... XN^EXPONN\n" );
  printf ( "  and compare to the exact integral.\n" );
  printf ( "\n" );
  printf ( "  The precision of each formula is known, and we only use\n" );
  printf ( "  a formula if its precision indicates it should be able to\n" );
  printf ( "  produce an exact result.\n" );

  for ( n = 1; n <= 7; n++ )
  {
    expon_dim = i4_max ( n, 2 );
    expon = ( int * ) malloc ( expon_dim * sizeof ( int ) );

    for ( i = 0; i < n; i++ )
    {
      expon[i] = 0;
    }
    en_r2_test ( n, expon );

    for ( i = 0; i < n; i++ )
    {
      expon[i] = 0;
    }
    expon[0] = 2;
    en_r2_test ( n, expon );

    for ( i = 0; i < n; i++ )
    {
      expon[i] = 0;
    }
    expon[1] = 4;
    en_r2_test ( n, expon );

    for ( i = 0; i < n; i++ )
    {
      expon[i] = 0;
    }
    i = 3 % n;
    expon[i] = 6;
    en_r2_test ( n, expon );

    for ( i = 0; i < n; i++ )
    {
      expon[i] = 0;
    }
    expon[0] = 2;
    expon[1] = 4;
    en_r2_test ( n, expon );

    for ( i = 0; i < n; i++ )
    {
      expon[i] = 0;
    }
    i = 4 % n;
    expon[i] = 8;
    en_r2_test ( n, expon );

    for ( i = 0; i < n; i++ )
    {
      expon[i] = 0;
    }
    i = 5 % n;
    expon[i] = 10;
    en_r2_test ( n, expon );

    for ( i = 0; i < n; i++ )
    {
      expon[i] = i + 1;
    }
    en_r2_test ( n, expon );

    for ( i = 0; i < n; i++ )
    {
      expon[i] = 2;
    }
    en_r2_test ( n, expon );

    free ( expon );
  }

  return;
}
/*****************************************************************************80*/

void en_r2_test ( int n, int expon[] )

/*****************************************************************************80

  Purpose:

    EN_R2_TEST tests the Stroud EN_R2 rules on a monomial.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 January 2010

  Author:

    John Burkardt
*/
{
  double c1;
  int d;
  double delta0;
  double err;
  double exact;
  double gamma0;
  int i;
  int o;
  int option;
  int p;
  double pi = 3.141592653589793;
  double quad;
  double *v;
  double volume_1d;
  double *w;
  double *x;

  printf ( "\n" );
  printf ( "  N = %d\n", n );
  printf ( "  EXPON = " );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %4d", expon[i] );
  }
  printf ( "\n" );
  d = i4vec_sum ( n, expon );
  printf ( "  Degree = %d\n", d );
  printf ( "\n" );

  exact = en_r2_monomial_integral ( n, expon );

  p = 1;

  if ( d <= p )
  {
    o =  en_r2_01_1_size ( n );
    x = ( double * ) malloc ( n * o * sizeof ( double ) );
    w = ( double * ) malloc ( o * sizeof ( double ) );
    en_r2_01_1 ( n, o, x, w );
    v = monomial_value ( n, o, x, expon );
    quad = r8vec_dot_product ( o, w, v );
    err = r8_abs ( quad - exact );
    printf ( "  EN_R2_01_1:      %6d  %14f  %14f\n", o, quad, err );
    free ( v );
    free ( w );
    free ( x );
  }

  p = 2;

  if ( d <= p )
  {
    o =  en_r2_02_xiu_size ( n );
    x = ( double * ) malloc ( n * o * sizeof ( double ) );
    w = ( double * ) malloc ( o * sizeof ( double ) );
    en_r2_02_xiu ( n, o, x, w );
    v = monomial_value ( n, o, x, expon );
    quad = r8vec_dot_product ( o, w, v );
    err = r8_abs ( quad - exact );
    printf ( "  EN_R2_02_XIU:    %6d  %14f  %14f\n", o, quad, err );
    free ( v );
    free ( w );
    free ( x );

    o = gw_02_xiu_size ( n );
    gamma0 = 2.0;
    delta0 = 0.0;
    c1 = 1.0;
    volume_1d = sqrt ( pi );
    x = ( double * ) malloc ( n * o * sizeof ( double ) );
    w = ( double * ) malloc ( o * sizeof ( double ) );
    gw_02_xiu ( n, o, gamma0, delta0, c1, volume_1d, x, w );
    v = monomial_value ( n, o, x, expon );
    quad = r8vec_dot_product ( o, w, v );
    err = r8_abs ( quad - exact );
    printf ( "  GW_02_XIU:       %6d  %14f  %14f\n", o, quad, err );
    free ( v );
    free ( w );
    free ( x );
  }

  p = 3;

  if ( d <= p )
  {
    o =  en_r2_03_1_size ( n );
    x = ( double * ) malloc ( n * o * sizeof ( double ) );
    w = ( double * ) malloc ( o * sizeof ( double ) );
    en_r2_03_1 ( n, o, x, w );
    v = monomial_value ( n, o, x, expon );
    quad = r8vec_dot_product ( o, w, v );
    err = r8_abs ( quad - exact );
    printf ( "  EN_R2_03_1:      %6d  %14f  %14f\n", o, quad, err );
    free ( v );
    free ( w );
    free ( x );

    o =  en_r2_03_2_size ( n );
    x = ( double * ) malloc ( n * o * sizeof ( double ) );
    w = ( double * ) malloc ( o * sizeof ( double ) );
    en_r2_03_2 ( n, o, x, w );
    v = monomial_value ( n, o, x, expon );
    quad = r8vec_dot_product ( o, w, v );
    err = r8_abs ( quad - exact );
    printf ( "  EN_R2_03_2:      %6d  %14f  %14f\n", o, quad, err );
    free ( v );
    free ( w );
    free ( x );

    o =  en_r2_03_xiu_size ( n );
    x = ( double * ) malloc ( n * o * sizeof ( double ) );
    w = ( double * ) malloc ( o * sizeof ( double ) );
    en_r2_03_xiu ( n, o, x, w );
    v = monomial_value ( n, o, x, expon );
    quad = r8vec_dot_product ( o, w, v );
    err = r8_abs ( quad - exact );
    printf ( "  EN_R2_03_XIU:    %6d  %14f  %14f\n", o, quad, err );
    free ( v );
    free ( w );
    free ( x );
  }

  p = 5;

  if ( d <= p )
  {
    if ( 2 <= n && n <= 7 )
    {
      option = 1;
      o =  en_r2_05_1_size ( n );
      x = ( double * ) malloc ( n * o * sizeof ( double ) );
      w = ( double * ) malloc ( o * sizeof ( double ) );
      en_r2_05_1 ( n, option, o, x, w );
      v = monomial_value ( n, o, x, expon );
      quad = r8vec_dot_product ( o, w, v );
      err = r8_abs ( quad - exact );
      printf ( "  EN_R2_05_1(1):   %6d  %14f  %14f\n", o, quad, err );
      free ( v );
      free ( w );
      free ( x );
    }

    if ( n == 3 || n == 5 || n == 6 )
    {
      option = 2;
      o =  en_r2_05_1_size ( n );
      x = ( double * ) malloc ( n * o * sizeof ( double ) );
      w = ( double * ) malloc ( o * sizeof ( double ) );
      en_r2_05_1 ( n, option, o, x, w );
      v = monomial_value ( n, o, x, expon );
      quad = r8vec_dot_product ( o, w, v );
      err = r8_abs ( quad - exact );
      printf ( "  EN_R2_05_1(2):   %6d  %14f  %14f\n", o, quad, err );
      free ( v );
      free ( w );
      free ( x );
    }

    o =  en_r2_05_2_size ( n );
    x = ( double * ) malloc ( n * o * sizeof ( double ) );
    w = ( double * ) malloc ( o * sizeof ( double ) );
    en_r2_05_2 ( n, o, x, w );
    v = monomial_value ( n, o, x, expon );
    quad = r8vec_dot_product ( o, w, v );
    err = r8_abs ( quad - exact );
    printf ( "  EN_R2_05_2:      %6d  %14f  %14f\n", o, quad, err );
    free ( v );
    free ( w );
    free ( x );

    if ( 3 <= n )
    {
      o =  en_r2_05_3_size ( n );
      x = ( double * ) malloc ( n * o * sizeof ( double ) );
      w = ( double * ) malloc ( o * sizeof ( double ) );
      en_r2_05_3 ( n, o, x, w );
      v = monomial_value ( n, o, x, expon );
      quad = r8vec_dot_product ( o, w, v );
      err = r8_abs ( quad - exact );
      printf ( "  EN_R2_05_3:      %6d  %14f  %14f\n", o, quad, err );
      free ( v );
      free ( w );
      free ( x );
    }

    o =  en_r2_05_4_size ( n );
    x = ( double * ) malloc ( n * o * sizeof ( double ) );
    w = ( double * ) malloc ( o * sizeof ( double ) );
    en_r2_05_4 ( n, o, x, w );
    v = monomial_value ( n, o, x, expon );
    quad = r8vec_dot_product ( o, w, v );
    err = r8_abs ( quad - exact );
    printf ( "  EN_R2_05_4:      %6d  %14f  %14f\n", o, quad, err );
    free ( v );
    free ( w );
    free ( x );

    o =  en_r2_05_5_size ( n );
    x = ( double * ) malloc ( n * o * sizeof ( double ) );
    w = ( double * ) malloc ( o * sizeof ( double ) );
    en_r2_05_5 ( n, o, x, w );
    v = monomial_value ( n, o, x, expon );
    quad = r8vec_dot_product ( o, w, v );
    err = r8_abs ( quad - exact );
    printf ( "  EN_R2_05_5:      %6d  %14f  %14f\n", o, quad, err );
    free ( v );
    free ( w );
    free ( x );

    if ( 5 <= n )
    {
      o =  en_r2_05_6_size ( n );
      x = ( double * ) malloc ( n * o * sizeof ( double ) );
      w = ( double * ) malloc ( o * sizeof ( double ) );
      en_r2_05_6 ( n, o, x, w );
      v = monomial_value ( n, o, x, expon );
      quad = r8vec_dot_product ( o, w, v );
      err = r8_abs ( quad - exact );
      printf ( "  EN_R2_05_6:      %6d  %14f  %14f\n", o, quad, err );
      free ( v );
      free ( w );
      free ( x );
    }
  }
  p = 7;

  if ( d <= p )
  {
    if ( n == 3 || n == 4 || n == 6 || n == 7 )
    {
      option = 1;
      o =  en_r2_07_1_size ( n );
      x = ( double * ) malloc ( n * o * sizeof ( double ) );
      w = ( double * ) malloc ( o * sizeof ( double ) );
      en_r2_07_1 ( n, option, o, x, w );
      v = monomial_value ( n, o, x, expon );
      quad = r8vec_dot_product ( o, w, v );
      err = r8_abs ( quad - exact );
      printf ( "  EN_R2_07_1(1):   %6d  %14f  %14f\n", o, quad, err );
      free ( v );
      free ( w );
      free ( x );
    }
    if ( n == 3 || n == 4 )
    {
      option = 2;
      o =  en_r2_07_1_size ( n );
      x = ( double * ) malloc ( n * o * sizeof ( double ) );
      w = ( double * ) malloc ( o * sizeof ( double ) );
      en_r2_07_1 ( n, option, o, x, w );
      v = monomial_value ( n, o, x, expon );
      quad = r8vec_dot_product ( o, w, v );
      err = r8_abs ( quad - exact );
      printf ( "  EN_R2_07_1(2):   %6d  %14f  %14f\n", o, quad, err );
      free ( v );
      free ( w );
      free ( x );
    }
    if ( 3 <= n )
    {
      o =  en_r2_07_2_size ( n );
      x = ( double * ) malloc ( n * o * sizeof ( double ) );
      w = ( double * ) malloc ( o * sizeof ( double ) );
      en_r2_07_2 ( n, o, x, w );
      v = monomial_value ( n, o, x, expon );
      quad = r8vec_dot_product ( o, w, v );
      err = r8_abs ( quad - exact );
      printf ( "  EN_R2_07_2:      %6d  %14f  %14f\n", o, quad, err );
      free ( v );
      free ( w );
      free ( x );
    }

    if ( 3 <= n && n <= 6 )
    {
      option = 1;
      o =  en_r2_07_3_size ( n );
      x = ( double * ) malloc ( n * o * sizeof ( double ) );
      w = ( double * ) malloc ( o * sizeof ( double ) );
      en_r2_07_3 ( n, option, o, x, w );
      v = monomial_value ( n, o, x, expon );
      quad = r8vec_dot_product ( o, w, v );
      err = r8_abs ( quad - exact );
      printf ( "  EN_R2_07_3(1):   %6d  %14f  %14f\n", o, quad, err );
      free ( v );
      free ( w );
      free ( x );
    }
    if ( n == 3 || n == 4 )
    {
      option = 2;
      o =  en_r2_07_3_size ( n );
      x = ( double * ) malloc ( n * o * sizeof ( double ) );
      w = ( double * ) malloc ( o * sizeof ( double ) );
      en_r2_07_3 ( n, option, o, x, w );
      v = monomial_value ( n, o, x, expon );
      quad = r8vec_dot_product ( o, w, v );
      err = r8_abs ( quad - exact );
      printf ( "  EN_R2_07_3(2):   %6d  %14f  %14f\n", o, quad, err );
      free ( v );
      free ( w );
      free ( x );
    }
  }

  p = 9;

  if ( d <= p )
  {
    if ( 3 <= n && n <= 6 )
    {
      option = 1;
      o =  en_r2_09_1_size ( n );
      x = ( double * ) malloc ( n * o * sizeof ( double ) );
      w = ( double * ) malloc ( o * sizeof ( double ) );
      en_r2_09_1 ( n, option, o, x, w );
      v = monomial_value ( n, o, x, expon );
      quad = r8vec_dot_product ( o, w, v );
      err = r8_abs ( quad - exact );
      printf ( "  EN_R2_09_1(1):   %6d  %14f  %14f\n", o, quad, err );
      free ( v );
      free ( w );
      free ( x );

      option = 2;
      o =  en_r2_09_1_size ( n );
      x = ( double * ) malloc ( n * o * sizeof ( double ) );
      w = ( double * ) malloc ( o * sizeof ( double ) );
      en_r2_09_1 ( n, option, o, x, w );
      v = monomial_value ( n, o, x, expon );
      quad = r8vec_dot_product ( o, w, v );
      err = r8_abs ( quad - exact );
      printf ( "  EN_R2_09_1(2):   %6d  %14f  %14f\n", o, quad, err );
      free ( v );
      free ( w );
      free ( x );
    }
  }

  p = 11;

  if ( d <= p )
  {
    if ( 3 <= n && n <= 5 )
    {
      option = 1;
      o =  en_r2_11_1_size ( n );
      x = ( double * ) malloc ( n * o * sizeof ( double ) );
      w = ( double * ) malloc ( o * sizeof ( double ) );
      en_r2_11_1 ( n, option, o, x, w );
      v = monomial_value ( n, o, x, expon );
      quad = r8vec_dot_product ( o, w, v );
      err = r8_abs ( quad - exact );
      printf ( "  EN_R2_11_1(1):   %6d  %14f  %14f\n", o, quad, err );
      free ( v );
      free ( w );
      free ( x );

      option = 2;
      o =  en_r2_11_1_size ( n );
      x = ( double * ) malloc ( n * o * sizeof ( double ) );
      w = ( double * ) malloc ( o * sizeof ( double ) );
      en_r2_11_1 ( n, option, o, x, w );
      v = monomial_value ( n, o, x, expon );
      quad = r8vec_dot_product ( o, w, v );
      err = r8_abs ( quad - exact );
      printf ( "  EN_R2_11_1(2):   %6d  %14f  %14f\n", o, quad, err );
      free ( v );
      free ( w );
      free ( x );
    }
  }
  printf ( "  EXACT                    %14f\n", exact );

  return;
}
/*****************************************************************************80*/

void test2075 ( )

/*****************************************************************************80

  Purpose:

    TEST2075 tests the rules for EPN with GLG on monomials.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    29 January 2010

  Author:

    John Burkardt
*/
{
# define TEST_NUM 5

  double alpha;
  double alpha_test[TEST_NUM] = { -0.5, 0.0, 0.5, 1.0, 2.0 };
  int *expon;
  int i;
  int n;
  int test;

  printf ( "\n" );
  printf ( "TEST2075\n" );
  printf ( "  Demonstrate the use of quadrature rules for the region\n" );
  printf ( "  EPN_GLG, that is, the positive half space [0,+oo)^N, with the\n" );
  printf ( "  weight W(ALPHA;X) = product ( 1 <= I <= N ) X(I)^ALPHA exp ( -X(I) )\n" );
  printf ( "\n" );
  printf ( "  We use the formulas to integrate various monomials of\n" );
  printf ( "  the form X(1)^E(1) * X(2)^E(2) * ... X(N)^E(N)\n" );
  printf ( "  and compare to the exact integral.\n" );
  printf ( "\n" );
  printf ( "  The precision of each formula is known, and we only use\n" );
  printf ( "  a formula if its precision indicates it should be able to\n" );
  printf ( "  produce an exact result.\n" );

  for ( n = 1; n <= 6; n++ )
  {
    expon = ( int * ) malloc ( n * sizeof ( int ) );

    i4vec_zero ( n, expon );
    for ( test = 0; test < TEST_NUM; test++ )
    {
      alpha = alpha_test[test];
      epn_glg_test ( n, expon, alpha );
    }

    i4vec_zero ( n, expon );
    expon[n-1] = 1;
    for ( test = 0; test < TEST_NUM; test++ )
    {
      alpha = alpha_test[test];
      epn_glg_test ( n, expon, alpha );
    }

    if ( 2 <= n )
    {
      i4vec_zero ( n, expon );
      expon[0] = 1;
      expon[1] = 1;
      for ( test = 0; test < TEST_NUM; test++ )
      {
        alpha = alpha_test[test];
        epn_glg_test ( n, expon, alpha );
      }
    }

    i4vec_zero ( n, expon );
    expon[0] = 2;
    for ( test = 0; test < TEST_NUM; test++ )
    {
      alpha = alpha_test[test];
      epn_glg_test ( n, expon, alpha );
    }
    free ( expon );
  }

  return;
# undef TEST_NUM
}
/*****************************************************************************80*/

void epn_glg_test ( int n, int expon[], double alpha )

/*****************************************************************************80

  Purpose:

    EPN_GLG_TEST tests the rules for EPN with GLG weight on a monomial.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    27 January 2010

  Author:

    John Burkardt
*/
{
  double c1;
  int d;
  double delta0;
  double err;
  double exact;
  double gamma0;
  int i;
  int o;
  int option;
  int p;
  double quad;
  double *v;
  double volume_1d;
  double *w;
  double *x;

  printf ( "\n" );
  printf ( "  N = %d\n", n );
  printf ( "  ALPHA = %f\n", alpha );
  printf ( "  EXPON = " );
  for ( i = 0; i < n; i++ )
  {
    printf ( "%4d", expon[i] );
  }
  printf ( "\n" );
  d = i4vec_sum ( n, expon );
  printf ( "  Degree = %d\n", d );
  printf ( "\n" );

  exact = epn_glg_monomial_integral ( n, expon, alpha );

  p = 0;

  if ( d <= p )
  {
    o =  epn_glg_00_1_size ( n, alpha );
    x = ( double * ) malloc ( n * o * sizeof ( double ) );
    w = ( double * ) malloc ( o * sizeof ( double ) );
    epn_glg_00_1 ( n, alpha, o, x, w );
    v = monomial_value ( n, o, x, expon );
    quad = r8vec_dot_product ( o, w, v );
    err = r8_abs ( quad - exact );
    printf ( "  EPN_GLG_00_1:     %6d  %14f  %14f\n", o, quad, err );
    free ( v );
    free ( w );
    free ( x );
  }

  p = 1;

  if ( d <= p )
  {
    o =  epn_glg_01_1_size ( n, alpha );
    x = ( double * ) malloc ( n * o * sizeof ( double ) );
    w = ( double * ) malloc ( o * sizeof ( double ) );
    epn_glg_01_1 ( n, alpha, o, x, w );
    v = monomial_value ( n, o, x, expon );
    quad = r8vec_dot_product ( o, w, v );
    err = r8_abs ( quad - exact );
    printf ( "  EPN_GLG_01_1:     %6d  %14f  %14f\n", o, quad, err );
    free ( v );
    free ( w );
    free ( x );
  }

  p = 2;

  if ( d <= p )
  {
    o =  epn_glg_02_xiu_size ( n, alpha );
    x = ( double * ) malloc ( n * o * sizeof ( double ) );
    w = ( double * ) malloc ( o * sizeof ( double ) );
    epn_glg_02_xiu ( n, alpha, o, x, w );
    v = monomial_value ( n, o, x, expon );
    quad = r8vec_dot_product ( o, w, v );
    err = r8_abs ( quad - exact );
    printf ( "  EPN_GLG_02_XIU:   %6d  %14f  %14f\n", o, quad, err );
    free ( v );
    free ( w );
    free ( x );

    o = gw_02_xiu_size ( n );
    gamma0 = - 1.0;
    delta0 = alpha + 1.0;
    c1 = - alpha - 1.0;
    volume_1d = r8_gamma ( 1.0 + alpha );
    x = ( double * ) malloc ( n * o * sizeof ( double ) );
    w = ( double * ) malloc ( o * sizeof ( double ) );
    gw_02_xiu ( n, o, gamma0, delta0, c1, volume_1d, x, w );
    v = monomial_value ( n, o, x, expon );
    quad = r8vec_dot_product ( o, w, v );
    err = r8_abs ( quad - exact );
    printf ( "  GW_02_XIU:        %6d  %14f  %14f\n", o, quad, err );
    free ( v );
    free ( w );
    free ( x );
  }
  printf ( "  EXACT                     %14f\n", exact );

  return;
}
/*****************************************************************************80*/

void test208 ( )

/*****************************************************************************80

  Purpose:

    TEST208 tests the rules for EPN with Laguerre weight on monomials.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    28 January 2010

  Author:

    John Burkardt
*/
{
  int *expon;
  int i;
  int n;
  int test;

  printf ( "\n" );
  printf ( "TEST208\n" );
  printf ( "  Demonstrate the use of quadrature rules for the region\n" );
  printf ( "  EPN_LAG, that is, the positive half space [0,+oo)^N, with the\n" );
  printf ( "  weight W(X) = product ( 1 <= I <= N ) exp ( -X(I) )\n" );
  printf ( "\n" );
  printf ( "  We use the formulas to integrate various monomials of\n" );
  printf ( "  the form X(1)^E(1) * X(2)^E(2) * ... X(N)^E(N)\n" );
  printf ( "  and compare to the exact integral.\n" );
  printf ( "\n" );
  printf ( "  The precision of each formula is known, and we only use\n" );
  printf ( "  a formula if its precision indicates it should be able to\n" );
  printf ( "  produce an exact result.\n" );

  for ( n = 1; n <= 6; n++ )
  {
    expon = ( int * ) malloc ( n * sizeof ( int ) );

    i4vec_zero ( n, expon );
    epn_lag_test ( n, expon );

    i4vec_zero ( n, expon );
    expon[n-1] = 1;
    epn_lag_test ( n, expon );

    if ( 2 <= n )
    {
      i4vec_zero ( n, expon );
      expon[0] = 1;
      expon[1] = 1;
      epn_lag_test ( n, expon );
    }

    i4vec_zero ( n, expon );
    expon[0] = 2;
    epn_lag_test ( n, expon );

    free ( expon );
  }

  return;
}
/*****************************************************************************80*/

void epn_lag_test ( int n, int expon[] )

/*****************************************************************************80

  Purpose:

    EPN_LAG_TEST tests the rules for EPN with Laguerre weight on a monomial.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    27 January 2010

  Author:

    John Burkardt
*/
{
  double c1;
  int d;
  double delta0;
  double err;
  double exact;
  double gamma0;
  int i;
  int o;
  int option;
  int p;
  double quad;
  double *v;
  double volume_1d;
  double *w;
  double *x;

  printf ( "\n" );
  printf ( "  N = %d\n", n );
  printf ( "  EXPON = " );
  for ( i = 0; i < n; i++ )
  {
    printf ( "%4d", expon[i] );
  }
  printf ( "\n" );
  d = i4vec_sum ( n, expon );
  printf ( "  Degree = %4d\n", d );
  printf ( "\n" );

  exact = epn_lag_monomial_integral ( n, expon );

  p = 0;

  if ( d <= p )
  {
    o =  epn_lag_00_1_size ( n );
    x = ( double * ) malloc ( n * o * sizeof ( double ) );
    w = ( double * ) malloc ( o * sizeof ( double ) );
    epn_lag_00_1 ( n, o, x, w );
    v = monomial_value ( n, o, x, expon );
    quad = r8vec_dot_product ( o, w, v );
    err = r8_abs ( quad - exact );
    printf ( "  EPN_LAG_00_1:     %6d  %14f  %14f\n", o, quad, err );
    free ( v );
    free ( w );
    free ( x );
  }

  p = 1;

  if ( d <= p )
  {
    o =  epn_lag_01_1_size ( n );
    x = ( double * ) malloc ( n * o * sizeof ( double ) );
    w = ( double * ) malloc ( o * sizeof ( double ) );
    epn_lag_01_1 ( n, o, x, w );
    v = monomial_value ( n, o, x, expon );
    quad = r8vec_dot_product ( o, w, v );
    err = r8_abs ( quad - exact );
    printf ( "  EPN_LAG_01_1:     %6d  %14f  %14f\n", o, quad, err );
    free ( v );
    free ( w );
    free ( x );
  }

  p = 2;

  if ( d <= p )
  {
    o =  epn_lag_02_xiu_size ( n );
    x = ( double * ) malloc ( n * o * sizeof ( double ) );
    w = ( double * ) malloc ( o * sizeof ( double ) );
    epn_lag_02_xiu ( n, o, x, w );
    v = monomial_value ( n, o, x, expon );
    quad = r8vec_dot_product ( o, w, v );
    err = r8_abs ( quad - exact );
    printf ( "  EPN_LAG_02_XIU:   %6d  %14f  %14f\n", o, quad, err );
    free ( v );
    free ( w );
    free ( x );

    o = gw_02_xiu_size ( n );
    gamma0 = - 1.0;
    delta0 = 1.0;
    c1 = - 1.0;
    volume_1d = 1.0;
    x = ( double * ) malloc ( n * o * sizeof ( double ) );
    w = ( double * ) malloc ( o * sizeof ( double ) );
    gw_02_xiu ( n, o, gamma0, delta0, c1, volume_1d, x, w );
    v = monomial_value ( n, o, x, expon );
    quad = r8vec_dot_product ( o, w, v );
    err = r8_abs ( quad - exact );
    printf ( "  GW_02_XIU:        %6d  %14f  %14f\n", o, quad, err );
    free ( v );
    free ( w );
    free ( x );
  }

  printf ( "  EXACT                     %14f\n", exact );

  return;
}
/*****************************************************************************80*/

void test21 ( )

/*****************************************************************************80

  Purpose:

   TEST21 tests HEXAGON_UNIT_SET and HEXAGON_SUM.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 April 2008

  Author:

    John Burkardt
*/
{
  double center[2] = { 0.0, 0.0 };
  int dim_num = 2;
  int i;
  int ihi;
  int ilo;
  char name[8];
  int num;
  int order;
  double rad;
  double result;
  int rule;
  int rule_max = 4;
  double *weight;
  double *xtab;
  double *ytab;

  rad = 2.0;

  printf ( "\n" );
  printf ( "TEST21\n" );
  printf ( "  HEXAGON_UNIT_SET sets a quadrature rule for the\n" );
  printf ( "    unit hexagon.\n" );
  printf ( "  HEXAGON_SUM evaluates the quadrature rule\n" );
  printf ( "    in an arbitrary hexagon.\n" );
  printf ( "\n" );
  printf ( "  We use a radius %f\n", rad );
  printf ( "  and center:\n" );
  printf ( "  CENTER = (%f, %f)\n", center[0], center[1] );
  printf ( "\n" );

  for ( ilo = 1; ilo <= rule_max; ilo = ilo + 5 )
  {
    ihi = i4_min ( ilo + 4, rule_max );

    printf ( "\n" );
    printf ( "  Rule:   " );
    for ( rule = ilo; rule <= ihi; rule++ )
    {
      printf ( "%6d", rule );
    }
    printf ( "\n" );
    printf ( "  Function\n" );
    printf ( "\n" );

    num = function_2d_num ( );

    for ( i = 1; i <= num; i++ )
    {
      function_2d_index = i;
      function_2d_name ( name );
      printf ( "  %8s", name );

      for ( rule = ilo; rule <= ihi; rule++ )
      {
        order = hexagon_unit_size ( rule );

        xtab = ( double * ) malloc ( order * sizeof ( double ) );
        ytab = ( double * ) malloc ( order * sizeof ( double ) );
        weight = ( double * ) malloc ( order * sizeof ( double ) );

        hexagon_unit_set ( rule, order, xtab, ytab, weight );

        result = hexagon_sum ( function_2d, center, rad, order, xtab, ytab,
          weight );

        printf ( "%14f", result );

        free ( xtab );
        free ( ytab );
        free ( weight );
      }
      printf ( "\n" );
    }
  }
  return;
}
/*****************************************************************************80*/

void test215 ( )

/*****************************************************************************80

  Purpose:

    TEST215 tests LENS_HALF_2D.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 April 2008

  Author:

    John Burkardt
*/
{
  double area;
  double center[2] = { 0.0, 0.0 };
  int dim_num = 2;
  int i;
  int i_max = 8;
  int order;
  double pi = 3.141592653589793;
  double r;
  double theta1;
  double theta2;
  double value;

  r = 1.0;

  printf ( "\n" );
  printf ( "TEST215\n" );
  printf ( "  LENS_HALF_2D approximates an integral within a\n" );
  printf ( "    circular half lens, defined by joining the endpoints\n" );
  printf ( "    of a circular arc.\n" );
  printf ( "\n" );
  printf ( "  Integrate F(X,Y) = 1\n" );
  printf ( "\n" );
  printf ( "      R            Theta1      Theta2        Area        Order Integral\n" );
  printf ( "\n" );

  for ( i = 0; i <= i_max; i++ )
  {
    theta1 = 0.0;
    theta2 = ( double ) ( i ) * 2.0 * pi / ( double ) ( i_max );

    area = lens_half_area_2d ( r, theta1, theta2 );

    printf ( "\n" );

    for ( order = 2; order <= 16; order = order + 2 )
    {
      value = lens_half_2d ( f_1_2d, center, r, theta1, theta2, order );
      printf ( "  %12f  %12f  %12f  %12f  %6d  %12f\n",
        r, theta1, theta2, area, order, value );
    }
  }

  printf ( "\n" );
  printf ( "\n" );
  printf ( "  Integrate F(X,Y) = X\n" );
  printf ( "\n" );
  printf ( "      R            Theta1      Theta2        Area        Order Integral\n" );
  printf ( "\n" );

  for ( i = 0; i <= i_max; i++ )
  {
    theta1 = 0.0;
    theta2 = ( double ) ( i ) * 2.0 * pi / ( double ) ( i_max );

    area = lens_half_area_2d ( r, theta1, theta2 );

    printf ( "\n" );

    for ( order = 2; order <= 16; order = order + 2 )
    {
      value = lens_half_2d ( f_x_2d, center, r, theta1, theta2, order );
      printf ( "  %12f  %12f  %12f  %12f  %6d  %12f\n",
        r, theta1, theta2, area, order, value );
    }
  }
  printf ( "\n" );
  printf ( "\n" );
  printf ( "  Integrate F(X,Y) = R\n" );
  printf ( "\n" );
  printf ( "      R            Theta1      Theta2        Area        Order Integral\n" );
  printf ( "\n" );

  for ( i = 0; i <= i_max; i++ )
  {
    theta1 = 0.0;
    theta2 = ( double ) ( i ) * 2.0 * pi / ( double ) ( i_max );

    area = lens_half_area_2d ( r, theta1, theta2 );

    printf ( "\n" );

    for ( order = 2; order <= 16; order = order + 2 )
    {
      value = lens_half_2d ( f_r_2d, center, r, theta1, theta2, order );
      printf ( "  %12f  %12f  %12f  %12f  %6d  %12f\n",
        r, theta1, theta2, area, order, value );
    }
  }

  return;
}
/*****************************************************************************80*/

void test22 ( )

/*****************************************************************************80

  Purpose:

    TEST22 tests OCTAHEDRON_UNIT_ND.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 April 2008

  Author:

    John Burkardt
*/
{
  int i;
  int n;
  int n_max = 3;
  char name[8];
  int num;
  double result;

  printf ( "\n" );
  printf ( "TEST22\n" );
  printf ( "  OCTAHEDRON_UNIT_ND approximates integrals in a unit\n" );
  printf ( "    octahedron in N dimensions.\n" );
  printf ( "\n" );
  printf ( "\n" );
  printf ( "    F(X)    N = 1    N = 2   N = 3\n" );
  printf ( "\n" );

  num = function_nd_num ( );

  for ( i = 1; i <= num; i++ )
  {
    function_nd_index = i;
    function_nd_name ( name );
    printf ( "  %8s", name );

    for ( n = 1; n <= n_max; n++ )
    {
      result = octahedron_unit_nd ( function_nd, n );
      printf ( "%14f", result );
    }
    printf ( "\n" );
  }
  return;
}
/*****************************************************************************80*/

void test23 ( )

/*****************************************************************************80

  Purpose:

    TEST23 tests PARALLELIPIPED_VOLUME_ND.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 April 2008

  Author:

    John Burkardt
*/
{
  int i;
  int j;
  int n;
  double *v;
  double volume;

  printf ( "\n" );
  printf ( "TEST23\n" );
  printf ( "  PARALLELIPIPED_VOLUME_ND computes the volume of a\n" );
  printf ( "    parallelipiped in N dimensions.\n" );
  printf ( "\n" );

  for ( n = 2; n <= 4; n++ )
  {
    printf ( "\n" );
    printf ( "  Spatial dimension N = %d\n", n );
/*
  Set the values of the parallelipiped.
*/
    v = setsim ( n );

    printf ( "\n" );
    printf ( "  Parallelipiped vertices:\n" );
    printf ( "\n" );

    for ( i = 0; i < n; i++ )
    {
      for ( j = 0; j < n + 1; j++ )
      {
        printf ( "%14f", v[i+j*n] );
      }
      printf ( "\n" );
    }
    volume = parallelipiped_volume_nd ( n, v );

    printf ( "\n" );
    printf ( "  Volume is %f\n", volume );

    free ( v );
  }
  return;
}
/*****************************************************************************80*/

void test24 ( )

/*****************************************************************************80

  Purpose:

    TEST24 tests POLYGON_**_2D;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 April 2008

  Author:

    John Burkardt
*/
{
  double result;
  int npts = 4;
  double x[4] = { 0.0, 1.0, 1.0, 0.0 };
  double y[4] = { 0.0, 0.0, 1.0, 1.0 };

  printf ( "\n" );
  printf ( "TEST24\n" );
  printf ( "  For a polygon in 2D:\n" );
  printf ( "  POLYGON_1_2D integrates 1\n" );
  printf ( "  POLYGON_X_2D integrates X\n" );
  printf ( "  POLYGON_Y_2D integrates Y\n" );
  printf ( "  POLYGON_XX_2D integrates X*X\n" );
  printf ( "  POLYGON_XY_2D integrates X*Y\n" );
  printf ( "  POLYGON_YY_2D integrates Y*Y\n" );
  printf ( "\n" );
  printf ( "  F(X,Y)    Integral\n" );
  printf ( "\n" );

  result = polygon_1_2d ( npts, x, y );
  printf ( "     1    %f\n", result );

  result = polygon_x_2d ( npts, x, y );
  printf ( "     X    %f\n", result );

  result = polygon_y_2d ( npts, x, y );
  printf ( "     Y    %f\n", result );

  result = polygon_xx_2d ( npts, x, y );
  printf ( "   X*X    %f\n", result );

  result = polygon_xy_2d ( npts, x, y );
  printf ( "   X*Y    %f\n", result );

  result = polygon_yy_2d ( npts, x, y );
  printf ( "   Y*Y    %f\n", result );

  return;
}
/*****************************************************************************80*/

void test25 ( )

/*****************************************************************************80

  Purpose:

    TEST25 tests PYRAMID_UNIT_O**_3D, PYRAMID_VOLUME_3D.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 April 2008

  Author:

    John Burkardt
*/
{
  int i;
  int j;
  int jhi;
  int jlo;
  char name[8];
  int num;
  int order[10] = { 1, 5, 6, 8, 8, 9, 13, 18, 27, 48 };
  double result;

  printf ( "\n" );
  printf ( "TEST25\n" );
  printf ( "  For the unit pyramid, we approximate integrals with:\n" );
  printf ( "  PYRAMID_UNIT_O01_3D, a 1 point rule.\n" );
  printf ( "  PYRAMID_UNIT_O05_3D, a 5 point rule.\n" );
  printf ( "  PYRAMID_UNIT_O06_3D, a 6 point rule.\n" );
  printf ( "  PYRAMID_UNIT_O08_3D, an 8 point rule.\n" );
  printf ( "  PYRAMID_UNIT_O08b_3D, an 8 point rule.\n" );
  printf ( "  PYRAMID_UNIT_O09_3D, a 9 point rule.\n" );
  printf ( "  PYRAMID_UNIT_O13_3D, a 13 point rule.\n" );
  printf ( "  PYRAMID_UNIT_O18_3D, a 18 point rule.\n" );
  printf ( "  PYRAMID_UNIT_O27_3D, a 27 point rule.\n" );
  printf ( "  PYRAMID_UNIT_O48_3D, a 48 point rule.\n" );
  printf ( "\n" );
  printf ( "  PYRAMID_UNIT_VOLUME_3D computes the volume of a unit pyramid.\n" );
  printf ( "\n" );
  printf ( "  Volume = %f\n", pyramid_unit_volume_3d ( ) );
  printf ( "\n" );

  num = function_3d_num ( );

  for ( jlo = 1; jlo <= num; jlo = jlo + 5 )
  {
    jhi = i4_min ( jlo + 4, num );
    printf ( "\n" );
    printf ( "  Order   " );
    for ( j = jlo; j <= jhi; j++ )
    {
      function_3d_index = j;
      function_3d_name ( name );
      printf ( "   %8s    ", name );
    }
    printf ( "\n" );
    printf ( "\n" );
    printf ( "  %5d", order[0] );
    for ( j = jlo; j <= jhi; j++ )
    {
      function_3d_index = j;
      printf ( "%14f", pyramid_unit_o01_3d ( function_3d ) );
    }
    printf ( "\n" );
    printf ( "  %5d", order[1] );
    for ( j = jlo; j <= jhi; j++ )
    {
      function_3d_index = j;
      printf ( "%14f", pyramid_unit_o05_3d ( function_3d ) );
    }
    printf ( "\n" );
    printf ( "  %5d", order[2] );
    for ( j = jlo; j <= jhi; j++ )
    {
      function_3d_index = j;
      printf ( "%14f", pyramid_unit_o06_3d ( function_3d ) );
    }
    printf ( "\n" );
    printf ( "  %5d", order[3] );
    for ( j = jlo; j <= jhi; j++ )
    {
      function_3d_index = j;
      printf ( "%14f", pyramid_unit_o08_3d ( function_3d ) );
    }
    printf ( "\n" );
    printf ( "  %5d", order[4] );
    for ( j = jlo; j <= jhi; j++ )
    {
      function_3d_index = j;
      printf ( "%14f", pyramid_unit_o08b_3d ( function_3d ) );
    }
    printf ( "\n" );
    printf ( "  %5d", order[5] );
    for ( j = jlo; j <= jhi; j++ )
    {
      function_3d_index = j;
      printf ( "%14f", pyramid_unit_o09_3d ( function_3d ) );
    }
    printf ( "\n" );
    printf ( "  %5d", order[6] );
    for ( j = jlo; j <= jhi; j++ )
    {
      function_3d_index = j;
      printf ( "%14f", pyramid_unit_o13_3d ( function_3d ) );
    }
    printf ( "\n" );
    printf ( "  %5d", order[7] );
    for ( j = jlo; j <= jhi; j++ )
    {
      function_3d_index = j;
      printf ( "%14f", pyramid_unit_o18_3d ( function_3d ) );
    }
    printf ( "\n" );
    printf ( "  %5d", order[8] );
    for ( j = jlo; j <= jhi; j++ )
    {
      function_3d_index = j;
      printf ( "%14f", pyramid_unit_o27_3d ( function_3d ) );
    }
    printf ( "\n" );
    printf ( "  %5d", order[9] );
    for ( j = jlo; j <= jhi; j++ )
    {
      function_3d_index = j;
      printf ( "%14f", pyramid_unit_o48_3d ( function_3d ) );
    }
    printf ( "\n" );
  }
  return;
}
/*****************************************************************************80*/

void test255 ( )

/*****************************************************************************80

  Purpose:

    TEST255 tests PYRAMID_UNIT_MONOMIAL_3D.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 March 2008

  Author:

    John Burkardt
*/
{
  int alpha;
  int beta;
  int degree_max = 4;
  int gamma;
  double value;

  printf ( "\n" );
  printf ( "TEST255\n" );
  printf ( "  For the unit pyramid,\n" );
  printf ( "  PYRAMID_UNIT_MONOMIAL_3D returns the exact value of the\n" );
  printf ( " integral of X^ALPHA Y^BETA Z^GAMMA\n" );
  printf ( "\n" );
  printf ( "  Volume = %f\n", pyramid_unit_volume_3d ( ) );
  printf ( "\n" );
  printf ( "     ALPHA      BETA     GAMMA      INTEGRAL\n" );
  printf ( "\n" );

  for ( alpha = 0; alpha <= degree_max; alpha++ )
  {
    for ( beta = 0; beta <= degree_max - alpha; beta++ )
    {
      for ( gamma = 0; gamma <= degree_max - alpha - beta; gamma++ )
      {
        value = pyramid_unit_monomial_3d ( alpha, beta, gamma );

        printf ( "  %8d  %8d  %8d  %14f\n", alpha, beta, gamma, value );
      }
    }
  }
  return;
}
/*****************************************************************************80*/

void test26 ( )

/*****************************************************************************80

  Purpose:

    TEST26 tests QMULT_1D.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 April 2008

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  int i;
  char name[8];
  int num;
  double result;

  a = -1.0;
  b = 1.0;

  printf ( "\n" );
  printf ( "TEST26\n" );
  printf ( "  QMULT_1D approximates an integral on a\n" );
  printf ( "    one-dimensional interval.\n" );
  printf ( "\n" );
  printf ( "  We use the interval:\n" );
  printf ( "  A = %f\n", a );
  printf ( "  B = %f\n", b );
  printf ( "\n" );
  printf ( "    F(X)     QMULT_1D\n" );
  printf ( "\n" );

  num = function_1d_num ( );

  for ( i = 1; i <= num; i++ )
  {
    function_1d_index = i;
    function_1d_name ( name );

    result = qmult_1d ( function_1d, a, b );
    printf ( "  %8s  %14f\n", name, result );
  }

  return;
}
/*****************************************************************************80*/

void test27 ( )

/*****************************************************************************80

  Purpose:

    TEST27 tests SIMPLEX_ND.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    23 April 2013

  Author:

    John Burkardt
*/
{
  int i;
  int j;
  int n;
  char name[8];
  int num;
  double result;
  double *v;

  printf ( "\n" );
  printf ( "TEST27\n" );
  printf ( "  SIMPLEX_ND approximates integrals inside an\n" );
  printf ( "  arbitrary simplex in ND.\n" );
  printf ( "\n" );

  for ( n = 2; n <= 4; n++ )
  {
    printf ( "\n" );
    printf ( "  Spatial dimension N = %d\n", n );
/*
  Restore values of simplex.
*/
    v = setsim ( n );

    printf ( "\n" );
    printf ( "  Simplex vertices:\n" );
    printf ( "\n" );

    for ( i = 0; i < n; i++ )
    {
      for ( j = 0; j < n + 1; j++ )
      {
        printf ( "%4f", v[i+j*n] );
      }
      printf ( "\n" );
    }

    free ( v );

    printf ( "\n" );
    printf ( "  F(X)    SIMPLEX_ND\n" );
    printf ( "\n" );

    num = function_nd_num ( );

    for ( i = 1; i <= num; i++ )
    {
      v = setsim ( n );
      function_nd_index = i;
      function_nd_name ( name );

      result = simplex_nd ( function_nd, n, v );
      printf ( "  %8s  %14f\n", name, result );
      free ( v );
    }
  }
  return;
}
/*****************************************************************************80*/

void test28 ( )

/*****************************************************************************80

  Purpose:

    TEST28 tests SIMPLEX_VOLUME_ND.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    06 April 2008

  Author:

    John Burkardt
*/
{
  int i;
  int j;
  int n;
  double *v;
  double volume;

  printf ( "\n" );
  printf ( "TEST28\n" );
  printf ( "  SIMPLEX_VOLUME_ND computes the volume of a simplex\n" );
  printf ( "    in N dimensions.\n" );
  printf ( "\n" );

  for ( n = 2; n <= 4; n++ )
  {
    printf ( "\n" );
    printf ( "  Spatial dimension N = %d\n", n );
/*
  Set the values of the simplex.
*/
    v = setsim ( n );

    printf ( "\n" );
    printf ( "  Simplex vertices:\n" );
    printf ( "\n" );

    for ( i = 0; i < n; i++ )
    {
      for ( j = 0; j < n + 1; j++ )
      {
        printf ( "%14f", v[i+j*n] );
      }
      printf ( "\n" );
    }
    volume = simplex_volume_nd ( n, v );

    printf ( "\n" );
    printf ( "  Volume is %f\n", volume );

    free ( v );
  }
  return;
}
/*****************************************************************************80*/

void test29 ( )

/*****************************************************************************80

  Purpose:

    TEST29 tests SIMPLEX_UNIT_**_ND.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 April 2008

  Author:

    John Burkardt
*/
{
  int i;
  int n;
  char name[8];
  double result1;
  double result2;
  double result3;
  double result4;
  double volume;

  printf ( "\n" );
  printf ( "TEST29\n" );
  printf ( "  For integrals in the unit simplex in ND,\n" );
  printf ( "  SIMPLEX_UNIT_01_ND uses a formula of degree 1.\n" );
  printf ( "  SIMPLEX_UNIT_03_ND uses a formula of degree 3.\n" );
  printf ( "  SIMPLEX_UNIT_05_ND uses a formula of degree 5.\n" );
  printf ( "  SIMPLEX_UNIT_05_2_ND uses a formula of degree 5.\n" );

  for ( i = 1; i <= 6; i++ )
  {
    function_nd_index = i;
    function_nd_name ( name );

    printf ( "\n" );
    printf ( "  Check the integral of %s\n", name );
    printf ( "\n" );
    printf ( "  N     Volume         #1              #3              #5              #5.2\n" );
    printf ( "\n" );

    for ( n = 2; n <= 16; n++ )
    {
      result1 = simplex_unit_01_nd ( function_nd, n );
      result2 = simplex_unit_03_nd ( function_nd, n );
      result3 = simplex_unit_05_nd ( function_nd, n );
      result4 = simplex_unit_05_2_nd ( function_nd, n );

      volume = simplex_unit_volume_nd ( n );

      printf ( "  %2d  %13f  %13f  %13f  %13f  %13f\n",
        n, volume, result1, result2, result3, result4 );
    }
  }
  return;
}
/*****************************************************************************80*/

void test30 ( )

/*****************************************************************************80

  Purpose:

    TEST30 tests SPHERE_UNIT_**_3D.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 April 2008

  Author:

    John Burkardt
*/
{
  int i;
  char name[8];
  int num;
  double result1;
  double result2;
  double result3;
  double result4;

  printf ( "\n" );
  printf ( "TEST30\n" );
  printf ( "  For integrals on the unit sphere in 3D:\n" );
  printf ( "  SPHERE_UNIT_07_3D uses a formula of degree 7.\n" );
  printf ( "  SPHERE_UNIT_11_3D uses a formula of degree 11.\n" );
  printf ( "  SPHERE_UNIT_14_3D uses a formula of degree 14.\n" );
  printf ( "  SPHERE_UNIT_15_3D uses a formula of degree 15.\n" );
  printf ( "\n" );
  printf ( "  Unit sphere area = %f\n", sphere_unit_area_nd ( 3 ) );
  printf ( "\n" );
  printf ( "    F(X)    S3S07        S3S11         S3S14          S3S15\n" );
  printf ( "\n" );

  num = function_3d_num ( );

  for ( i = 1; i <= num; i++ )
  {
    function_3d_index = i;
    function_3d_name ( name );

    result1 = sphere_unit_07_3d ( function_3d );
    result2 = sphere_unit_11_3d ( function_3d );
    result3 = sphere_unit_14_3d ( function_3d );
    result4 = sphere_unit_15_3d ( function_3d );

    printf ( "  %8s  %14f  %14f  %14f  %14f\n",
      name, result1, result2, result3, result4 );
  }
  return;
}
/******************************************************************************/

void test31 ( )

/******************************************************************************/
/*
  Purpose:

    TEST31 tests SPHERE_UNIT_**_ND.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 April 2008

  Author:

    John Burkardt
*/
{
  int i;
  int n;
  char name[8];
  int num;
  double result1;
  double result2;
  double result3;
  double result4;
  double result5;
  double result6;

  printf ( "\n" );
  printf ( "TEST31\n" );
  printf ( "  For integrals on the unit sphere in ND:\n" );
  printf ( "  SPHERE_UNIT_03_ND uses a formula of degree 3;\n" );
  printf ( "  SPHERE_UNIT_04_ND uses a formula of degree 4;\n" );
  printf ( "  SPHERE_UNIT_05_ND uses a formula of degree 5.\n" );
  printf ( "  SPHERE_UNIT_07_1_ND uses a formula of degree 7.\n" );
  printf ( "  SPHERE_UNIT_07_2_ND uses a formula of degree 7.\n" );
  printf ( "  SPHERE_UNIT_11_ND uses a formula of degree 11.\n" );
  printf ( "\n" );

  for ( n = 3; n <= 10; n++ )
  {
    printf ( "\n" );
    printf ( "  Spatial dimension N = %d\n", n );
    printf ( "  Unit sphere area = %f\n", sphere_unit_area_nd ( n ) );
    printf ( "\n" );
    printf ( "    Rule:     #3            #4            #5\n" );
    printf ( "              #7.1          #7.2          #11\n" );
    printf ( "    Function\n" );
    printf ( "\n" );

    num = function_nd_num ( );

    for ( i = 1; i <= num; i++ )
    {
      function_nd_index = i;
      function_nd_name ( name );

      result1 = sphere_unit_03_nd ( function_nd, n );
      result2 = sphere_unit_04_nd ( function_nd, n );
      result3 = sphere_unit_05_nd ( function_nd, n );
      result4 = sphere_unit_07_1_nd ( function_nd, n );
      result5 = sphere_unit_07_2_nd ( function_nd, n );
      result6 = sphere_unit_11_nd ( function_nd, n );

      printf ( "  %8s  %14f  %14f  %14f\n",
        name, result1, result2, result3 );
      printf ( "            %14f  %14f  %14f\n", result4, result5, result6 );
    }
  }
  return;
}
/*****************************************************************************80*/

void test32 ( )

/*****************************************************************************80

  Purpose:

    TEST32 tests SPHERE_05_ND, SPHERE_07_1_ND.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 April 2008

  Author:

    John Burkardt
*/
{
  double *center;
  int dim;
  int i;
  int n;
  char name[8];
  int num;
  double r;
  double result1;
  double result2;

  printf ( "\n" );
  printf ( "TEST32\n" );
  printf ( "  For integrals on a sphere in ND:\n" );
  printf ( "  SPHERE_05_ND uses a formula of degree 5.\n" );
  printf ( "  SPHERE_07_1_ND uses a formula of degree 7.\n" );
  printf ( "\n" );
  r = 2.0;

  for ( n = 2; n <= 4; n++ )
  {
    center = ( double * ) malloc ( n * sizeof ( double ) );
    for ( dim = 0; dim < n; dim++ )
    {
      center[dim] = 1.0;
    }
    printf ( "\n" );
    printf ( "  Spatial dimension N = %d\n", n );
    printf ( "  Sphere center = \n" );
    for ( dim = 0; dim < n; dim++ )
    {
      printf ( "%14f", center[dim] );
    }
    printf ( "\n" );
    printf ( "  Sphere radius = %f\n", r );
    printf ( "  Sphere area = %f\n", sphere_area_nd ( n, r ) );
    printf ( "\n" );
    printf ( "    Rule:     #5           #7.1\n" );
    printf ( "    Function\n" );
    printf ( "\n" );

    num = function_nd_num ( );

    for ( i = 1; i <= num; i++ )
    {
      function_nd_index = i;
      function_nd_name ( name );

      result1 = sphere_05_nd ( function_nd, n, center, r );
      result2 = sphere_07_1_nd ( function_nd, n, center, r );

      printf ( "  %8s  %14f  %14f\n", name, result1, result2 );
    }
    free ( center );
  }
  return;
}
/*****************************************************************************80*/

void test322 ( )

/*****************************************************************************80

  Purpose:

    TEST322 tests SPHERE_CAP_AREA_3D, SPHERE_CAP_AREA_ND.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 April 2008

  Author:

    John Burkardt
*/
{
  double area1;
  double area2;
  double center[3] = { 0.0, 0.0, 0.0 };
  int dim_num = 3;
  double h;
  int i;
  int ntest = 12;
  double r = 1.0;

  printf ( "\n" );
  printf ( "TEST322\n" );
  printf ( "  SPHERE_CAP_AREA_3D computes the volume of a\n" );
  printf ( "    3D spherical cap, defined by a plane that cuts the\n" );
  printf ( "    sphere to a thickness of H units.\n" );
  printf ( "  SPHERE_CAP_AREA_ND computes the volume of an\n" );
  printf ( "    ND spherical cap, defined by a plane that cuts the\n" );
  printf ( "    sphere to a thickness of H units.\n" );

  area1 = sphere_area_3d ( r );

  printf ( "\n" );
  printf ( "  Area of the total sphere in 3D = %f\n", area1 );

  printf ( "\n" );
  printf ( "        R           H           Cap         Cap\n" );
  printf ( "                                area_3d     area_nd\n" );
  printf ( "\n" );

  for ( i= 0; i <= ntest + 1; i++ )
  {
    h = 2.0 * r * ( double ) ( i ) / ( double ) ( ntest );

    area1 = sphere_cap_area_3d ( r, h );

    area2 = sphere_cap_area_nd ( dim_num, r, h );

    printf ( "  %12f  %12f  %12f  %12f\n", r, h, area1, area2 );
  }
  return;
}
/*****************************************************************************80*/

void test324 ( )

/*****************************************************************************80

  Purpose:

    TEST324 tests SPHERE_CAP_VOLUME_2D, SPHERE_CAP_VOLUME_ND.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 April 2008

  Author:

    John Burkardt
*/
{
  double center[2] = { 0.0, 0.0 };
  int dim_num = 2;
  double h;
  int i;
  int ntest = 12;
  double pi = 3.141592653589793;
  double r = 1.0;
  double volume1;
  double volume2;

  printf ( "\n" );
  printf ( "TEST324\n" );
  printf ( "  SPHERE_CAP_VOLUME_2D computes the volume (area) of a\n" );
  printf ( "    spherical cap, defined by a plane that cuts the\n" );
  printf ( "    sphere to a thickness of H units.\n" );
  printf ( "  SPHERE_CAP_VOLUME_ND does the same operation,\n" );
  printf ( "    but in N dimensions.\n" );
  printf ( "\n" );
  printf ( "  Using a radius R = %f\n", r );

  volume1 = sphere_volume_2d ( r );

  printf ( "\n" );
  printf ( "  Volume of the total sphere in 2D = %f\n", volume1 );

  printf ( "\n" );
  printf ( "        H           Cap        Cap\n" );
  printf ( "                    vol_2d     vol_nd\n" );
  printf ( "\n" );

  for ( i = 0; i <= ntest + 1; i++ )
  {
    h = 2.0 * r * ( double ) ( i ) / ( double ) ( ntest );

    volume1 = sphere_cap_volume_2d ( r, h );

    volume2 = sphere_cap_volume_nd ( dim_num, r, h );

    printf ( "  %12f  %12f  %12f\n", h, volume1, volume2 );
  }
  return;
}
/*****************************************************************************80*/

void test326 ( )

/*****************************************************************************80

  Purpose:

    TEST326 tests SPHERE_CAP_VOLUME_3D, SPHERE_CAP_VOLUME_ND.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 April 2008

  Author:

    John Burkardt
*/
{
  double center[3] = { 0.0, 0.0, 0.0 };
  int dim_num = 3;
  double h;
  int i;
  int ntest = 12;
  double r = 1.0;
  double volume1;
  double volume2;

  printf ( "\n" );
  printf ( "TEST326\n" );
  printf ( "  SPHERE_CAP_VOLUME_3D computes the volume of a\n" );
  printf ( "    spherical cap, defined by a plane that cuts the\n" );
  printf ( "    sphere to a thickness of H units.\n" );
  printf ( "  SPHERE_CAP_VOLUME_ND does the same operation,\n" );
  printf ( "    but in N dimensions.\n" );
  printf ( "\n" );
  printf ( "  Using a radius R = %f\n", r );

  volume1 = sphere_volume_3d ( r );

  printf ( "\n" );
  printf ( "  Volume of the total sphere in 3D = %f\n", volume1 );

  printf ( "\n" );
  printf ( "        H           Cap        Cap\n" );
  printf ( "                    volume_3d  volume_nd\n" );
  printf ( "\n" );

  for ( i = 0; i <= ntest + 1; i++ )
  {
    h = 2.0 * r * ( double ) ( i ) / ( double ) ( ntest );

    volume1 = sphere_cap_volume_3d ( r, h );

    volume2 = sphere_cap_volume_nd ( dim_num, r, h );

    printf ( "  %12f  %12f  %12f\n", h, volume1, volume2 );
  }

  return;
}
/*****************************************************************************80*/

void test33 ( )

/*****************************************************************************80

  Purpose:

    TEST33 tests SPHERE_CAP_AREA_ND, SPHERE_CAP_VOLUME_ND.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    06 April 2008

  Author:

    John Burkardt
*/
{
  double area;
  double h;
  int i;
  int n = 12;
  int dim_num;
  double pi = 3.141592653589793;
  double r;
  double volume;

  printf ( "\n" );
  printf ( "TEST33\n" );
  printf ( "  For a sphere in ND:\n" );
  printf ( "  SPHERE_CAP_AREA_ND computes the area\n" );
  printf ( "    of a spherical cap.\n" );
  printf ( "  SPHERE_CAP_VOLUME_ND computes the volume\n" );
  printf ( "    of a spherical cap.\n" );
  printf ( "\n" );

  r = 1.0;

  for ( dim_num = 2; dim_num <= 5; dim_num++ )
  {
    printf ( "\n" );
    printf ( "  Spatial dimension N = %d\n", dim_num );
    printf ( "  Radius =       %f\n", r );
    printf ( "  Area =         %f\n", sphere_area_nd ( dim_num, r ) );
    volume = sphere_volume_nd ( dim_num, r );
    printf ( "  Volume =       %f\n", volume );

    printf ( "\n" );
    printf ( "                 Sphere         Sphere\n" );
    printf ( "                 cap            cap\n" );
    printf ( "      H          area           volume\n" );
    printf ( "\n" );

    for ( i = 0; i <= n + 1; i++ )
    {
      h = ( double ) ( 2 * i ) * r / ( double ) ( n );
      area = sphere_cap_area_nd ( dim_num, r, h );
      volume = sphere_cap_volume_nd ( dim_num, r, h );
      printf ( "  %8f  %14f  %14f\n", h, area, volume );
    }
  }

  return;
}
/*****************************************************************************80*/

void test335 ( )

/*****************************************************************************80

  Purpose:

    TEST335 tests SPHERE_SHELL_03_ND.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 April 2008

  Author:

    John Burkardt
*/
{
  double center[3];
  int i;
  int j;
  int n;
  int n_max = 3;
  char name[8];
  int num;
  double result1;
  double result3;
  double result4;
  double result5;
  double result6;
  double r1;
  double r2;

  printf ( "\n" );
  printf ( "TEST335\n" );
  printf ( "  For integrals inside a spherical shell in ND:\n" );
  printf ( "  SPHERE_SHELL_03_ND approximates the integral.\n" );
  printf ( "\n" );
  printf ( "  We compare these results with those computed by\n" );
  printf ( "  from the difference of two ball integrals:\n" );
  printf ( "\n" );
  printf ( "  BALL_F1_ND approximates the integral;\n" );
  printf ( "  BALL_F3_ND approximates the integral\n" );
  printf ( "\n" );

  for ( j = 1; j <= 2; j++ )
  {
    if ( j == 1 )
    {
      r1 = 0.0;
      r2 = 1.0;
      center[0] = 0.0;
      center[1] = 0.0;
      center[2] = 0.0;
    }
    else
    {
      r1 = 2.0;
      r2 = 3.0;
      center[0] =  1.0;
      center[1] = -1.0;
      center[2] =  2.0;
    }

    for ( n = 2; n <= n_max; n++ )
    {
      printf ( "\n" );
      printf ( "  Spatial dimension N = %d\n", n );
      printf ( "  Sphere center:\n" );
      for ( i = 0; i < n; i++ )
      {
        printf ( "%10f", center[i] );
      }
      printf ( "  Inner sphere radius = %f\n", r1 );
      printf ( "  Outer sphere radius = %f\n", r2 );
      printf ( "  Spherical shell volume = %f\n",
        sphere_shell_volume_nd ( n, r1, r2 ) );
      printf ( "\n" );
      printf ( "\n" );
      printf ( "    Rule:      #3       F1(R2)-F1(R1)  F3(R2)-F3(R1)\n" );
      printf ( "    F(X)\n" );
      printf ( "\n" );

      num = function_nd_num ( );

      for ( i = 1; i <= num; i++ )
      {
        function_nd_index = i;
        function_nd_name ( name );
        printf ( "  %8s", name );

        result1 = sphere_shell_03_nd ( function_nd, n, center, r1, r2 );

        result3 = ball_f1_nd ( function_nd, n, center, r1 );
        result4 = ball_f1_nd ( function_nd, n, center, r2 );

        result5 = ball_f3_nd ( function_nd, n, center, r1 );
        result6 = ball_f3_nd ( function_nd, n, center, r2 );

        printf ( "%14f  %14f  %14f\n",
          result1, result4 - result3, result6 - result5 );
      }
    }
  }
  return;
}
/*****************************************************************************80*/

void test34 ( )

/*****************************************************************************80

  Purpose:

    TEST34 tests SPHERE_UNIT_AREA_ND, SPHERE_UNIT_AREA_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 April 2008

  Author:

    John Burkardt
*/
{
  double area;
  double area2;
  int dim_num;
  int n_data;

  printf ( "\n" );
  printf ( "TEST34:\n" );
  printf ( "  SPHERE_UNIT_AREA_ND evaluates the area of the unit\n" );
  printf ( "  sphere in N dimensions.\n" );
  printf ( "  SPHERE_UNIT_AREA_VALUES returns some test values.\n" );
  printf ( "\n" );
  printf ( "     dim_num    Exact          Computed\n" );
  printf ( "                Area           Area\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    sphere_unit_area_values ( &n_data, &dim_num, &area );

    if ( n_data == 0 )
    {
      break;
    }
    area2 = sphere_unit_area_nd ( dim_num );

    printf ( "  %8d  %10f  %10f\n", dim_num, area, area2 );
  }
  return;
}
/*****************************************************************************80*/

void test345 ( )

/*****************************************************************************80

  Purpose:

    TEST345 tests SPHERE_UNIT_VOLUME_ND, SPHERE_UNIT_VOLUME_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    06 April 2008

  Author:

    John Burkardt
*/
{
  int dim_num;
  int n_data;
  double volume;
  double volume2;

  printf ( "\n" );
  printf ( "TEST345:\n" );
  printf ( "  SPHERE_UNIT_VOLUME_ND evaluates the area of the unit\n" );
  printf ( "  sphere in N dimensions.\n" );
  printf ( "  SPHERE_UNIT_VOLUME_VALUES returns some test values.\n" );
  printf ( "\n" );
  printf ( "     dim_num    Exact          Computed\n" );
  printf ( "                Volume         Volume\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    sphere_unit_volume_values ( &n_data, &dim_num, &volume );

    if ( n_data == 0 )
    {
      break;
    }

    volume2 = sphere_unit_volume_nd ( dim_num );

    printf ( "  %8d  %10f  %10f\n", dim_num, volume, volume2 );
  }
  return;
}
/*****************************************************************************80*/

void test35 ( )

/*****************************************************************************80

  Purpose:

    TEST35 tests SQUARE_UNIT_SET, RECTANGLE_SUB_2D.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 April 2008

  Author:

    John Burkardt
*/
{
  int i;
  int j;
  char name[8];
  int num;
  int order;
  int nsub[2];
  double result;
  int rule;
  double *weight;
  double *xtab;
  double xval[2] = { 1.0, 3.0 };
  double *ytab;
  double yval[2] = { 2.0, 3.0 };

  printf ( "\n" );
  printf ( "TEST35\n" );
  printf ( "  SQUARE_UNIT_SET sets up a quadrature rule \n" );
  printf ( "    on a unit square.\n" );
  printf ( "  RECTANGLE_SUB_2D applies it to subrectangles of an\n" );
  printf ( "    arbitrary rectangle.\n" );
  printf ( "\n" );

  printf ( "\n" );
  printf ( "  The corners of the rectangle are:\n" );
  printf ( "\n" );
  printf ( "%14f  %14f\n", xval[0], yval[0] );
  printf ( "%14f  %14f\n", xval[1], yval[1] );
/*
  Get the quadrature abscissas and weights for a unit square.
*/
  rule = 2;
  order = square_unit_size ( rule );

  xtab = ( double * ) malloc ( order * sizeof ( double ) );
  ytab = ( double * ) malloc ( order * sizeof ( double ) );
  weight = ( double * ) malloc ( order * sizeof ( double ) );

  square_unit_set ( rule, order, xtab, ytab, weight );

  printf ( "\n" );
  printf ( "  Using unit square integration rule number %d\n", rule );
  printf ( "  Order of rule is %d\n", order );
/*
  Set the function.
*/
  num = function_2d_num ( );

  for ( i = 1; i <= num; i++ )
  {
    function_2d_index = i;
    function_2d_name ( name );
/*
  Try an increasing number of subdivisions.
*/
    printf ( "\n" );
    printf ( "    Function  Subdivisions  Integral\n" );
    printf ( "\n" );

    for ( j = 1; j <= 5; j++ )
    {
      nsub[0] = j;
      nsub[1] = 2 * j;

      result = rectangle_sub_2d ( function_2d, xval, yval, nsub, order, xtab,
        ytab, weight );

      printf ( "  %8s  %4d  %4d  %14f\n", name, nsub[0], nsub[1], result );
    }
  }

  free ( weight );
  free ( xtab );
  free ( ytab );

  return;
}
/*****************************************************************************80*/

void test36 ( )

/*****************************************************************************80

  Purpose:

    TEST36 tests SQUARE_UNIT_SET and SQUARE_SUM.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 April 2008

  Author:

    John Burkardt
*/
{
  double center[2] = { 2.0, 2.0 };
  int dim_num = 2;
  int i;
  int ihi;
  int ilo;
  char name[8];
  int num;
  int order;
  double r;
  double result;
  int rule;
  int rule_max = 6;
  double *weight;
  double *xtab;
  double *ytab;

  r = 3.0;

  printf ( "\n" );
  printf ( "TEST36\n" );
  printf ( "  SQUARE_UNIT_SET sets up quadrature on the unit square;\n" );
  printf ( "  SQUARE_SUM carries it out on an arbitrary square.\n" );
  printf ( "\n" );
  printf ( "  Square center:\n" );
  printf ( "  CENTER = ( %14f, %14f )\n", center[0], center[1] );
  printf ( "  Square radius is %f\n", r );

  for ( ilo = 1; ilo <= rule_max; ilo = ilo + 5 )
  {
    ihi = i4_min ( ilo + 4, rule_max );

    printf ( "\n" );
    printf ( "  Rule:" );
    for ( rule = ilo; rule <= ihi; rule++ )
    {
      printf ( "       %6d", rule );
    }
    printf ( "\n" );
    printf ( "  Function \n" );
    printf ( "\n" );

    num = function_2d_num ( );

    for ( i = 1; i <= num; i++ )
    {
      function_2d_index = i;
      function_2d_name ( name );

      printf ( "  %8d", name );

      for ( rule = ilo; rule <= ihi; rule++ )
      {
        order = square_unit_size ( rule );

        xtab = ( double * ) malloc ( order * sizeof ( double ) );
        ytab = ( double * ) malloc ( order * sizeof ( double ) );
        weight = ( double * ) malloc ( order * sizeof ( double ) );

        square_unit_set ( rule, order, xtab, ytab, weight );

        result = square_sum ( function_2d, center, r, order, xtab, ytab,
          weight );

        printf ( "%13f", result );

        free ( weight );
        free ( xtab );
        free ( ytab );
      }
      printf ( "\n" );
    }
  }
  return;
}
/*****************************************************************************80*/

void test37 ( )

/*****************************************************************************80

  Purpose:

    TEST37 tests SQUARE_UNIT_SET and SQUARE_UNIT_SUM.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 April 2008

  Author:

    John Burkardt
*/
{
  int i;
  int ihi;
  int ilo;
  char name[8];
  int num;
  int order;
  double result;
  int rule;
  int rule_max = 6;
  double *weight;
  double *xtab;
  double *ytab;

  printf ( "\n" );
  printf ( "TEST37\n" );
  printf ( "  SQUARE_UNIT_SET sets up quadrature on the unit square;\n" );
  printf ( "  SQUARE_UNIT_SUM carries it out on the unit square.\n" );
  printf ( "\n" );

  for ( ilo = 1; ilo <= rule_max; ilo = ilo + 5 )
  {
    ihi = i4_min ( ilo + 4, rule_max );

    printf ( "\n" );
    printf ( "  Rule:   " );
    for ( rule = ilo; rule <= ihi; rule++ )
    {
      printf ( "%6d", rule );
    }
    printf ( "\n" );
    printf ( "  Function\n" );
    printf ( "\n" );

    num = function_2d_num ( );

    for ( i = 1; i <= num; i++ )
    {
      function_2d_index = i;
      function_2d_name ( name );
      printf ( "  %8s", name );

      for ( rule = ilo; rule <= ihi; rule++ )
      {
        order = square_unit_size ( rule );

        xtab = ( double * ) malloc ( order * sizeof ( double ) );
        ytab = ( double * ) malloc ( order * sizeof ( double ) );
        weight = ( double * ) malloc ( order * sizeof ( double ) );

        square_unit_set ( rule, order, xtab, ytab, weight );

        result = square_unit_sum ( function_2d, order, xtab, ytab, weight );

        printf ( "%13f", result );

        free ( weight );
        free ( xtab );
        free ( ytab );
      }
      printf ( "\n" );
    }
  }
  return;
}
/*****************************************************************************80*/

void test38 ( )

/*****************************************************************************80

  Purpose:

    TEST38 tests TETRA_07, TETRA_TPRODUCT.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 April 2008

  Author:

    John Burkardt
*/
{
  int i;
  char name[8];
  int num;
  int order;
  int order_hi;
  int order_lo;
  int order_max = 9;
  double result;
  double x[4] = { 1.0, 4.0, 1.0, 1.0 };
  double y[4] = { 2.0, 2.0, 3.0, 2.0 };
  double z[4] = { 6.0, 6.0, 6.0, 8.0 };

  printf ( "\n" );
  printf ( "TEST38\n" );
  printf ( "  For integrals inside an arbitrary tetrahedron:\n" );
  printf ( "  TETRA_07 uses a formula of degree 7;\n" );
  printf ( "  TETRA_TPRODUCT uses a triangular product formula\n" );
  printf ( "    of varying degree.\n" );
  printf ( "\n" );
  printf ( "  Tetrahedron vertices:\n" );
  printf ( "\n" );
  for ( i = 0; i < 4; i++ )
  {
    printf ( "  %4f  %4f  %4f\n", x[i], y[i], z[i] );
  }
  printf ( "\n" );
  printf ( "  Tetrahedron unit volume = %f\n", tetra_unit_volume ( ) );
  printf ( "  Tetrahedron Volume = %f\n", tetra_volume ( x, y, z ) );
  printf ( "\n" );
  printf ( "\n" );
  printf ( "  F(X)    TETRA_07\n" );
  printf ( "          TETRA_TPRODUCT(1:4)\n" );
  printf ( "          TETRA_TPRODUCT(5:8)\n" );
  printf ( "          TETRA_TPRODUCT(9)\n" );
  printf ( "\n" );

  num = function_3d_num ( );

  for ( i = 1; i <= num; i++ )
  {
    function_3d_index = i;
    function_3d_name ( name );
    printf ( "  %8s", name );
    result = tetra_07 ( function_3d, x, y, z );
    printf ( "%14f\n", result );

    for ( order_lo = 1; order_lo <= order_max; order_lo = order_lo + 4 )
    {
      order_hi = i4_min ( order_lo + 3, order_max );
      printf ( "          " );
      for ( order = order_lo; order <= order_hi; order++ )
      {
        result = tetra_tproduct ( function_3d, order, x, y, z );
        printf ( "%16f", result );
      }
      printf ( "\n" );
    }
    printf ( "\n" );
  }
  return;
}
/*****************************************************************************80*/

void test39 ( )

/*****************************************************************************80

  Purpose:

    TEST39 tests TETRA_UNIT_SET and TETRA_UNIT_SUM.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 April 2008

  Author:

    John Burkardt
*/
{
  int i;
  int ihi;
  int ilo;
  char name[8];
  int num;
  int order;
  double result;
  int rule;
  int rule_max = 8;
  double *weight;
  double *xtab;
  double *ytab;
  double *ztab;

  printf ( "\n" );
  printf ( "TEST39\n" );
  printf ( "  TETRA_UNIT_SET sets quadrature rules\n" );
  printf ( "    for the unit tetrahedron;\n" );
  printf ( "  TETRA_UNIT_SUM applies them to the unit tetrahedron.\n" );
  printf ( "\n" );

  for ( ilo = 1; ilo <= rule_max; ilo = ilo + 5 )
  {
    ihi = i4_min ( ilo +  4, rule_max );

    printf ( "\n" );
    printf ( "  Rule:   " );
    for ( rule = ilo; rule <= ihi; rule++ )
    {
      printf ( "%6d", rule );
    }
    printf ( "\n" );
    printf ( "Function\n" );
    printf ( "\n" );

    num = function_3d_num ( );

    for ( i = 1; i <= num; i++ )
    {
      function_3d_index = i;
      function_3d_name ( name );
      printf ( "  %8s", name );

      for ( rule = ilo; rule <= ihi; rule++ )
      {
        order = tetra_unit_size ( rule );

        weight = ( double * ) malloc ( order * sizeof ( double ) );
        xtab = ( double * ) malloc ( order * sizeof ( double ) );
        ytab = ( double * ) malloc ( order * sizeof ( double ) );
        ztab = ( double * ) malloc ( order * sizeof ( double ) );

        tetra_unit_set ( rule, order, xtab, ytab, ztab, weight );

        result = tetra_unit_sum ( function_3d, order, xtab, ytab, ztab,
          weight );

        printf ( "%14f", result );

        free ( weight );
        free ( xtab );
        free ( ytab );
        free ( ztab );
      }
      printf ( "\n" );
    }
  }
  return;
}
/*****************************************************************************80*/

void test40 ( )

/*****************************************************************************80

  Purpose:

    TEST40 tests TETRA_UNIT_SET and TETRA_SUM.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 April 2008

  Author:

    John Burkardt
*/
{
  int i;
  int ihi;
  int ilo;
  char name[8];
  int num;
  int order;
  double result;
  int rule;
  int rule_max = 8;
  double value;
  double *weight;
  double x[4] = { 1.0, 4.0, 1.0, 1.0 };
  double *xtab;
  double y[4] = { 2.0, 2.0, 3.0, 2.0 };
  double *ytab;
  double z[4] = { 6.0, 6.0, 6.0, 8.0 };
  double *ztab;

  printf ( "\n" );
  printf ( "TEST40\n" );
  printf ( "  TETRA_UNIT_SET sets quadrature rules\n" );
  printf ( "    for the unit tetrahedron;\n" );
  printf ( "  TETRA_SUM applies them to an arbitrary tetrahedron.\n" );
  printf ( "\n" );
  printf ( "  Tetrahedron vertices:\n" );
  printf ( "\n" );
  for ( i = 0; i < 4; i++ )
  {
    printf ( "  %6f  %6f  %6f\n", x[i], y[i], z[i] );
  }


  for ( ilo = 1; ilo <= rule_max; ilo = ilo + 5 )
  {
    ihi = i4_min ( ilo +  4, rule_max );

    printf ( "\n" );
    printf ( "  Rule:   " );
    for ( rule = ilo; rule <= ihi; rule++ )
    {
      printf ( "       %7d", rule );
    }
    printf ( "\n" );
    printf ( "  Function\n" );
    printf ( "\n" );

    num = function_3d_num ( );

    for ( i = 1; i <= num; i++ )
    {
      function_3d_index = i;
      function_3d_name ( name );
      printf ( "  %8s", name );

      for ( rule = ilo; rule <= ihi; rule++ )
      {
        order = tetra_unit_size ( rule );

        weight = ( double * ) malloc ( order * sizeof ( double ) );
        xtab = ( double * ) malloc ( order * sizeof ( double ) );
        ytab = ( double * ) malloc ( order * sizeof ( double ) );
        ztab = ( double * ) malloc ( order * sizeof ( double ) );

        tetra_unit_set ( rule, order, xtab, ytab, ztab, weight );

        result = tetra_sum ( function_3d, x, y, z, order, xtab, ytab, ztab,
          weight );

        printf ( "%14f", result );

        free ( weight );
        free ( xtab );
        free ( ytab );
        free ( ztab );
      }
      printf ( "\n" );
    }
  }
  return;
}
/*****************************************************************************80*/

void test41 ( )

/*****************************************************************************80

  Purpose:

    TEST41 tests TRIANGLE_UNIT_SET, TRIANGLE_SUB.

  Discussion:

    Break up the triangle into NSUB*NSUB equal subtriangles.  Approximate
    the integral over the triangle by the sum of the integrals over each
    subtriangle.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 April 2008

  Author:

    John Burkardt
*/
{
  int i;
  char name[8];
  int nsub;
  int num;
  int order;
  double result;
  int rule;
  double *weight;
  double *xtab;
  double xval[3] = { 0.0, 0.0, 1.0 };
  double *ytab;
  double yval[3] = { 0.0, 1.0, 0.0 };

  printf ( "\n" );
  printf ( "TEST41\n" );
  printf ( "  TRIANGLE_UNIT_SET sets up a quadrature rule\n" );
  printf ( "    on a triangle.\n" );
  printf ( "  TRIANGLE_SUB applies it to subtriangles of an\n" );
  printf ( "    arbitrary triangle.\n" );
  printf ( "\n" );
  printf ( "  Triangle vertices:\n" );
  printf ( "\n" );
  printf ( "  %12f  %12f\n", xval[0], yval[0] );
  printf ( "  %12f  %12f\n", xval[1], yval[1] );
  printf ( "  %12f  %12f\n", xval[2], yval[2] );
/*
  Get the quadrature abscissas and weights for a unit triangle.
*/
  rule = 3;
  order = triangle_unit_size ( rule );

  xtab = ( double * ) malloc ( order * sizeof ( double ) );
  ytab = ( double * ) malloc ( order * sizeof ( double ) );
  weight = ( double * ) malloc ( order * sizeof ( double ) );

  triangle_unit_set ( rule, order, xtab, ytab, weight );

  printf ( "\n" );
  printf ( "  Using unit triangle quadrature rule %d\n", rule );
  printf ( "  Rule order = %d\n", order );
  printf ( "\n" );
  printf ( "  Function Nsub  Result\n" );
  printf ( "\n" );
/*
  Set the function.
*/
  num = function_2d_num ( );

  for ( i = 1; i <= num; i++ )
  {
    function_2d_index = i;
    function_2d_name ( name );
/*
  Try an increasing number of subdivisions.
*/
    for ( nsub = 1; nsub <= 5; nsub++ )
    {
      result = triangle_sub ( function_2d, xval, yval, nsub, order, xtab,
        ytab,  weight );
      printf ( "  %8s  %4d  %14f\n", name, nsub, result );
    }
  }
  free ( xtab );
  free ( ytab );
  free ( weight );

  return;
}
/*****************************************************************************80*/

void test42 ( )

/*****************************************************************************80

  Purpose:

    TEST42 tests TRIANGLE_UNIT_SET and TRIANGLE_UNIT_SUM.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 April 2008

  Author:

    John Burkardt
*/
{
  int i;
  int ihi;
  int ilo;
  char name[8];
  int num;
  int order;
  double result;
  int rule;
  int rule_max = 20;
  double *weight;
  double *xtab;
  double *ytab;

  printf ( "\n" );
  printf ( "TEST42\n" );
  printf ( "  TRIANGLE_UNIT_SET sets up a quadrature\n" );
  printf ( "    in the unit triangle,\n" );
  printf ( "  TRIANGLE_UNIT_SUM applies it.\n" );
  printf ( "\n" );

  for ( ilo = 1; ilo <= rule_max; ilo = ilo + 5 )
  {
    ihi = i4_min ( ilo + 4, rule_max );

    printf ( "\n" );
    printf ( "  Rule:   " );
    for ( rule = ilo; rule <= ihi; rule++ )
    {
      printf ( "       %4d", rule );
    }
    printf ( "\n" );
    printf ( "Function\n" );
    printf ( "\n" );

    num = function_2d_num ( );

    for ( i = 1; i <= num; i++ )
    {
      function_2d_index = i;
      function_2d_name ( name );
      printf ( "  %8s", name );

      for ( rule = ilo; rule <= ihi; rule++ )
      {
        order = triangle_unit_size ( rule );

        xtab = ( double * ) malloc ( order * sizeof ( double ) );
        ytab = ( double * ) malloc ( order * sizeof ( double ) );
        weight = ( double * ) malloc ( order * sizeof ( double ) );

        triangle_unit_set ( rule, order, xtab, ytab, weight );

        result = triangle_unit_sum ( function_2d, order, xtab, ytab,
          weight );

        printf ( "%14f", result );

        free ( xtab );
        free ( ytab );
        free ( weight );
      }
      printf ( "\n" );
    }
  }
  return;
}
/*****************************************************************************80*/

void test425 ( )

/*****************************************************************************80

  Purpose:

    TEST425 tests TRIANGLE_UNIT_SET.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 April 2008

  Author:

    John Burkardt
*/
{
  int a;
  int b;
  double coef;
  double err;
  double exact;
  int i;
  int order;
  double quad;
  int rule;
  int rule_max = 20;
  double value;
  double *weight;
  double *xtab;
  double *ytab;

  printf ( "\n" );
  printf ( "TEST425\n" );
  printf ( "  TRIANGLE_UNIT_SET sets up a quadrature\n" );
  printf ( "    in the unit triangle,\n" );
  printf ( "\n" );
  printf ( "  Estimate integral of X^A * Y^B.\n" );

  for ( a = 0; a <= 10; a++ )
  {
    for ( b = 0; b <= 10 - a; b++ )
    {
      coef = ( double ) ( a + b + 2 ) * ( double ) ( a + b + 1 );
      for ( i = 1; i <= b; i++ )
      {
        coef = coef *( double ) ( a + i ) / ( double ) ( i );
      }

      printf ( "\n" );
      printf ( "  A = %d  B = %d\n", a, b );
      printf ( "\n" );
      printf ( "  Rule       QUAD           ERROR\n" );
      printf ( "\n" );

      for ( rule = 1; rule <= rule_max; rule++ )
      {
        order = triangle_unit_size ( rule );

        xtab = ( double * ) malloc ( order * sizeof ( double ) );
        ytab = ( double * ) malloc ( order * sizeof ( double ) );
        weight = ( double * ) malloc ( order * sizeof ( double ) );

        triangle_unit_set ( rule, order, xtab, ytab, weight );

        quad = 0.0;

        for ( i = 0; i < order; i++ )
        {
          if ( a == 0 && b == 0 )
          {
            value = coef;
          }
          else if ( a == 0 && b != 0 )
          {
            value = coef * pow ( ytab[i], b );
          }
          else if ( a != 0 && b == 0 )
          {
            value = coef * pow ( xtab[i], a );
          }
          else if ( a != 0 && b != 0 )
          {
            value = coef * pow ( xtab[i], a ) * pow ( ytab[i], b );
          }
          quad = quad + 0.5 * weight[i] * value;
        }
        exact = 1.0;
        err = r8_abs ( exact - quad );
        printf ( "  %4d  %14f  %11e\n", rule, quad, err );

        free ( xtab );
        free ( ytab );
        free ( weight );
      }
    }
  }
  return;
}
/*****************************************************************************80*/

void test43 ( )

/*****************************************************************************80

  Purpose:

    TEST43 tests TRIANGLE_UNIT_PRODUCT_SET and TRIANGLE_UNIT_SUM.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 April 2008

  Author:

    John Burkardt
*/
{
  int i;
  int ihi;
  int ilo;
  char name[8];
  int num;
  int order;
  double result;
  int rule;
  int rule_max = 8;
  double *weight;
  double *xtab;
  double *ytab;

  printf ( "\n" );
  printf ( "TEST43\n" );
  printf ( "  TRIANGLE_UNIT_PRODUCT_SET sets up a product quadrature\n" );
  printf ( "    rule in the unit triangle,\n" );
  printf ( "  TRIANGLE_UNIT_SUM applies it.\n" );
  printf ( "\n" );

  for ( ilo = 1; ilo <= rule_max; ilo = ilo + 5 )
  {
    ihi = i4_min ( ilo +  4, rule_max );

    printf ( "\n" );
    printf ( "  Rule Order: " );
    for ( rule = ilo; rule <= ihi; rule++ )
    {
      printf ( "%6d", rule );
    }
    printf ( "\n" );
    printf ( "Function\n" );
    printf ( "\n" );

    num = function_2d_num ( );

    for ( i = 1; i <= num; i++ )
    {
      function_2d_index = i;
      function_2d_name ( name );
      printf ( "  %8s", name );

      for ( rule = ilo; rule <= ihi; rule = rule++ )
      {
        order = triangle_unit_product_size ( rule );

        xtab = ( double * ) malloc ( order * sizeof ( double ) );
        ytab = ( double * ) malloc ( order * sizeof ( double ) );
        weight = ( double * ) malloc ( order * sizeof ( double ) );

        triangle_unit_product_set ( rule, order, xtab, ytab, weight );

        result = triangle_unit_sum ( function_2d, order, xtab, ytab,
          weight );

        printf ( "%14f", result );

        free ( xtab );
        free ( ytab );
        free ( weight );
      }
      printf ( "\n" );
    }
  }
  return;
}
/*****************************************************************************80*/

void test44 ( )

/*****************************************************************************80

  Purpose:

    TEST44 tests TRIANGLE_UNIT_SET and TRIANGLE_SUM.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 May 2007

  Author:

    John Burkardt
*/
{
  int i;
  int ihi;
  int ilo;
  char name[8];
  int num;
  int order;
  double result;
  int rule;
  int rule_max = 20;
  double *weight;
  double *xtab;
  double xval[3] = { 1.0, 3.0, 1.0 };
  double *ytab;
  double yval[3] = { 1.0, 1.0, 4.0 };

  printf ( "\n" );
  printf ( "TEST44\n" );
  printf ( "  TRIANGLE_UNIT_SET sets up quadrature\n" );
  printf ( "  in the unit triangle,\n" );
  printf ( "  TRIANGLE_SUM applies it to an arbitrary triangle.\n" );
  printf ( "\n" );

  for ( ilo = 1; ilo <= rule_max; ilo = ilo + 5 )
  {
    ihi = i4_min ( ilo + 4, rule_max );

    printf ( "\n" );
    printf ( "  Rule:   " );
    for ( rule = ilo; rule <= ihi; rule++ )
    {
      printf ( "%6d", rule );
    }
    printf ( "\n" );

    printf ( "Function\n" );
    printf ( "\n" );

    num = function_2d_num ( );

    for ( i = 1; i <= num; i++ )
    {
      function_2d_index = i;
      function_2d_name ( name );
      printf ( "  %8s", name );

      for ( rule = ilo; rule <= ihi; rule++ )
      {
        order = triangle_unit_size ( rule );

        xtab = ( double * ) malloc ( order * sizeof ( double ) );
        ytab = ( double * ) malloc ( order * sizeof ( double ) );
        weight = ( double * ) malloc ( order * sizeof ( double ) );

        triangle_unit_set ( rule, order, xtab, ytab, weight );

        result = triangle_sum ( function_2d, xval, yval, order, xtab, ytab,
          weight );

        printf ( "%14f", result );

        free ( xtab );
        free ( ytab );
        free ( weight );
      }
      printf ( "\n" );
    }
  }
  return;
}
/*****************************************************************************80*/

void test45 ( )

/*****************************************************************************80

  Purpose:

    TEST45 tests TORUS_1.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    03 April 2008

  Author:

    John Burkardt
*/
{
  int i;
  int j;
  int j2;
  int n;
  char name[8];
  int num;
  double result;
  double r1;
  double r2;

  r1 = 0.5;
  r2 = 1.0;
  n = 10;

  printf ( "\n" );
  printf ( "TEST45\n" );
  printf ( "  TORUS_1 approximates integrals on a torus.\n" );
  printf ( "\n" );
  printf ( "  The order N will be varied.\n" );
  printf ( "\n" );
  printf ( "  Inner radius = %f\n", r1 );
  printf ( "  Outer radius = %f\n", r2 );
  printf ( "  Area = %f\n", torus_area_3d ( r1, r2 ) );
  printf ( "\n" );
  printf ( "    F(X)  " );
  for ( j = 1; j <= 5; j++ )
  {
    j2 = 2 * ( j - 1 );
    printf ( "%14d", i4_power ( 2, j2 ) );
  }
  printf ( "\n" );
  printf ( "\n" );

  num = function_3d_num ( );

  for ( i = 1; i <= num; i++ )
  {
    function_3d_index = i;
    function_3d_name ( name );

    printf ( "  %8s", name );

    for ( j = 1; j <= 5; j++ )
    {
      j2 = 2 * ( j - 1 );
      n = i4_power ( 2, j2 );
      result = torus_1 ( function_3d, r1, r2, n );
      printf ( "  %14f", result );
    }
    printf ( "\n" );
   }
  return;
}
/*****************************************************************************80*/

void test46 ( )

/*****************************************************************************80

  Purpose:

    TEST46 tests TORUS_5S2, TORUS_6S2 and TORUS_14S.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    03 April 2008

  Author:

    John Burkardt
*/
{
  int i;
  char name[8];
  int num;
  double result1;
  double result2;
  double result3;
  double r1;
  double r2;

  r1 = 0.5;
  r2 = 1.0;

  printf ( "\n" );
  printf ( "TEST46\n" );
  printf ( "  For the interior of a torus,\n" );
  printf ( "  TORUS_5S2,\n" );
  printf ( "  TORUS_6S2, and\n" );
  printf ( "  TORUS_5S2 approximate integrals.\n" );
  printf ( "\n" );
  printf ( "  Inner radius = %f\n", r1 );
  printf ( "  Outer radius = %f\n", r2 );
  printf ( "  Volume = %f\n", torus_volume_3d ( r1, r2 ) );
  printf ( "\n" );
  printf ( "    Rule:        #5S2          #6S2          #14S\n" );
  printf ( "    F(X)\n" );
  printf ( "\n" );

  num = function_3d_num ( );

  for ( i = 1; i <= num; i++ )
  {
    function_3d_index = i;
    function_3d_name ( name );

    result1 = torus_5s2 ( function_3d, r1, r2 );
    result2 = torus_6s2 ( function_3d, r1, r2 );
    result3 = torus_14s ( function_3d, r1, r2 );

    printf ( "  %8s  %14f  %14f  %14f\n", name, result1, result2, result3 );
  }
  return;
}
/*****************************************************************************80*/

void test47 ( )

/*****************************************************************************80

  Purpose:

    TEST47 tests TORUS_SQUARE_5C2 and TORUS_SQUARE_14C.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    03 April 2008

  Author:

    John Burkardt
*/
{
  int i;
  char name[8];
  int num;
  double result1;
  double result2;
  double r1;
  double r2;

  r1 = 1.0;
  r2 = 0.125;

  printf ( "\n" );
  printf ( "TEST47\n" );
  printf ( "  For integrals inside a torus with square cross-section:\n" );
  printf ( "  TORUS_SQUARE_5C2 approximates the integral;\n" );
  printf ( "  TORUS_SQUARE_14C approximates the integral.\n" );
  printf ( "\n" );
  printf ( "  Inner radius = %f\n", r1 );
  printf ( "  Outer radius = %f\n", r2 );
  printf ( "  Volume = %f\n", torus_square_volume_3d ( r1, r2 ) );
  printf ( "\n" );
  printf ( "    F(X)    5C2           14C\n" );
  printf ( "\n" );

  num = function_3d_num ( );

  for ( i = 1; i <= num; i++ )
  {
    function_3d_index = i;
    function_3d_name ( name );

    result1 = torus_square_5c2 ( function_3d, r1, r2 );
    result2 = torus_square_14c ( function_3d, r1, r2 );

    printf ( "  %8s  %14f  %14f\n", name, result1, result2 );
  }
  return;
}
/*****************************************************************************80*/

void test48 ( )

/*****************************************************************************80

  Purpose:

    TEST48 tests TVEC_EVEN, TVEC_EVEN2 and TVEC_EVEN3.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 April 2008

  Author:

    John Burkardt
*/
{
  int nt;
  double *t;

  printf ( "\n" );
  printf ( "TEST48\n" );
  printf ( "  For evenly spaced angles between 0 and 2*PI:\n" );
  printf ( "  TVEC_EVEN\n" );
  printf ( "  TVEC_EVEN2\n" );
  printf ( "  TVEC_EVEN3\n" );

  nt = 4;
  t = tvec_even ( nt );
  r8vec_print ( nt, t, "  TVEC_EVEN:" );
  free ( t );

  nt = 4;
  t = tvec_even2 ( nt );
  r8vec_print ( nt, t, "  TVEC_EVEN2:" );
  free ( t );

  nt = 4;
  t = tvec_even3 ( nt );
  r8vec_print ( nt, t, "  TVEC_EVEN3:" );
  free ( t );

  return;
}
/*****************************************************************************80*/

void test49 ( )

/*****************************************************************************80

  Purpose:

    TEST49 tests TVEC_EVEN_BRACKET, TVEC_EVEN_BRACKET2 and TVEC_EVEN_BRACKET3.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    03 April 2008

  Author:

    John Burkardt
*/
{
  int nt;
  double *t;
  double theta1;
  double theta2;

  printf ( "\n" );
  printf ( "TEST49\n" );
  printf ( "  For evenly spaced angles between THETA1 and THETA2:\n" );
  printf ( "  TVEC_EVEN_BRACKET\n" );
  printf ( "  TVEC_EVEN_BRACKET2.\n" );
  printf ( "  TVEC_EVEN_BRACKET3.\n" );

  theta1 = 30.0;
  theta2 = 90.0;

  printf ( "\n" );
  printf ( "  THETA1 = %f\n", theta1 );
  printf ( "  THETA2 = %f\n", theta2 );

  nt = 4;
  t = tvec_even_bracket ( nt, theta1, theta2 );
  r8vec_print ( nt, t, "  TVEC_EVEN_BRACKET" );
  free ( t ) ;

  nt = 5;
  t = tvec_even_bracket2 ( nt, theta1, theta2 );
  r8vec_print ( nt, t, "  TVEC_EVEN_BRACKET2" );
  free ( t );

  nt = 3;
  t = tvec_even_bracket3 ( nt, theta1, theta2 );
  r8vec_print ( nt, t, "  TVEC_EVEN_BRACKET3" );
  free ( t );

  return;
}
/*****************************************************************************80*/

double fu18 ( double x )

/*****************************************************************************80

  Purpose:

    FU18 is the upper limit of integration for x.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 March 2008

  Author:

    John Burkardt
*/
{
  double value;

  value = 1.0;

  return value;
}
/*****************************************************************************80*/

double fl18 ( double x )

/*****************************************************************************80

  Purpose:

    FL18 is the lower limit of integration for x.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 March 2008

  Author:

    John Burkardt
*/
{
  double value;

  value = - 1.0;

  return value;
}
/*****************************************************************************80*/

double fu28 ( double x, double y )

/*****************************************************************************80

  Purpose:

    FU28 computes the upper limit of integration for (x,y).

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 March 2008

  Author:

    John Burkardt
*/
{
  double value;

  value = 1.0;

  return value;
}
/*****************************************************************************80*/

double fl28 ( double x, double y )

/*****************************************************************************80

  Purpose:

    FL28 computes the lower limit of integration for (x,y).

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 March 2008

  Author:

    John Burkardt
*/
{
  double value;

  value = -1.0;

  return value;
}
/*****************************************************************************80*/

double function_1d ( double x )

/*****************************************************************************80

  Purpose:

    FUNCTION_1D evaluates the current 1D function.

  Discussion:

    This routine assumes that the global variable FUNCTION_1D_INDEX has been
    set, and is accessible.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 April 2008

  Author:

    John Burkardt

  Parameters:

    Input, double X, the value of the variable.

    Output, double FUNCTION_1D, the value of the function.
*/
{
  double value;

  if ( function_1d_index == 1 )
  {
    value = 1.0;
  }
  else if ( function_1d_index == 2 )
  {
    value = x;
  }
  else if ( function_1d_index == 3 )
  {
    value = x * x;
  }
  else if ( function_1d_index == 4 )
  {
    value = x * x * x * x;
  }
  else if ( function_1d_index == 5 )
  {
    value = x * x * x * x;
  }
  else if ( function_1d_index == 6 )
  {
    value = x * x * x * x * x;
  }
  else if ( function_1d_index == 7 )
  {
    value = x * x * x * x * x * x;
  }
  else if ( function_1d_index == 8 )
  {
    value = r8_abs ( x );
  }
  else if ( function_1d_index == 9 )
  {
    value = sin ( x );
  }
  else if ( function_1d_index == 10 )
  {
    value = exp ( x );
  }
  else if ( function_1d_index == 11 )
  {
    value = 1.0 / ( 1.0 + r8_abs ( x ) );
  }
  else if ( function_1d_index == 12 )
  {
    value = sqrt ( r8_abs ( x ) );
  }
  else
  {
    value = 0.0;
  }

  return value;
}
/*****************************************************************************80*/

void function_1d_name ( char *name )

/*****************************************************************************80

  Purpose:

    FUNCTION_1D_NAME returns the name of the current 1D function.

  Discussion:

    This routine assumes that the global variable FUNCTION_1D_INDEX has been
    set, and is accessible.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    06 April 2008

  Author:

    John Burkardt

  Parameters:

    Output, char NAME[8], the name of the current 1D function.
*/
{
  if ( function_1d_index == 1 )
  {
    strcpy ( name, "      1" );
  }
  else if ( function_1d_index == 2 )
  {
    strcpy ( name, "      X" );
  }
  else if ( function_1d_index == 3 )
  {
    strcpy ( name, "    X^2" );
  }
  else if ( function_1d_index == 4 )
  {
    strcpy ( name, "    X^3" );
  }
  else if ( function_1d_index == 5 )
  {
    strcpy ( name, "    X^4" );
  }
  else if ( function_1d_index == 6 )
  {
    strcpy ( name, "    X^5" );
  }
  else if ( function_1d_index == 7 )
  {
    strcpy ( name, "    X^6" );
  }
  else if ( function_1d_index == 8 )
  {
    strcpy ( name, "      R" );
  }
  else if ( function_1d_index == 9 )
  {
    strcpy ( name, " SIN(X)" );
  }
  else if ( function_1d_index == 10 )
  {
    strcpy ( name, " EXP(X)" );
  }
  else if ( function_1d_index == 11 )
  {
    strcpy ( name, "1/(1+R)" );
  }
  else if ( function_1d_index == 12 )
  {
    strcpy ( name, "SQRT(R)" );
  }
  else
  {
    strcpy ( name, "???????" );
  }

  return;
}
/*****************************************************************************80*/

int function_1d_num ( )

/*****************************************************************************80

  Purpose:

    FUNCTION_1D_NUM returns the number of 1D functions.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    06 April 2008

  Author:

    John Burkardt

  Parameters:

    Output, int FUNCTION_1D_NUM, the number of 1D functions.
*/
{
  int value = 12;

  return value;
}
/*****************************************************************************80*/

double function_2d ( double x, double y )

/*****************************************************************************80

  Purpose:

    FUNCTION_2D evaluates the current 2D function.

  Discussion:

    This routine assumes that the global variable FUNCTION_2D_INDEX has been
    set, and is accessible.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 April 2008

  Author:

    John Burkardt

  Parameters:

    Input, double X, double Y, the value of the variables.

    Output, double FUNCTION_2D, the value of the function.
*/
{
  double value;

  if ( function_2d_index == 1 )
  {
    value = 1.0;
  }
  else if ( function_2d_index == 2 )
  {
    value = x;
  }
  else if ( function_2d_index == 3 )
  {
    value = x * x;
  }
  else if ( function_2d_index == 4 )
  {
    value = x * x * x;
  }
  else if ( function_2d_index == 5 )
  {
    value = x * x * x * x;
  }
  else if ( function_2d_index == 6 )
  {
    value = x * x * x * x * x;
  }
  else if ( function_2d_index == 7 )
  {
    value = x * x * x * x * x * x;
  }
  else if ( function_2d_index == 8 )
  {
    value = sqrt ( x * x + y * y );
  }
  else if ( function_2d_index == 9 )
  {
    value = sin ( x );
  }
  else if ( function_2d_index == 10 )
  {
    value = exp ( x );
  }
  else if ( function_2d_index == 11 )
  {
    value = 1.0 / ( 1.0 + sqrt ( x * x + y * y ) );
  }
  else if ( function_2d_index == 12 )
  {
    value = sqrt ( sqrt ( x * x + y * y ) );
  }
  else
  {
    value = 0.0;
  }

  return value;
}
/*****************************************************************************80*/

void function_2d_name ( char *name )

/*****************************************************************************80

  Purpose:

    FUNCTION_2D_NAME returns the name of the current 2D function.

  Discussion:

    This routine assumes that the global variable FUNCTION_2D_INDEX has been
    set, and is accessible.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 April 2008

  Author:

    John Burkardt

  Parameters:

    Output, char NAME[8], the name of the current 2D function.
*/
{
  if ( function_2d_index == 1 )
  {
    strcpy ( name, "      1" );
  }
  else if ( function_2d_index == 2 )
  {
    strcpy ( name, "      X" );
  }
  else if ( function_2d_index == 3 )
  {
    strcpy ( name, "    X^2" );
  }
  else if ( function_2d_index == 4 )
  {
    strcpy ( name, "    X^3" );
  }
  else if ( function_2d_index == 5 )
  {
    strcpy ( name, "    X^4" );
  }
  else if ( function_2d_index == 6 )
  {
    strcpy ( name, "    X^5" );
  }
  else if ( function_2d_index == 7 )
  {
    strcpy ( name, "    X^6" );
  }
  else if ( function_2d_index == 8 )
  {
    strcpy ( name, "      R" );
  }
  else if ( function_2d_index == 9 )
  {
    strcpy ( name, " SIN(X)" );
  }
  else if ( function_2d_index == 10 )
  {
    strcpy ( name, " EXP(X)" );
  }
  else if ( function_2d_index == 11 )
  {
    strcpy ( name, "1/(1+R)" );
  }
  else if ( function_2d_index == 12 )
  {
    strcpy ( name, "SQRT(R)" );
  }
  else
  {
    strcpy ( name, "???????" );
  }

  return;
}
/*****************************************************************************80*/

int function_2d_num ( )

/*****************************************************************************80

  Purpose:

    FUNCTION_2D_NUM returns the number of 2D functions.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    06 April 2008

  Author:

    John Burkardt

  Parameters:

    Output, int FUNCTION_2D_NUM, the number of 2D functions.
*/
{
  int value = 12;

  return value;
}
/*****************************************************************************80*/

double function_3d ( double x, double y, double z )

/*****************************************************************************80

  Purpose:

    FUNCTION_3D evaluates a function F(X,Y,Z) of 3 variables.

/  Discussion:

    This routine assumes that the global variable FUNCTION_3D_INDEX has been
    set, and is accessible.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 March 2008

  Author:

    John Burkardt

  Parameters:

    Input, double X, Y, Z, the value of the variables.

    Output, double FUNC_3D, the value of the function.
*/
{
  double value;

  if ( function_3d_index == 1 )
  {
    value = 1.0;
  }
  else if ( function_3d_index == 2 )
  {
    value = x;
  }
  else if ( function_3d_index == 3 )
  {
    value = y;
  }
  else if ( function_3d_index == 4 )
  {
    value = z;
  }
  else if ( function_3d_index == 5 )
  {
    value = x * x;
  }
  else if ( function_3d_index == 6 )
  {
    value = x * y;
  }
  else if ( function_3d_index == 7 )
  {
    value = x * z;
  }
  else if ( function_3d_index == 8 )
  {
    value = y * y;
  }
  else if ( function_3d_index == 9 )
  {
    value = y * z;
  }
  else if ( function_3d_index == 10 )
  {
    value = z * z;
  }
  else if ( function_3d_index == 11 )
  {
    value = x * x * x;
  }
  else if ( function_3d_index == 12 )
  {
    value = x * y * z;
  }
  else if ( function_3d_index == 13 )
  {
    value = z * z * z;
  }
  else if ( function_3d_index == 14 )
  {
    value = x * x * x * x;
  }
  else if ( function_3d_index == 15 )
  {
    value = x * x * z * z;
  }
  else if ( function_3d_index == 16 )
  {
    value = z * z * z * z;
  }
  else if ( function_3d_index == 17 )
  {
    value = x * x * x * x * x;
  }
  else if ( function_3d_index == 18 )
  {
    value = pow ( x, 6 );
  }
  else if ( function_3d_index == 19 )
  {
    value = sqrt ( x * x + y * y + z * z );
  }
  else if ( function_3d_index == 20 )
  {
    value = sin ( x );
  }
  else if ( function_3d_index == 21 )
  {
    value = exp ( x );
  }
  else if ( function_3d_index == 22 )
  {
    value = 1.0 / sqrt ( 1.0 + x * x + y * y + z * z );
  }
  else if ( function_3d_index == 23 )
  {
    value = sqrt ( sqrt ( x * x + y * y + z * z ) );
  }
  else
  {
    value = 0.0;
  }

  return value;
}
/*****************************************************************************80*/

void function_3d_name ( char *name )

/*****************************************************************************80

  Purpose:

    FUNCTION_3D_NAME returns the name of the current 3D function.

  Discussion:

    This routine assumes that the global variable FUNCTION_3D_INDEX has been
    set, and is accessible.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 March 2008

  Author:

    John Burkardt

  Parameters:

    Output, char NAME[8], the name of the current 3D function.
*/
{
  if ( function_3d_index == 1 )
  {
    strcpy ( name, "      1" );
  }
  else if ( function_3d_index == 2 )
  {
    strcpy ( name, "      X" );
  }
  else if ( function_3d_index == 3 )
  {
    strcpy ( name, "      Y" );
  }
  else if ( function_3d_index == 4 )
  {
    strcpy ( name, "      Z" );
  }
  else if ( function_3d_index == 5 )
  {
    strcpy ( name, "    X*X" );
  }
  else if ( function_3d_index == 6 )
  {
    strcpy ( name, "    X*Y" );
  }
  else if ( function_3d_index == 7 )
  {
    strcpy ( name, "    X*Z" );
  }
  else if ( function_3d_index == 8 )
  {
    strcpy ( name, "    Y*Y" );
  }
  else if ( function_3d_index == 9 )
  {
    strcpy ( name, "    Y*Z" );
  }
  else if ( function_3d_index == 10 )
  {
    strcpy ( name, "    Z*Z" );
  }
  else if ( function_3d_index == 11 )
  {
    strcpy ( name, "    X^3" );
  }
  else if ( function_3d_index == 12 )
  {
    strcpy ( name, "  X*Y*Z" );
  }
  else if ( function_3d_index == 13 )
  {
    strcpy ( name, "  Z*Z*Z" );
  }
  else if ( function_3d_index == 14 )
  {
    strcpy ( name, "    X^4" );
  }
  else if ( function_3d_index == 15 )
  {
    strcpy ( name, "X^2 Z^2" );
  }
  else if ( function_3d_index == 16 )
  {
    strcpy ( name, "    Z^4" );
  }
  else if ( function_3d_index == 17 )
  {
    strcpy ( name, "    X^5" );
  }
  else if ( function_3d_index == 18 )
  {
    strcpy ( name, "    X^6" );
  }
  else if ( function_3d_index == 19 )
  {
    strcpy ( name, "      R" );
  }
  else if ( function_3d_index == 20 )
  {
    strcpy ( name, " SIN(X)" );
  }
  else if ( function_3d_index == 21 )
  {
    strcpy ( name, " EXP(X)" );
  }
  else if ( function_3d_index == 22 )
  {
    strcpy ( name, "1/(1+R)" );
  }
  else if ( function_3d_index == 23 )
  {
    strcpy ( name, "SQRT(R)" );
  }
  else
  {
    strcpy ( name, "???????" );
  }

  return;
}
/*****************************************************************************80*/

int function_3d_num ( )

/*****************************************************************************80

  Purpose:

    FUNCTION_3D_NUM returns the number of 3D functions.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 April 2008

  Author:

    John Burkardt

  Parameters:

    Output, int FUNCTION_3D_NUM, the number of 3D functions.
*/
{
  int value = 23;

  return value;
}
/*****************************************************************************80*/

double function_nd ( int n, double x[] )

/*****************************************************************************80

  Purpose:

    FUNCTION_ND evaluates the current ND function.

  Discussion:

    This routine assumes that the global variable FUNCTION_ND_INDEX has been
    set, and is accessible.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 March 2008

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of variables.

    Input, double X[N], the value of the variables.

    Output, double FUNCTION_ND, the value of the function.
*/
{
  int i;
  double temp;
  double value;

  if ( function_nd_index == 1 )
  {
    value = 1.0;
  }
  else if ( function_nd_index == 2 )
  {
    value = x[0];
  }
  else if ( function_nd_index == 3 )
  {
    value = pow ( x[0], 2 );
  }
  else if ( function_nd_index == 4 )
  {
    value = pow ( x[0], 3 );
  }
  else if ( function_nd_index == 5 )
  {
    value = pow ( x[0], 4 );
  }
  else if ( function_nd_index == 6 )
  {
    value = pow ( x[0], 5 );
  }
  else if ( function_nd_index == 7 )
  {
    value = pow ( x[0], 6 );
  }
  else if ( function_nd_index == 8 )
  {
    temp = 0.0;
    for ( i = 0; i < n; i++ )
    {
      temp = temp + x[i] * x[i];
    }
    value = sqrt ( temp );
  }
  else if ( function_nd_index == 9 )
  {
    value = sin ( x[0] );
  }
  else if ( function_nd_index == 10 )
  {
    value = exp ( x[0] );
  }
  else if ( function_nd_index == 11 )
  {
    temp = 0.0;
    for ( i = 0; i < n; i++ )
    {
      temp = temp + x[i] * x[i];
    }
    value = 1.0 / ( 1.0 + sqrt ( temp ) );
  }
  else if ( function_nd_index == 12 )
  {
    temp = 0.0;
    for ( i = 0; i < n; i++ )
    {
      temp = temp + x[i] * x[i];
    }
    value = sqrt ( sqrt ( temp ) );
  }
  else
  {
    value = 0.0;
  }

  return value;
}
/*****************************************************************************80*/

void function_nd_name ( char *name )

/*****************************************************************************80

  Purpose:

    FUNCTION_ND_NAME returns the name of the current ND function.

  Discussion:

    This routine assumes that the global variable FUNCTION_ND_INDEX has been
    set, and is accessible.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 March 2008

  Author:

    John Burkardt

  Parameters:

    Output, char NAME[8], the name of the current ND function.
*/
{
  if ( function_nd_index == 1 )
  {
    strcpy ( name, "      1" );
  }
  else if ( function_nd_index == 2 )
  {
    strcpy ( name, "      X" );
  }
  else if ( function_nd_index == 3 )
  {
    strcpy ( name, "    X^2" );
  }
  else if ( function_nd_index == 4 )
  {
    strcpy ( name, "    X^3" );
  }
  else if ( function_nd_index == 5 )
  {
    strcpy ( name, "    X^4" );
  }
  else if ( function_nd_index == 6 )
  {
    strcpy ( name, "    X^5" );
  }
  else if ( function_nd_index == 7 )
  {
    strcpy ( name, "    X^6" );
  }
  else if ( function_nd_index == 8 )
  {
    strcpy ( name, "      R" );
  }
  else if ( function_nd_index == 9 )
  {
    strcpy ( name, " SIN(X)" );
  }
  else if ( function_nd_index == 10 )
  {
    strcpy ( name, " EXP(X)" );
  }
  else if ( function_nd_index == 11 )
  {
    strcpy ( name, "1/(1+R)" );
  }
  else if ( function_nd_index == 12 )
  {
    strcpy ( name, "SQRT(R)" );
  }
  else
  {
    strcpy ( name, "???????" );
  }

  return;
}
/*****************************************************************************80*/

int function_nd_num ( )

/*****************************************************************************80

  Purpose:

    FUNCTION_ND_NUM returns the number of ND functions.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    06 April 2008

  Author:

    John Burkardt

  Parameters:

    Output, int FUNCTION_ND_NUM, the number of ND functions.
*/
{
  int value = 12;

  return value;
}
/*****************************************************************************80*/

double f_1_2d ( double x, double y )

/*****************************************************************************80

  Purpose:

    F_1_2D evaluates the function 1 in 2D.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 May 2007

  Author:

    John Burkardt

  Parameters:

    Input, double X, Y, the arguments.

    Output, double F_1_2D, the value of the function.
*/
{
  double value;

  value = 1.0;

  return value;
}
/*****************************************************************************80*/

double f_x_2d ( double x, double y )

/*****************************************************************************80

  Purpose:

    F_X_2D evaluates the function X in 2D.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 May 2007

  Author:

    John Burkardt

  Parameters:

    Input, double X, Y, the arguments.

    Output, double F_X_2D, the value of the function.
*/
{
  double value;

  value = x;

  return value;
}
/*****************************************************************************80*/

double f_r_2d ( double x, double y )

/*****************************************************************************80

  Purpose:

    F_R_2D evaluates the function sqrt ( X**2 + Y**2 ) in 2D.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 May 2007

  Author:

    John Burkardt

  Parameters:

    Input, double X, Y, the arguments.

    Output, double F_R_2D, the value of the function.
*/
{
  double value;

  value = sqrt ( x * x + y * y );

  return value;
}
/*****************************************************************************80*/

double mono_000_3d ( int n, double x[] )

/*****************************************************************************80

  Purpose:

    MONO_000_3D evaluates X^0 Y^0 Z^0.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 March 2008

  Author:

    John Burkardt

  Parameters:

    Input, int N, the spatial dimension (which is 3 here).

    Input, double X[N], the evaluation point.

    Output, double MONO_000_3D, the value of the monomial at X.
*/
{
  double value;

  value = 1.0;

  return value;
}
/*****************************************************************************80*/

double mono_111_3d ( int n, double x[] )

/*****************************************************************************80

  Purpose:

    MONO_111_3D evaluates X^1 Y^1 Z^1.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 March 2008

  Author:

    John Burkardt

  Parameters:

    Input, int N, the spatial dimension (which is 3 here).

    Input, double X[N], the evaluation point.

    Output, double MONO_111_3D, the value of the monomial at X.
*/
{
  double value;

  value = pow ( x[0], 1 )
        * pow ( x[1], 1 )
        * pow ( x[2], 1 );

  return value;
}
/*****************************************************************************80*/

double mono_202_3d ( int n, double x[] )

/*****************************************************************************80

  Purpose:

    MONO_202_3D evaluates X^2 Y^0 Z^2.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 March 2008

  Author:

    John Burkardt

  Parameters:

    Input, int N, the spatial dimension (which is 3 here).

    Input, double X[N], the evaluation point.

    Output, double MONO_202_3D, the value of the monomial at X.
*/
{
  double value;

  value = pow ( x[0], 2 )
        * pow ( x[1], 0 )
        * pow ( x[2], 2 );

  return value;
}
/*****************************************************************************80*/

double mono_422_3d ( int n, double x[] )

/*****************************************************************************80

  Purpose:

    MONO_422_3D evaluates X^4 Y^2 Z^2.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 March 2008

  Author:

    John Burkardt

  Parameters:

    Input, int N, the spatial dimension (which is 3 here).

    Input, double X[N], the evaluation point.

    Output, double MONO_422_3D, the value of the monomial at X.
*/
{
  double value;

  value = pow ( x[0], 4 )
        * pow ( x[1], 2 )
        * pow ( x[2], 2 );

  return value;
}
/*****************************************************************************80*/

double *setsim ( int n )

/*****************************************************************************80

  Purpose:

    SETSIM defines a unit simplex.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    06 April 2008

  Author:

    John Burkardt

  Parameters:

    Input, int N, the spatial dimension.

    Output, double V[N*(N+1)], the coordinates of the N+1 vertices.
*/
{
  int i;
  int j;
  double *v;

  v = ( double * ) malloc ( n * ( n + 1 ) * sizeof ( double ) );

  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n + 1; j++ )
    {
      v[i+j*n] = 0.0;
    }
  }

  for ( i = 0; i < n; i++ )
  {
    j = i + 1;
    v[i+j*n] = 1.0;
  }

  return v;
}
