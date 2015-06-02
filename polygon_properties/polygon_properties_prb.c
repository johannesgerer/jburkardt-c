# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <string.h>

# include "polygon_properties.h"

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );
void test04 ( );
void test05 ( );
void test06 ( );
void test07 ( );
void test08 ( );
void test09 ( );
void test10 ( );
void test11 ( );
void test12 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for POLYGON_PROPERTIES_PRB.

  Discussion:

    POLYGON_PROPERTIES_PRB tests the POLYGON_PROPERTIES library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 May 2014

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "POLYGON_PROPERTIES_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the POLYGON_PROPERTIES library.\n" );

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
  test12 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "POLYGON_PROPERTIES_PRB\n" );
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

    TEST01 tests POLYGON_ANGLES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    08 May 2014

  Author:

    John Burkardt
*/
{
  double *angle;
  int i;
  int n = 6;
  double v[2*6] = {
    0.0, 0.0,
    1.0, 0.0,
    2.0, 1.0,
    3.0, 0.0,
    3.0, 2.0,
    1.0, 2.0 };

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  For a polygon:\n" );
  printf ( "  POLYGON_ANGLES computes the angles.\n" );

  printf ( "\n" );
  printf ( "  Number of polygonal vertices = %d\n", n );

  r8mat_transpose_print ( 2, n, v, "  The polygon vertices:" );

  angle = polygon_angles ( n, v );

  printf ( "\n" );
  printf ( "  Polygonal angles in degrees:\n" );
  printf ( "\n" );

  for ( i = 0; i < n; i++ )
  {
    printf ( "  %8d  %14.6g\n", i, r8_degrees ( angle[i] ) );
  }

  free ( angle );

  return;
}
/******************************************************************************/

void test02 ( )

/******************************************************************************/
/*
  Purpose:

    TEST02 tests POLYGON_AREA.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    08 May 2014

  Author:

    John Burkardt
*/
{
  double area;
  double area_exact1 = 2.0;
  double area_exact2 = 6.0;
  int n1 = 4;
  int n2 = 8;
  double v1[2*4] = {
    1.0, 0.0, 
    2.0, 1.0, 
    1.0, 2.0, 
    0.0, 1.0 };
  double v2[2*8] = {
        0.0, 0.0, 
        3.0, 0.0, 
        3.0, 3.0, 
        2.0, 3.0, 
        2.0, 1.0, 
        1.0, 1.0, 
        1.0, 2.0, 
        0.0, 2.0 };

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  For a polygon:\n" );
  printf ( "  POLYGON_AREA computes the area.\n" );

  printf ( "\n" );
  printf ( "  Number of polygonal vertices = %d\n", n1 );
  r8mat_transpose_print ( 2, n1, v1, "  The polygon vertices:" );
  area = polygon_area ( n1, v1 );
  printf ( "\n" );
  printf ( "  Exact area is        %g\n", area_exact1 );
  printf ( "  The computed area is %g\n", area );

  printf ( "\n" );
  printf ( "  Number of polygonal vertices = %d\n", n2 );
  r8mat_transpose_print ( 2, n2, v2, "  The polygon vertices:" );
  area = polygon_area ( n2, v2 );
  printf ( "\n" );
  printf ( "  Exact area is        %g\n", area_exact2 );
  printf ( "  The computed area is %g\n", area );

  return;
}
/******************************************************************************/

void test03 ( )

/******************************************************************************/
/*
  Purpose:

    TEST03 tests POLYGON_AREA_2.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    08 May 2014

  Author:

    John Burkardt
*/
{
  double area;
  double area_exact1 = 2.0;
  double area_exact2 = 6.0;
  int n1 = 4;
  int n2 = 8;
  double v1[2*4] = {
    1.0, 0.0, 
    2.0, 1.0, 
    1.0, 2.0, 
    0.0, 1.0 };
  double v2[2*8] = {
        0.0, 0.0, 
        3.0, 0.0, 
        3.0, 3.0, 
        2.0, 3.0, 
        2.0, 1.0, 
        1.0, 1.0, 
        1.0, 2.0, 
        0.0, 2.0 };

  printf ( "\n" );
  printf ( "TEST03\n" );
  printf ( "  For a polygon:\n" );
  printf ( "  POLYGON_AREA_2 computes the area.\n" );

  printf ( "\n" );
  printf ( "  Number of polygonal vertices = %d\n", n1 );
  r8mat_transpose_print ( 2, n1, v1, "  The polygon vertices:" );
  area = polygon_area_2 ( n1, v1 );
  printf ( "\n" );
  printf ( "  Exact area is        %g\n", area_exact1 );
  printf ( "  The computed area is %g\n", area );

  printf ( "\n" );
  printf ( "  Number of polygonal vertices = %d\n", n2 );
  r8mat_transpose_print ( 2, n2, v2, "  The polygon vertices:" );
  area = polygon_area_2 ( n2, v2 );
  printf ( "\n" );
  printf ( "  Exact area is        %g\n", area_exact2 );
  printf ( "  The computed area is %g\n", area );

  return;
}
/******************************************************************************/

void test04 ( )

/******************************************************************************/
/*
  Purpose:

    TEST04 tests POLYGON_CENTROID and POLYGON_CENTROID_2.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    08 May 2014

  Author:

    John Burkardt
*/
{
  double *centroid;
  int n = 4;
  double v[2*4] = {
    1.0, 0.0,
    2.0, 1.0,
    1.0, 2.0,
    0.0, 1.0 };

  printf ( "\n" );
  printf ( "TEST04\n" );
  printf ( "  For a polygon:\n" );
  printf ( "  POLYGON_CENTROID computes the centroid.\n" );
  printf ( "  POLYGON_CENTROID_2 computes the centroid.\n" );

  r8mat_transpose_print ( 2, n, v, "  The polygon vertices:" );

  centroid = polygon_centroid ( n, v );
  r8vec_print ( 2, centroid, "  POLYGON_CENTROID:" );
  free ( centroid );

  centroid = polygon_centroid_2 ( n, v );
  r8vec_print ( 2, centroid, "  POLYGON_CENTROID_2:" );
  free ( centroid );

  return;
}
/******************************************************************************/

void test05 ( )

/******************************************************************************/
/*
  Purpose:

    TEST05 tests POLYGON_CONTAINS_POINT and POLYGON_CONTAINS_POINT_2.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    08 May 2014

  Author:

    John Burkardt
*/
{
  int i;
  int inside1;
  int inside2;
  int n = 5;
  double p[2];
  double p_test[2*4] = {
    1.0,  1.0, 
    3.0,  4.0, 
    0.0,  2.0, 
    0.5, -0.25 };
  int test;
  int test_num = 4;
  double v[2*5] = {
    0.0, 0.0, 
    1.0, 0.0, 
    2.0, 1.0, 
    1.0, 2.0, 
    0.0, 2.0 };
 
  printf ( "\n" );
  printf ( "TEST05\n" );
  printf ( "  POLYGON_CONTAINS_POINT determines if\n" );
  printf ( "  a point is in a polygon.\n" );
  printf ( "  POLYGON_CONTAINS_POINT_2 determines if\n" );
  printf ( "  a point is in a polygon.\n" );

  r8mat_transpose_print ( 2, n, v, "  The polygon vertices:" );

  printf ( "\n" );
  printf ( "          P          In1  In2\n" );
  printf ( "\n" );

  for ( test = 0; test < test_num; test++ )
  {
    p[0] = p_test[0+test*2];
    p[1] = p_test[1+test*2];
 
    inside1 = polygon_contains_point ( n, v, p );

    inside2 = polygon_contains_point_2 ( n, v, p );

    printf ( "  %14.6g  %14.6g    %d    %d\n", p[0], p[1], inside1, inside2 );
  } 

  return;
}
/******************************************************************************/

void test06 ( )

/******************************************************************************/
/*
  Purpose:

    TEST06 tests POLYGON_DIAMETER;

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    08 May 2014

  Author:

    John Burkardt
*/
{
  double diameter;
  double diameter_exact = 2.0;
  int n = 4;
  double v[2*4] = {
    1.0, 0.0, 
    2.0, 1.0, 
    1.0, 2.0, 
    0.0, 1.0 };

  printf ( "\n" );
  printf ( "TEST06\n" );
  printf ( "  For a polygon:\n" );
  printf ( "  POLYGON_DIAMETER computes the diameter;\n" );

  r8mat_transpose_print ( 2, n, v, "  The polygon vertices:" );

  diameter = polygon_diameter ( n, v );

  printf ( "\n" );
  printf ( "  Diameter ( computed ) %g\n", diameter );
  printf ( "  Diameter ( exact )    %g\n", diameter_exact );
 
  return;
}
/******************************************************************************/

void test07 ( )

/******************************************************************************/
/*
  Purpose:

    TEST07 tests POLYGON_EXPAND;

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    08 May 2014

  Author:

    John Burkardt
*/
{
  double h;
  int n = 4;
  double v[2*4] = {
    1.0, 1.0, 
    5.0, 1.0, 
    2.0, 4.0, 
    1.0, 3.0 };
  double *w;

  printf ( "\n" );
  printf ( "TEST07\n" );
  printf ( "  For a polygon:\n" );
  printf ( "  POLYGON_EXPAND expands it by an amount H.\n" );

  h = 0.5;

  r8mat_transpose_print ( 2, n, v, "  The polygon vertices:" );

  printf ( "\n" );
  printf ( "  The expansion amount H = %g\n", h );

  w = polygon_expand ( n, v, h );

  r8mat_transpose_print ( 2, n, w, "  The expanded polygon:" );

  free ( w );

  return;
}
/******************************************************************************/

void test08 ( )

/******************************************************************************/
/*
  Purpose:

    TEST08 tests POLYGON_INRAD_DATA, POLYGON_OUTRAD_DATA, POLYGON_SIDE_DATA;

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    08 May 2014

  Author:

    John Burkardt
*/
{
  double area;
  int n;
  double radin;
  double radout;
  double side;

  printf ( "\n" );
  printf ( "TEST08\n" );
  printf ( "  For a REGULAR polygon:\n" );
  printf ( "  the inradius, outradius and side are related.\n" );
  printf ( "  POLYGON_INRAD_DATA uses the inradius;\n" );
  printf ( "  POLYGON_OUTRAD_DATA uses the inradius;\n" );
  printf ( "  POLYGON_SIDE_DATA uses the inradius;\n" );

  for ( n = 3; n <= 5; n++ )
  {
    printf ( "\n" );
    printf ( "  Number of polygonal sides = %d\n", n );
    side = 1.0;
    printf ( "\n" );
    printf ( "  Assuming SIDE = %g\n", side );
    polygon_side_data ( n, side, &area, &radin, &radout );
    printf ( "    AREA =   %g\n", area );
    printf ( "    RADIN =  %g\n", radin );
    printf ( "    RADOUT = %g\n", radout );
    printf ( "\n" );
    printf ( "  Assuming RADIN = %g\n", radin );
    polygon_inrad_data ( n, radin, &area, &radout, &side );
    printf ( "    AREA =   %g\n", area );
    printf ( "    RADOUT = %g\n", radout );
    printf ( "    SIDE =   %g\n", side );
    printf ( "\n" );
    printf ( "  Assuming RADOUT = %g\n", radout );
    polygon_outrad_data ( n, radout, &area, &radin, &side );
    printf ( "    AREA =   %g\n", area );
    printf ( "    RADIN =  %g\n", radin );
    printf ( "    SIDE =   %g\n", side );
  }
  return;
}
/******************************************************************************/

void test09 ( )

/******************************************************************************/
/*
  Purpose:

    TEST09 tests POLYGON_INTEGRAL_*.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    08 May 2014

  Author:

    John Burkardt
*/
{
  int n = 4;
  double result;
  double v[2*4] = {
    0.0, 0.0, 
    1.0, 0.0, 
    1.0, 1.0, 
    0.0, 1.0 };

  printf ( "\n" );
  printf ( "TEST09\n" );
  printf ( "  For a polygon:\n" );
  printf ( "  POLYGON_INTEGRAL_1 integrates 1\n" );
  printf ( "  POLYGON_INTEGRAL_X integrates X\n" );
  printf ( "  POLYGON_INTEGRAL_Y integrates Y\n" );
  printf ( "  POLYGON_INTEGRAL_XX integrates X*X\n" );
  printf ( "  POLYGON_INTEGRAL_XY integrates X*Y\n" );
  printf ( "  POLYGON_INTEGRAL_YY integrates Y*Y\n" );

  r8mat_transpose_print ( 2, n, v, "  The polygon vertices:" );

  printf ( "\n" );
  printf ( "  F(X,Y)    Integral\n" );
  printf ( "\n" );

  result = polygon_integral_1 ( n, v );
  printf ( "    1    %g\n", result );

  result = polygon_integral_x ( n, v );
  printf ( "    X    %g\n", result );

  result = polygon_integral_y ( n, v );
  printf ( "    Y    %g\n", result );

  result = polygon_integral_xx ( n, v );
  printf ( "  X*X    %g\n", result );

  result = polygon_integral_xy ( n, v );
  printf ( "  X*Y    %g\n", result );

  result = polygon_integral_yy ( n, v );
  printf ( "  Y*Y    %g\n", result );

  return;
}
/******************************************************************************/

void test10 ( )

/******************************************************************************/
/*
  Purpose:

    TEST10 tests POLYGON_IS_CONVEX.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    09 May 2014

  Author:

    John Burkardt
*/
{ 
  double angle;
  int i;
  int n;
  int n01 = 1;
  int n02 = 2;
  int n03 = 3;
  int n04 = 3;
  int n05 = 3;
  int n06 = 4;
  int n07 = 5;
  int n08 = 5;
  int n09 = 6;
  int n10 = 6;
  int n11 = 8;
  const double r8_pi = 3.141592653589793;
  int result;
  int test;
  int test_num = 11;
  char title[255];
  double *v;
  double v01[2*1] = {
    0.0, 0.0 };
  double v02[2*2] = {
    0.0, 0.0, 
    1.0, 2.0 };
  double v03[2*3] = {
    0.0, 0.0, 
    2.0, 0.0, 
    1.0, 0.0 };
  double v04[2*3] = {
    0.0, 0.0, 
    1.0, 0.0, 
    0.0, 2.0 };
  double v05[2*3] = {
    0.0, 0.0, 
    0.0, 2.0, 
    1.0, 0.0 };
  double v06[2*4] = {
    1.0, 0.0, 
    2.0, 0.0, 
    3.0, 1.0, 
    0.0, 1.0 };
  double v07[2*5] = {
    0.0, 0.0, 
    0.5, 0.5, 
    1.0, 0.0, 
    1.0, 1.0, 
    0.0, 1.0 };
  double *v08;
  double *v09;
  double v10[2*6] = {
    0.0, 0.0, 
    2.0, 0.0, 
    1.0, 1.0, 
    0.0, 0.0, 
    2.0, 0.0, 
    1.0, 1.0 };
  double v11[2*8] = { 
    1.0, 0.0, 
    3.0, 0.0, 
    3.0, 3.0, 
    0.0, 3.0, 
    0.0, 1.0, 
    2.0, 1.0, 
    2.0, 2.0, 
    1.0, 2.0 };

  printf ( "\n" );
  printf ( "TEST10\n" );
  printf ( "  POLYGON_IS_CONVEX determines if a polygon\n" );
  printf ( "  is convex.\n" );

  for ( test = 1; test <= test_num; test++ )
  {
    if ( test == 1 )
    {
      n = n01;
      v = v01;
      strcpy ( title, "  A point:" );
    }
    else if ( test == 2 )
    {
      n = n02;
      v = v02;
      strcpy ( title, "  A line:" );
    }
    else if ( test == 3 )
    {
      n = n03;
      v = v03;
      strcpy ( title, "  A triangle:" );
    }
    else if ( test == 4 )
    {
      n = n04;
      v = v04;
      strcpy ( title, "  A CCW triangle:" );
    }
    else if ( test == 5 )
    {
      n = n05;
      v = v05;
      strcpy ( title, "  A CW triangle:" );
    }
    else if ( test == 6 )
    {
      n = n06;
      v = v06;
      strcpy ( title, "  Polygon with large angle:" );
    }
    else if ( test == 7 )
    {
      n = n07;
      v = v07;
      strcpy ( title, "  Polygon with huge angle:" );
    }
    else if ( test == 8 )
    {
      n = n08;
      v08 = ( double * ) malloc ( 2 * n * sizeof ( double ) );
      for ( i = 0; i < n; i++ )
      {
        angle = ( double ) ( i ) * 4.0 * r8_pi / ( double ) ( n );
        v08[0+i*2] = cos ( angle );
        v08[1+i*2] = sin ( angle );
      }
      v = v08;
      strcpy ( title, "  A five-pointed star:" );
    }
    else if ( test == 9 )
    {
      n = n09;
      v09 = ( double * ) malloc ( 2 * n * sizeof ( double ) );
      for ( i = 0; i < n; i++ )
      {
        angle = ( double ) ( i ) * 2.0 * r8_pi / ( double ) ( n );
        v09[0+i*2] = cos ( angle );
        v09[1+i*2] = sin ( angle );
      }
      v = v09;
      strcpy ( title, "  A hexagon:" );
    }
    else if ( test == 10 )
    {
      n = n10;
      v = v10;
      strcpy ( title, "  A triangle twice:" );
    }
    else if ( test == 11 )
    {
      n = n11;
      v = v11;
      strcpy ( title, "  Square knot:" );
    }

    r8mat_transpose_print ( 2, n, v, title );
    result = polygon_is_convex ( n, v );

    if ( result == -1 )
    {
      printf ( "  The polygon is not convex.\n" );
    }
    else if ( result == 0 )
    {
      printf ( "  The polygon is degenerate and convex.\n" );
    }
    else if ( result == 1 )
    {
      printf ( "  The polygon is convex and counterclockwise.\n" );
    }
    else if ( result == 2 )
    {
      printf ( "  The polygon is convex and clockwise.\n" );
    }
  }

  free ( v08 );
  free ( v09 );

  return;
}
/******************************************************************************/

void test11 ( )

/******************************************************************************/
/*
  Purpose:

    TEST11 tests POLYGON_LATTICE_AREA.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    09 May 2014

  Author:

    John Burkardt
*/
{
  double area;
  int b;
  int i;

  printf ( "\n" );
  printf ( "TEST11\n" );
  printf ( "  POLYGON_LATTICE_AREA returns the area\n" );
  printf ( "  of a polygon, measured in lattice points.\n" );

  i = 5;
  b = 6;

  printf ( "\n" );
  printf ( "  Number of interior lattice points = %d\n", i );
  printf ( "  Number of boundary lattice points = %d\n", b );

  area = polygon_lattice_area ( i, b );

  printf ( "  Area of polygon is %g\n", area );

  return;
}
/******************************************************************************/

void test12 ( )

/******************************************************************************/
/*
  Purpose:

    TEST12 tests POLYGON_SAMPLE.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    08 May 2014

  Author:

    John Burkardt
*/
{
  int n = 20;
  int nv = 6;
  int seed;
  double v[2*6] = {
    0.0, 0.0, 
    2.0, 0.0, 
    2.0, 1.0, 
    1.0, 1.0, 
    1.0, 2.0, 
    0.0, 1.0 };
  double *x;

  printf ( "\n" );
  printf ( "TEST12\n" );
  printf ( "  POLYGON_SAMPLE samples a polygon.\n" );

  seed = 123456789;

  x = polygon_sample ( nv, v, n, &seed );

  r8mat_transpose_print ( 2, n, x, "  Sample points:" );

  free ( x );

  return;
}
