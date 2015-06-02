# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "geometry.h"

int main ( );
void test0005 ( );
void test001 ( );
void test002 ( );
void test0351 ( );
void test0352 ( );
void test0477 ( );
void test0478 ( );
void test0616 ( );
void test0617 ( );
void test171 ( );
void test1712 ( );
void test1835 ( );
void test1836 ( );
void test2101 ( );
void test21011 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for GEOMETRY_PRB.

  Discussion:

    GEOMETRY_PRB tests the GEOMETRY library.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    14 April 2013

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "GEOMETRY_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the GEOMETRY library.\n" );

  test0005 ( );
  test001 ( );
  test002 ( );

  test0351 ( );
  test0352 ( );

  test0477 ( );
  test0478 ( );

  test0616 ( );
  test0617 ( );

  test171 ( );
  test1712 ( );

  test1835 ( );
  test1836 ( );

  test2101 ( );
  test21011 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "GEOMETRY_PRB\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void test0005 ( )

/******************************************************************************/
/*
  Purpose:

    TEST0005 tests ANGLE_BOX_2D.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    03 October 2010

  Author:

    John Burkardt
*/
{
# define DIM_NUM 2

  double dist;
  double p1[DIM_NUM];
  double p2[DIM_NUM];
  double p3[DIM_NUM];
  double p4[DIM_NUM];
  double p5[DIM_NUM];

  printf ( "\n" );
  printf ( "TEST0005\n" );
  printf ( "  ANGLE_BOX_2D\n" );
  printf ( "\n" );
  printf ( "  Compute P4 and P5, normal to\n" );
  printf ( "  line through P1 and P2, and\n" );
  printf ( "  line through P2 and P3,\n" );
  printf ( "  and DIST units from P2.\n" );
/*
  These points define the lines
    y = 0
  and
    y = 2x-6
*/
  p1[0] = 0.0;
  p1[1] = 0.0;
  p2[0] = 3.0;
  p2[1] = 0.0;
  p3[0] = 4.0;
  p3[1] = 2.0;
  dist = 1.0;

  printf ( "\n" );
  printf ( "  DIST %14f\n", dist );
  printf ( "  P1:  %14f  %14f\n", p1[0], p1[1] );
  printf ( "  P2:  %14f  %14f\n", p2[0], p2[1] );
  printf ( "  P3:  %14f  %14f\n", p3[0], p3[1] );

  angle_box_2d ( dist, p1, p2, p3, p4, p5 );

  printf ( "  P4:  %14f  %14f\n", p4[0], p4[1] );
  printf ( "  P5:  %14f  %14f\n", p5[0], p5[1] );
/*
  These points define the lines
    y = 0
  and
    y = 2x-6
*/
  p1[0] = 0.0;
  p1[1] = 0.0;
  p2[0] = 3.0;
  p2[1] = 0.0;
  p3[0] = 2.0;
  p3[1] = -2.0;
  dist = 1.0;

  printf ( "\n" );
  printf ( "  DIST %14f\n", dist );
  printf ( "  P1:  %14f  %14f\n", p1[0], p1[1] );
  printf ( "  P2:  %14f  %14f\n", p2[0], p2[1] );
  printf ( "  P3:  %14f  %14f\n", p3[0], p3[1] );

  angle_box_2d ( dist, p1, p2, p3, p4, p5 );

  printf ( "  P4:  %14f  %14f\n", p4[0], p4[1] );
  printf ( "  P5:  %14f  %14f\n", p5[0], p5[1] );
/*
  By setting P1 = P2, we are asking to extend the line
    y = 2x-6
  from P3 to P2 through to the other side.
*/
  p1[0] = 3.0;
  p1[1] = 0.0;
  p2[0] = 3.0;
  p2[1] = 0.0;
  p3[0] = 2.0;
  p3[1] = -2.0;
  dist = 1.0;

  printf ( "\n" );
  printf ( "  DIST %14f\n", dist );
  printf ( "  P1:  %14f  %14f\n", p1[0], p1[1] );
  printf ( "  P2:  %14f  %14f\n", p2[0], p2[1] );
  printf ( "  P3:  %14f  %14f\n", p3[0], p3[1] );

  angle_box_2d ( dist, p1, p2, p3, p4, p5 );

  printf ( "  P4:  %14f  %14f\n", p4[0], p4[1] );
  printf ( "  P5:  %14f  %14f\n", p5[0], p5[1] );

  return;
# undef DIM_NUM
}
/******************************************************************************/

void test001 ( )

/******************************************************************************/
/*
  Purpose:

    TEST001 tests ANGLE_CONTAINS_RAY_2D.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    23 October 2012

  Author:

    John Burkardt
*/
{
# define DIM_NUM 2
# define TEST_NUM 6

  int angle;
  int angle_num = 12;
  double angle_rad;
  int inside;
  double p[DIM_NUM];
  double p1[DIM_NUM];
  double p2[DIM_NUM];
  double p3[DIM_NUM];
  double pi = 3.141592653589793;
  int test;

  printf ( "\n" );
  printf ( "TEST001\n" );
  printf ( "  ANGLE_CONTAINS_RAY_2D sees if a ray lies within an angle.\n" );

  for ( test = 0; test < TEST_NUM; test++ )
  {
    if ( test == 0 )
    {
      p1[0] = 1.0;
      p1[1] = 0.0;

      p2[0] = 0.0;
      p2[1] = 0.0;

      p3[0] = 1.0;
      p3[1] = 1.0;
    }
    else if ( test == 1 )
    {
      p1[0] = 1.0;
      p1[1] = 0.0;

      p2[0] = 0.0;
      p2[1] = 0.0;

      p3[0] = 0.0;
      p3[1] = 1.0;
    }
    else if ( test == 2 )
    {
      p1[0] =  1.0;
      p1[1] = -1.0;

      p2[0] = 0.0;
      p2[1] = 0.0;

      p3[0] = 0.0;
      p3[1] = 1.0;
    }
    else if ( test == 3 )
    {
      p1[0] = 1.0;
      p1[1] = 0.0;

      p2[0] = 0.0;
      p2[1] = 0.0;

      p3[0] = -1.0;
      p3[1] =  0.0;
    }
    else if ( test == 4 )
    {
      p1[0] = 1.0;
      p1[1] = 0.0;

      p2[0] = 0.0;
      p2[1] = 0.0;

      p3[0] =  0.0;
      p3[1] = -1.0;
    }
    else if ( test == 5 )
    {
      p1[0] = 1.0;
      p1[1] = 0.0;

      p2[0] = 0.0;
      p2[1] = 0.0;

      p3[0] =  1.0;
      p3[1] = -0.01;
    }

    r8vec_print ( DIM_NUM, p1, "  Vertex A" );
    r8vec_print ( DIM_NUM, p2, "  Vertex B" );
    r8vec_print ( DIM_NUM, p3, "  Vertex C" );

    printf ( "\n" );
    printf ( "       X            Y       Inside?\n" );
    printf ( "\n" );

    for ( angle = 0; angle <= angle_num; angle++ )
    {
      angle_rad = ( double ) ( angle ) * 2.0 * pi / ( double ) angle_num;

      p[0] = cos ( angle_rad );
      p[1] = sin ( angle_rad );

      inside = angle_contains_ray_2d ( p1, p2, p3, p );

      printf ( "  %12g  %12g  %d\n", p[0], p[1], inside );
    }

  }

  return;
# undef DIM_NUM
# undef TEST_NUM
}
/******************************************************************************/

void test002 ( )

/******************************************************************************/
/*
  Purpose:

    TEST002 tests ANGLE_DEG_2D and ANGLE_RAD_ND.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    14 April 2013

  Author:

    John Burkardt
*/
{
  int angle_num = 12;
  int i;
  double temp1;
  double temp2;
  double temp3;
  double thetad;
  double thetar;
  double v1[2] = { 1.0, 0.0 };
  double v2[2];
  double v3[2] = { 0.0, 0.0 };

  printf ( "\n" );
  printf ( "TEST002\n" );
  printf ( "  ANGLE_DEG_2D computes an angle,\n" );
  printf ( "  ANGLE_RAD_ND computes an angle.\n" );
  printf ( "\n" );
  printf ( "  X  Y  Theta  ATAN2(y, x), ANGLE_RAD_ND, ANGLE_DEG_2D\n" );
  printf ( "\n" );

  for ( i = 0; i <= angle_num; i++ )
  {
    thetad = ( double ) ( i ) * 360.0 / ( double ) ( angle_num );
    thetar = degrees_to_radians ( thetad );

    v2[0] = cos ( thetar );
    v2[1] = sin ( thetar );

    temp1 = radians_to_degrees ( atan2 ( v2[1], v2[0] ) );

    temp2 = angle_rad_nd ( 2, v1, v2 );

    temp3 = angle_deg_2d ( v1, v3, v2 );

    printf ( "  %8.3f  %8.3f  %8.3f  %8.3f  %8.3f  %8.3f\n",
    v2[0], v2[1], thetad, temp1, temp2, temp3 );
  } 
  return;
}
/******************************************************************************/

void test0351 ( )

/******************************************************************************/
/*
  Purpose:

    TEST0351 tests LINE_PAR_POINT_NEAR_2D and LINE_PAR_POINT_DIST_2D.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    14 April 2013

  Author:

    John Burkardt
*/
{
# define TEST_NUM 3

  double dist;
  double f;
  double g;
  int i;
  double p[2];
  double *pn;
  double p_test[2*TEST_NUM] = {
    0.0,  0.0,
    5.0, -1.0,
    5.0,  3.0 };
  double t;
  int test;
  int test_num = TEST_NUM;
  double x0;
  double y0;

  printf ( "\n" );
  printf ( "TEST0351\n" );
  printf ( "  LINE_PAR_POINT_NEAR_2D finds the point on\n" );
  printf ( "  a parametric line (X0,Y0,F,G) nearest a point P in 2D.\n" );

  x0 = 1.0;
  y0 = 3.0;
  f = +1.0;
  g = -1.0;

  printf ( "\n" );
  printf ( "  Parametric line:\n" );
  printf ( "  X(t) = %g + %g * t\n", x0, f );
  printf ( "  Y(t) = %g + %g * t\n", y0, g );

  for ( test = 0; test < test_num; test++ )
  {
    for ( i = 0; i < 2; i++ )
    {
      p[i] = p_test[i+test*2];
    }

    r8vec_print ( 2, p, "  The point P:" );

    dist = line_par_point_dist_2d ( f, g, x0, y0, p );

    printf ( "  Distance = %g\n", dist );

    pn = line_par_point_near_2d ( f, g, x0, y0, p );

    r8vec_print ( 2, pn, "  Nearest point PN:" );

    dist = r8vec_norm_affine ( 2, p, pn );

    printf ( "  Distance recomputed = %g\n", dist );

    free ( pn );
  }
  return;
# undef TEST_NUM
}
/******************************************************************************/

void test0352 ( )

/******************************************************************************/
/*
  Purpose:

    TEST0352 tests LINE_PAR_POINT_DIST_3D and LINE_PAR_POINT_NEAR_3D.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    14 April 2013

  Author:

    John Burkardt
*/
{
# define TEST_NUM 3

  double dist;
  double f;
  double g;
  double h;
  int i;
  double p[3];
  double *pn;
  double p_test[3*TEST_NUM] = {
    0.0,  0.0, 2.0, 
    5.0, -1.0, 1.0, 
    5.0,  3.0, 3.0 };
  int test;
  int test_num = TEST_NUM;
  double x0;
  double y0;
  double z0;

  printf ( "\n" );
  printf ( "TEST0352\n" );
  printf ( "  LINE_PAR_POINT_DIST_3D finds the distance\n" );
  printf ( "  from a parametric line to a point in 3D.\n" );

  x0 = 1.0;
  y0 = 3.0;
  z0 = 2.0;

  f = +3.0;
  g = -3.0;
  h = -1.0;

  printf ( "\n" );
  printf ( "  Parametric line:\n" );
  printf ( "  X(t) = %g + %g * t\n", x0, f );
  printf ( "  Y(t) = %g + %g * t\n", y0, g );
  printf ( "  Z(t) = %g + %g * t\n", z0, h );

  for ( test = 0; test < test_num; test++ )
  {
    for ( i = 0; i < 3; i++ )
    {
      p[i] = p_test[i+test*3];
    }

    r8vec_print ( 3, p, "  The point P:" );

    dist = line_par_point_dist_3d ( f, g, h, x0, y0, z0, p );

    printf ( "  Distance = %g\n", dist );

    pn = line_par_point_near_3d ( f, g, h, x0, y0, z0, p );

    r8vec_print ( 3, pn, "  Nearest point PN:" );

    dist = r8vec_norm_affine ( 3, p, pn );

    printf ( "  Distance recomputed = %g\n", dist );

    free ( pn );
  }
  return;
# undef TEST_NUM
}
/******************************************************************************/

void test0477 ( )

/******************************************************************************/
/*
  Purpose:

    TEST0477 tests PARALLELOGRAM_AREA_2D.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    09 May 2010

  Author:

    John Burkardt
*/
{
  double area;
  double p[2*4] = {
    2.0, 7.0, 
    5.0, 7.0, 
    6.0, 9.0, 
    3.0, 9.0 };

  printf ( "\n" );
  printf ( "TEST0477\n" );
  printf ( "  PARALLELOGRAM_AREA_2D finds the area of a\n" );
  printf ( "  parallelogram in 2D.\n" );

  r8mat_transpose_print ( 2, 4, p, "  Vertices:" );

  area = parallelogram_area_2d ( p );

  printf ( "\n" );
  printf ( "  AREA = %f\n", area );

  return;
}
/******************************************************************************/

void test0478 ( )

/******************************************************************************/
/*
  Purpose:

    TEST0478 tests PARALLELOGRAM_AREA_3D.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    09 May 2010

  Author:

    John Burkardt
*/
{
  double area;
  double p[3*4] = {
    1.0,       2.0,       3.0, 
    2.4142137, 3.4142137, 3.0, 
    1.7071068, 2.7071068, 4.0, 
    0.2928931, 0.2928931, 4.0 };

  printf ( "\n" );
  printf ( "TEST0478\n" );
  printf ( "  PARALLELOGRAM_AREA_3D finds the area of a\n" );
  printf ( "  parallelogram in 3D.\n" );

  r8mat_transpose_print ( 3, 4, p, "  Vertices:" );

  area = parallelogram_area_3d ( p );

  printf ( "\n" );
  printf ( "  AREA = %f\n", area );

  return;
}
/******************************************************************************/

void test0616 ( )

/******************************************************************************/
/*
  Purpose:

    TEST0616 tests PLANE_NORMAL_QR_TO_XYZ and PLANE_NORMAL_XYZ_TO_QR.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 November 2010

  Author:

    John Burkardt
*/
{
  double dif;
  int i;
  int j;
  int m = 3;
  int n = 5;
  double *normal;
  double *pp;
  double pq[3];
  double pr[3];
  double *qr1;
  double *qr2;
  int seed;
  double t;
  double *xyz;

  seed = 123456789;

  printf ( "\n" );
  printf ( "TEST0616\n" );
  printf ( "  For a normal plane, with point PP and NORMAL vector,\n" );
  printf ( "  and in-plane basis vectors PQ and PR,\n" );
  printf ( "  PLANE_NORMAL_QR_TO_XYZ converts QR to XYZ coordinates;\n" );
  printf ( "  PLANE_NORMAL_XYZ_TO_QR converts XYZ to QR coordinates.\n" );
//
//  Choose PP and NORMAL at random.
//
  pp = r8vec_uniform_01_new ( m, &seed );

  normal = r8vec_uniform_01_new ( m, &seed );
//
//  Compute in-plane basis vectors PQ and PR.
//
  plane_normal_basis_3d ( pp, normal, pq, pr );
//
//  Choose random Q, R coordinates.
//
  qr1 = r8mat_uniform_01_new ( m - 1, n, &seed );
//
//  Convert to XYZ.
//
  xyz = plane_normal_qr_to_xyz ( pp, normal, pq, pr, n, qr1 );
//
//  Convert XYZ to QR.
//
  qr2 = plane_normal_xyz_to_qr ( pp, normal, pq, pr, n, xyz );

  dif = 0.0;
  for ( j = 0; j < n; j++ )
  {
    t = 0.0;
    for ( i = 0; i < m - 1; i++ )
    {
      t = t + pow ( qr1[0+j*2] - qr2[0+j*2], 2 );
    }
    t = sqrt ( t );
    dif = r8_max ( dif, t );
  }

  printf ( "\n" );
  printf ( "  Maximum difference was %f\n", dif );

  free ( normal );
  free ( pp );
  free ( qr1 );
  free ( qr2 );
  free ( xyz );

  return;
}
/******************************************************************************/

void test0617 ( )

/******************************************************************************/
/*
  Purpose:

    TEST0617 tests PLANE_NORMAL_TETRAHEDRON_INTERSECT.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    24 June 2010

  Author:

    John Burkardt
*/
{
  int i;
  int j;
  int k;
  int l;
  int int_num;
  double normal[3];
  double pint[3*4];
  double pp[3];
  double t[3*4] = {
    0.0, 0.0, 0.0, 
    1.0, 0.0, 0.0, 
    0.0, 1.0, 0.0, 
    0.0, 0.0, 1.0 };

  fprintf ( stdout, "\n" );
  fprintf ( stdout, "TEST0617\n" );
  fprintf ( stdout, "  PLANE_NORMAL_TETRAHEDRON_INTERSECT determines\n" );
  fprintf ( stdout, "  the intersection of a plane and tetrahedron.\n" );

  for ( k = 1; k <= 2; k++ )
  {
    if ( k == 1 )
    {
      normal[0] = 0.0;
      normal[1] = 0.0;
      normal[2] = 1.0;
    }
    else
    {
      normal[0] = 1.0 / sqrt ( 2.0 );
      normal[1] = 1.0 / sqrt ( 2.0 );
      normal[2] = 0.0;
    }

    fprintf ( stdout, "\n" );
    fprintf ( stdout, "  Plane normal vector number %d\n", k );
    fprintf ( stdout, "\n" );
    for ( i = 0; i < 3; i++ )
    {
      fprintf ( stdout, "  %f", normal[i] );
    }
    fprintf ( stdout, "\n" );

    for ( l = 0; l <= 6; l++ )
    {
      for ( i = 0; i < 3; i++ )
      {
        pp[i] = normal[i] * ( double ) ( l ) / 5.0;
      }
      plane_normal_tetrahedron_intersect ( pp, normal, t, &int_num, pint );

      fprintf ( stdout, "\n" );
      fprintf ( stdout, "  Point on plane:\n" );
      fprintf ( stdout, "\n" );
      for ( i = 0; i < 3; i++ )
      {
        fprintf ( stdout, "  %f", pp[i] );
      }
      fprintf ( stdout, "\n" );
      fprintf ( stdout, "\n" );
      fprintf ( stdout, "  Number of intersection points = %d\n", int_num );
      fprintf ( stdout, "\n" );
      for ( j = 0; j < int_num; j++ )
      {
        fprintf ( stdout, "  %4d", j );
        for ( i = 0; i < 3; i++ )
        {
          fprintf ( stdout, "  %f", pint[i+j*3] );
        }
        fprintf ( stdout, "\n" );
      }
    }
  }
  return;
}
/******************************************************************************/

void test171 ( )

/******************************************************************************/
/*
  Purpose:

    TEST171 tests QUAD_AREA_2D and QUAD_AREA2_2D;

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    09 May 2010

  Author:

    John Burkardt
*/
{
# define DIM_NUM 2

  double area;
  double quad[DIM_NUM*4] = {
    0.0, 0.0,
    1.0, 0.0,
    1.0, 1.0,
    0.0, 1.0 };

  printf ( "\n" );
  printf ( "TEST171\n" );
  printf ( "  For a quadrilateral in 2D:\n" );
  printf ( "  QUAD_AREA_2D finds the area;\n" );
  printf ( "  QUAD_AREA2_2D finds the area;\n" );

  r8mat_transpose_print ( DIM_NUM, 4, quad, "  The vertices:" );

  area = quad_area_2d ( quad );

  printf ( "\n" );
  printf ( "  QUAD_AREA_2D area is  %f\n", area );

  area = quad_area2_2d ( quad );

  printf ( "  QUAD_AREA2_2D area is %f\n", area );

  return;
# undef DIM_NUM
}
/******************************************************************************/

void test1712 ( )

/******************************************************************************/
/*
  Purpose:

    TEST1712 tests QUAD_AREA_3D;

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    09 May 2010

  Author:

    John Burkardt
*/
{
  double area;
  double area1;
  double area2;
  int i;
  int j;
  double q[3*4] = {
    2.0, 2.0, 0.0, 
    0.0, 0.0, 0.0, 
    1.0, 1.0, 1.0, 
    3.0, 3.0, 1.0 };
  double t[3*3];

  printf ( "\n" );
  printf ( "TEST1712\n" );
  printf ( "  For a quadrilateral in 3D:\n" );
  printf ( "  QUAD_AREA_3D finds the area.\n" );

  r8mat_transpose_print ( 3, 4, q, "  The vertices:");

  area = quad_area_3d ( q );

  printf ( "\n" );
  printf ( "  QUAD_AREA_3D area is     %f\n", area );

  for ( j = 0; j < 3; j++ )
  {
    for ( i = 0; i < 3; i++ )
    {
      t[i+j*3] = q[i+j*3];
    }
  }
  area1 = triangle_area_3d ( t );
  for ( j = 0; j < 2; j++ )
  {
    for ( i = 0; i < 3; i++ )
    {
      t[i+j*3] = q[i+(j+2)*3];
    }
  }
  for ( i = 0; i < 3; i++ )
  {
    t[i+2*3] = q[i+0*3];
  }
  area2 = triangle_area_3d ( t );
  printf ( "  Sum of TRIANGLE_AREA_3D: %f\n", area1 + area2 );

  return;
}
/******************************************************************************/

void test1835 ( )

/******************************************************************************/
/*
  Purpose:

    TEST1835 tests SPHERE_EXP2IMP_3D and SPHERE_IMP2EXP_3D.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    26 July 2011

  Author:

    John Burkardt
*/
{
# define DIM_NUM 3

  int i;
  double pc[DIM_NUM] = {  1.0, 2.0, 3.0 };
  double p1[DIM_NUM] = {  4.0, 2.0, 3.0 };
  double p2[DIM_NUM] = {  1.0, 5.0, 3.0 };
  double p3[DIM_NUM] = {  1.0, 2.0, 6.0 };
  double p4[DIM_NUM] = { -2.0, 2.0, 3.0 };
  double r = 3.0;

  printf ( "\n" );
  printf ( "TEST1835\n" );
  printf ( "  SPHERE_EXP2IMP_3D: explicit sphere => implicit form;\n" );
  printf ( "  SPHERE_IMP2EXP_3D: implicit sphere => explicit form.\n" );

  printf ( "\n" );
  printf ( "  Initial form of explicit sphere:\n" );
  printf ( "\n" );
  printf ( "  P1:" );
  for ( i = 0; i < DIM_NUM; i++ )
  {
    printf ( "  %12f", p1[i] );
  }
  printf ( "\n" );
  printf ( "  P2:" );
  for ( i = 0; i < DIM_NUM; i++ )
  {
    printf ( "  %12f", p2[i] );
  }
  printf ( "\n" );
  printf ( "  P3:" );
  for ( i = 0; i < DIM_NUM; i++ )
  {
    printf ( "  %12f", p3[i] );
  }
  printf ( "\n" );
  printf ( "  P4:" );
  for ( i = 0; i < DIM_NUM; i++ )
  {
    printf ( "  %12f", p4[i] );
  }
  printf ( "\n" );

  sphere_exp2imp_3d ( p1, p2, p3, p4, &r, pc );

  printf ( "\n" );
  printf ( "  Computed form of implicit sphere:\n" );
  printf ( "\n" );
  printf ( "  Imputed radius = %f\n", r );

  r8vec_print ( DIM_NUM, pc, "  Imputed center:" );

  sphere_imp2exp_3d ( r, pc, p1, p2, p3, p4 );

  printf ( "\n" );
  printf ( "  Computed form of explicit sphere:\n" );
  printf ( "\n" );
  printf ( "  P1:" );
  for ( i = 0; i < DIM_NUM; i++ )
  {
    printf ( "  %12f", p1[i] );
  }
  printf ( "\n" );
  printf ( "  P2:" );
  for ( i = 0; i < DIM_NUM; i++ )
  {
    printf ( "  %12f", p2[i] );
  }
  printf ( "\n" );
  printf ( "  P3:" );
  for ( i = 0; i < DIM_NUM; i++ )
  {
    printf ( "  %12f", p3[i] );
  }
  printf ( "\n" );
  printf ( "  P4:" );
  for ( i = 0; i < DIM_NUM; i++ )
  {
    printf ( "  %12f", p4[i] );
  }
  printf ( "\n" );

  return;
# undef DIM_NUM
}
/******************************************************************************/

void test1836 ( )

/******************************************************************************/
/*
  Purpose:

    TEST1836 tests SPHERE_EXP2IMP_ND.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    26 July 2011

  Author:

    John Burkardt
*/
{
# define N 3

  int n = N;
  double p[N*(N+1)] = {
    4.0, 2.0, 3.0, 
    1.0, 5.0, 3.0, 
    1.0, 2.0, 6.0, 
   -2.0, 2.0, 3.0 };
  double pc[N];
  double pc_true[N] = { 1.0, 2.0, 3.0 };
  double r;
  double r_true = 3.0;

  printf ( "\n" );
  printf ( "TEST1836\n" );
  printf ( "  SPHERE_EXP2IMP_ND: explicit sphere => implicit form;\n" );

  r8mat_transpose_print ( n, n + 1, p, "  Initial form of explicit sphere:" );

  sphere_exp2imp_nd ( n, p, &r, pc );

  printf ( "\n" );
  printf ( "  Computed form of implicit sphere:\n" );
  printf ( "\n" );
  printf ( "  Imputed radius = %f\n", r );
  printf ( "  True radius =    %f\n", r_true );

  r8vec_print ( n, pc, "  Imputed center" );

  r8vec_print ( n, pc_true, "  True center" );

  return;
# undef N
}
/******************************************************************************/

void test2101 ( )

/******************************************************************************/
/*
  Purpose:

    TEST2101 tests TRIANGLE_CIRCUMCENTER_2D and others.

  Discussion:

    The functions tested include
    * TRIANGLE_CIRCUMCENTER_2D;
    * TRIANGLE_CIRCUMCENTER_2D_2;
    * TRIANGLE_CIRCUMCENTER.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    28 October 2010

  Author:

    John Burkardt
*/
{
# define M 2
# define TEST_NUM 4

  int i;
  int j;
  int m = M;
  double *pc;
  double t[M*3];
  double t_test[M*3*TEST_NUM] = {
         10.0,  5.0, 
         11.0,  5.0, 
         10.0,  6.0, 
         10.0,  5.0, 
         11.0,  5.0, 
         10.5,  5.86602539, 
         10.0,  5.0, 
         11.0,  5.0, 
         10.5, 15.0, 
         10.0,  5.0, 
         11.0,  5.0, 
        20.0,   7.0 };
  int test;
  int test_num = TEST_NUM;

  printf ( "\n" );
  printf ( "TEST2101\n" );
  printf ( "  For a triangle in 2D, the circumenter can be computed by:\n" );
  printf ( "  TRIANGLE_CIRCUMCENTER_2D;\n" );
  printf ( "  TRIANGLE_CIRCUMCENTER_2D_2;\n" );
  printf ( "  TRIANGLE_CIRCUMCENTER (any dimension);\n" );

  for ( test = 0; test < test_num; test++ )
  {
    for ( j = 0; j < 3; j++ )
    {
      for ( i = 0; i < m; i++ )
      {
        t[i+j*m] = t_test[i+j*m+test*m*3];
      }
    }
    r8mat_transpose_print ( m, 3, t, "  Triangle vertices:" );

    pc = triangle_circumcenter_2d ( t );
    r8vec_print ( m, pc, "  Circumcenter by TRIANGLE_CIRCUMCENTER_2D:" );
    free ( pc );

    pc = triangle_circumcenter_2d_2 ( t );
    r8vec_print ( m, pc, "  Circumcenter by TRIANGLE_CIRCUMCENTER_2D_2:" );
    free ( pc );

    pc = triangle_circumcenter ( m, t );
    r8vec_print ( m, pc, "  Circumcenter by TRIANGLE_CIRCUMCENTER:" );
    free ( pc );
  }
  return;
# undef M
# undef TEST_NUM
}
/******************************************************************************/

void test21011 ( )

/******************************************************************************/
/*
  Purpose:

    TEST21011 tests TRIANGLE_CIRCUMCENTER.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    28 October 2010

  Author:

    John Burkardt
*/
{
# define M1 2
# define TEST_NUM 4

  double *a12;
  int i;
  int j;
  int k;
  int m2;
  double *o1;
  double *o2;
  int m1 = M1;
  double pc1[M1];
  double *pc2;
  int seed;
  double t1[M1*3];
  double *t2;
  double t_test[M1*3*TEST_NUM] = {
         10.0,  5.0, 
         11.0,  5.0, 
         10.0,  6.0, 
         10.0,  5.0, 
         11.0,  5.0, 
         10.5,  5.86602539, 
         10.0,  5.0, 
         11.0,  5.0, 
         10.5, 15.0, 
         10.0,  5.0, 
         11.0,  5.0, 
        20.0,   7.0 };
  int test;
  int test_num = TEST_NUM;

  printf ( "\n" );
  printf ( "TEST21011\n" );
  printf ( "  For a triangle in M dimensions, the circumenter can be computed by:\n" );
  printf ( "  TRIANGLE_CIRCUMCENTER;\n" );
//
//  Vary the dimension.
//
  for ( m2 = 2; m2 <= 5; m2++ )
  {
    seed = 123456789;

    printf ( "\n" );
    printf ( "  M2 = %d\n", m2 );

    t2 = ( double * ) malloc ( m2 * 3 * sizeof ( double ) );
//
//  Randomly choose a mapping P2 = O2 + A12 * ( P1 - O1 )
//
    a12 = r8mat_uniform_01_new ( m2, m1, &seed );
    o1 = r8vec_uniform_01_new ( m1, &seed );
    o2 = r8vec_uniform_01_new ( m2, &seed );
//
//  Map each M1-dimensional triangle into M2 space.
//
    for ( test = 0; test < test_num; test++ )
    {
      for ( j = 0; j < 3; j++ )
      {
        for ( i = 0; i < m1; i++ )
        {
          t1[i+j*m1] = t_test[i+j*m1+test*m1*3];
        }
      }
      for ( j = 0; j < 3; j++ )
      {
        t1[i+j*m1] = t1[i+j*m1] - o1[i];
      }

      for ( j = 0; j < 3; j++ )
      {
        for ( i = 0; i < m2; i++ )
        {
          t2[i+j*m2] = 0.0;
          for ( k = 0; k < m1; k++ )
          {
            t2[i+j*m2] = t2[i+j*m2] + a12[i+k*m2] * t1[k+j*m1];
          }
        }
      }
      for ( j = 0; j < 3; j++ )
      {
        for ( i = 0; i < m2; i++ )
        {
          t2[i+j*m2] = t2[i+j*m2] + o2[i];
        }
      }

      pc2 = triangle_circumcenter ( m2, t2 );

      r8vec_print ( m2, pc2, "  Circumcenter by TRIANGLE_CIRCUMCENTER:" );
      printf ( "\n" );
      printf ( "  Distances from circumcenter to vertices:\n" );
      printf ( "\n" );
      for ( j = 0; j < 3; j++ )
      {
        printf ( "  %f\n", r8vec_norm_affine ( m2, pc2, t2+j*m2 ) );
      }
      free ( pc2 );
    }
    free ( a12 );
    free ( o1 );
    free ( o2 );
    free ( t2 );
  }
  return;
}
