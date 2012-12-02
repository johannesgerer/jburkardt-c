# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "pwl_interp_2d_scattered.h"
# include "test_interp_2d.h"
# include "r8lib.h"

int main ( );
void test01 ( );
void test02 ( );
void test03 ( int prob );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    PWL_INTERP_2D_SCATTERED_PRB tests PWL_INTERP_2D_SCATTERED_PRB.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    25 October 2012

  Author:

    John Burkardt
*/
{
  int prob;
  int prob_num;

  timestamp ( );
  printf ( "\n" );
  printf ( "PWL_INTERP_2D_SCATTERED_PRB:\n" );
  printf ( "  C version\n" );
  printf ( "  Test PWL_INTERP_2D_SCATTERED.\n" );
  printf ( "  The R8LIB library is needed.\n" );
  printf ( "  This test also needs the TEST_INTERP_2D library.\n" );

  test01 ( );
  test02 ( );
/*
  Numerical tests.
*/
  prob_num = f00_num ( );

  for ( prob = 1; prob <= prob_num; prob++ )
  {
    test03 ( prob );
  }
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "PWL_INTERP_2D_SCATTERED_PRB:\n" );
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

    TEST01 tests R8TRIS2.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    25 October 2012

  Author:

    John Burkardt
*/
{
  int element_neighbor[3*2*9];
  int element_num;
  int element_order = 3;
  int node_num = 9;
  double node_xy[2*9] = {
       0.0, 0.0, 
       0.0, 1.0, 
       0.2, 0.5, 
       0.3, 0.6, 
       0.4, 0.5, 
       0.6, 0.4, 
       0.6, 0.5, 
       1.0, 0.0, 
       1.0, 1.0 };
  int triangle[3*2*9];

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  R8TRIS2 computes the Delaunay triangulation of\n" );
  printf ( "  a set of nodes in 2D.\n" );
/*
  Set up the Delaunay triangulation.
*/
  r8tris2 ( node_num, node_xy, &element_num, triangle, element_neighbor );

  triangulation_order3_print ( node_num, element_num, node_xy, 
    triangle, element_neighbor );

  return;
}
/******************************************************************************/

void test02 ( )

/******************************************************************************/
/*
  Purpose:

    TEST02 tests PWL_INTERP_2D_SCATTERED_VALUE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    25 October 2012

  Author:

    John Burkardt
*/
{
  int element_neighbor[3*2*9];
  int element_num;
  int element_order = 3;
  int i;
  int j;
  int k;
  int ni = 25;
  int node_num = 9;
  double node_xy[2*9] = {
       0.0, 0.0, 
       0.0, 1.0, 
       0.2, 0.5, 
       0.3, 0.6, 
       0.4, 0.5, 
       0.6, 0.4, 
       0.6, 0.5,
       1.0, 0.0, 
       1.0, 1.0 };
  int triangle[3*2*9];
  double x;
  double xyi[2*25];
  double y;
  double zd[9];
  double ze;
  double *zi;

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  PWL_INTERP_2D_SCATTERED_VALUE evaluates a\n" );
  printf ( "  piecewise linear interpolant to scattered data.\n" );
/*
  Set up the Delaunay triangulation.
*/
  r8tris2 ( node_num, node_xy, &element_num, triangle, element_neighbor );

  for ( j = 0; j < element_num; j++ )
  {
    for ( i = 0; i < 3; i++ )
    {
      if ( 0 < element_neighbor[i+j*3] )
      {
        element_neighbor[i+j*3] = element_neighbor[i+j*3] - 1;
      }
    }
  }

  triangulation_order3_print ( node_num, element_num, node_xy, 
    triangle, element_neighbor );
/*
  Define the Z data.
*/
  for ( i = 0; i < node_num; i++ )
  {
    x = node_xy[0+i*2];
    y = node_xy[1+i*2];
    zd[i] = x + 2.0 * y;
  }
/*
  Define the interpolation points.
*/
  k = 0;
  for ( i = 0; i <= 4; i++ )
  {
    for ( j = 0; j <= 4; j++ )
    {
      xyi[0+k*2] = ( i - 1 ) / 4.0;
      xyi[1+k*2] = ( j - 1 ) / 4.0;
      k = k + 1;
    }
  }
/*
  Evaluate the interpolant.
*/
  zi = pwl_interp_2d_scattered_value ( node_num, node_xy, zd, element_num, 
    triangle, element_neighbor, ni, xyi );

  printf ( "\n" );
  printf ( "     K      Xi(K)       Yi(K)       Zi(K)       Z(X,Y)\n" );
  printf ( "\n" );
  for ( k = 0; k < ni; k++ )
  {
    ze = xyi[0+k*2] + 2.0 * xyi[1+k*2];
    printf ( "  %4d  %10.4f  %10.4f  %10.4f  %10.4f\n",
      k, xyi[0+k*2], xyi[1+k*2], zi[k], ze );
  }

  free ( zi );

  return;
}
/******************************************************************************/

void test03 ( int prob )

/******************************************************************************/
/*
  Purpose:

    TEST03 tests PWL_INTERP_2D_SCATTERED_VALUE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    25 October 2012

  Author:

    John Burkardt
*/
{
  int *element_neighbor;
  int element_num;
  int g;
  int i;
  int j;
  int k;
  int nd;
  int ni = 25;
  double rms;
  int *triangle;
  double x;
  double *xd;
  double xi[25];
  double *xyd;
  double xyi[2*25];
  double y;
  double *yd;
  double yi[25];
  double *zd;
  double *ze;
  double *zi;

  g = 2;
  nd = g00_size ( g );

  printf ( "\n" );
  printf ( "TEST03\n" );
  printf ( "  PWL_INTERP_2D_SCATTERED_VALUE evaluates a\n" );
  printf ( "  piecewise linear interpolant to scattered data.\n" );
  printf ( "  Here, we use grid number %d\n", g );
  printf ( "  with %d scattered points in the unit square\n", nd );
  printf ( "  on problem %d\n", prob );
/*
  Get the data points and evaluate the function there.
*/
  xd = ( double * ) malloc ( nd * sizeof ( double ) );
  yd = ( double * ) malloc ( nd * sizeof ( double ) );

  g00_xy ( g, nd, xd, yd );

  zd = ( double * ) malloc ( nd * sizeof ( double ) );
  f00_f0 ( prob, nd, xd, yd, zd );

  xyd = ( double * ) malloc ( 2 * nd * sizeof ( double ) );

  for ( i = 0; i < nd; i++ )
  {
    xyd[0+i*2] = xd[i];
    xyd[1+i*2] = yd[i];
  }
/*
  Set up the Delaunay triangulation.
*/
  element_neighbor = ( int * ) malloc ( 3 * 2 * nd * sizeof ( int ) );
  triangle = ( int * ) malloc ( 3 * 2 * nd * sizeof ( int ) );

  r8tris2 ( nd, xyd, &element_num, triangle, element_neighbor );

  for ( j = 0; j < element_num; j++ )
  {
    for ( i = 0; i < 3; i++ )
    {
      if ( 0 < element_neighbor[i+j*3] )
      {
        element_neighbor[i+j*3] = element_neighbor[i+j*3] - 1;
      }
    }
  }
/*
  Define the interpolation points.
*/
  k = 0;
  for ( i = 1; i <= 5; i++ )
  {
    for ( j = 1; j <= 5; j++ )
    {
      xyi[0+k*2] = ( 2 * i - 1 ) / 10.0;
      xyi[1+k*2] = ( 2 * j - 1 ) / 10.0;
      k = k + 1;
    }
  }

  for ( k = 0; k < ni; k++ )
  {
    xi[k] = xyi[0+k*2];
    yi[k] = xyi[1+k*2];
  }
  ze = ( double * ) malloc ( ni * sizeof ( double ) );
  f00_f0 ( prob, ni, xi, yi, ze );
/*
  Evaluate the interpolant.
*/
  zi = pwl_interp_2d_scattered_value ( nd, xyd, zd, element_num, 
    triangle, element_neighbor, ni, xyi );

  rms = 0.0;
  for ( k = 0; k < ni; k++ )
  {
    rms = rms + pow ( zi[k] - ze[k], 2 );
  }
  rms = sqrt ( rms / ( double ) ( ni ) );

  printf ( "\n" );
  printf ( "  RMS error is %g\n", rms );

  printf ( "\n" );
  printf ( "     K      Xi(K)       Yi(K)       Zi(K)       Z(X,Y)\n" );
  printf ( "\n" );

  for ( k = 0; k < ni; k++ )
  {
    printf ( "  %4d  %10.4f  %10.4f  %10.4f  %10.4f\n",
      k, xyi[0+k*2], xyi[1+k*2], zi[k], ze[k] );
  }

  free ( element_neighbor );
  free ( triangle );
  free ( xd );
  free ( xyd );
  free ( yd );
  free ( zd );
  free ( ze );

  return;
}
