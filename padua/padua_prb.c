# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <string.h>

# include "padua.h"

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );
void test04 ( );
void test05 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for PADUA_PRB.

  Discussion:

    PADUA_PRB tests the PADUA library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 June 2014

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "PADUA_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the PADUA library.\n" );

  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
  test05 ( );
/*
  Terminate.
*/
  printf ( " \n" );
  printf ( "PADUA_PRB\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( " \n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void test01 ( )

/******************************************************************************/
/*
  Purpose:

    TEST01 tests PADUA_ORDER.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 April 2014

  Author:

    John Burkardt
*/
{
  int l;
  int n;

  printf ( " \n" );
  printf ( "TEST01\n" );
  printf ( "  PADUA_ORDER converts the level L into the order N\n" );
  printf ( "  of any Padua rule.\n" );
  printf ( " \n" );
  printf ( "     L         N\n" );
  printf ( " \n" );

  for ( l = 0; l <= 10; l++ )
  {
    n = padua_order ( l );
    printf ( "  %4d  %8d\n", l, n );
  }

  return;
}
/******************************************************************************/

void test02 ( )

/******************************************************************************/
/*
  Purpose:

    TEST02 tests PADUA_POINTS.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 June 2014

  Author:

    John Burkardt
*/
{
  int l;
  char label[255];
  int n;
  double *xy;

  printf ( " \n" );
  printf ( "TEST02\n" );
  printf ( "  PADUA_POINTS returns the points of a Padua rule.\n" );

  for ( l = 0; l <= 10; l++ )
  {
    n = padua_order ( l );
    xy = padua_points ( l );
    sprintf ( label, "  Level %d Padua points:", l );
    r8mat_transpose_print ( 2, n, xy, label );
    free ( xy );
  }

  return;
}
/******************************************************************************/

void test03 ( )

/******************************************************************************/
/*
  Purpose:

    TEST03 tests PADUA_PLOT.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 June 2014

  Author:

    John Burkardt
*/
{
  char filename[255];
  int l;
 
  printf ( " \n" );
  printf ( "TEST03\n" );
  printf ( "  PADUA_PLOT plots the Padua points.\n" );

  strcpy ( filename, "padua_00" );

  for ( l = 0; l <= 10; l++ )
  {
    padua_plot ( l, filename );
    filename_inc ( filename );
  }

  return;
}
/******************************************************************************/

void test04 ( )

/******************************************************************************/
/*
  Purpose:

    TEST04 tests PADUA_POINTS and PADUA_POINTS_SET.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 June 2014

  Author:

    John Burkardt
*/
{
  int j;
  int l;
  int n;
  double *xy1;
  double *x2;
  double *y2;

  printf ( " \n" );
  printf ( "TEST04\n" );
  printf ( "  PADUA_POINTS computes the points of a Padua rule.\n" );
  printf ( "  PADUA_POINTS_SET looks them up in a table.\n" );

  for ( l = 3; l <= 4; l++ )
  {
    n = padua_order ( l );
    xy1 = padua_points ( l );
    x2 = ( double * ) malloc ( n * sizeof ( double ) );
    y2 = ( double * ) malloc ( n * sizeof ( double ) );
    padua_points_set ( l, x2, y2 );
    printf ( "\n" );
    printf ( "  Level %d Padua points\n", l );
    printf ( "\n" );
    for ( j = 0; j < n; j++ )
    {
      printf ( "  %4d  %14.6g  %14.6g\n", j, xy1[0+j*2], xy1[1+j*2] );
      printf ( "        %14.6g  %14.6g\n",   x2[j], y2[j] );
    }
    free ( xy1 );
    free ( x2 );
    free ( y2 );
  }

  return;
}
/******************************************************************************/

void test05 ( )

/******************************************************************************/
/*
  Purpose:

    TEST05 tests PADUA_WEIGHTS and PADUA_WEIGHTS_SET.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 June 2014

  Author:

    John Burkardt
*/
{
  double diff;
  int j;
  int l;
  int n;
  double *w1;
  double *w2;

  printf ( " \n" );
  printf ( "TEST05\n" );
  printf ( "  PADUA_WEIGHTS computes quadrature weights of a Padua rule.\n" );
  printf ( "  PADUA_WEIGHTS_SET looks them up in a table.\n" );

  for ( l = 3; l <= 4; l++ )
  {
    n = padua_order ( l );
    w1 = padua_weights ( l );
    w2 = padua_weights_set ( l );
    printf ( "\n" );
    printf ( "  Level %d Padua weights\n", l );
    printf ( "              Computed          Lookup\n" );
    printf ( "\n" );
    diff = 0.0;
    for ( j = 0; j < n; j++ )
    {
      diff = r8_max ( diff, fabs ( w1[j] - w2[j] ) );
      printf ( "  %4d  %14.6g  %14.6g\n", j, w1[j], w2[j] );
    }
    printf ( "\n" );
    printf ( "  Maximum difference = %g\n", diff );
    free ( w1 );
    free ( w2 );
  }

  return;
}
