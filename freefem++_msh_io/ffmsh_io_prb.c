# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <string.h>

# include "ffmsh_io.h"

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for FFMSH_IO_PRB.

  Discussion:

    FFMSH_IO_PRB tests the FFMSH_IO library.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    21 December 2014

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "FFMSH_IO_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the FFMSH_IO library.\n" );

  test01 ( );
  test02 ( );
  test03 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "FFMSH_IO_PRB\n" );
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

    TEST01 gets the example data and prints it.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    21 December 2014

  Author:

    John Burkardt
*/
{
  int *e_l;
  int e_num;
  int *e_v;
  int *t_l;
  int t_num;
  int *t_v;
  int *v_l;
  int v_num;
  double *v_xy;

  printf ( "\n" );
  printf ( "TEST01:\n" );
  printf ( "  Get the example 2D data and print it.\n" );
/*
  Get example sizes.
*/
  ffmsh_2d_size_example ( &v_num, &e_num, &t_num );
/*
  Print example sizes.
*/
  ffmsh_2d_size_print ( "  Example Sizes:", v_num, e_num, t_num );
/*
  Allocate memory.
*/
  v_xy = ( double * ) malloc ( 2 * v_num * sizeof ( double ) );
  v_l = ( int * ) malloc ( v_num * sizeof ( int ) );
  e_v = ( int * ) malloc ( 2 * e_num * sizeof ( int ) );
  e_l = ( int * ) malloc ( e_num * sizeof ( int ) );
  t_v = ( int * ) malloc ( 3 * t_num * sizeof ( int ) );
  t_l = ( int * ) malloc ( t_num * sizeof ( int ) );
/*
  Get example data.
*/
  ffmsh_2d_data_example ( v_num, e_num, t_num, v_xy, v_l, e_v, e_l, 
    t_v, t_l );
/*
  Print example data.
*/
  ffmsh_2d_data_print ( "  Example data:", v_num, e_num, t_num, v_xy, 
    v_l, e_v, e_l, t_v, t_l );
/*
  Free memory.
*/
  free ( e_l );
  free ( e_v );
  free ( t_l );
  free ( t_v );
  free ( v_l );
  free ( v_xy );

  return;
}
/******************************************************************************/

void test02 ( )

/******************************************************************************/
/*
  Purpose:

    TEST02 writes the example data to a file.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    22 December 2014

  Author:

    John Burkardt
*/
{
  int *e_l;
  int e_num;
  int *e_v;
  char filename[] = "output.msh";
  int *t_l;
  int t_num;
  int *t_v;
  int *v_l;
  int v_num;
  double *v_xy;

  printf ( "\n" );
  printf ( "TEST02:\n" );
  printf ( "  Get the example 2D data and print it.\n" );
/*
  Get example sizes.
*/
  ffmsh_2d_size_example ( &v_num, &e_num, &t_num );
/*
  Allocate memory.
*/
  v_xy = ( double * ) malloc ( 2 * v_num * sizeof ( double ) );
  v_l = ( int * ) malloc ( v_num * sizeof ( int ) );
  e_v = ( int * ) malloc ( 2 * e_num * sizeof ( int ) );
  e_l = ( int * ) malloc ( e_num * sizeof ( int ) );
  t_v = ( int * ) malloc ( 3 * t_num * sizeof ( int ) );
  t_l = ( int * ) malloc ( t_num * sizeof ( int ) );
/*
  Get example data.
*/
  ffmsh_2d_data_example ( v_num, e_num, t_num, v_xy, v_l, e_v, e_l, 
    t_v, t_l );
/*
  Write the data to a file.
*/
  ffmsh_2d_write ( filename, v_num, e_num, t_num, v_xy, 
    v_l, e_v, e_l, t_v, t_l );

  printf ( "\n" );
  printf ( "  The data was written to \"%s\"\n", filename );
/*
  Free memory.
*/
  free ( e_l );
  free ( e_v );
  free ( t_l );
  free ( t_v );
  free ( v_l );
  free ( v_xy );

  return;
}
/******************************************************************************/

void test03 ( )

/******************************************************************************/
/*
  Purpose:

    TEST03 gets the example data from a file.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    22 December 2014

  Author:

    John Burkardt
*/
{
  int *e_l;
  int e_num;
  int *e_v;
  char filename[] = "input.msh";
  int *t_l;
  int t_num;
  int *t_v;
  int *v_l;
  int v_num;
  double *v_xy;

  printf ( "\n" );
  printf ( "TEST03:\n" );
  printf ( "  Read the example 2D data from a file.\n" );
/*
  Read sizes.
*/
  ffmsh_2d_size_read ( filename, &v_num, &e_num, &t_num );
/*
  Print sizes.
*/
  ffmsh_2d_size_print ( "  Example Sizes:", v_num, e_num, t_num );
/*
  Allocate memory.
*/
  v_xy = ( double * ) malloc ( 2 * v_num * sizeof ( double ) );
  v_l = ( int * ) malloc ( v_num * sizeof ( int ) );
  e_v = ( int * ) malloc ( 2 * e_num * sizeof ( int ) );
  e_l = ( int * ) malloc ( e_num * sizeof ( int ) );
  t_v = ( int * ) malloc ( 3 * t_num * sizeof ( int ) );
  t_l = ( int * ) malloc ( t_num * sizeof ( int ) );
/*
  Read data.
*/
  ffmsh_2d_data_read ( filename, v_num, e_num, t_num, v_xy, v_l, e_v, e_l, 
    t_v, t_l );
/*
  Print data.
*/
  ffmsh_2d_data_print ( "  Example data:", v_num, e_num, t_num, v_xy, 
    v_l, e_v, e_l, t_v, t_l );
/*
  Free memory.
*/
  free ( e_l );
  free ( e_v );
  free ( t_l );
  free ( t_v );
  free ( v_l );
  free ( v_xy );

  return;
}
