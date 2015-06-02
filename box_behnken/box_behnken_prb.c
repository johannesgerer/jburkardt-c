# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "box_behnken.h"

int main ( );

void test01 ( );
void test02 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for BOX_BEHNKEN_PRB.

  Discussion:

    BOX_BEHNKEN_PRB tests the BOX_BEHNKEN library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 November 2013

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "BOX_BEHNKEN_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the BOX_BEHNKEN library.\n" );

  test01 ( );
  test02 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "BOX_BEHNKEN_PRB\n" );
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

    TEST01 tests BOX_BEHNKEN.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    25 October 2008

  Author:

    John Burkardt
*/
{
# define DIM_NUM 3

  int dim_num = DIM_NUM;
  double range[DIM_NUM*2] = {
    0.0, 10.0,  5.0,
    1.0, 11.0, 15.0 };
  int x_num;
  double *x;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  BOX_BEHNKEN computes a Box-Behnken dataset.\n" );

  r8mat_transpose_print ( dim_num, 2, range, "  The ranges:" );

  x_num = box_behnken_size ( dim_num );

  printf ( "\n" );
  printf ( "  For dimension DIM_NUM = %d\n", dim_num );
  printf ( "  the Box-Behnken design is of size %d\n", x_num );

  x = box_behnken ( dim_num, x_num, range );

  r8mat_transpose_print ( dim_num, x_num, x, "  The Box-Behnken design:" );

  free ( x );

  return;
# undef DIM_NUM
}
/******************************************************************************/

void test02 ( )

/******************************************************************************/
/*
  Purpose:

    TEST02 tests R8MAT_WRITE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 February 2012

  Author:

    John Burkardt
*/
{
# define DIM_NUM 4

  int dim_num = DIM_NUM;
  char file_out_name[] = "box_behnken_04_33.txt";
  double range[DIM_NUM*2] = {
    0.0, 0.0, 0.0, 0.0,
    1.0, 1.0, 1.0, 1.0 };
  int x_num;
  double *x;

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  R8MAT_WRITE writes a Box-Behnken dataset\n" );
  printf ( "  to a file.\n" );

  r8mat_transpose_print ( dim_num, 2, range, "  The ranges:" );

  x_num = box_behnken_size ( dim_num );

  printf ( "\n" );
  printf ( "  For dimension DIM_NUM = %d\n", dim_num );
  printf ( "  the Box-Behnken design is of size %d\n", x_num );

  x = box_behnken ( dim_num, x_num, range );

  r8mat_write ( file_out_name, dim_num, x_num, x );

  free ( x );

  printf ( "\n" );
  printf ( "  The data was written to the file \"%s\"\n", file_out_name  );

  return;
# undef DIM_NUM
}

