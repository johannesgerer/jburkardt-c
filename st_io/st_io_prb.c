# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <string.h>

# include "st_io.h"

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for ST_IO_PRB.

  Discussion:

    ST_IO_PRB tests the ST_IO library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    25 January 2014

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "ST_IO_PRB:\n" );
  printf ( "  C version\n" );
  printf ( "  Test the ST_IO library.\n" );

  test01 ( );
  test02 ( );
  test03 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "ST_IO_PRB:\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void  test01 ( )

/******************************************************************************/
/*
  Purpose:

    TEST01 tests ST_WRITE.

  Discussion:

    The matrix is:

      11  12   0   0  15
      21  22   0   0   0
       0   0  33   0  35
       0   0   0  44   0
      51   0  53   0  55

    The index vectors are 1 based, and so have to be converted to
    0-base before being written.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    21 July 2014

  Author:

    John Burkardt
*/
{
  double ast[11] = {
    51.0, 12.0, 11.0, 33.0, 15.0, 
    53.0, 55.0, 22.0, 35.0, 44.0, 
    21.0 };
  int i_max;
  int i_min;
  int ist[11] = {
     5, 1, 1, 3, 1, 5, 5, 2, 3, 4, 2 };
  int j_max;
  int j_min;
  int jst[11] = {
     1, 2, 1, 3, 5, 3, 5, 2, 5, 4, 1 };
  int n = 5;
  int nst = 11;
  int m = 5;
  char output_filename[] = "a5by5.st";

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  ST_WRITE writes an ST file.\n" );

  i4vec_dec ( nst, ist );
  i4vec_dec ( nst, jst );

  i_min = i4vec_min ( nst, ist );
  i_max = i4vec_max ( nst, ist );
  j_min = i4vec_min ( nst, jst );
  j_max = i4vec_max ( nst, jst );

  st_header_print ( i_min, i_max, j_min, j_max, m, n, nst );

  st_print ( m, n, nst, ist, jst, ast, 
    "  Sparse Triplet (ST) data:" );

  st_write ( output_filename, m, n, nst, ist, jst, ast );

  printf ( "\n" );
  printf ( "  Wrote the matrix data to '%s'\n", output_filename );

  return;
}
/******************************************************************************/

void test02 ( )

/******************************************************************************/
/*
  Purpose:

    TEST02 tests ST_HEADER_READ, ST_DATA_READ.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    22 July 2014

  Author:

    John Burkardt
*/
{
  double *ast;
  char input_filename[] = "kershaw.st";
  int i_max;
  int i_min;
  int *ist;
  int j_max;
  int j_min;
  int *jst;
  int n;
  int m;
  int nst;

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  ST_HEADER_READ reads the header from an ST file.\n" );
  printf ( "  ST_DATA_READ reads the data from an ST file.\n" );
  printf ( "\n" );
  printf ( "  Read the data from '%s'\n", input_filename );

  st_header_read ( input_filename, &i_min, &i_max, &j_min, &j_max, &m, &n, &nst );

  st_header_print ( i_min, i_max, j_min, j_max, m, n, nst );

  ast = ( double * ) malloc ( nst * sizeof ( double ) );
  ist = ( int * ) malloc ( nst * sizeof ( int ) );
  jst = ( int * ) malloc ( nst * sizeof ( int ) );

  st_data_read ( input_filename, m, n, nst, ist, jst, ast );

  st_print ( m, n, nst, ist, jst, ast, 
    "  TEST02 matrix data read from file" );

  free ( ast );
  free ( ist );
  free ( jst );

  return;
}
/******************************************************************************/

void test03 ( )

/******************************************************************************/
/*
  Purpose:

    TEST03 tests ST_SORT_A.

  Discussion:

    The matrix is:

      11  12   0   0  15
      21  22   0   0   0
       0   0  33   0  35
       0   0   0  44   0
      51   0  53   0  55

    The index vectors are 1 based, and so have to be converted to
    0-base before being written.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    22 July 2014

  Author:

    John Burkardt
*/
{
  double ast[11] = {
    51.0, 12.0, 11.0, 33.0, 15.0, 
    53.0, 55.0, 22.0, 35.0, 44.0, 
    21.0 };
  int i_max;
  int i_min;
  int ist[11] = {
     5, 1, 1, 3, 1, 5, 5, 2, 3, 4, 2 };
  int j_max;
  int j_min;
  int jst[11] = {
     1, 2, 1, 3, 5, 3, 5, 2, 5, 4, 1 };
  int n = 5;
  int nst = 11;
  int m = 5;

  printf ( "\n" );
  printf ( "TEST03\n" );
  printf ( "  ST_SORT_A sorts an ST matrix by columns.\n" );

  i_min = i4vec_min ( nst, ist );
  i_max = i4vec_max ( nst, ist );
  j_min = i4vec_min ( nst, jst );
  j_max = i4vec_max ( nst, jst );

  st_header_print ( i_min, i_max, j_min, j_max, m, n, nst );

  st_print ( m, n, nst, ist, jst, ast, 
    "  Matrix data before sorting:" );

  st_sort_a ( m, n, nst, ist, jst, ast );

  st_print ( m, n, nst, ist, jst, ast, 
    "  Matrix data sorted by column:" );

  return;
}
