# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <string.h>

# include "cc_io.h"

int main ( );
void test01 ( );
void test02 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for CC_IO_PRB.

  Discussion:

    CC_IO_PRB tests the CC_IO library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 July 2014

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "CC_IO_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the CC_IO library.\n" );

  test01 ( );
  test02 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "CC_IO_PRB\n" );
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

    TEST01 tests CC_WRITE using a tiny matrix.

  Discussion:

    This test uses a trivial matrix whose full representation is:

          2  3  0  0  0
          3  0  4  0  6
      A = 0 -1 -3  2  0
          0  0  1  0  0
          0  4  2  0  1

    The 1-based CC representation is

      #  ICC  CCC  ACC
     --  ---  ---  ---
      1    1    1    2
      2    2         3

      3    1    3    3
      4    3        -1
      5    5         4

      6    2    6    4
      7    3        -3
      8    4         1
      9    5         2

     10    3   10    2

     11    2   11    6
     12    5         1

     13    *   13

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 July 2014

  Author:

    John Burkardt
*/
{
# define N 5
# define NCC 12

  double acc[NCC] = {
    2.0,  3.0,
    3.0, -1.0,  4.0,
    4.0, -3.0,  1.0, 2.0,
    2.0,
    6.0, 1.0 };
  int ccc[N+1] = {
    1, 3, 6, 10, 11, 13 };
  int icc[NCC] = {
    1, 2,
    1, 3, 5,
    2, 3, 4, 5,
    3,
    2, 5 };
  int m = N;
  int n = N;
  int ncc = NCC;
  char prefix[] = "simple";

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  Write a sparse matrix in CC format to 3 files.\n" );
/*
  Full storage statistics
*/
  printf ( "\n" );
  printf ( "  Full rows    M = %d\n", m );
  printf ( "  Full columns N = %d\n", n );
  printf ( "  Full storage   = %d\n", m * n );
/*
  Decrement the 1-based data.
*/
  i4vec_dec ( n + 1, ccc );
  i4vec_dec ( ncc, icc );
/*
  Print the CC matrix.
*/
  cc_print ( m, n, ncc, icc, ccc, acc, "  The matrix in 0-based CC format:" );
/*
  Write the matrix to 3 files.
*/
  cc_write ( prefix, ncc, n, icc, ccc, acc );

  return;
# undef NCC
}
/******************************************************************************/

void test02 ( )

/******************************************************************************/
/*
  Purpose:

    TEST02 tests CC_HEADER_READ and CC_DATA_READ.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 July 2014

  Author:

    John Burkardt
*/
{
  double *acc;
  int *ccc;
  int *icc;
  int m;
  int n;
  int ncc;
  char prefix[] = "simple";

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  Read a sparse matrix in CC format from 3 files.\n" );
/*
  Read the header.
*/
  cc_header_read ( prefix, &ncc, &n );
/*
  Allocate space.
*/
  acc = ( double * ) malloc ( ncc * sizeof ( double ) );
  ccc = ( int * ) malloc ( ( n + 1 ) * sizeof ( int ) );
  icc = ( int * ) malloc ( ncc * sizeof ( int ) );
/*
  Read the matrix data.
*/
  cc_data_read ( prefix, ncc, n, icc, ccc, acc );
/*
  Print the CC matrix.
*/
  m = n;
  cc_print ( m, n, ncc, icc, ccc, acc, "  The matrix in 0-based CC format:" );
/*
  Free memory.
*/
  free ( acc );
  free ( ccc );
  free ( icc );

  return;
}