# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <string.h>

# include "cc_to_st.h"

int main ( );
void test01 ( );
void test02 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for CC_TO_ST_PRB.

  Discussion:

    CC_TO_ST_PRB tests the CC_TO_ST library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    23 July 2014

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "CC_TO_ST_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the CC_TO_ST library.\n" );

  test01 ( );
  test02 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "CC_TO_ST_PRB\n" );
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

    TEST01 tests CC_TO_ST using a 1-based matrix.

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

    18 July 2014

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
  double *ast;
  int ccc[N+1] = {
    1, 3, 6, 10, 11, 13 };
  int i;
  int icc[NCC] = {
    1, 2, 
    1, 3, 5, 
    2, 3, 4, 5, 
    3, 
    2, 5 };
  int *ist;
  int *jst;
  int m = 5;
  int n = N;
  int ncc = NCC;
  int nst;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  Convert a 1-based CC matrix to ST format.\n" );
/*
  Print the CC matrix.
*/
  cc_print ( m, n, ncc, icc, ccc, acc, "  The CC matrix:" );
/*
  Convert it.
*/
  ist = ( int * ) malloc ( ncc * sizeof ( int ) );
  jst = ( int * ) malloc ( ncc * sizeof ( int ) );
  ast = ( double * ) malloc ( ncc * sizeof ( double ) );

  cc_to_st ( m, n, ncc, icc, ccc, acc, &nst, ist, jst, ast );
/*
  Print the ST matrix.
*/
  st_print ( m, n, nst, ist, jst, ast, "  The ST matrix:" );
/*
  Free memory.
*/
  free ( ast );
  free ( ist );
  free ( jst );

  return;
# undef N
# undef NCC
}
/******************************************************************************/

void test02 ( )

/******************************************************************************/
/*
  Purpose:

    TEST02 tests CC_TO_ST using a 0-based matrix.

  Discussion:

    This test uses a trivial matrix whose full representation is:

          2  3  0  0  0
          3  0  4  0  6
      A = 0 -1 -3  2  0
          0  0  1  0  0
          0  4  2  0  1

    The 0-based CC representation is

      #  ICC  CCC  ACC
     --  ---  ---  ---
      0    0    0    2
      1    1         3

      2    0    2    3
      3    2        -1
      4    4         4

      5    1    5    4
      6    2        -3
      7    3         1
      8    4         2

      9    2    9    2

     10    1   10    6
     11    4         1

     12    *   12

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    23 July 2014

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
  double *ast;
  int ccc[N+1] = {
    0, 2, 5, 9, 10, 12 };
  int i;
  int icc[NCC] = {
    0, 1, 
    0, 2, 4, 
    1, 2, 3, 4, 
    2, 
    1, 4 };
  int *ist;
  int *jst;
  int m = 5;
  int n = N;
  int ncc = NCC;
  int nst;

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  Convert a 0-based CC matrix to ST format.\n" );
/*
  Print the CC matrix.
*/
  cc_print ( m, n, ncc, icc, ccc, acc, "  The CC matrix:" );
/*
  Convert it.
*/
  ist = ( int * ) malloc ( ncc * sizeof ( int ) );
  jst = ( int * ) malloc ( ncc * sizeof ( int ) );
  ast = ( double * ) malloc ( ncc * sizeof ( double ) );

  cc_to_st ( m, n, ncc, icc, ccc, acc, &nst, ist, jst, ast );
/*
  Print the ST matrix.
*/
  st_print ( m, n, nst, ist, jst, ast, "  The ST matrix:" );
/*
  Free memory.
*/
  free ( ast );
  free ( ist );
  free ( jst );

  return;
}
