# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <string.h>

# include "cell.h"

int main ( );
void test01 ( );
void test02 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for CELL_PRB.

  Discussion:

    CELL_PRB tests the CELL library.

    An R8CVV is a "cell vector of vectors" of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    03 December 2012

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "CELL_PRB:\n" );
  printf ( "  C version\n" );
  printf ( "  Test the CELL library.\n" );

  test01 ( );
  test02 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "CELL_PRB:\n" );
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

    TEST01 stores some of Pascal's triangle in an R8CVV.

  Discussion:

    An R8CVV is a "cell array vector of vectors" of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 November 2012

  Author:

    John Burkardt
*/
{
  double *a;
  double *ai;
  double aij;
  int col;
  int i;
  int in[4] = { 0, 1, 4, 4 };
  int j;
  int jn[4] = { 1, 2, 3, 7 };
  int m = 5;
  int mn;
  int nn;
  int nr[5] = { 4, 5, 6, 7, 8 };
  int nr_max;
  int nv;
  int *roff;
  int row;
  double *vn;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  Use a real cell array (vector of vectors) to store rows 3:7\n" );
  printf ( "  of Pascal''s triangle.\n" );

  i4vec_print ( m, nr, "  The row sizes:" );
/*
  From the NR information:
  * determine the total size, MN
*/
  mn = r8cvv_size ( m, nr );
  printf ( "\n" );
  printf ( "  The storage for the cell array is %d\n", mn );
/*
  Allocate the cell array.
*/
  a = ( double * ) malloc ( mn * sizeof ( double ) );
/*
  Zero out the cell array.
*/
  for ( i = 0; i < mn; i++ )
  {
    a[i] = 0.0;
  }
/*
  Allocate a vector big enough to hold any single row.
*/
  nr_max = i4vec_max ( m, nr );
  ai = ( double * ) malloc ( nr_max * sizeof ( double ) );
/*
  From the NR information:
  * determine the offsets.
*/
  roff = r8cvv_offset ( m, nr );
  i4vec_print ( m + 1, roff, "  The row offsets:" );
/*
  Rows 1 through 5 of A will contain rows 3 through 7 of Pascal's triangle.
  Set these values one row at a time.
*/
  ai[0] = 1.0;

  for ( row = 1; row <= 7; row++ )
  {
    col = row + 1;
    ai[col-1] = ai[col-2];
    for ( col = row - 1; 1 <= col; col-- )
    {
      ai[col] = ai[col] + ai[col-1];
    }

    if ( 3 <= row )
    {
      i = row - 3;
      r8cvv_rset ( mn, a, m, roff, i, ai );
    }
  }
/*
  Print the cell array.
*/
  r8cvv_print ( mn, a, m, roff, "  Rows 3:7 of Pascal's Triangle:" );
/*
  Retrieve the entry from cell array row 2, column 3:
*/
  i = 2;
  j = 3;
  aij = r8cvv_iget ( mn, a, m, roff, i, j );
  printf ( "\n" );
  printf ( "  A(%d,%d) = %g\n", i, j, aij );
/*
  Retrieve row 3:
*/
  i = 3;
  ai = r8cvv_rget_new ( mn, a, m, roff, i );
  nv = roff[i+1] - roff[i];
  r8vec_transpose_print ( nv, ai, "  A(3,*):" );
/*
  Retrieve a list of entries.
*/
  nn = 4;
  vn = r8cvv_nget_new ( mn, a, m, roff, nn, in, jn );
  printf ( "\n" );
  printf ( "  Retrieve an arbitrary list of items:\n" );
  printf ( "\n" );
  for ( i = 0; i < nn; i++ )
  {
    printf ( "  A(%d,%d) = %g\n", in[i], jn[i], vn[i] );
  }
/*
  Free memory.
*/
  free ( a );
  free ( ai );
  free ( roff );
  free ( vn );

  return;
}
/******************************************************************************/

void test02 ( )

/******************************************************************************/
/*
  Purpose:

    TEST02 stores some of Pascal's triangle in an I4CVV.

  Discussion:

    An I4CVV is a "cell array vector of vectors" of I4's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    03 December 2012

  Author:

    John Burkardt
*/
{
  int *a;
  int *ai;
  int aij;
  int col;
  int i;
  int in[4] = { 0, 1, 4, 4 };
  int j;
  int jn[4] = { 1, 2, 3, 7 };
  int m = 5;
  int mn;
  int nn;
  int nr[5] = { 4, 5, 6, 7, 8 };
  int nr_max;
  int nv;
  int *roff;
  int row;
  int *vn;

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  Use an integer cell array (vector of vectors) to store rows 3:7\n" );
  printf ( "  of Pascal''s triangle.\n" );

  i4vec_print ( m, nr, "  The row sizes:" );
/*
  From the NR information:
  * determine the total size, MN
*/
  mn = i4cvv_size ( m, nr );
  printf ( "\n" );
  printf ( "  The storage for the cell array is %d\n", mn );
/*
  Allocate the cell array.
*/
  a = ( int * ) malloc ( mn * sizeof ( int ) );
/*
  Zero out the cell array.
*/
  for ( i = 0; i < mn; i++ )
  {
    a[i] = 0;
  }
/*
  Allocate a vector big enough to hold any single row.
*/
  nr_max = i4vec_max ( m, nr );
  ai = ( int * ) malloc ( nr_max * sizeof ( int ) );
/*
  From the NR information:
  * determine the offsets.
*/
  roff = i4cvv_offset ( m, nr );
  i4vec_print ( m + 1, roff, "  The row offsets:" );
/*
  Rows 1 through 5 of A will contain rows 3 through 7 of Pascal's triangle.
  Set these values one row at a time.
*/
  ai[0] = 1;

  for ( row = 1; row <= 7; row++ )
  {
    col = row + 1;
    ai[col-1] = ai[col-2];
    for ( col = row - 1; 1 <= col; col-- )
    {
      ai[col] = ai[col] + ai[col-1];
    }

    if ( 3 <= row )
    {
      i = row - 3;
      i4cvv_rset ( mn, a, m, roff, i, ai );
    }
  }
/*
  Print the cell array.
*/
  i4cvv_print ( mn, a, m, roff, "  Rows 3:7 of Pascal's Triangle:" );
/*
  Retrieve the entry from cell array row 2, column 3:
*/
  i = 2;
  j = 3;
  aij = i4cvv_iget ( mn, a, m, roff, i, j );
  printf ( "\n" );
  printf ( "  A(%d,%d) = %g\n", i, j, aij );
/*
  Retrieve row 3:
*/
  i = 3;
  ai = i4cvv_rget_new ( mn, a, m, roff, i );
  nv = roff[i+1] - roff[i];
  i4vec_transpose_print ( nv, ai, "  A(3,*):" );
/*
  Retrieve a list of entries.
*/
  nn = 4;
  vn = i4cvv_nget_new ( mn, a, m, roff, nn, in, jn );
  printf ( "\n" );
  printf ( "  Retrieve an arbitrary list of items:\n" );
  printf ( "\n" );
  for ( i = 0; i < nn; i++ )
  {
    printf ( "  A(%d,%d) = %d\n", in[i], jn[i], vn[i] );
  }
/*
  Free memory.
*/
  free ( a );
  free ( ai );
  free ( roff );
  free ( vn );

  return;
}
