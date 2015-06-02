# include <stdlib.h>
# include <stdio.h>

# include "index.h"

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

    MAIN is the main program for INDEX_PRB.

  Discussion:

    INDEX_PRB tests the INDEX library.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    28 November 2012

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "INDEX_PRB:\n" );
  printf ( "  C version\n" );
  printf ( "  Test the INDEX library.\n" );

  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
  test05 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "INDEX_PRB:\n" );
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

    TEST01 tests INDEX0 and INDEX1.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    28 November 2012

  Author:

    John Burkardt
*/
{
  int i;
  int i_max;
  int i_min;
  int index_max;
  int index_min;
  int value;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  INDEX0 indexes a 1D array with zero base,\n" );
  printf ( "  INDEX1 indexes a 1D array with  unit base.\n" );
  printf ( "\n" );
  printf ( "             Min Index   Max\n" );
  printf ( "\n" );

  i_min = 1;
  i = 3;
  i_max = 5;
  printf ( "  1D Index  %4d  %4d  %4d\n", i_min,     i,     i_max );

  value = index0 ( i_min, i, i_max );
  index_min = 0;
  index_max = index_min + i_max - i_min;
  printf ( "  Index0    %4d  %4d  %4d\n",index_min, value, index_max );

  value = index1 ( i_min, i, i_max );
  index_min = 1;
  index_max = index_min + i_max - i_min; 
  printf ( "  Index1    %4d  %4d  %4d\n",index_min, value, index_max );

  return;
}
/******************************************************************************/

void test02 ( )

/******************************************************************************/
/*
  Purpose:

    TEST02 tests INDEX01, INDEX10, INDEX12 and INDEX21.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    28 November 2012

  Author:

    John Burkardt
*/
{
  int i;
  int i_max;
  int i_min;
  int index_max;
  int index_min;
  int j;
  int j_max;
  int j_min;
  int value;

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  For a 2D array,\n" );
  printf ( "  INDEX01 column indexes with zero base,\n" );
  printf ( "  INDEX10 row indexes with zero base,\n" );
  printf ( "  INDEX12 column indexes with unit base,\n" );
  printf ( "  INDEX21 row indexes with unit base.\n" );
  printf ( "\n" );
  printf ( "                Min   Index     Max\n" );
  printf ( "\n" );

  i_min = 1;
  i = 3;
  i_max = 5;
  j_min = 1;
  j = 2;
  j_max = 4;
  printf ( "  2D Index:  %3d%3d  %3d%3d  %3d%3d\n", i_min, j_min, i, j, i_max, j_max );

  value = index01 ( i_min, i, i_max, j_min, j, j_max );
  index_min = 0;
  index_max = index_min + ( i_max - i_min + 1 ) * ( j_max - j_min + 1 ) - 1;
  printf ( "  INDEX01:   %6d  %6d  %6d\n", index_min, value, index_max );

  value = index10 ( i_min, i, i_max, j_min, j, j_max );
  index_min = 0;
  index_max = index_min + ( i_max - i_min + 1 ) * ( j_max - j_min + 1 ) - 1;
  printf ( "  INDEX10:   %6d  %6d  %6d\n", index_min, value, index_max );

  value = index12 ( i_min, i, i_max, j_min, j, j_max );
  index_min = 1;
  index_max = index_min + ( i_max - i_min + 1 ) * ( j_max - j_min + 1 ) - 1;
  printf ( "  INDEX12:   %6d  %6d  %6d\n", index_min, value, index_max );

  value = index21 ( i_min, i, i_max, j_min, j, j_max );
  index_min = 1;
  index_max = index_min + ( i_max - i_min + 1 ) * ( j_max - j_min + 1 ) - 1;
  printf ( "  INDEX21:   %6d  %6d  %6d\n", index_min, value, index_max );

  return;
}
/******************************************************************************/

void test03 ( )

/******************************************************************************/
/*
  Purpose:

    TEST03 tests INDEX012, INDEX123, INDEX210, and INDEX321.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    28 November 2012

  Author:

    John Burkardt
*/
{
  int i;
  int i_max;
  int i_min;
  int index_max;
  int index_min;
  int j;
  int j_max;
  int j_min;
  int k;
  int k_max;
  int k_min;
  int m;
  int value;

  printf ( "\n" );
  printf ( "TEST03\n" );
  printf ( "  For a 3D array,\n" );
  printf ( "  INDEX012 column indexes with zero base,\n" );
  printf ( "  INDEX123 column indexes with unit base,\n" );
  printf ( "  INDEX210 row indexes with zero base.\n" );
  printf ( "  INDEX321 row indexes with unit base.\n" );
  printf ( "\n" );
  printf ( "                   Min      Index        Max\n" );
  printf ( "\n" );

  i_min = 1;
  i = 3;
  i_max = 5;
  j_min = 1;
  j = 2;
  j_max = 4;
  k_min = 1;
  k = 1;
  k_max = 3;

  m = ( i_max - i_min + 1 ) 
    * ( j_max - j_min + 1 ) 
    * ( k_max - k_min + 1 );

  printf ( "  3D Index:  %3d%3d%3d  %3d%3d%3d  %3d%3d%3d\n",
    i_min, j_min, k_min, i, j, k, i_max, j_max, k_max );

  value = index012 ( i_min, i, i_max, j_min, j, j_max, k_min, k, k_max );
  index_min = 0;
  index_max = index_min + m - 1;
  printf ( "  INDEX012:  %9d  %9d  %9d\n", index_min, value, index_max );

  value = index123 ( i_min, i, i_max, j_min, j, j_max, k_min, k, k_max );
  index_min = 1;
  index_max = index_min + m - 1;
  printf ( "  INDEX123:  %9d  %9d  %9d\n", index_min, value, index_max );

  value = index210 ( i_min, i, i_max, j_min, j, j_max, k_min, k, k_max );
  index_min = 0;
  index_max = index_min + m - 1;
  printf ( "  INDEX210:  %9d  %9d  %9d\n", index_min, value, index_max );

  value = index321 ( i_min, i, i_max, j_min, j, j_max, k_min, k, k_max );
  index_min = 1;
  index_max = index_min + m - 1;
  printf ( "  INDEX321:  %9d  %9d  %9d\n", index_min, value, index_max );

  return;
}
/******************************************************************************/

void test04 ( )

/******************************************************************************/
/*
  Purpose:

    TEST04 tests INDEX0123, INDEX1234, INDEX3210, and INDEX4321.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    28 November 2012

  Author:

    John Burkardt
*/
{
  int i;
  int i_max;
  int i_min;
  int index_max;
  int index_min;
  int j;
  int j_max;
  int j_min;
  int k;
  int k_max;
  int k_min;
  int l;
  int l_max;
  int l_min;
  int m;
  int value;

  printf ( "\n" );
  printf ( "TEST04\n" );
  printf ( "  For a 4D array,\n" );
  printf ( "  INDEX0123 column indexes with zero base,\n" );
  printf ( "  INDEX1234 column indexes with unit base,\n" );
  printf ( "  INDEX3210 row indexes with zero base,\n" );
  printf ( "  INDEX4321 row indexes with unit base.\n" );
  printf ( "\n" );
  printf ( "                       Min         Index           Max\n" );
  printf ( "\n" );

  i_min = 1;
  i = 3;
  i_max = 5;
  j_min = 1;
  j = 2;
  j_max = 4;
  k_min = 1;
  k = 1;
  k_max = 3;
  l_min = 1;
  l = 2;
  l_max = 2;

  m = ( i_max - i_min + 1 ) 
    * ( j_max - j_min + 1 ) 
    * ( k_max - k_min + 1 ) 
    * ( l_max - l_min + 1 );

  printf ( "  4D Index:    %3d%3d%3d%3d  %3d%3d%3d%3d  %3d%3d%3d%3d\n",
    i_min, j_min, k_min, l_min, i, j, k, l, i_max, j_max, k_max, l_max );

  value = index0123 ( i_min, i, i_max, j_min, j, j_max, k_min, k, k_max, l_min, l, l_max );
  index_min = 0;
  index_max = index_min + m - 1;
  printf ( "  INDEX0123:   %12d  %12d  %12d\n", index_min, value, index_max );

  value = index1234 ( i_min, i, i_max, j_min, j, j_max, k_min, k, k_max, l_min, l, l_max );
  index_min = 1;
  index_max = index_min + m - 1;
  printf ( "  INDEX1234:   %12d  %12d  %12d\n", index_min, value, index_max );

  value = index3210 ( i_min, i, i_max, j_min, j, j_max, k_min, k, k_max, l_min, l, l_max );
  index_min = 0;
  index_max = index_min + m - 1;
  printf ( "  INDEX3210:   %12d  %12d  %12d\n", index_min, value, index_max );

  value = index4321 ( i_min, i, i_max, j_min, j, j_max, k_min, k, k_max, l_min, l, l_max );
  index_min = 1;
  index_max = index_min + m - 1;
  printf ( "  INDEX4321:   %12d  %12d  %12d\n", index_min, value, index_max );

  return;
}
/******************************************************************************/

void test05 ( )

/******************************************************************************/
/*
  Purpose:

    TEST05 tests INDEX0N, INDEX1N, INDEXN0 and INDEXN1.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    28 November 2012

  Author:

    John Burkardt

*/
{
  int i[4] = {3,2,1,2};
  int i_max[4] = {5,4,3,2};
  int i_min[4] = {1,1,1,1};
  int index_max;
  int index_min;
  int j;
  int m;
  int n = 4;
  int value;

  printf ( "\n" );
  printf ( "TEST05\n" );
  printf ( "  For an N-dimensional array,\n" );
  printf ( "  INDEX0N column indexes with zero base,\n" );
  printf ( "  INDEX1N column indexes with unit base,\n" );
  printf ( "  INDEXN0 row indexes with zero base,\n" );
  printf ( "  INDEXN1 row indexes with unit base.\n" );
  printf ( "\n" );
  printf ( "                       Min         Index           Max\n" );

  m = 1;
  for ( j = 0; j < n; j++ )
  {
    m = m * ( i_max[j] - i_min[j] + 1 );
  }

  printf ( "  ND Index:   %3d%3d%3d%3d  %3d%3d%3d%3d  %3d%3d%3d%3d\n",
    i_min[0], i_min[1], i_min[2], i_min[3], 
    i[0], i[1], i[2], i[3], 
    i_max[0], i_max[1], i_max[2], i_max[3] );

  value = index0n ( n, i_min, i, i_max );
  index_min = 0;
  index_max = index_min + m - 1;
  printf ( "  INDEX0N:    %12d  %12d  %12d\n", index_min, value, index_max );

  value = index1n ( n, i_min, i, i_max );
  index_min = 1;
  index_max = index_min + m - 1;
  printf ( "  INDEX1N:    %12d  %12d  %12d\n", index_min, value, index_max );

  value = indexn0 ( n, i_min, i, i_max );
  index_min = 0;
  index_max = index_min + m - 1;
  printf ( "  INDEXN0:    %12d  %12d  %12d\n", index_min, value, index_max );

  value = indexn1 ( n, i_min, i, i_max );
  index_min = 1;
  index_max = index_min + m - 1;
  printf ( "  INDEXN1:    %12d  %12d  %12d\n", index_min, value, index_max );

  return;
}


