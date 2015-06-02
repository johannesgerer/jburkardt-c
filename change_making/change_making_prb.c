# include <stdlib.h>
# include <stdio.h>

# include "change_making.h"

int main ( );
void test01 ( );
void test02 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for CHANGE_MAKING_PRB.

  Discussion:

    CHANGE_MAKING_PRB tests the CHANGE_MAKING library.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    20 August 2014

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "CHANGE_MAKING_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the CHANGE_MAKING library.\n" );

  test01 ( );
  test02 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "CHANGE_MAKING_PRB\n" );
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

    CHANGE_MAKING_TEST01 lists the problem data.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 August 2014

  Author:

    John Burkardt
*/
{
  int coin_num;
  int coin_num_list[7] = {
    3, 
    5, 
    6, 
    7, 
    3, 
    6, 
    3 };
  int coin_offset;
  int coin_offset_list[7] = {
     0, 
     3, 
     8, 
    14, 
    21, 
    24, 
    30 };
  int *coin_value;
  int coin_value_list[33] = {
     5,  9, 13, 
     1,  4,  5,  8, 11, 
     1,  5, 10, 25, 50, 100, 
     1,  2,  6, 12, 24,  48,  60, 
     1,  3,  4, 
    16, 17, 23, 24, 39,  40, 
     6,  9, 20 };
  int coin_value_list_num = 33;
  int i;
  int target;
  int target_list[7] = {
    19, 
    29, 
    96, 
    96, 
     6, 
   100, 
    43 };
  int test;
  int test_num = 7;

  printf ( "\n" );
  printf ( "TEST01:\n" );
  printf ( "  List the problem data.\n" );

  for ( test = 0; test < test_num; test++ )
  {
    coin_num = coin_num_list[test];
    coin_offset = coin_offset_list[test];
    coin_value = ( int * ) malloc ( coin_num * sizeof ( int ) );
    for ( i = 0; i < coin_num; i++ )
    {
      coin_value[i] = coin_value_list[i+coin_offset];
    }
    target = target_list[test];

    printf ( "\n" );
    printf ( "  Test %d:\n", test );
    printf ( "  Number of coins = %d\n", coin_num );
    printf ( "  Values = \n" );
    for ( i = 0; i < coin_num; i++ )
    {
      printf ( "  %4d\n", coin_value[i] );
    }
    printf ( "\n" );
    printf ( "  Target = %d\n", target );

    free ( coin_value );
  }

  return;
}
/******************************************************************************/

void test02 ( )

/******************************************************************************/
/*
  Purpose:

    CHANGE_MAKING_TEST02 uses CHANGE_MAKING_LIST on the problems.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 August 2014

  Author:

    John Burkardt
*/
{
  int *a;
  int coin_num;
  int coin_num_list[7] = {
    3, 
    5, 
    6, 
    7, 
    3, 
    6, 
    3 };
  int coin_offset;
  int coin_offset_list[7] = {
     0, 
     3, 
     8, 
    14, 
    21, 
    24, 
    30 };
  int *coin_value;
  int coin_value_list[33] = {
     5,  9, 13, 
     1,  4,  5,  8, 11, 
     1,  5, 10, 25, 50, 100, 
     1,  2,  6, 12, 24,  48,  60, 
     1,  3,  4, 
    16, 17, 23, 24, 39,  40, 
     6,  9, 20 };
  int coin_value_list_num = 33;
  int i;
  int i4_huge = 2147483647;
  int target;
  int target_list[7] = {
    19, 
    29, 
    96, 
    96, 
     6, 
   100, 
    43 };
  int test;
  int test_num = 7;

  printf ( "\n" );
  printf ( "TEST02:\n" );
  printf ( "  CHANGE_MAKING LIST computes A(T), the smallest number\n" );
  printf ( "  of coins needed to form a given sum T, by computing\n" );
  printf ( "  the list A(0) through A(T).\n" );

  for ( test = 0; test < test_num; test++ )
  {
    coin_num = coin_num_list[test];
    coin_offset = coin_offset_list[test];
    coin_value = ( int * ) malloc ( coin_num * sizeof ( int ) );
    for ( i = 0; i < coin_num; i++ )
    {
      coin_value[i] = coin_value_list[i+coin_offset];
    }
    target = target_list[test];

    printf ( "\n" );
    printf ( "  Test %d:\n", test );
    printf ( "  Number of coins = %d\n", coin_num );
    printf ( "  Values = \n" );
    for ( i = 0; i < coin_num; i++ )
    {
      printf ( "  %4d\n", coin_value[i] );
    }
    printf ( "\n" );
    printf ( "  Target = %d\n", target );

    a = change_making_list ( coin_num, coin_value, target );

    printf ( "\n" );
    for ( i = 0; i <= target; i++ )
    {
      if ( a[i] == i4_huge )
      {
        printf ( "  %4d  Not possible\n", i );
      }
      else
      {
        printf ( "  %4d  %4d\n", i, a[i] );
      }
    }

    free ( a );
    free ( coin_value );
  }

  return;
}

