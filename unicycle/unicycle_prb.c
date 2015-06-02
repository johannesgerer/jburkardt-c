# include <stdlib.h>
# include <stdio.h>

# include "unicycle.h"

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );
void test04 ( );
void test05 ( );
void test06 ( );
void test07 ( );
void test08 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for UNICYCLE_PRB.

  Discussion:

    UNICYCLE_PRB tests the UNICYCLE library.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 June 2012

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "UNICYCLE_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the UNICYCLE library.\n" );

  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
  test05 ( );
  test06 ( );
  test07 ( );
  test08 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "UNICYCLE_PRB\n" );
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

    TEST01 tests PERM_IS_UNICYCLE.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 June 2012

  Author:

    John Burkardt
*/
{
  int n = 5;
  int *p;
  int seed;
  int test;
  int test_num = 10;
  int *u;
  int value;

  seed = 123456789;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  PERM_IS_UNICYCLE determines whether a permutation\n" );
  printf ( "  is a unicyle\n" );

  for ( test = 1; test <= test_num; test++ )
  {
    p = perm_random ( n, &seed );

    value = perm_is_unicycle ( n, p );

    if ( value )
    {
      perm_print ( n, p, "  This permutation is a unicycle" );
      u = unicycle_index_to_sequence ( n, p );
      unicycle_print ( n, u, "  The permutation in sequence form" );
      free ( u );
    }
    else
    {
      perm_print ( n, p, "  This permutation is NOT a unicycle" );
    }
  }
  return;
}
/******************************************************************************/

void test02 ( )

/******************************************************************************/
/*
  Purpose:

    TEST02 tests UNICYCLE_ENUM.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 June 2012

  Author:

    John Burkardt
*/
{
  int n;
  int n_max = 10;
  int num;

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  UNICYCLE_ENUM enumerates the unicycles of N objects.\n" );
  printf ( "\n" );
  printf ( "  N    Number\n" );
  printf ( "\n" );

  for ( n = 0; n <= n_max; n++ )
  {
    num = unicycle_enum ( n );
    printf ( "  %3d  %8d\n", n, num );
  }
  return;
}
/******************************************************************************/

void test03 ( )

/******************************************************************************/
/*
  Purpose:

    TEST03 tests UNICYCLE_INVERSE;

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    13 January 2012

  Author:

    John Burkardt
*/
{
  int n = 7;
  int u[7] = { 1, 7, 6, 2, 4, 3, 5 };
  int *u_inverse;

  printf ( "\n" );
  printf ( "TEST03\n" );
  printf ( "  UNICYCLE_INVERSE inverts a unicycle;\n" );

  unicycle_print ( n, u, "  The original unicycle:" );
 
  u_inverse = unicycle_inverse ( n, u );
 
  unicycle_print ( n, u_inverse, "  The inverse unicycle:" );
 
  free ( u_inverse );

  return;
}
/******************************************************************************/

void test04 ( )

/******************************************************************************/
/*
  Purpose:

    TEST04 tests UNICYCLE_NEXT.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 June 2012

  Author:

    John Burkardt
*/
{
  int i;
  int n = 5;
  int rank;
  int *u;

  u = ( int * ) malloc ( n * sizeof ( int ) );

  printf ( "\n" );
  printf ( "TEST04\n" );
  printf ( "  UNICYCLE_NEXT generates unicycles in lex order.\n" );
  printf ( "\n" );
  rank = -1;
 
  for ( ; ; )
  {
    unicycle_next ( n, u, &rank );

    if ( rank == -1 )
    {
      break;
    }

    printf ( "  %3d:", rank );
    for ( i = 0; i < n; i++ )
    {
      printf ( "  %2d", u[i] );
    }
    printf ( "\n" );
  }

  free ( u );

  return;
}
/******************************************************************************/

void test05 ( )

/******************************************************************************/
/*
  Purpose:

    TEST05 tests UNICYCLE_RANDOM;

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 June 2012

  Author:

    John Burkardt
*/
{
  int i;
  int n = 5;
  int seed;
  int *u;

  printf ( "\n" );
  printf ( "TEST05\n" );
  printf ( "  UNICYCLE_RANDOM produces a random unicyle\n" );;
  printf ( "  For this test, N = %d\n", n );
  printf ( "\n" );

  seed = 123456789;

  for ( i = 1; i <= 5; i++ )
  {
    u = unicycle_random ( n, &seed );
    unicycle_print ( n, u, " " );
    free ( u );
  }
  return;
}
/******************************************************************************/

void test06 ( )

/******************************************************************************/
/*
  Purpose:

    TEST06 tests UNICYCLE_RANK.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 June 2012

  Author:

    John Burkardt
*/
{
  int n = 5;
  int rank;
  int u[5] = { 1, 5, 2, 3, 4 };

  printf ( "\n" );
  printf ( "TEST06\n" );
  printf ( "  UNICYCLE_RANK ranks a unicycle.\n" );

  unicycle_print ( n, u, "  The unicycle:" );
 
  rank = unicycle_rank ( n, u );
 
  printf ( "\n" );
  printf ( "  The rank is: %d\n", rank );
 
  return;
}
/******************************************************************************/

void test07 ( )

/******************************************************************************/
/*
  Purpose:

    TEST07 tests UNICYCLE_UNRANK.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 June 2012

  Author:

    John Burkardt
*/
{
  int n = 5;
  int rank;
  int *u;

  printf ( "\n" );
  printf ( "TEST07\n" );
  printf ( "  UNICYCLE_UNRANK, given a rank, computes the\n" );
  printf ( "  corresponding unicycle.\n" );
  printf ( "\n" );
  rank = 6;
  printf ( "  The requested rank is %d\n", rank );
 
  u = unicycle_unrank ( n, rank );
 
  unicycle_print ( n, u, "  The unicycle:" );
 
  free ( u );

  return;
}
/******************************************************************************/

void test08 ( )

/******************************************************************************/
/*
  Purpose:

    TEST08 tests UNICYCLE_INDEX, UNICYCLE_INDEX_TO_SEQUENCE.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 June 2012

  Author:

    John Burkardt
*/
{
  int n = 6;
  int *u;
  int *u_index;
  int *u2;
  int seed;
  int test;
  int test_num = 5;

  seed = 123456789;

  printf ( "\n" );
  printf ( "TEST08\n" );
  printf ( "  UNICYCLE_INDEX converts a unicycle to index form.\n" );
  printf ( "  UNICYCLE_INDEX_TO_SEQUENCE converts an index to unicycle form.\n" );

  for ( test = 1; test <= test_num; test++ )
  {
    u = unicycle_random ( n, &seed );

    unicycle_print ( n, u, "  The unicycle:" );

    u_index = unicycle_index ( n, u );
    
    unicycle_index_print ( n, u_index, "  The index form:" );

    u2 = unicycle_index_to_sequence ( n, u_index );

    unicycle_print ( n, u2, "  The unicycle recovered:" );

    free ( u );
    free ( u_index );
    free ( u2 );
  }
  return;
}
