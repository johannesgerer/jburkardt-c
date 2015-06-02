# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "bvec.h"

int main ( );

void bvec_add_test ( );
void bvec_complement2_test ( );
void bvec_mul_test ( );
void bvec_next_test ( );
void bvec_next_grlex_test ( );
void bvec_print_test ( );
void bvec_sub_test ( );
void bvec_to_i4_test ( );
void bvec_uniform_new_test ( );
void i4_bclr_test ( );
void i4_bset_test ( );
void i4_btest_test ( );
void i4_to_bvec_test ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for BVEC_PRB.

  Discussion:

    BVEC_PRB tests the BVEC library.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    13 March 2015

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "BVEC_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the BVEC library.\n" );

  bvec_add_test ( );
  bvec_complement2_test ( );
  bvec_mul_test ( );
  bvec_next_test ( );
  bvec_next_grlex_test ( );
  bvec_print_test ( );
  bvec_sub_test ( );
  bvec_to_i4_test ( );
  bvec_uniform_new_test ( );
  i4_bclr_test ( );
  i4_bset_test ( );
  i4_btest_test ( );
  i4_to_bvec_test ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "BVEC_PRB\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void bvec_add_test ( )

/******************************************************************************/
/*
  Purpose:

    BVEC_ADD_TEST tests BVEC_ADD;

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    09 June 2011

  Author:

    John Burkardt
*/
{
# define N 10

  int bvec1[N];
  int bvec2[N];
  int bvec3[N];
  int bvec4[N];
  int i;
  int j;
  int k;
  int l;
  int n = N;
  int seed = 123456789;
  int test;
  int test_num = 10;

  printf ( "\n" );
  printf ( "BVEC_ADD_TEST\n" );
  printf ( "  BVEC_ADD adds binary vectors representing integers;\n" );
  printf ( "\n" );
  printf ( "        I        J        I + J  BVEC_ADD\n" );
  printf ( "\n" );

  for ( test = 1; test <= test_num; test++ )
  { 
    i = i4_uniform_ab ( -100, 100, &seed );
    j = i4_uniform_ab ( -100, 100, &seed );

    k = i + j;

    i4_to_bvec ( i, n, bvec1 );
    i4_to_bvec ( j, n, bvec2 );

    bvec_add ( n, bvec1, bvec2, bvec3 );
    l = bvec_to_i4 ( n, bvec3 );

    printf ( "  %8d  %8d  %8d  %8d\n", i, j, k, l );
  }
  return;
# undef N
}
/******************************************************************************/

void bvec_complement2_test ( )

/******************************************************************************/
/*
  Purpose:

    BVEC_COMPLEMENT2_TEST tests BVEC_COMPLEMENT2;

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    14 November 2012

  Author:

    John Burkardt
*/
{
# define N 10

  int bvec1[N];
  int bvec2[N];
  int i;
  int j;
  int seed = 123456789;
  int test;
  int test_num = 5;

  printf ( "\n" );
  printf ( "BVEC_COMPLEMENT2_TEST\n" );
  printf ( "  BVEC_COMPLEMENT2 returns the two's complement\n" );
  printf ( "  of a (signed) binary vector;\n" );
  printf ( "\n" );

  for ( test = 1; test <= test_num; test++ )
  {
    i = i4_uniform_ab ( -100, 100, &seed );

    i4_to_bvec ( i, N, bvec1 );

    bvec_complement2 ( N, bvec1, bvec2 );

    j = bvec_to_i4 ( N, bvec2 );

    printf ( "\n" );
    printf ( "  I = %d\n", i );
    printf ( "  J = %d\n", j );
    bvec_print ( N, bvec1, "" );
    bvec_print ( N, bvec2, "" );

  }

  return;
# undef N
}
/******************************************************************************/

void bvec_mul_test ( )

/******************************************************************************/
/*
  Purpose:

    BVEC_MUL_TEST tests BVEC_MUL;

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    14 November 2012

  Author:

    John Burkardt
*/
{
# define N 15

  int bvec1[N];
  int bvec2[N];
  int bvec3[N];
  int i;
  int j;
  int k;
  int l;
  int seed = 123456789;
  int test;
  int test_num = 10;

  printf ( "\n" );
  printf ( "BVEC_MUL_TEST\n" );
  printf ( "  BVEC_MUL multiplies binary vectors\n" );
  printf ( "  representing integers;\n" );
  printf ( "\n" );
  printf ( "        I        J        I * J  BVEC_MUL\n" );
  printf ( "\n" );

  for ( test = 1; test <= test_num; test++ )
  { 
    i = i4_uniform_ab ( -100, 100, &seed );
    j = i4_uniform_ab ( -100, 100, &seed );

    k = i * j;

    i4_to_bvec ( i, N, bvec1 );
    i4_to_bvec ( j, N, bvec2 );
    bvec_mul ( N, bvec1, bvec2, bvec3 );
    l = bvec_to_i4 ( N, bvec3 );

    printf ( "  %8d  %8d  %8d  %8d\n", i, j, k, l );
  }

  return;
# undef N
}
/******************************************************************************/

void bvec_next_test ( )

/******************************************************************************/
/*
  Purpose:

    BVEC_NEXT_TEST tests BVEC_NEXT.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    02 January 2015

  Author:

    John Burkardt
*/
{ 
  int *b;
  int i;
  int n = 4;

  printf ( "\n" );
  printf ( "BVEC_NEXT_TEST\n" );
  printf ( "  BVEC_NEXT computes the 'next' BVEC.\n" );
  printf ( "\n" );

  b = ( int * ) malloc ( n * sizeof ( int ) );

  for ( i = 0; i < n; i++ )
  {
    b[i] = 0;
  }

  for ( i = 0; i <= 16; i++ )
  {
    bvec_print ( n, b, "" );
    bvec_next ( n, b );
  }

  free ( b );

  return;
}
/******************************************************************************/

void bvec_next_grlex_test ( )

/******************************************************************************/
/*
  Purpose:

    BVEC_NEXT_GRLEX_TEST tests BVEC_NEXT_GRLEX.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    13 March 2015

  Author:

    John Burkardt
*/
{ 
  int *b;
  int i;
  int j;
  int n = 4;

  printf ( "\n" );
  printf ( "BVEC_NEXT_GRLEX_TEST\n" );
  printf ( "  BVEC_NEXT_GRLEX computes binary vectors in GRLEX order.\n" );
  printf ( "\n" );

  b = ( int * ) malloc ( n * sizeof ( int ) );

  for ( j = 0; j < n; j++ )
  {
    b[j] = 0;
  }

  for ( i = 0; i <= 16; i++ )
  {
    printf ( "  %2d:  ", i );
    for ( j = 0; j < n; j++ )
    {
      printf ( "%d", b[j] );
    }
    printf ( "\n" );
    bvec_next_grlex ( n, b );
  }

  free ( b );

  return;
}
/******************************************************************************/

void bvec_print_test ( )

/******************************************************************************/
/*
  Purpose:

    BVEC_PRINT_TEST tests BVEC_PRINT.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    24 December 2014

  Author:

    John Burkardt
*/
{
  int n = 10;
  int bvec[10] = { 1, 0, 0, 1, 0, 1, 1, 1, 0, 0 };

  printf ( "\n" );
  printf ( "BVEC_PRINT_TEST\n" );
  printf ( "  BVEC_PRINT prints a binary vector.\n" );

  bvec_print ( n, bvec, "  BVEC:" );

  return;
}
/******************************************************************************/

void bvec_sub_test ( )

/******************************************************************************/
/*
  Purpose:

    BVEC_SUB_TEST tests BVEC_SUB;

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    09 June 2011

  Author:

    John Burkardt
*/
{
# define N 10

  int bvec1[N];
  int bvec2[N];
  int bvec3[N];
  int bvec4[N];
  int i;
  int j;
  int k;
  int l;
  int n = N;
  int seed = 123456789;
  int test;
  int test_num = 10;

  printf ( "\n" );
  printf ( "BVEC_SUB_TEST\n" );
  printf ( "  BVEC_SUB subtracts binary vectors representing integers;\n" );
  printf ( "\n" );
  printf ( "        I        J        I - J  BVEC_SUB\n" );
  printf ( "\n" );

  for ( test = 1; test <= test_num; test++ )
  { 
    i = i4_uniform_ab ( -100, 100, &seed );
    j = i4_uniform_ab ( -100, 100, &seed );

    k = i - j;

    i4_to_bvec ( i, n, bvec1 );
    i4_to_bvec ( j, n, bvec2 );
    bvec_sub ( n, bvec1, bvec2, bvec4 );
    l = bvec_to_i4 ( n, bvec4 );

    printf ( "  %8d  %8d  %8d  %8d\n", i, j, k, l );
  }
  return;
# undef N
}
/******************************************************************************/

void bvec_to_i4_test ( )

/******************************************************************************/
/*
  Purpose:

    BVEC_TO_I4_TEST tests BVEC_TO_I4;

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    19 February 2014

  Author:

    John Burkardt
*/
{
# define N 10

  int bvec[N];
  int i;
  int i2;
  int j;

  printf ( "\n" );
  printf ( "BVEC_TO_I4_TEST\n" );
  printf ( "  BVEC_TO_I4 converts a signed binary vector\n" );
  printf ( "  to an integer;\n" );
  printf ( "\n" );
  printf ( "  I --> BVEC  -->  I\n" );
  printf ( "\n" );

  for ( i = -3; i <= 10; i++ )
  {
    i4_to_bvec ( i, N, bvec );
    i2 = bvec_to_i4 ( N, bvec );

    printf ( "%3d  ", i );
    for ( j = 0; j < N; j++ )
    {
      printf ( "%1d", bvec[j] );
    }
    printf ( "  %3d\n", i2 );
  }

  return;
# undef N
}
/******************************************************************************/

void bvec_uniform_new_test ( )

/******************************************************************************/
/*
  Purpose:

    BVEC_UNIFORM_NEW_TEST tests BVEC_UNIFORM_NEW;

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    25 December 2014

  Author:

    John Burkardt
*/
{
  int *bvec;
  int n = 10;
  int seed = 123456789;
  int test;

  printf ( "\n" );
  printf ( "BVEC_UNIFORM_NEW_TEST\n" );
  printf ( "  BVEC_UNIFORM_NEW randomly generates BVEC's;\n" );
  printf ( "\n" );

  for ( test = 1; test <= 10; test++ )
  { 
    bvec = bvec_uniform_new ( n, &seed );
    bvec_print ( n, bvec, "" );
    free ( bvec );
  }
  return;
}
/******************************************************************************/

void i4_bclr_test ( )

/******************************************************************************/
/*
  Purpose:

    I4_BCLR_TEST tests I4_BCLR.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    23 January 2007

  Author:

    John Burkardt
*/
{
# define TEST_NUM 2

  int i4;
  int i4_test[TEST_NUM] = { 101, -31 };
  int ivec[32];
  int j1;
  int pos;
  int test;

  printf ( "\n" );
  printf ( "I4_BCLR_TEST\n" );
  printf ( "  I4_BCLR sets a given bit to 0.\n" );

  for ( test = 0; test < TEST_NUM; test++ )
  {
    i4 = i4_test[test];

    i4_to_bvec ( i4, 32, ivec );

    printf ( "\n" );
    printf ( "  Working on I4 = %d\n", i4 );
    printf ( "\n" );
    printf ( "       Pos     Digit       I4_BCLR\n" );
    printf ( "\n" );

    for ( pos = 0; pos <= 31; pos++ )
    {
      j1 = i4_bclr ( i4, pos );

      printf ( "  %8d  %8d  %8d\n", pos, ivec[pos], j1 );
    }
  }

  return;
# undef TEST_NUM
}
/******************************************************************************/

void i4_bset_test ( )

/******************************************************************************/
/*
  Purpose:

    I4_BSET_TEST tests I4_BSET.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    23 January 2007

  Author:

    John Burkardt
*/
{
# define TEST_NUM 2

  int i4;
  int i4_test[TEST_NUM] = { 101, -31 };
  int ivec[32];
  int j1;
  int pos;
  int test;

  printf ( "\n" );
  printf ( "I4_BSET_TEST\n" );
  printf ( "  I4_BSET sets a given bit to 0.\n" );

  for ( test = 0; test < TEST_NUM; test++ )
  {
    i4 = i4_test[test];

    i4_to_bvec ( i4, 32, ivec );

    printf ( "\n" );
    printf ( "  Working on I4 = %d\n", i4 );
    printf ( "\n" );
    printf ( "       Pos     Digit       I4_BSET\n" );
    printf ( "\n" );

    for ( pos = 0; pos <= 31; pos++ )
    {
      j1 = i4_bset ( i4, pos );

      printf ( "  %8d  %8d  %8d\n", pos, ivec[pos], j1 );
    }
  }

  return;
# undef TEST_NUM
}
/******************************************************************************/

void i4_btest_test ( )

/******************************************************************************/
/*
  Purpose:

    I4_BTEST_TEST tests I4_BTEST.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    23 January 2007

  Author:

    John Burkardt
*/
{
# define TEST_NUM 2

  int i4;
  int i4_test[TEST_NUM] = { 101, -31 };
  int ivec[32];
  int j1;
  int pos;
  int test;

  printf ( "\n" );
  printf ( "I4_BTEST_TEST\n" );
  printf ( "  I4_BTEST reports whether a given bit is 0 or 1.\n" );

  for ( test = 0; test < TEST_NUM; test++ )
  {
    i4 = i4_test[test];

    i4_to_bvec ( i4, 32, ivec );

    printf ( "\n" );
    printf ( "  Analyze the integer I4 = %d\n", i4 );
    printf ( "\n" );
    printf ( "       Pos     Digit  I4_BTEST\n" );
    printf ( "\n" );

    for ( pos = 0; pos <= 31; pos++ )
    {
      j1 = i4_btest ( i4, pos );

      printf ( "  %8d  %8d  %8d\n", pos, ivec[pos], j1 );
    }
  }

  return;
# undef TEST_NUM
}
/******************************************************************************/

void i4_to_bvec_test ( )

/******************************************************************************/
/*
  Purpose:

    I4_TO_BVEC_TEST tests I4_TO_BVEC.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    19 February 2014

  Author:

    John Burkardt
*/
{
# define N 10

  int bvec[N];
  int i;
  int i2;
  int j;

  printf ( "\n" );
  printf ( "I4_TO_BVEC_TEST\n" );
  printf ( "  I4_TO_BVEC converts an integer to a \n" );
  printf ( "  signed binary vector;\n" );
  printf ( "\n" );
  printf ( "  I --> BVEC  -->  I\n" );
  printf ( "\n" );

  for ( i = -3; i <= 10; i++ )
  {
    i4_to_bvec ( i, N, bvec );
    i2 = bvec_to_i4 ( N, bvec );

    printf ( "%3d  ", i );
    for ( j = 0; j < N; j++ )
    {
      printf ( "%1d", bvec[j] );
    }
    printf ( "  %3d\n", i2 );
  }

  return;
# undef N
}

