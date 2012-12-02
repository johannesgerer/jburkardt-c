# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "subset.h"

int main ( void );
void test000 ( void );
void test001 ( void );
void test002 ( void );
void test003 ( void );
void test9001 ( void );
void test9002 ( void );
void test9003 ( void );
void test005 ( void );
void test006 ( void );
void test011 ( void );
void test012 ( void );
void test0322 ( void );
void test03225 ( void );
void test0323 ( void );
void test06225 ( void );
void test06895 ( void );
void test069 ( void );
void test07715 ( void );
void test0955 ( void );
void test097 ( void );
void test15696 ( void );
void test15698 ( void );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    SUBSET_PRB calls the SUBSET tests.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    14 November 2012

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "SUBSET_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the SUBSET library.\n" );

  test000 ( );
  test001 ( );
  test002 ( );
  test003 ( );
  test9001 ( );
  test9002 ( );
  test9003 ( );
  test005 ( );
  test006 ( );

  test011 ( );
  test012 ( );

  test0322 ( );
  test03225 ( );
  test0323 ( );

  test06225 ( );
  test06895 ( );
  test069 ( );

  test07715 ( );

  test0955 ( );
  test097 ( );

  test15696 ( );
  test15698 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "SUBSET_PRB\n" );
  printf ( "  Normal end of execution.\n" );

  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void test000 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST000 tests RANDOM_INITIALIZE.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    02 December 2006

  Author:

    John Burkardt
*/
{
  unsigned long seed;

  printf ( "\n" );
  printf ( "TEST000\n" );
  printf ( "  Call RANDOM_INITIALIZE to initialize the\n" );
  printf ( "  RANDOM random number generator.\n" );

  seed = 0;
  seed = random_initialize ( seed );

  printf ( "\n" );
  printf ( "  RANDOM_INITIALIZE returns SEED = %lu\n", seed );

  return;
}
/******************************************************************************/

void test001 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST001 tests ASM_ENUM.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    08 February 2007

  Author:

    John Burkardt
*/
{
  int n;

  printf ( "\n" );
  printf ( "TEST001\n" );
  printf ( "  ASM_ENUM returns the number of alternating sign\n" );
  printf ( "  matrices of a given order.\n" );
  printf ( "\n" );

  for ( n = 0; n <= 7; n++ )
  {
    printf ( "  %4d  %6d\n", n, asm_enum ( n ) );
  }

  return;
}
/******************************************************************************/

void test002 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST002 tests ASM_TRIANGLE.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    01 July 2010

  Author:

    John Burkardt
*/
{
# define N_MAX 7

  int a[N_MAX+1];
  int i;
  int n;

  printf ( "\n" );
  printf ( "TEST002\n" );
  printf ( "  ASM_TRIANGLE returns a row of the alternating sign\n" );
  printf ( "  matrix triangle.\n" );
  printf ( "\n" );

  for ( n = 0; n <= N_MAX; n++ )
  {
    asm_triangle ( n, a );
    printf ( "%4d  ", n );
    for ( i = 0; i <= n; i++ )
    {
      printf ( "%8d  ", a[i] );
    }
    printf ( "\n" );
  }

  return;
# undef N_MAX
}
/******************************************************************************/

void test003 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST003 tests BELL and BELL_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    03 August 2010

  Author:

    John Burkardt
*/
{
  int c;
  int *c2;
  int n;
  int n_data;

  printf ( "\n" );
  printf ( "TEST003\n" );
  printf ( "  BELL computes Bell numbers.\n" );
  printf ( "  BELL_VALUES returns some exact values.\n" );
  printf ( "\n" );
  printf ( "  N  exact C(I)  computed C(I)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    bell_values ( &n_data, &n, &c );

    if ( n_data == 0 )
    {
      break;
    }

    c2 = ( int * ) malloc ( ( n + 1 ) * sizeof ( int ) );

    bell ( n, c2 );

    printf ( "  %4d  %8d  %8d\n", n, c, c2[n] );

    free ( c2 );
  }

  return;
}
/******************************************************************************/

void test9001 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST9001 tests BVEC_ADD and BVEC_SUB;

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
  printf ( "TEST9001\n" );
  printf ( "  BVEC_ADD adds binary vectors representing integers;\n" );
  printf ( "  BVEC_SUB subtracts binary vectors representing integers;\n" );
  printf ( "\n" );
  printf ( "        I        J        K = I + J    L = I - J\n" );
  printf ( "\n" );

  for ( test = 1; test <= test_num; test++ )
  { 
    i = i4_uniform ( -100, 100, &seed );
    j = i4_uniform ( -100, 100, &seed );

    printf ( "\n" );
    printf ( "  %8d  %8d\n", i, j );

    k = i + j;
    l = i - j;

    printf ( "  Directly:           %8d  %8d\n", k, l );

    i4_to_bvec ( i, n, bvec1 );
    i4_to_bvec ( j, n, bvec2 );

    bvec_add ( n, bvec1, bvec2, bvec3 );
    k = bvec_to_i4 ( n, bvec3 );

    bvec_sub ( n, bvec1, bvec2, bvec4 );
    l = bvec_to_i4 ( n, bvec4 );

    printf ( "  BVEC_ADD, BVEC_SUB  %8d  %8d\n", k, l );
  }
  return;
# undef N
}
/******************************************************************************/

void test9002 ( )

/******************************************************************************/
/*
  Purpose:

    TEST9002 tests BVEC_COMPLEMENT2;

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
  printf ( "TEST9002\n" );
  printf ( "  BVEC_COMPLEMENT2 returns the two's complement\n" );
  printf ( "  of a (signed) binary vector;\n" );
  printf ( "\n" );

  for ( test = 1; test <= test_num; test++ )
  {
    i = i4_uniform ( -100, 100, &seed );

    i4_to_bvec ( i, N, bvec1 );

    bvec_complement2 ( N, bvec1, bvec2 );

    j = bvec_to_i4 ( N, bvec2 );

    printf ( "\n" );
    printf ( "  I = %d\n", i );
    printf ( "  J = %d\n", j );
    bvec_print ( N, bvec1, " " );
    bvec_print ( N, bvec2, " " );

  }

  return;
# undef N
}
/******************************************************************************/

void test9003 ( )

/******************************************************************************/
/*
  Purpose:

    TEST9003 tests BVEC_MUL;

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
  int seed = 123456789;
  int test;
  int test_num = 10;

  printf ( "\n" );
  printf ( "TEST9003\n" );
  printf ( "  BVEC_MUL multiplies binary vectors\n" );
  printf ( "  representing integers;\n" );
  printf ( "\n" );
  printf ( "        I        J        K = I * J\n" );
  printf ( "\n" );

  for ( test = 1; test <= test_num; test++ )
  { 
    i = i4_uniform ( -100, 100, &seed );
    j = i4_uniform ( -100, 100, &seed );

    printf ( "\n" );
    printf ( "  %8d  %8d\n", i, j );

    k = i * j;

    printf ( "  Directly:           %8d\n", k );

    i4_to_bvec ( i, N, bvec1 );
    i4_to_bvec ( j, N, bvec2 );

    bvec_mul ( N, bvec1, bvec2, bvec3 );
    k = bvec_to_i4 ( N, bvec3 );

    printf ( "  BVEC_MUL            %8d\n", k );
  }

  return;
# undef N
}
/******************************************************************************/

void test005 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST005 tests CATALAN and CATALAN_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    02 October 2010

  Author:

    John Burkardt
*/
{
  int c;
  int *c2;
  int n;
  int n_data;

  printf ( "\n" );
  printf ( "TEST005\n" );
  printf ( "  CATALAN computes Catalan numbers.\n" );
  printf ( "  CATALAN_VALUES returns some exact values.\n" );
  printf ( "\n" );
  printf ( "  N  exact C(I)  computed C(I)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    catalan_values ( &n_data, &n, &c );

    if ( n_data == 0 )
    {
      break;
    }

    c2 = ( int * ) malloc ( ( n + 1 ) * sizeof ( int ) );

    catalan ( n, c2 );

    printf ( "  %4d  %8d  %8d\n", n, c, c2[n] );

    free ( c2 );
  }

  return;
}
/******************************************************************************/

void test006 ( )

/******************************************************************************/
/*
  Purpose:

    TEST006 tests CATALAN_ROW_NEXT.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 October 2006

  Author:

    John Burkardt
*/
{
# define N_MAX 10

  int c[N_MAX+1];
  int i;
  int n;
  int next;

  printf ( "\n" );
  printf ( "TEST006\n" );
  printf ( "  CATALAN_ROW_NEXT computes a row of the Catalan triangle.\n" );
  printf ( "\n" );
  printf ( "  First, compute row 7:\n" );
  printf ( "\n" );

  next = 0;
  n = 7;
  catalan_row_next ( next, n, c );

  printf ( "%4d", n );
  for ( i = 0; i <= n; i++ )
  {
    printf ( "  %8d", c[i] );
  }
  printf ( "\n" );

  printf ( "\n" );
  printf ( "  Now compute rows consecutively, one at a time:\n" );
  printf ( "\n" );

  next = 0;

  for ( n = 0; n <= N_MAX; n++ )
  {
    catalan_row_next ( next, n, c );
    next = 1;

    printf ( "%4d", n );
    for ( i = 0; i <= n; i++ )
    {
      printf ( "  %8d", c[i] );
    }
    printf ( "\n" );
  }

  return;
# undef N_MAX
}
/******************************************************************************/

void test011 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST011 tests COMB_NEXT.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    10 April 2009

  Author:

    John Burkardt
*/
{
# define N 5

  int a[N];
  int done;
  int i;
  int j;
  int k;
  int n = N;

  printf ( "\n" );
  printf ( "TEST011\n" );
  printf ( "  COMB_NEXT produces combinations.\n" );

  for ( k = 1; k <= n; k++ )
  {
    printf ( "\n" );
    printf ( "  Combinations of size K = %d\n", k );
    printf ( "\n" );

    done = 1;

    for ( ; ; )
    {
      comb_next ( n, k, a, &done );

      if ( done )
      {
        break;
      }

      for ( i = 0; i < k; i++ )
      {
        printf ( "  %4d", a[i] );
      }
      printf ( "\n" );
    }
  }

  return;
# undef N
}
/******************************************************************************/

void test012 ( )

/******************************************************************************/
/*
  Purpose:

    TEST012 tests COMB_ROW.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    14 May 2011

  Author:

    John Burkardt
*/
{
# define N 10

  int c[N+1];
  int i;
  int j;
  int next;

  printf ( "\n" );
  printf ( "TEST012\n" );
  printf ( "  COMB_ROW computes a row of the Pascal triangle.\n" );
  printf ( "\n" );

  next = 0;

  for ( i = 0; i <= N; i++ )
  {
    comb_row ( next, i, c );
    next = 1;
    printf ( "  %2d  ", i );
    for ( j = 0; j <= i; j++ )
    {
      printf ( "%5d", c[j] );
    }
    printf ( "\n" );
  }

  return;
# undef N
}
/******************************************************************************/

void test0322 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST0322 tests I4_BCLR.

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
  printf ( "TEST0322\n" );
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

void test03225 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST03225 tests I4_BSET.

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
  printf ( "TEST03225\n" );
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

void test0323 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST0323 tests I4_BTEST.

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
  printf ( "TEST0323\n" );
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

void test06225 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST06225 tests I4_PARTITIONS_NEXT.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    03 August 2010

  Author:

    John Burkardt
*/
{
  int i;
  int j;
  int m[3];
  int msum;
  int s = 3;

  printf ( "\n" );
  printf ( "TEST06225\n" );
  printf ( "  I4_PARTITIONS_NEXT produces the next\n" );
  printf ( "  nondecreasing partitions of an integer, and\n" );
  printf ( "  if necessary, increments the integer to keep on going.\n" );

  i = 0;
  m[0] = 0;
  m[1] = 0;
  m[2] = 0;

  printf ( "\n" );
  printf ( "   I Sum    Partition\n" );
  printf ( "\n" );
  msum = i4vec_sum ( s, m );
  printf ( "  %2d  %2d    ", i, msum );
  for ( j = 0; j < s; j++ )
  {
    printf ( "%2d", m[j] );
  }
  printf ( "\n" );

  for ( i = 1; i <= 15; i++ )
  {
    i4_partitions_next ( s, m );
    msum = i4vec_sum ( s, m );
    printf ( "  %2d  %2d    ", i, msum );
    for ( j = 0; j < s; j++ )
    {
      printf ( "%2d", m[j] );
    }
    printf ( "\n" );
  }
  printf ( "\n" );
  printf ( "  You can start from any legal partition.\n" );
  printf ( "  Here, we restart at ( 2, 1, 0 ).\n" );

  i = 0;
  m[0] = 2;
  m[1] = 1;
  m[2] = 0;

  printf ( "\n" );
  printf ( "   I Sum    Partition\n" );
  printf ( "\n" );
  msum = i4vec_sum ( s, m );
  printf ( "  %2d  %2d    ", i, msum );
  for ( j = 0; j < s; j++ )
  {
    printf ( "%2d", m[j] );
  }
  printf ( "\n" );

  for ( i = 1; i <= 15; i++ )
  {
    i4_partitions_next ( s, m );
    msum = i4vec_sum ( s, m );
    printf ( "  %2d  %2d    ", i, msum );
    for ( j = 0; j < s; j++ )
    {
      printf ( "%2d", m[j] );
    }
    printf ( "\n" );
  }
  return;
}
/******************************************************************************/

void test06895 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST06895 tests INVERSE_MOD_N.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    14 May 2011

  Author:

    John Burkardt
*/
{
  int b;
  int n;
  int y;
  int z;

  printf ( "\n" );
  printf ( "TEST06895\n" );
  printf ( "  INVERSE_MOD_N seeks Y, the inverse of B mod N,\n" );
  printf ( "  so that mod ( B * Y, N ) = 1, but returns 0\n" );
  printf ( "  if the inverse does not exist.\n" );

  printf ( "\n" );
  printf ( "     B     N     Y     Z = ( ( B * Y ) %% N )\n" );
  printf ( "\n" );

  for ( n = 1; n <= 10; n++ )
  {
    for ( b = 1; b < n; b++ )
    {
      y = inverse_mod_n ( b, n );
      z = ( ( b * y ) % n );
      printf ( "  %4d  %4d  %4d  %4d\n", b, n, y, z );
    }
  }

  return;
}
/******************************************************************************/

void test069 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST069 tests JFRAC_TO_RFRAC and RFRAC_TO_JFRAC.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    02 October 2010

  Author:

    John Burkardt
*/
{
# define MAXM 10

  int i;
  int m;
  double p[MAXM];
  double q[MAXM];
  double r[MAXM];
  double s[MAXM];
  int seed;
/*
  Generate the data, but force Q(M+1) to be 1.  
  That will make it easier to see that the two operations are inverses
  of each other.  JFRAC_TO_RFRAC is free to scale its output, and chooses
  a scaling in which Q(M+1) is 1.
*/
  seed = 123456789;
  m = 6;
  r8vec_uniform_01 ( m, &seed, p );
  r8vec_uniform_01 ( m + 1, &seed, q );

  for ( i = 0; i < m; i++ )
  {
    q[i] = q[i] / q[m];
  }
  q[m] = 1.0;

  printf ( "\n" );
  printf ( "TEST069\n" );
  printf ( "  RFRAC_TO_JFRAC converts a rational polynomial\n" );
  printf ( "  fraction to a J fraction.\n" );
  printf ( "  JFRAC_TO_RFRAC converts a J fraction\n" );
  printf ( "  to a rational polynomial fraction.\n" );
  printf ( "\n" );
  printf ( "  The original rational polynomial coefficients:\n" );
  printf ( "\n" );

  for ( i = 0; i < m; i++ )
  {
    printf ( "%14f  ", p[i] );
  }
  printf ( "\n" );

  for ( i = 0; i < m + 1; i++ )
  {
    printf ( "%14f  ", q[i] );
  }
  printf ( "\n" );
 
  rfrac_to_jfrac ( m, p, q, r, s );
 
  printf ( "\n" );
  printf ( "  The J fraction coefficients:\n" );
  printf ( "\n" );

  for ( i = 0; i < m; i++ )
  {
    printf ( "%14f  ", r[i] );
  }
  printf ( "\n" );
  for ( i = 0; i < m; i++ )
  {
    printf ( "%14f  ", s[i] );
  }
  printf ( "\n" );
 
  jfrac_to_rfrac ( m, r, s, p, q );

  printf ( "\n" );
  printf ( "  The recovered rational polynomial:\n" );
  printf ( "\n" );

  for ( i = 0; i < m; i++ )
  {
    printf ( "%14f  ", p[i] );
  }
  printf ( "\n" );

  for ( i = 0; i < m + 1; i++ )
  {
    printf ( "%14f  ", q[i] );
  }
  printf ( "\n" );

  return;
# undef MAXM
}
/******************************************************************************/

void test07715 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST07715 tests KSUB_RANDOM5.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    09 June 2011

  Author:

    John Burkardt
*/
{
  int *a;
  int i;
  int j;
  int k = 5;
  int n = 52;
  int seed;

  printf ( "\n" );
  printf ( "TEST07715\n" );
  printf ( "  KSUB_RANDOM5 generates a random K subset of an N set.\n" );
  printf ( "  Set size is N =    %d\n", n );
  printf ( "  Subset size is K = %d\n", k );
  printf ( "\n" );

  seed = 123456789;

  for ( i = 1; i <= 5; i++ )
  {
    a = ksub_random5 ( n, k, &seed );
    for ( j = 0; j < k; j++ )
    {
      printf ( "  %3d", a[j] );
    }
    printf ( "\n" );
    free ( a );
  }
  return;
}
/******************************************************************************/

void test0955 ( )

/******************************************************************************/
/*
  Purpose:

    TEST0955 tests PERM_INVERSE3_NEW.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    14 May 2011

  Author:

    John Burkardt
*/
{
# define N 7

  int perm[N] = { 3, 2, 4, 0, 6, 5, 1 };
  int *perm_inv;

  printf ( "\n" );
  printf ( "TEST0955\n" );
  printf ( "  PERM_INVERSE3_NEW inverts a permutation.\n" );

  perm_print ( N, perm, "  The original permutation:" );
 
  perm_inv = perm_inverse3_new ( N, perm );
 
  perm_print ( N, perm_inv, "  The inverted permutation:" );
 
  free ( perm_inv );

  return;
# undef N
}
/******************************************************************************/

void test097 ( )

/******************************************************************************/
/*
  Purpose:

    TEST097 tests PERM_LEX_NEXT and PERM_SIGN.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    14 November 2012

  Author:

    John Burkardt
*/
{
# define N 4

  int i;
  int more;
  int p[N];
  int rank;
  int p_sign;

  printf ( "\n" );
  printf ( "TEST097\n" );
  printf ( "  PERM_LEX_NEXT generates permutations in order.\n" );
  printf ( "  PERM_SIGN computes the sign of a permutation.\n" );
  printf ( "\n" );
  printf ( "  RANK  SIGN  Permutation\n" );
  printf ( "\n" );

  more = 0;
  rank = 0; 

  for ( ; ; )
  {
    perm_lex_next ( N, p, &more );

    p_sign = perm_sign ( N, p );

    if ( !more )
    {
      break;
    }

    printf ( "%4d  %4d", rank, p_sign );

    for ( i = 0; i < N; i++ )
    {
      printf ( "  %4d", p[i] );
    }
    printf ( "\n" );

    rank = rank + 1;

  }
 
  return;
# undef N
}
/******************************************************************************/

void test15696 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST15696 tests VECTOR_CONSTRAINED_NEXT7.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    01 July 2010

  Author:

    John Burkardt
*/
{
# define N_MAX 3

  double alpha[N_MAX] = { 4.0, 3.0, 5.0 };
  int i;
  int j;
  int more;
  int n;
  double q_max = 20.0;
  double q_min = 16.0;
  double total;
  int x[N_MAX];
  int x_max[N_MAX] = { 2, 6, 4 };

  fprintf ( stdout, "\n" );
  fprintf ( stdout, "TEST15696\n" );
  fprintf ( stdout, "  VECTOR_CONSTRAINED_NEXT7:\n" );
  fprintf ( stdout, "  Consider vectors:\n" );
  fprintf ( stdout, "    0 <= X(1:N) <= X_MAX(1:N),\n" );
  fprintf ( stdout, "  Set\n" );
  fprintf ( stdout, "    TOTAL = sum ( ALPHA(1:N) * X(1:N) )\n" );
  fprintf ( stdout, "  Accept only vectors for which:\n" );
  fprintf ( stdout, "    Q_MIN <= TOTAL <= Q_MAX\n" );

  for ( n = 2; n <= N_MAX; n++ )
  {
    more = 0;

    fprintf ( stdout, "\n" );
    fprintf ( stdout, "  ALPHA:" );
    for ( j = 0; j < n; j++ )
    {
      fprintf ( stdout, "  %8f", alpha[j] );
    }
    fprintf ( stdout, "\n" );
    fprintf ( stdout, "  Q_MIN:" );
    fprintf ( stdout, "  %8f", q_min );
    fprintf ( stdout, "\n" );
    fprintf ( stdout, "  Q_MAX:" );
    fprintf ( stdout, "  %8f", q_max );
    fprintf ( stdout, "\n" );
    fprintf ( stdout, "  X_MAX:" );
    for ( j = 0; j < n; j++ )
    {
      fprintf ( stdout, "  %4d", x_max[j] );
    }
    fprintf ( stdout, "\n" );
    fprintf ( stdout, "\n" );

    i = 0;

    for ( ; ; )
    {
      vector_constrained_next7 ( n, alpha, x_max, x, q_min, 
        q_max, &more );

      if ( !more )
      {
        break;
      }

      total = 0.0;
      for ( j = 0; j < n; j++ )
      {
        total = total + alpha[j] * ( double ) x[j];
      }
      i = i + 1;
      fprintf ( stdout, "  %8d", i );
      fprintf ( stdout, "  %14f", total );
      for ( j = 0; j < n; j++ )
      {
        fprintf ( stdout, "  %8d", x[j] );
      }
      fprintf ( stdout, "\n" );
    }
  }

  return;
# undef N_MAX
}
/******************************************************************************/

void test15698 ( )

/******************************************************************************/
/*
  Purpose:

    TEST15698 tests VECTOR_NEXT.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    01 July 2010

  Author:

    John Burkardt
*/
{
# define N_MAX 3

  int i;
  int j;
  int more;
  int n;
  int x[N_MAX];
  int x_max[N_MAX] = { 2, 6, 4 };
  int x_min[N_MAX] = { 1, 4, 3 };

  fprintf ( stdout, "\n" );
  fprintf ( stdout, "TEST15698\n" );
  fprintf ( stdout, "  VECTOR_NEXT:\n" );
  fprintf ( stdout, "  Generate all vectors X such that:\n" );
  fprintf ( stdout, "    X_MIN(1:N) <= X(1:N) <= X_MAX(1:N),\n" );

  for ( n = 2; n <= N_MAX; n++ )
  {
    more = 0;

    fprintf ( stdout, "\n" );
    fprintf ( stdout, "    X_MIN:" );
    for ( j = 0; j < n; j++ )
    {
      fprintf ( stdout, "  %8d", x_min[j] );
    }
    fprintf ( stdout, "\n" );

    i = 0;

    for ( ; ; )
    {
      vector_next ( n, x_min, x_max, x, &more );

      if ( !more )
      {
        break;
      }

      i = i + 1;
      fprintf ( stdout, "  %8d", i );
      for ( j = 0; j < n; j++ )
      {
        fprintf ( stdout, "  %8d", x[j] );
      }
      fprintf ( stdout, "\n" );
    }
    fprintf ( stdout, "    X_MAX:" );
    for ( j = 0; j < n; j++ )
    {
      fprintf ( stdout, "  %8d", x_max[j] );
    }
    fprintf ( stdout, "\n" );
  }

  return;
# undef N_MAX
}
