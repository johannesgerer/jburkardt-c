# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "partition_problem.h"

int main ( void );
void test01 ( int n, int w[] );
void test02 ( int n, int w[] );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for PARTITION_PROBLEM_PRB.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 May 2012

  Author:

    John Burkardt
*/
{
  int n;
  int test;
  int test_num = 5;
  int *w;
  int w1[5] = { 19, 17, 13, 9, 6 };
  int w2[9] = { 484, 114, 205, 288, 506, 503, 201, 127, 410 };
  int w3[10] = { 771, 121, 281, 854, 885, 734, 486, 1003, 83, 62 };
  int w4[10] = { 2, 10, 3, 8, 5, 7, 9, 5, 3, 2 };
  int w5[9] = { 3, 4, 3, 1, 3, 2, 3, 2, 1 };

  timestamp ( );
  printf ( "\n" );
  printf ( "PARTITION_PROBLEM_PRB:\n" );
  printf ( "  C version\n" );
  printf ( "  Test the PARTITION_PROBLEM library.\n" );
/*
  Find individual solutions.
*/
  for ( test = 1; test <= test_num; test++ )
  {
    if ( test == 1 )
    {
      n = 5;
      w = i4vec_copy_new ( n, w1 );
    }
    else if ( test == 2 )
    {
      n = 9;
      w = i4vec_copy_new ( n, w2 );
    }
    else if ( test == 3 )
    {
      n = 10;
      w = i4vec_copy_new ( n, w3 );
    }
    else if ( test == 4 )
    {
      n = 10;
      w = i4vec_copy_new ( n, w4 );
    }
    else if ( test == 5 )
    {
      n = 9;
      w = i4vec_copy_new ( n, w5 );
    }

    test01 ( n, w );

    free ( w );
  }
/*
  Count solutions.
*/
  for ( test = 1; test <= test_num; test++ )
  {
    if ( test == 1 )
    {
      n = 5;
      w = i4vec_copy_new ( n, w1 );
    }
    else if ( test == 2 )
    {
      n = 9;
      w = i4vec_copy_new ( n, w2 );
    }
    else if ( test == 3 )
    {
      n = 10;
      w = i4vec_copy_new ( n, w3 );
    }
    else if ( test == 4 )
    {
      n = 10;
      w = i4vec_copy_new ( n, w4 );
    }
    else if ( test == 5 )
    {
      n = 9;
      w = i4vec_copy_new ( n, w5 );
    }

    test02 ( n, w );

    free ( w );
  }
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "PARTITION_PROBLEM_PRB\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void test01 ( int n, int w[] )

/******************************************************************************/
/*
  Purpose:

    TEST01 tests PARTITION_BRUTE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 May 2012

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of weights.

    Input, int W[N], a set of weights.
*/
{
  int *c;
  int discrepancy;
  int i;
  int w0_sum;
  int w1_sum;

  printf ( "\n" );
  printf ( "TEST01:\n" );
  printf ( "  Partition a set of N integers W so that the subsets\n" );
  printf ( "  have equal sums.\n" );

  c = ( int * ) malloc ( n * sizeof ( int ) );

  partition_brute ( n, w, c, &discrepancy );

  printf ( "\n" );
  printf ( "     I        W0        W1\n" );
  printf ( "\n" );
  w0_sum = 0;
  w1_sum = 0;
  for ( i = 0; i < n; i++ )
  {
    if ( c[i] == 0 )
    {
      w0_sum = w0_sum + w[i];
      printf ( "  %4d  %8d\n", i, w[i] );
    }
    else
    {
      w1_sum = w1_sum + w[i];
      printf ( "  %4d            %8d\n", i, w[i] );
    }
  }
  printf ( "        --------  --------\n" );
  printf ( "        %8d  %8d\n", w0_sum, w1_sum );
  printf ( "\n" );
  printf ( "  Discrepancy = %8d\n", discrepancy );

  free ( c );

  return;
}
/******************************************************************************/

void test02 ( int n, int w[] )

/******************************************************************************/
/*
  Purpose:

    TEST02 tests PARTITION_COUNT.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 May 2012

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of weights.

    Input, int W[N], a set of weights.
*/
{
  int count;
  int i;

  printf ( "\n" );
  printf ( "TEST02:\n" );
  printf ( "  PARTITION_COUNT counts the number of exact solutions\n" );
  printf ( "  of the partition problem.\n" );

  count = partition_count ( n, w );

  printf ( "\n" );
  printf ( "     I        W\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %4d  %8d\n", i, w[i] );
  }
  printf ( "\n" );
  printf ( "  Number of solutions = %d\n", count );

  return;
}
