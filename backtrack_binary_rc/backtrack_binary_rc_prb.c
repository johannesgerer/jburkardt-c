# include <stdlib.h>
# include <stdio.h>
# include <stdbool.h>

# include "backtrack_binary_rc.h"

int main ( );
void test01 ( );
void test02 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for BACKTRACK_BINARY_RC_PRB.

  Discussion:

    BACKTRACK_BINARY_RC_PRB tests BACKTRACK_BINARY_RC.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 January 2014

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "BACKTRACK_BINARY_RC_PRB:\n" );
  printf ( "  C version\n" );
  printf ( "  Test the BACKTRACK_BINARY_RC library.\n" );

  test01 ( );
  test02 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "BACKTRACK_BINARY_RC_PRB:\n" );
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

    TEST01 seeks a selection of binary powers that have a given sum.

  Discussion:

    We consider the binary powers 1, 2, 4, ... 2^(n-1).

    We wish to select some of these powers, so that the sum is equal
    to a given target value.  We are actually simply seeking the binary
    representation of an integer.

    A partial solution is acceptable if it is less than the target value.

    We list the powers in descending order, so that the bactracking
    procedure makes the most significant choices first, thus quickly
    eliminating many unsuitable choices.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 January 2014

  Author:

    John Burkardt
*/
{
  int call_num;
  int choice[8];
  int factor;
  int i;
  int n = 8;
  int n2;
  bool reject;
  int result;
  int target;
  int targets[3] = { 73, 299, -3 };
  int test;
  int test_num = 3;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  Use BACKBIN_RC to find the binary expansion of\n" );
  printf ( "  an integer between 0 and 255.\n" );
  printf ( "  The choices are 0/1 for the 8 digits.\n" );

  for ( test = 0; test < test_num; test++ )
  {
    target = targets[test];
    printf ( "\n" );
    printf ( "  TARGET = %d\n", target );
    call_num = 0;
    n2 = -1;

    for ( ; ; )
    {
      backbin_rc ( n, reject, &n2, choice );
      call_num = call_num + 1;

      if ( n2 == -1 )
      {
        printf ( "  Termination without solution.\n" );
        break;
      }
/*
  Evaluate the integer determined by the choices.
*/
      factor = 1;
      for ( i = n; n2 < i; i-- )
      {
        factor = factor * 2;
      }

      result = 0;
      for ( i = 0; i < n2; i++ )
      {
        result = result * 2 + choice[i];
      }

      result = result * factor;
/*
  If the integer is too big, then we reject it, and
  all the related integers formed by making additional choices.
*/
      reject = ( target < result );
/*
  If we hit the target, then in this case, we can exit because
  the solution is unique.
*/
      if ( result == target )
      {
        break;
      }

    }

    printf ( "  Number of calls = %d\n", call_num );
    printf ( "  Binary search space = %d\n", i4_power ( 2, n ) );
    printf ( "  " );
    for ( i = 0; i < n; i++ )
    {
      printf ( "%2d", choice[i] );
    }
    printf ( "\n" );
  }

  return;
}
/******************************************************************************/

void test02 ( )

/******************************************************************************/
/*
  Purpose:

    TEST02 seeks a subset of a set of numbers which add to a given sum.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 January 2014

  Author:

    John Burkardt
*/
{
  int call_num;
  int choice[8];
  int i;
  int n = 8;
  int n2;
  bool reject;
  int result;
  int target = 53;
  int test;
  int w[8] = { 15, 22, 14, 26, 32, 9, 16, 8 };

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  Use BACKBIN_RC to seek subsets of a set W\n" );
  printf ( "  that sum to a given target value.\n" );
  printf ( "  The choices are 0/1 to select each element of W.\n" );

  printf ( "\n" );
  printf ( "  TARGET = %d\n", target );
  printf ( "\n" );
  call_num = 0;
  n2 = -1;

  for ( ; ; )
  {
    backbin_rc ( n, reject, &n2, choice );
    call_num = call_num + 1;

    if ( n2 == -1 )
    {
      break;
    }
/*
  Evaluate the partial sum.
*/
    result = 0;
    for ( i = 0; i < n2; i++ )
    {
      result = result + choice[i] * w[i];
    }
/*
  If the sum is too big, then we reject it, and
  all the related sums formed by making additional choices.
*/
    reject = ( target < result );
/*
  If we hit the target, print out the information.
*/
    if ( result == target && n2 == n )
    {
      printf ( "  " );
      for ( i = 0; i < n; i++ )
      {
        printf ( "%2d", choice[i] );
      }
      printf ( "\n" );
    }
  }

  printf ( "\n" );
  printf ( "  Number of calls = %d\n", call_num );
  printf ( "  Binary search space = %d\n", i4_power ( 2, n ) );

  return;
}
