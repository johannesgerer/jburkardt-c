# include <stdlib.h>
# include <stdio.h>
# include <stdbool.h>
# include <math.h>

# include "combination_lock.h"

int main ( );
void test01 ( );
void test02 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for COMBINATION_LOCK.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    13 May 2012

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "COMBINATION_LOCK\n" );
  printf ( "  C version\n" );
  printf ( "  Test the COMBINATION_LOCK libary.\n" );

  test01 ( );
  test02 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "COMBINATION_LOCK\n" );
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

    TEST01 tests BICYCLE_LOCK.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    13 May 2012

  Author:

    John Burkardt
*/
{
  int c;
  int m = 3;
  int n = 10;
  int seed;
  int step;
  
  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  A bicycle combination lock consists of 3 dials,\n" );
  printf ( "  each having 10 symbols, 0 through 9.\n" );
  printf ( "  We seek to determine the combination C.\n" );
/*
  Report on the problem data.
*/
  printf ( "\n" );
  printf ( "  The number of dials is M = %d\n", m );
  printf ( "  The number of symbols is N = %d\n", n );
  printf ( "  The number of possible combinations is M^N = %d\n", i4_power ( n, m ) );

  seed = get_seed ( );
  c = i4_uniform ( 0, 999, &seed );

  printf ( "  The \"secret\" combination is %d\n", c );

  step = bicycle_lock ( c );

  if ( step == -1 )
  {
    printf ( "\n" );
    printf ( "  The combination was not found!\n" );
  }
  else
  {
    printf ( "\n" );
    printf ( "  The combination was found on step %d\n", step );
  }
  return;
}
/******************************************************************************/

void test02 ( )

/******************************************************************************/
/*
  Purpose:

    TEST02 tests COMBINATION_LOCK.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    13 May 2012

  Author:

    John Burkardt
*/
{
  int c[4] = { 1, 2, 3, 4 };
  int m = 4;
  int n = 5;
  int step;
  
  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  A combination lock consists of M dials,\n" );
  printf ( "  each having N symbols.\n" );
  printf ( "  We seek to determine the combination C.\n" );
/*
  Report on the problem data.
*/
  printf ( "\n" );
  printf ( "  The number of dials is M = %d\n", m );
  printf ( "  The number of symbols is N = %d\n", n );
  printf ( "  The number of possible combinations is M^N = %d\n", i4_power ( n, m ) );

  i4vec_print ( m, c, "  The \"secret\" combination:" );

  step = combination_lock ( m, n, c );

  if ( step == -1 )
  {
    printf ( "\n" );
    printf ( "  The combination was not found!\n" );
  }
  else
  {
    printf ( "\n" );
    printf ( "  The combination was found on step %d\n", step );
  }

  return;
}
