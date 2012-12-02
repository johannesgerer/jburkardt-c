# include <stdlib.h>
# include <stdio.h>

# include "cycle_floyd.h"

int main ( );

void test01 ( );
int f1 ( int i );
void test02 ( );
int f2 ( int i );
void test03 ( );
int f3 ( int i );
void test04 ( );
int f4 ( int i );
void test05 ( );
int f5 ( int i );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for CYCLE_FLOYD_PRB.

  Discussion:

    CYCLE_FLOYD_PRB tests the CYCLE_FLOYD library.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    18 June 2012

  Author:

    John Burkardt
*/
{
  timestamp ( );

  printf ( "\n" );
  printf ( "CYCLE_FLOYD_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the CYCLE_FLOYD library.\n" );

  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
  test05 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "CYCLE_FLOYD_PRB\n" );
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

    TEST01 tests CYCLE_FLOYD for a tiny example.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    18 June 2012

  Author:

    John Burkardt
*/
{
  int lam;
  int mu;
  int x0;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  Test CYCLE_FLOYD on F1().\n" );
  printf ( "  f1(0) = 6.\n" );
  printf ( "  f1(1) = 6.\n" );
  printf ( "  f1(2) = 0.\n" );
  printf ( "  f1(3) = 1.\n" );
  printf ( "  f1(4) = 4.\n" );
  printf ( "  f1(5) = 3.\n" );
  printf ( "  f1(6) = 3.\n" );
  printf ( "  f1(7) = 4.\n" );
  printf ( "  f1(8) = 0.\n" );

  x0 = 2;
  printf ( "\n" );
  printf ( "  Starting argument X0 = %d\n", x0 );

  cycle_floyd ( f1, x0, &lam, &mu );

  printf ( "\n" );
  printf ( "  Reported cycle length is %d\n", lam );
  printf ( "  Expected value is 3\n" );
  printf ( "\n" );
  printf ( "  Reported distance to first cycle element is %d\n", mu );
  printf ( "  Expected value is 2\n" );

  return;
}
/******************************************************************************/

int f1 ( int i )

/******************************************************************************/
/*
  Purpose:

    F1 is the iteration function for example 1.

  Discussion:

    This function has two cycles:

    6, 3, 1, of length 3
    4, of length 1

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    18 June 2012

  Author:

    John Burkardt

  Parameters:
*/
{
  static int f_table[9] = { 6, 6, 0, 1, 4, 3, 3, 4, 0 };
  int value;

  value = f_table[i];

  return value;
}
/******************************************************************************/

void test02 ( )

/******************************************************************************/
/*
  Purpose:

    TEST02 tests CYCLE_FLOYD for F2.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    14 June 2012

  Author:

    John Burkardt
*/
{
  int lam;
  int mu;
  int x0;

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  Test CYCLE_FLOYD for F2().\n" );
  printf ( "  f2(i) = mod ( 22 * i + 1, 72 ).\n" );

  x0 = 0;
  printf ( "\n" );
  printf ( "  Starting argument X0 = %d\n", x0 );

  cycle_floyd ( f2, x0, &lam, &mu );

  printf ( "\n" );
  printf ( "  Reported cycle length is %d\n", lam );
  printf ( "  Expected value is 9\n" );
  printf ( "\n" );
  printf ( "  Reported distance to first cycle element is %d\n", mu );
  printf ( "  Expected value is 3\n" );

  return;
}
/******************************************************************************/

int f2 ( int i )

/******************************************************************************/
/*
  Purpose:

    F2 is the iteration function for example 2.

  Discussion:

    This function has a cycle

    3, 67, 35, 51, 43, 11, 27, 19, 59, of length 9

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    14 June 2012

  Author:

    John Burkardt

  Parameters:
*/
{
  int value;

  value = ( 22 * i + 1 ) % 72;

  return value;
}
/******************************************************************************/

void test03 ( )

/******************************************************************************/
/*
  Purpose:

    TEST03 tests CYCLE_FLOYD for F3.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    14 June 2012

  Author:

    John Burkardt
*/
{
  int lam;
  int mu;
  int x0;

  printf ( "\n" );
  printf ( "TEST03\n" );
  printf ( "  Test CYCLE_FLOYD for F3().\n" );
  printf ( "  f3(i) = mod ( 123 * i + 456, 100000 ).\n" );

  x0 = 789;
  printf ( "\n" );
  printf ( "  Starting argument X0 = %d\n", x0 );

  cycle_floyd ( f3, x0, &lam, &mu );

  printf ( "\n" );
  printf ( "  Reported cycle length is %d\n", lam );
  printf ( "  Expected value is 50000\n" );
  printf ( "\n" );
  printf ( "  Reported distance to first cycle element is %d\n", mu );
  printf ( "  Expected value is 0\n" );

  return;
}
/******************************************************************************/

int f3 ( int i )

/******************************************************************************/
/*
  Purpose:

    F3 is the iteration function for example 3.

  Discussion:

    This function has a cycle of length 50000

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    18 June 2012

  Author:

    John Burkardt

  Parameters:
*/
{
  int value;

  value = ( 123 * i + 456 ) % 1000000;

  return value;
}
/******************************************************************************/

void test04 ( )

/******************************************************************************/
/*
  Purpose:

    TEST04 tests CYCLE_FLOYD for F4.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    18 June 2012

  Author:

    John Burkardt
*/
{
  int lam;
  int mu;
  int x0;

  printf ( "\n" );
  printf ( "TEST04\n" );
  printf ( "  Test CYCLE_FLOYD for F4().\n" );
  printf ( "  f4(i) = mod ( 31421 * i + 6927, 65536 ).\n" );

  x0 = 1;
  printf ( "\n" );
  printf ( "  Starting argument X0 = %d\n", x0 );

  cycle_floyd ( f4, x0, &lam, &mu );

  printf ( "\n" );
  printf ( "  Reported cycle length is %d\n", lam );
  printf ( "  Expected value is 65536\n" );
  printf ( "\n" );
  printf ( "  Reported distance to first cycle element is %d\n", mu );
  printf ( "  Expected value is 0\n" );

  return;
}
/******************************************************************************/

int f4 ( int i )

/******************************************************************************/
/*
  Purpose:

    F4 is the iteration function for example 4.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    18 June 2012

  Author:

    John Burkardt

  Parameters:
*/
{
  int value;

  value = ( 31421 * i + 6927 ) % 65536;

  return value;
}
/******************************************************************************/

void test05 ( )

/******************************************************************************/
/*
  Purpose:

    TEST05 tests CYCLE_FLOYD for F5.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    18 June 2012

  Author:

    John Burkardt
*/
{
  int i;
  int lam;
  int mu;
  int x0;

  printf ( "\n" );
  printf ( "TEST05\n" );
  printf ( "  Test CYCLE_FLOYD for F5().\n" );
  printf ( "  f5(i) = mod ( 16383 * i + 1, 65536 ).\n" );

  x0 = 1;
  printf ( "\n" );
  printf ( "  Starting argument X0 = %d\n", x0 );

  cycle_floyd ( f5, x0, &lam, &mu );

  printf ( "\n" );
  printf ( "  Reported cycle length is %d\n", lam );
  printf ( "  Expected value is 8\n" );
  printf ( "\n" );
  printf ( "  Reported distance to first cycle element is %d\n", mu );
  printf ( "  Expected value is 0\n" );

  i = 0;
  x0 = 1;
  printf ( "  %d  %d\n", i, x0 );
  for ( i = 1; i <= 10; i++ )
  {
    x0 = f5 ( x0 );
    printf ( "  %d  %d\n", i, x0 );
  }

  return;
}
/******************************************************************************/

int f5 ( int i )

/******************************************************************************/
/*
  Purpose:

    F5 is the iteration function for example 4.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    18 June 2012

  Author:

    John Burkardt

  Parameters:
*/
{
  int value;

  value = ( 16383 * i + 1 ) % 65536;

  return value;
}
