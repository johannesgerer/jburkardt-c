# include <stdlib.h>
# include <stdio.h>

# include "randlc.h"

int main ( void );
void test01 ( void );
void test02 ( void );
void test03 ( void );
void test04 ( void );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    RANDLC_PRB calls sample problems for the RANDLC library.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 March 2010

  Author:

    John Burkardt
*/
{
  timestamp ( );

  printf ( "\n" );
  printf ( "RANDLC_PRB\n" );
  printf ( "  C version:\n" );
  printf ( "  Test the RANDLC library.\n" );

  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "RANDLC_PRB\n" );
  printf ( "  Normal end of execution.\n" );

  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void test01 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST01 tests RANDLC.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    08 March 2010

  Author:

    John Burkardt
*/
{
  int i;
  double seed;
  double seed_init = 123456789.0;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  RANDLC computes pseudorandom values \n" );
  printf ( "  in the interval [0,1].\n" );

  seed = seed_init;

  printf ( "\n" );
  printf ( "  The initial seed is %14.0f\n", seed_init );
  printf ( "\n" );
  printf ( "         I          RANDLC\n" );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    printf ( "  %8d  %14f\n", i, randlc ( &seed ) );
  }

  return;
}
/******************************************************************************/

void test02 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST02 tests RANDLC;

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    08 March 2010

  Author:

    John Burkardt
*/
{
# define N 1000

  int i;
  double seed;
  double seed_in;
  double seed_out;
  double u[N];
  double u_avg;
  double u_var;

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  RANDLC computes a sequence of uniformly distributed\n" );
  printf ( "  pseudorandom numbers.\n" );

  seed = 123456789.0;

  printf ( "\n" );
  printf ( "  Initial SEED = %14.0f\n", seed );

  printf ( "\n" );
  printf ( "  First 10 values:\n" );
  printf ( "\n" );
  printf ( "       I           Input          Output      RANDLC\n" );
  printf ( "                    SEED            SEED\n" );
  printf ( "\n" );

  for ( i = 0; i < 10; i++ )
  {
    seed_in = seed;
    u[i] = randlc ( &seed );
    seed_out = seed;
    printf ( "  %6d  %14.0f  %14.0f  %10f\n", i + 1, seed_in, seed_out, u[i] );
  }

  printf ( "\n" );
  printf ( "  Now call RANDLC %d times.\n", N );

  u_avg = 0.0;
  for ( i = 0; i < N; i++ )
  {
    u[i] = randlc ( &seed );
    u_avg = u_avg + u[i];
  }

  u_avg = u_avg / ( ( double ) N );

  u_var = 0.0;
  for ( i = 0; i < N; i++ )
  {
    u_var = u_var + ( u[i] - u_avg ) * ( u[i] - u_avg );
  }
  u_var = u_var / ( ( double ) ( N - 1 ) );

  printf ( "\n" );
  printf ( "  Average value = %f\n", u_avg );
  printf ( "  Expecting       %f\n", 0.5 );

  printf ( "\n" );
  printf ( "  Variance =      %f\n", u_var );
  printf ( "  Expecting       %f\n", 1.0 / 12.0 );

  return;
# undef N
}
/******************************************************************************/

void test03 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST03 tests RANDLC.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    08 March 2010

  Author:

    John Burkardt
*/
{
  int i;
  double seed;
  double seed_in;
  double seed_out;
  double seed_save;
  double x;

  printf ( "\n" );
  printf ( "TEST03\n" );
  printf ( "  RANDLC computes a sequence of pseudorandom numbers\n" );
  printf ( "  but all computations depend on the seed value.\n" );
  printf ( "  In this test, we show how a sequence of \"random\"\n" );
  printf ( "  values can be manipulated by accessing the seed.\n" );

  seed = 1066.0;

  printf ( "\n" );
  printf ( "  Set SEED to %14.0f\n", seed );
  printf ( "\n" );
  printf ( "  Now call RANDLC 10 times, and watch SEED.\n" );
  printf ( "\n" );
  printf ( "       I           Input          Output      RANDLC\n" );
  printf ( "                    SEED            SEED\n" );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    seed_in = seed;

    if ( i == 5 )
    {
      seed_save = seed;
    }
    x = randlc ( &seed );
    seed_out = seed;
    printf ( "  %6d  %14.0f  %14.0f  %10f\n", i, seed_in, seed_out, x );
  }

  seed = seed_save;

  printf ( "\n" );
  printf ( "  Reset SEED to its value at step 5, = %14.0f\n", seed );
  printf ( "\n" );
  printf ( "  Now call RANDLC 10 times, and watch how SEED\n" );
  printf ( "  and RANDLC restart themselves.\n" );
  printf ( "\n" );
  printf ( "       I           Input          Output      RANDLC\n" );
  printf ( "                    SEED            SEED\n" );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    seed_in = seed;
    x = randlc ( &seed );
    seed_out = seed;
    printf ( "  %6d  %14.0f  %14.0f  %10f\n", i, seed_in, seed_out, x );
  }

  seed = 0.0;

  printf ( "\n" );
  printf ( "  What happens with an initial zero SEED?\n" );
  printf ( "\n" );
  printf ( "       I           Input          Output      RANDLC\n" );
  printf ( "                    SEED            SEED\n" );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    seed_in = seed;
    x = randlc ( &seed );
    seed_out = seed;
    printf ( "  %6d  %14.0f  %14.0f  %10f\n", i, seed_in, seed_out, x );
  }

  seed = -123456789.0;

  printf ( "\n" );
  printf ( "  What happens with an initial negative SEED?\n" );
  printf ( "\n" );
  printf ( "       I           Input          Output      RANDLC\n" );
  printf ( "                    SEED            SEED\n" );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    seed_in = seed;
    x = randlc ( &seed );
    seed_out = seed;
    printf ( "  %6d  %14.0f  %14.0f  %10f\n", i, seed_in, seed_out, x );
  }

  return;
}
/******************************************************************************/

void test04 ( void )

/******************************************************************************/
/*
  Purpose:

    RANDLC_TEST04 tests RANDLC_JUMP.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 March 2010

  Author:

    John Burkardt
*/
{
  int i;
  int k;
  int klog;
  double seed;
  double x1;
  double x2;

  printf ( "\n" );
  printf ( "RANDLC_TEST04\n" );
  printf ( "  RANDLC_JUMP jumps directly to the K-th value\n" );
  printf ( "  returned by RANDLC.\n" );
  printf ( "\n" );
  printf ( "         K X(hard way)     X(jump)\n" );
  printf ( "\n" );

  k = 1;

  for ( klog = 1; klog <= 10; klog++ )
  {
    seed = 123456789.0;
    for ( i = 1; i <= k; i++ )
    {
      x1 = randlc ( &seed );
    }

    seed = 123456789.0;
    x2 = randlc_jump ( seed, k );

    printf ( "  %8d  %10f  %10f\n", k, x1, x2 );

    k = k * 2;
  }

  return;
}
