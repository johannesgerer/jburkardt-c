# include <stdlib.h>
# include <stdio.h>
# include <time.h>

# include "rnglib.h"

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );
void test04 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for RNGLIB_PRB.

  Discussion:

    RNGLIB_PRB tests the RNGLIB library.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    28 March 2013

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "RNGLIB_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the RNGLIB library.\n" );
/*
  Call tests.
*/
  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "RNGLIB_PRB\n" );
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

    TEST01 calls I4_UNI 10 times, just to show how it is done.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    05 August 2013

  Author:

    John Burkardt
*/
{
  int g;
  int i;
  int j;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  I4_UNI ( ) returns a random positive integer\n" );
  printf ( "  using the current generator.\n" );
/*
  Initialize the package.
*/
  printf ( "\n" );
  printf ( "  INITIALIZE initializes the random number generator.\n" );
  printf ( "  It only needs to be called once before using the package.\n" );

  initialize ( );
/*
  Set the generator index.
*/
  g = 0;
  cgn_set ( g );
  printf ( "\n" );
  printf ( "  Current generator index = %d\n", g );

  printf ( "\n" );
  printf ( "   I     I4_UNI ( )\n" );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++)
  {
    j = i4_uni ( );
    printf ( "  %2d  %12d\n", i, j );
  }

  return;
}
/******************************************************************************/

void test02 ( )

/******************************************************************************/
/*
  Purpose:

    TEST02 calls R4_UNI_01 10 times, just to show how it is done.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    05 August 2013

  Author:

    John Burkardt
*/
{
  int g;
  int i;
  float u;

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  R4_UNI_01 ( ) returns a random real number\n" );
  printf ( "  in [0,1] using the current generator.\n" );
/*
  Initialize the package.
*/
  printf ( "\n" );
  printf ( "  INITIALIZE initializes the random number generator.\n" );
  printf ( "  It only needs to be called once before using the package.\n" );

  initialize ( );
/*
  Set the current generator index.
*/
  g = 0;
  cgn_set ( g );
  printf ( "\n" );
  printf ( "  Current generator index = %d\n", g );

  printf ( "\n" );
  printf ( "   I     R4_UNI_01 ( )\n" );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    u = r4_uni_01 ( );
    printf ( "  %2d  %14.6g\n", i, u );
  }

  return;
}
/******************************************************************************/

void test03 ( )

/******************************************************************************/
/*
  Purpose:

    TEST03 demonstrates how the seed can be reset to its initial or last value.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    05 August 2013

  Author:

    John Burkardt
*/
{
  int g;
  int i;
  float u;

  printf ( "\n" );
  printf ( "TEST03\n" );
  printf ( "  R4_UNI_01 ( ) returns a random real number\n" );
  printf ( "  in [0,1] using the current generator.\n" );
/*
  Initialize the package.
*/
  printf ( "\n" );
  printf ( "  INITIALIZE initializes the random number generator.\n" );
  printf ( "  It only needs to be called once before using the package.\n" );

  initialize ( );
/*
  Set the current generator index.
*/
  g = 16;
  cgn_set ( g );
  printf ( "\n" );
  printf ( "  Current generator index = %d\n", g );
/*
  Force the generator to start at the initial seed.
*/
  printf ( "\n" );
  printf ( "  INIT_GENERATOR can reset the seed to the initial value,\n" );
  printf ( "  the last (previous) value, or a new seed.\n" );

  printf ( "\n" );
  printf ( "  INIT_GENERATOR ( 0 ) starts at the initial seed.\n" );

  init_generator ( 0 );

  printf ( "\n" );
  printf ( "   I    R4_UNI_01 ( )\n" );
  printf ( "\n" );
  for ( i = 1; i <= 10; i++ )
  {
    u = r4_uni_01 ( );
    printf ( "  %2d  %14.6g\n", i, u );
  }

  printf ( "\n" );
  printf ( "  Calling INIT_GENERATOR ( 0 ) again restarts\n" );
  printf ( "  at the initial seed.\n" );

  init_generator ( 0 );

  printf ( "\n" );
  printf ( "   I    R4_UNI_01 ( )\n" );
  printf ( "\n" );
  for ( i = 1; i <= 10; i++ )
  {
    u = r4_uni_01 ( );
    printf ( "  %2d  %14.6g\n", i, u );
  }

  printf ( "\n" );
  printf ( "  Calling INIT_GENERATOR ( 2 ) restarts\n" );
  printf ( "  at a new \"far ahead\" seed.\n" );

  init_generator ( 2 );

  printf ( "\n" );
  printf ( "   I    R4_UNI_01 ( )\n" );
  printf ( "\n" );
  for ( i = 1; i <= 10; i++ )
  {
    u = r4_uni_01 ( );
    printf ( "  %2d  %14.6g\n", i, u );
  }

  printf ( "\n" );
  printf ( "  Calling INIT_GENERATOR ( 1 ) restarts\n" );
  printf ( "  at the last seed (in this case, the \"far ahead\"\n" );
  printf ( "  seed specified on the previous call.)\n" );

  printf ( "\n" );
  printf ( "   I    R4_UNI_01 ( )\n" );
  printf ( "\n" );
  for ( i = 1; i <= 10; i++ )
  {
    u = r4_uni_01 ( );
    printf ( "  %2d  %14.6g\n", i, u );
    if ( ( i % 3 ) == 0 )
    {
      init_generator ( 1 );
      printf ( "  (Reset to last seed)\n" );
    }
  }

  return;
}
/******************************************************************************/

void test04 ( )

/******************************************************************************/
/*
  Purpose:

    TEST04 demonstrates the use of multiple streams.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    05 August 2013

  Author:

    John Burkardt
*/
{
  int g[3];
  int i;
  int j;
  float u[3];

  printf ( "\n" );
  printf ( "TEST04\n" );
  printf ( "  R4_UNI_01 ( ) returns a random real number\n" );
  printf ( "  in [0,1] using the current generator.\n" );
/*
  Initialize the package.
*/
  printf ( "\n" );
  printf ( "  INITIALIZE initializes the random number generator.\n" );
  printf ( "  It only needs to be called once before using the package.\n" );

  initialize ( );

  printf ( "\n" );
  printf ( "  Let us call generators #2, #5 and #8.\n" );
/*
  Use three separate generators.
  Force them to start at their initial seeds.
*/
  g[0] = 2;
  g[1] = 5;
  g[2] = 8;
  printf ( "\n" );
  for ( i = 0; i < 3; i++ )
  {
    printf ( "  Initialize generator %d\n", g[i] );
    cgn_set ( g[i] );
    init_generator ( 0 );
  }
/*
  Call the generators in the order 2, 5, 8.
*/
  printf ( "\n" );
  printf ( "   I    R4_UNI_01 ( 2 )  R4_UNI_01 ( 5 )  R4_UNI_01 ( 8 )\n" );
  printf ( "\n" );
  for ( i = 1; i <= 10; i++ )
  {
    printf ( "  %2d", i );
    for ( j = 0; j < 3; j++ )
    {
      cgn_set ( g[j] );
      u[j] = r4_uni_01 ( );
      printf ( "  %14.6g", u[j] );
    }
    printf ( "\n" );
  }
/*
  Restart the generators at their initial seeds.
*/
  g[0] = 5;
  g[1] = 8;
  g[2] = 2;
  printf ( "\n" );
  for ( i = 0; i < 3; i++ )
  {
    printf ( "  Reinitialize generator %d\n", g[i] );
    cgn_set ( g[i] );
    init_generator ( 0 );
  }
/*
  Call the generators in the order 5, 8, 2.
*/
  printf ( "\n" );
  printf ( "   I    R4_UNI_01 ( 5 )  R4_UNI_01 ( 8 )  R4_UNI_01 ( 2 )\n" );
  printf ( "\n" );
  for ( i = 1; i <= 10; i++ )
  {
    printf ( "  %2d", i );
    for ( j = 0; j < 3; j++ )
    {
      cgn_set ( g[j] );
      u[j] = r4_uni_01 ( );
      printf ( "  %14.6g", u[j] );
    }
    printf ( "\n" );
  }
  return;
}
