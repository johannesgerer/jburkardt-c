# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "latin_cover.h"

int main ( void );

void test01 ( void );
void test02 ( void );
void test03 ( void );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for LATIN_COVER_PRB.

  Discussion:

    LATIN_COVER_PRB tests the LATIN_COVER library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 August 2012

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "LATIN_COVER_PRB:\n" );
  printf ( "  C version\n" );
  printf ( "  Test the LATIN_COVER library.\n" );

  test01 ( );
  test02 ( );
  test03 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "LATIN_COVER_PRB:\n" );
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

    TEST01 tests LATIN_COVER.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 August 2012

  Author:

    John Burkardt
*/
{
  int *a;
  int n;
  int *p;
  int seed;
  int test;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  LATIN_COVER:\n" );

  for ( n = 3; n <= 9; n = n + 2 )
  {
    seed = 123456789;

    for ( test = 1; test <= 3; test++ )
    {
      p = perm_uniform_new ( n, &seed );
 
      perm_print ( n, p, "  Permutation" );

      a = latin_cover ( n, p );

      i4mat_print ( n, n, a, "  Latin cover" );

      free ( a );
      free ( p );
    }
  }
  return;
}
/******************************************************************************/

void test02 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST02 tests LATIN_COVER_2D.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    24 June 2012

  Author:

    John Burkardt
*/
{
  int *a;
  int n;
  int *p1;
  int *p2;
  int seed;
  int test;

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  LATIN_COVER_2D:\n" );

  for ( n = 3; n <= 9; n = n + 2 )
  {
    seed = 123456789;
    for ( test = 1; test <= 3; test++ )
    {
      p1 = perm_uniform_new ( n, &seed );
      perm_print ( n, p1, "  Permutation 1" );

      p2 = perm_uniform_new ( n, &seed ); 
      perm_print ( n, p2, "  Permutation 2" );

      a = latin_cover_2d ( n, p1, p2 );
      i4mat_print ( n, n, a, "  Latin cover" );

      free ( a );
      free ( p1 );
      free ( p2 );
    }
  }
  return;
}
/******************************************************************************/

void test03 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST03 tests LATIN_COVER_3D.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 August 2012

  Author:

    John Burkardt
*/
{
  int *a;
  int n;
  int *p1;
  int *p2;
  int *p3;
  int seed;
  int test;

  printf ( "\n" );
  printf ( "TEST03\n" );
  printf ( "  LATIN_COVER_3D\n" );

  for ( n = 3; n <= 9; n = n + 2 )
  {
    seed = 123456789;
    for ( test = 1; test <= 3; test++ )
    {
      p1 = perm_uniform_new ( n, &seed );
      perm_print ( n, p1, "  Permutation 1" );

      p2 = perm_uniform_new ( n, &seed );
      perm_print ( n, p2, "  Permutation 2" );

      p3 = perm_uniform_new ( n, &seed );
      perm_print ( n, p3, "  Permutation 1" );

      a = latin_cover_3d ( n, p1, p2, p3 );

      i4block_print ( n, n, n, a, "  Latin cover" );

      free ( a );
      free ( p1 );
      free ( p2 );
      free ( p3 );
    }
  }
  return;
}
