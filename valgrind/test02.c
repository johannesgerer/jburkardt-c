# include <stdlib.h>
# include <stdio.h>

void junk_data ( void );
int main ( void );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for TEST01.

  Discussion:

    TEST02 has some uninitialized data.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    19 May 2011
*/
{
  printf ( "\n" );
  printf ( "TEST02:\n" );
  printf ( "  C version\n" );
  printf ( "  A sample code for analysis by VALGRIND.\n" );

  junk_data ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  Normal end of execution.\n" );

  return 0;
}
/******************************************************************************/

void junk_data ( void )

/******************************************************************************/
/*
  Purpose:

    JUNK_DATA has some uninitialized variables.

  Discussion:

    VALGRIND's MEMCHECK program monitors uninitialized variables, but does
    not complain unless such a variable is used in a way that means its
    value affects the program's results, that is, the value is printed,
    or computed with.  Simply copying the unitialized data to another variable
    is of no concern.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    18 May 2011
*/
{
  int i;
  int *x;

  x = ( int * ) malloc ( 10 * sizeof ( int ) );
/*
  X = { 0, 1, 2, 3, 4, ?a, ?b, ?c, ?d, ?e }.
*/
  for ( i = 0; i < 5; i++ )
  {
    x[i] = i;
  }
/*
  Copy some values.
  X = { 0, 1, ?c, 3, 4, ?b, ?b, ?c, ?d, ?e }.
*/
  x[2] = x[7];
  x[5] = x[6];
/*
  Modify some uninitialized entries.
  Memcheck doesn't seem to care about this.
*/
  for ( i = 0; i < 10; i++ )
  {
    x[i] = 2 * x[i];
  }
/*
  Print X.
*/
  for ( i = 0; i < 10; i++ )
  {
    printf ( "  %d  %d\n", i, x[i] );
  }

  free ( x );

  return;
}
