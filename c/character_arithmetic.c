# include <stdlib.h>
# include <stdio.h>

int main ( int argc, char *argv[] );

/******************************************************************************/

int main ( int argc, char *argv[] )

/******************************************************************************/
/*
  Purpose:

    CHARACTER_ARITHMETIC does some simple arithmetic with characters.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    22 May 2012

  Author:

    John Burkardt
*/
{
  char c1, c2, c3, c4, c5, c6, c7, c8;

  c1 = 'a';
  c2 = 'b';
  c3 = 'c';
  c4 = 'd';
  c5 = 'z';
  c6 = 'A';
  c7 = '&';
  c8 = '7';
  
  printf ( "\n" );
  printf ( "CHARACTER_ARITHMETIC:\n" );
  printf ( "  C version\n" );
  printf ( "\n" );
  printf ( "  c1 = %c (character) = %d (integer)\n", c1, c1 );
  printf ( "  c2 = %c (character) = %d (integer)\n", c2, c2 );
  printf ( "  c3 = %c (character) = %d (integer)\n", c3, c3 );
  printf ( "  c4 = %c (character) = %d (integer)\n", c4, c4 );
  printf ( "  c5 = %c (character) = %d (integer)\n", c5, c5 );
  printf ( "  c6 = %c (character) = %d (integer)\n", c6, c6 );
  printf ( "  c7 = %c (character) = %d (integer)\n", c7, c7 );
  printf ( "  c8 = %c (character) = %d (integer)\n", c8, c8 );
  
  printf ( "\n" );
  printf ( "We can capitalize a lower case letter by adding A-a to it!\n" );
  printf ( "\n" );
  printf ( "  %c+A-a = %c\n", c1, c1+'A'-'a' );
  printf ( "  %c+A-a = %c\n", c2, c2+'A'-'a' );
  printf ( "  %c+A-a = %c\n", c3, c3+'A'-'a' );
  printf ( "  %c+A-a = %c\n", c4, c4+'A'-'a' );

  printf ( "\n" );
  printf ( "We can average two characters\n" );
  printf ( "\n" );
  printf ( "  1/2(%c+%c)=%c\n", c1, c3, ( c1 + c3 ) / 2 );
  printf ( "  1/2(%c+%c)=%c\n", c1, c5, ( c1 + c5 ) / 2 );
  return 0;
}
