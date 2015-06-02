# include <stdlib.h>
# include <stdio.h>

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    SCANF_DEMO shows how scanf() can be used to read input.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    23 May 2012

  Author:

    John Burkardt
*/
{
  char c;
  int n;
  char *s;
  float x;
  double y;

  printf ( "\n" );
  printf ( "SCANF_DEMO\n" );
  printf ( "  C version\n" );
  printf ( "  Show how scanf(format,address) is used to read data from the user.\n" );

  printf ( "\n" );
  printf ( "Enter an integer value for int n:\n" );
  scanf ( "%i", &n );
  printf ( "Enter a real value for float x:\n" );
  scanf ( "%f", &x );
/*
  To read a double, you cannot use the %f format.
  You must use "%lf".
*/
  printf ( "Enter a real value for double y:\n" );
  scanf ( "%lf", &y );
/*
  To read a character with scanf, it's important that
  the format string " %c" includes a blank before %c.
  Otherwise, C will read the carriage return left over from 
  your previous input line.
*/
  printf ( "Enter 1 character for char c:\n" );
  scanf ( " %c", &c );
/*
  A string is a peculiar quantity.  
  It's really a pointer to a list of characters, the last of which is
  a NULL.  To store 20 "useful" characters requires us to reserve 21 spaces.
  Also, since s is "already" a pointer, we just pass "s" to scanf,
  not "&s".
  Finally, note that scanf will only read the string up to the first
  blank space encountered.  So if we enter "Hi, mom!", scanf will
  only grab "Hi," to copy into s.
*/
  s = ( char * ) malloc ( 21 * sizeof ( char ) );
  printf ( "Enter a string of 20 characters or less and NO INTERNAL SPACES:\n" );
  scanf ( "%20s", s );
/*
  Let's see what we got:
*/
  printf ( "\n" );
  printf ( "  Here is what scanf() got for us:\n" );
  printf ( "\n" );
  printf ( "  n = %i\n", n );
  printf ( "  x = %f\n", x );
  printf ( "  y = %f\n", y );
  printf ( "  c = %c\n", c );
  printf ( "  s = %s\n", s );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "SCANF_DEMO\n" );
  printf ( "  Normal end of execution.\n" );

  return 0;
}
