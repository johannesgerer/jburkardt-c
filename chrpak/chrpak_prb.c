# include <stdlib.h>
# include <stdio.h>
# include <string.h>

# include "chrpak.h"

int main ( void );
void test011 ( void );
void test1126 ( void );
void test119 ( void );
void test137 ( void );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for CHRPAK_PRB.

  Discussion:

    CHRPAK_PRB calls the CHRPAK tests.

  Modified:

    25 January 2010

  Author:

    John Burkardt
*/
{
  timestamp ( );

  printf ( "\n" );
  printf ( "CHRPAK_PRB:\n" );
  printf ( "  C version\n" );
  printf ( "  Test the CHRPAK library.\n" );

  test011 ( );
  test1126 ( );
  test119 ( );
  test137 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "CHRPAK_PRB:\n" );
  printf ( "  Normal end of execution.\n" );

  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void test011 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST011 tests CH_CAP.

  Modified:

    19 January 2007

  Author:

    John Burkardt
*/
{
  char c;

  printf ( "\n" );
  printf ( "TEST011\n" );
  printf ( "  CH_CAP uppercases a character.\n" );
  printf ( "\n" );
  printf ( "  C  CH_CAP(C)\n" );
  printf ( "\n" );

  c = 'F';
  printf ( "  %c  %c\n", c, ch_cap ( c ) );
  c = 'f';
  printf ( "  %c  %c\n", c, ch_cap ( c ) );
  c = '1';
  printf ( "  %c  %c\n", c, ch_cap ( c ) );
  c = 'b';
  printf ( "  %c  %c\n", c, ch_cap ( c ) );
  c = 'B';
  printf ( "  %c  %c\n", c, ch_cap ( c ) );

  return;
}
/******************************************************************************/

void test1126 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST1126 tests S_LEN_TRIM.

  Modified:

    19 January 2007

  Author:

    John Burkardt
*/
{
  int i;
  char s[11];
  int test;

  printf ( "\n" );
  printf ( "TEST1126\n" );
  printf ( "  S_LEN_TRIM reports the length of a string to the last nonblank.\n" );
  printf ( "\n" );
  printf ( "  Here are some strings, and their lengths:\n" );
  printf ( "\n" );

  for ( test = 0; test < 4; test++ )
  {
    if ( test == 0 )
    {
      strcpy ( s, "HELLO" );
    }
    else if ( test == 1 )
    {
      strcpy ( s, "  B la nk" );
    }
    else if ( test == 2 )
    {
      strcpy ( s, " ");
    }
    else if ( test == 3 )
    {
      strcpy ( s, "1234567890" );
    }
    printf ( " \"%s\"  %d\n", s, s_len_trim ( s )  );
  }

  return;
}
/******************************************************************************/

void test119 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST119 tests S_REVERSE.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    25 January 2010

  Author:

    John Burkardt
*/
{
  char *s = "A man, a plan, a canal, Panama!";
  char *s2;

  printf ( "\n" );
  printf ( "TEST119\n" );
  printf ( "  S_REVERSE reverses a string.\n" );
  printf ( "\n" );
  printf ( "  Before: \"%s\".\n", s );
  s2 = s_reverse ( s );
  printf ( "  After:  \"%s\".\n", s2 );
 
  return;
}
/******************************************************************************/

void test137 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST137 tests S_WORD_COUNT.

  Modified:

    19 January 2007

  Author:

    John Burkardt
*/
{
  char s[40];

  printf ( "\n" );
  printf ( "TEST137\n" );
  printf ( "  S_WORD_COUNT counts the words in a string\n" );
  printf ( "\n" );
  printf ( "  STRING                      Words\n" );
  printf ( "\n" );
 
  strcpy ( s, "?" );
  printf ( "  %32s  %23d\n", s, s_word_count ( s ) );

  strcpy ( s, "A man, a plan, a canal - Panama!" );
  printf ( "  %32s  %23d\n", s, s_word_count ( s ) );

  strcpy ( s, " justone!word,-@#$ " );
  printf ( "  %32s  %23d\n", s, s_word_count ( s ) );

  strcpy ( s, "How about a day in the park?" );
  printf ( "  %32s  %23d\n", s, s_word_count ( s ) );

  return;
}
