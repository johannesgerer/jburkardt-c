# include <stdlib.h>
# include <stdio.h>
# include <string.h>
# include <time.h>

int main ( int argc, char** argv );
void test01 ( void );
void test02 ( void );
void test03 ( void );
void timestamp ( void );

/******************************************************************************/

int main ( int argc, char** argv )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for MICE_READER.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    29 April 2013

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "MICE_READER:\n" );
  printf ( "  C version\n" );
  printf ( "  Read a file of strings, one line at a time.\n" );

  test01 ( );
  test02 ( );
  test03 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "MICE_READER:\n" );
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

    TEST01 reads the file using an unusual format string.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    29 April 2013

  Author:

    John Burkardt
*/
{
  char buffer[256];
  FILE *infile;
  int line_num;
  int read;

  printf ( "\n" );
  printf ( "TEST01:\n" );
  printf ( "  Read a line at a time, using fscanf() and a special format.\n" );
  printf ( "  Note that this format misses:\n" );
  printf ( "  Initial and trailing blanks in a string.\n" );
  printf ( "  Any blank lines.\n" );

  infile = fopen ( "mice_file.txt", "rt" );

  printf ( "\n" );
  line_num = 0;

  while ( 1 )
  {
    read = fscanf ( infile, "%[^\n]\n", buffer );
    if ( read < 1 )
    {
      break;
    }
    line_num = line_num + 1;
    printf ( "%d: \"%s\"\n", line_num, buffer );
  }

  fclose ( infile );

  return;
}
/******************************************************************************/

void test02 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST02 reads the file using FGETC().

  Discussion:

    Note that C, which we think of as a character, is actually stored as
    an integer, so that we can allow it to have the special value EOF.

    While FGETC is very reliable, it does read the file at a very low level.
    If we want to recognize words and strings and numbers in the file,
    FGETC requires us to put those together ourselves.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    29 April 2013

  Author:

    John Burkardt
*/
{
  int c;
  int c_old;
  FILE *infile;
  int read;
  int word_num;

  printf ( "\n" );
  printf ( "TEST02:\n" );
  printf ( "  Read a character at a time, using fgetc() and putc().\n" );

  infile = fopen ( "mice_file.txt", "rt" );

  printf ( "\n" );
  word_num = 1;
  c = '\n';

  while ( 1 )
  {
    c_old = c;

    if ( c_old == '\n' )
    {
      printf ( "%d: \"", word_num );
    }

    c = fgetc ( infile );
    
    if ( c == EOF )
    {
      printf ( "EOF\n" );
      break;
    }

    if ( c == '\n' )
    {
      putchar ( '"' );
      word_num = word_num + 1;
    }
    putchar ( c );
  }

  fclose ( infile );

  return;
}
/******************************************************************************/

void test03 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST03 reads the file using a simple string format.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    29 April 2013

  Author:

    John Burkardt
*/
{
  char buffer[256];
  FILE *infile;
  int line_num;
  int read;

  printf ( "\n" );
  printf ( "TEST03:\n" );
  printf ( "  Read using fscanf() and the %%s format.\n" );
  printf ( "  This reads the file one word at a time.\n" );
  printf ( "  It ignores single blanks, double blanks, and blank lines.\n" );

  infile = fopen ( "mice_file.txt", "rt" );

  printf ( "\n" );
  line_num = 0;

  while ( 1 )
  {
    read = fscanf ( infile, "%s", buffer );
    if ( read < 1 )
    {
      break;
    }
    line_num = line_num + 1;
    printf ( "%d: \"%s\"\n", line_num, buffer );
  }

  fclose ( infile );

  return;
}
/******************************************************************************/

void timestamp ( void )

/******************************************************************************/
/*
  Purpose:

    TIMESTAMP prints the current YMDHMS date as a time stamp.

  Example:

    31 May 2001 09:45:54 AM

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    24 September 2003

  Author:

    John Burkardt

  Parameters:

    None
*/
{
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct tm *tm;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  printf ( "%s\n", time_buffer );

  return;
# undef TIME_SIZE
}
