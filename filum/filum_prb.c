# include <stdlib.h>
# include <stdio.h>

# include "filum.h"

int main ( void );
void test03 ( void );
void test06 ( void );
void test14 ( void );
void test22 ( void );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for FILUM_PRB.

  Discussion:

    FILUM_PRB tests the FILUM library.

  Modified:

    22 November 2011

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "FILUM_PRB:\n" );
  printf ( "  C version\n" );
  printf ( "  Test the FILUM library.\n" );

  test03 ( );
  test06 ( );

  test14 ( );

  test22 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "FILUM_PRB:\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void test03 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST03 tests FILE_COLUMN_COUNT.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    28 November 2009

  Author:

    John Burkardt 
*/
{
  int column_num;
  char file_name[] = "filum_prb_4by5.txt";

  printf ( "\n" );
  printf ( "TEST03\n" );
  printf ( "  FILE_COLUMN_COUNT counts the columns in a file.\n" );
  printf ( "\n" );
  printf ( "  It is assumed that the file contains a number of lines,\n" );
  printf ( "  with each line containing the same number of words.\n" );
  printf ( "  The task is to determine the number of words in a line,\n" );
  printf ( "  that is, the number of \"columns\" of text.\n" );

  printf ( "\n" );
  printf ( "  Examining the file:\"%s\".\n", file_name );

  column_num = file_column_count ( file_name );

  printf ( "\n" );
  printf ( "  Number of columns: %d\n", column_num );

  return;
}
/******************************************************************************/

void test06 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST006 tests FILE_EXIST.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    28 November 2009

  Author:

    John Burkardt
*/
{
  char filename1[] = "filum_prb.c";
  char filename2[] = "filum.c";
  char filename3[] = "raisin.txt";
  char filename4[] = "make.money.fast";

  printf ( "\n" );
  printf ( "TEST06\n" );
  printf ( "  FILE_EXIST reports whether a file 'exists'.\n" );
  printf ( "\n" );
  printf ( "  Exist?   File_name\n" );
  printf ( "\n" );
  printf ( "       %d  %s\n", file_exist ( filename1 ), filename1 );
  printf ( "       %d  %s\n", file_exist ( filename2 ), filename2 );
  printf ( "       %d  %s\n", file_exist ( filename3 ), filename3 );
  printf ( "       %d  %s\n", file_exist ( filename4 ), filename4 );

  return;
}
/******************************************************************************/

void test14 ( )

/******************************************************************************/
/*
  Purpose:

    TEST14 tests FILENAME_INC.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    22 November 2011

  Author:

    John Burkardt
*/
{
  int i;
  int j;
  char filename[80];

  printf ( "\n" );
  printf ( "TEST14\n" );
  printf ( "  FILENAME_INC increments a string\n" );
  printf ( "\n" );
  printf ( "     Input             Output\n" );

  for ( i = 0; i < 4; i++ )
  {
    if ( i == 0 )
    {
      strcpy ( filename, "file???.dat" );
    }
    else if ( i == 1 ) 
    {
      strcpy ( filename, "file072.dat" );
    }
    else if ( i == 2 ) 
    {
      strcpy ( filename, "2cat9.dat  " );
    }
    else if ( i == 3 ) 
    {
      strcpy ( filename, "fred98.txt " );
    }
    printf ( "\n" );
    for ( j = 1; j <= 4; j++ )
    {
      printf ( "  %11s  ", filename );

      filename_inc ( filename );

      printf ( "  %11s\n", filename );

      if ( s_len_trim ( filename ) <= 0 )
      {
        printf ( "  (File name not incrementable.  Quit loop!)\n" );
        break;
      }
    }
  }

  return;
}
/******************************************************************************/

void test22 ( )

/******************************************************************************/
/*
  Purpose:

    TEST22 tests FILE_ROW_COUNT.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    28 November 2009

  Author:

    John Burkardt
*/
{
  char filename[] = "filum_prb_test.txt";

  printf ( "\n" );
  printf ( "TEST22\n" );
  printf ( "  FILE_ROW_COUNT counts the lines in a file.\n" );
  printf ( "\n" );
  printf ( "  Examining file \"%s\".\n", filename );
  printf ( "\n" );
  printf ( "  Number of lines: %d\n", file_row_count ( filename ) );

  return;
}
