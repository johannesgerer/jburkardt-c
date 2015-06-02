# include <stdlib.h>
# include <stdio.h>
# include <string.h>

# include "file_name_sequence.h"

int main ( );
void test02 ( char *prefix, char *suffix, int first, int last );
void test03 ( char *prefix, char *suffix, int first, int last );
void test04 ( char *filename, int filename_num );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for FILE_NAME_SEQUENCE_PRB.

  Discussion:

    FILE_NAME_SEQUENCE_PRB tests the FILE_NAME_SEQUENCE library.

    There are situations such as animations or parallel processing in which
    it is necessary to generate a sequence of file names which include
    an embedded index that increases.  A simple example might be

      "fred0.txt", "fred1.txt", "fred2.txt"

    A side issue arises when the number of files is large enough that the
    number of digits in the index will vary.  Thus, if we are going to have
    15 files, do we want to number them as

      "fred00.txt" through "fred14.txt"

    which means, for one thing, that they will alphabetize properly, or
    will we be satisfied with

      "fred0.txt" through "fred14.txt" ?

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    23 July 2012

  Author:

    John Burkardt
*/
{
  char filename[81] = { "frodo_01345_lives.txt" };

  timestamp ( );
  printf ( "\n" );
  printf ( "FILE_NAME_SEQUENCE_PRB:\n" );
  printf ( "  C version\n" );
  printf ( "  Test the FILE_NAME_SEQUENCE library.\n" );
  printf ( "  Demonstrate ways of generating a numeric sequence of file names.\n" );

  test02 ( "fred", ".txt", 0, 12 );
  test03 ( "frid", ".txt", 99, 105 );
  test04 ( filename, 10 );
//
//  Terminate.
//
  printf ( "\n" );
  printf ( "FILE_NAME_SEQUENCE_PRB:\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
//****************************************************************************80

void test02 ( char *prefix, char *suffix, int first, int last )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 uses SPRINTF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 July 2012
//
//  Author:
//
//    John Burkardt
//
{
  char filename[81];
  int i;

  printf ( "\n" );
  printf ( "TEST02:\n" );
  printf ( "  FILENAME(I) = SPRINTF ( PREFIX, I, SUFFIX )\n" );
  printf ( "  PREFIX = \"%s\"\n", prefix );
  printf ( "  SUFFIX = \"%s\"\n", suffix );
  printf ( "  %d <= I <= %d\n", first, last );
  printf ( "  Numbers do NOT include leading zeros.\n" );
  printf ( "\n" );

  for ( i = first; i <= last; i++ )
  {
    sprintf ( filename, "%s%d%s", prefix, i, suffix );
    printf ( "  %4d:  \"%s\"\n", i, filename );
  }

  return;
}
//****************************************************************************80

void test03 ( char *prefix, char *suffix, int first, int last )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 uses SPRINTF and leading 0's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 July 2012
//
//  Author:
//
//    John Burkardt
//
{
  char filename[81];
  int i;

  printf ( "\n" );
  printf ( "TEST03:\n" );
  printf ( "  FILENAME(I) = SPRINTF ( PREFIX, I, SUFFIX )\n" );
  printf ( "  PREFIX = \"%s\"\n", prefix );
  printf ( "  SUFFIX = \"%s\"\n", suffix );
  printf ( "  %d <= I <= %d\n", first, last );
  printf ( "  Numbers DO include leading zeros.\n" );
  printf ( "\n" );

  for ( i = first; i <= last; i++ )
  {
    sprintf ( filename, "%s%04d%s", prefix, i, suffix );
    printf ( "  %4d:  \"%s\"\n", i, filename );
  }
  return;
}
//****************************************************************************80

void test04 ( char *filename, int filename_num )

//****************************************************************************80
//
//  Purpose:
//
//    TEST04 uses FILENAME_INC.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 July 2012
//
//  Author:
//
//    John Burkardt
//
{
  int i;

  printf ( "\n" );
  printf ( "TEST04:\n" );
  printf ( "  FILENAME(I+1) = FILENAME_INC ( FILENAME(I) )\n" );
  printf ( "  First FILENAME = \"%s\"\n", filename );
  printf ( "  Number of filenames = %d\n", filename_num );
  printf ( "  Numbers may include leading zeros.\n" );
  printf ( "\n" );

  for ( i = 1; i <= filename_num; i++ )
  {
    printf ( "  %4d:  \"%s\"\n", i, filename );
    filename_inc ( filename );
  }

  return;
}
