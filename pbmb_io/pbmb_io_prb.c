# include <stdlib.h>
# include <stdio.h>
# include <stdbool.h>
# include <time.h>

# include "pbmb_io.h"

int main ( int argc, char *argv[] );
bool test01 ( );
bool test02 ( );
void timestamp ( );

/******************************************************************************/

int main ( int argc, char *argv[] )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for PBMB_IO_PRB.

  Discussion:

    PBMB_IO_PRB calls the PBMB_IO test routines.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 June 2012

  Author:

    John Burkardt
*/
{
  bool error;

  timestamp ( );
  printf ( "\n" );
  printf ( "PBMB_IO_PRB:\n" );
  printf ( "  C version\n" );
  printf ( "  Test the PBMB_IO library.\n" );

  error = test01 ( );

  if ( error )
  {
    printf ( "\n" );
    printf ( "PBMB_IO_PRB - Fatal error!\n" );
    printf ( "  TEST01 terminated with an error.\n" );
    return 1;
  }

  error = test02 ( );

  if ( error )
  {
    printf ( "\n" );
    printf ( "PBMB_IO_PRB - Fatal error!\n" );
    printf ( "  TEST02 terminated with an error.\n" );
    return 1;
  }
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "PBMB_IO_PRB:\n" );
  printf ( "  Normal end of execution.\n" );

  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

bool test01 ( )

/******************************************************************************/
/*
  Purpose:

    TEST01 tests PBMB_EXAMPLE, PBMB_WRITE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 June 2012

  Author:

    John Burkardt
*/
{
  int *b;
  bool error;
  char *file_out_name = "pbmb_io_prb_01.pbm";
  unsigned char *indexb;
  int xsize = 300;
  int ysize = 300;

  printf ( "\n" );
  printf ( "TEST01:\n" );
  printf ( "  PBMB_EXAMPLE sets up PBMB data.\n" );
  printf ( "  PBMB_WRITE writes a PBMB file.\n" );
  printf ( "\n" );
  printf ( "  Writing the file \"%s\".\n", file_out_name );

  b = ( int * ) malloc ( xsize * ysize * sizeof ( int ) );

  error = pbmb_example ( xsize, ysize, b );

  if ( error )
  {
    printf ( "\n" );
    printf ( "TEST01 - Fatal error!\n" );
    printf ( "  PBMB_EXAMPLE failed!\n" );
    return error;
  }
  else
  {
    printf ( "\n" );
    printf ( "  PBMB_EXAMPLE has set up the data.\n" );
  }

  error = pbmb_write ( file_out_name, xsize, ysize, b );

  if ( error )
  {
    printf ( "\n" );
    printf ( "TEST01 - Fatal error!\n" );
    printf ( "  PBMB_WRITE failed!\n" );
  }
  else
  {
    printf ( "\n" );
    printf (  "  PBMB_WRITE was successful.\n" );
  }

  free ( b );
/*
  Now have PBMB_READ_TEST look at the file we think we created.
*/
  error = pbmb_read_test ( file_out_name );

  return error;
}
/******************************************************************************/

bool test02 ( )

/******************************************************************************/
/*
  Purpose:

    TEST02 tests PBMB_READ_HEADER, PBMB_READ_DATA.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 July 2012

  Author:

    John Burkardt
*/
{
  int *b;
  bool error;
  char *file_in_name = "pbmb_io_prb_02.pbm";
  int xsize;
  int ysize;

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  PBMB_READ reads the header and data of a PBMB file.\n" );
  printf ( "\n" );
  printf ( "  Reading the file \"%s\".\n", file_in_name );
/*
  Create a data file to read.
*/
  error = pbmb_write_test ( file_in_name );

  if ( error )
  {
    printf ( "\n" );
    printf ( "  PBMB_WRITE_TEST failed!\n" );
  }
  else
  {
    printf ( "\n" );
    printf ( "  PBMB_WRITE_TEST created some test data.\n" );
  }
/*
  Now have PBMB_READ try to read it.
*/
  error = pbmb_read ( file_in_name, &xsize, &ysize, &b );

  if ( error )
  {
    printf ( "\n" );
    printf ( "  PBMB_READ failed!\n" );
  }
  else
  {
    printf ( "\n" );
    printf ( "  PBMB_READ read the test data successfully.\n" );
  }

  free ( b );

  return error;
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
  size_t len;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  len = strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  fprintf ( stdout, "%s\n", time_buffer );

  return;
# undef TIME_SIZE
}
