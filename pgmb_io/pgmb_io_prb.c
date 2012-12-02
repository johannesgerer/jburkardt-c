# include <stdlib.h>
# include <stdio.h>
# include <stdbool.h>
# include <time.h>

# include "pgmb_io.h"

int main ( int argc, char *argv[] );
bool test01 ( );
bool test02 ( );
void timestamp ( );

/******************************************************************************/

int main ( int argc, char *argv[] )

/******************************************************************************/
/*
  Purpose:

    PGMB_IO_PRB calls the PGMB_IO test routines.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 June 2012

  Author:

    John Burkardt
*/
{
  bool error;

  timestamp ( );
  printf ( "\n" );
  printf ( "PGMB_IO_PRB:\n" );
  printf ( "  C version\n" );
  printf ( "  Test the PGMB_IO library.\n" );

  error = test01 ( );

  if ( error ) 
  {
    printf ( "\n" );
    printf ( "PGMB_IO_PRB - Fatal error!\n" );
    printf ( "  TEST01 terminated with an error.\n" );
    return 1;
  }

  error = test02 ( );

  if ( error ) 
  {
    printf ( "\n" );
    printf ( "PGMB_IO_PRB - Fatal error!\n" );
    printf ( "  TEST02 terminated with an error.\n" );
    return 1;
  }
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "PGMB_IO_PRB:\n" );
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

    TEST01 tests PGMB_EXAMPLE, PGMB_WRITE.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 June 2012

  Author:

    John Burkardt
*/
{
  bool error;
  char *file_out_name = "pgmb_io_prb_01.pgm";
  unsigned char *g;
  int i;
  unsigned char *indexg;
  int j;
  unsigned char maxg;
  int xsize = 300;
  int ysize = 300;

  printf ( "\n" );
  printf ( "TEST01:\n" );
  printf ( "  PGMB_EXAMPLE sets up PGMB data.\n" );
  printf ( "  PGMB_WRITE writes a PGMB file.\n" );
  printf ( "\n" );
  printf ( "  Writing the file \"%s\".\n", file_out_name );

  g = ( unsigned char * ) malloc ( xsize * ysize * sizeof ( unsigned char ) );

  error = pgmb_example ( xsize, ysize, g );

  if ( error )
  {
    printf ( "\n" );
    printf ( "TEST01 - Fatal error!\n" );
    printf ( "  PGMB_EXAMPLE failed!\n" );
    return error;
  }
  else
  {
    printf ( "\n" );
    printf ( "  PGMB_EXAMPLE has set up the data.\n" );
  }

  maxg = 0;
  indexg = g;

  for ( i = 0; i < xsize; i++ )
  {
    for ( j = 0; j < ysize; j++ )
    {
      if ( maxg < *indexg )
      {
        maxg = *indexg;
      }
      indexg = indexg + 1;
    }
  }
  printf ( "\n" );
  printf ( "  Gray scale data has maximum value %d\n", maxg );

  error = pgmb_write ( file_out_name, xsize, ysize, g );

  if ( error )
  {
    printf ( "\n" );
    printf ( "TEST01 - Fatal error!\n" );
    printf ( "  PGMB_WRITE failed!\n" );
  }
  else
  {
    printf ( "\n" );
    printf (  "  PGMB_WRITE was successful.\n" );
  }

  free ( g );
/*
  Now have PGMB_READ_TEST look at the file we think we created.
*/
  error = pgmb_read_test ( file_out_name );

  return error;
}
/******************************************************************************/

bool test02 ( )

/******************************************************************************/
/*
  Purpose:

    TEST02 tests PGMB_READ_HEADER, PGMB_READ_DATA.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 June 2012

  Author:

    John Burkardt
*/
{
  bool error;
  char *file_in_name = "pgmb_io_prb_02.pgm";
  unsigned char *g;
  unsigned char maxg;
  int xsize;
  int ysize;

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  PGMB_READ reads the header and data of a PGMB file.\n" );
  printf ( "\n" );
  printf ( "  Reading the file \"%s\".\n", file_in_name );
/*
  Create a data file to read.
*/
  error = pgmb_write_test ( file_in_name );

  if ( error )
  {
    printf ( "\n" );
    printf ( "  PGMB_WRITE_TEST failed!\n" );
  }
  else
  {
    printf ( "\n" );
    printf ( "  PGMB_WRITE_TEST created some test data.\n" );
  }
/*
  Now have PGMB_READ try to read it.
*/
  error = pgmb_read ( file_in_name, &xsize, &ysize, &maxg, &g );

  if ( error )
  {
    printf ( "\n" );
    printf ( "  PGMB_READ failed!\n" );
  }
  else
  {
    printf ( "\n" );
    printf ( "  PGMB_READ read the test data successfully.\n" );
  }

  free ( g );

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
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  printf ( "%s\n", time_buffer );

  return;
# undef TIME_SIZE
}
