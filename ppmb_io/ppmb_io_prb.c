# include <stdlib.h>
# include <stdio.h>
# include <stdbool.h>
# include <time.h>

# include "ppmb_io.h"

int main ( int argc, char *argv[] );
bool test01 ( );
bool test02 ( );
void timestamp ( );

/******************************************************************************/

int main ( int argc, char *argv[] )

/******************************************************************************/
/*
  Purpose:

    MAIN calls the PPMB_IO test routines.

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
  printf ( "PPMB_IO_PRB:\n" );
  printf ( "  C version\n" );
  printf ( "  Test the PPMB_IO library.\n" );

  error = test01 ( );
  error = test02 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "PPMB_IO_PRB:\n" );
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

    TEST01 tests PPM_EXAMPLE, PPMB_WRITE.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 June 2012

  Author:

    John Burkardt
*/
{
  unsigned char *b;
  bool error;
  char *file_name = "ppmb_io_prb_01.ppm";
  unsigned char *g;
  unsigned char *r;
  int xsize = 300;
  int ysize = 300;

  printf ( "\n" );
  printf ( "TEST01:\n" );
  printf ( "  PPMB_EXAMPLE sets up PPM data.\n" );
  printf ( "  PPMB_WRITE writes a binary PPM file.\n" );

  r = ( unsigned char * ) malloc ( xsize * ysize * sizeof ( unsigned char ) );
  g = ( unsigned char * ) malloc ( xsize * ysize * sizeof ( unsigned char ) );
  b = ( unsigned char * ) malloc ( xsize * ysize * sizeof ( unsigned char ) );

  error = ppmb_example ( xsize, ysize, r, g, b );

  if ( error )
  {
    printf ( "\n" );
    printf ( "TEST01 - Fatal error!\n" );
    printf ( "  PPMB_EXAMPLE failed!\n" );
    return error;
  }
  else
  {
    printf ( "\n" );
    printf ( "  PGMB_EXAMPLE has set up the data.\n" );
  }

  error = ppmb_write ( file_name, xsize, ysize, r, g, b );

  if ( error )
  {
    printf ( "  PPMB_WRITE failed!\n" );
  }
  else
  {
    printf ( "  PPMB_WRITE was successful.\n" );
  }

  free ( r );
  free ( g );
  free ( b );

  return error;
}
/******************************************************************************/

bool test02 ( )

/******************************************************************************/
/*
  Purpose:

    TEST02 tests PPMB_READ_HEADER, PPMB_READ_DATA.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    22 July 2011

  Author:

    John Burkardt
*/
{
  unsigned char *b;
  bool error;
  FILE *file_in;
  char *file_in_name = "ppmb_io_prb_02.ppm";
  unsigned char *g;
  int maxrgb;
  unsigned char *r;
  int xsize;
  int ysize;

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  PPMB_READ reads the header and data of a PPMB file.\n" );
 
  printf ( "\n" );
  printf ( "  Reading the file %s.\n", file_in_name );
/*
  Create a data file to read.
*/
  error = ppmb_write_test ( file_in_name );

  if ( error )
  {
    printf ( "\n" );
    printf ( "  PPMB_WRITE_TEST failed!\n" );
    return error;
  }
  else
  {
    printf ( "\n" );
    printf ( "  PPMB_WRITE_TEST created some test data.\n" );
  }
/*
  Now have PPMB_READ try to read it.
*/
  error = ppmb_read ( file_in_name, &xsize, &ysize, &maxrgb, &r, &g, &b );

  if ( error )
  {
    printf ( "\n" );
    printf ( "  PPMB_READ failed!\n" );
    return error;
  }
  else
  {
    printf ( "\n" );
    printf ( "  PPMB_READ read the test data successfully.\n" );
  }

  printf ( "\n" );
  printf ( "  The header was read successfully.\n" );
  printf ( "  Number of rows of data    = %d\n", xsize );
  printf ( "  Number of columns of data = %d\n", ysize );
  printf ( "  Maximum data value =        %d\n", maxrgb );

  free ( r );
  free ( g );
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
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  printf ( "%s\n", time_buffer );

  return;
# undef TIME_SIZE
}
