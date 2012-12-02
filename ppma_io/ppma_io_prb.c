# include <stdlib.h>
# include <stdio.h>
# include <time.h>

# include "ppma_io.h"

int main ( int argc, char *argv[] );
void test01 ( void );
void test02 ( void );
void timestamp ( void );

/******************************************************************************/

int main ( int argc, char *argv[] )

/******************************************************************************/
/*
  Purpose:

    PPMA_IO_PRB calls the PPMA_IO test routines.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    27 May 2011

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "PPMA_IO_PRB:\n" );
  printf ( "  C version\n" );
  printf ( "  Test the PPMA_IO library.\n" );

  test01 ( );
  test02 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "PPMA_IO_PRB:\n" );
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

    TEST01 tests PPMA_EXAMPLE and PPMA_WRITE.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    27 May 2011

  Author:

    John Burkardt
*/
{
  int *b;
  int error;
  char *file_out_name = "test01.ascii.ppm";
  int *g;
  int *r;
  int xsize = 300;
  int ysize = 300;

  printf ( "\n" );
  printf ( "TEST01:\n" );
  printf ( "  PPMA_EXAMPLE sets up PPM data.\n" );
  printf ( "  PPMA_WRITE writes an ASCII PPM file.\n" );
  printf ( "\n" );
  printf ( "  Writing the file \"%s\".\n", file_out_name );

  r = ( int * ) malloc ( xsize * ysize * sizeof ( int ) );
  g = ( int * ) malloc ( xsize * ysize * sizeof ( int ) );
  b = ( int * ) malloc ( xsize * ysize * sizeof ( int ) );

  ppma_example ( xsize, ysize, r, g, b );

  printf ( "\n" );
  printf ( "  PPMA_EXAMPLE has set up the data.\n" );

  error = ppma_write ( file_out_name, xsize, ysize, r, g, b );

  if ( error == 1 )
  {
    printf ( "\n" );
    printf ( "TEST01 - Fatal error!\n" );
    printf ( "  PPMA_WRITE failed!\n" );
    return;
  }

  printf ( "\n" );
  printf (  "  PPMA_WRITE was successful.\n" );

  free ( r );
  free ( g );
  free ( b );
/*
  Now have PPMA_READ_TEST look at the file we think we created.
*/
  ppma_read_test ( file_out_name );

  if ( error == 1 )
  {
    printf ( "\n" );
    printf ( "TEST01 - Fatal error!\n" );
    printf ( "  PPMA_READ_TEST failed to read the file we wrote!\n" );
    return;
  }

  printf ( "\n" );
  printf (  "  PPMA_READ_TEST was able to read the file we wrote.\n" );

  return;
}
/******************************************************************************/

void test02 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST02 tests PPMA_READ.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    27 May 2011

  Author:

    John Burkardt
*/
{
  int *b;
  int error;
  char *file_in_name = "test02.ascii.ppm";
  int *g;
  int i;
  int j;
  int k;
  int *r;
  int rgb_max;
  int xsize;
  int ysize;

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  PPMA_READ reads the header and data of a PPMA file.\n" );
  printf ( "\n" );
  printf ( "  Reading the file \"%s\".\n", file_in_name );
/*
  Create a data file to read.
*/
  error = ppma_write_test ( file_in_name );

  if ( error == 1 )
  {
    printf ( "\n" );
    printf ( "TEST02\n" );
    printf ( "  PPMA_WRITE_TEST failed!\n" );
    return;
  }

  printf ( "\n" );
  printf ( "  PPMA_WRITE_TEST created some test data.\n" );
/*
  Now have PPMA_READ try to read it.
*/
  ppma_read ( file_in_name, &xsize, &ysize, &rgb_max, &r, &g, &b );

  printf ( "\n" );
  printf ( "  PPMA_READ read the test data successfully.\n" );
  printf ( "\n" );
  printf ( "  Ten sample values:\n" );
  printf ( "\n" );
  for ( k = 0; k < 10; k++ )
  {
    i = ( ( 9 - k ) * 0 + k * ( xsize - 1 ) ) / 9;
    j = ( ( 9 - k ) * 0 + k * ( ysize - 1 ) ) / 9;
    printf ( "  %i4  %i4  %i6  %i6  %i6\n", i, j, r[i*ysize+j], g[i*ysize+j], g[i*ysize+j] );
  }

  free ( r );
  free ( g );
  free ( b );

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
  size_t len;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  len = strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  printf ( "%s\n", time_buffer );

  return;
# undef TIME_SIZE
}
