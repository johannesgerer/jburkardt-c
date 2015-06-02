# include <stdlib.h>
# include <stdio.h>
# include <time.h>

# include "pbma_io.h"

int main ( int argc, char *argv[] );
void test01 ( void );
void test02 ( void );
void timestamp ( void );

/******************************************************************************/

int main ( int argc, char *argv[] )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for PBMA_IO_PRB.

  Discussion:

    PBMA_IO_PRB tests the PBMA_IO library.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    05 June 2010

  Author:

    John Burkardt
*/
{
  timestamp ( );
  fprintf ( stdout, "\n" );
  fprintf ( stdout, "PBMA_IO_PRB:\n" );
  fprintf ( stdout, "  C version\n" );
  fprintf ( stdout, "  Test the PBMA_IO library.\n" );

  test01 ( );
  test02 ( );
/*
  Terminate.
*/
  fprintf ( stdout, "\n" );
  fprintf ( stdout, "PBMA_IO_PRB:\n" );
  fprintf ( stdout, "  Normal end of execution.\n" );
  fprintf ( stdout, "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void test01 ( )

/******************************************************************************/
/*
  Purpose:

    TEST01 tests PBMA_EXAMPLE, PBMA_WRITE.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    05 June 2010

  Author:

    John Burkardt
*/
{
  int *b;
  char file_out_name[80] = "pbma_io_prb_01.ascii.pbm";
  int xsize = 300;
  int ysize = 300;

  fprintf ( stdout, "\n" );
  fprintf ( stdout, "TEST01:\n" );
  fprintf ( stdout, "  PBMA_EXAMPLE sets up ASCII PBM data.\n" );
  fprintf ( stdout, "  PBMA_WRITE writes an ASCII PBM file.\n" );
  fprintf ( stdout, "\n" );
  fprintf ( stdout, "  Writing the file \"%s\".\n", file_out_name );

  b = ( int * ) malloc ( xsize * ysize * sizeof ( int ) );

  pbma_example ( xsize, ysize, b );

  fprintf ( stdout, "\n" );
  fprintf ( stdout, "  PBMA_EXAMPLE has set up the data.\n" );

  pbma_write ( file_out_name, xsize, ysize, b );

  fprintf ( stdout, "\n" );
  fprintf ( stdout,  "  PBMA_WRITE was successful.\n" );

  free ( b );
/*
  Now have PBMA_READ_TEST look at the file we think we created.
*/
  pbma_read_test ( file_out_name );

  fprintf ( stdout, "\n" );
  fprintf ( stdout,  "  PBMA_READ_TEST was able to read the file.\n" );

  return;
}
/******************************************************************************/

void test02 ( )

/******************************************************************************/
/*
  Purpose:

    TEST02 tests PBMA_READ_HEADER, PBMA_READ_DATA.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    05 June 2010

  Author:

    John Burkardt
*/
{
  int *b;
  char file_in_name[80] = "pbma_io_prb_02.ascii.pbm";
  int i;
  int j;
  int k;
  int xsize;
  int ysize;

  fprintf ( stdout, "\n" );
  fprintf ( stdout, "TEST02\n" );
  fprintf ( stdout, "  PBMA_READ reads the header and data of an ASCII PBM file.\n" );
/*
  Create a data file to read.
*/
  pbma_write_test ( file_in_name );

  fprintf ( stdout, "\n" );
  fprintf ( stdout, "  PBMA_WRITE_TEST created some test data.\n" );
/*
  Now have PBMA_READ try to read it.
*/
  pbma_read ( file_in_name, &xsize, &ysize, &b );

  fprintf ( stdout, "\n" );
  fprintf ( stdout, "  PBMA_READ was able to read the file we created.\n" );
  fprintf ( stdout, "\n" );
  fprintf ( stdout, "  Sample data:\n" );
  fprintf ( stdout, "\n" );

  for ( k = 0; k <= 29; k++ )
  {
    i = ( ( 29 - k ) * 0 + k * ( xsize - 1 ) ) / 29;
    j = ( ( 29 - k ) * 0 + k * ( ysize - 1 ) ) / 29;
    fprintf ( stdout, "%4d  %4d  %6d\n", i, j, b[i*ysize+j] );
  }

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
