# include <stdlib.h>
# include <stdio.h>
# include <string.h>
# include <time.h>

# include "pgma_io.h"

int main ( int argc, char *argv[] );
void test01 ( void )  ;
void test02 ( void )  ;
void test03 ( void )  ;
void timestamp ( void )  ;

/******************************************************************************/

int main ( int argc, char *argv[] )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for PGMA_IO_PRB.

  Discussion:

    PGMA_IO_PRB tests the PGMA_IO library.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    05 June 2010

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "PGMA_IO_PRB:\n" );
  printf ( "  C version\n" );
  printf ( "  Test the PGMA_IO library.\n" );

  test01 ( );
  test02 ( );
  test03 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "PGMA_IO_PRB:\n" );
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

    TEST01 tests PGMA_EXAMPLE, PGMA_WRITE.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    04 June 2010

  Author:

    John Burkardt
*/
{
  char file_out_name[80] = "pgma_io_prb_01.ascii.pgm";
  int *g;
  int xsize = 300;
  int ysize = 300;

  fprintf ( stdout, "\n" );
  fprintf ( stdout, "TEST01:\n" );
  fprintf ( stdout, "  PGMA_EXAMPLE sets up ASCII PGM data.\n" );
  fprintf ( stdout, "  PGMA_WRITE writes an ASCII PGM file.\n" );
  fprintf ( stdout, "\n" );
  fprintf ( stdout, "  Writing the file \"%s\".\n", file_out_name );

  g = ( int * ) malloc ( xsize * ysize * sizeof ( int ) );

  pgma_example ( xsize, ysize, g );

  fprintf ( stdout, "\n" );
  fprintf ( stdout, "  PGMA_EXAMPLE has set up the data.\n" );

  pgma_write ( file_out_name, xsize, ysize, g );

  fprintf ( stdout, "\n" );
  fprintf ( stdout,  "  PGMA_WRITE was successful.\n" );

  free ( g );
/*
  Now have PGMA_READ_TEST look at the file we think we created.
*/
  pgma_read_test ( file_out_name );

  fprintf ( stdout, "\n" );
  fprintf ( stdout,  "  PGMA_READ_TEST was able to read our file.\n" );

  return;
}
/******************************************************************************/

void test02 ( void )  

/******************************************************************************/
/*
  Purpose:

    TEST02 tests PGMA_READ_HEADER, PGMA_READ_DATA.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    04 June 2010

  Author:

    John Burkardt
*/
{
  char file_in_name[80] = "pgma_io_prb_02.ascii.pgm";
  int *g;
  int i;
  int j;
  int k;
  int maxg;
  int xsize;
  int ysize;

  fprintf ( stdout, "\n" );
  fprintf ( stdout, "TEST02\n" );
  fprintf ( stdout, "  PGMA_READ reads the header and data of an ASCII PGM file.\n" );
  fprintf ( stdout, "\n" );
  fprintf ( stdout, "  Reading the file \"%s\".\n", file_in_name );
/*
  Create a data file to read.
*/
  pgma_write_test ( file_in_name );
  
  fprintf ( stdout, "\n" );
  fprintf ( stdout, "  PGMA_WRITE_TEST created some test data.\n" );
/*
  Now have PGMA_READ try to read it.
*/
  pgma_read ( file_in_name, &xsize, &ysize, &maxg, &g );

  fprintf ( stdout, "\n" );
  fprintf ( stdout, "  PGMA_READ read the test data successfully.\n" );

  fprintf ( stdout, "\n" );
  fprintf ( stdout, "  Sample data:\n" );
  fprintf ( stdout, "\n" );
  for ( k = 0; k <= 9; k++ )
  {
    i = ( ( 9 - k ) * 0 + k * ( xsize - 1 ) ) / 9;
    j = ( ( 9 - k ) * 0 + k * ( ysize - 1 ) ) / 9;
    fprintf ( stdout, "%4d  %4d  %6d\n", i, j, g[i*ysize+j] );
  }

  free ( g );

  return;
}
/******************************************************************************/

void test03 ( void )  

/******************************************************************************/
/*
  Purpose:

    TEST03 tests PGMA_WRITE.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    04 June 2010

  Author:

    John Burkardt
*/
{
# define NGRAY 11

  char file_out_name[80] = "pgma_io_prb_03.ascii.pgm";
  int *g;
  double gray[NGRAY] = { 
    0.000, 0.291, 0.434, 0.540, 0.629,
    0.706, 0.774, 0.837, 0.895, 0.949,
    1.000 };
  int i;
  int j;
  int k;
  int xsize = 300;
  int ysize = 300;

  fprintf ( stdout, "\n" );
  fprintf ( stdout, "TEST03:\n" );
  fprintf ( stdout, "  PGMA_WRITE writes an ASCII PGM file.\n" );
  fprintf ( stdout, "\n" );
  fprintf ( stdout, "  In this example, we make a sort of grayscale\n" );
  fprintf ( stdout, "  checkerboard.\n" );

  g = ( int * ) malloc ( xsize * ysize * sizeof ( int ) );

  for ( i = 0; i < xsize; i++ )
  {
    for ( j = 0; j < ysize; j++ )
    {
      k = ( i + j ) * NGRAY / i4_min ( xsize, ysize );
      k = k % NGRAY;
      g[i*ysize+j] = ( int ) ( 255.0E+00 * gray[k] );
    }
  }

  fprintf ( stdout, "  Writing the file \"%s\".\n", file_out_name );

  pgma_write ( file_out_name, xsize, ysize, g );

  fprintf ( stdout, "\n" );
  fprintf ( stdout,  "  PGMA_WRITE was successful.\n" );

  free ( g );

  return;
# undef NGRAY
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
