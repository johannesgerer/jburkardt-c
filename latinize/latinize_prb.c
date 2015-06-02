# include <stdlib.h>
# include <stdio.h>
# include <time.h>
# include <string.h>

# include "latinize.h"

int main ( int argc, char *argv[] );
void test01 ( char *input_filename );

/******************************************************************************/

int main ( int argc, char *argv[] )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for LATINIZE_PRB.

  Discussion:

    LATINIZE_PRB tests the LATINIZE library.

    The dataset is presumed to be an M by N array of real numbers,
    where M is the spatial dimension, and N is the number of sample points.

    The dataset is presumed to be stored in a file, with N records,
    one per each sample point.  (Comment records may be included, 
    which begin with '#'.)

    The program reads the data file, "latinizes" the data, and writes
    the latinized data to a new file.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    04 February 2012

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "LATINIZE_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the LATINIZE library.\n" );
  printf ( "\n" );
  printf ( "  Read a dataset of N points in M dimensions,\n" );
  printf ( "  modify it into a Latin hypercube,\n" );
  printf ( "  write the modified dataset to a file.\n" );

  test01 ( "cvt_02_00010.txt" );
  test01 ( "cvt_03_00007.txt" );
  test01 ( "cvt_03_00056.txt" );
  test01 ( "cvt_07_00100.txt" );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "LATINIZE_PRB\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void test01 ( char *input_filename )

/******************************************************************************/
/*
  Purpose:

    TEST01 latinizes a given file.

  Discussion:

    The dataset is presumed to be an M by N array of real numbers,
    where M is the spatial dimension, and N is the number of sample points.

    The dataset is presumed to be stored in a file, with N records,
    one per each sample point.  (Comment records may be included, 
    which begin with '#'.)

    The program reads the data file, "latinizes" the data, and writes
    the latinized data to a new file.

  Modified:

    08 October 2004

  Author:

    John Burkardt
*/
{
  int m;
  int n;
  FILE *output;
  char *output_filename;
  double *table;
/*
  Need to create the output file name from the input filename.
*/
  output_filename = file_name_ext_swap ( input_filename, "latin.txt" );

  r8mat_header_read ( input_filename, &m, &n );

  printf ( "\n" );
  printf ( "  Read the header of \"%s\".\n", input_filename );
  printf ( "\n" );
  printf ( "  Spatial dimension M = %d\n", m );
  printf ( "  Number of points N  = %d\n", n );

  table = r8mat_data_read ( input_filename, m, n );

  printf ( "\n" );
  printf ( "  Read the data in \"%s\".\n", input_filename );

  r8mat_transpose_print_some ( m, n, table, 1, 1, 5, 5, 
    "  Small portion of data read from file:" );

  r8mat_latinize ( m, n, table );

  printf ( "\n" );
  printf ( "  Latinized the data.\n" );

  r8mat_transpose_print_some ( m, n, table, 1, 1, 5, 5, 
    "  Small portion of Latinized data:" );
/*
  Write the data to a file.
*/
  r8mat_write ( output_filename, m, n, table );

  printf ( "\n" );
  printf ( "  Wrote the latinized data to \"%s\".\n", output_filename );

  free ( table );

  return;
}
