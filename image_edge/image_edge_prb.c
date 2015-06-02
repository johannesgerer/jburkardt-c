# include <stdlib.h>
# include <stdio.h>
# include <time.h>

# include "image_edge.h"

int main ( int argc, char *argv[] );

/******************************************************************************/

int main ( int argc, char *argv[] )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for IMAGE_EDGE_PRB.

  Discussion:

    IMAGE_EDGE_PRB tests the IMAGE_EDGE library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 June 2010

  Author:

    John Burkardt
*/
{
  int *e;
  int *g;
  int *g_histo;
  int g_max;
  int i;
  char input_filename[80] = "coins.ascii.pgm";
  FILE *input_unit;
  int m;
  int n;
  char output_filename[80] = "coin_edges.ascii.pbm";

  timestamp ( );
  printf ( "\n" );
  printf ( "IMAGE_EDGE_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the IMAGE_EDGE library.\n" );

  printf ( "\n" );
  printf ( "  The input file is \"%s\".\n", input_filename );
//
//  Open the input file and read the data.
//
  input_unit = fopen ( input_filename, "rt" );

  pgma_read_header ( input_unit, &m, &n, &g_max );

  printf ( "\n" );
  printf ( "  Number of rows =          %d\n", m );
  printf ( "  Number of columns =       %d\n", n );
  printf ( "  Maximum pixel intensity = %d\n", g_max );

  g = ( int * ) malloc ( m * n * sizeof ( int ) );

  pgma_read_data ( input_unit, m, n, g );

  fclose ( input_unit );

  g_histo = i4mat_histogram ( m, n, g, 255 );

  printf ( "\n" );
  printf ( " Gray     Count\n" );
  printf ( "\n" );
  for ( i = 0; i <= 255; i++ )
  {
    printf ( "  %3d  %8d\n", i, g_histo[i] );
  }

  free ( g_histo );

  e = news ( m, n, g );
//
//  Write the edge information as a portable BIT map (0/1).
//
  pbma_write ( output_filename, m, n, e );

  printf ( "\n" );
  printf ( "  Wrote edge information to \"%s\".\n", output_filename );

  free ( e );
  free ( g );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "IMAGE_EDGE_PRB\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
