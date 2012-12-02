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
  fprintf ( stdout, "\n" );
  fprintf ( stdout, "IMAGE_EDGE_PRB\n" );
  fprintf ( stdout, "  C version\n" );
  fprintf ( stdout, "  Demonstrate the NEWS stencil for edge detection\n" );
  fprintf ( stdout, "  in images.\n" );

  fprintf ( stdout, "\n" );
  fprintf ( stdout, "  The input file is \"%s\".\n", input_filename );
//
//  Open the input file and read the data.
//
  input_unit = fopen ( input_filename, "rt" );

  pgma_read_header ( input_unit, &m, &n, &g_max );

  fprintf ( stdout, "\n" );
  fprintf ( stdout, "  Number of rows =          %d\n", m );
  fprintf ( stdout, "  Number of columns =       %d\n", n );
  fprintf ( stdout, "  Maximum pixel intensity = %d\n", g_max );

  g = ( int * ) malloc ( m * n * sizeof ( int ) );

  pgma_read_data ( input_unit, m, n, g );

  fclose ( input_unit );

  g_histo = i4mat_histogram ( m, n, g, 255 );

  fprintf ( stdout, "\n" );
  fprintf ( stdout, " Gray     Count\n" );
  fprintf ( stdout, "\n" );
  for ( i = 0; i <= 255; i++ )
  {
    fprintf ( stdout, "  %3d  %8d\n", i, g_histo[i] );
  }

  free ( g_histo );

  e = news ( m, n, g );
//
//  Write the edge information as a portable BIT map (0/1).
//
  pbma_write ( output_filename, m, n, e );

  fprintf ( stdout, "\n" );
  fprintf ( stdout, "  Wrote edge information to \"%s\".\n", output_filename );

  free ( e );
  free ( g );
/*
  Terminate.
*/
  fprintf ( stdout, "\n" );
  fprintf ( stdout, "IMAGE_EDGE_PRB\n" );
  fprintf ( stdout, "  Normal end of execution.\n" );

  fprintf ( stdout, "\n" );
  timestamp ( );

  return 0;
}
