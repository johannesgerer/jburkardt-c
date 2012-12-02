# include <stdlib.h>
# include <stdio.h>

# include "image_denoise.h"

int main ( int argc, char *argv[] );
void test01 ( void );

/******************************************************************************/

int main ( int argc, char *argv[] )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for IMAGE_DENOISE_PRB.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 February 2012

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "IMAGE_DENOISE_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the IMAGE_DENOISE library.\n" );

  test01 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "IMAGE_DENOISE_PRB\n" );
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

    TEST01 tests GRAY_MEDIAN_NEWS.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 February 2012

  Author:

    John Burkardt
*/
{
  int *g;
  int g_max;
  int *g2;
  char input_filename[80] = "glassware_noisy.ascii.pgm";
  FILE *input_unit;
  int m;
  int n;
  char output_filename[80] = "glassware_median_news.ascii.pgm";

  printf ( "\n" );
  printf ( "TEST01:\n" );
  printf ( "  GRAY_MEDIAN_NEWS uses a NEWS median filter \n" );
  printf ( "  on a noisy grayscale image.\n" );

  printf ( "\n" );
  printf ( "  The input file is \"%s\".\n", input_filename );
/*
  Open the input file and read the data.
*/
  input_unit = fopen ( input_filename, "r" );

  if ( !input_unit )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "TEST01 - Fatal error!\n" );
    fprintf ( stderr, "  Could not open the file \"%s\"\n", input_filename );
    exit ( 1 );
  }

  pgma_read_header ( input_unit, &m, &n, &g_max );

  printf ( "\n" );
  printf ( "  Number of rows =          %d\n", m );
  printf ( "  Number of columns =       %d\n", n );
  printf ( "  Maximum pixel intensity = %d\n", g_max );

  g = ( int * ) malloc ( m * n * sizeof ( int ) );

  pgma_read_data ( input_unit, m, n, g );

  fclose ( input_unit );

  g2 = gray_median_news ( m, n, g );
/*
  Write the denoised images.
*/
  pgma_write ( output_filename, m, n, g2 );

  printf ( "\n" );
  printf ( "  Wrote denoised image to \"%s\".\n", output_filename );

  free ( g );
  free ( g2 );

  return;
}
