# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "dislin.h"

int main ( int argc, char *argv[] );

/******************************************************************************/

int main ( int argc, char *argv[] )

/******************************************************************************/
/*
  Purpose:

    QUICKPLOT_COLOR demonstrates the DISLIN quickplot command QPLCLR.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:
 
    13 May 2012

  Reference:

    Helmut Michels,
    The Data Plotting Software DISLIN - version 10.4,
    Shaker Media GmbH, January 2010,
    ISBN13: 978-3-86858-517-9.
*/
{
# define N 100

  float fpi = 3.1415927 / 180.0;
  int i;
  int j;
  int n = N;
  float step;
  float x;
  float y;
  float zmat[N][N];

  printf ( "\n" );
  printf ( "QUICKPLOT_CURVE:\n" );
  printf ( "  C version:\n" );
  printf ( "  Demonstrate the DISLIN 'quickplot' command QPLCLR\n" );
  printf ( "  to make a color plot of a matrix of data.\n" );
/*
  Set up the X and Y data for the plot.
*/
  step = 360.0 / ( float ) ( N - 1 );
  for ( i = 0; i < N; i++ )
  {
    x = ( float ) i * step;
    for ( j = 0; j < N; j++ )
    {
      y = ( float ) j * step;
      zmat[i][j] = 2.0 * sin ( x * fpi ) * sin ( y * fpi );
    }
  }
/*
  Specify the format of the output file.
*/
  metafl ( "png" );
/*
  Specify that if a file already exists of the given name,
  the new data should overwrite the old.
*/
  filmod ( "delete" );
/*
  Specify the name of the output graphics file.
*/
  setfil ( "quickplot_color.png" );
/*
  Choose the page size and orientation.
*/
  setpag ( "usal" );
/*
  For PNG output, reverse the default black background to white.
*/
  scrmod ( "reverse" );
/*
  Open DISLIN.
*/
  disini ( );
/*
  Label the axes and the plot.
*/
  name ( "X-axis", "x" );
  name ( "Y-axis", "y" );
  titlin ( "Quick plot by QPLCLR", 2 );
/*
  Draw the curve.
*/
  qplclr ( ( float * ) zmat, n, n );
/*
  Close DISLIN.
*/
  disfin ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "QUICKPLOT_COLOR:\n" );
  printf ( "  Normal end of execution.\n" );

  return 0;
}
