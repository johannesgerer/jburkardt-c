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

    QUICKPLOT_CURVE demonstrates the DISLIN quickplot command QPLOT.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:
 
    22 April 2011

  Reference:

    Helmut Michels,
    The Data Plotting Software DISLIN - version 10.4,
    Shaker Media GmbH, January 2010,
    ISBN13: 978-3-86858-517-9.
*/
{
  int i;
  int n = 100;
  float pi = 3.1415926;
  float *xray;
  float *yray;

  printf ( "\n" );
  printf ( "QUICKPLOT_CURVE:\n" );
  printf ( "  C version:\n" );
  printf ( "  Demonstrate the DISLIN 'quickplot' command QPLOT\n" );
  printf ( "  to plot a curve.\n" );
/*
  Set up the X and Y data for the plot.
*/
  xray = ( float * ) malloc ( n * sizeof ( float ) );
  yray = ( float * ) malloc ( n * sizeof ( float ) );

  for ( i = 0; i < n; i++ )
  {
    xray[i] = ( float ) ( i ) * 360.0 / ( float ) ( n - 1 );
  }
  for ( i = 0; i < n; i++ )
  {
    yray[i] = sin ( pi * xray[i] / 180.0 );
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
  setfil ( "quickplot_curve.png" );
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
  name ( "<-- Angle in Degrees -->", "X" );
  name ( "<-- Sine (angle) -->", "Y" );
  titlin ( "Quick plot by QPLOT", 2 );
/*
  Draw the curve.
*/
  qplot ( xray, yray, n );
/*
  Close DISLIN.
*/
  disfin ( );
/*
  Free memory.
*/
  free ( xray );
  free ( yray );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "QUICKPLOT_CURVE:\n" );
  printf ( "  Normal end of execution.\n" );

  return 0;
}
