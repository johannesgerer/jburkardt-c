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

    QUICKPLOT_PIE demonstrates the DISLIN quickplot command QPLPIE.

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
  int i;
  int n = 5;
  float xray[5] = { 10.0, 20.0, 15.0, 5.0, 50.0 };

  printf ( "\n" );
  printf ( "QUICKPLOT_PIE:\n" );
  printf ( "  C version:\n" );
  printf ( "  Demonstrate the DISLIN 'quickplot' command QPLPIE\n" );
  printf ( "  to plot a curve.\n" );
  printf ( "\n" );
  printf ( "  Here, we plot 10 percent luck, 20 percent skill,\n" );
  printf ( "  15 percent concentrated power of will, 5 percent pleasure,\n" );
  printf ( "  50 percent pain.\n" );
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
  setfil ( "quickplot_pie.png" );
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
  titlin ( "Quick plot by QPLPIE", 2 );
/*
  Draw the curve.
*/
  qplpie ( xray, n );
/*
  Close DISLIN.
*/
  disfin ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "QUICKPLOT_PIE:\n" );
  printf ( "  Normal end of execution.\n" );

  return 0;
}
