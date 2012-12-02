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

    MAIN demonstrates the creation of a map plot.

  Modified:
 
    09 April 2011

  Reference:

    Helmut Michels,
    The Data Plotting Software DISLIN - version 10.4,
    Shaker Media GmbH, January 2010,
    ISBN13: 978-3-86858-517-9.
*/
{
  printf ( "\n" );
  printf ( "DISLIN_EX13:\n" );
  printf ( "  C version\n" );
  printf ( "  Demonstrate the creation of a map plot.\n" );
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
  setfil ( "dislin_ex13.png" );
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
  Plot a border around the page.
*/
  pagera ( );
/*
  Use the COMPLEX font.
*/
  complx ( );

  frame ( 3 );
  axspos ( 400, 1850 );
  axslen ( 2400, 1400 );

  name ( "Longitude", "x" );
  name ( "Latitude", "y" );
  titlin ( "World Coastlines and Lakes", 3 );

  labels ( "map", "xy" );
  grafmp ( -180.0, 180.0, -180.0, 90.0, -90.0, 90.0, -90.0, 30.0 );

  gridmp ( 1, 1 );
  color ( "green" );
  world ( );
  color ( "fore" );

  height ( 50 );
  title ( );
/*
  Close DISLIN.
*/
  disfin ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "DISLIN_EX13:\n" );
  printf ( "  Normal end of execution.\n" );

  return 0;
}
