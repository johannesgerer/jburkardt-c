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

    MAIN demonstrates the use of the CURVE routine.

  Modified:
 
    09 April 2011

  Reference:

    Helmut Michels,
    The Data Plotting Software DISLIN - version 10.4,
    Shaker Media GmbH, January 2010,
    ISBN13: 978-3-86858-517-9.
*/
{
  static char cbuf[41];
  static char *ctit = "Pie Charts (PIEGRF)";
  int i;
  int nya = 2800;
  static float xray[5] = { 1.0, 2.5, 2.0, 2.7, 1.8 };

  printf ( "\n" );
  printf ( "DISLIN_EX06:\n" );
  printf ( "  C version:\n" );
  printf ( "  Demonstrate the use of the PIEGRF routine, for\n" );
  printf ( "  plotting piechart data.\n" );
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
  setfil ( "dislin_ex06.png" );
/*
  Choose the page size and orientation.
*/
  setpag ( "usap" );
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
  axslen ( 1600, 1000 );
  titlin ( ctit, 2 );
  chnpie ( "both" );

  legini ( cbuf, 5, 8 );
  leglin ( cbuf, "FIRST", 1 );
  leglin ( cbuf, "SECOND", 2 );
  leglin ( cbuf, "THIRD", 3 );
  leglin ( cbuf, "FOURTH", 4 );
  leglin ( cbuf, "FIFTH", 5 );

  patcyc ( 1, 7L );
  patcyc ( 2, 4L );
  patcyc ( 3, 13L );
  patcyc ( 4, 3L );
  patcyc ( 5, 5L );

  for ( i = 0; i < 2; i++ )
  {
    axspos ( 250, nya-i*1200 );
    if ( i == 1 )
    {
      labels ( "data", "pie" );
      labpos ( "external", "pie" );
    }

    piegrf ( cbuf, 1, xray, 5 );

    if ( i == 1 )
    {
      height ( 50 );
      title ( );
    }
    endgrf ( );
  }
/*
  Close DISLIN.
*/
  disfin ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "DISLIN_EX06:\n" );
  printf ( "  Normal end of execution.\n" );

  return 0;
}
