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
# define N 50

  static char *ctit1 = "Surface Plot (SURMAT)";
  static char *ctit2 = "F(X,Y) = 2*SIN(X)*SIN(Y)";
  double fpi = 3.1415927 / 180.0;
  int i;
  int j;
  double step;
  double x;
  double y;
  float zmat[N][N];

  printf ( "\n" );
  printf ( "DISLIN_EX10:\n" );
  printf ( "  C version:\n" );
  printf ( "  Demonstrate the use of color 3D graphics.\n" );

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
  setfil ( "dislin_ex10.png" );
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
  axspos ( 200, 2600 );
  axslen ( 1800, 1800 );

  name ( "X-axis", "x" );
  name ( "Y-axis", "y" );
  name ( "Z-axis", "z" );

  titlin ( ctit1, 2 );
  titlin ( ctit2, 4 );

  view3d ( -5.0, -5.0, 4.0, "abs" );
  graf3d ( 0.0, 360.0, 0.0, 90.0, 0.0, 360.0, 0.0, 90.0, -3.0, 3.0, -3.0, 1.0 );
  height ( 50 );
  title ( );

  color ( "green" );
  surmat ( (float *) zmat, N, N, 1, 1 );
/*
  Close DISLIN.
*/
  disfin ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "DISLIN_EX10:\n" );
  printf ( "  Normal end of execution.\n" );

  return 0;
# undef N
}
