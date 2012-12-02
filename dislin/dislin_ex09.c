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
# define N 100

  double fpi = 3.1415927 / 180.0;
  int i;
  int j;
  double step;
  double x;
  double y;
  float zmat[N][N];

  printf ( "\n" );
  printf ( "DISLIN_EX09:\n" );
  printf ( "  C version:\n" );
  printf ( "  Demonstrate the creation of a 3D color plot.\n" );

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
  setfil ( "dislin_ex09.png" );
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
  Use the HARDWARE font.
*/
  hwfont ( );

  titlin ( "3-D Color Plot of the Function", 2 );
  titlin ( "F(X,Y) = 2 * SIN(X) * SIN(Y)", 4 );

  name ( "X-axis", "x" );
  name ( "Y-axis", "y" );
  name ( "Z-axis", "z" );

  intax ( );
  autres ( N, N );
  axspos ( 300, 1850 );
  ax3len ( 2200, 1400, 1400 );

  graf3 ( 0.0, 360.0, 0.0, 90.0, 0.0, 360.0, 0.0, 90.0, -2.0, 2.0, -2.0, 1.0 );
  crvmat ( (float *) zmat, N, N, 1, 1 );
  
  height ( 50 );
  title ( );
  mpaepl ( 3 );
/*
  Close DISLIN.
*/
  disfin ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "DISLIN_EX09:\n" );
  printf ( "  Normal end of execution.\n" );

  return 0;
}
