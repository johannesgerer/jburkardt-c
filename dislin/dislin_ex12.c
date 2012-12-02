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

    MAIN demonstrates the creation of a shaded contour plot.

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

  int i;
  int j;
  float step;
  float x;
  float xray[N];
  float y;
  float yray[N];
  float zlev[12];
  float zmat[N][N];

  printf ( "\n" );
  printf ( "DISLIN_EX12:\n" );
  printf ( "  C version:\n" );
  printf ( "  Demonstrate a shaded contour plot.\n" );

  step = 1.6 / ( float ) ( N - 1 );
  for (i = 0; i < N; i++)
  {
    x = 0.0 + ( float ) i * step;
    xray[i] = x;
    for ( j = 0; j < N; j++)
    {
      y = 0.0 + ( float ) j * step;
      yray[j] = y;
      zmat[i][j] = ( x * x - 1.0 ) * ( x * x - 1.0) 
                 + ( y * y - 1.0 ) * ( y * y - 1.0 );
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
  setfil ( "dislin_ex12.png" );
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

  mixalf ( );
  titlin ( "Shaded Contour Plot", 1 );
  titlin ( "F(X,Y) = (X[2$ - 1)[2$ + (Y[2$ - 1)[2$", 3 );
  name ( "X-axis", "x" );
  name ( "Y-axis", "y" );

  shdmod ( "poly", "contur" );
  axspos ( 450, 2670 );
  graf ( 0.0, 1.6, 0.0, 0.2, 0.0, 1.6, 0.0, 0.2 );

  for ( i = 1; i <= 12; i++ )
  {
    zlev[12-i] = 0.1 + ( float ) ( i - 1 ) * 0.1;
  }

  conshd ( xray, N, yray, N, (float *) zmat, zlev, 12 );

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
  printf ( "DISLIN_EX12:\n" );
  printf ( "  Normal end of execution.\n" );

  return 0;

# undef N
}
