# include <stdio.h>
# include <math.h>
# include <stdlib.h>

int main ( void );
void color ( int red, int green, int blue );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    MANDELBROT_PPM creates an ASCII PPM image of the Mandelbrot set.

  Notice:

    mandelbrot_ppm.C   by Eric R. Weeks   written 9-28-96
    weeks@physics.emory.edu
    http://www.physics.emory.edu/~weeks/

    This program is public domain, but this header must be left intact
    and unchanged.

  Author:

    Eric Weeks

  Local Parameters:

    Local, int HXRES, the horizontal resolution (number of pixels);

    Local, int HYRES, the vertical resolution (number of pixels);

    Local, int ITERMAX, the number of iterations to carry out.

    Local, double MAGNIFY, the magnification factor.
*/
{
  double cx;
  double cy;
  int hx;
  int hxres = 500;
  int hy;
  int hyres = 500;
  int it2;
  int iteration;
  int itermax = 100;
  double magnify = 1.0;
  double x;
  double x_new;
  double y;
  double y_new;
/* 
  Write the PPM header.
*/
  printf ( "P6\n" );
  printf ( "# CREATOR: Eric R Weeks / mandel program\n");
  printf ( "%d %d\n255\n", hxres, hyres );

  for ( hy = 1; hy <= hyres; hy++ )
  {
    for ( hx = 1; hx <= hxres; hx++ )
    {
      cx = (((double)hx)/((double)hxres)-0.5)/magnify*3.0-0.7;
      cy = (((double)hy)/((double)hyres)-0.5)/magnify*3.0;
      x = 0.0;
      y = 0.0;

      it2 = itermax + 1;

      for ( iteration = 1; iteration <= itermax; iteration++ )
      {
        x_new = x*x-y*y+cx;
        y_new = 2.0*x*y+cy;

        x = x_new;
        y = y_new;

        if ( 100.0 < x * x + y * y ) 
        { 
          it2 = iteration; 
          break; 
        }

      }

      if ( it2 < itermax )
      {
        color ( 200+(55*it2)/100, (230*(100-it2))/100, (230*(100-it2))/100 );
      }
      else
      {
        color(0,255,255);
      }
    }
  }
  return 0;
}
/******************************************************************************/

void color ( int red, int green, int blue )

/******************************************************************************/
/*
  Purpose:

    COLOR writes the R, G and B colors for one pixel.

  Author:

    Eric Weeks

  Parameters:

    Input, int RED, GREEN, BLUE, the colors for the pixel.
*/
{
  fputc ( ( char ) red, stdout );
  fputc ( ( char ) green, stdout );
  fputc ( ( char ) blue, stdout );

  return;
}
