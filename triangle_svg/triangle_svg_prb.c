# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <string.h>

# include "triangle_svg.h"

int main ( );
void test01 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for TRIANGLE_SVG_PRB.

  Discussion:

    TRIANGLE_SVG_PRB tests the TRIANGLE_SVG library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    21 April 2014

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "TRIANGLE_SVG_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the TRIANGLE_SVG library.\n" );

  test01 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "TRIANGLE_SVG_PRB\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void test01 ( )

/******************************************************************************/
/*
  Purpose:

    TEST01 calls TRIANGLE_SVG to plot a triangle and some points.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    21 April 2014

  Author:

    John Burkardt
*/
{ 
  double p[2*13] = {
    0.333333333333333, 0.333333333333333,
    0.479308067841923, 0.260345966079038,
    0.260345966079038, 0.479308067841923,
    0.260345966079038, 0.260345966079038,
    0.869739794195568, 0.065130102902216,
    0.065130102902216, 0.869739794195568,
    0.065130102902216, 0.065130102902216,
    0.638444188569809, 0.312865496004875,
    0.638444188569809, 0.048690315425316,
    0.312865496004875, 0.638444188569809,
    0.312865496004875, 0.048690315425316,
    0.048690315425316, 0.638444188569809,
    0.048690315425316, 0.312865496004875 };
  int p_num = 13;
  char plot_filename[] = "test01.svg";
  double t[2*3] = {
    0.0, 0.0, 
    1.0, 0.0, 
    0.0, 1.0 };

  triangle_svg ( plot_filename, t, p_num, p );

  return;
}
