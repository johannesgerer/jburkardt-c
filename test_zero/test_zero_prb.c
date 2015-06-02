# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <complex.h>
# include <string.h>

# include "test_zero.h"

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for TEST_ZERO_PRB.

  Discussion:

    TEST_ZERO_PRB tests the TEST_ZERO library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 May 2011

  Author:

    John Burkardt
*/
{
  double fatol = 1.0E-06;
  double fx;
  double fxa;
  double fxb;
  double fxc;
  int i;
  int max_step = 25;
  int prob;
  int prob_num;
  double *range;
  int root_num;
  int start_num;
  char *title;
  double x;
  double xa;
  double xatol = 1.0E-06;
  double xb;
  double xc;
  double xmax;
  double xmin;
  double xrtol = 1.0E-06;

  timestamp ( );
  printf ( "\n" );
  printf ( "TEST_ZERO_PRB\n" );
  printf ( "  C++ version\n" );
  printf ( "  Test the TEST_ZERO library.\n" );
  printf ( "\n" );
  printf ( "  Function value tolerance = %g\n", fatol );
  printf ( "  Root absolute tolerance =  %g\n", xatol );
  printf ( "  Root relative tolerance =  %g\n", xrtol );
  printf ( "  Maximum number of steps =  %d\n", max_step );
//
//  Find out how many problems there are
//
  prob_num = p00_prob_num ( );
  printf ( "\n" );
  printf ( "  Number of problems available is %d\n", prob_num );

  for ( prob = 1; prob <= prob_num; prob++ )
  {
//
//  Print the problem title.
//
    title = p00_title ( prob );
    printf ( "\n" );
    printf ( "\n" );
    printf ( "  Problem number %d\n", prob );
    printf ( "  \"%s\"\n", title );

    if ( prob == 16 )
    {
//    p16_p_print ( );
    }
//
//  Get the problem interval.
//
    range = p00_range ( prob );
    xmin = range[0];
    xmax = range[1];
    printf ( "  We seek roots between %g and %g\n", range[0], range[1] );
//
//  Get the number of roots.
//
    root_num = p00_root_num ( prob );

    printf ( "\n" );
    printf ( "  Number of known roots = %d\n", root_num );
//
//  Get the roots.
//
    if ( 0 < root_num )
    {
      printf ( "\n" );
      printf ( "     I          X          F(X)\n" );
      printf ( "\n" );
      for ( i = 1; i <= root_num; i++ )
      {
        x = p00_root ( prob, i );
        fx = p00_fx ( prob, x );
        printf ( "  %4d  %10g  %10g\n", i, x, fx );
      }
    }
//
//  Get the number of starting points.
//
    start_num = p00_start_num ( prob );

    printf ( "\n" );
    printf ( "  Number of starting points = %d\n", start_num );
//
//  Get the starting points.
//
    printf ( "\n" );
    printf ( "     I          X          F(X)\n" );
    printf ( "\n" );
    for ( i = 1; i <= start_num; i++ )
    {
      x = p00_start ( prob, i );
      fx = p00_fx ( prob, x );
      printf ( "  %4d  %10g  %10g\n", i, x, fx );
    }
//
//  Bisection.
//
    xa = p00_start ( prob, 1 );
    fxa = p00_fx ( prob, xa );

    for ( i = 2; i <= start_num; i++ )
    {
      xb = p00_start ( prob, i );
      fxb = p00_fx ( prob, xb );

      if ( r8_sign ( fxa ) != r8_sign ( fxb ) )
      {
         bisection ( fatol, max_step, prob, xatol, &xa, &xb, &fxa, &fxb );
         break;
      }
    }
//
//  Brent's method.
//
    xa = p00_start ( prob, 1 );
    fxa = p00_fx ( prob, xa );

    for ( i = 2; i <= start_num; i++ )
    {
      xb = p00_start ( prob, i );
      fxb = p00_fx ( prob, xb );

      if ( r8_sign ( fxa ) != r8_sign ( fxb ) )
      {
         brent ( fatol, max_step, prob, xatol, xrtol, &xa, &xb, &fxa, &fxb );
         break;
      }
    }
//
//  Muller's method.
//
    if ( 3 <= p00_start_num ( prob ) )
    {
      xa = p00_start ( prob, 1 );
      fxa = p00_fx ( prob, xa );
      xb = p00_start ( prob, 2 );
      fxb = p00_fx ( prob, xb );
      xc = p00_start ( prob, 3 );
      fxc = p00_fx ( prob, xc );

      muller ( fatol, max_step, prob, xatol, xrtol, &xa, &xb, &xc, &fxa, &fxb, &fxc );
    }
//
//  Newton's method.
//
    for ( i = 1; i <= start_num; i++ )
    {
      xa = p00_start ( prob, i );
      fxa = p00_fx ( prob, xa );
      newton ( fatol, max_step, prob, xatol, xmin, xmax, &xa, &fxa );
    }
//
//  Regula Falsi
//
    xa = p00_start ( prob, 1 );
    fxa = p00_fx ( prob, xa );

    for ( i = 2; i <= start_num; i++ )
    {
      xb = p00_start ( prob, i );
      fxb = p00_fx ( prob, xb );

      if ( r8_sign ( fxa ) != r8_sign ( fxb ) )
      {
         regula_falsi ( fatol, max_step, prob, xatol, &xa, &xb, &fxa, &fxb );
         break;
      }
    }
//
//  Secant.
//
    for ( i = 1; i < start_num; i++ )
    {
      xa = p00_start ( prob, i );
      fxa = p00_fx ( prob, xa );

      xb = p00_start ( prob, i + 1 );
      fxb = p00_fx ( prob, xb );

      secant ( fatol, max_step, prob, xatol, xmin, xmax, &xa, &xb, &fxa, &fxb );
    }
    free ( range );
    free ( title );
  }
//
//  Terminate.
//
  printf ( "\n" );
  printf ( "TEST_ZERO_PRB\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
