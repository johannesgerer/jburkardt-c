# include <stdio.h>
# include <stdlib.h>

# include "gnuplot_i.h"

# define SLEEP_LGTH  1

int main ( int argc, char *argv[] );

/******************************************************************************/

int main ( int argc, char *argv[] ) 

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for ANIM.

  Modified:

    24 June 2011
*/
{
  gnuplot_ctrl *h1;
  double phase;
 
  printf ( "\n" );
  printf ( "ANIM:\n" );
  printf ( "  C version\n" );
  printf ( "  Demonstrate how a running C program can create plots\n" );
  printf ( "  during execution by calling gnuplot, using the\n" );
  printf ( "  gnuplot_i interface program.\n" );
/*
  Open a GNUPLOT process.
*/
  h1 = gnuplot_init ( );

  if ( h1 == NULL )
  {
    printf ( "\n" );
    printf ( "ANIM - Fatal error!\n" );
    printf ( "  The gnuplot command is not in your path.\n" );
    exit ( 1 );
  }

  for ( phase = 0.1; phase < 10.0; phase = phase + 0.1 )
  {
    gnuplot_resetplot ( h1 );
    gnuplot_cmd ( h1, "plot sin(x+%g)", phase );
  }

  for ( phase = 10.0; 0.0 <= phase; phase = phase - 0.1 ) 
  {
    gnuplot_resetplot ( h1 );
    gnuplot_cmd ( h1, "plot sin(x+%g)", phase );
  }
/*
  Close the GNUPLOT process.
*/    
  gnuplot_close ( h1 );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "ANIM:\n" );
  printf ( "  Normal end of execution.\n" );

  return 0;
}

