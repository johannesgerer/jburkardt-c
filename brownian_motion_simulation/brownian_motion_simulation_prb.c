# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <string.h>

# include "brownian_motion_simulation.h"

int main ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    BROWNIAN_MOTION_SIMULATION_PRB tests the BROWNIAN_MOTION_SIMULATION library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 October 2012

  Author:

    John Burkardt
*/
{
  double d;
  double *dsq;
  char header[80];
  int k;
  int m;
  int n;
  int seed;
  double t;
  double *x;

  timestamp ( );
  printf ( "\n" );
  printf ( "BROWNIAN_MOTION_SIMULATION_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the BROWNIAN_MOTION_SIMULATION library.\n" );
/*
  Compute the path of a particle undergoing Brownian motion.
*/
  for ( m = 1; m <= 2; m++ )
  {
    n = 1001;
    d = 10.0;
    t = 1.0;
    seed = 123456789;
    x = brownian_motion_simulation ( m, n, d, t, &seed );
    if ( m == 1 )
    {
      strcpy ( header, "motion_1d" );
    }
    else if ( m == 2 )
    {
      strcpy ( header, "motion_2d" );
    }
    brownian_motion_display ( m, n, x, header );
    free ( x );
  }
/*
  Estimate the average displacement of the particle from the origin
  as a function of time.
*/
  for ( m = 1; m <= 3; m++ )
  {
    k = 40;
    n = 1001;
    d = 10.0;
    t = 1.0;
    seed = 123456789;

    dsq = brownian_displacement_simulation ( k, n, m, d, t, &seed );
    if ( m == 1 )
    {
      strcpy ( header, "displacement_1d" );
    }
    else if ( m == 2 )
    {
      strcpy ( header, "displacement_2d" );
    }
    else if ( m == 3 )
    {
      strcpy ( header, "displacement_3d" );
    }
    brownian_displacement_display ( k, n, d, t, dsq, header );
    free ( dsq );
  }
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "BROWNIAN_MOTION_SIMULATION_PRB\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
