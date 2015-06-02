# include <stdlib.h>
# include <stdio.h>

# define DIM     2               /* Dimension of points */
typedef int     tPointi[DIM];   /* type integer point */
# define PMAX    10000           /* Max # of pts in polygon */
   
typedef tPointi tPolygoni[PMAX];/* type integer polygon */

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    COMB generates a comb polygon.

  Discussion:

    A comb polygon is generated, with N vertices, listed in 
    counterclockwise order.

    The comb polygon is a good test case for software that triangulates
    a polygon.

  Licensing:

    This code is Copyright 1997 by Joseph O'Rourke.  It may be freely
    redistributed in its entirety provided that this copyright notice is
    not removed.

  Modified:

    30 April 2007

  Author:

    Joseph O'Rourke, Min Xu

  Reference:

    Joseph ORourke,
    Computational Geometry,
    Second Edition,
    Cambridge, 1998,
    ISBN: 0521649765,
    LC: QA448.D38.
*/
{
# define X 0
# define Y 1

  int i;
  int j;
  int n;
  tPolygoni P;

  printf ( "n=" );
  scanf ( "%d", &n );

  for ( i = 0; i < n/2; i++ )
  {
    P[2*i][X] = n-2 - (2*i);
    P[2*i][Y] = 0;
    P[2*i+1][X] = n-2 - (2*i + 1);
    P[2*i+1][Y] = 10;
  }
  
  P[n-1][X] = ( n - 2 ) / 2;
  P[n-1][Y] = -2;

  printf ( "%d\n", n );
  for ( i = 0; i< n; i++)
  {
    printf ( "%d %d\n", P[i][X], P[i][Y] );
  }
  return 0;

# undef X
# undef Y
}
