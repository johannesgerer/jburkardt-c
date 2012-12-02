# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "simplex_coordinates.h"

int main ( );
void test01 ( int n );
void test02 ( int n );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    SIMPLEX_COORDINATES_PRB tests SIMPLEX_COORDINATES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 September 2010

  Author:

    John Burkardt
*/
{
  int n;

  timestamp ( );
  printf ( "\n" );
  printf ( "SIMPLEX_COORDINATES_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the SIMPLEX_COORDINATES library.\n" );

  n = 3;
  test01 ( n );
  test02 ( n );

  n = 4;
  test01 ( n );
  test02 ( n );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "SIMPLEX_COORDINATES_PRB\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );

  timestamp ( );

  return 0;
}
/******************************************************************************/

void test01 ( int n )

/******************************************************************************/
/*
  Purpose:

    TEST01 calls SIMPLEX_COORDINATES1.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 September 2010

  Author:

    John Burkardt

  Parameters:

    Input, int N, the spatial dimension.
*/
{
  int i;
  int j;
  int k;
  double side;
  double volume;
  double volume2;
  double *x;
  double *xtx;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  Call SIMPLEX_COORDINATES1\n" );

  x = simplex_coordinates1 ( n );

  r8mat_transpose_print ( n, n + 1, x, "  Simplex vertex coordinates:" );

  side = 0.0;
  for ( i = 0; i < n; i++ )
  {
    side = side + pow ( x[i+0*n] - x[i+1*n], 2 );
  }
  side = sqrt ( side );

  volume = simplex_volume ( n, x );

  volume2 = sqrt ( ( double ) ( n + 1 ) ) / r8_factorial ( n ) 
    / sqrt ( pow ( 2.0, n ) ) * pow ( side, n );

  printf ( "\n" );
  printf ( "  Side length =     %f\n", side );
  printf ( "  Volume =          %f\n", volume );
  printf ( "  Expected volume = %f\n", volume2 );

  xtx = ( double * ) malloc ( ( n + 1 ) * ( n + 1 ) * sizeof ( double ) );

  for ( j = 0; j < n + 1; j++ )
  {
    for ( i = 0; i < n + 1; i++ )
    {
      xtx[i+j*(n+1)] = 0.0;
      for ( k = 0; k < n; k++ )
      {
        xtx[i+j*(n+1)] = xtx[i+j*(n+1)] + x[k+i*n] * x[k+j*n];
      }
    }
  }

  r8mat_transpose_print ( n + 1, n + 1, xtx, "  Dot product matrix:" );

  free ( x );
  free ( xtx );

  return;
}
/******************************************************************************/

void test02 ( int n )

/******************************************************************************/
/*
  Purpose:

    TEST02 calls SIMPLEX_COORDINATES2.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 September 2010

  Author:

    John Burkardt

  Parameters:

    Input, int N, the spatial dimension.
*/
{
  int i;
  int j;
  int k;
  double side;
  double volume;
  double volume2;
  double *x;
  double *xtx;

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  Call SIMPLEX_COORDINATES2\n" );

  x = simplex_coordinates2 ( n );

  r8mat_transpose_print ( n, n + 1, x, "  Simplex vertex coordinates:" );

  side = 0.0;
  for ( i = 0; i < n; i++ )
  {
    side = side + pow ( x[i+0*n] - x[i+1*n], 2 );
  }
  side = sqrt ( side );

  volume = simplex_volume ( n, x );

  volume2 = sqrt ( ( double ) ( n + 1 ) ) / r8_factorial ( n ) 
    / sqrt ( pow ( 2.0, n ) ) * pow ( side, n );

  printf ( "\n" );
  printf ( "  Side length =     %f\n", side );
  printf ( "  Volume =          %f\n", volume );
  printf ( "  Expected volume = %f\n", volume2 );

  xtx = ( double * ) malloc ( ( n + 1 ) * ( n + 1 ) * sizeof ( double ) );

  for ( j = 0; j < n + 1; j++ )
  {
    for ( i = 0; i < n + 1; i++ )
    {
      xtx[i+j*(n+1)] = 0.0;
      for ( k = 0; k < n; k++ )
      {
        xtx[i+j*(n+1)] = xtx[i+j*(n+1)] + x[k+i*n] * x[k+j*n];
      }
    }
  }

  r8mat_transpose_print ( n + 1, n + 1, xtx, "  Dot product matrix:" );

  free ( x );
  free ( xtx );

  return;
}
