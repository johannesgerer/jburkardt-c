# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "asa136.h"

int main ( void );
void test01 ( void );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for ASA136_PRB.

  Discussion:

    ASA136_PRB calls the ASA136 routines.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    09 November 2010

  Author:

    John Burkardt
*/
{
  timestamp ( );

  printf ( "\n" );
  printf ( "ASA136_PRB:\n" );
  printf ( "  C version\n" );
  printf ( "  Test the ASA136 library.\n" );

  test01 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "ASA136_PRB:\n" );
  printf ( "  Normal end of execution.\n" );

  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void test01 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST01 tries out the ASA136 routine.
  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    09 November 2010

  Author:

    John Burkardt
*/
{
  double *a;
  double *c;
  int i;
  int *ic1;
  int ifault;
  FILE *input;
  int iter;
  int j;
  int k = 5;
  int m = 100;
  int n = 2;
  int *nc;
  int nc_sum;
  float value;
  double *wss;
  double wss_sum;

  a = ( double * ) malloc ( m * n * sizeof ( double ) );
  c = ( double * ) malloc ( k * n * sizeof ( double ) );
  ic1 = ( int * ) malloc ( m * sizeof ( int ) );
  nc = ( int * ) malloc ( k * sizeof ( int ) );
  wss = ( double * ) malloc ( k * sizeof ( double ) );

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  Test the KMNS algorithm,\n" );
  printf ( "  Applied Statistics Algorithm #136.\n" );
/*
  Read the data.
*/
  input = fopen ( "points_100.txt", "rt" );

  for ( i = 1; i <= m; i++ )
  {
    for ( j = 1; j <= n; j++ )
    {
      fscanf ( input, "%f", &value );
      a[i-1+(j-1)*m] = value;
    }
  }

  fclose ( input );
/*
  Print a few data values.
*/
  printf ( "\n" );
  printf ( "  First 5 data values:\n" );
  printf ( "\n" );

  for ( i = 1; i <= 5; i++ )
  {
    printf ( "  %8d", i );
    for ( j = 1; j <= n; j++ )
    {
      printf ( "  %14f", a[i-1+(j-1)*m] );
    }
    printf ( "\n" );
  }
/*
  Initialize the cluster centers.
  Here, we arbitrarily make the first K data points cluster centers.
*/
  for ( i = 1; i <= k; i++ )
  {
    for ( j = 1; j <= n; j++ )
    {
      c[i-1+(j-1)*k] = a[i-1+(j-1)*m];
    }
  }

  iter = 50;
/*
  Compute the clusters.
*/
  kmns ( a, m, n, c, k, ic1, nc, iter, wss, &ifault );

  if ( ifault != 0 )
  {
    printf ( "\n" );
    printf ( "TEST01 - Fatal error!\n" );
    printf ( "  KMNS returned IFAULT = %d\n", ifault );
    return;
  }

  printf ( "\n" );
  printf ( "  Cluster  Population  Energy\n" );
  printf ( "\n" );

  nc_sum = 0;
  wss_sum = 0.0;

  for ( i = 1; i <= k; i++ )
  {
    printf ( "  %8d  %8d  %14f\n", i, nc[i-1], wss[i-1] );
    nc_sum = nc_sum + nc[i-1];
    wss_sum = wss_sum + wss[i-1];
  }

  printf ( "\n" );
  printf ( "     Total  %8d  %14f\n", nc_sum, wss_sum );

  free ( a );
  free ( c );
  free ( ic1 );
  free ( nc );
  free ( wss );

  return;
}
