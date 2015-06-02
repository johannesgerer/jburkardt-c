# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "asa058.h"

int main ( void );
void test01 ( void );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for ASA058_PRB.

  Discussion:

    ASA058_PRB tests the ASA058 library.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    03 February 2008

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "ASA058_PRB:\n" );
  printf ( "  C version\n" );
  printf ( "  Test the ASA058 library.\n" );

  test01 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "ASA058_PRB:\n" );
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

    TEST01 tries out the ASA058 routine.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    04 February 2008

  Author:

    John Burkardt
*/
{
# define K 5
# define LINE_MAX 80
# define M 2
# define N 100

  int b[N];
  double d[K*M];
  double dev[K];
  double dev_sum;
  int e[K];
  int e_sum;
  double f[N];
  int i;
  FILE *input;
  char input_filename[] = "points_100.txt";
  int j;
  int k2;
  char line[LINE_MAX];
  int nz;
  float temp;
  double x[N*M];

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  Test the CLUSTR algorithm.\n" );
  printf ( "  Applied Statistics Algorithm 58\n" );
/*
  Read the data.
*/
  printf ( "\n" );
  printf ( "  Reading the data.\n" );

  input = fopen ( input_filename, "r" );

  if ( !input )
  {
    printf ( "\n" );
    printf ( "TEST01 - Fatal error!\n" );
    printf ( "  Could not open the input file: \"%s\".\n", input_filename );
    return;
  }
/*
  It is unbelievable idiocy that FSCANF must read a FLOAT...
  If TEMP were double, I would just get zeros...as I did for
  half an hour trying to figure this error out.
*/
  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < 2; j++ )
    {
      fscanf ( input, "%f", &temp );
      x[i+j*N] = temp;
      printf ( "  %f", temp );
    }
    printf ( "\n" );
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
    for ( j = 1; j <= M; j++ )
    {
      printf ( "  %14f", x[i-1+(j-1)*N] );
    }
    printf ( "\n" );
  }
/*
  Initialize the cluster centers arbitrarily.
*/
  for ( i = 1; i <= K; i++ )
  {
    for ( j = 1; j <= M; j++ )
    {
      d[i-1+(j-1)*K] = x[i-1+(j-1)*N];
    }
  }
/*
  Compute the clusters.
*/
  nz = 1;
  k2 = K;

  clustr ( x, d, dev, b, f, e, N, M, K, nz, k2 );

  printf ( "\n" );
  printf ( "  Cluster  Population  Energy\n" );
  printf ( "\n" );

  for ( i = 1; i <= K; i++ )
  {
    printf ( "  %8d  %8d  %14f\n", i, e[i-1], dev[i-1] );
  }

  e_sum = 0;
  dev_sum = 0.0;

  for ( i = 1; i <= K; i++ )
  {
    e_sum = e_sum + e[i-1];
    dev_sum = dev_sum + dev[i-1];
  }

  printf ( "\n" );
  printf ( "     Total  %8d  %14f\n", e_sum, dev_sum );

  return;
# undef K
# undef LINE_MAX
# undef M
# undef N
}
