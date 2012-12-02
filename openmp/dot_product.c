# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <omp.h>

int main ( int argc, char *argv[] );
double test01 ( int n, double x[], double y[] );
double test02 ( int n, double x[], double y[] );

/******************************************************************************/

int main ( int argc, char *argv[] )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for DOT_PRODUCT.

  Discussion:

    This program illustrates how a vector dot product could be set up
    in a C program using OpenMP.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    18 April 2009

  Author:

    John Burkardt
*/
{
  double factor;
  int i;
  int n;
  double wtime;
  double *x;
  double xdoty;
  double *y;

  printf ( "\n" );
  printf ( "DOT_PRODUCT\n" );
  printf ( "  C/OpenMP version\n" );
  printf ( "\n" );
  printf ( "  A program which computes a vector dot product.\n" );

  printf ( "\n" );
  printf ( "  Number of processors available = %d\n", omp_get_num_procs ( ) );
  printf ( "  Number of threads =              %d\n", omp_get_max_threads ( ) );
/*
  Set up the vector data.
  N may be increased to get better timing data.

  The value FACTOR is chosen so that the correct value of the dot product 
  of X and Y is N.
*/
  n = 100;

  while ( n < 1000000 )
  {
    n = n * 10;

    x = ( double * ) malloc ( n * sizeof ( double ) );
    y = ( double * ) malloc ( n * sizeof ( double ) );

    factor = ( double ) ( n );
    factor = 1.0 / sqrt ( 2.0 * factor * factor + 3 * factor + 1.0 );

    for ( i = 0; i < n; i++ )
    {
      x[i] = ( i + 1 ) * factor;
    }

    for ( i = 0; i < n; i++ )
    {
      y[i] = ( i + 1 ) * 6 * factor;
    }

    printf ( "\n" );
/*
  Test #1
*/
    wtime = omp_get_wtime ( );

    xdoty = test01 ( n, x, y );

    wtime = omp_get_wtime ( ) - wtime;

    printf ( "  Sequential  %8d  %14.6e  %15.10f\n", n, xdoty, wtime );
/*
  Test #2
*/
    wtime = omp_get_wtime ( );

    xdoty = test02 ( n, x, y );

    wtime = omp_get_wtime ( ) - wtime;

    printf ( "  Parallel    %8d  %14.6e  %15.10f\n", n, xdoty, wtime );
  
    free ( x );
    free ( y );
  }
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "DOT_PRODUCT\n" );
  printf ( "  Normal end of execution.\n" );

  return 0;
}
/******************************************************************************/

double test01 ( int n, double x[], double y[] )

/******************************************************************************/
/*
  Purpose:

    TEST01 computes the dot product with no parallel processing directives.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    05 April 2008

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of the vectors.

    Input, double X[N], Y[N], the vectors.

    Output, double TEST01, the dot product of X and Y.
*/
{
  int i;
  double xdoty;

  xdoty = 0.0;

  for ( i = 0; i < n; i++ )
  {
    xdoty = xdoty + x[i] * y[i];
  }

  return xdoty;
}
/******************************************************************************/

double test02 ( int n, double x[], double y[] )

/******************************************************************************/
/*
  Purpose:

    TEST02 computes the dot product with parallel processing directives.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    05 April 2008

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of the vectors.

    Input, double X[N], Y[N], the vectors.

    Output, double TEST02, the dot product of X and Y.
*/
{
  int i;
  double xdoty;

  xdoty = 0.0;

# pragma omp parallel \
  shared ( n, x, y ) \
  private ( i )

# pragma omp for reduction ( + : xdoty )

  for ( i = 0; i < n; i++ )
  {
    xdoty = xdoty + x[i] * y[i];
  }

  return xdoty;
}

