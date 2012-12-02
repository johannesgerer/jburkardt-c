# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <omp.h>

int main ( int argc, char *argv[] );
double f ( double x, double y );
double cpu_time ( void );

/******************************************************************************/

int main ( int argc, char *argv[] )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for QUAD2D_OPENMP.
*/
{
  double a;
  double b;
  double error;
  double exact;
  int i;
  int j;
  int n;
  int nx;
  int ny;
  double pi;
  double total;
  double wtime;
  double x;
  double y;

  a = 0.0;
  b = 1.0;
  nx = 32768;
  ny = 32768;
  n = nx * ny;
  pi = 3.141592653589793;
  exact = pi * pi / 6.0;

  printf ( "\n" );
  printf ( "QUAD2D_OPENMP:\n" );
  printf ( "  C/OpenMP version\n" );
  printf ( "  Estimate the integral of f(x,y) over [0,1]x[0,1].\n" );
  printf ( "  f(x,y) = 1 / ( 1 - x * y ).\n" );
  printf ( "\n" );
  printf ( "  A        = %f\n", a );
  printf ( "  B        = %f\n", b );
  printf ( "  NX       = %d\n", nx );
  printf ( "  NY       = %d\n", ny );
  printf ( "  N        = %d\n", n );
  printf ( "  Exact    = %24.16f\n", exact );

  wtime = omp_get_wtime ( );

  total = 0.0;
# pragma omp parallel shared ( a, b, n ) private ( i, j, x, y )

# pragma omp for reduction ( + : total )

  for ( i = 1; i <= nx; i++ )
  {
    x = ( ( 2 * nx - 2 * i + 1 ) * a + ( 2 * i - 1 ) * b ) / ( 2 * nx );
    for ( j = 1; j <= ny; j++ )
    {
      y = ( ( 2 * ny - 2 * j + 1 ) * a + ( 2 * j - 1 ) * b ) / ( 2 * ny );
      total = total + f ( x, y );
    }
  }

  wtime = omp_get_wtime ( ) - wtime;

  total = ( b - a ) * ( b - a ) * total / ( double ) ( nx ) / ( double ) ( ny );
  error = fabs ( total - exact );
 
  printf ( "\n" );
  printf ( "  Estimate = %24.16f\n", total );
  printf ( "  Error    = %e\n", error );
  printf ( "  Time     = %f\n", wtime );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "QUAD_OPENMP:\n" );
  printf ( "  Normal end of execution.\n" );

  return 0;
}
/******************************************************************************/

double f ( double x, double y )

/******************************************************************************/
/*
  Purpose:

    F evaluates the function.
*/
{
  double value;

  value = 1.0 / ( 1.0 - x * y );

  return value;
}
/******************************************************************************/

double cpu_time ( )

/******************************************************************************/
/*
  Purpose:
 
    CPU_TIME reports the elapsed CPU time.

  Modified:

    27 September 2005

  Author:

    John Burkardt

  Parameters:

    Output, double CPU_TIME, the current total elapsed CPU time in second.
*/
{
  double value;

  value = ( double ) clock ( ) / ( double ) CLOCKS_PER_SEC;

  return value;
}
