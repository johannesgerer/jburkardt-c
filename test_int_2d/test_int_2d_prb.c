# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "test_int_2d.h"

int main ( void );
void test01 ( void );
void test02 ( void );
void test03 ( void );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for TEST_INT_2D_PRB.

  Discussion:

    TEST_INT_2D_PRB tests the TEST_INT_2D library.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    09 February 2012

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "TEST_INT_2D_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the TEST_INT_2D library.\n" );
 
  test01 ( );
  test02 ( );
  test03 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "TEST_INT_2D_PRB\n" );
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

    TEST01 applies a Monte Carlo rule.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    09 February 2012

  Author:

    John Burkardt
*/
{
  double a[2];
  double b[2];
  int dim;
  double error;
  double exact;
  double *fx;
  int i;
  int j;
  int n;
  int problem;
  int problem_num;
  double quad;
  int seed;
  double volume;
  double *x;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  Use a Monte Carlo rule.\n" );
  printf ( "\n" );
  printf ( "  Repeatedly multiply the number of points by 4.\n" );

  problem_num = p00_problem_num ( );

  printf ( "\n" );
  printf ( "   Problem      Points         Approx            Error\n" );

  for ( problem = 1; problem <= problem_num; problem++ )
  {
    printf ( "\n" );
    n = 1;
    for ( i = 1; i <= 12; i++ )
    {
      seed = 123456789;

      x = ( double * ) malloc ( 2 * n * sizeof ( double ) );
      fx = ( double * ) malloc ( n * sizeof ( double ) );

      x = r8mat_uniform_01 ( 2, n, &seed );

      p00_lim ( problem, a, b );

      for ( dim = 0; dim < 2; dim++ )
      {
        for ( j = 0; j < n; j++ )
        {
          x[dim+j*2] = ( 1.0 - x[dim+j*2] ) * a[dim]
                     +         x[dim+j*2]   * b[dim];
        }
      }
      volume = ( b[1] - a[1] ) * ( b[0] - a[0] );

      p00_fun ( problem, n, x, fx );

      quad = volume * r8vec_sum ( n, fx ) / ( double ) ( n );
   
      exact = p00_exact ( problem );

      error = r8_abs ( quad - exact );

      printf ( "  %8d  %10d  %14g  %14g\n", problem, n, quad, error );

      free ( fx );
      free ( x );

      n = n * 4;
    }
    printf ( "  %8d       Exact  %14g\n", problem, exact );
  }
  return;
}
/******************************************************************************/

void test02 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST02 applies a product of composite midpoint rules.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    09 February 2012

  Author:

    John Burkardt
*/
{
  double a[2];
  double b[2];
  int dim;
  double error;
  double exact;
  double *fx;
  int i;
  int ix;
  int iy;
  int k;
  int n;
  int nx;
  int ny;
  int problem;
  int problem_num;
  double quad;
  double volume;
  double *x;
  double xval;
  double yval;

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  Use a product of composite midpoint rules..\n" );
  printf ( "  Repeatedly multiply the number of points by 4.\n" );

  problem_num = p00_problem_num ( );

  printf ( "\n" );
  printf ( "   Problem      Points         Approx            Error\n" );

  for ( problem = 1; problem <= problem_num; problem++ )
  {
    printf ( "\n" );
    nx = 1;
    ny = 1;

    for ( i = 1; i <= 12; i++ )
    {
      n = nx * ny;

      x = ( double * ) malloc ( 2 * n * sizeof ( double ) );
      fx = ( double * ) malloc ( n * sizeof ( double ) );

      p00_lim ( problem, a, b );

      k = 0;

      for ( ix = 1; ix <= nx; ix++ )
      {
        xval = ( ( double ) ( 2 * nx - 2 * ix + 1 ) * a[0]   
               + ( double ) (          2 * ix - 1 ) * b[0] ) 
               / ( double ) ( 2 * nx              );

        for ( iy = 1; iy <= ny; iy++ )
        {
          yval = ( ( double ) ( 2 * ny - 2 * iy + 1 ) * a[1]
                 + ( double ) (          2 * iy - 1 ) * b[1] ) 
                 / ( double ) ( 2 * ny              );

          x[0+k*2] = xval;
          x[1+k*2] = yval;
          k = k + 1;
        }
      }

      volume = ( b[1] - a[1] ) * ( b[0] - a[0] );

      p00_fun ( problem, n, x, fx );

      quad = volume * r8vec_sum ( n, fx ) / ( double ) ( n );
   
      exact = p00_exact ( problem );

      error = r8_abs ( quad - exact );

      printf ( "  %8d  %10d  %14g  %14g\n", problem, n, quad, error );

      free ( fx );
      free ( x );

      nx = nx * 2;
      ny = ny * 2;
    }
    printf ( "  %8d       Exact  %14g\n", problem, exact );
  }
  return;
}
/******************************************************************************/

void test03 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST03 applies a product of Gauss-Legendre rules.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 February 2012

  Author:

    John Burkardt
*/
{
  double a[2];
  double b[2];
  int dim;
  double error;
  double exact;
  double  *fxy;
  int i;
  int ix;
  int iy;
  int j;
  int k;
  int nx;
  int nxy;
  int ny;
  int problem;
  int problem_num;
  double quad;
  double volume;
  double *w;
  double *wxy;
  double *x;
  double *xy;
  double xval;
  double yval;

  printf ( "\n" );
  printf ( "TEST03\n" );
  printf ( "  Use a product of Gauss-Legendre rules.\n" );
  printf ( "  The 1D rules essentially double in order.\n" );

  problem_num = p00_problem_num ( );

  printf ( "\n" );
  printf ( "   Problem      Points       Approx         Error\n" );

  for ( problem = 1; problem <= problem_num; problem++ )
  {
    printf ( "\n" );

    nx = 1;
    ny = 1;

    for ( i = 1; i <= 8; i++ )
    {
      x = ( double * ) malloc ( nx * sizeof ( double ) );
      w = ( double * ) malloc ( nx * sizeof ( double ) );

      legendre_dr_compute ( nx, x, w );

      nxy = nx * ny;

      wxy = ( double * ) malloc ( nxy * sizeof ( double ) );
      xy = ( double * ) malloc ( 2 * nxy * sizeof ( double ) );
      fxy = ( double * ) malloc ( nxy * sizeof ( double ) );

      p00_lim ( problem, a, b );

      k = 0;

      for ( ix = 0; ix < nx; ix++ )
      {
        xval = ( ( 1.0 + x[ix] ) * a[0]   
               + ( 1.0 - x[ix] ) * b[0] ) 
               /   2.0;
        for ( iy = 0; iy < ny; iy++ )
        {
          yval = ( ( 1.0 + x[iy] ) * a[1]   
                 + ( 1.0 - x[iy] ) * b[1] ) 
                 /   2.0;
          xy[0+k*2] = xval;
          xy[1+k*2] = yval;
          wxy[k] = w[ix] * w[iy];
          k = k + 1;
        }
      }
      volume = ( b[0] - a[0] ) * ( b[1] - a[1] );

      p00_fun ( problem, nxy, xy, fxy );

      quad = 0.0;
      for ( j = 0; j < nxy; j++ )
      {
        quad = quad + wxy[j] * fxy[j];
      }
      quad = quad * volume / 4.0;

      exact = p00_exact ( problem );

      error = r8_abs ( quad - exact );

      printf ( "  %8d  %10d  %14g  %14g\n", problem, nxy, quad, error );

      free ( fxy );
      free ( w );
      free ( wxy );
      free ( x );
      free ( xy );

      nx = 2 * nx + 1;
      ny = nx;
    }
    printf ( "  %8d       Exact  %14g\n", problem, exact );
  }
  return;
}
