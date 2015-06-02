# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "power_method.h"

int main ( void );
void test01 ( void );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for POWER_METHOD_PRB.

  Discussion:

    POWER_METHOD_PRB tests the POWER_METHOD library.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    28 May 2008

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "POWER_METHOD_PRB\n" );
  printf ( "  C version:\n" );
  printf ( "  Test the POWER_METHOD library.\n" );

  test01 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "POWER_METHOD_PRB:\n" );
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

    TEST01 uses POWER_METHOD on the Fibonacci2 matrix.

  Discussion:

    This matrix, despite having a single dominant eigenvalue, will generally
    converge only very slowly under the power method.  This has to do with
    the fact that the matrix has only 3 eigenvectors.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    20 July 2008

  Author:

    John Burkardt
*/
{
  double *a;
  double cos_x1x2;
  double ctime;
  double ctime1;
  double ctime2;
  int i;
  int it_max;
  int it_num;
  double lambda;
  int n = 50;
  double norm;
  double phi;
  int seed;
  double sin_x1x2;
  double tol;
  double *x;
  double *x2;

  a = fibonacci2 ( n );

  seed = 123456789;
  x = r8vec_uniform_01 ( n, &seed );

  it_max = 300;
  tol = 0.000001;

  phi = ( 1.0 + sqrt ( 5.0 ) ) / 2.0;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  Use the power method on the Fibonacci2 matrix.\n" );
  printf ( "\n" );
  printf ( "  Matrix order N       = %d\n", n );
  printf ( "  Maximum iterations   = %d\n", it_max );
  printf ( "  Error tolerance      = %e\n", tol );

  ctime1 = cpu_time ( );

  power_method ( n, a, x, it_max, tol, &lambda, &it_num );

  ctime2 = cpu_time ( );
  ctime = ctime2 - ctime1;

  printf ( "\n" );
  printf ( "  Number of iterations = %d\n", it_num );
  printf ( "  CPU time             = %f\n", ctime );
  printf ( "  Estimated eigenvalue = %24.16f\n", lambda );
  printf ( "  Correct value        = %24.16f\n", phi );
  printf ( "  Error                = %e\n", r8_abs ( lambda - phi ) );
/*
  X2 is the exact eigenvector.
*/
  x2 = ( double * ) malloc ( n * sizeof ( double ) );

  x2[0] = 1.0;
  for ( i = 1; i < n; i++ )
  {
    x2[i] = phi * x2[i-1];
  }
  norm = r8vec_norm_l2 ( n, x2 );
  for ( i = 0; i < n; i++ )
  {
    x2[i] = x2[i] / norm;
  }
/*
  The sine of the angle between X and X2 is a measure of error.
*/
  cos_x1x2 = r8vec_dot ( n, x, x2 );
  sin_x1x2 = sqrt ( ( 1.0 - cos_x1x2 ) * ( 1.0 + cos_x1x2 ) );

  printf ( "\n" );
  printf ( "  Sine of angle between true and estimated vectors = %e\n", sin_x1x2 );

  free ( a );
  free ( x );
  free ( x2 );

  return;
}
