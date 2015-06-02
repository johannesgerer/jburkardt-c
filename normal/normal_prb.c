# include <complex.h>
# include <stdlib.h>
# include <stdio.h>

# include "normal.h"

int main ( );
void c4_normal_01_test ( );
void c8_normal_01_test ( );
void i4_normal_ab_test ( );
void i8_normal_ab_test ( );
void r4_normal_01_test ( );
void r4_normal_ab_test ( );
void r8_normal_01_test ( );
void r8_normal_ab_test ( );
void r8mat_normal_01_new_test ( );
void r8vec_normal_01_new_test ( );

/**********************************************************************/

int main ( )

/**********************************************************************/
/*
  Purpose:

    MAIN is the main program for NORMAL_PRB.

  Discussion:

    NORMAL_PRB tests the NORMAL library.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    03 March 2015

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "NORMAL_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the NORMAL library.\n" );

  c4_normal_01_test ( );
  c8_normal_01_test ( );
  i4_normal_ab_test ( );
  i8_normal_ab_test ( );
  r4_normal_01_test ( );
  r4_normal_ab_test ( );
  r8_normal_01_test ( );
  r8_normal_ab_test ( );
  r8mat_normal_01_new_test ( );
  r8vec_normal_01_new_test ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "NORMAL_PRB\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/**********************************************************************/

void c4_normal_01_test ( )

/**********************************************************************/
/*
  Purpose:

    C4_NORMAL_01_TEST tests C4_NORMAL_01.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    03 March 2015

  Author:

    John Burkardt
*/
{
  int i;
  float complex r;
  int seed;

  printf ( "\n" );
  printf ( "C4_NORMAL_01_TEST\n" );
  printf ( "  C4_NORMAL_01 computes pseudonormal float complex values\n" );
  printf ( "  with mean 0.0 and standard deviation 1.0.\n" );

  seed = 123456789;

  printf ( "\n" );
  printf ( "  SEED = %d\n", seed );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    r = c4_normal_01 ( &seed );
    printf ( "  %6d  %14g,%14g\n", i, creal ( r ), cimag ( r ) );
  }

  return;
}
/**********************************************************************/

void c8_normal_01_test ( )

/**********************************************************************/
/*
  Purpose:

    C8_NORMAL_01_TEST tests C8_NORMAL_01.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    04 March 2015

  Author:

    John Burkardt
*/
{
  int i;
  double complex r;
  int seed;

  printf ( "\n" );
  printf ( "C8_NORMAL_01_TEST\n" );
  printf ( "  C8_NORMAL_01 computes pseudonormal double complex values\n" );
  printf ( "  with mean 0.0 and standard deviation 1.0.\n" );

  seed = 123456789;

  printf ( "\n" );
  printf ( "  SEED = %d\n", seed );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    r = c8_normal_01 ( &seed );
    printf ( "  %6d  %14g,%14g\n", i, creal ( r ), cimag ( r ) );
  }

  return;
}
/**********************************************************************/

void i4_normal_ab_test ( )

/**********************************************************************/
/*
  Purpose:

    I4_NORMAL_AB_TEST tests I4_NORMAL_AB.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    02 July 2006

  Author:

    John Burkardt
*/
{
  int i;
  float mu;
  int r;
  int seed;
  float sigma;

  printf ( "\n" );
  printf ( "I4_NORMAL_AB_TEST\n" );
  printf ( "  I4_NORMAL_AB computes pseudonormal integers\n" );
  printf ( "  with mean MU and standard deviation SIGMA.\n" );

  mu = 70.0;
  sigma = 10.0;
  seed = 123456789;

  printf ( "\n" );
  printf ( "  The mean = %f\n", mu );
  printf ( "  The standard deviation = %f\n", sigma );
  printf ( "  SEED = %d\n", seed );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    r = i4_normal_ab ( mu, sigma, &seed );
    printf ( "  %8d  %8d\n", i, r );
  }

  return;
}
/**********************************************************************/

void i8_normal_ab_test ( )

/**********************************************************************/
/*
  Purpose:

    I8_NORMAL_AB_TEST tests I8_NORMAL_AB.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    02 July 2006

  Author:

    John Burkardt
*/
{
  int i;
  double mu;
  long int r;
  long int seed;
  double sigma;

  printf ( "\n" );
  printf ( "I8_NORMAL_AB_TEST\n" );
  printf ( "  I8_NORMAL_AB computes pseudonormal integers\n" );
  printf ( "  with mean MU and standard deviation SIGMA.\n" );

  mu = 70.0;
  sigma = 10.0;
  seed = 123456789;

  printf ( "\n" );
  printf ( "  The mean = %f\n", mu );
  printf ( "  The standard deviation = %f\n", sigma );
  printf ( "  SEED = %d\n", seed );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    r = i8_normal_ab ( mu, sigma, &seed );
    printf ( "  %8d  %8d\n", i, r );
  }

  return;
}
/**********************************************************************/

void r4_normal_01_test ( )

/**********************************************************************/
/*
  Purpose:

    R4_NORMAL_01_TEST tests R4_NORMAL_01.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    29 June 2006

  Author:

    John Burkardt
*/
{
  int i;
  float r;
  int seed;

  printf ( "\n" );
  printf ( "R4_NORMAL_01_TEST\n" );
  printf ( "  R4_NORMAL_01 computes pseudonormal values\n" );
  printf ( "  with mean 0.0 and standard deviation 1.0.\n" );

  seed = 123456789;

  printf ( "\n" );
  printf ( "  SEED = %d\n", seed );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    r = r4_normal_01 ( &seed );
    printf ( "  %6d  %14f\n", i, r );
  }

  return;
}
/**********************************************************************/

void r4_normal_ab_test ( )

/**********************************************************************/
/*
  Purpose:

    R4_NORMAL_AB_TEST tests R4_NORMAL_AB.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    29 June 2006

  Author:

    John Burkardt
*/
{
  int i;
  float mu;
  float r;
  int seed;
  float sigma;

  printf ( "\n" );
  printf ( "R4_NORMAL_AB_TEST\n" );
  printf ( "  R4_NORMAL_AB computes pseudonormal values\n" );
  printf ( "  with mean MU and standard deviation SIGMA.\n" );

  mu = 10.0;
  sigma = 2.0;
  seed = 123456789;

  printf ( "\n" );
  printf ( "  The mean = %f\n", mu );
  printf ( "  The standard deviation = %f\n", sigma );
  printf ( "  SEED = %d\n", seed );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    r = r4_normal_ab ( mu, sigma, &seed );
    printf ( "  %6d  %14f\n", i, r );
  }

  return;
}
/**********************************************************************/

void r8_normal_01_test ( )

/**********************************************************************/
/*
  Purpose:

    R8_NORMAL_01_TEST tests R8_NORMAL_01.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    29 June 2006

  Author:

    John Burkardt
*/
{
  int i;
  double r;
  int seed;

  printf ( "\n" );
  printf ( "R8_NORMAL_01_TEST\n" );
  printf ( "  R8_NORMAL_01 computes pseudonormal values\n" );
  printf ( "  with mean 0.0 and standard deviation 1.0.\n" );

  seed = 123456789;

  printf ( "\n" );
  printf ( "  SEED = %d\n", seed );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    r = r8_normal_01 ( &seed );
    printf ( "  %6d  %14f\n", i, r );
  }

  return;
}
/**********************************************************************/

void r8_normal_ab_test ( )

/**********************************************************************/
/*
  Purpose:

    R8_NORMAL_AB_TEST tests R8_NORMAL_AB.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    29 June 2006

  Author:

    John Burkardt
*/
{
  int i;
  double mu;
  double r;
  int seed;
  double sigma;

  printf ( "\n" );
  printf ( "R8_NORMAL_AB_TEST\n" );
  printf ( "  R8_NORMAL_AB computes pseudonormal values\n" );
  printf ( "  with mean MU and standard deviation SIGMA.\n" );

  mu = 10.0;
  sigma = 2.0;
  seed = 123456789;

  printf ( "\n" );
  printf ( "  The mean = %f\n", mu );
  printf ( "  The standard deviation = %f\n", sigma );
  printf ( "  SEED = %d\n", seed );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    r = r8_normal_ab ( mu, sigma, &seed );
    printf ( "  %6d  %14f\n", i, r );
  }

  return;
}
/**********************************************************************/

void r8mat_normal_01_new_test ( )

/**********************************************************************/
/*
  Purpose:

    R8MAT_NORMAL_01_NEW_TEST tests R8MAT_NORMAL_01_NEW.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    29 June 2006

  Author:

    John Burkardt
*/
{
  int m = 5;
  int n = 4;
  double *r;
  int seed;

  printf ( "\n" );
  printf ( "R8MAT_NORMAL_01_NEW_TEST\n" );
  printf ( "  R8MAT_NORMAL_01_NEW computes a matrix of Normal 01 values.\n" );

  seed = 123456789;
  printf ( "\n" );
  printf ( "  SEED = %d\n", seed );

  r = r8mat_normal_01_new ( m, n, &seed );

  r8mat_print ( m, n, r, "  Matrix:" );
  
  free ( r );

  return;
}
/**********************************************************************/

void r8vec_normal_01_new_test ( )

/**********************************************************************/
/*
  Purpose:

   R8VEC_NORMAL_01_NEW_TEST tests R8VEC_NORMAL_01_NEW.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    29 June 2006

  Author:

    John Burkardt
*/
{
  int n = 10;
  double *r;
  int seed;

  printf ( "\n" );
  printf ( "R8VEC_NORMAL_01_NEW_TEST:\n" );
  printf ( "  R8VEC_NORMAL_01_NEW computes a vector of Normal 01 values.\n" );

  seed = 123456789;
  printf ( "\n" );
  printf ( "  SEED = %d\n", seed );

  r = r8vec_normal_01_new ( n, &seed );

  r8vec_print ( n, r, "  Vector:" );
  
  free ( r );

  return;
}
