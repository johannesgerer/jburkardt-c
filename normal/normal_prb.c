# include <stdlib.h>
# include <stdio.h>

# include "normal.h"

int main ( void );
void test03 ( void );
void test04 ( void );
void test05 ( void );
void test06 ( void );
void test07 ( void );
void test08 ( void );
void test09 ( void );
void test10 ( void );
void test11 ( void );
void test12 ( void );

/**********************************************************************/

int main ( void )

/**********************************************************************/
/*
  Purpose:

    NORMAL_PRB calls sample problems for the NORMAL library.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    08 January 2007

  Author:

    John Burkardt
*/
{
  timestamp ( );

  printf ( "\n" );
  printf ( "NORMAL_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the NORMAL library.\n" );

  test03 ( );
  test04 ( );
  test05 ( );
  test06 ( );
  test07 ( );
  test08 ( );
  test09 ( );

  test10 ( );
  test11 ( );
  test12 ( );
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

void test03 ( void )

/**********************************************************************/
/*
  Purpose:

    TEST03 tests I4_NORMAL.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    02 July 2006

  Author:

    John Burkardt
*/
{
  float a;
  float b;
  int i;
  int seed;
  int seed_init = 123456789;

  printf ( "\n" );
  printf ( "TEST03\n" );
  printf ( "  I4_NORMAL computes pseudonormal integers\n" );
  printf ( "  with mean A and standard deviation B.\n" );

  a = 70.0;
  b = 10.0;
  seed = seed_init;

  printf ( "\n" );
  printf ( "  The mean A = %f\n", a );
  printf ( "  The standard deviation B = %f\n", b );
  printf ( "  The initial seed is %d\n", seed_init );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    printf ( "  %8d  %8d\n", i, i4_normal ( a, b, &seed ) );
  }

  return;
}
/**********************************************************************/

void test04 ( void )

/**********************************************************************/
/*
  Purpose:

    TEST04 tests I8_NORMAL.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    02 July 2006

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  long int i;
  long int seed;
  long int seed_init = 123456789;

  printf ( "\n" );
  printf ( "TEST04\n" );
  printf ( "  I8_NORMAL computes pseudonormal integers\n" );
  printf ( "  with mean A and standard deviation B.\n" );

  a = 70.0;
  b = 10.0;
  seed = seed_init;

  printf ( "\n" );
  printf ( "  The mean A = %f\n", a );
  printf ( "  The standard deviation B = %f\n", b );
  printf ( "  The initial seed is %d\n", seed_init );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    printf ( "  %8d  %8d\n", i, i8_normal ( a, b, &seed ) );
  }

  return;
}
/**********************************************************************/

void test05 ( void )

/**********************************************************************/
/*
  Purpose:

    TEST05 tests R4_NORMAL.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    29 June 2006

  Author:

    John Burkardt
*/
{
  float a;
  float b;
  int i;
  int seed;
  int seed_init = 123456789;

  printf ( "\n" );
  printf ( "TEST05\n" );
  printf ( "  R4_NORMAL computes pseudonormal values\n" );
  printf ( "  with mean A and standard deviation B.\n" );

  a = 10.0;
  b = 2.0;
  seed = seed_init;

  printf ( "\n" );
  printf ( "  The mean A = %f\n", a );
  printf ( "  The standard deviation B = %f\n", b );
  printf ( "  The initial seed is %d\n", seed_init );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    printf ( "  %6d  %14f\n", i, r4_normal ( a, b, &seed ) );
  }

  return;
}
/**********************************************************************/

void test06 ( void )

/**********************************************************************/
/*
  Purpose:

    TEST06 tests R4_NORMAL_01.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    29 June 2006

  Author:

    John Burkardt
*/
{
  int i;
  int seed;
  int seed_init = 123456789;

  printf ( "\n" );
  printf ( "TEST06\n" );
  printf ( "  R4_NORMAL_01 computes pseudonormal values\n" );
  printf ( "  with mean 0 and standard deviation 1.\n" );

  seed = seed_init;

  printf ( "\n" );
  printf ( "  The initial seed is %d\n", seed_init );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    printf ( "  %6d  %14f\n", i, r4_normal_01 ( &seed ) );
  }

  return;
}
/**********************************************************************/

void test07 ( void )

/**********************************************************************/
/*
  Purpose:

    TEST07 tests R8_NORMAL.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    29 June 2006

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  int i;
  int seed;
  int seed_init = 123456789;

  printf ( "\n" );
  printf ( "TEST07\n" );
  printf ( "  R8_NORMAL computes pseudonormal values\n" );
  printf ( "  with mean A and standard deviation B.\n" );

  a = 10.0;
  b = 2.0;
  seed = seed_init;

  printf ( "\n" );
  printf ( "  The mean A = %f\n", a );
  printf ( "  The standard deviation B = %f\n", b );
  printf ( "  The initial seed is %d\n", seed_init );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    printf ( "  %6d  %14f\n", i, r8_normal ( a, b, &seed ) );
  }

  return;
}
/**********************************************************************/

void test08 ( void )

/**********************************************************************/
/*
  Purpose:

    TEST08 tests R8_NORMAL_01.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    29 June 2006

  Author:

    John Burkardt
*/
{
  int i;
  int seed;
  int seed_init = 123456789;

  printf ( "\n" );
  printf ( "TEST08\n" );
  printf ( "  R8_NORMAL_01 computes pseudonormal values\n" );
  printf ( "  with mean 0 and standard deviation 1.\n" );

  seed = seed_init;

  printf ( "\n" );
  printf ( "  The initial seed is %d\n", seed_init );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    printf ( "  %6d  %14f\n", i, r8_normal_01 ( &seed ) );
  }

  return;
}
/******************************************************************************/

void test09 ( void )

/******************************************************************************/
/*

  Purpose:

    TEST09 tests R8_NORMAL_01.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    08 January 2007

  Author:

    John Burkardt
*/
{
  int i;
  int seed;
  int seed_init = 123456789;
  int seed_input;
  int seed_output;
  double value;

  printf ( "\n" );
  printf ( "TEST09\n" );
  printf ( "  R8_NORMAL_01 computes pseudonormal values\n" );
  printf ( "  with mean 0 and standard deviation 1.\n" );
  printf ( "\n" );
  printf ( "  Verify that we can change the seed\n" );
  printf ( "  and get the desired results.\n" );
  printf ( "\n" );
  printf ( "  The initial seed is %d\n", seed_init );

  seed = seed_init;

  printf ( "\n" );
  printf ( "         I    Seed(in)   Seed(out)   R8_NORMAL_01\n" );
  printf ( "\n" );

  for ( i = 1; i <= 5; i++ )
  {
    seed_input = seed;
    value = r8_normal_01 ( &seed );
    seed_output = seed;

    printf ( "  %8d  %12d  %12d  %14f\n",
      i, seed_input, seed_output, value );
  }

  seed = seed_init;

  printf ( "\n" );
  printf ( "  Resetting seed to repeat, after an ODD number of steps.\n" );
  printf ( "\n" );

  for ( i = 6; i <= 10; i++ )
  {
    seed_input = seed;
    value = r8_normal_01 ( &seed );
    seed_output = seed;

    printf ( "  %8d  %12d  %12d  %14f\n",
      i, seed_input, seed_output, value );
  }

  seed = seed_init;

  printf ( "\n" );
  printf ( "  Resetting seed to repeat, after an EVEN number of steps.\n" );
  printf ( "\n" );


  for ( i = 11; i <= 15; i++ )
  {
    seed_input = seed;
    value = r8_normal_01 ( &seed );
    seed_output = seed;

    printf ( "  %8d  %12d  %12d  %14f\n",
      i, seed_input, seed_output, value );
  }

  return;
}
/**********************************************************************/

void test10 ( void )

/**********************************************************************/
/*
  Purpose:

   TEST10 tests R8_NORMAL_01;

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    29 June 2006

  Author:

    John Burkardt
*/
{
# define N 1000

  int i;
  int seed;
  int seed_in;
  int seed_out;
  double u[N];
  double u_avg;
  double u_var;

  printf ( "\n" );
  printf ( "TEST10\n" );
  printf ( "  R8_NORMAL_01 computes a sequence of normally distributed\n" );
  printf ( "  pseudorandom numbers.\n" );

  seed = 12345;

  printf ( "\n" );
  printf ( "  Initial SEED = %d\n", seed );

  printf ( "\n" );
  printf ( "  First 10 values:\n" );
  printf ( "\n" );
  printf ( "       I         Input        Output   R8_NORMAL_01\n" );
  printf ( "                  SEED          SEED\n" );
  printf ( "\n" );

  for ( i = 0; i < 10; i++ )
  {
    seed_in = seed;
    u[i] = r8_normal_01 ( &seed );
    seed_out = seed;

    printf ( "  %6d  %12d  %12d  %10f\n", i + 1, seed_in, seed_out, u[i] );
  }

  printf ( "\n" );
  printf ( "  Now call R8_NORMAL_01 %d times.\n", N );

  u_avg = 0.0;
  for ( i = 0; i < N; i++ )
  {
    u[i] = r8_normal_01 ( &seed );
    u_avg = u_avg + u[i];
  }

  u_avg = u_avg / ( ( double ) N );

  u_var = 0.0;
  for ( i = 0; i < N; i++ )
  {
    u_var = u_var + ( u[i] - u_avg ) * ( u[i] - u_avg );
  }
  u_var = u_var / ( ( double ) ( N - 1 ) );

  printf ( "\n" );
  printf ( "  Average value = %f\n", u_avg );
  printf ( "  Expecting       %f\n", 0.0 );

  printf ( "\n" );
  printf ( "  Variance =      %f\n", u_var );
  printf ( "  Expecting       %f\n", 1.0 );

  return;
# undef N
}
/**********************************************************************/

void test11 ( void )

/**********************************************************************/
/*
  Purpose:

    TEST11 tests R8_NORMAL_01 and R8MAT_NORMAL_01_NEW.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    29 June 2006

  Author:

    John Burkardt
*/
{
# define M 100
# define N 10

  double a[M*N];
  double *b;
  int i;
  int j;
  int k;
  int seed;
  int seed_init = 123456789;

  printf ( "\n" );
  printf ( "TEST11\n" );
  printf ( "  R8_NORMAL_01 computes pseudonormal values one at a time.\n" );
  printf ( "  R8MAT_NORMAL_01_NEW computes a matrix of values.\n" );
  printf ( "\n" );
  printf ( "  For the same initial seed, the results should be identical,\n" );
  printf ( "  but R8MAT_NORMAL_01_NEW might be faster.\n" );
  printf ( "\n" );
  printf ( "  Initial seed is %d\n", seed_init );

  seed = seed_init;
  for ( j = 0; j < N; j++ )
  {
    for ( i = 0; i < M; i++ )
    {
      a[i+j*M] = r8_normal_01 ( &seed );
    }
  }

  seed = seed_init;
  b = r8mat_normal_01_new ( M, N, &seed );

  printf ( "\n" );
  printf ( "      I       J      A[I,J]        B[I,J]\n" );
  printf ( "                 (R8_NORMAL_01)  (R8MAT_NORMAL_01)\n" );
  printf ( "\n" );

  for ( k = 0; k < 11; k++ )
  {
    i = ( k * ( M - 1 ) ) / 10;
    j = ( k * ( N - 1 ) ) / 10;

    printf ( "  %6d  %6d  %12f  %12f\n", i, j, a[i+j*M], b[i+j*M] );
  }
  
  free ( b );

  return;
# undef M
# undef N
}
/**********************************************************************/

void test12 ( void )

/**********************************************************************/
/*
  Purpose:

   TEST12 tests R8_NORMAL_01 and R8VEC_NORMAL_01_NEW.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    29 June 2006

  Author:

    John Burkardt
*/
{
# define N 10

  double a[N];
  double *b;
  int i;
  int j;
  int k;
  int seed;
  int seed_init = 123456789;

  printf ( "\n" );
  printf ( "TEST12\n" );
  printf ( "  R8_NORMAL_01 computes pseudeonormal values one at a time.\n" );
  printf ( "  R8VEC_NORMAL_01_NEW computes a vector of values.\n" );
  printf ( "\n" );
  printf ( "  For the same initial seed, the results should be identical,\n" );
  printf ( "  but R8VEC_NORMAL_01_NEW might be faster.\n" );
  printf ( "\n" );
  printf ( "  Initial seed is %d\n", seed_init );

  seed = seed_init;
  for ( i = 0; i < N; i++ )
  {
    a[i] = r8_normal_01 ( &seed );
  }

  seed = seed_init;
  b = r8vec_normal_01_new ( N, &seed );

  printf ( "\n" );
  printf ( "      I      A[I]          B[I]\n" );
  printf ( "         (R8_NORMAL_01)  (R8VEC_NORMAL_01_NEW)\n" );
  printf ( "\n" );

  for ( i = 1; i < N; i++ )
  {
    printf ( "  %6d  %12f  %12f\n", i, a[i], b[i] );
  }
  
  free ( b );

  return;
# undef N
}
