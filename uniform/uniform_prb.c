# include <stdlib.h>
# include <stdio.h>

# include "uniform.h"

int main ( void );
void test01 ( void );
void test02 ( void );
void test03 ( void );
void test04 ( void );
void test05 ( void );
void test06 ( void );
void test065 ( void );
void test07 ( void );
void test08 ( void );
void test09 ( void );

void test10 ( void );
void test11 ( void );
void test111 ( void );
void test112 ( void );
void test118 ( void );
void test119 ( void );
void test12 ( void );
void test13 ( void );
void test14 ( void );
void test15 ( void );
void test16 ( void );
void test17 ( void );
void test18 ( void );
void test19 ( void );
void test20 ( void );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    UNIFORM_PRB calls sample problems for the UNIFORM library.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    21 May 2008

  Author:

    John Burkardt
*/
{
  timestamp ( );

  printf ( "\n" );
  printf ( "UNIFORM_PRB\n" );
  printf ( "  C version:\n" );
  printf ( "  Test the UNIFORM library.\n" );

  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
  test05 ( );
  test06 ( );
  test065 ( );
  test07 ( );
  test08 ( );
  test09 ( );

  test10 ( );
  test11 ( );
  test111 ( );
  test112 ( );
  test118 ( );
  test119 ( );
  test12 ( );
  test13 ( );
  test14 ( );
  test15 ( );
  test16 ( );
  test17 ( );
  test18 ( );
  test19 ( );

  test20 ( );

  printf ( "\n" );
  printf ( "UNIFORM_PRB\n" );
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

    TEST01 tests C4_UNIFORM_01.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 January 2007

  Author:

    John Burkardt
*/
{
  int i;
  int seed;
  int seed_init = 123456789;
  complex value;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  C4_UNIFORM_01 computes pseudorandom complex values\n" );
  printf ( "  uniformly distributed in the unit circle.\n" );

  seed = seed_init;

  printf ( "\n" );
  printf ( "  The initial seed is %d\n", seed_init );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    value = c4_uniform_01 ( &seed );

    printf ( "  %6d  %12f  %12f\n", i, value.real, value.imag );
  }

  return;
}
/******************************************************************************/

void test02 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST02 tests C4VEC_UNIFORM_01_NEW.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 January 2007

  Author:

    John Burkardt
*/
{
  complex *cvec;
  int i;
  int n;
  int seed;
  int seed_init = 123456789;

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  C4VEC_UNIFORM_01_NEW computes pseudorandom complex values\n" );
  printf ( "  uniformly distributed in the unit circle.\n" );

  seed = seed_init;

  printf ( "\n" );
  printf ( "  The initial seed is %d\n", seed_init );
  printf ( "\n" );

  n = 10;
  cvec = c4vec_uniform_01_new ( n, &seed );

  for ( i = 0; i < n; i++ )
  {
    printf ( "  %8d  %12f  %12f\n", i, cvec[i].real, cvec[i].imag );
  }

  free ( cvec );

  return;
}
/******************************************************************************/

void test03 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST03 tests C8_UNIFORM_01.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 January 2007

  Author:

    John Burkardt
*/
{
  int i;
  int seed;
  int seed_init = 123456789;
  doublecomplex value;

  printf ( "\n" );
  printf ( "TEST03\n" );
  printf ( "  C8_UNIFORM_01 computes pseudorandom C8 values\n" );
  printf ( "  uniformly distributed in the unit circle.\n" );

  seed = seed_init;

  printf ( "\n" );
  printf ( "  The initial seed is %d\n", seed_init );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    value = c8_uniform_01 ( &seed );

    printf ( "  %6d  %12f  %12f\n", i, value.real, value.imag );
  }

  return;
}
/******************************************************************************/

void test04 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST04 tests C8VEC_UNIFORM_01_NEW.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 January 2007

  Author:

    John Burkardt
*/
{
  doublecomplex *cvec;
  int i;
  int n;
  int seed;
  int seed_init = 123456789;

  printf ( "\n" );
  printf ( "TEST04\n" );
  printf ( "  C8VEC_UNIFORM_01_NEW computes pseudorandom C8 values\n" );
  printf ( "  uniformly distributed in the unit circle.\n" );

  seed = seed_init;

  printf ( "\n" );
  printf ( "  The initial seed is %d\n", seed_init );
  printf ( "\n" );

  n = 10;
  cvec = c8vec_uniform_01_new ( n, &seed );

  for ( i = 0; i < n; i++ )
  {
    printf ( "  %8d  %12f  %12f\n", i, cvec[i].real, cvec[i].imag );
  }

  free ( cvec );

  return;
}
/******************************************************************************/

void test05 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST05 tests CH_UNIFORM_AB.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    29 June 2006

  Author:

    John Burkardt
*/
{
  char chi;
  char clo;
  int i;
  int seed;
  int seed_init = 123456789;

  printf ( "\n" );
  printf ( "TEST05\n" );
  printf ( "  CH_UNIFORM_AB computes pseudorandom characters\n" );
  printf ( "  in an interval [CLO,CHI].\n" );

  clo = 'A';
  chi = 'J';
  seed = seed_init;

  printf ( "\n" );
  printf ( "  The lower endpoint CLO = '%c'.\n", clo );
  printf ( "  The upper endpoint CHI = '%c'.\n", chi );
  printf ( "  The initial seed is %d\n", seed_init );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    printf ( "  %6d  '%c'\n", i, ch_uniform_ab ( clo, chi, &seed ) );
  }

  return;
}
/******************************************************************************/

void test06 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST06 tests GET_SEED.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    29 June 2006

  Author:

    John Burkardt
*/
{
  int i;
  int j;
  int oops;
  int oops_max = 10000;
  int seed;
  int seed_old;

  printf ( "\n" );
  printf ( "TEST06\n" );
  printf ( "  GET_SEED picks an initial seed value for R8_UNIFORM_01.\n" );
  printf ( "  The value chosen should vary over time, because\n" );
  printf ( "  the seed is based on reading the clock.\n" );
  printf ( "\n" );
  printf ( "  This is just the \"calendar\" clock, which does\n" );
  printf ( "  not change very fast, so calling GET_SEED several\n" );
  printf ( "  times in a row may result in the same value.\n" );

  oops = 0;
  seed = 12345678;
  seed_old = seed;

  printf ( "\n" );
  printf ( "  Initial seed is %d\n", seed );
  printf ( "\n" );
  printf ( "  Next 3 values of R8_UNIFORM_01:\n" );
  printf ( "\n" );

  for ( j = 1; j <= 3; j++ )
  {
    printf ( "  %10f\n", r8_uniform_01 ( &seed ) );
  }

  for ( i = 1; i <= 4; i++ )
  {
    while ( 1 )
    {
      seed = get_seed ( );

      if ( seed = seed_old )
      {
        seed_old = seed;
        break;
      }
      oops = oops + 1;
      if ( oops_max < oops ) 
      {
        printf ( "\n" );
        printf ( "  Oops\n" );
        printf ( "  Same seed returned for %d calls to GET_SEED\n", oops_max );
        printf ( "  Could be a bad algorithm, slow clock, or fast machine\n" );
        printf ( "  To avoid infinite loops, we take what we have now.\n" );
        oops = 0;
        oops_max = oops_max * 10;
        break;
      }
    }

    printf ( "\n" );
    printf ( "  New seed from GET_SEED is %d\n", seed );
    printf ( "\n" );
    printf ( "  Next 3 values of R8_UNIFORM_01:\n" );
    printf ( "\n" );

    for ( j = 1; j <= 3; j++ )
    {
      printf ( "  %10f\n", r8_uniform_01 ( &seed ) );
    }

  }
  return;
}
/******************************************************************************/

void test065 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST065 tests I4_SEED_ADVANCE.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    21 May 2008

  Author:

    John Burkardt
*/
{
  int seed;
  int seed_new;
  int step;

  printf ( "\n" );
  printf ( "TEST065\n" );
  printf ( "  I4_SEED_ADVANCE advances the seed.\n" );
  printf ( "\n" );
  printf ( "  Step        SEED input       SEED output\n" );
  printf ( "\n" );

  seed_new = 12345;

  for ( step = 1; step <= 10; step++ )
  {
    seed = seed_new;
    seed_new = i4_seed_advance ( seed );

    printf ( "  %4d  %16d  %16d\n", step, seed, seed_new );
  }

  return;
}
/******************************************************************************/

void test07 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST07 tests I4_UNIFORM_AB.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    29 June 2006

  Author:

    John Burkardt
*/
{
# define A 6
# define B 10

  int a = A;
  int b = B;
  int freq[B+1-A];
  int i;
  int j;
  int seed;
  int seed_init = 123456789;

  printf ( "\n" );
  printf ( "TEST07\n" );
  printf ( "  I4_UNIFORM_AB computes pseudorandom values\n" );
  printf ( "  in an interval [A,B].\n" );

  seed = seed_init;

  printf ( "\n" );
  printf ( "  The lower endpoint A = %d\n", a );
  printf ( "  The upper endpoint B = %d\n", b );
  printf ( "  The initial seed is %d\n", seed_init );
  printf ( "\n" );

  for ( i = a; i <= b; i++ )
  {
    freq[i-a] = 0;
  }
  for ( i = 1; i <= 10000; i++ )
  {
    j = i4_uniform_ab ( a, b, &seed );
    if ( j < a ) 
    {
      printf ( "  Illegal value J = %d\n", j );
    }
    else if ( j <= b )
    {
      freq[j-a] = freq[j-a] + 1;
    }
    else
    {
      printf ( "  Illegal value J = %d\n", j );
    }
  }

  printf ( "\n" );
  printf ( "         I    Frequency\n" );
  printf ( "\n" );
  for ( i = a; i <= b; i++ )
  {
    printf ( "  %8d  %8d\n", i, freq[i-a] );
  }

  return;
# undef A
# undef B
}
/******************************************************************************/

void test08 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST08 tests I4_UNIFORM_AB.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 January 2007

  Author:

    John Burkardt
*/
{
  int a = -100;
  int b = 200;
  int i;
  int j;
  int seed;
  int seed_init = 123456789;

  printf ( "\n" );
  printf ( "TEST08\n" );
  printf ( "  I4_UNIFORM_AB computes pseudorandom values\n" );
  printf ( "  in an interval [A,B].\n" );

  seed = seed_init;

  printf ( "\n" );
  printf ( "  The lower endpoint A = %d\n", a );
  printf ( "  The upper endpoint B = %d\n", b );
  printf ( "  The initial seed is %d\n", seed_init );
  printf ( "\n" );

  for ( i = 1; i <= 20; i++ )
  {
    j = i4_uniform_ab ( a, b, &seed );

    printf ( "  %8d  %d\n", i, j );
  }

  return;
}
/******************************************************************************/

void test09 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST09 tests I4_UNIFORM_0I

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
  float mean;
  int seed;
  float variance;
  int x[N];

  printf ( "\n" );
  printf ( "TEST09\n" );
  printf ( "  I4_UNIFORM_0I samples a uniform random\n" );
  printf ( "  integer distribution in [0,2**31-1].\n" );

  seed = 123456789;

  printf ( "\n" );
  printf ( "  Starting with seed = %d\n", seed );

  for ( i = 0; i < N; i++ )
  {
    x[i] = i4_uniform_0i ( &seed );
  }

  printf ( "\n" );
  printf ( "  First few values:\n" );
  printf ( "\n" );
  for ( i = 0; i < 5; i++ )
  {
    printf ( "  %6d  %6d\n", i, x[i] );
  }

  mean = i4vec_mean ( N, x );

  variance = i4vec_variance ( N, x );

  printf ( "\n" );
  printf ( "  Number of values computed was N = %d\n", N );
  printf ( "  Average value was %f\n", mean );
  printf ( "  Minimum value was %d\n", i4vec_min ( N, x ) );
  printf ( "  Maximum value was %d\n", i4vec_max ( N, x ) );
  printf ( "  Variance was %f\n", variance );

  return;
# undef N
}
/******************************************************************************/

void test10 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST10 tests I4VEC_UNIFORM_AB_NEW.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    29 June 2006

  Author:

    John Burkardt
*/
{
# define A 6
# define B 10
# define N 10000

  int a = A;
  int b = B;
  int freq[B+1-A];
  int i;
  int *i4vec;
  int j;
  int seed;
  int seed_init = 123456789;

  printf ( "\n" );
  printf ( "TEST10\n" );
  printf ( "  I4VEC_UNIFORM_AB_NEW computes pseudorandom values\n" );
  printf ( "  in an interval [A,B].\n" );

  seed = seed_init;

  printf ( "\n" );
  printf ( "  The lower endpoint A = %d\n", a );
  printf ( "  The upper endpoint B = %d\n", b );
  printf ( "  The initial seed is %d\n", seed_init );
  printf ( "\n" );

  for ( i = a; i <= b; i++ )
  {
    freq[i-a] = 0;
  }

  i4vec = i4vec_uniform_ab_new ( N, a, b, &seed );

  for ( i = 0; i < N; i++ )
  {
    j = i4vec[i];
    if ( j < a ) 
    {
      printf ( "  Illegal value J = %d\n", j );
    }
    else if ( j <= b )
    {
      freq[j-a] = freq[j-a] + 1;
    }
    else
    {
      printf ( "  Illegal value J = %d\n", j );
    }
  }

  printf ( "\n" );
  printf ( "         I    Frequency\n" );
  printf ( "\n" );
  for ( i = a; i <= b; i++ )
  {
    printf ( "  %8d  %8d\n", i, freq[i-a] );
  }

  free ( i4vec );

  return;
# undef A
# undef B
# undef N
}
/******************************************************************************/

void test11 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST11 tests I8_UNIFORM_AB.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    31 May 2007

  Author:

    John Burkardt
*/
{
  long long int a;
  long long int b;
  long long int i;
  long long int seed;
  long long int seed_init = 123456789LL;

  printf ( "\n" );
  printf ( "TEST11\n" );
  printf ( "  I8_UNIFORM_AB computes pseudorandom values\n" );
  printf ( "  in an interval [A,B].\n" );

  a = 100000LL;
  b = 800000LL;

  seed = seed_init;

  printf ( "\n" );
  printf ( "  The lower endpoint A = %d\n", a );
  printf ( "  The upper endpoint B = %d\n", b );
  printf ( "  The initial seed is    %d\n", seed_init );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    printf ( "  %8d  %24d\n", i, i8_uniform_ab ( a, b, &seed ) );
  }

  return;
}
/******************************************************************************/

void test111 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST111 tests L_UNIFORM.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    01 December 2007

  Author:

    John Burkardt
*/
{
  int i;
  int seed;
  int seed_init = 123456789;

  printf ( "\n" );
  printf ( "TEST111\n" );
  printf ( "  L_UNIFORM computes pseudorandom logical values.\n" );

  seed = seed_init;

  printf ( "\n" );
  printf ( "  The initial seed is %d\n", seed_init );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    printf ( "  %8d  %d\n", i, l_uniform ( &seed ) );
  }

  return;
}
/******************************************************************************/

void test112 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST112 tests LVEC_UNIFORM_NEW.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    01 December 2007

  Author:

    John Burkardt
*/
{
  int i;
  int *lvec;
  int n = 10;
  int seed;
  int seed_init = 123456789;

  printf ( "\n" );
  printf ( "TEST112\n" );
  printf ( "  LVEC_UNIFORM_NEW computes a vector of\n" );
  printf ( "  pseudorandom logical values.\n" );

  seed = seed_init;

  printf ( "\n" );
  printf ( "  The initial seed is %d\n", seed_init );

  lvec = lvec_uniform_new ( n, &seed );

  printf ( "\n" );

  for ( i = 0; i < n; i++ )
  {
    printf ( "  %8d  %d\n", i + 1, lvec[i] );
  }

  free ( lvec );

  return;
}
/******************************************************************************/

void test118 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST118 tests LCRG_ANBN.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    22 April 2008

  Author:

    John Burkardt
*/
{
  int a;
  int an;
  int b;
  int bn;
  int c;
  int n;

  printf ( "\n" );
  printf ( "TEST118\n" );
  printf ( "  LCRG_ANBN determines a linear congruential random\n" );
  printf ( "  number generator equivalent to N steps of a given one.\n" );
/*
  These parameters define the old (1969) IBM 360 random number generator:
*/
  a = 16807;
  b = 0;
  c = 2147483647;

  printf ( "\n" );
  printf ( "  LCRG parameters:\n" );
  printf ( "\n" );
  printf ( "  A = %12d\n", a );
  printf ( "  B = %12d\n", b );
  printf ( "  C = %12d\n", c );
  printf ( "\n" );
  printf ( "             N             A             B\n" );
  printf ( "\n" );

  for ( n = 0; n <= 10; n++ )
  {
    lcrg_anbn ( a, b, c, n, &an, &bn );
    printf ( "  %12d  %12d  %12d\n", n, an, bn );
  }

  return;
}
/******************************************************************************/

void test119 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST119 tests LCRG_ANBN.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    22 April 2008

  Author:

    John Burkardt
*/
{
  int a;
  int an;
  int b;
  int bn;
  int c;
  int j;
  int k;
  int n;
  int u;
  int v;
  int *x;
  int *y;

  printf ( "\n" );
  printf ( "TEST119\n" );
  printf ( "  LCRG_ANBN determines a linear congruential random\n" );
  printf ( "  number generator equivalent to N steps of a given one.\n" );
/*
  These parameters define the old (1969) IBM 360 random number generator:
*/
  a = 16807;
  b = 0;
  c = 2147483647;

  printf ( "\n" );
  printf ( "  LCRG parameters:\n" );
  printf ( "\n" );
  printf ( "  A  = %12d\n", a );
  printf ( "  B  = %12d\n", b );
  printf ( "  C  = %12d\n", c );
  printf ( "\n" );
  printf ( "                           N            In           Out\n" );
  printf ( "\n" );

  k = 0;
  u = 12345;
  printf ( "                %12d                %12d\n", k, u );
  for ( k = 1; k <= 11; k++ )
  {
    v = lcrg_evaluate ( a, b, c, u );
    printf ( "                %12d  %12d  %12d\n", k, u, v );
    u = v;
  }
/*
  Now try to replicate these results using N procesors.
*/
  n = 4;
  x = ( int * ) malloc ( n * sizeof ( int ) );
  y = ( int * ) malloc ( n * sizeof ( int ) );

  lcrg_anbn ( a, b, c, n, &an, &bn );

  printf ( "\n" );
  printf ( "  LCRG parameters:\n" );
  printf ( "\n" );
  printf ( "  AN = %12d\n", an );
  printf ( "  BN = %12d\n", bn );
  printf ( "  C  = %12d\n", c );
  printf ( "\n" );
  printf ( "             J             N            In           Out\n" );
  printf ( "\n" );

  x[0] = 12345;
  for ( j = 1; j < n; j++ )
  {
    x[j] = lcrg_evaluate ( a, b, c, x[j-1] );
  }

  for ( j = 0; j < n; j++ )
  {
    printf ( "  %12d  %12d                %12d\n", j+1, j, x[j] );
  }

  for ( k = n + 1; k <= 12; k = k + n )
  {
    for ( j = 0; j < n; j++ )
    {
      y[j] = lcrg_evaluate ( an, bn, c, x[j] );
      printf ( "  %12d  %12d  %12d  %12d\n", j+1, k+j-1, x[j], y[j] );
      x[j] = y[j];
    }
  }

  free ( x );
  free ( y );

  return;
}
/******************************************************************************/

void test12 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST12 tests LCRG_SEED

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    29 June 2006

  Author:

    John Burkardt
*/
{
  int a;
  int b;
  int c;
  int i;
  int seed;
  int seed_in;
  int seed_lcrg;
  int seed_out;
  int seed_start;
  float u;

  printf ( "\n" );
  printf ( "TEST12\n" );
  printf ( "  LCRG_SEED directly computes the updated value of a\n" );
  printf ( "  seed used by an linear congruential random number\n" );
  printf ( "  generator.\n" );
  printf ( "\n" );
  printf ( "       I          SEED          SEED          SEED    U\n" );
  printf ( "                 Input        Output          LCRG\n" );
  printf ( "\n" );
/*
  These parameters define the old (1969) IBM 360 random number generator:
*/
  a = 16807;
  b = 0;
  c = 2147483647;
/*
  This seed value was used in Pierre L'Ecuyer's article.
*/
  seed_start = 12345;

  seed = seed_start;
/*
  Compute 1000 random numbers "the hard way", that is, sequentially.
  Every now and then, call LCRG_SEED to compute SEED directly.
*/
  for ( i = 1; i <= 1000; i++ )
  {
    seed_in = seed;
    u = r4_uniform_01 ( &seed );
    seed_out = seed;

    if ( i <= 10 || i == 100 || i == 1000 )
    {
      seed_lcrg = lcrg_seed ( a, b, c, i, seed_start );

      printf ( "  %6d  %12d  %12d  %12d  %14f\n", 
        i, seed_in, seed_out, seed_lcrg, u );
    }
  }

  return;
}
/******************************************************************************/

void test13 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST13 tests R4_UNIFORM_AB.

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
  printf ( "TEST13\n" );
  printf ( "  R4_UNIFORM_AB computes pseudorandom values\n" );
  printf ( "  in an interval [A,B].\n" );

  a = 5.0;
  b = 10.0;
  seed = seed_init;

  printf ( "\n" );
  printf ( "  The lower endpoint A = %f\n", a );
  printf ( "  The upper endpoint B = %f\n", b );
  printf ( "  The initial seed is %d\n", seed_init );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    printf ( "  %6d  %14f\n", i, r4_uniform_ab ( a, b, &seed ) );
  }

  return;
}
/******************************************************************************/

void test14 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST14 tests R4_UNIFORM_01.

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
  printf ( "TEST14\n" );
  printf ( "  R4_UNIFORM_01 computes pseudorandom values \n" );
  printf ( "  in the interval [0,1].\n" );

  seed = seed_init;

  printf ( "\n" );
  printf ( "  The initial seed is %d\n", seed_init );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    printf ( "  %6d  %14f\n", i, r4_uniform_01 ( &seed ) );
  }

  return;
}
/******************************************************************************/

void test15 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST15 tests R8_UNIFORM_AB.

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
  printf ( "TEST15\n" );
  printf ( "  R8_UNIFORM_AB computes pseudorandom values\n" );
  printf ( "  in an interval [A,B].\n" );

  a = 5.0;
  b = 10.0;
  seed = seed_init;

  printf ( "\n" );
  printf ( "  The lower endpoint A = %f\n", a );
  printf ( "  The upper endpoint B = %f\n", b );
  printf ( "  The initial seed is %d\n", seed_init );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    printf ( "  %6d  %14f\n", i, r8_uniform_ab ( a, b, &seed ) );
  }

  return;
}
/******************************************************************************/

void test16 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST16 tests R8_UNIFORM_01.

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
  printf ( "TEST16\n" );
  printf ( "  R8_UNIFORM_01 computes pseudorandom values \n" );
  printf ( "  in the interval [0,1].\n" );

  seed = seed_init;

  printf ( "\n" );
  printf ( "  The initial seed is %d\n", seed_init );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    printf ( "  %6d  %14f\n", i, r8_uniform_01 ( &seed ) );
  }

  return;
}
/******************************************************************************/

void test17 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST17 tests R8_UNIFORM_01;

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
  printf ( "TEST17\n" );
  printf ( "  R8_UNIFORM_01 computes a sequence of uniformly distributed\n" );
  printf ( "  pseudorandom numbers.\n" );

  seed = 12345;

  printf ( "\n" );
  printf ( "  Initial SEED = %d\n", seed );

  printf ( "\n" );
  printf ( "  First 10 values:\n" );
  printf ( "\n" );
  printf ( "       I         Input        Output   R8_UNIFORM_01\n" );
  printf ( "                  SEED          SEED\n" );
  printf ( "\n" );

  for ( i = 0; i < 10; i++ )
  {
    seed_in = seed;
    u[i] = r8_uniform_01 ( &seed );
    seed_out = seed;
    printf ( "  %6d  %12d  %12d  %10f\n", i + 1, seed_in, seed_out, u[i] );
  }

  printf ( "\n" );
  printf ( "  Now call R8_UNIFORM_01 %d times.\n", N );

  u_avg = 0.0;
  for ( i = 0; i < N; i++ )
  {
    u[i] = r8_uniform_01 ( &seed );
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
  printf ( "  Expecting       %f\n", 0.5 );

  printf ( "\n" );
  printf ( "  Variance =      %f\n", u_var );
  printf ( "  Expecting       %f\n", 1.0 / 12.0 );

  return;
# undef N
}
/******************************************************************************/

void test18 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST18 tests R8_UNIFORM_01.

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
  int seed_in;
  int seed_out;
  int seed_save;
  double x;

  printf ( "\n" );
  printf ( "TEST18\n" );
  printf ( "  R8_UNIFORM_01 computes a sequence of pseudorandom numbers\n" );
  printf ( "  but all computations depend on the seed value.\n" );
  printf ( "  In this test, we show how a sequence of \"random\"\n" );
  printf ( "  values can be manipulated by accessing the seed.\n" );

  seed = 1066;

  printf ( "\n" );
  printf ( "  Set SEED to %d\n", seed );
  printf ( "\n" );
  printf ( "  Now call R8_UNIFORM_01 10 times, and watch SEED.\n" );
  printf ( "\n" );
  printf ( "       I         Input        Output   R8_UNIFORM_01\n" );
  printf ( "                  SEED          SEED\n" );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    seed_in = seed;

    if ( i == 5 )
    {
      seed_save = seed;
    }
    x = r8_uniform_01 ( &seed );
    seed_out = seed;
    printf ( "  %6d  %12d  %12d  %10f\n", i, seed_in, seed_out, x );
  }

  seed = seed_save;

  printf ( "\n" );
  printf ( "  Reset SEED to its value at step 5, = %d\n", seed );
  printf ( "\n" );
  printf ( "  Now call R8_UNIFORM_01 10 times, and watch how SEED\n" );
  printf ( "  and R8_UNIFORM_01 restart themselves.\n" );
  printf ( "\n" );
  printf ( "       I         Input        Output   R8_UNIFORM_01\n" );
  printf ( "                  SEED          SEED\n" );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    seed_in = seed;
    x = r8_uniform_01 ( &seed );
    seed_out = seed;
    printf ( "  %6d  %12d  %12d  %10f\n", i, seed_in, seed_out, x );
  }

  seed = -12345678;

  printf ( "\n" );
  printf ( "  What happens with an initial negative SEED?\n" );
  printf ( "\n" );
  printf ( "       I         Input        Output   R8_UNIFORM_01\n" );
  printf ( "                  SEED          SEED\n" );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    seed_in = seed;
    x = r8_uniform_01 ( &seed );
    seed_out = seed;
    printf ( "  %6d  %12d  %12d  %10f\n", i, seed_in, seed_out, x );
  }

  return;
}
/******************************************************************************/

void test19 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST19 tests R8_UNIFORM_01 and R8MAT_UNIFORM_01_NEW.

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
  printf ( "TEST19\n" );
  printf ( "  R8_UNIFORM_01 computes pseudorandom values one at a time.\n" );
  printf ( "  R8MAT_UNIFORM_01_NEW computes a matrix of values.\n" );
  printf ( "\n" );
  printf ( "  For the same initial seed, the results should be identical,\n" );
  printf ( "  but R8MAT_UNIFORM_01_NEW might be faster.\n" );
  printf ( "\n" );
  printf ( "  Initial seed is %d\n", seed_init );

  seed = seed_init;
  for ( j = 0; j < N; j++ )
  {
    for ( i = 0; i < M; i++ )
    {
      a[i+j*M] = r8_uniform_01 ( &seed );
    }
  }

  seed = seed_init;
  b = r8mat_uniform_01_new ( M, N, &seed );

  printf ( "\n" );
  printf ( "      I       J      A[I,J]        B[I,J]\n" );
  printf ( "                 (R8_UNIFORM_01)  (R8MAT_UNIFORM_01_NEW)\n" );
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
/******************************************************************************/

void test20 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST20 tests R8_UNIFORM_01 and R8VEC_UNIFORM_01_NEW.

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
  printf ( "TEST20\n" );
  printf ( "  R8_UNIFORM_01 computes pseudeorandom values one at a time.\n" );
  printf ( "  R8VEC_UNIFORM_01_NEW computes a vector of values.\n" );
  printf ( "\n" );
  printf ( "  For the same initial seed, the results should be identical,\n" );
  printf ( "  but R8VEC_UNIFORM_01_NEW might be faster.\n" );
  printf ( "\n" );
  printf ( "  Initial seed is %d\n", seed_init );

  seed = seed_init;
  for ( i = 0; i < N; i++ )
  {
    a[i] = r8_uniform_01 ( &seed );
  }

  seed = seed_init;
  b = r8vec_uniform_01_new ( N, &seed );

  printf ( "\n" );
  printf ( "      I      A[I]          B[I]\n" );
  printf ( "         (R8_UNIFORM_01)  (R8VEC_UNIFORM_01)\n" );
  printf ( "\n" );

  for ( i = 1; i < N; i++ )
  {
    printf ( "  %6d  %12f  %12f\n", i, a[i], b[i] );
  }
  
  free ( b );

  return;
# undef N
}
