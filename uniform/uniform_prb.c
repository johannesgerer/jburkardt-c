# include <complex.h>
# include <stdlib.h>
# include <stdio.h>

# include "uniform.h"

int main ( );
void bvec_uniform_new_test ( );
void c4_uniform_01_test ( );
void c4mat_uniform_01_new_test ( );
void c4vec_uniform_01_new_test ( );
void c8_uniform_01_test ( );
void c8mat_uniform_01_new_test ( );
void c8vec_uniform_01_new_test ( );
void ch_uniform_ab_test ( );
void get_seed_test ( );
void i4_seed_advance_test ( );
void i4_uniform_0i_test ( );
void i4_uniform_ab_test ( );
void i4mat_uniform_ab_new_test ( );
void i4vec_uniform_ab_new_test ( );
void l4_uniform_test ( );
void l4mat_uniform_new_test ( );
void l4vec_uniform_new_test ( );
void lcrg_anbn_test ( );
void lcrg_seed_test ( );
void r4_uniform_01_test ( );
void r4_uniform_ab_test ( );
void r4mat_uniform_ab_new_test ( );
void r4vec_uniform_ab_new_test ( );
void r8_uniform_01_test ( );
void r8_uniform_ab_test ( );
void r8col_uniform_abvec_new_test ( );
void r8mat_uniform_ab_new_test ( );
void r8row_uniform_abvec_new_test ( );
void r8vec_uniform_01_new_test ( );
void r8vec_uniform_ab_new_test ( );
void r8vec_uniform_abvec_new_test ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for UNIFORM_PRB.

  Discussion:

    UNIFORM_PRB tests the UNIFORM library.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    28 December 2014

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "UNIFORM_PRB\n" );
  printf ( "  C version:\n" );
  printf ( "  Test the UNIFORM library.\n" );

  bvec_uniform_new_test ( );

  c4_uniform_01_test ( );
  c4mat_uniform_01_new_test ( );
  c4vec_uniform_01_new_test ( );

  c8_uniform_01_test ( );
  c8mat_uniform_01_new_test ( );
  c8vec_uniform_01_new_test ( );

  ch_uniform_ab_test ( );

  get_seed_test ( );

  i4_seed_advance_test ( );

  i4_uniform_0i_test ( );
  i4_uniform_ab_test ( );
  i4mat_uniform_ab_new_test ( );
  i4vec_uniform_ab_new_test ( );

  l4_uniform_test ( );
  l4mat_uniform_new_test ( );
  l4vec_uniform_new_test ( );

  lcrg_anbn_test ( );
  lcrg_seed_test ( );

  r4_uniform_01_test ( );
  r4_uniform_ab_test ( );
  r4mat_uniform_ab_new_test ( );
  r4vec_uniform_ab_new_test ( );

  r8_uniform_01_test ( );
  r8_uniform_ab_test ( );
  r8mat_uniform_ab_new_test ( );
  r8vec_uniform_01_new_test ( );
  r8vec_uniform_ab_new_test ( );

  r8col_uniform_abvec_new_test ( );
  r8row_uniform_abvec_new_test ( );
  r8vec_uniform_abvec_new_test ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "UNIFORM_PRB\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void bvec_uniform_new_test ( )

/******************************************************************************/
/*
  Purpose:

    BVEC_UNIFORM_NEW_TEST tests BVEC_UNIFORM_NEW.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    26 December 2014

  Author:

    John Burkardt
*/
{
  int *b;
  int i;
  int n = 10;
  int seed;

  seed = 123456789;

  printf ( "\n" );
  printf ( "BVEC_UNIFORM_NEW_TEST\n" );
  printf ( "  BVEC_UNIFORM_NEW computes a binary vector.\n" );
  printf ( "\n" );
  printf ( "  The initial seed is %d\n", seed );
  printf ( "\n" );

  for ( i = 0; i < 10; i++ )
  {
    b = bvec_uniform_new ( n, &seed );
    bvec_print ( n, b, "" );
    free ( b );
  }

  return;
}
/******************************************************************************/

void c4_uniform_01_test ( )

/******************************************************************************/
/*
  Purpose:

    C4_UNIFORM_01_TEST tests C4_UNIFORM_01.

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
  float complex value;

  seed = 123456789;

  printf ( "\n" );
  printf ( "C4_UNIFORM_01_TEST\n" );
  printf ( "  C4_UNIFORM_01 computes pseudorandom complex values\n" );
  printf ( "  uniformly distributed in the unit circle.\n" );
  printf ( "\n" );
  printf ( "  The initial seed is %d\n", seed );

  printf ( "\n" );
  for ( i = 1; i <= 10; i++ )
  {
    value = c4_uniform_01 ( &seed );

    printf ( "  %6d  %12f  %12f\n", i, creal ( value ), cimag ( value ) );
  }

  return;
}
/******************************************************************************/

void c4mat_uniform_01_new_test ( )

/******************************************************************************/
/*
  Purpose:

    C4MAT_UNIFORM_01_NEW_TEST tests C4MAT_UNIFORM_01_NEW.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    14 December 2014

  Author:

    John Burkardt
*/
{
  float complex *c;
  int m;
  int n;
  int seed;

  m = 5;
  n = 2;
  seed = 123456789;

  printf ( "\n" );
  printf ( "C4MAT_UNIFORM_01_NEW_TEST\n" );
  printf ( "  C4MAT_UNIFORM_01_NEW computes pseudorandom complex values\n" );
  printf ( "  uniformly distributed in the unit circle.\n" );
  printf ( "\n" );
  printf ( "  The initial seed is %d\n", seed );

  c = c4mat_uniform_01_new ( m, n, &seed );

  c4mat_print ( m, n, c, "  Uniform C4MAT:" );

  free ( c );

  return;
}
/******************************************************************************/

void c4vec_uniform_01_new_test ( )

/******************************************************************************/
/*
  Purpose:

    C4VEC_UNIFORM_01_NEW_TEST tests C4VEC_UNIFORM_01_NEW.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 January 2007

  Author:

    John Burkardt
*/
{
  float complex *c;
  int n;
  int seed;

  n = 10;
  seed = 123456789;

  printf ( "\n" );
  printf ( "C4VEC_UNIFORM_01_NEW_TEST\n" );
  printf ( "  C4VEC_UNIFORM_01_NEW computes pseudorandom complex values\n" );
  printf ( "  uniformly distributed in the unit circle.\n" );
  printf ( "\n" );
  printf ( "  The initial seed is %d\n", seed );

  c = c4vec_uniform_01_new ( n, &seed );

  c4vec_print ( n, c, "  Uniform C4VEC:" );

  free ( c );

  return;
}
/******************************************************************************/

void c8_uniform_01_test ( )

/******************************************************************************/
/*
  Purpose:

    C8_UNIFORM_01_TEST tests C8_UNIFORM_01.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 December 2014

  Author:

    John Burkardt
*/
{
  int i;
  int seed;
  double complex value;

  seed = 123456789;

  printf ( "\n" );
  printf ( "C8_UNIFORM_01_TEST\n" );
  printf ( "  C8_UNIFORM_01 computes pseudorandom complex values\n" );
  printf ( "  uniformly distributed in the unit circle.\n" );
  printf ( "\n" );
  printf ( "  The initial seed is %d\n", seed );

  printf ( "\n" );
  for ( i = 1; i <= 10; i++ )
  {
    value = c8_uniform_01 ( &seed );

    printf ( "  %6d  %12f  %12f\n", i, creal ( value ), cimag ( value ) );
  }

  return;
}
/******************************************************************************/

void c8mat_uniform_01_new_test ( )

/******************************************************************************/
/*
  Purpose:

    C8MAT_UNIFORM_01_NEW_TEST tests C8MAT_UNIFORM_01_NEW.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 December 2014

  Author:

    John Burkardt
*/
{
  double complex *c;
  int m;
  int n;
  int seed;

  m = 5;
  n = 2;
  seed = 123456789;

  printf ( "\n" );
  printf ( "C8MAT_UNIFORM_01_NEW_TEST\n" );
  printf ( "  C8MAT_UNIFORM_01_NEW computes pseudorandom complex values\n" );
  printf ( "  uniformly distributed in the unit circle.\n" );
  printf ( "\n" );
  printf ( "  The initial seed is %d\n", seed );

  c = c8mat_uniform_01_new ( m, n, &seed );

  c8mat_print ( m, n, c, "  Uniform C4MAT:" );

  free ( c );

  return;
}
/******************************************************************************/

void c8vec_uniform_01_new_test ( )

/******************************************************************************/
/*
  Purpose:

    C8VEC_UNIFORM_01_NEW_TEST tests C8VEC_UNIFORM_01_NEW.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    16 December 2014

  Author:

    John Burkardt
*/
{
  double complex *c;
  int n;
  int seed;

  n = 10;
  seed = 123456789;

  printf ( "\n" );
  printf ( "C8VEC_UNIFORM_01_NEW_TEST\n" );
  printf ( "  C8VEC_UNIFORM_01_NEW computes pseudorandom complex values\n" );
  printf ( "  uniformly distributed in the unit circle.\n" );
  printf ( "\n" );
  printf ( "  The initial seed is %d\n", seed );

  c = c8vec_uniform_01_new ( n, &seed );

  c8vec_print ( n, c, "  Uniform C4VEC:" );

  free ( c );

  return;
}
/******************************************************************************/

void ch_uniform_ab_test ( )

/******************************************************************************/
/*
  Purpose:

    CH_UNIFORM_AB_TEST tests CH_UNIFORM_AB.

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

  clo = 'A';
  chi = 'J';
  seed = 123456789;

  printf ( "\n" );
  printf ( "CH_UNIFORM_AB_TEST\n" );
  printf ( "  CH_UNIFORM_AB computes pseudorandom characters\n" );
  printf ( "  in an interval [CLO,CHI].\n" );
  printf ( "\n" );
  printf ( "  The lower endpoint CLO = '%c'.\n", clo );
  printf ( "  The upper endpoint CHI = '%c'.\n", chi );
  printf ( "  The initial seed is %d\n", seed );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    printf ( "  %6d  '%c'\n", i, ch_uniform_ab ( clo, chi, &seed ) );
  }

  return;
}
/******************************************************************************/

void get_seed_test ( )

/******************************************************************************/
/*
  Purpose:

    GET_SEED_TEST tests GET_SEED.

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

  oops = 0;
  seed = 12345678;
  seed_old = seed;

  printf ( "\n" );
  printf ( "GET_SEED_TEST\n" );
  printf ( "  GET_SEED picks an initial seed value for R8_UNIFORM_01.\n" );
  printf ( "  The value chosen should vary over time, because\n" );
  printf ( "  the seed is based on reading the clock.\n" );
  printf ( "\n" );
  printf ( "  This is just the \"calendar\" clock, which does\n" );
  printf ( "  not change very fast, so calling GET_SEED several\n" );
  printf ( "  times in a row may result in the same value.\n" );
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

void i4_seed_advance_test ( )

/******************************************************************************/
/*
  Purpose:

    I4_SEED_ADVANCE_TEST tests I4_SEED_ADVANCE.

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

  seed_new = 12345;

  printf ( "\n" );
  printf ( "I4_SEED_ADVANCE_TEST\n" );
  printf ( "  I4_SEED_ADVANCE advances the seed.\n" );
  printf ( "\n" );
  printf ( "  Step        SEED input       SEED output\n" );
  printf ( "\n" );

  for ( step = 1; step <= 10; step++ )
  {
    seed = seed_new;
    seed_new = i4_seed_advance ( seed );

    printf ( "  %4d  %16d  %16d\n", step, seed, seed_new );
  }

  return;
}
/******************************************************************************/

void i4_uniform_0i_test ( )

/******************************************************************************/
/*
  Purpose:

    I4_UNIFORM_0I_TEST tests I4_UNIFORM_0I

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

  seed = 123456789;

  printf ( "\n" );
  printf ( "I4_UNIFORM_0I_TEST\n" );
  printf ( "  I4_UNIFORM_0I samples a uniform random\n" );
  printf ( "  integer distribution in [0,2^31-1].\n" );
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

void i4_uniform_ab_test ( )

/******************************************************************************/
/*
  Purpose:

    I4_UNIFORM_AB_TEST tests I4_UNIFORM_AB.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    27 October 2014

  Author:

    John Burkardt
*/
{
  int a = -100;
  int b = 200;
  int i;
  int j;
  int seed = 123456789;

  printf ( "\n" );
  printf ( "I4_UNIFORM_AB_TEST\n" );
  printf ( "  I4_UNIFORM_AB computes pseudorandom values\n" );
  printf ( "  in an interval [A,B].\n" );
  printf ( "\n" );
  printf ( "  The lower endpoint A = %d\n", a );
  printf ( "  The upper endpoint B = %d\n", b );
  printf ( "  The initial seed is %d\n", seed );
  printf ( "\n" );

  for ( i = 1; i <= 20; i++ )
  {
    j = i4_uniform_ab ( a, b, &seed );
    printf ( "  %8d  %d\n", i, j );
  }

  return;
}
/******************************************************************************/

void i4mat_uniform_ab_new_test ( )

/******************************************************************************/
/*
  Purpose:

    I4MAT_UNIFORM_AB_NEW_TEST tests I4MAT_UNIFORM_AB_NEW.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    18 December 2014

  Author:

    John Burkardt
*/
{
  int a = -100;
  int b = 200;
  int m = 5;
  int n = 4;
  int seed = 123456789;
  int *v;

  printf ( "\n" );
  printf ( "I4MAT_UNIFORM_AB_NEW_TEST\n" );
  printf ( "  I4MAT_UNIFORM_AB_NEW computes pseudorandom values\n" );
  printf ( "  in an interval [A,B].\n" );
  printf ( "\n" );
  printf ( "  The lower endpoint A = %d\n", a );
  printf ( "  The upper endpoint B = %d\n", b );
  printf ( "  The initial seed is %d\n", seed );

  v = i4mat_uniform_ab_new ( m, n, a, b, &seed );

  i4mat_print ( m, n, v, "  Uniform I4MAT:" );

  free ( v );

  return;
}
/******************************************************************************/

void i4vec_uniform_ab_new_test ( )

/******************************************************************************/
/*
  Purpose:

    I4VEC_UNIFORM_AB_NEW_TEST tests I4VEC_UNIFORM_AB_NEW.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    27 October 2014

  Author:

    John Burkardt
*/
{
  int a = -100;
  int b = 200;
  int n = 20;
  int seed = 123456789;
  int *v;

  printf ( "\n" );
  printf ( "I4VEC_UNIFORM_AB_NEW_TEST\n" );
  printf ( "  I4VEC_UNIFORM_AB_NEW computes pseudorandom values\n" );
  printf ( "  in an interval [A,B].\n" );
  printf ( "\n" );
  printf ( "  The lower endpoint A = %d\n", a );
  printf ( "  The upper endpoint B = %d\n", b );
  printf ( "  The initial seed is %d\n", seed );

  v = i4vec_uniform_ab_new ( n, a, b, &seed );

  i4vec_print ( n, v, "  Uniform I4VEC:" );

  free ( v );

  return;
}
/******************************************************************************/

void l4_uniform_test ( )

/******************************************************************************/
/*
  Purpose:

    L4_UNIFORM_TEST tests L4_UNIFORM.

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

  seed = 123456789;

  printf ( "\n" );
  printf ( "L4_UNIFORM_TEST\n" );
  printf ( "  L4_UNIFORM computes pseudorandom logical values.\n" );
  printf ( "\n" );
  printf ( "  The initial seed is %d\n", seed );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    printf ( "  %8d  %d\n", i, l4_uniform ( &seed ) );
  }

  return;
}
/******************************************************************************/

void l4mat_uniform_new_test ( )

/******************************************************************************/
/*
  Purpose:

    L4MAT_UNIFORM_NEW_TEST tests L4MAT_UNIFORM_NEW.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    23 December 2014

  Author:

    John Burkardt
*/
{
  int *l;
  int m = 5;
  int n = 4;
  int seed;

  seed = 123456789;

  printf ( "\n" );
  printf ( "L4MAT_UNIFORM_NEW_TEST\n" );
  printf ( "  L4MAT_UNIFORM_NEW computes a vector of\n" );
  printf ( "  pseudorandom logical values.\n" );
  printf ( "\n" );
  printf ( "  The initial seed is %d\n", seed );

  l = l4mat_uniform_new ( m, n, &seed );

  l4mat_print ( m, n, l, "  Uniform L4MAT:" );

  free ( l );

  return;
}
/******************************************************************************/

void l4vec_uniform_new_test ( )

/******************************************************************************/
/*
  Purpose:

    L4VEC_UNIFORM_NEW_TEST tests L4VEC_UNIFORM_NEW.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    20 December 2014

  Author:

    John Burkardt
*/
{
  int *l;
  int n = 10;
  int seed;

  seed = 123456789;

  printf ( "\n" );
  printf ( "L4VEC_UNIFORM_NEW_TEST\n" );
  printf ( "  L4VEC_UNIFORM_NEW computes a vector of\n" );
  printf ( "  pseudorandom logical values.\n" );
  printf ( "\n" );
  printf ( "  The initial seed is %d\n", seed );

  l = l4vec_uniform_new ( n, &seed );

  l4vec_print ( n, l, "  Uniform L4VEC:" );

  free ( l );

  return;
}
/******************************************************************************/

void lcrg_anbn_test ( )

/******************************************************************************/
/*
  Purpose:

    LCRG_ANBN_TEST tests LCRG_ANBN.

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
/*
  These parameters define the old (1969) IBM 360 random number generator:
*/
  a = 16807;
  b = 0;
  c = 2147483647;

  printf ( "\n" );
  printf ( "LCRG_ANBN_TEST\n" );
  printf ( "  LCRG_ANBN determines a linear congruential random\n" );
  printf ( "  number generator equivalent to N steps of a given one.\n" );

  printf ( "\n" );
  printf ( "  LCRG parameters:\n" );
  printf ( "\n" );
  printf ( "  A  = %12d\n", a );
  printf ( "  B  = %12d\n", b );
  printf ( "  C  = %12d\n", c );

  printf ( "\n" );
  printf ( "             N             A             B\n" );
  printf ( "\n" );

  for ( n = 0; n <= 10; n++ )
  {
    lcrg_anbn ( a, b, c, n, &an, &bn );
    printf ( "  %12d  %12d  %12d\n", n, an, bn );
  }

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

void lcrg_seed_test ( )

/******************************************************************************/
/*
  Purpose:

    LCRG_SEED_TEST tests LCRG_SEED

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
/*
  These parameters define the old (1969) IBM 360 random number generator:
*/
  a = 16807;
  b = 0;
  c = 2147483647;

  printf ( "\n" );
  printf ( "LCRG_SEED_TEST\n" );
  printf ( "  LCRG_SEED directly computes the updated value of a\n" );
  printf ( "  seed used by an linear congruential random number\n" );
  printf ( "  generator.\n" );
  printf ( "\n" );
  printf ( "       I          SEED          SEED          SEED    U\n" );
  printf ( "                 Input        Output          LCRG\n" );
  printf ( "\n" );
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

void r4_uniform_01_test ( )

/******************************************************************************/
/*
  Purpose:

    R4_UNIFORM_01_TEST tests R4_UNIFORM_01.

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

  seed = 123456789;

  printf ( "\n" );
  printf ( "R4_UNIFORM_01_TEST\n" );
  printf ( "  R4_UNIFORM_01 computes pseudorandom values \n" );
  printf ( "  in the interval [0,1].\n" );
  printf ( "\n" );
  printf ( "  The initial seed is %d\n", seed );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    printf ( "  %6d  %14f\n", i, r4_uniform_01 ( &seed ) );
  }

  return;
}
/******************************************************************************/

void r4_uniform_ab_test ( )

/******************************************************************************/
/*
  Purpose:

    R4_UNIFORM_AB_TEST tests R4_UNIFORM_AB.

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

  a = 5.0;
  b = 10.0;
  seed = 123456789;

  printf ( "\n" );
  printf ( "R4_UNIFORM_AB_TEST\n" );
  printf ( "  R4_UNIFORM_AB computes pseudorandom values\n" );
  printf ( "  in an interval [A,B].\n" );
  printf ( "\n" );
  printf ( "  The lower endpoint A = %f\n", a );
  printf ( "  The upper endpoint B = %f\n", b );
  printf ( "  The initial seed is %d\n", seed );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    printf ( "  %6d  %14f\n", i, r4_uniform_ab ( a, b, &seed ) );
  }

  return;
}
/******************************************************************************/

void r4mat_uniform_ab_new_test ( )

/******************************************************************************/
/*
  Purpose:

    R4MAT_UNIFORM_AB_NEW_TEST tests R4MAT_UNIFORM_AB_NEW.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    24 December 2014

  Author:

    John Burkardt
*/
{
  float a = -5.0;
  float b = +10.0;
  int m = 5;
  int n = 4;
  int seed = 123456789;
  float *v;

  printf ( "\n" );
  printf ( "R4MAT_UNIFORM_AB_NEW_TEST\n" );
  printf ( "  R4MAT_UNIFORM_AB_NEW computes pseudorandom values\n" );
  printf ( "  in an interval [A,B].\n" );
  printf ( "\n" );
  printf ( "  The lower endpoint A = %g\n", a );
  printf ( "  The upper endpoint B = %g\n", b );
  printf ( "  The initial seed is %d\n", seed );

  v = r4mat_uniform_ab_new ( m, n, a, b, &seed );

  r4mat_print ( m, n, v, "  Uniform R4MAT:" );

  free ( v );

  return;
}
/******************************************************************************/

void r4vec_uniform_ab_new_test ( )

/******************************************************************************/
/*
  Purpose:

    R4VEC_UNIFORM_AB_NEW_TEST tests R4VEC_UNIFORM_AB_NEW.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    24 December 2014

  Author:

    John Burkardt
*/
{
  float a = -5.0;
  float b = +10.0;
  int n = 20;
  int seed = 123456789;
  float *v;

  printf ( "\n" );
  printf ( "R4VEC_UNIFORM_AB_NEW_TEST\n" );
  printf ( "  R4VEC_UNIFORM_AB_NEW computes pseudorandom values\n" );
  printf ( "  in an interval [A,B].\n" );
  printf ( "\n" );
  printf ( "  The lower endpoint A = %g\n", a );
  printf ( "  The upper endpoint B = %g\n", b );
  printf ( "  The initial seed is %d\n", seed );

  v = r4vec_uniform_ab_new ( n, a, b, &seed );

  r4vec_print ( n, v, "  Uniform R4VEC:" );

  free ( v );

  return;
}
/******************************************************************************/

void r8_uniform_01_test ( )

/******************************************************************************/
/*
  Purpose:

    R8_UNIFORM_01_TEST tests R8_UNIFORM_01.

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

  seed = 123456789;

  printf ( "\n" );
  printf ( "R8_UNIFORM_01_TEST\n" );
  printf ( "  R8_UNIFORM_01 computes pseudorandom values \n" );
  printf ( "  in the interval [0,1].\n" );
  printf ( "\n" );
  printf ( "  The initial seed is %d\n", seed );

  printf ( "\n" );
  for ( i = 1; i <= 10; i++ )
  {
    printf ( "  %6d  %14f\n", i, r8_uniform_01 ( &seed ) );
  }

  return;
}
/******************************************************************************/

void r8_uniform_ab_test ( )

/******************************************************************************/
/*
  Purpose:

    R8_UNIFORM_AB_TEST tests R8_UNIFORM_AB.

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

  a = 5.0;
  b = 10.0;
  seed = 123456789;

  printf ( "\n" );
  printf ( "R8_UNIFORM_AB_TEST\n" );
  printf ( "  R8_UNIFORM_AB computes pseudorandom values\n" );
  printf ( "  in an interval [A,B].\n" );
  printf ( "\n" );
  printf ( "  The lower endpoint A = %f\n", a );
  printf ( "  The upper endpoint B = %f\n", b );
  printf ( "  The initial seed is %d\n", seed );
  printf ( "\n" );

  for ( i = 1; i <= 10; i++ )
  {
    printf ( "  %6d  %14f\n", i, r8_uniform_ab ( a, b, &seed ) );
  }

  return;
}
/******************************************************************************/

void r8mat_uniform_ab_new_test ( )

/******************************************************************************/
/*
  Purpose:

    R8MAT_UNIFORM_AB_NEW_TEST tests R8MAT_UNIFORM_AB_NEW.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    25 December 2014

  Author:

    John Burkardt
*/
{
  double a = -5.0;
  double b = +10.0;
  int m = 5;
  int n = 4;
  int seed = 123456789;
  double *v;

  printf ( "\n" );
  printf ( "R8MAT_UNIFORM_AB_NEW_TEST\n" );
  printf ( "  R8MAT_UNIFORM_AB_NEW computes pseudorandom values\n" );
  printf ( "  in an interval [A,B].\n" );
  printf ( "\n" );
  printf ( "  The lower endpoint A = %g\n", a );
  printf ( "  The upper endpoint B = %g\n", b );
  printf ( "  The initial seed is %d\n", seed );

  v = r8mat_uniform_ab_new ( m, n, a, b, &seed );

  r8mat_print ( m, n, v, "  Uniform R8MAT:" );

  free ( v );

  return;
}
/******************************************************************************/

void r8vec_uniform_01_new_test ( )

/******************************************************************************/
/*
  Purpose:

    R8VEC_UNIFORM_01_NEW_TEST tests R8VEC_UNIFORM_01_NEW.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    29 October 2014

  Author:

    John Burkardt
*/
{
  int n = 10;
  int seed;
  double *v;

  seed = 123456789;

  printf ( "\n" );
  printf ( "R8VEC_UNIFORM_01_NEW_TEST\n" );
  printf ( "  R8VEC_UNIFORM_01_NEW computes a random R8VEC.\n" );
  printf ( "\n" );
  printf ( "  Initial seed is %d\n", seed );

  v = r8vec_uniform_01_new ( n, &seed );

  r8vec_print ( n, v, "  Uniform R8VEC:" );
  
  free ( v );

  return;
}
/******************************************************************************/

void r8vec_uniform_ab_new_test ( )

/******************************************************************************/
/*
  Purpose:

    R8VEC_UNIFORM_AB_NEW_TEST tests R8VEC_UNIFORM_AB_NEW.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    29 October 2014

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  int n = 10;
  int seed;
  double *v;

  a = -1.0;
  b = 5.0;
  seed = 123456789;

  printf ( "\n" );
  printf ( "R8VEC_UNIFORM_AB_NEW_TEST\n" );
  printf ( "  R8VEC_UNIFORM_AB_NEW computes a random R8VEC.\n" );
  printf ( "\n" );
  printf ( "  %g <= X <= %g\n", a, b );
  printf ( "  Initial seed is %d\n", seed );

  v = r8vec_uniform_ab_new ( n, a, b, &seed );

  r8vec_print ( n, v, "  Uniform R8VEC:" );
  
  free ( v );

  return;
}
/******************************************************************************/

void r8col_uniform_abvec_new_test ( )

/******************************************************************************/
/*
  Purpose:

    R8COL_UNIFORM_ABVEC_NEW_TEST tests R8COL_UNIFORM_ABVEC_NEW.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    28 December 2014

  Author:

    John Burkardt
*/
{
  double a[5] = { 0.0, 0.20, 10.0, 52.0, -1.0 };
  double b[5] = { 1.0, 0.25, 20.0, 54.0, +1.0 };
  int i;
  int j;
  int m = 5;
  int n = 4;
  int seed;
  double *v;

  seed = 123456789;

  printf ( "\n" );
  printf ( "R8COL_UNIFORM_ABVEC_NEW_TEST\n" );
  printf ( "  R8COL_UNIFORM_ABVEC_NEW computes a random R8COL.\n" );
  printf ( "\n" );
  printf ( "  Initial seed is %d\n", seed );

  v = r8col_uniform_abvec_new ( m, n, a, b, &seed );

  printf ( "\n" );
  for ( i = 0; i < m; i++ )
  {
    printf ( "  %4d  %8.4f:  ", i, a[i] );
    for ( j = 0; j < n; j++ )
    {
      printf ( "  %8.4f", v[i+j*m] );
    }
    printf ( "    :%8.4f\n", b[i] );
  }
  
  free ( v );

  return;
}
/******************************************************************************/

void r8row_uniform_abvec_new_test ( )

/******************************************************************************/
/*
  Purpose:

    R8ROW_UNIFORM_ABVEC_NEW_TEST tests R8ROWL_UNIFORM_ABVEC_NEW.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    28 December 2014

  Author:

    John Burkardt
*/
{
  double a[5] = { 0.0, 0.20, 10.0, 52.0, -1.0 };
  double b[5] = { 1.0, 0.25, 20.0, 54.0, +1.0 };
  int i;
  int j;
  int m = 4;
  int n = 5;
  int seed;
  double *v;

  seed = 123456789;

  printf ( "\n" );
  printf ( "R8ROW_UNIFORM_ABVEC_NEW_TEST\n" );
  printf ( "  R8ROW_UNIFORM_ABVEC_NEW computes a random R8ROW.\n" );
  printf ( "\n" );
  printf ( "  Initial seed is %d\n", seed );
  printf ( "\n" );

  v = r8row_uniform_abvec_new ( m, n, a, b, &seed );

  for ( j = 0; j < n; j++ )
  {
    printf ( "  %8.4f", b[j] );
  }
  printf ( "\n" );
  printf ( "\n" );

  for ( i = 0; i < m; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      printf ( "  %8.4f", v[i+j*m] );
    }
    printf ( "\n" );
  }

  printf ( "\n" );
  for ( j = 0; j < n; j++ )
  {
    printf ( "  %8.4f", a[j] );
  }
  printf ( "\n" );

  free ( v );

  return;
}
/******************************************************************************/

void r8vec_uniform_abvec_new_test ( )

/******************************************************************************/
/*
  Purpose:

    R8VEC_UNIFORM_ABVEC_NEW_TEST tests R8VEC_UNIFORM_ABVEC_NEW.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    28 December 2014

  Author:

    John Burkardt
*/
{
  double a[5] = { 0.0, 0.20, 10.0, 52.0, -1.0 };
  double b[5] = { 1.0, 0.25, 20.0, 54.0, +1.0 };
  int i;
  int n = 5;
  int seed;
  double *v;

  seed = 123456789;

  printf ( "\n" );
  printf ( "R8VEC_UNIFORM_ABVEC_NEW_TEST\n" );
  printf ( "  R8VEC_UNIFORM_ABVEC_NEW computes a random R8VEC.\n" );
  printf ( "\n" );
  printf ( "  Initial seed is %d\n", seed );

  v = r8vec_uniform_abvec_new ( n, a, b, &seed );

  printf ( "\n" );
  printf ( "   I         A         X         B\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %4d  %8.4f  %8.4f  %8.4f\n", i, a[i], v[i], b[i] );
  }
  
  free ( v );

  return;
}
