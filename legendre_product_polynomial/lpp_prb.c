# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <string.h>
# include <time.h>

# include "lpp.h"

int main ( );

void comp_enum_test ( );
void comp_next_grlex_test ( );
void comp_random_grlex_test ( );
void comp_rank_grlex_test ( );
void comp_unrank_grlex_test ( );

void i4_choose_test ( );
void i4_uniform_ab_test ( );

void i4vec_permute_test ( );
void i4vec_print_test ( );
void i4vec_sort_heap_index_a_test ( );
void i4vec_sum_test ( );
void i4vec_uniform_ab_new_test ( );

void lp_coefficients_test ( );
void lp_value_test ( );
void lp_values_test ( );

void lpp_to_polynomial_test ( );
void lpp_value_test ( );

void mono_next_grlex_test ( );
void mono_print_test ( );
void mono_rank_grlex_test ( );
void mono_unrank_grlex_test ( );
void mono_upto_enum_test ( );
void mono_upto_next_grlex_test ( );
void mono_upto_random_test ( );
void mono_value_test ( );

void perm_uniform_test ( );

void polynomial_compress_test ( );
void polynomial_print_test ( );
void polynomial_sort_test ( );
void polynomial_value_test ( );

void r8mat_print_test ( );
void r8mat_print_some_test ( );
void r8mat_uniform_ab_new_test ( );

void r8vec_permute_test ( );
void r8vec_print_test ( );
void r8vec_uniform_ab_new_test ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for LPP_PRB.

  Discussion:

    LPP_PRB tests the LEGENDRE_PRODUCT_POLYNOMIAL library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    03 November 2014

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "LPP_PRB:\n" );
  printf ( "  C version\n" );
  printf ( "  Test the LEGENDRE_PRODUCT_POLYNOMIAL library.\n" );

  i4_choose_test ( );
  i4_uniform_ab_test ( );

  i4vec_permute_test ( );
  i4vec_print_test ( );
  i4vec_sort_heap_index_a_test ( );
  i4vec_sum_test ( );
  i4vec_uniform_ab_new_test ( );

  r8vec_permute_test ( );
  r8vec_print_test ( );
  r8vec_uniform_ab_new_test ( );

  r8mat_print_test ( );
  r8mat_print_some_test ( );
  r8mat_uniform_ab_new_test ( );

  perm_uniform_test ( );

  comp_enum_test ( );
  comp_next_grlex_test ( );
  comp_random_grlex_test ( );
  comp_rank_grlex_test ( );
  comp_unrank_grlex_test ( );

  mono_next_grlex_test ( );
  mono_print_test ( );
  mono_rank_grlex_test ( );
  mono_unrank_grlex_test ( );
  mono_upto_enum_test ( );
  mono_upto_next_grlex_test ( );
  mono_upto_random_test ( );
  mono_value_test ( );

  polynomial_compress_test ( );
  polynomial_print_test ( );
  polynomial_sort_test ( );
  polynomial_value_test ( );

  lp_coefficients_test ( );
  lp_value_test ( );
  lp_values_test ( );

  lpp_to_polynomial_test ( );
  lpp_value_test ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "LPP_PRB:\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void comp_enum_test ( )

/******************************************************************************/
/*
  Purpose:

    COMP_ENUM_TEST tests COMP_ENUM;

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    29 October 2014

  Author:

    John Burkardt
*/
{
  int num;
  int k;
  int n;

  printf ( "\n" );
  printf ( "COMP_ENUM_TEST\n" );
  printf ( "  COMP_ENUM counts compositions;\n" );
  printf ( "\n" );
  for ( n = 0; n <= 10; n++ )
  {
    for ( k = 1; k <= 10; k++ )
    {
      num = comp_enum ( n, k );
      printf ( "  %6d", num );
    }
    printf ( "\n" );
  }

  return;
}
/******************************************************************************/

void comp_next_grlex_test ( )

/******************************************************************************/
/*
  Purpose:

    COMP_NEXT_GRLEX_TEST tests COMP_NEXT_GRLEX.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    29 October 2014

  Author:

    John Burkardt
*/
{
  int j;
  int kc = 3;
  int nc;
  int rank;
  int xc[3];

  printf ( "\n" );
  printf ( "COMP_NEXT_GRLEX_TEST\n" );
  printf ( "  A COMP is a composition of an integer N into K parts.\n" );
  printf ( "  Each part is nonnegative.  The order matters.\n" );
  printf ( "  COMP_NEXT_GRLEX determines the next COMP in\n" );
  printf ( "  graded lexicographic (grlex) order.\n" );
  
  printf ( "\n" );
  printf ( "  Rank:     NC       COMP\n" );
  printf ( "  ----:     --   ------------\n" );

  for ( rank = 1; rank <= 71; rank++ )
  {
    if ( rank == 1 )
    {
      for ( j = 0; j < kc; j++ )
      {
        xc[j] = 0;
      }
    }
    else
    {
      comp_next_grlex ( kc, xc );
    }

    nc = i4vec_sum ( kc, xc );

    printf ( "   %3d: ", rank );
    printf ( "    %2d = ", nc );
    for ( j = 0; j < kc - 1; j++ )
    {
      printf ( "%2d + ", xc[j] );
    }
    printf ( "%2d\n", xc[kc-1] );
/*
  When XC(1) == NC, we have completed the compositions associated with
  a particular integer, and are about to advance to the next integer.
*/
    if ( xc[0] == nc )
    {
      printf ( "  ----:     --   ------------\n" );
    }
  }

  return;
}
/******************************************************************************/

void comp_random_grlex_test ( )

/******************************************************************************/
/*
  Purpose:

    COMP_RANDOM_GRLEX_TEST tests COMP_RANDOM_GRLEX.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    29 October 2014

  Author:

    John Burkardt
*/
{
  int j;
  int kc;
  int nc;
  int rank;
  int rank1;
  int rank2;
  int seed;
  int test;
  int *xc;

  printf ( "\n" );
  printf ( "COMP_RANDOM_GRLEX_TEST\n" );
  printf ( "  A COMP is a composition of an integer N into K parts.\n" );
  printf ( "  Each part is nonnegative.  The order matters.\n" );
  printf ( "  COMP_RANDOM_GRLEX selects a random COMP in\n" );
  printf ( "  graded lexicographic (grlex) order between indices RANK1 and RANK2.\n" );
  printf ( "\n" );

  kc = 3;
  rank1 = 20;
  rank2 = 60;
  seed = 123456789;

  for ( test = 1; test <= 5; test++ )
  {
    xc = comp_random_grlex ( kc, rank1, rank2, &seed, &rank );
    nc = i4vec_sum ( kc, xc );

    printf ( "   %3d: ", rank );
    printf ( "    %2d = ", nc );
    for ( j = 0; j < kc - 1; j++ )
    {
      printf ( "%2d + ", xc[j] );
    }
    printf ( "%2d\n", xc[kc-1] );
    free ( xc );
  }

  return;
}
/******************************************************************************/

void comp_rank_grlex_test ( )

/******************************************************************************/
/*
  Purpose:

    COMP_RANK_GRLEX_TEST tests COMP_RANK_GRLEX.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 October 2014

  Author:

    John Burkardt
*/
{
  int kc;
  int nc;
  int rank1;
  int rank2;
  int rank3;
  int rank4;
  int seed;
  int test;
  int *xc;

  printf ( "\n" );
  printf ( "COMP_RANK_GRLEX_TEST\n" );
  printf ( "  A COMP is a composition of an integer N into K parts.\n" );
  printf ( "  Each part is nonnegative.  The order matters.\n" );
  printf ( "  COMP_RANK_GRLEX determines the rank of a COMP\n" );
  printf ( "  from its parts.\n" );
  printf ( "\n" );
  printf ( "        Actual  Inferred\n" );
  printf ( "  Test    Rank      Rank\n" );
  printf ( "\n" );

  kc = 3;
  rank1 = 20;
  rank2 = 60;
  seed = 123456789;

  for ( test = 1; test <= 5; test++ )
  {
    xc = comp_random_grlex ( kc, rank1, rank2, &seed, &rank3 );
    rank4 = comp_rank_grlex ( kc, xc );
    printf ( "  %4d  %6d  %8d\n", test, rank3, rank4 );
    free ( xc );
  }
  return;
}
/******************************************************************************/

void comp_unrank_grlex_test ( )

/******************************************************************************/
/*
  Purpose:

    COMP_UNRANK_GRLEX_TEST tests COMP_UNRANK_GRLEX.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 December 2013

  Author:

    John Burkardt
*/
{
  int j;
  int kc = 3;
  int nc;
  int rank1;
  int rank2;
  int *xc;

  printf ( "\n" );
  printf ( "COMP_UNRANK_GRLEX_TEST\n" );
  printf ( "  A COMP is a composition of an integer N into K parts.\n" );
  printf ( "  Each part is nonnegative.  The order matters.\n" );
  printf ( "  COMP_UNRANK_GRLEX determines the parts\n" );
  printf ( "  of a COMP from its rank.\n" );
 
  printf ( "\n" );
  printf ( "  Rank: ->  NC       COMP    \n" );
  printf ( "  ----:     --   ------------ \n" );

  for ( rank1 = 1; rank1 <= 71; rank1++ )
  {
    xc = comp_unrank_grlex ( kc, rank1 );
    nc = i4vec_sum ( kc, xc );
    printf ( "   %3d: ", rank1 );
    printf ( "    %2d = ", nc );
    for ( j = 0; j < kc - 1; j++ )
    {
      printf ( "%2d + ", xc[j] );
    }
    printf ( "%2d\n", xc[kc-1] ); 
/*
  When XC(1) == NC, we have completed the compositions associated with
  a particular integer, and are about to advance to the next integer.
*/
    if ( xc[0] == nc )
    {
      printf ( "  ----:     --   ------------\n" );
    }
    free ( xc );
  }
  return;
}
/******************************************************************************/

void i4_choose_test ( )

/******************************************************************************/
/*
  Purpose:

    I4_CHOOSE_TEST tests I4_CHOOSE.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    27 October 2014

  Author:

    John Burkardt
*/
{
  int cnk;
  int k;
  int n;

  printf ( "\n" );
  printf ( "I4_CHOOSE_TEST\n" );
  printf ( "  I4_CHOOSE evaluates C(N,K).\n" );
  printf ( "\n" );
  printf ( "       N       K     CNK\n" );

  for ( n = 0; n <= 4; n++ )
  {
    printf ( "\n" );
    for ( k = 0; k <= n; k++ )
    {
      cnk = i4_choose ( n, k );

      printf ( "  %6d  %6d  %6d\n", n, k, cnk );
    }
  }

  return;
}
/******************************************************************************/

void i4_uniform_ab_test ( )

/******************************************************************************/
/*
  Purpose:

    I4_UNIFORM_TEST tests I4_UNIFORM_AB.

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
  printf ( "I4_UNIFORM_TEST\n" );
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

void i4vec_permute_test ( )

/******************************************************************************/
/*
  Purpose:

    I4VEC_PERMUTE_TEST tests I4VEC_PERMUTE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    25 October 2014

  Author:

    John Burkardt
*/
{
  int *a;
  int b;
  int c;
  int n = 12;
  int *p;
  int seed;

  printf ( "\n" );
  printf ( "I4VEC_PERMUTE_TEST\n" );
  printf ( "  I4VEC_PERMUTE reorders an I4VEC\n" );
  printf ( "  according to a given permutation.\n" );

  b = 0;
  c = n;
  seed = 123456789;
  a = i4vec_uniform_ab_new ( n, b, c, &seed );

  i4vec_print ( n, a, "  A[*], before rearrangement:" );

  p = perm_uniform_new ( n, &seed );

  i4vec_print ( n, p, "  Permutation vector P[*]:" );

  i4vec_permute ( n, p, a );

  i4vec_print ( n, a, "  A[P[*]]:" );

  free ( a );
  free ( p );

  return;
}
/******************************************************************************/

void i4vec_print_test ( )

/******************************************************************************/
/*
  Purpose:

    I4VEC_PRINT_TEST tests I4VEC_PRINT.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 October 2014

  Author:

    John Burkardt
*/
{
  int n = 4;
  int v[4] = { 91, 92, 93, 94 };

  printf ( "\n" );
  printf ( "I4VEC_PRINT_TEST\n" );
  printf ( "  I4VEC_PRINT prints an I4VEC\n" );

  i4vec_print ( n, v, "  Here is the I4VEC:" );

  return;
}
/******************************************************************************/

void i4vec_sort_heap_index_a_test ( )

/******************************************************************************/
/*
  Purpose:

    I4VEC_SORT_HEAP_INDEX_A_TEST tests I4VEC_SORT_HEAP_INDEX_A.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    27 October 2014

  Author:

    John Burkardt
*/
{
  int *a;
  int b;
  int c;
  int i;
  int *indx;
  int n = 20;
  int seed;

  printf ( "\n" );
  printf ( "I4VEC_SORT_HEAP_INDEX_A_TEST\n" );
  printf ( "  I4VEC_SORT_HEAP_INDEX_A creates an ascending\n" );
  printf ( "  sort index for an I4VEC.\n" );

  b = 0;
  c = 3 * n;
  seed = 123456789;

  a = i4vec_uniform_ab_new ( n, b, c, &seed );

  i4vec_print ( n, a, "  Unsorted array A:" );

  indx = i4vec_sort_heap_index_a ( n, a );

  i4vec_print ( n, indx, "  Sort vector INDX:" );

  printf ( "\n" );
  printf ( "       I   INDX(I)  A(INDX(I))\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %8d  %8d  %8d\n", i, indx[i], a[indx[i]] );
  }

  free ( a );
  free ( indx );

  return;
}
/******************************************************************************/

void i4vec_sum_test ( )

/******************************************************************************/
/*
  Purpose:

    I4VEC_SUM_TEST tests I4VEC_SUM.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    27 October 2014

  Author:

    John Burkardt
*/
{
  int *a;
  int hi;
  int lo;
  int n;
  int s;
  int seed;

  printf ( "\n" );
  printf ( "I4VEC_SUM_TEST\n" );
  printf ( "  I4VEC_SUM sums the entries of an I4VEC.\n" );

  n = 5;
  lo = 0;
  hi = 10;
  seed = 123456789;

  a = i4vec_uniform_ab_new ( n, lo, hi, &seed );
  i4vec_print ( n, a, "  The vector:" );

  s = i4vec_sum ( n, a );
  printf ( "\n" );
  printf ( "  The vector entries sum to %4d\n", s );

  free ( a );

  return;
}
/******************************************************************************/

void i4vec_uniform_ab_new_test ( )

/******************************************************************************/
/*
  Purpose:

    I4VEC_UNIFORM_AB_NEW_TEST tests I4_UNIFORM_AB.

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
  printf ( "\n" );

  v = i4vec_uniform_ab_new ( n, a, b, &seed );

  i4vec_print ( n, v, "  The vector:" );

  free ( v );

  return;
}
/******************************************************************************/

void lp_coefficients_test ( )

/******************************************************************************/
/*
  Purpose:

    LP_COEFFICIENTS_TEST tests LP_COEFFICIENTS.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 October 2014

  Author:

    John Burkardt
*/
{
  double *c;
  int *e;
  int *f;
  int i;
  char label[255];
  int m = 1;
  int n;
  int n_max = 10;
  int o;

  printf ( "\n" );
  printf ( "LP_COEFFICIENTS_TEST\n" );
  printf ( "  LP_COEFFICIENTS: coefficients of Legendre polynomial P(n,x).\n" );
  printf ( "\n" );

  for ( n = 0; n <= n_max; n++ )
  {
    c = ( double * ) malloc ( ( n + 1 ) * sizeof ( double ) );
    f = ( int * ) malloc ( ( n + 1 ) * sizeof ( int ) );

    lp_coefficients ( n, &o, c, f );

    e = ( int * ) malloc ( o * sizeof ( int ) );
    for ( i = 0; i < o; i++ )
    {
      e[i] = f[i] + 1;
    }
    sprintf ( label, "  P(%d,x) = ", n );
    polynomial_print ( m, o, c, e, label );

    free ( c );
    free ( e );
    free ( f );
   }

  return;
}
/******************************************************************************/

void lp_value_test ( )

/******************************************************************************/
/*
  Purpose:

    LP_VALUE_TEST tests LP_VALUE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 September 2014

  Author:

    John Burkardt
*/
{
  double e;
  int n;
  int n_data;
  int o;
  double x;
  double xvec[1];
  double fx1;
  double *fx2;

  n = 1;

  printf ( "\n" );
  printf ( "LP_VALUE_TEST:\n" );
  printf ( "  LP_VALUE evaluates a Legendre polynomial.\n" );
  printf ( "\n" );
  printf ( "                        Tabulated                 Computed\n" );
  printf ( "     O        X           L(O,X)                    L(O,X)" );
  printf ( "                   Error\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    lp_values ( &n_data, &o, &x, &fx1 );

    if ( n_data == 0 )
    {
      break;
    }
    xvec[0] = x;

    fx2 = lp_value ( n, o, xvec );

    e = fx1 - fx2[0];

    printf ( "  %4d  %12.8f  %24.16g  %24.16g  %8.2g\n",
      o, x, fx1, fx2[0], e );

    free ( fx2 );
  }

  return;
}
/******************************************************************************/

void lp_values_test ( )

/******************************************************************************/
/*
  Purpose:

    LP_VALUES_TEST tests LP_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 October 2014

  Author:

    John Burkardt
*/
{
  int n_data;
  int o;
  double x;
  double fx;

  printf ( "\n" );
  printf ( "LP_VALUES_TEST:\n" );
  printf ( "  LP_VALUES stores values of\n" );
  printf ( "  the Legendre polynomial P(o,x).\n" );
  printf ( "\n" );
  printf ( "                        Tabulated\n" );
  printf ( "     O        X           L(O,X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    lp_values ( &n_data, &o, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %4d  %12.8f  %24.16g \n", o, x, fx );
  }

  return;
}
/******************************************************************************/

void lpp_to_polynomial_test ( )

/******************************************************************************/
/*
  Purpose:

    LPP_TO_POLYNOMIAL_TEST tests LPP_TO_POLYNOMIAL.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 September 2014

  Author:

    John Burkardt
*/
{
  double *c;
  int *e;
  int i;
  int *l;
  char label[255];
  int m = 2;
  int o;
  int o_max;
  int rank;

  printf ( "\n" );
  printf ( "LPP_TO_POLYNOMIAL_TEST:\n" );
  printf ( "  LPP_TO_POLYNOMIAL is given a Legendre product polynomial\n" );
  printf ( "  and determines its polynomial representation.\n" );

  printf ( "\n" );
  printf ( "  Using spatial dimension M = %d\n", m );

  for ( rank = 1; rank <= 11; rank++ )
  {
    l = comp_unrank_grlex ( m, rank );

    o_max = 1;
    for ( i = 0; i < m; i++ )
    {
      o_max = o_max * ( l[i] + 2 ) / 2;
    }

    c = ( double * ) malloc ( o_max * sizeof ( double ) );
    e = ( int * ) malloc ( o_max * sizeof ( int ) );

    lpp_to_polynomial ( m, l, o_max, &o, c, e );

    sprintf ( label, "  LPP #%d = L(%d,X)*L(%d,Y) = ", rank, l[0], l[1] );

    printf ( "\n" );
    polynomial_print ( m, o, c, e, label );

    free ( c );
    free ( e );
    free ( l );
  }

  return;
}
/******************************************************************************/

void lpp_value_test ( )

/******************************************************************************/
/*
  Purpose:

    LPP_VALUE_TEST tests LPP_VALUE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 September 2014

  Author:

    John Burkardt
*/
{
  double *c;
  int *e;
  int i;
  int *l;
  int m = 3;
  int n = 1;
  int o;
  int o_max;
  int rank;
  int seed;
  double *v1;
  double *v2;
  double *x;
  double xhi;
  double xlo;

  printf ( "\n" );
  printf ( "LPP_VALUE_TEST:\n" );
  printf ( "  LPP_VALUE evaluates a Legendre product polynomial.\n" );

  xlo = -1.0;
  xhi = +1.0;
  seed = 123456789;
  x = r8vec_uniform_ab_new ( m, xlo, xhi, &seed );

  printf ( "\n" );
  printf ( "  Evaluate at X = " );
  for ( i = 0; i < m; i++ )
  {
    printf ( "%g  ", x[i+0*m] );
  }
  printf ( "\n" );
  printf ( "\n" );
  printf ( "  Rank  I1  I2  I3:  L(I1,X1)*L(I2,X2)*L(I3,X3)    P(X1,X2,X3)\n" );
  printf ( "\n" );

  for ( rank = 1; rank <= 20; rank++ )
  {
    l = comp_unrank_grlex ( m, rank );
/*
  Evaluate the LPP directly.
*/
    v1 = lpp_value ( m, n, l, x );
/*
  Convert the LPP to a polynomial.
*/
    o_max = 1;
    for ( i = 0; i < m; i++ )
    {
      o_max = o_max * ( l[i] + 2 ) / 2;
    }

    c = ( double * ) malloc ( o_max * sizeof ( double ) );
    e = ( int * ) malloc ( o_max * sizeof ( int ) );
 
    lpp_to_polynomial ( m, l, o_max, &o, c, e );
/*
  Evaluate the polynomial.
*/
    v2 = polynomial_value ( m, o, c, e, n, x );
/*
  Compare results.
*/
    printf ( "  %4d  %2d  %2d  %2d  %14.6g  %14.6g\n",
      rank, l[0], l[1], l[2], v1[0], v2[0] );
 
    free ( c );
    free ( e );
    free ( l );
    free ( v1 );
    free ( v2 );
  }

  free ( x );

  return;
}
/******************************************************************************/

void mono_next_grlex_test ( )

/******************************************************************************/
/*
  Purpose:

    MONO_NEXT_GRLEX_TEST tests MONO_NEXT_GRLEX.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 December 2013

  Author:

    John Burkardt
*/
{
  int a;
  int b;
  int m = 4;
  int i;
  int j;
  int k;
  int seed;
  int *x;

  printf ( "\n" );
  printf ( "MONO_NEXT_GRLEX_TEST\n" );
  printf ( "  MONO_NEXT_GRLEX computes the next monomial\n" );
  printf ( "  in M variables, in grlex order.\n" );
  printf ( "\n" );
  printf ( "  Let M =  %d\n", m );

  a = 0;
  b = 3;
  seed = 123456789;

  for ( i = 1; i <= 10; i++ )
  {
    x = i4vec_uniform_ab_new ( m, a, b, &seed );
    printf ( "\n" );
    printf ( "  " );
    for ( k = 0; k < m; k++ )
    {
      printf ( "%2d", x[k] );
    }
    printf ( "\n" );

    for ( j = 1; j <= 5; j++ )
    {
      mono_next_grlex ( m, x );
      printf ( "  " );
      for ( k = 0; k < m; k++ )
      {
        printf ( "%2d", x[k] );
      }
      printf ( "\n" );
    }
    free ( x );
  }

  return;
}
/******************************************************************************/

void mono_print_test ( )

/******************************************************************************/
/*
  Purpose:

    MONO_PRINT_TEST tests MONO_PRINT.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    03 November 2014

  Author:

    John Burkardt
*/
{
  int f1[1] = { 5 };
  int f2[1] = { -5 };
  int f3[4] = { 2, 1, 0, 3 };
  int f4[3] = { 17, -3, 199 };
  int m;

  printf ( "\n" );
  printf ( "MONO_PRINT_TEST\n" );
  printf ( "  MONO_PRINT can print out a monomial.\n" );
  printf ( "\n" );

  m = 1;
  mono_print ( m, f1, "  Monomial [5]:" );

  m = 1;
  mono_print ( m, f2, "  Monomial [5]:" );

  m = 4;
  mono_print ( m, f3, "  Monomial [2,1,0,3]:" );

  m = 3;
  mono_print ( m, f4, "  Monomial [17,-3,199]:" );

  return;
}
/******************************************************************************/

void mono_rank_grlex_test ( )

/******************************************************************************/
/*
  Purpose:

    MONO_RANK_GRLEX_TEST tests MONO_RANK_GRLEX.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 December 2013

  Author:

    John Burkardt
*/
{
  int m = 3;
  int i;
  int j;
  int n;
  int rank;
  int test;
  int test_num = 8;
  int x[3];
  int x_test[3*8] = {
    0, 0, 0, 
    1, 0, 0, 
    0, 0, 1, 
    0, 2, 0, 
    1, 0, 2, 
    0, 3, 1, 
    3, 2, 1, 
    5, 2, 1 };

  printf ( "\n" );
  printf ( "MONO_RANK_GRLEX_TEST\n" );
  printf ( "  MONO_RANK_GRLEX returns the rank of a monomial in the sequence\n" );
  printf ( "  of all monomials in M dimensions, in grlex order.\n" );

  printf ( "\n" );
  printf ( "  Print a monomial sequence with ranks assigned.\n" );

  n = 4;

  printf ( "\n" );
  printf ( "  Let M = %d\n", m );
  printf ( "      N = %d\n", n );
  printf ( "\n" );

  x[0] = 0;
  x[1] = 0;
  x[2] = 0; 

  i = 1;

  for ( ; ; )
  {
    printf ( "  %2d    ", i );
    for ( j = 0; j < m; j++ )
    {
      printf ( "%2d", x[j] );
    }
    printf ( "\n" );

    if ( x[0] == n ) 
    {
      break;
    }

    mono_upto_next_grlex ( m, n, x );
    i = i + 1;
  }

  printf ( "\n" );
  printf ( "  Now, given a monomial, retrieve its rank in the sequence:\n" );
  printf ( "\n" );

  for ( test = 0; test < test_num; test++ )
  {
    for ( j = 0; j < m; j++ )
    {
      x[j] = x_test[j+test*m];
    }
    rank = mono_rank_grlex ( m, x );

    printf ( "  %3d    ", rank );
    for ( j = 0; j < m; j++ )
    {
      printf ( "%2d", x[j] );
    }
    printf ( "\n" );
  }

  return;
}
/******************************************************************************/

void mono_unrank_grlex_test ( )

/******************************************************************************/
/*
  Purpose:

    MONO_UNRANK_GRLEX_TEST tests MONO_UNRANK_GRLEX.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 December 2013

  Author:

    John Burkardt
*/
{
  int m = 3;
  int i;
  int j;
  int n;
  int rank;
  int rank_max;
  int seed;
  int test;
  int test_num;
  int *x;

  printf ( "\n" );
  printf ( "MONO_UNRANK_GRLEX_TEST\n" );
  printf ( "  MONO_UNRANK_GRLEX is given a rank, and returns the corresponding\n" );
  printf ( "  monomial in the sequence of all monomials in M dimensions\n" );
  printf ( "  in grlex order.\n" );

  printf ( "\n" );
  printf ( "  For reference, print a monomial sequence with ranks.\n" );

  n = 4;
  rank_max = mono_upto_enum ( m, n );

  printf ( "\n" );
  printf ( "  Let M = %d\n", m );
  printf ( "      N = %d\n", n );
  printf ( "\n" );

  x = ( int * ) malloc ( m * sizeof ( int ) );
  for ( i = 0; i < m; i++ )
  {
    x[i] = 0;
  }

  i = 1;

  for ( ; ; )
  {
    printf ( "  %2d    ", i );
    for ( j = 0; j < m; j++ )
    {
      printf ( "%2d", x[j] );
    }
    printf ( "\n" );

    if ( x[0] == n ) 
    {
      break;
    }

    mono_upto_next_grlex ( m, n, x );
    i = i + 1;
  }

  printf ( "\n" );
  printf ( "  Now choose random ranks between 1 and %d\n", rank_max );
  printf ( "\n" );

  seed = 123456789;
  test_num = 5;

  for ( test = 1; test <= test_num; test++ )
  {
    rank = i4_uniform_ab ( 1, rank_max, &seed );   
    x = mono_unrank_grlex ( m, rank );
    printf ( "  %2d    ", rank );
    for ( j = 0; j < m; j++ )
    {
      printf ( "%2d", x[j] );
    }
    printf ( "\n" );
    free ( x );
  }

  return;
}
/******************************************************************************/

void mono_upto_enum_test ( )

/******************************************************************************/
/*
  Purpose:

    MONO_UPTO_ENUM_TEST tests MONO_UPTO_ENUM.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    18 November 2013

  Author:

    John Burkardt
*/
{
  int m;
  int n;
  int v;

  printf ( "\n" );
  printf ( "MONO_UPTO_ENUM_TEST\n" );
  printf ( "  MONO_UPTO_ENUM can enumerate the number of monomials\n" );
  printf ( "  in M variables, of total degree 0 up to N.\n" );

  printf ( "\n" );
  printf ( "    N:\n" );
  for ( n = 0; n <= 8; n++ )
  {
    printf ( "  %4d", n );
  }
  printf ( "\n" );
  printf ( "   M +------------------------------------------------------\n" );
  for ( m = 1; m <= 8; m++ )
  {
    printf ( "  %2d  |", m );
    for ( n = 0; n <= 8; n++ )
    {
      v = mono_upto_enum ( m, n );
      printf ( " %5d", v );
    }
    printf ( "\n" );
  }
  return;
}
/******************************************************************************/

void mono_upto_random_test ( )

/******************************************************************************/
/*
  Purpose:

    MONO_UPTO_RANDOM_TEST tests MONO_UPTO_RANDOM.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    21 November 2013

  Author:

    John Burkardt
*/
{
  int m = 3;
  int j;
  int n;
  int rank;
  int seed;
  int test;
  int test_num;
  int *x;

  printf ( "\n" );
  printf ( "MONO_UPTO_RANDOM_TEST\n" );
  printf ( "  MONO_UPTO_RANDOM selects at random a monomial\n" );
  printf ( "  in M dimensions of total degree no greater than N.\n" );

  n = 4;

  printf ( "\n" );
  printf ( "  Let M = %d\n", m );
  printf ( "      N = %d\n", n );
  printf ( "\n" );

  seed = 123456789;
  test_num = 5;

  for ( test = 1; test <= test_num; test++ )
  {
    x = mono_upto_random ( m, n, &seed, &rank );
    printf ( "  %2d    ", rank );
    for ( j = 0; j < m; j++ )
    {
      printf ( "%2d", x[j] );
    }
    printf ( "\n" );
    free ( x );
  }

  return;
}
/******************************************************************************/

void mono_upto_next_grlex_test ( )

/******************************************************************************/
/*
  Purpose:

    MONO_UPTO_NEXT_GRLEX_TEST tests MONO_UPTO_NEXT_GRLEX.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 December 2013

  Author:

    John Burkardt
*/
{
  int m = 3;
  int i;
  int j;
  int n;
  int x[3];

  printf ( "\n" );
  printf ( "MONO_UPTO_NEXT_GRLEX_TEST\n" );
  printf ( "  MONO_UPTO_NEXT_GRLEX can list the monomials\n" );
  printf ( "  in M variables, of total degree up to N,\n" );
  printf ( "  in grlex order, one at a time.\n" );
  printf ( "\n" );
  printf ( "  We start the process with (0,0,..0,0).\n" );
  printf ( "  The process ends with (N,0,...,0,0)\n" );

  n = 4;

  printf ( "\n" );
  printf ( "  Let M = %d\n", m );
  printf ( "      N = %d\n", n );
  printf ( "\n" );

  x[0] = 0;
  x[1] = 0;
  x[2] = 0; 

  i = 1;

  for ( ; ; )
  {
    printf ( "  %2d    ", i );
    for ( j = 0; j < m; j++ )
    {
      printf ( "%2d", x[j] );
    }
    printf ( "\n" );

    if ( x[0] == n ) 
    {
      break;
    }

    mono_upto_next_grlex ( m, n, x );
    i = i + 1;

  }

  return;
}
/******************************************************************************/

void mono_value_test ( )

/******************************************************************************/
/*
  Purpose:

    MONO_VALUE_TEST tests MONO_VALUE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    10 December 2013

  Author:

    John Burkardt
*/
{
  int m = 3;
  int *f;
  int j;
  int n;
  int nx = 2;
  int rank;
  int seed;
  int test;
  int test_num;
  double *v;
  double x[3*2] = {
     1.0, 2.0, 3.0, 
    -2.0, 4.0, 1.0 };

  printf ( "\n" );
  printf ( "MONO_VALUE_TEST\n" );
  printf ( "  MONO_VALUE evaluates a monomial.\n" );

  n = 6;

  printf ( "\n" );
  printf ( "  Let M = %d\n", m );
  printf ( "      N = %d\n", n );

  seed = 123456789;
  test_num = 5;

  for ( test = 1; test <= test_num; test++ )
  {
    f = mono_upto_random ( m, n, &seed, &rank );
    printf ( "\n" );
    mono_print ( m, f, "  M(X) = " );
    v = mono_value ( m, nx, f, x );
    for ( j = 0; j < nx; j++ )
    {
      printf ( "  M(%g,%g,%g) = %g\n", x[0+j*m], x[1+j*m], x[2+j*m], v[j] );
    }
    free ( f );
    free ( v );
  }

  return;
}
/******************************************************************************/

void perm_uniform_test ( )

/******************************************************************************/
/*
  Purpose:

    PERM_UNIFORM_TEST tests PERM_UNIFORM.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    25 October 2014

  Author:

    John Burkardt
*/
{
  int i;
  int n = 10;
  int *p;
  int seed;
  int test;

  printf ( "\n" );
  printf ( "PERM_UNIFORM_TEST\n" );
  printf ( "  PERM_UNIFORM randomly selects a permutation.\n" );
  printf ( "\n" );

  seed = 123456789;

  for ( test = 1; test <= 5; test++ )
  {
    p = perm_uniform_new ( n, &seed );
    printf ( "  " );
    for ( i = 0; i < n; i++ )
    {
      printf ( "%4d", p[i] );
    }
    printf ( "\n" );
    free ( p );
  }
  return;
}
/******************************************************************************/

void polynomial_compress_test ( )

/******************************************************************************/
/*
  Purpose:

    POLYNOMIAL_COMPRESS_TEST tests POLYNOMIAL_COMPRESS.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    27 October 2014

  Author:

    John Burkardt
*/
{
  double c[10] = { 7.0, - 5.0, 5.0, 9.0, 11.0, 3.0, 6.0, 0.0, - 13.0, 1.0E-20 };
  double c2[10];
  int m = 3;
  int e[10] = { 1, 2, 2, 4, 5, 5, 5, 12, 33, 35 }; 
  int e2[10];
  int j;
  int nx = 2;
  int o = 10;
  int o2;
  char title[255];

  printf ( "\n" );
  printf ( "POLYNOMIAL_COMPRESS_TEST\n" );
  printf ( "  POLYNOMIAL_COMPRESS compresses a polynomial.\n" );

  printf ( "\n" );
  strcpy ( title, "  Uncompressed P(X) = " );
  polynomial_print ( m, o, c, e, title );

  polynomial_compress ( o, c, e, &o2, c2, e2 );

  printf ( "\n" );
  strcpy ( title, "  Compressed P(X) = " );
  polynomial_print ( m, o2, c2, e2, title );

  return;
}
/******************************************************************************/

void polynomial_print_test ( )

/******************************************************************************/
/*
  Purpose:

    POLYNOMIAL_PRINT_TEST tests POLYNOMIAL_PRINT.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 December 2013

  Author:

    John Burkardt
*/
{
  double c[6] = { 7.0, - 5.0, 9.0, 11.0, 0.0, - 13.0 };
  int m = 3;
  int e[6] = { 1, 2, 4, 5, 12, 33 };
  int o = 6;
  char title[] = "  P1(X) =";

  printf ( "\n" );
  printf ( "POLYNOMIAL_PRINT_TEST\n" );
  printf ( "  POLYNOMIAL_PRINT prints a polynomial.\n" );

  printf ( "\n" );
  polynomial_print ( m, o, c, e, title );

  return;
}
/******************************************************************************/

void polynomial_sort_test ( )

/******************************************************************************/
/*
  Purpose:

    POLYNOMIAL_SORT_TEST tests POLYNOMIAL_SORT.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    27 October 2014

  Author:

    John Burkardt
*/
{
  double c[6] = { 0.0, 9.0, -5.0, - 13.0, 7.0, 11.0 };
  int m = 3;
  int e[6] = { 12, 4, 2, 33, 1, 5 };
  int o = 6;
  char title[255];

  printf ( "\n" );
  printf ( "POLYNOMIAL_SORT_TEST\n" );
  printf ( "  POLYNOMIAL_SORT sorts a polynomial by exponent index.\n" );

  printf ( "\n" );
  strcpy ( title, "  Unsorted polynomial:" );
  polynomial_print ( m, o, c, e, title );

  polynomial_sort ( o, c, e );

  printf ( "\n" );
  strcpy ( title, "  Sorted polynomial:" );
  polynomial_print ( m, o, c, e, title );

  return;
}
/******************************************************************************/

void polynomial_value_test ( )

/******************************************************************************/
/*
  Purpose:

    POLYNOMIAL_VALUE_TEST tests POLYNOMIAL_VALUE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    27 October 2014

  Author:

    John Burkardt
*/
{
  double c[6] = { 7.0, - 5.0, 9.0, 11.0, 0.0, - 13.0 };
  int m = 3;
  int e[6] = { 1, 2, 4, 5, 12, 33 }; 
  int j;
  int nx = 2;
  int o = 6;
  double *p;
  char title[] = "  P(X) =";
  double x[3*2] = {
     1.0, 2.0, 3.0,
    -2.0, 4.0, 1.0 };

  printf ( "\n" );
  printf ( "POLYNOMIAL_VALUE_TEST\n" );
  printf ( "  POLYNOMIAL_VALUE evaluates a polynomial.\n" );

  printf ( "\n" );
  polynomial_print ( m, o, c, e, title );

  p = polynomial_value ( m, o, c, e, nx, x );

  printf ( "\n" );
  for ( j = 0; j < nx; j++ )
  {
    printf ( "  P(%f,%f,%f) = %g\n",x[0+j*m], x[1+j*m], x[2+j*m], p[j] );
  }

  free ( p );

  return;
}
/******************************************************************************/

void r8mat_print_test ( )

/******************************************************************************/
/*
  Purpose:

    R8MAT_PRINT_TEST tests R8MAT_PRINT.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 August 2014

  Author:

    John Burkardt
*/
{
# define M 6
# define N 4

  double a[M*N];
  int i;
  int j;
  int m = M;
  int n = N;

  printf ( "\n" );
  printf ( "R8MAT_PRINT_TEST\n" );
  printf ( "  R8MAT_PRINT prints an R8MAT.\n" );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      a[i+j*m] = ( double ) ( ( i + 1 ) * 10 + ( j + 1 ) );
    }
  }
  r8mat_print ( m, n, a, "  The matrix:" );

  return;
# undef M
# undef N
}
/******************************************************************************/

void r8mat_print_some_test ( )

/******************************************************************************/
/*
  Purpose:

    R8MAT_PRINT_SOME_TEST tests R8MAT_PRINT_SOME.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 August 2014

  Author:

    John Burkardt
*/
{
# define M 6
# define N 4

  double a[M*N];
  int i;
  int j;
  int m = M;
  int n = N;

  printf ( "\n" );
  printf ( "R8MAT_PRINT_SOME_TEST\n" );
  printf ( "  R8MAT_PRINT_SOME prints some of an R8MAT.\n" );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      a[i+j*m] = ( double ) ( ( i + 1 ) * 10 + ( j + 1 ) );
    }
  }
  r8mat_print_some ( m, n, a, 2, 1, 4, 2, "  Rows 2:4, Cols 1:2:" );

  return;
# undef M
# undef N
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

    03 October 2005

  Author:

    John Burkardt
*/
{
# define M 5
# define N 4

  double *a;
  double b = 2.0E+00;
  double c = 10.0E+00;
  int seed = 123456789;

  printf ( "\n" );
  printf ( "R8MAT_UNIFORM_AB_NEW_TEST\n" );
  printf ( "  R8MAT_UNIFORM_AB_NEW sets an R8MAT to random values in [A,B].\n" );
  printf ( "\n" );

  a = r8mat_uniform_ab_new ( M, N, b, c, &seed );

  r8mat_print ( M, N, a, "  The random R8MAT:" );

  free ( a );

  return;
# undef M
# undef N
}
/******************************************************************************/

void r8vec_permute_test ( )

/******************************************************************************/
/*
  Purpose:

    R8VEC_PERMUTE_TEST tests R8VEC_PERMUTE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 October 2014

  Author:

    John Burkardt
*/
{
  int i;
  int n = 5;
  int p[5] = { 1, 3, 4, 0, 2 };
  double x[5] = { 1.1, 2.2, 3.3, 4.4, 5.5 };

  printf ( "\n" );
  printf ( "R8VEC_PERMUTE_TEST\n" );
  printf ( "  R8VEC_PERMUTE permutes an R8VEC.\n" );

  r8vec_print ( n, x, "  Original array X[]:" );

  i4vec_print ( n, p, "  Permutation vector P[]:" );

  r8vec_permute ( n, p, x );

  r8vec_print ( n, x, "  Permuted array X[P[]]:" );

  return;
}
/******************************************************************************/

void r8vec_print_test ( )

/******************************************************************************/
/*
  Purpose:

    R8VEC_PRINT_TEST tests R8VEC_PRINT.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 August 2014

  Author:

    John Burkardt
*/
{
  double a[4] = { 123.456, 0.000005, -1.0E+06, 3.14159265 };
  int n = 4;

  printf ( "\n" );
  printf ( "R8VEC_PRINT_TEST\n" );
  printf ( "  R8VEC_PRINT prints an R8VEC.\n" );

  r8vec_print ( n, a, "  The R8VEC:" );

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

    19 June 2012

  Author:

    John Burkardt
*/
{
# define N 10

  double a = 10.0;
  double b = 20.0;
  int j;
  double *r;
  int seed;

  printf ( "\n" );
  printf ( "R8VEC_UNIFORM_AB_NEW_TEST\n" );
  printf ( "  R8VEC_UNIFORM_AB_NEW returns a random R8VEC\n" );
  printf ( "  with entries in a given range [ A, B ]\n" );
  printf ( "\n" );
  printf ( "  For this problem:\n" );
  printf ( "  A = %g\n", a );
  printf ( "  B = %g\n", b );
  printf ( "\n" );

  seed = 123456789;

  for ( j = 1; j <= 3; j++ )
  {
    printf ( "\n" );
    printf ( "  Input SEED = %d\n", seed );
    printf ( "\n" );

    r = r8vec_uniform_ab_new ( N, a, b, &seed );

    r8vec_print ( N, r,  "  Random R8VEC:" );

    free ( r );
  }

  return;
# undef N
}
