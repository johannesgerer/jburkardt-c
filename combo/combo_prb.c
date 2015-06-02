# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "combo.h"

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );
void test04 ( );
void test05 ( );
void test06 ( );
void test07 ( );
void test08 ( );
void test09 ( );
void test10 ( );
void test11 ( );
void test12 ( );
void test13 ( );
void test14 ( );
void test15 ( );
void test16 ( );
void test17 ( );
void test18 ( );
void test19 ( );
void test20 ( );
void test21 ( );
void test22 ( );
void test23 ( );
void test24 ( );
void test25 ( );
void test26 ( );
void test27 ( );
void test28 ( );
void test29 ( );
void test30 ( );
void test31 ( );
void test32 ( );
void test33 ( );
void test34 ( );
void test35 ( );
void test36 ( );
void test37 ( );
void test38 ( );
void test39 ( );
void test40 ( );
void test41 ( );
void test42 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for COMBO_PRB.
//
//  Discussion:
//
//    COMBO_PRB tests the COMBO library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 August 2012
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  printf ( "\n" );
  printf ( "COMBO_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the COMBO library.\n" );

  test01 ( );
  test02 ( );
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
  test13 ( );
  test14 ( );
  test15 ( );
  test16 ( );
  test17 ( );
  test18 ( );
  test19 ( );

  test20 ( );
  test21 ( );
  test22 ( );
  test23 ( );
  test24 ( );
  test25 ( );
  test26 ( );
  test27 ( );
  test28 ( );
  test29 ( );

  test30 ( );
  test31 ( );
  test32 ( );
  test33 ( );
  test34 ( );
  test35 ( );
  test36 ( );
  test37 ( );
  test38 ( );
  test39 ( );

  test40 ( );
  test41 ( );
  test42 ( );
//
//  Terminate.
//
  printf ( "\n" );
  printf ( "COMBO_PRB\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
//****************************************************************************80

void test01 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST01 tests BAL_SEQ_ENUM, BAL_SEQ_RANK, BAL_SEQ_SUCCESSOR, BAL_SEQ_UNRANK.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 July 2011
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n = 5;
  int nseq;
  int rank;
  int rank_old;
  int *t;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  Balanced sequences:\n" );
  printf ( "\n" );
  printf ( "  BAL_SEQ_ENUM enumerates,\n" );
  printf ( "  BAL_SEQ_RANK ranks,\n" );
  printf ( "  BAL_SEQ_SUCCESSOR lists,\n" );
  printf ( "  BAL_SEQ_UNRANK unranks.\n" );
//
//  Enumerate.
//
  nseq = bal_seq_enum ( n );

  printf ( "\n" );
  printf ( "  For N = %d\n", n );
  printf ( "  the number of balanced sequences is %d\n", nseq );
  printf ( "\n" );
//
//  List.
//
  t = ( int * ) malloc ( 2 * n * sizeof ( int ) );

  rank = -1;

  for ( ; ; )
  {
    rank_old = rank;

    bal_seq_successor ( n, t, &rank );

    if ( rank <= rank_old )
    {
      break;
    }
    printf ( "  %4d", rank );
    for ( i = 0; i < 2*n; i++ )
    {
      printf ( "  %4d", t[i] );
    }
    printf ( "\n" );
  }
//
//  Unrank.
//
  rank = nseq / 2;

  free ( t );

  t = bal_seq_unrank ( rank, n );
  printf ( "\n" );
  printf ( "  The element of rank %d\n", rank );
  printf ( "\n" );
  i4vec_transpose_print ( 2 * n, t, " " );
//
//  Rank.
//
  rank = bal_seq_rank ( n, t );

  i4vec_transpose_print ( 2 * n, t, "  Element to be ranked:" );
  printf ( "\n" );
  printf ( "  Rank is computed as: %d\n", rank );

  free ( t );

  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests BAL_SEQ_TO_TABLEAU, TABLEAU_TO_BAL_SEQ.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 July 2011
//
//  Author:
//
//    John Burkardt
//
{
  int n = 4;
  int rank;
  int *t;
  int *tab;

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  BAL_SEQ_TO_TABLEAU converts a balanced\n" );
  printf ( "  sequence to a tableau;\n" );
  printf ( "  TABLEAU_TO_BAL_SEQ converts a tableau\n" );
  printf ( "  to a balanced sequence.\n" );
//
//  Pick a random balanced sequence.
//
  rank = 7;

  t = bal_seq_unrank ( rank, n );

  printf ( "\n" );
  printf ( "  Random balanced sequence:\n" );
  printf ( "\n" );
  i4vec_transpose_print ( 2 * n, t, " " );
//
//  Convert to a tableau.
//
  tab = bal_seq_to_tableau ( n, t );

  i4mat_print ( 2, n, tab, "  Corresponding tableau" );
//
//  Convert to a balanced sequence.
//
  free ( t );

  t = tableau_to_bal_seq ( n, tab );

  i4vec_transpose_print ( 2 * n, t, "  Corresponding balanced sequence:" );

  free ( t );
  free ( tab );

  return;
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 tests BELL_NUMBERS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 July 2011
//
//  Author:
//
//    John Burkardt
//
{
  int *b;
  int bn;
  int n;
  int n_data;

  n_data = 0;

  printf ( "\n" );
  printf ( "TEST03\n" );
  printf ( "  BELL_NUMBERS computes Bell numbers.\n" );
  printf ( "\n" );
  printf ( "     N          BELL(N)      BELL_NUMBERS(N)\n" );
  printf ( "\n" );
  for ( ; ; )
  {
    bell_values ( &n_data, &n, &bn );

    if ( n_data == 0 )
    {
      break;
    }
    b = bell_numbers ( n );
    printf ( "  %8d  %12d  %12d\n", n, bn, b[n] );
    free ( b );
  }

  return;
}
//****************************************************************************80

void test04 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST04 tests I4_CHOOSE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 July 2011
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int j;

  printf ( "\n" );
  printf ( "TEST04\n" );
  printf ( "  I4_CHOOSE computes binomial coefficients.\n" );

  for ( i = -1; i <= 5; i++ )
  {
    for ( j = - 1; j <= 5; j++ )
    {
      printf ( "  %4d  %4d  %12d\n", i, j, i4_choose ( i, j ) );
    }
  }
  return;
}
//****************************************************************************80

void test05 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST05 tests CYCLE_TO_PERM, PERM_TO_CYCLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 July 2011
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int j;
  int jlo;
  int *index;
  int n = 7;
  int ncycle;
  int nperm;
  int *p;
  int rank;
  int *t;

  printf ( "\n" );
  printf ( "TEST05\n" );
  printf ( "  CYCLE_TO_PERM converts a permutation from\n" );
  printf ( "  cycle to array form;\n" );
  printf ( "  PERM_TO_CYCLE converts a permutation from\n" );
  printf ( "  array to cycle form.\n" );
//
//  Enumerate.
//
  nperm = perm_enum ( n );
//
//  Choose a "random" permutation.
//
  rank = nperm / 2;

  p = perm_lex_unrank ( rank, n );

  perm_print ( n, p, "  Random permutation:" );
//
//  Convert the permutation to cycle form.
//
  t = ( int * ) malloc ( n * sizeof ( int ) );
  index = ( int * ) malloc ( n * sizeof ( int ) );

  perm_to_cycle ( n, p, &ncycle, t, index );

  printf ( "\n" );
  printf ( "  Corresponding cycle form:\n" );
  printf ( "  Number of cycles is %d\n", ncycle );
  printf ( "\n" );
  jlo = 0;
  for ( i = 1; i <= ncycle; i++ )
  {
    for ( j = jlo + 1; j <= jlo + index[i-1]; j++ )
    {
      printf ( "  %4d", t[j-1] );
    }
    printf ( "\n" );
    jlo = jlo + index[i-1];
  }
//
//  Convert the set partition back to an RGF.
//
  free ( p );

  p = cycle_to_perm ( n, ncycle, t, index );

  perm_print ( n, p, "  Corresponding permutation:" );

  free ( index );
  free ( p );
  free ( t );

  return;
}
//****************************************************************************80

void test06 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST06 tests DIST_ENUM and DIST_NEXT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 July 2011
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int idist;
  int k = 3;
  int m;
  int more;
  int num_dist;
  int *q;

  k = 3;
  m = 5;

  num_dist = dist_enum ( k, m );

  printf ( "\n" );
  printf ( "TEST06\n" );
  printf ( "  For a distribution of M indistinguishable\n" );
  printf ( "  objects among K distinguishable slots:\n" );
  printf ( "\n" );
  printf ( "  DIST_ENUM enumerates them;\n" );
  printf ( "  DIST_NEXT produces the \"next\" one.\n" );
  printf ( "\n" );
  printf ( "  Number of:\n" );
  printf ( "    indistinguishable objects = %d\n", m );
  printf ( "    distinguishable slots =     %d\n", k );
  printf ( "    distributions is            %d\n", num_dist );
  printf ( "\n" );

  idist = 0;
  more = 0;
  q = ( int * ) malloc ( k * sizeof ( int ) );

  for ( ; ; )
  {
    dist_next ( k, m, q, &more );

    if ( !more )
    {
      break;
    }

    idist = idist + 1;
    printf ( "  %4d", idist );
    for ( i = 0; i < k; i++ )
    {
      printf ( "  %2d", q[i] );
    }
    printf ( "\n" );
  }

  free ( q );

  return;
}
//****************************************************************************80

void test07 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST07 tests I4_FACTORIAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 July 2011
//
//  Author:
//
//    John Burkardt
//
{
  int fx;
  int fx2;
  int n;
  int x;

  printf ( "\n" );
  printf ( "TEST07:\n" );
  printf ( "  I4_FACTORIAL evaluates the factorial function.\n" );
  printf ( "\n" );
  printf ( "     X       Exact F       FACTORIAL(X)\n" );
  printf ( "\n" );

  n = 0;

  for ( ; ; )
  {
    i4_factorial_values ( &n, &x, &fx );

    if ( n == 0 )
    {
      break;
    }

    if ( x <= 0.0 )
    {
      continue;
    }

    fx2 = i4_factorial ( x );

    printf ( "  %4d  %12d  %12d\n", x, fx, fx2 );
  }
  return;
}
//****************************************************************************80

void test08 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST08 tests GRAY_CODE_*.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 July 2011
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n = 5;
  int ngray;
  int rank;
  int rank_old;
  int *t;

  t = ( int * ) malloc ( n * sizeof ( int ) );

  printf ( "\n" );
  printf ( "TEST08\n" );
  printf ( "  Gray codes:\n" );
  printf ( "\n" );
  printf ( "  GRAY_CODE_ENUM enumerates,\n" );
  printf ( "  GRAY_CODE_RANK ranks,\n" );
  printf ( "  GRAY_CODE_SUCCESSOR lists,\n" );
  printf ( "  GRAY_CODE_UNRANK unranks.\n" );
//
//  Enumerate.
//
  ngray = gray_code_enum ( n );

  printf ( "\n" );
  printf ( "  For N = %d\n", n );
  printf ( "  the number of Gray code elements is %d\n", ngray );
  printf ( "\n" );
//
//  List
//
  rank = -1;

  for ( ; ; )
  {
    rank_old = rank;

    gray_code_successor ( n, t, &rank );

    if ( rank <= rank_old )
    {
      break;
    }

    printf ( "  %4d", rank );
    for ( i = 0; i < n; i++ )
    {
      printf ( "  %4d", t[i] );
    }
    printf ( "\n" );
  }
//
//  Unrank.
//
  rank = ngray / 2;

  free ( t );

  t = gray_code_unrank ( rank, n );

  printf ( "\n" );
  printf ( "  The element of rank %d\n", rank );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %4d", t[i] );
  }
  printf ( "\n" );
//
//  Rank.
//
  rank = gray_code_rank ( n, t );

  i4vec_transpose_print ( n, t, "  Element to be ranked:" );

  printf ( "\n" );
  printf ( "  Computed rank: %d\n", rank );

  free ( t );

  return;
}
//****************************************************************************80

void test09 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST09 tests I4VEC_SEARCH_BINARY_A and I4VEC_SORT_INSERT_A.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 July 2011
//
//  Author:
//
//    John Burkardt
//
{
  int a[10] = { 6, 7, 1, 0, 4, 3, 2, 1, 5, 8 };
  int b;
  int index;
  int n = 10;

  printf ( "\n" );
  printf ( "TEST09\n" );
  printf ( "  Integer vectors:\n" );
  printf ( "\n" );
  printf ( "  I4VEC_SORT_INSERT_A ascending sorts;\n" );
  printf ( "  I4VEC_SEARCH_BINARY_A searches a ascending sorted vector.\n" );

  i4vec_print ( n, a, "  Before ascending sort:" );

  i4vec_sort_insert_a ( n, a );

  i4vec_print ( n, a, "  After ascending sort:" );

  b = 5;

  printf ( "\n" );
  printf ( "  Now search for an instance of the value %d\n", b );

  index = i4vec_search_binary_a ( n, a, b );

  printf ( "\n" );
  if ( index == 0 )
  {
    printf ( "  The value does not occur.\n" );
  }
  else
  {
    printf ( "  The value occurs at index = %d\n", index );
  }
  return;
}
//****************************************************************************80

void test10 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST10 tests I4VEC_SEARCH_BINARY_D and I4VEC_SORT_INSERT_D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 July 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 10

  int a[N] = { 6, 7, 1, 0, 4, 3, 2, 1, 5, 8 };
  int b;
  int index;
  int n = N;

  printf ( "\n" );
  printf ( "TEST10\n" );
  printf ( "  Integer vectors:\n" );
  printf ( "\n" );
  printf ( "  I4VEC_SORT_INSERT_D descending sorts;\n" );
  printf ( "  I4VEC_SEARCH_BINARY_D searches a descending \n" );
  printf ( " sorted vector.\n" );

  i4vec_print ( n, a, "  Before descending sort:" );

  i4vec_sort_insert_d ( n, a );

  i4vec_print ( n, a, "  After descending sort:" );

  b = 5;

  printf ( "\n" );
  printf ( "  Now search for an instance of the value %d\n", b );

  index = i4vec_search_binary_d ( n, a, b );

  printf ( "\n" );
  if ( index == 0 )
  {
    printf ( "  The value does not occur.\n" );
  }
  else
  {
    printf ( "  The value occurs at index = %d\n", index );
  }

  return;
# undef N
}
//****************************************************************************80

void test11 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST11 tests KNAPSACK_REORDER and KNAPSACK_01.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 July 2011
//
//  Author:
//
//    John Burkardt
//
{
# define N 5

  int i;
  double mass;
  double mass_limit = 26.0;
  int n = N;
  double p[N] = { 24.0, 13.0, 23.0, 15.0, 16.0 };
  double profit;
  double w[N] = { 12.0,  7.0, 11.0,  8.0,  9.0 };
  double x[N];

  printf ( "\n" );
  printf ( "TEST11\n" );
  printf ( "  KNAPSACK_REORDER reorders the knapsack data.\n" );
  printf ( "  KNAPSACK_01 solves the 0/1 knapsack problem.\n" );

  printf ( "\n" );
  printf ( "  Object, Profit, Mass, Profit Density\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %6d  %7g  %7g  %7g\n", i, p[i], w[i], p[i] / w[i] );
  }

  knapsack_reorder ( n, p, w );

  printf ( "\n" );
  printf ( "  After reordering by Profit Density:\n" );
  printf ( "\n" );
  printf ( "  Object, Profit, Mass, Profit Density\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %6d  %7g  %7g  %7g\n", i, p[i], w[i], p[i] / w[i] );
  }

  printf ( "\n" );
  printf ( "  Total mass restriction is %g\n", mass_limit );

  knapsack_01 ( n, mass_limit, p, w, x, &mass, &profit );

  printf ( "\n" );
  printf ( "  Object, Density, Choice, Profit, Mass\n" );
  printf ( "\n" );

  for ( i = 0; i < n; i++ )
  {
    printf ( "  %6d  %7g  %7g  %7g  %7g\n", 
      i, p[i] / w[i], x[i], x[i] * p[i], x[i] * w[i] );
  }

  printf ( "\n" );
  printf ( "  Total:            %g  %g", profit, mass );

  return;
# undef N
}
//****************************************************************************80

void test12 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST12 tests KNAPSACK_REORDER and KNAPSACK_RATIONAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 July 2011
//
//  Author:
//
//    John Burkardt
//
{
# define N 5

  int i;
  double mass;
  double mass_limit = 26.0;
  int n = N;
  double p[N] = { 24.0, 13.0, 23.0, 15.0, 16.0 };
  double profit;
  double w[N] = { 12.0,  7.0, 11.0,  8.0,  9.0 };
  double x[N];

  printf ( "\n" );
  printf ( "TEST12\n" );
  printf ( "  KNAPSACK_REORDER reorders the knapsack data.\n" );
  printf ( "  KNAPSACK_RATIONAL solves the rational knapsack problem.\n" );

  printf ( "\n" );
  printf ( "  Object, Profit, Mass, Profit Density\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %4d  %7g  %7g  %7g\n", i + 1, p[i], w[i], p[i] / w[i] );
  }

  knapsack_reorder ( n, p, w );

  printf ( "\n" );
  printf ( "  After reordering by Profit Density:\n" );
  printf ( "\n" );
  printf ( "  Object, Profit, Mass, Profit Density\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %4d  %7g  %7g  %7g\n", i + 1, p[i], w[i], p[i] / w[i] );
  }

  printf ( "\n" );
  printf ( "  Total mass restriction is %g\n", mass_limit );

  knapsack_rational ( n, mass_limit, p, w, x, &mass, &profit );

  printf ( "\n" );
  printf ( "  Object, Density, Choice, Profit, Mass\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %4d  %7g  %7g  %7g\n",
      i + 1, p[i] / w[i], x[i] * p[i], x[i] * w[i] );
  }

  printf ( "\n" );
  printf ( "  Total:            %g  %g\n", profit, mass );

  return;
# undef N
}
//****************************************************************************80

void test13 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST13 tests KSUBSET_COLEX_RANK, _SUCCESSOR, _UNRANK, _ENUM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 July 2011
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int k = 3;
  int n = 5;
  int nksub;
  int rank;
  int rank_old;
  int *t;

  printf ( "\n" );
  printf ( "TEST13\n" );
  printf ( "  K-subsets of an N set,\n" );
  printf ( "  using the colexicographic ordering:\n" );
  printf ( "\n" );
  printf ( "  KSUBSET_COLEX_RANK ranks,\n" );
  printf ( "  KSUBSET_COLEX_SUCCESSOR lists,\n" );
  printf ( "  KSUBSET_COLEX_UNRANK unranks.\n" );
  printf ( "  KSUBSET_ENUM enumerates,\n" );
//
//  Enumerate.
//
  nksub = ksubset_enum ( k, n );

  printf ( "\n" );
  printf ( "  For N = %d\n", n );
  printf ( "  the number of K subsets is %d\n", nksub );
  printf ( "\n" );
//
//  List
//
  t = ( int * ) malloc ( k * sizeof ( int ) );

  rank = -1;

  for ( ; ; )
  {
    rank_old = rank;

    ksubset_colex_successor ( k, n, t, &rank );

    if ( rank <= rank_old )
    {
      break;
    }

    printf ( "  %4d", rank );
    for ( i = 0; i < k; i++ )
    {
      printf ( "  %4d", t[i] );
    }
    printf ( "\n" );
  }
//
//  Unrank.
//
  rank = nksub / 2;

  free ( t );

  t = ksubset_colex_unrank ( rank, k, n );

  printf ( "\n" );
  printf ( "  The element of rank %d\n", rank );
  printf ( "\n" );
  i4vec_transpose_print ( k, t, " " );
//
//  Rank.
//
  rank = ksubset_colex_rank ( k, n, t );

  i4vec_transpose_print ( k, t, "  Element to be ranked:" );
  printf ( "\n" );
  printf ( "  Rank is computed as %d\n", rank );

  free ( t );

  return;
}
//****************************************************************************80

void test14 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST14 tests KSUBSET_ENUM, _LEX_RANK, _LEX_SUCCESSOR, _LEX_UNRANK.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 July 2011
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int k = 3;
  int n = 5;
  int nksub;
  int rank;
  int rank_old;
  int *t;

  printf ( "\n" );
  printf ( "TEST14\n" );
  printf ( "  K-subsets of an N set,\n" );
  printf ( "  using the lexicographic ordering:\n" );
  printf ( "\n" );
  printf ( "  KSUBSET_ENUM enumerates,\n" );
  printf ( "  KSUBSET_LEX_RANK ranks,\n" );
  printf ( "  KSUBSET_LEX_SUCCESSOR lists,\n" );
  printf ( "  KSUBSET_LEX_UNRANK unranks.\n" );
//
//  Enumerate.
//
  nksub = ksubset_enum ( k, n );

  printf ( "\n" );
  printf ( "  For N = %d\n", n );
  printf ( "  the number of K subsets is %d\n", nksub );
  printf ( "\n" );
//
//  List
//
  t = ( int * ) malloc ( k * sizeof ( int ) );

  rank = -1;

  for ( ; ; )
  {
    rank_old = rank;

    ksubset_lex_successor ( k, n, t, &rank );

    if ( rank <= rank_old )
    {
      break;
    }

    printf ( "  %4d", rank );
    for ( i = 0; i < k; i++ )
    {
      printf ( "  %4d", t[i] );
    }
    printf ( "\n" );
  }
//
//  Unrank.
//
  rank = nksub / 2;

  free ( t );

  t = ksubset_lex_unrank ( rank, k, n );

  printf ( "\n" );
  printf ( "  The element of rank %d\n", rank );
  printf ( "\n" );
  i4vec_transpose_print ( k, t, " " );
//
//  Rank.
//
  rank = ksubset_lex_rank ( k, n, t );

  i4vec_transpose_print ( k, t, "  Element to be ranked:" );
  printf ( "\n" );
  printf ( "  Rank is computed as %d\n", rank );

  free ( t );

  return;
}
//****************************************************************************80

void test15 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST15 tests KSUBSET_ENUM, _REVDOOR_RANK, _REVDOOR_SUCCESSOR, _REVDOOR_UNRANK.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 July 2011
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int k = 3;
  int n = 5;
  int nksub;
  int rank;
  int rank_old;
  int *t;

  printf ( "\n" );
  printf ( "TEST15\n" );
  printf ( "  K-subsets of an N set,\n" );
  printf ( "  using the revolving door ordering:\n" );
  printf ( "\n" );
  printf ( "  KSUBSET_ENUM enumerates,\n" );
  printf ( "  KSUBSET_REVDOOR_RANK ranks,\n" );
  printf ( "  KSUBSET_REVDOOR_SUCCESSOR lists,\n" );
  printf ( "  KSUBSET_REVDOOR_UNRANK unranks.\n" );
//
//  Enumerate.
//
  nksub = ksubset_enum ( k, n );

  printf ( "\n" );
  printf ( "  For N = %d\n", n );
  printf ( "  the number of K subsets is %d\n", nksub );
  printf ( "\n" );
//
//  List
//
  t = ( int * ) malloc ( k * sizeof ( int ) );

  rank = -1;

  for ( ; ; )
  {
    rank_old = rank;

    ksubset_revdoor_successor ( k, n, t, &rank );

    if ( rank <= rank_old )
    {
      break;
    }

    printf ( "  %4d", rank );
    for ( i = 0; i < k; i++ )
    {
      printf ( "  %4d", t[i] );
    }
    printf ( "\n" );
  }
//
//  Unrank.
//
  rank = nksub / 2;

  free ( t );

  t = ksubset_revdoor_unrank ( rank, k, n );

  printf ( "\n" );
  printf ( "  The element of rank %d\n", rank );
  printf ( "\n" );
  i4vec_transpose_print ( k, t, " " );
//
//  Rank.
//
  rank = ksubset_revdoor_rank ( k, n, t );

  i4vec_transpose_print ( k, t, "  Element to be ranked:" );
  printf ( "\n" );
  printf ( "  Rank is computed as %d\n", rank );

  free ( t );

  return;
}
//****************************************************************************80

void test16 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST16 tests MARRIAGE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 July 2011
//
//  Author:
//
//    John Burkardt
//
{
# define N 5

  int *fiancee;
  int i;
  int n = N;
  int *next;
  int prefer[N*N] = {
    2, 1, 2, 1, 5, 
    5, 2, 3, 3, 3, 
    1, 3, 5, 2, 2, 
    3, 4, 4, 4, 1, 
    4, 5, 1, 5, 4 };
  int rank[N*N] = {
    2, 4, 1, 4, 5, 
    4, 3, 3, 2, 2, 
    5, 5, 4, 1, 3, 
    3, 1, 2, 3, 1, 
    1, 2, 5, 5, 4 };

  printf ( "\n" );
  printf ( "TEST16\n" );
  printf ( "  MARRIAGE arranges a set of stable marriages\n" );
  printf ( "  given a set of preferences.\n" );

  fiancee = ( int * ) malloc ( n * sizeof ( int ) );
  next = ( int * ) malloc ( n * sizeof ( int ) );

  marriage ( n, prefer, rank, fiancee, next );

  printf ( "\n" );
  printf ( "  Man, Wife's rank, Wife\n" );
  printf ( "\n" );
  for ( i = 1; i <= n; i++ )
  {
    printf ( "  %4d  %4d  %4d\n", i, next[i-1], prefer[i-1+(next[i-1]-1)*n] );
  }

  printf ( "\n" );
  printf ( "  Woman, Husband's rank, Husband\n" );
  printf ( "\n" );
  for ( i = 1; i <= n; i++ )
  {
    printf ( "  %4d  %4d  %4d\n", i, rank[i-1+(fiancee[i-1]-1)*n], fiancee[i-1] );
  }

  printf ( "\n" );
  printf ( "  Correct result:\n" );
  printf ( "\n" );
  printf ( "  M:W 1  2  3  4  5\n" );
  printf ( "   1  +  .  .  .  .\n" );
  printf ( "   2  .  .  .  +  .\n" );
  printf ( "   3  .  .  .  .  +\n" );
  printf ( "   4  .  .  +  .  .\n" );
  printf ( "   5  .  +  .  .  .\n" );

  free ( fiancee );
  free ( next );

  return;
# undef N
}
//****************************************************************************80

void test17 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST17 tests MOUNTAIN.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 July 2011
//
//  Author:
//
//    John Burkardt
//
{
  int n = 5;
  int x;
  int y;

  printf ( "\n" );
  printf ( "TEST17\n" );
  printf ( "  MOUNTAIN computes mountain numbers.\n" );
  printf ( "\n" );
  printf ( "  Y  MXY\n" );
  printf ( "\n" );

  for ( y = 0; y <= n; y++ )
  {
    printf ( "  %2d   ", y );

    for ( x = 0; x <= 2 * n; x++ )
    {
      printf ( "  %4d", mountain ( n, x, y ) );
    }
    printf ( "\n" );
  }
  return;
}
//****************************************************************************80

void test18 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST18 tests NPART_ENUM, _RSF_LEX_RANK, _RSF_LEX_SUCCESSOR, _RSF_LEX_UNRANK.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 July 2011
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n;
  int npart = 3;
  int npartitions;
  int rank;
  int rank_old;
  int *t;

  n = 12;

  printf ( "\n" );
  printf ( "TEST18\n" );
  printf ( "  Partitions of N with NPART parts\n" );
  printf ( "  in reverse standard form:\n" );
  printf ( "\n" );
  printf ( "  NPART_ENUM enumerates,\n" );
  printf ( "  NPART_RSF_LEX_RANK ranks,\n" );
  printf ( "  NPART_RSF_LEX_SUCCESSOR lists;\n" );
  printf ( "  NPART_RSF_LEX_UNRANK unranks.\n" );
//
//  Enumerate.
//
  npartitions = npart_enum ( n, npart );

  printf ( "\n" );
  printf ( "  For N = %d\n", n );
  printf ( "  and NPART = %d\n", npart );
  printf ( "  the number of partitions is %d\n", npartitions );
  printf ( "\n" );
//
//  List.
//
  t = ( int * ) malloc ( npart * sizeof ( int ) );

  rank = -1;

  for ( ; ; )
  {
    rank_old = rank;

    npart_rsf_lex_successor ( n, npart, t, &rank );

    if ( rank <= rank_old )
    {
      break;
    }

    printf ( "  %4d", rank );
    for ( i = 0; i < npart; i++ )
    {
      printf ( "  %4d", t[i] );
    }
    printf ( "\n" );
  }
//
//  Unrank.
//
  rank = npartitions / 3;

  free ( t );

  t = npart_rsf_lex_unrank ( rank, n, npart );

  printf ( "\n" );
  printf ( "  The element of rank %d\n", rank );
  printf ( "\n" );
  i4vec_transpose_print ( npart, t, " " );
//
//  Rank.
//
  rank = npart_rsf_lex_rank ( n, npart, t );

  i4vec_transpose_print ( npart, t, "  Element to be ranked:" );
  printf ( "\n" );
  printf ( "  Rank is computed as %d\n", rank );

  free ( t );

  return;
}
//****************************************************************************80

void test19 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST19 tests NPART_RSF_LEX_RANDOM;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 July 2011
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n = 12;
  int npart = 3;
  int seed = 123456789;
  int *t;

  printf ( "\n" );
  printf ( "TEST19\n" );
  printf ( "  Partitions of N with NPART parts\n" );
  printf ( "  in reverse standard form:\n" );
  printf ( "\n" );
  printf ( "  NPART_RSF_LEX_RANDOM produces random examples.\n" );

  for ( i = 1; i <= 10; i++ )
  {
    t = npart_rsf_lex_random ( n, npart, &seed );
    i4vec_transpose_print ( npart, t, " " );
    free ( t );
  }

  return;
}
//****************************************************************************80

void test20 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST20 tests NPART_ENUM and NPART_SF_SUCCESSOR;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 July 2011
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n = 12;
  int npart = 3;
  int npartitions;
  int rank;
  int rank_old;
  int *t;

  printf ( "\n" );
  printf ( "TEST20\n" );
  printf ( "  Partitions of N with NPART parts\n" );
  printf ( "  in standard form:\n" );
  printf ( "\n" );
  printf ( "  NPART_ENUM enumerates,\n" );
  printf ( "  NPART_SF_LEX_SUCCESSOR lists.\n" );
//
//  Enumerate.
//
  npartitions = npart_enum ( n, npart );

  printf ( "\n" );
  printf ( "  For N = %d\n", n );
  printf ( "  and NPART = %d\n", npart );
  printf ( "  the number of partitions is %d\n", npartitions );
  printf ( "\n" );
//
//  List.
//
  t = ( int * ) malloc ( npart * sizeof ( int ) );

  rank = -1;

  for ( ; ; )
  {
    rank_old = rank;

    npart_sf_lex_successor ( n, npart, t, &rank );

    if ( rank <= rank_old )
    {
      break;
    }

    printf ( "  %4d", rank );
    for ( i = 0; i < npart; i++ )
    {
      printf ( "  %4d", t[i] );
    }
    printf ( "\n" );
  }
  return;
}
//****************************************************************************80

void test21 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST21 tests NPART_TABLE and PART_TABLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 July 2011
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int j;
  int maxn = 10;
  int maxpart = 5;
  int *p;
  int *p2;

  printf ( "\n" );
  printf ( "TEST21\n" );
  printf ( "  NPART_TABLE tabulates partitions\n" );
  printf ( "  of N with NPART parts;\n" );
  printf ( "  PART_TABLE tabulates partitions of N.\n" );

  p = npart_table ( maxn, maxpart );

  p2 = part_table ( maxn );

  printf ( "\n" );
  printf ( "    I P(I)  P(I,0) P(I,1) P(I,2) P(I,3) P(I,4) P(I,5)\n" );
  printf ( "\n" );

  for ( i = 0; i <= maxn; i++ )
  {
    printf ( "  %2d  %4d", i, p2[i] );
    for ( j = 0; j <= maxpart; j++ )
    {
      printf ( "  %4d", p[i+j*(maxn+1)] );
    }
    printf ( "\n" );
  }

  free ( p );
  free ( p2 );

  return;
}
//****************************************************************************80

void test22 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST22 tests PART_ENUM and PART_SUCCESSOR.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 July 2011
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n = 8;
  int npart;
  int npartitions;
  int rank;
  int rank_old;
  int *t;

  printf ( "\n" );
  printf ( "TEST22\n" );
  printf ( "  PART_SUCCESSOR produces partitions of N,\n" );
  printf ( "  PART_ENUM enumerates.\n" );
  printf ( "\n" );
  printf ( "  Partitions of N = %d\n", n );
//
//  Enumerate.
//
  npartitions = part_enum ( n );

  printf ( "\n" );
  printf ( "  For N = %d\n", n );
  printf ( "  the number of partitions is %d\n", npartitions );
  printf ( "\n" );
//
//  List.
//
  t = ( int * ) malloc ( n * sizeof ( int ) );

  rank = -1;

  for ( ; ; )
  {
    rank_old = rank;

    part_successor ( n, &npart, t, &rank );

    if ( rank <= rank_old )
    {
      break;
    }

    printf ( "  %4d", rank );
    for ( i = 0; i < npart; i++ )
    {
      printf ( "  %4d", t[i] );
    }
    printf ( "\n" );
  }

  free ( t );

  return;
}
//****************************************************************************80

void test23 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST23 tests PART_SUCCESSOR and PART_SF_CONJUGATE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 July 2011
//
//  Author:
//
//    John Burkardt
//
{
  int *b;
  int i;
  int n = 8;
  int npart;
  int npartb;
  int rank;
  int rank_old;
  int *t;

  printf ( "\n" );
  printf ( "TEST23\n" );
  printf ( "  PART_SUCCESSOR produces partitions of N,\n" );
  printf ( "  PART_SF_CONJUGATE produces the conjugate of a partition.\n" );
  printf ( "\n" );
  printf ( "  Partitions of N = %d\n", n );
//
//  List.
//
  t = ( int * ) malloc ( n * sizeof ( int ) );

  rank = -1;

  for ( ; ; )
  {
    rank_old = rank;

    part_successor ( n, &npart, t, &rank );

    if ( rank <= rank_old )
    {
      break;
    }

    printf ( "  %4d", rank );
    for ( i = 0; i < npart; i++ )
    {
      printf ( "  %4d", t[i] );
    }
    printf ( "\n" );

    b = part_sf_conjugate ( n, npart, t, &npartb );
    i4vec_transpose_print ( npartb, b, "  Con:" );
    free ( b );
  }

  free ( t );

  return;
}
//****************************************************************************80

void test24 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST24 tests PART_SF_MAJORIZE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 July 2011
//
//  Author:
//
//    John Burkardt
//
{
# define N 8

  int a[N] = { 2, 2, 2, 1, 1, 0, 0, 0 };
  int b[N] = { 3, 1, 1, 1, 1, 1, 0, 0 };
  int c[N] = { 2, 2, 1, 1, 1, 1, 0, 0 };
  int n = N;
  int nparta = 5;
  int npartb = 6;
  int npartc = 6;
  int result;

  printf ( "\n" );
  printf ( "TEST24\n" );
  printf ( "  PART_SF_MAJORIZE determines if one partition\n" );
  printf ( "  majorizes another.\n" );
  printf ( "\n" );
  printf ( "  Partitions of N = %d\n", n );
  printf ( "\n" );
  i4vec_transpose_print ( nparta, a, "  A: " );
  i4vec_transpose_print ( npartb, b, "  B: " );
  i4vec_transpose_print ( npartc, c, "  C: " );

  result = part_sf_majorize ( n, nparta, a, npartb, b );
  printf ( "\n" );
  printf ( "  A compare B: %d\n", result );
  result = part_sf_majorize ( n, npartb, b, npartc, c );
  printf ( "  B compare C: %d\n", result );
  result = part_sf_majorize ( n, npartc, c, nparta, a );
  printf ( "  C compare A: %d\n", result );
  result = part_sf_majorize ( n, npartc, c, npartc, c );
  printf ( "  C compare C: %d\n", result );

  return;
# undef N
}
//****************************************************************************80

void test25 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST25 tests PARTITION_GREEDY.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 July 2011
//
//  Author:
//
//    John Burkardt
//
{
# define N 10

  int a1[N] = { 2, 10, 3, 8, 5, 7, 9, 5, 3, 2 };
  int a2[N] = { 771, 121, 281, 854, 885, 734, 486, 1003, 83, 62 };
  int i;
  int *indx;
  int n = N;
  int sums[2];

  printf ( "\n" );
  printf ( "TEST25\n" );
  printf ( "  PARTITION_GREEDY partitions an integer vector into\n" );
  printf ( "  two subsets with nearly equal sum.\n" );
  printf ( "\n" );

  indx = partition_greedy ( n, a1 );

  printf ( "\n" );
  printf ( "\n" );
  printf ( "Data set #1 partitioned:\n" );
  printf ( "\n" );
  sums[0] = 0;
  sums[1] = 0;

  for ( i = 0; i < n; i++ )
  {
    if ( indx[i] == 1 ) 
    {
      sums[0] = sums[0] + a1[i];
      printf ( "  %4d\n", a1[i] );
    }
    else
    {
      sums[1] = sums[1] + a1[i];
      printf ( "        %4d\n", a1[i] );
    }
  }

  printf ( "\n" );
  printf ( "Sums:\n" );
  printf ( "\n" );
  printf ( "  %4d  %4d\n", sums[0], sums[1] );

  free ( indx );

  indx = partition_greedy ( n, a2 );

  printf ( "\n" );
  printf ( "\n" );
  printf ( "Data set #2 partitioned:\n" );
  printf ( "\n" );

  sums[0] = 0;
  sums[1] = 0;

  for ( i = 0; i < n; i++ )
  {
    if ( indx[i] == 1 ) 
    {
      sums[0] = sums[0] + a2[i];
      printf ( "  %4d\n", a2[i] );
    }
    else
    {
      sums[1] = sums[1] + a2[i];
      printf ( "        %4d\n", a2[i] );
    }
  }

  printf ( "\n" );
  printf ( "Sums:\n" );
  printf ( "\n" );
  printf ( "  %4d  %4d\n", sums[0], sums[1] );

  free ( indx );

  return;
# undef N
}
//****************************************************************************80

void test26 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST26 tests PARTN_ENUM, PARTN_SUCCESSOR and PART_SF_CONJUGATE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 July 2011
//
//  Author:
//
//    John Burkardt
//
{
  int *b;
  int i;
  int n = 11;
  int nmax;
  int npart;
  int npart2;
  int npartitions;
  int rank;
  int rank_old;
  int *t;

  printf ( "\n" );
  printf ( "TEST26\n" );
  printf ( "  Partitions of N with maximum element NMAX:\n" );
  printf ( "\n" );
  printf ( "  PARTN_SUCCESSOR lists;\n" );
  printf ( "  PARTN_ENUM enumerates.\n" );

  nmax = 4;
//
//  Enumerate.
//
  npartitions = partn_enum ( n, nmax );

  printf ( "\n" );
  printf ( "  For N = %d\n", n );
  printf ( "  and NMAX = %d\n", nmax );
  printf ( "  the number of partitions is %d\n", npartitions );
  printf ( "\n" );
//
//  List.
//
  t = ( int * ) malloc ( n * sizeof ( int ) );

  rank = -1;

  for ( ; ; )
  {
    rank_old = rank;

    partn_successor ( n, nmax, &npart, t, &rank );

    if ( rank <= rank_old )
    {
      break;
    }

    printf ( "  %4d", rank );
    for ( i = 0; i < npart; i++ )
    {
      printf ( "  %4d", t[i] );
    }
    printf ( "\n" );
  }
//
//  List conjugates.
//
  printf ( "\n" );
  printf ( "  Repeat, but list RSF conjugated partitions.\n" );
  printf ( "\n" );

  rank = -1;

  for ( ; ; )
  {
    rank_old = rank;

    partn_successor ( n, nmax, &npart, t, &rank );

    if ( rank <= rank_old )
    {
      break;
    }

    b = part_sf_conjugate ( n, npart, t, &npart2 );

    i4vec_reverse ( npart2, b );

    printf ( "  %4d", rank );
    for ( i = 0; i < npart2; i++ )
    {
      printf ( "  %4d", b[i] );
    }
    printf ( "\n" );
    free ( b );
  }

  free ( t );

  return;
}
//****************************************************************************80

void test27 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST27 tests PERM_INV and PERM_MUL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 July 2011
//
//  Author:
//
//    John Burkardt
//
{
  int n = 4;
  int nperm;
  int *p;
  int *q;
  int *r;
  int rank;

  printf ( "\n" );
  printf ( "TEST27\n" );
  printf ( "  Permutations of the integers:\n" );
  printf ( "\n" );
  printf ( "  PERM_INV computes an inverse permutation,\n" );
  printf ( "  PERM_MUL multiplies two permutations.\n" );
//
//  Enumerate.
//
  nperm = perm_enum ( n );
//
//  Unrank.
//
  rank = nperm / 2;

  p = perm_lex_unrank ( rank, n );

  perm_print ( n, p, "  The permutation P:" );
//
//  Invert.
//
  q = perm_inv ( n, p );

  perm_print ( n, q, "  The inverse permutation Q:" );
//
//  Multiply.
//
  r = perm_mul ( n, p, q );

  perm_print ( n, r, "  The product R = P * Q:" );

  free ( p );
  free ( q );
  free ( r );

  return;
}
//****************************************************************************80

void test28 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST28 tests PERM_ENUM, _LEX_RANK, _LEX_SUCCESSOR, _LEX_UNRANK.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 July 2011
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n = 4;
  int nperm;
  int *pi;
  int rank;
  int rank_old;

  printf ( "\n" );
  printf ( "TEST28\n" );
  printf ( "  Permutations of the integers,\n" );
  printf ( "  using the lexicographic ordering:\n" );
  printf ( "\n" );
  printf ( "  PERM_ENUM enumerates,\n" );
  printf ( "  PERM_LEX_RANK ranks,\n" );
  printf ( "  PERM_LEX_SUCCESSOR lists,\n" );
  printf ( "  PERM_LEX_UNRANK unranks.\n" );
//
//  Enumerate.
//
  nperm = perm_enum ( n );

  printf ( "\n" );
  printf ( "  For N = %d\n", n );
  printf ( "  the number of permutations is %d\n", nperm );
  printf ( "\n" );
//
//  List
//
  pi = ( int * ) malloc ( n * sizeof ( int ) );

  rank = -1;

  for ( ; ; )
  {
    rank_old = rank;

    perm_lex_successor ( n, pi, &rank );

    if ( rank <= rank_old )
    {
      break;
    }

    printf ( "  %4d", rank );
    for ( i = 0; i < n; i++ )
    {
      printf ( "  %4d", pi[i] );
    }
    printf ( "\n" );
  }
//
//  Unrank.
//
  rank = nperm / 2;

  free ( pi );
  pi = perm_lex_unrank ( rank, n );

  printf ( "\n" );
  printf ( "  The element of rank %d\n", rank );

  perm_print ( n, pi, " " );
//
//  Rank.
//
  rank = perm_lex_rank ( n, pi );

  perm_print ( n, pi, "  Element to be ranked:" );
  printf ( "\n" );
  printf ( "  Rank is computed as %d\n", rank );

  free ( pi );

  return;
}
//****************************************************************************80

void test29 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST29 tests PERM_TJ_ENUM, _TJ_RANK, _TJ_SUCCESSOR, _TJ_UNRANK.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 July 2011
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n = 4;
  int nperm;
  int *pi;
  int rank;
  int rank_old;

  printf ( "\n" );
  printf ( "TEST29\n" );
  printf ( "  Permutations of the integers\n" );
  printf ( "  using the Trotter-Johnson ordering:\n" );
  printf ( "\n" );
  printf ( "  PERM_ENUM enumerates,\n" );
  printf ( "  PERM_TJ_RANK ranks,\n" );
  printf ( "  PERM_TJ_SUCCESSOR lists,\n" );
  printf ( "  PERM_TJ_UNRANK unranks.\n" );
//
//  Enumerate.
//
  nperm = perm_enum ( n );

  printf ( "\n" );
  printf ( "  For N = %d\n", n );
  printf ( "  the number of permutations is %d\n", nperm );
  printf ( "\n" );
//
//  List
//
  pi = ( int * ) malloc ( n * sizeof ( int ) );

  rank = -1;

  for ( ; ; )
  {
    rank_old = rank;

    perm_tj_successor ( n, pi, &rank );

    if ( rank <= rank_old )
    {
      break;
    }

    printf ( "  %4d", rank );
    for ( i = 0; i < n; i++ )
    {
      printf ( "  %4d", pi[i] );
    }
    printf ( "\n" );
  }
//
//  Unrank.
//
  rank = nperm / 2;

  free ( pi );

  pi = perm_tj_unrank ( rank, n );

  printf ( "\n" );
  printf ( "  The element of rank %d\n", rank );

  perm_print ( n, pi, " " );
//
//  Rank.
//
  rank = perm_tj_rank ( n, pi );

  perm_print ( n, pi, "  The element to be ranked:" );
  printf ( "\n" );
  printf ( "  Rank is computed as %d\n", rank );

  free ( pi );

  return;
}
//****************************************************************************80

void test30 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST30 tests PRUEFER_ENUM, PRUEFER_RANK, PRUEFER_SUCCESSOR, PRUEFER_UNRANK.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 July 2011
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n = 4;
  int ncode;
  int *p;
  int rank;
  int rank_old;

  printf ( "\n" );
  printf ( "TEST30\n" );
  printf ( "  Pruefer codes:\n" );
  printf ( "\n" );
  printf ( "  PRUEFER_ENUM enumerates,\n" );
  printf ( "  PRUEFER_RANK ranks,\n" );
  printf ( "  PRUEFER_SUCCESSOR lists,\n" );
  printf ( "  PRUEFER_UNRANK unranks.\n" );
//
//  Enumerate.
//
  ncode = pruefer_enum ( n );

  printf ( "\n" );
  printf ( "  For N = %d\n", n );
  printf ( "  the number of Pruefer codes is %d\n", ncode );
  printf ( "\n" );
//
//  List
//
  p = ( int * ) malloc ( ( n - 2 ) * sizeof ( int ) );

  rank = -1;

  for ( ; ; )
  {
    rank_old = rank;

    pruefer_successor ( n, p, &rank );

    if ( rank <= rank_old )
    {
      break;
    }

    printf ( "  %4d", rank );
    for ( i = 0; i < n - 2; i++ )
    {
      printf ( "  %4d", p[i] );
    }
    printf ( "\n" );
  }
//
//  Unrank.
//
  rank = ncode / 2;

  free ( p );

  p = pruefer_unrank ( rank, n );

  printf ( "\n" );
  printf ( "  The element of rank %d\n", rank );
  printf ( "\n" );
  i4vec_transpose_print ( n - 2, p, " " );
//
//  Rank.
//
  rank = pruefer_rank ( n, p );

  i4vec_transpose_print ( n - 2, p, "  Element to be ranked:" );
  printf ( "\n" );
  printf ( "  Rank is computed as %d\n", rank );

  free ( p );

  return;
}
//****************************************************************************80

void test31 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST31 tests PRUEFER_TO_TREE and TREE_TO_PRUEFER.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 July 2011
//
//  Author:
//
//    John Burkardt
//
{
  int i4_hi;
  int i4_lo;
  int j;
  int n = 5;
  int *p;
  int pruefer_num;
  int rank;
  int seed = 123456789;
  int *t;
  int test;
  int test_num = 5;

  printf ( "\n" );
  printf ( "TEST31\n" );
  printf ( "  PRUEFER_TO_TREE converts a Pruefer code to a tree;\n" );
  printf ( "  TREE_TO_PRUEFER converts a tree to a Pruefer code.\n" );

  pruefer_num = pruefer_enum ( n );

  i4_lo = 0;
  i4_hi = pruefer_num - 1;

  for ( test = 1; test <= test_num; test++ )
  {
//
//  Pick a "random" Pruefer code.
//
    rank = i4_uniform ( i4_lo, i4_hi, &seed );

    p = pruefer_unrank ( rank, n );

    printf ( "\n" );
    printf ( "  Random Pruefer code of rank %d\n", rank );
    i4vec_transpose_print ( n - 2, p, " " );
//
//  Convert the Pruefer code to a tree.
//
    t = pruefer_to_tree_new ( n, p );

    printf ( "\n" );
    printf ( "  Edge list for the corresponding tree:\n" );
    printf ( "\n" );
    for ( j = 0; j < n - 1; j++ )
    {
      printf ( "  %2d  %4d  %4d\n", j, t[0+j*2], t[1+j*2] );
    }
//
//  Convert the tree to a Pruefer code.
//
    free ( p );

    p = tree_to_pruefer ( n, t );

    printf ( "\n" );
    i4vec_transpose_print ( n - 2, p, "  Pruefer code:" );

    free ( p );
    free ( t );
  }
  return;
}
//****************************************************************************80

void test32 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST32 tests QUEENS and BACKTRACK.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 July 2011
//
//  Author:
//
//    John Burkardt
//
{
# define N 8

  int iarray[N];
  int indx;
  int istack[N*N];
  int k;
  int n = N;
  int maxstack = N * N;
  int nstack = 0;

  printf ( "\n" );
  printf ( "TEST32\n" );
  printf ( "  QUEENS produces nonattacking queens\n" );
  printf ( "  on a chessboard.\n" );
  printf ( "  BACKTRACK supervises a backtrack search.\n" );
  printf ( "\n" );

  indx = 0;

  for ( ; ; )
  {
    backtrack ( n, iarray, &indx, &k, &nstack, istack, maxstack );

    if ( indx == 1 )
    {
      i4vec_transpose_print ( n, iarray, " " );
    }
    else if ( indx == 2 )
    {
      queens ( n, iarray, k, &nstack, istack, maxstack );
    }
    else
    {
      break;
    }
  }

  return;
# undef N
}
//****************************************************************************80

void test33 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST33 tests RGF_G_TABLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 July 2011
//
//  Author:
//
//    John Burkardt
//
{
  int *d;
  int i;
  int j;
  int m = 6;

  printf ( "\n" );
  printf ( "TEST33\n" );
  printf ( "  RGF_G_TABLE tabulates generalized restricted\n" );
  printf ( "  growth functions.\n" );
  printf ( "\n" );

  d = rgf_g_table ( m );

  for ( i = 0; i <= m; i++ )
  {
    for ( j = 0; j <= m - i; j++ )
    {
      printf ( "  %4d", d[i+j*(m+1)] );
    }
    printf ( "\n" );
  }

  free ( d );

  return;
}
//****************************************************************************80

void test34 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST34 tests RGF_ENUM, RGF_RANK, RGF_SUCCESSOR, RGF_UNRANK.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 July 2011
//
//  Author:
//
//    John Burkardt
//
{
  int *f;
  int i;
  int m = 4;
  int nrgf;
  int rank;
  int rank_old;

  printf ( "\n" );
  printf ( "TEST34\n" );
  printf ( "  Restricted growth functions:\n" );
  printf ( "\n" );
  printf ( "  RGF_ENUM enumerates,\n" );
  printf ( "  RGF_RANK ranks,\n" );
  printf ( "  RGF_SUCCESSOR lists;\n" );
  printf ( "  RGF_UNRANK unranks.\n" );
//
//  Enumerate.
//
  nrgf = rgf_enum ( m );

  printf ( "\n" );
  printf ( "  For M = %d\n", m );
  printf ( "  the number of RGF's is %d\n", nrgf );
  printf ( "\n" );
//
//  List.
//
  f = ( int * ) malloc ( m * sizeof ( int ) );

  rank = -1;

  for ( ; ; )
  {
    rank_old = rank;

    rgf_successor ( m, f, &rank );

    if ( rank <= rank_old )
    {
      break;
    }

    printf ( "  %4d", rank );
    for ( i = 0; i < m; i++ )
    {
      printf ( "  %4d", f[i] );
    }
    printf ( "\n" );
  }
//
//  Unrank.
//
  rank = nrgf / 2;

  free ( f );
  f = rgf_unrank ( rank, m );

  printf ( "\n" );
  printf ( "  The element of rank %d\n", rank );
  printf ( "\n" );
  i4vec_transpose_print ( m, f, " " );
//
//  Rank.
//
  rank = rgf_rank ( m, f );

  i4vec_transpose_print ( m, f, "  Element to be ranked:" );
  printf ( "\n" );
  printf ( "  Rank is computed as %d\n", rank );

  free ( f );

  return;
}
//****************************************************************************80

void test35 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST35 tests RGF_TO_SETPART and SETPART_TO_RGF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 July 2011
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int j;
  int jlo;
  int *f;
  int *index;
  int m = 8;
  int nsub;
  int rank;
  int *s;

  printf ( "\n" );
  printf ( "TEST35\n" );
  printf ( "  RGF_TO_SETPART converts a balanced\n" );
  printf ( "  sequence to a restricted growth function;\n" );
  printf ( "  SETPART_TO_RGF converts a restricted growth\n" );
  printf ( "  function to a balanced sequence.\n" );
//
//  Choose a "random" RGF.
//
  rank = 7;
  f = rgf_unrank ( rank, m );

  printf ( "\n" );
  printf ( "  Random restricted growth function:\n" );
  printf ( "\n" );
  i4vec_transpose_print ( m, f, " " );
//
//  Convert the RGF to a set partition.
//
  s = ( int * ) malloc ( m * sizeof ( int ) );
  index = ( int * ) malloc ( m * sizeof ( int ) );

  rgf_to_setpart ( m, f, &nsub, s, index );

  printf ( "\n" );
  printf ( "  Corresponding set partition\n" );
  printf ( "\n" );
  jlo = 1;
  for ( i = 1; i <= nsub; i++ )
  {
    for ( j = jlo; j <= index[i-1]; j++ )
    {
      printf ( "  %4d", s[j-1] );
    }
    printf ( "\n" );
    jlo = index[i-1] + 1;
  }
//
//  Convert the set partition back to an RGF.
//
  free ( f );

  f = setpart_to_rgf ( m, nsub, s, index );

  i4vec_transpose_print ( m, f, "  Corresponding RGF:" );

  free ( f );
  free ( index );
  free ( s );

  return;
}
//****************************************************************************80

void test36 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST36 tests SETPART_ENUM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 July 2011
//
//  Author:
//
//    John Burkardt
//
{
  int n;
  int npart;

  printf ( "\n" );
  printf ( "TEST36\n" );
  printf ( "  Set partitions:\n" );
  printf ( "\n" );
  printf ( "  SETPART_ENUM enumerates.\n" );
  printf ( "\n" );
//
//  Enumerate.
//
  for ( n = 1; n <= 6; n++ )
  {
    npart = setpart_enum ( n );
    printf ( "  %4d  %4d\n", n, npart );
  }

  return;
}
//****************************************************************************80

void test37 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST37 tests STIRLING_NUMBERS1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 July 2011
//
//  Author:
//
//    John Burkardt
//
{
  int maxm = 6;
  int maxn = 6;
  int *s;

  printf ( "\n" );
  printf ( "TEST37\n" );
  printf ( "  STIRLING_NUMBERS1 computes a table of Stirling\n" );
  printf ( "  numbers of the first kind.\n" );

  s = stirling_numbers1 ( maxm, maxn );

  i4mat_print ( maxm + 1, maxn + 1, s, "  Stirling number of first kind" ); 

  free ( s );

  return;
}
//****************************************************************************80

void test38 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST38 tests STIRLING_NUMBERS2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 July 2011
//
//  Author:
//
//    John Burkardt
//
{
  int maxm = 6;
  int maxn = 6;
  int *s;

  printf ( "\n" );
  printf ( "TEST38\n" );
  printf ( "  STIRLING_NUMBERS2 computes a table of Stirling\n" );
  printf ( "  numbers of the second kind.\n" );

  s = stirling_numbers2 ( maxm, maxn );

  i4mat_print ( maxm + 1, maxn + 1, s, "  Stirling number of second kind" ); 

  free ( s );

  return;
}
//****************************************************************************80

void test39 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST39 tests SUBSET_COLEX_RANK, _COLEX_SUCCESSOR, _COLEX_UNRANK, _ENUM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 July 2011
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n = 5;
  int nsub;
  int rank;
  int rank_old;
  int *t;

  printf ( "\n" );
  printf ( "TEST39\n" );
  printf ( "  All subsets of a set,\n" );
  printf ( "  using the colexicographic ordering:\n" );
  printf ( "\n" );
  printf ( "  SUBSET_COLEX_RANK ranks,\n" );
  printf ( "  SUBSET_COLEX_SUCCESSOR lists,\n" );
  printf ( "  SUBSET_COLEX_UNRANK unranks.\n" );
  printf ( "  SUBSET_ENUM enumerates.\n" );
//
//  Enumerate.
//
  nsub = subset_enum ( n );

  printf ( "\n" );
  printf ( "  For N = %d\n", n );
  printf ( "  the number of subsets is %d\n", nsub );
  printf ( "\n" );
//
//  List
//
  t = ( int * ) malloc ( n * sizeof ( int ) );

  rank = -1;

  for ( ; ; )
  {
    rank_old = rank;

    subset_colex_successor ( n, t, &rank );

    if ( rank <= rank_old )
    {
      break;
    }

    printf ( "  %4d", rank );
    for ( i = 0; i < n; i++ )
    {
      printf ( "  %4d", t[i] );
    }
    printf ( "\n" );
  }
//
//  Unrank.
//
  rank = nsub / 3;

  free ( t );
  t = subset_colex_unrank ( rank, n );

  printf ( "\n" );
  printf ( "  The element of rank %d\n", rank );
  printf ( "\n" );
  i4vec_transpose_print ( n, t, " " );
//
//  Rank.
//
  rank = subset_colex_rank ( n, t );

  i4vec_transpose_print ( n, t, "  Element to be ranked:" );
  printf ( "\n" );
  printf ( "  Rank is computed as %d\n", rank );

  free ( t );

  return;
}
//****************************************************************************80

void test40 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST40 tests SUBSET_ENUM, _LEX_RANK, _LEX_SUCCESSOR, _LEX_UNRANK.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 July 2011
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n = 5;
  int nsub;
  int rank;
  int rank_old;
  int *t;

  printf ( "\n" );
  printf ( "TEST40\n" );
  printf ( "  All subsets of a set,\n" );
  printf ( "  using the lexicographic ordering:\n" );
  printf ( "\n" );
  printf ( "  SUBSET_ENUM enumerates,\n" );
  printf ( "  SUBSET_LEX_RANK ranks,\n" );
  printf ( "  SUBSET_LEX_SUCCESSOR lists,\n" );
  printf ( "  SUBSET_LEX_UNRANK unranks.\n" );
//
//  Enumerate.
//
  nsub = subset_enum ( n );

  printf ( "\n" );
  printf ( "  For N = %d\n", n );
  printf ( "  the number of subsets is %d\n", nsub );
  printf ( "\n" );
//
//  List
//
  t = ( int * ) malloc ( n * sizeof ( int ) );

  rank = -1;

  for ( ; ; )
  {
    rank_old = rank;

    subset_lex_successor ( n, t, &rank );

    if ( rank <= rank_old )
    {
      break;
    }

    printf ( "  %4d", rank );
    for ( i = 0; i < n; i++ )
    {
      printf ( "  %4d", t[i] );
    }
    printf ( "\n" );
  }
//
//  Unrank.
//
  rank = nsub / 3;

  free ( t );

  t = subset_lex_unrank ( rank, n );

  printf ( "\n" );
  printf ( "  The element of rank %d\n", rank );
  printf ( "\n" );
  i4vec_transpose_print ( n, t, " " );
//
//  Rank.
//
  rank = subset_lex_rank ( n, t );

  i4vec_transpose_print ( n, t, "  Element to be ranked:" );
  printf ( "\n" );
  printf ( "  Rank is computed as %d\n", rank );

  free ( t );

  return;
}
//****************************************************************************80

void test41 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST41 tests SUBSETSUM_SWAP.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 July 2011
//
//  Author:
//
//    John Burkardt
//
{
# define N 7

  int a[N] = { 12, 8, 11, 30, 8, 3, 7 };
  int i;
  int index[N];
  int n = N;
  int sum_achieved;
  int sum_desired = 17;

  printf ( "\n" );
  printf ( "TEST41\n" );
  printf ( "  SUBSETSUM_SWAP seeks a solution of the subset\n" );
  printf ( "  sum problem using pair swapping.\n" );
  printf ( "\n" );
  printf ( "  The desired sum is %d\n", sum_desired );

  sum_achieved = subsetsum_swap ( n, a, sum_desired, index );

  printf ( "\n" );
  printf ( "    A(I), INDEX(I)\n" );
  printf ( "\n" );

  for ( i = 0; i < n; i++ )
  {
    printf ( "  %5d  %5d\n", a[i], index[i] );
  }

  printf ( "\n" );
  printf ( "  The achieved sum is %d\n", sum_achieved );

  return;
# undef N
}
//****************************************************************************80

void test42 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST42 tests TREE_ENUM, TREE_RANK, TREE_SUCCESSOR, TREE_UNRANK.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 July 2011
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int j;
  int n = 4;
  int rank;
  int rank_old;
  int *t;
  int tree_num;

  printf ( "\n" );
  printf ( "TEST42\n" );
  printf ( "  Trees:\n" );
  printf ( "\n" );
  printf ( "  TREE_ENUM enumerates,\n" );
  printf ( "  TREE_RANK ranks,\n" );
  printf ( "  TREE_SUCCESSOR lists,\n" );
  printf ( "  TREE_UNRANK unranks.\n" );
//
//  Enumerate.
//
  tree_num = tree_enum ( n );

  printf ( "\n" );
  printf ( "  For N = %d\n", n );
  printf ( "  the number of trees is %d\n", tree_num );
  printf ( "\n" );
//
//  List
//
  t = ( int * ) malloc ( 2 * ( n - 1 ) * sizeof ( int ) );

  rank = -1;

  for ( ; ; )
  {
    rank_old = rank;

    tree_successor ( n, t, &rank );

    if ( rank <= rank_old )
    {
      break;
    }
    printf ( "  %4d", rank );
    for ( j = 0; j < n - 1; j++ )
    {
      printf ( "  %4d", t[0+j*2] );
    }
    printf ( "\n" );
    printf ( "      " );
    for ( j = 0; j < n - 1; j++ )
    {
      printf ( "  %4d", t[1+j*2] );
    }
    printf ( "\n" );
  }
//
//  Unrank.
//
  rank = tree_num / 2;

  free ( t );

  t = tree_unrank ( rank, n );

  printf ( "\n" );
  printf ( "  The element of rank %d\n", rank );

  i4mat_print ( 2, n - 1, t, " " );
//
//  Rank.
//
  rank = tree_rank ( n, t );

  i4mat_print ( 2, n - 1, t, "  Element to be ranked:" );
  printf ( "\n" );
  printf ( "  Rank is computed as %d\n", rank );

  free ( t );

  return;
}
