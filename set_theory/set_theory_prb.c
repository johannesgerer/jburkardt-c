# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <stdbool.h>

# include "set_theory.h"

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );
void test04 ( );
void test05 ( );
void test06 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for SET_THEORY_PRB.

  Discussion:

    SET_THEORY_PRB tests the SET_THEORY library.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    20 September 2011

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "SET_THEORY_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the SET_THEORY library.\n" );
  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
  test05 ( );
  test06 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "SET_THEORY_PRB\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void test01 ( )

/******************************************************************************/
/*
  Purpose:

    TEST01 tests the B4SET routines.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 September 2011

  Author:

    John Burkardt
*/
{
  int a;
  int a_num;
  int *a_numeric;
  int b;
  int b_num = 16;
  int b_numeric[16] = {
    3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 33, 36, 39, 42, 45, 48 };
  int c;
  int d;
  int e;
  int f;
  int g;
  int h;
  int i;
  int n = 32;
  int u;
  int *u_numeric;
  int w;
  int w_num = 5;
  int w_numeric[5] = { 1, 11, 21, 31, 41 };
  int x;
  int y;
  int y_num = 4;
  int y_numeric[4] = { 4, 5, 6, 7 };

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  Test the set theory functions\n" );
  printf ( "  with the B4SET representation of a set.\n" );
/*
  Define the universal set.
*/
  u_numeric = ( int * ) malloc ( n * sizeof ( int ) );

  for ( i = 0; i < n; i++ )
  {
    u_numeric[i] = i + 1;
  }
  u = i4vec_to_b4set ( n, u_numeric, n );
/*
  Define the set A by a numeric property.
*/
  a_numeric = ( int * ) malloc ( n * sizeof ( int ) );
  a_num = 0;
  for ( i = 1; i <= n; i++ )
  {
    if ( ( i % 5 ) == 0 )
    {
      a_numeric[a_num] = i;
      a_num = a_num + 1;
    }
  }
  a = i4vec_to_b4set ( a_num, a_numeric, n );
  b4set_transpose_print ( n, a, "  A: " );
/*
  Define the set by starting with a numeric list of entries.
*/
  b = i4vec_to_b4set ( b_num, b_numeric, n );
  b4set_transpose_print ( n, b, "  B: " );
/*
  C is the complement of B (with respect to the universal set).
*/
  c = b4set_complement ( n, b );
  b4set_transpose_print ( n, c, "  C = ~ B:" );
/*
  D is the intersection of A and B.
*/
  d = b4set_intersect ( n, a, b );
  b4set_transpose_print ( n, d, "  D = A intersect B:" );
/*
  E is the intersection of A and B.
*/
  e = b4set_union ( n, a, b );
  b4set_transpose_print ( n, e, "  E = A union B:" );
/*
  F is the symmetric difference of A and B.
*/
  f = b4set_xor ( n, a, b );
  b4set_transpose_print ( n, f, "  F = A xor B:" );
/*
  G is the complement of B with respect to A.
  H is the complement of A with respect to B.
*/
  g = b4set_complement_relative ( n, a, b );
  b4set_transpose_print ( n, g, "  G = A ~ B:" );

  h = b4set_complement_relative ( n, b, a );
  b4set_transpose_print ( n, h, "  H = B ~ A:" );
/*
  B4SET_IS_MEMBER checks if an element is in a set.
*/
  printf ( "\n" );
  printf ( "  B4SET_IS_MEMBER ( i, A ) reports whether i is a member of A\n" );
  printf ( "\n" );

  for ( i = 10; i <= 20; i++ )
  {
    if ( b4set_is_member ( n, i, a ) )
    {
      printf ( "  %d is a member of A.\n", i );
    }
    else
    {
      printf ( "  %d is not a member of A.\n", i );
    }
  }
/*
  B4SET_IS_SUBSET checks whether a set is a subset.
*/
  printf ( "\n" );
  printf ( "  B4SET_IS_SUBSET ( D, A ) reports whether D is a subset of A\n" );
  printf ( "\n" );

  d = b4set_intersect ( n, a, b );

  if ( b4set_is_subset ( n, d, a ) )
  {
    printf ( "  ( A intersect B ) is a subset of A.\n" );
  }
  else
  {
    printf ( "  ( A intersect B)  is not a subset of A.\n" );
  }
/*
  B4SET_INSERT adds an item to a set.
*/
  w = i4vec_to_b4set ( w_num, w_numeric, n );
  b4set_transpose_print ( n, w, "  W: " );

  x = 6;
  w = b4set_insert ( n, x, w );
  b4set_transpose_print ( n, w, "  W := W + 6:" );

  x = 31;
  w = b4set_delete ( n, x, w );
  b4set_transpose_print ( n, w, "  W := W - 31:" );

  y = i4vec_to_b4set ( y_num, y_numeric, n );
  w = b4set_union ( n, w, y );

  b4set_transpose_print ( n, w, "  W := W union [ 4, 5, 6, 7 ]:" );

  free ( a_numeric );
  free ( u_numeric );

  return;
}
/******************************************************************************/

void test02 ( )

/******************************************************************************/
/*
  Purpose:

    TEST02 tests B4SET_COLEX_RANK, _COLEX_SUCCESSOR, _COLEX_UNRANK, _ENUM.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 September 2011

  Author:

    John Burkardt
*/
{
  int n = 5;
  int nsub;
  int rank;
  int rank_old;
  char s[80];
  int t;
  char title[80];

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  All subsets of a set,\n" );
  printf ( "  using the colexicographic ordering\n" );
  printf ( "  with the B4SET representation of a set.\n" );
  printf ( "\n" );
  printf ( "  B4SET_COLEX_RANK ranks,\n" );
  printf ( "  B4SET_COLEX_SUCCESSOR lists,\n" );
  printf ( "  B4SET_COLEX_UNRANK unranks.\n" );
  printf ( "  B4SET_ENUM enumerates.\n" );
/*
  Enumerate.
*/
  nsub = b4set_enum ( n );

  printf ( "\n" );
  printf ( "  For N = %d\n", n );
  printf ( "  the number of subsets is %d\n", nsub );
  printf ( "\n" );
/*
  List
*/
  rank = -1;

  for ( ; ; )
  {
    rank_old = rank;

    b4set_colex_successor ( n, &t, &rank );

    if ( rank <= rank_old )
    {
      break;
    }
    sprintf ( title, "Rank: %d", rank );
    b4set_transpose_print ( n, t, title );
  }
/*
  Unrank.
*/
  rank = nsub / 3;

  t = b4set_colex_unrank ( rank, n );

  sprintf ( title, "The element of rank ", rank );
  b4set_transpose_print ( n, t, title );
/*
  Rank.
*/
  rank = b4set_colex_rank ( n, t );

  printf ( "\n" );
  printf ( "  The rank of this element is computed as %d\n", rank );

  return;
}
/******************************************************************************/

void test03 ( )

/******************************************************************************/
/*
  Purpose:

    TEST03 tests B4SET_LEX_RANK, _LEX_SUCCESSOR, _LEX_UNRANK, _ENUM.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 September 2011

  Author:

    John Burkardt
*/
{
  int n = 5;
  int nsub;
  int rank;
  int rank_old;
  char s[80];
  int t;
  char title[80];

  printf ( "\n" );
  printf ( "TEST03\n" );
  printf ( "  All subsets of a set,\n" );
  printf ( "  using the lexicographic ordering,\n" );
  printf ( "  with the B4SET representation of a set.\n" );
  printf ( "\n" );
  printf ( "  B4SET_LEX_RANK ranks,\n" );
  printf ( "  B4SET_LEX_SUCCESSOR lists,\n" );
  printf ( "  B4SET_LEX_UNRANK unranks.\n" );
  printf ( "  B4SET_ENUM enumerates.\n" );
/*
  Enumerate.
*/
  nsub = b4set_enum ( n );

  printf ( "\n" );
  printf ( "  For N = %d\n", n );
  printf ( "  the number of subsets is %d\n", nsub );
  printf ( "\n" );
/*
  List
*/
  rank = -1;

  for ( ; ; )
  {
    rank_old = rank;

    b4set_lex_successor ( n, &t, &rank );

    if ( rank <= rank_old )
    {
      break;
    }
    sprintf ( title, "Rank: %d", rank );
    b4set_transpose_print ( n, t, title );
  }
/*
  Unrank.
*/
  rank = nsub / 3;

  t = b4set_lex_unrank ( rank, n );

  sprintf ( title, "The element of rank ", rank );
  b4set_transpose_print ( n, t, title );
/*
  Rank.
*/
  rank = b4set_lex_rank ( n, t );

  printf ( "\n" );
  printf ( "  The rank of this element is computed as %d\n", rank );

  return;
}
/******************************************************************************/

void test04 ( )

/******************************************************************************/
/*
  Purpose:

    TEST04 tests the LSET routines.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 August 2011

  Author:

    John Burkardt
*/
{
  bool *a;
  bool *b;
  int b_num = 16;
  int b_numeric[16] = { 
    3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 33, 36, 39, 42, 45, 48 };
  bool *c;
  bool *d;
  bool *e;
  bool *f;
  bool *g;
  bool *h;
  int i;
  int n = 50;
  bool *u;
  int *u_numeric;
  bool *w;
  int w_num = 5;
  int w_numeric[5] = { 1, 11, 21, 31, 41 };
  int x;
  bool *y;
  int y_num = 4;
  int y_numeric[4] = { 16, 26, 36, 46 };

  printf ( "\n" );
  printf ( "TEST04\n" );
  printf ( "  Test the set theory functions\n" );
  printf ( "  with the LSET representation of a set.\n" );
/*
  Define the universal set.
*/
  u_numeric = ( int * ) malloc ( n * sizeof ( int ) );
  u = ( bool * ) malloc ( n * sizeof ( bool ) );

  for ( i = 0; i < n; i++ )
  {
    u_numeric[i] = i + 1;
    u[i] = true;
  }
/*
  Define the set A by a numeric property.
*/
  a = ( bool * ) malloc ( n * sizeof ( bool ) );
  for ( i = 0; i < n; i++ )
  {
    a[i] = ( ( u_numeric[i] % 5 ) == 0 );
  }
  lset_transpose_print ( n, a, "  A: " );
/*
  Define the set by starting with a numeric list of entries.
*/
  b = i4vec_to_lset ( b_num, b_numeric, n );
  lset_transpose_print ( n, b, "  B: " );
/*
  C is the complement of B (with respect to the universal set).
*/
  c = lset_complement ( n, b );
  lset_transpose_print ( n, c, "  C = ~ B:" );
/*
  D is the intersection of A and B.
*/
  d = lset_intersect ( n, a, b );
  lset_transpose_print ( n, d, "  D = A intersect B:" );
/*
  E is the intersection of A and B.
*/
  e = lset_union ( n, a, b );
  lset_transpose_print ( n, e, "  E = A union B:" );
/*
  F is the symmetric difference of A and B.
*/
  f = lset_xor ( n, a, b );
  lset_transpose_print ( n, f, "  F = A xor B:" );
/*
  G is the complement of B with respect to A.
  H is the complement of A with respect to B.
*/
  g = lset_complement_relative ( n, a, b );
  lset_transpose_print ( n, g, "  G = A ~ B:" );

  h = lset_complement_relative ( n, b, a );
  lset_transpose_print ( n, h, "  H = B ~ A:" );
/*
  LSET_IS_MEMBER checks if an element is in a set.
*/
  printf ( "\n" );
  printf ( "  LSET_IS_MEMBER ( i, A ) reports whether i is a member of A\n" );
  printf ( "\n" );

  for ( i = 10; i <= 20; i++ )
  {
    if ( lset_is_member ( n, i, a ) )
    {
      printf ( "  %d is a member of A.\n", i );
    }
    else
    {
      printf ( "  %d is not a member of A.\n", i );
    }
  }
/*
  LSET_IS_SUBSET checks whether a set is a subset.
*/
  printf ( "\n" );
  printf ( "  LSET_IS_SUBSET ( D, A ) reports whether D is a subset of A\n" );
  printf ( "\n" );

  free ( d );
  d = lset_intersect ( n, a, b );

  if ( lset_is_subset ( n, d, a ) )
  {
    printf ( "  ( A intersect B ) is a subset of A.\n" );
  }
  else
  {
    printf ( "  ( A intersect B)  is not a subset of A.\n" );
  }
/*
  LSET_INSERT adds an item to a set.
*/
  w = i4vec_to_lset ( w_num, w_numeric, n );
  lset_transpose_print ( n, w, "  W:" );

  x = 6;
  lset_insert ( n, x, w );
  lset_transpose_print ( n, w, "  W := W + 6:" );

  x = 31;
  lset_delete ( n, x, w );
  lset_transpose_print ( n, w, "  W := W - 31:" );
  free ( w );

  y = i4vec_to_lset ( y_num, y_numeric, n );
  w = lset_union ( n, w, y );

  lset_transpose_print ( n, w, "  W := W union [16, 26, 36, 46]:" );

  free ( a );
  free ( b );
  free ( c );
  free ( d );
  free ( e );
  free ( f );
  free ( g );
  free ( h );
  free ( u );
  free ( u_numeric );
  free ( w );
  free ( y );

  return;
}
/******************************************************************************/

void test05 ( )

/******************************************************************************/
/*
  Purpose:

    TEST05 tests LSET_COLEX_RANK, _COLEX_SUCCESSOR, _COLEX_UNRANK, _ENUM.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 August 2011

  Author:

    John Burkardt
*/
{
  int i;
  int n = 5;
  int nsub;
  int rank;
  int rank_old;
  char s[80];
  bool *t;
  char title[80];

  printf ( "\n" );
  printf ( "TEST05\n" );
  printf ( "  All subsets of a set,\n" );
  printf ( "  using the colexicographic ordering\n" );
  printf ( "  with the LSET representation of a set.\n" );
  printf ( "\n" );
  printf ( "  LSET_COLEX_RANK ranks,\n" );
  printf ( "  LSET_COLEX_SUCCESSOR lists,\n" );
  printf ( "  LSET_COLEX_UNRANK unranks.\n" );
  printf ( "  LSET_ENUM enumerates.\n" );

  t = ( bool * ) malloc ( n * sizeof ( bool ) );
/*
  Enumerate.
*/
  nsub = lset_enum ( n );

  printf ( "\n" );
  printf ( "  For N = %d\n", n );
  printf ( "  the number of subsets is %d\n", nsub );
  printf ( "\n" );
/*
  List
*/
  rank = -1;

  while ( true )
  {
    rank_old = rank;

    lset_colex_successor ( n, t, &rank );

    if ( rank <= rank_old )
    {
      break;
    }
    sprintf ( title, "Rank: %d", rank );
    lset_transpose_print ( n, t, title );
  }
/*
  Unrank.
*/
  rank = nsub / 3;

  free ( t );

  t = lset_colex_unrank ( rank, n );

  sprintf ( title, "  The element of rank ", rank );
  lset_transpose_print ( n, t, title );
/*
  Rank.
*/
  rank = lset_colex_rank ( n, t );

  printf ( "\n" );
  printf ( "  The rank of this element is computed as %d\n", rank );

  free ( t );

  return;
}
/******************************************************************************/

void test06 ( )

/******************************************************************************/
/*
  Purpose:

    TEST06 tests LSET_LEX_RANK, _LEX_SUCCESSOR, _LEX_UNRANK, _ENUM.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 September 2011

  Author:

    John Burkardt
*/
{
  int n = 5;
  int nsub;
  int rank;
  int rank_old;
  char s[80];
  bool *t;
  char title[80];

  printf ( "\n" );
  printf ( "TEST06\n" );
  printf ( "  All subsets of a set,\n" );
  printf ( "  using the lexicographic ordering,\n" );
  printf ( "  with the LSET representation of a set.\n" );
  printf ( "\n" );
  printf ( "  LSET_LEX_RANK ranks,\n" );
  printf ( "  LSET_LEX_SUCCESSOR lists,\n" );
  printf ( "  LSET_LEX_UNRANK unranks.'\n" );
  printf ( "  LSET_ENUM enumerates.\n" );

  t = ( bool * ) malloc ( n * sizeof ( bool ) );
/*
  Enumerate.
*/
  nsub = lset_enum ( n );

  printf ( "\n" );
  printf ( "  For N = %d\n", n );
  printf ( "  the number of subsets is %d\n", nsub );
  printf ( "\n" );
/*
  List
*/
  rank = -1;

  while ( 1 )
  {
    rank_old = rank;

    lset_lex_successor ( n, t, &rank );

    if ( rank <= rank_old )
    {
      break;
    }
    sprintf ( title, "Rank: %d", rank );
    lset_transpose_print ( n, t, title );
  }
/*
  Unrank.
*/
  rank = nsub / 3;

  free ( t );
  t = lset_lex_unrank ( rank, n );

  sprintf ( title, "  The element of rank %d", rank );
  lset_transpose_print ( n, t, title );
/*
  Rank.
*/
  rank = lset_lex_rank ( n, t );

  printf ( "\n" );
  printf ( "  The rank of this element is computed as %d\n", rank );

  free ( t );

  return;
}
