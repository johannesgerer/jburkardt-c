# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <string.h>
# include <stdbool.h>

# include "set_theory.h"

/******************************************************************************/

int b4set_colex_rank ( int n, int t )

/******************************************************************************/
/*
  Purpose:

    B4SET_COLEX_RANK computes the colexicographic rank of a B4SET.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 September 2011

  Author:

    John Burkardt

  Reference:

    Donald Kreher, Douglas Simpson,
    Combinatorial Algorithms,
    CRC Press, 1998,
    ISBN: 0-8493-3988-X,
    LC: QA164.K73.

  Parameters:

    Input, int N, the order of the master set.

    Input, int T, the set.

    Output, int RANK, the rank of the set.
*/
{
  int i;
  int rank;

  rank = 0;

  for ( i = 0; i < n; i++ )
  {
    if ( i4_btest ( t, i ) )
    {
      rank = rank + i4_power ( 2, i );
    }
  }
  return rank;
}
/******************************************************************************/

void b4set_colex_successor ( int n, int *t, int *rank )

/******************************************************************************/
/*
  Purpose:

    B4SET_COLEX_SUCCESSOR computes the colexicographic successor of a B4SET.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 September 2011

  Author:

    John Burkardt

  Reference:

    Donald Kreher, Douglas Simpson,
    Combinatorial Algorithms,
    CRC Press, 1998,
    ISBN: 0-8493-3988-X,
    LC: QA164.K73.

  Parameters:

    Input, int N, the order of the master set.

    Input/output, int *T, describes a set.  
    On input, T describes a set.
    On output, T describes the next set in the ordering.
    If the input T was the last in the ordering, then the output T
    will be the first.

    Input/output, int *RANK, the rank.
    If RANK = -1 on input, then the routine understands that this is
    the first call, and that the user wishes the routine to supply
    the first element in the ordering, which has RANK = 0.
    In general, the input value of RANK is increased by 1 for output,
    unless the very last element of the ordering was input, in which
    case the output value of RANK is 0.
*/
{
  int i;
/*
  Return the first element.
*/
  if ( *rank == -1 )
  {
    *t = 0;
    *rank = 0;
    return;
  }
  for ( i = 0; i < n; i++ )
  {
    if ( ! i4_btest ( *t, i ) )
    {
      *t = i4_bset ( *t, i );
      rank = rank + 1;
      return;
    }
    else
    {
      *t = i4_bclr ( *t, i );
    }
  }
  *rank = 0;
  return;
}
/******************************************************************************/

int b4set_colex_unrank ( int rank, int n )

/******************************************************************************/
/*
  Purpose:

    B4SET_COLEX_UNRANK computes the B4SET of given colexicographic rank.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 September 2011

  Author:

    John Burkardt

  Reference:

    Donald Kreher, Douglas Simpson,
    Combinatorial Algorithms,
    CRC Press, 1998,
    ISBN: 0-8493-3988-X,
    LC: QA164.K73.

  Parameters:

    Input, int RANK, the rank of the set.

    Input, int N, the order of the master set.

    Output, int B4SET_COLEX_UNRANK, the set of the given rank.
*/
{
  int i;
  int rank_copy;
  int sub_num;
  int t;
/*
  Check.
*/
  if ( n < 1 )
  {
    fprintf ( stderr,  "\n" );
    fprintf ( stderr,  "B4SET_COLEX_UNRANK - Fatal error!\n" );
    fprintf ( stderr,  "  Input N is illegal.\n" );
    exit ( 1 );
  }

  sub_num = b4set_enum ( n );

  if ( rank < 0 || sub_num < rank )
  {
    fprintf ( stderr,  "\n" );
    fprintf ( stderr,  "B4SET_COLEX_UNRANK - Fatal error!\n" );
    fprintf ( stderr,  "  The input rank is illegal.\n" );
    exit ( 1 );
  }

  rank_copy = rank;

  t = 0;

  for ( i = 0; i < n; i++ )
  {
    if ( ( rank_copy % 2 ) == 1 )
    {
      t = i4_bset ( t, i );
    }
    else
    {
      t = i4_bclr ( t, i );
    }
    rank_copy = rank_copy / 2;
  }
  return t;
}
/******************************************************************************/

int b4set_complement ( int n, int a )

/******************************************************************************/
/*
  Purpose:

    B4SET_COMPLEMENT computes the complement of a B4SET.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 September 2011

  Author:

    John Burkardt

  Reference:

    Donald Kreher, Douglas Simpson,
    Combinatorial Algorithms,
    CRC Press, 1998,
    ISBN: 0-8493-3988-X,
    LC: QA164.K73.

  Parameters:

    Input, int N, the order of the master set.

    Input, int A, the set.

    Output, int B4SET_COMPLEMENT, the complement of A.
*/
{
  int b;
  int i;

  b = 0;

  for ( i = 0; i < n; i++ )
  {
    if ( ! i4_btest ( a, i ) )
    {
      b = i4_bset ( b, i );
    }
  }
  return b;
}
/******************************************************************************/

int b4set_complement_relative ( int n, int a, int b )

/******************************************************************************/
/*
  Purpose:

    B4SET_COMPLEMENT_RELATIVE computes the relative complement of a B4SET.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 September 2011

  Author:

    John Burkardt

  Reference:

    Donald Kreher, Douglas Simpson,
    Combinatorial Algorithms,
    CRC Press, 1998,
    ISBN: 0-8493-3988-X,
    LC: QA164.K73.

  Parameters:

    Input, int N, the order of the master set.

    Input, int A, the set.

    Input, int B, the set with respect to which 
    the complement is taken.

    Output, int B4SET_COMPLEMENT, the complement of A with respect to B.
*/
{
  int c;
  int i;

  c = 0;

  for ( i = 0; i < n; i++ )
  {
    if (  i4_btest ( a, i ) && ! i4_btest ( b, i ) )
    {
      c = i4_bset ( c, i );
    }
  }
  return c;
}
/******************************************************************************/

int b4set_delete ( int n, int a, int t )

/******************************************************************************/
/*
  Purpose:

    B4SET_DELETE deletes an element from a B4SET.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 September 2011

  Author:

    John Burkardt

  Reference:

    Donald Kreher, Douglas Simpson,
    Combinatorial Algorithms,
    CRC Press, 1998,
    ISBN: 0-8493-3988-X,
    LC: QA164.K73.

  Parameters:

    Input, int N, the order of the master set.

    Input, int A, an item.

    Input, int T, a set.

    Output, int B4SET, the set with A deleted.
*/
{
  int value;

  if ( a < 1 || n < a )
  {
    fprintf ( stderr,  "\n" );
    fprintf ( stderr,  "B4SET_DELETE - Fatal error!\n" );
    fprintf ( stderr,  "  1 <= A <= N fails.\n" );
    exit ( 1 );
  }

  value = i4_bclr ( t, a - 1 );

  return value;
}
/******************************************************************************/

int b4set_distance ( int n, int t1, int t2 )

/******************************************************************************/
/*
  Purpose:

    B4SET_DISTANCE computes the Hamming distance between two B4SET's.

  Discussion:

    The sets T1 and T2 are assumed to be subsets of a set of N elements.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 September 2011

  Author:

    John Burkardt

  Reference:

    Donald Kreher, Douglas Simpson,
    Combinatorial Algorithms,
    CRC Press, 1998,
    ISBN: 0-8493-3988-X,
    LC: QA164.K73.

  Parameters:

    Input, int N, the order of the master set.

    Input, int T1, T2, two sets.

    Output, int B4SET_DISTANCE, the Hamming distance between T1 and T2,
    defined as the number of elements of the master set which are
    in either T1 or T2 but not both.
*/
{
  int dist;
  int i;

  dist = 0;

  for ( i = 0; i < n; i++ )
  {
    if ( (   i4_btest ( t1, i ) && ! i4_btest ( t2, i ) ) ||
         ( ! i4_btest ( t1, i ) &&   i4_btest ( t2, i ) ) )
    {
      dist = dist + 1;
    }
  }
  return dist;
}
/******************************************************************************/

int b4set_enum ( int n )

/******************************************************************************/
/*
  Purpose:

    B4SET_ENUM enumerates the B4SET's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 September 2011

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of the master set.

    Output, int B4SET_ENUM, the number of distinct sets.
*/
{
  int value;

  value = i4_power ( 2, n );

  return value;
}
/******************************************************************************/

int b4set_index ( int n, int a, int t )

/******************************************************************************/
/*
  Purpose:

    B4SET_INDEX returns the index of an element of a B4SET.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 September 2011

  Author:

    John Burkardt

  Reference:

    Donald Kreher, Douglas Simpson,
    Combinatorial Algorithms,
    CRC Press, 1998,
    ISBN: 0-8493-3988-X,
    LC: QA164.K73.

  Parameters:

    Input, int N, the order of the master set.

    Input, int A, the item.

    Input, int T, a set.

    Output, int B4SET_INDEX, the index of the item in the set,
    or -1 if the item is not in the set.
*/
{
  int i;
  int value;

  if ( a < 1 || n < a )
  {
    value = -1;
  }
  else
  {
    value = 0;
    for ( i = 1; i <= a; i++ )
    {
      if ( i4_btest ( t, i - 1 ) )
      {
        value = value + 1;
      }
    }
  }
  return value;
}
/******************************************************************************/

int b4set_insert ( int n, int a, int t )

/******************************************************************************/
/*
  Purpose:

    B4SET_INSERT inserts an item into a B4SET.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 September 2011

  Author:

    John Burkardt

  Reference:

    Donald Kreher, Douglas Simpson,
    Combinatorial Algorithms,
    CRC Press, 1998,
    ISBN: 0-8493-3988-X,
    LC: QA164.K73.

  Parameters:

    Input, int N, the order of the master set.

    Input, int A, the item.
    1 <= A <= N.

    Input, int T, a set.

    Output, B4SET_INSERT, the modified set.
*/
{
  int value;

  if ( a < 1 || n < a )
  {
    fprintf ( stderr,  "\n" );
    fprintf ( stderr,  "B4SET_INSERT - Fatal error!\n" );
    fprintf ( stderr,  "  1 <= A <= N fails.\n" );
    exit ( 1 );
  }

  value = i4_bset ( t, a - 1 );

  return value;
}
/******************************************************************************/

int b4set_intersect ( int n, int a, int b )

/******************************************************************************/
/*
  Purpose:

    B4SET_INTERSECT computes the intersection of two B4SET's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 September 2011

  Author:

    John Burkardt

  Reference:

    Donald Kreher, Douglas Simpson,
    Combinatorial Algorithms,
    CRC Press, 1998,
    ISBN: 0-8493-3988-X,
    LC: QA164.K73.

  Parameters:

    Input, int N, the order of the master set.

    Input, int A, B, two sets.

    Output, int B4SET_INTERSECT, the intersection of A and B.
*/
{
  int i;
  int value;

  value = 0;

  for ( i = 0; i < n; i++ )
  {
    if ( i4_btest ( a, i ) && i4_btest ( b, i ) )
    {
      value = i4_bset ( value, i );
    }
  }
  return value;
}
/******************************************************************************/

bool b4set_is_empty ( int n, int t )

/******************************************************************************/
/*
  Purpose:

    B4SET_IS_EMPTY determines if a B4SET is empty.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 August 2011

  Author:

    John Burkardt

  Reference:

    Donald Kreher, Douglas Simpson,
    Combinatorial Algorithms,
    CRC Press, 1998,
    ISBN: 0-8493-3988-X,
    LC: QA164.K73.

  Parameters:

    Input, int N, the order of the master set.

    Input, int T, a set.

    Output, bool B4SET_IS_EMPTY is TRUE if T is empty.
*/
{
  bool value;

  value = ( t == 0 );

  return value;
}
/******************************************************************************/

bool b4set_is_equal ( int n, int t1, int t2 )

/******************************************************************************/
/*
  Purpose:

    B4SET_IS_EQUAL determines if two B4SET's are equal.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 September 2011

  Author:

    John Burkardt

  Reference:

    Donald Kreher, Douglas Simpson,
    Combinatorial Algorithms,
    CRC Press, 1998,
    ISBN: 0-8493-3988-X,
    LC: QA164.K73.

  Parameters:

    Input, int N, the order of the master set.

    Input, int T1, T2, two sets.

    Output, bool B4SET_IS_EQUAL, is TRUE if T1 equals T2.
*/
{
  bool value;

  value = ( t1 == t2 );

  return value;
}
/******************************************************************************/

bool b4set_is_member ( int n, int a, int t )

/******************************************************************************/
/*
  Purpose:

    B4SET_IS_MEMBER determines if an item is a member of a B4SET.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 September 2011

  Author:

    John Burkardt

  Reference:

    Donald Kreher, Douglas Simpson,
    Combinatorial Algorithms,
    CRC Press, 1998,
    ISBN: 0-8493-3988-X,
    LC: QA164.K73.

  Parameters:

    Input, int N, the order of the master set.

    Input, int A, an item.

    Input, int T, a set.

    Output, bool B4SET_IS_MEMBER, is TRUE if A is an element of T.
*/
{
  bool value;

  if ( 1 <= a && a <= n )
  {
    value = i4_btest ( t, a - 1 );
  }
  else
  {
    value = false;
  }
  return value;
}
/******************************************************************************/

bool b4set_is_subset ( int n, int t1, int t2 )

/******************************************************************************/
/*
  Purpose:

    B4SET_IS_SUBSET determines if one B4SET is a subset of another.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 September 2011

  Author:

    John Burkardt

  Reference:

    Donald Kreher, Douglas Simpson,
    Combinatorial Algorithms,
    CRC Press, 1998,
    ISBN: 0-8493-3988-X,
    LC: QA164.K73.

  Parameters:

    Input, int N, the order of the master set.

    Input, int T1, T2, two sets.

    Output, bool B4SET_IS_SUBSET, is TRUE if T1 is a subset of T2.
*/
{
  int i;
  int value;

  value = true;
 
  for ( i = 0; i < n; i++ )
  {
    if ( i4_btest ( t1, i ) && ! i4_btest ( t2, i ) )
    {
      value = false;
      return value;
    }
  }
  return value;
}
/******************************************************************************/

int b4set_lex_rank ( int n, int t )

/******************************************************************************/
/*
  Purpose:

    B4SET_LEX_RANK computes the lexicographic rank of a B4SET.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 September 2011

  Author:

    John Burkardt

  Reference:

    Donald Kreher, Douglas Simpson,
    Combinatorial Algorithms,
    CRC Press, 1998,
    ISBN: 0-8493-3988-X,
    LC: QA164.K73.

  Parameters:

    Input, int N, the order of the master set.

    Input, int T, the set.

    Output, int B4SET_LEX_RANK, the rank of the set.
*/
{
  int i;
  int rank;

  rank = 0;
  for ( i = 0; i < n; i++ )
  {
    if ( i4_btest ( t, i ) )
    {
      rank = rank + i4_power ( 2, n - i - 1 );
    }
  }
  return rank;
}
/******************************************************************************/

void b4set_lex_successor ( int n, int *t, int *rank )

/******************************************************************************/
/*
  Purpose:

    B4SET_LEX_SUCCESSOR computes the lexicographic successor of a B4SET.

  Discussion:

    In the original code, there is a last element with no successor.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 September 2011

  Author:

    John Burkardt

  Reference:

    Donald Kreher, Douglas Simpson,
    Combinatorial Algorithms,
    CRC Press, 1998,
    ISBN: 0-8493-3988-X,
    LC: QA164.K73.

  Parameters:

    Input, int N, the order of the master set.

    Input/output, int *T, describes a set.
    On input, T describes a set.
    On output, T describes the next set in the ordering.
    If the input T was the last in the ordering, then the output T
    will be the first.

    Input/output, int *RANK, the rank.
    If RANK = -1 on input, then the routine understands that this is
    the first call, and that the user wishes the routine to supply
    the first element in the ordering, which has RANK = 0.
    In general, the input value of RANK is increased by 1 for output,
    unless the very last element of the ordering was input, in which
    case the output value of RANK is 0.
*/
{
  int i;
/*
  Return the first element.
*/
  if ( *rank == -1 )
  {
    *t = 0;
    *rank = 0;
    return;
  }

  for ( i = n - 1; 0 <= i; i-- )
  {
    if ( ! i4_btest ( *t, i ) )
    {
      *t = i4_bset ( *t, i );
      *rank = *rank + 1;
      return;
    }
    else
    {
      *t = i4_bclr ( *t, i );
    }
  }
  *rank = 0;
  return;
}
/******************************************************************************/

int b4set_lex_unrank ( int rank, int n )

/******************************************************************************/
/*
  Purpose:

    B4SET_LEX_UNRANK computes the B4SET of given lexicographic rank.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 September 2011

  Author:

    John Burkardt

  Reference:

    Donald Kreher, Douglas Simpson,
    Combinatorial Algorithms,
    CRC Press, 1998,
    ISBN: 0-8493-3988-X,
    LC: QA164.K73.

  Parameters:

    Input, int RANK, the rank of the set.

    Input, int N, the order of the master set.

    Output, int B4SET_LEX_UNRANK, the set of the given rank.
*/
{
  int i;
  int nsub;
  int rank_copy;
  int set_num;
  int t;
/*
  Check.
*/
  if ( n < 1 )
  {
    fprintf ( stderr,  "\n" );
    fprintf ( stderr,  "B4SET_LEX_UNRANK - Fatal error!\n" );
    fprintf ( stderr,  "  Input N is illegal.\n" );
    exit ( 1 );
  }

  set_num = b4set_enum ( n );

  if ( rank < 0 || set_num < rank )
  {
    fprintf ( stderr,  "\n" );
    fprintf ( stderr,  "B4SET_LEX_UNRANK - Fatal error!\n" );
    fprintf ( stderr,  "  The input rank is illegal.\n" );
    exit ( 1 );
  }

  rank_copy = rank;

  t = 0;

  for ( i = n - 1; 0 <= i; i-- )
  {
    if ( ( rank_copy % 2 ) == 1 )
    {
      t = i4_bset ( t, i );
    }
    else
    {
      t = i4_bclr ( t, i );
    }
    rank_copy = rank_copy / 2;
  }
  return t;
}
/******************************************************************************/

int b4set_random ( int n, int *seed )

/******************************************************************************/
/*
  Purpose:

    B4SET_RANDOM sets a rondom B4SET.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 September 2011

  Author:

    John Burkardt

  Reference:

    Donald Kreher, Douglas Simpson,
    Combinatorial Algorithms,
    CRC Press, 1998,
    ISBN: 0-8493-3988-X,
    LC: QA164.K73.

  Parameters:

    Input, int N, the order of the master set.

    Input/output, int *SEED, a seed for the random 
    number generator.

    Output, int A, the random B4SET.
*/
{
  int a;
  bool *a_log;

  a_log = lset_random ( n, seed );

  a = lset_to_b4set ( n, a_log );

  free ( a_log );

  return a;
}
/******************************************************************************/

bool *b4set_to_lset ( int n, int t )

/******************************************************************************/
/*
  Purpose:

    B4SET_TO_LSET converts a B4SET to an LSET.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 September 2011

  Author:

    John Burkardt

  Reference:

    Donald Kreher, Douglas Simpson,
    Combinatorial Algorithms,
    CRC Press, 1998,
    ISBN: 0-8493-3988-X,
    LC: QA164.K73.

  Parameters:

    Input, int N, the order of the master set.

    Input, int T, the set.

    Input, bool B4SET_TO_LSET[N], the LSET version of the set.
*/
{
  bool *a;
  int i;

  a = ( bool * ) malloc ( n * sizeof ( bool ) );

  for ( i = 0; i < n; i++ )
  {
    a[i] = i4_btest ( t, i );
  }

  return a;
}
/******************************************************************************/

void b4set_transpose_print ( int n, int t, char *title )

/******************************************************************************/
/*
  Purpose:

    B4SET_TRANSPOSE_PRINT prints a B4SET "transposed".

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 September 2011

  Author:

    John Burkardt

  Reference:

    Donald Kreher, Douglas Simpson,
    Combinatorial Algorithms,
    CRC Press, 1998,
    ISBN: 0-8493-3988-X,
    LC: QA164.K73.

  Parameters:

    Input, int N, the order of the master set.

    Input, int T, the set.

    Input, char *TITLE, a title.
*/
{
  int i;
  int s;

  printf ( "\n" );
  printf ( "%s\n", title );

  if ( t == 0 )
  {
    printf ( "  (Empty set)\n" );
  }
  else
  {
    s = 0;
    for ( i = 0; i < n; i++ )
    {
      if ( i4_btest ( t, i ) )
      {
        printf ( "  %2d", i + 1 );
        s = s + 4;
      }
      if ( 76 < s || ( 0 < s && i == n - 1 ) )
      {
        printf ( "\n" );
        s = 0;
      }
    }
  }
  return;
}
/******************************************************************************/

int b4set_union ( int n, int a, int b )

/******************************************************************************/
/*
  Purpose:

    B4SET_UNION computes the union of two B4SET's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 September 2011

  Author:

    John Burkardt

  Reference:

    Donald Kreher, Douglas Simpson,
    Combinatorial Algorithms,
    CRC Press, 1998,
    ISBN: 0-8493-3988-X,
    LC: QA164.K73.

  Parameters:

    Input, int N, the order of the master set.

    Input, int A, B, two sets.

    Output, int B4SET_UNION, the union of A and B.
*/
{
  int c;
  int i;

  c = 0;

  for ( i = 0; i < n; i++ )
  {
    if ( i4_btest ( a, i ) || i4_btest ( b, i ) )
    {
      c = i4_bset ( c, i );
    }
  }

  return c;
}
/******************************************************************************/

int b4set_weight ( int n, int t )

/******************************************************************************/
/*
  Purpose:

    B4SET_WEIGHT computes the Hamming weight of a B4SET.

  Discussion:

    The Hamming weight is simply the number of elements in the set.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 September 2011

  Author:

    John Burkardt

  Reference:

    Donald Kreher, Douglas Simpson,
    Combinatorial Algorithms,
    CRC Press, 1998,
    ISBN: 0-8493-3988-X,
    LC: QA164.K73.

  Parameters:

    Input, int N, the order of the master set..

    Input, int T, the set.

    Output, int B4SET_WEIGHT, the Hamming weight of the set T.
*/
{
  int i;
  int weight;

  weight = 0;

  for ( i = 0; i < n; i++ )
  {
    if ( i4_btest ( t, i ) )
    {
      weight = weight + 1;
    }
  }

  return weight;
}
/******************************************************************************/

int b4set_xor ( int n, int a, int b )

/******************************************************************************/
/*
  Purpose:

    B4SET_XOR computes the symmetric difference of two B4SET's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 September 2011

  Author:

    John Burkardt

  Reference:

    Donald Kreher, Douglas Simpson,
    Combinatorial Algorithms,
    CRC Press, 1998,
    ISBN: 0-8493-3988-X,
    LC: QA164.K73.

  Parameters:

    Input, int N, the order of the master set.

    Input, int A, B, two sets.

    Output, int B4SET_XOR, the symmetric difference of A and B.
*/
{
  int c;
  int i;

  c = 0;

  for ( i = 0; i < n; i++ )
  {
    if ( (   i4_btest ( a, i ) && ! i4_btest ( b, i ) ) ||
         ( ! i4_btest ( a, i ) &&   i4_btest ( b, i ) ) )
    {
      c = i4_bset ( c, i );
    }
  }
  return c;
}

/******************************************************************************/

char digit_to_ch ( int i )

/******************************************************************************/
/*
  Purpose:

    DIGIT_TO_CH returns the base 10 digit character corresponding to a digit.

  Example:

     I     C
   -----  ---
     0    '0'
     1    '1'
   ...    ...
     9    '9'
    10    '*'
   -83    '*'

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 June 2003

  Author:

    John Burkardt

  Parameters:

    Input, int I, the digit, which should be between 0 and 9.

    Output, char DIGIT_TO_CH, the appropriate character '0'
    through '9' or '*'.
*/
{
  char c;

  if ( 0 <= i && i <= 9 )
  {
    c = '0' + i;
  }
  else
  {
    c = '*';
  }

  return c;
}
/******************************************************************************/

int i4_bclr ( int i, int p )

/******************************************************************************/
/*
  Purpose:

    I4_BCLR clears a bit of an I4.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 September 2011

  Author:

    John Burkardt

  Reference:

    Donald Kreher, Douglas Simpson,
    Combinatorial Algorithms,
    CRC Press, 1998,
    ISBN: 0-8493-3988-X,
    LC: QA164.K73.

  Parameters:

    Input, int I, the I4 whose bit is to be cleared.

    Input, int P, the position of the bit.
    0 <= P < 32.

    Output, int I4_BCLR, the modified value.
*/
{
  int value;

  value = i & ~( 1 << p );

  return value;
}
/******************************************************************************/

int i4_bset ( int i, int p )

/******************************************************************************/
/*
  Purpose:

    I4_BSET sets a bit of an I4.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 September 2011

  Author:

    John Burkardt

  Reference:

    Donald Kreher, Douglas Simpson,
    Combinatorial Algorithms,
    CRC Press, 1998,
    ISBN: 0-8493-3988-X,
    LC: QA164.K73.

  Parameters:

    Input, int I, the I4 whose bit is to be set.

    Input, int P, the position of the bit.
    0 <= P < 32.

    Output, int I4_BSET, the modified value.
*/
{
  int value;

  value = i | ( 1 << p );

  return value;
}
/******************************************************************************/

bool i4_btest ( int i, int p )

/******************************************************************************/
/*
  Purpose:

    I4_BTEST returns a bit of an I4.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 September 2011

  Author:

    John Burkardt

  Reference:

    Donald Kreher, Douglas Simpson,
    Combinatorial Algorithms,
    CRC Press, 1998,
    ISBN: 0-8493-3988-X,
    LC: QA164.K73.

  Parameters:

    Input, int I, the I4 whose bit is to be returned.

    Input, int P, the position of the bit.
    0 <= P < 32.

    Output, bool I4_BTEST, the value of the bit.
*/
{
  bool value;

  value = i & ( 1 << p );

  return value;
}
/******************************************************************************/

int i4_log_10 ( int i )

/******************************************************************************/
/*
  Purpose:

    I4_LOG_10 returns the integer part of the logarithm base 10 of an I4.

  Example:

        I  I4_LOG_10
    -----  --------
        0    0
        1    0
        2    0
        9    0
       10    1
       11    1
       99    1
      100    2
      101    2
      999    2
     1000    3
     1001    3
     9999    3
    10000    4

  Discussion:

    I4_LOG_10 ( I ) + 1 is the number of decimal digits in I.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    04 January 2004

  Author:

    John Burkardt

  Parameters:

    Input, int I, the number whose logarithm base 10 is desired.

    Output, int I4_LOG_10, the integer part of the logarithm base 10 of
    the absolute value of X.
*/
{
  int i_abs;
  int ten_pow;
  int value;

  if ( i == 0 )
  {
    value = 0;
  }
  else
  {
    value = 0;
    ten_pow = 10;

    i_abs = abs ( i );

    while ( ten_pow <= i_abs )
    {
      value = value + 1;
      ten_pow = ten_pow * 10;
    }

  }

  return value;
}
/******************************************************************************/

int i4_max ( int i1, int i2 )

/******************************************************************************/
/*
  Purpose:

    I4_MAX returns the maximum of two I4's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 October 1998

  Author:

    John Burkardt

  Parameters:

    Input, int I1, I2, are two integers to be compared.

    Output, int I4_MAX, the larger of I1 and I2.
*/
{
  int value;

  if ( i2 < i1 )
  {
    value = i1;
  }
  else
  {
    value = i2;
  }
  return value;
}
/******************************************************************************/

int i4_min ( int i1, int i2 )

/******************************************************************************/
/*
  Purpose:

    I4_MIN returns the minimum of two I4's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 October 1998

  Author:

    John Burkardt

  Parameters:

    Input, int I1, I2, two integers to be compared.

    Output, int I4_MIN, the smaller of I1 and I2.
*/
{
  int value;

  if ( i1 < i2 )
  {
    value = i1;
  }
  else
  {
    value = i2;
  }
  return value;
}
/******************************************************************************/

int i4_power ( int i, int j )

/******************************************************************************/
/*
  Purpose:

    I4_POWER returns the value of I^J.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 April 2004

  Author:

    John Burkardt

  Parameters:

    Input, int I, J, the base and the power.  J should be nonnegative.

    Output, int I4_POWER, the value of I^J.
*/
{
  int k;
  int value;

  if ( j < 0 )
  {
    if ( i == 1 )
    {
      value = 1;
    }
    else if ( i == 0 )
    {
      fprintf ( stderr,  "\n" );
      fprintf ( stderr,  "I4_POWER - Fatal error!\n" );
      fprintf ( stderr,  "  I^J requested, with I = 0 and J negative.\n" );
      exit ( 1 );
    }
    else
    {
      value = 0;
    }
  }
  else if ( j == 0 )
  {
    if ( i == 0 )
    {
      fprintf ( stderr,  "\n" );
      fprintf ( stderr,  "I4_POWER - Fatal error!\n" );
      fprintf ( stderr,  "  I^J requested, with I = 0 and J = 0.\n" );
      exit ( 1 );
    }
    else
    {
      value = 1;
    }
  }
  else if ( j == 1 )
  {
    value = i;
  }
  else
  {
    value = 1;
    for ( k = 1; k <= j; k++ )
    {
      value = value * i;
    }
  }
  return value;
}
/******************************************************************************/

char *i4_to_s ( int i )

/******************************************************************************/
/*
  Purpose:

    I4_TO_S converts an I4 to a string.

  Example:

    INTVAL  S

         1  1
        -1  -1
         0  0
      1952  1952
    123456  123456
   1234567  1234567

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 May 2011

  Author:

    John Burkardt

  Parameters:

    Input, int I, an integer to be converted.

    Output, char *I4_TO_S, the representation of the integer.
*/
{
  int digit;
  int j;
  int length;
  int ten_power;
  char *s;
  static double ten = 10.0;

  length = i4_log_10 ( i );

  ten_power = ( int ) ( pow ( ten, length ) );

  if ( i < 0 )
  {
    length = length + 1;
  }
/*
  Add one position for the trailing null.
*/
  length = length + 1;

  s = ( char * ) malloc ( length * sizeof ( char ) );

  if ( i == 0 )
  {
    s[0] = '0';
    s[1] = '\0';
    return s;
  }
/*
  Now take care of the sign.
*/
  j = 0;
  if ( i < 0 )
  {
    s[j] = '-';
    j = j + 1;
    i = abs ( i );
  }
/*
  Find the leading digit of I, strip it off, and stick it into the string.
*/
  while ( 0 < ten_power )
  {
    digit = i / ten_power;
    s[j] = digit_to_ch ( digit );
    j = j + 1;
    i = i - digit * ten_power;
    ten_power = ten_power / 10;
  }
/*
  Tack on the trailing NULL.
*/
  s[j] = '\0';
  j = j + 1;

  return s;
}
/******************************************************************************/

int i4_width ( int i )

/******************************************************************************/
/*
  Purpose:

    I4_WIDTH returns the "width" of an I4.

  Example:

        I  I4_WIDTH
    -----  -------
    -1234    5
     -123    4
      -12    3
       -1    2
        0    1
        1    1
       12    2
      123    3
     1234    4
    12345    5

  Discussion:

    The width of an integer is the number of characters necessary to print it.

    The width of an integer can be useful when setting the appropriate output
    format for a vector or array of values.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 September 2005

  Author:

    John Burkardt

  Parameters:

    Input, int I, the number whose width is desired.

    Output, int I4_WIDTH, the number of characters necessary to represent
    the integer in base 10, including a negative sign if necessary.
*/
{
  int width;

  if ( 0 <= i )
  {
    width = i4_log_10 ( i ) + 1;
  }
  else
  {
    width = i4_log_10 ( i ) + 2;
  }

  return width;
}
/******************************************************************************/

int i4_xor ( int i, int p )

/******************************************************************************/
/*
  Purpose:

    I4_XOR xor's a bit of an I4.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 September 2011

  Author:

    John Burkardt

  Reference:

    Donald Kreher, Douglas Simpson,
    Combinatorial Algorithms,
    CRC Press, 1998,
    ISBN: 0-8493-3988-X,
    LC: QA164.K73.

  Parameters:

    Input, int I, the I4 whose bit is to be xor'ed.

    Input, int P, the position of the bit.
    0 <= P < 32.

    Output, int I4_BCLR, the modified value.
*/
{
  int value;

  value = i ^ ( 1 << p );

  return value;
}
/******************************************************************************/

int i4vec_to_b4set ( int n_num, int a_num[], int n )

/******************************************************************************/
/*
  Purpose:

    I4VEC_TO_B4SET converts an I4VEC to a B4SET.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 September 2011

  Author:

    John Burkardt

  Reference:

    Donald Kreher, Douglas Simpson,
    Combinatorial Algorithms,
    CRC Press, 1998,
    ISBN: 0-8493-3988-X,
    LC: QA164.K73.

  Parameters:

    Input, int N_NUM, the number of numeric entries.

    Input, int A_NUM[N_NUM], the numeric vector.
    Entries of A_NUM should be between 1 and 32.

    Input, int N, the order of the master set.
    N <= 32.

    Output, int I4VEC_TO_B4SET, the corresponding B4SET.
*/
{
  int a;
  int i;
  int pos;
  int pos_max;

  a = 0;
  pos_max = i4_min ( 32, n );

  for ( i = 0; i < n_num; i++ )
  {
    pos = a_num[i];
    if ( 1 <= pos && pos <= pos_max )
    {
      a = i4_bset ( a, pos - 1 );
    }
  }

  return a;
}
/******************************************************************************/

bool *i4vec_to_lset ( int n_num, int a_num[], int n )

/******************************************************************************/
/*
  Purpose:

    I4VEC_TO_LSET converts an I4VEC to an LSET.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 August 2011

  Author:

    John Burkardt

  Reference:

    Donald Kreher, Douglas Simpson,
    Combinatorial Algorithms,
    CRC Press, 1998,
    ISBN: 0-8493-3988-X,
    LC: QA164.K73.

  Parameters:

    Input, int N_NUM, the number of numeric entries.

    Input, int A_NUM[N_NUM], the numeric vector.

    Input, int N, the order of the master set.

    Output, bool I4VEC_TO_LSET[N], the corresponding LSET.
*/
{
  bool *a;
  int i;

  a = ( bool * ) malloc ( n * sizeof ( bool ) );

  for ( i = 0; i < n; i++ )
  {
    a[i] = false;
  }

  for ( i = 0; i < n_num; i++ )
  {
    lset_insert ( n, a_num[i], a );
  }

  return a;
}
/******************************************************************************/

int *i4vec_uniform_ab_new ( int n, int a, int b, int *seed )

/******************************************************************************/
/*
  Purpose:

    I4VEC_UNIFORM_AB_NEW returns a scaled pseudorandom I4VEC.

  Discussion:

    The pseudorandom numbers should be uniformly distributed
    between A and B.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    24 May 2012

  Author:

    John Burkardt

  Reference:

    Paul Bratley, Bennett Fox, Linus Schrage,
    A Guide to Simulation,
    Second Edition,
    Springer, 1987,
    ISBN: 0387964673,
    LC: QA76.9.C65.B73.

    Bennett Fox,
    Algorithm 647:
    Implementation and Relative Efficiency of Quasirandom
    Sequence Generators,
    ACM Transactions on Mathematical Software,
    Volume 12, Number 4, December 1986, pages 362-376.

    Pierre L'Ecuyer,
    Random Number Generation,
    in Handbook of Simulation,
    edited by Jerry Banks,
    Wiley, 1998,
    ISBN: 0471134031,
    LC: T57.62.H37.

    Peter Lewis, Allen Goodman, James Miller,
    A Pseudo-Random Number Generator for the System/360,
    IBM Systems Journal,
    Volume 8, Number 2, 1969, pages 136-143.

  Parameters:

    Input, integer N, the dimension of the vector.

    Input, int A, B, the limits of the interval.

    Input/output, int *SEED, the "seed" value, which should NOT be 0.
    On output, SEED has been updated.

    Output, int I4VEC_UNIFORM_AB_NEW[N], a vector of random values 
    between A and B.
*/
{
  int c;
  int i;
  int i4_huge = 2147483647;
  int k;
  float r;
  int value;
  int *x;

  if ( *seed == 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "I4VEC_UNIFORM_AB_NEW - Fatal error!\n" );
    fprintf ( stderr, "  Input value of SEED = 0.\n" );
    exit ( 1 );
  }
/*
  Guaranteee A <= B.
*/
  if ( b < a )
  {
    c = a;
    a = b;
    b = c;
  }

  x = ( int * ) malloc ( n * sizeof ( int ) );

  for ( i = 0; i < n; i++ )
  {
    k = *seed / 127773;

    *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

    if ( *seed < 0 )
    {
      *seed = *seed + i4_huge;
    }

    r = ( float ) ( *seed ) * 4.656612875E-10;
/*
  Scale R to lie between A-0.5 and B+0.5.
*/
    r = ( 1.0 - r ) * ( ( float ) a - 0.5 ) 
      +         r   * ( ( float ) b + 0.5 );
/*
  Use rounding to convert R to an integer between A and B.
*/
    value = round ( r );
/*
  Guarantee A <= VALUE <= B.
*/
    if ( value < a )
    {
      value = a;
    }
    if ( b < value )
    {
      value = b;
    }

    x[i] = value;
  }

  return x;
}
/******************************************************************************/

int lset_colex_rank ( int n, bool t[] )

/******************************************************************************/
/*
  Purpose:

    LSET_COLEX_RANK computes the colexicographic rank of an LSET.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 September 2011

  Author:

    John Burkardt

  Reference:

    Donald Kreher, Douglas Simpson,
    Combinatorial Algorithms,
    CRC Press, 1998,
    ISBN: 0-8493-3988-X,
    LC: QA164.K73.

  Parameters:

    Input, int N, the order of the master set.

    Input, bool T[N], the set.

    Output, int LSET_COLEX_RANK, the rank of the set.
*/
{
  int i;
  int rank;

  rank = 0;

  for ( i = 0; i < n; i++ )
  {
    if ( t[i] )
    {
      rank = rank + i4_power ( 2, i );
    }
  }
  return rank;
}
/******************************************************************************/

void lset_colex_successor ( int n, bool t[], int *rank )

/******************************************************************************/
/*
  Purpose:

    LSET_COLEX_SUCCESSOR computes the colexicographic successor of an LSET.

  Discussion:

    In the original code, there is a last element with no successor.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 September 2011

  Author:

    John Burkardt

  Reference:

    Donald Kreher, Douglas Simpson,
    Combinatorial Algorithms,
    CRC Press, 1998,
    ISBN: 0-8493-3988-X,
    LC: QA164.K73.

  Parameters:

    Input, int N, the order of the master set.

    Input/output, bool T[N], describes a set.  
    On input, T describes a set.
    On output, T describes the next set in the ordering.
    If the input T was the last in the ordering, then the output T
    will be the first.

    Input/output, int *RANK, the rank.
    If RANK = -1 on input, then the routine understands that this is
    the first call, and that the user wishes the routine to supply
    the first element in the ordering, which has RANK = 0.
    In general, the input value of RANK is increased by 1 for output,
    unless the very last element of the ordering was input, in which
    case the output value of RANK is 0.
*/
{
  int i;
/*
  Return the first element.
*/
  if ( *rank == -1 )
  {
    for ( i = 0; i < n; i++ )
    {
      t[i] = false;
    }
    *rank = 0;
    return;
  }

  for ( i = 0; i < n; i++ )
  {
    if ( ! t[i] )
    {
      t[i] = true;
      *rank = *rank + 1;
      return;
    }
    else
    {
      t[i] = false;
    }
  }
  *rank = 0;

  return;
}
/******************************************************************************/

bool *lset_colex_unrank ( int rank, int n )

/******************************************************************************/
/*
  Purpose:

    LSET_COLEX_UNRANK computes the LSET of given colexicographic rank.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 September 2011

  Author:

    John Burkardt

  Reference:

    Donald Kreher, Douglas Simpson,
    Combinatorial Algorithms,
    CRC Press, 1998,
    ISBN: 0-8493-3988-X,
    LC: QA164.K73.

  Parameters:

    Input, int RANK, the rank of the set.

    Input, int N, the order of the master set.

    Output, bool LSET_COLEX_UNRANK[N], the set of the given rank.
*/
{
  int i;
  int rank_copy;
  int sub_num;
  bool *t;
/*
  Check.
*/
  if ( n < 1 )
  {
    fprintf ( stderr,  "\n" );
    fprintf ( stderr,  "LSET_COLEX_UNRANK - Fatal error!\n" );
    fprintf ( stderr,  "  Input N is illegal.\n" );
    exit ( 1 );
  }

  sub_num = lset_enum ( n );

  if ( rank < 0 || sub_num < rank )
  {
    fprintf ( stderr,  "\n" );
    fprintf ( stderr,  "LSET_COLEX_UNRANK - Fatal error!\n" );
    fprintf ( stderr,  "  The input rank is illegal.\n" );
    exit ( 1 );
  }

  rank_copy = rank;

  t = ( bool * ) malloc ( n * sizeof ( bool ) );

  for ( i = 0; i < n; i++ )
  {
    t[i] = ( ( rank_copy % 2 ) == 1 );
    rank_copy = rank_copy / 2;
  }
  return t;
}
/******************************************************************************/

bool *lset_complement ( int n, bool a[] )

/******************************************************************************/
/*
  Purpose:

    LSET_COMPLEMENT computes the complement of an LSET.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 September 2011

  Author:

    John Burkardt

  Reference:

    Donald Kreher, Douglas Simpson,
    Combinatorial Algorithms,
    CRC Press, 1998,
    ISBN: 0-8493-3988-X,
    LC: QA164.K73.

  Parameters:

    Input, int N, the order of the master set.

    Input, logical A[N], the set.

    Output, logical B[N], the complement of A.
*/
{
  bool *b;
  int i;
  
  b = ( bool * ) malloc ( n * sizeof ( bool ) );

  for ( i = 0; i < n; i++ )
  {
    b[i] = ! a[i];
  }
  return b;
}
/******************************************************************************/

bool *lset_complement_relative ( int n, bool a[], bool b[] )

/******************************************************************************/
/*
  Purpose:

    LSET_COMPLEMENT_RELATIVE computes the relative complement of an LSET.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 September 2011

  Author:

    John Burkardt

  Reference:

    Donald Kreher, Douglas Simpson,
    Combinatorial Algorithms,
    CRC Press, 1998,
    ISBN: 0-8493-3988-X,
    LC: QA164.K73.

  Parameters:

    Input, int N, the order of the master set.

    Input, bool A[N], the set.

    Input, bool B[N], the set with respect to which the complement is taken.

    Output, bool LSET_COMPLEMENT_RELATIVE[N], the complement of A with 
    respect to B.
*/
{
  bool *c;
  int i;

  c = ( bool * ) malloc ( n * sizeof ( bool ) );

  for ( i = 0; i < n; i++ )
  {
    c[i] = ( a[i] && ( ! b[i] ) );
  }
  return c;
}
/******************************************************************************/

void lset_delete ( int n, int a, bool t[] )

/******************************************************************************/
/*
  Purpose:

    LSET_DELETE deletes an element from an LSET.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 September 2011

  Author:

    John Burkardt

  Reference:

    Donald Kreher, Douglas Simpson,
    Combinatorial Algorithms,
    CRC Press, 1998,
    ISBN: 0-8493-3988-X,
    LC: QA164.K73.

  Parameters:

    Input, int N, the order of the master set.

    Input, int A, an item, between 1 and N.

    Input/output, bool T[N], a set.
    On output, T[A] = FALSE.
*/
{
  if ( a < 1 || n < a )
  {
    fprintf ( stderr,  "\n" );
    fprintf ( stderr,  "LSET_DELETE - Fatal error!\n" );
    fprintf ( stderr,  "  1 <= A <= N fails.\n" );
    exit ( 1 );
  }

  t[a-1] = false;

  return;
}
/******************************************************************************/

int lset_distance ( int n, bool t1[], int t2[] )

/******************************************************************************/
/*
  Purpose:

    LSET_DISTANCE computes the Hamming distance between two LSET's.

  Discussion:

    The sets T1 and T2 are assumed to be subsets of a set of N elements.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 September 2011

  Author:

    John Burkardt

  Reference:

    Donald Kreher, Douglas Simpson,
    Combinatorial Algorithms,
    CRC Press, 1998,
    ISBN: 0-8493-3988-X,
    LC: QA164.K73.

  Parameters:

    Input, int N, the order of the master set.

    Input, bool T1[N], T2[N], two sets.

    Output, int LSET_DISTANCE, the Hamming distance between T1 and T2,
    defined as the number of elements of the master set which are
    in either T1 or T2 but not both.
*/
{
  int dist;
  int i;

  dist = 0;

  for ( i = 0; i < n; i++ )
  {
    if ( (    t1[i]   && ( !t2[i] ) ) ||
         ( ( !t1[i] ) &&    t2[i]   )  )
    {
      dist = dist + 1;
    }
  }
  return dist;
}
/******************************************************************************/

int lset_enum ( int n )

/******************************************************************************/
/*
  Purpose:

    LSET_ENUM enumerates the LSET's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 September 2011

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of the master set.

    Output, int LSET_ENUM, the number of distinct sets.
*/
{
  int set_num;

  set_num = i4_power ( 2, n );

  return set_num;
}
/******************************************************************************/

int lset_index ( int n, int a, bool t[] )

/******************************************************************************/
/*
  Purpose:

    LSET_INDEX returns the index of an element of an LSET.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 August 2011

  Author:

    John Burkardt

  Reference:

    Donald Kreher, Douglas Simpson,
    Combinatorial Algorithms,
    CRC Press, 1998,
    ISBN: 0-8493-3988-X,
    LC: QA164.K73.

  Parameters:

    Input, int N, the order of the master set.

    Input, int A, the item.

    Input, bool T[N], a set.

    Output, int LSET_INDEX, the index of the item in the set,
    or -1 if the item is not in the set.
*/
{
  int i;
  int value;

  if ( a < 1 || n < a )
  {
    value = -1;
  }
  else
  {
    value = 0;
    for ( i = 0; i < a; i++ )
    {
      if ( t[i] )
      {
        value = value + 1;
      }
    }
  }
  return value;
}
/******************************************************************************/

void lset_insert ( int n, int a, bool t[] )

/******************************************************************************/
/*
  Purpose:

    LSET_INSERT inserts an item into an LSET.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 September 2011

  Author:

    John Burkardt

  Reference:

    Donald Kreher, Douglas Simpson,
    Combinatorial Algorithms,
    CRC Press, 1998,
    ISBN: 0-8493-3988-X,
    LC: QA164.K73.

  Parameters:

    Input, int N, the order of the master set.

    Input, int A, the item.
    1 <= A <= N.

    Input/output, bool T[N], a set.
    On output, T[A-1] = TRUE.
*/
{
  if ( a < 1 || n < a )
  {
    fprintf ( stderr,  "\n" );
    fprintf ( stderr,  "LSET_INSERT - Fatal error!\n" );
    fprintf ( stderr,  "  1 <= A <= N fails.\n" );
    exit ( 1 );
  }

  t[a-1] = true;

  return;
}
/******************************************************************************/

bool *lset_intersect ( int n, bool a[], bool b[] )

/******************************************************************************/
/*
  Purpose:

    LSET_INTERSECT computes the intersection of two LSET's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 September 2011

  Author:

    John Burkardt

  Reference:

    Donald Kreher, Douglas Simpson,
    Combinatorial Algorithms,
    CRC Press, 1998,
    ISBN: 0-8493-3988-X,
    LC: QA164.K73.

  Parameters:

    Input, int N, the order of the master set.

    Input, bool A[N], B[N], two sets.

    Output, bool LSET_INTERSECT[N], the intersection of A and B.
*/
{
  bool *c;
  int i;

  c = ( bool * ) malloc ( n * sizeof ( bool ) );

  for ( i = 0; i < n; i++ )
  {
    c[i] = ( a[i] && b[i] );
  }
  return c;
}
/******************************************************************************/

bool lset_is_empty ( int n, bool t[] )

/******************************************************************************/
/*
  Purpose:

    LSET_IS_EMPTY determines if an LSET is empty.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 September 2011

  Author:

    John Burkardt

  Reference:

    Donald Kreher, Douglas Simpson,
    Combinatorial Algorithms,
    CRC Press, 1998,
    ISBN: 0-8493-3988-X,
    LC: QA164.K73.

  Parameters:

    Input, int N, the order of the master set.

    Input, bool T[N], a set.

    Output, bool LSET_IS_EMPTY is TRUE if T is empty.
*/
{
  int i;
  bool value;

  value = false;

  for ( i = 0; i < n; i++ )
  {
    if ( t[i] ) 
    {
      return value;
    }
  }
  value = true;
  return value;
}
/******************************************************************************/

bool lset_is_equal ( int n, bool t1[], bool t2[] )

/******************************************************************************/
/*
  Purpose:

    LSET_IS_EQUAL determines if two LSET's are equal.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 September 2011

  Author:

    John Burkardt

  Reference:

    Donald Kreher, Douglas Simpson,
    Combinatorial Algorithms,
    CRC Press, 1998,
    ISBN: 0-8493-3988-X,
    LC: QA164.K73.

  Parameters:

    Input, int N, the order of the master set.

    Input, bool T1[N], T2[N], two sets.

    Output, bool LSET_IS_EQUAL, is TRUE if T1 equals T2.
*/
{
  int i;
  bool value;

  value = false;
  for ( i = 0; i < n; i++ )
  {
    if ( t1[i] != t2[i] )
    {
      return value;
    }
  }
  value = true;
  return value;
}
/******************************************************************************/

bool lset_is_member ( int n, int a, bool t[] )

/******************************************************************************/
/*
  Purpose:

    LSET_IS_MEMBER determines if an item is a member of an LSET.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 September 2011

  Author:

    John Burkardt

  Reference:

    Donald Kreher, Douglas Simpson,
    Combinatorial Algorithms,
    CRC Press, 1998,
    ISBN: 0-8493-3988-X,
    LC: QA164.K73.

  Parameters:

    Input, int N, the order of the master set.

    Input, int A, an item.

    Input, bool T[N], a set.

    Output, bool LSET_IS_MEMBER, is TRUE if A is an element of T.
*/
{
  bool value;

  if ( 1 <= a && a <= n )
  {
    value = t[a-1];
  }
  else
  {
    value = false;
  }
  return value;
}
/******************************************************************************/

bool lset_is_subset ( int n, bool t1[], bool t2[] )

/******************************************************************************/
/*
  Purpose:

    LSET_IS_SUBSET determines if one LSET is a subset of another.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 September 2011

  Author:

    John Burkardt

  Reference:

    Donald Kreher, Douglas Simpson,
    Combinatorial Algorithms,
    CRC Press, 1998,
    ISBN: 0-8493-3988-X,
    LC: QA164.K73.

  Parameters:

    Input, int N, the order of the master set.

    Input, bool T1[N], T2[N], two sets.

    Output, bool LSET_IS_SUBSET, is TRUE if T1 is a subset of T2.
*/
{
  int i;
  bool value;

  value = true;

  for ( i = 0; i < n; i++ )
  {
    if ( t1[i] && ! t2[i] )
    {
      value = false;
      return value;
    }
  }
  return value;
}
/******************************************************************************/

int lset_lex_rank ( int n, bool t[] )

/******************************************************************************/
/*
  Purpose:

    LSET_LEX_RANK computes the lexicographic rank of an LSET.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 September 2011

  Author:

    John Burkardt

  Reference:

    Donald Kreher, Douglas Simpson,
    Combinatorial Algorithms,
    CRC Press, 1998,
    ISBN: 0-8493-3988-X,
    LC: QA164.K73.

  Parameters:

    Input, int N, the order of the master set.

    Input, bool T[N], the set.

    Output, int LSET_LEX_RANK, the rank of the set.
*/
{
  int i;
  int rank;

  rank = 0;

  for ( i = 0; i < n; i++ )
  {
    if ( t[i] )
    {
      rank = rank + i4_power ( 2, n - i - 1 );
    }
  }
  return rank;
}
/******************************************************************************/

void lset_lex_successor ( int n, bool t[], int *rank )

/******************************************************************************/
/*
  Purpose:

    LSET_LEX_SUCCESSOR computes the lexicographic successor of an LSET.

  Discussion:

    In the original code, there is a last element with no successor.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 September 2011

  Author:

    John Burkardt

  Reference:

    Donald Kreher, Douglas Simpson,
    Combinatorial Algorithms,
    CRC Press, 1998,
    ISBN: 0-8493-3988-X,
    LC: QA164.K73.

  Parameters:

    Input, int N, the order of the master set.

    Input/output, bool T[N], describes a set.
    On input, T describes a set.
    On output, T describes the next set in the ordering.
    If the input T was the last in the ordering, then the output T
    will be the first.

    Input/output, int *RANK, the rank.
    If RANK = -1 on input, then the routine understands that this is
    the first call, and that the user wishes the routine to supply
    the first element in the ordering, which has RANK = 0.
    In general, the input value of RANK is increased by 1 for output,
    unless the very last element of the ordering was input, in which
    case the output value of RANK is 0.
*/
{
  int i;
/*
  Return the first element.
*/
  if ( *rank == -1 )
  {
    for ( i = 0; i < n; i++ )
    {
      t[i] = false;
    }
    *rank = 0;
    return;
  }

  for ( i = n - 1; 0 <= i; i-- )
  {
    if ( ! t[i] )
    {
      t[i] = true;
      *rank = *rank + 1;
      return;
    }
    else
    {
      t[i] = false;
    }
  }

  *rank = 0;

  return;
}
/******************************************************************************/

bool *lset_lex_unrank ( int rank, int n )

/******************************************************************************/
/*
  Purpose:

    LSET_LEX_UNRANK computes the LSET of given lexicographic rank.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 September 2011

  Author:

    John Burkardt

  Reference:

    Donald Kreher, Douglas Simpson,
    Combinatorial Algorithms,
    CRC Press, 1998,
    ISBN: 0-8493-3988-X,
    LC: QA164.K73.

  Parameters:

    Input, int RANK, the rank of the set.

    Input, int N, the order of the master set.

    Output, bool LSET_LEX_UNRANK[N], the set of the given rank.
*/
{
  int i;
  int nsub;
  int rank_copy;
  int set_num;
  bool *t;
/*
  Check.
*/
  if ( n < 1 )
  {
    fprintf ( stderr,  "\n" );
    fprintf ( stderr,  "LSET_LEX_UNRANK - Fatal error!\n" );
    fprintf ( stderr,  "  Input N is illegal.\n" );
    exit ( 1 );
  }

  set_num = lset_enum ( n );

  if ( rank < 0 || set_num < rank )
  {
    fprintf ( stderr,  "\n" );
    fprintf ( stderr,  "LSET_LEX_UNRANK - Fatal error!\n" );
    fprintf ( stderr,  "  The input rank is illegal.\n" );
    exit ( 1 );
  }

  rank_copy = rank;

  t = ( bool * ) malloc ( n * sizeof ( bool ) );

  for ( i = n - 1; 0 <= i; i-- )
  {
    if ( ( rank_copy % 2 ) == 1 )
    {
      t[i] = true;
    }
    else
    {
      t[i] = false;
    }
    rank_copy = rank_copy / 2;
  }
  return t;
}
/******************************************************************************/

bool *lset_random ( int n, int *seed )

/******************************************************************************/
/*
  Purpose:

    LSET_RANDOM sets a rondom LSET.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 September 2011

  Author:

    John Burkardt

  Reference:

    Donald Kreher, Douglas Simpson,
    Combinatorial Algorithms,
    CRC Press, 1998,
    ISBN: 0-8493-3988-X,
    LC: QA164.K73.

  Parameters:

    Input, int N, the order of the master set.

    Input/output, int *SEED, a seed for the random 
    number generator.

    Output, bool LSET_RANDOM[N].
*/
{
  bool *a;
  static int i4_huge = 2147483647;
  static int i4_huge_half = 1073741823;
  int i;
  int k;


  if ( seed == 0 )
  {
    fprintf ( stderr,  "\n" );
    fprintf ( stderr,  "LSET_RANDOM - Fatal error!\n" );
    fprintf ( stderr,  "  Input value of SEED = 0.\n" );
    exit ( 1 );
  }

  a = ( bool * ) malloc ( n * sizeof ( bool ) );

  for ( i = 0; i < n; i++ )
  {
    k = *seed / 127773;

    *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

    if ( seed < 0 )
    {
      *seed = *seed + i4_huge;
    }

    a[i] = ( i4_huge_half < *seed );
  }

  return a;

  return a;
}
/******************************************************************************/

int lset_to_b4set ( int n, bool a_log[] )

/******************************************************************************/
/*
  Purpose:

    LSET_TO_B4SET converts an I4VEC to a B4SET.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 September 2011

  Author:

    John Burkardt

  Reference:

    Donald Kreher, Douglas Simpson,
    Combinatorial Algorithms,
    CRC Press, 1998,
    ISBN: 0-8493-3988-X,
    LC: QA164.K73.

  Parameters:

    Input, int N, the order of the master set.
    N <= 32.

    Input, bool A_LOG[N], the logical representation of the set.

    Output, int A, the corresponding B4SET.
*/
{
  int a;
  int i;

  a = 0;

  for ( i = 0; i < n; i++ )
  {
    if ( a_log[i] )
    {
      a = i4_bset ( a, i );
    }
  }

  return a;
}
/******************************************************************************/

void lset_transpose_print ( int n, bool t[], char *title )

/******************************************************************************/
/*
  Purpose:

    LSET_TRANSPOSE_PRINT prints an LSET "transposed".

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 September 2011

  Author:

    John Burkardt

  Reference:

    Donald Kreher, Douglas Simpson,
    Combinatorial Algorithms,
    CRC Press, 1998,
    ISBN: 0-8493-3988-X,
    LC: QA164.K73.

  Parameters:

    Input, int N, the order of the master set.

    Input, bool T[N], the set.

    Input, char *TITLE, a title.
*/
{
  int i;
  int p;
  int s;

  printf ( "\n" );
  printf ( "%s\n", title );

  p = 0;

  if ( lset_is_empty ( n, t ) )
  {
    printf ( "  (Empty set)\n" );
  }
  else
  {
    for ( i = 0; i < n; i++ )
    {
      if ( t[i] )
      {
        printf ( "  %d", i + 1 );
        s = 2 + i4_width ( i + 1 );
        p = p + s;
        if ( 80 <= p )
        {
          printf ( "\n" );
          p = 0;
        }        
      }
    }
  }

  if ( 0 < p )
  {
    printf ( "\n" );
  }

  return;
}
/******************************************************************************/

bool *lset_union ( int n, bool a[], bool b[] )

/******************************************************************************/
/*
  Purpose:

    LSET_UNION computes the union of two LSET's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 September 2011

  Author:

    John Burkardt

  Reference:

    Donald Kreher, Douglas Simpson,
    Combinatorial Algorithms,
    CRC Press, 1998,
    ISBN: 0-8493-3988-X,
    LC: QA164.K73.

  Parameters:

    Input, int N, the order of the master set.

    Input, bool A[N], B[N], two sets.

    Output, bool LSET_UNION[N], the union of A and B.
*/
{
  bool *c;
  int i;

  c = ( bool * ) malloc ( n * sizeof ( bool ) );

  for ( i = 0; i < n; i++ )
  {
    c[i] = ( a[i] || b[i] );
  }
  return c;
}
/******************************************************************************/

int lset_weight ( int n, bool t[] )

/******************************************************************************/
/*
  Purpose:

    LSET_WEIGHT computes the Hamming weight of an LSET.

  Discussion:

    The Hamming weight is simply the number of elements in the set.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 September 2011

  Author:

    John Burkardt

  Reference:

    Donald Kreher, Douglas Simpson,
    Combinatorial Algorithms,
    CRC Press, 1998,
    ISBN: 0-8493-3988-X,
    LC: QA164.K73.

  Parameters:

    Input, int N, the order of the master set..

    Input, bool T[N], the set.

    Output, int LSET_WEIGHT, the Hamming weight of the set T.
*/
{
  int i;
  int weight;

  weight = 0;

  for ( i = 0; i < n; i++ )
  {
    if ( t[i] )
    {
      weight = weight + 1;
    }
  }
  return weight;
}
/******************************************************************************/

bool *lset_xor ( int n, bool a[], bool b[] )

/******************************************************************************/
/*
  Purpose:

    LSET_XOR computes the symmetric difference of two LSET's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 September 2011

  Author:

    John Burkardt

  Reference:

    Donald Kreher, Douglas Simpson,
    Combinatorial Algorithms,
    CRC Press, 1998,
    ISBN: 0-8493-3988-X,
    LC: QA164.K73.

  Parameters:

    Input, int N, the order of the master set.

    Input, bool A[N], B[N], two sets.

    Output, bool LSET_XOR[N], the symmetric difference of A and B.
*/
{
  bool *c;
  int i;

  c = ( bool * ) malloc ( n * sizeof ( bool ) );

  for ( i = 0; i < n; i++ )
  {
    c[i] =  (       a[i]   && ( ! b[i] ) ) ||
              ( ( ! a[i] ) &&     b[i]   );
  }
  return c;
}
/******************************************************************************/

void lvec_transpose_print ( int n, bool t[], char *title )

/******************************************************************************/
/*
  Purpose:

    LVEC_TRANSPOSE_PRINT prints an LVEC "transposed".

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 September 2011

  Author:

    John Burkardt

  Reference:

    Donald Kreher, Douglas Simpson,
    Combinatorial Algorithms,
    CRC Press, 1998,
    ISBN: 0-8493-3988-X,
    LC: QA164.K73.

  Parameters:

    Input, int N, the order of the master set.

    Input, bool T[N], the set.

    Input, char *TITLE, a title.
*/
{
  int i;
  int ihi;
  int ilo;

  fprintf ( stderr,  "\n" );
  fprintf ( stderr,  "%s\n", title );
  fprintf ( stderr,  "\n" );

  for ( ilo = 0; ilo < n; ilo = ilo + 80 )
  {
    ihi = i4_min ( ilo + 80, n );
    for ( i = ilo; i < ihi; i++ )
    {
      printf ( "%d", t[i] );
    }
    printf ( "\n" );
  }
  return;
}
/******************************************************************************/

void timestamp ( void )

/******************************************************************************/
/*
  Purpose:

    TIMESTAMP prints the current YMDHMS date as a time stamp.

  Example:

    31 May 2001 09:45:54 AM

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    24 September 2003

  Author:

    John Burkardt

  Parameters:

    None
*/
{
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct tm *tm;
  size_t len;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  len = strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  fprintf ( stdout, "%s\n", time_buffer );

  return;
# undef TIME_SIZE
}
