# include <stdlib.h>
# include <stdio.h>
# include <time.h>

# include "monomial.h"

int main ( );
void mono_between_enum_test ( );
void mono_between_next_grevlex_test ( );
void mono_between_next_grlex_test ( );
void mono_between_random_test ( );
void mono_next_grevlex_test ( );
void mono_next_grlex_test ( );
void mono_print_test ( );
void mono_rank_grlex_test ( );
void mono_total_enum_test ( );
void mono_total_next_grevlex_test ( );
void mono_total_next_grlex_test ( );
void mono_total_random_test ( );
void mono_unrank_grlex_test ( );
void mono_upto_enum_test ( );
void mono_upto_next_grevlex_test ( );
void mono_upto_next_grlex_test ( );
void mono_upto_random_test ( );
void mono_value_test ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for MONOMIAL_PRB.

  Discussion:

    MONOMIAL_PRB tests the MONOMIAL library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 February 2015

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "MONOMIAL_PRB\n" );
  printf ( "  C version.\n" );
  printf ( "  Test the MONOMIAL library.\n" );

  mono_between_enum_test ( );
  mono_between_next_grevlex_test ( );
  mono_between_next_grlex_test ( );
  mono_between_random_test ( );

  mono_next_grevlex_test ( );
  mono_next_grlex_test ( );
  mono_print_test ( );
  mono_rank_grlex_test ( );

  mono_total_enum_test ( );
  mono_total_next_grevlex_test ( );
  mono_total_next_grlex_test ( );
  mono_total_random_test ( );

  mono_unrank_grlex_test ( );

  mono_upto_enum_test ( );
  mono_upto_next_grevlex_test ( );
  mono_upto_next_grlex_test ( );
  mono_upto_random_test ( );

  mono_value_test ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "MONOMIAL_PRB\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void mono_between_enum_test ( )

/******************************************************************************/
/*
  Purpose:

    MONO_BETWEEN_ENUM_TEST tests MONO_BETWEEN_ENUM.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    18 November 2013

  Author:

    John Burkardt
*/
{
  int m;
  int n1;
  int n2;
  int v;

  printf ( "\n" );
  printf ( "MONO_BETWEEN_ENUM_TEST\n" );
  printf ( "  MONO_BETWEEN_ENUM can enumerate the number of monomials\n" );
  printf ( "  in M variables, of total degree between N1 and N2.\n" );

  m = 3;
  printf ( "\n" );
  printf ( "  Using spatial dimension M = %d\n", m );
  printf ( "\n" );
  printf ( "   N2:" );
  for ( n2 = 0; n2 <= 8; n2++ )
  {
    printf ( "  %4d", n2 );
  }
  printf ( "\n" );
  printf ( "  N1 +------------------------------------------------------\n" );
  for ( n1 = 0; n1 <= 8; n1++ )
  {
    printf ( "  %2d |", n1 );
    for ( n2 = 0; n2 <= 8; n2++ )
    {
      v = mono_between_enum ( m, n1, n2 );
      printf ( "  %4d", v );
    }
    printf ( "\n" );
  }
  return;
}
/******************************************************************************/

void mono_between_next_grevlex_test ( )

/******************************************************************************/
/*
  Purpose:

    MONO_BETWEEN_NEXT_GREVLEX_TEST tests MONO_BETWEEN_NEXT_GREVLEX.

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
  int n1;
  int n2;
  int x[3];

  printf ( "\n" );
  printf ( "MONO_BETWEEN_NEXT_GREVLEX_TEST\n" );
  printf ( "  MONO_BETWEEN_NEXT_GREVLEX can list the monomials\n" );
  printf ( "  in M variables, of total degree N between N1 and N2,\n" );
  printf ( "  in grevlex order, one at a time.\n" );
  printf ( "\n" );
  printf ( "  We start the process with (0,0,...,0,N1).\n" );
  printf ( "  The process ends with (N2,0,...,0,0)\n" );

  n1 = 2;
  n2 = 3;

  printf ( "\n" );
  printf ( "  Let M =  %d\n", m );
  printf ( "      N1 = %d\n", n1 );
  printf ( "      N2 = %d\n", n2 );
  printf ( "\n" );

  x[0] = 0;
  x[1] = 0;
  x[2] = n1; 

  i = 1;

  for ( ; ; )
  {
    printf ( "  %2d    ", i );
    for ( j = 0; j < m; j++ )
    {
      printf ( "%2d", x[j] );
    }
    printf ( "\n" );

    if ( x[0] == n2 )
    {
      break;
    }
 
    mono_between_next_grevlex ( m, n1, n2, x );
    i = i + 1;
  }

  return;
}
/******************************************************************************/

void mono_between_next_grlex_test ( )

/******************************************************************************/
/*
  Purpose:

    MONO_BETWEEN_NEXT_GRLEX_TEST tests MONO_BETWEEN_NEXT_GRLEX.

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
  int n1;
  int n2;
  int x[3];

  printf ( "\n" );
  printf ( "MONO_BETWEEN_NEXT_GRLEX_TEST\n" );
  printf ( "  MONO_BETWEEN_NEXT_GRLEX can list the monomials\n" );
  printf ( "  in M variables, of total degree N between N1 and N2,\n" );
  printf ( "  in grlex order, one at a time.\n" );
  printf ( "\n" );
  printf ( "  We start the process with (0,0,...,0,N1).\n" );
  printf ( "  The process ends with (N2,0,...,0,0)\n" );

  n1 = 2;
  n2 = 3;

  printf ( "\n" );
  printf ( "  Let M =  %d\n", m );
  printf ( "      N1 = %d\n", n1 );
  printf ( "      N2 = %d\n", n2 );
  printf ( "\n" );

  x[0] = 0;
  x[1] = 0;
  x[2] = n1; 

  i = 1;

  for ( ; ; )
  {
    printf ( "  %2d    ", i );
    for ( j = 0; j < m; j++ )
    {
      printf ( "%2d", x[j] );
    }
    printf ( "\n" );

    if ( x[0] == n2 )
    {
      break;
    }
 
    mono_between_next_grlex ( m, n1, n2, x );
    i = i + 1;
  }

  return;
}
/******************************************************************************/

void mono_between_random_test ( )

/******************************************************************************/
/*
  Purpose:

    MONO_BETWEEN_RANDOM_TEST tests MONO_BETWEEN_RANDOM.

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
  int n1;
  int n2;
  int rank;
  int seed;
  int test;
  int test_num;
  int *x;

  printf ( "\n" );
  printf ( "MONO_BETWEEN_RANDOM_TEST\n" );
  printf ( "  MONO_BETWEEN_RANDOM selects at random a monomial\n" );
  printf ( "  in M dimensions of total degree between N1 and N2.\n" );

  n1 = 2;
  n2 = 3;

  printf ( "\n" );
  printf ( "  Let M =  %d\n", m );
  printf ( "      N1 = %d\n", n1 );
  printf ( "      N2 = %d\n", n2 );
  printf ( "\n" );

  seed = 123456789;
  test_num = 5;

  for ( test = 1; test <= test_num; test++ )
  {
    x = mono_between_random ( m, n1, n2, &seed, &rank );
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

void mono_next_grevlex_test ( )

/******************************************************************************/
/*
  Purpose:

    MONO_NEXT_GREVLEX_TEST tests MONO_NEXT_GREVLEX.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 February 2015

  Author:

    John Burkardt
*/
{
  int d;
  int i;
  int k;
  int m = 4;
  int *x;

  printf ( "\n" );
  printf ( "MONO_NEXT_GREVLEX_TEST\n" );
  printf ( "  MONO_NEXT_GREVLEX computes the next monomial\n" );
  printf ( "  in M variables in grevlex order.\n" );
  printf ( "\n" );
  printf ( "  Let M =  %d\n", m );

  x = ( int * ) malloc ( m * sizeof ( int ) );
  k = 0;
  for ( i = 0; i < m; i++ )
  {
    x[i] = 0;
  }

  for ( ; ; )
  {
    d = i4vec_sum ( m, x );
    printf ( "  %2d  %2d  |", k, d );
    for ( i = 0; i < m; i++ )
    {
      printf ( "%2d", x[i] );
    }
    printf ( "\n" );

    if ( x[0] == 3 )
    {
      break;
    }
    mono_next_grevlex ( m, x );
    k = k + 1;
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

    05 February 2015

  Author:

    John Burkardt
*/
{
  int d;
  int i;
  int k;
  int m = 4;
  int *x;

  printf ( "\n" );
  printf ( "MONO_NEXT_GRLEX_TEST\n" );
  printf ( "  MONO_NEXT_GRLEX computes the next monomial\n" );
  printf ( "  in M variables in grlex order.\n" );
  printf ( "\n" );
  printf ( "  Let M =  %d\n", m );

  x = ( int * ) malloc ( m * sizeof ( int ) );
  k = 0;
  for ( i = 0; i < m; i++ )
  {
    x[i] = 0;
  }

  for ( ; ; )
  {
    d = i4vec_sum ( m, x );
    printf ( "  %2d  %2d  |", k, d );
    for ( i = 0; i < m; i++ )
    {
      printf ( "%2d", x[i] );
    }
    printf ( "\n" );

    if ( x[0] == 3 )
    {
      break;
    }
    mono_next_grlex ( m, x );
    k = k + 1;
  }

  free ( x );

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

void mono_total_enum_test ( )

/******************************************************************************/
/*
  Purpose:

    MONO_TOTAL_ENUM_TEST tests MONO_TOTAL_ENUM.

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
  printf ( "MONO_TOTAL_ENUM_TEST\n" );
  printf ( "  MONO_TOTAL_ENUM can enumerate the number of monomials\n" );
  printf ( "  in M variables, of total degree N.\n" );

  printf ( "\n" );
  printf ( "    N:" );
  for ( n = 0; n <= 8; n++ )
  {
    printf ( "  %4d", n );
  }
  printf ( "\n" );
  printf ( "   M +------------------------------------------------------\n" );
  for ( m = 1; m <= 8; m++ )
  {
    printf ( "  %2d |", m );
    for ( n = 0; n <= 8; n++ )
    {
      v = mono_total_enum ( m, n );
      printf ( "  %4d", v );
    }
    printf ( "\n" );
  }
  return;
}
/******************************************************************************/

void mono_total_next_grevlex_test ( )

/******************************************************************************/
/*
  Purpose:

    MONO_TOTAL_NEXT_GREVLEX_TEST tests MONO_TOTAL_NEXT_GREVLEX.

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
  printf ( "MONO_TOTAL_NEXT_GREVLEX_TEST\n" );
  printf ( "  MONO_TOTAL_NEXT_GREVLEX can list the monomials\n" );
  printf ( "  in M variables, of total degree N,\n" );
  printf ( "  in grevlex order, one at a time.\n" );
  printf ( "\n" );
  printf ( "  We start the process with (0,0,...,0,N).\n" );
  printf ( "  The process ends with (N,0,...,0,0)\n" );

  n = 3;

  printf ( "\n" );
  printf ( "  Let M = %d\n", m );
  printf ( "      N = %d\n", n );
  printf ( "\n" );

  x[0] = 0;
  x[1] = 0;
  x[2] = n; 

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

    mono_total_next_grevlex ( m, n, x );
    i = i + 1;
  }

  return;
}
/******************************************************************************/

void mono_total_next_grlex_test ( )

/******************************************************************************/
/*
  Purpose:

    MONO_TOTAL_NEXT_GRLEX_TEST tests MONO_TOTAL_NEXT_GRLEX.

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
  printf ( "MONO_TOTAL_NEXT_GRLEX_TEST\n" );
  printf ( "  MONO_TOTAL_NEXT_GRLEX can list the monomials\n" );
  printf ( "  in M variables, of total degree N,\n" );
  printf ( "  in grlex order, one at a time.\n" );
  printf ( "\n" );
  printf ( "  We start the process with (0,0,...,0,N).\n" );
  printf ( "  The process ends with (N,0,...,0,0)\n" );

  n = 3;

  printf ( "\n" );
  printf ( "  Let M = %d\n", m );
  printf ( "      N = %d\n", n );
  printf ( "\n" );

  x[0] = 0;
  x[1] = 0;
  x[2] = n; 

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

    mono_total_next_grlex ( m, n, x );
    i = i + 1;
  }

  return;
}
/******************************************************************************/

void mono_total_random_test ( )

/******************************************************************************/
/*
  Purpose:

    MONO_TOTAL_RANDOM_TEST tests MONO_TOTAL_RANDOM.

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
  printf ( "MONO_TOTAL_RANDOM_TEST\n" );
  printf ( "  MONO_TOTAL_RANDOM selects at random a monomial\n" );
  printf ( "  in M dimensions of total degree N.\n" );

  n = 4;

  printf ( "\n" );
  printf ( "  Let M = %d\n", m );
  printf ( "      N = %d\n", n );
  printf ( "\n" );

  seed = 123456789;
  test_num = 5;

  for ( test = 1; test <= test_num; test++ )
  {
    x = mono_total_random ( m, n, &seed, &rank );
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

void mono_upto_next_grevlex_test ( )

/******************************************************************************/
/*
  Purpose:

    MONO_UPTO_NEXT_GREVLEX_TEST tests MONO_UPTO_NEXT_GREVLEX.

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
  printf ( "MONO_UPTO_NEXT_GREVLEX_TEST\n" );
  printf ( "  MONO_UPTO_NEXT_GREVLEX can list the monomials\n" );
  printf ( "  in M variables, of total degree up to N,\n" );
  printf ( "  in grevlex order, one at a time.\n" );
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

    mono_upto_next_grevlex ( m, n, x );
    i = i + 1;

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

